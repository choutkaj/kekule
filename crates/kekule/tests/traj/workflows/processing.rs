//! Samples processed as they are read: each step alone and all of them
//! equal the in-memory functions and the streaming adapters, smoothing
//! averages a window of processed samples imaged onto the sample, and a set
//! read from a file processes like the same set in memory.

use std::collections::HashMap;
use std::panic::AssertUnwindSafe;
use std::sync::{Arc, Mutex};

use kekule::diagnostics::UnrepresentablePolicy;
use kekule::geometry::{Matrix3, PeriodicCell, PeriodicGeometry, Point3, RigidTransform, Vector3};
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{
    Distribution, Forces, Model, Order, Positions, Sample, SampleRead, SampleReadError, SampleSet,
    SampleView, Trajectory, TrajectoryFrame, Velocities,
};
use kekule::topology::{AtomSelection, Topology, TopologyAtomIndex, TopologyBuilder};
use kekule::traj::analysis::FrameSuperposer;
use kekule::traj::io::trr::{TrrLambdaPolicy, TrrWriteOptions};
use kekule::traj::io::{open_indexed_sample_set, write_trajectory, TrajectoryWriteOptions};
use kekule::traj::periodic::{self, MoleculeImager};
use kekule::traj::{
    FitReference, FrameBuffer, Processed, ProcessingError, ProcessingStep, ProcessingSteps,
    TrajectoryFormat,
};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};

use crate::support::linear_carbon_topology;

/// Seven frames of two three-atom chains crossing the faces of a 1 nm cell,
/// turning a little, with velocities, forces, and a property.
fn source() -> Trajectory {
    let mut graph = MoleculeEditor::new();
    let chain = ["C", "C", "O"].map(|symbol| {
        graph
            .add_atom(Atom::new(Element::from_symbol(symbol).unwrap()))
            .unwrap()
    });
    for pair in chain.windows(2) {
        graph.add_bond(pair[0], pair[1], BondOrder::Single).unwrap();
    }
    let molecule = graph.finish().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(molecule).unwrap();
    builder.add_instance(definition).unwrap();
    builder.add_instance(definition).unwrap();
    let topology = Arc::new(builder.build().unwrap());
    let mut trajectory = Trajectory::new(topology);
    for index in 0..7 {
        let t = index as f64;
        let x = (0.85 + 0.07 * t) % 1.0;
        let y = (0.1 + 0.9 - 0.05 * t) % 1.0;
        let turn = 0.02 * t;
        let wrap = |p: Point3| Point3::new(p.x.rem_euclid(1.0), p.y.rem_euclid(1.0), p.z);
        let points = [
            Point3::new(x, 0.3, 0.2),
            wrap(Point3::new(x + 0.14, 0.3 + turn, 0.2)),
            wrap(Point3::new(x + 0.21, 0.42 + turn, 0.25)),
            Point3::new(0.6, y, 0.7),
            wrap(Point3::new(0.6 + turn, y + 0.15, 0.7)),
            wrap(Point3::new(0.72, y + 0.24, 0.75 - turn)),
        ];
        let mut frame =
            TrajectoryFrame::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap())
                .with_cell(
                    PeriodicCell::orthorhombic(
                        Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER),
                        [true; 3],
                    )
                    .unwrap(),
                );
        frame.set_time(Some(Quantity::new(t, PICOSECOND))).unwrap();
        frame.set_step(Some(index * 10));
        frame.set_velocities(Some(Velocities::zeros(6))).unwrap();
        frame.set_forces(Some(Forces::zeros(6))).unwrap();
        frame
            .properties_mut()
            .owner_mut()
            .insert(
                PropertyKey::new("label").unwrap(),
                PropertyValue::Int(index as i64),
            )
            .unwrap();
        trajectory.push(frame).unwrap();
    }
    trajectory
}

fn set_of(trajectory: &Trajectory) -> SampleSet {
    SampleSet::from_trajectory(trajectory.clone()).unwrap()
}

fn atoms(trajectory: &Trajectory, indices: &[u32]) -> AtomSelection {
    AtomSelection::from_indices(
        trajectory.shared_topology(),
        indices.iter().map(|&index| TopologyAtomIndex::new(index)),
    )
    .unwrap()
}

/// Checks that `sample` holds everything `frame` does.
fn assert_sample(sample: &Sample, frame: &TrajectoryFrame, context: &str) {
    assert_eq!(sample.positions(), frame.positions(), "{context}");
    assert_eq!(sample.cell(), frame.cell(), "{context}");
    assert_eq!(sample.velocities(), frame.velocities(), "{context}");
    assert_eq!(sample.forces(), frame.forces(), "{context}");
    assert_eq!(sample.properties(), frame.properties(), "{context}");
}

/// Checks that every sample of `reader` is the frame of `expected`.
fn assert_reads(reader: &dyn SampleRead, expected: &Trajectory, context: &str) {
    assert_eq!(reader.len(), expected.len(), "{context}");
    for (index, frame) in expected.iter().enumerate() {
        let sample = reader.read(index).unwrap();
        assert_sample(
            &sample,
            frame.payload(),
            &format!("{context}: sample {index}"),
        );
    }
}

/// `trajectory` with every molecule made whole.
fn whole(trajectory: &Trajectory) -> Trajectory {
    let mut whole = trajectory.clone();
    periodic::make_molecules_whole(&mut whole).unwrap();
    whole
}

/// `trajectory` with the centroid of `atoms` moved to the cell's center.
fn centered(trajectory: &Trajectory, atoms: &AtomSelection) -> Trajectory {
    let mut centered = trajectory.clone();
    for index in 0..centered.len() {
        let mut frame = centered.get_mut(index).unwrap();
        let [a, b, c] = frame.cell().unwrap().vectors().into_value();
        let target = (a + b + c) * 0.5;
        let values = frame.positions().values();
        let points = values.value();
        let mut sum = Vector3::zero();
        for atom in atoms.indices() {
            sum += points[atom.index()] - Point3::origin();
        }
        let shift = target - sum / atoms.indices().len() as f64;
        let transform =
            RigidTransform::new(Matrix3::identity(), Quantity::new(shift, NANOMETER)).unwrap();
        frame.apply_rigid_transform(&transform).unwrap();
    }
    centered
}

/// `trajectory` superposed onto its frame 0 through `atoms`.
fn fitted(trajectory: &Trajectory, atoms: &AtomSelection) -> Trajectory {
    let mut fitted = trajectory.clone();
    fitted.superpose(0, atoms, Default::default()).unwrap();
    fitted
}

#[test]
fn each_step_alone_equals_the_in_memory_functions_and_the_streaming_adapters() {
    let trajectory = source();
    let topology = trajectory.shared_topology();
    let set = set_of(&trajectory);
    let all = AtomSelection::all(Arc::clone(&topology));
    let first = atoms(&trajectory, &[0, 1, 2]);
    let imager = MoleculeImager::new(Arc::clone(&topology));

    // Making whole.
    let expected = whole(&trajectory);
    let processed = Processed::new(&set, ProcessingSteps::default().with_whole(true)).unwrap();
    assert_reads(&processed, &expected, "whole");
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    for (index, frame) in trajectory.iter().enumerate() {
        buffer.copy_from(frame).unwrap();
        imager.make_whole(index, &mut buffer).unwrap();
        assert_sample(&processed.read(index).unwrap(), &buffer, "streamed whole");
    }

    // Wrapping around the first chain.
    let mut expected = trajectory.clone();
    periodic::image_molecules(&mut expected, &first).unwrap();
    let processed =
        Processed::new(&set, ProcessingSteps::default().with_wrap(first.clone())).unwrap();
    assert_reads(&processed, &expected, "wrap");
    for (index, frame) in trajectory.iter().enumerate() {
        buffer.copy_from(frame).unwrap();
        imager.image(index, &mut buffer, &first).unwrap();
        assert_sample(&processed.read(index).unwrap(), &buffer, "streamed wrap");
    }
    // Wrapping makes whole too.
    let both = ProcessingSteps::default()
        .with_whole(true)
        .with_wrap(first.clone());
    assert_reads(
        &Processed::new(&set, both).unwrap(),
        &expected,
        "whole and wrap",
    );

    // Centering.
    let expected = centered(&trajectory, &first);
    let processed =
        Processed::new(&set, ProcessingSteps::default().with_center(first.clone())).unwrap();
    assert_reads(&processed, &expected, "center");

    // Fitting onto sample 0 and onto given positions.
    let expected = fitted(&trajectory, &all);
    let processed = Processed::new(
        &set,
        ProcessingSteps::default().with_fit(all.clone(), FitReference::Sample(0)),
    )
    .unwrap();
    assert_reads(&processed, &expected, "fit");
    let reference = trajectory.get(0).unwrap();
    let superposer = FrameSuperposer::new(&reference, &all, Default::default());
    for (index, frame) in trajectory.iter().enumerate() {
        buffer.copy_from(frame).unwrap();
        superposer.superpose(index, &mut buffer).unwrap();
        assert_sample(&processed.read(index).unwrap(), &buffer, "streamed fit");
    }
    let positions = trajectory.get(3).unwrap().positions().clone();
    let processed = Processed::new(
        &set,
        ProcessingSteps::default().with_fit(first.clone(), FitReference::Positions(positions)),
    )
    .unwrap();
    let mut expected = trajectory.clone();
    let model = Model::new(
        Arc::clone(&topology),
        trajectory.get(3).unwrap().positions().clone(),
    )
    .unwrap();
    expected
        .superpose(&model, &first, Default::default())
        .unwrap();
    assert_reads(&processed, &expected, "fit onto positions");

    // The subset.
    let expected = set.subset(&first).unwrap();
    let processed =
        Processed::new(&set, ProcessingSteps::default().with_subset(first.clone())).unwrap();
    assert_eq!(processed.topology().atom_count(), 3);
    assert_eq!(expected.topology().atom_count(), 3);
    for index in 0..set.len() {
        assert_eq!(
            &processed.read(index).unwrap(),
            expected.get(index).unwrap().payload()
        );
    }

    // No step reads the samples as stored.
    let processed = Processed::new(&set, ProcessingSteps::default()).unwrap();
    assert!(processed.steps().is_empty());
    assert_reads(&processed, &trajectory, "no step");
    assert_eq!(processed.order(), set.order());
}

#[test]
fn every_step_together_equals_the_functions_in_order() {
    let trajectory = source();
    let set = set_of(&trajectory);
    let first = atoms(&trajectory, &[0, 1, 2]);
    let second = atoms(&trajectory, &[3, 4, 5]);
    let fit = atoms(&trajectory, &[0, 1, 2, 4]);

    let mut expected = trajectory.clone();
    periodic::image_molecules(&mut expected, &second).unwrap();
    let expected = fitted(&centered(&expected, &first), &fit);
    let expected = SampleSet::from_trajectory(expected)
        .unwrap()
        .subset(&second)
        .unwrap();

    let steps = ProcessingSteps::default()
        .with_whole(true)
        .with_wrap(second.clone())
        .with_center(first)
        .with_fit(fit, FitReference::Sample(0))
        .with_subset(second);
    let processed = Processed::new(&set, steps.clone()).unwrap();
    for index in 0..set.len() {
        let sample = processed.read(index).unwrap();
        assert_eq!(&sample, expected.get(index).unwrap().payload(), "{index}");
    }
    // Without smoothing, the steps read as before.
    let smoothed = steps.with_smoothing(3).unwrap();
    assert_eq!(smoothed.smoothing(), Some(3));
    let plain = Processed::new(&set, smoothed.without_smoothing()).unwrap();
    for index in 0..set.len() {
        assert_eq!(plain.read(index).unwrap(), processed.read(index).unwrap());
    }
}

/// The centroid of each molecule of `sample`, in fractional coordinates of
/// its cell from the cell's center.
fn molecule_centers_from_cell_center(topology: &Topology, sample: &Sample) -> Vec<[f64; 3]> {
    let cell = sample.cell().unwrap();
    let [a, b, c] = cell.vectors().into_value();
    let middle = (a + b + c) * 0.5;
    let lattice = PeriodicGeometry::new(*cell).unwrap();
    let values = sample.positions().values();
    let points = values.value();
    let count = topology.atom_count() / 2;
    (0..2)
        .map(|molecule| {
            let mut sum = Vector3::zero();
            for atom in molecule * count..(molecule + 1) * count {
                sum += points[atom] - Point3::origin();
            }
            let offset = sum / count as f64 - middle;
            lattice
                .fractional(Quantity::new(offset, NANOMETER))
                .unwrap()
        })
        .collect()
}

#[test]
fn wrapping_then_centering_on_the_anchors_keeps_every_molecule_in_the_cell() {
    let trajectory = source();
    let topology = trajectory.shared_topology();
    let set = set_of(&trajectory);
    for anchors in [
        atoms(&trajectory, &[0, 1, 2]),
        atoms(&trajectory, &[3, 4, 5]),
        AtomSelection::all(Arc::clone(&topology)),
    ] {
        let steps = ProcessingSteps::default()
            .with_wrap(anchors.clone())
            .with_center(anchors);
        let processed = Processed::new(&set, steps).unwrap();
        for index in 0..set.len() {
            let sample = processed.read(index).unwrap();
            for center in molecule_centers_from_cell_center(&topology, &sample) {
                assert!(
                    center.iter().all(|value| value.abs() <= 0.5 + 1.0e-12),
                    "sample {index}: {center:?}"
                );
            }
        }
    }
}

/// What smoothing over `window` samples gives for `stepped`, the samples
/// after the other steps, as the host computed it: the window imaged onto
/// sample `index` and averaged.
fn smoothed_by_hand(stepped: &Trajectory, window: usize, index: usize) -> Vec<Point3> {
    let half = window / 2;
    let imager = MoleculeImager::new(stepped.shared_topology());
    let center = stepped.get(index).unwrap();
    let mut sum = center
        .positions()
        .values()
        .value()
        .iter()
        .map(|point| *point - Point3::origin())
        .collect::<Vec<_>>();
    let mut count = 1.0;
    let mut buffer = FrameBuffer::new(stepped.shared_topology());
    let last = stepped.len() - 1;
    for other in index.saturating_sub(half)..=(index + half).min(last) {
        if other == index {
            continue;
        }
        buffer.copy_from(stepped.get(other).unwrap()).unwrap();
        imager.image_onto(other, &mut buffer, &center).unwrap();
        for (total, point) in sum.iter_mut().zip(*buffer.positions().values().value()) {
            *total += *point - Point3::origin();
        }
        count += 1.0;
    }
    sum.iter().map(|v| Point3::origin() + *v / count).collect()
}

#[test]
fn smoothing_averages_the_window_imaged_onto_each_sample() {
    let trajectory = source();
    let set = set_of(&trajectory);
    let all = AtomSelection::all(trajectory.shared_topology());
    let stepped = fitted(&whole(&trajectory), &all);
    for window in [1, 3, 5, 15] {
        let steps = ProcessingSteps::default()
            .with_whole(true)
            .with_fit(all.clone(), FitReference::Sample(0))
            .with_smoothing(window)
            .unwrap();
        let processed = Processed::new(&set, steps.clone()).unwrap();
        let plain = Processed::new(&set, steps.without_smoothing()).unwrap();
        for index in 0..set.len() {
            let sample = processed.read(index).unwrap();
            let frame = stepped.get(index).unwrap();
            assert_eq!(
                *sample.positions().values().value(),
                &smoothed_by_hand(&stepped, window, index)[..],
                "window {window}, sample {index}"
            );
            // Everything else is the sample's own.
            assert_eq!(sample.cell(), frame.cell());
            assert_eq!(sample.velocities(), frame.velocities());
            assert_eq!(sample.properties(), frame.properties());
            assert_eq!(plain.read(index).unwrap().positions(), frame.positions());
        }
    }
    // Smoothing alone averages the samples as stored, imaged.
    let processed =
        Processed::new(&set, ProcessingSteps::default().with_smoothing(3).unwrap()).unwrap();
    for index in 0..set.len() {
        assert_eq!(
            *processed.read(index).unwrap().positions().values().value(),
            &smoothed_by_hand(&trajectory, 3, index)[..]
        );
    }
}

#[test]
fn smoothing_then_the_subset_keeps_the_smoothed_atoms_of_the_subset() {
    let trajectory = source();
    let set = set_of(&trajectory);
    let second = atoms(&trajectory, &[3, 4, 5]);
    let steps = ProcessingSteps::default()
        .with_whole(true)
        .with_smoothing(3)
        .unwrap()
        .with_subset(second);
    let processed = Processed::new(&set, steps).unwrap();
    assert_eq!(processed.topology().atom_count(), 3);
    let stepped = whole(&trajectory);
    for index in 0..set.len() {
        assert_eq!(
            *processed.read(index).unwrap().positions().values().value(),
            &smoothed_by_hand(&stepped, 3, index)[3..6],
            "sample {index}"
        );
    }
}

#[test]
fn smoothing_images_each_sample_of_the_window_in_its_own_cell() {
    let mut trajectory = source();
    for index in 0..trajectory.len() {
        let edge = 1.0 + 0.03 * index as f64;
        let cell = PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(edge, edge, 1.0), NANOMETER),
            [true; 3],
        )
        .unwrap();
        trajectory.get_mut(index).unwrap().set_cell(Some(cell));
    }
    let set = set_of(&trajectory);
    for window in [3, 5] {
        let steps = ProcessingSteps::default()
            .with_whole(true)
            .with_smoothing(window)
            .unwrap();
        let processed = Processed::new(&set, steps).unwrap();
        let stepped = whole(&trajectory);
        for index in 0..set.len() {
            let sample = processed.read(index).unwrap();
            assert_eq!(
                *sample.positions().values().value(),
                &smoothed_by_hand(&stepped, window, index)[..],
                "window {window}, sample {index}"
            );
            assert_eq!(sample.cell(), trajectory.get(index).unwrap().cell());
        }
    }
}

#[test]
fn samples_without_a_cell_skip_imaging_and_center_where_sample_0_was() {
    let mut trajectory = source();
    for index in 0..trajectory.len() {
        trajectory.get_mut(index).unwrap().set_cell(None);
    }
    let set = set_of(&trajectory);
    let first = atoms(&trajectory, &[0, 1, 2]);
    let steps = ProcessingSteps::default()
        .with_whole(true)
        .with_wrap(first.clone())
        .with_center(first.clone())
        .with_smoothing(3)
        .unwrap();
    // Smoothing reads the samples, so construction checks it too.
    Processed::new(&set, steps.clone()).unwrap();
    let processed = Processed::new(&set, steps.without_smoothing()).unwrap();
    let centroid = |sample: &Sample| {
        let values = sample.positions().values();
        let points = values.value();
        let mut sum = Vector3::zero();
        for atom in first.indices() {
            sum += points[atom.index()] - Point3::origin();
        }
        sum / 3.0
    };
    let rest = centroid(&set.read(0).unwrap());
    for index in 0..set.len() {
        let sample = processed.read(index).unwrap();
        assert!((centroid(&sample) - rest).norm() < 1.0e-12, "{index}");
        assert!(sample.cell().is_none());
    }
}

#[test]
fn a_set_read_from_a_file_processes_like_the_set_in_memory() {
    let trajectory = source();
    let directory = std::env::temp_dir().join(format!("kekule-processing-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("run.trr");
    write_trajectory(
        &path,
        &trajectory,
        TrajectoryWriteOptions::new(TrajectoryFormat::Trr).with_trr_options(
            TrrWriteOptions::default()
                .with_lambda_policy(TrrLambdaPolicy::RequireZero)
                .with_frame_properties(UnrepresentablePolicy::Omit),
        ),
    )
    .unwrap();
    let file =
        open_indexed_sample_set(&path, trajectory.shared_topology(), Default::default()).unwrap();
    let memory = file.load(0..file.len()).unwrap();
    let first = atoms(&trajectory, &[0, 1, 2]);
    let steps = ProcessingSteps::default()
        .with_wrap(first.clone())
        .with_fit(first, FitReference::Sample(2))
        .with_smoothing(3)
        .unwrap();
    let from_file = Processed::new(Arc::new(file), steps.clone()).unwrap();
    let from_memory = Processed::new(memory, steps).unwrap();
    for index in 0..from_file.len() {
        assert_eq!(
            from_file.read(index).unwrap(),
            from_memory.read(index).unwrap()
        );
    }
    assert_eq!(from_file.order(), from_memory.order());
    let (inner, _) = from_file.into_parts();
    drop(inner);
    let _ = std::fs::remove_dir_all(&directory);
}

#[test]
fn steps_that_cannot_process_the_samples_fail_on_construction() {
    let trajectory = source();
    let set = set_of(&trajectory);
    let foreign = AtomSelection::all(linear_carbon_topology(6));
    let none = AtomSelection::empty(trajectory.shared_topology());
    let all = AtomSelection::all(trajectory.shared_topology());

    assert_eq!(
        ProcessingSteps::default().with_smoothing(4).unwrap_err(),
        ProcessingError::InvalidWindow { window: 4 }
    );
    for (steps, step) in [
        (
            ProcessingSteps::default().with_wrap(foreign.clone()),
            ProcessingStep::Wrap,
        ),
        (
            ProcessingSteps::default().with_center(foreign.clone()),
            ProcessingStep::Center,
        ),
        (
            ProcessingSteps::default().with_fit(foreign.clone(), FitReference::Sample(0)),
            ProcessingStep::Fit,
        ),
        (
            ProcessingSteps::default().with_subset(foreign),
            ProcessingStep::Subset,
        ),
    ] {
        assert_eq!(
            Processed::new(&set, steps).unwrap_err(),
            ProcessingError::SelectionTopologyMismatch { step }
        );
    }
    assert_eq!(
        Processed::new(&set, ProcessingSteps::default().with_center(none)).unwrap_err(),
        ProcessingError::EmptySelection {
            step: ProcessingStep::Center
        }
    );
    assert_eq!(
        Processed::new(
            &set,
            ProcessingSteps::default().with_fit(all.clone(), FitReference::Sample(7))
        )
        .unwrap_err(),
        ProcessingError::ReferenceOutOfRange { index: 7, len: 7 }
    );
    let short = Positions::new(Quantity::new([Point3::origin(); 2], NANOMETER)).unwrap();
    assert!(matches!(
        Processed::new(
            &set,
            ProcessingSteps::default().with_fit(all.clone(), FitReference::Positions(short))
        ),
        Err(ProcessingError::Reference(_))
    ));
    // A collapsed reference cannot be fitted onto.
    let collapsed = Positions::new(Quantity::new([Point3::origin(); 6], NANOMETER)).unwrap();
    assert!(matches!(
        Processed::new(
            &set,
            ProcessingSteps::default().with_fit(all, FitReference::Positions(collapsed))
        ),
        Err(ProcessingError::Sample(SampleReadError::Alignment {
            sample: 0,
            ..
        }))
    ));
    // A sample past the last is out of range.
    let processed = Processed::new(&set, ProcessingSteps::default().with_whole(true)).unwrap();
    assert_eq!(
        processed.read(7),
        Err(SampleReadError::IndexOutOfRange { index: 7, len: 7 })
    );
}

/// A reader of its own whose samples have one atom too few.
#[derive(Debug)]
struct Short(SampleSet);

impl SampleRead for Short {
    fn topology(&self) -> &Arc<Topology> {
        SampleRead::topology(&self.0)
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn order(&self) -> Option<&Order> {
        self.0.order()
    }

    fn distribution(&self) -> Option<&Distribution> {
        None
    }

    fn peek(&self, _index: usize) -> Option<SampleView<'_, dyn SampleRead>> {
        None
    }

    fn read(&self, index: usize) -> Result<Sample, SampleReadError> {
        let points = self.0.read(index)?.positions().values().value()[1..].to_vec();
        Ok(Sample::new(
            Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap(),
        ))
    }
}

/// A reader of its own that panics reading sample 4.
#[derive(Debug)]
struct Panicking(SampleSet);

impl SampleRead for Panicking {
    fn topology(&self) -> &Arc<Topology> {
        SampleRead::topology(&self.0)
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn order(&self) -> Option<&Order> {
        self.0.order()
    }

    fn distribution(&self) -> Option<&Distribution> {
        None
    }

    fn peek(&self, _index: usize) -> Option<SampleView<'_, dyn SampleRead>> {
        None
    }

    fn read(&self, index: usize) -> Result<Sample, SampleReadError> {
        assert_ne!(index, 4, "the host's reader fails on sample 4");
        self.0.read(index)
    }
}

#[test]
fn a_panic_in_the_reader_below_leaves_the_other_samples_readable() {
    let trajectory = source();
    let set = set_of(&trajectory);
    let panicking = Panicking(set_of(&trajectory));
    let first = atoms(&trajectory, &[0, 1, 2]);
    let steps = ProcessingSteps::default()
        .with_wrap(first.clone())
        .with_center(first.clone())
        .with_fit(first, FitReference::Sample(0))
        .with_smoothing(3)
        .unwrap();
    let processed = Processed::new(&panicking, steps.clone()).unwrap();
    let expected = Processed::new(&set, steps).unwrap();
    for _ in 0..2 {
        // Samples 3 and 5 smooth over sample 4.
        for index in [3, 4, 5] {
            let read = std::panic::catch_unwind(AssertUnwindSafe(|| processed.read(index)));
            assert!(read.is_err(), "sample {index}");
        }
        for index in [0, 1, 2, 6] {
            assert_eq!(
                processed.read(index).unwrap(),
                expected.read(index).unwrap(),
                "sample {index}"
            );
        }
    }
}

/// A reader of its own that returns samples nothing else holds, as a set in
/// a file does, and keeps the address of each sample's positions.
#[derive(Debug)]
struct Unique {
    set: SampleSet,
    addresses: Mutex<HashMap<usize, usize>>,
}

fn address(sample: &Sample) -> usize {
    sample.positions().values().value().as_ptr() as usize
}

impl SampleRead for Unique {
    fn topology(&self) -> &Arc<Topology> {
        SampleRead::topology(&self.set)
    }

    fn len(&self) -> usize {
        self.set.len()
    }

    fn order(&self) -> Option<&Order> {
        self.set.order()
    }

    fn distribution(&self) -> Option<&Distribution> {
        None
    }

    fn peek(&self, _index: usize) -> Option<SampleView<'_, dyn SampleRead>> {
        None
    }

    fn read(&self, index: usize) -> Result<Sample, SampleReadError> {
        let stored = self.set.read(index)?;
        let points = stored.positions().values().value().to_vec();
        let positions = Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap();
        let model = Model::new(Arc::clone(SampleRead::topology(&self.set)), positions).unwrap();
        let sample = Sample::from_model(model).with_cell(*stored.cell().unwrap());
        self.addresses
            .lock()
            .unwrap()
            .insert(index, address(&sample));
        Ok(sample)
    }
}

#[test]
fn the_steps_change_a_unique_sample_in_place_and_copy_a_shared_one() {
    let trajectory = source();
    let topology = trajectory.shared_topology();
    let set = set_of(&trajectory);
    let first = atoms(&trajectory, &[0, 1, 2]);
    let steps = ProcessingSteps::default()
        .with_whole(true)
        .with_wrap(first.clone())
        .with_center(first)
        .with_fit(AtomSelection::all(topology), FitReference::Sample(0));
    let unique = Unique {
        set: set_of(&trajectory),
        addresses: Mutex::new(HashMap::new()),
    };
    let from_unique = Processed::new(&unique, steps.clone()).unwrap();
    let from_shared = Processed::new(&set, steps).unwrap();
    for index in 0..set.len() {
        // A unique sample is processed in the allocation it was read into.
        let sample = from_unique.read(index).unwrap();
        assert_eq!(
            address(&sample),
            unique.addresses.lock().unwrap()[&index],
            "sample {index}"
        );
        let shared = from_shared.read(index).unwrap();
        assert_eq!(sample.positions(), shared.positions(), "sample {index}");
        // A shared sample is copied, and the set's sample is unchanged.
        let stored = set.read(index).unwrap();
        assert!(!shared.shares_structure_with(&stored), "sample {index}");
        assert_ne!(address(&shared), address(&stored), "sample {index}");
        assert_sample(&stored, trajectory.get(index).unwrap().payload(), "stored");
    }
}

#[test]
fn samples_of_another_layout_fail_instead_of_being_processed() {
    let trajectory = source();
    let first = atoms(&trajectory, &[0, 1, 2]);
    let short = Short(set_of(&trajectory));
    assert!(matches!(
        Processed::new(&short, ProcessingSteps::default().with_center(first)),
        Err(ProcessingError::Sample(SampleReadError::Structure {
            sample: 0,
            ..
        }))
    ));
}

/// A reader of a set in memory that counts the reads of each sample and
/// lends none.
#[derive(Debug)]
struct Counting {
    set: SampleSet,
    reads: Mutex<HashMap<usize, usize>>,
}

impl Counting {
    fn reads(&self, index: usize) -> usize {
        self.reads.lock().unwrap().get(&index).copied().unwrap_or(0)
    }
}

impl SampleRead for Counting {
    fn topology(&self) -> &Arc<Topology> {
        SampleRead::topology(&self.set)
    }

    fn len(&self) -> usize {
        self.set.len()
    }

    fn order(&self) -> Option<&Order> {
        self.set.order()
    }

    fn distribution(&self) -> Option<&Distribution> {
        None
    }

    fn peek(&self, _index: usize) -> Option<SampleView<'_, dyn SampleRead>> {
        None
    }

    fn read(&self, index: usize) -> Result<Sample, SampleReadError> {
        *self.reads.lock().unwrap().entry(index).or_default() += 1;
        self.set.read(index)
    }
}

#[test]
fn the_reader_lends_sample_0_it_processed_on_construction() {
    let trajectory = source();
    let first = atoms(&trajectory, &[0, 1, 2]);
    let steps = ProcessingSteps::default()
        .with_wrap(first.clone())
        .with_center(first.clone())
        .with_fit(
            first,
            FitReference::Positions(trajectory.get(3).unwrap().payload().positions().clone()),
        )
        .with_smoothing(3)
        .unwrap()
        .with_subset(atoms(&trajectory, &[0, 1, 2, 3]));
    let counting = Counting {
        set: set_of(&trajectory),
        reads: Mutex::new(HashMap::new()),
    };
    let processed = Processed::new(&counting, steps.clone()).unwrap();
    let expected = Processed::new(set_of(&trajectory), steps).unwrap();
    let reads = counting.reads(0);
    assert!(reads > 0);

    let lent = processed.peek(0).unwrap();
    assert_eq!(lent.index(), 0);
    assert_eq!(lent.payload(), &expected.read(0).unwrap());
    let read = processed.read(0).unwrap();
    assert!(read.shares_structure_with(lent.payload()));
    assert_eq!(counting.reads(0), reads);
    assert!(processed.peek(1).is_none());
    assert_eq!(processed.read(1).unwrap(), expected.read(1).unwrap());

    let empty = SampleSet::new(trajectory.shared_topology());
    let processed = Processed::new(&empty, ProcessingSteps::default().with_whole(true)).unwrap();
    assert!(processed.peek(0).is_none());
}
