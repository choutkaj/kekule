//! Samples processed as they are read: each step alone and all of them
//! equal the in-memory functions and the streaming adapters, smoothing
//! averages a window of processed samples imaged onto the sample, and a set
//! read from a file processes like the same set in memory.

use std::sync::Arc;

use kekule::diagnostics::UnrepresentablePolicy;
use kekule::geometry::{Matrix3, PeriodicCell, Point3, RigidTransform, Vector3};
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{
    Distribution, Forces, Order, Positions, Sample, SampleRead, SampleReadError, SampleSet,
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
    let model = kekule::structure::Model::new(
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
    let unsmoothed = processed.unsmoothed();
    for index in 0..set.len() {
        let sample = processed.read(index).unwrap();
        assert_eq!(&sample, expected.get(index).unwrap().payload(), "{index}");
        assert_eq!(unsmoothed.read(index).unwrap(), sample);
    }

    // With smoothing, the layer below is the same steps without it.
    let smoothed = Processed::new(&set, steps.clone().with_smoothing(3).unwrap()).unwrap();
    let plain = Processed::new(&set, steps.without_smoothing()).unwrap();
    for index in 0..set.len() {
        assert_eq!(
            smoothed.unsmoothed().read(index).unwrap(),
            plain.read(index).unwrap()
        );
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
        let processed = Processed::new(&set, steps).unwrap();
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
            assert_eq!(
                processed.unsmoothed().read(index).unwrap().positions(),
                frame.positions()
            );
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
    let processed = Processed::new(&set, steps).unwrap();
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
        let sample = processed.unsmoothed().read(index).unwrap();
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

    fn peek(&self, _index: usize) -> Option<SampleView<'_>> {
        None
    }

    fn read(&self, index: usize) -> Result<Sample, SampleReadError> {
        let points = self.0.read(index)?.positions().values().value()[1..].to_vec();
        Ok(Sample::new(
            Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap(),
        ))
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
