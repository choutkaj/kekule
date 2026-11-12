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
    SampleView, Velocities,
};
use kekule::topology::{AtomSelection, Topology, TopologyAtomIndex, TopologyBuilder};
use kekule::traj::analysis::FrameSuperposer;
use kekule::traj::io::trr::{TrrLambdaPolicy, TrrWriteOptions};
use kekule::traj::io::{open_indexed_sample_set, write_sample_set, TrajectoryWriteOptions};
use kekule::traj::periodic::{self, MoleculeImager};
use kekule::traj::{
    FitReference, FrameBuffer, FrameView, Processed, ProcessingError, ProcessingStep,
    ProcessingSteps, TrajectoryFormat,
};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};

use crate::support::{
    contents, fingerprint_bytes, linear_carbon_topology, reads, time_ordered, Frame, Recorded,
};

/// Seven frames of two three-atom chains crossing the faces of a 1 nm cell,
/// turning a little, with velocities, forces, and a property, in a time
/// order.
fn source() -> SampleSet {
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
    let mut frames = Vec::new();
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
        let mut frame = Frame::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap())
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
        frames.push(frame);
    }
    time_ordered(topology, frames)
}

fn atoms(set: &SampleSet, indices: &[u32]) -> AtomSelection {
    AtomSelection::from_indices(
        set.shared_topology(),
        indices.iter().map(|&index| TopologyAtomIndex::new(index)),
    )
    .unwrap()
}

/// Checks that `sample` holds everything `frame` does.
fn assert_sample(sample: &Sample, frame: &Sample, context: &str) {
    assert_eq!(sample.positions(), frame.positions(), "{context}");
    assert_eq!(sample.cell(), frame.cell(), "{context}");
    assert_eq!(sample.velocities(), frame.velocities(), "{context}");
    assert_eq!(sample.forces(), frame.forces(), "{context}");
    assert_eq!(sample.properties(), frame.properties(), "{context}");
}

/// Checks that every sample of `reader` is the sample of `expected`, and
/// that the samples read are the next recorded ones.
fn assert_reads(
    reader: &dyn SampleRead,
    expected: &SampleSet,
    context: &str,
    recorded: &mut Recorded,
) {
    recorded.check(&reads(reader));
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

/// `set` with every molecule made whole.
fn whole(set: &SampleSet) -> SampleSet {
    let mut whole = set.clone();
    periodic::make_molecules_whole(&mut whole).unwrap();
    whole
}

/// `set` with the centroid of `atoms` moved to the cell's center.
fn centered(set: &SampleSet, atoms: &AtomSelection) -> SampleSet {
    let mut centered = set.clone();
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

/// `set` superposed onto its frame 0 through `atoms`.
fn fitted(set: &SampleSet, atoms: &AtomSelection) -> SampleSet {
    let mut fitted = set.clone();
    fitted.superpose(0, atoms, Default::default()).unwrap();
    fitted
}

#[test]
fn each_step_alone_equals_the_in_memory_functions_and_the_streaming_adapters() {
    let mut recorded = Recorded::new(&[
        0xa549_77de_911b_bb78, // the source
        0xb309_664f_9f45_1e52, // whole
        0xdb65_e88b_76ff_1106, // wrap
        0xdb65_e88b_76ff_1106, // whole and wrap
        0xb079_5215_fe05_ca8c, // center
        0xbd59_005d_be9e_672f, // fit
        0xd87d_000f_0650_7b67, // fit onto positions
        0x5d79_7569_13b9_f2f7, // subset
        0xcf9e_8746_f544_05c8, // no step
    ]);
    let set = source();
    recorded.check(&contents(&set));
    let topology = set.shared_topology();
    let all = AtomSelection::all(Arc::clone(&topology));
    let first = atoms(&set, &[0, 1, 2]);
    let imager = MoleculeImager::new(Arc::clone(&topology));

    // Making whole.
    let expected = whole(&set);
    let processed = Processed::new(&set, ProcessingSteps::default().with_whole(true)).unwrap();
    assert_reads(&processed, &expected, "whole", &mut recorded);
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    for (index, frame) in set.iter().enumerate() {
        buffer
            .copy_from(FrameView::from_sample(frame).unwrap())
            .unwrap();
        imager.make_whole(index, &mut buffer).unwrap();
        assert_sample(&processed.read(index).unwrap(), &buffer, "streamed whole");
    }

    // Wrapping around the first chain.
    let mut expected = set.clone();
    periodic::image_molecules(&mut expected, &first).unwrap();
    let processed =
        Processed::new(&set, ProcessingSteps::default().with_wrap(first.clone())).unwrap();
    assert_reads(&processed, &expected, "wrap", &mut recorded);
    for (index, frame) in set.iter().enumerate() {
        buffer
            .copy_from(FrameView::from_sample(frame).unwrap())
            .unwrap();
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
        &mut recorded,
    );

    // Centering.
    let expected = centered(&set, &first);
    let processed =
        Processed::new(&set, ProcessingSteps::default().with_center(first.clone())).unwrap();
    assert_reads(&processed, &expected, "center", &mut recorded);

    // Fitting onto sample 0 and onto given positions.
    let expected = fitted(&set, &all);
    let processed = Processed::new(
        &set,
        ProcessingSteps::default().with_fit(all.clone(), FitReference::Sample(0)),
    )
    .unwrap();
    assert_reads(&processed, &expected, "fit", &mut recorded);
    let reference = set.get(0).unwrap();
    let superposer = FrameSuperposer::new(&reference, &all, Default::default());
    for (index, frame) in set.iter().enumerate() {
        buffer
            .copy_from(FrameView::from_sample(frame).unwrap())
            .unwrap();
        superposer.superpose(index, &mut buffer).unwrap();
        assert_sample(&processed.read(index).unwrap(), &buffer, "streamed fit");
    }
    let positions = set.get(3).unwrap().positions().clone();
    let processed = Processed::new(
        &set,
        ProcessingSteps::default().with_fit(first.clone(), FitReference::Positions(positions)),
    )
    .unwrap();
    let mut expected = set.clone();
    let model = Model::new(
        Arc::clone(&topology),
        set.get(3).unwrap().positions().clone(),
    )
    .unwrap();
    expected
        .superpose(&model, &first, Default::default())
        .unwrap();
    assert_reads(&processed, &expected, "fit onto positions", &mut recorded);

    // The subset.
    let expected = set.subset(&first).unwrap();
    let processed =
        Processed::new(&set, ProcessingSteps::default().with_subset(first.clone())).unwrap();
    recorded.check(&reads(&processed));
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
    assert_reads(&processed, &set, "no step", &mut recorded);
    assert_eq!(processed.order(), set.order());
}

#[test]
fn every_step_together_equals_the_functions_in_order() {
    let mut recorded = Recorded::new(&[
        0xa549_77de_911b_bb78, // the source
        0x1585_a1a0_e321_46c4, // the functions in order
        0x2c30_0947_d5eb_c322, // processed
    ]);
    let set = source();
    recorded.check(&contents(&set));
    let first = atoms(&set, &[0, 1, 2]);
    let second = atoms(&set, &[3, 4, 5]);
    let fit = atoms(&set, &[0, 1, 2, 4]);

    let mut expected = set.clone();
    periodic::image_molecules(&mut expected, &second).unwrap();
    let expected = fitted(&centered(&expected, &first), &fit);
    let expected = expected.subset(&second).unwrap();
    recorded.check(&contents(&expected));

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
    recorded.check(&reads(&processed));
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
    let mut recorded = Recorded::new(&[
        0xa549_77de_911b_bb78, // the source
        0xf586_3c36_0dbd_08a7, // the first chain
        0x5282_d536_197f_f0bd, // the second chain
        0x2f4b_73d8_061b_8ff1, // every atom
    ]);
    let set = source();
    recorded.check(&contents(&set));
    let topology = set.shared_topology();
    for anchors in [
        atoms(&set, &[0, 1, 2]),
        atoms(&set, &[3, 4, 5]),
        AtomSelection::all(Arc::clone(&topology)),
    ] {
        let steps = ProcessingSteps::default()
            .with_wrap(anchors.clone())
            .with_center(anchors);
        let processed = Processed::new(&set, steps).unwrap();
        recorded.check(&reads(&processed));
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
fn smoothed_by_hand(stepped: &SampleSet, window: usize, index: usize) -> Vec<Point3> {
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
        buffer
            .copy_from(FrameView::from_sample(stepped.get(other).unwrap()).unwrap())
            .unwrap();
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
    let mut recorded = Recorded::new(&[
        0xa549_77de_911b_bb78, // the source
        0x738f_0850_1fc8_470d, // window 1
        0x738f_0850_1fc8_470d, // window 1 unsmoothed
        0xd0f2_4d0b_fe3a_89dd, // window 3
        0x738f_0850_1fc8_470d, // window 3 unsmoothed
        0xf805_325f_129d_8c29, // window 5
        0x738f_0850_1fc8_470d, // window 5 unsmoothed
        0xc14d_ba4b_c065_ca6b, // window 15
        0x738f_0850_1fc8_470d, // window 15 unsmoothed
        0x1c81_cd72_f135_b300, // smoothing alone
    ]);
    let set = source();
    recorded.check(&contents(&set));
    let all = AtomSelection::all(set.shared_topology());
    let stepped = fitted(&whole(&set), &all);
    for window in [1, 3, 5, 15] {
        let steps = ProcessingSteps::default()
            .with_whole(true)
            .with_fit(all.clone(), FitReference::Sample(0))
            .with_smoothing(window)
            .unwrap();
        let processed = Processed::new(&set, steps.clone()).unwrap();
        let plain = Processed::new(&set, steps.without_smoothing()).unwrap();
        recorded.check(&reads(&processed));
        recorded.check(&reads(&plain));
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
    recorded.check(&reads(&processed));
    for index in 0..set.len() {
        assert_eq!(
            *processed.read(index).unwrap().positions().values().value(),
            &smoothed_by_hand(&set, 3, index)[..]
        );
    }
}

#[test]
fn smoothing_then_the_subset_keeps_the_smoothed_atoms_of_the_subset() {
    let mut recorded = Recorded::new(&[
        0xa549_77de_911b_bb78, // the source
        0x50f9_dbed_b335_a8a2, // processed
    ]);
    let set = source();
    recorded.check(&contents(&set));
    let second = atoms(&set, &[3, 4, 5]);
    let steps = ProcessingSteps::default()
        .with_whole(true)
        .with_smoothing(3)
        .unwrap()
        .with_subset(second);
    let processed = Processed::new(&set, steps).unwrap();
    recorded.check(&reads(&processed));
    assert_eq!(processed.topology().atom_count(), 3);
    let stepped = whole(&set);
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
    let mut recorded = Recorded::new(&[
        0x9e0a_7a50_c2a0_9c34, // the source in growing cells
        0xf831_92cf_db4e_6b20, // window 3
        0x4311_3f0f_d9d6_1209, // window 5
    ]);
    let mut set = source();
    for index in 0..set.len() {
        let edge = 1.0 + 0.03 * index as f64;
        let cell = PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(edge, edge, 1.0), NANOMETER),
            [true; 3],
        )
        .unwrap();
        set.get_mut(index).unwrap().set_cell(Some(cell));
    }
    recorded.check(&contents(&set));
    for window in [3, 5] {
        let steps = ProcessingSteps::default()
            .with_whole(true)
            .with_smoothing(window)
            .unwrap();
        let processed = Processed::new(&set, steps).unwrap();
        recorded.check(&reads(&processed));
        let stepped = whole(&set);
        for index in 0..set.len() {
            let sample = processed.read(index).unwrap();
            assert_eq!(
                *sample.positions().values().value(),
                &smoothed_by_hand(&stepped, window, index)[..],
                "window {window}, sample {index}"
            );
            assert_eq!(sample.cell(), set.get(index).unwrap().cell());
        }
    }
}

#[test]
fn samples_without_a_cell_skip_imaging_and_center_where_sample_0_was() {
    let mut recorded = Recorded::new(&[
        0x45cf_5136_e607_3f68, // the source without cells
        0x4c6d_9a00_90c7_bb02, // processed
    ]);
    let mut set = source();
    for index in 0..set.len() {
        set.get_mut(index).unwrap().set_cell(None);
    }
    recorded.check(&contents(&set));
    let first = atoms(&set, &[0, 1, 2]);
    let steps = ProcessingSteps::default()
        .with_whole(true)
        .with_wrap(first.clone())
        .with_center(first.clone())
        .with_smoothing(3)
        .unwrap();
    // Smoothing reads the samples, so construction checks it too.
    Processed::new(&set, steps.clone()).unwrap();
    let processed = Processed::new(&set, steps.without_smoothing()).unwrap();
    recorded.check(&reads(&processed));
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
    let mut recorded = Recorded::new(&[
        0x4b7f_d20b_0ee3_d80c, // processed from the file
    ]);
    let set = source();
    let directory = std::env::temp_dir().join(format!("kekule-processing-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("run.trr");
    write_sample_set(
        &path,
        &set,
        TrajectoryWriteOptions::new(TrajectoryFormat::Trr).with_trr_options(
            TrrWriteOptions::default()
                .with_lambda_policy(TrrLambdaPolicy::RequireZero)
                .with_frame_properties(UnrepresentablePolicy::Omit),
        ),
    )
    .unwrap();
    // The bytes the trajectory of these frames wrote.
    assert_eq!(
        fingerprint_bytes(&std::fs::read(&path).unwrap()),
        0xcb70_4f4e_30e5_b7f1
    );
    let file = open_indexed_sample_set(&path, set.shared_topology(), Default::default()).unwrap();
    let memory = file.load(0..file.len()).unwrap();
    let first = atoms(&set, &[0, 1, 2]);
    let steps = ProcessingSteps::default()
        .with_wrap(first.clone())
        .with_fit(first, FitReference::Sample(2))
        .with_smoothing(3)
        .unwrap();
    let from_file = Processed::new(Arc::new(file), steps.clone()).unwrap();
    let from_memory = Processed::new(memory, steps).unwrap();
    recorded.check(&reads(&from_file));
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
    let set = source();
    let foreign = AtomSelection::all(linear_carbon_topology(6));
    let none = AtomSelection::empty(set.shared_topology());
    let all = AtomSelection::all(set.shared_topology());

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
    let mut recorded = Recorded::new(&[
        0x8131_df28_bb39_9479, // processed from the set
    ]);
    let set = source();
    let panicking = Panicking(set.clone());
    let first = atoms(&set, &[0, 1, 2]);
    let steps = ProcessingSteps::default()
        .with_wrap(first.clone())
        .with_center(first.clone())
        .with_fit(first, FitReference::Sample(0))
        .with_smoothing(3)
        .unwrap();
    let processed = Processed::new(&panicking, steps.clone()).unwrap();
    let expected = Processed::new(&set, steps).unwrap();
    recorded.check(&reads(&expected));
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
    let mut recorded = Recorded::new(&[
        0x5a14_dfd8_4689_5cb5, // processed from the set
    ]);
    let set = source();
    let topology = set.shared_topology();
    let first = atoms(&set, &[0, 1, 2]);
    let steps = ProcessingSteps::default()
        .with_whole(true)
        .with_wrap(first.clone())
        .with_center(first)
        .with_fit(AtomSelection::all(topology), FitReference::Sample(0));
    let unique = Unique {
        set: set.clone(),
        addresses: Mutex::new(HashMap::new()),
    };
    let from_unique = Processed::new(&unique, steps.clone()).unwrap();
    let from_shared = Processed::new(&set, steps).unwrap();
    recorded.check(&reads(&from_shared));
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
        assert_sample(&stored, set.get(index).unwrap().payload(), "stored");
    }
}

#[test]
fn samples_of_another_layout_fail_instead_of_being_processed() {
    let set = source();
    let first = atoms(&set, &[0, 1, 2]);
    let short = Short(set.clone());
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
    let mut recorded = Recorded::new(&[
        0x2242_d7ac_7251_1755, // processed from the set
    ]);
    let set = source();
    let first = atoms(&set, &[0, 1, 2]);
    let steps = ProcessingSteps::default()
        .with_wrap(first.clone())
        .with_center(first.clone())
        .with_fit(
            first,
            FitReference::Positions(set.get(3).unwrap().payload().positions().clone()),
        )
        .with_smoothing(3)
        .unwrap()
        .with_subset(atoms(&set, &[0, 1, 2, 3]));
    let counting = Counting {
        set: set.clone(),
        reads: Mutex::new(HashMap::new()),
    };
    let processed = Processed::new(&counting, steps.clone()).unwrap();
    let expected = Processed::new(set.clone(), steps).unwrap();
    recorded.check(&reads(&expected));
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

    let empty = SampleSet::new(set.shared_topology());
    let processed = Processed::new(&empty, ProcessingSteps::default().with_whole(true)).unwrap();
    assert!(processed.peek(0).is_none());
}
