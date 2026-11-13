//! Geometries processed as they are read: each step alone and all of them
//! equal the in-memory functions and the streaming adapters, smoothing
//! averages a window of processed geometries imaged onto the geometry, and a set
//! read from a file processes like the same set in memory.

use std::collections::HashMap;
use std::panic::AssertUnwindSafe;
use std::sync::{Arc, Mutex};

use kekule::diagnostics::UnrepresentablePolicy;
use kekule::math::{Lattice, Matrix3, PeriodicCell, Point3, RigidTransform, Vector3};
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{
    Distribution, Forces, Geometry, GeometryView, Order, Positions, Structure, StructureSet,
    StructureSetRead, StructureSetReadError, Velocities,
};
use kekule::topology::{AtomSelection, Topology, TopologyAtomIndex, TopologyBuilder};
use kekule::traj::analysis::FrameSuperposer;
use kekule::traj::io::trr::{TrrLambdaPolicy, TrrWriteOptions};
use kekule::traj::io::{open_indexed_structure_set, write_structure_set, TrajectoryWriteOptions};
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
fn source() -> StructureSet {
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

fn atoms(set: &StructureSet, indices: &[u32]) -> AtomSelection {
    AtomSelection::from_indices(
        set.shared_topology(),
        indices.iter().map(|&index| TopologyAtomIndex::new(index)),
    )
    .unwrap()
}

/// Checks that `geometry` holds everything `frame` does.
fn assert_geometry(geometry: &Geometry, frame: &Geometry, context: &str) {
    assert_eq!(geometry.positions(), frame.positions(), "{context}");
    assert_eq!(geometry.cell(), frame.cell(), "{context}");
    assert_eq!(geometry.velocities(), frame.velocities(), "{context}");
    assert_eq!(geometry.forces(), frame.forces(), "{context}");
    assert_eq!(geometry.properties(), frame.properties(), "{context}");
}

/// Checks that every geometry of `reader` is the geometry of `expected`, and
/// that the geometries read are the next recorded ones.
fn assert_reads(
    reader: &dyn StructureSetRead,
    expected: &StructureSet,
    context: &str,
    recorded: &mut Recorded,
) {
    recorded.check(&reads(reader));
    assert_eq!(reader.len(), expected.len(), "{context}");
    for (index, frame) in expected.iter().enumerate() {
        let geometry = reader.read(index).unwrap();
        assert_geometry(
            &geometry,
            frame.payload(),
            &format!("{context}: geometry {index}"),
        );
    }
}

/// `set` with every molecule made whole.
fn whole(set: &StructureSet) -> StructureSet {
    let mut whole = set.clone();
    periodic::make_molecules_whole(&mut whole).unwrap();
    whole
}

/// `set` with the centroid of `atoms` moved to the cell's center.
fn centered(set: &StructureSet, atoms: &AtomSelection) -> StructureSet {
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
fn fitted(set: &StructureSet, atoms: &AtomSelection) -> StructureSet {
    let mut fitted = set.clone();
    fitted.superpose(0, atoms, Default::default()).unwrap();
    fitted
}

#[test]
fn each_step_alone_equals_the_in_memory_functions_and_the_streaming_adapters() {
    let mut recorded = Recorded::new(&[
        0x7324_51a7_a937_6e3d, // the source
        0x688c_a6df_25cd_509b, // whole
        0x023c_5470_ede1_9063, // wrap
        0x023c_5470_ede1_9063, // whole and wrap
        0xcc67_9e97_d1dc_f223, // center
        0xc3ff_c9a2_11e5_1140, // fit
        0x5ab2_8aee_a0d0_68b8, // fit onto positions
        0x9bd7_b6f6_1f55_88c2, // subset
        0x90c1_dfd5_4a90_a6eb, // no step
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
            .copy_from(FrameView::from_geometry(frame).unwrap())
            .unwrap();
        imager.make_whole(index, &mut buffer).unwrap();
        assert_geometry(&processed.read(index).unwrap(), &buffer, "streamed whole");
    }

    // Wrapping around the first chain.
    let mut expected = set.clone();
    periodic::image_molecules(&mut expected, &first).unwrap();
    let processed =
        Processed::new(&set, ProcessingSteps::default().with_wrap(first.clone())).unwrap();
    assert_reads(&processed, &expected, "wrap", &mut recorded);
    for (index, frame) in set.iter().enumerate() {
        buffer
            .copy_from(FrameView::from_geometry(frame).unwrap())
            .unwrap();
        imager.image(index, &mut buffer, &first).unwrap();
        assert_geometry(&processed.read(index).unwrap(), &buffer, "streamed wrap");
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

    // Fitting onto geometry 0 and onto given positions.
    let expected = fitted(&set, &all);
    let processed = Processed::new(
        &set,
        ProcessingSteps::default().with_fit(all.clone(), FitReference::Geometry(0)),
    )
    .unwrap();
    assert_reads(&processed, &expected, "fit", &mut recorded);
    let reference = set.get(0).unwrap();
    let superposer = FrameSuperposer::new(&reference, &all, Default::default());
    for (index, frame) in set.iter().enumerate() {
        buffer
            .copy_from(FrameView::from_geometry(frame).unwrap())
            .unwrap();
        superposer.superpose(index, &mut buffer).unwrap();
        assert_geometry(&processed.read(index).unwrap(), &buffer, "streamed fit");
    }
    let positions = set.get(3).unwrap().positions().clone();
    let processed = Processed::new(
        &set,
        ProcessingSteps::default().with_fit(first.clone(), FitReference::Positions(positions)),
    )
    .unwrap();
    let mut expected = set.clone();
    let structure = Structure::new(
        Arc::clone(&topology),
        set.get(3).unwrap().positions().clone(),
    )
    .unwrap();
    expected
        .superpose(&structure, &first, Default::default())
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

    // No step reads the geometries as stored.
    let processed = Processed::new(&set, ProcessingSteps::default()).unwrap();
    assert!(processed.steps().is_empty());
    assert_reads(&processed, &set, "no step", &mut recorded);
    assert_eq!(processed.order(), set.order());
}

#[test]
fn every_step_together_equals_the_functions_in_order() {
    let mut recorded = Recorded::new(&[
        0x7324_51a7_a937_6e3d, // the source
        0xb000_d3c6_d255_1797, // the functions in order
        0x90cf_cd04_baa7_ea17, // processed
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
        .with_fit(fit, FitReference::Geometry(0))
        .with_subset(second);
    let processed = Processed::new(&set, steps.clone()).unwrap();
    for index in 0..set.len() {
        let geometry = processed.read(index).unwrap();
        assert_eq!(&geometry, expected.get(index).unwrap().payload(), "{index}");
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

/// The centroid of each molecule of `geometry`, in fractional coordinates of
/// its cell from the cell's center.
fn molecule_centers_from_cell_center(topology: &Topology, geometry: &Geometry) -> Vec<[f64; 3]> {
    let cell = geometry.cell().unwrap();
    let [a, b, c] = cell.vectors().into_value();
    let middle = (a + b + c) * 0.5;
    let lattice = Lattice::new(*cell).unwrap();
    let values = geometry.positions().values();
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
        0x7324_51a7_a937_6e3d, // the source
        0x002e_84b2_93ff_0174, // the first chain
        0x7f51_eedb_a255_2e7e, // the second chain
        0x8893_61fe_4886_31ea, // every atom
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
            let geometry = processed.read(index).unwrap();
            for center in molecule_centers_from_cell_center(&topology, &geometry) {
                assert!(
                    center.iter().all(|value| value.abs() <= 0.5 + 1.0e-12),
                    "geometry {index}: {center:?}"
                );
            }
        }
    }
}

/// What smoothing over `window` geometries gives for `stepped`, the geometries
/// after the other steps, as the host computed it: the window imaged onto
/// geometry `index` and averaged.
fn smoothed_by_hand(stepped: &StructureSet, window: usize, index: usize) -> Vec<Point3> {
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
            .copy_from(FrameView::from_geometry(stepped.get(other).unwrap()).unwrap())
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
        0x7324_51a7_a937_6e3d, // the source
        0x6c96_d2a6_8979_a782, // window 1
        0x6c96_d2a6_8979_a782, // window 1 unsmoothed
        0x9d03_93be_c02b_7034, // window 3
        0x6c96_d2a6_8979_a782, // window 3 unsmoothed
        0x3cc8_320c_447c_faf2, // window 5
        0x6c96_d2a6_8979_a782, // window 5 unsmoothed
        0x904b_595b_f2dd_24b0, // window 15
        0x6c96_d2a6_8979_a782, // window 15 unsmoothed
        0x24d0_413e_e6bd_d1c9, // smoothing alone
    ]);
    let set = source();
    recorded.check(&contents(&set));
    let all = AtomSelection::all(set.shared_topology());
    let stepped = fitted(&whole(&set), &all);
    for window in [1, 3, 5, 15] {
        let steps = ProcessingSteps::default()
            .with_whole(true)
            .with_fit(all.clone(), FitReference::Geometry(0))
            .with_smoothing(window)
            .unwrap();
        let processed = Processed::new(&set, steps.clone()).unwrap();
        let plain = Processed::new(&set, steps.without_smoothing()).unwrap();
        recorded.check(&reads(&processed));
        recorded.check(&reads(&plain));
        for index in 0..set.len() {
            let geometry = processed.read(index).unwrap();
            let frame = stepped.get(index).unwrap();
            assert_eq!(
                *geometry.positions().values().value(),
                &smoothed_by_hand(&stepped, window, index)[..],
                "window {window}, geometry {index}"
            );
            // Everything else is the geometry's own.
            assert_eq!(geometry.cell(), frame.cell());
            assert_eq!(geometry.velocities(), frame.velocities());
            assert_eq!(geometry.properties(), frame.properties());
            assert_eq!(plain.read(index).unwrap().positions(), frame.positions());
        }
    }
    // Smoothing alone averages the geometries as stored, imaged.
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
        0x7324_51a7_a937_6e3d, // the source
        0x5df4_afd9_7e57_3cc5, // processed
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
            "geometry {index}"
        );
    }
}

#[test]
fn smoothing_images_each_sample_of_the_window_in_its_own_cell() {
    let mut recorded = Recorded::new(&[
        0xaa06_98fc_bcf7_6889, // the source in growing cells
        0x6a28_bcd9_942d_a727, // window 3
        0x75b5_caac_7b62_5ea0, // window 5
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
            let geometry = processed.read(index).unwrap();
            assert_eq!(
                *geometry.positions().values().value(),
                &smoothed_by_hand(&stepped, window, index)[..],
                "window {window}, geometry {index}"
            );
            assert_eq!(geometry.cell(), set.get(index).unwrap().cell());
        }
    }
}

#[test]
fn samples_without_a_cell_skip_imaging_and_center_where_sample_0_was() {
    let mut recorded = Recorded::new(&[
        0xda94_1fe8_2f88_f055, // the source without cells
        0x304b_3c10_03f7_2a17, // processed
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
    // Smoothing reads the geometries, so construction checks it too.
    Processed::new(&set, steps.clone()).unwrap();
    let processed = Processed::new(&set, steps.without_smoothing()).unwrap();
    recorded.check(&reads(&processed));
    let centroid = |geometry: &Geometry| {
        let values = geometry.positions().values();
        let points = values.value();
        let mut sum = Vector3::zero();
        for atom in first.indices() {
            sum += points[atom.index()] - Point3::origin();
        }
        sum / 3.0
    };
    let rest = centroid(&set.read(0).unwrap());
    for index in 0..set.len() {
        let geometry = processed.read(index).unwrap();
        assert!((centroid(&geometry) - rest).norm() < 1.0e-12, "{index}");
        assert!(geometry.cell().is_none());
    }
}

#[test]
fn a_set_read_from_a_file_processes_like_the_set_in_memory() {
    let mut recorded = Recorded::new(&[
        0x4374_27bf_dbd5_f4b7, // processed from the file
    ]);
    let set = source();
    let directory = std::env::temp_dir().join(format!("kekule-processing-{}", std::process::id()));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("run.trr");
    write_structure_set(
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
    let file =
        open_indexed_structure_set(&path, set.shared_topology(), Default::default()).unwrap();
    let memory = file.load(0..file.len()).unwrap();
    let first = atoms(&set, &[0, 1, 2]);
    let steps = ProcessingSteps::default()
        .with_wrap(first.clone())
        .with_fit(first, FitReference::Geometry(2))
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
            ProcessingSteps::default().with_fit(foreign.clone(), FitReference::Geometry(0)),
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
            ProcessingSteps::default().with_fit(all.clone(), FitReference::Geometry(7))
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
        Err(ProcessingError::Geometry(
            StructureSetReadError::Alignment { geometry: 0, .. }
        ))
    ));
    // A geometry past the last is out of range.
    let processed = Processed::new(&set, ProcessingSteps::default().with_whole(true)).unwrap();
    assert_eq!(
        processed.read(7),
        Err(StructureSetReadError::IndexOutOfRange { index: 7, len: 7 })
    );
}

/// A reader of its own whose geometries have one atom too few.
#[derive(Debug)]
struct Short(StructureSet);

impl StructureSetRead for Short {
    fn topology(&self) -> &Arc<Topology> {
        StructureSetRead::topology(&self.0)
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

    fn peek(&self, _index: usize) -> Option<GeometryView<'_, dyn StructureSetRead>> {
        None
    }

    fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        let points = self.0.read(index)?.positions().values().value()[1..].to_vec();
        Ok(Geometry::new(
            Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap(),
        ))
    }
}

/// A reader of its own that panics reading geometry 4.
#[derive(Debug)]
struct Panicking(StructureSet);

impl StructureSetRead for Panicking {
    fn topology(&self) -> &Arc<Topology> {
        StructureSetRead::topology(&self.0)
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

    fn peek(&self, _index: usize) -> Option<GeometryView<'_, dyn StructureSetRead>> {
        None
    }

    fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        assert_ne!(index, 4, "the host's reader fails on geometry 4");
        self.0.read(index)
    }
}

#[test]
fn a_panic_in_the_reader_below_leaves_the_other_samples_readable() {
    let mut recorded = Recorded::new(&[
        0x4bf5_8692_bcc2_a3a0, // processed from the set
    ]);
    let set = source();
    let panicking = Panicking(set.clone());
    let first = atoms(&set, &[0, 1, 2]);
    let steps = ProcessingSteps::default()
        .with_wrap(first.clone())
        .with_center(first.clone())
        .with_fit(first, FitReference::Geometry(0))
        .with_smoothing(3)
        .unwrap();
    let processed = Processed::new(&panicking, steps.clone()).unwrap();
    let expected = Processed::new(&set, steps).unwrap();
    recorded.check(&reads(&expected));
    for _ in 0..2 {
        // Geometries 3 and 5 smooth over geometry 4.
        for index in [3, 4, 5] {
            let read = std::panic::catch_unwind(AssertUnwindSafe(|| processed.read(index)));
            assert!(read.is_err(), "geometry {index}");
        }
        for index in [0, 1, 2, 6] {
            assert_eq!(
                processed.read(index).unwrap(),
                expected.read(index).unwrap(),
                "geometry {index}"
            );
        }
    }
}

/// A reader of its own that returns geometries nothing else holds, as a set in
/// a file does, and keeps the address of each geometry's positions.
#[derive(Debug)]
struct Unique {
    set: StructureSet,
    addresses: Mutex<HashMap<usize, usize>>,
}

fn address(geometry: &Geometry) -> usize {
    geometry.positions().values().value().as_ptr() as usize
}

impl StructureSetRead for Unique {
    fn topology(&self) -> &Arc<Topology> {
        StructureSetRead::topology(&self.set)
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

    fn peek(&self, _index: usize) -> Option<GeometryView<'_, dyn StructureSetRead>> {
        None
    }

    fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        let stored = self.set.read(index)?;
        let points = stored.positions().values().value().to_vec();
        let positions = Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap();
        let structure =
            Structure::new(Arc::clone(StructureSetRead::topology(&self.set)), positions).unwrap();
        let geometry = Geometry::from_structure(structure).with_cell(*stored.cell().unwrap());
        self.addresses
            .lock()
            .unwrap()
            .insert(index, address(&geometry));
        Ok(geometry)
    }
}

#[test]
fn the_steps_change_a_unique_sample_in_place_and_copy_a_shared_one() {
    let mut recorded = Recorded::new(&[
        0xb33d_d873_02fa_717a, // processed from the set
    ]);
    let set = source();
    let topology = set.shared_topology();
    let first = atoms(&set, &[0, 1, 2]);
    let steps = ProcessingSteps::default()
        .with_whole(true)
        .with_wrap(first.clone())
        .with_center(first)
        .with_fit(AtomSelection::all(topology), FitReference::Geometry(0));
    let unique = Unique {
        set: set.clone(),
        addresses: Mutex::new(HashMap::new()),
    };
    let from_unique = Processed::new(&unique, steps.clone()).unwrap();
    let from_shared = Processed::new(&set, steps).unwrap();
    recorded.check(&reads(&from_shared));
    for index in 0..set.len() {
        // A unique geometry is processed in the allocation it was read into.
        let geometry = from_unique.read(index).unwrap();
        assert_eq!(
            address(&geometry),
            unique.addresses.lock().unwrap()[&index],
            "geometry {index}"
        );
        let shared = from_shared.read(index).unwrap();
        assert_eq!(geometry.positions(), shared.positions(), "geometry {index}");
        // A shared geometry is copied, and the set's geometry is unchanged.
        let stored = set.read(index).unwrap();
        assert!(!shared.shares_structure_with(&stored), "geometry {index}");
        assert_ne!(address(&shared), address(&stored), "geometry {index}");
        assert_geometry(&stored, set.get(index).unwrap().payload(), "stored");
    }
}

#[test]
fn samples_of_another_layout_fail_instead_of_being_processed() {
    let set = source();
    let first = atoms(&set, &[0, 1, 2]);
    let short = Short(set.clone());
    assert!(matches!(
        Processed::new(&short, ProcessingSteps::default().with_center(first)),
        Err(ProcessingError::Geometry(StructureSetReadError::Geometry {
            geometry: 0,
            ..
        }))
    ));
}

/// A reader of a set in memory that counts the reads of each geometry and
/// lends none.
#[derive(Debug)]
struct Counting {
    set: StructureSet,
    reads: Mutex<HashMap<usize, usize>>,
}

impl Counting {
    fn reads(&self, index: usize) -> usize {
        self.reads.lock().unwrap().get(&index).copied().unwrap_or(0)
    }
}

impl StructureSetRead for Counting {
    fn topology(&self) -> &Arc<Topology> {
        StructureSetRead::topology(&self.set)
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

    fn peek(&self, _index: usize) -> Option<GeometryView<'_, dyn StructureSetRead>> {
        None
    }

    fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        *self.reads.lock().unwrap().entry(index).or_default() += 1;
        self.set.read(index)
    }
}

#[test]
fn the_reader_lends_sample_0_it_processed_on_construction() {
    let mut recorded = Recorded::new(&[
        0x91c1_dcde_1993_d380, // processed from the set
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

    let empty = StructureSet::new(set.shared_topology());
    let processed = Processed::new(&empty, ProcessingSteps::default().with_whole(true)).unwrap();
    assert!(processed.peek(0).is_none());
}
