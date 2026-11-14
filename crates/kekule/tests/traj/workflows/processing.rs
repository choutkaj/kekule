//! Geometries processed as they are read: each step alone and all of them
//! equal the in-memory functions and the streaming adapters, smoothing
//! averages a window of processed geometries imaged onto the geometry, and a set
//! read from a file processes like the same set in memory.

use kekule::diagnostics::Countable;
use kekule::diagnostics::IndexOutOfRange;
use std::collections::HashMap;
use std::panic::AssertUnwindSafe;
use std::sync::{Arc, Mutex};

use kekule::diagnostics::UnrepresentablePolicy;
use kekule::math::{Lattice, Matrix3, PeriodicCell, Point3, RigidTransform, Vector3};
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{
    Forces, Geometry, GeometryView, Order, Positions, Structure, StructureSet, StructureSetRead,
    StructureSetReadError, Velocities,
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

use crate::support::{linear_carbon_topology, time_ordered, Frame};

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

/// Checks that every geometry of `reader` is the geometry of `expected`.
fn assert_reads(reader: &dyn StructureSetRead, expected: &StructureSet, context: &str) {
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

/// The bonds of [`source`]'s two chains.
const BONDS: [(usize, usize); 4] = [(0, 1), (1, 2), (3, 4), (4, 5)];

const TOLERANCE: f64 = 1.0e-12;

fn points(geometry: &Geometry) -> Vec<Point3> {
    geometry.positions().values().value().to_vec()
}

/// The shortest image of `vector` in [`source`]'s unit cubic cell.
fn minimum_image(vector: Vector3) -> Vector3 {
    Vector3::new(
        vector.x - vector.x.round(),
        vector.y - vector.y.round(),
        vector.z - vector.z.round(),
    )
}

fn centroid(points: &[Point3], atoms: &[usize]) -> Point3 {
    let mut sum = Vector3::zero();
    for &atom in atoms {
        sum += points[atom] - Point3::origin();
    }
    Point3::origin() + sum / atoms.len() as f64
}

/// Checks that every bond of `processed` is the shortest image of that bond
/// in `stored`: each chain is whole, a translated copy of its shortest-image
/// shape.
fn assert_whole(processed: &[Point3], stored: &[Point3], context: &str) {
    for (a, b) in BONDS {
        let bond = processed[b] - processed[a];
        let expected = minimum_image(stored[b] - stored[a]);
        assert!(
            (bond - expected).norm() < TOLERANCE,
            "{context}: bond {a}-{b}: {bond:?}, expected {expected:?}"
        );
    }
}

/// Checks that every atom of `processed` is an image of the atom in `stored`,
/// moved by whole cell vectors.
fn assert_images(processed: &[Point3], stored: &[Point3], context: &str) {
    for (atom, (moved, at)) in processed.iter().zip(stored).enumerate() {
        let shift = *moved - *at;
        for value in [shift.x, shift.y, shift.z] {
            assert!(
                (value - value.round()).abs() < TOLERANCE,
                "{context}: atom {atom} moved by {shift:?}"
            );
        }
    }
}

/// Checks that `processed` is `stored` moved by one translation.
fn assert_translated(processed: &[Point3], stored: &[Point3], context: &str) {
    let shift = processed[0] - stored[0];
    for (atom, (moved, at)) in processed.iter().zip(stored).enumerate() {
        assert!(
            ((*moved - *at) - shift).norm() < TOLERANCE,
            "{context}: atom {atom}"
        );
    }
}

/// Checks that `processed` is `stored` moved rigidly: every distance between
/// two atoms is kept.
fn assert_rigid(processed: &[Point3], stored: &[Point3], context: &str) {
    for a in 0..stored.len() {
        for b in a + 1..stored.len() {
            let kept = (processed[b] - processed[a]).norm();
            let before = (stored[b] - stored[a]).norm();
            assert!(
                (kept - before).abs() < TOLERANCE,
                "{context}: atoms {a} and {b}: {kept} against {before}"
            );
        }
    }
}

fn rmsd(points: &[Point3], reference: &[Point3]) -> f64 {
    let sum = points
        .iter()
        .zip(reference)
        .map(|(point, reference)| (*point - *reference).norm().powi(2))
        .sum::<f64>();
    (sum / points.len() as f64).sqrt()
}

/// The center of [`source`]'s unit cell.
fn cell_center() -> Point3 {
    Point3::new(0.5, 0.5, 0.5)
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
    let set = source();
    let topology = set.shared_topology();
    let all = AtomSelection::all(Arc::clone(&topology));
    let first = atoms(&set, &[0, 1, 2]);
    let imager = MoleculeImager::new(Arc::clone(&topology));

    // Making whole.
    let expected = whole(&set);
    let processed = Processed::new(&set, ProcessingSteps::default().with_whole(true)).unwrap();
    assert_reads(&processed, &expected, "whole");
    // Each chain is whole, its atoms moved by whole cell vectors.
    for (index, stored) in set.iter().enumerate() {
        let read = points(&processed.read(index).unwrap());
        let stored = points(stored.payload());
        let context = format!("whole: geometry {index}");
        assert_whole(&read, &stored, &context);
        assert_images(&read, &stored, &context);
    }
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
    assert_reads(&processed, &expected, "wrap");
    // Each chain is whole, the first centered in the cell and the second at
    // its image nearest the first.
    for (index, stored) in set.iter().enumerate() {
        let read = points(&processed.read(index).unwrap());
        let stored = points(stored.payload());
        let context = format!("wrap: geometry {index}");
        assert_whole(&read, &stored, &context);
        let first = centroid(&read, &[0, 1, 2]);
        assert!((first - cell_center()).norm() < TOLERANCE, "{context}");
        let second = centroid(&read, &[3, 4, 5]) - first;
        assert!(
            [second.x, second.y, second.z]
                .iter()
                .all(|value| value.abs() <= 0.5),
            "{context}: {second:?}"
        );
        // The first chain's first atom moved by the centering translation;
        // the second chain moved by it and whole cell vectors.
        let shift = read[0] - stored[0];
        let moved = stored
            .iter()
            .map(|point| *point + shift)
            .collect::<Vec<_>>();
        assert_images(&read[3..], &moved[3..], &context);
    }
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
    );

    // Centering.
    let expected = centered(&set, &first);
    let processed =
        Processed::new(&set, ProcessingSteps::default().with_center(first.clone())).unwrap();
    assert_reads(&processed, &expected, "center");
    // Every atom moved by one translation, which put the first chain's
    // centroid at the cell's center.
    for (index, stored) in set.iter().enumerate() {
        let read = points(&processed.read(index).unwrap());
        let context = format!("center: geometry {index}");
        assert_translated(&read, &points(stored.payload()), &context);
        assert!(
            (centroid(&read, &[0, 1, 2]) - cell_center()).norm() < TOLERANCE,
            "{context}"
        );
    }

    // Fitting onto geometry 0 and onto given positions.
    let expected = fitted(&set, &all);
    let processed = Processed::new(
        &set,
        ProcessingSteps::default().with_fit(all.clone(), FitReference::Geometry(0)),
    )
    .unwrap();
    assert_reads(&processed, &expected, "fit");
    // Each geometry moved rigidly, no farther from geometry 0 than before,
    // and geometry 0, fitted onto itself, stayed.
    let reference_points = points(set.get(0).unwrap().payload());
    for (index, stored) in set.iter().enumerate() {
        let read = points(&processed.read(index).unwrap());
        let stored = points(stored.payload());
        let context = format!("fit: geometry {index}");
        assert_rigid(&read, &stored, &context);
        assert!(
            rmsd(&read, &reference_points) <= rmsd(&stored, &reference_points) + TOLERANCE,
            "{context}"
        );
    }
    let stayed = points(&processed.read(0).unwrap());
    for (atom, (read, stored)) in stayed.iter().zip(&reference_points).enumerate() {
        assert!((*read - *stored).norm() < TOLERANCE, "fit: atom {atom}");
    }
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
    assert_reads(&processed, &expected, "fit onto positions");
    // Each geometry moved rigidly, and geometry 3, fitted onto its own
    // positions, stayed.
    for (index, stored) in set.iter().enumerate() {
        let read = points(&processed.read(index).unwrap());
        let context = format!("fit onto positions: geometry {index}");
        assert_rigid(&read, &points(stored.payload()), &context);
    }
    let stayed = points(&processed.read(3).unwrap());
    let stored = points(set.get(3).unwrap().payload());
    for (atom, (read, stored)) in stayed.iter().zip(&stored).enumerate() {
        assert!(
            (*read - *stored).norm() < TOLERANCE,
            "fit onto positions: atom {atom}"
        );
    }

    // The subset.
    let expected = set.subset(&first).unwrap();
    let processed =
        Processed::new(&set, ProcessingSteps::default().with_subset(first.clone())).unwrap();
    assert_eq!(processed.topology().atom_count(), 3);
    assert_eq!(expected.topology().atom_count(), 3);
    for index in 0..set.len() {
        let read = processed.read(index).unwrap();
        assert_eq!(&read, expected.get(index).unwrap().payload());
        // The first chain's atoms, as stored.
        let stored = set.get(index).unwrap();
        assert_eq!(points(&read), points(stored.payload())[..3]);
        assert_eq!(read.cell(), stored.cell());
    }

    // No step reads the geometries as stored.
    let processed = Processed::new(&set, ProcessingSteps::default()).unwrap();
    assert!(processed.steps().is_empty());
    assert_reads(&processed, &set, "no step");
    assert_eq!(processed.order(), set.order());
}

#[test]
fn every_step_together_equals_the_functions_in_order() {
    let set = source();
    let first = atoms(&set, &[0, 1, 2]);
    let second = atoms(&set, &[3, 4, 5]);
    let fit = atoms(&set, &[0, 1, 2, 4]);

    let mut expected = set.clone();
    periodic::image_molecules(&mut expected, &second).unwrap();
    let expected = fitted(&centered(&expected, &first), &fit);
    let expected = expected.subset(&second).unwrap();

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
        // The second chain, whole and moved rigidly.
        let read = points(&geometry);
        let stored = points(set.get(index).unwrap().payload());
        let chain = [
            stored[3],
            stored[3] + minimum_image(stored[4] - stored[3]),
            stored[3] + minimum_image(stored[5] - stored[3]),
        ];
        assert_rigid(&read, &chain, &format!("geometry {index}"));
    }
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
    let set = source();
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
        for index in 0..set.len() {
            let geometry = processed.read(index).unwrap();
            assert_whole(
                &points(&geometry),
                &points(set.get(index).unwrap().payload()),
                &format!("geometry {index}"),
            );
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
fn smoothing_averages_the_window_imaged_onto_each_geometry() {
    let set = source();
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
        for index in 0..set.len() {
            let geometry = processed.read(index).unwrap();
            let frame = stepped.get(index).unwrap();
            // A window of one geometry is the geometry alone.
            if window == 1 {
                assert_eq!(geometry.positions(), frame.positions());
            }
            // A window wider than the set averages every geometry.
            if window == 15 {
                assert_eq!(
                    *geometry.positions().values().value(),
                    &smoothed_by_hand(&stepped, 2 * set.len() + 1, index)[..]
                );
            }
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
    for index in 0..set.len() {
        assert_eq!(
            *processed.read(index).unwrap().positions().values().value(),
            &smoothed_by_hand(&set, 3, index)[..]
        );
    }
}

#[test]
fn smoothing_then_the_subset_keeps_the_smoothed_atoms_of_the_subset() {
    let set = source();
    let second = atoms(&set, &[3, 4, 5]);
    let steps = ProcessingSteps::default()
        .with_whole(true)
        .with_smoothing(3)
        .unwrap()
        .with_subset(second);
    let processed = Processed::new(&set, steps).unwrap();
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
fn smoothing_images_each_geometry_of_the_window_in_its_own_cell() {
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
    for window in [3, 5] {
        let steps = ProcessingSteps::default()
            .with_whole(true)
            .with_smoothing(window)
            .unwrap();
        let processed = Processed::new(&set, steps).unwrap();
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
fn geometries_without_a_cell_skip_imaging_and_center_where_geometry_0_was() {
    let mut set = source();
    for index in 0..set.len() {
        set.get_mut(index).unwrap().set_cell(None);
    }
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
        // Without imaging, every atom moved by the one centering translation.
        assert_translated(
            &points(&geometry),
            &points(set.get(index).unwrap().payload()),
            &format!("geometry {index}"),
        );
    }
}

#[test]
fn a_set_read_from_a_file_processes_like_the_set_in_memory() {
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
                .with_properties(UnrepresentablePolicy::Omit),
        ),
    )
    .unwrap();
    let file =
        open_indexed_structure_set(&path, set.shared_topology(), Default::default()).unwrap();
    let memory = file.load(0..file.len()).unwrap();
    // The file holds the frames: their times, steps, cells, and positions to
    // TRR's single precision.
    assert_eq!(memory.order(), set.order());
    assert_eq!(memory.len(), set.len());
    for (index, (read, stored)) in memory.iter().zip(&set).enumerate() {
        assert_eq!(read.cell(), stored.cell(), "geometry {index}");
        for (read, stored) in points(read.payload()).iter().zip(points(stored.payload())) {
            assert!((*read - stored).norm() < 1.0e-6, "geometry {index}");
        }
    }
    let first = atoms(&set, &[0, 1, 2]);
    let steps = ProcessingSteps::default()
        .with_wrap(first.clone())
        .with_fit(first, FitReference::Geometry(2))
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
fn steps_that_cannot_process_the_geometries_fail_on_construction() {
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
        ProcessingError::ReferenceOutOfRange(IndexOutOfRange::new(Countable::Geometry, 7, 7))
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
        Err(StructureSetReadError::IndexOutOfRange(
            IndexOutOfRange::new(Countable::Geometry, 7, 7)
        ))
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

    fn peek(&self, _index: usize) -> Option<GeometryView<'_, dyn StructureSetRead>> {
        None
    }

    fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        assert_ne!(index, 4, "the host's reader fails on geometry 4");
        self.0.read(index)
    }
}

#[test]
fn a_panic_in_the_reader_below_leaves_the_other_geometries_readable() {
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
fn the_steps_change_a_unique_geometry_in_place_and_copy_a_shared_one() {
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
        assert!(!shared.shares_storage_with(&stored), "geometry {index}");
        assert_ne!(address(&shared), address(&stored), "geometry {index}");
        assert_geometry(&stored, set.get(index).unwrap().payload(), "stored");
    }
}

#[test]
fn geometries_of_another_layout_fail_instead_of_being_processed() {
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

    fn peek(&self, _index: usize) -> Option<GeometryView<'_, dyn StructureSetRead>> {
        None
    }

    fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        *self.reads.lock().unwrap().entry(index).or_default() += 1;
        self.set.read(index)
    }
}

#[test]
fn the_reader_lends_geometry_0_it_processed_on_construction() {
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
    let reads = counting.reads(0);
    assert!(reads > 0);

    let lent = processed.peek(0).unwrap();
    assert_eq!(lent.index(), 0);
    assert_eq!(lent.payload(), &expected.read(0).unwrap());
    let read = processed.read(0).unwrap();
    assert!(read.shares_storage_with(lent.payload()));
    assert_eq!(counting.reads(0), reads);
    assert!(processed.peek(1).is_none());
    assert_eq!(processed.read(1).unwrap(), expected.read(1).unwrap());

    let empty = StructureSet::new(set.shared_topology());
    let processed = Processed::new(&empty, ProcessingSteps::default().with_whole(true)).unwrap();
    assert!(processed.peek(0).is_none());
}
