use kekule::diagnostics::TopologyMismatch;
use std::sync::Arc;

use kekule::math::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{Forces, Place, Positions, StructureSet, Velocities};
use kekule::topology::{AtomSelection, Topology, TopologyAtomIndex, TopologyBuilder};
use kekule::traj::periodic::{self, PeriodicError};
use kekule::units::{
    Quantity, CANONICAL_FORCE_UNIT, CANONICAL_VELOCITY_UNIT, NANOMETER, PICOSECOND,
};

use crate::support::{linear_carbon_topology, time_ordered, topology, Frame};

fn cell(length: f64) -> PeriodicCell {
    PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(length, length, length), NANOMETER),
        [true; 3],
    )
    .unwrap()
}

fn frame(points: &[[f64; 3]], cell: Option<PeriodicCell>) -> Frame {
    let points = points
        .iter()
        .map(|p| Point3::new(p[0], p[1], p[2]))
        .collect::<Vec<_>>();
    let mut frame = Frame::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap());
    frame.set_cell(cell);
    frame
}

fn xyz(set: &StructureSet, frame: usize) -> Vec<Point3> {
    set.get(frame)
        .unwrap()
        .positions()
        .values()
        .value()
        .to_vec()
}

/// Asserts that `copy` has the geometries of `set` with only their positions
/// moved: the same count, and every cell, velocity, and force unchanged.
fn assert_only_positions_moved(set: &StructureSet, copy: &StructureSet) {
    assert_eq!(copy.len(), set.len());
    for (original, moved) in set.iter().zip(copy.iter()) {
        assert_eq!(moved.positions().len(), original.positions().len());
        assert_eq!(moved.cell(), original.cell());
        assert_eq!(moved.velocities(), original.velocities());
        assert_eq!(moved.forces(), original.forces());
    }
}

/// Copying forms of the in-place periodic operations, each checking that
/// the operation moved positions alone.
fn whole(set: &StructureSet) -> Result<StructureSet, PeriodicError> {
    let mut copy = set.clone();
    periodic::make_molecules_whole(&mut copy)?;
    assert_only_positions_moved(set, &copy);
    Ok(copy)
}

fn imaged(set: &StructureSet, anchors: &AtomSelection) -> Result<StructureSet, PeriodicError> {
    let mut copy = set.clone();
    periodic::image_molecules(&mut copy, anchors)?;
    assert_only_positions_moved(set, &copy);
    Ok(copy)
}

fn unwrapped(set: &StructureSet) -> Result<StructureSet, PeriodicError> {
    let mut copy = set.clone();
    periodic::unwrap_trajectory(&mut copy)?;
    assert_only_positions_moved(set, &copy);
    Ok(copy)
}

fn close(actual: Point3, expected: [f64; 3]) {
    let expected = Point3::new(expected[0], expected[1], expected[2]);
    assert!(
        (actual - expected).norm() < 1.0e-10,
        "{actual:?} != {expected:?}"
    );
}

#[test]
fn making_whole_preserves_complete_metadata_source_and_topology_and_supports_in_place() {
    let topology = linear_carbon_topology(2);
    let mut source_frame = frame(&[[0.9, 0.2, 0.3], [0.1, 0.2, 0.3]], Some(cell(1.0)));
    source_frame
        .set_time(Some(Quantity::new(2.5, PICOSECOND)))
        .unwrap();
    source_frame.set_step(Some(25));
    source_frame
        .set_velocities(Some(
            Velocities::new(Quantity::new(
                [Vector3::new(1.0, 2.0, 3.0); 2],
                CANONICAL_VELOCITY_UNIT,
            ))
            .unwrap(),
        ))
        .unwrap();
    source_frame
        .set_forces(Some(
            Forces::new(Quantity::new(
                [Vector3::new(4.0, 5.0, 6.0); 2],
                CANONICAL_FORCE_UNIT,
            ))
            .unwrap(),
        ))
        .unwrap();
    source_frame
        .properties_mut()
        .owner_mut()
        .insert(PropertyKey::new("frame").unwrap(), PropertyValue::Int(7))
        .unwrap();
    source_frame
        .properties_mut()
        .atoms_mut()
        .insert(
            PropertyKey::new("atom").unwrap(),
            PropertyColumn::Int(vec![Some(1), Some(2)]),
        )
        .unwrap();
    let mut source = time_ordered(topology.clone(), [source_frame]);
    source
        .get_mut(0)
        .unwrap()
        .properties_mut()
        .bonds_mut()
        .insert(
            PropertyKey::new("bond").unwrap(),
            PropertyColumn::Int(vec![Some(3)]),
        )
        .unwrap();
    source
        .properties_mut()
        .insert(
            PropertyKey::new("run").unwrap(),
            PropertyValue::String("original".into()),
        )
        .unwrap();
    let before = format!("{source:?}");
    let whole = whole(&source).unwrap();
    close(xyz(&whole, 0)[0], [0.9, 0.2, 0.3]);
    close(xyz(&whole, 0)[1], [1.1, 0.2, 0.3]);
    assert_eq!(format!("{source:?}"), before);
    assert!(Arc::ptr_eq(&topology, &whole.shared_topology()));
    assert_eq!(whole.properties(), source.properties());
    let original = source.get(0).unwrap();
    let transformed = whole.get(0).unwrap();
    assert_eq!(transformed.cell(), original.cell());
    assert_eq!(transformed.velocities(), original.velocities());
    assert_eq!(transformed.forces(), original.forces());
    assert_eq!(transformed.place(), original.place());
    assert_eq!(
        transformed.place(),
        Some(Place::Time {
            time: Some(Quantity::new(2.5, PICOSECOND)),
            step: Some(25)
        })
    );
    assert_eq!(transformed.properties(), original.properties());
    let velocity_pointer = source
        .get(0)
        .unwrap()
        .velocities()
        .unwrap()
        .values()
        .value()
        .as_ptr();
    periodic::make_molecules_whole(&mut source).unwrap();
    assert_eq!(format!("{source:?}"), format!("{whole:?}"));
    assert_eq!(
        source
            .get(0)
            .unwrap()
            .velocities()
            .unwrap()
            .values()
            .value()
            .as_ptr(),
        velocity_pointer
    );
}

#[test]
fn making_whole_uses_true_shortest_images_for_skewed_rotated_and_partial_cells() {
    let topology = linear_carbon_topology(2);
    for rotate in [false, true] {
        let rotate_vector = |v: Vector3| {
            if rotate {
                Vector3::new(-v.y, v.x, v.z)
            } else {
                v
            }
        };
        for (axes, expected) in [
            ([true; 3], [0.031, -0.102, 0.0]),
            ([true, false, true], [-0.069, 0.098, 0.0]),
        ] {
            let basis = [
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.9, 0.2, 0.0),
                Vector3::new(0.0, 0.0, 1.0),
            ]
            .map(rotate_vector);
            let cell = PeriodicCell::new(Quantity::new(basis, NANOMETER), axes).unwrap();
            let delta = rotate_vector(Vector3::new(0.931, 0.098, 0.0));
            let set = time_ordered(
                topology.clone(),
                [frame(&[[0.0; 3], [delta.x, delta.y, delta.z]], Some(cell))],
            );
            let whole = whole(&set).unwrap();
            let expected = rotate_vector(Vector3::new(expected[0], expected[1], expected[2]));
            // The first atom anchors the molecule; the second moves to the
            // shortest image of the bond.
            close(xyz(&whole, 0)[0], [0.0; 3]);
            close(xyz(&whole, 0)[1], [expected.x, expected.y, expected.z]);
        }
    }
}

#[test]
fn reconstruction_checks_ring_closures_and_rolls_back_late_failures() {
    let topology = topology(&["C", "C", "C"], &[(0, 1), (1, 2), (2, 0)]);
    let good = frame(
        &[[0.9, 0.0, 0.0], [0.1, 0.0, 0.0], [0.0, 0.1, 0.0]],
        Some(cell(1.0)),
    );
    let bad = frame(
        &[[0.1, 0.0, 0.0], [0.4, 0.0, 0.0], [0.8, 0.0, 0.0]],
        Some(cell(1.0)),
    );
    let mut set = time_ordered(topology, [good, bad]);
    let before = format!("{set:?}");
    assert!(matches!(
        periodic::make_molecules_whole(&mut set),
        Err(PeriodicError::InconsistentBondImages { frame: 1, .. })
    ));
    assert_eq!(format!("{set:?}"), before);
    assert!(matches!(
        whole(&set),
        Err(PeriodicError::InconsistentBondImages { frame: 1, .. })
    ));
}

#[test]
fn ring_closure_validation_is_independent_of_cell_aspect_ratio() {
    let topology = topology(&["C", "C", "C"], &[(0, 1), (1, 2), (2, 0)]);
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(1.0e-10, 10.0, 10.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    let set = time_ordered(
        topology,
        [frame(
            &[
                [0.1e-10, 0.0, 0.0],
                [0.4e-10, 0.0, 0.0],
                [0.8e-10, 0.0, 0.0],
            ],
            Some(cell),
        )],
    );
    assert!(matches!(
        whole(&set),
        Err(PeriodicError::InconsistentBondImages { frame: 0, .. })
    ));
}

#[test]
fn reconstruction_reuses_bond_image_ties_when_traversal_reverses_an_edge() {
    let topology = topology(&["C", "C", "C"], &[(0, 2), (1, 2)]);
    let cell = PeriodicCell::new(
        Quantity::new(
            [
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.9, 0.2, 0.0),
                Vector3::new(0.0, 0.0, 1.0),
            ],
            NANOMETER,
        ),
        [true; 3],
    )
    .unwrap();
    let set = time_ordered(
        topology,
        [frame(
            &[[0.0; 3], [0.51, 0.0, 0.0], [0.01, 0.0, 0.0]],
            Some(cell),
        )],
    );
    let whole = whole(&set).unwrap();
    let points = xyz(&whole, 0);
    close(points[0], [0.0; 3]);
    close(points[2], [0.01, 0.0, 0.0]);
    // There are two equally short images of the 1--2 bond, whose raw
    // vector (0.5, 0, 0) shifts by -b or by a - b to (-0.4, -0.2, 0) or
    // (0.4, 0.2, 0), both of squared length 0.2.
    assert!(((points[1] - points[2]).norm_squared() - 0.2).abs() < 1.0e-14);
    let bond = points[1] - points[2];
    assert!(
        (bond - Vector3::new(-0.4, -0.2, 0.0)).norm() < 1.0e-10
            || (bond - Vector3::new(0.4, 0.2, 0.0)).norm() < 1.0e-10,
        "{bond:?}"
    );
}

fn imaging_topology() -> Arc<Topology> {
    let molecule = kekule::smiles::to_molecules("CC").unwrap().pop().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(molecule.clone()).unwrap();
    builder.add_instance(definition).unwrap();
    builder.add_instance(definition).unwrap();
    builder
        .add_molecule((kekule::smiles::to_molecules("O").unwrap().pop().unwrap()).clone())
        .unwrap();
    Arc::new(builder.build().unwrap())
}

#[test]
fn imaging_expands_anchor_atoms_to_molecules_and_centers_complete_groups() {
    let topology = imaging_topology();
    let source = frame(
        &[
            [0.9, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [0.2, 0.0, 0.0],
            [0.3, 0.0, 0.0],
            [0.95, 0.0, 0.0],
        ],
        Some(cell(1.0)),
    );
    let mut set = time_ordered(topology.clone(), [source]);
    let anchors = AtomSelection::from_indices(
        Arc::clone(&topology),
        [TopologyAtomIndex::new(0), TopologyAtomIndex::new(3)],
    )
    .unwrap();
    let before = format!("{set:?}");
    let imaged = imaged(&set, &anchors).unwrap();
    for (actual, x) in xyz(&imaged, 0)
        .into_iter()
        .zip([0.275, 0.475, 0.575, 0.675, 0.325])
    {
        close(actual, [x, 0.5, 0.5]);
    }
    assert_eq!(format!("{set:?}"), before);
    assert!(Arc::ptr_eq(&topology, &imaged.shared_topology()));
    periodic::image_molecules(&mut set, &anchors).unwrap();
    assert_eq!(format!("{set:?}"), format!("{imaged:?}"));
    assert_eq!(imaged.get(0).unwrap().cell(), set.get(0).unwrap().cell());
}

#[test]
fn imaging_rejects_empty_or_foreign_anchors_and_preserves_nonperiodic_coordinates() {
    let topology = linear_carbon_topology(2);
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER),
        [true, false, false],
    )
    .unwrap();
    let set = time_ordered(
        topology.clone(),
        [frame(&[[0.9, 2.0, 3.0], [0.1, 2.1, 3.2]], Some(cell))],
    );
    let empty = AtomSelection::from_atoms(Arc::clone(&topology), []).unwrap();
    assert!(matches!(
        imaged(&set, &empty),
        Err(PeriodicError::EmptyAnchors)
    ));
    let foreign = AtomSelection::all(Arc::clone(&linear_carbon_topology(2)));
    assert!(matches!(
        imaged(&set, &foreign),
        Err(PeriodicError::TopologyMismatch(TopologyMismatch::Selection))
    ));
    let imaged = imaged(&set, &AtomSelection::all(Arc::clone(&topology))).unwrap();
    close(xyz(&imaged, 0)[0], [0.4, 2.0, 3.0]);
    close(xyz(&imaged, 0)[1], [0.6, 2.1, 3.2]);
}

#[test]
fn unwrapping_tracks_multiple_crossings_and_preserves_the_first_frame() {
    let topology = linear_carbon_topology(1);
    let frames = [0.9, 0.1, 0.4, 0.8, 0.2].map(|x| frame(&[[x, 0.0, 0.0]], Some(cell(1.0))));
    let mut set = time_ordered(topology.clone(), frames);
    let before = format!("{set:?}");
    let unwrapped = unwrapped(&set).unwrap();
    assert_eq!(xyz(&unwrapped, 0), xyz(&set, 0));
    for (frame, expected) in [0.9, 1.1, 1.4, 1.8, 2.2].into_iter().enumerate() {
        close(xyz(&unwrapped, frame)[0], [expected, 0.0, 0.0]);
    }
    assert!(Arc::ptr_eq(&topology, &unwrapped.shared_topology()));
    assert_eq!(format!("{set:?}"), before);
    periodic::unwrap_trajectory(&mut set).unwrap();
    assert_eq!(format!("{set:?}"), format!("{unwrapped:?}"));
}

#[test]
fn unwrapping_uses_current_triclinic_cells_and_only_periodic_fractional_axes() {
    let topology = linear_carbon_topology(1);
    let first = PeriodicCell::new(
        Quantity::new(
            [
                Vector3::new(1.0, 0.0, 0.0),
                Vector3::new(0.4, 1.0, 0.0),
                Vector3::new(0.0, 0.0, 1.0),
            ],
            NANOMETER,
        ),
        [true, false, true],
    )
    .unwrap();
    let second = PeriodicCell::new(
        Quantity::new(
            [
                Vector3::new(2.0, 0.0, 0.0),
                Vector3::new(0.6, 1.5, 0.0),
                Vector3::new(0.0, 0.0, 1.0),
            ],
            NANOMETER,
        ),
        [true, false, true],
    )
    .unwrap();
    // Fractional positions (0.9, 0.1, 0.0) -> (0.1, 1.4, 0.0).
    let set = time_ordered(
        topology,
        [
            frame(&[[0.94, 0.1, 0.0]], Some(first)),
            frame(&[[1.04, 2.1, 0.0]], Some(second)),
        ],
    );
    let unwrapped = unwrapped(&set).unwrap();
    // The first frame is the reference and stays.
    close(xyz(&unwrapped, 0)[0], [0.94, 0.1, 0.0]);
    // Unwrapped fractional position (1.1, 1.4, 0.0) in the second cell.
    close(xyz(&unwrapped, 1)[0], [3.04, 2.1, 0.0]);
    assert_eq!(unwrapped.get(1).unwrap().cell(), Some(&second));
}

#[test]
fn periodic_failures_are_transactional_and_include_frame_context() {
    let topology = linear_carbon_topology(1);
    let mut missing = time_ordered(
        topology.clone(),
        [
            frame(&[[0.9, 0.0, 0.0]], Some(cell(1.0))),
            frame(&[[0.1, 0.0, 0.0]], None),
        ],
    );
    let before = format!("{missing:?}");
    assert!(matches!(
        periodic::unwrap_trajectory(&mut missing),
        Err(PeriodicError::MissingCell { frame: 1 })
    ));
    assert!(matches!(
        periodic::make_molecules_whole(&mut missing),
        Err(PeriodicError::MissingCell { frame: 1 })
    ));
    assert!(matches!(
        periodic::image_molecules(&mut missing, &AtomSelection::all(Arc::clone(&topology))),
        Err(PeriodicError::MissingCell { frame: 1 })
    ));
    assert_eq!(format!("{missing:?}"), before);

    let mut ambiguous = time_ordered(
        topology.clone(),
        [
            frame(&[[0.9, 0.0, 0.0]], Some(cell(1.0))),
            frame(&[[0.4, 0.0, 0.0]], Some(cell(1.0))),
        ],
    );
    let before = format!("{ambiguous:?}");
    assert!(matches!(
        periodic::unwrap_trajectory(&mut ambiguous),
        Err(PeriodicError::AmbiguousDisplacement {
            frame: 1,
            axis: 0,
            ..
        })
    ));
    assert_eq!(format!("{ambiguous:?}"), before);

    let partial = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER),
        [true, false, false],
    )
    .unwrap();
    let changed = time_ordered(
        topology,
        [
            frame(&[[0.0; 3]], Some(cell(1.0))),
            frame(&[[0.1; 3]], Some(partial)),
        ],
    );
    assert!(matches!(
        unwrapped(&changed),
        Err(PeriodicError::PeriodicAxesChanged { frame: 1 })
    ));
}

#[test]
fn empty_trajectories_and_single_atom_molecules_need_no_special_cases() {
    let topology = linear_carbon_topology(1);
    let empty = StructureSet::new(topology.clone());
    assert!(whole(&empty).unwrap().is_empty());
    assert!(imaged(&empty, &AtomSelection::all(Arc::clone(&topology)))
        .unwrap()
        .is_empty());
    assert!(unwrapped(&empty).unwrap().is_empty());
    let single = time_ordered(topology, [frame(&[[3.0, 4.0, 5.0]], Some(cell(1.0)))]);
    assert_eq!(xyz(&whole(&single).unwrap(), 0), xyz(&single, 0));
    assert_eq!(xyz(&unwrapped(&single).unwrap(), 0), xyz(&single, 0));
}
