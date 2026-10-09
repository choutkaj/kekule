//! Frame selection, replacement, editing, and extraction of in-memory
//! trajectories.

use std::sync::Arc;

use kekule::core::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{
    Forces, Positions, RealizationError, RealizationView, Trajectory, TrajectoryFrame, Velocities,
};
use kekule::topology::{AtomSelection, Topology, TopologyAtomIndex, TopologyBondIndex};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};

fn chain(atoms: usize) -> Arc<Topology> {
    let mut editor = MoleculeEditor::new();
    let mut previous = None;
    for _ in 0..atoms {
        let atom = editor
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .unwrap();
        if let Some(previous) = previous {
            editor.add_bond(previous, atom, BondOrder::Single).unwrap();
        }
        previous = Some(atom);
    }
    Arc::new(Topology::from_molecule(editor.finish().unwrap()).unwrap())
}

/// Four frames carrying every kind of frame state.
fn annotated() -> Trajectory {
    let mut trajectory = Trajectory::new(chain(3));
    trajectory
        .properties_mut()
        .insert(PropertyKey::new("run").unwrap(), PropertyValue::Int(7))
        .unwrap();
    for index in 0..4 {
        let mut frame = TrajectoryFrame::new(
            Positions::new(Quantity::new(
                [
                    Point3::new(index as f64, 0.0, 0.0),
                    Point3::new(0.1, 0.0, 0.0),
                    Point3::new(0.0, 0.1, 0.0),
                ],
                NANOMETER,
            ))
            .unwrap(),
        );
        frame
            .set_time(Some(Quantity::new(index as f64 * 0.5, PICOSECOND)))
            .unwrap();
        frame.set_step(Some(index * 10));
        frame.set_velocities(Some(Velocities::zeros(3))).unwrap();
        frame.set_forces(Some(Forces::zeros(3))).unwrap();
        let mut conformation = frame.conformation_mut();
        conformation.set_cell(Some(
            PeriodicCell::orthorhombic(
                Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
                [true; 3],
            )
            .unwrap(),
        ));
        let mut properties = conformation.properties_mut();
        properties
            .owner_mut()
            .insert(
                PropertyKey::new("frame").unwrap(),
                PropertyValue::Int(index as i64),
            )
            .unwrap();
        properties
            .atoms_mut()
            .insert(
                PropertyKey::new("atom").unwrap(),
                PropertyColumn::Int(vec![Some(1), None, Some(3)]),
            )
            .unwrap();
        trajectory.push(frame).unwrap();
        // Bond rows exist once the frame is bound to the topology.
        trajectory
            .get_mut(index as usize)
            .unwrap()
            .conformation_mut()
            .properties_mut()
            .bonds_mut()
            .insert(
                PropertyKey::new("bond").unwrap(),
                PropertyColumn::Int(vec![Some(4), Some(5)]),
            )
            .unwrap();
    }
    trajectory
}

#[test]
fn frame_selection_retains_exact_state_order_duplicates_and_empty_topology_binding() {
    let original = annotated();
    let before = format!("{original:?}");
    let selected = original.select([3, 1, 1, 0]).unwrap();
    assert!(Arc::ptr_eq(
        &selected.shared_topology(),
        &original.shared_topology()
    ));
    assert_eq!(selected.properties(), original.properties());
    for (actual, index) in selected.iter().zip([3, 1, 1, 0]) {
        assert_eq!(actual.payload(), original.get(index).unwrap().payload());
    }
    assert_eq!(
        selected.validate_monotonic_time(true),
        Err(RealizationError::NonMonotonicTime { frame: 1 })
    );
    let strided = original.select((0..4).step_by(2)).unwrap();
    assert_eq!(strided.get(1).unwrap().step(), Some(20));
    assert_eq!(
        strided.get(1).unwrap().time(),
        Some(Quantity::new(1.0, PICOSECOND))
    );
    let empty = original.select([]).unwrap();
    assert!(empty.is_empty());
    assert!(Arc::ptr_eq(
        &empty.shared_topology(),
        &original.shared_topology()
    ));
    assert_eq!(empty.properties(), original.properties());
    assert_eq!(
        original.select([0, 4]).unwrap_err(),
        RealizationError::IndexOutOfRange { index: 4, len: 4 }
    );
    assert_eq!(format!("{original:?}"), before);
    let atoms = AtomSelection::all(&original.shared_topology());
    assert_eq!(original.subset(&atoms).unwrap().len(), 4);
}

#[test]
fn owned_frames_are_complete_and_independent_and_replacement_is_transactional() {
    let mut original = annotated();
    let old = original.get(1).unwrap().payload().clone();
    let mut owned = old.clone();
    owned
        .conformation_mut()
        .set_positions(Quantity::new([Point3::origin(); 3], NANOMETER))
        .unwrap();
    owned.set_step(Some(999));
    assert_eq!(original.get(1).unwrap().payload(), &old);
    assert_eq!(owned.properties(), old.properties());
    assert_eq!(owned.velocities(), old.velocities());
    assert_eq!(owned.forces(), old.forces());
    assert_eq!(owned.cell(), old.cell());
    assert_eq!(original.replace(1, owned.clone()).unwrap(), old);
    assert_eq!(original.get(1).unwrap().payload(), &owned);
    let before = format!("{original:?}");
    let invalid = TrajectoryFrame::new(Positions::zeros(2));
    assert_eq!(
        original.replace(4, invalid.clone()).unwrap_err(),
        RealizationError::IndexOutOfRange { index: 4, len: 4 }
    );
    assert!(original.replace(1, invalid).is_err());
    assert_eq!(format!("{original:?}"), before);
}

#[test]
#[allow(clippy::forget_non_drop)] // Regression: forgetting an editor must never bypass validation.
fn stored_editor_keeps_dimensions_even_after_columns_are_removed_or_editor_is_forgotten() {
    let mut trajectory = annotated();
    let key = PropertyKey::new("bond").unwrap();
    let before = trajectory.get(0).unwrap().payload().clone();
    let smaller = TrajectoryFrame::new(Positions::zeros(2))
        .properties()
        .clone();
    {
        let mut frame = trajectory.get_mut(0).unwrap();
        assert!(frame.set_velocities(Some(Velocities::zeros(4))).is_err());
        assert!(frame.set_forces(Some(Forces::zeros(4))).is_err());
        assert!(frame
            .set_time(Some(Quantity::new(f64::NAN, PICOSECOND)))
            .is_err());
        let mut conformation = frame.conformation_mut();
        assert!(conformation
            .set_positions(Quantity::new([Point3::origin(); 2], NANOMETER))
            .is_err());
        assert!(conformation
            .set_occupancy(TopologyAtomIndex::new(0), Some(f64::NAN))
            .is_err());
        assert!(conformation.set_properties(smaller).is_err());
    }
    assert_eq!(trajectory.get(0).unwrap().payload(), &before);
    {
        let mut frame = trajectory.get_mut(0).unwrap();
        let mut conformation = frame.conformation_mut();
        let mut properties = conformation.properties_mut();
        properties.bonds_mut().remove(&key);
        assert!(properties
            .bonds_mut()
            .insert(key.clone(), PropertyColumn::Int(vec![None; 1]))
            .is_err());
        assert_eq!(properties.bonds().len(), 2);
        properties
            .bonds_mut()
            .insert(key.clone(), PropertyColumn::Int(vec![None; 2]))
            .unwrap();
        properties
            .bonds_mut()
            .set_value(key, TopologyBondIndex::new(1), Some(PropertyValue::Int(9)))
            .unwrap();
        conformation
            .set_positions(Quantity::new([Point3::origin(); 3], NANOMETER))
            .unwrap();
        frame.set_step(Some(500));
        // Invariants do not depend on Drop.
        std::mem::forget(frame);
    }
    assert_eq!(trajectory.get(0).unwrap().step(), Some(500));
    let topology = trajectory.shared_topology();
    let edited = trajectory.get(0).unwrap().payload().clone();
    RealizationView::new(&topology, &edited).unwrap();
    assert!(trajectory.get_mut(10).is_none());
}

#[test]
fn consuming_extraction_moves_complete_frames_and_preserves_collection_context() {
    let trajectory = annotated();
    let topology = trajectory.shared_topology();
    let properties = trajectory.properties().clone();
    let expected: Vec<_> = trajectory
        .iter()
        .map(|frame| frame.payload().clone())
        .collect();
    let pointers: Vec<_> = trajectory
        .iter()
        .map(|frame| {
            (
                frame.positions().values().value().as_ptr(),
                frame.velocities().unwrap().values().value().as_ptr(),
                frame.forces().unwrap().values().value().as_ptr(),
            )
        })
        .collect();
    let (actual_topology, actual_properties, frames) = trajectory.into_parts();
    assert!(Arc::ptr_eq(&actual_topology, &topology));
    assert_eq!(actual_properties, properties);
    assert_eq!(frames, expected);
    for (frame, (positions, velocities, forces)) in frames.iter().zip(pointers) {
        assert_eq!(frame.positions().values().value().as_ptr(), positions);
        assert_eq!(
            frame.velocities().unwrap().values().value().as_ptr(),
            velocities
        );
        assert_eq!(frame.forces().unwrap().values().value().as_ptr(), forces);
    }
    let trajectory = Trajectory::from_items(actual_topology, frames).unwrap();
    let pointer = trajectory
        .get(0)
        .unwrap()
        .positions()
        .values()
        .value()
        .as_ptr();
    let frames = trajectory.into_items();
    assert_eq!(frames, expected);
    assert_eq!(frames[0].positions().values().value().as_ptr(), pointer);
    assert!(Trajectory::new(topology).into_items().is_empty());
}
