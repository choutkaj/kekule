//! Frame selection, replacement, editing, and extraction of sample sets in
//! a time order.

use std::sync::Arc;

use kekule::math::{PeriodicCell, Point3, Vector3};
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{
    Forces, Order, Place, Positions, Sample, SampleSet, SampleSetError, TimeAxis, Velocities,
};
use kekule::topology::{AtomSelection, Topology, TopologyBondIndex};
use kekule::traj::FrameView;
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

fn timed(picoseconds: f64, step: u64) -> Option<Place> {
    Some(Place::Time {
        time: Some(Quantity::new(picoseconds, PICOSECOND)),
        step: Some(step),
    })
}

/// A time order with times and steps, without values yet.
fn time_order() -> Order {
    let times = Quantity::new(Vec::new(), PICOSECOND);
    Order::Time(TimeAxis::new(Some(times), Some(Vec::new())).unwrap())
}

/// Four frames carrying every kind of frame state.
fn annotated() -> SampleSet {
    let mut set = SampleSet::new(chain(3));
    set.set_order(Some(time_order())).unwrap();
    set.properties_mut()
        .insert(PropertyKey::new("run").unwrap(), PropertyValue::Int(7))
        .unwrap();
    for index in 0..4 {
        let mut sample = Sample::new(
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
        sample.set_velocities(Some(Velocities::zeros(3))).unwrap();
        sample.set_forces(Some(Forces::zeros(3))).unwrap();
        sample.set_cell(Some(
            PeriodicCell::orthorhombic(
                Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
                [true; 3],
            )
            .unwrap(),
        ));
        let mut properties = sample.properties_mut();
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
        set.push(sample, timed(index as f64 * 0.5, index * 10))
            .unwrap();
        // Bond rows exist once the sample is bound to the topology.
        set.get_mut(index as usize)
            .unwrap()
            .properties_mut()
            .bonds_mut()
            .insert(
                PropertyKey::new("bond").unwrap(),
                PropertyColumn::Int(vec![Some(4), Some(5)]),
            )
            .unwrap();
    }
    set
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
        let expected = original.get(index).unwrap();
        assert_eq!(actual.payload(), expected.payload());
        assert_eq!(actual.place(), expected.place());
    }
    assert_eq!(
        selected.validate_monotonic_time(),
        Err(SampleSetError::NonMonotonicTime { sample: 1 })
    );
    let strided = original.select((0..4).step_by(2)).unwrap();
    assert_eq!(strided.get(1).unwrap().place(), timed(1.0, 20));
    let empty = original.select([]).unwrap();
    assert!(empty.is_empty());
    assert!(Arc::ptr_eq(
        &empty.shared_topology(),
        &original.shared_topology()
    ));
    assert_eq!(empty.properties(), original.properties());
    assert_eq!(
        original.select([0, 4]).unwrap_err(),
        SampleSetError::IndexOutOfRange { index: 4, len: 4 }
    );
    assert_eq!(format!("{original:?}"), before);
    let atoms = AtomSelection::all(original.shared_topology());
    assert_eq!(original.subset(&atoms).unwrap().len(), 4);
}

#[test]
fn owned_frames_are_complete_and_independent_and_replacement_is_transactional() {
    let mut original = annotated();
    let old = original.get(1).unwrap().payload().clone();
    let mut owned = old.clone();
    owned
        .set_positions(Quantity::new([Point3::origin(); 3], NANOMETER))
        .unwrap();
    assert_eq!(original.get(1).unwrap().payload(), &old);
    assert_eq!(owned.properties(), old.properties());
    assert_eq!(owned.velocities(), old.velocities());
    assert_eq!(owned.forces(), old.forces());
    assert_eq!(owned.cell(), old.cell());
    assert_eq!(original.replace(1, owned.clone()).unwrap(), old);
    assert_eq!(original.get(1).unwrap().payload(), &owned);
    // The replacement takes the replaced frame's time and step.
    assert_eq!(original.get(1).unwrap().place(), timed(0.5, 10));
    let before = format!("{original:?}");
    let invalid = Sample::new(Positions::zeros(2));
    assert_eq!(
        original.replace(4, invalid.clone()).unwrap_err(),
        SampleSetError::IndexOutOfRange { index: 4, len: 4 }
    );
    assert!(original.replace(1, invalid).is_err());
    assert_eq!(format!("{original:?}"), before);
}

#[test]
#[allow(clippy::forget_non_drop)] // Regression: forgetting an editor must never bypass validation.
fn stored_editor_keeps_dimensions_even_after_columns_are_removed_or_editor_is_forgotten() {
    let mut set = annotated();
    let key = PropertyKey::new("bond").unwrap();
    let before = set.get(0).unwrap().payload().clone();
    let smaller = Sample::new(Positions::zeros(2)).properties().clone();
    {
        let mut sample = set.get_mut(0).unwrap();
        assert!(sample.set_velocities(Some(Velocities::zeros(4))).is_err());
        assert!(sample.set_forces(Some(Forces::zeros(4))).is_err());
        assert!(sample
            .set_positions(Quantity::new([Point3::origin(); 2], NANOMETER))
            .is_err());
        assert!(sample.set_properties(smaller).is_err());
    }
    assert_eq!(
        set.set_place(0, timed(f64::NAN, 0)),
        Err(SampleSetError::NonFiniteOrderValue { sample: 0 })
    );
    assert_eq!(set.get(0).unwrap().payload(), &before);
    assert_eq!(set.get(0).unwrap().place(), timed(0.0, 0));
    {
        let mut sample = set.get_mut(0).unwrap();
        let mut properties = sample.properties_mut();
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
        sample
            .set_positions(Quantity::new([Point3::origin(); 3], NANOMETER))
            .unwrap();
        // Invariants do not depend on Drop.
        std::mem::forget(sample);
    }
    set.set_place(0, timed(0.0, 500)).unwrap();
    assert_eq!(set.get(0).unwrap().place(), timed(0.0, 500));
    let topology = set.shared_topology();
    let edited = set.get(0).unwrap().payload().clone();
    FrameView::new(&topology, &edited, None, None).unwrap();
    assert!(set.get_mut(10).is_none());
}

#[test]
fn consuming_extraction_moves_complete_frames() {
    let set = annotated();
    let topology = set.shared_topology();
    let order = set.order().cloned();
    let expected: Vec<_> = set.iter().map(|frame| frame.payload().clone()).collect();
    let pointers: Vec<_> = set
        .iter()
        .map(|frame| {
            (
                frame.positions().values().value().as_ptr(),
                frame.velocities().unwrap().values().value().as_ptr(),
                frame.forces().unwrap().values().value().as_ptr(),
            )
        })
        .collect();
    let frames = set.into_samples();
    assert_eq!(frames, expected);
    for (frame, (positions, velocities, forces)) in frames.iter().zip(pointers) {
        assert_eq!(frame.positions().values().value().as_ptr(), positions);
        assert_eq!(
            frame.velocities().unwrap().values().value().as_ptr(),
            velocities
        );
        assert_eq!(frame.forces().unwrap().values().value().as_ptr(), forces);
    }
    let mut set = SampleSet::from_samples(Arc::clone(&topology), frames).unwrap();
    set.set_order(order).unwrap();
    assert_eq!(set.get(3).unwrap().place(), timed(1.5, 30));
    let pointer = set.get(0).unwrap().positions().values().value().as_ptr();
    let frames = set.into_samples();
    assert_eq!(frames, expected);
    assert_eq!(frames[0].positions().values().value().as_ptr(), pointer);
    assert!(SampleSet::new(topology).into_samples().is_empty());
}
