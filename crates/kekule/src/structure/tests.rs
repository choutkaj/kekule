use super::*;
use crate::core::{Atom, BondOrder, Element, MoleculeEditor};
use crate::geometry::Point3;
use crate::properties::{Properties, PropertyKey, PropertyValue};
use crate::topology::TopologyBuilder;
use crate::units::{Quantity, NANOMETER, SQUARE_ANGSTROM, SQUARE_NANOMETER};
use std::sync::Arc;

fn model_fixture() -> (
    Model,
    crate::topology::InstanceAtomId,
    crate::topology::InstanceBondId,
) {
    let mut editor = MoleculeEditor::new();
    let carbon = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let oxygen = editor
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    let bond = editor.add_bond(carbon, oxygen, BondOrder::Single).unwrap();
    let molecule = editor.finish().unwrap();
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(&molecule).unwrap();
    let topology = Arc::new(builder.build().unwrap());
    let atom = crate::topology::InstanceAtomId::new(instance, carbon);
    let positions = Positions::new(Quantity::new(
        [Point3::origin(), Point3::new(1.0, 0.0, 0.0)],
        NANOMETER,
    ))
    .unwrap();
    (
        Model::new(topology, positions).unwrap(),
        atom,
        crate::topology::InstanceBondId::new(instance, bond),
    )
}

#[test]
fn model_properties_and_canonical_atom_fields_share_one_store() {
    let (mut model, atom, bond) = model_fixture();
    model
        .properties_mut()
        .insert(
            PropertyKey::new("method").unwrap(),
            PropertyValue::String("test".into()),
        )
        .unwrap();
    model.set_occupancy(atom, Some(0.75)).unwrap();
    model
        .set_b_factor(atom, Some(Quantity::new(12.5, SQUARE_ANGSTROM)))
        .unwrap();
    model
        .set_bond_property_value(
            PropertyKey::new("bond_tag").unwrap(),
            bond,
            Some(PropertyValue::Int(2)),
        )
        .unwrap();
    assert_eq!(model.occupancy(atom).unwrap(), Some(0.75));
    assert_eq!(
        model.b_factor(atom).unwrap().unwrap().unit(),
        SQUARE_NANOMETER
    );
    assert!(model
        .properties()
        .atoms()
        .get(&PropertyKey::new("occupancy").unwrap())
        .is_some());
    assert!(matches!(
        model.properties_mut().atoms_mut().set_value(
            PropertyKey::new("occupancy").unwrap(),
            0,
            Some(PropertyValue::Int(1)),
        ),
        Err(crate::properties::PropertyError::ReservedKey(_))
    ));
    assert_eq!(
        model
            .bond_property_value(&PropertyKey::new("bond_tag").unwrap(), bond)
            .unwrap(),
        Some(PropertyValue::Int(2))
    );
}

#[test]
fn model_slice_projects_entity_properties_and_drops_owner_properties() {
    let (mut model, atom, bond) = model_fixture();
    model
        .properties_mut()
        .insert(PropertyKey::new("energy").unwrap(), PropertyValue::Int(3))
        .unwrap();
    model
        .set_atom_property_value(
            PropertyKey::new("selected").unwrap(),
            atom,
            Some(PropertyValue::Bool(true)),
        )
        .unwrap();
    model
        .set_bond_property_value(
            PropertyKey::new("bond_selected").unwrap(),
            bond,
            Some(PropertyValue::String("yes".into())),
        )
        .unwrap();
    let selection =
        crate::topology::AtomSelection::from_atoms(&model.shared_topology(), [atom]).unwrap();
    let sliced = model.slice(&selection).unwrap();
    assert!(sliced.properties().owner_is_empty());
    assert_eq!(
        sliced
            .properties()
            .atoms()
            .value(&PropertyKey::new("selected").unwrap(), 0)
            .unwrap(),
        Some(PropertyValue::Bool(true))
    );
    let all_atoms = model.topology().atom_ids().to_vec();
    let all_selection =
        crate::topology::AtomSelection::from_atoms(&model.shared_topology(), all_atoms).unwrap();
    let fully_retained = model.slice(&all_selection).unwrap();
    assert_eq!(
        fully_retained
            .properties()
            .bonds()
            .value(&PropertyKey::new("bond_selected").unwrap(), 0)
            .unwrap(),
        Some(PropertyValue::String("yes".into()))
    );
}

#[test]
fn ensemble_collection_and_member_properties_are_separate() {
    let (mut model, atom, bond) = model_fixture();
    model
        .set_atom_property_value(
            PropertyKey::new("member_atom").unwrap(),
            atom,
            Some(PropertyValue::Int(1)),
        )
        .unwrap();
    model
        .set_bond_property_value(
            PropertyKey::new("member_bond").unwrap(),
            bond,
            Some(PropertyValue::Int(2)),
        )
        .unwrap();
    let mut ensemble = Ensemble::from_models(&[model]).unwrap();
    ensemble
        .properties_mut()
        .insert(
            PropertyKey::new("collection").unwrap(),
            PropertyValue::Bool(true),
        )
        .unwrap();
    ensemble
        .member_mut(0)
        .unwrap()
        .properties_mut()
        .insert(PropertyKey::new("member").unwrap(), PropertyValue::Int(1))
        .unwrap();
    assert!(!ensemble.properties().owner_is_empty());
    assert!(!ensemble.member(0).unwrap().properties().owner_is_empty());
    assert!(ensemble.member(0).unwrap().properties().atoms().has_data());
    assert!(ensemble.member(0).unwrap().properties().bonds().has_data());
}

#[test]
fn model_rejects_property_dimension_mismatch() {
    let (model, _, _) = model_fixture();
    let wrong = Properties::realization(1, model.topology().bond_count());
    assert!(matches!(
        Model::with_properties(
            model.shared_topology(),
            model.positions().clone(),
            None,
            wrong
        ),
        Err(ModelError::AtomPropertyCountMismatch { .. })
    ));
}
