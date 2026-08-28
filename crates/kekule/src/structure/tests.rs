use super::*;
use crate::core::{Atom, BondOrder, Element, MoleculeEditor};
use crate::geometry::Point3;
use crate::properties::{Properties, PropertyColumn, PropertyKey, PropertyValue};
use crate::topology::TopologyBuilder;
use crate::units::{Quantity, KELVIN, NANOMETER, SQUARE_ANGSTROM, SQUARE_NANOMETER};
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
        .insert_property(
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
        .atom_properties()
        .get(&PropertyKey::new("occupancy").unwrap())
        .is_some());
    assert!(matches!(
        model.set_atom_property_value(
            PropertyKey::new("occupancy").unwrap(),
            atom,
            Some(PropertyValue::Int(1)),
        ),
        Err(ModelError::Property(
            crate::properties::PropertyError::ReservedKey(_)
        ))
    ));
    assert_eq!(
        model
            .bond_property_value(&PropertyKey::new("bond_tag").unwrap(), bond)
            .unwrap(),
        Some(PropertyValue::Int(2))
    );
    assert!(model.set_occupancy(atom, Some(f64::NAN)).is_err());
    assert!(model
        .set_b_factor(atom, Some(Quantity::new(1.0, KELVIN)))
        .is_err());

    let mut malformed = Properties::realization(model.atom_count(), model.topology().bond_count());
    malformed
        .atoms_mut()
        .insert(
            PropertyKey::new("occupancy").unwrap(),
            PropertyColumn::Int(vec![Some(1), None]),
        )
        .unwrap();
    assert!(matches!(
        model.set_properties(malformed),
        Err(ModelError::Property(
            crate::properties::PropertyError::InvalidCanonicalProperty(_)
        ))
    ));
}

#[test]
fn model_slice_projects_entity_properties_and_drops_owner_properties() {
    let (mut model, atom, bond) = model_fixture();
    model
        .insert_property(PropertyKey::new("energy").unwrap(), PropertyValue::Int(3))
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
            .atom_properties()
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
            .bond_properties()
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
        .insert_property(
            PropertyKey::new("collection").unwrap(),
            PropertyValue::Bool(true),
        )
        .unwrap();
    ensemble
        .member_mut(0)
        .unwrap()
        .insert_property(PropertyKey::new("member").unwrap(), PropertyValue::Int(1))
        .unwrap();
    assert!(!ensemble.properties().owner_is_empty());
    assert!(!ensemble.member(0).unwrap().properties().owner_is_empty());
    assert!(ensemble.member(0).unwrap().atom_properties().has_data());
    assert!(ensemble.member(0).unwrap().bond_properties().has_data());

    let member = ensemble.member_mut(0).unwrap();
    assert!(matches!(
        member.set_atom_property_value(
            PropertyKey::new("b_factor").unwrap(),
            0,
            Some(PropertyValue::Int(2)),
        ),
        Err(EnsembleError::Property(error))
            if matches!(*error, crate::properties::PropertyError::ReservedKey(_))
    ));
    member.set_occupancy_at(0, Some(0.5)).unwrap();
    member
        .set_b_factor_at(0, Some(Quantity::new(10.0, SQUARE_ANGSTROM)))
        .unwrap();
    assert_eq!(member.occupancy_at(0).unwrap(), Some(0.5));
    assert_eq!(
        member.b_factor_at(0).unwrap().unwrap().unit(),
        SQUARE_NANOMETER
    );
    assert!(member.set_occupancy_at(0, Some(f64::INFINITY)).is_err());
    assert!(member
        .set_b_factor_at(0, Some(Quantity::new(1.0, KELVIN)))
        .is_err());

    let mut malformed =
        Properties::realization(member.positions().len(), member.bond_properties().len());
    malformed
        .atoms_mut()
        .insert(
            PropertyKey::new("b_factor").unwrap(),
            PropertyColumn::String(vec![Some("bad".into()), None]),
        )
        .unwrap();
    assert!(matches!(
        member.set_properties(malformed),
        Err(EnsembleError::Property(error))
            if matches!(
                *error,
                crate::properties::PropertyError::InvalidCanonicalProperty(_)
            )
    ));

    let cloned = ensemble.clone();
    assert_eq!(cloned.properties(), ensemble.properties());
    let selection =
        crate::topology::AtomSelection::from_atoms(&ensemble.shared_topology(), [atom]).unwrap();
    let sliced = ensemble.slice(&selection).unwrap();
    assert!(sliced.properties().owner_is_empty());
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
