use super::*;
use crate::core::{Atom, BondOrder, Element, MoleculeEditor};
use crate::geometry::{PeriodicCell, Point3, Vector3};
use crate::properties::{Properties, PropertyColumn, PropertyKey, PropertyValue};
use crate::topology::TopologyBuilder;
use crate::units::{Quantity, ANGSTROM, KELVIN, NANOMETER, SQUARE_ANGSTROM, SQUARE_NANOMETER};
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

fn single_atom_topology() -> crate::topology::Topology {
    let mut editor = MoleculeEditor::new();
    editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    crate::topology::Topology::from_molecule(&editor.finish().unwrap()).unwrap()
}

fn single_position(x: f64) -> Positions {
    Positions::new(Quantity::new([Point3::new(x, 0.0, 0.0)], NANOMETER)).unwrap()
}

#[test]
fn canonical_model_constructors_accept_owned_and_shared_topology() {
    let owned = Model::new(single_atom_topology(), single_position(1.0)).unwrap();
    let shared = owned.shared_topology();
    let shared_model = Model::new(Arc::clone(&shared), single_position(2.0)).unwrap();
    assert!(Arc::ptr_eq(&shared_model.shared_topology(), &shared));

    let owned_complete = Model::with_properties(
        single_atom_topology(),
        single_position(3.0),
        None,
        Properties::realization(1, 0),
    )
    .unwrap();
    let shared_complete = Model::with_properties(
        owned_complete.shared_topology(),
        single_position(4.0),
        None,
        Properties::realization(1, 0),
    )
    .unwrap();
    assert!(owned_complete
        .topology()
        .same_layout(shared_complete.topology()));
}

#[test]
fn canonical_ensemble_constructors_accept_owned_and_shared_topology() {
    let owned_new = Ensemble::new(single_atom_topology());
    let shared = owned_new.shared_topology();
    let shared_new = Ensemble::new(Arc::clone(&shared));
    assert!(Arc::ptr_eq(&shared_new.shared_topology(), &shared));

    let owned_members = Ensemble::from_members(
        single_atom_topology(),
        [EnsembleMember::new(single_position(1.0), 0)],
    )
    .unwrap();
    let shared_members = Ensemble::from_members(
        owned_members.shared_topology(),
        [EnsembleMember::new(single_position(2.0), 0)],
    )
    .unwrap();
    assert!(owned_members
        .topology()
        .same_layout(shared_members.topology()));
}

#[test]
fn ensemble_member_views_project_borrowed_and_owned_models_in_stable_order() {
    let topology = Arc::new(single_atom_topology());
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(10.0, 10.0, 10.0), ANGSTROM),
        [true; 3],
    )
    .unwrap();
    let key = PropertyKey::new("score").unwrap();
    let mut first = EnsembleMember::new(single_position(1.0), 0);
    first.set_cell(Some(cell));
    first
        .insert_property(
            PropertyKey::new("method").unwrap(),
            PropertyValue::String("sampled".into()),
        )
        .unwrap();
    first
        .set_atom_property(0, key.clone(), Some(PropertyValue::Int(7)))
        .unwrap();
    first.set_weight(Some(0.25)).unwrap();
    let mut second = EnsembleMember::new(single_position(2.0), 0);
    second.set_weight(Some(0.75)).unwrap();
    let mut ensemble = Ensemble::from_members(Arc::clone(&topology), [first, second]).unwrap();

    assert_eq!(
        ensemble
            .members()
            .map(|member| member.positions().values().value()[0].x)
            .collect::<Vec<_>>(),
        [1.0, 2.0]
    );
    assert_eq!(
        ensemble
            .members()
            .map(EnsembleMemberView::weight)
            .collect::<Vec<_>>(),
        [Some(0.25), Some(0.75)]
    );

    let member = ensemble.member(0).unwrap();
    assert!(Arc::ptr_eq(&member.shared_topology(), &topology));
    let borrowed = member.as_model();
    assert!(Arc::ptr_eq(&borrowed.shared_topology(), &topology));
    assert_eq!(
        borrowed.positions().values().value().as_ptr(),
        member.positions().values().value().as_ptr()
    );
    let owned = member.to_model();
    assert!(Arc::ptr_eq(&owned.shared_topology(), &topology));
    assert_eq!(owned.positions(), member.positions());
    assert_eq!(owned.cell(), member.cell());
    assert_eq!(owned.properties(), member.properties());

    ensemble
        .member_mut(0)
        .unwrap()
        .set_atom_property(0, key.clone(), Some(PropertyValue::Int(9)))
        .unwrap();
    assert_eq!(
        owned.atom_properties().value(&key, 0).unwrap(),
        Some(PropertyValue::Int(7))
    );
}

#[test]
fn model_view_materialization_clones_realization_and_shares_topology() {
    let mut source = Model::new(single_atom_topology(), single_position(1.0)).unwrap();
    let topology = source.shared_topology();
    let key = PropertyKey::new("score").unwrap();
    let atom = topology.atom_ids()[0];
    source
        .set_atom_property(atom, key.clone(), Some(PropertyValue::Int(4)))
        .unwrap();
    let owned = source.view().to_model();

    assert!(Arc::ptr_eq(&owned.shared_topology(), &topology));
    assert_eq!(owned.positions(), source.positions());
    assert_eq!(owned.properties(), source.properties());
    source
        .set_atom_property(atom, key.clone(), Some(PropertyValue::Int(8)))
        .unwrap();
    source
        .set_position(atom, Quantity::new(Point3::new(9.0, 0.0, 0.0), NANOMETER))
        .unwrap();
    assert_eq!(
        owned.atom_property(atom, &key).unwrap(),
        Some(PropertyValue::Int(4))
    );
    assert_eq!(owned.position(atom).unwrap().value().x, 1.0);
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
        .set_bond_property(
            bond,
            PropertyKey::new("bond_tag").unwrap(),
            Some(PropertyValue::Int(2)),
        )
        .unwrap();
    model
        .set_atom_property(
            atom,
            PropertyKey::new("atom_tag").unwrap(),
            Some(PropertyValue::String("carbon".into())),
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
        model.set_atom_property(
            atom,
            PropertyKey::new("occupancy").unwrap(),
            Some(PropertyValue::Int(1)),
        ),
        Err(ModelError::Property(
            crate::properties::PropertyError::ReservedKey(_)
        ))
    ));
    assert_eq!(
        model
            .atom_property(atom, &PropertyKey::new("atom_tag").unwrap())
            .unwrap(),
        Some(PropertyValue::String("carbon".into()))
    );
    assert_eq!(
        model
            .bond_property(bond, &PropertyKey::new("bond_tag").unwrap())
            .unwrap(),
        Some(PropertyValue::Int(2))
    );
    let view = model.view();
    assert_eq!(
        view.atom_property(atom, &PropertyKey::new("atom_tag").unwrap())
            .unwrap(),
        Some(PropertyValue::String("carbon".into()))
    );
    assert_eq!(
        view.bond_property(bond, &PropertyKey::new("bond_tag").unwrap())
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
        .set_atom_property(
            atom,
            PropertyKey::new("selected").unwrap(),
            Some(PropertyValue::Bool(true)),
        )
        .unwrap();
    model
        .set_bond_property(
            bond,
            PropertyKey::new("bond_selected").unwrap(),
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
        .set_atom_property(
            atom,
            PropertyKey::new("member_atom").unwrap(),
            Some(PropertyValue::Int(1)),
        )
        .unwrap();
    model
        .set_bond_property(
            bond,
            PropertyKey::new("member_bond").unwrap(),
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
    assert_eq!(
        ensemble
            .member(0)
            .unwrap()
            .atom_property(0, &PropertyKey::new("member_atom").unwrap())
            .unwrap(),
        Some(PropertyValue::Int(1))
    );
    assert_eq!(
        ensemble
            .member(0)
            .unwrap()
            .bond_property(0, &PropertyKey::new("member_bond").unwrap())
            .unwrap(),
        Some(PropertyValue::Int(2))
    );

    let member = ensemble.member_mut(0).unwrap();
    member
        .set_atom_property(
            0,
            PropertyKey::new("member_local").unwrap(),
            Some(PropertyValue::Int(3)),
        )
        .unwrap();
    assert_eq!(
        member
            .atom_property(0, &PropertyKey::new("member_local").unwrap())
            .unwrap(),
        Some(PropertyValue::Int(3))
    );
    assert!(matches!(
        member.set_atom_property(
            0,
            PropertyKey::new("b_factor").unwrap(),
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
