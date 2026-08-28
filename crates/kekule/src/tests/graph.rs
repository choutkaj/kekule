use super::*;
use crate::properties::{PropertyKey, PropertyValue};

#[test]
fn empty_molecule_has_no_atoms_or_bonds() {
    let mol = crate::core::MoleculeEditor::new();

    assert_eq!(mol.atom_count(), 0);
    assert_eq!(mol.bond_count(), 0);
    assert_eq!(mol.formal_charge(), 0);
    assert!(mol.atoms().next().is_none());
    assert!(mol.bonds().next().is_none());
}

#[test]
fn formal_charge_sums_only_live_atom_payloads() {
    let mut mol = crate::core::MoleculeEditor::new();
    let positive = mol
        .add_atom(charged_atom("N", 3))
        .expect("atom identifier capacity");
    mol.add_atom(charged_atom("O", -1))
        .expect("atom identifier capacity");
    let deleted = mol
        .add_atom(charged_atom("Cl", -2))
        .expect("atom identifier capacity");

    assert_eq!(mol.formal_charge(), 0);
    mol.delete_atom(deleted)
        .expect("charged atom should delete");
    assert_eq!(mol.formal_charge(), 2);

    mol.atom_mut(positive)
        .expect("nitrogen should exist")
        .formal_charge = 1;
    assert_eq!(mol.formal_charge(), 0);
}

#[test]
fn atom_insertion_assigns_stable_typed_ids() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(oxygen()).expect("atom identifier capacity");

    assert_eq!(a.raw(), 0);
    assert_eq!(b.raw(), 1);
    assert_eq!(mol.atom_count(), 2);
    assert_eq!(
        mol.atom(a).expect("first atom exists").element.symbol(),
        "C"
    );
    assert_eq!(
        mol.atom(b).expect("second atom exists").element.symbol(),
        "O"
    );
    assert_eq!(mol.atom_ids().collect::<Vec<_>>(), vec![a, b]);
}

#[test]
fn bond_insertion_assigns_stable_typed_ids() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let bond = mol
        .add_bond(a, b, BondOrder::Single)
        .expect("bond should be valid");

    assert_eq!(bond.raw(), 0);
    assert_eq!(mol.bond_count(), 1);
    assert_eq!(
        mol.bond(bond).expect("bond should exist").endpoints(),
        (a, b)
    );
    assert_eq!(mol.bond_ids().collect::<Vec<_>>(), vec![bond]);
}

#[test]
fn invalid_atom_ids_are_rejected() {
    let mut mol = crate::core::MoleculeEditor::new();
    let atom = mol.add_atom(carbon()).expect("atom identifier capacity");

    assert_eq!(
        mol.atom(AtomId::new(99))
            .expect_err("missing atom should fail"),
        MoleculeError::InvalidAtomId(AtomId::new(99))
    );
    mol.delete_atom(atom).expect("atom should delete");
    assert_eq!(
        mol.atom(atom).expect_err("deleted atom should fail"),
        MoleculeError::InvalidAtomId(atom)
    );
    assert_eq!(
        mol.add_bond(atom, AtomId::new(99), BondOrder::Single)
            .expect_err("deleted endpoint should fail"),
        MoleculeError::InvalidAtomId(atom)
    );
}

#[test]
fn invalid_bond_ids_are_rejected() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let bond = mol
        .add_bond(a, b, BondOrder::Single)
        .expect("bond should be valid");

    assert_eq!(
        mol.bond(BondId::new(99))
            .expect_err("missing bond should fail"),
        MoleculeError::InvalidBondId(BondId::new(99))
    );
    mol.delete_bond(bond).expect("bond should delete");
    assert_eq!(
        mol.bond(bond).expect_err("deleted bond should fail"),
        MoleculeError::InvalidBondId(bond)
    );
    assert_eq!(
        mol.delete_bond(bond)
            .expect_err("deleting bond twice should fail"),
        MoleculeError::InvalidBondId(bond)
    );
}

#[test]
fn self_bonds_are_rejected() {
    let mut mol = crate::core::MoleculeEditor::new();
    let atom = mol.add_atom(carbon()).expect("atom identifier capacity");

    let err = mol
        .add_bond(atom, atom, BondOrder::Single)
        .expect_err("self-bond should fail");
    assert_eq!(err, MoleculeError::SelfBond(atom));
}

#[test]
fn duplicate_bond_is_rejected() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    mol.add_bond(a, b, BondOrder::Single)
        .expect("first bond should be valid");

    let err = mol
        .add_bond(a, b, BondOrder::Double)
        .expect_err("duplicate should fail");
    assert_eq!(err, MoleculeError::DuplicateBond { a, b });

    let reverse_err = mol
        .add_bond(b, a, BondOrder::Double)
        .expect_err("reverse duplicate should fail");
    assert_eq!(reverse_err, MoleculeError::DuplicateBond { a: b, b: a });
}

#[test]
fn neighbor_iteration_reports_live_adjacent_atoms() {
    let mut mol = crate::core::MoleculeEditor::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let isolated = mol.add_atom(carbon()).expect("atom identifier capacity");
    mol.add_bond(center, left, BondOrder::Single)
        .expect("left bond should be valid");
    mol.add_bond(center, right, BondOrder::Double)
        .expect("right bond should be valid");

    assert_eq!(
        sorted_atom_ids(mol.neighbors(center).expect("center exists")),
        vec![left, right]
    );
    assert_eq!(
        mol.neighbors(isolated)
            .expect("isolated atom exists")
            .collect::<Vec<_>>(),
        Vec::<AtomId>::new()
    );
    match mol.neighbors(AtomId::new(99)) {
        Ok(_) => panic!("missing atom should fail"),
        Err(err) => assert_eq!(err, MoleculeError::InvalidAtomId(AtomId::new(99))),
    };
}

#[test]
fn incident_bond_iteration_reports_live_bonds() {
    let mut mol = crate::core::MoleculeEditor::new();
    let center = mol.add_atom(carbon()).expect("atom identifier capacity");
    let left = mol.add_atom(carbon()).expect("atom identifier capacity");
    let right = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let left_bond = mol
        .add_bond(center, left, BondOrder::Single)
        .expect("left bond should be valid");
    let right_bond = mol
        .add_bond(center, right, BondOrder::Double)
        .expect("right bond should be valid");

    assert_eq!(
        sorted_bond_ids(
            mol.incident_bonds(center)
                .expect("center exists")
                .map(|(id, _)| id)
        ),
        vec![left_bond, right_bond]
    );

    mol.delete_bond(left_bond).expect("left bond should delete");
    assert_eq!(
        mol.incident_bonds(center)
            .expect("center still exists")
            .map(|(id, _)| id)
            .collect::<Vec<_>>(),
        vec![right_bond]
    );
    match mol.incident_bonds(AtomId::new(99)) {
        Ok(_) => panic!("missing atom should fail"),
        Err(err) => assert_eq!(err, MoleculeError::InvalidAtomId(AtomId::new(99))),
    };
}

#[test]
fn bond_between_finds_live_undirected_bonds() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let c = mol.add_atom(carbon()).expect("atom identifier capacity");
    let bond = mol
        .add_bond(a, b, BondOrder::Single)
        .expect("bond should be valid");

    assert_eq!(mol.bond_between(a, b).expect("atoms exist"), Some(bond));
    assert_eq!(mol.bond_between(b, a).expect("atoms exist"), Some(bond));
    assert_eq!(mol.bond_between(a, c).expect("atoms exist"), None);
}

#[test]
fn bond_deletion_preserves_remaining_ids_and_counts() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let first = mol
        .add_bond(a, b, BondOrder::Single)
        .expect("first bond should be valid");
    let second = mol
        .add_bond(b, c, BondOrder::Double)
        .expect("second bond should be valid");

    let removed = mol.delete_bond(first).expect("first bond should delete");

    assert_eq!(removed.a(), a);
    assert_eq!(mol.bond_count(), 1);
    assert_eq!(mol.bond(first), Err(MoleculeError::InvalidBondId(first)));
    assert_eq!(
        mol.bond(second).expect("second bond remains").order,
        BondOrder::Double
    );
    assert_eq!(mol.bond_ids().collect::<Vec<_>>(), vec![second]);
    assert_eq!(
        mol.neighbors(b)
            .expect("middle atom exists")
            .collect::<Vec<_>>(),
        vec![c]
    );
}

#[test]
fn atom_deletion_removes_incident_bonds_and_preserves_remaining_ids() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let c = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let first = mol
        .add_bond(a, b, BondOrder::Single)
        .expect("first bond should be valid");
    let second = mol
        .add_bond(b, c, BondOrder::Double)
        .expect("second bond should be valid");

    let removed = mol.delete_atom(b).expect("middle atom should delete");

    assert_eq!(removed.element.symbol(), "C");
    assert_eq!(mol.atom_count(), 2);
    assert_eq!(mol.bond_count(), 0);
    assert_eq!(mol.atom(b), Err(MoleculeError::InvalidAtomId(b)));
    assert_eq!(
        mol.atom(a).expect("first atom remains").element.symbol(),
        "C"
    );
    assert_eq!(
        mol.atom(c).expect("third atom remains").element.symbol(),
        "O"
    );
    assert_eq!(mol.bond(first), Err(MoleculeError::InvalidBondId(first)));
    assert_eq!(mol.bond(second), Err(MoleculeError::InvalidBondId(second)));
    assert_eq!(mol.atom_ids().collect::<Vec<_>>(), vec![a, c]);
    assert_eq!(
        mol.neighbors(a)
            .expect("first atom exists")
            .collect::<Vec<_>>(),
        Vec::<AtomId>::new()
    );
}

#[test]
fn adding_after_deletion_allocates_new_ids() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(carbon()).expect("atom identifier capacity");
    let first_bond = mol
        .add_bond(a, b, BondOrder::Single)
        .expect("bond should be valid");
    mol.delete_bond(first_bond).expect("bond should delete");
    mol.delete_atom(a).expect("atom should delete");

    let c = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let second_bond = mol
        .add_bond(b, c, BondOrder::Double)
        .expect("new bond should be valid");

    assert_eq!(c.raw(), 2);
    assert_eq!(second_bond.raw(), 1);
    assert_eq!(mol.atom_ids().collect::<Vec<_>>(), vec![b, c]);
    assert_eq!(mol.bond_ids().collect::<Vec<_>>(), vec![second_bond]);
}

#[test]
fn every_topology_mutation_invalidates_fresh_perception() {
    let mut mol = crate::core::MoleculeEditor::new();
    mark_all_fresh(&mut mol);
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    assert_all_stale(&mol);

    mark_all_fresh(&mut mol);
    let b = mol.add_atom(oxygen()).expect("atom identifier capacity");
    assert_all_stale(&mol);

    mark_all_fresh(&mut mol);
    let bond = mol
        .add_bond(a, b, BondOrder::Single)
        .expect("bond should be valid");
    assert_all_stale(&mol);

    mark_all_fresh(&mut mol);
    mol.delete_bond(bond).expect("bond should delete");
    assert_all_stale(&mol);

    mark_all_fresh(&mut mol);
    mol.delete_atom(a).expect("atom should delete");
    assert_all_stale(&mol);
}

#[test]
fn absent_perception_remains_absent_after_topology_mutation() {
    let mut mol = crate::core::MoleculeEditor::new();

    mol.add_atom(carbon()).expect("atom identifier capacity");

    assert!(!mol.perception().has_valence());
    assert!(!mol.perception().has_rings());
    assert!(!mol.perception().has_aromaticity());
    assert!(!mol.perception().has_stereo());
}

#[test]
fn properties_can_be_mutated_without_topology_changes() {
    let mut mol = crate::core::MoleculeEditor::new();
    let a = mol.add_atom(carbon()).expect("atom identifier capacity");
    let b = mol.add_atom(oxygen()).expect("atom identifier capacity");
    let bond = mol
        .add_bond(a, b, BondOrder::Single)
        .expect("bond should be valid");
    let name = PropertyKey::new("name").unwrap();
    let role = PropertyKey::new("role").unwrap();
    let source = PropertyKey::new("source").unwrap();
    mol.insert_property(
        name.clone(),
        PropertyValue::String("carbon monoxide".to_owned()),
    )
    .unwrap();
    mol.set_atom_property(
        a,
        role.clone(),
        Some(PropertyValue::String("donor".to_owned())),
    )
    .unwrap();
    mol.set_bond_property(bond, source.clone(), Some(PropertyValue::Bool(true)))
        .unwrap();

    assert_eq!(mol.atom_count(), 2);
    assert_eq!(mol.bond_count(), 1);
    assert_eq!(
        mol.properties().get(&name),
        Some(&PropertyValue::String("carbon monoxide".to_owned()))
    );
    assert_eq!(
        mol.atom_property(a, &role).unwrap(),
        Some(PropertyValue::String("donor".to_owned()))
    );
    assert_eq!(
        mol.bond_property(bond, &source).unwrap(),
        Some(PropertyValue::Bool(true))
    );
}

#[test]
fn atom_and_bond_properties_follow_stable_ids_across_tombstones() {
    let mut molecule = crate::core::MoleculeEditor::new();
    let first = molecule.add_atom(carbon()).unwrap();
    let removed = molecule.add_atom(carbon()).unwrap();
    let last = molecule.add_atom(oxygen()).unwrap();
    let removed_bond = molecule
        .add_bond(first, removed, BondOrder::Single)
        .unwrap();
    let retained_bond = molecule.add_bond(first, last, BondOrder::Double).unwrap();
    let tag = PropertyKey::new("tag").unwrap();
    molecule
        .set_atom_property(last, tag.clone(), Some(PropertyValue::Int(3)))
        .unwrap();
    molecule
        .set_bond_property(
            retained_bond,
            tag.clone(),
            Some(PropertyValue::String("retained".into())),
        )
        .unwrap();
    molecule.delete_bond(removed_bond).unwrap();
    molecule.delete_atom(removed).unwrap();

    assert!(molecule
        .set_atom_property(removed, tag.clone(), Some(PropertyValue::Int(9)))
        .is_err());
    assert!(molecule
        .set_bond_property(removed_bond, tag.clone(), Some(PropertyValue::Int(9)))
        .is_err());
    assert!(!molecule
        .atom_properties()
        .row_has_data(removed.index())
        .unwrap());
    assert!(!molecule
        .bond_properties()
        .row_has_data(removed_bond.index())
        .unwrap());

    assert_eq!(
        molecule.atom_property(last, &tag).unwrap(),
        Some(PropertyValue::Int(3))
    );
    assert_eq!(
        molecule.bond_property(retained_bond, &tag).unwrap(),
        Some(PropertyValue::String("retained".into()))
    );
}

#[test]
fn property_and_coordinate_edits_preserve_computed_state() {
    let (mut mol, atoms, bonds) = ring_molecule(
        &["C", "C", "C"],
        &[BondOrder::Single, BondOrder::Single, BondOrder::Single],
    );
    rings_api::perceive_ring_set(&mut mol).expect("ring perception should succeed");
    let _ = valence_api::perceive_valence(&mut mol, ValenceModel::RdkitLike);
    mol.begin_aromaticity(AromaticityModel::RdkitLike);
    let before = mol.perception().clone();

    mol.set_atom_property(
        atoms[0],
        PropertyKey::new("label").unwrap(),
        Some(PropertyValue::String("a".to_owned())),
    )
    .unwrap();
    mol.set_bond_property(
        bonds[0],
        PropertyKey::new("score").unwrap(),
        Some(PropertyValue::Int(1)),
    )
    .unwrap();
    mol.insert_property(
        PropertyKey::new("name").unwrap(),
        PropertyValue::String("triangle".to_owned()),
    )
    .unwrap();
    assert_eq!(mol.perception(), &before);
    assert!(mol.ring_membership().is_some());
    assert!(mol.ring_set().is_some());
}
