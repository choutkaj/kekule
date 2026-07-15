use super::*;

#[test]
fn empty_molecule_has_no_atoms_or_bonds() {
    let mol = Molecule::new();

    assert_eq!(mol.atom_count(), 0);
    assert_eq!(mol.bond_count(), 0);
    assert!(mol.atoms().next().is_none());
    assert!(mol.bonds().next().is_none());
}

#[test]
fn atom_insertion_assigns_stable_typed_ids() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(oxygen());

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
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(oxygen());
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
    let mut mol = Molecule::new();
    let atom = mol.add_atom(carbon());

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
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(oxygen());
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
    let mut mol = Molecule::new();
    let atom = mol.add_atom(carbon());

    let err = mol
        .add_bond(atom, atom, BondOrder::Single)
        .expect_err("self-bond should fail");
    assert_eq!(err, MoleculeError::SelfBond(atom));
}

#[test]
fn duplicate_bond_is_rejected() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(carbon());
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
    let mut mol = Molecule::new();
    let center = mol.add_atom(carbon());
    let left = mol.add_atom(carbon());
    let right = mol.add_atom(oxygen());
    let isolated = mol.add_atom(carbon());
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
    let mut mol = Molecule::new();
    let center = mol.add_atom(carbon());
    let left = mol.add_atom(carbon());
    let right = mol.add_atom(oxygen());
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
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(oxygen());
    let c = mol.add_atom(carbon());
    let bond = mol
        .add_bond(a, b, BondOrder::Single)
        .expect("bond should be valid");

    assert_eq!(mol.bond_between(a, b).expect("atoms exist"), Some(bond));
    assert_eq!(mol.bond_between(b, a).expect("atoms exist"), Some(bond));
    assert_eq!(mol.bond_between(a, c).expect("atoms exist"), None);
}

#[test]
fn bond_deletion_preserves_remaining_ids_and_counts() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(carbon());
    let c = mol.add_atom(oxygen());
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
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(carbon());
    let c = mol.add_atom(oxygen());
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
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(carbon());
    let first_bond = mol
        .add_bond(a, b, BondOrder::Single)
        .expect("bond should be valid");
    mol.delete_bond(first_bond).expect("bond should delete");
    mol.delete_atom(a).expect("atom should delete");

    let c = mol.add_atom(oxygen());
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
    let mut mol = Molecule::new();
    mark_all_fresh(&mut mol);
    let a = mol.add_atom(carbon());
    assert_all_stale(&mol);

    mark_all_fresh(&mut mol);
    let b = mol.add_atom(oxygen());
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
    let mut mol = Molecule::new();

    mol.add_atom(carbon());

    assert!(!mol.perception().has_valence());
    assert!(!mol.perception().has_rings());
    assert!(!mol.perception().has_aromaticity());
    assert!(!mol.perception().has_cip_descriptors());
}

#[test]
fn property_maps_can_be_mutated_without_topology_changes() {
    let mut mol = Molecule::new();
    let a = mol.add_atom(carbon());
    let b = mol.add_atom(oxygen());
    let bond = mol
        .add_bond(a, b, BondOrder::Single)
        .expect("bond should be valid");
    mol.props_mut().insert(
        "name".to_owned(),
        PropValue::String("carbon monoxide".to_owned()),
    );
    mol.atom_mut(a)
        .expect("atom exists")
        .props
        .insert("role".to_owned(), PropValue::String("donor".to_owned()));
    mol.bond_mut(bond)
        .expect("bond exists")
        .props
        .insert("source".to_owned(), PropValue::Bool(true));

    assert_eq!(mol.atom_count(), 2);
    assert_eq!(mol.bond_count(), 1);
    assert_eq!(
        mol.props().get("name"),
        Some(&PropValue::String("carbon monoxide".to_owned()))
    );
    assert_eq!(
        mol.atom(a).expect("atom exists").props.get("role"),
        Some(&PropValue::String("donor".to_owned()))
    );
    assert_eq!(
        mol.bond(bond).expect("bond exists").props.get("source"),
        Some(&PropValue::Bool(true))
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
    mol.begin_aromaticity(AromaticityProvenance::Imported);
    let before = mol.perception().clone();

    mol.atom_mut(atoms[0])
        .expect("atom exists")
        .props
        .insert("label".to_owned(), PropValue::String("a".to_owned()));
    mol.bond_mut(bonds[0])
        .expect("bond exists")
        .props
        .insert("score".to_owned(), PropValue::Int(1));
    mol.props_mut()
        .insert("name".to_owned(), PropValue::String("triangle".to_owned()));
    let mut conformer = Conformer::new();
    conformer.set_position(atoms[0], Point3::new(0.0, 0.0, 0.0));
    let conformer_id = mol.add_conformer(conformer).expect("valid conformer");
    mol.conformer_mut(conformer_id)
        .expect("conformer exists")
        .set_position(atoms[1], Point3::new(1.0, 0.0, 0.0));

    assert_eq!(mol.perception(), &before);
    assert!(mol.ring_membership().is_some());
    assert!(mol.ring_set().is_some());
}

#[test]
fn conformer_attachment_rejects_coordinates_for_non_live_atoms_transactionally() {
    let mut mol = Molecule::new();
    let deleted = mol.add_atom(carbon());
    let live = mol.add_atom(oxygen());
    mol.delete_atom(deleted).expect("delete atom");
    let mut conformer = Conformer::new();
    conformer.set_position(deleted, Point3::new(0.0, 0.0, 0.0));
    conformer.set_position(live, Point3::new(1.0, 0.0, 0.0));

    assert!(matches!(
        mol.add_conformer(conformer),
        Err(MoleculeError::InvalidAtomId(id)) if id == deleted
    ));
    assert!(mol.conformers().next().is_none());
}
