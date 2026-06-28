use super::*;

#[test]
fn canonical_ranking_groups_symmetric_atoms() {
    let mut molecule = read_smiles_str("CC(C)C", SmilesParseOptions).expect("isobutane parses");
    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("isobutane sanitizes");

    let ranking = canonical_atom_ranking(&molecule.mol);

    assert_eq!(ranking.rank_count(), 2);
    assert_eq!(
        ranking.rank_of(AtomId::new(0)),
        ranking.rank_of(AtomId::new(2))
    );
    assert_eq!(
        ranking.rank_of(AtomId::new(0)),
        ranking.rank_of(AtomId::new(3))
    );
    assert_ne!(
        ranking.rank_of(AtomId::new(0)),
        ranking.rank_of(AtomId::new(1))
    );
}

#[test]
fn canonical_ranking_is_stable_across_atom_order_for_path_roles() {
    let mut first = SmallMolecule {
        mol: Molecule::new(),
    };
    let first_terminal_a = first.mol.add_atom(carbon());
    let first_center = first.mol.add_atom(carbon());
    let first_terminal_b = first.mol.add_atom(carbon());
    first
        .mol
        .add_bond(first_terminal_a, first_center, BondOrder::Single)
        .expect("bond should be valid");
    first
        .mol
        .add_bond(first_center, first_terminal_b, BondOrder::Single)
        .expect("bond should be valid");
    sanitize_small_molecule(&mut first, SanitizeOptions::default()).expect("propane sanitizes");

    let mut second = SmallMolecule {
        mol: Molecule::new(),
    };
    let second_center = second.mol.add_atom(carbon());
    let second_terminal_a = second.mol.add_atom(carbon());
    let second_terminal_b = second.mol.add_atom(carbon());
    second
        .mol
        .add_bond(second_center, second_terminal_a, BondOrder::Single)
        .expect("bond should be valid");
    second
        .mol
        .add_bond(second_center, second_terminal_b, BondOrder::Single)
        .expect("bond should be valid");
    sanitize_small_molecule(&mut second, SanitizeOptions::default()).expect("propane sanitizes");

    let first_ranking = canonical_atom_ranking(&first.mol);
    let second_ranking = canonical_atom_ranking(&second.mol);

    assert_eq!(
        first_ranking.rank_of(first_center),
        second_ranking.rank_of(second_center)
    );
    assert_eq!(
        first_ranking.rank_of(first_terminal_a),
        second_ranking.rank_of(second_terminal_a)
    );
    assert_eq!(
        first_ranking.rank_of(first_terminal_b),
        second_ranking.rank_of(second_terminal_b)
    );
}

#[test]
fn canonical_ranking_uses_isotope_hydrogens_and_atom_maps() {
    let mut molecule = read_smiles_str("[13CH3:7][CH3]", SmilesParseOptions)
        .expect("mapped isotope molecule parses");
    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())
        .expect("mapped isotope molecule sanitizes");

    let ranking = canonical_atom_ranking(&molecule.mol);

    assert_ne!(
        ranking.rank_of(AtomId::new(0)),
        ranking.rank_of(AtomId::new(1))
    );
}
