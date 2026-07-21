use super::*;

#[test]
fn canonical_ranking_groups_symmetric_atoms() {
    let mut molecule = read_smiles("CC(C)C").expect("isobutane parses");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("isobutane sanitizes");

    let ranking = canon::atom_ranking(molecule.graph());

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
    let mut first = SmallMolecule::new();
    let first_terminal_a = first.graph_mut().add_atom(carbon());
    let first_center = first.graph_mut().add_atom(carbon());
    let first_terminal_b = first.graph_mut().add_atom(carbon());
    first
        .graph_mut()
        .add_bond(first_terminal_a, first_center, BondOrder::Single)
        .expect("bond should be valid");
    first
        .graph_mut()
        .add_bond(first_center, first_terminal_b, BondOrder::Single)
        .expect("bond should be valid");
    perception_api::sanitize_with_options(&mut first, SanitizeOptions::default())
        .expect("propane sanitizes");

    let mut second = SmallMolecule::new();
    let second_center = second.graph_mut().add_atom(carbon());
    let second_terminal_a = second.graph_mut().add_atom(carbon());
    let second_terminal_b = second.graph_mut().add_atom(carbon());
    second
        .graph_mut()
        .add_bond(second_center, second_terminal_a, BondOrder::Single)
        .expect("bond should be valid");
    second
        .graph_mut()
        .add_bond(second_center, second_terminal_b, BondOrder::Single)
        .expect("bond should be valid");
    perception_api::sanitize_with_options(&mut second, SanitizeOptions::default())
        .expect("propane sanitizes");

    let first_ranking = canon::atom_ranking(first.graph());
    let second_ranking = canon::atom_ranking(second.graph());

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
    let mut molecule = read_smiles("[13CH3:7][CH3]").expect("mapped isotope molecule parses");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("mapped isotope molecule sanitizes");

    let ranking = canon::atom_ranking(molecule.graph());

    assert_ne!(
        ranking.rank_of(AtomId::new(0)),
        ranking.rank_of(AtomId::new(1))
    );
}

#[test]
fn canonical_ranking_ignores_kekule_choice_for_perceived_aromatic_bonds() {
    let mut molecule = read_smiles("c1ccc2ccccc2c1").expect("naphthalene parses");
    perception_api::sanitize_with_options(&mut molecule, SanitizeOptions::default())
        .expect("naphthalene sanitizes");

    let ranking = canon::atom_ranking(molecule.graph());

    assert_eq!(ranking.rank_count(), 3);
}
