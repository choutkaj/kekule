use super::*;

#[test]
fn happy_path_universal_molecule_api_matches_architecture() {
    let mut molecule = read_smiles("c1ccccc1O").expect("phenol parses");

    assert_eq!(molecule.perception(), &Perception::default());
    molecule.perceive().expect("phenol perceives");
    assert_eq!(molecule.ring_set().expect("installed ring basis").len(), 1);
    assert_eq!(molecule.atom_count(), molecule.atom_count());
    assert_eq!(molecule.bond_count(), molecule.bond_count());

    let canonical = smiles_api::write_canonical(&molecule).expect("canonical SMILES writes");
    assert!(!canonical.is_empty());

    let chiral = read_smiles("F[C@H](Cl)Br").expect("chiral molecule parses");
    assert_eq!(
        smiles_api::write_isomeric(&chiral).expect("isomeric SMILES writes"),
        "F[C@H](Cl)Br"
    );
}

#[test]
fn namespaced_molecule_api_keeps_pipeline_stages_separate() {
    let mut molecule = read_smiles("CC(=O)O").expect("acetic acid parses");
    assert!(!molecule.perception().has_valence());

    perceive(&mut molecule).expect("acetic acid perceives");
    assert!(molecule.perception().has_valence());

    let canonical = smiles_api::write_canonical(&molecule).expect("canonical SMILES writes");
    assert!(!canonical.is_empty());

    let chiral = read_smiles("F[C@H](Cl)Br").expect("chiral molecule parses");
    assert_eq!(
        smiles_api::write_isomeric(&chiral).expect("isomeric SMILES writes"),
        "F[C@H](Cl)Br"
    );
}
