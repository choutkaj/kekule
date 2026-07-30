use super::*;

#[test]
fn happy_path_small_molecule_api_matches_architecture() {
    let mut molecule = SmallMolecule::from_smiles("c1ccccc1O").expect("phenol parses");

    let report = molecule.sanitize().expect("phenol sanitizes");
    assert_eq!(report.ring_count, Some(1));
    assert_eq!(molecule.atom_count(), molecule.graph().atom_count());
    assert_eq!(molecule.bond_count(), molecule.graph().bond_count());

    let canonical = molecule
        .to_canonical_smiles()
        .expect("canonical SMILES writes");
    assert!(!canonical.is_empty());

    let chiral = SmallMolecule::from_smiles("F[C@H](Cl)Br").expect("chiral molecule parses");
    assert_eq!(
        chiral.to_isomeric_smiles().expect("isomeric SMILES writes"),
        "F[C@H](Cl)Br"
    );
}

#[test]
fn namespaced_small_molecule_api_keeps_parsing_and_sanitization_separate() {
    let mut molecule = read_smiles("CC(=O)O").expect("acetic acid parses");
    assert!(!molecule.graph().perception().has_valence());

    perception_api::sanitize(&mut molecule).expect("acetic acid sanitizes");
    assert!(molecule.graph().perception().has_valence());

    let canonical = smiles_api::write_canonical(&molecule).expect("canonical SMILES writes");
    assert!(!canonical.is_empty());

    let chiral = read_smiles("F[C@H](Cl)Br").expect("chiral molecule parses");
    assert_eq!(
        smiles_api::write_isomeric(&chiral).expect("isomeric SMILES writes"),
        "F[C@H](Cl)Br"
    );
}
