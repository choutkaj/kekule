use super::*;

#[test]
fn happy_path_small_molecule_api_matches_architecture() {
    let mut molecule = SmallMolecule::from_smiles("c1ccccc1O").expect("phenol parses");

    let normalization = molecule.normalize().expect("phenol normalizes");
    assert!(normalization.warnings.is_empty());
    molecule.perceive().expect("phenol perceives");
    assert_eq!(
        molecule
            .graph()
            .ring_set()
            .expect("installed ring basis")
            .len(),
        1
    );
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
fn namespaced_small_molecule_api_keeps_pipeline_stages_separate() {
    let mut molecule = read_smiles("CC(=O)O").expect("acetic acid parses");
    assert!(!molecule.graph().perception().has_valence());

    normalize_and_perceive(&mut molecule).expect("acetic acid normalizes_and_perceives");
    assert!(molecule.graph().perception().has_valence());

    let canonical = smiles_api::write_canonical(&molecule).expect("canonical SMILES writes");
    assert!(!canonical.is_empty());

    let chiral = read_smiles("F[C@H](Cl)Br").expect("chiral molecule parses");
    assert_eq!(
        smiles_api::write_isomeric(&chiral).expect("isomeric SMILES writes"),
        "F[C@H](Cl)Br"
    );
}
