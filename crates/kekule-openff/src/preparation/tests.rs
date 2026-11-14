// Expected graphs are RDKit's sanitized forms, which upstream identifies.
fn same_graph(prepared: &kekule::molecule::Molecule, expected: &str) -> bool {
    let mut expected = kekule::smiles::to_molecules(expected).unwrap().remove(0);
    expected.perceive().unwrap();
    crate::identity::same_graph(prepared, &expected).unwrap()
}

#[test]
fn neutral_valence_five_nitrogen_matches_charge_separated_input() {
    for (s, expected) in [
        ("[H:4][N:2](=[O:1])=[O:3]", "[H][N+](=O)[O-]"),
        ("[H:4][N:2](#[N:3])[F:1]", "[H][N+](=[N-])F"),
    ] {
        let input = kekule::smiles::to_molecules(s).unwrap().remove(0);
        let prepared = crate::explicit(&input).unwrap();
        assert!(same_graph(&prepared, expected), "{s}");
        assert!(input.atoms().all(|(_, a)| a.formal_charge() == 0));
    }
}
#[test]
fn ash_phosphorus_entry_matches_rdkit_without_mutating_input() {
    let input = kekule::smiles::to_molecules("[H:5][N:3]=[P:2](=[O:1])[H:4]")
        .unwrap()
        .remove(0);
    let normalized = crate::explicit(&input).unwrap();
    assert!(same_graph(&normalized, "[H]N=[P+]([H])[O-]"));
    assert!(input.atoms().all(|(_, a)| a.formal_charge() == 0));
    assert_eq!(
        input.atom_ids().collect::<Vec<_>>(),
        normalized.atom_ids().collect::<Vec<_>>()
    );
    assert_eq!(input.formal_charge(), normalized.formal_charge());
}
