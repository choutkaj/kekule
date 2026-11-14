use super::*;

#[test]
fn aromatic_localization_limit_is_structured_and_transactional() {
    let mut molecule = crate::molecule::MoleculeEditor::new();
    let atoms = (0..6)
        .map(|_| {
            molecule
                .add_atom(Atom::new(
                    crate::molecule::Element::from_symbol("C").unwrap(),
                ))
                .unwrap()
        })
        .collect::<Vec<_>>();
    let aromatic_bonds = (0..6)
        .map(|index| {
            molecule
                .add_bond(atoms[index], atoms[(index + 1) % 6], BondOrder::Single)
                .unwrap()
        })
        .collect::<BTreeSet<_>>();
    let component = source_aromatic_bond_components(&molecule, &aromatic_bonds)
        .into_iter()
        .next()
        .expect("imported aromatic component");
    let before = molecule.clone();

    let error =
        try_localize_aromatic_component_with_limit(&mut molecule, &component, &aromatic_bonds, 0)
            .expect_err("zero matching budget should fail structurally");

    assert!(matches!(
        error,
        NormalizationError::AromaticLocalizationLimit { limit, .. }
            if limit == LimitExceeded::new(Resource::SearchStates, 0, 1)
    ));
    assert_eq!(molecule, before);
}
