use super::*;

fn prepared(smiles: &str) -> Molecule {
    crate::explicit(&kekule::smiles::to_molecules(smiles).unwrap().remove(0)).unwrap()
}

#[test]
fn selection_requires_the_exact_entry_molecule() {
    let index = LookupIndex::new([
        (prepared("[H:4][N:2]=[N:1][F:3]"), 7),
        (prepared("[H]O[H]"), 3),
    ])
    .unwrap();
    let select = |s: &str| index.select(&prepared(s)).unwrap();
    assert_eq!(select("FN=N[H]"), Some(7));
    assert_eq!(select("[H]O[H]"), Some(3));
    // Another charge placement, stereo, or constitution is not the entry.
    assert_eq!(select("[H][N-]N=[F+]"), None);
    assert_eq!(select("[H]/N=N/F"), None);
    assert_eq!(select("[H]N([H])F"), None);
}

#[test]
fn mapping_relaxes_charges_and_orders_then_stereo() {
    let entry = prepared("[H:4][N:2]=[N:1][F:3]");
    let query = prepared("[F+:3]=[N:1][N-:2][H:4]");
    let map = mapping(&query, &entry).unwrap();
    let elements = |m: &Molecule, ids: Vec<AtomId>| {
        ids.into_iter()
            .map(|a| m.atom(a).unwrap().element())
            .collect::<Vec<_>>()
    };
    assert_eq!(
        elements(&query, query.atom_ids().collect()),
        elements(&entry, map)
    );
    let other = prepared("[H]N([H])F");
    assert_eq!(
        mapping(&other, &entry).unwrap_err().kind(),
        ErrorKind::Identity
    );
}
