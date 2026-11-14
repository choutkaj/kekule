fn explicit_hydrogens(smiles: &str) -> kekule::molecule::Molecule {
    let mut molecule = kekule::smiles::to_molecules(smiles).unwrap().remove(0);
    molecule.perceive().unwrap();
    molecule.add_hydrogens(Default::default()).unwrap();
    molecule
}

// Regression: features had a 4096-atom cap, and resonance work grew with
// the whole molecule for every conjugated fragment.
#[test]
fn resonance_averages_scale_past_former_size_cap() {
    // Polyglycine: one amide resonance fragment per residue.
    let molecule = explicit_hydrogens(&format!("N{}CC(=O)O", "CC(=O)N".repeat(800)));
    assert!(molecule.atom_count() > 4096);
    let averages = super::average_charges(&molecule).unwrap();
    assert_eq!(averages.len(), molecule.atom_count());
    // Neutral amides have one dominant form; nothing becomes charged.
    assert!(averages.iter().all(|q| *q == 0.0));
}

// Regression: normalization stopped after 200 applications of one rule,
// as upstream NAGL does, so molecules with more matching groups failed.
#[test]
fn normalization_converges_past_two_hundred_matching_groups() {
    // Each zwitterionic imidic group normalizes to a neutral amide.
    let molecule = explicit_hydrogens(&format!("C{}", "C(C([O-])=[NH2+])".repeat(250)));
    assert_eq!(
        molecule
            .atoms()
            .filter(|(_, a)| a.formal_charge() != 0)
            .count(),
        500
    );
    let normalized = super::normalize(&molecule).unwrap();
    assert!(normalized.atoms().all(|(_, a)| a.formal_charge() == 0));
    assert_eq!(normalized.atom_count(), molecule.atom_count());
}

#[test]
fn complete_ash_features_match_pinned_external_reference() {
    let mut count = 0;
    for r in crate::reference_records() {
        if r["features"]["status"] != "ok" {
            continue;
        }
        let m = kekule::smiles::to_molecules(r["mapped_smiles"].as_str().unwrap())
            .unwrap()
            .remove(0);
        let m = crate::explicit(&m).unwrap();
        let specs = r["features"]["values"]
            .as_array()
            .unwrap()
            .iter()
            .map(|f| serde_json::from_value(f["config"].clone()).unwrap())
            .collect::<Vec<_>>();
        let actual = super::features(&m, &specs).unwrap();
        for (i, (_, atom)) in m.atoms().enumerate() {
            let index = atom.atom_map().unwrap().get() as usize - 1;
            let expected = r["features"]["values"]
                .as_array()
                .unwrap()
                .iter()
                .flat_map(|f| f["values"][index].as_array().unwrap())
                .map(|v| v.as_f64().unwrap())
                .collect::<Vec<_>>();
            assert_eq!(expected.len(), 22);
            for (j, &expected) in expected.iter().enumerate() {
                assert!(
                    (f64::from(actual[(i, j)]) - expected).abs() <= 1e-6,
                    "{} atom {} column {j}",
                    r["input"]["id"],
                    index + 1
                );
            }
        }
        count += 1;
    }
    assert_eq!(count, 22);
}
