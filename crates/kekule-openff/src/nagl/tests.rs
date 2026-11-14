use super::*;

#[test]
fn lookup_preparation_drops_isotopes_without_mutating_input() {
    // PubChem CID 143783 from the locked pubchem-100k corpus. The OpenFF
    // Toolkit's charge representation retains no isotopic masses.
    let mut input = kekule::smiles::to_molecules("[2H]C([2H])([2H])C([2H])([2H])SCC")
        .unwrap()
        .remove(0);
    input.perceive().unwrap();
    input.add_hydrogens(Default::default()).unwrap();
    let before = input.clone();
    let prepared = lookup_molecule(&input).unwrap();
    assert!(prepared.atoms().all(|(_, a)| a.isotope().is_none()));
    assert_eq!(input, before);
    assert_eq!(
        prepared.atom_ids().collect::<Vec<_>>(),
        input.atom_ids().collect::<Vec<_>>()
    );
}

#[cfg(feature = "ash")]
#[test]
fn lookup_selects_exact_entry_molecules_ignoring_isotopes() {
    // Heavy water selects the same Ash entry as water, as upstream does.
    let model = NaglModel::ash().unwrap();
    let select = |smiles: &str| {
        let m = kekule::smiles::to_molecules(smiles).unwrap().remove(0);
        model.lookup_key(&m).unwrap().map(str::to_owned)
    };
    assert_eq!(select("[H]O[H]").as_deref(), Some("InChI=1/H2O/h1H2"));
    assert_eq!(select("[2H]O[2H]"), select("[H]O[H]"));
    assert_eq!(
        select("[H]C([H])([H])[H]").as_deref(),
        Some("InChI=1/CH4/h1H4")
    );
}

#[test]
fn invalid_charge_outputs_fail_without_publishing() {
    let source = ChargeSource::Inference {
        model: ModelIdentity::new("test".into(), "0".repeat(64)).unwrap(),
    };
    for values in [vec![], vec![f64::NAN], vec![f64::INFINITY]] {
        assert!(normalize_charges(values, 0, source.clone()).is_err());
    }
    let q = normalize_charges(vec![0.2, 0.3], 1, source).unwrap();
    assert!((q.charges.value().iter().sum::<f64>() - 1.0).abs() < 1e-14);
}
