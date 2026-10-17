//! Native NAGL charges, and the lookup and network paths behind them, against
//! the external reference audit.

use kekule::{molecule::Molecule, smiles};

use crate::ForceField;

fn molecule(s: &str) -> Molecule {
    smiles::to_molecules(s).unwrap().remove(0)
}

#[test]
fn native_charges_match_external_reference() {
    let model = crate::NaglModel::ash().unwrap();
    assert_eq!(model.lookup_entry_count(), 13944);
    // Isotope masses do not change the reference toolkit's charge lookup.
    let mut isotope = molecule("[2H]OC");
    isotope.perceive().unwrap();
    isotope.add_hydrogens(Default::default()).unwrap();
    let mut ordinary = molecule("[H]OC");
    ordinary.perceive().unwrap();
    ordinary.add_hydrogens(Default::default()).unwrap();
    let q = model.assign_charges(&isotope).unwrap();
    let expected = model.assign_charges(&ordinary).unwrap();
    assert_eq!(q.source, expected.source);
    assert_eq!(q.charges, expected.charges);
    let mut large = molecule(&"C".repeat(350));
    large.perceive().unwrap();
    large.add_hydrogens(Default::default()).unwrap();
    assert!(large.atom_count() > 1024);
    let assignment = model.assign_charges(&large).unwrap();
    assert!(matches!(
        assignment.source,
        crate::ChargeSource::Inference { .. }
    ));
    assert_eq!(
        assignment.charges,
        model.infer_charges(&large).unwrap().charges
    );
    assert!(assignment.charges.value().iter().sum::<f64>().abs() < 1e-10);
    // Regression: inference had a 4096-atom cap.
    let mut oversized = molecule(&"C".repeat(1500));
    oversized.perceive().unwrap();
    oversized.add_hydrogens(Default::default()).unwrap();
    assert!(oversized.atom_count() > 4096);
    let assignment = model.assign_charges(&oversized).unwrap();
    assert_eq!(assignment.charges.value().len(), oversized.atom_count());
    assert!(assignment.charges.value().iter().sum::<f64>().abs() < 1e-10);
    assert_eq!(
        model.atom_features(&oversized).unwrap().len(),
        oversized.atom_count()
    );
    assert!(matches!(q.source, crate::ChargeSource::Lookup { .. }));
    let report: serde_json::Value = serde_json::from_reader(flate2::read::GzDecoder::new(
        &include_bytes!("../../tests/fixtures/audit.json.gz")[..],
    ))
    .unwrap();
    let ff = ForceField::rosemary().unwrap();
    for r in report["records"].as_array().unwrap() {
        let m = molecule(r["mapped_smiles"].as_str().unwrap());
        let p = ff.parameterize_molecule(m.clone(), &model).unwrap();
        let expected = &r["system"]["charges"];
        for (i, (_, a)) in m.atoms().enumerate() {
            let q = expected[a.atom_map().unwrap().get() as usize - 1]["charge"]
                .as_f64()
                .unwrap();
            assert!(
                (p.charges().value()[i] - q).abs() <= 5e-5,
                "{}",
                r["input"]["id"]
            );
        }
        assert!((p.charges().value().iter().sum::<f64>() - m.formal_charge() as f64).abs() < 1e-10);
    }
    // Sharing the immutable engine/model must produce the same lookup and
    // inference results as serial execution, without mutating the inputs.
    let cases: Vec<_> = report["records"]
        .as_array()
        .unwrap()
        .iter()
        .take(6)
        .map(|r| molecule(r["mapped_smiles"].as_str().unwrap()))
        .collect();
    let expected: Vec<_> = cases
        .iter()
        .map(|m| {
            ff.parameterize_molecule(m.clone(), &model)
                .unwrap()
                .charges()
                .clone()
        })
        .collect();
    std::thread::scope(|scope| {
        for _ in 0..4 {
            scope.spawn(|| {
                for (m, q) in cases.iter().zip(&expected) {
                    let before = m.clone();
                    assert_eq!(
                        ff.parameterize_molecule(m.clone(), &model)
                            .unwrap()
                            .charges(),
                        q
                    );
                    assert_eq!(*m, before);
                }
            });
        }
    });
}
