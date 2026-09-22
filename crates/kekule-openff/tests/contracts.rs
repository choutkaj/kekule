use kekule::{core::Molecule, smiles};
use kekule_openff::ForceField;

fn molecule(s: &str) -> Molecule {
    smiles::to_molecules(s).unwrap().remove(0)
}
#[test]
fn rosemary_compiles_and_last_rule_wins() {
    let ff = ForceField::rosemary().unwrap();
    let labels = ff.label_molecule(&molecule("[H]C([H])([H])O[H]")).unwrap();
    assert_eq!(labels.bonds.len(), 5);
    assert_eq!(labels.angles.len(), 7);
    assert_eq!(labels.constraints.len(), 4);
    assert_eq!(labels.vdw.len(), 6);
    let xml = include_str!("../data/rosemary.offxml").replace("</Bonds>",
        "<Bond id=\"override\" smirks=\"[*:1]~[*:2]\" length=\"1 * angstrom\" k=\"1 * kilocalorie_per_mole * angstrom ** -2\"/></Bonds>");
    let overridden = ForceField::from_offxml(&xml)
        .unwrap()
        .label_molecule(&molecule("[H]C([H])([H])O[H]"))
        .unwrap();
    assert_eq!(overridden.bonds.len(), 5);
    assert!(overridden.bonds.values().all(|p| p.id == "override"));
    assert!(ff
        .label_molecule(&molecule("CO"))
        .unwrap_err()
        .to_string()
        .contains("implicit hydrogens"));
}
#[test]
fn unsupported_functional_forms_fail() {
    let xml = include_str!("../data/rosemary.offxml");
    for changed in [
        xml.replace("OEAroModel_MDL", "other"),
        xml.replace("Lennard-Jones-12-6", "Buckingham"),
        xml.replace("<Author>", "<VirtualSites/><Author>"),
        xml.replace("angstrom ** 1", "kelvin ** 1"),
        xml.replace("<Bonds version", "<Bonds unknown=\"true\" version"),
        xml.replace("<NAGLCharges version=\"0.3\"", "<NAGLCharges version=\"999\""),
        xml.replace("phase1=\"180.0 * degree ** 1\"", "phase1=\"30.0 * degree ** 1\""),
        xml.replace("</Bonds>", "<Bond id=\"bad\" smirks=\"[*:1]~[*:3]\" length=\"1 * angstrom\" k=\"1 * kilocalorie_per_mole * angstrom ** -2\"/></Bonds>"),
    ] {
        assert!(ForceField::from_offxml(&changed).is_err());
    }
}

#[test]
fn all_reference_rule_labels_match_in_atom_map_space() {
    use std::collections::BTreeSet;
    let report: serde_json::Value = serde_json::from_reader(flate2::read::GzDecoder::new(
        &include_bytes!("fixtures/audit.json.gz")[..],
    ))
    .unwrap();
    let ff = ForceField::rosemary().unwrap();
    for r in report["records"].as_array().unwrap() {
        let m = molecule(r["mapped_smiles"].as_str().unwrap());
        let before = m.clone();
        let labels = ff.label_molecule(&m).unwrap();
        assert_eq!(m, before);
        let check = |handler: &str, observed: Vec<(Vec<u32>, String, String)>| {
            let canonical = |mut atoms: Vec<u32>| {
                if handler == "ImproperTorsions" {
                    let mut outer = [atoms[0], atoms[2], atoms[3]];
                    outer.sort();
                    atoms = vec![outer[0], atoms[1], outer[1], outer[2]];
                } else {
                    let reverse = atoms.iter().rev().copied().collect::<Vec<_>>();
                    if reverse < atoms {
                        atoms = reverse;
                    }
                }
                atoms
            };
            let expected = r["assigned_parameters"][handler]
                .as_array()
                .unwrap()
                .iter()
                .map(|p| {
                    let maps = p["atoms"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .map(|i| {
                            r["charge_atom_maps"][i.as_u64().unwrap() as usize]
                                .as_u64()
                                .unwrap() as u32
                        })
                        .collect();
                    (
                        canonical(maps),
                        p["id"].as_str().unwrap().to_owned(),
                        p["smirks"].as_str().unwrap().to_owned(),
                    )
                })
                .collect::<BTreeSet<_>>();
            let observed = observed
                .into_iter()
                .map(|(a, id, s)| (canonical(a), id, s))
                .collect::<BTreeSet<_>>();
            assert_eq!(observed, expected, "{} {handler}", r["input"]["id"]);
        };
        macro_rules! check {
            ($handler:literal,$field:ident) => {
                check(
                    $handler,
                    labels
                        .$field
                        .iter()
                        .map(|(a, p)| {
                            (
                                a.iter()
                                    .map(|&a| m.atom(a).unwrap().atom_map.unwrap())
                                    .collect(),
                                p.id.clone(),
                                p.smirks.clone(),
                            )
                        })
                        .collect(),
                );
            };
        }
        check!("Bonds", bonds);
        check!("Angles", angles);
        check!("ProperTorsions", proper_torsions);
        check!("ImproperTorsions", improper_torsions);
        check!("Constraints", constraints);
        check!("vdW", vdw);
    }
}

#[test]
#[ignore = "requires the separately exported, checksum-pinned Ash model; set KEKULE_OPENFF_MODEL"]
fn native_charges_match_external_reference() {
    let model = kekule_openff::NaglModel::load(
        std::env::var("KEKULE_OPENFF_MODEL")
            .expect("set KEKULE_OPENFF_MODEL to the exported bundle"),
    )
    .unwrap();
    assert_eq!(model.lookup_entry_count(), 13944);
    // Isotope masses do not change the reference toolkit's charge lookup.
    let mut isotope = molecule("[2H]OC");
    isotope.perceive().unwrap();
    isotope.add_hydrogens().unwrap();
    let mut ordinary = molecule("[H]OC");
    ordinary.perceive().unwrap();
    ordinary.add_hydrogens().unwrap();
    let q = model.assign_charges(&isotope).unwrap();
    let expected = model.assign_charges(&ordinary).unwrap();
    assert_eq!(q.source, expected.source);
    assert_eq!(q.charges, expected.charges);
    let mut large = molecule(&"C".repeat(350));
    large.perceive().unwrap();
    large.add_hydrogens().unwrap();
    assert!(large.atom_count() > 1024);
    let assignment = model.assign_charges(&large).unwrap();
    assert!(matches!(
        assignment.source,
        kekule_openff::ChargeSource::Inference { .. }
    ));
    assert_eq!(
        assignment.charges,
        model.infer_charges(&large).unwrap().charges
    );
    assert!(assignment.charges.value().iter().sum::<f64>().abs() < 1e-10);
    let oversized = molecule(&"C".repeat(4097));
    for error in [
        model.assign_charges(&oversized).unwrap_err(),
        model.infer_charges(&oversized).unwrap_err(),
        model.atom_features(&oversized).unwrap_err(),
    ] {
        assert!(error.to_string().contains("4096"));
    }
    assert!(matches!(
        q.source,
        kekule_openff::ChargeSource::Lookup { .. }
    ));
    let report: serde_json::Value = serde_json::from_reader(flate2::read::GzDecoder::new(
        &include_bytes!("fixtures/audit.json.gz")[..],
    ))
    .unwrap();
    let ff = ForceField::rosemary().unwrap();
    for r in report["records"].as_array().unwrap() {
        let m = molecule(r["mapped_smiles"].as_str().unwrap());
        let p = ff.parameterize_molecule(&m, &model).unwrap();
        let expected = &r["system"]["charges"];
        for (i, (_, a)) in m.atoms().enumerate() {
            let q = expected[a.atom_map.unwrap() as usize - 1]["charge"]
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
            ff.parameterize_molecule(m, &model)
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
                    assert_eq!(ff.parameterize_molecule(m, &model).unwrap().charges(), q);
                    assert_eq!(*m, before);
                }
            });
        }
    });
}
#[test]
fn unassigned_chemistry_fails() {
    assert!(ForceField::rosemary()
        .unwrap()
        .label_molecule(&molecule("[Xe]"))
        .is_err());
}
