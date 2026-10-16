use kekule::{molecule::Molecule, smiles};
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
    let error = ff.label_molecule(&molecule("CO")).unwrap_err();
    assert_eq!(error.kind(), kekule_openff::ErrorKind::UnsupportedMolecule);
    assert!(error.to_string().contains("implicit hydrogens"));
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
fn unassigned_chemistry_fails() {
    assert!(ForceField::rosemary()
        .unwrap()
        .label_molecule(&molecule("[Xe]"))
        .is_err());
}
