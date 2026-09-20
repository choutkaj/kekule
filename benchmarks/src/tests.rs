use crate::{
    compare::{differences, normalize_benchmark_for_comparison_in_place},
    dataset::{Dataset, Input},
    features::evaluate,
};
use kekule::{molfile, smiles};
use serde_json::{json, Value};
use std::{collections::BTreeSet, fs, path::Path};
fn output(feature: &str, path: &str, text: &str) -> Value {
    evaluate(
        feature,
        &Input {
            path: path.into(),
            text: text.into(),
        },
    )
    .unwrap()
}
fn compare(feature: &str, mut a: Value, mut b: Value) -> bool {
    normalize_benchmark_for_comparison_in_place(feature, &mut a);
    normalize_benchmark_for_comparison_in_place(feature, &mut b);
    differences(feature, &a, &b).exact()
}

#[test]
fn cip_sequence_descriptors_use_reference_spelling_without_losing_pseudoasymmetry() {
    for (source, label, ordinary) in [
        (r"C/C=C/1\CC[C@H](C)CC1", "z", "Z"),
        (r"C/C=C\1/CC[C@H](C)CC1", "e", "E"),
    ] {
        let value = output("stereo.cip", "input.smi", source);
        assert_eq!(
            value["records"][0]["bond_descriptors"],
            json!([{"begin_atom_index":1,"end_atom_index":2,"descriptor":label}])
        );
        let mut wrong = value.clone();
        wrong["records"][0]["bond_descriptors"][0]["descriptor"] = json!(ordinary);
        assert!(!compare("stereo.cip", value, wrong));
    }
}

#[test]
fn same_local_atom_environments_do_not_hide_different_graphs() {
    let a = output("io.smiles.parse", "input.smi", "C1CCCCC1");
    let b = output("io.smiles.parse", "input.smi", "C1CC1.C1CC1");
    assert!(!compare("io.smiles.parse", a, b));
}
#[test]
fn aromatic_observations_preserve_triple_bond_order() {
    for source in ["C1=CC=CC#C1", "C1=C(C(=CC#C1)Cl)Cl"] {
        let value = output("io.smiles.parse", "input.smi", source);
        let bonds = value["records"][0]["components"][0]["bonds"]
            .as_array()
            .unwrap();
        let triple = bonds
            .iter()
            .filter(|bond| bond["bond_type"] == "TRIPLE")
            .collect::<Vec<_>>();
        assert_eq!(triple.len(), 1, "{source}");
        assert_eq!(triple[0]["is_aromatic"], true, "{source}");
        assert_eq!(
            bonds
                .iter()
                .filter(|bond| bond["bond_type"] == "AROMATIC")
                .count(),
            5
        );
        let mut wrong_order = value.clone();
        for bond in wrong_order["records"][0]["components"][0]["bonds"]
            .as_array_mut()
            .unwrap()
        {
            if bond["bond_type"] == "TRIPLE" {
                bond["bond_type"] = json!("AROMATIC");
            }
        }
        assert!(!compare("io.smiles.parse", value, wrong_order));
    }
}

#[test]
fn query_records_separate_the_query_from_the_source_title() {
    let value = output("query.smarts", "input.smi", "[#6]-[#8]  external query");
    assert_eq!(value["records"][0]["smarts"], "[#6]-[#8]");
    assert_eq!(value["records"][0]["title"], "external query");
    assert_eq!(value["records"][0]["atom_count"], 2);
}

#[test]
fn smiles_records_use_library_record_names_and_preserve_record_order() {
    let value = output(
        "io.smiles.parse",
        "input.smi",
        "# source records\n\t[Na+].[Cl-] ||  sůl | sample \n\n CC\tsecond\tname\n",
    );
    let records = value["records"].as_array().unwrap();
    assert_eq!(records.len(), 2);
    assert_eq!(records[0]["record_index"], 0);
    assert_eq!(records[0]["title"], "sůl | sample");
    assert_eq!(records[0]["components"].as_array().unwrap().len(), 2);
    assert_eq!(records[1]["record_index"], 1);
    assert_eq!(records[1]["title"], "second\tname");
}

#[test]
fn cxsmiles_supported_chemistry_remains_asserted_in_benchmark_observations() {
    let value = output(
        "io.smiles.parse",
        "input.smi",
        "F[C@H](Cl)Br |&1:1,r| grouped",
    );
    assert_eq!(value["records"][0]["title"], "grouped");
    assert_eq!(
        value["records"][0]["components"][0]["groups"],
        json!([
            {"kind":"and", "members":[{"type":"tetrahedral","focus":[1]}]}
        ])
    );
    let relative = output("io.smiles.parse", "input.smi", "F[C@H](Cl)Br |r|");
    assert_eq!(
        relative["records"][0]["components"][0]["groups"][0]["kind"],
        "relative"
    );
    let radical = output("io.smiles.parse", "input.smi", "[CH2] |^3:0|");
    let atom = &radical["records"][0]["components"][0]["atoms"][0];
    assert_eq!(atom["radical_electrons"], 2);
    assert_eq!(atom["spin_multiplicity"], 1);
    let mut omitted = value.clone();
    omitted["records"][0]["components"][0]["groups"] = json!([]);
    assert!(!compare("io.smiles.parse", value, omitted));
}

#[test]
fn cxsmiles_is_never_silently_reduced_to_plain_smiles() {
    for feature in [
        "io.smiles.parse",
        "stereo.representation",
        "io.smiles.isomeric",
        "descriptor.molecular",
    ] {
        let input = Input {
            path: "input.smi".into(),
            text: "F[C@H](Cl)Br |unsupported:1| name".into(),
        };
        assert!(evaluate(feature, &input)
            .unwrap_err()
            .to_string()
            .contains("CXSMILES"));
    }
}

#[test]
fn observation_schema_requires_measured_nullable_fields_and_rejects_extras() {
    let original = output("io.smiles.parse", "input.smi", "CC");
    crate::observation::validate("io.smiles.parse", &original).unwrap();
    for field in [
        "isotope",
        "atom_map",
        "spin_multiplicity",
        "radical_electrons",
    ] {
        let mut value = original.clone();
        value["records"][0]["components"][0]["atoms"][0]
            .as_object_mut()
            .unwrap()
            .remove(field);
        assert!(
            crate::observation::validate("io.smiles.parse", &value).is_err(),
            "{field}"
        );
    }
    let mut value = original;
    value["records"][0]["unexpected"] = json!(true);
    assert!(crate::observation::validate("io.smiles.parse", &value).is_err());
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(crate::observation::finite(value).is_err());
    }
}
#[test]
fn radical_observations_assert_electron_count_and_source_spin_separately() {
    use kekule::core::{AtomRadical, HydrogenDeclaration};
    for spin in [None, Some(1), Some(3)] {
        let molecule = smiles::to_molecules("[CH2]").unwrap().pop().unwrap();
        let id = molecule.atom_ids().next().unwrap();
        let mut editor = molecule.into_editor();
        {
            let mut atom = editor.atom_mut(id).unwrap();
            atom.hydrogens = HydrogenDeclaration::Fixed(2);
            atom.radical = AtomRadical::new(2, spin);
        }
        let molecule = editor.finish().unwrap();
        let value = if spin.is_none() {
            output("io.smiles.parse", "input.smi", "[CH2]")
        } else {
            output(
                "io.mol.parse",
                "input.mol",
                &molfile::write_v3000(&molecule).unwrap(),
            )
        };
        let atom = &value["records"][0]["components"][0]["atoms"][0];
        assert_eq!(atom["radical_electrons"], 2);
        assert_eq!(atom["spin_multiplicity"], json!(spin));
        for (field, wrong) in [
            ("radical_electrons", json!(1)),
            ("spin_multiplicity", json!(2)),
        ] {
            let mut changed = value.clone();
            changed["records"][0]["components"][0]["atoms"][0][field] = wrong;
            assert!(!compare("io.mol.parse", value.clone(), changed));
        }
        let mut invalid = value;
        invalid["records"][0]["components"][0]["atoms"][0]["spin_multiplicity"] = json!(0);
        assert!(crate::observation::validate("io.mol.parse", &invalid).is_err());
    }
}

#[test]
fn stereo_and_coordinate_units_match_the_reference_conventions() {
    let value = output("io.smiles.parse", "input.smi", "F[C@H](Cl)Br");
    assert_eq!(
        value["records"][0]["components"][0]["stereo"][0]["parity"],
        1
    );
    let source = simple_sdf_record("test").replacen("    0.0000", "   10.0000", 1);
    let value = output("io.mol.parse", "input.sdf", &source);
    assert!(
        (value["records"][0]["components"][0]["atoms"][0]["coord"][0]
            .as_f64()
            .unwrap()
            - 10.0)
            .abs()
            < 1e-12
    );
}
#[test]
fn moving_stereo_between_similar_centers_changes_the_observation() {
    let a = output("io.smiles.parse", "input.smi", "F[C@H](Cl)C[C@@H](F)Br");
    let b = output("io.smiles.parse", "input.smi", "F[C@@H](Cl)C[C@H](F)Br");
    assert!(!compare("io.smiles.parse", a, b));
}
#[test]
fn dative_direction_is_significant() {
    let a = json!({"bonds":[{"begin_atom_index":0,"end_atom_index":1,"bond_type":"DATIVE"}]});
    let b = json!({"bonds":[{"begin_atom_index":1,"end_atom_index":0,"bond_type":"DATIVE"}]});
    assert!(!compare("io.mol.parse", a, b));
}
#[test]
fn mol_parse_observes_bond_order_and_coordinates() {
    let molecule = smiles::to_molecules("CCO").unwrap().remove(0);
    let original = molfile::write_v2000(&molecule).unwrap() + "$$$$\n";
    let bond_changed = original.replacen("  1  2  1", "  1  2  2", 1);
    let coord_changed = original.replacen("    0.0000", "   10.0000", 1);
    assert_ne!(original, bond_changed);
    assert_ne!(original, coord_changed);
    let a = output("io.mol.parse", "input.sdf", &original);
    assert!(!compare(
        "io.mol.parse",
        a.clone(),
        output("io.mol.parse", "input.sdf", &bond_changed)
    ));
    assert!(!compare(
        "io.mol.parse",
        a,
        output("io.mol.parse", "input.sdf", &coord_changed)
    ));
}
#[test]
fn sdf_feature_uses_sdf_records_independently_of_file_suffix() {
    let first = simple_sdf_record("first").replace("$$$$", ">  <_Private>\nkept\n\n$$$$");
    let source = first + &simple_sdf_record("second");
    let expected = output("io.sdf.parse", "input.sdf", &source);
    assert_eq!(expected["records"].as_array().unwrap().len(), 2);
    assert_eq!(expected["records"][0]["title"], "first");
    assert_eq!(expected["records"][1]["title"], "second");
    assert_eq!(
        expected["records"][0]["properties"],
        json!([{"name":"_Private","value":"kept"}])
    );
    assert_eq!(expected["records"][1]["properties"], json!([]));
    for path in ["input.mol", "input.mdl"] {
        assert_eq!(output("io.sdf.parse", path, &source), expected);
    }
    let single = simple_sdf_record("standalone");
    let without_delimiter = single.strip_suffix("$$$$\n").unwrap();
    assert_eq!(
        output("io.sdf.parse", "input.mol", without_delimiter),
        output("io.sdf.parse", "input.sdf", &single)
    );
    assert!(evaluate(
        "io.sdf.parse",
        &Input {
            path: "input.sdf".into(),
            text: without_delimiter.into(),
        }
    )
    .is_err());
}

#[test]
fn sdf_writer_preserves_record_boundaries_and_delimiter_like_metadata() {
    let first = simple_sdf_record("cost $$$$ sample").replace(
        "$$$$\n",
        ">  <PRICE$$$$>\nfirst $$$$ line\n>  <not_a_header>\nlast line\n\n$$$$\n",
    );
    let source = first + &simple_sdf_record("");
    let value = output("io.sdf.v2000.write", "input.sdf", &source);
    let text = value["written"][0]["text"].as_str().unwrap();
    let expected = kekule::sdf::parse_str(&source).unwrap();
    let actual = kekule::sdf::parse_str(text).unwrap();
    assert_eq!(actual.records().len(), 2);
    for (actual, expected) in actual.records().iter().zip(expected.records()) {
        assert_eq!(actual.title(), expected.title());
        assert_eq!(
            actual
                .data_fields()
                .iter()
                .map(|f| (f.name(), f.value()))
                .collect::<Vec<_>>(),
            expected
                .data_fields()
                .iter()
                .map(|f| (f.name(), f.value()))
                .collect::<Vec<_>>()
        );
    }
}

#[test]
fn duplicate_sdf_fields_are_not_overwritten() {
    let original =
        simple_sdf_record("molecule").replace("$$$$", ">  <ID>\nfirst\n\n>  <ID>\nlast\n\n$$$$");
    let changed = original.replace("first", "changed");
    let a = output("io.sdf.parse", "input.sdf", &original);
    assert_eq!(
        a["records"][0]["properties"],
        json!([{"name":"ID","value":"first"},{"name":"ID","value":"last"}])
    );
    assert!(!compare(
        "io.sdf.parse",
        a,
        output("io.sdf.parse", "input.sdf", &changed)
    ));
}
#[test]
fn selected_ring_comparison_preserves_cycle_edges() {
    let ring_set = |rings: Value| {
        json!({"records":[{
            "record_index":0,"status":"ok","title":"","rings":rings
        }]})
    };
    let original = ring_set(json!([[0, 1, 2, 3], [4, 5, 6]]));
    for equivalent in [
        json!([[5, 6, 4], [2, 3, 0, 1]]),
        json!([[0, 3, 2, 1], [6, 5, 4]]),
    ] {
        assert!(compare(
            "algo.rings.sssr",
            original.clone(),
            ring_set(equivalent)
        ));
    }
    // The same four atoms can trace a different cycle in a dense graph.
    assert!(!compare(
        "algo.rings.sssr",
        original.clone(),
        ring_set(json!([[0, 1, 3, 2], [4, 5, 6]]))
    ));
    assert!(!compare(
        "algo.rings.sssr",
        original,
        ring_set(json!([[0, 1, 2, 3]]))
    ));
}

#[test]
fn disconnected_inputs_reach_molecular_algorithms() {
    for feature in [
        "io.mol.parse",
        "io.sdf.parse",
        "algo.rings.fast",
        "descriptor.molecular",
    ] {
        let value = output(feature, "input.sdf", &disconnected_sdf_record());
        if feature.starts_with("io.") {
            assert_eq!(
                value["records"][0]["components"].as_array().unwrap().len(),
                2
            );
        } else {
            assert_eq!(value["records"].as_array().unwrap().len(), 2);
        }
    }
}
#[test]
fn writers_emit_source_coordinates_for_independent_reading() {
    let source = simple_sdf_record("test").replacen("    0.0000", "   10.0000", 1);
    let value = output("io.sdf.v2000.write", "input.sdf", &source);
    let text = value["written"][0]["text"].as_str().unwrap();
    assert!(text.contains("10.0000"));
    assert!(text.starts_with("test\n"));
}
#[test]
fn valence_uses_runtime_values_without_aromatic_nitrogen_rewrites() {
    let value = output("algo.valence.rdkit-like", "input.smi", "c1cc[nH]c1");
    assert_eq!(value["records"][0]["atoms"][0]["explicit_valence"], 3);
    let n = &value["records"][0]["atoms"][3];
    assert_eq!(n["explicit_hydrogens"], 1);
}
#[test]
fn mmcif_compares_all_categories_and_distinguishes_missing_tokens() {
    let source = "data_test\n_custom.value ?\nloop_\n_atom_site.id\n_atom_site.label_entity_id\n_atom_site.pdbx_formal_charge\n1 9 2\n";
    let a = output("io.mmcif.parse", "input.cif", source);
    assert_eq!(
        a["blocks"][0]["values"]["_atom_site.label_entity_id"],
        json!(["9"])
    );
    assert!(!compare(
        "io.mmcif.parse",
        a.clone(),
        output("io.mmcif.parse", "input.cif", &source.replace(" ?", " ."))
    ));
    assert!(!compare(
        "io.mmcif.parse",
        a,
        output(
            "io.mmcif.parse",
            "input.cif",
            &source.replace("1 9 2", "1 9 -1")
        )
    ));
}
#[test]
fn mmcif_observations_preserve_quotes_and_source_whitespace() {
    let source = "data_x\n_x.quote 'don't stop'\n_x.text\n;  first  \n  second\t\n;\n";
    let actual = output("io.mmcif.parse", "input.cif", source);
    assert_eq!(
        actual["blocks"][0]["values"]["_x.quote"],
        json!(["don't stop"])
    );
    assert_eq!(
        actual["blocks"][0]["values"]["_x.text"],
        json!(["  first  \n  second\t"])
    );
    let mut trimmed = actual.clone();
    trimmed["blocks"][0]["values"]["_x.text"] = json!(["  first\n  second"]);
    assert!(!compare("io.mmcif.parse", actual, trimmed));
}

#[test]
fn source_membership_includes_every_enamine_sdf_and_every_pdb() {
    for (dataset, feature, count) in [
        ("enamine-diversity", "io.sdf.parse", 50240),
        ("pdb-1000", "bio.secondary-structure.dssp", 1000),
    ] {
        let data = Dataset::open(dataset).unwrap();
        let mut ids = Vec::new();
        for path in data.fixtures(feature) {
            ids.extend(data.members(&path).unwrap());
        }
        assert_eq!(ids.len(), count);
        assert_eq!(
            ids.into_iter().collect::<BTreeSet<_>>(),
            data.selection(usize::MAX).unwrap()
        );
    }
}

fn disconnected_sdf_record() -> String {
    // Two components with interleaved atoms and bonds in source order.
    "disconnected
  benchmark-regression

  5  3  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
  3  5  1  0  0  0  0
  2  4  1  0  0  0  0
  1  3  1  0  0  0  0
M  CHG  1   2  -1
M  END
>  <ID>
salt-regression

$$$$
"
    .to_owned()
}

pub(crate) fn simple_sdf_record(title: &str) -> String {
    format!(
        "{title}
  xtask-test

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"
    )
}

#[test]
fn package_metadata_is_release_consistent() {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let workspace_manifest: toml::Value = toml::from_str(
        &fs::read_to_string(workspace_root.join("Cargo.toml"))
            .expect("workspace manifest should read"),
    )
    .expect("workspace manifest should parse");
    let workspace_package = &workspace_manifest["workspace"]["package"];
    let workspace_version = workspace_package["version"]
        .as_str()
        .expect("workspace package version should be a string");
    assert!(!workspace_version.is_empty());
    assert_eq!(
        workspace_package["repository"].as_str(),
        Some("https://github.com/choutkaj/kekule")
    );

    for (relative_path, package_name, publish) in [
        ("crates/kekule/Cargo.toml", "kekule", true),
        (
            "crates/kekule-potentials/Cargo.toml",
            "kekule-potentials",
            true,
        ),
        ("crates/kekule-traj/Cargo.toml", "kekule-traj", true),
        ("benchmarks/Cargo.toml", "kekule-bench", false),
    ] {
        let manifest: toml::Value = toml::from_str(
            &fs::read_to_string(workspace_root.join(relative_path))
                .expect("package manifest should read"),
        )
        .expect("package manifest should parse");
        assert_eq!(manifest["package"]["name"].as_str(), Some(package_name));
        assert_eq!(
            manifest["package"]["version"]["workspace"].as_bool(),
            Some(true)
        );
        assert_eq!(manifest["package"]["publish"].as_bool(), Some(publish));
    }

    let fuzz_manifest: toml::Value = toml::from_str(
        &fs::read_to_string(workspace_root.join("fuzz/Cargo.toml"))
            .expect("fuzz manifest should read"),
    )
    .expect("fuzz manifest should parse");
    assert_eq!(
        fuzz_manifest["package"]["name"].as_str(),
        Some("kekule-fuzz")
    );
    assert_eq!(fuzz_manifest["package"]["version"].as_str(), Some("0.0.0"));
    assert_eq!(fuzz_manifest["package"]["publish"].as_bool(), Some(false));

    for (relative_path, section, dependency) in [
        ("crates/kekule-traj/Cargo.toml", "dependencies", "kekule"),
        (
            "crates/kekule-potentials/Cargo.toml",
            "dependencies",
            "kekule",
        ),
        (
            "crates/kekule-potentials/Cargo.toml",
            "dev-dependencies",
            "kekule-traj",
        ),
    ] {
        let manifest: toml::Value = toml::from_str(
            &fs::read_to_string(workspace_root.join(relative_path))
                .expect("package manifest should read"),
        )
        .expect("package manifest should parse");
        assert_eq!(
            manifest[section][dependency]["version"].as_str(),
            Some(workspace_version),
            "{relative_path} must require the workspace release version of {dependency}"
        );
    }
}
