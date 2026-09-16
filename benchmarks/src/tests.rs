use super::*;
fn output(feature: &str, path: &str, text: &str) -> Value {
    evaluate(
        feature,
        "regression",
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
    first_json_diff("$", &a, &b).is_none()
}

#[test]
fn same_local_atom_environments_do_not_hide_different_graphs() {
    let a = output("io.smiles.parse", "input.smi", "C1CCCCC1");
    let b = output("io.smiles.parse", "input.smi", "C1CC1.C1CC1");
    assert!(!compare("io.smiles.parse", a, b));
}
#[test]
fn query_records_separate_the_query_from_the_source_title() {
    let value = output("query.smarts", "input.smi", "[#6]-[#8]  external query");
    assert_eq!(value["records"][0]["smarts"], "[#6]-[#8]");
    assert_eq!(value["records"][0]["title"], "external query");
    assert_eq!(value["records"][0]["atom_count"], 2);
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

fn simple_sdf_record(title: &str) -> String {
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
