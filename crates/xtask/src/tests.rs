use super::*;

fn perceive(molecule: &mut Molecule) {
    molecule.perceive().expect("molecule should be perceived");
}

fn one_smiles(input: &str) -> Result<Molecule, String> {
    let mut molecules = kekule::smiles::to_molecules(input).map_err(|error| error.to_string())?;
    if molecules.len() != 1 {
        return Err(format!(
            "expected one connected molecule, found {}",
            molecules.len()
        ));
    }
    Ok(molecules.pop().expect("component count was checked"))
}

#[test]
fn corpus_data_is_required_only_when_requested() {
    assert!(!corpus_requires_data("pubchem-1k", false));
    assert!(corpus_requires_data("pubchem-1k", true));
    assert!(!corpus_requires_data("pdb-1000", false));
    assert!(corpus_requires_data("pdb-1000", true));
}

#[test]
fn tracked_derived_corpora_have_matching_selection_ids_and_exact_prefixes() {
    let mut locks = BTreeMap::new();
    for corpus in ["pubchem-100", "pubchem-1k", "pdb-10", "pdb-100", "pdb-1000"] {
        let descriptor =
            read_tracked_corpus_descriptor(corpus).expect("tracked corpus descriptor should parse");
        let lock = read_source_lock(corpus).expect("tracked source lock should parse");
        check_corpus_lock(&descriptor, &lock)
            .expect("tracked descriptor and source lock should agree");
        locks.insert(corpus.to_owned(), lock);
    }
    check_nested_corpora(&locks).expect("derived corpora should be exact prefixes");
}
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};

#[test]
fn value_after_flag_finds_following_value() {
    let args = vec![
        "benchmark".to_owned(),
        "--benchmark".to_owned(),
        "core.graph".to_owned(),
    ];

    assert_eq!(value_after_flag(&args, "--benchmark"), Some("core.graph"));
}

#[test]
fn package_metadata_is_release_consistent() {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
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
        ("crates/xtask/Cargo.toml", "xtask", false),
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

    for legacy_path in [
        "crates/molecular",
        "crates/molecular-dreiding",
        "crates/molecular-trajectory-io",
        "crates/kekule-dreiding",
        "crates/kekule-trajectory-io",
    ] {
        assert!(
            !workspace_root.join(legacy_path).exists(),
            "legacy package directory still exists: {legacy_path}"
        );
    }
}

#[test]
fn benchmark_jobs_uses_a_memory_safe_default_and_accepts_override() {
    let default_jobs = benchmark_jobs(&[]).expect("default worker count should resolve");
    assert!((1..=4).contains(&default_jobs));
    assert_eq!(
        benchmark_jobs(&["--jobs".to_owned(), "2".to_owned()]).expect("explicit jobs should parse"),
        2
    );
    assert!(benchmark_jobs(&["--jobs".to_owned(), "0".to_owned()]).is_err());
    assert!(benchmark_jobs(&["--jobs".to_owned(), "many".to_owned()]).is_err());
    assert!(benchmark_args(&[
        "--benchmark".to_owned(),
        "all".to_owned(),
        "--jobs".to_owned()
    ])
    .is_err());
    assert!(benchmark_args(&[
        "--benchmark".to_owned(),
        "io.smiles.canonical".to_owned(),
        "--corpus".to_owned(),
        "pubchem-100k".to_owned(),
        "--fixture".to_owned(),
        "data/packs/pack_001.smi".to_owned(),
    ])
    .is_ok());
    assert!(benchmark_args(&[
        "--benchmark".to_owned(),
        "stereo.perception".to_owned(),
        "--corpus".to_owned(),
        "pubchem-100k".to_owned(),
        "--accept-implementation-goldens".to_owned(),
    ])
    .is_ok());
    assert!(benchmark_args(&["--update".to_owned()]).is_err());
    assert!(!include_str!("cli.rs").contains("Some(\"validate\")"));
}

#[test]
fn implementation_golden_acceptance_is_limited_to_manual_semantic_references() {
    let root = temp_workspace_root("accept-implementation-goldens");
    let corpus_root = root.join("benchmarks/corpora/smoke");
    let manifest_path = corpus_root.join("features/stereo.perception.toml");
    fs::create_dir_all(manifest_path.parent().expect("manifest parent"))
        .expect("manifest directory");
    fs::create_dir_all(corpus_root.join("data")).expect("data directory");
    fs::write(corpus_root.join("data/example.smi"), "CC CID:1\n").expect("fixture should write");
    let mut manifest = BenchmarkManifest {
        benchmark_id: "stereo.perception".to_owned(),
        corpus_id: "smoke".to_owned(),
        reference_tool: "rdkit".to_owned(),
        reference_version: "RDKit 2026.03.3".to_owned(),
        comparison_mode: COMPARISON_MODE_IMPLEMENTATION_GOLDEN.to_owned(),
        fixtures: vec!["data/example.smi".to_owned()],
        _notes: Vec::new(),
    };
    assert!(accept_implementation_goldens(&manifest_path, &manifest, 2).is_err());

    manifest.reference_tool = "pubchem-manual-semantic".to_owned();
    manifest.reference_version = "PubChem PUG REST 2026-07-05".to_owned();
    accept_implementation_goldens(&manifest_path, &manifest, 2)
        .expect("manual semantic golden should be accepted");
    let golden_path = corpus_root.join("golden/stereo.perception/data_example.smi.json.gz");
    let golden: Value =
        serde_json::from_str(&read_gzip_string(&golden_path).expect("golden should decompress"))
            .expect("golden should be JSON");
    assert_eq!(golden["feature_id"], "stereo.perception");
    assert_eq!(golden["reference"]["runtime_dependency"], false);
    assert!(golden["expected"]["records"].is_array());

    fs::remove_dir_all(root).ok();
}

#[test]
fn progress_bars_are_compact_and_deterministic() {
    assert_eq!(progress_bar(0, 4), "[------------------------] 0/4   0%");
    assert_eq!(progress_bar(2, 4), "[############------------] 2/4  50%");
    assert_eq!(progress_bar(4, 4), "[########################] 4/4 100%");
    assert_eq!(benchmark_worker_count(16, 3), 3);
    assert_eq!(benchmark_worker_count(0, 3), 1);
}

#[test]
fn benchmark_discovery_is_manifest_only_and_deterministic() {
    let root = temp_workspace_root("benchmark-discovery");
    write_test_benchmark_manifest(&root, "zeta", "corpus-b");
    write_test_benchmark_manifest(&root, "alpha", "corpus-b");
    write_test_benchmark_manifest(&root, "alpha", "corpus-a");

    let targets =
        discover_benchmark_targets_from(&root).expect("benchmark manifests should be discovered");
    let identities = targets
        .iter()
        .map(|target| (target.benchmark_id.as_str(), target.corpus_id.as_str()))
        .collect::<Vec<_>>();
    assert_eq!(
        identities,
        vec![
            ("alpha", "corpus-a"),
            ("alpha", "corpus-b"),
            ("zeta", "corpus-b"),
        ]
    );
    fs::remove_dir_all(root).ok();
}

#[test]
fn selecting_one_benchmark_and_corpus_uses_manifest_directly() {
    let root = temp_workspace_root("benchmark-selection");
    write_test_benchmark_manifest(&root, "io.smiles.parse", "smoke");
    write_test_benchmark_manifest(&root, "io.smiles.write", "smoke");

    let targets = select_benchmark_targets_from(&root, "io.smiles.parse", "smoke")
        .expect("one benchmark and corpus should select");
    assert_eq!(targets.len(), 1);
    assert_eq!(targets[0].benchmark_id, "io.smiles.parse");
    assert_eq!(targets[0].corpus_id, "smoke");
    fs::remove_dir_all(root).ok();
}

#[test]
fn unknown_benchmark_corpus_and_missing_pair_errors_are_clear() {
    let root = temp_workspace_root("benchmark-selection-errors");
    write_test_benchmark_manifest(&root, "io.smiles.parse", "smoke");
    write_test_benchmark_manifest(&root, "io.smiles.write", "other");

    let error = select_benchmark_targets_from(&root, "missing", "smoke")
        .expect_err("unknown benchmark should fail");
    assert!(error.to_string().contains("unknown benchmark: missing"));

    let error = select_benchmark_targets_from(&root, "io.smiles.parse", "missing")
        .expect_err("unknown corpus should fail");
    assert!(error.to_string().contains("unknown corpus: missing"));

    let error = select_benchmark_targets_from(&root, "io.smiles.parse", "other")
        .expect_err("missing benchmark/corpus pair should fail");
    assert!(error
        .to_string()
        .contains("no benchmark manifest for benchmark `io.smiles.parse` and corpus `other`"));
    fs::remove_dir_all(root).ok();
}

#[test]
fn corpus_builders_reuse_locked_membership_unless_reselection_is_explicit() {
    let rdkit_small = include_str!("../../../benchmarks/reference/rdkit/build_corpus.py");
    let rdkit_large = include_str!("../../../benchmarks/reference/rdkit/build_large_corpora.py");
    let biopython = include_str!("../../../benchmarks/reference/biopython/build_corpus.py");

    for builder in [rdkit_small, rdkit_large, biopython] {
        assert!(builder.contains("\"--reselect\""));
        assert!(builder.contains("[\"selection_id\"]"));
        assert!(!builder.contains("selection_seed"));
    }
    assert!(rdkit_small.contains("load_locked_entries("));
    assert!(rdkit_small.contains("check_locked_prefix("));
    assert!(rdkit_large.contains("load_locked_pubchem_candidates("));
    assert!(biopython.contains("load_locked_entries("));
    assert!(biopython.contains("checked-in PDB source locks violate nested-prefix membership"));
}

#[test]
fn implementation_dispatch_uses_current_molfile_benchmark_ids() {
    let root = temp_workspace_root("mol-benchmark-dispatch");
    let fixture = root.join("fixture.sdf");
    fs::write(&fixture, simple_sdf_record("methane")).expect("fixture should write");

    for benchmark_id in [
        "io.mol.v2000.parse",
        "io.mol.v2000.write",
        "io.mol.v3000.parse",
        "io.mol.v3000.write",
    ] {
        let expected = implementation_expected(benchmark_id, "pubchem-1k", &fixture)
            .expect("benchmark should compare");
        assert_eq!(expected["records"][0]["status"], "ok");
    }

    fs::remove_dir_all(root).ok();
}

#[test]
fn implementation_dispatch_supports_mmcif_document_rows() {
    let root = temp_workspace_root("mmcif-document-dispatch");
    let fixture = root.join("fixture.cif");
    fs::write(
        &fixture,
        r#"data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.auth_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.auth_comp_id
_atom_site.label_asym_id
_atom_site.auth_asym_id
_atom_site.label_seq_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA CA . ALA ALA A A 1 1 ? 1.00 10.00 1.0 2.0 3.0 1
"#,
    )
    .expect("fixture should write");

    let expected = implementation_expected("io.mmcif.parse", "pdb-100", &fixture)
        .expect("mmCIF document benchmark should compare");
    let atom_site = &expected["atom_site_rows"];
    assert_eq!(atom_site["status"], "ok");
    assert_eq!(atom_site["row_count"], 1);
    assert_eq!(atom_site["rows"][0]["id"], "1");
    assert_eq!(atom_site["rows"][0]["label_alt_id"], Value::Null);
    assert_eq!(atom_site["rows"][0]["pdbx_PDB_ins_code"], Value::Null);
    assert_eq!(atom_site["rows"][0]["Cartn_z"], "3.0");

    fs::remove_dir_all(root).ok();
}
#[test]
fn implementation_dispatch_supports_hydrogen_transforms() {
    let root = temp_workspace_root("hydrogen-transforms-dispatch");
    let fixture = root.join("fixture.sdf");
    fs::write(&fixture, simple_sdf_record("methane")).expect("fixture should write");

    let expected = implementation_expected("chem.hydrogen-transforms", "pubchem-1k", &fixture)
        .expect("benchmark should compare");
    let record = &expected["records"][0];

    assert_eq!(record["status"], "ok");
    assert_eq!(record["atom_count_after_add"], 5);
    assert_eq!(
        record["added_hydrogens_by_parent"],
        json!([{ "parent_atom_index": 0, "count": 4 }])
    );
    assert_eq!(record["round_trip"]["status"], "ok");

    fs::remove_dir_all(root).ok();
}

#[test]
fn implementation_dispatch_supports_query_benchmark() {
    let root = temp_workspace_root("query-benchmark-dispatch");
    let smarts_fixture = root.join("fixture.smi");
    fs::write(&smarts_fixture, "CCO\nC1=CC=CC=C1\n").expect("fixture should write");

    let parsed = implementation_expected("query.smarts", "pubchem-1k", &smarts_fixture)
        .expect("SMARTS benchmark should compare");
    assert_eq!(parsed["records"][0]["status"], "ok");
    assert_eq!(parsed["records"][0]["atom_count"], 3);
    assert_eq!(parsed["records"][1]["bond_count"], 6);

    let molecule_fixture = root.join("fixture.sdf");
    fs::write(&molecule_fixture, simple_sdf_record("methane")).expect("fixture should write");
    let matched = implementation_expected("algo.substructure.vf2", "pubchem-1k", &molecule_fixture)
        .expect("substructure benchmark should compare");
    assert_eq!(matched["records"][0]["status"], "ok");
    assert_eq!(matched["records"][0]["queries"][0]["smarts"], "[#6]");
    assert_eq!(matched["records"][0]["queries"][0]["matches"], json!([[0]]));

    fs::remove_dir_all(root).ok();
}

#[test]
fn implementation_dispatch_uses_current_isomeric_smiles_benchmark_id() {
    let root = temp_workspace_root("isomeric-smiles-benchmark-dispatch");
    let fixture = root.join("fixture.smi");
    fs::write(
        &fixture,
        [
            "CCO CID:plain",
            "C[C@@H](C(=O)O)N CID:tetrahedral",
            "C(=C\\F)\\F CID:double-bond",
        ]
        .join("\n"),
    )
    .expect("fixture should write");

    let expected = implementation_expected("io.smiles.isomeric", "pubchem-1k", &fixture)
        .expect("benchmark should compare");
    let records = expected["records"]
        .as_array()
        .expect("records should be an array");

    assert_eq!(records.len(), 2);
    assert!(records.iter().all(|record| record["status"] == "ok"));
    assert!(!records[0]["stereo"]["atom_descriptors"]
        .as_array()
        .expect("atom descriptors should be an array")
        .is_empty());
    assert!(!records[1]["stereo"]["bond_descriptors"]
        .as_array()
        .expect("bond descriptors should be an array")
        .is_empty());

    fs::remove_dir_all(root).ok();
}

#[test]
fn nonisomeric_smiles_benchmark_excludes_stereo_syntax() {
    for smiles in ["C[C@H](N)C", "C/C=C/C", "C\\C=C\\C", "C*"] {
        assert_eq!(
            smiles_unsupported_subset_reason(smiles),
            Some("unsupported"),
            "{smiles}"
        );
    }
    assert_eq!(smiles_unsupported_subset_reason("CCO"), None);
}

#[test]
fn stereo_and_nonisomeric_benchmark_use_distinct_smiles_subsets() {
    let root = temp_workspace_root("smiles-benchmark-subsets");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "C[C@H](N)C CID:stereo\n").expect("fixture should write");

    let stereo_records = read_smiles_records(&fixture).expect("stereo records");
    assert_eq!(stereo_records[0].status, "ok");
    assert!(stereo_records[0].molecule.is_some());

    let nonisomeric_records =
        read_nonisomeric_smiles_records(&fixture).expect("nonisomeric records");
    assert_eq!(nonisomeric_records[0].status, "unsupported");
    assert!(nonisomeric_records[0].molecule.is_none());

    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_cip_benchmark_compares_only_descriptor_bearing_records() {
    let root = temp_workspace_root("stereo-cip-descriptor-filter");
    let fixture = root.join("fixture.smi");
    fs::write(
        &fixture,
        [
            "CC CID:no-stereo",
            "C(#N)[Hg-2](C#N)(C#N)C#N.[K+].[K+] CID:unsupported-no-stereo",
            "C[C@H](N)C(=O)O CID:stereo",
        ]
        .join("\n"),
    )
    .expect("fixture should write");

    let expected = implementation_expected("stereo.cip", "pubchem-1k", &fixture)
        .expect("benchmark should compare");
    let records = expected["records"]
        .as_array()
        .expect("records should be an array");

    assert_eq!(records.len(), 1);
    assert_eq!(records[0]["title"], "CID:stereo");
    assert!(!records[0]["atom_descriptors"]
        .as_array()
        .expect("atom descriptors should be an array")
        .is_empty());

    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_cip_benchmark_uses_rdkit_default_hydrogen_indexing() {
    let root = temp_workspace_root("stereo-cip-rdkit-h-index");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "[H][C@](F)(Cl)Br CID:explicit-h\n").expect("fixture should write");

    let expected = implementation_expected("stereo.cip", "pubchem-1k", &fixture)
        .expect("benchmark should compare");
    let records = expected["records"]
        .as_array()
        .expect("records should be an array");

    assert_eq!(records.len(), 1);
    assert_eq!(records[0]["atom_count"], 4);
    assert_eq!(records[0]["bond_count"], 3);
    assert_eq!(records[0]["atom_descriptors"][0]["atom_index"], 0);

    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_cip_benchmark_reads_all_sdf_pack_records() {
    let root = temp_workspace_root("stereo-cip-sdf-pack");
    let fixture = root.join("fixture.sdf");
    fs::write(
        &fixture,
        [
            chiral_wedge_sdf_record("first"),
            chiral_wedge_sdf_record("second"),
        ]
        .join(""),
    )
    .expect("fixture should write");

    let expected = implementation_expected("stereo.cip", "pubchem-1k", &fixture)
        .expect("benchmark should compare");
    let records = expected["records"]
        .as_array()
        .expect("records should be an array");

    assert_eq!(records.len(), 2);
    assert_eq!(records[0]["title"], "first");
    assert_eq!(records[1]["title"], "second");
    assert!(records
        .iter()
        .all(|record| record["atom_count"].as_u64() == Some(5)));
    assert!(records.iter().all(|record| !record["atom_descriptors"]
        .as_array()
        .expect("atom descriptors should be an array")
        .is_empty()));

    fs::remove_dir_all(root).ok();
}

#[test]
fn pack_members_support_custom_sdf_property_and_smiles_title_prefix() {
    let root = temp_workspace_root("pack-members");
    let sdf_path = root.join("pack.sdf");
    fs::write(
        &sdf_path,
        [
            simple_sdf_record_with_property("first", "Catalog ID", "Z111"),
            simple_sdf_record_with_property("second", "Catalog ID", "Z222"),
        ]
        .join(""),
    )
    .expect("sdf pack should write");
    let sdf_pack = SourcePack {
        path: "pack.sdf".to_owned(),
        format: "sdf-v2000".to_owned(),
        count: 2,
        members: vec!["Z111".to_owned(), "Z222".to_owned()],
        sha256: "0".repeat(64),
        member_id_property: Some("Catalog ID".to_owned()),
        member_title_prefix: None,
    };
    assert_eq!(
        read_pack_members(&sdf_path, &sdf_pack).expect("sdf members should read"),
        sdf_pack.members
    );

    let smiles_path = root.join("pack.smi");
    fs::write(&smiles_path, "CC ID:Z111\nCO ID:Z222\n").expect("smiles pack should write");
    let smiles_pack = SourcePack {
        path: "pack.smi".to_owned(),
        format: "smiles".to_owned(),
        count: 2,
        members: vec!["Z111".to_owned(), "Z222".to_owned()],
        sha256: "0".repeat(64),
        member_id_property: None,
        member_title_prefix: Some("ID:".to_owned()),
    };
    assert_eq!(
        read_pack_members(&smiles_path, &smiles_pack).expect("smiles members should read"),
        smiles_pack.members
    );

    fs::remove_dir_all(root).ok();
}

#[test]
fn unsupported_comparison_mode_is_rejected() {
    let manifest = BenchmarkManifest {
        benchmark_id: "example".to_owned(),
        corpus_id: "smoke".to_owned(),
        reference_tool: "rdkit".to_owned(),
        reference_version: "RDKit test".to_owned(),
        comparison_mode: "planned".to_owned(),
        fixtures: vec!["data/example.sdf".to_owned()],
        _notes: Vec::new(),
    };
    assert!(check_comparison_mode(Path::new("example.toml"), &manifest).is_err());
}

#[test]
fn benchmark_comparison_detects_matches_and_differences() {
    let root = temp_workspace_root("comparison-match-and-difference");
    let corpus_root = root.join("benchmarks").join("corpora").join("smoke");
    let manifest_dir = corpus_root.join("features");
    let data_dir = corpus_root.join("data");
    let golden_dir = corpus_root.join("golden").join("io.smiles.parse");
    fs::create_dir_all(&manifest_dir).expect("manifest dir should create");
    fs::create_dir_all(&data_dir).expect("data dir should create");
    fs::create_dir_all(&golden_dir).expect("golden dir should create");
    let manifest_path = manifest_dir.join("io.smiles.parse.toml");
    fs::write(
        &manifest_path,
        "feature_id = \"io.smiles.parse\"\ncorpus_id = \"smoke\"\nreference_tool = \"rdkit\"\nreference_version = \"RDKit test\"\ncomparison_mode = \"implementation-golden\"\nfixtures = [\"data/one.smi\", \"data/two.smi\"]\n",
    )
    .expect("manifest should write");
    for (fixture, text) in [("data/one.smi", "C CID:1\n"), ("data/two.smi", "O CID:2\n")] {
        fs::write(corpus_root.join(fixture), text).expect("fixture should write");
        let expected = if fixture == "data/one.smi" {
            implementation_expected("io.smiles.parse", "smoke", &corpus_root.join(fixture))
                .expect("implementation output should serialize")
        } else {
            json!({
                "records": [{
                    "record_index": 999,
                    "status": "intentionally_wrong",
                }]
            })
        };
        let golden = json!({
            "schema_version": GOLDEN_SCHEMA_VERSION,
            "feature_id": "io.smiles.parse",
            "corpus_id": "smoke",
            "fixture_path": fixture,
            "input_sha256": hash_file(&corpus_root.join(fixture)).expect("fixture should hash"),
            "reference": {
                "tool": "rdkit",
                "version": "RDKit test",
                "runtime_dependency": false,
            },
            "expected": expected,
        });
        write_gzip_json(
            &golden_dir.join(format!("{}.json.gz", slugify_fixture(fixture))),
            &golden,
        );
    }
    let manifest = read_benchmark_manifest(&manifest_path).expect("manifest should read");

    let comparison = compare_golden_outputs(&manifest_path, &manifest, 1, None)
        .expect("comparison should complete");

    assert_eq!(comparison.match_count, 1);
    assert_eq!(comparison.difference_count, 1);
    assert!(comparison
        .first_difference
        .as_deref()
        .is_some_and(|failure| failure.contains("data/two.smi")));
    fs::remove_dir_all(root).ok();
}

#[test]
fn benchmark_run_does_not_rewrite_dashboard_or_results_metadata() {
    let root = temp_workspace_root("benchmark-does-not-rewrite-metadata");
    let manifest_path = write_test_benchmark_manifest(&root, "io.smiles.parse", "smoke");
    let corpus_root = root.join("benchmarks/corpora/smoke");
    let fixture = corpus_root.join("data/example.smi");
    fs::create_dir_all(fixture.parent().expect("fixture parent")).expect("fixture directory");
    fs::write(&fixture, "C CID:1\n").expect("fixture should write");
    fs::write(
        corpus_root.join("sources.lock.json"),
        serde_json::to_vec_pretty(&json!({
            "schema_version": 1,
            "corpus_id": "smoke",
            "source": "test",
            "selection_id": "test",
            "entries": [{
                "id": "1",
                "category": "test",
                "files": [{
                    "path": "data/example.smi",
                    "url": "https://example.invalid/example.smi",
                    "sha256": hash_file(&fixture).expect("fixture should hash"),
                }]
            }],
            "packs": [],
        }))
        .expect("source lock should serialize"),
    )
    .expect("source lock should write");

    let expected = implementation_expected("io.smiles.parse", "smoke", &fixture)
        .expect("implementation output should serialize");
    let golden_dir = corpus_root.join("golden/io.smiles.parse");
    fs::create_dir_all(&golden_dir).expect("golden directory");
    write_gzip_json(
        &golden_dir.join("data_example.smi.json.gz"),
        &json!({
            "schema_version": GOLDEN_SCHEMA_VERSION,
            "feature_id": "io.smiles.parse",
            "corpus_id": "smoke",
            "fixture_path": "data/example.smi",
            "input_sha256": hash_file(&fixture).expect("fixture should hash"),
            "reference": {
                "tool": "rdkit",
                "version": "RDKit test",
                "runtime_dependency": false,
            },
            "expected": expected,
        }),
    );

    let dashboard_path = root.join("features/DASHBOARD.html");
    fs::create_dir_all(dashboard_path.parent().expect("dashboard parent"))
        .expect("dashboard directory");
    fs::write(&dashboard_path, "sentinel dashboard\n").expect("dashboard sentinel");
    let results_path = corpus_root.join("results.toml");
    fs::write(&results_path, "sentinel results\n").expect("results sentinel");

    let target = BenchmarkTarget {
        benchmark_id: "io.smiles.parse".to_owned(),
        corpus_id: "smoke".to_owned(),
        manifest_path,
    };
    let mut progress = BenchmarkProgress::start(1, 1);
    let comparison =
        run_target(&target, None, false, 1, &mut progress).expect("benchmark target should match");
    assert_eq!(comparison.match_count, 1);
    assert_eq!(
        fs::read_to_string(dashboard_path).expect("dashboard sentinel should read"),
        "sentinel dashboard\n"
    );
    assert_eq!(
        fs::read_to_string(results_path).expect("results sentinel should read"),
        "sentinel results\n"
    );
    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_perception_benchmark_records_reference_preparation_errors_per_record() {
    let molecule = one_smiles("C(C)(C)(C)(C)C")
        .expect("pentavalent neutral carbon should remain an interpretation-valid graph");
    let mut record = IndexedSmallRecord {
        record_index: 0,
        title: "invalid neutral-carbon valence".to_owned(),
        molecule,
        sdf_fields: BTreeMap::new(),
    };

    let value = stereo_perception_record_json(&mut record);

    assert_eq!(
        value.get("status").and_then(Value::as_str),
        Some("perception_error")
    );
    assert!(value.get("report").is_none());
}

#[test]
fn smiles_component_benchmarks_preserve_source_record_cardinality() {
    let root = temp_workspace_root("smiles-component-benchmark-cardinality");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "CC.Cl.Cl multi\nC=C connected\n").expect("fixture should write");

    let parsed = implementation_expected("io.smiles.parse", "pubchem-1k", &fixture)
        .expect("parse benchmark should serialize");
    assert_eq!(parsed["records"].as_array().map(Vec::len), Some(2));
    assert_eq!(parsed["records"][0]["status"], "ok");
    assert_eq!(parsed["records"][0]["raw"]["atom_count"], 4);
    assert_eq!(parsed["records"][0]["raw"]["bond_count"], 1);
    assert!(parsed["records"][0].get("normalized_perceived").is_some());
    assert!(parsed["records"][0].get("write_round_trip").is_some());

    let written = implementation_expected("io.smiles.write", "pubchem-1k", &fixture)
        .expect("write benchmark should serialize");
    assert_eq!(written["records"].as_array().map(Vec::len), Some(2));
    assert_eq!(written["records"][0]["status"], "ok");
    assert_eq!(
        written["records"][0]["normalized_perceived"]["atom_count"],
        4
    );
    assert_eq!(
        written["records"][0]["normalized_perceived"]["bond_count"],
        1
    );

    let records = read_nonisomeric_smiles_records(&fixture).expect("fixture should interpret");
    let reparsed = records[0]
        .components
        .iter()
        .map(|molecule| {
            let text = smiles::write(molecule).expect("component should write");
            let document = smiles::parse_str(&text).expect("written component should parse");
            smiles::interpret(&document)
                .expect("written component should interpret")
                .to_molecule()
                .expect("written component should remain connected")
        })
        .collect::<Vec<_>>();
    assert_eq!(
        written["records"][0]["normalized_perceived"],
        smiles_components_perceived_semantic_json(&reparsed)
    );
    let connected_written = smiles::write(
        records[1]
            .molecule
            .as_ref()
            .expect("connected record should have one molecule"),
    )
    .expect("connected record should write");
    let connected_document =
        smiles::parse_str(&connected_written).expect("written component should parse");
    let connected_reparsed = smiles::interpret(&connected_document)
        .expect("written component should interpret")
        .to_molecule()
        .expect("written component should remain connected");
    assert_eq!(
        written["records"][1]["normalized_perceived"],
        smiles_perceived_semantic_json(connected_reparsed)
    );

    let skipped = IndexedSmilesRecord {
        record_index: 0,
        status: "ok".to_owned(),
        title: "missing components".to_owned(),
        input_smiles: "CC.Cl.Cl".to_owned(),
        molecule: None,
        components: Vec::new(),
    };
    let skipped = smiles_write_record_json(&skipped).expect("error record should serialize");
    assert_ne!(skipped["status"], "ok");

    let stereo = implementation_expected("stereo.perception", "pubchem-1k", &fixture)
        .expect("stereo benchmark should serialize");
    assert_eq!(stereo["records"].as_array().map(Vec::len), Some(2));
    assert_eq!(stereo["records"][0]["status"], "ok");
    assert_eq!(stereo["records"][0]["atom_count"], 4);
    assert_eq!(stereo["records"][0]["bond_count"], 1);
    assert!(stereo["records"][0]["report"].get("candidates").is_some());
    assert!(stereo["records"][0]
        .get("source_stereo_element_indices")
        .is_none());

    let stereo_fixture = root.join("stereo.smi");
    fs::write(&stereo_fixture, "F/C=C/F.F/C=C/F directional\n")
        .expect("stereo fixture should write");
    let stereo = implementation_expected("stereo.perception", "pubchem-1k", &stereo_fixture)
        .expect("component stereo benchmark should serialize");
    assert_eq!(stereo["records"].as_array().map(Vec::len), Some(1));
    assert_eq!(
        stereo["records"][0]["report"]["assembled_elements"]
            .as_array()
            .map(Vec::len),
        Some(2)
    );
    assert_eq!(
        stereo["records"][0]["report"]["assembled_elements"][0]["index"],
        0
    );
    assert_eq!(
        stereo["records"][0]["report"]["assembled_elements"][1]["index"],
        1
    );

    fs::remove_dir_all(root).ok();
}

#[test]
fn dssp_comparison_matches_residues_by_source_identity_not_container_order() {
    let mut expected = json!({
        "status": "ok",
        "residues": [
            {"chain_id": "B", "sequence_id": 1, "insertion_code": null, "label_chain_id": "B", "label_sequence_id": 1, "residue_name": "ALA", "sheet": 4, "strand": 8, "ladders": [19, null]},
            {"chain_id": "D", "sequence_id": 1, "insertion_code": null, "label_chain_id": "D", "label_sequence_id": 1, "residue_name": "GLY", "sheet": 7, "strand": 9, "ladders": [21, 19]}
        ]
    });
    let mut actual = json!({
        "status": "ok",
        "residues": [
            {"chain_id": "D", "sequence_id": 1, "insertion_code": null, "label_chain_id": "D", "label_sequence_id": 1, "residue_name": "GLY", "sheet": 12, "strand": 16, "ladders": [31, 30]},
            {"chain_id": "B", "sequence_id": 1, "insertion_code": null, "label_chain_id": "B", "label_sequence_id": 1, "residue_name": "ALA", "sheet": 10, "strand": 15, "ladders": [30, null]}
        ]
    });
    normalize_benchmark_for_comparison_in_place("bio.secondary-structure.dssp", &mut expected);
    normalize_benchmark_for_comparison_in_place("bio.secondary-structure.dssp", &mut actual);
    assert_eq!(expected, actual);
}

#[test]
fn comparison_normalizes_undirected_bonds_and_ring_order() {
    let expected = json!({
        "records": [{
            "bonds": [
                {"index": 0, "begin_atom_index": 5, "end_atom_index": 0, "bond_type": "SINGLE", "stereo": "STEREONONE"}
            ],
            "rings": [[5, 3, 1]]
        }]
    });
    let actual = json!({
        "records": [{
            "bonds": [
                {"index": 7, "begin_atom_index": 0, "end_atom_index": 5, "bond_type": "SINGLE", "stereo": "STEREONONE"}
            ],
            "rings": [[1, 3, 5]]
        }]
    });

    assert_eq!(
        normalize_for_comparison(&expected),
        normalize_for_comparison(&actual)
    );
}

#[test]
fn smiles_semantic_records_assert_topology_and_atom_identity() {
    let single = one_smiles("CC").expect("single bond should parse");
    let double = one_smiles("C=C").expect("double bond should parse");
    assert_ne!(
        smiles_perceived_bonds_json(&single),
        smiles_perceived_bonds_json(&double)
    );

    let aromatic = one_smiles("c1ccccc1").expect("benzene should parse");
    let mut perceived_aromatic = aromatic.clone();
    perceive(&mut perceived_aromatic);
    assert_eq!(
        explicit_valence_json(&perceived_aromatic, AtomId::new(0)),
        3
    );
    let mut aromatic_cyclohexyne = one_smiles("C1=CC#CC=C1").expect("cyclohexyne parses");
    perceive(&mut aromatic_cyclohexyne);
    let alkyne_atoms = aromatic_cyclohexyne
        .bonds()
        .find_map(|(id, bond)| {
            (aromatic_cyclohexyne.bond_is_aromatic(id).ok().flatten() == Some(true)
                && bond.order == BondOrder::Triple)
                .then_some(bond.endpoints())
        })
        .expect("aromaticized triple bond is retained");
    assert_eq!(
        explicit_valence_json(&aromatic_cyclohexyne, alkyne_atoms.0),
        4
    );
    assert_eq!(
        explicit_valence_json(&aromatic_cyclohexyne, alkyne_atoms.1),
        4
    );
    let mut thiophene = one_smiles("c1ccsc1").expect("thiophene parses");
    perceive(&mut thiophene);
    let sulfur_id = thiophene
        .atoms()
        .find_map(|(id, atom)| (atom.element.symbol() == "S").then_some(id))
        .expect("sulfur atom");
    assert_eq!(explicit_valence_json(&thiophene, sulfur_id), 2);
    let mut phosphorus_ring = one_smiles("C(F)(F)(F)P1P(P(P(P1C(F)(F)F)C(F)(F)F)C(F)(F)F)C(F)(F)F")
        .expect("phosphorus ring parses");
    perceive(&mut phosphorus_ring);
    for (phosphorus_id, _phosphorus) in phosphorus_ring
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "P")
    {
        assert_eq!(
            phosphorus_ring.atom_is_aromatic(phosphorus_id).unwrap(),
            Some(true)
        );
        assert_eq!(explicit_valence_json(&phosphorus_ring, phosphorus_id), 3);
    }
    let mut phosphinine = one_smiles("C1=CC=PC=C1").expect("phosphinine parses");
    perceive(&mut phosphinine);
    let phosphinine_phosphorus = phosphinine
        .atoms()
        .find_map(|(id, atom)| (atom.element.symbol() == "P").then_some(id))
        .expect("phosphinine phosphorus");
    assert_eq!(
        explicit_valence_json(&phosphinine, phosphinine_phosphorus),
        3
    );
    let document = kekule::smiles::parse_str("CN(C)CCO.C1=CC=C2C(=C1)C3=NC4=C5C=CC=CC5=C([N-]4)N=C6C7=CC=CC=C7C(=N6)N=C8C9=CC=CC=C9C(=N8)N=C2[N-]3.[Cu+2]")
        .expect("anionic macrocycle mixture parses");
    let mut anionic_macrocycle = kekule::smiles::interpret(&document)
        .expect("anionic macrocycle mixture interprets")
        .to_molecules()
        .swap_remove(1);
    perceive(&mut anionic_macrocycle);
    let anionic_nitrogen = anionic_macrocycle
        .atoms()
        .find_map(|(id, atom)| {
            (atom.element.symbol() == "N"
                && atom.formal_charge < 0
                && anionic_macrocycle.atom_is_aromatic(id).ok().flatten() == Some(true))
            .then_some(id)
        })
        .expect("anionic aromatic nitrogen");
    assert_eq!(
        explicit_valence_json(&anionic_macrocycle, anionic_nitrogen),
        2
    );
    let mut cyclopentadienyl =
        one_smiles("[CH-]1[C-]=[C-][C-]=[C-]1").expect("cyclopentadienyl anion parses");
    perceive(&mut cyclopentadienyl);
    let anionic_carbon_with_h = cyclopentadienyl
        .atoms()
        .find_map(|(id, atom)| {
            (atom.element.symbol() == "C"
                && atom.formal_charge < 0
                && cyclopentadienyl.atom_is_aromatic(id).ok().flatten() == Some(true)
                && atom.hydrogens.explicit_count() > 0)
                .then_some(id)
        })
        .expect("anionic aromatic carbon with explicit hydrogen");
    let anionic_carbon = cyclopentadienyl
        .atom(anionic_carbon_with_h)
        .expect("anionic carbon should exist");
    assert_eq!(
        explicit_valence_json(&cyclopentadienyl, anionic_carbon_with_h)
            + anionic_carbon.hydrogens.explicit_count(),
        3
    );
    let mut substituted_cyclopentadienyl =
        one_smiles("C[C-]1[C-]=[C-][C-]=[C-]1").expect("substituted cyclopentadienyl parses");
    perceive(&mut substituted_cyclopentadienyl);
    let substituted_anionic_carbon = substituted_cyclopentadienyl
        .atoms()
        .find_map(|(id, atom)| {
            let degree = substituted_cyclopentadienyl
                .incident_bonds(id)
                .ok()?
                .count();
            (atom.element.symbol() == "C"
                && atom.formal_charge < 0
                && substituted_cyclopentadienyl
                    .atom_is_aromatic(id)
                    .ok()
                    .flatten()
                    == Some(true)
                && degree == 3)
                .then_some(id)
        })
        .expect("substituted anionic carbon");
    assert_eq!(
        explicit_valence_json(&substituted_cyclopentadienyl, substituted_anionic_carbon,),
        3
    );
    let mut fused_triazine =
        one_smiles("O=[N+]([O-])c2cc(-c1nn5c(=O)c(C=Cc3c(O)ccc4c3cccc4)nnc5s1)ccc2")
            .expect("fused triazine should parse");
    perceive(&mut fused_triazine);
    let tricoordinate_aromatic_nitrogen = fused_triazine
        .atoms()
        .find_map(|(id, atom)| {
            let aromatic_degree = fused_triazine
                .incident_bonds(id)
                .ok()?
                .filter(|(bond, _)| {
                    fused_triazine.bond_is_aromatic(*bond).ok().flatten() == Some(true)
                })
                .count();
            (atom.element.symbol() == "N"
                && fused_triazine.atom_is_aromatic(id).ok().flatten() == Some(true)
                && aromatic_degree >= 3)
                .then_some(id)
        })
        .expect("tri-coordinate aromatic nitrogen");
    assert_eq!(
        explicit_valence_json(&fused_triazine, tricoordinate_aromatic_nitrogen),
        3
    );
    let localized_bonds = smiles_perceived_bonds_json(&aromatic);
    assert_eq!(
        localized_bonds
            .iter()
            .filter(|bond| bond["bond_type"] == "SINGLE" && bond["is_aromatic"] == false)
            .count(),
        3
    );
    assert_eq!(
        localized_bonds
            .iter()
            .filter(|bond| bond["bond_type"] == "DOUBLE" && bond["is_aromatic"] == false)
            .count(),
        3
    );
    assert!(perceived_aromatic
        .bonds()
        .all(|(_, bond)| matches!(bond.order, BondOrder::Single | BondOrder::Double)));
    assert!(smiles_perceived_bonds_json(&perceived_aromatic)
        .iter()
        .all(|bond| bond["bond_type"] == "AROMATIC" && bond["is_aromatic"] == true));

    let labeled = one_smiles("[13CH3:7]C").expect("labeled carbon should parse");
    let atoms = smiles_perceived_atoms_json(&labeled);
    assert!(atoms
        .iter()
        .any(|atom| atom["isotope"] == 13 && atom["atom_map"] == 7));
    assert!(atoms.iter().all(|atom| atom["neighbors"].is_array()));
}

#[test]
fn canonical_smiles_records_do_not_prefilter_unsupported_categories() {
    let root = temp_workspace_root("canonical-no-prefilter");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "* CID:example\n").expect("fixture should write");

    let records = read_canonical_smiles_records(&fixture).expect("records should load");

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].record_index, 0);
    assert_eq!(records[0].status, "parse_error");
    assert_eq!(records[0].input_smiles, "*");
    assert!(records[0].molecule.is_none());
}

#[test]
fn canonical_smiles_benchmark_perceives_before_writing() {
    let root = temp_workspace_root("canonical-perceive-before-write");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "C1=CC=CC=C1 CID:benzene\n").expect("fixture should write");

    let records = read_canonical_smiles_records(&fixture).expect("records should load");
    let item =
        canonical_smiles_record_json(&records[0], true).expect("canonical record should render");

    assert_eq!(item["status"], "ok");
    assert_eq!(item["canonical_smiles"], "c1ccccc1");
}

#[test]
fn canonical_smiles_benchmark_matches_rdkit_parse_status_for_invalid_input() {
    let root = temp_workspace_root("canonical-invalid-input");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "[Cl-](Br)Br CID:invalid\n").expect("fixture should write");

    let records = read_canonical_smiles_records(&fixture).expect("records should load");
    let item =
        canonical_smiles_record_json(&records[0], false).expect("canonical record should render");

    assert_eq!(item["status"], "parse_error");
}

#[test]
fn smiles_semantics_match_rdkit_aromatic_carbonyl_valence() {
    let molecule =
        one_smiles("CCCCCCCc1cc2c(=O)ccn(O)c2cc1").expect("aromatic carbonyl SMILES should parse");

    let item = smiles_perceived_semantic_json(molecule);
    let atoms = item["atoms"]
        .as_array()
        .expect("perceived atoms should be an array");

    assert!(atoms.iter().any(|atom| {
        atom["symbol"] == "C"
            && atom["aromatic"] == true
            && atom["explicit_valence"] == 4
            && atom["neighbors"].as_array().is_some_and(|neighbors| {
                neighbors.iter().any(|neighbor| {
                    neighbor["bond_type"] == "DOUBLE"
                        && neighbor["atom"]
                            .as_str()
                            .is_some_and(|key| key.starts_with("008|O|0|0|0|0|false|2|"))
                })
            })
    }));
    assert!(!atoms.iter().any(|atom| {
        atom["symbol"] == "C" && atom["aromatic"] == true && atom["explicit_valence"] == 5
    }));
    assert!(atoms.iter().any(|atom| {
        atom["symbol"] == "N"
            && atom["aromatic"] == true
            && atom["explicit_valence"] == 3
            && atom["neighbors"].as_array().is_some_and(|neighbors| {
                neighbors.iter().any(|neighbor| {
                    neighbor["bond_type"] == "SINGLE"
                        && neighbor["atom"]
                            .as_str()
                            .is_some_and(|key| key.starts_with("008|O|0|0|0|1|false|1|"))
                })
            })
    }));
    assert!(!atoms.iter().any(|atom| {
        atom["symbol"] == "N" && atom["aromatic"] == true && atom["explicit_valence"] == 4
    }));
}

#[test]
fn smiles_semantics_match_rdkit_aromatic_nh_no_implicit_flag() {
    let molecule = one_smiles("[nH]1cccc1").expect("aromatic nH SMILES should parse");

    let item = smiles_perceived_semantic_json(molecule);
    let atoms = item["atoms"]
        .as_array()
        .expect("perceived atoms should be an array");

    assert!(atoms.iter().any(|atom| {
        atom["symbol"] == "N"
            && atom["aromatic"] == true
            && atom["explicit_hydrogens"] == 1
            && atom["implicit_hydrogens"] == 0
            && atom["no_implicit_hydrogens"] == false
    }));
}

#[test]
fn smiles_semantics_derive_promoted_aromatic_nh_valence_without_feedback() {
    let molecule = one_smiles("CCOC(=O)C1=C(C(=C(N1)C)C(=O)OC(C)(C)C)C")
        .expect("substituted pyrrole SMILES should parse");

    let item = smiles_perceived_semantic_json(molecule);
    let atoms = item["atoms"]
        .as_array()
        .expect("perceived atoms should be an array");

    assert!(
        atoms.iter().any(|atom| {
            atom["symbol"] == "N"
                && atom["aromatic"] == true
                && atom["explicit_hydrogens"] == 1
                && atom["implicit_hydrogens"] == 0
                && atom["no_implicit_hydrogens"] == false
                && atom["explicit_valence"] == 3
        }),
        "{atoms:#?}"
    );
    assert!(!atoms.iter().any(|atom| {
        atom["symbol"] == "N" && atom["aromatic"] == true && atom["explicit_valence"] == 4
    }));
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

fn simple_sdf_record_with_property(title: &str, property: &str, value: &str) -> String {
    let mut record = simple_sdf_record(title);
    let marker = "M  END\n";
    let replacement = format!("M  END\n>  <{property}>  (1)\n{value}\n\n");
    record = record.replacen(marker, &replacement, 1);
    record
}

fn chiral_wedge_sdf_record(title: &str) -> String {
    format!(
        "{title}
  xtask-test

  5  4  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 F   0  0  0  0  0  0  0  0  0  0  0  0
   -1.0000    0.0000    0.0000 Cl  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.0000    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -1.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  1  0  0  0
  1  3  1  0  0  0  0
  1  4  1  0  0  0  0
  1  5  1  0  0  0  0
M  END
$$$$
"
    )
}

fn temp_workspace_root(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time should be available")
        .as_nanos();
    let root = env::temp_dir().join(format!("kekule-xtask-{label}-{}-{nonce}", process::id()));
    fs::create_dir_all(&root).expect("temporary workspace root should create");
    root
}

fn write_test_benchmark_manifest(root: &Path, benchmark: &str, corpus: &str) -> PathBuf {
    let corpus_root = root.join("benchmarks").join("corpora").join(corpus);
    let manifest_path = corpus_root
        .join("features")
        .join(format!("{benchmark}.toml"));
    fs::create_dir_all(manifest_path.parent().expect("manifest parent"))
        .expect("manifest directory should create");
    fs::write(
        corpus_root.join("corpus.toml"),
        format!(
            "id = \"{corpus}\"\ntitle = \"Test corpus\"\nkind = \"small-molecule\"\nready = true\nexpected_count = 1\nlocal_only = true\nselection_id = \"test\"\nformats = [\"smiles\"]\nbuild_command = \"test\"\n"
        ),
    )
    .expect("corpus descriptor should write");
    fs::write(
        &manifest_path,
        format!(
            "feature_id = \"{benchmark}\"\ncorpus_id = \"{corpus}\"\nreference_tool = \"rdkit\"\nreference_version = \"RDKit test\"\ncomparison_mode = \"implementation-golden\"\nfixtures = [\"data/example.smi\"]\n"
        ),
    )
    .expect("benchmark manifest should write");
    manifest_path
}

fn write_gzip_json(path: &Path, value: &Value) {
    let file = fs::File::create(path).expect("gzip file should create");
    let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    encoder
        .write_all(value.to_string().as_bytes())
        .expect("gzip json should write");
    encoder.finish().expect("gzip json should finish");
}
