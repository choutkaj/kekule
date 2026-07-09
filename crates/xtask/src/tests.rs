use super::*;
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};

#[test]
fn value_after_flag_finds_following_value() {
    let args = vec![
        "validate".to_owned(),
        "--feature".to_owned(),
        "core.graph".to_owned(),
    ];

    assert_eq!(value_after_flag(&args, "--feature"), Some("core.graph"));
}

#[test]
fn read_feature_parses_typed_metadata() {
    let root = temp_feature_root("read-feature");
    write_feature(
        &root,
        "example.feature",
        r#"id = "example.feature"
title = "Example"
area = "infrastructure"
version = 2
implemented = false
validated = true
description = "Example feature."
depends_on = ["core.graph"]
validation_required = []
"#,
    );

    let feature = read_feature(&root.join("example.feature").join("feature.toml"))
        .expect("feature should parse");

    assert_eq!(feature.id, "example.feature");
    assert_eq!(feature.version, 2);
    assert!(!feature.implemented);
    assert!(feature.validated);
    assert_eq!(feature.depends_on, vec!["core.graph"]);
    fs::remove_dir_all(root).ok();
}

#[test]
fn read_feature_rejects_bad_boolean_deprecated_keys_missing_docs_and_directory_mismatch() {
    let root = temp_feature_root("bad-feature");
    write_feature(
        &root,
        "bad.bool",
        r#"id = "bad.bool"
title = "Bad"
area = "infrastructure"
version = 1
implemented = maybe
validated = false
description = "Bad feature."
depends_on = []
validation_required = []
"#,
    );
    assert!(read_feature(&root.join("bad.bool").join("feature.toml")).is_err());

    write_feature(
        &root,
        "bad.deprecated",
        r#"id = "bad.deprecated"
title = "Bad"
area = "infrastructure"
version = 1
priority = "P0"
implemented = false
validated = false
description = "Bad feature."
depends_on = []
validation_required = []
"#,
    );
    assert!(read_feature(&root.join("bad.deprecated").join("feature.toml")).is_err());

    write_feature(
        &root,
        "bad.version",
        r#"id = "bad.version"
title = "Bad"
area = "infrastructure"
version = 0
implemented = false
validated = false
description = "Bad feature."
depends_on = []
validation_required = []
"#,
    );
    assert!(read_feature(&root.join("bad.version").join("feature.toml")).is_err());

    write_feature_without_doc(
        &root,
        "missing.doc",
        r#"id = "missing.doc"
title = "Bad"
area = "infrastructure"
version = 1
implemented = false
validated = false
description = "Bad feature."
depends_on = []
validation_required = []
"#,
    );
    assert!(read_feature(&root.join("missing.doc").join("feature.toml")).is_err());

    write_feature(
        &root,
        "real.id",
        r#"id = "wrong.id"
title = "Bad"
area = "infrastructure"
version = 1
implemented = false
validated = false
description = "Bad feature."
depends_on = []
validation_required = []
"#,
    );
    assert!(read_feature(&root.join("real.id").join("feature.toml")).is_err());
    fs::remove_dir_all(root).ok();
}

#[test]
fn read_features_sorts_skips_templates_and_validates_dependencies() {
    let root = temp_feature_root("feature-set");
    write_feature(
        &root,
        "z.feature",
        r#"id = "z.feature"
title = "Zed"
area = "core"
version = 1
implemented = true
validated = false
description = "Z feature."
depends_on = ["a.feature"]
validation_required = []
"#,
    );
    write_feature(
        &root,
        "a.feature",
        r#"id = "a.feature"
title = "Aye"
area = "core"
version = 1
implemented = false
validated = false
description = "A feature."
depends_on = []
validation_required = []
"#,
    );
    fs::create_dir_all(root.join("_template")).expect("template dir should create");
    fs::write(root.join("_template").join("feature.toml"), "not = valid")
        .expect("template metadata should write");

    let features = read_features_from(&root).expect("feature set should parse");

    assert_eq!(
        features
            .iter()
            .map(|feature| feature.id.as_str())
            .collect::<Vec<_>>(),
        vec!["a.feature", "z.feature"]
    );

    write_feature(
        &root,
        "bad.dependency",
        r#"id = "bad.dependency"
title = "Bad"
area = "core"
version = 1
implemented = false
validated = false
description = "Bad dependency."
depends_on = ["missing.feature"]
validation_required = []
"#,
    );
    assert!(read_features_from(&root).is_err());
    fs::remove_dir_all(root).ok();
}

#[test]
fn render_dashboard_is_stable_and_uses_compact_validation_cells() {
    let features = vec![
        Feature {
            id: "a.feature".to_owned(),
            title: "Aye".to_owned(),
            area: "core".to_owned(),
            version: 1,
            implemented: false,
            validated: false,
            description: "A feature.".to_owned(),
            depends_on: Vec::new(),
            validation_required: Vec::new(),
        },
        Feature {
            id: "z.feature".to_owned(),
            title: "Zed".to_owned(),
            area: "io".to_owned(),
            version: 3,
            implemented: true,
            validated: false,
            description: "Z feature.".to_owned(),
            depends_on: vec!["a.feature".to_owned()],
            validation_required: Vec::new(),
        },
        Feature {
            id: "failing.feature".to_owned(),
            title: "Failing".to_owned(),
            area: "validation".to_owned(),
            version: 1,
            implemented: true,
            validated: false,
            description: "Feature with counted failures.".to_owned(),
            depends_on: Vec::new(),
            validation_required: vec!["smoke".to_owned()],
        },
        Feature {
            id: "missing.feature".to_owned(),
            title: "Missing".to_owned(),
            area: "validation".to_owned(),
            version: 1,
            implemented: true,
            validated: false,
            description: "Feature without recorded status.".to_owned(),
            depends_on: Vec::new(),
            validation_required: vec!["pubchem-100".to_owned()],
        },
    ];
    let statuses = BTreeMap::from([(
        "failing.feature".to_owned(),
        ValidationStatus {
            feature_id: "failing.feature".to_owned(),
            corpora: BTreeMap::from([(
                "smoke".to_owned(),
                CorpusStatus::from_failed_run(FailedValidationRun {
                    fixture_count: 7,
                    compared_count: 4,
                    failed_count: 3,
                    first_failure: "fixture `data/bad.sdf` differs".to_owned(),
                    reference_tool: "rdkit".to_owned(),
                    reference_version: "RDKit test".to_owned(),
                    manifest_hash: "0".repeat(64),
                })
                .expect("failed status should build"),
            )]),
        },
    )]);
    let corpus_info = BTreeMap::from([
        (
            "smoke".to_owned(),
            CorpusDashboardInfo {
                id: "smoke".to_owned(),
                label: "smoke".to_owned(),
                title: "Checked-in external smoke corpus".to_owned(),
                expected_count: 7,
            },
        ),
        (
            "pubchem-1k".to_owned(),
            CorpusDashboardInfo {
                id: "pubchem-1k".to_owned(),
                label: "PubChem 1k".to_owned(),
                title: "PubChem deterministic 1000-compound corpus".to_owned(),
                expected_count: 1000,
            },
        ),
    ]);

    let dashboard = render_dashboard(&features, &statuses, &corpus_info);

    assert!(dashboard.starts_with("<!doctype html>\n"));
    assert!(dashboard.contains("<table id=\"feature-dashboard\">"));
    assert!(dashboard.contains("th.area, td.area { text-align: left; }"));
    assert!(dashboard.contains("<th class=\"compact area\" data-sort-type=\"text\" title=\"Area\"><button class=\"sort\" type=\"button\" aria-label=\"Sort by Area\">Area</button></th>"));
    assert!(dashboard.contains("<td class=\"compact area\" data-sort-value=\"core\">core</td>"));
    assert!(!dashboard
        .contains("aria-label=\"Sort by Area\"><span class=\"rotated-label\">Area</span>"));
    assert!(dashboard.contains(
        "<span class=\"rotated-label\"><span class=\"rotated-name\">Implemented</span></span>"
    ));
    assert!(dashboard.contains("height: 168px"));
    assert!(dashboard.contains("left: calc(50% + 23px)"));
    assert!(dashboard.contains("bottom: 12px"));
    assert!(dashboard.contains("width: 144px"));
    assert!(dashboard.contains("height: 46px"));
    assert!(dashboard.contains("display: flex"));
    assert!(dashboard.contains("rotate(-90deg)"));
    assert!(dashboard.contains("transform-origin: left bottom"));
    assert!(dashboard.contains("overflow: hidden"));
    assert!(dashboard.contains("white-space: nowrap"));
    assert!(!dashboard.contains("Validated"));
    assert!(dashboard.contains(
        "<span class=\"rotated-name\">smoke</span><br><span class=\"rotated-count\">(n=7)</span>"
    ));
    assert!(dashboard.contains("<span class=\"rotated-name\">pubchem-1k</span><br><span class=\"rotated-count\">(n=1000)</span>"));
    assert!(dashboard.contains("<code>a.feature</code>"));
    assert!(dashboard.contains("data-sort-value=\"0\""));
    assert!(dashboard.contains("<code>z.feature</code>"));
    assert!(dashboard.contains("data-sort-value=\"1\""));
    assert!(dashboard.contains("aria-label=\"failed: 3 non-passing case(s)\""));
    assert!(dashboard.contains("<span class=\"count\">3</span>"));
    assert!(dashboard.contains("<span class=\"unknown\">?</span>unknown"));
    assert!(dashboard.contains(
        "<span class=\"unknown\" aria-label=\"unknown\" title=\"no recorded validation status\">?</span>"
    ));
    assert!(!dashboard.contains(
        "<span class=\"bad\" aria-label=\"failed\" title=\"no recorded validation status\">"
    ));
    assert!(dashboard.contains("button.addEventListener('click'"));
    assert!(dashboard.ends_with('\n'));
}

#[test]
fn skill_metadata_parser_and_check_validate_repo_skill_contract() {
    let root = temp_feature_root("skills-check");
    write_skill(
        &root,
        "feature-work",
        r#"---
name: feature-work
description: Builder skill.
---
# Feature Work
add -> optional research -> plan -> implement
Use feature.md. Set implemented = true with evidence and declare validation_required.
Molecular validation fixtures must be externally supplied.
Run cargo xtask dashboard --check and cargo xtask validate --feature <feature-id> --corpus smoke.
"#,
    );
    write_skill(
        &root,
        "feature-review",
        r#"---
name: feature-review
description: Review skill.
---
# Feature Review
Independent audit for architecture and validation claims.
Read feature.md. Run cargo test --workspace and cargo xtask validate --feature <feature-id> --corpus smoke.
"#,
    );

    check_skills(&root).expect("skills should pass");

    fs::write(
        root.join("feature-review").join("SKILL.md"),
        "# Missing frontmatter",
    )
    .expect("skill should rewrite");
    assert!(check_skills(&root).is_err());
    fs::remove_dir_all(root).ok();
}

#[test]
fn validation_manifest_path_is_feature_scoped() {
    assert_eq!(
        validation_manifest_path("core.graph", "smoke"),
        PathBuf::from("validation/corpora/smoke/features/core.graph.toml")
    );
}

#[test]
fn validate_jobs_defaults_to_available_parallelism_and_accepts_override() {
    let default_jobs = validation_jobs(&[]).expect("default worker count should resolve");
    assert!(default_jobs >= 1);
    assert_eq!(
        validation_jobs(&["--jobs".to_owned(), "2".to_owned()])
            .expect("explicit jobs should parse"),
        2
    );
    assert!(validation_jobs(&["--jobs".to_owned(), "0".to_owned()]).is_err());
    assert!(validation_jobs(&["--jobs".to_owned(), "many".to_owned()]).is_err());
    assert!(validate_args(&[
        "--feature".to_owned(),
        "all".to_owned(),
        "--jobs".to_owned()
    ])
    .is_err());
}

#[test]
fn progress_bars_are_compact_and_deterministic() {
    assert_eq!(progress_bar(0, 4), "[------------------------] 0/4   0%");
    assert_eq!(progress_bar(2, 4), "[############------------] 2/4  50%");
    assert_eq!(progress_bar(4, 4), "[########################] 4/4 100%");
    assert_eq!(validation_worker_count(16, 3), 3);
    assert_eq!(validation_worker_count(0, 3), 1);
}

#[test]
fn all_selectors_expand_only_applicable_feature_corpus_pairs() {
    let features = vec![
        Feature {
            id: "small".to_owned(),
            title: "Small".to_owned(),
            area: "io".to_owned(),
            version: 1,
            implemented: true,
            validated: false,
            description: "Small feature.".to_owned(),
            depends_on: Vec::new(),
            validation_required: vec!["smoke".to_owned(), "pubchem-100".to_owned()],
        },
        Feature {
            id: "macro".to_owned(),
            title: "Macro".to_owned(),
            area: "bio".to_owned(),
            version: 1,
            implemented: true,
            validated: false,
            description: "Macro feature.".to_owned(),
            depends_on: Vec::new(),
            validation_required: vec!["smoke".to_owned(), "pdb-10".to_owned()],
        },
    ];

    assert_eq!(
        validation_targets(&features, "all", "pubchem-100")
            .into_iter()
            .map(|(feature, corpus)| (feature.id.as_str(), corpus))
            .collect::<Vec<_>>(),
        vec![("small", "pubchem-100".to_owned())]
    );
    assert_eq!(validation_targets(&features, "small", "all").len(), 2);
    assert_eq!(
        validation_targets(&features, "macro", "pubchem-100").len(),
        0
    );
}

#[test]
fn implementation_dispatch_uses_current_molfile_feature_ids() {
    let root = temp_feature_root("mol-feature-dispatch");
    let fixture = root.join("fixture.sdf");
    fs::write(&fixture, simple_sdf_record("methane")).expect("fixture should write");

    for feature in [
        "io.mol.v2000.parse",
        "io.mol.v2000.write",
        "io.mol.v3000.parse",
        "io.mol.v3000.write",
    ] {
        let expected =
            implementation_expected(feature, "smoke", &fixture).expect("feature should compare");
        assert_eq!(expected["records"][0]["status"], "ok");
    }

    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_cip_validation_compares_only_descriptor_bearing_records() {
    let root = temp_feature_root("stereo-cip-descriptor-filter");
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

    let expected =
        implementation_expected("stereo.cip", "smoke", &fixture).expect("feature should compare");
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
fn stereo_cip_validation_uses_rdkit_default_hydrogen_indexing() {
    let root = temp_feature_root("stereo-cip-rdkit-h-index");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "[H][C@](F)(Cl)Br CID:explicit-h\n").expect("fixture should write");

    let expected =
        implementation_expected("stereo.cip", "smoke", &fixture).expect("feature should compare");
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
fn stereo_cip_validation_reads_all_sdf_pack_records() {
    let root = temp_feature_root("stereo-cip-sdf-pack");
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

    let expected =
        implementation_expected("stereo.cip", "smoke", &fixture).expect("feature should compare");
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
    let root = temp_feature_root("pack-members");
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
fn current_status_drives_overall_validation_and_metadata_sync() {
    let root = temp_feature_root("status-sync");
    let (features_root, validation_root, manifest_path) = write_evidence_test_repo(&root);
    let metadata_path = features_root.join("example").join("feature.toml");

    let feature = Feature {
        id: "example".to_owned(),
        title: "Example".to_owned(),
        area: "io".to_owned(),
        version: 1,
        implemented: true,
        validated: false,
        description: "Example feature.".to_owned(),
        depends_on: Vec::new(),
        validation_required: vec!["smoke".to_owned()],
    };
    let manifest = read_validation_manifest(&manifest_path).expect("manifest should read");
    let evidence =
        build_validation_evidence(&root, &manifest_path, &manifest).expect("evidence should build");
    let corpus_status = CorpusStatus {
        passed: true,
        fixture_count: 1,
        compared_count: 1,
        reference_tool: "rdkit".to_owned(),
        reference_version: "RDKit test".to_owned(),
        manifest_hash: hash_file(&manifest_path).expect("manifest should hash"),
        failed_count: 0,
        first_failure: None,
        evidence_schema_version: Some(VALIDATION_EVIDENCE_SCHEMA_VERSION),
        evidence_hash: Some(evidence.sha256),
        evidence_inputs: evidence.inputs,
        validated_at_unix: 1,
    };
    let status = ValidationStatus {
        feature_id: feature.id.clone(),
        corpora: BTreeMap::from([("smoke".to_owned(), corpus_status)]),
    };
    let statuses = BTreeMap::from([(feature.id.clone(), status.clone())]);

    assert!(overall_validated_at(
        &feature,
        Some(&status),
        &validation_root
    ));
    sync_feature_validation_flags_at(
        std::slice::from_ref(&feature),
        &statuses,
        &features_root,
        &validation_root,
    )
    .expect("metadata should sync");
    assert!(fs::read_to_string(&metadata_path)
        .expect("metadata should read")
        .contains("validated = true"));

    fs::write(&manifest_path, "changed = true\n").expect("manifest should change");
    assert!(!overall_validated_at(
        &feature,
        Some(&status),
        &validation_root
    ));
    fs::remove_dir_all(root).ok();
}

#[test]
fn evidence_changes_after_material_input_changes() {
    let root = temp_feature_root("evidence-change");
    let (_, _, manifest_path) = write_evidence_test_repo(&root);
    let manifest = read_validation_manifest(&manifest_path).expect("manifest should read");
    let original =
        build_validation_evidence(&root, &manifest_path, &manifest).expect("evidence should build");

    fs::write(root.join("crates/molecules/src/lib.rs"), "changed source\n")
        .expect("source should mutate");
    let source_changed =
        build_validation_evidence(&root, &manifest_path, &manifest).expect("evidence should build");
    assert_ne!(original.sha256, source_changed.sha256);

    fs::write(
        root.join("validation/corpora/smoke/data/example.sdf"),
        "changed fixture\n",
    )
    .expect("fixture should mutate");
    let fixture_changed =
        build_validation_evidence(&root, &manifest_path, &manifest).expect("evidence should build");
    assert_ne!(source_changed.sha256, fixture_changed.sha256);

    fs::write(
        root.join("validation/corpora/smoke/golden/example/data_example.sdf.json.gz"),
        "changed golden\n",
    )
    .expect("golden should mutate");
    let golden_changed =
        build_validation_evidence(&root, &manifest_path, &manifest).expect("evidence should build");
    assert_ne!(fixture_changed.sha256, golden_changed.sha256);
    fs::remove_dir_all(root).ok();
}

#[test]
fn manual_semantic_reference_evidence_does_not_require_generator_files() {
    let root = temp_feature_root("manual-reference-evidence");
    let (_, _, manifest_path) = write_evidence_test_repo(&root);
    fs::write(
        &manifest_path,
        "feature_id = \"example\"\ncorpus_id = \"smoke\"\nreference_tool = \"enamine-manual-semantic\"\nreference_version = \"Enamine Discovery Diversity Set 2026-07-05\"\ncomparison_mode = \"implementation-golden\"\nfixtures = [\"data/example.sdf\"]\n",
    )
    .expect("manual manifest should write");
    fs::remove_dir_all(root.join("validation/reference")).ok();

    let manifest = read_validation_manifest(&manifest_path).expect("manifest should read");
    let evidence =
        build_validation_evidence(&root, &manifest_path, &manifest).expect("evidence should build");

    assert!(evidence
        .inputs
        .iter()
        .all(|input| !input.path.starts_with("validation/reference/")));
    fs::remove_dir_all(root).ok();
}

#[test]
fn evidence_hash_normalizes_text_line_endings() {
    let root = temp_feature_root("evidence-line-endings");
    let path = root.join("source.rs");
    fs::write(&path, "first\nsecond\n").expect("LF source should write");
    let lf_hash = hash_evidence_file(&path).expect("LF evidence should hash");

    fs::write(&path, "first\r\nsecond\r\n").expect("CRLF source should write");
    let crlf_hash = hash_evidence_file(&path).expect("CRLF evidence should hash");

    assert_eq!(lf_hash, crlf_hash);
    fs::remove_dir_all(root).ok();
}

#[test]
fn current_status_requires_known_nonempty_evidence() {
    let root = temp_feature_root("status-rejects-stale");
    let (_, validation_root, manifest_path) = write_evidence_test_repo(&root);
    let manifest = read_validation_manifest(&manifest_path).expect("manifest should read");
    let evidence =
        build_validation_evidence(&root, &manifest_path, &manifest).expect("evidence should build");
    let feature = Feature {
        id: "example".to_owned(),
        title: "Example".to_owned(),
        area: "io".to_owned(),
        version: 1,
        implemented: true,
        validated: false,
        description: "Example feature.".to_owned(),
        depends_on: Vec::new(),
        validation_required: vec!["smoke".to_owned()],
    };
    let mut corpus_status = CorpusStatus {
        passed: true,
        fixture_count: 1,
        compared_count: 1,
        reference_tool: "rdkit".to_owned(),
        reference_version: "RDKit test".to_owned(),
        manifest_hash: hash_file(&manifest_path).expect("manifest should hash"),
        failed_count: 0,
        first_failure: None,
        evidence_schema_version: Some(VALIDATION_EVIDENCE_SCHEMA_VERSION),
        evidence_hash: Some(evidence.sha256),
        evidence_inputs: evidence.inputs,
        validated_at_unix: 1,
    };
    let status = ValidationStatus {
        feature_id: feature.id.clone(),
        corpora: BTreeMap::from([("smoke".to_owned(), corpus_status.clone())]),
    };
    assert!(overall_validated_at(
        &feature,
        Some(&status),
        &validation_root
    ));

    corpus_status.evidence_schema_version = Some(999);
    let status = ValidationStatus {
        feature_id: feature.id.clone(),
        corpora: BTreeMap::from([("smoke".to_owned(), corpus_status.clone())]),
    };
    assert!(!overall_validated_at(
        &feature,
        Some(&status),
        &validation_root
    ));

    corpus_status.evidence_schema_version = Some(VALIDATION_EVIDENCE_SCHEMA_VERSION);
    corpus_status.compared_count = 0;
    let status = ValidationStatus {
        feature_id: feature.id.clone(),
        corpora: BTreeMap::from([("smoke".to_owned(), corpus_status)]),
    };
    assert!(!overall_validated_at(
        &feature,
        Some(&status),
        &validation_root
    ));
    fs::remove_dir_all(root).ok();
}

#[test]
fn dashboard_text_comparison_ignores_platform_line_endings() {
    assert_eq!(
        normalize_text_line_endings("one\r\ntwo\rthree\n"),
        "one\ntwo\nthree\n"
    );
}

#[test]
fn recorded_dashboard_status_is_portable_but_current_status_is_content_addressed() {
    let feature = Feature {
        id: "portable.feature".to_owned(),
        title: "Portable".to_owned(),
        area: "infrastructure".to_owned(),
        version: 1,
        implemented: true,
        validated: true,
        description: "Portable dashboard evidence.".to_owned(),
        depends_on: Vec::new(),
        validation_required: vec!["smoke".to_owned()],
    };
    let status = ValidationStatus {
        feature_id: feature.id.clone(),
        corpora: BTreeMap::from([(
            "smoke".to_owned(),
            CorpusStatus {
                passed: true,
                fixture_count: 1,
                compared_count: 1,
                reference_tool: "rdkit".to_owned(),
                reference_version: "test".to_owned(),
                manifest_hash: "0".repeat(64),
                failed_count: 0,
                first_failure: None,
                evidence_schema_version: Some(VALIDATION_EVIDENCE_SCHEMA_VERSION),
                evidence_hash: Some("1".repeat(64)),
                evidence_inputs: vec![EvidenceInput {
                    path: "missing.fixture".to_owned(),
                    sha256: "2".repeat(64),
                }],
                validated_at_unix: 1,
            },
        )]),
    };

    assert!(recorded_overall_validated(&feature, Some(&status)));
    assert!(!overall_validated_at(
        &feature,
        Some(&status),
        Path::new("definitely-missing-validation-root")
    ));
}

#[test]
fn old_status_toml_defaults_failure_summary_fields() {
    let status: CorpusStatusFile = toml::from_str(
        r#"
corpus_id = "smoke"

[features.example]
passed = true
fixture_count = 1
compared_count = 1
reference_tool = "rdkit"
reference_version = "RDKit test"
manifest_hash = "0000000000000000000000000000000000000000000000000000000000000000"
evidence_schema_version = 2
evidence_hash = "1111111111111111111111111111111111111111111111111111111111111111"
evidence_inputs = []
validated_at_unix = 1
"#,
    )
    .expect("old status shape should deserialize");
    let corpus_status = status
        .features
        .get("example")
        .expect("feature status should exist");

    assert_eq!(corpus_status.failed_count, 0);
    assert_eq!(corpus_status.first_failure, None);
}

#[test]
fn unsupported_comparison_mode_is_rejected() {
    let manifest = ValidationManifest {
        feature_id: "example".to_owned(),
        corpus_id: "smoke".to_owned(),
        reference_tool: "rdkit".to_owned(),
        reference_version: "RDKit test".to_owned(),
        comparison_mode: "planned".to_owned(),
        fixtures: vec!["data/example.sdf".to_owned()],
        _notes: Vec::new(),
    };
    assert!(validate_comparison_mode(Path::new("example.toml"), &manifest).is_err());
}

#[test]
fn validation_comparison_counts_multiple_fixture_failures() {
    let root = temp_feature_root("comparison-counts-failures");
    let corpus_root = root.join("validation").join("corpora").join("smoke");
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
            "expected": {
                "records": [{
                    "record_index": 999,
                    "status": "intentionally_wrong",
                }]
            },
        });
        write_gzip_json(
            &golden_dir.join(format!("{}.json.gz", slugify_fixture(fixture))),
            &golden,
        );
    }
    let manifest = read_validation_manifest(&manifest_path).expect("manifest should read");

    let comparison = validate_golden_outputs(&manifest_path, &manifest, 1, None)
        .expect("comparison should complete");

    assert_eq!(comparison.compared_count, 0);
    assert_eq!(comparison.failed_count, 2);
    assert!(comparison
        .first_failure
        .as_deref()
        .is_some_and(|failure| failure.contains("data/one.smi")));
    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_perception_validation_records_sanitize_errors_per_record() {
    let molecule = smiles::read_str("C(#N)[Hg-2](C#N)(C#N)C#N.[K+].[K+]")
        .expect("unsupported valence molecule should parse");
    let mut record = IndexedSmallRecord {
        record_index: 0,
        title: "unsupported element".to_owned(),
        molecule,
    };

    let value = stereo_perception_record_json(&mut record);

    assert_eq!(
        value.get("status").and_then(Value::as_str),
        Some("sanitize_error")
    );
    assert!(value.get("report").is_none());
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
    let single = smiles::read_str_with_options("CC", SmilesParseOptions::default())
        .expect("single bond should parse");
    let double = smiles::read_str_with_options("C=C", SmilesParseOptions::default())
        .expect("double bond should parse");
    assert_ne!(
        smiles_sanitized_bonds_json(single.graph()),
        smiles_sanitized_bonds_json(double.graph())
    );

    let aromatic = smiles::read_str_with_options("c1ccccc1", SmilesParseOptions::default())
        .expect("benzene should parse");
    let mut sanitized_aromatic = aromatic.clone();
    perception::sanitize_with_options(&mut sanitized_aromatic, SanitizeOptions::default())
        .expect("benzene should sanitize");
    assert_eq!(
        explicit_valence_json(sanitized_aromatic.graph(), AtomId::new(0)),
        3
    );
    let mut thiophene = smiles::read_str_with_options("c1ccsc1", SmilesParseOptions::default())
        .expect("thiophene parses");
    perception::sanitize_with_options(&mut thiophene, SanitizeOptions::default())
        .expect("thiophene should sanitize");
    let sulfur_id = thiophene
        .graph()
        .atoms()
        .find_map(|(id, atom)| (atom.element.symbol() == "S").then_some(id))
        .expect("sulfur atom");
    assert_eq!(explicit_valence_json(thiophene.graph(), sulfur_id), 2);
    let mut anionic_macrocycle = smiles::read_str_with_options(
        "CN(C)CCO.C1=CC=C2C(=C1)C3=NC4=C5C=CC=CC5=C([N-]4)N=C6C7=CC=CC=C7C(=N6)N=C8C9=CC=CC=C9C(=N8)N=C2[N-]3.[Cu+2]",
        SmilesParseOptions::default(),
    )
    .expect("anionic macrocycle parses");
    perception::sanitize_with_options(&mut anionic_macrocycle, SanitizeOptions::default())
        .expect("anionic macrocycle should sanitize");
    let anionic_nitrogen = anionic_macrocycle
        .graph()
        .atoms()
        .find_map(|(id, atom)| {
            (atom.element.symbol() == "N" && atom.formal_charge < 0 && atom.aromatic).then_some(id)
        })
        .expect("anionic aromatic nitrogen");
    assert_eq!(
        explicit_valence_json(anionic_macrocycle.graph(), anionic_nitrogen),
        2
    );
    let mut cyclopentadienyl =
        smiles::read_str_with_options("[CH-]1[C-]=[C-][C-]=[C-]1", SmilesParseOptions::default())
            .expect("cyclopentadienyl anion parses");
    perception::sanitize_with_options(&mut cyclopentadienyl, SanitizeOptions::default())
        .expect("cyclopentadienyl anion should sanitize");
    let anionic_carbon_with_h = cyclopentadienyl
        .graph()
        .atoms()
        .find_map(|(id, atom)| {
            (atom.element.symbol() == "C"
                && atom.formal_charge < 0
                && atom.aromatic
                && atom.explicit_hydrogens > 0)
                .then_some(id)
        })
        .expect("anionic aromatic carbon with explicit hydrogen");
    let anionic_carbon = cyclopentadienyl
        .graph()
        .atom(anionic_carbon_with_h)
        .expect("anionic carbon should exist");
    assert_eq!(
        explicit_valence_json(cyclopentadienyl.graph(), anionic_carbon_with_h)
            + anionic_carbon.explicit_hydrogens,
        3
    );
    let mut fused_triazine = smiles::read_str_with_options(
        "O=[N+]([O-])c2cc(-c1nn5c(=O)c(C=Cc3c(O)ccc4c3cccc4)nnc5s1)ccc2",
        SmilesParseOptions::default(),
    )
    .expect("fused triazine should parse");
    perception::sanitize_with_options(&mut fused_triazine, SanitizeOptions::default())
        .expect("fused triazine should sanitize");
    let tricoordinate_aromatic_nitrogen = fused_triazine
        .graph()
        .atoms()
        .find_map(|(id, atom)| {
            let aromatic_degree = fused_triazine
                .graph()
                .incident_bonds(id)
                .ok()?
                .filter(|(_, bond)| bond.aromatic)
                .count();
            (atom.element.symbol() == "N" && atom.aromatic && aromatic_degree >= 3).then_some(id)
        })
        .expect("tri-coordinate aromatic nitrogen");
    assert_eq!(
        explicit_valence_json(fused_triazine.graph(), tricoordinate_aromatic_nitrogen),
        3
    );
    assert!(smiles_sanitized_bonds_json(aromatic.graph())
        .iter()
        .all(|bond| bond["bond_type"] == "AROMATIC" && bond["is_aromatic"] == true));

    let labeled = smiles::read_str_with_options("[13CH3:7]C", SmilesParseOptions::default())
        .expect("labeled carbon should parse");
    let atoms = smiles_sanitized_atoms_json(labeled.graph());
    assert!(atoms
        .iter()
        .any(|atom| atom["isotope"] == 13 && atom["atom_map"] == 7));
    assert!(atoms.iter().all(|atom| atom["neighbors"].is_array()));
}

#[test]
fn canonical_smiles_records_do_not_prefilter_unsupported_categories() {
    let root = temp_feature_root("canonical-no-prefilter");
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
fn canonical_smiles_validation_sanitizes_before_writing() {
    let root = temp_feature_root("canonical-sanitize-before-write");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "C1=CC=CC=C1 CID:benzene\n").expect("fixture should write");

    let records = read_canonical_smiles_records(&fixture).expect("records should load");
    let item =
        canonical_smiles_record_json(&records[0], true).expect("canonical record should render");

    assert_eq!(item["status"], "ok");
    assert_eq!(item["canonical_smiles"], "c1ccccc1");
}

#[test]
fn smiles_semantics_match_rdkit_aromatic_carbonyl_valence() {
    let molecule = smiles::read_str_with_options(
        "CCCCCCCc1cc2c(=O)ccn(O)c2cc1",
        SmilesParseOptions::default(),
    )
    .expect("aromatic carbonyl SMILES should parse");

    let item = smiles_sanitized_semantic_json(molecule);
    let atoms = item["atoms"]
        .as_array()
        .expect("sanitized atoms should be an array");

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
    let molecule = smiles::read_str_with_options("[nH]1cccc1", SmilesParseOptions::default())
        .expect("aromatic nH SMILES should parse");

    let item = smiles_sanitized_semantic_json(molecule);
    let atoms = item["atoms"]
        .as_array()
        .expect("sanitized atoms should be an array");

    assert!(atoms.iter().any(|atom| {
        atom["symbol"] == "N"
            && atom["aromatic"] == true
            && atom["explicit_hydrogens"] == 1
            && atom["implicit_hydrogens"] == 0
            && atom["no_implicit_hydrogens"] == false
    }));
}

#[test]
fn smiles_semantics_match_rdkit_promoted_aromatic_nh_valence() {
    let molecule = smiles::read_str_with_options(
        "CCOC(=O)C1=C(C(=C(N1)C)C(=O)OC(C)(C)C)C",
        SmilesParseOptions::default(),
    )
    .expect("substituted pyrrole SMILES should parse");

    let item = smiles_sanitized_semantic_json(molecule);
    let atoms = item["atoms"]
        .as_array()
        .expect("sanitized atoms should be an array");

    assert!(atoms.iter().any(|atom| {
        atom["symbol"] == "N"
            && atom["aromatic"] == true
            && atom["explicit_hydrogens"] == 1
            && atom["implicit_hydrogens"] == 0
            && atom["no_implicit_hydrogens"] == false
            && atom["explicit_valence"] == 3
    }));
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

fn temp_feature_root(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time should be available")
        .as_nanos();
    let root = env::temp_dir().join(format!("molecules-xtask-{label}-{}-{nonce}", process::id()));
    fs::create_dir_all(&root).expect("temp feature root should create");
    root
}

fn write_evidence_test_repo(root: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let features_root = root.join("features");
    let validation_root = root.join("validation");
    let feature_dir = features_root.join("example");
    let corpus_root = validation_root.join("corpora").join("smoke");
    let manifest_dir = corpus_root.join("features");
    fs::create_dir_all(&feature_dir).expect("feature dir should create");
    fs::create_dir_all(&manifest_dir).expect("manifest dir should create");
    fs::write(
        feature_dir.join("feature.toml"),
        "id = \"example\"\nvalidated = false\n",
    )
    .expect("metadata should write");
    fs::write(feature_dir.join("feature.md"), "# Example\n").expect("feature doc should write");
    let manifest_path = manifest_dir.join("example.toml");
    fs::write(
            &manifest_path,
            "feature_id = \"example\"\ncorpus_id = \"smoke\"\nreference_tool = \"rdkit\"\nreference_version = \"RDKit test\"\ncomparison_mode = \"implementation-golden\"\nfixtures = [\"data/example.sdf\"]\n",
        )
        .expect("manifest should write");
    fs::create_dir_all(corpus_root.join("data")).expect("data dir should create");
    fs::create_dir_all(corpus_root.join("golden").join("example"))
        .expect("golden dir should create");
    fs::write(corpus_root.join("corpus.toml"), "id = \"smoke\"\n")
        .expect("corpus descriptor should write");
    fs::write(corpus_root.join("sources.lock.json"), "{}\n").expect("source lock should write");
    fs::write(corpus_root.join("data").join("example.sdf"), "fixture\n")
        .expect("fixture should write");
    fs::write(
        corpus_root
            .join("golden")
            .join("example")
            .join("data_example.sdf.json.gz"),
        "golden\n",
    )
    .expect("golden should write");
    fs::write(root.join("Cargo.toml"), "[workspace]\n").expect("cargo toml should write");
    fs::write(root.join("Cargo.lock"), "# lock\n").expect("cargo lock should write");
    for path in [
        "crates/molecules/Cargo.toml",
        "crates/xtask/Cargo.toml",
        "crates/molecules/src/lib.rs",
        "crates/xtask/src/main.rs",
        "validation/reference/rdkit/run_feature.py",
        "validation/reference/rdkit/environment.yml",
    ] {
        let path = root.join(path);
        fs::create_dir_all(path.parent().expect("test path should have parent"))
            .expect("test parent should create");
        fs::write(path, "test\n").expect("test evidence file should write");
    }
    (features_root, validation_root, manifest_path)
}

fn write_feature(root: &Path, id: &str, metadata: &str) {
    let dir = root.join(id);
    fs::create_dir_all(&dir).expect("feature dir should create");
    fs::write(dir.join("feature.toml"), metadata).expect("feature metadata should write");
    fs::write(dir.join("feature.md"), "# Feature\n").expect("feature doc should write");
}

fn write_feature_without_doc(root: &Path, id: &str, metadata: &str) {
    let dir = root.join(id);
    fs::create_dir_all(&dir).expect("feature dir should create");
    fs::write(dir.join("feature.toml"), metadata).expect("feature metadata should write");
}

fn write_skill(root: &Path, name: &str, text: &str) {
    let dir = root.join(name);
    fs::create_dir_all(&dir).expect("skill dir should create");
    fs::write(dir.join("SKILL.md"), text).expect("skill should write");
}

fn write_gzip_json(path: &Path, value: &Value) {
    let file = fs::File::create(path).expect("gzip file should create");
    let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    encoder
        .write_all(value.to_string().as_bytes())
        .expect("gzip json should write");
    encoder.finish().expect("gzip json should finish");
}
