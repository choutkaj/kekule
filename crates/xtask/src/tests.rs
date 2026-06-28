use super::*;
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
fn render_dashboard_is_stable_and_uses_boolean_labels() {
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
    ];

    let dashboard = render_dashboard(&features, &BTreeMap::new());

    assert!(dashboard.starts_with("<!doctype html>\n"));
    assert!(dashboard.contains("<table id=\"feature-dashboard\">"));
    assert!(dashboard.contains("<span>Implemented</span>"));
    assert!(!dashboard.contains("Validated"));
    assert!(dashboard.contains("<span>Tiny</span>"));
    assert!(dashboard.contains("<code>a.feature</code>"));
    assert!(dashboard.contains("data-sort-value=\"0\""));
    assert!(dashboard.contains("<code>z.feature</code>"));
    assert!(dashboard.contains("data-sort-value=\"1\""));
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
Run cargo xtask dashboard --check and cargo xtask validate --feature <feature-id> --corpus tiny.
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
Read feature.md. Run cargo test --workspace and cargo xtask validate --feature <feature-id> --corpus tiny.
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
        validation_manifest_path("core.graph", "tiny"),
        PathBuf::from("validation/corpora/tiny/features/core.graph.toml")
    );
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
            validation_required: vec!["tiny".to_owned(), "pubchem-100".to_owned()],
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
            validation_required: vec!["tiny".to_owned(), "pdb-10".to_owned()],
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
        validation_required: vec!["tiny".to_owned()],
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
        evidence_schema_version: Some(VALIDATION_EVIDENCE_SCHEMA_VERSION),
        evidence_hash: Some(evidence.sha256),
        evidence_inputs: evidence.inputs,
        validated_at_unix: 1,
    };
    let status = ValidationStatus {
        feature_id: feature.id.clone(),
        corpora: BTreeMap::from([("tiny".to_owned(), corpus_status)]),
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
        root.join("validation/corpora/tiny/data/example.sdf"),
        "changed fixture\n",
    )
    .expect("fixture should mutate");
    let fixture_changed =
        build_validation_evidence(&root, &manifest_path, &manifest).expect("evidence should build");
    assert_ne!(source_changed.sha256, fixture_changed.sha256);

    fs::write(
        root.join("validation/corpora/tiny/golden/example/data_example.sdf.json.gz"),
        "changed golden\n",
    )
    .expect("golden should mutate");
    let golden_changed =
        build_validation_evidence(&root, &manifest_path, &manifest).expect("evidence should build");
    assert_ne!(fixture_changed.sha256, golden_changed.sha256);
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
        validation_required: vec!["tiny".to_owned()],
    };
    let mut corpus_status = CorpusStatus {
        passed: true,
        fixture_count: 1,
        compared_count: 1,
        reference_tool: "rdkit".to_owned(),
        reference_version: "RDKit test".to_owned(),
        manifest_hash: hash_file(&manifest_path).expect("manifest should hash"),
        evidence_schema_version: Some(VALIDATION_EVIDENCE_SCHEMA_VERSION),
        evidence_hash: Some(evidence.sha256),
        evidence_inputs: evidence.inputs,
        validated_at_unix: 1,
    };
    let status = ValidationStatus {
        feature_id: feature.id.clone(),
        corpora: BTreeMap::from([("tiny".to_owned(), corpus_status.clone())]),
    };
    assert!(overall_validated_at(
        &feature,
        Some(&status),
        &validation_root
    ));

    corpus_status.evidence_schema_version = Some(999);
    let status = ValidationStatus {
        feature_id: feature.id.clone(),
        corpora: BTreeMap::from([("tiny".to_owned(), corpus_status.clone())]),
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
        corpora: BTreeMap::from([("tiny".to_owned(), corpus_status)]),
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
        validation_required: vec!["tiny".to_owned()],
    };
    let status = ValidationStatus {
        feature_id: feature.id.clone(),
        corpora: BTreeMap::from([(
            "tiny".to_owned(),
            CorpusStatus {
                passed: true,
                fixture_count: 1,
                compared_count: 1,
                reference_tool: "rdkit".to_owned(),
                reference_version: "test".to_owned(),
                manifest_hash: "0".repeat(64),
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
fn unsupported_comparison_mode_is_rejected() {
    let manifest = ValidationManifest {
        feature_id: "example".to_owned(),
        corpus_id: "tiny".to_owned(),
        reference_tool: "rdkit".to_owned(),
        reference_version: "RDKit test".to_owned(),
        comparison_mode: "planned".to_owned(),
        fixtures: vec!["data/example.sdf".to_owned()],
        _notes: Vec::new(),
    };
    assert!(validate_comparison_mode(Path::new("example.toml"), &manifest).is_err());
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
    let single = read_smiles_str("CC", SmilesParseOptions).expect("single bond should parse");
    let double = read_smiles_str("C=C", SmilesParseOptions).expect("double bond should parse");
    assert_ne!(
        smiles_sanitized_bonds_json(&single.mol),
        smiles_sanitized_bonds_json(&double.mol)
    );

    let aromatic = read_smiles_str("c1ccccc1", SmilesParseOptions).expect("benzene should parse");
    let mut sanitized_aromatic = aromatic.clone();
    sanitize_small_molecule(&mut sanitized_aromatic, SanitizeOptions::default())
        .expect("benzene should sanitize");
    assert_eq!(
        explicit_valence_json(&sanitized_aromatic.mol, AtomId::new(0)),
        3
    );
    let mut thiophene = read_smiles_str("c1ccsc1", SmilesParseOptions).expect("thiophene parses");
    sanitize_small_molecule(&mut thiophene, SanitizeOptions::default())
        .expect("thiophene should sanitize");
    let sulfur_id = thiophene
        .mol
        .atoms()
        .find_map(|(id, atom)| (atom.element.symbol() == "S").then_some(id))
        .expect("sulfur atom");
    assert_eq!(explicit_valence_json(&thiophene.mol, sulfur_id), 2);
    assert!(smiles_sanitized_bonds_json(&aromatic.mol)
        .iter()
        .all(|bond| bond["bond_type"] == "AROMATIC" && bond["is_aromatic"] == true));

    let labeled =
        read_smiles_str("[13CH3:7]C", SmilesParseOptions).expect("labeled carbon should parse");
    let atoms = smiles_sanitized_atoms_json(&labeled.mol);
    assert!(atoms
        .iter()
        .any(|atom| atom["isotope"] == 13 && atom["atom_map"] == 7));
    assert!(atoms.iter().all(|atom| atom["neighbors"].is_array()));
}

#[test]
fn canonical_smiles_records_do_not_prefilter_unsupported_categories() {
    let root = temp_feature_root("canonical-no-prefilter");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "C/C=C\\C CID:example\n").expect("fixture should write");

    let records = read_canonical_smiles_records(&fixture).expect("records should load");

    assert_eq!(records.len(), 1);
    assert_eq!(records[0].record_index, 0);
    assert_eq!(records[0].status, "parse_error");
    assert_eq!(records[0].input_smiles, "C/C=C\\C");
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
    let corpus_root = validation_root.join("corpora").join("tiny");
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
            "feature_id = \"example\"\ncorpus_id = \"tiny\"\nreference_tool = \"rdkit\"\nreference_version = \"RDKit test\"\ncomparison_mode = \"implementation-golden\"\nfixtures = [\"data/example.sdf\"]\n",
        )
        .expect("manifest should write");
    fs::create_dir_all(corpus_root.join("data")).expect("data dir should create");
    fs::create_dir_all(corpus_root.join("golden").join("example"))
        .expect("golden dir should create");
    fs::write(corpus_root.join("corpus.toml"), "id = \"tiny\"\n")
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
