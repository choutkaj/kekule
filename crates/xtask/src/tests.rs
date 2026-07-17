use super::*;

#[test]
fn corpus_data_is_required_only_when_requested() {
    assert!(!corpus_requires_data("pubchem-1k", false));
    assert!(corpus_requires_data("pubchem-1k", true));
    assert!(!corpus_requires_data("pdb-1000", false));
    assert!(corpus_requires_data("pdb-1000", true));
}
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
domains = ["infrastructure"]
version = 2
status = "planned"
description = "Example feature."
depends_on = ["core.graph"]
validation_required = []
"#,
    );

    let feature = read_feature(&root.join("example.feature").join("feature.toml"))
        .expect("feature should parse");

    assert_eq!(feature.id, "example.feature");
    assert_eq!(feature.version, 2);
    assert_eq!(feature.status, FeatureStatus::Planned);
    assert!(!feature.status.has_implementation());
    assert_eq!(feature.domains, vec![FeatureDomain::Infrastructure]);
    assert_eq!(feature.depends_on, vec!["core.graph"]);
    fs::remove_dir_all(root).ok();
}

#[test]
fn feature_status_parses_the_release_vocabulary() {
    #[derive(Deserialize)]
    struct StatusOnly {
        status: FeatureStatus,
    }

    for (name, expected, has_implementation) in [
        ("planned", FeatureStatus::Planned, false),
        ("experimental", FeatureStatus::Experimental, true),
        ("supported", FeatureStatus::Supported, true),
        ("deprecated", FeatureStatus::Deprecated, true),
    ] {
        let parsed: StatusOnly =
            toml::from_str(&format!("status = \"{name}\"")).expect("release status should parse");
        assert_eq!(parsed.status, expected);
        assert_eq!(parsed.status.has_implementation(), has_implementation);
    }
}

#[test]
fn local_only_corpus_descriptors_match_the_registry() {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    for corpus_id in [
        "pubchem-1k",
        "pubchem-100k",
        "pl-rex",
        "enamine-diversity",
        "pdb-100",
        "pdb-1000",
    ] {
        let path = workspace_root
            .join("validation/corpora")
            .join(corpus_id)
            .join("corpus.toml");
        let text = fs::read_to_string(&path).expect("corpus descriptor should read");
        let descriptor: CorpusDescriptor =
            toml::from_str(&text).expect("corpus descriptor should parse");
        let registered = validation_corpus(corpus_id).expect("corpus should be registered");
        assert_eq!(descriptor.id, registered.id);
        assert_eq!(descriptor.local_only, registered.local_only);
    }
}

#[test]
fn read_feature_accepts_local_only_required_corpora() {
    let root = temp_feature_root("local-only-required-corpus");
    write_feature(
        &root,
        "valid.requirement",
        r#"id = "valid.requirement"
title = "Valid requirement"
area = "infrastructure"
domains = ["infrastructure"]
version = 1
status = "supported"
description = "Valid local-only validation requirement."
depends_on = []
validation_required = ["pubchem-1k"]
"#,
    );

    let feature = read_feature(&root.join("valid.requirement").join("feature.toml"))
        .expect("known local-only corpora may be required");
    assert_eq!(feature.validation_required, vec!["pubchem-1k"]);
    assert!(VALIDATION_CORPORA.iter().all(|corpus| corpus.local_only));
    assert!(validation_corpus("smoke").is_none());
    assert!(validation_corpus("pubchem-100").is_none());
    assert!(validation_corpus("pdb-10").is_none());
    fs::remove_dir_all(root).ok();
}

#[test]
fn required_validation_assignments_need_manifests() {
    let root = temp_feature_root("required-validation-manifest");
    let feature = Feature {
        id: "required.feature".to_owned(),
        title: "Required".to_owned(),
        area: "validation".to_owned(),
        domains: vec![FeatureDomain::SmallMolecule],
        version: 1,
        status: FeatureStatus::Supported,
        description: "Required validation feature.".to_owned(),
        depends_on: Vec::new(),
        validation_required: vec!["pubchem-1k".to_owned()],
    };

    let error = validate_required_manifests_from(&root, std::slice::from_ref(&feature))
        .expect_err("missing required manifest should fail");
    assert!(error
        .to_string()
        .contains("requires validation corpus `pubchem-1k`"));

    let manifest = validation_manifest_path_from(&root, &feature.id, "pubchem-1k");
    fs::create_dir_all(manifest.parent().expect("manifest parent"))
        .expect("manifest directory should create");
    fs::write(&manifest, "").expect("manifest marker should write");
    validate_required_manifests_from(&root, &[feature])
        .expect("present required manifest should pass");
    fs::remove_dir_all(root).ok();
}
#[test]
fn read_feature_rejects_required_validation_for_planned_features() {
    let root = temp_feature_root("planned-required-validation");
    write_feature(
        &root,
        "planned.validation",
        r#"id = "planned.validation"
title = "Planned validation"
area = "infrastructure"
domains = ["infrastructure"]
version = 1
status = "planned"
description = "Planned features cannot require validation."
depends_on = []
validation_required = ["pubchem-1k"]
"#,
    );

    let error = read_feature(&root.join("planned.validation").join("feature.toml"))
        .expect_err("planned features should not require validation");
    assert!(error
        .to_string()
        .contains("status `planned` but declares required validation"));
    fs::remove_dir_all(root).ok();
}

#[test]
fn read_feature_rejects_unknown_status_removed_keys_and_shape_errors() {
    let root = temp_feature_root("bad-feature");
    write_feature(
        &root,
        "bad.bool",
        r#"id = "bad.bool"
title = "Bad"
area = "infrastructure"
domains = ["infrastructure"]
version = 1
status = "unknown"
description = "Bad feature."
depends_on = []
validation_required = []
"#,
    );
    assert!(read_feature(&root.join("bad.bool").join("feature.toml")).is_err());

    write_feature(
        &root,
        "bad.implemented",
        r#"id = "bad.implemented"
title = "Bad"
area = "infrastructure"
domains = ["infrastructure"]
version = 1
implemented = false
description = "Removed metadata field."
depends_on = []
validation_required = []
"#,
    );
    assert!(read_feature(&root.join("bad.implemented").join("feature.toml")).is_err());

    write_feature(
        &root,
        "bad.deprecated",
        r#"id = "bad.deprecated"
title = "Bad"
area = "infrastructure"
domains = ["infrastructure"]
version = 1
priority = "P0"
status = "planned"
description = "Bad feature."
depends_on = []
validation_required = []
"#,
    );
    assert!(read_feature(&root.join("bad.deprecated").join("feature.toml")).is_err());

    write_feature(
        &root,
        "bad.removed",
        r#"id = "bad.removed"
title = "Bad"
area = "infrastructure"
domains = ["infrastructure"]
version = 1
status = "planned"
validated = false
description = "Removed metadata field."
depends_on = []
validation_required = []
"#,
    );
    assert!(read_feature(&root.join("bad.removed").join("feature.toml")).is_err());

    write_feature(
        &root,
        "bad.version",
        r#"id = "bad.version"
title = "Bad"
area = "infrastructure"
domains = ["infrastructure"]
version = 0
status = "planned"
description = "Bad feature."
depends_on = []
validation_required = []
"#,
    );
    assert!(read_feature(&root.join("bad.version").join("feature.toml")).is_err());

    write_feature(
        &root,
        "bad.domains",
        r#"id = "bad.domains"
title = "Bad"
area = "infrastructure"
domains = ["infrastructure", "small-molecule"]
version = 1
status = "planned"
description = "Bad feature domains."
depends_on = []
validation_required = []
"#,
    );
    assert!(read_feature(&root.join("bad.domains").join("feature.toml")).is_err());

    write_feature_without_doc(
        &root,
        "missing.doc",
        r#"id = "missing.doc"
title = "Bad"
area = "infrastructure"
domains = ["infrastructure"]
version = 1
status = "planned"
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
domains = ["infrastructure"]
version = 1
status = "planned"
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
domains = ["small-molecule", "macromolecule"]
version = 1
status = "experimental"
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
domains = ["small-molecule", "macromolecule"]
version = 1
status = "experimental"
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
domains = ["small-molecule"]
version = 1
status = "planned"
description = "Bad dependency."
depends_on = ["missing.feature"]
validation_required = []
"#,
    );
    assert!(read_features_from(&root).is_err());
    fs::remove_dir_all(root).ok();
}

#[test]
fn feature_graph_rejects_duplicate_self_cyclic_and_incompatible_dependencies() {
    let base = feature_for_test("base", FeatureStatus::Supported, &[]);

    let duplicate = feature_for_test("duplicate", FeatureStatus::Experimental, &["base", "base"]);
    let error = validate_feature_set(&[base.clone(), duplicate])
        .expect_err("duplicate dependencies should be rejected");
    assert!(error.to_string().contains("more than once"));

    let self_dependent = feature_for_test("self", FeatureStatus::Planned, &["self"]);
    let error =
        validate_feature_set(&[self_dependent]).expect_err("self dependencies should be rejected");
    assert!(error.to_string().contains("depends on itself"));

    let cycle_a = feature_for_test("cycle.a", FeatureStatus::Planned, &["cycle.b"]);
    let cycle_b = feature_for_test("cycle.b", FeatureStatus::Planned, &["cycle.a"]);
    let error = validate_feature_set(&[cycle_a, cycle_b])
        .expect_err("dependency cycles should be rejected");
    assert!(error
        .to_string()
        .contains("feature dependency graph contains a cycle: cycle.a -> cycle.b -> cycle.a"));

    let experimental = feature_for_test("experimental", FeatureStatus::Experimental, &[]);
    let supported = feature_for_test("supported", FeatureStatus::Supported, &["experimental"]);
    let error = validate_feature_set(&[experimental, supported])
        .expect_err("supported features should require supported dependencies");
    assert!(error
        .to_string()
        .contains("`supported` features may depend only on `supported` features"));

    let planned = feature_for_test("planned", FeatureStatus::Planned, &[]);
    let experimental = feature_for_test("experimental", FeatureStatus::Experimental, &["planned"]);
    let error = validate_feature_set(&[planned, experimental])
        .expect_err("experimental features should not depend on planned work");
    assert!(error.to_string().contains(
        "`experimental` features may depend only on `experimental` or `supported` features"
    ));

    let supported = feature_for_test("supported", FeatureStatus::Supported, &[]);
    let experimental =
        feature_for_test("experimental", FeatureStatus::Experimental, &["supported"]);
    let deprecated = feature_for_test("deprecated", FeatureStatus::Deprecated, &["experimental"]);
    let planned = feature_for_test("planned", FeatureStatus::Planned, &["deprecated"]);
    validate_feature_set(&[supported, experimental, deprecated, planned])
        .expect("each status should accept its documented dependency maturity");
}

#[test]
fn feature_dependency_layers_are_deterministic() {
    let features = vec![
        feature_for_test("leaf", FeatureStatus::Supported, &["right", "left"]),
        feature_for_test("right", FeatureStatus::Supported, &["root"]),
        feature_for_test("root", FeatureStatus::Supported, &[]),
        feature_for_test("left", FeatureStatus::Supported, &["root"]),
    ];

    validate_feature_set(&features).expect("feature graph should be valid");
    let layers = feature_dependency_layers(&features).expect("layers should resolve");
    assert_eq!(
        layers
            .iter()
            .map(|layer| {
                layer
                    .iter()
                    .map(|feature| feature.id.as_str())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>(),
        vec![vec!["root"], vec!["left", "right"], vec!["leaf"]]
    );
}

#[test]
fn render_dashboard_is_stable_and_uses_compact_validation_cells() {
    let features = vec![
        Feature {
            id: "a.feature".to_owned(),
            title: "Aye".to_owned(),
            area: "core".to_owned(),
            domains: vec![FeatureDomain::SmallMolecule, FeatureDomain::Macromolecule],
            version: 1,
            status: FeatureStatus::Supported,
            description: "A feature.".to_owned(),
            depends_on: Vec::new(),
            validation_required: Vec::new(),
        },
        Feature {
            id: "z.feature".to_owned(),
            title: "Zed".to_owned(),
            area: "io".to_owned(),
            domains: vec![FeatureDomain::SmallMolecule],
            version: 3,
            status: FeatureStatus::Supported,
            description: "Z feature.".to_owned(),
            depends_on: vec!["a.feature".to_owned()],
            validation_required: Vec::new(),
        },
        Feature {
            id: "failing.feature".to_owned(),
            title: "Failing".to_owned(),
            area: "validation".to_owned(),
            domains: vec![FeatureDomain::SmallMolecule],
            version: 1,
            status: FeatureStatus::Deprecated,
            description: "Feature with counted failures.".to_owned(),
            depends_on: Vec::new(),
            validation_required: vec!["pubchem-1k".to_owned()],
        },
        Feature {
            id: "missing.feature".to_owned(),
            title: "Missing".to_owned(),
            area: "validation".to_owned(),
            domains: vec![FeatureDomain::Macromolecule],
            version: 1,
            status: FeatureStatus::Experimental,
            description: "Feature without recorded status.".to_owned(),
            depends_on: Vec::new(),
            validation_required: vec!["pdb-100".to_owned()],
        },
        Feature {
            id: "harness.feature".to_owned(),
            title: "Harness".to_owned(),
            area: "infrastructure".to_owned(),
            domains: vec![FeatureDomain::Infrastructure],
            version: 2,
            status: FeatureStatus::Supported,
            description: "Infrastructure feature.".to_owned(),
            depends_on: Vec::new(),
            validation_required: Vec::new(),
        },
    ];
    let statuses = BTreeMap::from([(
        "failing.feature".to_owned(),
        ValidationStatus {
            feature_id: "failing.feature".to_owned(),
            corpora: BTreeMap::from([(
                "pubchem-1k".to_owned(),
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
                kind: CorpusKind::Mixed,
                expected_count: 7,
                features: BTreeMap::from([
                    (
                        "failing.feature".to_owned(),
                        CorpusFeatureDashboardInfo {
                            reference_tool: "rdkit".to_owned(),
                            reference_version: "RDKit 2026.03.3".to_owned(),
                        },
                    ),
                    (
                        "missing.feature".to_owned(),
                        CorpusFeatureDashboardInfo {
                            reference_tool: "biopython".to_owned(),
                            reference_version: "Biopython 1.87 / mkdssp version 4.6.1".to_owned(),
                        },
                    ),
                ]),
            },
        ),
        (
            "pubchem-1k".to_owned(),
            CorpusDashboardInfo {
                id: "pubchem-1k".to_owned(),
                label: "PubChem 1k".to_owned(),
                title: "PubChem deterministic 1000-compound corpus".to_owned(),
                kind: CorpusKind::SmallMolecule,
                expected_count: 1000,
                features: BTreeMap::from([
                    (
                        "a.feature".to_owned(),
                        CorpusFeatureDashboardInfo {
                            reference_tool: "rdkit".to_owned(),
                            reference_version: "RDKit 2026.03.3".to_owned(),
                        },
                    ),
                    (
                        "failing.feature".to_owned(),
                        CorpusFeatureDashboardInfo {
                            reference_tool: "rdkit".to_owned(),
                            reference_version: "RDKit 2026.03.3".to_owned(),
                        },
                    ),
                ]),
            },
        ),
        (
            "pdb-100".to_owned(),
            CorpusDashboardInfo {
                id: "pdb-100".to_owned(),
                label: "PDB 100".to_owned(),
                title: "PDB deterministic 100-entry corpus".to_owned(),
                kind: CorpusKind::Macromolecule,
                expected_count: 100,
                features: BTreeMap::from([(
                    "missing.feature".to_owned(),
                    CorpusFeatureDashboardInfo {
                        reference_tool: "biopython".to_owned(),
                        reference_version: "Biopython 1.87 / mkdssp version 4.6.1".to_owned(),
                    },
                )]),
            },
        ),
    ]);

    let dashboard = render_dashboard(&features, &statuses, &corpus_info);

    assert!(dashboard.starts_with("<!doctype html>\n"));
    assert!(
        dashboard.contains("<table id=\"small-molecules-dashboard\" class=\"feature-dashboard\">")
    );
    assert!(
        dashboard.contains("<table id=\"macromolecules-dashboard\" class=\"feature-dashboard\">")
    );
    assert!(dashboard.contains(
        "<table id=\"infrastructure-dashboard\" class=\"feature-dashboard infrastructure-table\">"
    ));
    assert!(dashboard.contains("<h2>Small molecules</h2>"));
    assert!(dashboard.contains("<h2>Macromolecules</h2>"));
    assert!(dashboard.contains("<h2>Infrastructure and harness</h2>"));
    assert!(dashboard.contains("<h2>Feature dependency graph</h2>"));
    let infrastructure_position = dashboard
        .find("<h2>Infrastructure and harness</h2>")
        .expect("infrastructure section should be present");
    let graph_position = dashboard
        .find("<h2>Feature dependency graph</h2>")
        .expect("dependency graph should be present");
    assert!(
        infrastructure_position < graph_position,
        "all feature tables should precede the dependency graph"
    );
    assert!(dashboard.contains("class=\"feature-graph\""));
    assert!(dashboard.contains("marker-end=\"url(#feature-graph-arrow)\""));
    assert!(dashboard.contains("<a href=\"./z.feature/feature.md\">"));
    assert!(dashboard.contains("layer 0"));
    assert!(dashboard.contains("layer 1"));
    assert!(dashboard.contains("<strong>Reference codebase:</strong> RDKit v2026.03.3"));
    assert!(dashboard.contains("<strong>Reference codebase:</strong> Biopython v1.87"));
    assert!(dashboard.contains("<strong>DSSP executable:</strong> mkdssp v4.6.1"));
    assert!(dashboard.contains("th.area, td.area { text-align: left; }"));
    assert!(dashboard.contains("<th class=\"compact area\" data-sort-type=\"text\" title=\"Area\"><button class=\"sort\" type=\"button\" aria-label=\"Sort by Area\">Area</button></th>"));
    assert!(dashboard.contains("<td class=\"compact area\" data-sort-value=\"core\">core</td>"));
    assert!(!dashboard
        .contains("aria-label=\"Sort by Area\"><span class=\"rotated-label\">Area</span>"));
    assert!(dashboard.contains("aria-label=\"Sort by Status\">Status</button>"));
    assert!(dashboard.contains("<span class=\"feature-status status-supported\">supported</span>"));
    assert!(dashboard
        .contains("<span class=\"feature-status status-experimental\">experimental</span>"));
    assert!(
        dashboard.contains("<span class=\"feature-status status-deprecated\">deprecated</span>")
    );
    assert!(!dashboard.contains(">Implemented<"));
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
    assert!(!dashboard.contains("<span class=\"rotated-name\">smoke</span>"));
    assert!(dashboard.contains(
        "<span class=\"rotated-name\">pubchem-1k</span><br><span class=\"rotated-count\">(n=1000)</span>"
    ));
    assert!(dashboard.contains(
        "<span class=\"rotated-name\">pdb-100</span><br><span class=\"rotated-count\">(n=100)</span>"
    ));
    assert_eq!(dashboard.matches("<code>a.feature</code>").count(), 2);
    assert!(dashboard.contains("data-sort-value=\"0\""));
    assert!(dashboard.contains("<code>z.feature</code>"));
    assert!(dashboard.contains("<code>harness.feature</code>"));
    assert!(dashboard.contains("data-sort-value=\"1\""));
    assert!(dashboard.contains("aria-label=\"failed: 3 non-passing case(s)\""));
    assert!(dashboard.contains("<span class=\"count\">3</span>"));
    assert!(dashboard.contains("<span class=\"unknown\">?</span>unknown"));
    assert!(dashboard.contains(
        "<span class=\"unknown\" aria-label=\"unknown\" title=\"no recorded validation status; reference: Biopython v1.87"
    ));
    assert!(!dashboard.contains(
        "<span class=\"bad\" aria-label=\"failed\" title=\"no recorded validation status\">"
    ));
    assert!(dashboard.contains("document.querySelectorAll('table.feature-dashboard')"));
    assert!(dashboard.contains("button.addEventListener('click'"));
    assert!(dashboard.ends_with('\n'));
}

#[test]
fn dashboard_corpus_cells_show_optional_manifest_evidence() {
    let feature = Feature {
        id: "optional.feature".to_owned(),
        title: "Optional".to_owned(),
        area: "validation".to_owned(),
        domains: vec![FeatureDomain::SmallMolecule],
        version: 1,
        status: FeatureStatus::Supported,
        description: "Feature with optional corpus evidence.".to_owned(),
        depends_on: Vec::new(),
        validation_required: Vec::new(),
    };
    let passed = CorpusStatus {
        passed: true,
        fixture_count: 1,
        compared_count: 1,
        reference_tool: "rdkit".to_owned(),
        reference_version: "RDKit test".to_owned(),
        manifest_hash: "a".repeat(64),
        failed_count: 0,
        first_failure: None,
        evidence_schema_version: Some(VALIDATION_EVIDENCE_SCHEMA_VERSION),
        evidence_hash: Some("b".repeat(64)),
        evidence_input_count: 1,
        legacy_evidence_inputs: Vec::new(),
        validated_at_unix: 1,
    };
    let status = ValidationStatus {
        feature_id: feature.id.clone(),
        corpora: BTreeMap::from([("pubchem-1k".to_owned(), passed)]),
    };
    let reference = CorpusFeatureDashboardInfo {
        reference_tool: "rdkit".to_owned(),
        reference_version: "RDKit test".to_owned(),
    };

    assert!(dashboard_corpus_cell(
        &feature,
        Some(&status),
        "pubchem-1k",
        Some(&reference),
        true,
    )
    .contains("aria-label=\"passed\""));
    assert!(
        dashboard_corpus_cell(&feature, None, "pubchem-1k", Some(&reference), true)
            .contains("title=\"no recorded validation status; reference: RDKit vtest\"")
    );
    assert!(
        dashboard_corpus_cell(&feature, None, "pubchem-1k", None, true)
            .contains("aria-label=\"not required\"")
    );
    assert!(
        dashboard_corpus_cell(&feature, Some(&status), "pubchem-1k", None, true)
            .contains("aria-label=\"not required\"")
    );
    assert!(dashboard_corpus_cell(
        &feature,
        Some(&status),
        "pubchem-1k",
        Some(&reference),
        false
    )
    .contains("aria-label=\"not required\""));
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
Use feature.md. Set status = "supported" with evidence, declare depends_on, and declare validation_required.
Molecular validation fixtures must be externally supplied.
Run cargo xtask dashboard --check and cargo xtask validate --feature <feature-id> --corpus <corpus-id>.
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
Read feature.md. Run cargo test --workspace and cargo xtask validate --feature <feature-id> --corpus <corpus-id>.
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
fn validate_jobs_uses_a_memory_safe_default_and_accepts_override() {
    let default_jobs = validation_jobs(&[]).expect("default worker count should resolve");
    assert!((1..=4).contains(&default_jobs));
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
    assert!(validate_args(&[
        "--feature".to_owned(),
        "io.smiles.canonical".to_owned(),
        "--corpus".to_owned(),
        "pubchem-100k".to_owned(),
        "--fixture".to_owned(),
        "data/packs/pack_001.smi".to_owned(),
    ])
    .is_ok());
    assert!(validate_args(&[
        "--feature".to_owned(),
        "stereo.perception".to_owned(),
        "--corpus".to_owned(),
        "pubchem-100k".to_owned(),
        "--accept-implementation-goldens".to_owned(),
    ])
    .is_ok());
}

#[test]
fn implementation_golden_acceptance_is_limited_to_manual_semantic_references() {
    let root = temp_feature_root("accept-implementation-goldens");
    let corpus_root = root.join("validation/corpora/smoke");
    let manifest_path = corpus_root.join("features/stereo.perception.toml");
    fs::create_dir_all(manifest_path.parent().expect("manifest parent"))
        .expect("features directory");
    fs::create_dir_all(corpus_root.join("data")).expect("data directory");
    fs::write(corpus_root.join("data/example.smi"), "CC CID:1\n").expect("fixture should write");
    let mut manifest = ValidationManifest {
        feature_id: "stereo.perception".to_owned(),
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
    assert_eq!(validation_worker_count(16, 3), 3);
    assert_eq!(validation_worker_count(0, 3), 1);
}

#[test]
fn validation_defaults_to_all_and_all_includes_available_manifest_backed_features() {
    assert_eq!(validation_corpus_selector(&[]), "all");
    assert_eq!(
        validation_corpus_selector(&[
            "--feature".to_owned(),
            "all".to_owned(),
            "--corpus".to_owned(),
            "pubchem-1k".to_owned(),
        ]),
        "pubchem-1k"
    );

    let root = temp_feature_root("all-validation-corpora");
    let features = vec![
        Feature {
            id: "small".to_owned(),
            title: "Small".to_owned(),
            area: "io".to_owned(),
            domains: vec![FeatureDomain::SmallMolecule],
            version: 1,
            status: FeatureStatus::Supported,
            description: "Small feature.".to_owned(),
            depends_on: Vec::new(),
            validation_required: vec!["pubchem-1k".to_owned()],
        },
        Feature {
            id: "macro".to_owned(),
            title: "Macro".to_owned(),
            area: "bio".to_owned(),
            domains: vec![FeatureDomain::Macromolecule],
            version: 1,
            status: FeatureStatus::Experimental,
            description: "Macro feature.".to_owned(),
            depends_on: Vec::new(),
            validation_required: vec!["pdb-100".to_owned()],
        },
        Feature {
            id: "planned".to_owned(),
            title: "Planned".to_owned(),
            area: "descriptors".to_owned(),
            domains: vec![FeatureDomain::SmallMolecule],
            version: 1,
            status: FeatureStatus::Planned,
            description: "Planned feature.".to_owned(),
            depends_on: Vec::new(),
            validation_required: Vec::new(),
        },
        Feature {
            id: "deprecated".to_owned(),
            title: "Deprecated".to_owned(),
            area: "descriptors".to_owned(),
            domains: vec![FeatureDomain::SmallMolecule],
            version: 1,
            status: FeatureStatus::Deprecated,
            description: "Deprecated feature.".to_owned(),
            depends_on: Vec::new(),
            validation_required: Vec::new(),
        },
    ];
    for (feature, corpus) in [
        ("small", "pubchem-1k"),
        ("small", "pubchem-100k"),
        ("small", "enamine-diversity"),
        ("macro", "pdb-100"),
        ("macro", "pdb-1000"),
        ("planned", "pubchem-1k"),
        ("deprecated", "pubchem-100k"),
    ] {
        let path = validation_manifest_path_from(&root, feature, corpus);
        fs::create_dir_all(path.parent().expect("manifest parent"))
            .expect("manifest directory should create");
        fs::write(path, "").expect("manifest marker should write");
    }

    assert_eq!(
        validation_targets_from(&root, &features, "all", "pubchem-1k")
            .into_iter()
            .map(|(feature, corpus)| (feature.id.as_str(), corpus))
            .collect::<Vec<_>>(),
        vec![("small", "pubchem-1k".to_owned())]
    );
    assert_eq!(
        validation_targets_from(&root, &features, "all", "all")
            .into_iter()
            .map(|(feature, corpus)| (feature.id.as_str(), corpus))
            .collect::<Vec<_>>(),
        vec![
            ("small", "pubchem-1k".to_owned()),
            ("small", "pubchem-100k".to_owned()),
            ("small", "enamine-diversity".to_owned()),
            ("macro", "pdb-100".to_owned()),
            ("macro", "pdb-1000".to_owned()),
            ("deprecated", "pubchem-100k".to_owned()),
        ]
    );
    assert_eq!(
        validation_targets_from(&root, &features, "small", "all")
            .into_iter()
            .map(|(feature, corpus)| (feature.id.as_str(), corpus))
            .collect::<Vec<_>>(),
        vec![
            ("small", "pubchem-1k".to_owned()),
            ("small", "pubchem-100k".to_owned()),
            ("small", "enamine-diversity".to_owned()),
        ]
    );
    assert_eq!(
        validation_targets_from(&root, &features, "small", "pubchem-100k")
            .into_iter()
            .map(|(feature, corpus)| (feature.id.as_str(), corpus))
            .collect::<Vec<_>>(),
        vec![("small", "pubchem-100k".to_owned())]
    );
    assert_eq!(
        validation_targets_from(&root, &features, "macro", "pubchem-1k")
            .into_iter()
            .map(|(feature, corpus)| (feature.id.as_str(), corpus))
            .collect::<Vec<_>>(),
        vec![("macro", "pubchem-1k".to_owned())]
    );
    fs::remove_dir_all(root).ok();
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
        let expected = implementation_expected(feature, "pubchem-1k", &fixture)
            .expect("feature should compare");
        assert_eq!(expected["records"][0]["status"], "ok");
    }

    fs::remove_dir_all(root).ok();
}

#[test]
fn implementation_dispatch_supports_mmcif_document_rows() {
    let root = temp_feature_root("mmcif-document-dispatch");
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
        .expect("mmCIF document feature should compare");
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
fn implementation_dispatch_supports_hydrogen_normalization() {
    let root = temp_feature_root("hydrogen-normalization-dispatch");
    let fixture = root.join("fixture.sdf");
    fs::write(&fixture, simple_sdf_record("methane")).expect("fixture should write");

    let expected = implementation_expected("chem.hydrogen-normalization", "pubchem-1k", &fixture)
        .expect("feature should compare");
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
fn implementation_dispatch_supports_query_validation() {
    let root = temp_feature_root("query-validation-dispatch");
    let smarts_fixture = root.join("fixture.smi");
    fs::write(&smarts_fixture, "CCO\nC1=CC=CC=C1\n").expect("fixture should write");

    let parsed = implementation_expected("query.smarts", "pubchem-1k", &smarts_fixture)
        .expect("SMARTS feature should compare");
    assert_eq!(parsed["records"][0]["status"], "ok");
    assert_eq!(parsed["records"][0]["atom_count"], 3);
    assert_eq!(parsed["records"][1]["bond_count"], 6);

    let molecule_fixture = root.join("fixture.sdf");
    fs::write(&molecule_fixture, simple_sdf_record("methane")).expect("fixture should write");
    let matched = implementation_expected("algo.substructure.vf2", "pubchem-1k", &molecule_fixture)
        .expect("substructure feature should compare");
    assert_eq!(matched["records"][0]["status"], "ok");
    assert_eq!(matched["records"][0]["queries"][0]["smarts"], "[#6]");
    assert_eq!(matched["records"][0]["queries"][0]["matches"], json!([[0]]));

    fs::remove_dir_all(root).ok();
}

#[test]
fn implementation_dispatch_uses_current_isomeric_smiles_feature_id() {
    let root = temp_feature_root("isomeric-smiles-feature-dispatch");
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
        .expect("feature should compare");
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
fn nonisomeric_smiles_validation_excludes_stereo_syntax() {
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
fn stereo_and_nonisomeric_validation_use_distinct_smiles_subsets() {
    let root = temp_feature_root("smiles-validation-subsets");
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

    let expected = implementation_expected("stereo.cip", "pubchem-1k", &fixture)
        .expect("feature should compare");
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

    let expected = implementation_expected("stereo.cip", "pubchem-1k", &fixture)
        .expect("feature should compare");
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

    let expected = implementation_expected("stereo.cip", "pubchem-1k", &fixture)
        .expect("feature should compare");
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
fn recorded_corpus_status_requires_known_nonempty_evidence() {
    let root = temp_feature_root("status-rejects-stale");
    let (_, _, manifest_path) = write_evidence_test_repo(&root);
    let manifest = read_validation_manifest(&manifest_path).expect("manifest should read");
    let evidence =
        build_validation_evidence(&root, &manifest_path, &manifest).expect("evidence should build");
    let mut corpus_status = CorpusStatus {
        passed: true,
        fixture_count: 1,
        compared_count: 1,
        reference_tool: "rdkit".to_owned(),
        reference_version: "RDKit test".to_owned(),
        manifest_hash: hash_evidence_file(&manifest_path).expect("manifest should hash"),
        failed_count: 0,
        first_failure: None,
        evidence_schema_version: Some(VALIDATION_EVIDENCE_SCHEMA_VERSION),
        evidence_hash: Some(evidence.sha256),
        evidence_input_count: evidence.inputs.len(),
        legacy_evidence_inputs: Vec::new(),
        validated_at_unix: 1,
    };
    assert!(recorded_corpus_status_passed(Some(&corpus_status)));

    corpus_status.evidence_schema_version = Some(999);
    assert!(!recorded_corpus_status_passed(Some(&corpus_status)));

    corpus_status.evidence_schema_version = Some(VALIDATION_EVIDENCE_SCHEMA_VERSION);
    corpus_status.compared_count = 0;
    assert!(!recorded_corpus_status_passed(Some(&corpus_status)));

    corpus_status.compared_count = corpus_status.fixture_count;
    corpus_status.evidence_input_count = 0;
    assert!(!recorded_corpus_status_passed(Some(&corpus_status)));
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
fn dashboard_status_requires_a_current_manifest() {
    let feature = Feature {
        id: "portable.feature".to_owned(),
        title: "Portable".to_owned(),
        area: "infrastructure".to_owned(),
        domains: vec![FeatureDomain::Infrastructure],
        version: 1,
        status: FeatureStatus::Supported,
        description: "Portable dashboard evidence.".to_owned(),
        depends_on: Vec::new(),
        validation_required: vec!["pubchem-1k".to_owned()],
    };
    let status = ValidationStatus {
        feature_id: feature.id.clone(),
        corpora: BTreeMap::from([(
            "pubchem-1k".to_owned(),
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
                evidence_input_count: 1,
                legacy_evidence_inputs: Vec::new(),
                validated_at_unix: 1,
            },
        )]),
    };
    let reference = CorpusFeatureDashboardInfo {
        reference_tool: "rdkit".to_owned(),
        reference_version: "test".to_owned(),
    };

    assert!(
        dashboard_corpus_cell(&feature, Some(&status), "pubchem-1k", None, true)
            .contains("no recorded validation status")
    );
    assert!(dashboard_corpus_cell(
        &feature,
        Some(&status),
        "pubchem-1k",
        Some(&reference),
        true,
    )
    .contains("aria-label=\"passed\""));
}
#[test]
fn status_writer_prunes_entries_without_manifests() {
    let root = temp_feature_root("status-manifest-pruning");
    let status_path = validation_status_path_from(&root, "pubchem-1k");
    fs::create_dir_all(status_path.parent().expect("status parent"))
        .expect("status directory should create");
    fs::write(&status_path, "stale").expect("stale status should write");

    let feature_status = CorpusStatus::from_failed_run(FailedValidationRun {
        fixture_count: 1,
        compared_count: 0,
        failed_count: 1,
        first_failure: "fixture differs".to_owned(),
        reference_tool: "rdkit".to_owned(),
        reference_version: "test".to_owned(),
        manifest_hash: "0".repeat(64),
    })
    .expect("failed status should build");
    let statuses = BTreeMap::from([(
        "example.feature".to_owned(),
        ValidationStatus {
            feature_id: "example.feature".to_owned(),
            corpora: BTreeMap::from([("pubchem-1k".to_owned(), feature_status)]),
        },
    )]);
    let selected = BTreeSet::from(["pubchem-1k".to_owned()]);

    write_validation_statuses_from(&root, &statuses, &selected)
        .expect("status pruning should succeed");
    assert!(!status_path.exists());

    let manifest = validation_manifest_path_from(&root, "example.feature", "pubchem-1k");
    fs::create_dir_all(manifest.parent().expect("manifest parent"))
        .expect("manifest directory should create");
    fs::write(&manifest, "").expect("manifest marker should write");
    write_validation_statuses_from(&root, &statuses, &selected)
        .expect("manifest-backed status should write");
    let stored = read_corpus_status(&status_path).expect("written status should parse");
    assert!(stored.features.contains_key("example.feature"));

    fs::remove_dir_all(root).ok();
}
#[test]
fn old_status_toml_defaults_failure_summary_fields() {
    let root = temp_feature_root("legacy-status");
    let path = root.join("status.toml");
    fs::write(
        &path,
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
validated_at_unix = 1

[[features.example.evidence_inputs]]
path = "validation/corpora/smoke/data/example.sdf"
sha256 = "2222222222222222222222222222222222222222222222222222222222222222"
"#,
    )
    .expect("legacy status fixture should write");
    let status = read_corpus_status(&path).expect("old status shape should deserialize");
    let corpus_status = status
        .features
        .get("example")
        .expect("feature status should exist");

    assert_eq!(corpus_status.failed_count, 0);
    assert_eq!(corpus_status.first_failure, None);
    assert_eq!(corpus_status.evidence_input_count, 1);
    assert_eq!(corpus_status.legacy_evidence_inputs.len(), 1);

    let compact = toml::to_string(&status).expect("legacy status should reserialize");
    assert!(compact.contains("evidence_input_count = 1"));
    assert!(!compact.contains("evidence_inputs"));
    fs::remove_dir_all(root).ok();
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
    let molecule =
        SmallMolecule::from_smiles("c1cccc1").expect("invalid aromatic molecule should parse");
    let mut record = IndexedSmallRecord {
        record_index: 0,
        title: "invalid aromatic representation".to_owned(),
        molecule,
        sdf_fields: BTreeMap::new(),
    };

    let value = stereo_perception_record_json(&mut record);

    assert_eq!(
        value.get("status").and_then(Value::as_str),
        Some("sanitize_error")
    );
    assert!(value.get("report").is_none());
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
    normalize_feature_for_comparison_in_place("bio.secondary-structure.dssp", &mut expected);
    normalize_feature_for_comparison_in_place("bio.secondary-structure.dssp", &mut actual);
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
    let single = SmallMolecule::from_smiles("CC").expect("single bond should parse");
    let double = SmallMolecule::from_smiles("C=C").expect("double bond should parse");
    assert_ne!(
        smiles_sanitized_bonds_json(single.graph()),
        smiles_sanitized_bonds_json(double.graph())
    );

    let aromatic = SmallMolecule::from_smiles("c1ccccc1").expect("benzene should parse");
    let mut sanitized_aromatic = aromatic.clone();
    perception::sanitize_with_options(&mut sanitized_aromatic, SanitizeOptions::default())
        .expect("benzene should sanitize");
    assert_eq!(
        explicit_valence_json(sanitized_aromatic.graph(), AtomId::new(0)),
        3
    );
    let mut aromatic_cyclohexyne =
        SmallMolecule::from_smiles("C1=CC#CC=C1").expect("cyclohexyne parses");
    perception::sanitize_with_options(&mut aromatic_cyclohexyne, SanitizeOptions::default())
        .expect("cyclohexyne should sanitize");
    let alkyne_atoms = aromatic_cyclohexyne
        .graph()
        .bonds()
        .find_map(|(id, bond)| {
            (aromatic_cyclohexyne
                .graph()
                .bond_is_aromatic(id)
                .ok()
                .flatten()
                == Some(true)
                && bond.order == BondOrder::Triple)
                .then_some(bond.endpoints())
        })
        .expect("aromaticized triple bond is retained");
    assert_eq!(
        explicit_valence_json(aromatic_cyclohexyne.graph(), alkyne_atoms.0),
        4
    );
    assert_eq!(
        explicit_valence_json(aromatic_cyclohexyne.graph(), alkyne_atoms.1),
        4
    );
    let mut thiophene = SmallMolecule::from_smiles("c1ccsc1").expect("thiophene parses");
    perception::sanitize_with_options(&mut thiophene, SanitizeOptions::default())
        .expect("thiophene should sanitize");
    let sulfur_id = thiophene
        .graph()
        .atoms()
        .find_map(|(id, atom)| (atom.element.symbol() == "S").then_some(id))
        .expect("sulfur atom");
    assert_eq!(explicit_valence_json(thiophene.graph(), sulfur_id), 2);
    let mut phosphorus_ring =
        SmallMolecule::from_smiles("C(F)(F)(F)P1P(P(P(P1C(F)(F)F)C(F)(F)F)C(F)(F)F)C(F)(F)F")
            .expect("phosphorus ring parses");
    perception::sanitize_with_options(&mut phosphorus_ring, SanitizeOptions::default())
        .expect("phosphorus ring should sanitize");
    for (phosphorus_id, _phosphorus) in phosphorus_ring
        .graph()
        .atoms()
        .filter(|(_, atom)| atom.element.symbol() == "P")
    {
        assert_eq!(
            phosphorus_ring
                .graph()
                .atom_is_aromatic(phosphorus_id)
                .unwrap(),
            Some(true)
        );
        assert_eq!(
            explicit_valence_json(phosphorus_ring.graph(), phosphorus_id),
            3
        );
    }
    let mut phosphinine = SmallMolecule::from_smiles("C1=CC=PC=C1").expect("phosphinine parses");
    perception::sanitize_with_options(&mut phosphinine, SanitizeOptions::default())
        .expect("phosphinine should sanitize");
    let phosphinine_phosphorus = phosphinine
        .graph()
        .atoms()
        .find_map(|(id, atom)| (atom.element.symbol() == "P").then_some(id))
        .expect("phosphinine phosphorus");
    assert_eq!(
        explicit_valence_json(phosphinine.graph(), phosphinine_phosphorus),
        3
    );
    let mut anionic_macrocycle = SmallMolecule::from_smiles("CN(C)CCO.C1=CC=C2C(=C1)C3=NC4=C5C=CC=CC5=C([N-]4)N=C6C7=CC=CC=C7C(=N6)N=C8C9=CC=CC=C9C(=N8)N=C2[N-]3.[Cu+2]")
    .expect("anionic macrocycle parses");
    perception::sanitize_with_options(&mut anionic_macrocycle, SanitizeOptions::default())
        .expect("anionic macrocycle should sanitize");
    let anionic_nitrogen = anionic_macrocycle
        .graph()
        .atoms()
        .find_map(|(id, atom)| {
            (atom.element.symbol() == "N"
                && atom.formal_charge < 0
                && anionic_macrocycle
                    .graph()
                    .atom_is_aromatic(id)
                    .ok()
                    .flatten()
                    == Some(true))
            .then_some(id)
        })
        .expect("anionic aromatic nitrogen");
    assert_eq!(
        explicit_valence_json(anionic_macrocycle.graph(), anionic_nitrogen),
        2
    );
    let mut cyclopentadienyl = SmallMolecule::from_smiles("[CH-]1[C-]=[C-][C-]=[C-]1")
        .expect("cyclopentadienyl anion parses");
    perception::sanitize_with_options(&mut cyclopentadienyl, SanitizeOptions::default())
        .expect("cyclopentadienyl anion should sanitize");
    let anionic_carbon_with_h = cyclopentadienyl
        .graph()
        .atoms()
        .find_map(|(id, atom)| {
            (atom.element.symbol() == "C"
                && atom.formal_charge < 0
                && cyclopentadienyl.graph().atom_is_aromatic(id).ok().flatten() == Some(true)
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
    let mut substituted_cyclopentadienyl = SmallMolecule::from_smiles("C[C-]1[C-]=[C-][C-]=[C-]1")
        .expect("substituted cyclopentadienyl parses");
    perception::sanitize_with_options(
        &mut substituted_cyclopentadienyl,
        SanitizeOptions::default(),
    )
    .expect("substituted cyclopentadienyl should sanitize");
    let substituted_anionic_carbon = substituted_cyclopentadienyl
        .graph()
        .atoms()
        .find_map(|(id, atom)| {
            let degree = substituted_cyclopentadienyl
                .graph()
                .incident_bonds(id)
                .ok()?
                .count();
            (atom.element.symbol() == "C"
                && atom.formal_charge < 0
                && substituted_cyclopentadienyl
                    .graph()
                    .atom_is_aromatic(id)
                    .ok()
                    .flatten()
                    == Some(true)
                && degree == 3)
                .then_some(id)
        })
        .expect("substituted anionic carbon");
    assert_eq!(
        explicit_valence_json(
            substituted_cyclopentadienyl.graph(),
            substituted_anionic_carbon,
        ),
        3
    );
    let mut fused_triazine = SmallMolecule::from_smiles(
        "O=[N+]([O-])c2cc(-c1nn5c(=O)c(C=Cc3c(O)ccc4c3cccc4)nnc5s1)ccc2",
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
                .filter(|(bond, _)| {
                    fused_triazine
                        .graph()
                        .bond_is_aromatic(*bond)
                        .ok()
                        .flatten()
                        == Some(true)
                })
                .count();
            (atom.element.symbol() == "N"
                && fused_triazine.graph().atom_is_aromatic(id).ok().flatten() == Some(true)
                && aromatic_degree >= 3)
                .then_some(id)
        })
        .expect("tri-coordinate aromatic nitrogen");
    assert_eq!(
        explicit_valence_json(fused_triazine.graph(), tricoordinate_aromatic_nitrogen),
        3
    );
    assert!(smiles_sanitized_bonds_json(aromatic.graph())
        .iter()
        .all(|bond| bond["bond_type"] == "AROMATIC" && bond["is_aromatic"] == true));

    let labeled = SmallMolecule::from_smiles("[13CH3:7]C").expect("labeled carbon should parse");
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
fn canonical_smiles_validation_matches_rdkit_parse_status_for_unsanitizable_input() {
    let root = temp_feature_root("canonical-unsanitizable-input");
    let fixture = root.join("fixture.smi");
    fs::write(&fixture, "[Cl-](Br)Br CID:invalid\n").expect("fixture should write");

    let records = read_canonical_smiles_records(&fixture).expect("records should load");
    let item =
        canonical_smiles_record_json(&records[0], false).expect("canonical record should render");

    assert_eq!(item["status"], "parse_error");
}

#[test]
fn smiles_semantics_match_rdkit_aromatic_carbonyl_valence() {
    let molecule = SmallMolecule::from_smiles("CCCCCCCc1cc2c(=O)ccn(O)c2cc1")
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
    let molecule =
        SmallMolecule::from_smiles("[nH]1cccc1").expect("aromatic nH SMILES should parse");

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
    let molecule = SmallMolecule::from_smiles("CCOC(=O)C1=C(C(=C(N1)C)C(=O)OC(C)(C)C)C")
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

fn feature_for_test(id: &str, status: FeatureStatus, depends_on: &[&str]) -> Feature {
    Feature {
        id: id.to_owned(),
        title: id.to_owned(),
        area: "test".to_owned(),
        domains: vec![FeatureDomain::Infrastructure],
        version: 1,
        status,
        description: format!("Test feature {id}."),
        depends_on: depends_on
            .iter()
            .map(|dependency| (*dependency).to_owned())
            .collect(),
        validation_required: Vec::new(),
    }
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
        "id = \"example\"\ntitle = \"Example\"\narea = \"test\"\ndomains = [\"small-molecule\"]\nversion = 1\nstatus = \"supported\"\ndescription = \"Example feature.\"\ndepends_on = []\nvalidation_required = [\"smoke\"]\n",
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
