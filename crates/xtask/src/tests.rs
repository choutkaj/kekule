use super::*;

fn normalize_and_perceive(molecule: &mut SmallMolecule) {
    molecule.normalize().expect("molecule should normalize");
    molecule.perceive().expect("molecule should be perceived");
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
            .join("benchmarks/corpora")
            .join(corpus_id)
            .join("corpus.toml");
        let text = fs::read_to_string(&path).expect("corpus descriptor should read");
        let descriptor: CorpusDescriptor =
            toml::from_str(&text).expect("corpus descriptor should parse");
        let registered = benchmark_corpus(corpus_id).expect("corpus should be registered");
        assert_eq!(descriptor.id, registered.id);
        assert_eq!(descriptor.local_only, registered.local_only);
    }
}

#[test]
fn kekule_package_metadata_uses_the_initial_release_contract() {
    let workspace_root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let workspace_manifest: toml::Value = toml::from_str(
        &fs::read_to_string(workspace_root.join("Cargo.toml"))
            .expect("workspace manifest should read"),
    )
    .expect("workspace manifest should parse");
    let workspace_package = &workspace_manifest["workspace"]["package"];
    assert_eq!(workspace_package["version"].as_str(), Some("0.1.0"));
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
            Some("0.1.0")
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
fn feature_metadata_does_not_require_benchmark_manifests() {
    let root = temp_feature_root("optional-benchmark-manifest");
    write_feature(
        &root,
        "valid.feature",
        r#"id = "valid.feature"
title = "Valid feature"
area = "infrastructure"
domains = ["infrastructure"]
version = 1
status = "supported"
description = "Feature without a benchmark manifest."
depends_on = []
"#,
    );

    let feature = read_feature(&root.join("valid.feature").join("feature.toml"))
        .expect("feature should not require a benchmark manifest");
    assert!(benchmark_targets_from(&root, std::slice::from_ref(&feature), "all", "all").is_empty());
    assert!(BENCHMARK_CORPORA.iter().all(|corpus| corpus.local_only));
    assert!(benchmark_corpus("smoke").is_none());
    assert!(benchmark_corpus("pubchem-100").is_none());
    assert!(benchmark_corpus("pdb-10").is_none());
    fs::remove_dir_all(root).ok();
}

#[test]
fn feature_schema_rejects_removed_validation_and_benchmark_requirement_keys() {
    let root = temp_feature_root("removed-benchmark-requirement");
    for removed_key in ["validation_required", "benchmark_required"] {
        write_feature(
            &root,
            removed_key,
            &format!(
                r#"id = "{removed_key}"
title = "Removed key"
area = "infrastructure"
domains = ["infrastructure"]
version = 1
status = "supported"
description = "Removed schema key."
depends_on = []
{removed_key} = ["pubchem-1k"]
"#
            ),
        );
        let error = read_feature(&root.join(removed_key).join("feature.toml"))
            .expect_err("removed schema key should be rejected");
        assert!(error.to_string().contains("unknown field"));
    }
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
fn render_dashboard_is_stable_and_uses_compact_benchmark_cells() {
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
        },
        Feature {
            id: "failing.feature".to_owned(),
            title: "Failing".to_owned(),
            area: "benchmark".to_owned(),
            domains: vec![FeatureDomain::SmallMolecule],
            version: 1,
            status: FeatureStatus::Deprecated,
            description: "Feature with counted failures.".to_owned(),
            depends_on: Vec::new(),
        },
        Feature {
            id: "missing.feature".to_owned(),
            title: "Missing".to_owned(),
            area: "benchmark".to_owned(),
            domains: vec![FeatureDomain::Macromolecule],
            version: 1,
            status: FeatureStatus::Experimental,
            description: "Feature without recorded status.".to_owned(),
            depends_on: Vec::new(),
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
        },
    ];
    let results = BTreeMap::from([(
        "failing.feature".to_owned(),
        BenchmarkResults {
            feature_id: "failing.feature".to_owned(),
            corpora: BTreeMap::from([(
                "pubchem-1k".to_owned(),
                BenchmarkResult {
                    outcome: BenchmarkResultOutcome::Differences,
                    scope: "full".to_owned(),
                    fixture_count: 7,
                    compared_count: 4,
                    difference_count: 3,
                    first_detail: Some("fixture `data/bad.sdf` differs".to_owned()),
                    reference_tool: Some("rdkit".to_owned()),
                    reference_version: Some("RDKit test".to_owned()),
                    manifest_digest: Some("0".repeat(64)),
                    input_digest_schema_version: Some(BENCHMARK_INPUT_DIGEST_SCHEMA_VERSION),
                    input_digest: Some("1".repeat(64)),
                    input_count: 1,
                    legacy_source: None,
                    benchmarked_at_unix: 1,
                },
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

    let dashboard = render_dashboard(&features, &results, &corpus_info);

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
    assert!(dashboard.contains("aria-label=\"last differences\""));
    assert!(dashboard.contains("<span class=\"count\">3</span>"));
    assert!(dashboard.contains("<span class=\"available\">A</span>available"));
    assert!(dashboard.contains(
        "<span class=\"available\" aria-label=\"available\" title=\"benchmark available; no recorded result; reference: Biopython v1.87"
    ));
    assert!(dashboard.contains("never affect feature status or release health"));
    assert!(dashboard.contains("document.querySelectorAll('table.feature-dashboard')"));
    assert!(dashboard.contains("button.addEventListener('click'"));
    assert!(dashboard.ends_with('\n'));
}

#[test]
fn dashboard_corpus_cells_show_optional_benchmark_observations() {
    let feature = Feature {
        id: "optional.feature".to_owned(),
        title: "Optional".to_owned(),
        area: "benchmark".to_owned(),
        domains: vec![FeatureDomain::SmallMolecule],
        version: 1,
        status: FeatureStatus::Supported,
        description: "Feature with optional corpus evidence.".to_owned(),
        depends_on: Vec::new(),
    };
    let matched = BenchmarkResult {
        outcome: BenchmarkResultOutcome::Match,
        scope: "full".to_owned(),
        fixture_count: 1,
        compared_count: 1,
        difference_count: 0,
        first_detail: None,
        reference_tool: Some("rdkit".to_owned()),
        reference_version: Some("RDKit test".to_owned()),
        manifest_digest: Some("a".repeat(64)),
        input_digest_schema_version: Some(BENCHMARK_INPUT_DIGEST_SCHEMA_VERSION),
        input_digest: Some("b".repeat(64)),
        input_count: 1,
        legacy_source: None,
        benchmarked_at_unix: 1,
    };
    let results = BenchmarkResults {
        feature_id: feature.id.clone(),
        corpora: BTreeMap::from([("pubchem-1k".to_owned(), matched)]),
    };
    let reference = CorpusFeatureDashboardInfo {
        reference_tool: "rdkit".to_owned(),
        reference_version: "RDKit test".to_owned(),
    };

    assert!(dashboard_corpus_cell(
        &feature,
        Some(&results),
        "pubchem-1k",
        Some(&reference),
        true,
    )
    .contains("aria-label=\"last match\""));
    assert!(
        dashboard_corpus_cell(&feature, None, "pubchem-1k", Some(&reference), true)
            .contains("title=\"benchmark available; no recorded result; reference: RDKit vtest\"")
    );
    assert!(
        dashboard_corpus_cell(&feature, None, "pubchem-1k", None, true)
            .contains("aria-label=\"not available\"")
    );
    assert!(
        dashboard_corpus_cell(&feature, Some(&results), "pubchem-1k", None, true)
            .contains("aria-label=\"not available\"")
    );
    assert!(dashboard_corpus_cell(
        &feature,
        Some(&results),
        "pubchem-1k",
        Some(&reference),
        false
    )
    .contains("aria-label=\"not available\""));
}

#[test]
fn benchmark_manifest_path_is_feature_scoped() {
    assert_eq!(
        benchmark_manifest_path("core.graph", "smoke"),
        PathBuf::from("benchmarks/corpora/smoke/features/core.graph.toml")
    );
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
        "--feature".to_owned(),
        "all".to_owned(),
        "--jobs".to_owned()
    ])
    .is_err());
    assert!(benchmark_args(&[
        "--feature".to_owned(),
        "io.smiles.canonical".to_owned(),
        "--corpus".to_owned(),
        "pubchem-100k".to_owned(),
        "--fixture".to_owned(),
        "data/packs/pack_001.smi".to_owned(),
    ])
    .is_ok());
    assert!(benchmark_args(&[
        "--feature".to_owned(),
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
    let root = temp_feature_root("accept-implementation-goldens");
    let corpus_root = root.join("benchmarks/corpora/smoke");
    let manifest_path = corpus_root.join("features/stereo.perception.toml");
    fs::create_dir_all(manifest_path.parent().expect("manifest parent"))
        .expect("features directory");
    fs::create_dir_all(corpus_root.join("data")).expect("data directory");
    fs::write(corpus_root.join("data/example.smi"), "CC CID:1\n").expect("fixture should write");
    let mut manifest = BenchmarkManifest {
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
    assert_eq!(benchmark_worker_count(16, 3), 3);
    assert_eq!(benchmark_worker_count(0, 3), 1);
}

#[test]
fn benchmark_defaults_to_small_corpora_and_all_includes_broad_manifest_backed_features() {
    assert_eq!(benchmark_corpus_selector(&[]), "baseline");
    assert_eq!(
        benchmark_corpus_selector(&[
            "--feature".to_owned(),
            "all".to_owned(),
            "--corpus".to_owned(),
            "pubchem-1k".to_owned(),
        ]),
        "pubchem-1k"
    );

    let root = temp_feature_root("all-benchmark-corpora");
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
        let path = benchmark_manifest_path_from(&root, feature, corpus);
        fs::create_dir_all(path.parent().expect("manifest parent"))
            .expect("manifest directory should create");
        fs::write(path, "").expect("manifest marker should write");
    }

    assert_eq!(
        benchmark_targets_from(&root, &features, "all", "baseline")
            .into_iter()
            .map(|(feature, corpus)| (feature.id.as_str(), corpus))
            .collect::<Vec<_>>(),
        vec![
            ("small", "pubchem-1k".to_owned()),
            ("macro", "pdb-100".to_owned()),
        ]
    );
    assert_eq!(
        benchmark_targets_from(&root, &features, "all", "pubchem-1k")
            .into_iter()
            .map(|(feature, corpus)| (feature.id.as_str(), corpus))
            .collect::<Vec<_>>(),
        vec![("small", "pubchem-1k".to_owned())]
    );
    assert_eq!(
        benchmark_targets_from(&root, &features, "all", "all")
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
        benchmark_targets_from(&root, &features, "small", "all")
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
        benchmark_targets_from(&root, &features, "small", "pubchem-100k")
            .into_iter()
            .map(|(feature, corpus)| (feature.id.as_str(), corpus))
            .collect::<Vec<_>>(),
        vec![("small", "pubchem-100k".to_owned())]
    );
    assert_eq!(
        benchmark_targets_from(&root, &features, "macro", "pubchem-1k")
            .into_iter()
            .map(|(feature, corpus)| (feature.id.as_str(), corpus))
            .collect::<Vec<_>>(),
        Vec::<(&str, String)>::new()
    );
    fs::remove_dir_all(root).ok();
}

#[test]
fn concrete_missing_manifest_error_is_clear_while_all_selectors_skip() {
    assert_eq!(
        concrete_missing_manifest_error("io.smiles.parse", "pdb-100").as_deref(),
        Some("no benchmark manifest for feature `io.smiles.parse` and corpus `pdb-100`")
    );
    assert_eq!(concrete_missing_manifest_error("all", "pdb-100"), None);
    assert_eq!(
        concrete_missing_manifest_error("io.smiles.parse", "all"),
        None
    );
    assert_eq!(
        concrete_missing_manifest_error("io.smiles.parse", "baseline"),
        None
    );
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
fn implementation_dispatch_supports_hydrogen_transforms() {
    let root = temp_feature_root("hydrogen-transforms-dispatch");
    let fixture = root.join("fixture.sdf");
    fs::write(&fixture, simple_sdf_record("methane")).expect("fixture should write");

    let expected = implementation_expected("chem.hydrogen-transforms", "pubchem-1k", &fixture)
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
fn implementation_dispatch_supports_query_benchmark() {
    let root = temp_feature_root("query-benchmark-dispatch");
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
    let root = temp_feature_root("smiles-benchmark-subsets");
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
fn stereo_cip_benchmark_uses_rdkit_default_hydrogen_indexing() {
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
fn stereo_cip_benchmark_reads_all_sdf_pack_records() {
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
fn benchmark_digest_changes_after_material_input_changes() {
    let root = temp_feature_root("digest-change");
    let (_, _, manifest_path) = write_digest_test_repo(&root);
    let manifest = read_benchmark_manifest(&manifest_path).expect("manifest should read");
    let original = build_benchmark_input_digest(&root, &manifest_path, &manifest)
        .expect("digest should build");

    fs::write(root.join("crates/kekule/src/lib.rs"), "changed source\n")
        .expect("source should mutate");
    let source_changed = build_benchmark_input_digest(&root, &manifest_path, &manifest)
        .expect("digest should build");
    assert_ne!(original.sha256, source_changed.sha256);

    fs::write(
        root.join("benchmarks/corpora/smoke/data/example.sdf"),
        "changed fixture\n",
    )
    .expect("fixture should mutate");
    let fixture_changed = build_benchmark_input_digest(&root, &manifest_path, &manifest)
        .expect("digest should build");
    assert_ne!(source_changed.sha256, fixture_changed.sha256);

    fs::write(
        root.join("benchmarks/corpora/smoke/golden/example/data_example.sdf.json.gz"),
        "changed golden\n",
    )
    .expect("golden should mutate");
    let golden_changed = build_benchmark_input_digest(&root, &manifest_path, &manifest)
        .expect("digest should build");
    assert_ne!(fixture_changed.sha256, golden_changed.sha256);

    fs::write(
        root.join("benchmarks/reference/rdkit/run_feature.py"),
        "changed generator\n",
    )
    .expect("generator should mutate");
    let generator_changed = build_benchmark_input_digest(&root, &manifest_path, &manifest)
        .expect("digest should build");
    assert_ne!(golden_changed.sha256, generator_changed.sha256);

    fs::write(
        &manifest_path,
        "feature_id = \"example\"\ncorpus_id = \"smoke\"\nreference_tool = \"rdkit\"\nreference_version = \"RDKit changed\"\ncomparison_mode = \"implementation-golden\"\nfixtures = [\"data/example.sdf\"]\n",
    )
    .expect("manifest should mutate");
    let changed_manifest =
        read_benchmark_manifest(&manifest_path).expect("changed manifest should read");
    let manifest_changed = build_benchmark_input_digest(&root, &manifest_path, &changed_manifest)
        .expect("digest should build");
    assert_ne!(generator_changed.sha256, manifest_changed.sha256);

    write_digest_lock(
        &root.join("Cargo.lock"),
        "registry+https://example.invalid/index",
    );
    let dependency_changed = build_benchmark_input_digest(&root, &manifest_path, &changed_manifest)
        .expect("digest should build");
    assert_ne!(manifest_changed.sha256, dependency_changed.sha256);
    fs::remove_dir_all(root).ok();
}

#[test]
fn manual_semantic_reference_digest_does_not_require_generator_files() {
    let root = temp_feature_root("manual-reference-digest");
    let (_, _, manifest_path) = write_digest_test_repo(&root);
    fs::write(
        &manifest_path,
        "feature_id = \"example\"\ncorpus_id = \"smoke\"\nreference_tool = \"enamine-manual-semantic\"\nreference_version = \"Enamine Discovery Diversity Set 2026-07-05\"\ncomparison_mode = \"implementation-golden\"\nfixtures = [\"data/example.sdf\"]\n",
    )
    .expect("manual manifest should write");
    fs::remove_dir_all(root.join("benchmarks/reference")).ok();

    let manifest = read_benchmark_manifest(&manifest_path).expect("manifest should read");
    let digest = build_benchmark_input_digest(&root, &manifest_path, &manifest)
        .expect("digest should build");
    assert!(digest.input_count > 0);
    fs::remove_dir_all(root).ok();
}

#[test]
fn benchmark_hash_normalizes_text_line_endings() {
    let root = temp_feature_root("digest-line-endings");
    let path = root.join("source.rs");
    fs::write(&path, "first\nsecond\n").expect("LF source should write");
    let lf_hash = hash_normalized_file(&path).expect("LF input should hash");

    fs::write(&path, "first\r\nsecond\r\n").expect("CRLF source should write");
    let crlf_hash = hash_normalized_file(&path).expect("CRLF input should hash");

    assert_eq!(lf_hash, crlf_hash);
    fs::remove_dir_all(root).ok();
}

#[test]
fn benchmark_digest_ignores_checkout_and_workspace_identity() {
    let first = temp_feature_root("digest-identity-first");
    let second = temp_feature_root("digest-identity-second");
    let (_, _, first_manifest_path) = write_digest_test_repo(&first);
    let (_, _, second_manifest_path) = write_digest_test_repo(&second);
    let first_manifest =
        read_benchmark_manifest(&first_manifest_path).expect("manifest should read");
    let second_manifest =
        read_benchmark_manifest(&second_manifest_path).expect("manifest should read");

    let first_digest = build_benchmark_input_digest(&first, &first_manifest_path, &first_manifest)
        .expect("first digest should build");
    let second_digest =
        build_benchmark_input_digest(&second, &second_manifest_path, &second_manifest)
            .expect("second digest should build");
    assert_eq!(first_digest.sha256, second_digest.sha256);

    fs::write(
        second.join("Cargo.toml"),
        "[workspace]\n[workspace.package]\nrepository = \"https://example.invalid/renamed-repository\"\n",
    )
    .expect("workspace identity should mutate");
    write_digest_lock_with_local_name(&second.join("Cargo.lock"), "renamed-workspace-package");
    fs::write(
        second.join("features/example/feature.md"),
        "# Renamed repository documentation\n",
    )
    .expect("feature documentation should mutate");
    let renamed_digest =
        build_benchmark_input_digest(&second, &second_manifest_path, &second_manifest)
            .expect("renamed digest should build");
    assert_eq!(first_digest.sha256, renamed_digest.sha256);

    fs::remove_dir_all(first).ok();
    fs::remove_dir_all(second).ok();
}

#[test]
fn benchmark_digest_ignores_core_package_directory_identity() {
    let first = temp_feature_root("digest-package-identity-first");
    let second = temp_feature_root("digest-package-identity-second");
    let (_, _, first_manifest_path) = write_digest_test_repo(&first);
    let (_, _, second_manifest_path) = write_digest_test_repo(&second);
    let first_manifest =
        read_benchmark_manifest(&first_manifest_path).expect("manifest should read");
    let second_manifest =
        read_benchmark_manifest(&second_manifest_path).expect("manifest should read");

    let first_core_source = first.join("crates/kekule/src");
    let renamed_package = second.join("crates/renamed-core");
    fs::rename(second.join("crates/kekule"), &renamed_package)
        .expect("core package directory should rename");
    let second_core_source = renamed_package.join("src");

    let first_digest = build_benchmark_input_digest_with_core_source_root(
        &first,
        &first_manifest_path,
        &first_manifest,
        &first_core_source,
    )
    .expect("first digest should build");
    let second_digest = build_benchmark_input_digest_with_core_source_root(
        &second,
        &second_manifest_path,
        &second_manifest,
        &second_core_source,
    )
    .expect("renamed-package digest should build");

    assert_eq!(first_digest.sha256, second_digest.sha256);
    fs::remove_dir_all(first).ok();
    fs::remove_dir_all(second).ok();
}

#[test]
fn dashboard_text_comparison_ignores_platform_line_endings() {
    assert_eq!(
        normalize_text_line_endings("one\r\ntwo\rthree\n"),
        "one\ntwo\nthree\n"
    );
}

#[test]
fn dashboard_result_requires_a_current_manifest() {
    let feature = Feature {
        id: "portable.feature".to_owned(),
        title: "Portable".to_owned(),
        area: "infrastructure".to_owned(),
        domains: vec![FeatureDomain::Infrastructure],
        version: 1,
        status: FeatureStatus::Supported,
        description: "Portable dashboard evidence.".to_owned(),
        depends_on: Vec::new(),
    };
    let results = BenchmarkResults {
        feature_id: feature.id.clone(),
        corpora: BTreeMap::from([(
            "pubchem-1k".to_owned(),
            BenchmarkResult {
                outcome: BenchmarkResultOutcome::Match,
                scope: "full".to_owned(),
                fixture_count: 1,
                compared_count: 1,
                difference_count: 0,
                first_detail: None,
                reference_tool: Some("rdkit".to_owned()),
                reference_version: Some("test".to_owned()),
                manifest_digest: Some("0".repeat(64)),
                input_digest_schema_version: Some(BENCHMARK_INPUT_DIGEST_SCHEMA_VERSION),
                input_digest: Some("1".repeat(64)),
                input_count: 1,
                legacy_source: None,
                benchmarked_at_unix: 1,
            },
        )]),
    };
    let reference = CorpusFeatureDashboardInfo {
        reference_tool: "rdkit".to_owned(),
        reference_version: "test".to_owned(),
    };

    assert!(
        dashboard_corpus_cell(&feature, Some(&results), "pubchem-1k", None, true)
            .contains("no benchmark manifest")
    );
    assert!(dashboard_corpus_cell(
        &feature,
        Some(&results),
        "pubchem-1k",
        Some(&reference),
        true,
    )
    .contains("aria-label=\"last match\""));
}
#[test]
fn result_writer_prunes_entries_without_manifests_and_records_errors() {
    let root = temp_feature_root("result-manifest-pruning");
    let results_path = benchmark_results_path_from(&root, "pubchem-1k");
    fs::create_dir_all(results_path.parent().expect("results parent"))
        .expect("results directory should create");
    fs::write(&results_path, "stale").expect("stale results should write");

    let feature_result = BenchmarkResult::from_run(BenchmarkRun {
        outcome: BenchmarkResultOutcome::Error,
        scope: "fixture:data/example.sdf".to_owned(),
        fixture_count: 1,
        compared_count: 0,
        difference_count: 0,
        first_detail: Some("fixture could not be read".to_owned()),
        reference_tool: Some("rdkit".to_owned()),
        reference_version: Some("test".to_owned()),
        manifest_digest: Some("0".repeat(64)),
        input_digest: None,
    })
    .expect("error result should build");
    let results = BTreeMap::from([(
        "example.feature".to_owned(),
        BenchmarkResults {
            feature_id: "example.feature".to_owned(),
            corpora: BTreeMap::from([("pubchem-1k".to_owned(), feature_result)]),
        },
    )]);
    let selected = BTreeSet::from(["pubchem-1k".to_owned()]);

    write_benchmark_results_from(&root, &results, &selected)
        .expect("result pruning should succeed");
    assert!(!results_path.exists());

    let manifest = benchmark_manifest_path_from(&root, "example.feature", "pubchem-1k");
    fs::create_dir_all(manifest.parent().expect("manifest parent"))
        .expect("manifest directory should create");
    fs::write(&manifest, "").expect("manifest marker should write");
    write_benchmark_results_from(&root, &results, &selected)
        .expect("manifest-backed result should write");
    let stored = read_corpus_results(&results_path).expect("written results should parse");
    let stored_result = stored
        .features
        .get("example.feature")
        .expect("feature result should exist");
    assert_eq!(stored_result.outcome, BenchmarkResultOutcome::Error);
    assert_eq!(stored_result.scope, "fixture:data/example.sdf");
    assert_eq!(
        stored_result.first_detail.as_deref(),
        Some("fixture could not be read")
    );

    fs::remove_dir_all(root).ok();
}
#[test]
fn historical_legacy_result_is_replaced_by_a_new_target_snapshot() {
    let root = temp_feature_root("legacy-result");
    let path = benchmark_results_path_from(&root, "pubchem-1k");
    fs::create_dir_all(path.parent().expect("result parent")).expect("result directory");
    fs::write(
        &path,
        r#"
corpus_id = "pubchem-1k"

[features.example]
outcome = "match"
scope = "full"
fixture_count = 1
compared_count = 1
manifest_digest = "0000000000000000000000000000000000000000000000000000000000000000"
input_digest_schema_version = 2
input_digest = "1111111111111111111111111111111111111111111111111111111111111111"
input_count = 1
legacy_source = "validation-status-v2"
benchmarked_at_unix = 1
"#,
    )
    .expect("legacy result fixture should write");
    let legacy = read_corpus_results(&path).expect("legacy result should deserialize");
    assert_eq!(
        legacy.features["example"].legacy_source.as_deref(),
        Some("validation-status-v2")
    );

    let manifest = benchmark_manifest_path_from(&root, "example", "pubchem-1k");
    fs::create_dir_all(manifest.parent().expect("manifest parent")).expect("manifest directory");
    fs::write(manifest, "").expect("manifest marker");
    let current = BenchmarkResult::from_run(BenchmarkRun {
        outcome: BenchmarkResultOutcome::Differences,
        scope: "fixture:data/example.sdf".to_owned(),
        fixture_count: 1,
        compared_count: 1,
        difference_count: 1,
        first_detail: Some("changed".to_owned()),
        reference_tool: None,
        reference_version: None,
        manifest_digest: None,
        input_digest: None,
    })
    .expect("current result should build");
    let results = BTreeMap::from([(
        "example".to_owned(),
        BenchmarkResults {
            feature_id: "example".to_owned(),
            corpora: BTreeMap::from([("pubchem-1k".to_owned(), current)]),
        },
    )]);
    write_benchmark_results_from(&root, &results, &BTreeSet::from(["pubchem-1k".to_owned()]))
        .expect("current result should replace legacy result");
    let replaced = read_corpus_results(&path).expect("replacement should parse");
    let result = replaced
        .features
        .get("example")
        .expect("feature result should exist");
    assert_eq!(result.outcome, BenchmarkResultOutcome::Differences);
    assert_eq!(result.legacy_source, None);
    fs::remove_dir_all(root).ok();
}

#[test]
fn unsupported_comparison_mode_is_rejected() {
    let manifest = BenchmarkManifest {
        feature_id: "example".to_owned(),
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
fn benchmark_comparison_counts_multiple_fixture_failures() {
    let root = temp_feature_root("comparison-counts-failures");
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
    let manifest = read_benchmark_manifest(&manifest_path).expect("manifest should read");

    let comparison = compare_golden_outputs(&manifest_path, &manifest, 1, None)
        .expect("comparison should complete");

    assert_eq!(comparison.compared_count, 0);
    assert_eq!(comparison.difference_count, 2);
    assert!(comparison
        .first_difference
        .as_deref()
        .is_some_and(|failure| failure.contains("data/one.smi")));
    fs::remove_dir_all(root).ok();
}

#[test]
fn stereo_perception_benchmark_records_reference_preparation_errors_per_record() {
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
        Some("normalization_or_perception_error")
    );
    assert!(value.get("report").is_none());
}

#[test]
fn smiles_component_benchmarks_preserve_source_record_cardinality() {
    let root = temp_feature_root("smiles-component-benchmark-cardinality");
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
                .into_molecule()
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
        .into_molecule()
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
    assert!(stereo["records"][0].get("normalization_report").is_none());

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
        smiles_perceived_bonds_json(single.graph()),
        smiles_perceived_bonds_json(double.graph())
    );

    let aromatic = SmallMolecule::from_smiles("c1ccccc1").expect("benzene should parse");
    let mut perceived_aromatic = aromatic.clone();
    normalize_and_perceive(&mut perceived_aromatic);
    assert_eq!(
        explicit_valence_json(perceived_aromatic.graph(), AtomId::new(0)),
        3
    );
    let mut aromatic_cyclohexyne =
        SmallMolecule::from_smiles("C1=CC#CC=C1").expect("cyclohexyne parses");
    normalize_and_perceive(&mut aromatic_cyclohexyne);
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
    normalize_and_perceive(&mut thiophene);
    let sulfur_id = thiophene
        .graph()
        .atoms()
        .find_map(|(id, atom)| (atom.element.symbol() == "S").then_some(id))
        .expect("sulfur atom");
    assert_eq!(explicit_valence_json(thiophene.graph(), sulfur_id), 2);
    let mut phosphorus_ring =
        SmallMolecule::from_smiles("C(F)(F)(F)P1P(P(P(P1C(F)(F)F)C(F)(F)F)C(F)(F)F)C(F)(F)F")
            .expect("phosphorus ring parses");
    normalize_and_perceive(&mut phosphorus_ring);
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
    normalize_and_perceive(&mut phosphinine);
    let phosphinine_phosphorus = phosphinine
        .graph()
        .atoms()
        .find_map(|(id, atom)| (atom.element.symbol() == "P").then_some(id))
        .expect("phosphinine phosphorus");
    assert_eq!(
        explicit_valence_json(phosphinine.graph(), phosphinine_phosphorus),
        3
    );
    let document = kekule::smiles::parse_str("CN(C)CCO.C1=CC=C2C(=C1)C3=NC4=C5C=CC=CC5=C([N-]4)N=C6C7=CC=CC=C7C(=N6)N=C8C9=CC=CC=C9C(=N8)N=C2[N-]3.[Cu+2]")
        .expect("anionic macrocycle mixture parses");
    let mut anionic_macrocycle = kekule::smiles::interpret(&document)
        .expect("anionic macrocycle mixture interprets")
        .into_molecules()
        .swap_remove(1);
    normalize_and_perceive(&mut anionic_macrocycle);
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
    normalize_and_perceive(&mut cyclopentadienyl);
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
    normalize_and_perceive(&mut substituted_cyclopentadienyl);
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
    normalize_and_perceive(&mut fused_triazine);
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
    assert!(smiles_perceived_bonds_json(aromatic.graph())
        .iter()
        .all(|bond| bond["bond_type"] == "AROMATIC" && bond["is_aromatic"] == false));
    assert!(perceived_aromatic
        .graph()
        .bonds()
        .all(|(_, bond)| bond.order != BondOrder::Aromatic));
    assert!(smiles_perceived_bonds_json(perceived_aromatic.graph())
        .iter()
        .all(|bond| bond["bond_type"] == "AROMATIC" && bond["is_aromatic"] == true));

    let labeled = SmallMolecule::from_smiles("[13CH3:7]C").expect("labeled carbon should parse");
    let atoms = smiles_perceived_atoms_json(labeled.graph());
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
fn canonical_smiles_benchmark_normalizes_and_perceives_before_writing() {
    let root = temp_feature_root("canonical-perceive-before-write");
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
    let root = temp_feature_root("canonical-invalid-input");
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
    let molecule =
        SmallMolecule::from_smiles("[nH]1cccc1").expect("aromatic nH SMILES should parse");

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
    let molecule = SmallMolecule::from_smiles("CCOC(=O)C1=C(C(=C(N1)C)C(=O)OC(C)(C)C)C")
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

fn temp_feature_root(label: &str) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("time should be available")
        .as_nanos();
    let root = env::temp_dir().join(format!("kekule-xtask-{label}-{}-{nonce}", process::id()));
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
    }
}

fn write_digest_test_repo(root: &Path) -> (PathBuf, PathBuf, PathBuf) {
    let features_root = root.join("features");
    let benchmark_root = root.join("benchmarks");
    let feature_dir = features_root.join("example");
    let corpus_root = benchmark_root.join("corpora").join("smoke");
    let manifest_dir = corpus_root.join("features");
    fs::create_dir_all(&feature_dir).expect("feature dir should create");
    fs::create_dir_all(&manifest_dir).expect("manifest dir should create");
    fs::write(
        feature_dir.join("feature.toml"),
        "id = \"example\"\ntitle = \"Example\"\narea = \"test\"\ndomains = [\"small-molecule\"]\nversion = 1\nstatus = \"supported\"\ndescription = \"Example feature.\"\ndepends_on = []\n",
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
    write_digest_lock(
        &root.join("Cargo.lock"),
        "registry+https://github.com/rust-lang/crates.io-index",
    );
    for path in [
        "crates/kekule/Cargo.toml",
        "crates/xtask/Cargo.toml",
        "crates/kekule/src/lib.rs",
        "crates/xtask/src/main.rs",
        "benchmarks/reference/rdkit/run_feature.py",
        "benchmarks/reference/rdkit/environment.yml",
    ] {
        let path = root.join(path);
        fs::create_dir_all(path.parent().expect("test path should have parent"))
            .expect("test parent should create");
        fs::write(path, "test\n").expect("test evidence file should write");
    }
    (features_root, benchmark_root, manifest_path)
}

fn write_digest_lock(path: &Path, external_source: &str) {
    write_digest_lock_contents(path, "workspace-package", external_source);
}

fn write_digest_lock_with_local_name(path: &Path, local_name: &str) {
    write_digest_lock_contents(
        path,
        local_name,
        "registry+https://github.com/rust-lang/crates.io-index",
    );
}

fn write_digest_lock_contents(path: &Path, local_name: &str, external_source: &str) {
    fs::write(
        path,
        format!(
            r#"version = 3

[[package]]
name = "{local_name}"
version = "0.1.0"

[[package]]
name = "serde"
version = "1.0.0"
source = "{external_source}"
checksum = "0000000000000000000000000000000000000000000000000000000000000000"
dependencies = []
"#
        ),
    )
    .expect("digest lock should write");
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

fn write_gzip_json(path: &Path, value: &Value) {
    let file = fs::File::create(path).expect("gzip file should create");
    let mut encoder = flate2::write::GzEncoder::new(file, flate2::Compression::default());
    encoder
        .write_all(value.to_string().as_bytes())
        .expect("gzip json should write");
    encoder.finish().expect("gzip json should finish");
}
