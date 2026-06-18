use std::collections::BTreeMap;
use std::env;
use std::error::Error;
use std::fs;
use std::path::{Path, PathBuf};
use std::process;

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let mut args = env::args().skip(1);
    match args.next().as_deref() {
        Some("dashboard") => dashboard(args.collect()),
        Some("validate") => validate(args.collect()),
        Some("features") => list_features(),
        Some("skills") => skills(args.collect()),
        _ => {
            print_help();
            Ok(())
        }
    }
}

fn print_help() {
    eprintln!(
        "usage:\n  cargo xtask dashboard [--check]\n  cargo xtask validate --feature FEATURE_ID\n  cargo xtask skills --check\n  cargo xtask features"
    );
}

fn dashboard(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    let check = args.iter().any(|arg| arg == "--check");
    let features = read_features()?;
    let rendered = render_dashboard(&features);
    let path = Path::new("features/DASHBOARD.md");

    if check {
        let existing = fs::read_to_string(path)?;
        if existing != rendered {
            return Err(boxed_error(
                "features/DASHBOARD.md is out of date; run `cargo xtask dashboard`",
            ));
        }
    } else {
        fs::write(path, rendered)?;
    }
    Ok(())
}

fn validate(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    let feature = value_after_flag(&args, "--feature")
        .ok_or_else(|| boxed_error("missing required flag: --feature FEATURE_ID"))?;
    let features = read_features()?;
    if !features.iter().any(|candidate| candidate.id == feature) {
        return Err(boxed_error(format!("unknown feature: {feature}")));
    }

    println!("validation harness found feature `{feature}`");
    let manifest_path = validation_manifest_path(feature);
    if manifest_path.exists() {
        let manifest = read_validation_manifest(&manifest_path)?;
        if manifest.feature_id != feature {
            return Err(boxed_error(format!(
                "{} declares feature_id `{}`, expected `{feature}`",
                manifest_path.display(),
                manifest.feature_id
            )));
        }
        println!(
            "validation manifest uses {} {}",
            manifest.reference_tool, manifest.reference_version
        );
        validate_manifest_paths(&manifest_path, &manifest)?;
        println!(
            "validation manifest lists {} fixture(s)",
            manifest.fixtures.len()
        );
    } else {
        println!("no reference validation manifest configured for `{feature}`");
    }
    Ok(())
}

fn list_features() -> Result<(), Box<dyn Error>> {
    for feature in read_features()? {
        println!(
            "{}\t{}\tv{}\timplemented={}\tvalidated={}",
            feature.id, feature.area, feature.version, feature.implemented, feature.validated
        );
    }
    Ok(())
}

fn skills(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    if args.iter().any(|arg| arg != "--check") {
        return Err(boxed_error("usage: cargo xtask skills --check"));
    }
    check_skills(Path::new(".codex/skills"))?;
    println!("repo-local feature skills are in sync");
    Ok(())
}

fn value_after_flag<'a>(args: &'a [String], flag: &str) -> Option<&'a str> {
    args.windows(2)
        .find(|window| window[0] == flag)
        .map(|window| window[1].as_str())
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Feature {
    id: String,
    title: String,
    area: String,
    version: u32,
    implemented: bool,
    validated: bool,
    description: String,
    depends_on: Vec<String>,
}

fn read_features() -> Result<Vec<Feature>, Box<dyn Error>> {
    read_features_from(Path::new("features"))
}

fn read_features_from(root: &Path) -> Result<Vec<Feature>, Box<dyn Error>> {
    let mut features = Vec::new();
    for entry in fs::read_dir(root)? {
        let entry = entry?;
        let path = entry.path();
        if !path.is_dir() || is_hidden_or_template(&path) {
            continue;
        }
        let metadata_path = path.join("feature.toml");
        if metadata_path.exists() {
            features.push(read_feature(&metadata_path)?);
        }
    }
    features.sort_by(|a, b| a.id.cmp(&b.id));
    validate_feature_set(&features)?;
    Ok(features)
}

fn is_hidden_or_template(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .map(|name| name.starts_with('_'))
        .unwrap_or(false)
}

fn read_feature(path: &Path) -> Result<Feature, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    let map = parse_simple_toml(&text);
    reject_deprecated_feature_keys(&map, path)?;
    let feature = Feature {
        id: required(&map, "id", path)?,
        title: required(&map, "title", path)?,
        area: required(&map, "area", path)?,
        version: required_u32(&map, "version", path)?,
        implemented: required_bool(&map, "implemented", path)?,
        validated: required_bool(&map, "validated", path)?,
        description: required(&map, "description", path)?,
        depends_on: required_string_array(&map, "depends_on", path)?,
    };
    validate_feature(&feature, path)?;
    Ok(feature)
}

fn required(
    map: &BTreeMap<String, String>,
    key: &str,
    path: &Path,
) -> Result<String, Box<dyn Error>> {
    map.get(key)
        .cloned()
        .ok_or_else(|| boxed_error(format!("{} is missing `{key}`", path.display())))
}

fn reject_deprecated_feature_keys(
    map: &BTreeMap<String, String>,
    path: &Path,
) -> Result<(), Box<dyn Error>> {
    for key in ["priority", "status", "last_ai_review"] {
        if map.contains_key(key) {
            return Err(boxed_error(format!(
                "{} uses deprecated feature metadata key `{key}`",
                path.display()
            )));
        }
    }
    Ok(())
}

fn parse_simple_toml(text: &str) -> BTreeMap<String, String> {
    let mut map = BTreeMap::new();
    let mut pending: Option<(String, String)> = None;
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some((key, mut value)) = pending.take() {
            value.push(' ');
            value.push_str(line);
            if line.ends_with(']') {
                map.insert(key, normalize_value(&value));
            } else {
                pending = Some((key, value));
            }
            continue;
        }
        if let Some((key, value)) = line.split_once('=') {
            let key = key.trim().to_owned();
            let value = value.trim();
            if value.starts_with('[') && !value.ends_with(']') {
                pending = Some((key, value.to_owned()));
            } else {
                map.insert(key, normalize_value(value));
            }
        }
    }
    map
}

fn required_bool(
    map: &BTreeMap<String, String>,
    key: &str,
    path: &Path,
) -> Result<bool, Box<dyn Error>> {
    match required(map, key, path)?.as_str() {
        "true" => Ok(true),
        "false" => Ok(false),
        value => Err(boxed_error(format!(
            "{} has invalid boolean `{key}` value `{value}`",
            path.display()
        ))),
    }
}

fn required_u32(
    map: &BTreeMap<String, String>,
    key: &str,
    path: &Path,
) -> Result<u32, Box<dyn Error>> {
    let value = required(map, key, path)?;
    let parsed = value.parse::<u32>().map_err(|_| {
        boxed_error(format!(
            "{} has invalid integer `{key}` value `{value}`",
            path.display()
        ))
    })?;
    if parsed == 0 {
        return Err(boxed_error(format!(
            "{} has invalid zero `{key}` value",
            path.display()
        )));
    }
    Ok(parsed)
}

fn required_string_array(
    map: &BTreeMap<String, String>,
    key: &str,
    path: &Path,
) -> Result<Vec<String>, Box<dyn Error>> {
    let value = required(map, key, path)?;
    parse_string_array(&value).ok_or_else(|| {
        boxed_error(format!(
            "{} has invalid string array `{key}` value `{value}`",
            path.display()
        ))
    })
}

fn parse_string_array(value: &str) -> Option<Vec<String>> {
    let value = value.trim();
    let inner = value.strip_prefix('[')?.strip_suffix(']')?.trim();
    if inner.is_empty() {
        return Some(Vec::new());
    }
    inner
        .split(',')
        .filter(|item| !item.trim().is_empty())
        .map(|item| {
            let item = item.trim();
            if item.starts_with('"') && item.ends_with('"') && item.len() >= 2 {
                Some(item[1..item.len() - 1].to_owned())
            } else {
                None
            }
        })
        .collect()
}

fn validate_feature(feature: &Feature, path: &Path) -> Result<(), Box<dyn Error>> {
    let expected_id = path
        .parent()
        .and_then(Path::file_name)
        .and_then(|name| name.to_str())
        .ok_or_else(|| {
            boxed_error(format!(
                "cannot determine feature directory for {}",
                path.display()
            ))
        })?;
    if feature.id != expected_id {
        return Err(boxed_error(format!(
            "{} declares id `{}`, expected `{expected_id}`",
            path.display(),
            feature.id
        )));
    }
    for (key, value) in [
        ("title", feature.title.as_str()),
        ("area", feature.area.as_str()),
        ("description", feature.description.as_str()),
    ] {
        if value.trim().is_empty() {
            return Err(boxed_error(format!(
                "{} has empty required field `{key}`",
                path.display()
            )));
        }
    }
    let feature_doc = path.with_file_name("feature.md");
    if !feature_doc.exists() {
        return Err(boxed_error(format!(
            "{} is missing required feature.md",
            path.parent().unwrap_or_else(|| Path::new(".")).display()
        )));
    }
    Ok(())
}

fn validate_feature_set(features: &[Feature]) -> Result<(), Box<dyn Error>> {
    let mut seen = BTreeMap::<&str, ()>::new();
    for feature in features {
        if seen.insert(feature.id.as_str(), ()).is_some() {
            return Err(boxed_error(format!(
                "duplicate feature id `{}`",
                feature.id
            )));
        }
    }
    for feature in features {
        for dependency in &feature.depends_on {
            if !seen.contains_key(dependency.as_str()) {
                return Err(boxed_error(format!(
                    "feature `{}` depends on unknown feature `{dependency}`",
                    feature.id
                )));
            }
        }
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ValidationManifest {
    feature_id: String,
    reference_tool: String,
    reference_version: String,
    fixtures: Vec<String>,
}

fn validation_manifest_path(feature: &str) -> PathBuf {
    Path::new("validation")
        .join("features")
        .join(feature)
        .join("validation.toml")
}

fn read_validation_manifest(path: &Path) -> Result<ValidationManifest, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    let map = parse_simple_toml(&text);
    Ok(ValidationManifest {
        feature_id: required(&map, "feature_id", path)?,
        reference_tool: required(&map, "reference_tool", path)?,
        reference_version: required(&map, "reference_version", path)?,
        fixtures: optional_string_array(&map, "fixtures", path)?,
    })
}

fn optional_string_array(
    map: &BTreeMap<String, String>,
    key: &str,
    path: &Path,
) -> Result<Vec<String>, Box<dyn Error>> {
    match map.get(key) {
        Some(value) => parse_string_array(value).ok_or_else(|| {
            boxed_error(format!(
                "{} has invalid string array `{key}` value `{value}`",
                path.display()
            ))
        }),
        None => Ok(Vec::new()),
    }
}

fn validate_manifest_paths(
    manifest_path: &Path,
    manifest: &ValidationManifest,
) -> Result<(), Box<dyn Error>> {
    let base = manifest_path.parent().ok_or_else(|| {
        boxed_error(format!(
            "{} has no parent directory",
            manifest_path.display()
        ))
    })?;
    for fixture in &manifest.fixtures {
        let path = base.join(fixture);
        if !path.exists() {
            return Err(boxed_error(format!(
                "{} references missing fixture `{fixture}`",
                manifest_path.display()
            )));
        }
    }
    Ok(())
}

fn normalize_value(value: &str) -> String {
    let value = value.trim().trim_end_matches(',').trim();
    if value.starts_with('"') && value.ends_with('"') && value.len() >= 2 {
        value[1..value.len() - 1].to_owned()
    } else {
        value.to_owned()
    }
}

fn render_dashboard(features: &[Feature]) -> String {
    let mut out = String::new();
    out.push_str("# Feature Dashboard\n\n");
    out.push_str("Generated from `features/*/feature.toml`. Do not hand-edit this file.\n\n");
    out.push_str("| Feature | Title | Area | Version | Implemented | Validated |\n");
    out.push_str("|---|---|---|---:|:---:|:---:|\n");
    for feature in features {
        out.push_str(&format!(
            "| `{}` | {} | {} | {} | {} | {} |\n",
            feature.id,
            feature.title,
            feature.area,
            feature.version,
            checkmark(feature.implemented),
            checkmark(feature.validated)
        ));
    }
    out
}

fn checkmark(value: bool) -> &'static str {
    if value {
        "✅"
    } else {
        "❌"
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SkillMetadata {
    name: String,
    description: String,
}

fn check_skills(root: &Path) -> Result<(), Box<dyn Error>> {
    for expected in expected_skills() {
        let path = root.join(expected.name).join("SKILL.md");
        if !path.exists() {
            return Err(boxed_error(format!(
                "missing repo-local skill `{}` at {}",
                expected.name,
                path.display()
            )));
        }
        let text = fs::read_to_string(&path)?;
        let metadata = parse_skill_metadata(&text, &path)?;
        if metadata.name != expected.name {
            return Err(boxed_error(format!(
                "{} declares skill name `{}`, expected `{}`",
                path.display(),
                metadata.name,
                expected.name
            )));
        }
        let lower = text.to_lowercase();
        for required in expected.required_phrases {
            if !lower.contains(&required.to_lowercase()) {
                return Err(boxed_error(format!(
                    "{} is missing required phrase `{required}`",
                    path.display()
                )));
            }
        }
    }
    Ok(())
}

struct ExpectedSkill {
    name: &'static str,
    required_phrases: &'static [&'static str],
}

fn expected_skills() -> &'static [ExpectedSkill] {
    &[
        ExpectedSkill {
            name: "feature-work",
            required_phrases: &[
                "add -> optional research -> plan -> implement",
                "feature.md",
                "implemented = true",
                "validated = true",
                "cargo xtask dashboard --check",
                "cargo xtask validate --feature",
            ],
        },
        ExpectedSkill {
            name: "feature-review",
            required_phrases: &[
                "independent audit",
                "architecture",
                "validation claims",
                "feature.md",
                "cargo test --workspace",
                "cargo xtask validate --feature",
            ],
        },
    ]
}

fn parse_skill_metadata(text: &str, path: &Path) -> Result<SkillMetadata, Box<dyn Error>> {
    let mut lines = text.lines();
    if lines.next() != Some("---") {
        return Err(boxed_error(format!(
            "{} is missing YAML frontmatter",
            path.display()
        )));
    }
    let mut fields = BTreeMap::<String, String>::new();
    for line in lines.by_ref() {
        if line == "---" {
            let name = fields
                .get("name")
                .cloned()
                .ok_or_else(|| boxed_error(format!("{} is missing `name`", path.display())))?;
            let description = fields.get("description").cloned().ok_or_else(|| {
                boxed_error(format!("{} is missing `description`", path.display()))
            })?;
            if name.trim().is_empty() || description.trim().is_empty() {
                return Err(boxed_error(format!(
                    "{} has empty skill frontmatter",
                    path.display()
                )));
            }
            return Ok(SkillMetadata { name, description });
        }
        if let Some((key, value)) = line.split_once(':') {
            fields.insert(key.trim().to_owned(), value.trim().to_owned());
        }
    }
    Err(boxed_error(format!(
        "{} has unterminated YAML frontmatter",
        path.display()
    )))
}

fn boxed_error(message: impl Into<String>) -> Box<dyn Error> {
    std::io::Error::other(message.into()).into()
}

#[cfg(test)]
mod tests {
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
    fn simple_toml_parser_normalizes_strings_booleans_and_arrays() {
        let parsed = parse_simple_toml(
            r#"
            id = "core.graph"
            implemented = true
            depends_on = []
            "#,
        );

        assert_eq!(parsed.get("id"), Some(&"core.graph".to_owned()));
        assert_eq!(parsed.get("implemented"), Some(&"true".to_owned()));
        assert_eq!(parsed.get("depends_on"), Some(&"[]".to_owned()));
    }

    #[test]
    fn string_array_parser_accepts_empty_and_string_lists() {
        assert_eq!(parse_string_array("[]"), Some(Vec::new()));
        assert_eq!(
            parse_string_array(r#"["core.graph", "validation.harness"]"#),
            Some(vec![
                "core.graph".to_owned(),
                "validation.harness".to_owned()
            ])
        );
        assert_eq!(parse_string_array("[core.graph]"), None);
    }

    #[test]
    fn simple_toml_parser_accepts_multiline_string_arrays() {
        let parsed = parse_simple_toml(
            r#"
            feature_id = "io.sdf.v2000.parse"
            fixtures = [
              "fixtures/a.sdf",
              "fixtures/b.sdf",
            ]
            "#,
        );

        assert_eq!(
            parse_string_array(parsed.get("fixtures").expect("fixtures should parse")),
            Some(vec![
                "fixtures/a.sdf".to_owned(),
                "fixtures/b.sdf".to_owned()
            ])
        );
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
            },
        ];

        let dashboard = render_dashboard(&features);

        assert!(
            dashboard.contains("| Feature | Title | Area | Version | Implemented | Validated |")
        );
        assert!(dashboard.contains("| `a.feature` | Aye | core | 1 | ❌ | ❌ |"));
        assert!(dashboard.contains("| `z.feature` | Zed | io | 3 | ✅ | ❌ |"));
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
Use feature.md. Only set implemented = true and validated = true with evidence.
Run cargo xtask dashboard --check and cargo xtask validate --feature <feature-id>.
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
Read feature.md. Run cargo test --workspace and cargo xtask validate --feature <feature-id>.
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
            validation_manifest_path("core.graph"),
            PathBuf::from("validation/features/core.graph/validation.toml")
        );
    }

    #[test]
    fn validation_manifest_reads_and_checks_fixture_paths() {
        let root = temp_feature_root("validation-manifest");
        let feature_dir = root.join("validation").join("features").join("example");
        let fixture_dir = feature_dir.join("fixtures");
        fs::create_dir_all(&fixture_dir).expect("fixture dir should create");
        fs::write(fixture_dir.join("ok.txt"), "{}").expect("fixture should write");
        let manifest_path = feature_dir.join("validation.toml");
        fs::write(
            &manifest_path,
            r#"feature_id = "example"
reference_tool = "manual-fixtures"
reference_version = "test"
fixtures = [
  "fixtures/ok.txt",
]
"#,
        )
        .expect("manifest should write");

        let manifest = read_validation_manifest(&manifest_path).expect("manifest should parse");
        assert_eq!(manifest.fixtures, vec!["fixtures/ok.txt"]);
        validate_manifest_paths(&manifest_path, &manifest).expect("fixture should exist");

        fs::write(
            &manifest_path,
            r#"feature_id = "example"
reference_tool = "manual-fixtures"
reference_version = "test"
fixtures = [
  "fixtures/missing.txt",
]
"#,
        )
        .expect("manifest should rewrite");
        let manifest = read_validation_manifest(&manifest_path).expect("manifest should parse");
        assert!(validate_manifest_paths(&manifest_path, &manifest).is_err());
        fs::remove_dir_all(root).ok();
    }

    fn temp_feature_root(label: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time should be available")
            .as_nanos();
        let root =
            env::temp_dir().join(format!("molecules-xtask-{label}-{}-{nonce}", process::id()));
        fs::create_dir_all(&root).expect("temp feature root should create");
        root
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
}
