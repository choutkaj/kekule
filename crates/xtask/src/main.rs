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
        _ => {
            print_help();
            Ok(())
        }
    }
}

fn print_help() {
    eprintln!("usage:\n  cargo xtask dashboard [--check]\n  cargo xtask validate --feature FEATURE_ID\n  cargo xtask features");
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
    } else {
        println!("no reference validation manifest configured for `{feature}`");
    }
    Ok(())
}

fn list_features() -> Result<(), Box<dyn Error>> {
    for feature in read_features()? {
        println!("{}\t{}\t{}", feature.id, feature.priority, feature.status);
    }
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
    priority: String,
    status: String,
    implemented: bool,
    validated: bool,
    last_ai_review: String,
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
    let feature = Feature {
        id: required(&map, "id", path)?,
        title: required(&map, "title", path)?,
        area: required(&map, "area", path)?,
        priority: required(&map, "priority", path)?,
        status: required(&map, "status", path)?,
        implemented: required_bool(&map, "implemented", path)?,
        validated: required_bool(&map, "validated", path)?,
        last_ai_review: required(&map, "last_ai_review", path)?,
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

fn parse_simple_toml(text: &str) -> BTreeMap<String, String> {
    let mut map = BTreeMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some((key, value)) = line.split_once('=') {
            map.insert(key.trim().to_owned(), normalize_value(value));
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
        ("last_ai_review", feature.last_ai_review.as_str()),
    ] {
        if value.trim().is_empty() {
            return Err(boxed_error(format!(
                "{} has empty required field `{key}`",
                path.display()
            )));
        }
    }
    if !matches!(feature.priority.as_str(), "P0" | "P1" | "P2" | "P3") {
        return Err(boxed_error(format!(
            "{} has invalid priority `{}`",
            path.display(),
            feature.priority
        )));
    }
    if !matches!(
        feature.status.as_str(),
        "planned" | "implemented" | "validated" | "deferred" | "blocked"
    ) {
        return Err(boxed_error(format!(
            "{} has invalid status `{}`",
            path.display(),
            feature.status
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
    })
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
    out.push_str("| Feature | Title | Area | Priority | Status | Implemented | Validated | Last AI review |\n");
    out.push_str("|---|---|---|---|---|---|---|---|\n");
    for feature in features {
        out.push_str(&format!(
            "| `{}` | {} | {} | {} | {} | {} | {} | {} |\n",
            feature.id,
            feature.title,
            feature.area,
            feature.priority,
            feature.status,
            yes_no(feature.implemented),
            yes_no(feature.validated),
            feature.last_ai_review
        ));
    }
    out
}

fn yes_no(value: bool) -> &'static str {
    if value {
        "yes"
    } else {
        "no"
    }
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
    fn read_feature_parses_typed_metadata() {
        let root = temp_feature_root("read-feature");
        write_feature(
            &root,
            "example.feature",
            r#"id = "example.feature"
title = "Example"
area = "infrastructure"
priority = "P1"
status = "planned"
implemented = false
validated = true
last_ai_review = "2026-06-17"
description = "Example feature."
depends_on = ["core.graph"]
"#,
        );

        let feature = read_feature(&root.join("example.feature").join("feature.toml"))
            .expect("feature should parse");

        assert_eq!(feature.id, "example.feature");
        assert_eq!(feature.priority, "P1");
        assert!(!feature.implemented);
        assert!(feature.validated);
        assert_eq!(feature.depends_on, vec!["core.graph"]);
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn read_feature_rejects_bad_boolean_priority_status_and_directory_mismatch() {
        let root = temp_feature_root("bad-feature");
        write_feature(
            &root,
            "bad.bool",
            r#"id = "bad.bool"
title = "Bad"
area = "infrastructure"
priority = "P0"
status = "planned"
implemented = maybe
validated = false
last_ai_review = "2026-06-17"
description = "Bad feature."
depends_on = []
"#,
        );
        assert!(read_feature(&root.join("bad.bool").join("feature.toml")).is_err());

        write_feature(
            &root,
            "bad.priority",
            r#"id = "bad.priority"
title = "Bad"
area = "infrastructure"
priority = "P9"
status = "planned"
implemented = false
validated = false
last_ai_review = "2026-06-17"
description = "Bad feature."
depends_on = []
"#,
        );
        assert!(read_feature(&root.join("bad.priority").join("feature.toml")).is_err());

        write_feature(
            &root,
            "bad.status",
            r#"id = "bad.status"
title = "Bad"
area = "infrastructure"
priority = "P0"
status = "maybe"
implemented = false
validated = false
last_ai_review = "2026-06-17"
description = "Bad feature."
depends_on = []
"#,
        );
        assert!(read_feature(&root.join("bad.status").join("feature.toml")).is_err());

        write_feature(
            &root,
            "real.id",
            r#"id = "wrong.id"
title = "Bad"
area = "infrastructure"
priority = "P0"
status = "planned"
implemented = false
validated = false
last_ai_review = "2026-06-17"
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
priority = "P0"
status = "implemented"
implemented = true
validated = false
last_ai_review = "2026-06-17"
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
priority = "P0"
status = "planned"
implemented = false
validated = false
last_ai_review = "2026-06-17"
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
priority = "P0"
status = "planned"
implemented = false
validated = false
last_ai_review = "2026-06-17"
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
                priority: "P0".to_owned(),
                status: "planned".to_owned(),
                implemented: false,
                validated: false,
                last_ai_review: "2026-06-17".to_owned(),
                description: "A feature.".to_owned(),
                depends_on: Vec::new(),
            },
            Feature {
                id: "z.feature".to_owned(),
                title: "Zed".to_owned(),
                area: "io".to_owned(),
                priority: "P1".to_owned(),
                status: "implemented".to_owned(),
                implemented: true,
                validated: false,
                last_ai_review: "2026-06-18".to_owned(),
                description: "Z feature.".to_owned(),
                depends_on: vec!["a.feature".to_owned()],
            },
        ];

        let dashboard = render_dashboard(&features);

        assert!(dashboard
            .contains("| `a.feature` | Aye | core | P0 | planned | no | no | 2026-06-17 |"));
        assert!(dashboard
            .contains("| `z.feature` | Zed | io | P1 | implemented | yes | no | 2026-06-18 |"));
        assert!(dashboard.ends_with('\n'));
    }

    #[test]
    fn validation_manifest_path_is_feature_scoped() {
        assert_eq!(
            validation_manifest_path("core.graph"),
            PathBuf::from("validation/features/core.graph/validation.toml")
        );
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
    }
}
