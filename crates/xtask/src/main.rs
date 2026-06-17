use std::collections::BTreeMap;
use std::env;
use std::error::Error;
use std::fs;
use std::path::Path;
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
    let path = Path::new("features").join(feature);
    if !path.join("feature.toml").exists() {
        return Err(boxed_error(format!("unknown feature: {feature}")));
    }

    println!("validation harness found feature `{feature}`");
    println!("reference validation is feature-specific and may not be implemented yet");
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

#[derive(Debug, Clone, Default)]
struct Feature {
    id: String,
    title: String,
    area: String,
    priority: String,
    status: String,
    implemented: String,
    validated: String,
    last_ai_review: String,
}

fn read_features() -> Result<Vec<Feature>, Box<dyn Error>> {
    let mut features = Vec::new();
    for entry in fs::read_dir("features")? {
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
    Ok(Feature {
        id: required(&map, "id", path)?,
        title: required(&map, "title", path)?,
        area: required(&map, "area", path)?,
        priority: required(&map, "priority", path)?,
        status: required(&map, "status", path)?,
        implemented: required(&map, "implemented", path)?,
        validated: required(&map, "validated", path)?,
        last_ai_review: required(&map, "last_ai_review", path)?,
    })
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
            yes_no(&feature.implemented),
            yes_no(&feature.validated),
            feature.last_ai_review
        ));
    }
    out
}

fn yes_no(value: &str) -> &'static str {
    match value {
        "true" => "yes",
        "false" => "no",
        _ => "unknown",
    }
}

fn boxed_error(message: impl Into<String>) -> Box<dyn Error> {
    std::io::Error::other(message.into()).into()
}
