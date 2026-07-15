use crate::*;

pub(crate) fn list_features() -> Result<(), Box<dyn Error>> {
    let features = read_features()?;
    let statuses = read_validation_statuses(&features)?;
    for feature in &features {
        println!(
            "{}\t{}\tv{}\timplemented={}\tvalidated={}",
            feature.id,
            feature.area,
            feature.version,
            feature.implemented,
            overall_validated(feature, statuses.get(&feature.id))
        );
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Feature {
    pub(crate) id: String,
    pub(crate) title: String,
    pub(crate) area: String,
    pub(crate) version: u32,
    pub(crate) implemented: bool,
    pub(crate) validated: bool,
    pub(crate) description: String,
    pub(crate) depends_on: Vec<String>,
    pub(crate) validation_required: Vec<String>,
}

pub(crate) fn read_features() -> Result<Vec<Feature>, Box<dyn Error>> {
    read_features_from(Path::new("features"))
}

pub(crate) fn read_features_from(root: &Path) -> Result<Vec<Feature>, Box<dyn Error>> {
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

pub(crate) fn is_hidden_or_template(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .map(|name| name.starts_with('_'))
        .unwrap_or(false)
}

pub(crate) fn read_feature(path: &Path) -> Result<Feature, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    let feature: Feature = toml::from_str(&text)
        .map_err(|error| boxed_error(format!("{}: {error}", path.display())))?;
    validate_feature(&feature, path)?;
    Ok(feature)
}

pub(crate) fn validate_feature(feature: &Feature, path: &Path) -> Result<(), Box<dyn Error>> {
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
    if feature.version == 0 {
        return Err(boxed_error(format!(
            "{} has invalid zero `version` value",
            path.display()
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
    let mut seen_corpora = BTreeSet::new();
    for corpus in &feature.validation_required {
        if !is_known_corpus(corpus) {
            return Err(boxed_error(format!(
                "{} requires unknown validation corpus `{corpus}`",
                path.display()
            )));
        }
        if validation_corpus(corpus).is_some_and(|registered| registered.local_only) {
            return Err(boxed_error(format!(
                "{} lists local-only validation corpus `{corpus}` in `validation_required`; local-only corpora may be run explicitly but cannot determine repository-wide validation state",
                path.display()
            )));
        }
        if !seen_corpora.insert(corpus) {
            return Err(boxed_error(format!(
                "{} lists validation corpus `{corpus}` more than once",
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

pub(crate) fn validate_feature_set(features: &[Feature]) -> Result<(), Box<dyn Error>> {
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
