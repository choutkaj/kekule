use crate::*;

pub(crate) fn list_features() -> Result<(), Box<dyn Error>> {
    let features = read_features()?;
    validate_required_manifests(&features)?;
    for feature in &features {
        println!(
            "{}\t{}\tv{}\tstatus={}",
            feature.id,
            feature.area,
            feature.version,
            feature.status.as_str()
        );
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deserialize, Serialize)]
#[serde(rename_all = "kebab-case")]
pub(crate) enum FeatureDomain {
    SmallMolecule,
    Macromolecule,
    Infrastructure,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deserialize, Serialize)]
#[serde(rename_all = "kebab-case")]
pub(crate) enum FeatureStatus {
    Planned,
    Experimental,
    Supported,
    Deprecated,
}

impl FeatureStatus {
    pub(crate) const fn as_str(self) -> &'static str {
        match self {
            Self::Planned => "planned",
            Self::Experimental => "experimental",
            Self::Supported => "supported",
            Self::Deprecated => "deprecated",
        }
    }

    pub(crate) const fn sort_value(self) -> u8 {
        match self {
            Self::Planned => 0,
            Self::Experimental => 1,
            Self::Supported => 2,
            Self::Deprecated => 3,
        }
    }

    pub(crate) const fn has_implementation(self) -> bool {
        !matches!(self, Self::Planned)
    }

    pub(crate) const fn permits_dependency(self, dependency: Self) -> bool {
        match self {
            Self::Planned => true,
            Self::Experimental => {
                matches!(dependency, Self::Experimental | Self::Supported)
            }
            Self::Supported => matches!(dependency, Self::Supported),
            Self::Deprecated => dependency.has_implementation(),
        }
    }

    pub(crate) const fn dependency_contract(self) -> &'static str {
        match self {
            Self::Planned => "any registered status",
            Self::Experimental => "`experimental` or `supported` features",
            Self::Supported => "`supported` features",
            Self::Deprecated => "implemented features",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Feature {
    pub(crate) id: String,
    pub(crate) title: String,
    pub(crate) area: String,
    pub(crate) domains: Vec<FeatureDomain>,
    pub(crate) version: u32,
    pub(crate) status: FeatureStatus,
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

pub(crate) fn validate_required_manifests(features: &[Feature]) -> Result<(), Box<dyn Error>> {
    validate_required_manifests_from(Path::new("."), features)
}

pub(crate) fn validate_required_manifests_from(
    root: &Path,
    features: &[Feature],
) -> Result<(), Box<dyn Error>> {
    for feature in features {
        for corpus in &feature.validation_required {
            let path = validation_manifest_path_from(root, &feature.id, corpus);
            if !path.exists() {
                return Err(boxed_error(format!(
                    "{} requires validation corpus `{corpus}` but manifest {} is missing",
                    feature.id,
                    path.display()
                )));
            }
        }
    }
    Ok(())
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
    if feature.domains.is_empty() {
        return Err(boxed_error(format!(
            "{} has empty required field `domains`",
            path.display()
        )));
    }
    let mut seen_domains = BTreeSet::new();
    for domain in &feature.domains {
        if !seen_domains.insert(domain) {
            return Err(boxed_error(format!(
                "{} lists feature domain `{domain:?}` more than once",
                path.display()
            )));
        }
    }
    if feature.domains.contains(&FeatureDomain::Infrastructure) && feature.domains.len() != 1 {
        return Err(boxed_error(format!(
            "{} combines the `infrastructure` domain with a chemistry domain",
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
    let mut seen_dependencies = BTreeSet::new();
    for dependency in &feature.depends_on {
        if dependency.trim().is_empty() {
            return Err(boxed_error(format!(
                "{} has an empty feature ID in `depends_on`",
                path.display()
            )));
        }
        if dependency == &feature.id {
            return Err(boxed_error(format!(
                "{} declares a self-dependency in `depends_on`",
                path.display()
            )));
        }
        if !seen_dependencies.insert(dependency) {
            return Err(boxed_error(format!(
                "{} lists dependency `{dependency}` more than once",
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
        if !seen_corpora.insert(corpus) {
            return Err(boxed_error(format!(
                "{} lists validation corpus `{corpus}` more than once",
                path.display()
            )));
        }
    }
    if feature.status == FeatureStatus::Planned && !feature.validation_required.is_empty() {
        return Err(boxed_error(format!(
            "{} has status `planned` but declares required validation; planned features have no implementation to validate",
            path.display()
        )));
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
    let mut seen = BTreeMap::<&str, &Feature>::new();
    for feature in features {
        if seen.insert(feature.id.as_str(), feature).is_some() {
            return Err(boxed_error(format!(
                "duplicate feature id `{}`",
                feature.id
            )));
        }
    }
    for feature in features {
        let mut dependencies = BTreeSet::new();
        for dependency in &feature.depends_on {
            if dependency == &feature.id {
                return Err(boxed_error(format!(
                    "feature `{}` depends on itself",
                    feature.id
                )));
            }
            if !dependencies.insert(dependency.as_str()) {
                return Err(boxed_error(format!(
                    "feature `{}` lists dependency `{dependency}` more than once",
                    feature.id
                )));
            }
            let Some(dependency_feature) = seen.get(dependency.as_str()) else {
                return Err(boxed_error(format!(
                    "feature `{}` depends on unknown feature `{dependency}`",
                    feature.id
                )));
            };
            if !feature.status.permits_dependency(dependency_feature.status) {
                return Err(boxed_error(format!(
                    "feature `{}` has status `{}` but depends on `{dependency}` with status `{}`; `{}` features may depend only on {}",
                    feature.id,
                    feature.status.as_str(),
                    dependency_feature.status.as_str(),
                    feature.status.as_str(),
                    feature.status.dependency_contract()
                )));
            }
        }
    }
    feature_topological_order(features)?;
    Ok(())
}

pub(crate) fn feature_dependency_layers(
    features: &[Feature],
) -> Result<Vec<Vec<&Feature>>, Box<dyn Error>> {
    let order = feature_topological_order(features)?;
    let indexes = features
        .iter()
        .enumerate()
        .map(|(index, feature)| (feature.id.as_str(), index))
        .collect::<BTreeMap<_, _>>();
    let mut depths = vec![0usize; features.len()];
    for index in order {
        depths[index] = features[index]
            .depends_on
            .iter()
            .filter_map(|dependency| indexes.get(dependency.as_str()))
            .map(|dependency| depths[*dependency] + 1)
            .max()
            .unwrap_or(0);
    }
    let mut layers = vec![Vec::new(); depths.iter().copied().max().unwrap_or(0) + 1];
    for (feature, depth) in features.iter().zip(depths) {
        layers[depth].push(feature);
    }
    for layer in &mut layers {
        layer.sort_by(|left, right| left.id.cmp(&right.id));
    }
    Ok(layers)
}

fn feature_topological_order(features: &[Feature]) -> Result<Vec<usize>, Box<dyn Error>> {
    let indexes = features
        .iter()
        .enumerate()
        .map(|(index, feature)| (feature.id.as_str(), index))
        .collect::<BTreeMap<_, _>>();
    let mut indegrees = features
        .iter()
        .map(|feature| feature.depends_on.len())
        .collect::<Vec<_>>();
    let mut dependents = vec![Vec::new(); features.len()];
    for (dependent, feature) in features.iter().enumerate() {
        for dependency in &feature.depends_on {
            let Some(dependency) = indexes.get(dependency.as_str()).copied() else {
                return Err(boxed_error(format!(
                    "feature `{}` depends on unknown feature `{dependency}`",
                    feature.id
                )));
            };
            dependents[dependency].push(dependent);
        }
    }
    for values in &mut dependents {
        values.sort_by(|left, right| features[*left].id.cmp(&features[*right].id));
    }
    let mut ready = features
        .iter()
        .enumerate()
        .filter(|(index, _)| indegrees[*index] == 0)
        .map(|(_, feature)| feature.id.as_str())
        .collect::<BTreeSet<_>>();
    let mut order = Vec::with_capacity(features.len());
    while let Some(id) = ready.pop_first() {
        let index = indexes[id];
        order.push(index);
        for dependent in &dependents[index] {
            indegrees[*dependent] -= 1;
            if indegrees[*dependent] == 0 {
                ready.insert(features[*dependent].id.as_str());
            }
        }
    }
    if order.len() != features.len() {
        let cycle = find_feature_dependency_cycle(features, &indexes)
            .unwrap_or_else(|| vec!["unknown cycle".to_owned()]);
        return Err(boxed_error(format!(
            "feature dependency graph contains a cycle: {}",
            cycle.join(" -> ")
        )));
    }
    Ok(order)
}

fn find_feature_dependency_cycle(
    features: &[Feature],
    indexes: &BTreeMap<&str, usize>,
) -> Option<Vec<String>> {
    fn visit(
        index: usize,
        features: &[Feature],
        indexes: &BTreeMap<&str, usize>,
        states: &mut [u8],
        stack: &mut Vec<usize>,
    ) -> Option<Vec<String>> {
        states[index] = 1;
        stack.push(index);
        let mut dependencies = features[index]
            .depends_on
            .iter()
            .filter_map(|dependency| indexes.get(dependency.as_str()).copied())
            .collect::<Vec<_>>();
        dependencies.sort_by(|left, right| features[*left].id.cmp(&features[*right].id));
        for dependency in dependencies {
            if states[dependency] == 0 {
                if let Some(cycle) = visit(dependency, features, indexes, states, stack) {
                    return Some(cycle);
                }
            } else if states[dependency] == 1 {
                let start = stack
                    .iter()
                    .position(|candidate| *candidate == dependency)
                    .unwrap_or(0);
                let mut cycle = stack[start..]
                    .iter()
                    .map(|candidate| features[*candidate].id.clone())
                    .collect::<Vec<_>>();
                cycle.push(features[dependency].id.clone());
                return Some(cycle);
            }
        }
        stack.pop();
        states[index] = 2;
        None
    }

    let mut states = vec![0u8; features.len()];
    let mut stack = Vec::new();
    for index in 0..features.len() {
        if states[index] == 0 {
            if let Some(cycle) = visit(index, features, indexes, &mut states, &mut stack) {
                return Some(cycle);
            }
        }
    }
    None
}
