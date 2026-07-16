use crate::*;

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ValidationManifest {
    pub(crate) feature_id: String,
    pub(crate) corpus_id: String,
    pub(crate) reference_tool: String,
    pub(crate) reference_version: String,
    pub(crate) comparison_mode: String,
    #[serde(default)]
    pub(crate) fixtures: Vec<String>,
    #[serde(default, rename = "notes")]
    pub(crate) _notes: Vec<String>,
}

#[derive(Debug)]
pub(crate) struct ValidationRun {
    pub(crate) fixture_count: usize,
    pub(crate) compared_count: usize,
    pub(crate) reference_tool: String,
    pub(crate) reference_version: String,
    pub(crate) manifest_hash: String,
    pub(crate) evidence: ValidationEvidence,
}

#[derive(Debug)]
pub(crate) struct FailedValidationRun {
    pub(crate) fixture_count: usize,
    pub(crate) compared_count: usize,
    pub(crate) failed_count: usize,
    pub(crate) first_failure: String,
    pub(crate) reference_tool: String,
    pub(crate) reference_version: String,
    pub(crate) manifest_hash: String,
}

#[derive(Debug)]
pub(crate) enum ValidationOutcome {
    Passed(ValidationRun),
    Failed(FailedValidationRun),
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ValidationEvidence {
    pub(crate) schema_version: u32,
    pub(crate) comparison_mode: String,
    pub(crate) inputs: Vec<EvidenceInput>,
    pub(crate) sha256: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct EvidenceInput {
    pub(crate) path: String,
    pub(crate) sha256: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct CorpusStatus {
    pub(crate) passed: bool,
    pub(crate) fixture_count: usize,
    pub(crate) compared_count: usize,
    pub(crate) reference_tool: String,
    pub(crate) reference_version: String,
    pub(crate) manifest_hash: String,
    #[serde(default)]
    pub(crate) failed_count: usize,
    #[serde(default)]
    pub(crate) first_failure: Option<String>,
    #[serde(default)]
    pub(crate) evidence_schema_version: Option<u32>,
    #[serde(default)]
    pub(crate) evidence_hash: Option<String>,
    #[serde(default)]
    pub(crate) evidence_inputs: Vec<EvidenceInput>,
    pub(crate) validated_at_unix: u64,
}

impl CorpusStatus {
    pub(crate) fn from_run(
        run: ValidationRun,
        existing: Option<&CorpusStatus>,
    ) -> Result<Self, Box<dyn Error>> {
        let unchanged_evidence = existing.and_then(|status| status.evidence_hash.as_deref())
            == Some(run.evidence.sha256.as_str());
        let validated_at_unix = if unchanged_evidence {
            existing
                .map(|status| status.validated_at_unix)
                .unwrap_or(SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs())
        } else {
            SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs()
        };
        Ok(Self {
            passed: true,
            fixture_count: run.fixture_count,
            compared_count: run.compared_count,
            reference_tool: run.reference_tool,
            reference_version: run.reference_version,
            manifest_hash: run.manifest_hash,
            failed_count: 0,
            first_failure: None,
            evidence_schema_version: Some(run.evidence.schema_version),
            evidence_hash: Some(run.evidence.sha256),
            evidence_inputs: run.evidence.inputs,
            validated_at_unix,
        })
    }

    pub(crate) fn from_failed_run(run: FailedValidationRun) -> Result<Self, Box<dyn Error>> {
        Ok(Self {
            passed: false,
            fixture_count: run.fixture_count,
            compared_count: run.compared_count,
            reference_tool: run.reference_tool,
            reference_version: run.reference_version,
            manifest_hash: run.manifest_hash,
            failed_count: run.failed_count,
            first_failure: Some(run.first_failure),
            evidence_schema_version: None,
            evidence_hash: None,
            evidence_inputs: Vec::new(),
            validated_at_unix: SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs(),
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ValidationStatus {
    pub(crate) feature_id: String,
    pub(crate) corpora: BTreeMap<String, CorpusStatus>,
}

impl ValidationStatus {
    pub(crate) fn new(feature_id: &str) -> Self {
        Self {
            feature_id: feature_id.to_owned(),
            corpora: BTreeMap::new(),
        }
    }
}

pub(crate) fn is_known_corpus(corpus: &str) -> bool {
    validation_corpus(corpus).is_some()
}

pub(crate) fn validation_corpus(corpus: &str) -> Option<&'static ValidationCorpus> {
    VALIDATION_CORPORA
        .iter()
        .find(|candidate| candidate.id == corpus)
}

pub(crate) fn validation_manifest_path(feature: &str, corpus: &str) -> PathBuf {
    Path::new("validation")
        .join("corpora")
        .join(corpus)
        .join("features")
        .join(format!("{feature}.toml"))
}

pub(crate) fn validation_manifest_path_from(root: &Path, feature: &str, corpus: &str) -> PathBuf {
    root.join(validation_manifest_path(feature, corpus))
}

pub(crate) fn read_validation_manifest(path: &Path) -> Result<ValidationManifest, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    toml::from_str(&text).map_err(|error| boxed_error(format!("{}: {error}", path.display())))
}

pub(crate) fn validate_comparison_mode(
    manifest_path: &Path,
    manifest: &ValidationManifest,
) -> Result<(), Box<dyn Error>> {
    if manifest.comparison_mode != COMPARISON_MODE_IMPLEMENTATION_GOLDEN {
        return Err(boxed_error(format!(
            "{} uses unsupported comparison_mode `{}`",
            manifest_path.display(),
            manifest.comparison_mode
        )));
    }
    Ok(())
}

pub(crate) fn validation_targets<'a>(
    features: &'a [Feature],
    feature_selector: &str,
    corpus_selector: &str,
) -> Vec<(&'a Feature, String)> {
    validation_targets_from(Path::new("."), features, feature_selector, corpus_selector)
}

pub(crate) fn validation_targets_from<'a>(
    root: &Path,
    features: &'a [Feature],
    feature_selector: &str,
    corpus_selector: &str,
) -> Vec<(&'a Feature, String)> {
    if corpus_selector == "all" {
        let mut targets = Vec::new();
        for feature in features {
            if feature_selector != "all" && feature.id != feature_selector {
                continue;
            }
            for corpus in VALIDATION_CORPORA {
                let required = feature
                    .validation_required
                    .iter()
                    .any(|required| required == corpus.id);
                let manifest_exists =
                    validation_manifest_path_from(root, &feature.id, corpus.id).exists();
                if required || (feature.status.has_implementation() && manifest_exists) {
                    targets.push((feature, corpus.id.to_owned()));
                }
            }
        }
        return targets;
    }

    if corpus_selector != "all"
        && validation_corpus(corpus_selector).is_some_and(|corpus| corpus.local_only)
    {
        return features
            .iter()
            .filter(|feature| feature_selector == "all" || feature.id == feature_selector)
            .filter(|feature| {
                feature_selector != "all"
                    || validation_manifest_path_from(root, &feature.id, corpus_selector).exists()
            })
            .map(|feature| (feature, corpus_selector.to_owned()))
            .collect();
    }

    let mut targets = Vec::new();
    for feature in features {
        if feature_selector != "all" && feature.id != feature_selector {
            continue;
        }
        for corpus in &feature.validation_required {
            if corpus_selector == "all" || corpus == corpus_selector {
                targets.push((feature, corpus.clone()));
            }
        }
    }
    targets
}
