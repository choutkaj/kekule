use crate::*;

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct BenchmarkManifest {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "kebab-case")]
pub(crate) enum BenchmarkResultOutcome {
    Match,
    Differences,
    Error,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct BenchmarkInputDigest {
    pub(crate) schema_version: u32,
    pub(crate) input_count: usize,
    pub(crate) sha256: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct BenchmarkRun {
    pub(crate) outcome: BenchmarkResultOutcome,
    pub(crate) scope: String,
    pub(crate) fixture_count: usize,
    pub(crate) compared_count: usize,
    pub(crate) difference_count: usize,
    pub(crate) first_detail: Option<String>,
    pub(crate) reference_tool: Option<String>,
    pub(crate) reference_version: Option<String>,
    pub(crate) manifest_digest: Option<String>,
    pub(crate) input_digest: Option<BenchmarkInputDigest>,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct BenchmarkResult {
    pub(crate) outcome: BenchmarkResultOutcome,
    pub(crate) scope: String,
    pub(crate) fixture_count: usize,
    pub(crate) compared_count: usize,
    #[serde(default)]
    pub(crate) difference_count: usize,
    #[serde(default)]
    pub(crate) first_detail: Option<String>,
    #[serde(default)]
    pub(crate) reference_tool: Option<String>,
    #[serde(default)]
    pub(crate) reference_version: Option<String>,
    #[serde(default)]
    pub(crate) manifest_digest: Option<String>,
    #[serde(default)]
    pub(crate) input_digest_schema_version: Option<u32>,
    #[serde(default)]
    pub(crate) input_digest: Option<String>,
    #[serde(default)]
    pub(crate) input_count: usize,
    #[serde(default)]
    pub(crate) legacy_source: Option<String>,
    pub(crate) benchmarked_at_unix: u64,
}

impl BenchmarkResult {
    pub(crate) fn from_run(run: BenchmarkRun) -> Result<Self, Box<dyn Error>> {
        let (input_digest_schema_version, input_digest, input_count) =
            run.input_digest.map_or((None, None, 0), |digest| {
                (
                    Some(digest.schema_version),
                    Some(digest.sha256),
                    digest.input_count,
                )
            });
        Ok(Self {
            outcome: run.outcome,
            scope: run.scope,
            fixture_count: run.fixture_count,
            compared_count: run.compared_count,
            difference_count: run.difference_count,
            first_detail: run.first_detail,
            reference_tool: run.reference_tool,
            reference_version: run.reference_version,
            manifest_digest: run.manifest_digest,
            input_digest_schema_version,
            input_digest,
            input_count,
            legacy_source: None,
            benchmarked_at_unix: SystemTime::now().duration_since(UNIX_EPOCH)?.as_secs(),
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct BenchmarkResults {
    pub(crate) feature_id: String,
    pub(crate) corpora: BTreeMap<String, BenchmarkResult>,
}

impl BenchmarkResults {
    pub(crate) fn new(feature_id: &str) -> Self {
        Self {
            feature_id: feature_id.to_owned(),
            corpora: BTreeMap::new(),
        }
    }
}

pub(crate) fn is_known_corpus(corpus: &str) -> bool {
    benchmark_corpus(corpus).is_some()
}

pub(crate) fn benchmark_corpus(corpus: &str) -> Option<&'static BenchmarkCorpus> {
    BENCHMARK_CORPORA
        .iter()
        .find(|candidate| candidate.id == corpus)
}

pub(crate) fn benchmark_manifest_path(feature: &str, corpus: &str) -> PathBuf {
    Path::new("benchmarks")
        .join("corpora")
        .join(corpus)
        .join("features")
        .join(format!("{feature}.toml"))
}

pub(crate) fn benchmark_manifest_path_from(root: &Path, feature: &str, corpus: &str) -> PathBuf {
    root.join(benchmark_manifest_path(feature, corpus))
}

pub(crate) fn read_benchmark_manifest(path: &Path) -> Result<BenchmarkManifest, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    toml::from_str(&text).map_err(|error| boxed_error(format!("{}: {error}", path.display())))
}

pub(crate) fn check_comparison_mode(
    manifest_path: &Path,
    manifest: &BenchmarkManifest,
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

pub(crate) fn benchmark_targets<'a>(
    features: &'a [Feature],
    feature_selector: &str,
    corpus_selector: &str,
) -> Vec<(&'a Feature, String)> {
    benchmark_targets_from(Path::new("."), features, feature_selector, corpus_selector)
}

pub(crate) fn benchmark_targets_from<'a>(
    root: &Path,
    features: &'a [Feature],
    feature_selector: &str,
    corpus_selector: &str,
) -> Vec<(&'a Feature, String)> {
    let corpora = BENCHMARK_CORPORA
        .iter()
        .filter(|corpus| match corpus_selector {
            "all" => true,
            "baseline" => corpus.default,
            concrete => corpus.id == concrete,
        })
        .collect::<Vec<_>>();
    let mut targets = Vec::new();
    for feature in features {
        if (feature_selector == "all" || feature.id == feature_selector)
            && feature.status.has_implementation()
        {
            for corpus in &corpora {
                if benchmark_manifest_path_from(root, &feature.id, corpus.id).exists() {
                    targets.push((feature, corpus.id.to_owned()));
                }
            }
        }
    }
    targets
}
