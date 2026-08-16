use crate::*;

#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct BenchmarkManifest {
    // `feature_id` is retained only as legacy on-disk benchmark-schema vocabulary.
    #[serde(rename = "feature_id")]
    pub(crate) benchmark_id: String,
    pub(crate) corpus_id: String,
    pub(crate) reference_tool: String,
    pub(crate) reference_version: String,
    pub(crate) comparison_mode: String,
    #[serde(default)]
    pub(crate) fixtures: Vec<String>,
    #[serde(default, rename = "notes")]
    pub(crate) _notes: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct BenchmarkTarget {
    pub(crate) benchmark_id: String,
    pub(crate) corpus_id: String,
    pub(crate) manifest_path: PathBuf,
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

pub(crate) fn discover_benchmark_targets_from(
    root: &Path,
) -> Result<Vec<BenchmarkTarget>, Box<dyn Error>> {
    let corpora_root = root.join("benchmarks").join("corpora");
    let mut corpus_dirs = fs::read_dir(&corpora_root)?
        .map(|entry| entry.map(|entry| entry.path()))
        .collect::<Result<Vec<_>, _>>()?;
    corpus_dirs.sort();

    let mut targets = Vec::new();
    for corpus_dir in corpus_dirs {
        if !corpus_dir.is_dir() || !corpus_dir.join("corpus.toml").is_file() {
            continue;
        }
        let corpus_id = utf8_file_name(&corpus_dir, "corpus directory")?;
        let manifest_dir = corpus_dir.join("features");
        if !manifest_dir.is_dir() {
            continue;
        }
        let mut manifest_paths = fs::read_dir(&manifest_dir)?
            .map(|entry| entry.map(|entry| entry.path()))
            .collect::<Result<Vec<_>, _>>()?;
        manifest_paths.sort();
        for manifest_path in manifest_paths {
            if manifest_path
                .extension()
                .and_then(|extension| extension.to_str())
                != Some("toml")
            {
                continue;
            }
            let benchmark_id = manifest_path
                .file_stem()
                .and_then(|stem| stem.to_str())
                .ok_or_else(|| {
                    boxed_error(format!(
                        "{} has a non-UTF-8 benchmark manifest name",
                        manifest_path.display()
                    ))
                })?
                .to_owned();
            let manifest = read_benchmark_manifest(&manifest_path)?;
            if manifest.benchmark_id != benchmark_id {
                return Err(boxed_error(format!(
                    "{} declares legacy feature_id `{}`, expected benchmark ID `{benchmark_id}`",
                    manifest_path.display(),
                    manifest.benchmark_id
                )));
            }
            if manifest.corpus_id != corpus_id {
                return Err(boxed_error(format!(
                    "{} declares corpus_id `{}`, expected `{corpus_id}`",
                    manifest_path.display(),
                    manifest.corpus_id
                )));
            }
            targets.push(BenchmarkTarget {
                benchmark_id,
                corpus_id: corpus_id.clone(),
                manifest_path,
            });
        }
    }
    targets.sort_by(|left, right| {
        (&left.benchmark_id, &left.corpus_id).cmp(&(&right.benchmark_id, &right.corpus_id))
    });
    Ok(targets)
}

pub(crate) fn select_benchmark_targets(
    benchmark_selector: &str,
    corpus_selector: &str,
) -> Result<Vec<BenchmarkTarget>, Box<dyn Error>> {
    select_benchmark_targets_from(Path::new("."), benchmark_selector, corpus_selector)
}

pub(crate) fn select_benchmark_targets_from(
    root: &Path,
    benchmark_selector: &str,
    corpus_selector: &str,
) -> Result<Vec<BenchmarkTarget>, Box<dyn Error>> {
    let targets = discover_benchmark_targets_from(root)?;
    if benchmark_selector != "all"
        && !targets
            .iter()
            .any(|target| target.benchmark_id == benchmark_selector)
    {
        return Err(boxed_error(format!(
            "unknown benchmark: {benchmark_selector}"
        )));
    }
    let corpora = discover_corpus_ids_from(root)?;
    if corpus_selector != "all" && !corpora.iter().any(|corpus| corpus == corpus_selector) {
        return Err(boxed_error(format!("unknown corpus: {corpus_selector}")));
    }

    let selected = targets
        .into_iter()
        .filter(|target| {
            (benchmark_selector == "all" || target.benchmark_id == benchmark_selector)
                && (corpus_selector == "all" || target.corpus_id == corpus_selector)
        })
        .collect::<Vec<_>>();
    if selected.is_empty() && benchmark_selector != "all" && corpus_selector != "all" {
        return Err(boxed_error(format!(
            "no benchmark manifest for benchmark `{benchmark_selector}` and corpus `{corpus_selector}`"
        )));
    }
    Ok(selected)
}

fn utf8_file_name(path: &Path, kind: &str) -> Result<String, Box<dyn Error>> {
    path.file_name()
        .and_then(|name| name.to_str())
        .map(str::to_owned)
        .ok_or_else(|| boxed_error(format!("{} has a non-UTF-8 {kind} name", path.display())))
}
