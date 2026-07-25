use crate::*;

mod compare;
mod digest;
mod implementation;
mod manifest;
mod progress;
mod results;

pub(crate) use compare::*;
pub(crate) use digest::*;
pub(crate) use implementation::*;
pub(crate) use manifest::*;
pub(crate) use progress::*;
pub(crate) use results::*;

pub(crate) fn benchmark(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    benchmark_args(&args)?;
    let feature_selector = value_after_flag(&args, "--feature")
        .ok_or_else(|| boxed_error("missing required flag: --feature FEATURE_ID|all"))?;
    let corpus_selector = benchmark_corpus_selector(&args);
    let fixture_selector = value_after_flag(&args, "--fixture");
    let accept_goldens = args
        .iter()
        .any(|arg| arg == "--accept-implementation-goldens");
    if fixture_selector.is_some()
        && (feature_selector == "all" || corpus_selector == "all" || corpus_selector == "baseline")
    {
        return Err(boxed_error(
            "--fixture requires one concrete --feature and explicit --corpus",
        ));
    }
    if accept_goldens
        && (feature_selector == "all" || corpus_selector == "all" || corpus_selector == "baseline")
    {
        return Err(boxed_error(
            "--accept-implementation-goldens requires one concrete --feature and explicit --corpus",
        ));
    }

    let features = read_features()?;
    if feature_selector != "all"
        && !features
            .iter()
            .any(|candidate| candidate.id == feature_selector)
    {
        return Err(boxed_error(format!("unknown feature: {feature_selector}")));
    }
    if corpus_selector != "all"
        && corpus_selector != "baseline"
        && !is_known_corpus(corpus_selector)
    {
        return Err(boxed_error(format!("unknown corpus: {corpus_selector}")));
    }

    let targets = benchmark_targets(&features, feature_selector, corpus_selector);
    if targets.is_empty() {
        if let Some(error) = concrete_missing_manifest_error(feature_selector, corpus_selector) {
            return Err(boxed_error(error));
        }
        println!(
            "no applicable benchmark targets for feature `{feature_selector}` and corpus `{corpus_selector}`"
        );
        return Ok(());
    }

    let jobs = benchmark_jobs(&args)?;
    let mut progress = BenchmarkProgress::start(targets.len(), jobs);
    let mut results = read_benchmark_results(&features)?;
    let mut selected_corpora = BTreeSet::new();
    let mut differences = Vec::new();
    let mut matched = 0;
    let hash_cache = BenchmarkHashCache::default();

    for (target_index, (feature, corpus)) in targets.into_iter().enumerate() {
        progress.target_start(target_index + 1, &feature.id, &corpus);
        selected_corpora.insert(corpus.clone());
        let manifest_path = benchmark_manifest_path(&feature.id, &corpus);
        let scope = fixture_selector
            .map(|fixture| format!("fixture:{fixture}"))
            .unwrap_or_else(|| "full".to_owned());
        let manifest_digest = hash_normalized_file(&manifest_path).ok();
        let result = run_target(
            feature,
            &corpus,
            &manifest_path,
            fixture_selector,
            accept_goldens,
            jobs,
            &scope,
            manifest_digest.clone(),
            &hash_cache,
            &mut progress,
        )
        .unwrap_or_else(|error| BenchmarkRun {
            outcome: BenchmarkResultOutcome::Error,
            scope: scope.clone(),
            fixture_count: 0,
            compared_count: 0,
            difference_count: 0,
            first_detail: Some(error.to_string()),
            reference_tool: None,
            reference_version: None,
            manifest_digest,
            input_digest: None,
        });

        match result.outcome {
            BenchmarkResultOutcome::Match => {
                matched += 1;
                progress.target_match(result.compared_count, result.fixture_count);
            }
            BenchmarkResultOutcome::Differences => {
                progress.target_differences(
                    result.difference_count,
                    result.compared_count,
                    result.fixture_count,
                );
                differences.push(format!(
                    "{} [{corpus}]: {} difference(s); first: {}",
                    feature.id,
                    result.difference_count,
                    result.first_detail.as_deref().unwrap_or("not recorded")
                ));
            }
            BenchmarkResultOutcome::Error => {
                progress.target_error();
                differences.push(format!(
                    "{} [{corpus}]: benchmark error: {}",
                    feature.id,
                    result.first_detail.as_deref().unwrap_or("not recorded")
                ));
            }
        }

        let snapshot = BenchmarkResult::from_run(result)?;
        results
            .entry(feature.id.clone())
            .or_insert_with(|| BenchmarkResults::new(&feature.id))
            .corpora
            .insert(corpus, snapshot);
    }

    write_benchmark_results(&results, &selected_corpora)?;
    let corpus_info = read_dashboard_corpus_info()?;
    let rendered = render_dashboard(&features, &results, &corpus_info);
    write_atomic_text(Path::new(DASHBOARD_PATH), &rendered)?;
    println!("recorded benchmark results and refreshed dashboard");
    println!("benchmark matched {matched} target(s)");

    if !differences.is_empty() {
        for difference in &differences {
            eprintln!("benchmark observation: {difference}");
        }
        return Err(boxed_error(format!(
            "{} benchmark target(s) reported differences or errors",
            differences.len()
        )));
    }
    Ok(())
}

pub(crate) fn concrete_missing_manifest_error(
    feature_selector: &str,
    corpus_selector: &str,
) -> Option<String> {
    (feature_selector != "all"
        && corpus_selector != "all"
        && corpus_selector != "baseline")
        .then(|| {
            format!(
                "no benchmark manifest for feature `{feature_selector}` and corpus `{corpus_selector}`"
            )
        })
}

#[allow(clippy::too_many_arguments)]
fn run_target(
    feature: &Feature,
    corpus: &str,
    manifest_path: &Path,
    fixture_selector: Option<&str>,
    accept_goldens: bool,
    jobs: usize,
    scope: &str,
    manifest_digest: Option<String>,
    hash_cache: &BenchmarkHashCache,
    progress: &mut BenchmarkProgress,
) -> Result<BenchmarkRun, Box<dyn Error>> {
    let mut manifest = read_benchmark_manifest(manifest_path)?;
    if manifest.feature_id != feature.id {
        return Err(boxed_error(format!(
            "{} declares feature_id `{}`, expected `{}`",
            manifest_path.display(),
            manifest.feature_id,
            feature.id
        )));
    }
    if manifest.corpus_id != corpus {
        return Err(boxed_error(format!(
            "{} declares corpus_id `{}`, expected `{corpus}`",
            manifest_path.display(),
            manifest.corpus_id
        )));
    }
    if let Some(fixture) = fixture_selector {
        if !manifest
            .fixtures
            .iter()
            .any(|candidate| candidate == fixture)
        {
            return Err(boxed_error(format!(
                "{} does not declare fixture `{fixture}`",
                manifest_path.display()
            )));
        }
        manifest.fixtures.retain(|candidate| candidate == fixture);
    }
    progress.manifest(&manifest.reference_tool, &manifest.reference_version);
    check_comparison_mode(manifest_path, &manifest)?;
    if manifest.fixtures.is_empty() {
        return Err(boxed_error(format!(
            "{} must list at least one benchmark fixture",
            manifest_path.display()
        )));
    }
    check_manifest_paths(manifest_path, &manifest)?;
    if accept_goldens {
        accept_implementation_goldens(manifest_path, &manifest, jobs)?;
        println!(
            "  accepted {} reviewed implementation golden(s)",
            manifest.fixtures.len()
        );
    }

    let input_digest =
        build_benchmark_input_digest_cached(Path::new("."), manifest_path, &manifest, hash_cache)?;
    let worker_count = benchmark_worker_count(jobs, manifest.fixtures.len());
    let fixture_progress = FixtureProgress::start(manifest.fixtures.len(), worker_count);
    let comparison_result = compare_golden_outputs_cached(
        manifest_path,
        &manifest,
        jobs,
        Some(&fixture_progress),
        hash_cache,
    );
    fixture_progress.finish();
    let comparison = comparison_result?;
    let outcome = if comparison.difference_count == 0 {
        BenchmarkResultOutcome::Match
    } else {
        BenchmarkResultOutcome::Differences
    };
    Ok(BenchmarkRun {
        outcome,
        scope: scope.to_owned(),
        fixture_count: manifest.fixtures.len(),
        compared_count: comparison.compared_count,
        difference_count: comparison.difference_count,
        first_detail: comparison.first_difference,
        reference_tool: Some(manifest.reference_tool),
        reference_version: Some(manifest.reference_version),
        manifest_digest,
        input_digest: Some(input_digest),
    })
}

pub(crate) fn benchmark_corpus_selector(args: &[String]) -> &str {
    value_after_flag(args, "--corpus").unwrap_or("baseline")
}
