use crate::*;

mod compare;
mod golden;
mod implementation;
mod manifest;
mod progress;

pub(crate) use compare::*;
pub(crate) use golden::*;
pub(crate) use implementation::*;
pub(crate) use manifest::*;
pub(crate) use progress::*;

pub(crate) fn benchmark(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    benchmark_args(&args)?;
    let benchmark_selector = value_after_flag(&args, "--benchmark")
        .ok_or_else(|| boxed_error("missing required flag: --benchmark BENCHMARK_ID|all"))?;
    let corpus_selector = value_after_flag(&args, "--corpus")
        .ok_or_else(|| boxed_error("missing required flag: --corpus CORPUS_ID|all"))?;
    let fixture_selector = value_after_flag(&args, "--fixture");
    let accept_goldens = args
        .iter()
        .any(|arg| arg == "--accept-implementation-goldens");
    if fixture_selector.is_some() && (benchmark_selector == "all" || corpus_selector == "all") {
        return Err(boxed_error(
            "--fixture requires one concrete --benchmark and --corpus",
        ));
    }
    if accept_goldens && (benchmark_selector == "all" || corpus_selector == "all") {
        return Err(boxed_error(
            "--accept-implementation-goldens requires one concrete --benchmark and --corpus",
        ));
    }

    let targets = select_benchmark_targets(benchmark_selector, corpus_selector)?;
    if targets.is_empty() {
        println!(
            "no benchmark manifests selected for benchmark `{benchmark_selector}` and corpus `{corpus_selector}`"
        );
        return Ok(());
    }

    let jobs = benchmark_jobs(&args)?;
    let target_count = targets.len();
    let mut progress = BenchmarkProgress::start(target_count, jobs);
    let mut failures = Vec::new();
    let mut matched_targets = 0usize;
    let mut matched_fixtures = 0usize;

    for (target_index, target) in targets.iter().enumerate() {
        progress.target_start(target_index + 1, &target.benchmark_id, &target.corpus_id);
        match run_target(
            target,
            fixture_selector,
            accept_goldens,
            jobs,
            &mut progress,
        ) {
            Ok(comparison) if comparison.difference_count == 0 => {
                matched_targets += 1;
                matched_fixtures += comparison.match_count;
                progress.target_match(comparison.match_count, comparison.match_count);
            }
            Ok(comparison) => {
                matched_fixtures += comparison.match_count;
                progress.target_differences(
                    comparison.difference_count,
                    comparison.match_count,
                    comparison.match_count + comparison.difference_count,
                );
                failures.push(format!(
                    "{} [{}]: {} difference(s); first: {}",
                    target.benchmark_id,
                    target.corpus_id,
                    comparison.difference_count,
                    comparison
                        .first_difference
                        .as_deref()
                        .unwrap_or("not recorded")
                ));
            }
            Err(error) => {
                progress.target_error();
                failures.push(format!(
                    "{} [{}]: benchmark error: {error}",
                    target.benchmark_id, target.corpus_id
                ));
            }
        }
    }

    println!("benchmark summary: {matched_fixtures} fixture matches; {matched_targets}/{target_count} targets fully matched");
    if !failures.is_empty() {
        for failure in &failures {
            eprintln!("benchmark result: {failure}");
        }
        return Err(boxed_error(format!(
            "{} benchmark target(s) reported differences or errors",
            failures.len()
        )));
    }
    Ok(())
}

pub(crate) fn run_target(
    target: &BenchmarkTarget,
    fixture_selector: Option<&str>,
    accept_goldens: bool,
    jobs: usize,
    progress: &mut BenchmarkProgress,
) -> Result<BenchmarkComparison, Box<dyn Error>> {
    let mut manifest = read_benchmark_manifest(&target.manifest_path)?;
    if manifest.benchmark_id != target.benchmark_id {
        return Err(boxed_error(format!(
            "{} declares legacy feature_id `{}`, expected benchmark ID `{}`",
            target.manifest_path.display(),
            manifest.benchmark_id,
            target.benchmark_id
        )));
    }
    if manifest.corpus_id != target.corpus_id {
        return Err(boxed_error(format!(
            "{} declares corpus_id `{}`, expected `{}`",
            target.manifest_path.display(),
            manifest.corpus_id,
            target.corpus_id
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
                target.manifest_path.display()
            )));
        }
        manifest.fixtures.retain(|candidate| candidate == fixture);
    }
    progress.manifest(&manifest.reference_tool, &manifest.reference_version);
    check_comparison_mode(&target.manifest_path, &manifest)?;
    if manifest.fixtures.is_empty() {
        return Err(boxed_error(format!(
            "{} must list at least one benchmark fixture",
            target.manifest_path.display()
        )));
    }
    check_manifest_paths(&target.manifest_path, &manifest)?;
    if accept_goldens {
        accept_implementation_goldens(&target.manifest_path, &manifest, jobs)?;
        println!(
            "  accepted {} reviewed implementation golden(s)",
            manifest.fixtures.len()
        );
    }

    let worker_count = benchmark_worker_count(jobs, manifest.fixtures.len());
    let fixture_progress = FixtureProgress::start(manifest.fixtures.len(), worker_count);
    let comparison = compare_golden_outputs(
        &target.manifest_path,
        &manifest,
        jobs,
        Some(&fixture_progress),
    );
    fixture_progress.finish();
    comparison
}
