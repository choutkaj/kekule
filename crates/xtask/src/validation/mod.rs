use crate::*;

mod compare;
mod evidence;
mod implementation;
mod manifest;
mod status;

pub(crate) use compare::*;
pub(crate) use evidence::*;
pub(crate) use implementation::*;
pub(crate) use manifest::*;
pub(crate) use status::*;

pub(crate) fn validate(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    validate_args(&args)?;
    let feature_selector = value_after_flag(&args, "--feature")
        .ok_or_else(|| boxed_error("missing required flag: --feature FEATURE_ID"))?;
    let corpus_selector = value_after_flag(&args, "--corpus").unwrap_or("smoke");
    let update = args.iter().any(|arg| arg == "--update");
    let features = read_features()?;
    if feature_selector != "all"
        && !features
            .iter()
            .any(|candidate| candidate.id == feature_selector)
    {
        return Err(boxed_error(format!("unknown feature: {feature_selector}")));
    }
    if corpus_selector != "all" && !is_known_corpus(corpus_selector) {
        return Err(boxed_error(format!("unknown corpus: {corpus_selector}")));
    }

    let targets = validation_targets(&features, feature_selector, corpus_selector);
    if targets.is_empty() {
        println!(
            "no applicable validation targets for feature `{feature_selector}` and corpus `{corpus_selector}`"
        );
        return Ok(());
    }
    let jobs = validation_jobs(&args)?;
    println!("validation worker count: {jobs}");

    let mut statuses = read_validation_statuses(&features)?;
    let mut failures = Vec::new();
    let mut passed = 0;
    let mut update_corpora = BTreeSet::new();
    for (feature, corpus) in targets {
        println!("validating `{}` against `{corpus}`", feature.id);
        if update {
            update_corpora.insert(corpus.clone());
        }
        let manifest_path = validation_manifest_path(&feature.id, &corpus);
        if !manifest_path.exists() {
            failures.push(format!(
                "{} is missing required manifest {}",
                feature.id,
                manifest_path.display()
            ));
            continue;
        }
        if update {
            statuses
                .entry(feature.id.clone())
                .or_insert_with(|| ValidationStatus::new(&feature.id))
                .corpora
                .remove(&corpus);
        }
        let result = (|| -> Result<ValidationOutcome, Box<dyn Error>> {
            let manifest = read_validation_manifest(&manifest_path)?;
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
            println!(
                "validation manifest uses {} {}",
                manifest.reference_tool, manifest.reference_version
            );
            validate_comparison_mode(&manifest_path, &manifest)?;
            if manifest.fixtures.is_empty() {
                return Err(boxed_error(format!(
                    "{} must list at least one fixture for required validation",
                    manifest_path.display()
                )));
            }
            validate_manifest_paths(&manifest_path, &manifest)?;
            println!(
                "validation manifest lists {} fixture(s)",
                manifest.fixtures.len()
            );
            let comparison = validate_golden_outputs(&manifest_path, &manifest, jobs)?;
            if comparison.compared_count > 0 {
                println!(
                    "validation compared {} golden file(s)",
                    comparison.compared_count
                );
            }
            if comparison.failed_count > 0 {
                println!(
                    "validation found {} non-passing fixture(s)",
                    comparison.failed_count
                );
                return Ok(ValidationOutcome::Failed(FailedValidationRun {
                    fixture_count: manifest.fixtures.len(),
                    compared_count: comparison.compared_count,
                    failed_count: comparison.failed_count,
                    first_failure: comparison
                        .first_failure
                        .unwrap_or_else(|| "fixture comparison failed".to_owned()),
                    reference_tool: manifest.reference_tool,
                    reference_version: manifest.reference_version,
                    manifest_hash: hash_file(&manifest_path)?,
                }));
            }
            if comparison.compared_count != manifest.fixtures.len() {
                return Err(boxed_error(format!(
                    "{} compared {} fixture(s), expected {}",
                    manifest_path.display(),
                    comparison.compared_count,
                    manifest.fixtures.len()
                )));
            }
            let evidence = build_validation_evidence(Path::new("."), &manifest_path, &manifest)?;
            Ok(ValidationOutcome::Passed(ValidationRun {
                fixture_count: manifest.fixtures.len(),
                compared_count: comparison.compared_count,
                reference_tool: manifest.reference_tool,
                reference_version: manifest.reference_version,
                manifest_hash: hash_file(&manifest_path)?,
                evidence,
            }))
        })();

        match result {
            Ok(ValidationOutcome::Passed(run)) => {
                passed += 1;
                if update {
                    let existing = statuses
                        .get(&feature.id)
                        .and_then(|status| status.corpora.get(&corpus));
                    let updated = CorpusStatus::from_run(run, existing)?;
                    statuses
                        .entry(feature.id.clone())
                        .or_insert_with(|| ValidationStatus::new(&feature.id))
                        .corpora
                        .insert(corpus, updated);
                }
            }
            Ok(ValidationOutcome::Failed(run)) => {
                let failure = format!(
                    "{} [{corpus}]: {} non-passing fixture(s); first failure: {}",
                    feature.id, run.failed_count, run.first_failure
                );
                if update {
                    let updated = CorpusStatus::from_failed_run(run)?;
                    statuses
                        .entry(feature.id.clone())
                        .or_insert_with(|| ValidationStatus::new(&feature.id))
                        .corpora
                        .insert(corpus, updated);
                }
                failures.push(failure);
            }
            Err(error) => failures.push(format!("{} [{corpus}]: {error}", feature.id)),
        }
    }

    if update {
        write_validation_statuses(&statuses, &update_corpora)?;
        sync_feature_validation_flags(&features, &statuses)?;
        let refreshed_features = read_features()?;
        let corpus_info = read_dashboard_corpus_info()?;
        let rendered = render_dashboard(&refreshed_features, &statuses, &corpus_info);
        write_atomic_text(Path::new(DASHBOARD_PATH), &rendered)?;
        println!("updated validation status and dashboard");
    }

    println!("validation passed {passed} target(s)");
    if !failures.is_empty() {
        for failure in &failures {
            eprintln!("validation failure: {failure}");
        }
        return Err(boxed_error(format!(
            "{} validation target(s) failed",
            failures.len()
        )));
    }
    Ok(())
}
