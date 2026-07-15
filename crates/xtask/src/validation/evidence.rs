use crate::*;

pub(crate) fn hash_file(path: &Path) -> Result<String, Box<dyn Error>> {
    let mut hasher = Sha256::new();
    hasher.update(fs::read(path)?);
    Ok(format!("{:x}", hasher.finalize()))
}

pub(crate) fn build_validation_evidence(
    repo_root: &Path,
    manifest_path: &Path,
    manifest: &ValidationManifest,
) -> Result<ValidationEvidence, Box<dyn Error>> {
    let corpus_root = manifest_path
        .parent()
        .and_then(Path::parent)
        .ok_or_else(|| boxed_error(format!("{} has no corpus root", manifest_path.display())))?;
    let mut paths = BTreeSet::<PathBuf>::new();
    paths.insert(manifest_path.to_path_buf());
    paths.insert(corpus_root.join("corpus.toml"));
    paths.insert(corpus_root.join("sources.lock.json"));
    paths.insert(
        repo_root
            .join("features")
            .join(&manifest.feature_id)
            .join("feature.toml"),
    );
    paths.insert(
        repo_root
            .join("features")
            .join(&manifest.feature_id)
            .join("feature.md"),
    );
    paths.insert(repo_root.join("Cargo.toml"));
    paths.insert(repo_root.join("Cargo.lock"));
    paths.insert(repo_root.join("crates/molecules/Cargo.toml"));
    paths.insert(repo_root.join("crates/xtask/Cargo.toml"));

    collect_files(&repo_root.join("crates/molecules/src"), &mut paths)?;
    collect_files(&repo_root.join("crates/xtask/src"), &mut paths)?;

    match manifest.reference_tool.as_str() {
        "rdkit" => {
            let reference_root = repo_root.join("validation/reference/rdkit");
            paths.insert(reference_root.join("run_feature.py"));
            paths.insert(reference_root.join("environment.yml"));
        }
        "biopython" => {
            let reference_root = repo_root.join("validation/reference/biopython");
            paths.insert(reference_root.join("run_feature.py"));
            paths.insert(reference_root.join("environment.yml"));
        }
        value if is_manual_semantic_reference_tool(value) => {}
        value => {
            return Err(boxed_error(format!(
                "{} uses unsupported reference_tool `{value}`",
                manifest_path.display()
            )))
        }
    }

    for fixture in &manifest.fixtures {
        paths.insert(corpus_root.join(fixture));
        paths.insert(
            corpus_root
                .join("golden")
                .join(&manifest.feature_id)
                .join(format!("{}.json.gz", slugify_fixture(fixture))),
        );
    }

    let mut inputs = Vec::new();
    for path in paths {
        if !path.exists() {
            return Err(boxed_error(format!(
                "validation evidence input is missing: {}",
                path.display()
            )));
        }
        if !path.is_file() {
            continue;
        }
        inputs.push(EvidenceInput {
            path: relative_path(repo_root, &path)?,
            sha256: hash_evidence_file(&path)?,
        });
    }
    inputs.sort_by(|left, right| left.path.cmp(&right.path));
    let evidence_document = json!({
        "schema_version": VALIDATION_EVIDENCE_SCHEMA_VERSION,
        "comparison_mode": manifest.comparison_mode,
        "inputs": inputs,
    });
    let mut hasher = Sha256::new();
    hasher.update(serde_json::to_vec(&evidence_document)?);
    let sha256 = format!("{:x}", hasher.finalize());
    let inputs = serde_json::from_value(
        evidence_document
            .get("inputs")
            .cloned()
            .ok_or_else(|| boxed_error("evidence document has no inputs"))?,
    )?;
    Ok(ValidationEvidence {
        schema_version: VALIDATION_EVIDENCE_SCHEMA_VERSION,
        comparison_mode: manifest.comparison_mode.clone(),
        inputs,
        sha256,
    })
}

pub(crate) fn collect_files(
    root: &Path,
    paths: &mut BTreeSet<PathBuf>,
) -> Result<(), Box<dyn Error>> {
    if !root.exists() {
        return Err(boxed_error(format!("{} does not exist", root.display())));
    }
    for entry in fs::read_dir(root)? {
        let path = entry?.path();
        if path.is_dir() {
            collect_files(&path, paths)?;
        } else {
            paths.insert(path);
        }
    }
    Ok(())
}

pub(crate) fn hash_evidence_file(path: &Path) -> Result<String, Box<dyn Error>> {
    let raw = fs::read(path)?;
    let normalized_text = String::from_utf8(raw.clone())
        .ok()
        .map(|text| text.replace("\r\n", "\n").replace('\r', "\n"));
    let bytes = if path.file_name().and_then(|name| name.to_str()) == Some("feature.toml")
        && path
            .components()
            .any(|component| component.as_os_str() == "features")
    {
        normalized_text
            .ok_or_else(|| boxed_error(format!("{} is not UTF-8", path.display())))?
            .lines()
            .map(|line| {
                if line.trim_start().starts_with("validated =") {
                    "validated = <generated>"
                } else {
                    line
                }
            })
            .collect::<Vec<_>>()
            .join("\n")
            .into_bytes()
    } else if let Some(text) = normalized_text {
        text.into_bytes()
    } else {
        raw
    };
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    Ok(format!("{:x}", hasher.finalize()))
}

pub(crate) fn relative_path(repo_root: &Path, path: &Path) -> Result<String, Box<dyn Error>> {
    let relative = if path.is_absolute() {
        path.strip_prefix(repo_root).unwrap_or(path)
    } else {
        path
    };
    relative
        .to_str()
        .map(|value| value.replace('\\', "/").trim_start_matches("./").to_owned())
        .ok_or_else(|| boxed_error(format!("{} is not valid UTF-8", path.display())))
}

pub(crate) fn read_gzip_string(path: &Path) -> Result<String, Box<dyn Error>> {
    let file = fs::File::open(path)?;
    let mut decoder = GzDecoder::new(file);
    let mut text = String::new();
    decoder.read_to_string(&mut text)?;
    Ok(text)
}

pub(crate) fn accept_implementation_goldens(
    manifest_path: &Path,
    manifest: &ValidationManifest,
    jobs: usize,
) -> Result<(), Box<dyn Error>> {
    if !is_manual_semantic_reference_tool(&manifest.reference_tool) {
        return Err(boxed_error(format!(
            "{} uses generator-backed reference tool `{}`; only *-manual-semantic implementation goldens can be accepted from the Rust implementation",
            manifest_path.display(),
            manifest.reference_tool
        )));
    }
    let corpus_root = manifest_path
        .parent()
        .and_then(Path::parent)
        .ok_or_else(|| boxed_error(format!("{} has no corpus root", manifest_path.display())))?;
    let worker_count = validation_worker_count(jobs, manifest.fixtures.len());
    if worker_count == 1 {
        for fixture in &manifest.fixtures {
            accept_one_implementation_golden(corpus_root, manifest, fixture)?;
        }
        return Ok(());
    }

    let next_fixture = std::sync::Mutex::new(0usize);
    let results = std::sync::Mutex::new(vec![None; manifest.fixtures.len()]);
    std::thread::scope(|scope| {
        for _ in 0..worker_count {
            scope.spawn(|| loop {
                let index = {
                    let mut next = next_fixture
                        .lock()
                        .expect("implementation golden queue lock should not be poisoned");
                    if *next >= manifest.fixtures.len() {
                        break;
                    }
                    let index = *next;
                    *next += 1;
                    index
                };
                let result = accept_one_implementation_golden(
                    corpus_root,
                    manifest,
                    &manifest.fixtures[index],
                )
                .map_err(|error| error.to_string());
                results
                    .lock()
                    .expect("implementation golden result lock should not be poisoned")[index] =
                    Some(result);
            });
        }
    });
    for result in results
        .into_inner()
        .expect("implementation golden result lock should not be poisoned")
    {
        result
            .ok_or_else(|| boxed_error("implementation golden worker recorded no result"))?
            .map_err(boxed_error)?;
    }
    Ok(())
}

fn accept_one_implementation_golden(
    corpus_root: &Path,
    manifest: &ValidationManifest,
    fixture: &str,
) -> Result<(), Box<dyn Error>> {
    let fixture_path = corpus_root.join(fixture);
    let expected =
        implementation_expected(&manifest.feature_id, &manifest.corpus_id, &fixture_path)?;
    let document = json!({
        "schema_version": GOLDEN_SCHEMA_VERSION,
        "feature_id": manifest.feature_id,
        "corpus_id": manifest.corpus_id,
        "fixture_id": slugify_fixture(fixture),
        "fixture_path": fixture,
        "input_sha256": hash_file(&fixture_path)?,
        "reference": {
            "tool": manifest.reference_tool,
            "version": manifest.reference_version,
            "runtime_dependency": false,
        },
        "expected": expected,
    });
    let mut encoder = GzBuilder::new()
        .mtime(0)
        .write(Vec::new(), Compression::default());
    serde_json::to_writer_pretty(&mut encoder, &document)?;
    encoder.write_all(b"\n")?;
    let compressed = encoder.finish()?;
    let golden_path = corpus_root
        .join("golden")
        .join(&manifest.feature_id)
        .join(format!("{}.json.gz", slugify_fixture(fixture)));
    write_atomic_bytes(&golden_path, &compressed)
}

pub(crate) fn validate_manifest_paths(
    manifest_path: &Path,
    manifest: &ValidationManifest,
) -> Result<(), Box<dyn Error>> {
    let base = manifest_path
        .parent()
        .and_then(Path::parent)
        .ok_or_else(|| {
            boxed_error(format!(
                "{} has no parent directory",
                manifest_path.display()
            ))
        })?;
    for fixture in &manifest.fixtures {
        let path = base.join(fixture);
        if !path.exists() {
            return Err(boxed_error(format!(
                "{} references missing fixture `{fixture}`",
                manifest_path.display()
            )));
        }
    }
    let lock = read_source_lock(&manifest.corpus_id)?;
    let pinned_paths = lock
        .entries
        .iter()
        .flat_map(|entry| entry.files.iter().map(|file| file.path.as_str()))
        .chain(lock.packs.iter().map(|pack| pack.path.as_str()))
        .collect::<BTreeSet<_>>();
    for fixture in &manifest.fixtures {
        if !pinned_paths.contains(fixture.as_str()) {
            return Err(boxed_error(format!(
                "{} fixture `{fixture}` is not pinned by sources.lock.json",
                manifest_path.display()
            )));
        }
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct ValidationComparison {
    pub(crate) compared_count: usize,
    pub(crate) failed_count: usize,
    pub(crate) first_failure: Option<String>,
}

pub(crate) fn validate_golden_outputs(
    manifest_path: &Path,
    manifest: &ValidationManifest,
    jobs: usize,
    progress: Option<&FixtureProgress>,
) -> Result<ValidationComparison, Box<dyn Error>> {
    if manifest.fixtures.is_empty() {
        return Ok(ValidationComparison {
            compared_count: 0,
            failed_count: 0,
            first_failure: None,
        });
    }
    let base = manifest_path
        .parent()
        .and_then(Path::parent)
        .ok_or_else(|| {
            boxed_error(format!(
                "{} has no parent directory",
                manifest_path.display()
            ))
        })?;
    let worker_count = validation_worker_count(jobs, manifest.fixtures.len());
    if worker_count == 1 {
        return validate_golden_outputs_serial(manifest_path, manifest, base, progress);
    }

    let next_fixture = std::sync::Mutex::new(0usize);
    let results = std::sync::Mutex::new(
        (0..manifest.fixtures.len())
            .map(|_| None)
            .collect::<Vec<Option<Result<FixtureComparison, String>>>>(),
    );
    std::thread::scope(|scope| {
        for _ in 0..worker_count {
            scope.spawn(|| loop {
                let index = {
                    let mut next = next_fixture
                        .lock()
                        .expect("validation fixture queue lock should not be poisoned");
                    if *next >= manifest.fixtures.len() {
                        None
                    } else {
                        let index = *next;
                        *next += 1;
                        Some(index)
                    }
                };
                let Some(index) = index else {
                    break;
                };
                let fixture = &manifest.fixtures[index];
                let result = compare_one_golden(manifest_path, base, manifest, fixture)
                    .map_err(|error| error.to_string());
                if let Some(progress) = progress {
                    progress.fixture_finished();
                }
                results
                    .lock()
                    .expect("validation result lock should not be poisoned")[index] = Some(result);
            });
        }
    });

    let results = results
        .into_inner()
        .expect("validation result lock should not be poisoned");
    let mut comparison = ValidationComparison {
        compared_count: 0,
        failed_count: 0,
        first_failure: None,
    };
    for result in results {
        let result = result
            .ok_or_else(|| boxed_error("validation worker did not record a fixture result"))?;
        record_fixture_comparison(&mut comparison, result.map_err(boxed_error)?);
    }
    Ok(comparison)
}

fn validate_golden_outputs_serial(
    manifest_path: &Path,
    manifest: &ValidationManifest,
    base: &Path,
    progress: Option<&FixtureProgress>,
) -> Result<ValidationComparison, Box<dyn Error>> {
    let mut comparison = ValidationComparison {
        compared_count: 0,
        failed_count: 0,
        first_failure: None,
    };
    for fixture in &manifest.fixtures {
        let result = compare_one_golden(manifest_path, base, manifest, fixture);
        if let Some(progress) = progress {
            progress.fixture_finished();
        }
        let result = result?;
        record_fixture_comparison(&mut comparison, result);
    }
    Ok(comparison)
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum FixtureComparison {
    Passed,
    Failed(String),
}

fn record_fixture_comparison(
    comparison: &mut ValidationComparison,
    fixture_result: FixtureComparison,
) {
    match fixture_result {
        FixtureComparison::Passed => comparison.compared_count += 1,
        FixtureComparison::Failed(failure) => {
            eprintln!("fixture comparison failure: {failure}");
            comparison.failed_count += 1;
            comparison.first_failure.get_or_insert(failure);
        }
    }
}

fn compare_one_golden(
    manifest_path: &Path,
    base: &Path,
    manifest: &ValidationManifest,
    fixture: &str,
) -> Result<FixtureComparison, Box<dyn Error>> {
    let fixture_path = base.join(fixture);
    let golden_path = base
        .join("golden")
        .join(&manifest.feature_id)
        .join(format!("{}.json.gz", slugify_fixture(fixture)));
    if !golden_path.exists() {
        return Err(boxed_error(format!(
            "{} is missing golden file for fixture `{fixture}`",
            manifest_path.display()
        )));
    }
    let golden: Value = serde_json::from_str(&read_gzip_string(&golden_path)?)?;
    validate_golden_metadata(&golden_path, &golden, manifest, fixture, &fixture_path)?;
    let expected = golden
        .get("expected")
        .ok_or_else(|| boxed_error(format!("{} is missing `expected`", golden_path.display())))?;
    let actual =
        match implementation_expected(&manifest.feature_id, &manifest.corpus_id, &fixture_path) {
            Ok(actual) => actual,
            Err(error) => {
                return Ok(FixtureComparison::Failed(format!(
                    "fixture `{fixture}` implementation output failed: {error}"
                )))
            }
        };
    let expected = normalize_for_comparison(expected);
    let actual = normalize_for_comparison(&actual);
    if let Some(diff) = first_json_diff(&manifest.feature_id, "$", &expected, &actual) {
        return Ok(FixtureComparison::Failed(format!(
            "{} differs from implementation output for fixture `{fixture}`: {diff}",
            golden_path.display()
        )));
    }
    Ok(FixtureComparison::Passed)
}

pub(crate) fn validate_golden_metadata(
    golden_path: &Path,
    golden: &Value,
    manifest: &ValidationManifest,
    fixture: &str,
    fixture_path: &Path,
) -> Result<(), Box<dyn Error>> {
    if golden.get("schema_version") != Some(&json!(GOLDEN_SCHEMA_VERSION)) {
        return Err(boxed_error(format!(
            "{} has unsupported schema_version",
            golden_path.display()
        )));
    }
    if golden.get("feature_id").and_then(Value::as_str) != Some(manifest.feature_id.as_str()) {
        return Err(boxed_error(format!(
            "{} feature_id does not match manifest",
            golden_path.display()
        )));
    }
    if golden.get("corpus_id").and_then(Value::as_str) != Some(manifest.corpus_id.as_str()) {
        return Err(boxed_error(format!(
            "{} corpus_id does not match manifest",
            golden_path.display()
        )));
    }
    if golden.get("fixture_path").and_then(Value::as_str) != Some(fixture) {
        return Err(boxed_error(format!(
            "{} fixture_path does not match manifest",
            golden_path.display()
        )));
    }
    let input_sha256 = golden
        .get("input_sha256")
        .and_then(Value::as_str)
        .ok_or_else(|| boxed_error(format!("{} is missing input_sha256", golden_path.display())))?;
    let fixture_hash = hash_file(fixture_path)?;
    if input_sha256 != fixture_hash {
        return Err(boxed_error(format!(
            "{} input_sha256 does not match current fixture `{fixture}`",
            golden_path.display()
        )));
    }
    let reference = golden
        .get("reference")
        .and_then(Value::as_object)
        .ok_or_else(|| boxed_error(format!("{} is missing reference", golden_path.display())))?;
    if reference.get("tool").and_then(Value::as_str) != Some(manifest.reference_tool.as_str()) {
        return Err(boxed_error(format!(
            "{} reference.tool does not match manifest",
            golden_path.display()
        )));
    }
    let golden_version = reference
        .get("version")
        .and_then(Value::as_str)
        .ok_or_else(|| {
            boxed_error(format!(
                "{} reference.version is missing",
                golden_path.display()
            ))
        })?;
    if reference_version_label(&manifest.reference_tool, golden_version)
        != manifest.reference_version
    {
        return Err(boxed_error(format!(
            "{} reference.version does not match manifest",
            golden_path.display()
        )));
    }
    if reference.get("runtime_dependency").and_then(Value::as_bool) != Some(false) {
        return Err(boxed_error(format!(
            "{} must record reference.runtime_dependency=false",
            golden_path.display()
        )));
    }
    Ok(())
}

pub(crate) fn is_manual_semantic_reference_tool(tool: &str) -> bool {
    tool.ends_with("-manual-semantic")
}

pub(crate) fn reference_version_label(tool: &str, version: &str) -> String {
    match tool {
        "rdkit" if !version.starts_with("RDKit ") => format!("RDKit {version}"),
        "biopython" if !version.starts_with("Biopython ") => format!("Biopython {version}"),
        _ => version.to_owned(),
    }
}
