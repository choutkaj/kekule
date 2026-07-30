use crate::*;

#[derive(Default)]
pub(crate) struct BenchmarkHashCache {
    exact: std::sync::Mutex<BTreeMap<PathBuf, String>>,
    normalized: std::sync::Mutex<BTreeMap<PathBuf, String>>,
}

impl BenchmarkHashCache {
    pub(crate) fn exact_hash(&self, path: &Path) -> Result<String, Box<dyn Error>> {
        if let Some(hash) = self
            .exact
            .lock()
            .map_err(|_| boxed_error("benchmark exact-hash cache lock was poisoned"))?
            .get(path)
            .cloned()
        {
            return Ok(hash);
        }
        let hash = hash_file(path)?;
        let mut cache = self
            .exact
            .lock()
            .map_err(|_| boxed_error("benchmark exact-hash cache lock was poisoned"))?;
        Ok(cache.entry(path.to_path_buf()).or_insert(hash).clone())
    }

    pub(crate) fn normalized_hash(&self, path: &Path) -> Result<String, Box<dyn Error>> {
        if let Some(hash) = self
            .normalized
            .lock()
            .map_err(|_| boxed_error("benchmark normalized-hash cache lock was poisoned"))?
            .get(path)
            .cloned()
        {
            return Ok(hash);
        }
        let hash = hash_normalized_file(path)?;
        let mut cache = self
            .normalized
            .lock()
            .map_err(|_| boxed_error("benchmark normalized-hash cache lock was poisoned"))?;
        Ok(cache.entry(path.to_path_buf()).or_insert(hash).clone())
    }
}

fn digest_hex(digest: impl AsRef<[u8]>) -> String {
    digest
        .as_ref()
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect()
}

pub(crate) fn hash_file(path: &Path) -> Result<String, Box<dyn Error>> {
    let file = fs::File::open(path)?;
    let mut reader = std::io::BufReader::new(file);
    let mut buffer = vec![0u8; 64 * 1024];
    let mut hasher = Sha256::new();
    loop {
        let read = reader.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
    }
    Ok(digest_hex(hasher.finalize()))
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
struct DigestInput {
    path: String,
    sha256: String,
}

#[cfg(test)]
pub(crate) fn build_benchmark_input_digest(
    repo_root: &Path,
    manifest_path: &Path,
    manifest: &BenchmarkManifest,
) -> Result<BenchmarkInputDigest, Box<dyn Error>> {
    build_benchmark_input_digest_cached(
        repo_root,
        manifest_path,
        manifest,
        &BenchmarkHashCache::default(),
    )
}

#[cfg(test)]
pub(crate) fn build_benchmark_input_digest_with_core_source_root(
    repo_root: &Path,
    manifest_path: &Path,
    manifest: &BenchmarkManifest,
    core_source_root: &Path,
) -> Result<BenchmarkInputDigest, Box<dyn Error>> {
    build_benchmark_input_digest_cached_with_core_source_root(
        repo_root,
        manifest_path,
        manifest,
        &BenchmarkHashCache::default(),
        core_source_root,
    )
}

pub(crate) fn build_benchmark_input_digest_cached(
    repo_root: &Path,
    manifest_path: &Path,
    manifest: &BenchmarkManifest,
    hash_cache: &BenchmarkHashCache,
) -> Result<BenchmarkInputDigest, Box<dyn Error>> {
    build_benchmark_input_digest_cached_with_core_source_root(
        repo_root,
        manifest_path,
        manifest,
        hash_cache,
        &repo_root.join("crates/kekule/src"),
    )
}

fn build_benchmark_input_digest_cached_with_core_source_root(
    repo_root: &Path,
    manifest_path: &Path,
    manifest: &BenchmarkManifest,
    hash_cache: &BenchmarkHashCache,
    core_source_root: &Path,
) -> Result<BenchmarkInputDigest, Box<dyn Error>> {
    let corpus_root = manifest_path
        .parent()
        .and_then(Path::parent)
        .ok_or_else(|| boxed_error(format!("{} has no corpus root", manifest_path.display())))?;
    let mut paths = BTreeSet::<PathBuf>::new();
    paths.insert(manifest_path.to_path_buf());
    paths.insert(corpus_root.join("corpus.toml"));
    paths.insert(corpus_root.join("sources.lock.json"));
    collect_files(core_source_root, &mut paths)?;
    collect_files(&repo_root.join("crates/xtask/src"), &mut paths)?;

    match manifest.reference_tool.as_str() {
        "rdkit" => {
            let reference_root = repo_root.join("benchmarks/reference/rdkit");
            paths.insert(reference_root.join("run_feature.py"));
            paths.insert(reference_root.join("environment.yml"));
        }
        "biopython" => {
            let reference_root = repo_root.join("benchmarks/reference/biopython");
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
                "benchmark input is missing: {}",
                path.display()
            )));
        }
        if !path.is_file() {
            continue;
        }
        inputs.push(DigestInput {
            path: benchmark_digest_input_path(repo_root, core_source_root, &path)?,
            sha256: hash_cache.normalized_hash(&path)?,
        });
    }
    inputs.push(DigestInput {
        path: "Cargo.lock#external-dependencies".to_owned(),
        sha256: external_dependency_hash(&repo_root.join("Cargo.lock"))?,
    });
    inputs.sort_by(|left, right| left.path.cmp(&right.path));
    let input_count = inputs.len();
    let digest_document = json!({
        "schema_version": BENCHMARK_INPUT_DIGEST_SCHEMA_VERSION,
        "comparison_mode": manifest.comparison_mode,
        "inputs": inputs,
    });
    let mut hasher = Sha256::new();
    hasher.update(serde_json::to_vec(&digest_document)?);
    let sha256 = digest_hex(hasher.finalize());
    Ok(BenchmarkInputDigest {
        schema_version: BENCHMARK_INPUT_DIGEST_SCHEMA_VERSION,
        input_count,
        sha256,
    })
}

fn benchmark_digest_input_path(
    repo_root: &Path,
    core_source_root: &Path,
    path: &Path,
) -> Result<String, Box<dyn Error>> {
    if path.strip_prefix(core_source_root).is_ok() {
        return Ok(format!(
            "implementation/core/src/{}",
            relative_path(core_source_root, path)?
        ));
    }
    relative_path(repo_root, path)
}

fn external_dependency_hash(path: &Path) -> Result<String, Box<dyn Error>> {
    let lock: toml::Value = toml::from_str(&fs::read_to_string(path)?)?;
    let mut packages = lock
        .get("package")
        .and_then(toml::Value::as_array)
        .ok_or_else(|| boxed_error(format!("{} has no package array", path.display())))?
        .iter()
        .filter(|package| package.get("source").is_some())
        .map(|package| {
            let mut dependencies = package
                .get("dependencies")
                .and_then(toml::Value::as_array)
                .map(|values| {
                    values
                        .iter()
                        .filter_map(toml::Value::as_str)
                        .map(str::to_owned)
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default();
            dependencies.sort();
            json!({
                "name": package.get("name").and_then(toml::Value::as_str),
                "version": package.get("version").and_then(toml::Value::as_str),
                "source": package.get("source").and_then(toml::Value::as_str),
                "checksum": package.get("checksum").and_then(toml::Value::as_str),
                "dependencies": dependencies,
            })
        })
        .collect::<Vec<_>>();
    packages.sort_by_key(|package| package.to_string());
    let mut hasher = Sha256::new();
    hasher.update(serde_json::to_vec(&packages)?);
    Ok(digest_hex(hasher.finalize()))
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

pub(crate) fn hash_normalized_file(path: &Path) -> Result<String, Box<dyn Error>> {
    let mut hasher = Sha256::new();
    match String::from_utf8(fs::read(path)?) {
        Ok(text) if text.contains('\r') => {
            hasher.update(text.replace("\r\n", "\n").replace('\r', "\n"));
        }
        Ok(text) => hasher.update(text),
        Err(error) => hasher.update(error.into_bytes()),
    }
    Ok(digest_hex(hasher.finalize()))
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
    manifest: &BenchmarkManifest,
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
    let worker_count = benchmark_worker_count(jobs, manifest.fixtures.len());
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
    manifest: &BenchmarkManifest,
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

pub(crate) fn check_manifest_paths(
    manifest_path: &Path,
    manifest: &BenchmarkManifest,
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
pub(crate) struct BenchmarkComparison {
    pub(crate) compared_count: usize,
    pub(crate) difference_count: usize,
    pub(crate) first_difference: Option<String>,
}

#[cfg(test)]
pub(crate) fn compare_golden_outputs(
    manifest_path: &Path,
    manifest: &BenchmarkManifest,
    jobs: usize,
    progress: Option<&FixtureProgress>,
) -> Result<BenchmarkComparison, Box<dyn Error>> {
    compare_golden_outputs_cached(
        manifest_path,
        manifest,
        jobs,
        progress,
        &BenchmarkHashCache::default(),
    )
}

pub(crate) fn compare_golden_outputs_cached(
    manifest_path: &Path,
    manifest: &BenchmarkManifest,
    jobs: usize,
    progress: Option<&FixtureProgress>,
    hash_cache: &BenchmarkHashCache,
) -> Result<BenchmarkComparison, Box<dyn Error>> {
    if manifest.fixtures.is_empty() {
        return Ok(BenchmarkComparison {
            compared_count: 0,
            difference_count: 0,
            first_difference: None,
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
    let worker_count = benchmark_worker_count(jobs, manifest.fixtures.len());
    if worker_count == 1 {
        return compare_golden_outputs_serial(manifest_path, manifest, base, progress, hash_cache);
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
                        .expect("benchmark fixture queue lock should not be poisoned");
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
                let result = compare_one_golden(manifest_path, base, manifest, fixture, hash_cache)
                    .map_err(|error| error.to_string());
                if let Some(progress) = progress {
                    progress.fixture_finished();
                }
                results
                    .lock()
                    .expect("benchmark result lock should not be poisoned")[index] = Some(result);
            });
        }
    });

    let results = results
        .into_inner()
        .expect("benchmark result lock should not be poisoned");
    let mut comparison = BenchmarkComparison {
        compared_count: 0,
        difference_count: 0,
        first_difference: None,
    };
    for result in results {
        let result = result
            .ok_or_else(|| boxed_error("benchmark worker did not record a fixture result"))?;
        record_fixture_comparison(&mut comparison, result.map_err(boxed_error)?);
    }
    Ok(comparison)
}

fn compare_golden_outputs_serial(
    manifest_path: &Path,
    manifest: &BenchmarkManifest,
    base: &Path,
    progress: Option<&FixtureProgress>,
    hash_cache: &BenchmarkHashCache,
) -> Result<BenchmarkComparison, Box<dyn Error>> {
    let mut comparison = BenchmarkComparison {
        compared_count: 0,
        difference_count: 0,
        first_difference: None,
    };
    for fixture in &manifest.fixtures {
        let result = compare_one_golden(manifest_path, base, manifest, fixture, hash_cache);
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
    Match,
    Difference(String),
}

fn record_fixture_comparison(
    comparison: &mut BenchmarkComparison,
    fixture_result: FixtureComparison,
) {
    match fixture_result {
        FixtureComparison::Match => comparison.compared_count += 1,
        FixtureComparison::Difference(difference) => {
            eprintln!("fixture comparison difference: {difference}");
            comparison.difference_count += 1;
            comparison.first_difference.get_or_insert(difference);
        }
    }
}

fn compare_one_golden(
    manifest_path: &Path,
    base: &Path,
    manifest: &BenchmarkManifest,
    fixture: &str,
    hash_cache: &BenchmarkHashCache,
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
    let mut golden: Value = serde_json::from_str(&read_gzip_string(&golden_path)?)?;
    check_golden_metadata(
        &golden_path,
        &golden,
        manifest,
        fixture,
        &fixture_path,
        hash_cache,
    )?;
    let expected = golden
        .get_mut("expected")
        .ok_or_else(|| boxed_error(format!("{} is missing `expected`", golden_path.display())))?;
    let mut actual =
        match implementation_expected(&manifest.feature_id, &manifest.corpus_id, &fixture_path) {
            Ok(actual) => actual,
            Err(error) => {
                return Ok(FixtureComparison::Difference(format!(
                    "fixture `{fixture}` implementation output failed: {error}"
                )))
            }
        };
    normalize_feature_for_comparison_in_place(&manifest.feature_id, expected);
    normalize_feature_for_comparison_in_place(&manifest.feature_id, &mut actual);
    if let Some(diff) = first_json_diff(&manifest.feature_id, "$", expected, &actual) {
        return Ok(FixtureComparison::Difference(format!(
            "{} differs from implementation output for fixture `{fixture}`: {diff}",
            golden_path.display()
        )));
    }
    Ok(FixtureComparison::Match)
}

pub(crate) fn check_golden_metadata(
    golden_path: &Path,
    golden: &Value,
    manifest: &BenchmarkManifest,
    fixture: &str,
    fixture_path: &Path,
    hash_cache: &BenchmarkHashCache,
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
    let fixture_hash = hash_cache.exact_hash(fixture_path)?;
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
