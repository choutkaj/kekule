use crate::*;

pub(crate) fn corpus(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    if args.first().map(String::as_str) != Some("check") {
        return Err(boxed_error(
            "usage: cargo xtask corpus check --corpus CORPUS_ID|all [--require-data]",
        ));
    }
    let args = &args[1..];
    let selector = value_after_flag(args, "--corpus")
        .ok_or_else(|| boxed_error("missing required flag: --corpus CORPUS_ID|all"))?;
    let require_data = args.iter().any(|arg| arg == "--require-data");
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--corpus" => index += 2,
            "--require-data" => index += 1,
            arg => return Err(boxed_error(format!("unknown corpus check argument: {arg}"))),
        }
    }
    let available_corpora = discover_corpus_ids()?;
    if selector != "all" && !available_corpora.iter().any(|corpus| corpus == selector) {
        return Err(boxed_error(format!("unknown corpus: {selector}")));
    }

    let corpora = available_corpora
        .iter()
        .map(String::as_str)
        .filter(|id| selector == "all" || selector == *id)
        .collect::<Vec<_>>();
    let mut locks = BTreeMap::new();
    for corpus_id in &corpora {
        let descriptor = read_corpus_descriptor(corpus_id)?;
        if !descriptor.ready {
            println!("corpus `{corpus_id}` is declared but not built; skipping integrity checks");
            continue;
        }
        let lock = read_source_lock(corpus_id)?;
        check_corpus_lock(&descriptor, &lock)?;
        check_corpus_artifacts(
            corpus_id,
            &lock,
            corpus_requires_data(corpus_id, require_data),
            &descriptor.build_command,
        )?;
        println!(
            "corpus `{corpus_id}` has {} pinned entries; integrity checks completed",
            lock.entries.len()
        );
        locks.insert((*corpus_id).to_owned(), lock);
    }
    for corpus_id in ["pubchem-100", "pubchem-1k", "pdb-10", "pdb-100", "pdb-1000"] {
        if !locks.contains_key(corpus_id) {
            let descriptor = read_tracked_corpus_descriptor(corpus_id)?;
            let lock = read_source_lock(corpus_id)?;
            check_corpus_lock(&descriptor, &lock)?;
            locks.insert(corpus_id.to_owned(), lock);
        }
    }
    check_nested_corpora(&locks)?;
    Ok(())
}

pub(crate) fn corpus_requires_data(_corpus: &str, requested: bool) -> bool {
    requested
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub(crate) enum CorpusKind {
    SmallMolecule,
    Macromolecule,
    Mixed,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct CorpusDescriptor {
    pub(crate) id: String,
    pub(crate) title: String,
    #[serde(rename = "kind")]
    pub(crate) _kind: CorpusKind,
    pub(crate) ready: bool,
    pub(crate) expected_count: usize,
    #[serde(default)]
    #[serde(rename = "local_only")]
    pub(crate) _local_only: bool,
    #[serde(default)]
    pub(crate) parent: Option<String>,
    #[serde(default)]
    pub(crate) selection_id: Option<String>,
    #[serde(default)]
    pub(crate) formats: Vec<String>,
    #[serde(default)]
    pub(crate) categories: BTreeMap<String, usize>,
    #[serde(default, rename = "notes")]
    pub(crate) _notes: Vec<String>,
    pub(crate) build_command: String,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct SourceLock {
    pub(crate) schema_version: u32,
    pub(crate) corpus_id: String,
    pub(crate) source: String,
    pub(crate) selection_id: String,
    pub(crate) entries: Vec<SourceEntry>,
    #[serde(default)]
    pub(crate) packs: Vec<SourcePack>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct SourceEntry {
    pub(crate) id: String,
    pub(crate) category: String,
    pub(crate) files: Vec<SourceFile>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct SourceFile {
    pub(crate) path: String,
    pub(crate) url: String,
    pub(crate) sha256: String,
    #[serde(default, rename = "record_type")]
    pub(crate) _record_type: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct SourcePack {
    pub(crate) path: String,
    pub(crate) format: String,
    pub(crate) count: usize,
    pub(crate) members: Vec<String>,
    pub(crate) sha256: String,
    #[serde(default)]
    pub(crate) member_id_property: Option<String>,
    #[serde(default)]
    pub(crate) member_title_prefix: Option<String>,
}

pub(crate) fn corpus_root(corpus: &str) -> PathBuf {
    workspace_root()
        .join("benchmarks")
        .join("corpora")
        .join(corpus)
}

pub(crate) fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

pub(crate) fn discover_corpus_ids() -> Result<Vec<String>, Box<dyn Error>> {
    discover_corpus_ids_from(&workspace_root())
}

pub(crate) fn discover_corpus_ids_from(root: &Path) -> Result<Vec<String>, Box<dyn Error>> {
    let corpora_root = root.join("benchmarks").join("corpora");
    let mut ids = Vec::new();
    for entry in fs::read_dir(&corpora_root)? {
        let path = entry?.path();
        if !path.is_dir() || !path.join("corpus.toml").is_file() {
            continue;
        }
        let id = path
            .file_name()
            .and_then(|name| name.to_str())
            .ok_or_else(|| {
                boxed_error(format!("{} has a non-UTF-8 corpus name", path.display()))
            })?;
        ids.push(id.to_owned());
    }
    ids.sort();
    Ok(ids)
}

pub(crate) fn is_known_corpus(corpus: &str) -> bool {
    corpus_descriptor_path(corpus).is_file()
}

pub(crate) fn corpus_descriptor_path(corpus: &str) -> PathBuf {
    corpus_root(corpus).join("corpus.toml")
}

pub(crate) fn read_corpus_descriptor(corpus: &str) -> Result<CorpusDescriptor, Box<dyn Error>> {
    read_tracked_corpus_descriptor(corpus)
}

pub(crate) fn read_tracked_corpus_descriptor(
    corpus: &str,
) -> Result<CorpusDescriptor, Box<dyn Error>> {
    let path = corpus_descriptor_path(corpus);
    let text = fs::read_to_string(&path)?;
    let descriptor: CorpusDescriptor = toml::from_str(&text)
        .map_err(|error| boxed_error(format!("{}: {error}", path.display())))?;
    if descriptor.id != corpus {
        return Err(boxed_error(format!(
            "{} declares id `{}`, expected `{}`",
            path.display(),
            descriptor.id,
            corpus
        )));
    }
    Ok(descriptor)
}

pub(crate) fn read_source_lock(corpus: &str) -> Result<SourceLock, Box<dyn Error>> {
    let path = corpus_root(corpus).join("sources.lock.json");
    read_source_lock_path(&path)
}

pub(crate) fn read_source_lock_path(path: &Path) -> Result<SourceLock, Box<dyn Error>> {
    let text = fs::read_to_string(path)
        .map_err(|error| boxed_error(format!("{} is unavailable: {error}", path.display())))?;
    serde_json::from_str(&text).map_err(|error| boxed_error(format!("{}: {error}", path.display())))
}

pub(crate) fn check_corpus_lock(
    descriptor: &CorpusDescriptor,
    lock: &SourceLock,
) -> Result<(), Box<dyn Error>> {
    if descriptor.title.trim().is_empty()
        || descriptor.formats.is_empty()
        || lock.source.trim().is_empty()
    {
        return Err(boxed_error(format!(
            "{} has incomplete corpus metadata",
            descriptor.id
        )));
    }
    if let Some(parent) = &descriptor.parent {
        if !is_known_corpus(parent) {
            return Err(boxed_error(format!(
                "{} names unknown parent corpus `{parent}`",
                descriptor.id
            )));
        }
    }
    if lock.schema_version != 1 || lock.corpus_id != descriptor.id {
        return Err(boxed_error(format!(
            "{} has incompatible source lock metadata",
            descriptor.id
        )));
    }
    if descriptor.selection_id.as_deref() != Some(lock.selection_id.as_str()) {
        return Err(boxed_error(format!(
            "{} selection ID does not match corpus.toml",
            descriptor.id
        )));
    }
    if lock.entries.len() != descriptor.expected_count {
        return Err(boxed_error(format!(
            "{} contains {} entries, expected {}",
            descriptor.id,
            lock.entries.len(),
            descriptor.expected_count
        )));
    }
    let mut ids = BTreeSet::new();
    let mut categories = BTreeMap::<String, usize>::new();
    for entry in &lock.entries {
        if !ids.insert(entry.id.as_str()) {
            return Err(boxed_error(format!(
                "{} repeats source id `{}`",
                descriptor.id, entry.id
            )));
        }
        *categories.entry(entry.category.clone()).or_default() += 1;
        for file in &entry.files {
            if !file.url.starts_with("https://") || !is_sha256(&file.sha256) {
                return Err(boxed_error(format!(
                    "{} entry `{}` has invalid source provenance",
                    descriptor.id, entry.id
                )));
            }
        }
    }
    if !descriptor.categories.is_empty() && categories != descriptor.categories {
        return Err(boxed_error(format!(
            "{} category counts differ: expected {:?}, found {:?}",
            descriptor.id, descriptor.categories, categories
        )));
    }
    Ok(())
}

pub(crate) fn check_corpus_artifacts(
    corpus: &str,
    lock: &SourceLock,
    require_data: bool,
    build_command: &str,
) -> Result<(), Box<dyn Error>> {
    let root = corpus_root(corpus);
    let manifest_dir = root.join("features");
    if !manifest_dir.exists() {
        return Err(boxed_error(format!(
            "{} has no benchmark manifests",
            root.display()
        )));
    }
    let mut golden_paths = Vec::new();
    for entry in fs::read_dir(&manifest_dir)? {
        let path = entry?.path();
        if path.extension().and_then(|ext| ext.to_str()) != Some("toml") {
            continue;
        }
        let manifest = read_benchmark_manifest(&path)?;
        let benchmark_id = path
            .file_stem()
            .and_then(|stem| stem.to_str())
            .ok_or_else(|| {
                boxed_error(format!(
                    "{} has a non-UTF-8 benchmark manifest name",
                    path.display()
                ))
            })?;
        if manifest.benchmark_id != benchmark_id {
            return Err(boxed_error(format!(
                "{} declares legacy feature_id `{}`, expected benchmark ID `{benchmark_id}`",
                path.display(),
                manifest.benchmark_id
            )));
        }
        if manifest.corpus_id != corpus {
            return Err(boxed_error(format!(
                "{} declares corpus `{}`",
                path.display(),
                manifest.corpus_id
            )));
        }
        for fixture in &manifest.fixtures {
            let golden = root
                .join("golden")
                .join(benchmark_id)
                .join(format!("{}.json.gz", slugify_fixture(fixture)));
            if !golden.exists() {
                return Err(boxed_error(format!(
                    "{} is missing golden {}",
                    corpus,
                    golden.display()
                )));
            }
            golden_paths.push(golden);
        }
    }
    golden_paths.sort();
    validate_gzip_json_files(&golden_paths)?;
    if !require_data {
        return Ok(());
    }
    for entry in &lock.entries {
        for file in &entry.files {
            check_data_file(&root, &file.path, &file.sha256, build_command)?;
        }
    }
    for pack in &lock.packs {
        if pack.count != pack.members.len() {
            return Err(boxed_error(format!(
                "{} pack `{}` count does not match members",
                corpus, pack.path
            )));
        }
        check_data_file(&root, &pack.path, &pack.sha256, build_command)?;
        let actual_members = read_pack_members(&root.join(&pack.path), pack)?;
        if actual_members != pack.members {
            return Err(boxed_error(format!(
                "{} pack `{}` member order differs from sources.lock.json",
                corpus, pack.path
            )));
        }
    }
    Ok(())
}

fn validate_gzip_json_files(paths: &[PathBuf]) -> Result<(), Box<dyn Error>> {
    let worker_count = std::thread::available_parallelism()
        .map(usize::from)
        .unwrap_or(1)
        .min(paths.len().max(1));
    if worker_count == 1 {
        for path in paths {
            validate_gzip_json(path)?;
        }
        return Ok(());
    }

    let next_path = std::sync::atomic::AtomicUsize::new(0);
    let results = std::sync::Mutex::new(vec![None; paths.len()]);
    std::thread::scope(|scope| {
        for _ in 0..worker_count {
            scope.spawn(|| loop {
                let index = next_path.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                if index >= paths.len() {
                    break;
                }
                let result = validate_gzip_json(&paths[index]).map_err(|error| error.to_string());
                results
                    .lock()
                    .expect("corpus golden result lock should not be poisoned")[index] =
                    Some(result);
            });
        }
    });
    for result in results
        .into_inner()
        .expect("corpus golden result lock should not be poisoned")
    {
        result
            .ok_or_else(|| boxed_error("corpus golden worker recorded no result"))?
            .map_err(boxed_error)?;
    }
    Ok(())
}

fn validate_gzip_json(path: &Path) -> Result<(), Box<dyn Error>> {
    let file = fs::File::open(path)?;
    let decoder = GzDecoder::new(file);
    let mut deserializer = serde_json::Deserializer::from_reader(decoder);
    serde::de::IgnoredAny::deserialize(&mut deserializer)?;
    deserializer.end()?;
    Ok(())
}

pub(crate) fn read_pack_members(
    path: &Path,
    pack: &SourcePack,
) -> Result<Vec<String>, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    match pack.format.as_str() {
        "smiles" => text
            .lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                let prefix = pack.member_title_prefix.as_deref().unwrap_or("CID:");
                line.split_whitespace()
                    .last()
                    .and_then(|title| title.strip_prefix(prefix))
                    .map(str::to_owned)
                    .ok_or_else(|| {
                        boxed_error(format!(
                            "{} contains a SMILES row without a `{prefix}` title",
                            path.display(),
                        ))
                    })
            })
            .collect(),
        "sdf-v2000" => {
            let property = pack
                .member_id_property
                .as_deref()
                .unwrap_or("PUBCHEM_COMPOUND_CID");
            let mut members = Vec::new();
            for record in text.split("$$$$") {
                if record.trim().is_empty() {
                    continue;
                }
                let mut lines = record.lines();
                let id = loop {
                    let Some(line) = lines.next() else {
                        return Err(boxed_error(format!(
                            "{} contains an SDF record without `{property}`",
                            path.display(),
                        )));
                    };
                    if sdf_data_header_name(line) != Some(property) {
                        continue;
                    }
                    break lines
                        .by_ref()
                        .find(|value| !value.trim().is_empty())
                        .ok_or_else(|| boxed_error(format!("missing `{property}` value")))?;
                };
                members.push(id.trim().to_owned());
            }
            Ok(members)
        }
        value => Err(boxed_error(format!(
            "{} uses unsupported pack format `{value}`",
            path.display()
        ))),
    }
}

fn sdf_data_header_name(line: &str) -> Option<&str> {
    let trimmed = line.trim_start();
    if !trimmed.starts_with('>') {
        return None;
    }
    let start = trimmed.find('<')?;
    let end = trimmed[start + 1..].find('>')? + start + 1;
    Some(&trimmed[start + 1..end])
}

pub(crate) fn check_data_file(
    corpus_root: &Path,
    relative: &str,
    expected_hash: &str,
    build_command: &str,
) -> Result<(), Box<dyn Error>> {
    let path = corpus_root.join(relative);
    if !path.exists() {
        return Err(boxed_error(format!(
            "{} is missing; build it with `{build_command}`",
            path.display()
        )));
    }
    let actual = hash_file(&path)?;
    if actual != expected_hash {
        return Err(boxed_error(format!(
            "{} checksum differs: expected {expected_hash}, found {actual}",
            path.display()
        )));
    }
    Ok(())
}

pub(crate) fn check_nested_corpora(
    locks: &BTreeMap<String, SourceLock>,
) -> Result<(), Box<dyn Error>> {
    for (child, parent) in [
        ("pubchem-100", "pubchem-1k"),
        ("pdb-10", "pdb-100"),
        ("pdb-100", "pdb-1000"),
    ] {
        let (Some(child_lock), Some(parent_lock)) = (locks.get(child), locks.get(parent)) else {
            continue;
        };
        let child_ids = child_lock
            .entries
            .iter()
            .map(|entry| entry.id.as_str())
            .collect::<Vec<_>>();
        let parent_ids = parent_lock
            .entries
            .iter()
            .take(child_ids.len())
            .map(|entry| entry.id.as_str())
            .collect::<Vec<_>>();
        if child_ids != parent_ids {
            return Err(boxed_error(format!(
                "{child} is not an exact prefix of {parent}"
            )));
        }
    }
    Ok(())
}
