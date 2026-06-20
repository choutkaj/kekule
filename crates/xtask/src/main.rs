use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::error::Error;
use std::fs;
use std::io::Read;
use std::path::{Path, PathBuf};
use std::process;
use std::time::{SystemTime, UNIX_EPOCH};

use flate2::read::GzDecoder;
use molecules::prelude::{
    perceive_aromaticity, perceive_ring_membership, perceive_ring_set, perceive_valence,
    read_mmcif_str, read_mol_v2000_str, read_smiles_str, sanitize_small_molecule, write_mol_v2000,
    write_sdf_v2000, write_smiles, AromaticityModel, Atom, AtomId, AtomRadical, Bond, BondOrder,
    BondStereo, MacroMolecule, MmcifParseOptions, Molecule, PropValue, SanitizeOptions,
    SdfParseOptions, SdfRecord, SmallMolecule, SmilesParseOptions, SmilesWriteOptions,
    ValenceModel,
};
use molecules::read_sdf_v2000_records;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};

const VALIDATION_CORPORA: &[(&str, &str)] = &[
    ("tiny", "Tiny"),
    ("pubchem-100", "PubChem 100"),
    ("pubchem-1000", "PubChem 1000"),
    ("pl-rex", "PL-REX"),
    ("enamine-diversity", "Enamine diversity"),
    ("pdb-10", "PDB 10"),
    ("pdb-100", "PDB 100"),
];
const DASHBOARD_PATH: &str = "features/DASHBOARD.html";
const VALIDATION_EVIDENCE_SCHEMA_VERSION: u32 = 2;
const GOLDEN_SCHEMA_VERSION: u32 = 1;
const COMPARISON_MODE_IMPLEMENTATION_GOLDEN: &str = "implementation-golden";

fn main() {
    if let Err(error) = run() {
        eprintln!("error: {error}");
        process::exit(1);
    }
}

fn run() -> Result<(), Box<dyn Error>> {
    let mut args = env::args().skip(1);
    match args.next().as_deref() {
        Some("dashboard") => dashboard(args.collect()),
        Some("validate") => validate(args.collect()),
        Some("corpus") => corpus(args.collect()),
        Some("features") => list_features(),
        Some("skills") => skills(args.collect()),
        _ => {
            print_help();
            Ok(())
        }
    }
}

fn print_help() {
    eprintln!(
        "usage:\n  cargo xtask dashboard [--check]\n  cargo xtask validate --feature FEATURE_ID|all [--corpus CORPUS_ID|all] [--update]\n  cargo xtask corpus check --corpus CORPUS_ID|all [--require-data]\n  cargo xtask skills --check\n  cargo xtask features"
    );
}

fn dashboard(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    let check = args.iter().any(|arg| arg == "--check");
    let features = read_features()?;
    let statuses = read_validation_statuses(&features)?;
    ensure_validation_flags_synced(&features, &statuses)?;
    let rendered = render_dashboard(&features, &statuses);
    let path = Path::new(DASHBOARD_PATH);

    if check {
        let existing = fs::read_to_string(path)?;
        if existing != rendered {
            return Err(boxed_error(
                "features/DASHBOARD.html is out of date; run `cargo xtask dashboard`",
            ));
        }
    } else {
        write_atomic_text(path, &rendered)?;
    }
    Ok(())
}

fn validate(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    validate_args(&args)?;
    let feature_selector = value_after_flag(&args, "--feature")
        .ok_or_else(|| boxed_error("missing required flag: --feature FEATURE_ID"))?;
    let corpus_selector = value_after_flag(&args, "--corpus").unwrap_or("tiny");
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
        let result = (|| -> Result<ValidationRun, Box<dyn Error>> {
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
            let compared = validate_golden_outputs(&manifest_path, &manifest)?;
            if compared > 0 {
                println!("validation compared {compared} golden file(s)");
            }
            if compared != manifest.fixtures.len() {
                return Err(boxed_error(format!(
                    "{} compared {compared} fixture(s), expected {}",
                    manifest_path.display(),
                    manifest.fixtures.len()
                )));
            }
            let evidence = build_validation_evidence(Path::new("."), &manifest_path, &manifest)?;
            Ok(ValidationRun {
                fixture_count: manifest.fixtures.len(),
                compared_count: compared,
                reference_tool: manifest.reference_tool,
                reference_version: manifest.reference_version,
                manifest_hash: hash_file(&manifest_path)?,
                evidence,
            })
        })();

        match result {
            Ok(run) => {
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
            Err(error) => failures.push(format!("{} [{corpus}]: {error}", feature.id)),
        }
    }

    if update {
        write_validation_statuses(&statuses, &update_corpora)?;
        sync_feature_validation_flags(&features, &statuses)?;
        let refreshed_features = read_features()?;
        let rendered = render_dashboard(&refreshed_features, &statuses);
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

fn list_features() -> Result<(), Box<dyn Error>> {
    let features = read_features()?;
    let statuses = read_validation_statuses(&features)?;
    for feature in &features {
        println!(
            "{}\t{}\tv{}\timplemented={}\tvalidated={}",
            feature.id,
            feature.area,
            feature.version,
            feature.implemented,
            overall_validated(feature, statuses.get(&feature.id))
        );
    }
    Ok(())
}

fn skills(args: Vec<String>) -> Result<(), Box<dyn Error>> {
    if args.iter().any(|arg| arg != "--check") {
        return Err(boxed_error("usage: cargo xtask skills --check"));
    }
    check_skills(Path::new(".codex/skills"))?;
    println!("repo-local feature skills are in sync");
    Ok(())
}

fn value_after_flag<'a>(args: &'a [String], flag: &str) -> Option<&'a str> {
    args.windows(2)
        .find(|window| window[0] == flag)
        .map(|window| window[1].as_str())
}

fn validate_args(args: &[String]) -> Result<(), Box<dyn Error>> {
    let mut index = 0;
    while index < args.len() {
        match args[index].as_str() {
            "--feature" | "--corpus" => {
                if index + 1 >= args.len() {
                    return Err(boxed_error(format!("missing value after {}", args[index])));
                }
                index += 2;
            }
            "--update" => index += 1,
            arg => return Err(boxed_error(format!("unknown validate argument: {arg}"))),
        }
    }
    Ok(())
}

fn corpus(args: Vec<String>) -> Result<(), Box<dyn Error>> {
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
    if selector != "all" && !is_known_corpus(selector) {
        return Err(boxed_error(format!("unknown corpus: {selector}")));
    }

    let corpora = VALIDATION_CORPORA
        .iter()
        .map(|(id, _)| *id)
        .filter(|id| selector == "all" || selector == *id)
        .collect::<Vec<_>>();
    let mut locks = BTreeMap::new();
    for corpus_id in &corpora {
        let descriptor = read_corpus_descriptor(corpus_id)?;
        if descriptor.id != *corpus_id {
            return Err(boxed_error(format!(
                "{} declares id `{}`, expected `{corpus_id}`",
                corpus_descriptor_path(corpus_id).display(),
                descriptor.id
            )));
        }
        if !descriptor.ready {
            println!("corpus `{corpus_id}` is declared but not built; skipping integrity checks");
            continue;
        }
        let lock = read_source_lock(corpus_id)?;
        check_corpus_lock(&descriptor, &lock)?;
        check_corpus_artifacts(corpus_id, &lock, require_data, &descriptor.build_command)?;
        println!(
            "corpus `{corpus_id}` has {} pinned entries and passed integrity checks",
            lock.entries.len()
        );
        locks.insert((*corpus_id).to_owned(), lock);
    }
    check_nested_corpora(&locks)?;
    Ok(())
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
struct CorpusDescriptor {
    id: String,
    title: String,
    kind: String,
    ready: bool,
    expected_count: usize,
    #[serde(default)]
    parent: Option<String>,
    #[serde(default)]
    seed: Option<String>,
    #[serde(default)]
    formats: Vec<String>,
    #[serde(default)]
    categories: BTreeMap<String, usize>,
    #[serde(default, rename = "notes")]
    _notes: Vec<String>,
    build_command: String,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct SourceLock {
    schema_version: u32,
    corpus_id: String,
    source: String,
    selection_seed: String,
    entries: Vec<SourceEntry>,
    #[serde(default)]
    packs: Vec<SourcePack>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct SourceEntry {
    id: String,
    category: String,
    files: Vec<SourceFile>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct SourceFile {
    path: String,
    url: String,
    sha256: String,
    #[serde(default)]
    record_type: Option<String>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct SourcePack {
    path: String,
    format: String,
    count: usize,
    members: Vec<String>,
    sha256: String,
}

fn corpus_root(corpus: &str) -> PathBuf {
    Path::new("validation").join("corpora").join(corpus)
}

fn corpus_descriptor_path(corpus: &str) -> PathBuf {
    corpus_root(corpus).join("corpus.toml")
}

fn read_corpus_descriptor(corpus: &str) -> Result<CorpusDescriptor, Box<dyn Error>> {
    let path = corpus_descriptor_path(corpus);
    let text = fs::read_to_string(&path)?;
    toml::from_str(&text).map_err(|error| boxed_error(format!("{}: {error}", path.display())))
}

fn read_source_lock(corpus: &str) -> Result<SourceLock, Box<dyn Error>> {
    let path = corpus_root(corpus).join("sources.lock.json");
    let text = fs::read_to_string(&path)
        .map_err(|error| boxed_error(format!("{} is unavailable: {error}", path.display())))?;
    serde_json::from_str(&text).map_err(|error| boxed_error(format!("{}: {error}", path.display())))
}

fn check_corpus_lock(
    descriptor: &CorpusDescriptor,
    lock: &SourceLock,
) -> Result<(), Box<dyn Error>> {
    if descriptor.title.trim().is_empty()
        || descriptor.kind.trim().is_empty()
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
    if descriptor.seed.as_deref() != Some(lock.selection_seed.as_str()) {
        return Err(boxed_error(format!(
            "{} selection seed does not match corpus.toml",
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

fn check_corpus_artifacts(
    corpus: &str,
    lock: &SourceLock,
    require_data: bool,
    build_command: &str,
) -> Result<(), Box<dyn Error>> {
    let root = corpus_root(corpus);
    let features_dir = root.join("features");
    if !features_dir.exists() {
        return Err(boxed_error(format!(
            "{} has no feature manifests",
            root.display()
        )));
    }
    for entry in fs::read_dir(&features_dir)? {
        let path = entry?.path();
        if path.extension().and_then(|ext| ext.to_str()) != Some("toml") {
            continue;
        }
        let manifest = read_validation_manifest(&path)?;
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
                .join(&manifest.feature_id)
                .join(format!("{}.json.gz", slugify_fixture(fixture)));
            if !golden.exists() {
                return Err(boxed_error(format!(
                    "{} is missing golden {}",
                    corpus,
                    golden.display()
                )));
            }
            let _: Value = serde_json::from_str(&read_gzip_string(&golden)?)?;
        }
    }
    if validation_status_path(corpus).exists() {
        read_corpus_status(&validation_status_path(corpus))?;
    }
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
        let actual_members = read_pack_members(&root.join(&pack.path), &pack.format)?;
        if actual_members != pack.members {
            return Err(boxed_error(format!(
                "{} pack `{}` member order differs from sources.lock.json",
                corpus, pack.path
            )));
        }
    }
    Ok(())
}

fn read_pack_members(path: &Path, format: &str) -> Result<Vec<String>, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    match format {
        "smiles" => text
            .lines()
            .filter(|line| !line.trim().is_empty())
            .map(|line| {
                line.split_whitespace()
                    .last()
                    .and_then(|title| title.strip_prefix("CID:"))
                    .map(str::to_owned)
                    .ok_or_else(|| {
                        boxed_error(format!(
                            "{} contains a SMILES row without a CID title",
                            path.display()
                        ))
                    })
            })
            .collect(),
        "sdf-v2000" => {
            let marker = "> <PUBCHEM_COMPOUND_CID>";
            let mut members = Vec::new();
            for record in text.split("$$$$") {
                if record.trim().is_empty() {
                    continue;
                }
                let position = record.find(marker).ok_or_else(|| {
                    boxed_error(format!(
                        "{} contains an SDF record without PUBCHEM_COMPOUND_CID",
                        path.display()
                    ))
                })?;
                let cid = record[position + marker.len()..]
                    .trim_start_matches(['\r', '\n'])
                    .lines()
                    .next()
                    .ok_or_else(|| boxed_error("missing PubChem CID value"))?;
                members.push(cid.trim().to_owned());
            }
            Ok(members)
        }
        value => Err(boxed_error(format!(
            "{} uses unsupported pack format `{value}`",
            path.display()
        ))),
    }
}

fn check_data_file(
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

fn check_nested_corpora(locks: &BTreeMap<String, SourceLock>) -> Result<(), Box<dyn Error>> {
    for (child, parent) in [("pubchem-100", "pubchem-1000"), ("pdb-10", "pdb-100")] {
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

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Feature {
    id: String,
    title: String,
    area: String,
    version: u32,
    implemented: bool,
    validated: bool,
    description: String,
    depends_on: Vec<String>,
    validation_required: Vec<String>,
}

fn read_features() -> Result<Vec<Feature>, Box<dyn Error>> {
    read_features_from(Path::new("features"))
}

fn read_features_from(root: &Path) -> Result<Vec<Feature>, Box<dyn Error>> {
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

fn is_hidden_or_template(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .map(|name| name.starts_with('_'))
        .unwrap_or(false)
}

fn read_feature(path: &Path) -> Result<Feature, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    let feature: Feature = toml::from_str(&text)
        .map_err(|error| boxed_error(format!("{}: {error}", path.display())))?;
    validate_feature(&feature, path)?;
    Ok(feature)
}

fn validate_feature(feature: &Feature, path: &Path) -> Result<(), Box<dyn Error>> {
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
    let feature_doc = path.with_file_name("feature.md");
    if !feature_doc.exists() {
        return Err(boxed_error(format!(
            "{} is missing required feature.md",
            path.parent().unwrap_or_else(|| Path::new(".")).display()
        )));
    }
    Ok(())
}

fn validate_feature_set(features: &[Feature]) -> Result<(), Box<dyn Error>> {
    let mut seen = BTreeMap::<&str, ()>::new();
    for feature in features {
        if seen.insert(feature.id.as_str(), ()).is_some() {
            return Err(boxed_error(format!(
                "duplicate feature id `{}`",
                feature.id
            )));
        }
    }
    for feature in features {
        for dependency in &feature.depends_on {
            if !seen.contains_key(dependency.as_str()) {
                return Err(boxed_error(format!(
                    "feature `{}` depends on unknown feature `{dependency}`",
                    feature.id
                )));
            }
        }
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct ValidationManifest {
    feature_id: String,
    corpus_id: String,
    reference_tool: String,
    reference_version: String,
    comparison_mode: String,
    #[serde(default)]
    fixtures: Vec<String>,
    #[serde(default, rename = "notes")]
    _notes: Vec<String>,
}

#[derive(Debug)]
struct ValidationRun {
    fixture_count: usize,
    compared_count: usize,
    reference_tool: String,
    reference_version: String,
    manifest_hash: String,
    evidence: ValidationEvidence,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct ValidationEvidence {
    schema_version: u32,
    comparison_mode: String,
    inputs: Vec<EvidenceInput>,
    sha256: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct EvidenceInput {
    path: String,
    sha256: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct CorpusStatus {
    passed: bool,
    fixture_count: usize,
    compared_count: usize,
    reference_tool: String,
    reference_version: String,
    manifest_hash: String,
    #[serde(default)]
    evidence_schema_version: Option<u32>,
    #[serde(default)]
    evidence_hash: Option<String>,
    #[serde(default)]
    evidence_inputs: Vec<EvidenceInput>,
    validated_at_unix: u64,
}

impl CorpusStatus {
    fn from_run(
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
            evidence_schema_version: Some(run.evidence.schema_version),
            evidence_hash: Some(run.evidence.sha256),
            evidence_inputs: run.evidence.inputs,
            validated_at_unix,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ValidationStatus {
    feature_id: String,
    corpora: BTreeMap<String, CorpusStatus>,
}

impl ValidationStatus {
    fn new(feature_id: &str) -> Self {
        Self {
            feature_id: feature_id.to_owned(),
            corpora: BTreeMap::new(),
        }
    }
}

fn is_known_corpus(corpus: &str) -> bool {
    VALIDATION_CORPORA
        .iter()
        .any(|(candidate, _)| *candidate == corpus)
}

fn validation_manifest_path(feature: &str, corpus: &str) -> PathBuf {
    Path::new("validation")
        .join("corpora")
        .join(corpus)
        .join("features")
        .join(format!("{feature}.toml"))
}

fn read_validation_manifest(path: &Path) -> Result<ValidationManifest, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    toml::from_str(&text).map_err(|error| boxed_error(format!("{}: {error}", path.display())))
}

fn validate_comparison_mode(
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

fn validation_targets<'a>(
    features: &'a [Feature],
    feature_selector: &str,
    corpus_selector: &str,
) -> Vec<(&'a Feature, String)> {
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

fn hash_file(path: &Path) -> Result<String, Box<dyn Error>> {
    let mut hasher = Sha256::new();
    hasher.update(fs::read(path)?);
    Ok(format!("{:x}", hasher.finalize()))
}

fn build_validation_evidence(
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

    let reference_root = match manifest.reference_tool.as_str() {
        "rdkit" => repo_root.join("validation/reference/rdkit"),
        "biopython" => repo_root.join("validation/reference/biopython"),
        value => {
            return Err(boxed_error(format!(
                "{} uses unsupported reference_tool `{value}`",
                manifest_path.display()
            )))
        }
    };
    paths.insert(reference_root.join("run_feature.py"));
    paths.insert(reference_root.join("environment.yml"));

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

fn collect_files(root: &Path, paths: &mut BTreeSet<PathBuf>) -> Result<(), Box<dyn Error>> {
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

fn hash_evidence_file(path: &Path) -> Result<String, Box<dyn Error>> {
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

fn relative_path(repo_root: &Path, path: &Path) -> Result<String, Box<dyn Error>> {
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

fn read_gzip_string(path: &Path) -> Result<String, Box<dyn Error>> {
    let file = fs::File::open(path)?;
    let mut decoder = GzDecoder::new(file);
    let mut text = String::new();
    decoder.read_to_string(&mut text)?;
    Ok(text)
}

fn validate_manifest_paths(
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

fn validate_golden_outputs(
    manifest_path: &Path,
    manifest: &ValidationManifest,
) -> Result<usize, Box<dyn Error>> {
    if manifest.fixtures.is_empty() {
        return Ok(0);
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
    let mut compared = 0;
    for fixture in &manifest.fixtures {
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
        let expected = golden.get("expected").ok_or_else(|| {
            boxed_error(format!("{} is missing `expected`", golden_path.display()))
        })?;
        let actual = implementation_expected(&manifest.feature_id, &fixture_path)?;
        let expected = normalize_for_comparison(expected);
        let actual = normalize_for_comparison(&actual);
        if let Some(diff) = first_json_diff("$", &expected, &actual) {
            return Err(boxed_error(format!(
                "{} differs from implementation output for fixture `{fixture}`: {diff}",
                golden_path.display()
            )));
        }
        compared += 1;
    }
    Ok(compared)
}

fn validate_golden_metadata(
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

fn reference_version_label(tool: &str, version: &str) -> String {
    match tool {
        "rdkit" if !version.starts_with("RDKit ") => format!("RDKit {version}"),
        "biopython" if !version.starts_with("Biopython ") => format!("Biopython {version}"),
        _ => version.to_owned(),
    }
}

fn implementation_expected(feature: &str, fixture_path: &Path) -> Result<Value, Box<dyn Error>> {
    match feature {
        "io.sdf.v2000.parse" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({ "records": records.iter().map(sdf_record_json).collect::<Vec<_>>() }))
        }
        "io.sdf.v2000.write" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            let molecules = records
                .into_iter()
                .map(|record| record.molecule)
                .collect::<Vec<_>>();
            let written = write_sdf_v2000(&molecules)?;
            let records = read_sdf_v2000_records(&written, SdfParseOptions::default())?
                .into_iter()
                .enumerate()
                .map(|(index, record)| IndexedSmallRecord {
                    record_index: index,
                    title: record.title,
                    molecule: record.molecule,
                })
                .collect::<Vec<_>>();
            Ok(json!({ "records": records.iter().map(sdf_record_json).collect::<Vec<_>>() }))
        }
        "io.mol.v2000.parse" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({ "records": records.iter().map(mol_parse_record_json).collect::<Vec<_>>() }))
        }
        "core.conformers" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({ "records": records.iter().map(conformer_record_json).collect::<Vec<_>>() }))
        }
        "io.mol.v2000.write" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            let records = records
                .into_iter()
                .enumerate()
                .map(|(index, record)| {
                    let written = write_mol_v2000(&record.molecule)?;
                    let molecule = read_mol_v2000_str(&written)?;
                    Ok(IndexedSmallRecord {
                        record_index: index,
                        title: molecule_title(&molecule.mol),
                        molecule,
                    })
                })
                .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
            Ok(json!({ "records": records.iter().map(mol_record_json).collect::<Vec<_>>() }))
        }
        "io.smiles.parse" => {
            let records = read_smiles_records(fixture_path)?;
            Ok(
                json!({ "records": records.iter().map(smiles_parse_record_json).collect::<Vec<_>>() }),
            )
        }
        "io.smiles.write" => {
            let records = read_smiles_records(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter()
                    .map(smiles_write_record_json)
                    .collect::<Result<Vec<_>, Box<dyn Error>>>()?
            }))
        }
        "algo.rings.fast" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(ring_membership_record_json).collect::<Vec<_>>() }),
            )
        }
        "algo.rings.sssr" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(ring_set_record_json).collect::<Vec<_>>() }),
            )
        }
        "algo.valence.rdkit-like" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(valence_record_json).collect::<Vec<_>>() }),
            )
        }
        "chem.sanitize.rdkit-like" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(sanitized_atom_record_json).collect::<Vec<_>>() }),
            )
        }
        "algo.aromaticity.rdkit-like" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(aromaticity_record_json).collect::<Vec<_>>() }),
            )
        }
        "io.mmcif.parse" | "bio.hierarchy.smcra" => {
            let input = fs::read_to_string(fixture_path)?;
            let molecule = read_mmcif_str(&input, MmcifParseOptions::default())?;
            Ok(mmcif_expected_json(&molecule))
        }
        _ => Err(boxed_error(format!(
            "no implementation comparison configured for feature `{feature}`"
        ))),
    }
}

#[derive(Debug, Clone)]
struct IndexedSmallRecord {
    record_index: usize,
    title: String,
    molecule: SmallMolecule,
}

#[derive(Debug, Clone)]
struct IndexedSmilesRecord {
    record_index: usize,
    status: String,
    title: String,
    input_smiles: String,
    molecule: Option<SmallMolecule>,
}

fn read_small_records_by_suffix(path: &Path) -> Result<Vec<IndexedSmallRecord>, Box<dyn Error>> {
    let input = fs::read_to_string(path)?;
    if matches!(
        path.extension().and_then(|ext| ext.to_str()),
        Some("mol" | "mdl")
    ) {
        let molecule = read_mol_v2000_str(&input)?;
        return Ok(vec![IndexedSmallRecord {
            record_index: 0,
            title: molecule_title(&molecule.mol),
            molecule,
        }]);
    }
    Ok(read_sdf_v2000_records(&input, SdfParseOptions::default())?
        .into_iter()
        .enumerate()
        .map(|(index, record)| small_record(index, record))
        .collect())
}

fn small_record(index: usize, record: SdfRecord) -> IndexedSmallRecord {
    IndexedSmallRecord {
        record_index: index,
        title: record.title,
        molecule: record.molecule,
    }
}

fn read_smiles_records(path: &Path) -> Result<Vec<IndexedSmilesRecord>, Box<dyn Error>> {
    let mut records = Vec::new();
    for (index, raw_line) in fs::read_to_string(path)?.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut parts = line.splitn(2, char::is_whitespace);
        let smiles = parts.next().unwrap_or_default().to_owned();
        let title = parts.next().unwrap_or_default().trim().to_owned();
        if smiles_unsupported_subset_reason(&smiles).is_some() {
            records.push(IndexedSmilesRecord {
                record_index: index,
                status: "unsupported".to_owned(),
                title,
                input_smiles: smiles,
                molecule: None,
            });
            continue;
        }
        let (status, molecule) = match read_smiles_str(&smiles, SmilesParseOptions) {
            Ok(molecule) => ("ok".to_owned(), Some(molecule)),
            Err(_) => ("parse_error".to_owned(), None),
        };
        records.push(IndexedSmilesRecord {
            record_index: index,
            status,
            title,
            input_smiles: smiles,
            molecule,
        });
    }
    Ok(records)
}

fn smiles_unsupported_subset_reason(smiles: &str) -> Option<&'static str> {
    smiles
        .chars()
        .any(|ch| matches!(ch, '@' | '/' | '\\' | '*'))
        .then_some("unsupported")
}

fn sdf_record_json(record: &IndexedSmallRecord) -> Value {
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "atoms": atoms_json(mol),
        "bonds": bonds_json(mol),
        "properties": sdf_properties_json(mol),
    })
}

fn mol_record_json(record: &IndexedSmallRecord) -> Value {
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "atoms": atoms_json(mol),
        "bonds": bonds_json(mol),
    })
}

fn conformer_record_json(record: &IndexedSmallRecord) -> Value {
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "conformers": mol.conformers().map(|(_, conformer)| {
            mol.atom_ids()
                .filter_map(|atom_id| {
                    conformer.position(atom_id).map(|point| json!({
                        "atom_index": atom_id.raw(),
                        "x": point.x,
                        "y": point.y,
                        "z": point.z,
                    }))
                })
                .collect::<Vec<_>>()
        }).collect::<Vec<_>>(),
        "atoms": mol.atoms().map(|(id, atom)| conformer_atom_json(id, atom)).collect::<Vec<_>>(),
    })
}

fn mol_parse_record_json(record: &IndexedSmallRecord) -> Value {
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_count": mol.atom_count(),
        "conformers": conformers_json(mol),
        "atoms": atoms_json(mol),
    })
}

fn conformers_json(mol: &Molecule) -> Vec<Vec<Value>> {
    mol.conformers()
        .map(|(_, conformer)| {
            mol.atom_ids()
                .filter_map(|atom_id| {
                    conformer.position(atom_id).map(|point| {
                        json!({
                            "atom_index": atom_id.raw(),
                            "x": point.x,
                            "y": point.y,
                            "z": point.z,
                        })
                    })
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>()
}

fn conformer_atom_json(id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "isotope": atom.isotope,
        "explicit_hydrogens": atom.explicit_hydrogens,
        "atom_map": atom.atom_map,
        "aromatic": atom.aromatic,
    })
}

fn ring_membership_record_json(record: &mut IndexedSmallRecord) -> Value {
    let membership = perceive_ring_membership(&mut record.molecule.mol);
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_in_ring": mol.atom_ids().map(|id| membership.atom_in_ring(id)).collect::<Vec<_>>(),
        "bond_in_ring": mol.bond_ids().map(|id| membership.bond_in_ring(id)).collect::<Vec<_>>(),
    })
}

fn ring_set_record_json(record: &mut IndexedSmallRecord) -> Value {
    match perceive_ring_set(&mut record.molecule.mol) {
        Ok(ring_set) => json!({
            "record_index": record.record_index,
            "status": "ok",
            "title": record.title,
            "rings": ring_set
                .rings()
                .iter()
                .map(|ring| ring.atoms.iter().map(|atom| atom.raw()).collect::<Vec<_>>())
                .collect::<Vec<_>>(),
        }),
        Err(_) => json!({
            "record_index": record.record_index,
            "status": "resource_error",
            "title": record.title,
        }),
    }
}

fn sanitized_atom_record_json(record: &mut IndexedSmallRecord) -> Value {
    match sanitize_small_molecule(&mut record.molecule, SanitizeOptions::default()) {
        Ok(_) => json!({
            "record_index": record.record_index,
            "status": "ok",
            "title": record.title,
            "atoms": basic_atoms_json(&record.molecule.mol),
        }),
        Err(_) => json!({
            "record_index": record.record_index,
            "status": "sanitize_error",
            "title": record.title,
        }),
    }
}

fn valence_record_json(record: &mut IndexedSmallRecord) -> Value {
    let report = perceive_valence(&mut record.molecule.mol, ValenceModel::RdkitLike);
    if !report.is_ok() {
        return json!({
            "record_index": record.record_index,
            "status": "valence_error",
            "title": record.title,
        });
    }
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atoms": record
            .molecule
            .mol
            .atoms()
            .map(|(id, atom)| valence_atom_json(&record.molecule.mol, id, atom))
            .collect::<Vec<_>>(),
    })
}

fn aromaticity_record_json(record: &mut IndexedSmallRecord) -> Value {
    let status = sanitize_small_molecule(
        &mut record.molecule,
        SanitizeOptions {
            perceive_valence: true,
            perceive_rings: true,
            perceive_aromaticity: false,
        },
    )
    .and_then(|_| {
        perceive_aromaticity(&mut record.molecule.mol, AromaticityModel::RdkitLike)
            .map_err(molecules::prelude::SanitizeError::Aromaticity)
    });
    if status.is_err() {
        return json!({
            "record_index": record.record_index,
            "status": "sanitize_error",
            "title": record.title,
        });
    }
    let mol = &record.molecule.mol;
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "atom_aromatic": mol.atoms().map(|(_, atom)| atom.aromatic).collect::<Vec<_>>(),
        "bond_aromatic": mol.bonds().map(|(_, bond)| bond.aromatic).collect::<Vec<_>>(),
    })
}

fn smiles_write_record_json(record: &IndexedSmilesRecord) -> Result<Value, Box<dyn Error>> {
    let Some(molecule) = &record.molecule else {
        return Ok(smiles_error_record_json(record));
    };
    let written = write_smiles(molecule, SmilesWriteOptions)?;
    let reparsed = match read_smiles_str(&written, SmilesParseOptions) {
        Ok(reparsed) => reparsed,
        Err(_) => {
            return Ok(json!({
                "record_index": record.record_index,
                "status": "write_reparse_error",
                "title": record.title,
                "input_smiles": record.input_smiles,
            }));
        }
    };
    Ok(json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "input_smiles": record.input_smiles,
        "sanitized": smiles_sanitized_semantic_json(reparsed),
    }))
}

fn smiles_parse_record_json(record: &IndexedSmilesRecord) -> Value {
    let Some(molecule) = &record.molecule else {
        return smiles_error_record_json(record);
    };
    let written = write_smiles(molecule, SmilesWriteOptions);
    let round_trip = match written
        .as_ref()
        .map_err(|_| ())
        .and_then(|text| read_smiles_str(text, SmilesParseOptions).map_err(|_| ()))
    {
        Ok(reparsed) => smiles_sanitized_semantic_json(reparsed),
        Err(_) => json!({ "status": "write_reparse_error" }),
    };
    json!({
        "record_index": record.record_index,
        "status": "ok",
        "title": record.title,
        "input_smiles": record.input_smiles,
        "raw": smiles_raw_semantic_json(molecule),
        "sanitized": smiles_sanitized_semantic_json(molecule.clone()),
        "write_round_trip": round_trip,
    })
}

fn smiles_error_record_json(record: &IndexedSmilesRecord) -> Value {
    json!({
        "record_index": record.record_index,
        "status": record.status,
        "title": record.title,
        "input_smiles": record.input_smiles,
    })
}

fn smiles_raw_semantic_json(molecule: &SmallMolecule) -> Value {
    let mol = &molecule.mol;
    json!({
        "atom_count": mol.atom_count(),
        "bond_count": mol.bond_count(),
        "atoms": basic_atoms_json(mol),
        "bonds": basic_bonds_json(mol),
    })
}

fn smiles_sanitized_semantic_json(mut molecule: SmallMolecule) -> Value {
    match sanitize_small_molecule(&mut molecule, SanitizeOptions::default()) {
        Ok(_) => {
            let mol = &molecule.mol;
            json!({
                "status": "ok",
                "atom_count": mol.atom_count(),
                "bond_count": mol.bond_count(),
                "atoms": smiles_sanitized_atoms_json(mol),
                "bonds": smiles_sanitized_bonds_json(mol),
            })
        }
        Err(_) => json!({ "status": "sanitize_error" }),
    }
}

fn smiles_sanitized_bonds_json(mol: &Molecule) -> Vec<Value> {
    let mut bonds = mol
        .bonds()
        .map(|(_, bond)| {
            let left = mol.atom(bond.a()).expect("bond endpoint should exist");
            let right = mol.atom(bond.b()).expect("bond endpoint should exist");
            let mut endpoints = [
                smiles_sanitized_atom_key(mol, bond.a(), left),
                smiles_sanitized_atom_key(mol, bond.b(), right),
            ];
            endpoints.sort();
            json!({
                "endpoint_atoms": endpoints,
                "bond_type": smiles_semantic_bond_type(bond),
                "is_aromatic": bond.aromatic,
            })
        })
        .collect::<Vec<_>>();
    bonds.sort_by_key(|value| value.to_string());
    bonds
}

fn smiles_sanitized_atoms_json(mol: &Molecule) -> Vec<Value> {
    let mut atoms = mol
        .atoms()
        .map(|(id, atom)| {
            let (explicit_hydrogens, implicit_hydrogens) = smiles_effective_hydrogens(atom);
            let explicit_valence = explicit_valence_json(mol, id) + explicit_hydrogens;
            let mut neighbors = mol
                .incident_bonds(id)
                .expect("atom should exist")
                .map(|(_, bond)| {
                    let neighbor_id = if bond.a() == id { bond.b() } else { bond.a() };
                    let neighbor = mol.atom(neighbor_id).expect("bond endpoint should exist");
                    json!({
                        "atom": smiles_sanitized_atom_key(mol, neighbor_id, neighbor),
                        "bond_type": smiles_semantic_bond_type(bond),
                        "is_aromatic": bond.aromatic,
                    })
                })
                .collect::<Vec<_>>();
            neighbors.sort_by_key(|value| value.to_string());
            (
                smiles_sanitized_atom_key(mol, id, atom),
                json!({
                    "atomic_number": atom.element.atomic_number(),
                    "symbol": atom.element.symbol(),
                    "formal_charge": atom.formal_charge,
                    "isotope": atom.isotope,
                    "explicit_hydrogens": explicit_hydrogens,
                    "implicit_hydrogens": implicit_hydrogens,
                    "no_implicit_hydrogens": atom.no_implicit_hydrogens,
                    "explicit_valence": explicit_valence,
                    "atom_map": atom.atom_map,
                    "aromatic": atom.aromatic,
                    "neighbors": neighbors,
                }),
            )
        })
        .collect::<Vec<_>>();
    atoms.sort_by(|left, right| {
        left.0
            .cmp(&right.0)
            .then_with(|| left.1.to_string().cmp(&right.1.to_string()))
    });
    atoms.into_iter().map(|(_, value)| value).collect()
}

fn smiles_sanitized_atom_key(mol: &Molecule, id: AtomId, atom: &Atom) -> String {
    let (explicit_hydrogens, implicit_hydrogens) = smiles_effective_hydrogens(atom);
    let explicit_valence = explicit_valence_json(mol, id) + explicit_hydrogens;
    format!(
        "{:03}|{}|{}|{}|{}|{}|{}|{}|{}|{}",
        atom.element.atomic_number(),
        atom.element.symbol(),
        atom.formal_charge,
        atom.isotope.unwrap_or(0),
        explicit_hydrogens,
        implicit_hydrogens,
        atom.no_implicit_hydrogens,
        explicit_valence,
        atom.atom_map.unwrap_or(0),
        atom.aromatic
    )
}

fn smiles_semantic_bond_type(bond: &Bond) -> &'static str {
    if bond.aromatic {
        "AROMATIC"
    } else {
        bond_order_json(bond.order)
    }
}

fn smiles_effective_hydrogens(atom: &Atom) -> (u8, u8) {
    if atom.element.symbol() == "N"
        && atom.aromatic
        && atom.explicit_hydrogens == 0
        && atom.implicit_hydrogens == Some(1)
    {
        (1, 0)
    } else {
        (
            atom.explicit_hydrogens,
            atom.implicit_hydrogens.unwrap_or(0),
        )
    }
}

fn atoms_json(mol: &Molecule) -> Vec<Value> {
    mol.atoms()
        .map(|(id, atom)| atom_json(id, atom))
        .collect::<Vec<_>>()
}

fn atom_json(id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "isotope": atom.isotope,
        "explicit_hydrogens": atom.explicit_hydrogens,
        "atom_map": atom.atom_map,
        "radical": atom.radical.map(radical_json),
        "unpaired_electrons": atom.radical.map(AtomRadical::unpaired_electron_count).unwrap_or(0),
        "aromatic": atom.aromatic,
    })
}

fn basic_atoms_json(mol: &Molecule) -> Vec<Value> {
    mol.atoms()
        .map(|(id, atom)| basic_atom_json(id, atom))
        .collect::<Vec<_>>()
}

fn basic_atom_json(id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "isotope": atom.isotope,
        "explicit_hydrogens": atom.explicit_hydrogens,
        "atom_map": atom.atom_map,
        "aromatic": atom.aromatic,
    })
}

fn valence_atom_json(mol: &Molecule, id: AtomId, atom: &Atom) -> Value {
    json!({
        "index": id.raw(),
        "atomic_number": atom.element.atomic_number(),
        "symbol": atom.element.symbol(),
        "formal_charge": atom.formal_charge,
        "explicit_hydrogens": atom.explicit_hydrogens,
        "implicit_hydrogens": atom.implicit_hydrogens.unwrap_or(0),
        "explicit_valence": explicit_valence_json(mol, id) + atom.explicit_hydrogens,
    })
}

fn explicit_valence_json(mol: &Molecule, atom: AtomId) -> u8 {
    mol.incident_bonds(atom)
        .ok()
        .into_iter()
        .flatten()
        .map(|(_, bond)| match bond.order {
            BondOrder::Zero | BondOrder::Dative => 0,
            BondOrder::Single | BondOrder::Aromatic => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
        })
        .sum()
}

fn bonds_json(mol: &Molecule) -> Vec<Value> {
    mol.bonds()
        .map(|(id, bond)| bond_json(id.raw(), bond))
        .collect::<Vec<_>>()
}

fn bond_json(index: u32, bond: &Bond) -> Value {
    json!({
        "index": index,
        "begin_atom_index": bond.a().raw(),
        "end_atom_index": bond.b().raw(),
        "bond_type": bond_order_json(bond.order),
        "is_aromatic": bond.aromatic,
        "stereo": bond_stereo_json(bond.order, bond.stereo),
        "bond_direction": bond_direction_json(bond.order, bond.stereo),
    })
}

fn basic_bonds_json(mol: &Molecule) -> Vec<Value> {
    mol.bonds()
        .map(|(id, bond)| basic_bond_json(id.raw(), bond))
        .collect::<Vec<_>>()
}

fn basic_bond_json(index: u32, bond: &Bond) -> Value {
    json!({
        "index": index,
        "begin_atom_index": bond.a().raw(),
        "end_atom_index": bond.b().raw(),
        "bond_type": bond_order_json(bond.order),
        "is_aromatic": bond.aromatic,
        "stereo": legacy_bond_stereo_json(bond.stereo),
    })
}

fn legacy_bond_stereo_json(stereo: Option<BondStereo>) -> &'static str {
    match stereo {
        None | Some(BondStereo::Unspecified) => "STEREONONE",
        Some(BondStereo::E) => "STEREOE",
        Some(BondStereo::Z) => "STEREOZ",
        Some(BondStereo::Up) | Some(BondStereo::Down) | Some(BondStereo::Any) => "STEREOANY",
    }
}

fn radical_json(radical: AtomRadical) -> &'static str {
    match radical {
        AtomRadical::Singlet => "SINGLET",
        AtomRadical::Doublet => "DOUBLET",
        AtomRadical::Triplet => "TRIPLET",
    }
}

fn bond_order_json(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Zero => "ZERO",
        BondOrder::Single => "SINGLE",
        BondOrder::Double => "DOUBLE",
        BondOrder::Triple => "TRIPLE",
        BondOrder::Quadruple => "QUADRUPLE",
        BondOrder::Aromatic => "AROMATIC",
        BondOrder::Dative => "DATIVE",
    }
}

fn bond_stereo_json(order: BondOrder, stereo: Option<BondStereo>) -> &'static str {
    match (order, stereo) {
        (_, None | Some(BondStereo::Unspecified)) => "STEREONONE",
        (_, Some(BondStereo::E)) => "STEREOE",
        (_, Some(BondStereo::Z)) => "STEREOZ",
        (BondOrder::Double, Some(BondStereo::Any)) => "STEREOANY",
        _ => "STEREONONE",
    }
}

fn bond_direction_json(order: BondOrder, stereo: Option<BondStereo>) -> &'static str {
    match (order, stereo) {
        (BondOrder::Single, Some(BondStereo::Up)) => "BEGINWEDGE",
        (BondOrder::Single, Some(BondStereo::Down)) => "BEGINDASH",
        (BondOrder::Single, Some(BondStereo::Any)) => "UNKNOWN",
        _ => "NONE",
    }
}

fn sdf_properties_json(mol: &Molecule) -> Value {
    let mut props = serde_json::Map::new();
    for (key, value) in mol.props() {
        let Some(name) = key.strip_prefix("sdf.field.") else {
            continue;
        };
        if let PropValue::String(text) = value {
            props.insert(name.to_owned(), json!(text));
        }
    }
    Value::Object(props)
}

fn molecule_title(mol: &Molecule) -> String {
    match mol.props().get("sdf.title") {
        Some(PropValue::String(title)) => title.clone(),
        _ => String::new(),
    }
}

fn mmcif_expected_json(molecule: &MacroMolecule) -> Value {
    json!({
        "atom_site_rows": atom_site_rows_json(molecule),
        "structure": structure_json(molecule),
    })
}

fn atom_site_rows_json(molecule: &MacroMolecule) -> Value {
    let rows = molecule
        .hierarchy
        .atom_sites()
        .map(|(site_id, site)| {
            let residue = molecule
                .hierarchy
                .residue(site.residue)
                .expect("residue exists");
            let chain = molecule
                .hierarchy
                .chain(residue.chain)
                .expect("chain exists");
            let model = molecule.hierarchy.model(chain.model).expect("model exists");
            let atom = molecule.mol.atom(site.atom).expect("atom exists");
            let point = first_conformer_point(molecule, site.atom);
            json!({
                "group_PDB": site.metadata.group_pdb.clone(),
                "id": site.metadata.atom_site_id.clone().unwrap_or_else(|| (site_id.raw() + 1).to_string()),
                "type_symbol": site.metadata.type_symbol.clone().unwrap_or_else(|| atom.element.symbol().to_owned()),
                "label_atom_id": site.metadata.label_atom_id,
                "auth_atom_id": site.metadata.auth_atom_id,
                "label_alt_id": site.metadata.label_alt_id,
                "label_comp_id": residue.label_comp_id,
                "auth_comp_id": residue.author_comp_id,
                "label_asym_id": site.metadata.label_asym_id.clone().unwrap_or_else(|| chain.label_id.clone()),
                "auth_asym_id": site.metadata.auth_asym_id.clone().or_else(|| chain.author_id.clone()),
                "label_seq_id": residue.label_seq_id.map(|value| value.to_string()),
                "auth_seq_id": residue.author_seq_id,
                "pdbx_PDB_ins_code": residue.insertion_code,
                "occupancy": site.metadata.occupancy_raw.clone().or_else(|| site.metadata.occupancy.map(|value| format!("{value:.2}"))),
                "B_iso_or_equiv": site.metadata.b_factor_raw.clone().or_else(|| site.metadata.b_factor.map(|value| format!("{value:.2}"))),
                "Cartn_x": site.metadata.cartn_x_raw.clone().or_else(|| point.map(|point| format!("{:.3}", point.x))),
                "Cartn_y": site.metadata.cartn_y_raw.clone().or_else(|| point.map(|point| format!("{:.3}", point.y))),
                "Cartn_z": site.metadata.cartn_z_raw.clone().or_else(|| point.map(|point| format!("{:.3}", point.z))),
                "pdbx_PDB_model_num": model.model_id,
            })
        })
        .collect::<Vec<_>>();
    json!({
        "status": "ok",
        "row_count": rows.len(),
        "rows": rows,
    })
}

fn structure_json(molecule: &MacroMolecule) -> Value {
    json!({
        "status": "ok",
        "models": molecule.hierarchy.models().map(|(model_id, model)| {
            json!({
                "id": model_id.raw(),
                "chains": model.chains.iter().map(|chain_id| {
                    let chain = molecule.hierarchy.chain(*chain_id).expect("chain exists");
                    json!({
                        "id": chain.author_id.clone().unwrap_or_else(|| chain.label_id.clone()),
                        "residues": chain.residues.iter().map(|residue_id| {
                            let residue = molecule.hierarchy.residue(*residue_id).expect("residue exists");
                            json!({
                                "name": residue.name,
                                "hetflag": residue_hetflag_json(molecule, residue),
                                "sequence_id": residue_sequence_json(residue),
                                "insertion_code": residue.insertion_code,
                                "atoms": residue.atom_sites.iter().map(|site_id| {
                                    let site = molecule.hierarchy.atom_site(*site_id).expect("site exists");
                                    let atom = molecule.mol.atom(site.atom).expect("atom exists");
                                    let name = site
                                        .metadata
                                        .label_atom_id
                                        .clone()
                                        .unwrap_or_else(|| atom.element.symbol().to_owned());
                                    let coord = first_conformer_point(molecule, site.atom)
                                        .map(|point| json!([point.x, point.y, point.z]))
                                        .unwrap_or(Value::Null);
                                    json!({
                                        "name": name,
                                        "full_name": name,
                                        "altloc": site.metadata.label_alt_id,
                                        "element": site.metadata.type_symbol.clone().unwrap_or_else(|| atom.element.symbol().to_owned()),
                                        "occupancy": site.metadata.occupancy,
                                        "bfactor": site.metadata.b_factor,
                                        "coord": coord,
                                    })
                                }).collect::<Vec<_>>(),
                            })
                        }).collect::<Vec<_>>(),
                    })
                }).collect::<Vec<_>>(),
            })
        }).collect::<Vec<_>>(),
    })
}

fn first_conformer_point(
    molecule: &MacroMolecule,
    atom: AtomId,
) -> Option<molecules::prelude::Point3> {
    molecule
        .mol
        .first_conformer()
        .and_then(|(_, conformer)| conformer.position(atom))
}

fn residue_sequence_json(residue: &molecules::prelude::Residue) -> Value {
    if let Some(author_seq_id) = &residue.author_seq_id {
        return author_seq_id
            .parse::<i32>()
            .map(Value::from)
            .unwrap_or_else(|_| json!(author_seq_id));
    }
    residue.label_seq_id.map(Value::from).unwrap_or(Value::Null)
}

fn residue_hetflag_json(molecule: &MacroMolecule, residue: &molecules::prelude::Residue) -> Value {
    let is_hetatm = residue.atom_sites.iter().any(|site_id| {
        molecule
            .hierarchy
            .atom_site(*site_id)
            .ok()
            .and_then(|site| site.metadata.group_pdb.as_deref())
            == Some("HETATM")
    });
    if is_hetatm {
        if residue.name == "HOH" {
            json!("W")
        } else {
            json!(format!("H_{}", residue.name))
        }
    } else {
        Value::Null
    }
}

fn first_json_diff(path: &str, expected: &Value, actual: &Value) -> Option<String> {
    match (expected, actual) {
        (Value::Object(expected), Value::Object(actual)) => {
            for key in expected.keys() {
                let next = format!("{path}.{key}");
                let Some(actual_value) = actual.get(key) else {
                    return Some(format!("{next} missing from actual output"));
                };
                if let Some(diff) = first_json_diff(&next, &expected[key], actual_value) {
                    return Some(diff);
                }
            }
            for key in actual.keys() {
                if !expected.contains_key(key) {
                    return Some(format!("{path}.{key} present only in actual output"));
                }
            }
            None
        }
        (Value::Array(expected), Value::Array(actual)) => {
            if expected.len() != actual.len() {
                return Some(format!(
                    "{path} length differs: expected {}, actual {}",
                    expected.len(),
                    actual.len()
                ));
            }
            for (index, (expected_value, actual_value)) in expected.iter().zip(actual).enumerate() {
                if let Some(diff) =
                    first_json_diff(&format!("{path}[{index}]"), expected_value, actual_value)
                {
                    return Some(diff);
                }
            }
            None
        }
        (Value::Number(expected), Value::Number(actual))
            if path.contains(".coord[")
                && expected
                    .as_f64()
                    .zip(actual.as_f64())
                    .map(|(expected, actual)| (expected - actual).abs() <= 0.0015)
                    .unwrap_or(false) =>
        {
            None
        }
        _ if expected == actual => None,
        _ => Some(format!(
            "{path} differs: expected {}, actual {}",
            expected, actual
        )),
    }
}

fn normalize_for_comparison(value: &Value) -> Value {
    match value {
        Value::Array(items) => Value::Array(
            items
                .iter()
                .map(normalize_for_comparison)
                .collect::<Vec<_>>(),
        ),
        Value::Object(object) => {
            let mut normalized = serde_json::Map::new();
            for (key, value) in object {
                normalized.insert(key.clone(), normalize_for_comparison(value));
            }
            normalize_undirected_bond_object(&mut normalized);
            normalize_bond_array_object(&mut normalized);
            normalize_ring_set_object(&mut normalized);
            normalize_coord_object(&mut normalized);
            Value::Object(normalized)
        }
        _ => value.clone(),
    }
}

fn normalize_coord_object(object: &mut serde_json::Map<String, Value>) {
    let Some(Value::Array(coord)) = object.get_mut("coord") else {
        return;
    };
    for value in coord.iter_mut() {
        if let Some(number) = value.as_f64() {
            *value = json!((number * 1000.0).round() / 1000.0);
        }
    }
}

fn normalize_undirected_bond_object(object: &mut serde_json::Map<String, Value>) {
    let Some(begin) = object.get("begin_atom_index").and_then(Value::as_u64) else {
        return;
    };
    let Some(end) = object.get("end_atom_index").and_then(Value::as_u64) else {
        return;
    };
    if begin > end {
        object.insert("begin_atom_index".to_owned(), json!(end));
        object.insert("end_atom_index".to_owned(), json!(begin));
    }
}

fn normalize_bond_array_object(object: &mut serde_json::Map<String, Value>) {
    let Some(Value::Array(bonds)) = object.get_mut("bonds") else {
        return;
    };
    for bond in bonds.iter_mut() {
        if let Value::Object(bond) = bond {
            bond.remove("index");
        }
    }
    bonds.sort_by_key(bond_sort_key);
    for (index, bond) in bonds.iter_mut().enumerate() {
        if let Value::Object(bond) = bond {
            bond.insert("index".to_owned(), json!(index));
        }
    }
}

fn bond_sort_key(value: &Value) -> (u64, u64, String, String) {
    let Some(object) = value.as_object() else {
        return (u64::MAX, u64::MAX, String::new(), String::new());
    };
    (
        object
            .get("begin_atom_index")
            .and_then(Value::as_u64)
            .unwrap_or(u64::MAX),
        object
            .get("end_atom_index")
            .and_then(Value::as_u64)
            .unwrap_or(u64::MAX),
        object
            .get("bond_type")
            .and_then(Value::as_str)
            .unwrap_or_default()
            .to_owned(),
        object
            .get("stereo")
            .and_then(Value::as_str)
            .unwrap_or_default()
            .to_owned(),
    )
}

fn normalize_ring_set_object(object: &mut serde_json::Map<String, Value>) {
    let Some(Value::Array(rings)) = object.get_mut("rings") else {
        return;
    };
    for ring in rings.iter_mut() {
        let Value::Array(atoms) = ring else {
            continue;
        };
        atoms.sort_by_key(|value| value.as_u64().unwrap_or(u64::MAX));
    }
    rings.sort_by(|left, right| {
        let left = left
            .as_array()
            .map(|items| items.iter().filter_map(Value::as_u64).collect::<Vec<_>>())
            .unwrap_or_default();
        let right = right
            .as_array()
            .map(|items| items.iter().filter_map(Value::as_u64).collect::<Vec<_>>())
            .unwrap_or_default();
        left.cmp(&right)
    });
}

fn slugify_fixture(fixture: &str) -> String {
    let mut slug = String::new();
    let mut previous_was_separator = false;
    for ch in fixture.chars() {
        if ch.is_ascii_alphanumeric() || matches!(ch, '.' | '-' | '_') {
            slug.push(ch);
            previous_was_separator = false;
        } else if !previous_was_separator {
            slug.push('_');
            previous_was_separator = true;
        }
    }
    slug.trim_matches(['.', '_', '-']).to_owned()
}

fn is_sha256(value: &str) -> bool {
    value.len() == 64 && value.bytes().all(|byte| byte.is_ascii_hexdigit())
}

fn read_validation_statuses(
    features: &[Feature],
) -> Result<BTreeMap<String, ValidationStatus>, Box<dyn Error>> {
    let mut statuses = BTreeMap::new();
    for (corpus, _) in VALIDATION_CORPORA {
        let path = validation_status_path(corpus);
        if path.exists() {
            let status = read_corpus_status(&path)?;
            if status.corpus_id != *corpus {
                return Err(boxed_error(format!(
                    "{} declares corpus_id `{}`, expected `{corpus}`",
                    path.display(),
                    status.corpus_id
                )));
            }
            for (feature_id, feature_status) in status.features {
                if !features.iter().any(|feature| feature.id == feature_id) {
                    return Err(boxed_error(format!(
                        "{} records unknown feature `{feature_id}`",
                        path.display()
                    )));
                }
                statuses
                    .entry(feature_id.clone())
                    .or_insert_with(|| ValidationStatus::new(&feature_id))
                    .corpora
                    .insert((*corpus).to_owned(), feature_status);
            }
        }
    }
    Ok(statuses)
}

fn validation_status_path(corpus: &str) -> PathBuf {
    Path::new("validation")
        .join("corpora")
        .join(corpus)
        .join("status.toml")
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct CorpusStatusFile {
    corpus_id: String,
    #[serde(default)]
    features: BTreeMap<String, CorpusStatus>,
}

fn read_corpus_status(path: &Path) -> Result<CorpusStatusFile, Box<dyn Error>> {
    let text = fs::read_to_string(path)?;
    toml::from_str(&text).map_err(|error| boxed_error(format!("{}: {error}", path.display())))
}

fn write_validation_statuses(
    statuses: &BTreeMap<String, ValidationStatus>,
    selected_corpora: &BTreeSet<String>,
) -> Result<(), Box<dyn Error>> {
    for (corpus, _) in VALIDATION_CORPORA {
        if !selected_corpora.contains(*corpus) {
            continue;
        }
        let mut corpus_status = CorpusStatusFile {
            corpus_id: (*corpus).to_owned(),
            features: BTreeMap::new(),
        };
        for (feature_id, status) in statuses {
            if let Some(feature_status) = status.corpora.get(*corpus) {
                corpus_status
                    .features
                    .insert(feature_id.clone(), feature_status.clone());
            }
        }
        if corpus_status.features.is_empty() {
            continue;
        }
        let path = validation_status_path(corpus);
        fs::create_dir_all(
            path.parent()
                .ok_or_else(|| boxed_error("status path has no parent"))?,
        )?;
        let text = toml::to_string_pretty(&corpus_status)?;
        write_atomic_text(
            &path,
            &format!("# Generated by `cargo xtask validate --update`. Do not hand-edit.\n{text}"),
        )?;
    }
    Ok(())
}

fn corpus_passed(feature: &Feature, status: Option<&ValidationStatus>, corpus: &str) -> bool {
    corpus_passed_at(feature, status, corpus, Path::new("validation"))
}

fn corpus_passed_at(
    feature: &Feature,
    status: Option<&ValidationStatus>,
    corpus: &str,
    validation_root: &Path,
) -> bool {
    if !feature
        .validation_required
        .iter()
        .any(|item| item == corpus)
    {
        return false;
    }
    let Some(corpus_status) = status.and_then(|status| status.corpora.get(corpus)) else {
        return false;
    };
    if !corpus_status.passed {
        return false;
    }
    if corpus_status.fixture_count == 0
        || corpus_status.compared_count != corpus_status.fixture_count
    {
        return false;
    }
    if corpus_status.evidence_schema_version != Some(VALIDATION_EVIDENCE_SCHEMA_VERSION) {
        return false;
    }
    let Some(recorded_evidence_hash) = corpus_status.evidence_hash.as_deref() else {
        return false;
    };
    let manifest_path = validation_root
        .join("corpora")
        .join(corpus)
        .join("features")
        .join(format!("{}.toml", feature.id));
    if !manifest_path.exists()
        || hash_file(&manifest_path)
            .map(|hash| hash != corpus_status.manifest_hash)
            .unwrap_or(true)
    {
        return false;
    }
    let repo_root = validation_root.parent().unwrap_or_else(|| Path::new("."));
    read_validation_manifest(&manifest_path)
        .and_then(|manifest| build_validation_evidence(repo_root, &manifest_path, &manifest))
        .map(|evidence| {
            evidence.schema_version == VALIDATION_EVIDENCE_SCHEMA_VERSION
                && evidence.sha256 == recorded_evidence_hash
                && !evidence.inputs.is_empty()
        })
        .unwrap_or(false)
}

fn overall_validated(feature: &Feature, status: Option<&ValidationStatus>) -> bool {
    overall_validated_at(feature, status, Path::new("validation"))
}

fn overall_validated_at(
    feature: &Feature,
    status: Option<&ValidationStatus>,
    validation_root: &Path,
) -> bool {
    feature.implemented
        && !feature.validation_required.is_empty()
        && feature
            .validation_required
            .iter()
            .all(|corpus| corpus_passed_at(feature, status, corpus, validation_root))
}

fn sync_feature_validation_flags(
    features: &[Feature],
    statuses: &BTreeMap<String, ValidationStatus>,
) -> Result<(), Box<dyn Error>> {
    sync_feature_validation_flags_at(
        features,
        statuses,
        Path::new("features"),
        Path::new("validation"),
    )
}

fn sync_feature_validation_flags_at(
    features: &[Feature],
    statuses: &BTreeMap<String, ValidationStatus>,
    features_root: &Path,
    validation_root: &Path,
) -> Result<(), Box<dyn Error>> {
    for feature in features {
        let validated = overall_validated_at(feature, statuses.get(&feature.id), validation_root);
        if validated == feature.validated {
            continue;
        }
        let path = features_root.join(&feature.id).join("feature.toml");
        let text = fs::read_to_string(&path)?;
        let replacement = format!("validated = {validated}");
        let mut replaced = false;
        let rewritten = text
            .lines()
            .map(|line| {
                if line.trim_start().starts_with("validated =") {
                    replaced = true;
                    replacement.as_str()
                } else {
                    line
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        if !replaced {
            return Err(boxed_error(format!(
                "{} is missing `validated`",
                path.display()
            )));
        }
        write_atomic_text(&path, &format!("{rewritten}\n"))?;
    }
    Ok(())
}

fn write_atomic_text(path: &Path, text: &str) -> Result<(), Box<dyn Error>> {
    let parent = path
        .parent()
        .ok_or_else(|| boxed_error(format!("{} has no parent", path.display())))?;
    fs::create_dir_all(parent)?;
    let file_name = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| boxed_error(format!("{} has no file name", path.display())))?;
    let nonce = SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos();
    let tmp = parent.join(format!(".{file_name}.{nonce}.tmp"));
    fs::write(&tmp, text)?;
    match fs::rename(&tmp, path) {
        Ok(()) => Ok(()),
        Err(error) if path.exists() => {
            fs::remove_file(path)?;
            fs::rename(&tmp, path).map_err(|rename_error| {
                boxed_error(format!(
                    "failed to replace {} after initial rename error {error}: {rename_error}",
                    path.display()
                ))
            })
        }
        Err(error) => Err(boxed_error(format!(
            "failed to replace {}: {error}",
            path.display()
        ))),
    }
}

fn ensure_validation_flags_synced(
    features: &[Feature],
    statuses: &BTreeMap<String, ValidationStatus>,
) -> Result<(), Box<dyn Error>> {
    for feature in features {
        let derived = overall_validated(feature, statuses.get(&feature.id));
        if feature.validated != derived {
            return Err(boxed_error(format!(
                "feature `{}` has validated={}, but corpus evidence derives validated={derived}; run validation with --update",
                feature.id, feature.validated
            )));
        }
    }
    Ok(())
}

fn render_dashboard(features: &[Feature], statuses: &BTreeMap<String, ValidationStatus>) -> String {
    let mut out = String::new();
    out.push_str("<!doctype html>\n");
    out.push_str("<html lang=\"en\">\n");
    out.push_str("<head>\n");
    out.push_str("<meta charset=\"utf-8\">\n");
    out.push_str("<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n");
    out.push_str("<title>Feature Dashboard</title>\n");
    out.push_str("<style>\n");
    out.push_str(":root { color-scheme: light dark; --border: #d0d7de; --head: #f6f8fa; --ok: #1a7f37; --bad: #cf222e; --muted: #656d76; --text: #24292f; --bg: #ffffff; }\n");
    out.push_str("@media (prefers-color-scheme: dark) { :root { --border: #30363d; --head: #161b22; --ok: #3fb950; --bad: #ff7b72; --muted: #8b949e; --text: #c9d1d9; --bg: #0d1117; } }\n");
    out.push_str("body { margin: 24px; background: var(--bg); color: var(--text); font: 14px/1.4 system-ui, -apple-system, Segoe UI, sans-serif; }\n");
    out.push_str("h1 { margin: 0 0 4px; font-size: 24px; }\n");
    out.push_str("p { margin: 0 0 18px; color: var(--muted); }\n");
    out.push_str(".dashboard-wrap { overflow-x: auto; }\n");
    out.push_str("table { border-collapse: collapse; width: 100%; min-width: 980px; }\n");
    out.push_str(
        "th, td { border: 1px solid var(--border); padding: 6px 8px; vertical-align: middle; }\n",
    );
    out.push_str("thead th { position: sticky; top: 0; z-index: 1; height: 112px; background: var(--head); white-space: nowrap; }\n");
    out.push_str("tbody tr:nth-child(even) { background: color-mix(in srgb, var(--head) 45%, transparent); }\n");
    out.push_str("th.text, td.text { text-align: left; }\n");
    out.push_str("th.compact, td.compact, th.rotated, td.marker { text-align: center; }\n");
    out.push_str(
        "th.rotated { width: 42px; min-width: 42px; padding: 0; vertical-align: bottom; }\n",
    );
    out.push_str("th.rotated button { height: 108px; width: 42px; padding: 0; display: flex; align-items: flex-end; justify-content: center; }\n");
    out.push_str("th.rotated span { display: inline-block; transform: rotate(-60deg); transform-origin: bottom left; width: 96px; text-align: left; }\n");
    out.push_str(
        "button.sort { all: unset; cursor: pointer; color: inherit; font-weight: 650; }\n",
    );
    out.push_str(
        "button.sort:focus-visible { outline: 2px solid Highlight; outline-offset: 2px; }\n",
    );
    out.push_str("th[aria-sort=\"ascending\"] button.sort::after { content: \" \\25B2\"; font-size: 10px; color: var(--muted); }\n");
    out.push_str("th[aria-sort=\"descending\"] button.sort::after { content: \" \\25BC\"; font-size: 10px; color: var(--muted); }\n");
    out.push_str(".ok { color: var(--ok); font-weight: 700; }\n");
    out.push_str(".bad { color: var(--bad); font-weight: 700; }\n");
    out.push_str(".na { color: var(--muted); }\n");
    out.push_str("code { font-family: ui-monospace, SFMono-Regular, Consolas, monospace; font-size: 13px; }\n");
    out.push_str("</style>\n");
    out.push_str("</head>\n");
    out.push_str("<body>\n");
    out.push_str("<h1>Feature Dashboard</h1>\n");
    out.push_str("<p>Generated from feature metadata and validation status. Do not hand-edit this file.</p>\n");
    out.push_str("<div class=\"dashboard-wrap\">\n");
    out.push_str("<table id=\"feature-dashboard\">\n");
    out.push_str("<thead>\n<tr>");
    out.push_str("<th class=\"text\" data-sort-type=\"text\"><button class=\"sort\" type=\"button\">Feature</button></th>");
    out.push_str("<th class=\"text\" data-sort-type=\"text\"><button class=\"sort\" type=\"button\">Title</button></th>");
    out.push_str("<th class=\"compact rotated\" data-sort-type=\"text\"><button class=\"sort\" type=\"button\"><span>Area</span></button></th>");
    out.push_str("<th class=\"compact rotated\" data-sort-type=\"number\"><button class=\"sort\" type=\"button\"><span>Version</span></button></th>");
    out.push_str("<th class=\"compact rotated\" data-sort-type=\"number\"><button class=\"sort\" type=\"button\"><span>Implemented</span></button></th>");
    for (_, label) in VALIDATION_CORPORA {
        out.push_str(&format!(
            "<th class=\"rotated\" data-sort-type=\"number\"><button class=\"sort\" type=\"button\"><span>{}</span></button></th>",
            escape_html(label)
        ));
    }
    out.push_str("</tr>\n</thead>\n<tbody>\n");
    for feature in features {
        let status = statuses.get(&feature.id);
        out.push_str(&format!(
            "<tr><td class=\"text\" data-sort-value=\"{0}\"><code>{0}</code></td>",
            escape_html(&feature.id)
        ));
        out.push_str(&format!(
            "<td class=\"text\" data-sort-value=\"{}\">{}</td>",
            escape_html(&feature.title),
            escape_html(&feature.title)
        ));
        out.push_str(&format!(
            "<td class=\"compact\" data-sort-value=\"{}\">{}</td>",
            escape_html(&feature.area),
            escape_html(&feature.area)
        ));
        out.push_str(&format!(
            "<td class=\"compact\" data-sort-value=\"{}\">{}</td>",
            feature.version, feature.version
        ));
        out.push_str(&format!(
            "<td class=\"marker\" data-sort-value=\"{}\">{}</td>",
            bool_sort_value(feature.implemented),
            dashboard_marker(Some(feature.implemented))
        ));
        for (corpus, _) in VALIDATION_CORPORA {
            let marker = if feature
                .validation_required
                .iter()
                .any(|required| required == corpus)
            {
                Some(corpus_passed(feature, status, corpus))
            } else {
                None
            };
            out.push_str(&format!(
                "<td class=\"marker\" data-sort-value=\"{}\">{}</td>",
                optional_bool_sort_value(marker),
                dashboard_marker(marker)
            ));
        }
        out.push_str("</tr>\n");
    }
    out.push_str("</tbody>\n</table>\n</div>\n");
    out.push_str("<script>\n");
    out.push_str("(() => {\n");
    out.push_str("  const table = document.getElementById('feature-dashboard');\n");
    out.push_str("  const tbody = table.tBodies[0];\n");
    out.push_str("  const headers = Array.from(table.tHead.rows[0].cells);\n");
    out.push_str("  const value = (row, index, type) => {\n");
    out.push_str("    const raw = row.cells[index].dataset.sortValue || row.cells[index].textContent.trim();\n");
    out.push_str("    return type === 'number' ? Number(raw) : raw.toLocaleLowerCase();\n");
    out.push_str("  };\n");
    out.push_str("  headers.forEach((header, index) => {\n");
    out.push_str("    const button = header.querySelector('button.sort');\n");
    out.push_str("    if (!button) return;\n");
    out.push_str("    button.addEventListener('click', () => {\n");
    out.push_str("      const ascending = header.getAttribute('aria-sort') !== 'ascending';\n");
    out.push_str("      headers.forEach(other => other.removeAttribute('aria-sort'));\n");
    out.push_str(
        "      header.setAttribute('aria-sort', ascending ? 'ascending' : 'descending');\n",
    );
    out.push_str("      const type = header.dataset.sortType || 'text';\n");
    out.push_str("      const rows = Array.from(tbody.rows);\n");
    out.push_str("      rows.sort((left, right) => {\n");
    out.push_str("        const a = value(left, index, type);\n");
    out.push_str("        const b = value(right, index, type);\n");
    out.push_str("        if (a < b) return ascending ? -1 : 1;\n");
    out.push_str("        if (a > b) return ascending ? 1 : -1;\n");
    out.push_str("        return value(left, 0, 'text').localeCompare(value(right, 0, 'text'));\n");
    out.push_str("      });\n");
    out.push_str("      rows.forEach(row => tbody.appendChild(row));\n");
    out.push_str("    });\n");
    out.push_str("  });\n");
    out.push_str("})();\n");
    out.push_str("</script>\n");
    out.push_str("</body>\n</html>\n");
    out
}

fn dashboard_marker(value: Option<bool>) -> &'static str {
    match value {
        Some(true) => "<span class=\"ok\" aria-label=\"passed\">&#10003;</span>",
        Some(false) => "<span class=\"bad\" aria-label=\"failed\">&#10007;</span>",
        None => "<span class=\"na\" aria-label=\"not required\">-</span>",
    }
}

fn bool_sort_value(value: bool) -> &'static str {
    if value {
        "1"
    } else {
        "0"
    }
}

fn optional_bool_sort_value(value: Option<bool>) -> &'static str {
    match value {
        Some(true) => "1",
        Some(false) => "0",
        None => "-1",
    }
}

fn escape_html(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
#[derive(Debug, Clone, PartialEq, Eq)]
struct SkillMetadata {
    name: String,
    description: String,
}

fn check_skills(root: &Path) -> Result<(), Box<dyn Error>> {
    for expected in expected_skills() {
        let path = root.join(expected.name).join("SKILL.md");
        if !path.exists() {
            return Err(boxed_error(format!(
                "missing repo-local skill `{}` at {}",
                expected.name,
                path.display()
            )));
        }
        let text = fs::read_to_string(&path)?;
        let metadata = parse_skill_metadata(&text, &path)?;
        if metadata.name != expected.name {
            return Err(boxed_error(format!(
                "{} declares skill name `{}`, expected `{}`",
                path.display(),
                metadata.name,
                expected.name
            )));
        }
        let lower = text.to_lowercase();
        for required in expected.required_phrases {
            if !lower.contains(&required.to_lowercase()) {
                return Err(boxed_error(format!(
                    "{} is missing required phrase `{required}`",
                    path.display()
                )));
            }
        }
    }
    Ok(())
}

struct ExpectedSkill {
    name: &'static str,
    required_phrases: &'static [&'static str],
}

fn expected_skills() -> &'static [ExpectedSkill] {
    &[
        ExpectedSkill {
            name: "feature-work",
            required_phrases: &[
                "add -> optional research -> plan -> implement",
                "feature.md",
                "implemented = true",
                "validation_required",
                "externally supplied",
                "cargo xtask dashboard --check",
                "cargo xtask validate --feature",
                "--corpus",
            ],
        },
        ExpectedSkill {
            name: "feature-review",
            required_phrases: &[
                "independent audit",
                "architecture",
                "validation claims",
                "feature.md",
                "cargo test --workspace",
                "cargo xtask validate --feature",
                "--corpus",
            ],
        },
    ]
}

fn parse_skill_metadata(text: &str, path: &Path) -> Result<SkillMetadata, Box<dyn Error>> {
    let mut lines = text.lines();
    if lines.next() != Some("---") {
        return Err(boxed_error(format!(
            "{} is missing YAML frontmatter",
            path.display()
        )));
    }
    let mut fields = BTreeMap::<String, String>::new();
    for line in lines.by_ref() {
        if line == "---" {
            let name = fields
                .get("name")
                .cloned()
                .ok_or_else(|| boxed_error(format!("{} is missing `name`", path.display())))?;
            let description = fields.get("description").cloned().ok_or_else(|| {
                boxed_error(format!("{} is missing `description`", path.display()))
            })?;
            if name.trim().is_empty() || description.trim().is_empty() {
                return Err(boxed_error(format!(
                    "{} has empty skill frontmatter",
                    path.display()
                )));
            }
            return Ok(SkillMetadata { name, description });
        }
        if let Some((key, value)) = line.split_once(':') {
            fields.insert(key.trim().to_owned(), value.trim().to_owned());
        }
    }
    Err(boxed_error(format!(
        "{} has unterminated YAML frontmatter",
        path.display()
    )))
}

fn boxed_error(message: impl Into<String>) -> Box<dyn Error> {
    std::io::Error::other(message.into()).into()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn value_after_flag_finds_following_value() {
        let args = vec![
            "validate".to_owned(),
            "--feature".to_owned(),
            "core.graph".to_owned(),
        ];

        assert_eq!(value_after_flag(&args, "--feature"), Some("core.graph"));
    }

    #[test]
    fn read_feature_parses_typed_metadata() {
        let root = temp_feature_root("read-feature");
        write_feature(
            &root,
            "example.feature",
            r#"id = "example.feature"
title = "Example"
area = "infrastructure"
version = 2
implemented = false
validated = true
description = "Example feature."
depends_on = ["core.graph"]
validation_required = []
"#,
        );

        let feature = read_feature(&root.join("example.feature").join("feature.toml"))
            .expect("feature should parse");

        assert_eq!(feature.id, "example.feature");
        assert_eq!(feature.version, 2);
        assert!(!feature.implemented);
        assert!(feature.validated);
        assert_eq!(feature.depends_on, vec!["core.graph"]);
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn read_feature_rejects_bad_boolean_deprecated_keys_missing_docs_and_directory_mismatch() {
        let root = temp_feature_root("bad-feature");
        write_feature(
            &root,
            "bad.bool",
            r#"id = "bad.bool"
title = "Bad"
area = "infrastructure"
version = 1
implemented = maybe
validated = false
description = "Bad feature."
depends_on = []
validation_required = []
"#,
        );
        assert!(read_feature(&root.join("bad.bool").join("feature.toml")).is_err());

        write_feature(
            &root,
            "bad.deprecated",
            r#"id = "bad.deprecated"
title = "Bad"
area = "infrastructure"
version = 1
priority = "P0"
implemented = false
validated = false
description = "Bad feature."
depends_on = []
validation_required = []
"#,
        );
        assert!(read_feature(&root.join("bad.deprecated").join("feature.toml")).is_err());

        write_feature(
            &root,
            "bad.version",
            r#"id = "bad.version"
title = "Bad"
area = "infrastructure"
version = 0
implemented = false
validated = false
description = "Bad feature."
depends_on = []
validation_required = []
"#,
        );
        assert!(read_feature(&root.join("bad.version").join("feature.toml")).is_err());

        write_feature_without_doc(
            &root,
            "missing.doc",
            r#"id = "missing.doc"
title = "Bad"
area = "infrastructure"
version = 1
implemented = false
validated = false
description = "Bad feature."
depends_on = []
validation_required = []
"#,
        );
        assert!(read_feature(&root.join("missing.doc").join("feature.toml")).is_err());

        write_feature(
            &root,
            "real.id",
            r#"id = "wrong.id"
title = "Bad"
area = "infrastructure"
version = 1
implemented = false
validated = false
description = "Bad feature."
depends_on = []
validation_required = []
"#,
        );
        assert!(read_feature(&root.join("real.id").join("feature.toml")).is_err());
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn read_features_sorts_skips_templates_and_validates_dependencies() {
        let root = temp_feature_root("feature-set");
        write_feature(
            &root,
            "z.feature",
            r#"id = "z.feature"
title = "Zed"
area = "core"
version = 1
implemented = true
validated = false
description = "Z feature."
depends_on = ["a.feature"]
validation_required = []
"#,
        );
        write_feature(
            &root,
            "a.feature",
            r#"id = "a.feature"
title = "Aye"
area = "core"
version = 1
implemented = false
validated = false
description = "A feature."
depends_on = []
validation_required = []
"#,
        );
        fs::create_dir_all(root.join("_template")).expect("template dir should create");
        fs::write(root.join("_template").join("feature.toml"), "not = valid")
            .expect("template metadata should write");

        let features = read_features_from(&root).expect("feature set should parse");

        assert_eq!(
            features
                .iter()
                .map(|feature| feature.id.as_str())
                .collect::<Vec<_>>(),
            vec!["a.feature", "z.feature"]
        );

        write_feature(
            &root,
            "bad.dependency",
            r#"id = "bad.dependency"
title = "Bad"
area = "core"
version = 1
implemented = false
validated = false
description = "Bad dependency."
depends_on = ["missing.feature"]
validation_required = []
"#,
        );
        assert!(read_features_from(&root).is_err());
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn render_dashboard_is_stable_and_uses_boolean_labels() {
        let features = vec![
            Feature {
                id: "a.feature".to_owned(),
                title: "Aye".to_owned(),
                area: "core".to_owned(),
                version: 1,
                implemented: false,
                validated: false,
                description: "A feature.".to_owned(),
                depends_on: Vec::new(),
                validation_required: Vec::new(),
            },
            Feature {
                id: "z.feature".to_owned(),
                title: "Zed".to_owned(),
                area: "io".to_owned(),
                version: 3,
                implemented: true,
                validated: false,
                description: "Z feature.".to_owned(),
                depends_on: vec!["a.feature".to_owned()],
                validation_required: Vec::new(),
            },
        ];

        let dashboard = render_dashboard(&features, &BTreeMap::new());

        assert!(dashboard.starts_with("<!doctype html>\n"));
        assert!(dashboard.contains("<table id=\"feature-dashboard\">"));
        assert!(dashboard.contains("<span>Implemented</span>"));
        assert!(!dashboard.contains("Validated"));
        assert!(dashboard.contains("<span>Tiny</span>"));
        assert!(dashboard.contains("<code>a.feature</code>"));
        assert!(dashboard.contains("data-sort-value=\"0\""));
        assert!(dashboard.contains("<code>z.feature</code>"));
        assert!(dashboard.contains("data-sort-value=\"1\""));
        assert!(dashboard.contains("button.addEventListener('click'"));
        assert!(dashboard.ends_with('\n'));
    }

    #[test]
    fn skill_metadata_parser_and_check_validate_repo_skill_contract() {
        let root = temp_feature_root("skills-check");
        write_skill(
            &root,
            "feature-work",
            r#"---
name: feature-work
description: Builder skill.
---
# Feature Work
add -> optional research -> plan -> implement
Use feature.md. Set implemented = true with evidence and declare validation_required.
Molecular validation fixtures must be externally supplied.
Run cargo xtask dashboard --check and cargo xtask validate --feature <feature-id> --corpus tiny.
"#,
        );
        write_skill(
            &root,
            "feature-review",
            r#"---
name: feature-review
description: Review skill.
---
# Feature Review
Independent audit for architecture and validation claims.
Read feature.md. Run cargo test --workspace and cargo xtask validate --feature <feature-id> --corpus tiny.
"#,
        );

        check_skills(&root).expect("skills should pass");

        fs::write(
            root.join("feature-review").join("SKILL.md"),
            "# Missing frontmatter",
        )
        .expect("skill should rewrite");
        assert!(check_skills(&root).is_err());
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn validation_manifest_path_is_feature_scoped() {
        assert_eq!(
            validation_manifest_path("core.graph", "tiny"),
            PathBuf::from("validation/corpora/tiny/features/core.graph.toml")
        );
    }

    #[test]
    fn all_selectors_expand_only_applicable_feature_corpus_pairs() {
        let features = vec![
            Feature {
                id: "small".to_owned(),
                title: "Small".to_owned(),
                area: "io".to_owned(),
                version: 1,
                implemented: true,
                validated: false,
                description: "Small feature.".to_owned(),
                depends_on: Vec::new(),
                validation_required: vec!["tiny".to_owned(), "pubchem-100".to_owned()],
            },
            Feature {
                id: "macro".to_owned(),
                title: "Macro".to_owned(),
                area: "bio".to_owned(),
                version: 1,
                implemented: true,
                validated: false,
                description: "Macro feature.".to_owned(),
                depends_on: Vec::new(),
                validation_required: vec!["tiny".to_owned(), "pdb-10".to_owned()],
            },
        ];

        assert_eq!(
            validation_targets(&features, "all", "pubchem-100")
                .into_iter()
                .map(|(feature, corpus)| (feature.id.as_str(), corpus))
                .collect::<Vec<_>>(),
            vec![("small", "pubchem-100".to_owned())]
        );
        assert_eq!(validation_targets(&features, "small", "all").len(), 2);
        assert_eq!(
            validation_targets(&features, "macro", "pubchem-100").len(),
            0
        );
    }

    #[test]
    fn current_status_drives_overall_validation_and_metadata_sync() {
        let root = temp_feature_root("status-sync");
        let (features_root, validation_root, manifest_path) = write_evidence_test_repo(&root);
        let metadata_path = features_root.join("example").join("feature.toml");

        let feature = Feature {
            id: "example".to_owned(),
            title: "Example".to_owned(),
            area: "io".to_owned(),
            version: 1,
            implemented: true,
            validated: false,
            description: "Example feature.".to_owned(),
            depends_on: Vec::new(),
            validation_required: vec!["tiny".to_owned()],
        };
        let manifest = read_validation_manifest(&manifest_path).expect("manifest should read");
        let evidence = build_validation_evidence(&root, &manifest_path, &manifest)
            .expect("evidence should build");
        let corpus_status = CorpusStatus {
            passed: true,
            fixture_count: 1,
            compared_count: 1,
            reference_tool: "rdkit".to_owned(),
            reference_version: "RDKit test".to_owned(),
            manifest_hash: hash_file(&manifest_path).expect("manifest should hash"),
            evidence_schema_version: Some(VALIDATION_EVIDENCE_SCHEMA_VERSION),
            evidence_hash: Some(evidence.sha256),
            evidence_inputs: evidence.inputs,
            validated_at_unix: 1,
        };
        let status = ValidationStatus {
            feature_id: feature.id.clone(),
            corpora: BTreeMap::from([("tiny".to_owned(), corpus_status)]),
        };
        let statuses = BTreeMap::from([(feature.id.clone(), status.clone())]);

        assert!(overall_validated_at(
            &feature,
            Some(&status),
            &validation_root
        ));
        sync_feature_validation_flags_at(
            std::slice::from_ref(&feature),
            &statuses,
            &features_root,
            &validation_root,
        )
        .expect("metadata should sync");
        assert!(fs::read_to_string(&metadata_path)
            .expect("metadata should read")
            .contains("validated = true"));

        fs::write(&manifest_path, "changed = true\n").expect("manifest should change");
        assert!(!overall_validated_at(
            &feature,
            Some(&status),
            &validation_root
        ));
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn evidence_changes_after_material_input_changes() {
        let root = temp_feature_root("evidence-change");
        let (_, _, manifest_path) = write_evidence_test_repo(&root);
        let manifest = read_validation_manifest(&manifest_path).expect("manifest should read");
        let original = build_validation_evidence(&root, &manifest_path, &manifest)
            .expect("evidence should build");

        fs::write(root.join("crates/molecules/src/lib.rs"), "changed source\n")
            .expect("source should mutate");
        let source_changed = build_validation_evidence(&root, &manifest_path, &manifest)
            .expect("evidence should build");
        assert_ne!(original.sha256, source_changed.sha256);

        fs::write(
            root.join("validation/corpora/tiny/data/example.sdf"),
            "changed fixture\n",
        )
        .expect("fixture should mutate");
        let fixture_changed = build_validation_evidence(&root, &manifest_path, &manifest)
            .expect("evidence should build");
        assert_ne!(source_changed.sha256, fixture_changed.sha256);

        fs::write(
            root.join("validation/corpora/tiny/golden/example/data_example.sdf.json.gz"),
            "changed golden\n",
        )
        .expect("golden should mutate");
        let golden_changed = build_validation_evidence(&root, &manifest_path, &manifest)
            .expect("evidence should build");
        assert_ne!(fixture_changed.sha256, golden_changed.sha256);
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn evidence_hash_normalizes_text_line_endings() {
        let root = temp_feature_root("evidence-line-endings");
        let path = root.join("source.rs");
        fs::write(&path, "first\nsecond\n").expect("LF source should write");
        let lf_hash = hash_evidence_file(&path).expect("LF evidence should hash");

        fs::write(&path, "first\r\nsecond\r\n").expect("CRLF source should write");
        let crlf_hash = hash_evidence_file(&path).expect("CRLF evidence should hash");

        assert_eq!(lf_hash, crlf_hash);
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn current_status_requires_known_nonempty_evidence() {
        let root = temp_feature_root("status-rejects-stale");
        let (_, validation_root, manifest_path) = write_evidence_test_repo(&root);
        let manifest = read_validation_manifest(&manifest_path).expect("manifest should read");
        let evidence = build_validation_evidence(&root, &manifest_path, &manifest)
            .expect("evidence should build");
        let feature = Feature {
            id: "example".to_owned(),
            title: "Example".to_owned(),
            area: "io".to_owned(),
            version: 1,
            implemented: true,
            validated: false,
            description: "Example feature.".to_owned(),
            depends_on: Vec::new(),
            validation_required: vec!["tiny".to_owned()],
        };
        let mut corpus_status = CorpusStatus {
            passed: true,
            fixture_count: 1,
            compared_count: 1,
            reference_tool: "rdkit".to_owned(),
            reference_version: "RDKit test".to_owned(),
            manifest_hash: hash_file(&manifest_path).expect("manifest should hash"),
            evidence_schema_version: Some(VALIDATION_EVIDENCE_SCHEMA_VERSION),
            evidence_hash: Some(evidence.sha256),
            evidence_inputs: evidence.inputs,
            validated_at_unix: 1,
        };
        let status = ValidationStatus {
            feature_id: feature.id.clone(),
            corpora: BTreeMap::from([("tiny".to_owned(), corpus_status.clone())]),
        };
        assert!(overall_validated_at(
            &feature,
            Some(&status),
            &validation_root
        ));

        corpus_status.evidence_schema_version = Some(999);
        let status = ValidationStatus {
            feature_id: feature.id.clone(),
            corpora: BTreeMap::from([("tiny".to_owned(), corpus_status.clone())]),
        };
        assert!(!overall_validated_at(
            &feature,
            Some(&status),
            &validation_root
        ));

        corpus_status.evidence_schema_version = Some(VALIDATION_EVIDENCE_SCHEMA_VERSION);
        corpus_status.compared_count = 0;
        let status = ValidationStatus {
            feature_id: feature.id.clone(),
            corpora: BTreeMap::from([("tiny".to_owned(), corpus_status)]),
        };
        assert!(!overall_validated_at(
            &feature,
            Some(&status),
            &validation_root
        ));
        fs::remove_dir_all(root).ok();
    }

    #[test]
    fn unsupported_comparison_mode_is_rejected() {
        let manifest = ValidationManifest {
            feature_id: "example".to_owned(),
            corpus_id: "tiny".to_owned(),
            reference_tool: "rdkit".to_owned(),
            reference_version: "RDKit test".to_owned(),
            comparison_mode: "planned".to_owned(),
            fixtures: vec!["data/example.sdf".to_owned()],
            _notes: Vec::new(),
        };
        assert!(validate_comparison_mode(Path::new("example.toml"), &manifest).is_err());
    }

    #[test]
    fn comparison_normalizes_undirected_bonds_and_ring_order() {
        let expected = json!({
            "records": [{
                "bonds": [
                    {"index": 0, "begin_atom_index": 5, "end_atom_index": 0, "bond_type": "SINGLE", "stereo": "STEREONONE"}
                ],
                "rings": [[5, 3, 1]]
            }]
        });
        let actual = json!({
            "records": [{
                "bonds": [
                    {"index": 7, "begin_atom_index": 0, "end_atom_index": 5, "bond_type": "SINGLE", "stereo": "STEREONONE"}
                ],
                "rings": [[1, 3, 5]]
            }]
        });

        assert_eq!(
            normalize_for_comparison(&expected),
            normalize_for_comparison(&actual)
        );
    }

    #[test]
    fn smiles_semantic_records_assert_topology_and_atom_identity() {
        let single = read_smiles_str("CC", SmilesParseOptions).expect("single bond should parse");
        let double = read_smiles_str("C=C", SmilesParseOptions).expect("double bond should parse");
        assert_ne!(
            smiles_sanitized_bonds_json(&single.mol),
            smiles_sanitized_bonds_json(&double.mol)
        );

        let aromatic =
            read_smiles_str("c1ccccc1", SmilesParseOptions).expect("benzene should parse");
        assert!(smiles_sanitized_bonds_json(&aromatic.mol)
            .iter()
            .all(|bond| bond["bond_type"] == "AROMATIC" && bond["is_aromatic"] == true));

        let labeled =
            read_smiles_str("[13CH3:7]C", SmilesParseOptions).expect("labeled carbon should parse");
        let atoms = smiles_sanitized_atoms_json(&labeled.mol);
        assert!(atoms
            .iter()
            .any(|atom| atom["isotope"] == 13 && atom["atom_map"] == 7));
        assert!(atoms.iter().all(|atom| atom["neighbors"].is_array()));
    }

    fn temp_feature_root(label: &str) -> PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("time should be available")
            .as_nanos();
        let root =
            env::temp_dir().join(format!("molecules-xtask-{label}-{}-{nonce}", process::id()));
        fs::create_dir_all(&root).expect("temp feature root should create");
        root
    }

    fn write_evidence_test_repo(root: &Path) -> (PathBuf, PathBuf, PathBuf) {
        let features_root = root.join("features");
        let validation_root = root.join("validation");
        let feature_dir = features_root.join("example");
        let corpus_root = validation_root.join("corpora").join("tiny");
        let manifest_dir = corpus_root.join("features");
        fs::create_dir_all(&feature_dir).expect("feature dir should create");
        fs::create_dir_all(&manifest_dir).expect("manifest dir should create");
        fs::write(
            feature_dir.join("feature.toml"),
            "id = \"example\"\nvalidated = false\n",
        )
        .expect("metadata should write");
        fs::write(feature_dir.join("feature.md"), "# Example\n").expect("feature doc should write");
        let manifest_path = manifest_dir.join("example.toml");
        fs::write(
            &manifest_path,
            "feature_id = \"example\"\ncorpus_id = \"tiny\"\nreference_tool = \"rdkit\"\nreference_version = \"RDKit test\"\ncomparison_mode = \"implementation-golden\"\nfixtures = [\"data/example.sdf\"]\n",
        )
        .expect("manifest should write");
        fs::create_dir_all(corpus_root.join("data")).expect("data dir should create");
        fs::create_dir_all(corpus_root.join("golden").join("example"))
            .expect("golden dir should create");
        fs::write(corpus_root.join("corpus.toml"), "id = \"tiny\"\n")
            .expect("corpus descriptor should write");
        fs::write(corpus_root.join("sources.lock.json"), "{}\n").expect("source lock should write");
        fs::write(corpus_root.join("data").join("example.sdf"), "fixture\n")
            .expect("fixture should write");
        fs::write(
            corpus_root
                .join("golden")
                .join("example")
                .join("data_example.sdf.json.gz"),
            "golden\n",
        )
        .expect("golden should write");
        fs::write(root.join("Cargo.toml"), "[workspace]\n").expect("cargo toml should write");
        fs::write(root.join("Cargo.lock"), "# lock\n").expect("cargo lock should write");
        for path in [
            "crates/molecules/Cargo.toml",
            "crates/xtask/Cargo.toml",
            "crates/molecules/src/lib.rs",
            "crates/xtask/src/main.rs",
            "validation/reference/rdkit/run_feature.py",
            "validation/reference/rdkit/environment.yml",
        ] {
            let path = root.join(path);
            fs::create_dir_all(path.parent().expect("test path should have parent"))
                .expect("test parent should create");
            fs::write(path, "test\n").expect("test evidence file should write");
        }
        (features_root, validation_root, manifest_path)
    }

    fn write_feature(root: &Path, id: &str, metadata: &str) {
        let dir = root.join(id);
        fs::create_dir_all(&dir).expect("feature dir should create");
        fs::write(dir.join("feature.toml"), metadata).expect("feature metadata should write");
        fs::write(dir.join("feature.md"), "# Feature\n").expect("feature doc should write");
    }

    fn write_feature_without_doc(root: &Path, id: &str, metadata: &str) {
        let dir = root.join(id);
        fs::create_dir_all(&dir).expect("feature dir should create");
        fs::write(dir.join("feature.toml"), metadata).expect("feature metadata should write");
    }

    fn write_skill(root: &Path, name: &str, text: &str) {
        let dir = root.join(name);
        fs::create_dir_all(&dir).expect("skill dir should create");
        fs::write(dir.join("SKILL.md"), text).expect("skill should write");
    }
}
