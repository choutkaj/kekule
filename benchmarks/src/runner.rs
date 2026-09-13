use crate::*;
use indicatif::{ProgressBar, ProgressFinish, ProgressStyle};
use rayon::prelude::*;
use std::hint::black_box;
use std::time::{Instant, SystemTime, UNIX_EPOCH};

const FEATURES: &[&str] = &[
    "io.smiles.parse",
    "io.smiles.write",
    "io.smiles.canonical",
    "io.smiles.isomeric",
    "io.mol.v2000.parse",
    "io.mol.v2000.write",
    "io.mol.v3000.parse",
    "io.mol.v3000.write",
    "io.sdf.v2000.parse",
    "io.sdf.v2000.write",
    "io.mmcif.parse",
    "algo.rings.fast",
    "algo.rings.sssr",
    "algo.valence.rdkit-like",
    "algo.aromaticity.rdkit-like",
    "algo.canonical-ranking",
    "algo.substructure.vf2",
    "query.smarts",
    "chem.perception.default",
    "chem.hydrogen-transforms",
    "descriptor.molecular",
    "descriptor.rotatable-bonds.rdkit-strict",
    "stereo.representation",
    "stereo.perception",
    "stereo.cip",
    "bio.secondary-structure.dssp",
];

struct Options {
    feature: String,
    dataset: String,
    limit: Option<usize>,
    jobs: usize,
    output: PathBuf,
    reference_python: Option<PathBuf>,
}

fn options(args: &[String]) -> Result<Options, Box<dyn Error>> {
    let mut values = BTreeMap::new();
    if !args.len().is_multiple_of(2) {
        return Err(boxed_error("each option requires a value; use --help"));
    }
    for pair in args.chunks_exact(2) {
        if ![
            "--feature",
            "--dataset",
            "--limit",
            "--jobs",
            "--output",
            "--reference-python",
        ]
        .contains(&pair[0].as_str())
        {
            return Err(boxed_error(format!("unknown option: {}", pair[0])));
        }
        if values.insert(pair[0].as_str(), pair[1].as_str()).is_some() {
            return Err(boxed_error(format!("duplicate option: {}", pair[0])));
        }
    }
    let positive = |flag| -> Result<Option<usize>, Box<dyn Error>> {
        values
            .get(flag)
            .map(|value| {
                value
                    .parse::<usize>()
                    .ok()
                    .filter(|n| *n > 0)
                    .ok_or_else(|| boxed_error(format!("{flag} must be a positive integer")))
            })
            .transpose()
    };
    Ok(Options {
        feature: values
            .get("--feature")
            .ok_or("--feature is required")?
            .to_string(),
        dataset: values
            .get("--dataset")
            .ok_or("--dataset is required")?
            .to_string(),
        limit: positive("--limit")?,
        jobs: positive("--jobs")?
            .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, usize::from)),
        output: values
            .get("--output")
            .map(PathBuf::from)
            .unwrap_or_else(|| {
                Path::new(env!("CARGO_MANIFEST_DIR"))
                    .join("../target/benchmarks")
                    .join(format!(
                        "run-{}-{}.json",
                        process::id(),
                        SystemTime::now()
                            .duration_since(UNIX_EPOCH)
                            .unwrap_or_default()
                            .as_nanos()
                    ))
            }),
        reference_python: values.get("--reference-python").map(PathBuf::from),
    })
}

#[derive(Default, Serialize)]
struct Counts {
    selected: usize,
    compared: usize,
    agrees: usize,
    disagrees: usize,
    unsupported: usize,
    errors: usize,
}

#[derive(Default, Serialize)]
struct NumericErrors {
    values: usize,
    mean_absolute_error: f64,
    maximum_absolute_error: f64,
}

fn numeric_errors(
    path: &str,
    expected: &Value,
    actual: &Value,
    metrics: &mut BTreeMap<String, NumericErrors>,
) {
    match (expected, actual) {
        (Value::Object(e), Value::Object(a)) => {
            for (key, value) in e {
                if let Some(other) = a.get(key) {
                    numeric_errors(&format!("{path}.{key}"), value, other, metrics);
                }
            }
        }
        (Value::Array(e), Value::Array(a)) => {
            for (value, other) in e.iter().zip(a) {
                numeric_errors(&format!("{path}[]"), value, other, metrics);
            }
        }
        (Value::Number(e), Value::Number(a)) if e.is_f64() || a.is_f64() => {
            if let Some((e, a)) = e.as_f64().zip(a.as_f64()) {
                let error = (e - a).abs();
                let metric = metrics.entry(path.to_owned()).or_default();
                metric.values += 1;
                metric.mean_absolute_error +=
                    (error - metric.mean_absolute_error) / metric.values as f64;
                metric.maximum_absolute_error = metric.maximum_absolute_error.max(error);
            }
        }
        _ => (),
    }
}

#[derive(Serialize)]
struct ResultRow {
    dataset: String,
    feature: String,
    reference_tool: String,
    reference_version: String,
    reference_kind: &'static str,
    preparation: Vec<String>,
    source_lock_sha256: String,
    selection_id: Value,
    selection_sha256: String,
    counts: Counts,
    numeric_errors: BTreeMap<String, NumericErrors>,
    time_ms: Option<f64>,
    fixtures: Vec<Value>,
    differences: Vec<Value>,
    execution_errors: Vec<Value>,
    live_reference_disagrees: usize,
    live_references: Vec<Value>,
}

fn expected_record(expected: &Value, index: usize) -> Result<Value, Box<dyn Error>> {
    if let Some(records) = expected.get("records").and_then(Value::as_array) {
        let mut value = Value::Object(
            expected
                .as_object()
                .ok_or("reference output must be an object")?
                .iter()
                .filter(|(key, _)| key.as_str() != "records")
                .map(|(key, value)| (key.clone(), value.clone()))
                .collect(),
        );
        let selected = records
            .iter()
            .filter(|record| record["record_index"].as_u64() == Some(index as u64))
            .cloned()
            .collect::<Vec<_>>();
        if selected.len() > 1 {
            return Err(boxed_error("duplicate reference record index"));
        }
        value["records"] = json!(selected);
        Ok(value)
    } else {
        Ok(expected.clone())
    }
}

fn restore_index(value: &mut Value, index: usize) {
    if let Some(records) = value.get_mut("records").and_then(Value::as_array_mut) {
        for record in records {
            record["record_index"] = json!(index);
        }
    }
}

fn output_status(value: &Value) -> &str {
    if let Some(records) = value["records"].as_array() {
        if records.is_empty() {
            return "unsupported";
        }
        for record in records {
            let status = record["status"].as_str().unwrap_or("ok");
            if status != "ok" {
                return status;
            }
        }
    }
    nested_failure_status(value).unwrap_or("ok")
}

fn nested_failure_status(value: &Value) -> Option<&str> {
    match value {
        Value::Object(fields) => fields
            .get("status")
            .and_then(Value::as_str)
            .filter(|status| *status != "ok")
            .or_else(|| fields.values().find_map(nested_failure_status)),
        Value::Array(values) => values.iter().find_map(nested_failure_status),
        _ => None,
    }
}

fn workload_sizes(value: &Value, sizes: &mut BTreeMap<String, u64>) {
    if let Some(count) = value["atom_site_rows"]["row_count"].as_u64() {
        *sizes.entry("atom_site_rows".into()).or_default() += count;
    }
    if let Some(residues) = value["residues"].as_array() {
        *sizes.entry("compared_residues".into()).or_default() += residues.len() as u64;
    }
    for record in value["records"].as_array().into_iter().flatten() {
        if let Some(count) = record["atom_count"]
            .as_u64()
            .or_else(|| record["normalized_perceived"]["atom_count"].as_u64())
        {
            *sizes.entry("atoms".into()).or_default() += count;
        }
    }
}

fn compare_record(
    row: &mut ResultRow,
    id: &str,
    fixture: &str,
    index: usize,
    mut expected: Value,
    mut actual: Value,
) {
    restore_index(&mut actual, index);
    normalize_benchmark_for_comparison_in_place(&row.feature, &mut expected);
    normalize_benchmark_for_comparison_in_place(&row.feature, &mut actual);
    let difference = first_json_diff(&row.feature, "$", &expected, &actual);
    let es = output_status(&expected);
    let a = output_status(&actual);
    // Keep status equality assertions, but matching failures never count as agreement.
    if let Some(diff) = difference {
        row.counts.compared += 1;
        row.counts.disagrees += 1;
        row.differences.push(json!({"id":id,"fixture":fixture,"record_index":index,"status":"disagrees","difference":diff,"reference_status":es,"kekule_status":a}));
    } else if es == "unsupported" || es == "no_analyzable_residues" {
        row.counts.unsupported += 1;
        row.differences.push(json!({"id":id,"fixture":fixture,"record_index":index,"status":"unsupported","reason":es}));
    } else if es != "ok" || a != "ok" {
        row.counts.errors += 1;
        row.differences.push(json!({"id":id,"fixture":fixture,"record_index":index,"status":"error","reference_status":es,"kekule_status":a}));
    } else {
        row.counts.compared += 1;
        row.counts.agrees += 1;
    }
    numeric_errors("$", &expected, &actual, &mut row.numeric_errors);
}

// Keep memory bounded while allowing records from separate files to share a batch.
const BATCH_SIZE: usize = 256;

struct Work {
    fixture: usize,
    id: String,
    record_index: usize,
    input: Input,
    expected: Value,
}

fn evaluate_batch(
    row: &mut ResultRow,
    pending: &mut Vec<Work>,
    pool: &rayon::ThreadPool,
    evaluate: &(impl Fn(&str, &str, &Input) -> Result<Value, Box<dyn Error>> + Sync),
    opts: &Options,
    reference: &Reference,
    progress: &ProgressBar,
) {
    if pending.is_empty() {
        return;
    }
    let start = Instant::now();
    let results: Vec<_> = pool.install(|| {
        pending
            .par_iter()
            .map(|work| {
                evaluate(&row.feature, &row.dataset, black_box(&work.input))
                    .map_err(|error| error.to_string())
            })
            .collect()
    });
    *row.time_ms.get_or_insert(0.0) += start.elapsed().as_secs_f64() * 1000.0;

    let mut inputs = Vec::new();
    let mut outputs = Vec::new();
    let mut ids = Vec::new();
    // Indexed parallel collection preserves source order, including failures.
    for (work, result) in pending.drain(..).zip(results) {
        let fixture = &mut row.fixtures[work.fixture];
        let path = fixture["path"]
            .as_str()
            .expect("fixture path is recorded")
            .to_owned();
        fixture["evaluated_ids"]
            .as_array_mut()
            .unwrap()
            .push(json!(work.id));
        fixture["evaluated_records"] = json!(fixture["evaluated_ids"].as_array().unwrap().len());
        match result {
            Ok(actual) => {
                let mut sizes = BTreeMap::new();
                workload_sizes(&actual, &mut sizes);
                for (key, count) in sizes {
                    let total = fixture["workload"][&key].as_u64().unwrap_or(0) + count;
                    fixture["workload"][key] = json!(total);
                }
                if opts.reference_python.is_some() {
                    inputs.push(work.input);
                    outputs.push(actual.clone());
                    ids.push(work.id.clone());
                }
                compare_record(
                    row,
                    &work.id,
                    &path,
                    work.record_index,
                    work.expected,
                    actual,
                );
            }
            Err(error) => {
                row.counts.errors += 1;
                row.differences.push(json!({"id":work.id,"fixture":path,"record_index":work.record_index,"status":"error","message":error}));
            }
        }
    }
    if let Some(python) = &opts.reference_python {
        match reference_run(python, &row.feature, &inputs, &outputs, reference) {
            Ok(result) if !result.is_null() => {
                row.live_reference_disagrees += result["disagrees"].as_u64().unwrap_or(0) as usize;
                row.live_references.push(json!({"ids":ids,"result":result}));
            }
            Ok(_) => (),
            Err(error) => row.execution_errors.push(json!({"status":"error","stage":"live reference","ids":ids,"message":error.to_string()})),
        }
    }
    update_progress(&row.counts, progress);
}

fn update_progress(counts: &Counts, progress: &ProgressBar) {
    progress.set_position(
        (counts.agrees + counts.disagrees + counts.unsupported + counts.errors) as u64,
    );
}

fn run_feature(
    dataset_id: &str,
    dataset: &Dataset,
    feature: &str,
    opts: &Options,
    progress: &ProgressBar,
    evaluate: impl Fn(&str, &str, &Input) -> Result<Value, Box<dyn Error>> + Sync,
) -> Result<ResultRow, Box<dyn Error>> {
    let reference = dataset
        .references
        .get(feature)
        .ok_or("no reference evidence for this feature/dataset")?;
    let selected = dataset.selection(opts.limit.unwrap_or(usize::MAX))?;
    progress.set_length(selected.len() as u64);
    let mut row = ResultRow {
        dataset: dataset_id.into(),
        feature: feature.into(),
        reference_tool: reference.tool.clone(),
        reference_version: reference.version.clone(),
        reference_kind: if reference.tool.ends_with("-manual-semantic") {
            "implementation snapshot; not independent validation"
        } else {
            "external reference"
        },
        preparation: reference.notes.clone(),
        source_lock_sha256: dataset.lock_sha256.clone(),
        selection_id: dataset.lock["selection_id"].clone(),
        selection_sha256: sha256(serde_json::to_string(&selected)?.as_bytes()),
        counts: Counts::default(),
        numeric_errors: BTreeMap::new(),
        time_ms: None,
        fixtures: vec![],
        differences: vec![],
        execution_errors: vec![],
        live_reference_disagrees: 0,
        live_references: vec![],
    };
    let mut seen = BTreeSet::new();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opts.jobs)
        .build()?;
    let mut pending = Vec::with_capacity(BATCH_SIZE);
    for fixture in &reference.fixtures {
        let members = dataset.members(fixture)?;
        if !members.iter().any(|id| selected.contains(id)) {
            continue;
        }
        let path = safe_join(&dataset.root, fixture)?;
        let result = (|| -> Result<(), Box<dyn Error>> {
            let source = Input::read(&path)?;
            let expected = dataset.reference(feature, fixture, &source, reference)?;
            let texts = split_records(&source, members.len())?;
            let fixture_index = row.fixtures.len();
            row.fixtures.push(json!({"path":fixture,"sha256":sha256(source.text.as_bytes()),"selected_ids":members.iter().filter(|id| selected.contains(*id)).collect::<Vec<_>>(),"workload":{},"evaluated_ids":[],"evaluated_records":0}));
            for (index, (id, text)) in members.iter().zip(texts).enumerate() {
                if !selected.contains(id) {
                    continue;
                }
                let input = Input {
                    path: path.clone(),
                    text,
                };
                let expected = expected_record(&expected, index)?;
                seen.insert(id.clone());
                row.counts.selected += 1;
                pending.push(Work {
                    fixture: fixture_index,
                    id: id.clone(),
                    record_index: index,
                    input,
                    expected,
                });
                if pending.len() == BATCH_SIZE {
                    evaluate_batch(
                        &mut row,
                        &mut pending,
                        &pool,
                        &evaluate,
                        opts,
                        reference,
                        progress,
                    );
                }
            }
            Ok(())
        })();
        if let Err(error) = result {
            let unhandled = members
                .iter()
                .filter(|id| selected.contains(*id) && !seen.contains(*id))
                .cloned()
                .collect::<Vec<_>>();
            for id in &unhandled {
                seen.insert(id.clone());
            }
            row.counts.selected += unhandled.len();
            row.counts.errors += unhandled.len();
            // Infrastructure errors are retained even if evaluation already finished.
            row.execution_errors.push(json!({"fixture":fixture,"status":"error","ids":unhandled,"message":error.to_string()}));
            update_progress(&row.counts, progress);
        }
    }
    evaluate_batch(
        &mut row,
        &mut pending,
        &pool,
        &evaluate,
        opts,
        reference,
        progress,
    );
    for id in selected.difference(&seen) {
        row.counts.selected += 1;
        row.counts.unsupported += 1;
        row.differences.push(json!({"id":id,"status":"unsupported","reason":"no input/reference for this feature in the pinned evidence"}));
    }
    update_progress(&row.counts, progress);
    progress.finish();
    Ok(row)
}

fn reference_run(
    python: &Path,
    feature: &str,
    inputs: &[Input],
    outputs: &[Value],
    reference: &Reference,
) -> Result<Value, Box<dyn Error>> {
    if reference.tool.ends_with("-manual-semantic") || inputs.is_empty() {
        return Ok(Value::Null);
    }
    let mut child = process::Command::new(python)
        .arg(Path::new(env!("CARGO_MANIFEST_DIR")).join("reference/run.py"))
        .stdin(process::Stdio::piped())
        .stdout(process::Stdio::piped())
        .stderr(process::Stdio::piped())
        .spawn()?;
    let request = json!({"feature":feature,"inputs":inputs.iter().map(|input| json!({"path":input.path,"text":input.text})).collect::<Vec<_>>()});
    child
        .stdin
        .take()
        .ok_or("missing reference stdin")?
        .write_all(serde_json::to_string(&request)?.as_bytes())?;
    let output = child.wait_with_output()?;
    if !output.status.success() {
        return Err(boxed_error(format!(
            "reference process failed: {}",
            String::from_utf8_lossy(&output.stderr)
        )));
    }
    let mut response: Value = serde_json::from_slice(&output.stdout)?;
    let expected = response
        .as_object_mut()
        .ok_or("reference response must be an object")?
        .remove("expected")
        .ok_or("reference response has no expected outputs")?;
    let expected = expected
        .as_array()
        .ok_or("reference outputs must be an array")?;
    if expected.len() != inputs.len() || outputs.len() != inputs.len() {
        return Err(boxed_error("live reference input/output count mismatch"));
    }
    let mut differences = Vec::new();
    for (index, (expected, actual)) in expected.iter().zip(outputs).enumerate() {
        let mut expected = expected.clone();
        let mut actual = actual.clone();
        normalize_benchmark_for_comparison_in_place(feature, &mut expected);
        normalize_benchmark_for_comparison_in_place(feature, &mut actual);
        if let Some(diff) = first_json_diff(feature, "$", &expected, &actual) {
            differences.push(json!({"selected_index":index,"difference":diff}));
        }
    }
    response["disagrees"] = json!(differences.len());
    response["differences"] = json!(differences);
    Ok(response)
}

pub(crate) fn run() -> Result<(), Box<dyn Error>> {
    let args = env::args().skip(1).collect::<Vec<_>>();
    if args.is_empty() || args == ["--help"] {
        println!("cargo benchmark --feature FEATURE|all --dataset DATASET|all [--limit N] [--jobs N] [--output FILE] [--reference-python PYTHON]\ncargo benchmark --list\nSee benchmarks/GUIDE.md for data preparation and measurement scope.");
        return Ok(());
    }
    if args == ["--list"] {
        for id in DATASETS {
            let dataset = Dataset::open(id)?;
            println!(
                "{id}: {}",
                dataset
                    .references
                    .keys()
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(", ")
            );
        }
        return Ok(());
    }
    let opts = options(&args)?;
    if opts.feature != "all" && !FEATURES.contains(&opts.feature.as_str()) {
        return Err(boxed_error(format!("unknown feature: {}", opts.feature)));
    }
    if opts.dataset != "all" && !DATASETS.contains(&opts.dataset.as_str()) {
        return Err(boxed_error(format!("unknown dataset: {}", opts.dataset)));
    }
    let mut rows = Vec::new();
    for id in DATASETS
        .iter()
        .filter(|id| opts.dataset == "all" || **id == opts.dataset)
    {
        let dataset = Dataset::open(id)?;
        for feature in FEATURES.iter().filter(|f| {
            (opts.feature == "all" || **f == opts.feature) && dataset.references.contains_key(**f)
        }) {
            let progress = ProgressBar::new(0)
                .with_style(
                    ProgressStyle::with_template(
                        "{prefix} {bar:20} {pos}/{len} {percent:>3}% {msg}",
                    )?
                    .progress_chars("██░"),
                )
                .with_prefix(format!("{id} {feature}"))
                .with_finish(ProgressFinish::AbandonWithMessage("stopped".into()));
            let row = run_feature(id, &dataset, feature, &opts, &progress, evaluate)?;
            println!(
                "{id} {feature}: {}/{} agree; {} disagree; {} unsupported; {} errors; Time {} ms",
                row.counts.agrees,
                row.counts.compared,
                row.counts.disagrees,
                row.counts.unsupported,
                row.counts.errors,
                row.time_ms
                    .map(|time| format!("{time:.3}"))
                    .unwrap_or_else(|| "n/a".into())
            );
            if row.reference_kind.starts_with("implementation") {
                println!("  snapshot agreement only; this feature has no independent reference");
            }
            if !row.execution_errors.is_empty() || row.live_reference_disagrees > 0 {
                println!(
                    "  {} execution errors; {} live-reference differences",
                    row.execution_errors.len(),
                    row.live_reference_disagrees
                );
            }
            rows.push(row);
        }
    }
    if rows.is_empty() {
        return Err(boxed_error(
            "no reference evidence for the requested feature/dataset",
        ));
    }
    let git = |args: &[&str]| {
        process::Command::new("git")
            .args(args)
            .current_dir(env!("CARGO_MANIFEST_DIR"))
            .output()
            .ok()
            .filter(|o| o.status.success())
            .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_owned())
    };
    let failed = rows.iter().any(|row| {
        row.counts.errors > 0
            || row.counts.disagrees > 0
            || !row.execution_errors.is_empty()
            || row.live_reference_disagrees > 0
    });
    let report = json!({"schema_version":3,"revision":git(&["rev-parse","HEAD"]),"working_tree_status":git(&["status","--porcelain"]),"os":env::consts::OS,"architecture":env::consts::ARCH,"machine":env::var("COMPUTERNAME").or_else(|_|env::var("HOSTNAME")).ok(),"cpu":env::var("PROCESSOR_IDENTIFIER").ok(),"debug_assertions":cfg!(debug_assertions),"rustc":process::Command::new("rustc").arg("--version").output().ok().map(|o|String::from_utf8_lossy(&o.stdout).trim().to_owned()),"timing_scope":"one parallel evaluation per record: parse + feature + result materialization; time_ms sums elapsed batch times, including scheduling and failed evaluations; preloaded input; excludes pool creation, correctness comparison, result drop and JSON encoding; no warmup or repetitions","results":rows});
    if let Some(parent) = opts.output.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    let mut file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&opts.output)?;
    file.write_all(serde_json::to_string_pretty(&report)?.as_bytes())?;
    println!("Report: {}", opts.output.display());
    if failed {
        return Err(boxed_error(
            "benchmark reported differences or errors; see report",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn row() -> ResultRow {
        ResultRow {
            dataset: "test".into(),
            feature: "io.smiles.parse".into(),
            reference_tool: "test".into(),
            reference_version: "1".into(),
            reference_kind: "external reference",
            preparation: vec![],
            source_lock_sha256: String::new(),
            selection_id: Value::Null,
            selection_sha256: String::new(),
            counts: Counts::default(),
            numeric_errors: BTreeMap::new(),
            time_ms: None,
            fixtures: vec![],
            differences: vec![],
            execution_errors: vec![],
            live_reference_disagrees: 0,
            live_references: vec![],
        }
    }
    #[test]
    fn matching_failures_and_missing_records_are_not_agreement() {
        let mut row = row();
        let failed = json!({"records":[{"record_index":0,"status":"parse_error"}]});
        compare_record(&mut row, "bad", "input", 0, failed.clone(), failed);
        compare_record(
            &mut row,
            "no-stereo",
            "input",
            0,
            json!({"records":[]}),
            json!({"records":[]}),
        );
        assert_eq!(row.counts.agrees, 0);
        assert_eq!(row.counts.errors, 1);
        assert_eq!(row.counts.unsupported, 1);
        let nested = json!({"records":[{"record_index":0,"status":"ok","normalized_perceived":{"status":"perception_error"}}]});
        compare_record(
            &mut row,
            "failed-perception",
            "input",
            0,
            nested.clone(),
            nested,
        );
        assert_eq!(row.counts.agrees, 0);
        assert_eq!(row.counts.errors, 2);
    }

    #[test]
    fn explicit_empty_cip_descriptors_agree_but_new_assignments_disagree() {
        let mut row = row();
        row.feature = "stereo.cip".into();
        let empty = json!({"records":[{"record_index":0,"status":"ok","atom_count":2,"bond_count":1,"atom_descriptors":[],"bond_descriptors":[]}]});
        compare_record(&mut row, "plain", "input", 0, empty.clone(), empty.clone());
        let mut assigned = empty.clone();
        assigned["records"][0]["atom_descriptors"] = json!([{"atom_index":0,"descriptor":"R"}]);
        compare_record(&mut row, "false-positive", "input", 0, empty, assigned);
        assert_eq!(row.counts.agrees, 1);
        assert_eq!(row.counts.disagrees, 1);
        assert_eq!(row.counts.unsupported, 0);
    }

    #[test]
    fn progress_counts_each_completed_outcome_and_leaves_pending_inputs_unfinished() {
        let mut row = row();
        row.counts.selected = 5;
        for (id, status, expected_value, actual_value) in [
            ("match", "ok", 1, 1),
            ("difference", "ok", 1, 2),
            ("unsupported", "unsupported", 1, 1),
            ("failed", "parse_error", 1, 1),
        ] {
            compare_record(
                &mut row,
                id,
                "input",
                0,
                json!({"records":[{"record_index":0,"status":status,"value":expected_value}]}),
                json!({"records":[{"record_index":0,"status":status,"value":actual_value}]}),
            );
        }
        let progress = ProgressBar::hidden();
        progress.set_length(5);
        update_progress(&row.counts, &progress);
        assert_eq!(progress.position(), 4);
        assert_eq!(progress.length(), Some(5));
        assert!(!progress.is_finished());
    }
    #[test]
    fn every_failing_record_is_retained_and_extra_fields_are_asserted() {
        let mut row = row();
        for index in [2, 5] {
            compare_record(
                &mut row,
                &index.to_string(),
                "input",
                index,
                json!({"records":[{"record_index":index,"status":"ok","charge":0}]}),
                json!({"records":[{"record_index":0,"status":"ok","charge":0,"extra":true}]}),
            );
        }
        assert_eq!(row.counts.disagrees, 2);
        assert_eq!(row.differences.len(), 2);
        assert_eq!(row.differences[1]["record_index"], 5);
    }
    #[test]
    fn missing_cip_reference_cannot_hide_a_new_assignment() {
        let mut row = row();
        compare_record(
            &mut row,
            "stereo",
            "input",
            0,
            json!({"records":[]}),
            json!({"records":[{"record_index":0,"status":"ok","cip":"R"}]}),
        );
        assert_eq!(row.counts.disagrees, 1);
        assert_eq!(row.counts.unsupported, 0);
    }
    #[test]
    fn cli_rejects_old_flags_duplicates_repetitions_and_zero_limit() {
        for args in [
            vec!["--benchmark", "all"],
            vec!["--feature", "x", "--feature", "y"],
            vec!["--feature", "x", "--dataset", "y", "--samples", "1"],
            vec!["--feature", "x", "--dataset", "y", "--limit", "0"],
            vec!["--feature", "x", "--dataset", "y", "--jobs", "0"],
        ] {
            assert!(options(&args.into_iter().map(str::to_owned).collect::<Vec<_>>()).is_err());
        }
    }
    #[test]
    fn numeric_errors_measure_bias_not_correlation() {
        let mut metrics = BTreeMap::new();
        numeric_errors(
            "$",
            &json!([1.0, 2.0, 3.0]),
            &json!([2.0, 4.0, 6.0]),
            &mut metrics,
        );
        assert_eq!(metrics["$[]"].mean_absolute_error, 2.0);
        assert_eq!(metrics["$[]"].maximum_absolute_error, 3.0);
    }
    #[test]
    fn selection_preserves_original_reference_index() {
        let expected = json!({"policy":"fixed","records":[{"record_index":2,"value":1},{"record_index":9,"value":2}]});
        assert_eq!(
            expected_record(&expected, 9).unwrap(),
            json!({"policy":"fixed","records":[{"record_index":9,"value":2}]})
        );
        assert_eq!(
            expected_record(&expected, 7).unwrap(),
            json!({"policy":"fixed","records":[]})
        );
    }
    #[test]
    fn each_record_is_evaluated_once_for_timing_and_correctness_even_after_error() {
        let dataset = Dataset::open("smoke").unwrap();
        let mut opts =
            options(&["--feature", "io.sdf.v2000.parse", "--dataset", "smoke"].map(str::to_owned))
                .unwrap();
        for fail_first in [false, true] {
            let mut baseline = None;
            for jobs in [1, 2] {
                opts.jobs = jobs;
                let calls = std::sync::Mutex::new(BTreeMap::new());
                let first_fixture = &dataset.references["io.sdf.v2000.parse"].fixtures[0];
                let progress = ProgressBar::hidden();
                let row = run_feature(
                    "smoke",
                    &dataset,
                    "io.sdf.v2000.parse",
                    &opts,
                    &progress,
                    |feature, dataset, input| {
                        *calls
                            .lock()
                            .unwrap()
                            .entry((input.path.clone(), input.text.clone()))
                            .or_insert(0) += 1;
                        if fail_first && input.path.ends_with(first_fixture) {
                            Err(boxed_error("injected evaluation failure"))
                        } else {
                            evaluate(feature, dataset, input)
                        }
                    },
                )
                .unwrap();
                let calls = calls.into_inner().unwrap();
                assert_eq!(calls.len(), 4);
                assert!(calls.values().all(|count| *count == 1));
                assert_eq!(row.counts.agrees, if fail_first { 3 } else { 4 });
                assert_eq!(row.counts.errors, usize::from(fail_first));
                assert_eq!(row.counts.disagrees, 0);
                assert!(row.execution_errors.is_empty());
                assert_eq!(progress.position(), row.counts.selected as u64);
                assert_eq!(progress.length(), Some(20));
                assert!(progress.is_finished());
                let evaluated_records: u64 = row
                    .fixtures
                    .iter()
                    .map(|fixture| fixture["evaluated_records"].as_u64().unwrap())
                    .sum();
                assert_eq!(evaluated_records, 4);
                assert!(row.time_ms.unwrap().is_finite());
                assert!(row.time_ms.unwrap() >= 0.0);
                let result = json!({"counts":row.counts,"differences":row.differences,"fixtures":row.fixtures});
                if let Some(expected) = &baseline {
                    assert_eq!(&result, expected);
                } else {
                    baseline = Some(result);
                }
            }
        }
    }

    #[test]
    fn records_from_separate_files_execute_concurrently() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let dataset = Dataset::open("smoke").unwrap();
        let opts = options(
            &[
                "--feature",
                "io.sdf.v2000.parse",
                "--dataset",
                "smoke",
                "--jobs",
                "2",
            ]
            .map(str::to_owned),
        )
        .unwrap();
        let entered = AtomicUsize::new(0);
        let row = run_feature(
            "smoke",
            &dataset,
            "io.sdf.v2000.parse",
            &opts,
            &ProgressBar::hidden(),
            |feature, dataset, input| {
                entered.fetch_add(1, Ordering::SeqCst);
                // A serial regression fails within five seconds instead of deadlocking.
                let deadline = Instant::now() + std::time::Duration::from_secs(5);
                while entered.load(Ordering::SeqCst) < 2 {
                    if Instant::now() >= deadline {
                        return Err(boxed_error("record evaluations did not overlap"));
                    }
                    std::thread::yield_now();
                }
                evaluate(feature, dataset, input)
            },
        )
        .unwrap();
        assert_eq!(entered.load(Ordering::SeqCst), 4);
        assert_eq!(row.counts.agrees, 4);
        assert_eq!(row.counts.errors, 0);
        assert_eq!(row.fixtures.len(), 4);
    }

    #[test]
    fn jobs_default_to_all_available_cpus() {
        let opts = options(&["--feature", "all", "--dataset", "smoke"].map(str::to_owned)).unwrap();
        assert_eq!(
            opts.jobs,
            std::thread::available_parallelism().map_or(1, usize::from)
        );
    }

    #[test]
    fn whole_dataset_is_selected_unless_limit_is_explicit() {
        // Metadata-only datasets cross every former automatic cap. Missing
        // reference inputs must stay visible in the selected/unsupported counts.
        for (id, feature, size) in [
            ("pdb-1000", "bio.secondary-structure.dssp", 1_000),
            ("pubchem-100k", "stereo.cip", 1_001),
            ("pubchem-100k", "algo.substructure.vf2", 1_001),
            ("pubchem-100k", "io.smiles.canonical", 1_001),
            ("pubchem-100k", "algo.canonical-ranking", 1_001),
            ("pubchem-100k", "io.smiles.parse", 100_001),
        ] {
            let dataset = Dataset {
                root: PathBuf::new(),
                lock: json!({"corpus_id":id,"entries":(0..size).map(|n| json!({"id":n.to_string()})).collect::<Vec<_>>()}),
                references: BTreeMap::from([(
                    feature.into(),
                    Reference {
                        tool: "test".into(),
                        version: "1".into(),
                        fixtures: vec![],
                        notes: vec![],
                    },
                )]),
                lock_sha256: String::new(),
            };
            let mut opts =
                options(&["--feature", feature, "--dataset", id, "--jobs", "1"].map(str::to_owned))
                    .unwrap();
            assert_eq!(opts.limit, None);
            for (limit, expected_count) in [(None, size), (Some(17), 17)] {
                opts.limit = limit;
                let progress = ProgressBar::hidden();
                let row = run_feature(id, &dataset, feature, &opts, &progress, |_, _, _| {
                    panic!("metadata-only dataset has no evaluation inputs")
                })
                .unwrap();
                assert_eq!(row.counts.selected, expected_count, "{id} {feature}");
                assert_eq!(row.counts.unsupported, expected_count);
                assert_eq!(row.counts.errors, 0);
                assert!(row.execution_errors.is_empty());
                assert_eq!(progress.position(), expected_count as u64);
                assert_eq!(progress.length(), Some(expected_count as u64));
                assert!(progress.is_finished());
            }
        }
    }
}
