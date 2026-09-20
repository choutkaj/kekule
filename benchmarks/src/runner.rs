use crate::{
    boxed_error,
    compare::{differences, normalize_benchmark_for_comparison_in_place},
    dataset::{safe_join, sha256, split_records, Dataset, Input, DATASETS},
    features::{self, evaluate},
};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    env,
    error::Error,
    fs,
    io::Write,
    path::{Path, PathBuf},
    process,
};
mod dashboard;
#[cfg(test)]
mod storage_tests;
mod stored;
use indicatif::{ProgressBar, ProgressFinish, ProgressStyle};
use rayon::prelude::*;
use std::io::{Read, Seek, SeekFrom};
use std::time::{Instant, SystemTime, UNIX_EPOCH};
use stored::{GeneratedGoldens, Metadata, StoredGoldens};

const FEATURES: &[&str] = &[
    "io.smiles.parse",
    "io.smiles.write",
    "io.smiles.canonical",
    "io.smiles.isomeric",
    "io.mol.parse",
    "io.mol.v2000.write",
    "io.mol.v3000.write",
    "io.sdf.parse",
    "io.sdf.v2000.write",
    "io.mmcif.parse",
    "algo.rings.fast",
    "algo.rings.sssr",
    "algo.valence.rdkit-like",
    "algo.aromaticity.rdkit-like",
    "algo.aromaticity.mdl",
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

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "status", rename_all = "snake_case", deny_unknown_fields)]
enum Outcome {
    Ok { value: Value },
    Error { message: String },
}

impl Outcome {
    fn from_result(result: Result<Value, Box<dyn Error>>) -> Self {
        match result {
            Ok(value) => match failure(&value) {
                Some(message) => Self::Error { message },
                None => Self::Ok { value },
            },
            Err(error) => Self::Error {
                message: error.to_string(),
            },
        }
    }
}

// Neither empty records nor matching errors establish correctness.
fn failure(value: &Value) -> Option<String> {
    match value {
        Value::Object(fields) => {
            if let Some(status) = fields.get("status").and_then(Value::as_str) {
                if status != "ok" {
                    return Some(
                        fields
                            .get("message")
                            .and_then(Value::as_str)
                            .unwrap_or(status)
                            .to_owned(),
                    );
                }
            }
            if fields
                .get("records")
                .is_some_and(|records| records.as_array().is_none_or(Vec::is_empty))
            {
                return Some("adapter returned no records".into());
            }
            fields.values().find_map(failure)
        }
        Value::Array(values) => values.iter().find_map(failure),
        _ => None,
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
struct Reference {
    tool: String,
    version: String,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ReferenceResponse {
    reference: Reference,
    results: Vec<Outcome>,
    time_ms: f64,
}

fn reference_batch(
    python: &Path,
    request: &Value,
    count: usize,
) -> Result<ReferenceResponse, Box<dyn Error>> {
    let mut child = process::Command::new(python)
        .arg(Path::new(env!("CARGO_MANIFEST_DIR")).join("reference/run.py"))
        .stdin(process::Stdio::piped())
        .stdout(process::Stdio::piped())
        .stderr(process::Stdio::piped())
        .spawn()?;
    let mut input = child.stdin.take().ok_or("missing reference stdin")?;
    let request = serde_json::to_vec(request)?;
    let input_thread = std::thread::spawn(move || input.write_all(&request));
    let mut stdout = child.stdout.take().ok_or("missing reference stdout")?;
    let mut stderr = child.stderr.take().ok_or("missing reference stderr")?;
    let out_thread = std::thread::spawn(move || {
        let mut bytes = Vec::new();
        stdout.read_to_end(&mut bytes).map(|_| bytes)
    });
    let err_thread = std::thread::spawn(move || {
        let mut bytes = Vec::new();
        stderr.read_to_end(&mut bytes).map(|_| bytes)
    });
    let start = Instant::now();
    let mut timed_out = false;
    let status = loop {
        if let Some(status) = child.try_wait()? {
            break status;
        }
        if start.elapsed() > std::time::Duration::from_secs(300) {
            timed_out = true;
            child.kill()?;
            break child.wait()?;
        }
        std::thread::sleep(std::time::Duration::from_millis(10));
    };
    if timed_out {
        return Err(boxed_error(
            "reference exceeded its 300-second resource limit",
        ));
    }
    let stdout = out_thread
        .join()
        .map_err(|_| "reference stdout reader panicked")??;
    let stderr = err_thread
        .join()
        .map_err(|_| "reference stderr reader panicked")??;
    let input_result = input_thread
        .join()
        .map_err(|_| "reference stdin writer panicked")?;
    let output = process::Output {
        status,
        stdout,
        stderr,
    };
    if output.status.success() {
        input_result?;
    }
    if !output.status.success() {
        return Err(boxed_error(format!(
            "reference process failed: {}",
            String::from_utf8_lossy(&output.stderr)
        )));
    }
    let response: ReferenceResponse = serde_json::from_slice(&output.stdout)?;
    if response.results.len() != count || !response.time_ms.is_finite() || response.time_ms < 0.0 {
        return Err(boxed_error("invalid reference result count or time"));
    }
    Ok(response)
}

/// Retry a failed process batch one case at a time so a crashing reference
/// cannot erase the observations of unrelated cases.
fn reference_run(
    python: &Path,
    request: &Value,
    count: usize,
) -> Result<ReferenceResponse, Box<dyn Error>> {
    isolate_reference(request, count, |request, count| {
        reference_batch(python, request, count)
    })
}
fn isolate_reference(
    request: &Value,
    count: usize,
    mut run: impl FnMut(&Value, usize) -> Result<ReferenceResponse, Box<dyn Error>>,
) -> Result<ReferenceResponse, Box<dyn Error>> {
    match run(request, count) {
        Ok(response) => Ok(response),
        Err(batch_error) => {
            let description = run(&json!({"feature":request["feature"],"describe":true}), 0)?;
            let field = if request.get("written").is_some() {
                "written"
            } else {
                "inputs"
            };
            let mut response = ReferenceResponse {
                reference: description.reference,
                results: Vec::new(),
                time_ms: 0.0,
            };
            for item in request[field]
                .as_array()
                .ok_or("invalid reference request")?
            {
                if count == 1 {
                    response.results.push(Outcome::Error {
                        message: batch_error.to_string(),
                    });
                } else {
                    match run(&json!({"feature":request["feature"],field:[item]}), 1) {
                        Ok(single) => {
                            if single.reference != response.reference {
                                return Err(boxed_error(
                                    "reference version changed during isolation",
                                ));
                            }
                            response.time_ms += single.time_ms;
                            response.results.extend(single.results);
                        }
                        Err(error) => response.results.push(Outcome::Error {
                            message: error.to_string(),
                        }),
                    }
                }
            }
            Ok(response)
        }
    }
}

struct Options {
    feature: String,
    dataset: String,
    limit: usize,
    jobs: usize,
    generate: bool,
    python: Option<PathBuf>,
    goldens: PathBuf,
    output: PathBuf,
    runs_dir: PathBuf,
    started_at_unix_ms: u64,
}

fn options(args: &[String]) -> Result<Options, Box<dyn Error>> {
    let generate = args.first().is_some_and(|arg| arg == "generate");
    let args = if generate { &args[1..] } else { args };
    let mut values = BTreeMap::new();
    if !args.len().is_multiple_of(2) {
        return Err(boxed_error("each option requires a value"));
    }
    for pair in args.chunks_exact(2) {
        if ![
            "--feature",
            "--dataset",
            "--limit",
            "--jobs",
            "--writer-python",
            "--python",
            "--goldens",
            "--output",
        ]
        .contains(&pair[0].as_str())
        {
            return Err(boxed_error(format!("unknown option: {}", pair[0])));
        }
        if values.insert(pair[0].as_str(), pair[1].as_str()).is_some() {
            return Err(boxed_error("duplicate option"));
        }
    }
    let positive = |flag, default| -> Result<usize, Box<dyn Error>> {
        match values.get(flag) {
            None => Ok(default),
            Some(value) => value
                .parse::<usize>()
                .ok()
                .filter(|n| *n > 0)
                .ok_or_else(|| boxed_error(format!("{flag} must be a positive integer"))),
        }
    };
    if !generate && values.contains_key("--python") {
        return Err(boxed_error("--python is only used by the generate command"));
    }
    if generate && values.contains_key("--writer-python") {
        return Err(boxed_error("use --python for golden generation"));
    }
    let started_at_unix_ms = SystemTime::now().duration_since(UNIX_EPOCH)?.as_millis() as u64;
    let runs_dir = env::var_os("KEKULE_BENCHMARK_RUNS_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| Path::new(env!("CARGO_MANIFEST_DIR")).join("runs"));
    Ok(Options {
        feature: values
            .get("--feature")
            .ok_or("--feature is required")?
            .to_string(),
        dataset: values
            .get("--dataset")
            .ok_or("--dataset is required")?
            .to_string(),
        limit: positive("--limit", usize::MAX)?,
        jobs: positive(
            "--jobs",
            std::thread::available_parallelism().map_or(1, usize::from),
        )?,
        generate,
        python: if generate {
            Some(
                values
                    .get("--python")
                    .map(PathBuf::from)
                    .unwrap_or_else(|| "python".into()),
            )
        } else {
            values.get("--writer-python").map(PathBuf::from)
        },
        goldens: values
            .get("--goldens")
            .map(PathBuf::from)
            .unwrap_or_else(|| Path::new(env!("CARGO_MANIFEST_DIR")).join("goldens")),
        output: values
            .get("--output")
            .map(PathBuf::from)
            .unwrap_or_else(|| {
                runs_dir.join(format!("run-{started_at_unix_ms}-{}.json", process::id()))
            }),
        runs_dir,
        started_at_unix_ms,
    })
}

#[derive(Default, Serialize)]
struct Summary {
    dataset: String,
    feature: String,
    source_ids: usize,
    cases: usize,
    agrees: usize,
    disagrees: usize,
    errors: usize,
    not_applicable: usize,
    input_errors: usize,
    reference_errors: usize,
    kekule_errors: usize,
    writer_validation_errors: usize,
    observation_errors: usize,
    exact_agrees: usize,
    structural_differences: usize,
    numerical_differences: usize,
    golden: Option<Metadata>,
    kekule_ms: f64,
    reference_ms: f64,
}

fn comparison(feature: &str, expected: &Outcome, actual: &Outcome) -> (&'static str, Value) {
    match (expected, actual) {
        (Outcome::Ok { value: expected }, Outcome::Ok { value: actual }) => {
            if let Err(error) = crate::observation::validate(feature, expected)
                .and_then(|()| crate::observation::validate(feature, actual))
            {
                return ("error", json!({"cause":error.to_string()}));
            }
            let mut expected = expected.clone();
            let mut actual = actual.clone();
            // The public model writer has no title or SDF-property input.
            // Its source reference remains unmodified in the case report.
            if feature.starts_with("io.mol.") && features::is_writer(feature) {
                for record in expected["records"].as_array_mut().unwrap() {
                    record["title"] = json!("");
                    record["properties"] = json!([]);
                }
            }
            normalize_benchmark_for_comparison_in_place(feature, &mut expected);
            normalize_benchmark_for_comparison_in_place(feature, &mut actual);
            let diff = differences(feature, &expected, &actual);
            (
                if diff.agrees() { "agrees" } else { "disagrees" },
                json!({"exact":diff.exact(),"differences":diff}),
            )
        }
        _ => (
            "error",
            json!({"reference_failed":matches!(expected,Outcome::Error{..}),"kekule_failed":matches!(actual,Outcome::Error{..})}),
        ),
    }
}

struct Case {
    id: String,
    fixture: String,
    index: usize,
    input: Input,
}

fn json_line<W: Write + ?Sized>(file: &mut W, value: &Value) -> Result<(), Box<dyn Error>> {
    let mut line = serde_json::to_vec(value)?;
    line.push(b'\n');
    file.write_all(&line)?;
    Ok(())
}

type ReferenceRunner = fn(&Path, &Value, usize) -> Result<ReferenceResponse, Box<dyn Error>>;
type Evaluator = fn(&str, &Input) -> Result<Value, Box<dyn Error>>;
enum RunMode {
    Compare(StoredGoldens),
    Generate {
        writer: GeneratedGoldens,
        metadata: Metadata,
    },
}
struct BatchRun<'a> {
    opts: &'a Options,
    pool: &'a rayon::ThreadPool,
    progress: &'a ProgressBar,
    mode: RunMode,
    results: &'a mut (dyn Write + Send + Sync),
    reference: ReferenceRunner,
    evaluate: Evaluator,
}

impl BatchRun<'_> {
    fn batch(&mut self, cases: &[Case], row: &mut Summary) -> Result<(), Box<dyn Error>> {
        if cases.is_empty() {
            return Ok(());
        }
        if let RunMode::Generate {
            writer: goldens,
            metadata,
        } = &mut self.mode
        {
            // Explicit generation invokes only the reference. Never evaluate Kekule.
            self.progress.set_message("reference");
            let request = json!({"feature": row.feature, "inputs": cases.iter()
                .map(|case| json!({"path":case.input.path,"text":case.input.text})).collect::<Vec<_>>()});
            let response = (self.reference)(
                self.opts
                    .python
                    .as_deref()
                    .ok_or("generation requires Python")?,
                &request,
                cases.len(),
            )?;
            if metadata
                .reference
                .as_ref()
                .is_some_and(|reference| reference != &response.reference)
            {
                return Err(boxed_error("mixed reference versions during generation"));
            }
            metadata.reference = Some(response.reference.clone());
            row.reference_ms += response.time_ms;
            for (case, expected) in cases.iter().zip(response.results) {
                if let Outcome::Ok { value } = &expected {
                    crate::observation::validate(&row.feature, value)?;
                }
                json_line(
                    goldens,
                    &json!({"dataset":row.dataset,"feature":row.feature,"id":case.id,
                    "fixture":case.fixture,"record_index":case.index,"input_sha256":sha256(case.input.text.as_bytes()),
                    "reference":response.reference,"expected":expected}),
                )?;
                row.cases += 1;
                if matches!(expected, Outcome::Error { .. }) {
                    row.errors += 1;
                }
                self.progress.inc(1);
            }
            return Ok(());
        }
        let RunMode::Compare(stored) = &mut self.mode else {
            unreachable!()
        };
        let expected = cases
            .iter()
            .map(|case| stored.expected(case))
            .collect::<Result<Vec<_>, _>>()?;
        self.progress.set_message("Kekule");
        let start = Instant::now();
        let mut actual: Vec<_> = self.pool.install(|| {
            cases
                .par_iter()
                .map(|case| {
                    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        (self.evaluate)(&row.feature, &case.input)
                    }));
                    let value = match result {
                        Ok(result) => Outcome::from_result(result),
                        Err(payload) => Outcome::Error {
                            message: format!(
                                "Kekule panicked: {}",
                                payload
                                    .downcast_ref::<String>()
                                    .map(String::as_str)
                                    .or_else(|| payload.downcast_ref::<&str>().copied())
                                    .unwrap_or("unknown panic")
                            ),
                        },
                    };
                    self.progress.inc(1);
                    value
                })
                .collect()
        });
        row.kekule_ms += start.elapsed().as_secs_f64() * 1000.0;
        let emitted = features::is_writer(&row.feature).then(|| actual.clone());
        if features::is_writer(&row.feature) {
            // This only reads newly written text. It never recalculates goldens.
            self.progress.set_message("validate writer output");
            let python = self
                .opts
                .python
                .clone()
                .or_else(|| env::var_os("KEKULE_WRITER_PYTHON").map(PathBuf::from))
                .unwrap_or_else(|| "python".into());
            let request = json!({"feature":row.feature,"written":actual});
            match (self.reference)(&python, &request, cases.len()) {
                Ok(response) => {
                    row.reference_ms += response.time_ms;
                    if stored.metadata.reference.as_ref()!=Some(&response.reference) {
                        return Err(boxed_error("writer reader version differs from stored reference"));
                    }
                    actual = response.results;
                },
                Err(error) => actual = vec![Outcome::Error { message: format!("writer validation failed: {error}; select its interpreter with --writer-python or KEKULE_WRITER_PYTHON") }; cases.len()],
            }
        }
        for (case_index, ((case, expected), actual)) in
            cases.iter().zip(&expected).zip(&actual).enumerate()
        {
            let (status, difference) = comparison(&row.feature, expected, actual);
            row.cases += 1;
            row.reference_errors += usize::from(matches!(expected, Outcome::Error { .. }));
            let implementation_outcome = emitted
                .as_ref()
                .map_or(actual, |values| &values[case_index]);
            row.kekule_errors +=
                usize::from(matches!(implementation_outcome, Outcome::Error { .. }));
            row.writer_validation_errors += usize::from(
                emitted.is_some()
                    && matches!(implementation_outcome, Outcome::Ok { .. })
                    && matches!(actual, Outcome::Error { .. }),
            );
            row.observation_errors += usize::from(difference.get("cause").is_some());
            row.exact_agrees += usize::from(difference["exact"].as_bool() == Some(true));
            row.numerical_differences +=
                difference["differences"]["numerical"].as_u64().unwrap_or(0) as usize;
            row.structural_differences += difference["differences"]["structural"]
                .as_u64()
                .unwrap_or(0) as usize;
            match status {
                "agrees" => row.agrees += 1,
                "disagrees" => row.disagrees += 1,
                _ => row.errors += 1,
            }
            json_line(
                self.results,
                &json!({"dataset":row.dataset,"feature":row.feature,"id":case.id,
                "fixture":case.fixture,"record_index":case.index,"status":status,"difference":difference,
                "expected":expected,"actual":actual,"written":emitted.as_ref().map(|values| &values[case_index])}),
            )?;
        }
        self.progress.set_message(format!(
            "{} agree, {} differ, {} errors",
            row.agrees, row.disagrees, row.errors
        ));
        Ok(())
    }

    fn input_error(
        &mut self,
        row: &mut Summary,
        id: &str,
        fixture: Option<&str>,
        message: &str,
    ) -> Result<(), Box<dyn Error>> {
        if matches!(self.mode, RunMode::Generate { .. }) && fixture.is_some() {
            return Err(boxed_error(format!(
                "cannot generate a reference from invalid input {id}: {message}"
            )));
        }
        if fixture.is_some() {
            row.cases += 1;
            row.errors += 1;
            row.input_errors += 1;
        } else {
            row.not_applicable += 1;
        }
        self.progress.inc(1);
        let outcome = Outcome::Error {
            message: message.into(),
        };
        if let RunMode::Generate {
            writer: goldens, ..
        } = &mut self.mode
        {
            json_line(
                goldens,
                &json!({"dataset":row.dataset,"feature":row.feature,"id":id,
                "fixture":fixture,"record_index":null,"input_sha256":null,"reference":null,"expected":outcome}),
            )
        } else {
            json_line(
                self.results,
                &json!({"dataset":row.dataset,"feature":row.feature,"id":id,
                "fixture":fixture,"record_index":null,"status":if fixture.is_none(){"not_applicable"}else{"input_error"},"difference":message,
                "expected":outcome,"actual":outcome,"written":null}),
            )
        }
    }
}

fn progress_length(
    dataset: &Dataset,
    selected: &BTreeSet<String>,
    fixtures: &[String],
) -> Result<u64, Box<dyn Error>> {
    let mut seen = BTreeSet::new();
    let mut count = 0;
    for fixture in fixtures {
        for id in dataset.members(fixture)? {
            if selected.contains(&id) {
                count += 1;
                seen.insert(id);
            }
        }
    }
    Ok(count + selected.difference(&seen).count() as u64)
}

pub(crate) fn run() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = env::args().skip(1).collect();
    if args.is_empty() || args == ["--help"] {
        println!("cargo benchmark --feature FEATURE|all --dataset DATASET|all [--limit N] [--jobs N] [--goldens DIR] [--output FILE] [--writer-python PYTHON]\ncargo benchmark generate --feature FEATURE|all --dataset DATASET|all [--python PYTHON] [--goldens DIR] [--limit N]\ncargo benchmark --list");
        return Ok(());
    }
    if args == ["--list"] {
        println!(
            "Datasets: {}\nFeatures:\n{}",
            DATASETS.join(", "),
            FEATURES.join("\n")
        );
        return Ok(());
    }
    let opts = options(&args)?;
    if opts.feature != "all" && !FEATURES.contains(&opts.feature.as_str()) {
        return Err(boxed_error("unknown feature"));
    }
    if opts.dataset != "all" && !DATASETS.contains(&opts.dataset.as_str()) {
        return Err(boxed_error("unknown dataset"));
    }
    if let Some(parent) = opts.output.parent().filter(|p| !p.as_os_str().is_empty()) {
        fs::create_dir_all(parent)?;
    }
    let create = |path: &Path| {
        fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)
    };
    let mut report = create(&opts.output)?;
    let result_path = opts.output.with_extension("cases.jsonl");
    let mut results: Box<dyn Write + Send + Sync> = if opts.generate {
        Box::new(std::io::sink())
    } else {
        Box::new(create(&result_path)?)
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opts.jobs)
        .build()?;
    let mut rows = Vec::new();
    let implementation = implementation_identity()?;
    write_report(
        &mut report,
        &opts,
        &rows,
        &implementation,
        false,
        false,
        None,
    )?;
    let execution = (|| -> Result<(), Box<dyn Error>> {
        for id in DATASETS
            .iter()
            .filter(|id| opts.dataset == "all" || **id == opts.dataset)
        {
            let dataset = Dataset::open(id)?;
            let selected = dataset.selection(opts.limit)?;
            for feature in FEATURES
                .iter()
                .filter(|feature| opts.feature == "all" || **feature == opts.feature)
            {
                let mut row = Summary {
                    dataset: id.to_string(),
                    feature: feature.to_string(),
                    source_ids: selected.len(),
                    ..Summary::default()
                };
                let fixtures = dataset.fixtures(feature);
                let progress = ProgressBar::new(progress_length(&dataset, &selected, &fixtures)?)
                    .with_style(
                        ProgressStyle::with_template(
                            "{prefix} {bar:20} {pos}/{len} {percent:>3}% {msg}",
                        )?
                        .progress_chars("██░"),
                    )
                    .with_prefix(format!("{id} {feature}"))
                    .with_finish(ProgressFinish::AbandonWithMessage("stopped".into()));
                progress.enable_steady_tick(std::time::Duration::from_millis(100));
                progress.set_message(if opts.generate {
                    "reference"
                } else {
                    "loading stored goldens"
                });
                let golden_path = opts.goldens.join(id).join(format!("{feature}.jsonl.gz"));
                let lock_hash =
                    stored::text_hash(&fs::read_to_string(dataset.root.join("sources.lock.json"))?);
                let mode = if opts.generate {
                    fs::create_dir_all(golden_path.parent().ok_or("golden parent missing")?)?;
                    RunMode::Generate {
                        writer: GeneratedGoldens::create(&golden_path)?,
                        metadata: Metadata {
                            schema: 2,
                            contract_sha256: stored::feature_contract_hash(feature),
                            dataset: id.to_string(),
                            feature: feature.to_string(),
                            input_lock_sha256: lock_hash,
                            sha256: String::new(),
                            reference: None,
                            reference_code_sha256: Some(stored::reference_code_hash()?),
                            origin: "independent generation".into(),
                            cases: 0,
                        },
                    }
                } else {
                    let stored = StoredGoldens::load(&golden_path, id, feature, &lock_hash)?;
                    row.golden = Some(stored.metadata.clone());
                    RunMode::Compare(stored)
                };
                let mut run = BatchRun {
                    opts: &opts,
                    pool: &pool,
                    progress: &progress,
                    mode,
                    results: &mut *results,
                    reference: reference_run,
                    evaluate,
                };
                let mut seen = BTreeSet::new();
                let mut pending = Vec::new();
                for fixture in fixtures {
                    let members = dataset.members(&fixture)?;
                    if !members.iter().any(|id| selected.contains(id)) {
                        continue;
                    }
                    let read = (|| -> Result<_, Box<dyn Error>> {
                        let path = safe_join(&dataset.root, &fixture)?;
                        let input = Input::read(&path)?;
                        dataset.verify(&fixture, &input)?;
                        Ok((path, split_records(&input, members.len())?))
                    })();
                    match read {
                        Ok((path, texts)) => {
                            if members.len() == 1 && texts.len() > 1 {
                                progress.inc_length((texts.len() - 1) as u64);
                            }
                            for (index, text) in texts.into_iter().enumerate() {
                                let id = if members.len() == 1 {
                                    &members[0]
                                } else {
                                    &members[index]
                                };
                                if !selected.contains(id) {
                                    continue;
                                }
                                seen.insert(id.clone());
                                pending.push(Case {
                                    id: id.clone(),
                                    fixture: fixture.clone(),
                                    index,
                                    input: Input {
                                        path: path.clone(),
                                        text,
                                    },
                                });
                                if pending.len() == 256 {
                                    run.batch(&pending, &mut row)?;
                                    pending.clear();
                                }
                            }
                        }
                        Err(error) => {
                            run.batch(&pending, &mut row)?;
                            pending.clear();
                            for id in members.iter().filter(|id| selected.contains(*id)) {
                                seen.insert(id.clone());
                                run.input_error(&mut row, id, Some(&fixture), &error.to_string())?;
                            }
                        }
                    }
                }
                run.batch(&pending, &mut row)?;
                for id in selected.difference(&seen) {
                    run.input_error(&mut row, id, None, "no source input in the required format")?;
                }
                match run.mode {
                    RunMode::Generate {
                        writer,
                        mut metadata,
                    } => {
                        metadata.cases = row.cases + row.not_applicable;
                        writer.finish(metadata)?;
                    }
                    RunMode::Compare(mut stored) => {
                        progress.set_message("validate remaining stored records");
                        stored.finish()?;
                    }
                }
                progress.finish_and_clear();
                if opts.generate {
                    println!(
                        "{} {}: stored {} reference outcomes ({} errors) in {}",
                        row.dataset,
                        row.feature,
                        row.cases,
                        row.errors,
                        golden_path.display()
                    );
                } else {
                    println!(
                        "{} {}: {}/{} agree; {} disagree; {} errors",
                        row.dataset, row.feature, row.agrees, row.cases, row.disagrees, row.errors
                    );
                }
                rows.push(row);
                write_report(
                    &mut report,
                    &opts,
                    &rows,
                    &implementation,
                    false,
                    false,
                    None,
                )?;
            }
        }
        Ok(())
    })();
    if let Err(error) = execution {
        write_report(
            &mut report,
            &opts,
            &rows,
            &implementation,
            false,
            false,
            Some(error.to_string()),
        )?;
        dashboard::publish(&opts);
        return Err(error);
    }
    let passed = !rows.is_empty()
        && rows.iter().any(|row| row.cases > 0)
        && rows.iter().all(|row| {
            if opts.generate {
                row.errors == 0
            } else {
                row.agrees == row.cases
            }
        });
    write_report(
        &mut report,
        &opts,
        &rows,
        &implementation,
        true,
        passed,
        None,
    )?;
    println!("Report: {}", opts.output.display());
    dashboard::publish(&opts);
    if !passed {
        return Err(boxed_error(if opts.generate {
            "reference generation retained errors; see report"
        } else {
            "benchmark comparison failed; see report"
        }));
    }
    Ok(())
}

fn implementation_identity() -> Result<Value, Box<dyn Error>> {
    let git = |args: &[&str]| -> Option<String> {
        let output = process::Command::new("git")
            .args(args)
            .current_dir(env!("CARGO_MANIFEST_DIR"))
            .output()
            .ok()?;
        if !output.status.success() {
            return None;
        }
        Some(String::from_utf8(output.stdout).ok()?.trim().into())
    };
    let executable = stored::file_hash(&env::current_exe()?)?;
    let status = git(&["status", "--porcelain"]);
    Ok(
        json!({"revision":git(&["rev-parse","HEAD"]),"dirty":status.as_ref().map(|status|!status.is_empty()),"working_tree_status_sha256":status.as_ref().map(|status|sha256(status.as_bytes())),
        "executable_sha256":executable,"reference_code_sha256":stored::reference_code_hash()?,"contract_sha256":stored::contract_hash(),"feature_contracts":{"query.smarts":stored::feature_contract_hash("query.smarts")}}),
    )
}
fn write_report(
    file: &mut fs::File,
    opts: &Options,
    rows: &[Summary],
    implementation: &Value,
    complete: bool,
    passed: bool,
    error: Option<String>,
) -> Result<(), Box<dyn Error>> {
    file.seek(SeekFrom::Start(0))?;
    file.set_len(0)?;
    serde_json::to_writer_pretty(
        &mut *file,
        &json!({"schema":2,"mode":if opts.generate{"generate"}else{"compare"},
        "started_at_unix_ms":opts.started_at_unix_ms,
        "complete":complete,"passed":passed,"error":error,"implementation":implementation,"goldens":opts.goldens,
        "cases":if opts.generate{None}else{Some(opts.output.with_extension("cases.jsonl"))},"results":rows}),
    )?;
    file.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn errors_empty_records_and_extra_fields_never_pass() {
        for value in [
            Value::Null,
            json!({}),
            json!({"records":[]}),
            json!({"status":"unsupported"}),
            json!({"status":"no_analyzable_residues"}),
        ] {
            let outcome = Outcome::Ok { value };
            assert_eq!(comparison("test", &outcome, &outcome).0, "error");
        }
        let error = Outcome::Error {
            message: "failed".into(),
        };
        assert_eq!(comparison("test", &error, &error).0, "error");
        let expected = Outcome::Ok {
            value: json!({"records":[{"x":1}]}),
        };
        let actual = Outcome::Ok {
            value: json!({"records":[{"x":1,"extra":0}]}),
        };
        assert_eq!(comparison("test", &expected, &actual).0, "error");
    }
    #[test]
    fn protocol_rejects_ambiguous_outcomes_and_extra_metadata() {
        assert!(serde_json::from_value::<Outcome>(
            json!({"status":"ok","value":{},"message":"error"})
        )
        .is_err());
        assert!(serde_json::from_value::<Reference>(
            json!({"tool":"snapshot","version":"1","notes":[]})
        )
        .is_err());
    }
}
