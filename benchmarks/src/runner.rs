use crate::*;
#[cfg(test)]
mod storage_tests;
mod stored;
use indicatif::{ProgressBar, ProgressFinish, ProgressStyle};
use rayon::prelude::*;
use std::io::{BufRead, BufReader};
use std::time::{Instant, SystemTime, UNIX_EPOCH};
use stored::{GeneratedGoldens, StoredGoldens};

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
                    return Some(format!("adapter returned {status}"));
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

#[derive(Deserialize, Serialize)]
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

fn reference_run(
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
    child
        .stdin
        .take()
        .ok_or("missing reference stdin")?
        .write_all(serde_json::to_string(request)?.as_bytes())?;
    let output = child.wait_with_output()?;
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

struct Options {
    feature: String,
    dataset: String,
    limit: usize,
    jobs: usize,
    generate: bool,
    python: Option<PathBuf>,
    goldens: PathBuf,
    output: PathBuf,
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
                Path::new(env!("CARGO_MANIFEST_DIR"))
                    .join("../target/benchmarks")
                    .join(format!(
                        "run-{}-{}.json",
                        process::id(),
                        SystemTime::now()
                            .duration_since(UNIX_EPOCH)
                            .unwrap()
                            .as_nanos()
                    ))
            }),
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
    kekule_ms: f64,
    reference_ms: f64,
}

fn comparison(
    feature: &str,
    expected: &Outcome,
    actual: &Outcome,
) -> (&'static str, Option<String>) {
    match (expected, actual) {
        (Outcome::Ok { value: expected }, Outcome::Ok { value: actual }) => {
            let collection = match feature {
                "io.mmcif.parse" => "blocks",
                "bio.secondary-structure.dssp" => "residues",
                _ => "records",
            };
            if [expected, actual].iter().any(|value| {
                value
                    .get(collection)
                    .and_then(Value::as_array)
                    .is_none_or(Vec::is_empty)
            }) {
                return ("error", Some(format!("missing or empty {collection}")));
            }
            if let Some(error) = failure(expected).or_else(|| failure(actual)) {
                return ("error", Some(error));
            }
            let mut expected = expected.clone();
            let mut actual = actual.clone();
            normalize_benchmark_for_comparison_in_place(feature, &mut expected);
            normalize_benchmark_for_comparison_in_place(feature, &mut actual);
            match first_json_diff("$", &expected, &actual) {
                None => ("agrees", None),
                Some(diff) => ("disagrees", Some(diff)),
            }
        }
        _ => (
            "error",
            Some("reference or Kekule evaluation failed".into()),
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
type Evaluator = fn(&str, &str, &Input) -> Result<Value, Box<dyn Error>>;
struct BatchRun<'a> {
    opts: &'a Options,
    pool: &'a rayon::ThreadPool,
    progress: &'a ProgressBar,
    stored: Result<StoredGoldens, String>,
    generated: Option<GeneratedGoldens>,
    results: &'a mut (dyn Write + Send + Sync),
    reference: ReferenceRunner,
    evaluate: Evaluator,
}

impl BatchRun<'_> {
    fn batch(&mut self, cases: &[Case], row: &mut Summary) -> Result<(), Box<dyn Error>> {
        if cases.is_empty() {
            return Ok(());
        }
        if let Some(goldens) = &mut self.generated {
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
            row.reference_ms += response.time_ms;
            for (case, expected) in cases.iter().zip(response.results) {
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
        let expected: Vec<_> = cases
            .iter()
            .map(|case| match &self.stored {
                Ok(goldens) => goldens.expected(case),
                Err(error) => Outcome::Error {
                    message: error.clone(),
                },
            })
            .collect();
        self.progress.set_message("Kekule");
        let start = Instant::now();
        let mut actual: Vec<_> = self.pool.install(|| {
            cases
                .par_iter()
                .map(|case| {
                    let value = Outcome::from_result((self.evaluate)(
                        &row.feature,
                        &row.dataset,
                        &case.input,
                    ));
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
        row.cases += 1;
        row.errors += 1;
        self.progress.inc(1);
        let outcome = Outcome::Error {
            message: message.into(),
        };
        if let Some(goldens) = &mut self.generated {
            json_line(
                goldens,
                &json!({"dataset":row.dataset,"feature":row.feature,"id":id,
                "fixture":fixture,"record_index":null,"input_sha256":null,"reference":null,"expected":outcome}),
            )
        } else {
            json_line(
                self.results,
                &json!({"dataset":row.dataset,"feature":row.feature,"id":id,
                "fixture":fixture,"record_index":null,"status":"error","difference":message,
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
            let (stored, generated) = if opts.generate {
                fs::create_dir_all(golden_path.parent().ok_or("golden path has no parent")?)?;
                (
                    Err("generation does not read goldens".into()),
                    Some(GeneratedGoldens::create(&golden_path)?),
                )
            } else {
                (
                    Ok(StoredGoldens::load(&golden_path, id, feature, &selected)?),
                    None,
                )
            };
            let mut run = BatchRun {
                opts: &opts,
                pool: &pool,
                progress: &progress,
                stored,
                generated,
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
            if let Some(generated) = run.generated.take() {
                generated.finish()?;
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
        }
    }
    let passed = !rows.is_empty()
        && rows.iter().all(|row| {
            row.cases > 0
                && if opts.generate {
                    row.errors == 0
                } else {
                    row.agrees == row.cases
                }
        });
    serde_json::to_writer_pretty(
        &mut report,
        &json!({"mode":if opts.generate { "generate" } else { "compare" },"passed":passed,"goldens":opts.goldens,"cases":if opts.generate { None } else { Some(result_path) },"results":rows}),
    )?;
    println!("Report: {}", opts.output.display());
    if !passed {
        return Err(boxed_error(if opts.generate {
            "reference generation retained errors; see report"
        } else {
            "benchmark comparison failed; see report"
        }));
    }
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
        assert_eq!(comparison("test", &expected, &actual).0, "disagrees");
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
