//! Kekule observer for the scientific benchmark.
//!
//! Reads one JSON request per stdin line, `{id, task, format, text, options}`,
//! and writes one response per line in input order:
//! `{id, status: ok|error|panic, value?, error?: {kind, message}}`.
//!
//! A value is a list of facts, each `[key..., value]`, keyed by source
//! positions: SMILES atom order, Molfile atom-block rows, or mmCIF
//! `_atom_site.id`. The reference observers emit the same facts, so the
//! comparison never depends on how either toolkit stores chemistry. Routing,
//! selection and comparison live in the Python orchestrator; this binary only
//! reports what Kekule's public API computes.

mod bio;
mod input;
mod small;

use std::io::{self, BufRead, Write};
use std::panic::{self, AssertUnwindSafe};
use std::process::ExitCode;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use serde_json::{json, Value};

/// A failure reported for one input, never for the whole run.
#[derive(Debug)]
pub struct Failure {
    pub kind: &'static str,
    pub message: String,
}

impl Failure {
    pub fn new(kind: &'static str, error: impl std::fmt::Display) -> Self {
        Self {
            kind,
            message: error.to_string(),
        }
    }
}

pub type Facts = Vec<Value>;

fn dispatch(task: &str, format: &str, text: &str, options: &Value) -> Result<Facts, Failure> {
    if let Some(task) = task.strip_prefix("small.") {
        return small::observe(task, format, text, options);
    }
    if let Some(task) = task.strip_prefix("bio.") {
        return bio::observe(task, format, text, options);
    }
    Err(Failure::new("request", format!("unknown task {task:?}")))
}

fn respond(line: &str) -> Value {
    let request: Value = match serde_json::from_str(line) {
        Ok(request) => request,
        Err(error) => {
            return json!({"id": null, "status": "error",
                "error": {"kind": "request", "message": error.to_string()}})
        }
    };
    let id = request["id"].clone();
    let field = |name: &str| request[name].as_str().unwrap_or_default().to_owned();
    let (task, format, text) = (field("task"), field("format"), field("text"));
    let options = request.get("options").cloned().unwrap_or(Value::Null);
    match panic::catch_unwind(AssertUnwindSafe(|| {
        dispatch(&task, &format, &text, &options)
    })) {
        Ok(Ok(facts)) => json!({"id": id, "status": "ok", "value": facts}),
        Ok(Err(failure)) => json!({"id": id, "status": "error",
            "error": {"kind": failure.kind, "message": failure.message}}),
        Err(payload) => {
            let message = payload
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| payload.downcast_ref::<&str>().map(|s| (*s).to_owned()))
                .unwrap_or_else(|| "panic".to_owned());
            json!({"id": id, "status": "panic", "error": {"kind": "panic", "message": message}})
        }
    }
}

fn jobs() -> Result<usize, String> {
    let mut arguments = std::env::args().skip(1);
    let mut jobs = std::thread::available_parallelism().map_or(1, usize::from);
    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--jobs" => {
                jobs = arguments
                    .next()
                    .and_then(|value| value.parse().ok())
                    .filter(|&value| value > 0)
                    .ok_or("--jobs needs a positive integer")?;
            }
            other => return Err(format!("unknown argument {other:?}")),
        }
    }
    Ok(jobs)
}

fn main() -> ExitCode {
    let jobs = match jobs() {
        Ok(jobs) => jobs,
        Err(error) => {
            eprintln!("kekule-observe: {error}");
            return ExitCode::FAILURE;
        }
    };
    // Panics are reported per input; keep the default hook from printing them.
    panic::set_hook(Box::new(|_| {}));
    let lines = match io::stdin().lock().lines().collect::<Result<Vec<_>, _>>() {
        Ok(lines) => lines,
        Err(error) => {
            eprintln!("kekule-observe: {error}");
            return ExitCode::FAILURE;
        }
    };
    let lines = lines
        .into_iter()
        .filter(|line| !line.trim().is_empty())
        .collect::<Vec<_>>();
    let responses = Mutex::new(vec![Value::Null; lines.len()]);
    let next = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        for _ in 0..jobs.min(lines.len().max(1)) {
            scope.spawn(|| loop {
                let index = next.fetch_add(1, Ordering::Relaxed);
                let Some(line) = lines.get(index) else { break };
                let response = respond(line);
                responses.lock().expect("no poisoned responses")[index] = response;
            });
        }
    });
    let mut out = io::BufWriter::new(io::stdout().lock());
    for response in responses.into_inner().expect("no poisoned responses") {
        if writeln!(out, "{response}").is_err() {
            return ExitCode::FAILURE;
        }
    }
    if out.flush().is_err() {
        return ExitCode::FAILURE;
    }
    ExitCode::SUCCESS
}
