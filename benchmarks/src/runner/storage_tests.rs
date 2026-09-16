use super::*;

fn opts(generate: bool) -> Options {
    Options {
        feature: "io.smiles.parse".into(),
        dataset: "test".into(),
        limit: 4,
        jobs: 2,
        generate,
        python: Some("must-not-run".into()),
        goldens: "unused".into(),
        output: "unused".into(),
    }
}
fn case(id: &str) -> Case {
    Case {
        id: id.into(),
        fixture: "input.smi".into(),
        index: 0,
        input: Input {
            path: "input.smi".into(),
            text: "CC".into(),
        },
    }
}
fn golden(id: &str, expected: Outcome) -> Value {
    json!({"dataset":"test","feature":"io.smiles.parse","id":id,"fixture":"input.smi","record_index":0,
        "input_sha256":sha256(b"CC"),"reference":{"tool":"rdkit","version":"test"},"expected":expected})
}
fn temporary_file() -> (PathBuf, fs::File) {
    let path = env::temp_dir().join(format!(
        "kekule-golden-test-{}-{}",
        process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&path)
        .unwrap();
    (path, file)
}

#[test]
fn ordinary_options_do_not_require_python_or_generation() {
    let args = ["--feature", "io.smiles.parse", "--dataset", "smoke"].map(str::to_owned);
    let normal = options(&args).unwrap();
    assert!(!normal.generate);
    assert!(normal.python.is_none());
    let mut generation = vec!["generate".into()];
    generation.extend(args);
    let generation = options(&generation).unwrap();
    assert!(generation.generate);
    assert_eq!(generation.python, Some(PathBuf::from("python")));
}

#[test]
fn stored_comparison_never_invokes_reference_and_counts_every_outcome() {
    let opts = opts(false);
    let cases = [
        case("match"),
        case("reference-error"),
        case("stale"),
        case("missing"),
    ];
    let value = evaluate("io.smiles.parse", "test", &cases[0].input).unwrap();
    let mut stale = golden(
        "stale",
        Outcome::Ok {
            value: value.clone(),
        },
    );
    stale["input_sha256"] = json!(sha256(b"different source"));
    let contents = [
        golden("match", Outcome::Ok { value }),
        golden(
            "reference-error",
            Outcome::Error {
                message: "reference rejected source".into(),
            },
        ),
        stale,
    ]
    .iter()
    .map(Value::to_string)
    .collect::<Vec<_>>()
    .join("\n");
    let selected = cases.iter().map(|case| case.id.clone()).collect();
    let stored =
        StoredGoldens::read(contents.as_bytes(), "test", "io.smiles.parse", &selected).unwrap();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap();
    let progress = ProgressBar::hidden();
    progress.set_length(4);
    let (path, mut results) = temporary_file();
    let mut run = BatchRun {
        opts: &opts,
        pool: &pool,
        progress: &progress,
        stored: Ok(stored),
        generated: None,
        results: &mut results,
        reference: |_, _, _| panic!("ordinary comparison must not run a reference"),
        evaluate,
    };
    let mut row = Summary {
        dataset: "test".into(),
        feature: "io.smiles.parse".into(),
        ..Default::default()
    };
    run.batch(&cases, &mut row).unwrap();
    assert_eq!(
        (row.cases, row.agrees, row.disagrees, row.errors),
        (4, 1, 0, 3)
    );
    assert_eq!(progress.position(), 4);
    assert_eq!(progress.length(), Some(4));
    assert!(!progress.is_finished());
    drop(run);
    drop(results);
    let output = fs::read_to_string(&path).unwrap();
    assert_eq!(output.lines().count(), 4);
    assert!(output.contains("checksum differs"));
    assert!(output.contains("missing stored golden"));
    fs::remove_file(path).unwrap();
}

#[test]
fn generation_calls_reference_only_and_stores_its_errors() {
    let opts = opts(true);
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let progress = ProgressBar::hidden();
    progress.set_length(2);
    let (golden_path, golden_file) = temporary_file();
    drop(golden_file);
    fs::remove_file(&golden_path).unwrap();
    let (path, mut results) = temporary_file();
    let mut run = BatchRun {
        opts: &opts,
        pool: &pool,
        progress: &progress,
        stored: Err("unused".into()),
        generated: Some(GeneratedGoldens::create(&golden_path).unwrap()),
        results: &mut results,
        reference: |_, request, count| {
            assert_eq!(count, 2);
            assert!(request.get("written").is_none());
            assert_eq!(request["inputs"].as_array().unwrap().len(), 2);
            Ok(ReferenceResponse {
                reference: Reference {
                    tool: "rdkit".into(),
                    version: "test".into(),
                },
                time_ms: 1.0,
                results: vec![
                    Outcome::Ok {
                        value: json!({"records":[{"reference":7}]}),
                    },
                    Outcome::Error {
                        message: "reference error".into(),
                    },
                ],
            })
        },
        evaluate: |_, _, _| panic!("golden generation must not evaluate Kekule"),
    };
    let mut row = Summary {
        dataset: "test".into(),
        feature: "io.smiles.parse".into(),
        ..Default::default()
    };
    run.batch(&[case("success"), case("failure")], &mut row)
        .unwrap();
    run.generated.take().unwrap().finish().unwrap();
    assert_eq!((row.cases, row.errors, row.kekule_ms), (2, 1, 0.0));
    assert_eq!(progress.position(), 2);
    drop(run);
    drop(results);
    let stored = StoredGoldens::load(
        &golden_path,
        "test",
        "io.smiles.parse",
        &BTreeSet::from(["success".into(), "failure".into()]),
    )
    .unwrap();
    assert!(matches!(
        stored.expected(&case("success")),
        Outcome::Ok { .. }
    ));
    assert!(matches!(
        stored.expected(&case("failure")),
        Outcome::Error { .. }
    ));
    fs::remove_file(golden_path).unwrap();
    fs::remove_file(path).unwrap();
}

#[test]
fn invalid_or_duplicate_stored_goldens_are_rejected() {
    let baseline = golden(
        "id",
        Outcome::Ok {
            value: json!({"records":[{}]}),
        },
    );
    let selected = BTreeSet::from(["id".into()]);
    for (field, value) in [
        ("dataset", json!("other")),
        ("feature", json!("other")),
        ("input_sha256", json!("")),
        ("reference", Value::Null),
        ("extra", json!(true)),
    ] {
        let mut invalid = baseline.clone();
        invalid[field] = value;
        assert!(StoredGoldens::read(
            invalid.to_string().as_bytes(),
            "test",
            "io.smiles.parse",
            &selected
        )
        .is_err());
    }
    let duplicate = format!("{baseline}\n{baseline}\n");
    assert!(
        StoredGoldens::read(duplicate.as_bytes(), "test", "io.smiles.parse", &selected).is_err()
    );
}

#[test]
fn progress_includes_multiple_formats_and_missing_inputs() {
    let dataset = Dataset {
        root: PathBuf::new(),
        lock: json!({"entries":[
        {"id":"a","files":[{"path":"a.smi"},{"path":"a.sdf"}]},{"id":"b","files":[]}]}),
    };
    assert_eq!(
        progress_length(
            &dataset,
            &BTreeSet::from(["a".into(), "b".into()]),
            &dataset.fixtures("algo.rings.fast")
        )
        .unwrap(),
        3
    );
}

#[test]
fn interrupted_generation_is_not_published_and_existing_goldens_are_not_replaced() {
    let (path, file) = temporary_file();
    drop(file);
    fs::remove_file(&path).unwrap();
    {
        let mut output = GeneratedGoldens::create(&path).unwrap();
        output.write_all(b"unfinished").unwrap();
        assert!(!path.exists());
    }
    assert!(!path.exists());
    let mut output = GeneratedGoldens::create(&path).unwrap();
    output.write_all(b"complete").unwrap();
    output.finish().unwrap();
    let bytes = fs::read(&path).unwrap();
    assert!(GeneratedGoldens::create(&path).is_err());
    assert_eq!(fs::read(&path).unwrap(), bytes);
    fs::remove_file(path).unwrap();
}

#[test]
fn writer_validation_sends_only_emitted_text_and_uses_stored_expectation() {
    let mut opts = opts(false);
    opts.feature = "io.smiles.write".into();
    let mut golden = golden(
        "writer",
        Outcome::Ok {
            value: json!({"records":[{"checked":1}]}),
        },
    );
    golden["feature"] = json!(opts.feature);
    let stored = StoredGoldens::read(
        golden.to_string().as_bytes(),
        "test",
        &opts.feature,
        &BTreeSet::from(["writer".into()]),
    )
    .unwrap();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(1)
        .build()
        .unwrap();
    let progress = ProgressBar::hidden();
    let (path, mut results) = temporary_file();
    let mut run = BatchRun {
        opts: &opts,
        pool: &pool,
        progress: &progress,
        stored: Ok(stored),
        generated: None,
        results: &mut results,
        evaluate,
        reference: |_, request, count| {
            assert!(request.get("inputs").is_none());
            assert_eq!(count, 1);
            assert_eq!(request["written"][0]["status"], "ok");
            assert!(request["written"][0]["value"]["written"][0]["text"]
                .as_str()
                .unwrap()
                .starts_with("CC"));
            Ok(ReferenceResponse {
                reference: Reference {
                    tool: "rdkit".into(),
                    version: "test".into(),
                },
                time_ms: 1.0,
                results: vec![Outcome::Ok {
                    value: json!({"records":[{"checked":1}]}),
                }],
            })
        },
    };
    let mut row = Summary {
        dataset: "test".into(),
        feature: opts.feature.clone(),
        ..Default::default()
    };
    run.batch(&[case("writer")], &mut row).unwrap();
    assert_eq!((row.cases, row.agrees, row.errors), (1, 1, 0));
    drop(run);
    drop(results);
    fs::remove_file(path).unwrap();
}
