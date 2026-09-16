use super::*;

struct Fixture(PathBuf);
impl Fixture {
    fn new() -> Self {
        let path = env::temp_dir().join(format!(
            "kekule-golden-test-{}-{}",
            process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn path(&self) -> PathBuf {
        self.0.join("reference.jsonl.gz")
    }
    fn store(&self, feature: &str, entries: &[Value]) -> StoredGoldens {
        let mut writer = GeneratedGoldens::create(&self.path()).unwrap();
        for entry in entries {
            json_line(&mut writer, entry).unwrap();
        }
        writer.finish(metadata(feature, entries.len())).unwrap();
        self.load(feature).unwrap()
    }
    fn load(&self, feature: &str) -> Result<StoredGoldens, Box<dyn Error>> {
        StoredGoldens::load(&self.path(), "test", feature, "lock")
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}
fn metadata(feature: &str, cases: usize) -> Metadata {
    Metadata {
        schema: 2,
        contract_sha256: stored::contract_hash(),
        dataset: "test".into(),
        feature: feature.into(),
        input_lock_sha256: "lock".into(),
        sha256: String::new(),
        reference: Some(reference()),
        reference_code_sha256: None,
        origin: "unit regression".into(),
        cases,
    }
}
fn reference() -> Reference {
    Reference {
        tool: "rdkit".into(),
        version: "test".into(),
    }
}
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
fn case(index: usize) -> Case {
    Case {
        id: index.to_string(),
        fixture: "input.smi".into(),
        index,
        input: Input {
            path: "input.smi".into(),
            text: "CC".into(),
        },
    }
}
fn observation() -> Value {
    evaluate("io.smiles.parse", &case(0).input).unwrap()
}
fn golden(feature: &str, index: usize, expected: Outcome) -> Value {
    json!({"dataset":"test","feature":feature,"id":index.to_string(),"fixture":"input.smi","record_index":index,"input_sha256":sha256(b"CC"),"reference":reference(),"expected":expected})
}
fn success() -> Outcome {
    Outcome::Ok {
        value: observation(),
    }
}
fn pool() -> rayon::ThreadPool {
    rayon::ThreadPoolBuilder::new()
        .num_threads(2)
        .build()
        .unwrap()
}

#[test]
fn ordinary_options_do_not_require_python_or_generation() {
    let args = ["--feature", "io.smiles.parse", "--dataset", "smoke"].map(str::to_owned);
    let normal = options(&args).unwrap();
    assert!(!normal.generate);
    assert!(normal.python.is_none());
    let mut args = vec!["generate".into()];
    args.extend(normal_args());
    let generate = options(&args).unwrap();
    assert!(generate.generate);
    assert_eq!(generate.python, Some("python".into()));
}

#[test]
fn repository_text_fingerprints_are_portable_but_input_bytes_remain_exact() {
    assert_eq!(
        stored::text_hash("first\r\nsecond\r\n"),
        stored::text_hash("first\nsecond\n")
    );
    assert_ne!(sha256(b"first\r\n"), sha256(b"first\n"));
}

#[test]
fn dssp_partner_identity_and_omega_are_required_observations() {
    use std::io::{BufRead, BufReader};
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("goldens/smoke/bio.secondary-structure.dssp.jsonl.gz");
    let reader = BufReader::new(flate2::read::GzDecoder::new(fs::File::open(path).unwrap()));
    let value = reader
        .lines()
        .map(|line| serde_json::from_str::<Value>(&line.unwrap()).unwrap())
        .find(|entry| entry["expected"]["status"] == "ok")
        .unwrap()["expected"]["value"]
        .clone();
    let expected = Outcome::Ok {
        value: value.clone(),
    };
    assert_eq!(
        comparison("bio.secondary-structure.dssp", &expected, &expected).0,
        "agrees"
    );
    let mut changed = value.clone();
    let partner = changed["residues"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .flat_map(|residue| residue["beta_partners"].as_array_mut().unwrap())
        .find(|partner| !partner.is_null())
        .unwrap();
    partner["partner_sequence_id"] = json!(99999);
    assert_eq!(
        comparison(
            "bio.secondary-structure.dssp",
            &expected,
            &Outcome::Ok { value: changed }
        )
        .0,
        "disagrees"
    );
    for field in ["omega_degrees", "beta_partners"] {
        let mut missing = value.clone();
        missing["residues"][0]
            .as_object_mut()
            .unwrap()
            .remove(field);
        assert_eq!(
            comparison(
                "bio.secondary-structure.dssp",
                &expected,
                &Outcome::Ok { value: missing }
            )
            .0,
            "error"
        );
    }
}
fn normal_args() -> Vec<String> {
    ["--feature", "io.smiles.parse", "--dataset", "smoke"]
        .map(str::to_owned)
        .into()
}

#[test]
fn stored_comparison_never_invokes_reference_and_counts_every_outcome() {
    let opts = opts(false);
    let fixture = Fixture::new();
    let mut stale = golden(&opts.feature, 2, success());
    stale["input_sha256"] = json!(sha256(b"changed"));
    let stored = fixture.store(
        &opts.feature,
        &[
            golden(&opts.feature, 0, success()),
            golden(
                &opts.feature,
                1,
                Outcome::Error {
                    message: "reference rejected source".into(),
                },
            ),
            stale,
        ],
    );
    let pool = pool();
    let progress = ProgressBar::hidden();
    progress.set_length(4);
    let mut results = Vec::new();
    let mut run = BatchRun {
        opts: &opts,
        pool: &pool,
        progress: &progress,
        mode: RunMode::Compare(stored),
        results: &mut results,
        reference: |_, _, _| panic!("ordinary comparison must not run a reference"),
        evaluate,
    };
    let mut row = Summary {
        dataset: "test".into(),
        feature: opts.feature.clone(),
        ..Default::default()
    };
    run.batch(&(0..4).map(case).collect::<Vec<_>>(), &mut row)
        .unwrap();
    assert_eq!(
        (row.cases, row.agrees, row.disagrees, row.errors),
        (4, 1, 0, 3)
    );
    assert_eq!(progress.position(), 4);
    drop(run);
    let output = String::from_utf8(results).unwrap();
    assert_eq!(output.lines().count(), 4);
    assert!(output.contains("checksum differs"));
    assert!(output.contains("missing stored golden"));
}

#[test]
fn generation_calls_reference_only_and_stores_its_errors() {
    let opts = opts(true);
    let fixture = Fixture::new();
    let pool = pool();
    let progress = ProgressBar::hidden();
    let mut results = Vec::new();
    let mut run = BatchRun {
        opts: &opts,
        pool: &pool,
        progress: &progress,
        mode: RunMode::Generate {
            writer: GeneratedGoldens::create(&fixture.path()).unwrap(),
            metadata: metadata(&opts.feature, 2),
        },
        results: &mut results,
        reference: |_, request, count| {
            assert_eq!(count, 2);
            assert!(request.get("written").is_none());
            Ok(ReferenceResponse {
                reference: reference(),
                time_ms: 1.0,
                results: vec![
                    success(),
                    Outcome::Error {
                        message: "reference error".into(),
                    },
                ],
            })
        },
        evaluate: |_, _| panic!("generation must not evaluate Kekule"),
    };
    let mut row = Summary {
        dataset: "test".into(),
        feature: opts.feature.clone(),
        ..Default::default()
    };
    run.batch(&[case(0), case(1)], &mut row).unwrap();
    let RunMode::Generate { writer, metadata } = run.mode else {
        panic!()
    };
    writer.finish(metadata).unwrap();
    assert_eq!((row.cases, row.errors, row.kekule_ms), (2, 1, 0.0));
    let mut stored = fixture.load(&opts.feature).unwrap();
    assert!(matches!(
        stored.expected(&case(0)).unwrap(),
        Outcome::Ok { .. }
    ));
    assert!(matches!(
        stored.expected(&case(1)).unwrap(),
        Outcome::Error { .. }
    ));
    stored.finish().unwrap();
}

#[test]
fn invalid_input_cannot_publish_a_reference_set() {
    let opts = opts(true);
    let fixture = Fixture::new();
    let pool = pool();
    let progress = ProgressBar::hidden();
    let mut results = Vec::new();
    let mut run = BatchRun {
        opts: &opts,
        pool: &pool,
        progress: &progress,
        mode: RunMode::Generate {
            writer: GeneratedGoldens::create(&fixture.path()).unwrap(),
            metadata: metadata(&opts.feature, 0),
        },
        results: &mut results,
        reference: |_, _, _| unreachable!(),
        evaluate: |_, _| unreachable!(),
    };
    assert!(run
        .input_error(
            &mut Summary::default(),
            "corrupt",
            Some("input.smi"),
            "checksum differs"
        )
        .is_err());
    drop(run);
    assert!(!fixture.path().exists());
    assert!(!fixture.path().with_extension("meta.json").exists());
}

#[test]
fn stream_rejects_duplicates_out_of_order_and_wrong_counts() {
    for entries in [
        vec![golden("io.smiles.parse", 0, success()); 2],
        vec![
            golden("io.smiles.parse", 1, success()),
            golden("io.smiles.parse", 0, success()),
        ],
    ] {
        let fixture = Fixture::new();
        let mut stored = fixture.store("io.smiles.parse", &entries);
        assert!(stored.finish().is_err());
    }
    let fixture = Fixture::new();
    let mut writer = GeneratedGoldens::create(&fixture.path()).unwrap();
    json_line(&mut writer, &golden("io.smiles.parse", 0, success())).unwrap();
    writer.finish(metadata("io.smiles.parse", 2)).unwrap();
    assert!(fixture.load("io.smiles.parse").unwrap().finish().is_err());
}

#[test]
fn invalid_metadata_checksums_and_observation_fields_are_rejected() {
    for (field, value) in [
        ("dataset", json!("other")),
        ("feature", json!("other")),
        ("schema", json!(0)),
        ("contract_sha256", json!("old")),
        ("input_lock_sha256", json!("old")),
        ("sha256", json!("wrong")),
    ] {
        let fixture = Fixture::new();
        drop(fixture.store(
            "io.smiles.parse",
            &[golden("io.smiles.parse", 0, success())],
        ));
        let path = fixture.path().with_extension("meta.json");
        let mut meta: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        meta[field] = value;
        fs::write(path, serde_json::to_vec(&meta).unwrap()).unwrap();
        assert!(fixture.load("io.smiles.parse").is_err(), "{field}");
    }
    for (field, value) in [
        ("dataset", json!("other")),
        ("input_sha256", json!("")),
        ("reference", Value::Null),
        ("extra", json!(true)),
        ("expected", json!({"status":"ok","value":{"records":[{}]}})),
    ] {
        let fixture = Fixture::new();
        let mut entry = golden("io.smiles.parse", 0, success());
        entry[field] = value;
        let mut writer = GeneratedGoldens::create(&fixture.path()).unwrap();
        json_line(&mut writer, &entry).unwrap();
        writer.finish(metadata("io.smiles.parse", 1)).unwrap();
        assert!(fixture.load("io.smiles.parse").is_err(), "{field}");
    }
}

#[test]
fn interrupted_generation_is_not_published_and_existing_goldens_are_not_replaced() {
    let fixture = Fixture::new();
    {
        let mut writer = GeneratedGoldens::create(&fixture.path()).unwrap();
        writer.write_all(b"unfinished").unwrap();
        assert!(!fixture.path().exists());
    }
    assert!(!fixture.path().exists());
    let stored = fixture.store(
        "io.smiles.parse",
        &[golden("io.smiles.parse", 0, success())],
    );
    drop(stored);
    let bytes = fs::read(fixture.path()).unwrap();
    assert!(GeneratedGoldens::create(&fixture.path()).is_err());
    assert_eq!(fs::read(fixture.path()).unwrap(), bytes);
}

fn identity() -> Outcome {
    Outcome::Ok {
        value: json!({"records":[{"record_index":0,"status":"ok","title":"","identity":"CC"}]}),
    }
}
#[test]
fn writer_validation_uses_only_emitted_text_and_rejects_version_drift() {
    for version_matches in [true, false] {
        let mut opts = opts(false);
        opts.feature = "io.smiles.write".into();
        let fixture = Fixture::new();
        let stored = fixture.store(&opts.feature, &[golden(&opts.feature, 0, identity())]);
        let pool = pool();
        let progress = ProgressBar::hidden();
        let mut results = Vec::new();
        let mut run = BatchRun {
            opts: &opts,
            pool: &pool,
            progress: &progress,
            mode: RunMode::Compare(stored),
            results: &mut results,
            evaluate,
            reference: if version_matches {
                |_, request, count| {
                    assert!(request.get("inputs").is_none());
                    assert_eq!(count, 1);
                    assert!(request["written"][0]["value"]["written"][0]["text"]
                        .as_str()
                        .unwrap()
                        .starts_with("CC"));
                    Ok(ReferenceResponse {
                        reference: reference(),
                        time_ms: 1.0,
                        results: vec![identity()],
                    })
                }
            } else {
                |_, _, _| {
                    Ok(ReferenceResponse {
                        reference: Reference {
                            tool: "rdkit".into(),
                            version: "different".into(),
                        },
                        time_ms: 1.0,
                        results: vec![identity()],
                    })
                }
            },
        };
        let mut row = Summary {
            dataset: "test".into(),
            feature: opts.feature.clone(),
            ..Default::default()
        };
        let result = run.batch(&[case(0)], &mut row);
        if version_matches {
            result.unwrap();
            assert_eq!((row.cases, row.agrees, row.errors), (1, 1, 0));
        } else {
            assert!(result.unwrap_err().to_string().contains("version differs"));
        }
    }
}

#[test]
fn writer_reader_failure_is_not_misattributed_to_kekule_execution() {
    let mut opts = opts(false);
    opts.feature = "io.smiles.write".into();
    let fixture = Fixture::new();
    let stored = fixture.store(&opts.feature, &[golden(&opts.feature, 0, identity())]);
    let pool = pool();
    let progress = ProgressBar::hidden();
    let mut results = Vec::new();
    let mut run = BatchRun {
        opts: &opts,
        pool: &pool,
        progress: &progress,
        mode: RunMode::Compare(stored),
        results: &mut results,
        evaluate,
        reference: |_, _, _| Err(boxed_error("reference interpreter unavailable")),
    };
    let mut row = Summary {
        dataset: "test".into(),
        feature: opts.feature.clone(),
        ..Default::default()
    };
    run.batch(&[case(0)], &mut row).unwrap();
    assert_eq!(
        (row.errors, row.kekule_errors, row.writer_validation_errors),
        (1, 0, 1)
    );
}

#[test]
fn panics_are_counted_and_missing_formats_do_not_inflate_cases() {
    let opts = opts(false);
    let fixture = Fixture::new();
    let stored = fixture.store(&opts.feature, &[golden(&opts.feature, 0, success())]);
    let pool = pool();
    let progress = ProgressBar::hidden();
    let mut results = Vec::new();
    let mut run = BatchRun {
        opts: &opts,
        pool: &pool,
        progress: &progress,
        mode: RunMode::Compare(stored),
        results: &mut results,
        reference: |_, _, _| unreachable!(),
        evaluate: |_, _| panic!("injected defect"),
    };
    let mut row = Summary {
        dataset: "test".into(),
        feature: opts.feature.clone(),
        ..Default::default()
    };
    run.batch(&[case(0)], &mut row).unwrap();
    run.input_error(&mut row, "absent", None, "no source format")
        .unwrap();
    run.input_error(&mut row, "corrupt", Some("corrupt.smi"), "checksum")
        .unwrap();
    assert_eq!(
        (
            row.cases,
            row.errors,
            row.not_applicable,
            row.input_errors,
            row.kekule_errors
        ),
        (2, 2, 1, 1, 1)
    );
    drop(run);
    assert!(String::from_utf8(results)
        .unwrap()
        .contains("injected defect"));
}

#[test]
fn progress_includes_multiple_formats_and_missing_inputs() {
    let dataset = Dataset {
        root: PathBuf::new(),
        lock: json!({"entries":[{"id":"a","files":[{"path":"a.smi"},{"path":"a.sdf"}]},{"id":"b","files":[]}]}),
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
fn a_crashing_reference_cannot_erase_the_rest_of_the_batch() {
    let request = json!({"feature":"io.smiles.parse","inputs":[{"text":"crash"},{"text":"CC"}]});
    let response = isolate_reference(&request, 2, |request, count| {
        if count == 2 || request["inputs"][0]["text"] == "crash" {
            return Err(boxed_error("injected process failure"));
        }
        if count == 0 {
            assert_eq!(request["describe"], true);
        }
        Ok(ReferenceResponse {
            reference: reference(),
            time_ms: 0.0,
            results: if count == 0 { vec![] } else { vec![success()] },
        })
    })
    .unwrap();
    assert_eq!(response.results.len(), 2);
    assert!(
        matches!(&response.results[0],Outcome::Error{message} if message.contains("injected process failure"))
    );
    assert!(matches!(&response.results[1], Outcome::Ok { .. }));
}

#[test]
fn mol_model_writer_has_an_explicit_title_contract_without_mutating_the_source() {
    let input = Input {
        path: "input.sdf".into(),
        text: crate::tests::simple_sdf_record("external title"),
    };
    let original = evaluate("io.mol.parse", &input).unwrap();
    let expected = Outcome::Ok {
        value: original.clone(),
    };
    let mut emitted = original.clone();
    emitted["records"][0]["title"] = json!("");
    let actual = Outcome::Ok { value: emitted };
    assert_eq!(
        comparison("io.mol.v2000.write", &expected, &actual).0,
        "agrees"
    );
    assert_eq!(
        comparison("io.mol.parse", &expected, &actual).0,
        "disagrees"
    );
    assert!(matches!(expected,Outcome::Ok{value} if value==original));
}
