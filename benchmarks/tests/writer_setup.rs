use serde_json::Value;
use std::{
    fs,
    io::Read,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};

#[test]
fn unavailable_writer_reader_stops_before_measuring_cases() {
    let root = std::env::temp_dir().join(format!(
        "kekule-writer-setup-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir_all(&root).unwrap();
    let report_path = root.join("report.json");
    let output = Command::new(env!("CARGO_BIN_EXE_kekule-bench"))
        .env("KEKULE_BENCHMARK_RUNS_DIR", root.join("runs"))
        .env("KEKULE_DASHBOARD_PYTHON", root.join("no-python"))
        .args([
            "--feature",
            "io.smiles.write",
            "--dataset",
            "smoke",
            "--writer-python",
        ])
        .arg(root.join("no-python"))
        .arg("--output")
        .arg(&report_path)
        .output()
        .unwrap();
    assert!(!output.status.success());
    let report: Value = serde_json::from_slice(&fs::read(&report_path).unwrap()).unwrap();
    assert_eq!(report["complete"], false);
    assert_eq!(report["passed"], false);
    assert!(report["error"]
        .as_str()
        .unwrap()
        .contains("writer reader unavailable"));
    assert!(report["results"].as_array().unwrap().is_empty());
    let mut cases = String::new();
    flate2::read::GzDecoder::new(
        fs::File::open(report_path.with_extension("cases.jsonl.gz")).unwrap(),
    )
    .read_to_string(&mut cases)
    .unwrap();
    assert!(cases.is_empty());
    fs::remove_dir_all(root).unwrap();
}
