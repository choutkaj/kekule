use serde_json::Value;
use std::{
    fs,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};

#[test]
fn unavailable_dashboard_python_preserves_reports_and_the_scientific_status() {
    let root = std::env::temp_dir().join(format!(
        "kekule-dashboard-unavailable-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let output = Command::new(env!("CARGO_BIN_EXE_kekule-bench"))
        .env("KEKULE_BENCHMARK_RUNS_DIR", &root)
        .env("KEKULE_DASHBOARD_PYTHON", root.join("no-python"))
        .args(["--feature", "io.smiles.parse", "--dataset", "smoke"])
        .output()
        .unwrap();
    assert!(output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("dashboard was not refreshed"));
    let reports: Vec<_> = fs::read_dir(&root)
        .unwrap()
        .map(|item| item.unwrap().path())
        .filter(|path| path.extension().is_some_and(|ext| ext == "json"))
        .collect();
    assert_eq!(reports.len(), 1);
    let report: Value = serde_json::from_slice(&fs::read(&reports[0]).unwrap()).unwrap();
    assert_eq!(report["passed"], true);
    assert!(reports[0].with_extension("cases.jsonl").exists());
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn comparisons_archive_and_refresh_even_when_they_disagree() {
    let root = std::env::temp_dir().join(format!(
        "kekule-runs-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    let runs = root.join("runs");
    fs::create_dir_all(&runs).unwrap();
    for (feature, passing) in [("io.smiles.parse", true), ("io.sdf.parse", false)] {
        let mut command = Command::new(env!("CARGO_BIN_EXE_kekule-bench"));
        command.env("KEKULE_BENCHMARK_RUNS_DIR", &runs).args([
            "--feature",
            feature,
            "--dataset",
            "smoke",
        ]);
        // A non-JSON output name is still discoverable through its archived summary.
        if !passing {
            command.arg("--output").arg(runs.join("custom.report"));
        }
        let output = command.output().unwrap();
        let stderr = String::from_utf8_lossy(&output.stderr);
        assert_eq!(output.status.success(), passing, "{stderr}");
        assert!(!stderr.contains("dashboard was not refreshed"), "{stderr}");
        let script = fs::read_to_string(runs.join("dashboard-data.js")).unwrap();
        let data: Value = serde_json::from_str(
            script
                .strip_prefix("window.updateKekuleBenchmarks(")
                .unwrap()
                .trim()
                .strip_suffix(");")
                .unwrap(),
        )
        .unwrap();
        let history = data["runs"].as_array().unwrap();
        assert_eq!(history.len(), if passing { 1 } else { 2 });
        assert_eq!(history[0]["complete"], true);
        assert_eq!(history[0]["passed"], passing);
        assert_eq!(history[0]["results"][0]["feature"], feature);
        assert!(history[0]["started_at_unix_ms"].as_u64().unwrap() > 0);
        assert!(runs.join("index.html").exists());
        if !passing {
            assert!(runs.join("custom.cases.jsonl").exists());
            assert!(
                history[0]["started_at_unix_ms"].as_u64()
                    >= history[1]["started_at_unix_ms"].as_u64()
            );
        }
    }
    fs::remove_dir_all(root).unwrap();
}
