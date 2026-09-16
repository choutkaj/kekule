use std::fs;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

#[test]
fn unavailable_goldens_stop_the_cli_without_reporting_case_failures() {
    let root = std::env::temp_dir().join(format!(
        "kekule-missing-goldens-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&root).unwrap();
    let report = root.join("report.json");
    let cases = report.with_extension("cases.jsonl");
    let golden = root.join("smoke").join("io.smiles.parse.jsonl.gz");
    for missing in [true, false] {
        if !missing {
            fs::create_dir(root.join("smoke")).unwrap();
            fs::write(&golden, b"invalid gzip file").unwrap();
        }
        let output = Command::new(env!("CARGO_BIN_EXE_kekule-bench"))
            .args([
                "--feature",
                "io.smiles.parse",
                "--dataset",
                "smoke",
                "--limit",
                "1",
                "--goldens",
            ])
            .arg(&root)
            .arg("--output")
            .arg(&report)
            .output()
            .unwrap();
        assert!(!output.status.success());
        let stderr = String::from_utf8(output.stderr).unwrap();
        let stdout = String::from_utf8(output.stdout).unwrap();
        assert!(stderr.contains(golden.to_str().unwrap()), "{stderr}");
        if missing {
            assert!(
                stderr
                    .contains("cargo benchmark generate --feature io.smiles.parse --dataset smoke"),
                "{stderr}"
            );
        } else {
            assert!(stderr.contains("cannot load stored goldens"), "{stderr}");
            assert!(!stderr.contains("Generate them once"), "{stderr}");
        }
        assert!(!stdout.contains("agree"), "{stdout}");
        assert!(!stdout.contains("errors"), "{stdout}");
        assert!(fs::read(&cases).unwrap().is_empty());
        let summary: serde_json::Value =
            serde_json::from_slice(&fs::read(&report).unwrap()).unwrap();
        assert_eq!(summary["complete"], false);
        assert_eq!(summary["passed"], false);
        assert!(summary["error"]
            .as_str()
            .unwrap()
            .contains("cannot load stored goldens"));
        assert_eq!(summary["results"], serde_json::json!([]));
        assert_eq!(golden.exists(), !missing);
        fs::remove_file(&report).unwrap();
        fs::remove_file(&cases).unwrap();
    }
    fs::remove_file(golden).unwrap();
    fs::remove_dir(root.join("smoke")).unwrap();
    fs::remove_dir(root).unwrap();
}
