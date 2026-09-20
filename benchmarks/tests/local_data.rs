use std::fs;
use std::path::Path;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

fn git(root: &Path, args: &[&str]) -> String {
    let output = Command::new("git")
        .current_dir(root)
        .args(args)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "git {args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    String::from_utf8(output.stdout).unwrap()
}

#[test]
fn staging_keeps_bulk_data_local_and_includes_bundled_inputs_and_provenance() {
    let root = std::env::temp_dir().join(format!(
        "kekule-local-data-{}-{}",
        std::process::id(),
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir(&root).unwrap();
    fs::copy(
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../.gitignore"),
        root.join(".gitignore"),
    )
    .unwrap();
    let mut expected = vec![".gitignore".to_owned()];
    for dataset in [
        "smoke",
        "rdkit-queries",
        "rdkit-structures",
        "pubchem-100k",
        "future-dataset",
    ] {
        let bundled = matches!(dataset, "smoke" | "rdkit-queries" | "rdkit-structures");
        for (suffix, tracked) in [
            ("data/input.sdf", bundled),
            ("data/nested/input.sdf", bundled),
            ("sources.lock.json", true),
        ] {
            let name = format!("benchmarks/corpora/{dataset}/{suffix}");
            let path = root.join(&name);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(path, b"fixture").unwrap();
            if tracked {
                expected.push(name);
            }
        }
        for (suffix, tracked) in [
            ("io.sdf.parse.jsonl.gz", bundled),
            ("io.sdf.parse.jsonl.meta.json", true),
            ("interrupted.tmp", false),
        ] {
            let name = format!("benchmarks/goldens/{dataset}/{suffix}");
            let path = root.join(&name);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(path, b"fixture").unwrap();
            if tracked {
                expected.push(name);
            }
        }
    }
    for dataset in ["smoke", "pubchem-100k"] {
        for (suffix, tracked) in [("jsonl.meta.json", true), ("jsonl.gz", dataset == "smoke")] {
            let name = format!(
                "benchmarks/goldens/legacy/query-smarts-v1/{dataset}/query.smarts.{suffix}"
            );
            let path = root.join(&name);
            fs::create_dir_all(path.parent().unwrap()).unwrap();
            fs::write(path, b"archived reference").unwrap();
            if tracked {
                expected.push(name);
            }
        }
    }
    fs::write(root.join("benchmarks/goldens/bundle.tar.gz"), b"fixture").unwrap();
    fs::create_dir(root.join("benchmarks/runs")).unwrap();
    for name in [
        "run.json",
        "run.cases.jsonl",
        "run.cases.jsonl.gz",
        "index.html",
        "dashboard-data.js",
    ] {
        fs::write(root.join("benchmarks/runs").join(name), b"local run").unwrap();
    }
    fs::write(
        root.join("benchmarks/.dashboard-python"),
        b"local executable",
    )
    .unwrap();
    git(&root, &["init", "--quiet"]);
    git(&root, &["add", "."]);
    let staged = git(&root, &["ls-files"]);
    let mut actual: Vec<_> = staged.lines().map(str::to_owned).collect();
    actual.sort();
    expected.sort();
    assert_eq!(actual, expected);
    assert!(root
        .join("benchmarks/goldens/pubchem-100k/io.sdf.parse.jsonl.gz")
        .exists());
    fs::remove_dir_all(root).unwrap();
}
