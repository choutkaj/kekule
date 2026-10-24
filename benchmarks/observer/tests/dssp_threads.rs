//! Parallel DSSP acceptance: on every curated benchmark entry, and on the
//! `pdb-1000` corpus when it is present, assignment on several threads
//! returns exactly the single-threaded result, or the same error.
//!
//! `pdb-1000` is read from `benchmarks/corpora/pdb-1000/data`, or from the
//! directory `KEKULE_PDB_1000` names; it is not part of the repository.

use std::io::Read;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use kekule::dssp::{self, DsspOptions};
use kekule::macromolecular::{ModelSelection, StructureInterpretOptions};
use kekule::mmcif;

const THREADS: [usize; 3] = [2, 3, 8];

fn entries(root: &Path, suffix: &str) -> Vec<PathBuf> {
    let mut files = std::fs::read_dir(root)
        .expect("readable dataset directory")
        .map(|entry| entry.expect("dataset entry").path())
        .filter(|path| path.to_string_lossy().ends_with(suffix))
        .collect::<Vec<_>>();
    files.sort();
    files
}

/// Compares every thread count with one thread; returns whether the entry
/// was analyzed. Entries Kekule does not interpret are skipped.
fn check(path: &Path, text: &str) -> bool {
    let Ok(document) = mmcif::parse_str(text, Default::default()) else {
        return false;
    };
    let Ok(interpretation) = document.interpret(
        StructureInterpretOptions::default().with_model_selection(ModelSelection::First),
    ) else {
        return false;
    };
    let model = interpretation.into_model();
    let serial = dssp::assign(&model, DsspOptions::default());
    for threads in THREADS {
        let options =
            DsspOptions::default().with_threads(NonZeroUsize::new(threads).expect("nonzero"));
        assert_eq!(
            dssp::assign(&model, options),
            serial,
            "{}: {threads} threads",
            path.display()
        );
    }
    serial.is_ok()
}

#[test]
fn dssp_threads_match_one_thread_on_the_curated_entries() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../datasets/bio/curated/pdb");
    let mut analyzed = 0;
    for path in entries(&root, ".cif.gz") {
        let mut text = String::new();
        flate2::read::GzDecoder::new(std::fs::File::open(&path).expect("dataset file"))
            .read_to_string(&mut text)
            .expect("gzip mmCIF");
        analyzed += usize::from(check(&path, &text));
    }
    assert!(analyzed >= 15, "{analyzed} curated entries analyzed");
}

#[test]
#[ignore = "needs the pdb-1000 corpus; run in release with --ignored"]
fn dssp_threads_match_one_thread_on_pdb_1000() {
    let root = std::env::var_os("KEKULE_PDB_1000").map_or_else(
        || Path::new(env!("CARGO_MANIFEST_DIR")).join("../corpora/pdb-1000/data"),
        PathBuf::from,
    );
    if !root.is_dir() {
        eprintln!("skipped: no pdb-1000 corpus at {}", root.display());
        return;
    }
    let mut analyzed = 0;
    for path in entries(&root, ".cif") {
        let text = std::fs::read_to_string(&path).expect("corpus file");
        analyzed += usize::from(check(&path, &text));
    }
    eprintln!("{analyzed} pdb-1000 entries analyzed");
    assert!(analyzed > 0);
}
