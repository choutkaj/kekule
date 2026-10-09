//! Fixture access shared by the integration test binaries.
#![allow(dead_code)]

use std::path::{Path, PathBuf};

use kekule::core::Molecule;
use kekule::{sdf, smiles};

/// Absolute path of a file under `tests/fixtures`.
pub fn fixture(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(relative)
}

/// Absolute path of a file under the repository's `benchmarks/corpora`.
pub fn corpus(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../benchmarks/corpora")
        .join(relative)
}

/// Interprets one connected molecule from SMILES or from the first record of
/// an SDF file, written as `sdf:<path>` relative to `directory`. Perception is
/// not run.
pub fn read_molecule(input: &str, directory: &Path) -> Molecule {
    let mut molecules = match input.strip_prefix("sdf:") {
        Some(file) => {
            let text = std::fs::read_to_string(directory.join(file))
                .unwrap_or_else(|error| panic!("{file}: {error}"));
            let document = sdf::parse_str(&text).unwrap_or_else(|error| panic!("{file}: {error}"));
            document.records()[0]
                .interpret()
                .unwrap_or_else(|error| panic!("{file}: {error}"))
                .into_molecules()
        }
        None => smiles::to_molecules(input).unwrap_or_else(|error| panic!("{input}: {error}")),
    };
    assert_eq!(molecules.len(), 1, "{input} must be one connected molecule");
    molecules.pop().unwrap()
}

/// One row of `fixtures/perception/aromaticity.tsv`.
pub struct AromaticityCase {
    pub label: String,
    pub input: String,
    pub aromatic_atoms: Vec<usize>,
    pub nonaromatic_bonds: Vec<(usize, usize)>,
}

impl AromaticityCase {
    pub fn molecule(&self) -> Molecule {
        read_molecule(&self.input, &fixture("perception"))
    }
}

pub fn aromaticity_cases() -> Vec<AromaticityCase> {
    let path = fixture("perception/aromaticity.tsv");
    std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()))
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            assert_eq!(fields.len(), 4, "malformed aromaticity row: {line}");
            AromaticityCase {
                label: fields[0].to_owned(),
                input: fields[1].to_owned(),
                aromatic_atoms: index_ranges(fields[2]),
                nonaromatic_bonds: index_pairs(fields[3]),
            }
        })
        .collect()
}

/// Parses `0-5,7` into `[0, 1, 2, 3, 4, 5, 7]`; `-` is empty.
fn index_ranges(field: &str) -> Vec<usize> {
    if field == "-" {
        return Vec::new();
    }
    field
        .split(',')
        .flat_map(|part| {
            let (start, end) = part.split_once('-').unwrap_or((part, part));
            start.parse::<usize>().unwrap()..=end.parse::<usize>().unwrap()
        })
        .collect()
}

/// Parses `3:12,4:6` into `[(3, 12), (4, 6)]`; `-` is empty.
fn index_pairs(field: &str) -> Vec<(usize, usize)> {
    if field == "-" {
        return Vec::new();
    }
    field
        .split(',')
        .map(|pair| {
            let (left, right) = pair.split_once(':').unwrap();
            (left.parse().unwrap(), right.parse().unwrap())
        })
        .collect()
}
