//! Externally supplied molecules shared by every invariant.
//!
//! Sources are the bundled benchmark smoke and RDKit structure corpora, the
//! CIP validation fixture, and the perception regression fixtures. None of
//! these is constructed for a particular test.

use std::path::Path;

use kekule::core::{Molecule, StereoElementKind};
use kekule::structure::Model;
use kekule::{molfile, sdf, smiles};

use crate::support::{aromaticity_cases, corpus, fixture};

/// One connected molecule, interpreted but not perceived.
pub struct Sample {
    pub label: String,
    pub molecule: Molecule,
    /// The source model when the molecule came from a single-component
    /// Molfile or SDF record, retaining its drawing or 3D coordinates.
    pub model: Option<Model>,
}

impl Sample {
    pub fn perceived(&self) -> Molecule {
        let mut molecule = self.molecule.clone();
        molecule
            .perceive()
            .unwrap_or_else(|error| panic!("{}: {error}", self.label));
        molecule
    }

    /// SMILES cannot encode atropisomeric axes yet; writers must reject them.
    pub fn has_axis_stereo(&self) -> bool {
        self.molecule
            .stereo_elements()
            .any(|(_, element)| matches!(element.kind, StereoElementKind::Axis(_)))
    }
}

/// Source records that are not molecules Kekule can represent, with the
/// reason. Each must keep failing interpretation; if one starts to succeed,
/// move it into the corpus.
pub const EXPECTED_REJECTIONS: &[(&str, &str)] = &[
    (
        "rdkit-structures/data/chebi_15469.v3k.mol",
        "R-group query atom",
    ),
    (
        "rdkit-structures/data/chebi_57262.v3k.mol",
        "R-group query atom",
    ),
    (
        "rdkit-structures/data/chebi_57262.v3k.2.mol",
        "R-group query atom",
    ),
    (
        "rdkit-structures/data/atropisomers/BMS-986142_atropBad2.mol",
        "V3000 atom configuration with fewer than four tetrahedral carriers",
    ),
];

/// A corpus sample that exposes a known, unfixed defect in one invariant.
pub struct KnownDefect {
    pub invariant: &'static str,
    pub sample: &'static str,
    /// Text that identifies the expected failure message. Only matching
    /// failures are exempt; any other failure of the sample is still reported.
    pub failure: &'static str,
    pub defect: &'static str,
}

const EXPLICIT_H_CANONICAL: &str = "canonical SMILES depends on explicit hydrogen atoms";
const EXPLICIT_H_WRITER: &str =
    "canonical writer demands hydrogen perception for a perceived explicit-H aromatic NH";

/// Each entry must keep failing exactly as described; when a fix makes it
/// pass, the invariant reports the entry so it can be deleted.
pub const KNOWN_DEFECTS: &[KnownDefect] = &[
    KnownDefect {
        invariant: "smiles_round_trip",
        sample: "molfile:cid_2244.sdf#0",
        failure: "left: \"CC(=O)Oc1ccccc1C(=O)O\"\n right: \"CC(=O)Oc1ccccc1C(O)=O\"",
        defect: EXPLICIT_H_CANONICAL,
    },
    KnownDefect {
        invariant: "hydrogen_round_trip",
        sample: "molfile:cid_2244.sdf#0",
        failure: r#"collapsed canonical Some("CC(=O)Oc1ccccc1C(=O)O") != Some("CC(=O)Oc1ccccc1C(O)=O")"#,
        defect: EXPLICIT_H_CANONICAL,
    },
    KnownDefect {
        invariant: "smiles_round_trip",
        sample: "molfile:pubchem-10250.sdf#0",
        failure: "canonical writer failed: canonical SMILES bracket atom a2 requires installed \
                  hydrogen perception",
        defect: EXPLICIT_H_WRITER,
    },
    KnownDefect {
        invariant: "hydrogen_round_trip",
        sample: "molfile:pubchem-10250.sdf#0",
        failure: r#"collapsed canonical Some("O=c1[nH]c(=O)c2nccnc2[nH]1") != None"#,
        defect: EXPLICIT_H_WRITER,
    },
];

/// Fails on every (sample label, message) failure that does not match a known
/// defect of `invariant`, and on every known defect of `invariant` that no
/// longer produces its expected failure.
pub fn assert_invariant(invariant: &str, failures: Vec<(String, String)>) {
    let known = KNOWN_DEFECTS
        .iter()
        .filter(|defect| defect.invariant == invariant)
        .collect::<Vec<_>>();
    let expected = |label: &str, message: &str, defect: &KnownDefect| {
        defect.sample == label && message.contains(defect.failure)
    };
    let mut report = failures
        .iter()
        .filter(|(label, message)| !known.iter().any(|defect| expected(label, message, defect)))
        .map(|(label, message)| format!("{label}: {message}"))
        .collect::<Vec<_>>();
    for defect in known {
        if !failures
            .iter()
            .any(|(label, message)| expected(label, message, defect))
        {
            report.push(format!(
                "{}: known defect ({}) no longer fails as expected; remove it from KNOWN_DEFECTS",
                defect.sample, defect.defect
            ));
        }
    }
    assert!(report.is_empty(), "{}", report.join("\n"));
}

pub fn molecules() -> Vec<Sample> {
    let mut samples = Vec::new();
    // SDF regression records are read below with their source models.
    for case in aromaticity_cases()
        .into_iter()
        .filter(|case| !case.input.starts_with("sdf:"))
    {
        samples.push(Sample {
            label: format!("aromaticity:{}", case.label),
            molecule: case.molecule(),
            model: None,
        });
    }
    for (label, components) in regression_mixtures() {
        for (index, molecule) in components.into_iter().enumerate() {
            samples.push(Sample {
                label: format!("regression:{label}#{index}"),
                molecule,
                model: None,
            });
        }
    }
    for path in files(&corpus("smoke/data/pubchem_smiles"), "txt") {
        let text = std::fs::read_to_string(&path).unwrap();
        for (index, molecule) in smiles::to_molecules(text.trim())
            .unwrap_or_else(|error| panic!("{}: {error}", path.display()))
            .into_iter()
            .enumerate()
        {
            samples.push(Sample {
                label: format!("smoke:{}#{index}", name(&path)),
                molecule,
                model: None,
            });
        }
    }
    let molfiles = [
        (corpus("smoke/data/pubchem"), "sdf"),
        (corpus("smoke/data/rdkit_atropisomers"), "mol"),
        (corpus("rdkit-structures/data"), "mol"),
        (corpus("rdkit-structures/data/atropisomers"), "mol"),
        (fixture("cip"), "sdf"),
        (fixture("perception"), "sdf"),
    ];
    for (directory, extension) in molfiles {
        for path in files(&directory, extension) {
            if EXPECTED_REJECTIONS
                .iter()
                .any(|(rejected, _)| path.ends_with(rejected))
            {
                continue;
            }
            let model =
                read_model(&path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
            let molecules = model
                .topology()
                .molecules()
                .map(|instance| instance.molecule().clone())
                .collect::<Vec<_>>();
            let single = molecules.len() == 1;
            for (index, molecule) in molecules.into_iter().enumerate() {
                samples.push(Sample {
                    label: format!("molfile:{}#{index}", name(&path)),
                    molecule,
                    model: single.then(|| model.clone()),
                });
            }
        }
    }
    samples
}

/// Multi-component regression records, interpreted component by component.
pub fn regression_mixtures() -> Vec<(String, Vec<Molecule>)> {
    let path = fixture("perception/regressions.smi");
    std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()))
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let (label, input) = line.split_once('\t').expect("label<TAB>SMILES");
            let molecules =
                smiles::to_molecules(input).unwrap_or_else(|error| panic!("{label}: {error}"));
            (label.to_owned(), molecules)
        })
        .collect()
}

/// Reads a single-record Molfile or the first SDF record as a model.
pub fn read_model(path: &Path) -> Result<Model, String> {
    let text = std::fs::read_to_string(path).map_err(|error| error.to_string())?;
    if path.extension().is_some_and(|extension| extension == "sdf") {
        let document = sdf::parse_str(&text).map_err(|error| error.to_string())?;
        document.records()[0]
            .to_model()
            .map_err(|error| error.to_string())
    } else {
        molfile::parse_str(&text)
            .map_err(|error| error.to_string())?
            .to_model()
            .map_err(|error| error.to_string())
    }
}

fn files(directory: &Path, extension: &str) -> Vec<std::path::PathBuf> {
    let mut paths = std::fs::read_dir(directory)
        .unwrap_or_else(|error| panic!("{}: {error}", directory.display()))
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.extension().is_some_and(|value| value == extension))
        .collect::<Vec<_>>();
    paths.sort();
    paths
}

fn name(path: &Path) -> String {
    path.file_name().unwrap().to_string_lossy().into_owned()
}

#[test]
fn corpus_reads_every_source_and_keeps_rejections_explicit() {
    let samples = molecules();
    for prefix in ["aromaticity:", "regression:", "smoke:", "molfile:"] {
        let count = samples
            .iter()
            .filter(|sample| sample.label.starts_with(prefix))
            .count();
        assert!(count >= 8, "only {count} samples from {prefix}");
    }
    assert!(
        samples
            .iter()
            .filter(|sample| sample.model.is_some())
            .count()
            >= 60
    );
    assert!(samples.iter().any(Sample::has_axis_stereo));
    for (path, reason) in EXPECTED_REJECTIONS {
        assert!(
            read_model(&corpus(path)).is_err(),
            "{path} ({reason}) is now accepted; move it into the corpus"
        );
    }
}
