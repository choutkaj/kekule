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

/// Samples that expose a known, unfixed defect, as (invariant, sample label,
/// defect). An entry must keep failing its invariant; when a fix makes it
/// pass, the invariant reports the entry so it can be deleted.
pub const KNOWN_DEFECTS: &[(&str, &str, &str)] = &[
    (
        "smiles_round_trip",
        "molfile:cid_2244.sdf#0",
        "canonical SMILES depends on explicit hydrogen atoms",
    ),
    (
        "hydrogen_round_trip",
        "molfile:cid_2244.sdf#0",
        "canonical SMILES depends on explicit hydrogen atoms",
    ),
    (
        "smiles_round_trip",
        "molfile:pubchem-10250.sdf#0",
        "canonical writer demands hydrogen perception for a perceived explicit-H aromatic NH",
    ),
    (
        "hydrogen_round_trip",
        "molfile:pubchem-10250.sdf#0",
        "canonical writer demands hydrogen perception for a perceived explicit-H aromatic NH",
    ),
    (
        "hydrogen_round_trip",
        "molfile:JDQ443_3d.mol#0",
        "CIP leaves an atropisomeric axis unresolved after explicit hydrogens are collapsed",
    ),
    (
        "hydrogen_round_trip",
        "molfile:Sotorasib_3d.mol#0",
        "CIP leaves an atropisomeric axis unresolved after explicit hydrogens are collapsed",
    ),
];

/// Fails on every (sample label, message) failure that is not a known defect
/// of `invariant`, and on every known defect of `invariant` that now passes.
pub fn assert_invariant(invariant: &str, failures: Vec<(String, String)>) {
    let known = KNOWN_DEFECTS
        .iter()
        .filter(|(name, ..)| *name == invariant)
        .collect::<Vec<_>>();
    let mut report = failures
        .iter()
        .filter(|(label, _)| !known.iter().any(|(_, known, _)| known == label))
        .map(|(label, message)| format!("{label}: {message}"))
        .collect::<Vec<_>>();
    for (_, label, defect) in known {
        if !failures.iter().any(|(failed, _)| failed == label) {
            report.push(format!(
                "{label}: known defect ({defect}) no longer reproduces; remove it from KNOWN_DEFECTS"
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
