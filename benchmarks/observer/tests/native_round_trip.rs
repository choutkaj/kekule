//! Native serialization acceptance: every curated benchmark input that Kekule
//! interprets round-trips through the native format exactly.
//!
//! Small molecules (Molfile, SDF, SMILES) and macromolecular entries (mmCIF,
//! as models and, for several coordinate models, ensembles) are checked for
//! molecular and layout equality, properties, classifications, perception,
//! and bitwise-identical coordinates. Inputs Kekule does not interpret are
//! the benchmark's business, not this test's, and are skipped.

use std::io::Read;
use std::path::{Path, PathBuf};

use kekule::geometry::Point3;
use kekule::molecule::Molecule;
use kekule::native::{self, NativeObject, NativeReadOptions, NativeSource, NativeWriteOptions};
use kekule::structure::{Ensemble, Model};
use kekule::topology::Topology;
use kekule::{mmcif, molfile, sdf, smiles};

fn dataset(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../datasets")
        .join(name)
        .join("curated")
}

fn files(root: &Path) -> Vec<PathBuf> {
    let mut found = Vec::new();
    let mut pending = vec![root.to_path_buf()];
    while let Some(directory) = pending.pop() {
        for entry in std::fs::read_dir(&directory).expect("readable dataset directory") {
            let path = entry.expect("dataset entry").path();
            if path.is_dir() {
                pending.push(path);
            } else {
                found.push(path);
            }
        }
    }
    found.sort();
    found
}

/// Encodes, decodes, and checks that the decoded object encodes to the same
/// bytes.
fn round_trip<'a>(source: impl Into<NativeSource<'a>>) -> NativeObject {
    let bytes = native::write(source, NativeWriteOptions::default()).expect("encodes");
    let decoded = native::read(bytes.as_slice(), NativeReadOptions::default()).expect("decodes");
    let again = match &decoded {
        NativeObject::Molecule(value) => {
            native::write(value.as_ref(), NativeWriteOptions::default())
        }
        NativeObject::Topology(value) => native::write(value, NativeWriteOptions::default()),
        NativeObject::Model(value) => native::write(value, NativeWriteOptions::default()),
        NativeObject::Ensemble(value) => native::write(value, NativeWriteOptions::default()),
        other => panic!("unexpected {:?}", other.kind()),
    }
    .expect("re-encodes");
    assert!(again == bytes, "re-encoding changed the bytes");
    decoded
}

fn assert_same_molecule(decoded: &Molecule, source: &Molecule) {
    assert_eq!(decoded, source);
    assert_eq!(decoded.perception(), source.perception());
    assert_eq!(decoded.properties(), source.properties());
}

fn assert_same_topology(decoded: &Topology, source: &Topology) {
    assert!(decoded.same_layout(source));
    assert!(!decoded.shares_layout(source));
    assert_eq!(decoded.atom_ids(), source.atom_ids());
    assert_eq!(decoded.properties(), source.properties());
    for (left, right) in decoded.definitions().zip(source.definitions()) {
        assert_same_molecule(left.molecule(), right.molecule());
    }
    for (left, right) in decoded.instances().zip(source.instances()) {
        assert_eq!(left.class_assignment(), right.class_assignment());
    }
    for (left, right) in decoded.residues().zip(source.residues()) {
        assert_eq!(left.class_assignment(), right.class_assignment());
        assert_eq!(
            (left.name(), left.author_seq_id(), left.insertion_code()),
            (right.name(), right.author_seq_id(), right.insertion_code())
        );
        assert_eq!(left.chain().author_id(), right.chain().author_id());
    }
}

fn coordinate_bits(points: &[Point3]) -> Vec<[u64; 3]> {
    points
        .iter()
        .map(|point| [point.x.to_bits(), point.y.to_bits(), point.z.to_bits()])
        .collect()
}

fn assert_same_model(decoded: &Model, source: &Model) {
    assert_same_topology(decoded.topology(), source.topology());
    assert_eq!(
        coordinate_bits(decoded.positions().values().value()),
        coordinate_bits(source.positions().values().value())
    );
    assert_eq!(decoded.cell(), source.cell());
    assert_eq!(decoded.occupancies(), source.occupancies());
    assert_eq!(decoded.b_factors(), source.b_factors());
    assert_eq!(decoded.properties(), source.properties());
}

fn assert_same_ensemble(decoded: &Ensemble, source: &Ensemble) {
    assert_same_topology(decoded.topology(), source.topology());
    assert_eq!(decoded.len(), source.len());
    for (left, right) in decoded.iter().zip(source.iter()) {
        assert_eq!(left.payload(), right.payload());
        assert_eq!(
            coordinate_bits(left.payload().positions().values().value()),
            coordinate_bits(right.payload().positions().values().value())
        );
    }
}

/// Checks a model with and without perception, so stored perception is
/// exercised wherever Kekule can perceive the input.
fn check_model(mut model: Model) {
    let decoded = round_trip(&model).into_model().expect("a model");
    assert_same_model(&decoded, &model);
    if model.perceive().is_ok() {
        let decoded = round_trip(&model).into_model().expect("a model");
        assert_same_model(&decoded, &model);
    }
}

#[test]
fn curated_small_molecules_round_trip_exactly() {
    let mut checked = 0;
    for path in files(&dataset("small")) {
        let Some(extension) = path.extension().and_then(|extension| extension.to_str()) else {
            continue;
        };
        let Ok(text) = std::fs::read_to_string(&path) else {
            continue;
        };
        match extension {
            // ChEMBL Molfiles carry SDF data fields without a final record
            // delimiter, so both read as SDF records; an input that still does
            // not falls back to the Molfile reader.
            "mol" | "sdf" => {
                let mut record_text = text.clone();
                if !text.trim_end().ends_with("$$$$") {
                    record_text.push_str(
                        "
$$$$
",
                    );
                }
                if let Ok(document) = sdf::parse_str(&record_text, Default::default()) {
                    for record in document.records() {
                        if let Ok(interpretation) = record.interpret(Default::default()) {
                            check_model(interpretation.model().clone());
                            checked += 1;
                        }
                    }
                } else if let Ok(interpretation) = molfile::parse_str(&text, Default::default())
                    .and_then(|document| document.interpret(Default::default()))
                {
                    check_model(interpretation.into_model());
                    checked += 1;
                }
            }
            "smi" => {
                let line = text.lines().next().unwrap_or_default();
                let Ok(molecules) =
                    smiles::to_molecules(line.split_whitespace().next().unwrap_or(""))
                else {
                    continue;
                };
                for mut molecule in molecules {
                    let _ = molecule.perceive();
                    let decoded = round_trip(&molecule).into_molecule().expect("a molecule");
                    assert_same_molecule(&decoded, &molecule);
                }
                checked += 1;
            }
            _ => {}
        }
    }
    assert!(
        checked > 600,
        "only {checked} small-molecule inputs were checked"
    );
}

#[test]
fn curated_biomolecules_round_trip_exactly() {
    let mut checked = 0;
    for path in files(&dataset("bio")) {
        if !path.to_string_lossy().ends_with(".cif.gz") {
            continue;
        }
        let mut text = String::new();
        flate2::read::GzDecoder::new(std::fs::File::open(&path).expect("dataset file"))
            .read_to_string(&mut text)
            .expect("gzip mmCIF");
        let Ok(document) = mmcif::parse_str(&text, Default::default()) else {
            continue;
        };
        if let Ok(interpretation) = document.interpret(Default::default()) {
            check_model(interpretation.into_model());
            checked += 1;
        }
        if let Ok(interpretation) = document.interpret_ensemble(Default::default()) {
            let ensemble = interpretation.into_ensemble();
            if ensemble.len() > 1 {
                let decoded = round_trip(&ensemble).into_ensemble().expect("an ensemble");
                assert_same_ensemble(&decoded, &ensemble);
            }
        }
    }
    assert!(
        checked >= 15,
        "only {checked} macromolecular entries were checked"
    );
}
