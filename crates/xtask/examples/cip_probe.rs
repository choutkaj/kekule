//! JSON-lines adapter for the optional external CIP comparison runner.
//!
//! Build with `cargo build -p xtask --example cip_probe --locked`.

use std::io::{self, BufRead};

use kekule::{core::StereoElementKind, smiles, stereo};
use serde_json::{json, Value};

fn failure(stage: &str, error: impl std::fmt::Debug) -> Value {
    json!({"status": "error", "stage": stage, "message": format!("{error:?}")})
}

fn probe(input: &str) -> Result<Value, Value> {
    let document = smiles::parse_str(input).map_err(|error| failure("parse", error))?;
    let interpretation =
        smiles::interpret(&document).map_err(|error| failure("interpret", error))?;
    let mut atoms = Vec::new();
    let mut bonds = Vec::new();
    let mut atom_count = 0;
    let mut bond_count = 0;
    for mut molecule in interpretation.into_molecules() {
        molecule
            .perceive()
            .map_err(|error| failure("perception", error))?;
        stereo::assign_cip_descriptors(&mut molecule).map_err(|error| failure("cip", error))?;
        for (id, element) in molecule.stereo_elements() {
            let Some(descriptor) = molecule
                .cip_descriptor(id)
                .map_err(|error| failure("descriptor", error))?
            else {
                continue;
            };
            let descriptor = format!("{descriptor:?}");
            match &element.kind {
                StereoElementKind::Tetrahedral(stereo) => {
                    atoms.push((atom_count + stereo.center.index(), descriptor));
                }
                StereoElementKind::DoubleBond(stereo) => bonds.push((
                    atom_count + stereo.left.index().min(stereo.right.index()),
                    atom_count + stereo.left.index().max(stereo.right.index()),
                    descriptor,
                )),
                StereoElementKind::Axis(stereo) => {
                    let (left, right) = molecule
                        .bond(stereo.axis)
                        .map_err(|error| failure("axis", error))?
                        .endpoints();
                    bonds.push((
                        atom_count + left.index().min(right.index()),
                        atom_count + left.index().max(right.index()),
                        descriptor,
                    ));
                }
            }
        }
        atom_count += molecule.atom_count();
        bond_count += molecule.bond_count();
    }
    atoms.sort();
    bonds.sort();
    Ok(json!({
        "status": "ok",
        "atom_count": atom_count,
        "bond_count": bond_count,
        "labels": {"atoms": atoms, "bonds": bonds},
    }))
}

fn main() {
    for line in io::stdin().lock().lines() {
        let value = match line {
            Ok(line) => probe(&line).unwrap_or_else(|error| error),
            Err(error) => {
                println!("{}", failure("stdin", error));
                break;
            }
        };
        println!("{value}");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn disconnected_input_preserves_indices_and_explicit_isotopic_hydrogen() {
        let value = probe("[Na+].[2H][C@](F)(Cl)Br").unwrap();
        assert_eq!(value["atom_count"], 6);
        assert_eq!(value["bond_count"], 4);
        assert_eq!(value["labels"]["atoms"], json!([[2, "S"]]));
        assert_eq!(value["labels"]["bonds"], json!([]));
    }

    #[test]
    fn invalid_input_is_an_error_record_instead_of_an_empty_assignment() {
        let value = probe("C(").unwrap_err();
        assert_eq!(value["status"], "error");
        assert_eq!(value["stage"], "parse");
        assert!(value["message"]
            .as_str()
            .unwrap()
            .contains("SmilesParseError"));
    }

    #[test]
    fn bond_labels_use_atom_endpoints_across_components_and_ring_closures() {
        let value = probe("C1CC2CCC1C2.F/C=C/Cl").unwrap();
        assert_eq!(value["atom_count"], 11);
        assert_eq!(value["bond_count"], 11);
        assert_eq!(value["labels"]["atoms"], json!([]));
        assert_eq!(value["labels"]["bonds"], json!([[8, 9, "E"]]));
    }
}
