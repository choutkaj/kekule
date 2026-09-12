//! JSON-lines adapter exposing emitted SMILES to optional external graph checks.

use std::io::{self, BufRead};

use kekule::smiles;
use serde_json::{json, Value};

fn probe(input: &str) -> Result<Value, Box<dyn std::error::Error>> {
    let mut molecules = smiles::parse_str(input)?.interpret()?.into_molecules();
    for molecule in &mut molecules {
        molecule.perceive()?;
    }
    let write = |mode| -> Result<String, String> {
        let mut parts = molecules
            .iter()
            .map(|mol| smiles::write_molecule(mol, smiles::SmilesWriteOptions { mode }))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|error| error.to_string())?;
        if mode == smiles::SmilesWriteMode::Canonical {
            parts.sort();
        }
        Ok(parts.join("."))
    };
    Ok(
        json!({"status": "ok", "isomeric": write(smiles::SmilesWriteMode::Isomeric), "canonical": write(smiles::SmilesWriteMode::Canonical)}),
    )
}

fn main() {
    for line in io::stdin().lock().lines() {
        let value = match line {
            Ok(line) => probe(&line).unwrap_or_else(
                |error| json!({"status": "error", "message": format!("{error:?}")}),
            ),
            Err(error) => {
                println!(
                    "{}",
                    json!({"status": "error", "message": error.to_string()})
                );
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
    fn emits_all_components_without_fixing_metal_neighbor_hydrogen_declarations() {
        let input = "O[Fe]=O.O[Fe]=O.[Fe]";
        assert_eq!(probe(input).unwrap()["isomeric"], json!({"Ok": input}));
        assert!(probe("C(").is_err());
    }
}
