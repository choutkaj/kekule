//! JSON-lines adapter exposing emitted SMILES to optional external graph checks.

use std::io::{self, BufRead};

use kekule::smiles;
use serde_json::{json, Value};

fn probe(input: &str) -> Result<Value, Box<dyn std::error::Error>> {
    let mut written = Vec::new();
    for mut molecule in smiles::parse_str(input)?.interpret()?.into_molecules() {
        molecule.perceive()?;
        written.push(smiles::write_isomeric(&molecule)?);
    }
    Ok(json!({"status": "ok", "smiles": written.join(".")}))
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
        assert_eq!(
            probe(input).unwrap(),
            json!({"status": "ok", "smiles": input})
        );
        assert!(probe("C(").is_err());
    }
}
