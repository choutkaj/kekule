#![no_main]

use libfuzzer_sys::fuzz_target;
use kekule::smiles::{interpret, parse_str, write};

#[path = "../../crates/kekule/tests/support/smiles_contract.rs"]
mod contract;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(document) = parse_str(input) {
        let Ok(interpreted) = interpret(&document) else {
            return;
        };
        for molecule in interpreted.molecules() {
            if let Ok(output) = write(molecule) {
                let document = parse_str(&output).expect("writer output must parse");
                interpret(&document).expect("writer output must interpret");
            }
            // Bound repeated canonical searches independently of parser limits.
            // Larger records still exercise the parsing/plain-writing path above.
            if molecule.atom_count() <= 32 && molecule.bond_count() <= 48 {
                let mut perceived = molecule.clone();
                if perceived.perceive().is_ok() {
                    contract::check_molecule(&perceived);
                }
            }
        }
    }
});
