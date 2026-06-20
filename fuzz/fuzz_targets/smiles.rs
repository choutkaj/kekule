#![no_main]

use libfuzzer_sys::fuzz_target;
use molecules::{read_smiles_str, write_smiles, SmilesParseOptions, SmilesWriteOptions};

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(molecule) = read_smiles_str(input, SmilesParseOptions) {
        if let Ok(output) = write_smiles(&molecule, SmilesWriteOptions) {
            let _ = read_smiles_str(&output, SmilesParseOptions);
        }
    }
});
