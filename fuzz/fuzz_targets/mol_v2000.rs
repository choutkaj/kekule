#![no_main]

use libfuzzer_sys::fuzz_target;
use molecules::{read_mol_v2000_str, write_mol_v2000};

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(molecule) = read_mol_v2000_str(input) {
        if let Ok(output) = write_mol_v2000(&molecule) {
            let _ = read_mol_v2000_str(&output);
        }
    }
});
