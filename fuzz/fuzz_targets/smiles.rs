#![no_main]

use libfuzzer_sys::fuzz_target;
use molecules::smiles::{
    read_str_with_options, write_with_options, SmilesParseOptions, SmilesWriteOptions,
};

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(molecule) = read_str_with_options(input, SmilesParseOptions::default()) {
        if let Ok(output) = write_with_options(&molecule, SmilesWriteOptions::default()) {
            let _ = read_str_with_options(&output, SmilesParseOptions::default());
        }
    }
});
