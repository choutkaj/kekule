#![no_main]

use libfuzzer_sys::fuzz_target;
use molecules::smiles::{interpret, parse_str, write_with_options, SmilesWriteOptions};

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(document) = parse_str(input) {
        let Ok(molecule) = interpret(&document) else { return };
        if let Ok(output) = write_with_options(&molecule, SmilesWriteOptions::default()) {
            if let Ok(document) = parse_str(&output) {
                let _ = interpret(&document);
            }
        }
    }
});
