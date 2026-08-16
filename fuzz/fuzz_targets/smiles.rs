#![no_main]

use libfuzzer_sys::fuzz_target;
use kekule::smiles::{interpret, parse_str, write};

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
                if let Ok(document) = parse_str(&output) {
                    let _ = interpret(&document);
                }
            }
        }
    }
});
