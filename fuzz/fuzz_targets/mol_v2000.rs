#![no_main]

use kekule::molfile::{interpret, parse_str, write, MolfileWriteOptions, MolfileWriteVersion};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(document) = parse_str(input) {
        let Ok(interpreted) = interpret(&document) else {
            return;
        };
        for molecule in interpreted.molecules() {
            let options = MolfileWriteOptions {
                version: MolfileWriteVersion::V2000,
            };
            if let Ok(output) = write(molecule, options) {
                if let Ok(document) = parse_str(&output) {
                    let _ = interpret(&document);
                }
            }
        }
    }
});
