#![no_main]

use kekule::molfile::{parse_str, write, MolfileWriteOptions, MolfileWriteVersion};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(document) = parse_str(input, Default::default()) {
        let Ok(interpreted) = document.interpret(Default::default()) else {
            return;
        };
        for molecule in interpreted.molecules() {
            let options = MolfileWriteOptions {
                version: MolfileWriteVersion::V2000,
                ..Default::default()
            };
            if let Ok(output) = write(molecule, options) {
                if let Ok(document) = parse_str(&output, Default::default()) {
                    let _ = document.interpret(Default::default());
                }
            }
        }
    }
});
