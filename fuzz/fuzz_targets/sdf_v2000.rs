#![no_main]

use kekule::sdf::{
    interpret, parse_str, parse_str_with_options, write, MolfileWriteVersion, SdfParseOptions,
    SdfWriteOptions,
};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(document) = parse_str_with_options(
        input,
        SdfParseOptions {
            allow_missing_final_delimiter: true,
            ..SdfParseOptions::default()
        },
    ) {
        let Ok(interpreted) = interpret(&document) else {
            return;
        };
        let options = SdfWriteOptions {
            version: MolfileWriteVersion::V2000,
        };
        if let Ok(output) = write(interpreted.records(), options) {
            if let Ok(document) = parse_str(&output) {
                let _ = interpret(&document);
            }
        }
    }
});
