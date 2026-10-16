#![no_main]

use kekule::molfile::MolfileWriteVersion;
use kekule::sdf::{parse_str, records_with_options, write, SdfParseOptions, SdfWriteOptions};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // Streaming accepts arbitrary bytes: non-UTF-8 records fail alone.
    let options = SdfParseOptions {
        allow_missing_final_delimiter: true,
        ..SdfParseOptions::default()
    };
    let interpreted = records_with_options(data, options)
        .filter_map(Result::ok)
        .filter_map(|record| record.interpret().ok())
        .collect::<Vec<_>>();
    let options = SdfWriteOptions {
        version: MolfileWriteVersion::V2000,
        ..Default::default()
    };
    if let Ok(output) = write(&interpreted, options) {
        if let Ok(document) = parse_str(&output) {
            for record in document.records() {
                let _ = record.interpret();
            }
        }
    }
});
