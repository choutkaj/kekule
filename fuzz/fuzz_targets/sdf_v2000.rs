#![no_main]

use kekule::molfile::MolfileWriteVersion;
use kekule::sdf::{parse_str, records, write, SdfParseOptions, SdfWriteOptions};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    // Streaming accepts arbitrary bytes: non-UTF-8 records fail alone.
    let options = SdfParseOptions::default().with_allow_missing_final_delimiter(true);
    let interpreted = records(data, options)
        .filter_map(Result::ok)
        .filter_map(|record| record.interpret(Default::default()).ok())
        .collect::<Vec<_>>();
    let options = SdfWriteOptions::default().with_version(MolfileWriteVersion::V2000);
    if let Ok(output) = write(&interpreted, options) {
        if let Ok(document) = parse_str(&output, Default::default()) {
            for record in document.records() {
                let _ = record.interpret(Default::default());
            }
        }
    }
});
