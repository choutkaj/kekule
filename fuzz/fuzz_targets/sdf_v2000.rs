#![no_main]

use libfuzzer_sys::fuzz_target;
use molecules::sdf::{interpret, parse_str, write_v2000, SdfParseOptions};

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(document) = parse_str(
        input,
        SdfParseOptions {
            allow_missing_final_delimiter: true,
        },
    ) {
        let Ok(interpreted) = interpret(&document) else {
            return;
        };
        if let Ok(output) = write_v2000(interpreted.records()) {
            if let Ok(document) = parse_str(&output, SdfParseOptions::default()) {
                let _ = interpret(&document);
            }
        }
    }
});
