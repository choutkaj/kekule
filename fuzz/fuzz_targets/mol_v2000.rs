#![no_main]

use libfuzzer_sys::fuzz_target;
use molecules::molfile::{interpret, parse_str, write_v2000};

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(document) = parse_str(input) {
        let Ok(molecule) = interpret(&document) else { return };
        if let Ok(output) = write_v2000(&molecule) {
            if let Ok(document) = parse_str(&output) {
                let _ = interpret(&document);
            }
        }
    }
});
