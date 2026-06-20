#![no_main]

use libfuzzer_sys::fuzz_target;
use molecules::{read_sdf_v2000_records, write_sdf_v2000, SdfParseOptions};

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(records) = read_sdf_v2000_records(
        input,
        SdfParseOptions {
            allow_missing_final_delimiter: true,
        },
    ) {
        let molecules = records
            .into_iter()
            .map(|record| record.molecule)
            .collect::<Vec<_>>();
        if let Ok(output) = write_sdf_v2000(&molecules) {
            let _ = read_sdf_v2000_records(&output, SdfParseOptions::default());
        }
    }
});
