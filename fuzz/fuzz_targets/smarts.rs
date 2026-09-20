#![no_main]
use kekule::query::{parse_smarts_with_options, SmartsParseOptions};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    let result = parse_smarts_with_options(
        input,
        SmartsParseOptions {
            max_input_bytes: 4096,
            max_atoms: 64,
            max_bonds: 128,
            max_recursive_depth: 8,
            max_total_expression_nodes: 1024,
            ..Default::default()
        },
    );
    if let Err(error) = result {
        let span = error.span();
        assert!(span.start <= span.end && span.end <= input.len());
    }
});
