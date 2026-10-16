#![no_main]
use kekule::{smarts, smarts::SmartsParseOptions};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    let result = smarts::parse_str(
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
        let location = error.location().expect("SMARTS errors are located");
        let column = location
            .column()
            .expect("SMARTS errors carry a byte column");
        assert!(location.line() == 1 && (1..=input.len() + 1).contains(&column));
    }
});
