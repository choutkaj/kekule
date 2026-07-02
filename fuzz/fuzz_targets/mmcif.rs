#![no_main]

use libfuzzer_sys::fuzz_target;
use molecules::mmcif::{read_str, MmcifParseOptions};

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(molecule) = read_str(
        input,
        MmcifParseOptions {
            max_input_bytes: 64 * 1024,
            max_tokens: 16 * 1024,
            max_token_bytes: 16 * 1024,
            max_atom_site_rows: 4 * 1024,
            ..MmcifParseOptions::default()
        },
    ) {
        for atom in molecule.graph().atom_ids() {
            let _ = molecule.graph().atom(atom);
            let _ = molecule.hierarchy().atom_site_for_atom(atom);
        }
    }
});
