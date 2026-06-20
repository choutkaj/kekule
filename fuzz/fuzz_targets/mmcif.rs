#![no_main]

use libfuzzer_sys::fuzz_target;
use molecules::{read_mmcif_str, MmcifParseOptions};

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(molecule) = read_mmcif_str(
        input,
        MmcifParseOptions {
            max_input_bytes: 64 * 1024,
            max_tokens: 16 * 1024,
            max_token_bytes: 16 * 1024,
            max_atom_site_rows: 4 * 1024,
            ..MmcifParseOptions::default()
        },
    ) {
        for atom in molecule.mol.atom_ids() {
            let _ = molecule.mol.atom(atom);
            let _ = molecule.hierarchy.atom_site_for_atom(atom);
        }
    }
});
