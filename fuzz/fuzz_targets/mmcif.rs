#![no_main]

use libfuzzer_sys::fuzz_target;
use molecules::mmcif::{
    interpret, parse_str, MmcifInterpretOptions, MmcifModelSelection, MmcifParseOptions,
};

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(document) = parse_str(
        input,
        MmcifParseOptions {
            max_input_bytes: 64 * 1024,
            max_tokens: 16 * 1024,
            max_token_bytes: 16 * 1024,
            max_atom_site_rows: 4 * 1024,
            ..MmcifParseOptions::default()
        },
    ) {
        let options = MmcifInterpretOptions {
            model_selection: MmcifModelSelection::First,
            ..MmcifInterpretOptions::default()
        };
        if let Ok(interpreted) = interpret(&document, options) {
            for atom in interpreted.model().topology().atom_ids() {
                let _ = interpreted.model().topology().atom(*atom);
                let _ = interpreted.model().position(*atom);
            }
        }
    }
});
