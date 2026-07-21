#![no_main]

use libfuzzer_sys::fuzz_target;
use molecular::mmcif::{
    interpret, parse_str, MmcifEntry, MmcifInterpretOptions, MmcifModelSelection, MmcifParseOptions,
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
        for block in document.blocks() {
            let _ = document.block(block.name());
            for entry in block.entries() {
                match entry {
                    MmcifEntry::Item(item) => {
                        let _ = block.item(item.tag());
                        let _ = item.value().optional_text();
                    }
                    MmcifEntry::Loop(table) => {
                        for row in 0..table.row_count() {
                            let _ = table.row(row);
                            for tag in table.tags() {
                                let _ = table.value(row, tag);
                            }
                        }
                    }
                }
            }
        }

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
