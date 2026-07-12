#![no_main]

use libfuzzer_sys::fuzz_target;
use molecules::mmcif::{
    interpret, parse_str, MmcifEntry, MmcifInterpretOptions, MmcifParseOptions,
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

        if let Ok(interpreted) = interpret(&document, MmcifInterpretOptions::default()) {
            for molecule in interpreted
                .contents()
                .small_molecules()
                .chain(interpreted.contents().solvent().molecules())
            {
                for atom in molecule.graph().atom_ids() {
                    let _ = molecule.graph().atom(atom);
                }
            }
            for molecule in interpreted.contents().macromolecules() {
                for atom in molecule.graph().atom_ids() {
                    let _ = molecule.graph().atom(atom);
                    let _ = molecule.hierarchy().atom_site_for_atom(atom);
                }
            }
        }
    }
});
