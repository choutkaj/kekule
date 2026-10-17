#![no_main]

use kekule::macromolecular::{ModelSelection, StructureInterpretOptions};
use kekule::mmcif::{parse_str, MmcifEntry, MmcifParseOptions};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    if let Ok(document) = parse_str(
        input,
        MmcifParseOptions::default()
            .with_max_input_bytes(64 * 1024)
            .with_max_tokens(16 * 1024)
            .with_max_token_bytes(16 * 1024)
            .with_max_atom_site_rows(4 * 1024),
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
                    _ => {}
                }
            }
        }

        let options =
            StructureInterpretOptions::default().with_model_selection(ModelSelection::First);
        if let Ok(interpreted) = document.interpret(options) {
            for atom in interpreted.model().topology().atom_ids() {
                let _ = interpreted.model().topology().atom(*atom);
                let _ = interpreted.model().position(*atom);
            }
        }
        if let Ok(interpreted) = document.interpret_ensemble(Default::default()) {
            for view in interpreted
                .ensemble()
                .iter()
                .map(|member| member.as_model_view())
            {
                for atom in view.topology().atom_ids() {
                    let _ = view.topology().atom(*atom);
                    let _ = view.position(*atom);
                }
            }
        }
    }
});
