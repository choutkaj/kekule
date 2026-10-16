#![no_main]

use kekule::macromolecular::{
    ModelSelection, StructureEnsembleInterpretOptions, StructureInterpretOptions,
};
use kekule::pdb::{
    interpret_ensemble, interpret_ensemble_with_options, interpret_with_options,
    parse_str_with_options, PdbParseOptions,
};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    let Ok(document) = parse_str_with_options(
        input,
        PdbParseOptions {
            max_input_bytes: 64 * 1024,
            max_atom_records: 4 * 1024,
        },
    ) else {
        return;
    };
    for atom in document.atoms() {
        let _ = (atom.name(), atom.element(), atom.formal_charge());
    }
    let _ = document.alternate_locations();

    let options = StructureInterpretOptions {
        model_selection: ModelSelection::First,
        ..StructureInterpretOptions::default()
    };
    if let Ok(interpreted) = interpret_with_options(&document, options) {
        for atom in interpreted.model().topology().atom_ids() {
            let _ = interpreted.model().topology().atom(*atom);
            let _ = interpreted.model().position(*atom);
        }
    }
    if let Ok(interpreted) = interpret_ensemble(&document) {
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
    let _ = interpret_ensemble_with_options(
        &document,
        StructureEnsembleInterpretOptions {
            model_ids: Some(Vec::new()),
            ..StructureEnsembleInterpretOptions::default()
        },
    );
});
