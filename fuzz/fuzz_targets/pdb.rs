#![no_main]

use kekule::macromolecular::{ModelSelection, StructureInterpretOptions};
use kekule::pdb::{parse_str, PdbParseOptions};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    let Ok(document) = parse_str(
        input,
        PdbParseOptions::default()
            .with_max_input_bytes(64 * 1024)
            .with_max_atom_records(4 * 1024),
    ) else {
        return;
    };
    for atom in document.atoms() {
        let _ = (atom.name(), atom.element(), atom.formal_charge());
    }
    let _ = document.alternate_locations();

    let options = StructureInterpretOptions::default().with_model_selection(ModelSelection::First);
    if let Ok(interpreted) = document.interpret(options.clone()) {
        for atom in interpreted.model().topology().atom_ids() {
            let _ = interpreted.model().topology().atom(*atom);
            let _ = interpreted.model().position(*atom);
        }
        // Whatever the writer accepts reads back with the same atoms.
        if let Ok(written) = kekule::pdb::write(interpreted.model(), Default::default()) {
            let reread = parse_str(&written, PdbParseOptions::default())
                .and_then(|document| document.interpret(options))
                .expect("written PDB reads back");
            assert_eq!(
                reread.model().topology().atom_count(),
                interpreted.model().topology().atom_count()
            );
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
});
