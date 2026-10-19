#![no_main]

use kekule::diagnostics::UnsupportedContent;
use kekule::xyz::{
    parse_str, XyzEnsembleInterpretOptions, XyzFrameSelection, XyzInterpretOptions,
    XyzParseOptions,
};
use libfuzzer_sys::fuzz_target;

fuzz_target!(|data: &[u8]| {
    let Ok(input) = std::str::from_utf8(data) else {
        return;
    };
    let Ok(document) = parse_str(
        input,
        XyzParseOptions::default()
            .with_max_input_bytes(64 * 1024)
            .with_max_atom_records(4 * 1024),
    ) else {
        return;
    };
    let options = XyzInterpretOptions::default()
        .with_frame(XyzFrameSelection::First)
        .with_extra_columns(UnsupportedContent::Report);
    if let Ok(interpreted) = document.interpret(options) {
        let topology = interpreted.topology();
        assert_eq!(topology.atom_count(), document.frames()[0].atoms().len());
        assert_eq!(
            topology.bond_count(),
            interpreted.report().inferred_bonds().len()
        );
    }
    let options =
        XyzEnsembleInterpretOptions::default().with_extra_columns(UnsupportedContent::Report);
    if let Ok(interpreted) = document.interpret_ensemble(options) {
        assert_eq!(interpreted.ensemble().len(), document.frames().len());
        assert_eq!(interpreted.reports().len(), document.frames().len());
    }
});
