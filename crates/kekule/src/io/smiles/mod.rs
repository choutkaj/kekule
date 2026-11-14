//! SMILES and CXSMILES: record parsing, interpretation into connected
//! molecules, and ordinary, isomeric, and canonical writing.

mod canonical;
mod cx;
mod emit;
mod entry_points;
mod interpret;
mod parse;
mod write;

pub use canonical::{
    MAX_CANONICAL_GRAPH_SLOTS, MAX_CANONICAL_INPUT_COMPLEXITY, MAX_CANONICAL_PENDING_LABELS,
    MAX_CANONICAL_REFINEMENT_WORK, MAX_CANONICAL_SEARCH_STATES,
};
pub use cx::{SmilesCxField, SmilesCxFieldKind, SmilesCxStereoGroupKind};
pub use entry_points::{
    parse_reader, parse_str, records, to_molecules, to_topology, write, write_to, SmilesSource,
    SmilesWriteMode, SmilesWriteOptions,
};
pub use interpret::{
    SmilesAtomMapping, SmilesBondMapping, SmilesComponentInterpretation, SmilesInterpretOptions,
    SmilesInterpretation, SmilesInterpretationReport,
};
pub use parse::{
    SmilesDocument, SmilesDocumentToken, SmilesDocumentTokenKind, SmilesParseOptions, SmilesRecords,
};
pub use write::{MAX_DIRECTIONAL_COMBINATIONS, MAX_DIRECTIONAL_VISITS};
