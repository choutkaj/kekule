//! Molfile: V2000 and V3000 connection-table parsing, interpretation of
//! Molfile and SDF records, and writing.

pub(super) mod document;
pub(crate) mod entry_points;
pub(super) mod interpret;
mod v2000;
mod v3000;
pub(super) mod write;

pub(crate) use document::{parse_molfile_document, parse_molfile_reader};
pub use document::{
    MolfileAtom, MolfileAtomParity, MolfileBond, MolfileBondOrder, MolfileBondStereo,
    MolfileDocument, MolfileHeader, MolfileParseOptions, MolfileRadical, MolfileStereoGroup,
    MolfileStereoGroupKind, MolfileUnsupported, MolfileUnsupportedContent, MolfileVersion,
};
pub use interpret::{
    MolfileAtomMapping, MolfileBondMapping, MolfileInterpretOptions, MolfileInterpretation,
    MolfileInterpretationReport, MolfileInterpretationWarning,
};
pub(crate) use write::{write_molfile, write_molfile_to};
pub use write::{MolfileSource, MolfileWriteOptions, MolfileWriteVersion};
