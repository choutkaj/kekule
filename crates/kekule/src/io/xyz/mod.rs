//! XYZ: frame parsing and interpretation with inferred bonds.

mod document;
pub(crate) mod entry_points;
mod interpret;
pub(crate) mod write;

pub(crate) use document::{parse_xyz_reader, parse_xyz_str};
pub use document::{XyzAtom, XyzDocument, XyzFrame, XyzFrames, XyzParseOptions};
pub use interpret::{
    XyzFrameSelection, XyzInferredBond, XyzInterpretIssue, XyzInterpretOptions, XyzInterpretation,
    XyzInterpretationReport, XyzStructureSetInterpretOptions, XyzStructureSetInterpretation,
    DEFAULT_XYZ_MAX_ATOM_CONTACTS,
};
pub use write::{XyzSource, XyzWriteOptions};
