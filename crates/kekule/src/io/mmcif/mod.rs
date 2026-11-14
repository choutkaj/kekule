//! mmCIF: document parsing, category access, and writing. Interpretation is
//! shared with PDB in `io::macromolecular`.

pub(super) mod category;
mod document;
pub(crate) mod entry_points;
pub(super) mod write;

pub(crate) use document::{parse_mmcif_reader, parse_mmcif_str};
pub use document::{
    MmcifBlock, MmcifDocument, MmcifEntry, MmcifItem, MmcifLoopTable, MmcifParseOptions, MmcifValue,
};
pub(crate) use write::{write_mmcif, write_mmcif_to};
pub use write::{
    MmcifBlockSource, MmcifEntityClassifications, MmcifSharedCategory, MmcifWriteErrorKind,
    MmcifWriteOptions,
};
