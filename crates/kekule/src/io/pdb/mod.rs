//! PDB: record parsing and writing. Interpretation is shared with mmCIF in
//! `io::macromolecular`.

pub(super) mod document;
pub(crate) mod entry_points;
mod write;

pub(crate) use document::{parse_pdb_reader, parse_pdb_str};
pub use document::{
    PdbAtom, PdbAtomKind, PdbConect, PdbCryst1, PdbDocument, PdbLink, PdbModel, PdbParseOptions,
    PdbPartner, PdbRecord, PdbSsbond, PdbTer,
};
pub(crate) use write::{write_pdb, write_pdb_to};
pub use write::{
    PdbAtomField, PdbCellReason, PdbChainIdPolicy, PdbSource, PdbWriteErrorKind, PdbWriteOptions,
};
