//! SMILES parsing, interpretation, and molecule/topology writing.
//!
//! [`smiles::parse_str`](crate::smiles::parse_str) preserves source syntax in a
//! [`smiles::SmilesDocument`](crate::smiles::SmilesDocument), while
//! [`SmilesDocument::interpret`](crate::smiles::SmilesDocument::interpret) publishes
//! canonical connected molecules. [`smiles::to_molecules`](crate::smiles::to_molecules) is the concise
//! parse-and-interpret path. Dot-separated components remain separate molecules
//! and perception is never run implicitly. [`smiles::write`](crate::smiles::write) writes a molecule
//! or a whole topology as ordinary, isomeric, or canonical SMILES.
//!
//! ```
//! use kekule::smiles::{self, SmilesWriteOptions};
//!
//! let molecules = smiles::to_molecules("CCO.[Na+].[Cl-]")?;
//! assert_eq!(molecules.len(), 3);
//! assert_eq!(molecules[0].atom_count(), 3);
//! assert_eq!(smiles::write(&molecules[0], SmilesWriteOptions::ordinary())?, "CCO");
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub use crate::io::smiles::{
    parse_reader, parse_str, records, to_molecules, to_topology, write, write_to, SmilesSource,
    SmilesWriteMode, SmilesWriteOptions,
};
pub use crate::io::{
    SmilesAtomMapping, SmilesBondMapping, SmilesComponentInterpretation, SmilesCxField,
    SmilesCxFieldKind, SmilesCxStereoGroupKind, SmilesDocument, SmilesDocumentToken,
    SmilesDocumentTokenKind, SmilesInterpretOptions, SmilesInterpretation,
    SmilesInterpretationReport, SmilesParseOptions, SmilesRecords, MAX_CANONICAL_GRAPH_SLOTS,
    MAX_CANONICAL_INPUT_COMPLEXITY, MAX_CANONICAL_PENDING_LABELS, MAX_CANONICAL_REFINEMENT_WORK,
    MAX_CANONICAL_SEARCH_STATES, MAX_DIRECTIONAL_COMBINATIONS, MAX_DIRECTIONAL_VISITS,
};
