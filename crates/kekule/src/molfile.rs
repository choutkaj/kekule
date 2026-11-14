//! Molfile parsing, interpretation, and V2000/V3000 writing.
//!
//! A parsed [`molfile::MolfileDocument`](crate::molfile::MolfileDocument) is format state. Its interpretation can
//! project to connected molecules or to a geometry-bearing
//! [`crate::structure::Structure`]. Writers reject chemistry that the selected
//! Molfile version cannot represent.
//!
//! ```
//! use kekule::{molfile, smiles};
//!
//! let ethanol = smiles::to_molecules("CCO")?.remove(0);
//! let (text, _) = molfile::write(&ethanol, Default::default())?;
//! let document = molfile::parse_str(&text, Default::default())?;
//! let interpretation = document.interpret(Default::default())?;
//! assert_eq!(interpretation.structure().atom_count(), 3);
//! assert_eq!(interpretation.into_molecules().len(), 1);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub use crate::io::molfile::entry_points::{parse_reader, parse_str, write, write_to};
pub use crate::io::{
    MolfileAtom, MolfileAtomMapping, MolfileAtomParity, MolfileBond, MolfileBondMapping,
    MolfileBondOrder, MolfileBondStereo, MolfileDocument, MolfileHeader, MolfileInterpretIssue,
    MolfileInterpretOptions, MolfileInterpretation, MolfileInterpretationReport,
    MolfileParseOptions, MolfileRadical, MolfileSource, MolfileStereoGroup, MolfileStereoGroupKind,
    MolfileUnsupported, MolfileUnsupportedContent, MolfileVersion, MolfileWriteOptions,
    MolfileWriteVersion,
};
