//! PDB-format parsing, interpretation, and writing.
//!
//! A PDB file holds one entry, whose coordinate models share one
//! interpretation scope. Parsing types the unit-cell, coordinate,
//! chain-termination, model, and connection records and keeps every other
//! record verbatim;
//! interpretation runs the pipeline [`mmcif`](crate::mmcif) uses and returns the same
//! [`macromolecular`](crate::macromolecular) reports.
//!
//! The format has no entity or label records. A chain's polymer is the run of
//! residues its `TER` record ends, or, without one, the residues from its first
//! to its last `ATOM` record; water and single residues are never polymer.
//! Residues in such a run are adjacent in sequence. `LINK` and `SSBOND`
//! records are covalent connections. `CONECT` records list connectivity
//! without bond types: within a residue they supplement its template and
//! stand in for distance inference, and a partner listed two or three times
//! has a double or triple bond. A `LINK` or `CONECT` between a metal and a
//! non-metal is coordination: within a residue a dative bond, between
//! residues no bond but a record of the report
//! ([`StructureInterpretationReport::coordination`](crate::macromolecular::StructureInterpretationReport::coordination)). Atoms
//! without an element symbol take the one their name's alignment implies,
//! and the report lists them. A `CRYST1` cell is the structure's periodic cell
//! under the rule mmCIF `_cell` follows; otherwise the report says why not.
//!
//! Interpretation keeps coordinate records in source order as dense atom order,
//! and [`write()`](crate::pdb::write) emits records in dense atom order, so a written
//! entry interprets back to the same atoms in the same order.
//!
//! ```
//! use kekule::pdb;
//!
//! let text = "\
//! HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O
//! HETATM    2  H1  HOH A   1       0.960   0.000   0.000  1.00  0.00           H
//! HETATM    3  H2  HOH A   1      -0.240   0.930   0.000  1.00  0.00           H
//! END
//! ";
//! let document = pdb::parse_str(text, Default::default())?;
//! let interpretation = document.interpret(Default::default())?;
//! assert_eq!(interpretation.structure().atom_count(), 3);
//! assert_eq!(interpretation.topology().instance_count(), 1);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub use crate::io::pdb::entry_points::{parse_reader, parse_str, write, write_to};
pub use crate::io::{
    PdbAtom, PdbAtomKind, PdbChainIdPolicy, PdbConect, PdbCryst1, PdbDocument, PdbLink, PdbModel,
    PdbParseOptions, PdbPartner, PdbRecord, PdbSource, PdbSsbond, PdbTer, PdbWriteErrorKind,
    PdbWriteOptions,
};
