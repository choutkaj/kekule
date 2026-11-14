//! The atom selection language: text such as `protein and name CA` or
//! `within 5 of resname LIG and not water`, evaluated to an
//! [`AtomSelection`](crate::topology::AtomSelection).
//!
//! The keyword grammar follows VMD, PyMOL and MDAnalysis where they agree:
//! `and`, `or`, `not`, parentheses, value lists (`resname ALA GLY`), ranges
//! (`resid 10-50`), comparisons (`bfactor > 50`), quoted values, `within`
//! distances in ångströms unless a unit follows, expansion (`byres`,
//! `bychain`, `bymolecule`, `bonded to`), and macros such as `protein`,
//! `water`, `ligand`, `backbone`, and `sidechain`. Chain, residue, and atom
//! identifiers are author identifiers, falling back to label identifiers;
//! `label_*` terms match label identifiers only. The full grammar is in
//! `docs/selection-language.md`.
//!
//! [`select::parse`](crate::select::parse) returns a [`select::SelectionExpr`](crate::select::SelectionExpr) whose `Display` is
//! its canonical text. [`select::SelectionExpr::evaluate`](crate::select::SelectionExpr::evaluate) evaluates
//! topology terms against a topology;
//! [`select::SelectionExpr::evaluate_structure`](crate::select::SelectionExpr::evaluate_structure) also evaluates spatial and
//! site terms against a structure, geometry, or trajectory frame. Chemistry
//! terms read installed perception and never perceive. Applications add
//! `#...`, `$name`, and keyword terms through a [`select::SelectionResolver`](crate::select::SelectionResolver),
//! which may answer spatial references with another object's points and add
//! such points to any `within` reference, as in `within 5 of ligand`.
//!
//! ```
//! use kekule::select::{self, SelectEvaluateOptions, SelectParseOptions};
//! use kekule::smiles;
//!
//! let topology = smiles::to_topology("c1ccccc1O.[OH2]")?.perceived()?;
//! let expression = select::parse("aromatic or (element O and not water)", SelectParseOptions::default())?;
//! assert_eq!(expression.to_string(), "aromatic or element O and not water");
//! let atoms = expression.evaluate(topology, SelectEvaluateOptions::default())?;
//! assert_eq!(atoms.len(), 7);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub use crate::selection_language::{
    parse, HostRequest, HostRole, HostTerm, HostTermError, HostTermKind, HostValue,
    ReferencePoints, SelectEvaluateError, SelectEvaluateOptions, SelectParseOptions, SelectionExpr,
    SelectionResolver, DEFAULT_SELECT_MAX_DEPTH, DEFAULT_SELECT_MAX_INPUT_BYTES,
    DEFAULT_SELECT_MAX_SMARTS_MATCHES,
};
