//! Syntax-independent molecular queries and SMARTS parsing.
//!
//! [`QueryGraph`] is the canonical query representation. Use [`parse_smarts`]
//! to construct one from SMARTS syntax, then pass it to the
//! [`crate::substructure`] matching facade. Parsing a query and matching it
//! against a molecule are separate operations; matching does not perceive the
//! target implicitly.
//! Stereo constraints compare represented local configuration under the atom
//! mapping before match deduplication or limits. They do not compare CIP labels
//! or interpret enhanced stereo group relationships. An achiral query places no
//! restriction on the target's stereochemistry.
//!
//! ```
//! use kekule::{query::parse_smarts, smiles, substructure};
//!
//! let mut target = smiles::to_molecules("CCO")?.pop().unwrap();
//! target.perceive()?;
//! let query = parse_smarts("[#6]-[#8]")?;
//! let matched = substructure::find_substructure_match(&target, &query)?;
//! assert!(matched.is_some());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod expression;
mod graph;
mod smarts;
mod stereo;

pub use expression::*;
pub use graph::*;
pub use smarts::*;
pub use stereo::*;
