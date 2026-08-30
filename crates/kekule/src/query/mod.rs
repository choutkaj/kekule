//! Syntax-independent molecular queries and SMARTS parsing.
//!
//! [`QueryGraph`] is the canonical query representation. Use [`parse_smarts`]
//! to construct one from SMARTS syntax, then pass it to the
//! [`crate::substructure`] matching facade. Parsing a query and matching it
//! against a molecule are separate operations; matching does not perceive the
//! target implicitly.
//!
//! ```
//! use kekule::{query::parse_smarts, smiles, substructure};
//!
//! let target = smiles::to_molecules("CCO")?.pop().unwrap();
//! let query = parse_smarts("[#6]-[#8]")?;
//! let matched = substructure::find_substructure_match(&target, &query)?;
//! assert!(matched.is_some());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

mod expression;
mod graph;
mod smarts;

pub use expression::*;
pub use graph::*;
pub use smarts::*;
