//! Syntax-independent molecular queries.
//!
//! [`QueryGraph`] is the canonical query representation. Use
//! [`crate::smarts::parse_str`] to construct one from SMARTS syntax, then pass it to the
//! [`crate::substructure`] matching facade. Parsing a query and matching it
//! against a molecule are separate operations; matching does not perceive the
//! target implicitly.
//! Stereo constraints compare represented local configuration under the atom
//! mapping before match deduplication or limits. They do not compare CIP labels
//! or require CIP assignment. CXSMARTS `a`, `o`, `&` and `r` fields construct
//! enhanced groups. These enable correlated matching automatically; plain queries
//! may opt in through [`crate::substructure::SubstructureMatchOptions::use_enhanced_stereo`].
//! An achiral query places no restriction on the target's stereochemistry.
//!
//! ```
//! use kekule::{smarts, smiles, substructure};
//!
//! let mut target = smiles::to_molecules("CCO")?.pop().unwrap();
//! target.perceive()?;
//! let query = smarts::parse_str("[#6]-[#8]", Default::default())?;
//! let matched = substructure::PreparedTarget::new(&target).first(&query, Default::default())?;
//! assert!(matched.is_some());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

#![doc = include_str!("dialect.md")]

mod expression;
mod graph;
pub(crate) mod smarts;
mod stereo;

pub use expression::{
    AtomExpression, AtomPredicate, BondExpression, BondPredicate, QueryExpression,
    QueryExpressionError, MAX_QUERY_EXPRESSION_DEPTH, MAX_QUERY_EXPRESSION_NODES,
};
pub use graph::{
    QueryAtom, QueryAtomId, QueryBond, QueryBondId, QueryGraph, QueryGraphBuilder, QueryGraphError,
    QueryStereoIssue, MAX_QUERY_RECURSION_DEPTH, MAX_QUERY_TOTAL_ATOMS,
    MAX_QUERY_TOTAL_EXPRESSION_NODES,
};
pub use stereo::{QueryStereoConstraint, QueryStereoGroup};
