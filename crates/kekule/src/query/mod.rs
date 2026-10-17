//! Syntax-independent molecular queries.
//!
//! [`QueryGraph`] is the canonical query representation. Use
//! [`crate::smarts::parse_str`] to construct one from SMARTS syntax, then pass it to the
//! [`crate::substructure`] matching facade. Parsing a query and matching it
//! against a molecule are separate operations; matching does not perceive the
//! target implicitly.
//!
//! Each [`QueryAtom`] and [`QueryBond`] carries a boolean [`QueryExpression`]
//! ([`AtomExpression`], [`BondExpression`]) over primitive [`AtomPredicate`]s
//! and [`BondPredicate`]s, including recursive subqueries. Build a graph
//! programmatically with [`QueryGraphBuilder`], which rejects invalid edits with
//! [`QueryGraphError`]; expression and graph sizes are bounded by the `MAX_*`
//! constants. Local stereo requirements are [`QueryStereoConstraint`]s, and
//! enhanced stereo relationships are [`QueryStereoGroup`]s.
//!
//! Stereo constraints compare represented local configuration under the atom
//! mapping before match deduplication or limits. They do not compare CIP labels
//! or require CIP assignment. CXSMARTS `a`, `o`, `&` and `r` fields construct
//! enhanced groups. These enable correlated matching automatically; plain queries
//! may opt in through [`crate::substructure::SubstructureMatchOptions::use_enhanced_stereo`].
//! An achiral query places no restriction on the target's stereochemistry.
//!
//! The same query can be parsed from SMARTS or built programmatically:
//!
//! ```
//! use kekule::molecule::{AtomId, BondOrder, Element};
//! use kekule::query::{AtomExpression, AtomPredicate, BondExpression, BondPredicate, QueryGraph};
//! use kekule::{smarts, smiles, substructure};
//!
//! let mut target = smiles::to_molecules("CCO")?.pop().unwrap();
//! target.perceive()?;
//! let parsed = smarts::parse_str("[#6]-[#8]", Default::default())?;
//!
//! let mut builder = QueryGraph::builder();
//! let carbon = builder.add_atom(AtomExpression::predicate(AtomPredicate::Element(Element::C)))?;
//! let oxygen = builder.add_atom(AtomExpression::predicate(AtomPredicate::Element(Element::O)))?;
//! let single = BondExpression::predicate(BondPredicate::Order(BondOrder::Single));
//! builder.add_bond(carbon, oxygen, single)?;
//! let built = builder.build()?;
//!
//! let prepared = substructure::PreparedTarget::new(&target);
//! for query in [&parsed, &built] {
//!     let matched = prepared.first(query, Default::default())?.expect("C-O matches");
//!     assert_eq!(matched.atoms(), [AtomId::new(1), AtomId::new(2)]);
//! }
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
