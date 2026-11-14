//! Syntax-independent substructure matching algorithms.
//!
//! Matching consumes [`query::QueryGraph`](crate::query::QueryGraph) and current target perception state;
//! it never invokes parsing, interpretation, or perception implicitly.
//! [`PreparedTarget`] borrows one molecule and
//! answers repeated queries against it;
//! [`PreparedTopologyTarget`] does the
//! same for every molecule instance of a topology.
//!
//! ```
//! use kekule::substructure::{PreparedTarget, SubstructureMatchOptions};
//! use kekule::{smarts, smiles};
//!
//! let mut toluene = smiles::to_molecules("Cc1ccccc1")?.remove(0);
//! toluene.perceive()?;
//! let aromatic_carbon = smarts::parse_str("c", Default::default())?;
//! let target = PreparedTarget::new(&toluene);
//! let matches = target.all(&aromatic_carbon, SubstructureMatchOptions::default())?;
//! assert_eq!(matches.len(), 6);
//! let nitrogen = smarts::parse_str("N", Default::default())?;
//! assert!(target.first(&nitrogen, Default::default())?.is_none());
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

pub use crate::algorithms::{
    MatchCompletion, PreparedTarget, PreparedTopologyTarget, QueryMatch, QueryPerception,
    SubstructureMatchError, SubstructureMatchOptions, SubstructureMatchWork, TaggedMatchError,
    TaggedQuery, TopologyQueryMatch, TupleMode, MAX_SUBSTRUCTURE_QUERY_ATOMS,
};
