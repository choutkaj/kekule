use crate::diagnostics::LimitExceeded;
use std::fmt;
mod engine;
mod prepared;
mod stereo;
mod topology;
pub use prepared::{PreparedTarget, TaggedMatchError, TaggedQuery, TupleMode};
pub use topology::{PreparedTopologyTarget, TopologyQueryMatch};

use crate::molecule::{AtomId, Molecule};
use crate::query::{QueryAtomId, QueryGraph};

/// Absolute query-size ceiling for the recursive bounded matcher.
pub const MAX_SUBSTRUCTURE_QUERY_ATOMS: usize = 256;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SubstructureMatchOptions {
    /// Maximum matches. Complete enumeration (`all`) reports
    /// exceeding it as [`SubstructureMatchError::LimitExceeded`]; a visitor
    /// stops successfully once it is reached.
    pub max_matches: usize,
    /// Maximum candidate assignments visited by backtracking.
    pub max_search_states: usize,
    /// Maximum query size accepted by the recursive search.
    pub max_query_atoms: usize,
    /// Maximum query-atom by target-atom compatibility matrix size.
    pub max_candidate_pairs: usize,
    /// Collapse query-automorphism duplicates by target atom set.
    pub uniquify: bool,
    /// Honor target enhanced groups as correlated configuration choices.
    /// Query stereo groups always enable this behavior. Plain SMARTS retains
    /// literal local-parity matching by default for backward compatibility.
    pub use_enhanced_stereo: bool,
}

impl SubstructureMatchOptions {
    /// The [`Default`] bound of [`Self::max_matches`].
    pub const DEFAULT_MAX_MATCHES: usize = 1_000;
    /// The [`Default`] bound of [`Self::max_search_states`].
    pub const DEFAULT_MAX_SEARCH_STATES: usize = 1_000_000;
    /// The [`Default`] bound of [`Self::max_candidate_pairs`].
    pub const DEFAULT_MAX_CANDIDATE_PAIRS: usize = 2_000_000;
}

impl Default for SubstructureMatchOptions {
    fn default() -> Self {
        Self {
            max_matches: Self::DEFAULT_MAX_MATCHES,
            max_search_states: Self::DEFAULT_MAX_SEARCH_STATES,
            max_query_atoms: MAX_SUBSTRUCTURE_QUERY_ATOMS,
            max_candidate_pairs: Self::DEFAULT_MAX_CANDIDATE_PAIRS,
            uniquify: true,
            use_enhanced_stereo: false,
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SubstructureMatchWork {
    pub query_atoms: usize,
    pub target_atoms: usize,
    pub candidate_pairs: usize,
    pub search_states: usize,
    pub matches: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryMatch {
    atoms: Vec<AtomId>,
}

impl QueryMatch {
    pub fn atom(&self, query_atom: QueryAtomId) -> Option<AtomId> {
        self.atoms.get(query_atom.index()).copied()
    }

    /// Target atoms in query-atom order.
    pub fn atoms(&self) -> &[AtomId] {
        &self.atoms
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum QueryPerception {
    Valence,
    RingMembership,
    RingBasis,
    Aromaticity,
}

impl fmt::Display for QueryPerception {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Valence => f.write_str("valence"),
            Self::RingBasis => f.write_str("compatible SSSR ring basis"),
            Self::RingMembership => f.write_str("ring membership"),
            Self::Aromaticity => f.write_str("aromaticity"),
        }
    }
}

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SubstructureMatchError {
    InvalidOptions(&'static str),
    MissingPerception(QueryPerception),
    IncompatiblePerception(QueryPerception),
    /// A configured bound was exceeded; `work` is the work done before it.
    LimitExceeded {
        /// The exceeded bound.
        limit: LimitExceeded,
        /// The work done before the bound was exceeded.
        work: SubstructureMatchWork,
    },
}

impl SubstructureMatchError {
    pub fn work(&self) -> Option<SubstructureMatchWork> {
        match self {
            Self::LimitExceeded { work, .. } => Some(*work),
            Self::InvalidOptions(_)
            | Self::MissingPerception(_)
            | Self::IncompatiblePerception(_) => None,
        }
    }
}

impl fmt::Display for SubstructureMatchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidOptions(message) => write!(f, "invalid substructure options: {message}"),
            Self::IncompatiblePerception(perception) => {
                write!(f, "substructure query requires compatible {perception}")
            }
            Self::MissingPerception(perception) => write!(
                f,
                "substructure query requires current {perception} perception on the target"
            ),
            Self::LimitExceeded { limit, .. } => write!(f, "substructure search {limit}"),
        }
    }
}

impl std::error::Error for SubstructureMatchError {}

/// Whether a streaming traversal exhausted the search or was stopped by its visitor.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchCompletion {
    Complete,
    Stopped,
}

fn validate_options(options: SubstructureMatchOptions) -> Result<(), SubstructureMatchError> {
    for (name, value) in [
        ("max_matches must be greater than zero", options.max_matches),
        (
            "max_search_states must be greater than zero",
            options.max_search_states,
        ),
        (
            "max_query_atoms must be greater than zero",
            options.max_query_atoms,
        ),
        (
            "max_candidate_pairs must be greater than zero",
            options.max_candidate_pairs,
        ),
    ] {
        if value == 0 {
            return Err(SubstructureMatchError::InvalidOptions(name));
        }
    }
    if options.max_query_atoms > MAX_SUBSTRUCTURE_QUERY_ATOMS {
        return Err(SubstructureMatchError::InvalidOptions(
            "max_query_atoms exceeds the stack-safe matcher ceiling",
        ));
    }
    Ok(())
}
