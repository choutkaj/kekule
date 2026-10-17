use crate::diagnostics::{LimitExceeded, Resource};
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

/// Resource bounds and match semantics for substructure search.
///
/// [`Default`] gives fixed bounds, keeps one match per target atom set, and
/// matches stereo by local parity only.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SubstructureMatchOptions {
    max_matches: usize,
    max_search_states: usize,
    max_query_atoms: usize,
    max_candidate_pairs: usize,
    uniquify: bool,
    use_enhanced_stereo: bool,
}

impl SubstructureMatchOptions {
    /// Maximum matches. Enumeration with `all`, `visit`, or
    /// [`TaggedQuery::find_matches`] reports exceeding it as
    /// [`SubstructureMatchError::LimitExceeded`]; `first` stops at one match.
    pub const fn max_matches(&self) -> usize {
        self.max_matches
    }

    /// Sets [`Self::max_matches`].
    ///
    /// # Errors
    ///
    /// [`SubstructureMatchError::InvalidLimit`] when `max_matches` is zero.
    pub fn with_max_matches(mut self, max_matches: usize) -> Result<Self, SubstructureMatchError> {
        if max_matches == 0 {
            return Err(SubstructureMatchError::InvalidLimit {
                resource: Resource::Matches,
                limit: max_matches,
            });
        }
        self.max_matches = max_matches;
        Ok(self)
    }

    /// Maximum candidate assignments visited by backtracking.
    pub const fn max_search_states(&self) -> usize {
        self.max_search_states
    }

    /// Sets [`Self::max_search_states`].
    ///
    /// # Errors
    ///
    /// [`SubstructureMatchError::InvalidLimit`] when `max_search_states` is zero.
    pub fn with_max_search_states(
        mut self,
        max_search_states: usize,
    ) -> Result<Self, SubstructureMatchError> {
        if max_search_states == 0 {
            return Err(SubstructureMatchError::InvalidLimit {
                resource: Resource::SearchStates,
                limit: max_search_states,
            });
        }
        self.max_search_states = max_search_states;
        Ok(self)
    }

    /// Maximum query size accepted by the recursive search.
    pub const fn max_query_atoms(&self) -> usize {
        self.max_query_atoms
    }

    /// Sets [`Self::max_query_atoms`].
    ///
    /// # Errors
    ///
    /// [`SubstructureMatchError::InvalidLimit`] when `max_query_atoms` is zero
    /// or exceeds [`MAX_SUBSTRUCTURE_QUERY_ATOMS`].
    pub fn with_max_query_atoms(
        mut self,
        max_query_atoms: usize,
    ) -> Result<Self, SubstructureMatchError> {
        if max_query_atoms == 0 || max_query_atoms > MAX_SUBSTRUCTURE_QUERY_ATOMS {
            return Err(SubstructureMatchError::InvalidLimit {
                resource: Resource::Atoms,
                limit: max_query_atoms,
            });
        }
        self.max_query_atoms = max_query_atoms;
        Ok(self)
    }

    /// Maximum query-atom by target-atom compatibility matrix size.
    pub const fn max_candidate_pairs(&self) -> usize {
        self.max_candidate_pairs
    }

    /// Sets [`Self::max_candidate_pairs`].
    ///
    /// # Errors
    ///
    /// [`SubstructureMatchError::InvalidLimit`] when `max_candidate_pairs` is zero.
    pub fn with_max_candidate_pairs(
        mut self,
        max_candidate_pairs: usize,
    ) -> Result<Self, SubstructureMatchError> {
        if max_candidate_pairs == 0 {
            return Err(SubstructureMatchError::InvalidLimit {
                resource: Resource::CandidatePairs,
                limit: max_candidate_pairs,
            });
        }
        self.max_candidate_pairs = max_candidate_pairs;
        Ok(self)
    }

    /// Collapse query-automorphism duplicates by target atom set.
    pub const fn uniquify(&self) -> bool {
        self.uniquify
    }

    /// Sets [`Self::uniquify`].
    #[must_use]
    pub const fn with_uniquify(mut self, uniquify: bool) -> Self {
        self.uniquify = uniquify;
        self
    }

    /// Honor target enhanced groups as correlated configuration choices.
    /// Query stereo groups always enable this behavior. Plain SMARTS retains
    /// literal local-parity matching by default for backward compatibility.
    pub const fn use_enhanced_stereo(&self) -> bool {
        self.use_enhanced_stereo
    }

    /// Sets [`Self::use_enhanced_stereo`].
    #[must_use]
    pub const fn with_use_enhanced_stereo(mut self, use_enhanced_stereo: bool) -> Self {
        self.use_enhanced_stereo = use_enhanced_stereo;
        self
    }
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

/// Work done by a substructure search, reported with
/// [`SubstructureMatchError::LimitExceeded`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SubstructureMatchWork {
    /// Atoms in the query.
    pub query_atoms: usize,
    /// Atoms in the target.
    pub target_atoms: usize,
    /// Query-atom by target-atom compatibility pairs evaluated, including
    /// those of recursive subqueries.
    pub candidate_pairs: usize,
    /// Candidate assignments visited by backtracking.
    pub search_states: usize,
    /// Matches found.
    pub matches: usize,
}

/// One embedding of a query in a molecule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryMatch {
    atoms: Vec<AtomId>,
}

impl QueryMatch {
    /// The target atom matched by `query_atom`, or `None` when the query has
    /// no such atom.
    pub fn atom(&self, query_atom: QueryAtomId) -> Option<AtomId> {
        self.atoms.get(query_atom.index()).copied()
    }

    /// Target atoms in query-atom order.
    pub fn atoms(&self) -> &[AtomId] {
        &self.atoms
    }
}

/// Target perception a query predicate depends on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum QueryPerception {
    /// Valence perception, for hydrogen counts, connectivity, and valence.
    Valence,
    /// Ring membership, for ring-membership and ring-bond-count predicates.
    RingMembership,
    /// An SSSR-like ring basis, for ring-count and smallest-ring-size
    /// predicates.
    RingBasis,
    /// Aromaticity perception, for aromatic atom and bond predicates.
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

/// Substructure search could not run to completion.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SubstructureMatchError {
    /// A configured bound is zero, or above its ceiling
    /// ([`MAX_SUBSTRUCTURE_QUERY_ATOMS`] for query atoms).
    InvalidLimit {
        /// The bounded resource.
        resource: Resource,
        /// The configured bound.
        limit: usize,
    },
    /// The query needs this perception, which the target has not installed.
    MissingPerception(QueryPerception),
    /// The target's installed perception uses a model the query cannot use,
    /// such as a ring basis other than the SSSR-like one.
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
    /// The work done before a bound was exceeded, or `None` for any other
    /// error.
    pub fn work(&self) -> Option<SubstructureMatchWork> {
        match self {
            Self::LimitExceeded { work, .. } => Some(*work),
            Self::InvalidLimit { .. }
            | Self::MissingPerception(_)
            | Self::IncompatiblePerception(_) => None,
        }
    }
}

impl fmt::Display for SubstructureMatchError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidLimit { resource, limit } => {
                write!(f, "invalid substructure {resource} limit {limit}")
            }
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
#[non_exhaustive]
pub enum MatchCompletion {
    /// Every match was delivered.
    Complete,
    /// The visitor stopped the traversal.
    Stopped,
}
