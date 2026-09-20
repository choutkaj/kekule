use std::fmt;
mod engine;
mod prepared;
mod stereo;
mod topology;
pub use prepared::*;
pub use topology::*;

use crate::core::{AtomId, Molecule};
use crate::query::{QueryAtomId, QueryGraph};

/// Absolute query-size ceiling for the recursive bounded matcher.
pub const MAX_SUBSTRUCTURE_QUERY_ATOMS: usize = 256;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SubstructureMatchOptions {
    /// Maximum returned matches; reaching this cap stops successfully.
    pub max_matches: usize,
    /// Maximum candidate assignments visited by backtracking.
    pub max_search_states: usize,
    /// Maximum query size accepted by the recursive search.
    pub max_query_atoms: usize,
    /// Maximum query-atom by target-atom compatibility matrix size.
    pub max_candidate_pairs: usize,
    /// Collapse query-automorphism duplicates by target atom set.
    pub uniquify: bool,
}

impl Default for SubstructureMatchOptions {
    fn default() -> Self {
        Self {
            max_matches: 1_000,
            max_search_states: 1_000_000,
            max_query_atoms: MAX_SUBSTRUCTURE_QUERY_ATOMS,
            max_candidate_pairs: 2_000_000,
            uniquify: true,
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
    ResourceLimit {
        resource: &'static str,
        observed: usize,
        limit: usize,
        work: SubstructureMatchWork,
    },
}

impl SubstructureMatchError {
    pub fn work(&self) -> Option<SubstructureMatchWork> {
        match self {
            Self::ResourceLimit { work, .. } => Some(*work),
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
            Self::ResourceLimit {
                resource,
                observed,
                limit,
                ..
            } => write!(
                f,
                "substructure {resource} limit exceeded: observed {observed}, limit {limit}"
            ),
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

pub fn find_substructure_match(
    target: &Molecule,
    query: &QueryGraph,
) -> Result<Option<QueryMatch>, SubstructureMatchError> {
    let options = SubstructureMatchOptions {
        max_matches: 1,
        uniquify: false,
        ..Default::default()
    };
    Ok(
        find_substructure_matches_with_options(target, query, options)?
            .into_iter()
            .next(),
    )
}
pub fn find_substructure_matches(
    target: &Molecule,
    query: &QueryGraph,
) -> Result<Vec<QueryMatch>, SubstructureMatchError> {
    find_substructure_matches_with_options(target, query, SubstructureMatchOptions::default())
}
/// Legacy bounded collection: reaching `max_matches` stops successfully.
pub fn find_substructure_matches_with_options(
    target: &Molecule,
    query: &QueryGraph,
    options: SubstructureMatchOptions,
) -> Result<Vec<QueryMatch>, SubstructureMatchError> {
    let prepared = PreparedTarget::new(target);
    let mut results = Vec::new();
    prepared.visit(query, options, false, &mut |atoms| {
        results.push(QueryMatch {
            atoms: atoms.to_vec(),
        });
        true
    })?;
    Ok(results)
}
/// Complete collection, retaining every embedding unless `uniquify` is requested.
/// The match cap is a resource bound: exceeding it returns an error, never a partial list.
pub fn find_substructure_matches_complete(
    target: &Molecule,
    query: &QueryGraph,
    options: SubstructureMatchOptions,
) -> Result<Vec<QueryMatch>, SubstructureMatchError> {
    let mut matches = Vec::new();
    visit_substructure_matches(target, query, options, |m| {
        matches.push(m.clone());
        std::ops::ControlFlow::Continue(())
    })?;
    Ok(matches)
}
/// Streams matches. Delivered matches remain provisional until `Complete` is returned.
/// A resource error may follow earlier callbacks; no complete-result claim is then made.
pub fn visit_substructure_matches(
    target: &Molecule,
    query: &QueryGraph,
    options: SubstructureMatchOptions,
    mut visitor: impl FnMut(&QueryMatch) -> std::ops::ControlFlow<()>,
) -> Result<MatchCompletion, SubstructureMatchError> {
    PreparedTarget::new(target).visit(query, options, true, &mut |atoms| {
        visitor(&QueryMatch {
            atoms: atoms.to_vec(),
        })
        .is_continue()
    })
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
