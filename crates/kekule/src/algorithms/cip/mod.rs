use crate::diagnostics::LimitExceeded;
use std::fmt;

use crate::algorithms::StereoValidationIssue;
use crate::molecule::*;

mod ranking;

type CipResult<T> = std::result::Result<T, CipAssignmentIssue>;
type RankingResult<T> = std::result::Result<T, CipRankingError>;

/// Failure to compare rooted ligands without assigning a stereo descriptor.
/// A resource limit or unresolved comparison is never an equal-priority result.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum CipRankingError {
    /// Expansion exceeded [`CipAssignmentOptions::max_depth`]
    /// ([`Resource::DigraphDepth`](crate::diagnostics::Resource::DigraphDepth)) or [`CipAssignmentOptions::max_nodes`]
    /// ([`Resource::DigraphNodes`](crate::diagnostics::Resource::DigraphNodes)).
    LimitExceeded(LimitExceeded),
    UnresolvedPriority,
    InvalidStereo {
        issue: StereoValidationIssue,
    },
    /// Run valence perception explicitly before ranking ligands whose inferred
    /// hydrogen counts are unknown. No perception is installed by ranking.
    UnknownHydrogenCount {
        atom: AtomId,
    },
}

impl fmt::Display for CipRankingError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LimitExceeded(limit) => write!(f, "CIP ranking {limit}"),
            Self::UnresolvedPriority => f.write_str("CIP ligand priority could not be resolved"),
            Self::InvalidStereo { issue } => write!(f, "invalid stereo for CIP ranking: {issue:?}"),
            Self::UnknownHydrogenCount { atom } => write!(
                f,
                "CIP ranking requires known implicit hydrogen counts: {atom}"
            ),
        }
    }
}

impl std::error::Error for CipRankingError {}

impl CipRankingError {
    fn for_element(self, element: StereoElementId) -> CipAssignmentIssue {
        match self {
            Self::LimitExceeded(limit) => CipAssignmentIssue::LimitExceeded { element, limit },
            Self::InvalidStereo { issue } => CipAssignmentIssue::InvalidStereo { issue },
            Self::UnknownHydrogenCount { atom } => {
                CipAssignmentIssue::UnknownHydrogenCount { atom }
            }
            Self::UnresolvedPriority => CipAssignmentIssue::UnresolvedPriority { element },
        }
    }
}

pub(crate) use ranking::bond_reference_atoms;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Explicit bounds on rooted ligand expansion.
pub struct CipAssignmentOptions {
    /// Maximum number of bonds explored beyond a carrier atom.
    ///
    /// Constitutional comparisons retain their digraphs and compare at
    /// increasing depths, also checking the last available shell before a
    /// bound prevents expansion. Rule 1a can prove an ordering without a
    /// complete digraph. A tied, truncated digraph cannot advance to
    /// later sequence rules or establish that a center is nonstereogenic.
    /// Defaults to 64; the independent node bound also applies.
    pub max_depth: usize,
    /// Maximum nodes in one ligand expansion or shared auxiliary digraph.
    pub max_nodes: usize,
}

impl CipAssignmentOptions {
    /// The [`Default`] bound of [`Self::max_depth`].
    pub const DEFAULT_MAX_DEPTH: usize = 64;
    /// The [`Default`] bound of [`Self::max_nodes`].
    pub const DEFAULT_MAX_NODES: usize = 100_000;
}

impl Default for CipAssignmentOptions {
    fn default() -> Self {
        Self {
            max_depth: Self::DEFAULT_MAX_DEPTH,
            max_nodes: Self::DEFAULT_MAX_NODES,
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CipAssignmentReport {
    pub assigned: Vec<CipAssignment>,
    pub skipped: Vec<CipSkipped>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CipAssignmentError {
    pub issues: Vec<CipAssignmentIssue>,
}

impl fmt::Display for CipAssignmentError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "CIP assignment reported {} issue(s)",
            self.issues.len()
        )
    }
}

impl std::error::Error for CipAssignmentError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CipAssignment {
    pub element: StereoElementId,
    pub descriptor: StereoDescriptor,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CipSkipped {
    pub element: StereoElementId,
    pub reason: CipSkippedReason,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CipSkippedReason {
    UnknownConfiguration,
    NotStereogenic,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CipAssignmentIssue {
    InvalidStereo {
        issue: StereoValidationIssue,
    },
    UnresolvedPriority {
        element: StereoElementId,
    },
    /// Ranking the element's ligands exceeded a digraph bound.
    LimitExceeded {
        /// The stereo element being ranked.
        element: StereoElementId,
        /// The exceeded bound.
        limit: LimitExceeded,
    },
    /// The atom's implicit hydrogen count is inferred but not yet perceived.
    /// Missing valence perception is not zero hydrogens, so no ranking ran.
    UnknownHydrogenCount {
        atom: AtomId,
    },
}

/// Assigns CIP descriptors from represented local stereochemistry.
///
/// Uses ordered sequence rules 1a–6, path-dependent auxiliary descriptors, and
/// fractional atomic numbers for mancude systems. Existing CIP perception is
/// replaced only after every requested assignment succeeds. Unknown and
/// nonstereogenic configurations are reported as skipped.
/// Explicit bond configurations are ranked regardless of ring size,
/// aromaticity, or endpoint elements; those restrictions belong to stereo
/// perception rather than descriptor assignment.
///
/// Ligand expansion counts every implicit hydrogen, so every atom's implicit
/// hydrogen count must be known: fixed, or inferred by valence perception
/// (`Molecule::perceive` or the valence facade). Otherwise assignment reports
/// [`CipAssignmentIssue::UnknownHydrogenCount`] for each such atom, ranks
/// nothing, and leaves perception unchanged.
///
/// `options` bounds the expansion. Exceeding a bound returns an issue and
/// preserves the previous perception; no label is inferred from a truncated
/// tie.
pub fn assign_cip_descriptors(
    mol: &mut Molecule,
    options: CipAssignmentOptions,
) -> std::result::Result<CipAssignmentReport, CipAssignmentError> {
    ranking::assign_cip_descriptors(mol, options)
}
