use std::fmt;

use crate::algorithms::StereoValidationIssue;
use crate::core::*;

mod ranking;

type CipResult<T> = std::result::Result<T, CipAssignmentIssue>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Explicit bounds on rooted ligand expansion.
pub struct CipAssignmentOptions {
    /// Maximum number of bonds explored beyond a carrier atom.
    ///
    /// Constitutional comparisons expand progressively, stopping as soon as
    /// Rule 1a proves the ordering. A tied, truncated digraph cannot advance to
    /// later sequence rules or establish that a center is nonstereogenic.
    pub max_depth: usize,
    /// Maximum nodes in one ligand expansion or shared auxiliary digraph.
    pub max_nodes: usize,
}

impl Default for CipAssignmentOptions {
    fn default() -> Self {
        Self {
            max_depth: 32,
            max_nodes: 100_000,
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
    ResourceLimitExceeded {
        element: StereoElementId,
        max_nodes: usize,
    },
    DepthLimitExceeded {
        element: StereoElementId,
        max_depth: usize,
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
pub fn assign_cip_descriptors(
    mol: &mut Molecule,
) -> std::result::Result<CipAssignmentReport, CipAssignmentError> {
    assign_cip_descriptors_with_options(mol, CipAssignmentOptions::default())
}

/// Assigns CIP descriptors with explicit expansion bounds.
///
/// Exceeding a bound returns an issue and preserves the previous perception;
/// no label is inferred from a truncated tie.
pub fn assign_cip_descriptors_with_options(
    mol: &mut Molecule,
    options: CipAssignmentOptions,
) -> std::result::Result<CipAssignmentReport, CipAssignmentError> {
    ranking::assign_cip_descriptors_with_options(mol, options)
}
