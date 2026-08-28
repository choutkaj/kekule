use std::fmt;

use crate::algorithms::StereoValidationIssue;
use crate::core::*;

mod ranking;

type CipResult<T> = std::result::Result<T, CipAssignmentIssue>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CipAssignmentOptions {
    pub max_depth: usize,
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
}

pub fn assign_cip_descriptors(
    mol: &mut Molecule,
) -> std::result::Result<CipAssignmentReport, CipAssignmentError> {
    assign_cip_descriptors_with_options(mol, CipAssignmentOptions::default())
}

pub fn assign_cip_descriptors_with_options(
    mol: &mut Molecule,
    options: CipAssignmentOptions,
) -> std::result::Result<CipAssignmentReport, CipAssignmentError> {
    ranking::assign_cip_descriptors_with_options(mol, options)
}
