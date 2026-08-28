use std::fmt;

use crate::algorithms::{validate_stereo, StereoValidationIssue};
use crate::core::*;

mod ranking;

use ranking::{
    assign_cip_element, assign_deferred_tetrahedral_rule6, descriptor_is_absolute_tetrahedral,
    element_is_finally_nonstereogenic, set_stereo_descriptor, CipElementAssignment,
};

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
    if let Err(error) = validate_stereo(mol) {
        return Err(CipAssignmentError {
            issues: error
                .issues
                .into_iter()
                .map(|issue| CipAssignmentIssue::InvalidStereo { issue })
                .collect(),
        });
    }

    let previous_stereo = mol.replace_stereo_perception(Some(StereoPerception::default()));
    let mut report = CipAssignmentReport::default();
    let mut issues = Vec::new();

    let mut pending = mol
        .stereo_elements()
        .map(|(id, element)| (id, element.clone()))
        .collect::<Vec<_>>();

    while !pending.is_empty() {
        let round_mol = mol.clone();
        let mut next_pending = Vec::new();
        let mut round_assignments = Vec::new();
        let mut assigned_this_round = false;
        for (id, element) in pending {
            match assign_cip_element(&round_mol, id, &element, options) {
                CipElementAssignment::Assigned(descriptor) => {
                    round_assignments.push((id, descriptor));
                    assigned_this_round = true;
                }
                CipElementAssignment::Skipped(reason) => {
                    report.skipped.push(CipSkipped {
                        element: id,
                        reason,
                    });
                }
                CipElementAssignment::Deferred => next_pending.push((id, element)),
                CipElementAssignment::Issue(issue) => issues.push(issue),
            }
        }
        for (id, descriptor) in round_assignments {
            set_stereo_descriptor(mol, id, descriptor);
            report.assigned.push(CipAssignment {
                element: id,
                descriptor,
            });
        }
        if !assigned_this_round {
            match assign_deferred_tetrahedral_rule6(mol, &next_pending, options) {
                Ok(assignments) if !assignments.is_empty() => {
                    let has_absolute_assignment = assignments
                        .iter()
                        .any(|(_, descriptor)| descriptor_is_absolute_tetrahedral(*descriptor));
                    let assignments_to_apply = assignments
                        .into_iter()
                        .filter(|(_, descriptor)| {
                            !has_absolute_assignment
                                || descriptor_is_absolute_tetrahedral(*descriptor)
                        })
                        .collect::<Vec<_>>();
                    let assigned_ids = assignments_to_apply
                        .iter()
                        .map(|(id, _)| *id)
                        .collect::<Vec<_>>();
                    for (id, descriptor) in assignments_to_apply {
                        set_stereo_descriptor(mol, id, descriptor);
                        report.assigned.push(CipAssignment {
                            element: id,
                            descriptor,
                        });
                    }
                    pending = next_pending
                        .into_iter()
                        .filter(|(id, _)| !assigned_ids.contains(id))
                        .collect();
                    continue;
                }
                Ok(_) => {}
                Err(issue) => {
                    issues.push(issue);
                    break;
                }
            }
            for (id, element) in next_pending {
                match element_is_finally_nonstereogenic(mol, id, &element, options) {
                    Ok(true) => report.skipped.push(CipSkipped {
                        element: id,
                        reason: CipSkippedReason::NotStereogenic,
                    }),
                    Ok(false) => {
                        issues.push(CipAssignmentIssue::UnresolvedPriority { element: id })
                    }
                    Err(issue) => issues.push(issue),
                }
            }
            break;
        }
        pending = next_pending;
    }
    if issues.is_empty() {
        Ok(report)
    } else {
        drop(mol.replace_stereo_perception(previous_stereo));
        Err(CipAssignmentError { issues })
    }
}
