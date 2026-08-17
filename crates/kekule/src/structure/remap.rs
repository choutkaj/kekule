use std::fmt;
use std::sync::Arc;

use crate::topology::{
    InstanceAtomId, InstanceBondId, Topology, TopologyAtomIndex, TopologyBondIndex, TopologyMapping,
};

/// Validates that `mapping` transfers complete per-atom state from `source`
/// to `target`.
///
/// Added or otherwise unmapped target atoms are rejected because complete
/// dense state cannot invent values for them.
pub fn validate_complete_atom_mapping(
    source: &Arc<Topology>,
    target: &Arc<Topology>,
    mapping: &TopologyMapping,
) -> Result<(), TopologyRemapError> {
    if !mapping.is_source(source) {
        return Err(TopologyRemapError::MappingSourceMismatch);
    }
    if !mapping.is_target(target) {
        return Err(TopologyRemapError::MappingTargetMismatch);
    }
    if let Some(target_atom) = mapping.added_atoms().first().copied() {
        return Err(TopologyRemapError::AddedAtomsRequireState { target_atom });
    }
    if mapping.atom_index_pairs().len() != target.atom_count() {
        let mapped = mapping
            .atom_pairs()
            .map(|(_, target)| target)
            .collect::<std::collections::BTreeSet<_>>();
        let target_atom = target
            .atom_ids()
            .iter()
            .copied()
            .find(|atom| !mapped.contains(atom))
            .expect("incomplete target mapping has an unmapped target atom");
        return Err(TopologyRemapError::AddedAtomsRequireState { target_atom });
    }
    Ok(())
}

/// Remaps one complete dense atom array into the target topology's
/// authoritative dense order.
pub fn dense_atom_values<T: Clone>(
    source_values: &[T],
    source: &Arc<Topology>,
    target: &Arc<Topology>,
    mapping: &TopologyMapping,
) -> Result<Vec<T>, TopologyRemapError> {
    if source_values.len() != source.atom_count() {
        return Err(TopologyRemapError::SourceAtomCountMismatch {
            expected: source.atom_count(),
            actual: source_values.len(),
        });
    }
    validate_complete_atom_mapping(source, target, mapping)?;

    let mut values = std::iter::repeat_with(|| None)
        .take(target.atom_count())
        .collect::<Vec<Option<T>>>();
    for (source_index, target_index) in mapping.atom_index_pairs() {
        let slot = &mut values[target_index.index()];
        if slot.is_some() {
            return Err(TopologyRemapError::DuplicateTargetAssignment { target_index });
        }
        *slot = Some(source_values[source_index.index()].clone());
    }
    values
        .into_iter()
        .enumerate()
        .map(|(index, value)| {
            value.ok_or_else(|| TopologyRemapError::AddedAtomsRequireState {
                target_atom: target.atom_ids()[index],
            })
        })
        .collect()
}

/// Validates that `mapping` transfers complete per-bond state from
/// `source` to `target`.
pub fn validate_complete_bond_mapping(
    source: &Arc<Topology>,
    target: &Arc<Topology>,
    mapping: &TopologyMapping,
) -> Result<(), TopologyRemapError> {
    if !mapping.is_source(source) {
        return Err(TopologyRemapError::MappingSourceMismatch);
    }
    if !mapping.is_target(target) {
        return Err(TopologyRemapError::MappingTargetMismatch);
    }
    if let Some(target_bond) = mapping.added_bonds().first().copied() {
        return Err(TopologyRemapError::AddedBondsRequireState { target_bond });
    }
    if mapping.bond_index_pairs().len() != target.bond_count() {
        let mapped = mapping
            .bond_pairs()
            .map(|(_, target)| target)
            .collect::<std::collections::BTreeSet<_>>();
        let target_bond = target
            .bond_ids()
            .iter()
            .copied()
            .find(|bond| !mapped.contains(bond))
            .expect("incomplete target mapping has an unmapped target bond");
        return Err(TopologyRemapError::AddedBondsRequireState { target_bond });
    }
    Ok(())
}

/// Remaps one complete dense bond array into the target topology's
/// authoritative dense order.
pub fn dense_bond_values<T: Clone>(
    source_values: &[T],
    source: &Arc<Topology>,
    target: &Arc<Topology>,
    mapping: &TopologyMapping,
) -> Result<Vec<T>, TopologyRemapError> {
    if source_values.len() != source.bond_count() {
        return Err(TopologyRemapError::SourceBondCountMismatch {
            expected: source.bond_count(),
            actual: source_values.len(),
        });
    }
    validate_complete_bond_mapping(source, target, mapping)?;

    let mut values = std::iter::repeat_with(|| None)
        .take(target.bond_count())
        .collect::<Vec<Option<T>>>();
    for (source_index, target_index) in mapping.bond_index_pairs() {
        let slot = &mut values[target_index.index()];
        if slot.is_some() {
            return Err(TopologyRemapError::DuplicateTargetBondAssignment { target_index });
        }
        *slot = Some(source_values[source_index.index()].clone());
    }
    values
        .into_iter()
        .enumerate()
        .map(|(index, value)| {
            value.ok_or_else(|| TopologyRemapError::AddedBondsRequireState {
                target_bond: target.bond_ids()[index],
            })
        })
        .collect()
}

/// Failure to remap complete topology-bound structure state.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TopologyRemapError {
    /// The source state is not bound to the supplied source topology.
    SourceTopologyMismatch,
    /// The destination state is not bound to the supplied target topology.
    TargetTopologyMismatch,
    /// The mapping is not sourced from the supplied source topology.
    MappingSourceMismatch,
    /// The mapping does not target the supplied target topology.
    MappingTargetMismatch,
    /// A complete source dense array has an invalid length.
    SourceAtomCountMismatch { expected: usize, actual: usize },
    /// A target atom has no source state under this mapping.
    AddedAtomsRequireState { target_atom: InstanceAtomId },
    /// A complete source dense bond array has an invalid length.
    SourceBondCountMismatch { expected: usize, actual: usize },
    /// A target bond has no source state under this mapping.
    AddedBondsRequireState { target_bond: InstanceBondId },
    /// More than one source value was assigned to one target dense index.
    DuplicateTargetAssignment { target_index: TopologyAtomIndex },
    /// More than one source value was assigned to one target dense bond index.
    DuplicateTargetBondAssignment { target_index: TopologyBondIndex },
    /// One ensemble member could not be remapped.
    Member {
        member: usize,
        error: Box<TopologyRemapError>,
    },
}

impl fmt::Display for TopologyRemapError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SourceTopologyMismatch => {
                formatter.write_str("source state does not belong to the supplied topology")
            }
            Self::TargetTopologyMismatch => {
                formatter.write_str("target state does not belong to the supplied topology")
            }
            Self::MappingSourceMismatch => {
                formatter.write_str("topology mapping does not match the source topology")
            }
            Self::MappingTargetMismatch => {
                formatter.write_str("topology mapping does not match the target topology")
            }
            Self::SourceAtomCountMismatch { expected, actual } => write!(
                formatter,
                "source state requires {expected} atoms, but received {actual}"
            ),
            Self::AddedAtomsRequireState { target_atom } => write!(
                formatter,
                "target atom {target_atom} has no mapped source state"
            ),
            Self::SourceBondCountMismatch { expected, actual } => write!(
                formatter,
                "source state requires {expected} bonds, but received {actual}"
            ),
            Self::AddedBondsRequireState { target_bond } => write!(
                formatter,
                "target bond {target_bond} has no mapped source state"
            ),
            Self::DuplicateTargetAssignment { target_index } => {
                write!(formatter, "target {target_index} received duplicate state")
            }
            Self::DuplicateTargetBondAssignment { target_index } => {
                write!(formatter, "target {target_index} received duplicate state")
            }
            Self::Member { member, error } => {
                write!(formatter, "cannot remap ensemble member {member}: {error}")
            }
        }
    }
}

impl std::error::Error for TopologyRemapError {}
