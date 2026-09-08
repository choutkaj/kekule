use std::collections::BTreeSet;
use std::fmt;
use std::sync::Arc;

use crate::topology::{InstanceAtomId, Topology, TopologyAtomIndex};

/// Ordered, one-to-one atom pairs for comparing two exact topology snapshots.
///
/// Pairs associate moving atoms with reference atoms. They may cover a subset
/// of either system and preserve caller order, including the order of weights.
/// This describes an explicit geometric correspondence, not chemical equivalence
/// or a mapping for transferring topology, hierarchy, or annotations.
///
/// No coordinates are copied. Reuse this value across models or frames sharing
/// its two topology allocations. After an operation publishes a new topology
/// snapshot (including perception), construct a new correspondence for that snapshot.
#[derive(Debug, Clone)]
pub struct AtomCorrespondence {
    moving: Arc<Topology>,
    reference: Arc<Topology>,
    pairs: Vec<(TopologyAtomIndex, TopologyAtomIndex)>,
}

impl AtomCorrespondence {
    /// Validates explicit moving/reference atom pairs without sorting them.
    ///
    /// Each atom must resolve on its supplied side and occur at most once on
    /// that side. Elements and bonding are not compared: the caller chooses the
    /// scientific correspondence. Empty or short correspondences are permitted;
    /// each calculation enforces its own minimum size and geometry requirements.
    pub fn from_pairs(
        moving: &Arc<Topology>,
        reference: &Arc<Topology>,
        pairs: impl IntoIterator<Item = (InstanceAtomId, InstanceAtomId)>,
    ) -> Result<Self, AtomCorrespondenceError> {
        let mut moving_seen = BTreeSet::new();
        let mut reference_seen = BTreeSet::new();
        let mut indices = Vec::new();
        for (pair_index, (moving_atom, reference_atom)) in pairs.into_iter().enumerate() {
            let resolve = |side, topology: &Topology, atom, seen: &mut BTreeSet<_>| {
                let index =
                    topology
                        .atom_index(atom)
                        .ok_or(AtomCorrespondenceError::InvalidAtom {
                            side,
                            pair_index,
                            atom,
                        })?;
                if !seen.insert(index) {
                    return Err(AtomCorrespondenceError::DuplicateAtom {
                        side,
                        pair_index,
                        atom,
                    });
                }
                Ok(index)
            };
            let left = resolve(
                CorrespondenceSide::Moving,
                moving,
                moving_atom,
                &mut moving_seen,
            )?;
            let right = resolve(
                CorrespondenceSide::Reference,
                reference,
                reference_atom,
                &mut reference_seen,
            )?;
            indices.push((left, right));
        }
        Ok(Self {
            moving: Arc::clone(moving),
            reference: Arc::clone(reference),
            pairs: indices,
        })
    }

    /// Pairs all atoms by dense order after checking [`Topology::same_layout`].
    ///
    /// This is useful for independent imports of the same record. It does not
    /// infer graph isomorphism, reorder atoms, or resolve symmetry. Perception,
    /// properties, and coordinate differences do not prevent layout matching.
    pub fn from_same_layout(
        moving: &Arc<Topology>,
        reference: &Arc<Topology>,
    ) -> Result<Self, AtomCorrespondenceError> {
        if !moving.same_layout(reference) {
            return Err(AtomCorrespondenceError::LayoutMismatch);
        }
        Ok(Self {
            moving: Arc::clone(moving),
            reference: Arc::clone(reference),
            pairs: moving
                .atom_ids()
                .iter()
                .map(|&atom| {
                    let index = moving.atom_index(atom).expect("published atom index");
                    (index, index)
                })
                .collect(),
        })
    }

    pub fn moving_topology(&self) -> &Topology {
        &self.moving
    }

    pub fn reference_topology(&self) -> &Topology {
        &self.reference
    }

    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// Validated dense pairs in caller order, bound to this value's topologies.
    /// Check [`Self::ensure_compatible`] before indexing another owner's arrays.
    pub fn index_pairs(&self) -> &[(TopologyAtomIndex, TopologyAtomIndex)] {
        &self.pairs
    }

    /// Semantic moving/reference atom pairs in caller order.
    pub fn atom_pairs(
        &self,
    ) -> impl ExactSizeIterator<Item = (InstanceAtomId, InstanceAtomId)> + '_ {
        self.pairs.iter().map(|&(moving, reference)| {
            (
                self.moving.atom_id(moving).expect("validated moving atom"),
                self.reference
                    .atom_id(reference)
                    .expect("validated reference atom"),
            )
        })
    }

    /// Checks exact allocation identity independently on both sides.
    pub fn ensure_compatible(
        &self,
        moving: &Topology,
        reference: &Topology,
    ) -> Result<(), AtomCorrespondenceError> {
        for (side, expected, actual) in [
            (CorrespondenceSide::Moving, self.moving.as_ref(), moving),
            (
                CorrespondenceSide::Reference,
                self.reference.as_ref(),
                reference,
            ),
        ] {
            if !std::ptr::eq(expected, actual) {
                return Err(AtomCorrespondenceError::TopologyMismatch { side });
            }
        }
        Ok(())
    }
}

/// The side of a moving-to-reference correspondence implicated in an error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorrespondenceSide {
    Moving,
    Reference,
}

impl fmt::Display for CorrespondenceSide {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Moving => "moving",
            Self::Reference => "reference",
        })
    }
}

/// Invalid atom pairs or use with different topology snapshots.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum AtomCorrespondenceError {
    LayoutMismatch,
    InvalidAtom {
        side: CorrespondenceSide,
        pair_index: usize,
        atom: InstanceAtomId,
    },
    DuplicateAtom {
        side: CorrespondenceSide,
        pair_index: usize,
        atom: InstanceAtomId,
    },
    TopologyMismatch {
        side: CorrespondenceSide,
    },
}

impl fmt::Display for AtomCorrespondenceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::LayoutMismatch => formatter
                .write_str("topology layouts differ; supply explicit moving/reference atom pairs"),
            Self::InvalidAtom {
                side,
                pair_index,
                atom,
            } => write!(
                formatter,
                "correspondence pair {pair_index} has invalid {side} atom {atom:?}",
            ),
            Self::DuplicateAtom {
                side,
                pair_index,
                atom,
            } => write!(
                formatter,
                "correspondence pair {pair_index} repeats {side} atom {atom:?}",
            ),
            Self::TopologyMismatch { side } => write!(
                formatter,
                "correspondence belongs to a different {side} topology allocation",
            ),
        }
    }
}

impl std::error::Error for AtomCorrespondenceError {}
