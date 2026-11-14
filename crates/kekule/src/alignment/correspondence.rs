use std::collections::BTreeSet;
use std::fmt;
use std::sync::Arc;

use crate::diagnostics::TopologyMismatch;
use crate::topology::{InstanceAtomId, Topology, TopologyAtomIndex};

use super::AlignmentError;

mod labels;

pub use labels::{AmbiguousLabel, AtomLabel, LabelCorrespondence, LabelScheme};

/// Ordered, one-to-one atom pairs for comparing two topology layouts.
///
/// Pairs associate moving atoms with reference atoms. They may cover a subset
/// of either system and preserve caller order, including the order of weights.
/// This describes an explicit geometric correspondence, not chemical equivalence
/// or a mapping for transferring topology, hierarchy, or annotations.
///
/// Pairs come from the caller ([`Self::from_pairs`]), from dense order when
/// both sides share a layout ([`Self::from_same_layout`]), or from equal
/// hierarchy labels: chain, residue number, insertion code, and atom name
/// ([`Self::from_labels`]).
///
/// No coordinates are copied. Reuse this value across structures or frames sharing
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
    ///
    /// # Errors
    ///
    /// [`AlignmentError::InvalidAtom`] when an atom is absent from its
    /// side's topology, and [`AlignmentError::DuplicateAtom`] when an
    /// atom repeats on its side; both name the first offending pair.
    pub fn from_pairs(
        moving: impl Into<Arc<Topology>>,
        reference: impl Into<Arc<Topology>>,
        pairs: impl IntoIterator<Item = (InstanceAtomId, InstanceAtomId)>,
    ) -> Result<Self, AlignmentError> {
        let moving = &moving.into();
        let reference = &reference.into();
        let mut moving_seen = BTreeSet::new();
        let mut reference_seen = BTreeSet::new();
        let mut indices = Vec::new();
        for (pair_index, (moving_atom, reference_atom)) in pairs.into_iter().enumerate() {
            let resolve = |side, topology: &Topology, atom, seen: &mut BTreeSet<_>| {
                let index = topology
                    .atom_index(atom)
                    .ok_or(AlignmentError::InvalidAtom {
                        side,
                        pair_index,
                        atom,
                    })?;
                if !seen.insert(index) {
                    return Err(AlignmentError::DuplicateAtom {
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
    ///
    /// # Errors
    ///
    /// [`AlignmentError::LayoutMismatch`] when the two topologies do
    /// not have the same layout.
    pub fn from_same_layout(
        moving: impl Into<Arc<Topology>>,
        reference: impl Into<Arc<Topology>>,
    ) -> Result<Self, AlignmentError> {
        let moving = &moving.into();
        let reference = &reference.into();
        if !moving.same_layout(reference) {
            return Err(AlignmentError::LayoutMismatch);
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

    /// The topology the moving atoms belong to.
    pub fn moving_topology(&self) -> &Topology {
        &self.moving
    }

    /// The topology the reference atoms belong to.
    pub fn reference_topology(&self) -> &Topology {
        &self.reference
    }

    /// Number of atom pairs.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    /// Whether there are no atom pairs.
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

    /// Checks, independently on each side, that a topology shares the layout
    /// of the one this correspondence was built for
    /// ([`Topology::shares_layout`]).
    ///
    /// # Errors
    ///
    /// [`AlignmentError::TopologyMismatch`] naming the first side,
    /// [`TopologyMismatch::Moving`] then [`TopologyMismatch::Reference`],
    /// whose topology does not share the layout.
    pub fn ensure_compatible(
        &self,
        moving: &Topology,
        reference: &Topology,
    ) -> Result<(), AlignmentError> {
        for (side, expected, actual) in [
            (TopologyMismatch::Moving, self.moving.as_ref(), moving),
            (
                TopologyMismatch::Reference,
                self.reference.as_ref(),
                reference,
            ),
        ] {
            if !expected.shares_layout(actual) {
                return Err(AlignmentError::TopologyMismatch(side));
            }
        }
        Ok(())
    }
}

/// The side of a moving-to-reference correspondence implicated in an error.
///
/// This enum is closed: a correspondence has exactly two sides.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CorrespondenceSide {
    /// The moving side, whose coordinates are fitted.
    Moving,
    /// The reference side, which stays fixed.
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
