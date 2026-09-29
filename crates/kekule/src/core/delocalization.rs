//! Reconstructible conjugation flags and connected resonance groups.
use std::collections::BTreeSet;

use super::{AtomId, BondId};

/// Rules used to assign conjugation, independently of aromaticity's model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConjugationModel {
    /// RDKit 2026.03.3 bond-conjugation rules.
    RdkitLike,
}

/// Installed bond conjugation. Aromatic bonds are also conjugated.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConjugationPerception {
    pub(crate) model: ConjugationModel,
    pub(crate) atoms: BTreeSet<AtomId>,
    pub(crate) bonds: BTreeSet<BondId>,
}

impl ConjugationPerception {
    /// Returns the conjugation model.
    pub const fn model(&self) -> ConjugationModel {
        self.model
    }
    /// Iterates the endpoints of conjugated bonds in stable ID order.
    pub fn atoms(&self) -> impl ExactSizeIterator<Item = AtomId> + DoubleEndedIterator + '_ {
        self.atoms.iter().copied()
    }
    /// Iterates conjugated bonds in stable ID order.
    pub fn bonds(&self) -> impl ExactSizeIterator<Item = BondId> + DoubleEndedIterator + '_ {
        self.bonds.iter().copied()
    }
}

/// One connected component of conjugated bonds, including aromatic bonds.
///
/// Membership does not imply equivalent bond orders, equal contributor weights,
/// or more than one contributor under a particular enumeration policy.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResonanceGroup {
    pub atoms: Vec<AtomId>,
    pub bonds: Vec<BondId>,
}

/// Installed conjugated-group partition. This contains no enumerated structures.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResonancePerception {
    pub(crate) groups: Vec<ResonanceGroup>,
}

impl ResonancePerception {
    /// Returns groups ordered by their smallest bond ID.
    pub fn groups(&self) -> &[ResonanceGroup] {
        &self.groups
    }
    /// Returns the group containing an atom, or `None` outside conjugation.
    pub fn atom_group(&self, atom: AtomId) -> Option<usize> {
        self.groups.iter().position(|g| g.atoms.contains(&atom))
    }
    /// Returns the group containing a bond, or `None` outside conjugation.
    pub fn bond_group(&self, bond: BondId) -> Option<usize> {
        self.groups.iter().position(|g| g.bonds.contains(&bond))
    }
}
