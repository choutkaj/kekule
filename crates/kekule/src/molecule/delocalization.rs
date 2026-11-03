//! Reconstructible conjugation flags and resonance delocalization.
use std::collections::{BTreeMap, BTreeSet};

use super::{AtomId, BondId};

/// Rules used to assign conjugation, independently of aromaticity's model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
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
    /// Whether one atom is an endpoint of a conjugated bond.
    pub fn atom_is_conjugated(&self, atom: AtomId) -> bool {
        self.atoms.contains(&atom)
    }
    /// Whether one bond is conjugated; aromatic bonds are.
    pub fn bond_is_conjugated(&self, bond: BondId) -> bool {
        self.bonds.contains(&bond)
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

/// Rules used to perceive resonance delocalization.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum ResonanceModel {
    /// Contributors that move one formal charge along an alternating path of
    /// non-aromatic bonds, keep every octet as it is, and add no charge; see
    /// [`perceive_resonance`](crate::perception::resonance::perceive_resonance).
    ChargeMoves,
}

/// Connected atoms and bonds: one delocalized group of installed resonance
/// perception, or one connected component of conjugated bonds that
/// contributor enumeration varies.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResonanceGroup {
    /// The endpoints of the group's bonds, in ascending ID order.
    pub atoms: Vec<AtomId>,
    /// The group's bonds, in ascending ID order.
    pub bonds: Vec<BondId>,
}

/// Installed resonance delocalization: the bonds whose order differs between
/// the contributors of their group, partitioned into connected groups.
///
/// Contributors are never enumerated or stored. Aromatic bonds are left to
/// aromaticity and never belong to a group.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResonancePerception {
    pub(crate) model: ResonanceModel,
    pub(crate) groups: Vec<ResonanceGroup>,
    pub(crate) atom_groups: BTreeMap<AtomId, usize>,
    pub(crate) bond_groups: BTreeMap<BondId, usize>,
}

impl ResonancePerception {
    /// Orders `groups` by smallest bond ID and indexes their members.
    pub(crate) fn new(model: ResonanceModel, mut groups: Vec<ResonanceGroup>) -> Self {
        groups.sort_by_key(|group| group.bonds.first().copied());
        let mut atom_groups = BTreeMap::new();
        let mut bond_groups = BTreeMap::new();
        for (index, group) in groups.iter().enumerate() {
            atom_groups.extend(group.atoms.iter().map(|&atom| (atom, index)));
            bond_groups.extend(group.bonds.iter().map(|&bond| (bond, index)));
        }
        Self {
            model,
            groups,
            atom_groups,
            bond_groups,
        }
    }
    /// Returns the resonance model.
    pub const fn model(&self) -> ResonanceModel {
        self.model
    }
    /// Returns the delocalized groups, ordered by their smallest bond ID.
    pub fn groups(&self) -> &[ResonanceGroup] {
        &self.groups
    }
    /// Whether resonance delocalizes one bond.
    pub fn bond_is_delocalized(&self, bond: BondId) -> bool {
        self.bond_groups.contains_key(&bond)
    }
    /// Iterates the delocalized bonds in ascending ID order.
    pub fn bonds(&self) -> impl ExactSizeIterator<Item = BondId> + DoubleEndedIterator + '_ {
        self.bond_groups.keys().copied()
    }
    /// Returns the index of the group containing an atom, or `None` when no
    /// delocalized bond ends at it.
    pub fn atom_group(&self, atom: AtomId) -> Option<usize> {
        self.atom_groups.get(&atom).copied()
    }
    /// Returns the index of the group containing a bond, or `None` when the
    /// bond is not delocalized.
    pub fn bond_group(&self, bond: BondId) -> Option<usize> {
        self.bond_groups.get(&bond).copied()
    }
}
