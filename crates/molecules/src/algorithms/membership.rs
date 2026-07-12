use crate::core::*;

use super::rings::compute_ring_membership;

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct RingMembership {
    pub(super) atom_flags: Vec<bool>,
    pub(super) bond_flags: Vec<bool>,
}

impl RingMembership {
    pub fn atom_in_ring(&self, atom: AtomId) -> bool {
        self.atom_flags.get(atom.index()).copied().unwrap_or(false)
    }

    pub fn bond_in_ring(&self, bond: BondId) -> bool {
        self.bond_flags.get(bond.index()).copied().unwrap_or(false)
    }

    pub fn ring_atom_ids(&self) -> impl Iterator<Item = AtomId> + '_ {
        self.atom_flags
            .iter()
            .enumerate()
            .filter_map(|(index, in_ring)| in_ring.then_some(AtomId::new(index as u32)))
    }

    pub fn ring_bond_ids(&self) -> impl Iterator<Item = BondId> + '_ {
        self.bond_flags
            .iter()
            .enumerate()
            .filter_map(|(index, in_ring)| in_ring.then_some(BondId::new(index as u32)))
    }
}

pub fn perceive_ring_membership(mol: &mut Molecule) -> RingMembership {
    let (membership, _) = compute_ring_membership(mol);
    mol.install_ring_membership(membership.clone());
    membership
}
