use crate::core::*;

use super::rings::compute_ring_membership;

pub fn perceive_ring_membership(mol: &mut Molecule) -> RingMembership {
    let (membership, _) = compute_ring_membership(mol);
    mol.install_ring_membership(membership.clone());
    membership
}
