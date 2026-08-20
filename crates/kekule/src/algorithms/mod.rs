mod aromaticity;
mod canonical;
mod cip;
mod hydrogens;
mod membership;
mod rings;
mod rotatable_bonds;
mod stereo;
mod substructure;
mod valence;

pub(crate) use crate::core::{RingMembership, ValenceModel};
pub use aromaticity::*;
pub use canonical::*;
pub use cip::*;
pub use hydrogens::*;
pub use membership::*;
pub use rings::*;
pub use rotatable_bonds::*;
pub use stereo::*;
pub(crate) use stereo::{
    atom_axis_carriers, atom_hydrogen_count, coordinates_are_planar,
    double_bond_between_aromatic_atoms, double_bond_endpoint_carriers,
    double_bond_has_noncarbon_endpoint, double_bond_is_in_ring,
    tetrahedral_orientation_from_points, tetrahedral_points,
};

pub(crate) fn compute_graph_ring_membership(
    molecule: &crate::core::Molecule,
) -> crate::core::RingMembership {
    rings::compute_ring_membership(molecule)
}

pub(crate) fn graph_bond_in_ring_smaller_than(
    molecule: &crate::core::Molecule,
    bond: crate::core::BondId,
    ring_size: usize,
) -> bool {
    rings::bond_in_ring_smaller_than(molecule, bond, ring_size)
}
pub use substructure::*;
pub use valence::*;
