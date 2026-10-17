//! Closed-shell nitrogen inversion eligibility under the RDKit-like convention.
//! This is a discrete stereo model, not a prediction of an inversion barrier.
use super::*;
use crate::algorithms::aromaticity::{count_rdkit_like_atom_pi_electrons, rdkit_outer_electrons};
use crate::algorithms::valence::rdkit_default_valence;

pub(super) fn classified(mol: &Molecule, center: AtomId) -> bool {
    let atom = mol.atom(center).expect("live stereo center");
    atom.element.symbol() == "N"
        && atom.formal_charge == 0
        && atom.radical.is_none()
        && mol
            .incident_bonds(center)
            .expect("live center")
            .all(|(_, bond)| {
                matches!(
                    bond.order,
                    BondOrder::Single | BondOrder::Double | BondOrder::Triple
                )
            })
}

pub(super) fn lone_pair_is_stereogenic(mol: &Molecule, center: AtomId) -> bool {
    if !classified(mol, center)
        || mol
            .incident_bonds(center)
            .expect("live center")
            .any(|(id, bond)| {
                bond.order != BondOrder::Single
                    || mol
                        .perception()
                        .aromaticity()
                        .map(|aromaticity| aromaticity.bond_is_aromatic(id))
                        == Some(true)
                    || conjugated_neighbor(mol, center, bond.other_atom(center))
            })
    {
        return false;
    }
    // A three-membered ring can be identified directly, independently of any
    // selected ring basis. Hydrogen storage does not change ligand geometry.
    let neighbors: Vec<_> = mol.neighbors(center).expect("live center").collect();
    if neighbors.iter().enumerate().any(|(i, &a)| {
        neighbors[i + 1..].iter().any(|&b| {
            mol.bond_between(a, b).is_some_and(|id| {
                !matches!(
                    mol.bond(id).expect("live bond").order,
                    BondOrder::Zero | BondOrder::Dative
                )
            })
        })
    }) {
        return true;
    }
    bridgehead(mol, center)
}

fn degree(mol: &Molecule, atom: AtomId) -> usize {
    mol.neighbors(atom).expect("live atom").count() + usize::from(atom_hydrogen_count(mol, atom))
}

fn conjugation_donor(mol: &Molecule, id: AtomId) -> bool {
    let atom = mol.atom(id).expect("live atom");
    let valence =
        crate::algorithms::explicit_valence(mol, id) + usize::from(atom_hydrogen_count(mol, id));
    if atom.formal_charge == 0
        && rdkit_default_valence(atom).is_some_and(|minimum| valence > usize::from(minimum))
    {
        return false;
    }
    let outer = rdkit_outer_electrons(atom);
    (atom.element.atomic_number() <= 10
        || !matches!(outer, 5 | 6)
        || outer == 6 && degree(mol, id) < 2)
        && count_rdkit_like_atom_pi_electrons(mol, id, atom).is_some_and(|n| n > 0)
}

fn conjugated_neighbor(mol: &Molecule, center: AtomId, neighbor: AtomId) -> bool {
    if !(2..=3).contains(&degree(mol, neighbor)) || !conjugation_donor(mol, neighbor) {
        return false;
    }
    mol.incident_bonds(neighbor)
        .expect("live neighbor")
        .any(|(id, bond)| {
            let other = bond.other_atom(neighbor);
            other != center
                && (matches!(bond.order, BondOrder::Double | BondOrder::Triple)
                    || mol
                        .perception()
                        .aromaticity()
                        .map(|aromaticity| aromaticity.bond_is_aromatic(id))
                        == Some(true))
                && conjugation_donor(mol, other)
        })
}

fn bridgehead(mol: &Molecule, center: AtomId) -> bool {
    let membership = mol
        .perception()
        .rings()
        .map(|rings| rings.membership())
        .expect("prepared ring membership");
    if mol
        .incident_bonds(center)
        .expect("live center")
        .filter(|(id, _)| membership.bond_in_ring(*id))
        .count()
        < 3
    {
        return false;
    }
    let rings: Vec<_> = mol
        .perception()
        .rings()
        .and_then(|rings| rings.ring_set())
        .expect("prepared ring set")
        .rings()
        .iter()
        .filter(|ring| ring.atoms.contains(&center))
        .collect();
    // RDKit's bridgehead convention: every selected ring at the center must
    // overlap another such ring by at least two bonds. Three incident ring
    // bonds alone would also accept ordinary fused-ring junctions.
    !rings.is_empty()
        && rings.iter().enumerate().all(|(i, ring)| {
            rings.iter().enumerate().any(|(j, other)| {
                i != j
                    && ring
                        .bonds
                        .iter()
                        .filter(|bond| other.bonds.contains(bond))
                        .take(2)
                        .count()
                        == 2
            })
        })
}
