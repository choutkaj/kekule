use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use super::*;
use crate::core::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AromaticityModel {
    RdkitLike,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AromaticityError {
    UnsupportedElement(AtomId),
    InvalidAromaticRepresentation(AtomId),
    RingPerception(RingPerceptionError),
}

impl fmt::Display for AromaticityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedElement(id) => {
                write!(f, "unsupported aromaticity element at atom {id}")
            }
            Self::InvalidAromaticRepresentation(id) => {
                write!(f, "invalid aromatic representation at atom {id}")
            }
            Self::RingPerception(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for AromaticityError {}

pub fn perceive_aromaticity(
    mol: &mut Molecule,
    model: AromaticityModel,
) -> std::result::Result<(), AromaticityError> {
    perceive_aromaticity_with_ring_options(mol, model, RingPerceptionOptions::default())
}

pub fn perceive_aromaticity_with_ring_options(
    mol: &mut Molecule,
    model: AromaticityModel,
    ring_options: RingPerceptionOptions,
) -> std::result::Result<(), AromaticityError> {
    match model {
        AromaticityModel::RdkitLike => perceive_rdkit_like_aromaticity(mol, ring_options),
    }
}

fn perceive_rdkit_like_aromaticity(
    mol: &mut Molecule,
    ring_options: RingPerceptionOptions,
) -> std::result::Result<(), AromaticityError> {
    mol.perception.aromaticity = invalidate(mol.perception.aromaticity);
    mol.perception.stereo = invalidate(mol.perception.stereo);
    let imported_aromatic_components = imported_aromatic_bond_components(mol);
    let imported_explicit_aromatic_singles = mol
        .bonds()
        .filter_map(|(bond_id, bond)| {
            (matches!(bond.order, BondOrder::Single)
                && mol.atom(bond.a()).is_ok_and(|atom| atom.aromatic)
                && mol.atom(bond.b()).is_ok_and(|atom| atom.aromatic))
            .then_some(bond_id)
        })
        .collect::<BTreeSet<_>>();
    let ring_set = perceive_ring_set_with_options(mol, ring_options)
        .map_err(AromaticityError::RingPerception)?;
    for atom in mol.atoms.iter_mut().flatten() {
        atom.aromatic = false;
    }
    for bond in mol.bonds.iter_mut().flatten() {
        bond.aromatic = false;
    }

    let ring_aromatic = ring_set
        .rings()
        .iter()
        .map(|ring| {
            if ring.atoms.len() > 7 {
                return Ok(false);
            }
            let electrons = aromatic_ring_pi_electrons(mol, ring)?;
            Ok(electrons >= 6 && (electrons - 2) % 4 == 0)
        })
        .collect::<std::result::Result<Vec<_>, AromaticityError>>()?;
    let non_aromatic_fusion_singles = ring_set
        .rings()
        .iter()
        .flat_map(|ring| ring.bonds.iter().copied())
        .filter(|bond_id| {
            mol.bond(*bond_id).is_ok_and(|bond| {
                matches!(bond.order, BondOrder::Single)
                    && mol
                        .atom(bond.a())
                        .is_ok_and(|atom| atom.element.symbol() == "C")
                    && mol
                        .atom(bond.b())
                        .is_ok_and(|atom| atom.element.symbol() == "C")
            })
        })
        .filter(|bond_id| {
            let containing_rings = ring_set
                .rings()
                .iter()
                .enumerate()
                .filter(|(_, ring)| ring.bonds.contains(bond_id))
                .collect::<Vec<_>>();
            if containing_rings
                .iter()
                .any(|(_, ring)| fused_component_is_all_carbon(mol, ring))
            {
                return false;
            }
            containing_rings.iter().any(|(index, ring)| {
                if ring_aromatic[*index] {
                    return false;
                }
                let non_donor_five_ring = ring.atoms.len() == 5
                    && !fused_component_is_all_carbon(mol, ring)
                    && !ring_has_chalcogen_donor(mol, ring)
                    && ring_hetero_donor_count(mol, ring) < 2
                    && !ring_has_nitrogen_lone_pair_donor(mol, ring);
                let multi_hetero_dione_ring = ring.atoms.len() == 6
                    && containing_rings.len() > 1
                    && ring_hetero_donor_count(mol, ring) >= 2
                    && ring_terminal_exocyclic_pi_bond_count(mol, ring) >= 2;
                non_donor_five_ring || multi_hetero_dione_ring
            })
        })
        .collect::<BTreeSet<_>>();
    let protected_non_aromatic_bonds = imported_explicit_aromatic_singles
        .union(&non_aromatic_fusion_singles)
        .copied()
        .collect::<BTreeSet<_>>();

    for (ring, aromatic) in ring_set.rings().iter().zip(ring_aromatic) {
        if aromatic {
            for atom_id in &ring.atoms {
                if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
                    atom.aromatic = true;
                }
            }
            for bond_id in &ring.bonds {
                if let Some(bond) = mol.bonds[bond_id.index()].as_mut() {
                    bond.aromatic = !protected_non_aromatic_bonds.contains(bond_id);
                }
            }
        }
    }
    perceive_fused_aromatic_components(mol, ring_set.rings(), &protected_non_aromatic_bonds)?;
    perceive_fused_single_exocyclic_carbon_rings(
        mol,
        ring_set.rings(),
        &protected_non_aromatic_bonds,
    );
    clear_terminal_chalcogen_oxo_ring_atoms(mol, ring_set.rings());
    clear_ring_oxo_chalcogen_atoms(mol, ring_set.rings());
    clear_fused_lactam_bridge_ring_atoms(mol, ring_set.rings());
    clear_imide_carbonyl_ring_atoms(mol, ring_set.rings());
    clear_fused_lactone_bridge_ring_atoms(mol, ring_set.rings());
    clear_saturated_tertiary_amine_ring_atoms(mol, ring_set.rings());
    clear_saturated_chalcogen_bridge_atoms(mol, ring_set.rings());
    clear_fused_carbonyl_bridge_atoms(mol, ring_set.rings());
    clear_saturated_fused_carbon_ring_atoms(mol, ring_set.rings());
    clear_terminal_aromatic_imine_fragments(mol);
    clear_aromatic_amidine_carbon_atoms(mol, ring_set.rings());
    clear_saturated_aromatic_carbon_atoms(mol);
    clear_orphan_aromatic_atoms(mol);
    for component in imported_aromatic_components {
        if !component.iter().any(|atom_id| {
            mol.atom(*atom_id)
                .map(|atom| atom.aromatic)
                .unwrap_or(false)
        }) {
            return Err(AromaticityError::InvalidAromaticRepresentation(
                component[0],
            ));
        }
        restore_imported_aromatic_component(mol, &component);
    }

    mol.perception.aromaticity = ComputedState::Fresh;
    Ok(())
}

fn restore_imported_aromatic_component(mol: &mut Molecule, component: &[AtomId]) {
    let atoms = component.iter().copied().collect::<BTreeSet<_>>();
    for atom_id in component {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = true;
        }
    }
    for bond in mol.bonds.iter_mut().flatten() {
        if matches!(bond.order, BondOrder::Aromatic)
            && atoms.contains(&bond.a())
            && atoms.contains(&bond.b())
        {
            bond.aromatic = true;
        }
    }
}

fn imported_aromatic_bond_components(mol: &Molecule) -> Vec<Vec<AtomId>> {
    let mut adjacency = BTreeMap::<AtomId, Vec<AtomId>>::new();
    for (_, bond) in mol
        .bonds()
        .filter(|(_, bond)| matches!(bond.order, BondOrder::Aromatic))
    {
        adjacency.entry(bond.a()).or_default().push(bond.b());
        adjacency.entry(bond.b()).or_default().push(bond.a());
    }

    let mut components = Vec::new();
    let mut visited = BTreeSet::new();
    for start in adjacency.keys().copied() {
        if !visited.insert(start) {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![start];
        while let Some(atom_id) = stack.pop() {
            component.push(atom_id);
            if let Some(neighbors) = adjacency.get(&atom_id) {
                for neighbor in neighbors.iter().rev().copied() {
                    if visited.insert(neighbor) {
                        stack.push(neighbor);
                    }
                }
            }
        }
        component.sort();
        components.push(component);
    }
    components
}

fn clear_saturated_chalcogen_bridge_atoms(mol: &mut Molecule, rings: &[Ring]) {
    let atoms_to_clear = rings
        .iter()
        .enumerate()
        .filter(|(index, ring)| {
            ring.atoms.len() >= 6
                && rings.iter().enumerate().any(|(other_index, other)| {
                    other_index != *index && rings_share_bond(ring, other)
                })
        })
        .flat_map(|(_, ring)| {
            ring.atoms
                .iter()
                .copied()
                .filter(|atom_id| is_saturated_chalcogen_bridge(mol, ring, *atom_id))
        })
        .collect::<BTreeSet<_>>();

    for atom_id in &atoms_to_clear {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = false;
        }
    }
    for bond in mol.bonds.iter_mut().flatten() {
        if atoms_to_clear.contains(&bond.a()) || atoms_to_clear.contains(&bond.b()) {
            bond.aromatic = false;
        }
    }
}

fn clear_orphan_aromatic_atoms(mol: &mut Molecule) {
    let atoms_to_clear = mol
        .atoms()
        .filter_map(|(atom_id, atom)| {
            (atom.aromatic
                && mol
                    .incident_bonds(atom_id)
                    .map_or(true, |mut bonds| !bonds.any(|(_, bond)| bond.aromatic)))
            .then_some(atom_id)
        })
        .collect::<BTreeSet<_>>();

    for atom_id in atoms_to_clear {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = false;
        }
    }
}

fn clear_terminal_aromatic_imine_fragments(mol: &mut Molecule) {
    let atoms_to_clear = mol
        .bonds()
        .filter_map(|(_, bond)| {
            if !bond.aromatic || !matches!(bond.order, BondOrder::Double) {
                return None;
            }
            let left = mol.atom(bond.a()).ok()?;
            let right = mol.atom(bond.b()).ok()?;
            let (nitrogen_id, carbon_id) = match (left.element.symbol(), right.element.symbol()) {
                ("N", "C") => (bond.a(), bond.b()),
                ("C", "N") => (bond.b(), bond.a()),
                _ => return None,
            };
            let carbon_aromatic_bonds = aromatic_incident_bond_count(mol, carbon_id);
            let nitrogen_aromatic_bonds = aromatic_incident_bond_count(mol, nitrogen_id);
            (carbon_aromatic_bonds <= 1 && nitrogen_aromatic_bonds <= 2)
                .then_some([nitrogen_id, carbon_id])
        })
        .flatten()
        .collect::<BTreeSet<_>>();

    for atom_id in &atoms_to_clear {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = false;
        }
    }
    for bond in mol.bonds.iter_mut().flatten() {
        if atoms_to_clear.contains(&bond.a()) || atoms_to_clear.contains(&bond.b()) {
            bond.aromatic = false;
        }
    }
}

fn clear_saturated_aromatic_carbon_atoms(mol: &mut Molecule) {
    let atoms_to_clear = mol
        .atoms()
        .filter_map(|(atom_id, atom)| {
            if atom.element.symbol() != "C" || !atom.aromatic {
                return None;
            }
            let bonds = mol
                .incident_bonds(atom_id)
                .ok()?
                .map(|(_, bond)| bond)
                .collect::<Vec<_>>();
            (!bonds.is_empty()
                && bonds
                    .iter()
                    .all(|bond| matches!(bond.order, BondOrder::Single))
                && bonds.iter().any(|bond| !bond.aromatic))
            .then_some(atom_id)
        })
        .collect::<BTreeSet<_>>();

    for atom_id in &atoms_to_clear {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = false;
        }
    }
    for bond in mol.bonds.iter_mut().flatten() {
        if atoms_to_clear.contains(&bond.a()) || atoms_to_clear.contains(&bond.b()) {
            bond.aromatic = false;
        }
    }
}

fn clear_aromatic_amidine_carbon_atoms(mol: &mut Molecule, rings: &[Ring]) {
    let amidine_carbons = mol
        .atoms()
        .filter_map(|(atom_id, atom)| {
            if atom.element.symbol() != "C" || !atom.aromatic {
                return None;
            }
            let bonds = mol
                .incident_bonds(atom_id)
                .ok()?
                .map(|(_, bond)| bond)
                .collect::<Vec<_>>();
            let has_imine_nitrogen = bonds.iter().any(|bond| {
                matches!(bond.order, BondOrder::Double)
                    && mol
                        .atom(bond.other_atom(atom_id))
                        .is_ok_and(|other| other.element.symbol() == "N")
            });
            let has_single_nitrogen = bonds.iter().any(|bond| {
                matches!(bond.order, BondOrder::Single)
                    && mol
                        .atom(bond.other_atom(atom_id))
                        .is_ok_and(|other| other.element.symbol() == "N")
            });
            let has_saturated_carbon_neighbor = bonds.iter().any(|bond| {
                let neighbor_id = bond.other_atom(atom_id);
                matches!(bond.order, BondOrder::Single)
                    && !bond.aromatic
                    && atoms_share_ring(rings, atom_id, neighbor_id)
                    && mol
                        .atom(neighbor_id)
                        .is_ok_and(|other| other.element.symbol() == "C" && !other.aromatic)
            });
            (has_imine_nitrogen && has_single_nitrogen && has_saturated_carbon_neighbor)
                .then_some(atom_id)
        })
        .collect::<BTreeSet<_>>();
    let mut atoms_to_clear = amidine_carbons.clone();
    for carbon_id in &amidine_carbons {
        atoms_to_clear.extend(
            mol.incident_bonds(*carbon_id)
                .ok()
                .into_iter()
                .flatten()
                .filter_map(|(_, bond)| {
                    let neighbor_id = bond.other_atom(*carbon_id);
                    mol.atom(neighbor_id)
                        .is_ok_and(|atom| atom.element.symbol() == "N")
                        .then_some(neighbor_id)
                }),
        );
    }

    for atom_id in &atoms_to_clear {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = false;
        }
    }
    for bond in mol.bonds.iter_mut().flatten() {
        if atoms_to_clear.contains(&bond.a()) || atoms_to_clear.contains(&bond.b()) {
            bond.aromatic = false;
        }
    }
}

fn atoms_share_ring(rings: &[Ring], left: AtomId, right: AtomId) -> bool {
    rings
        .iter()
        .any(|ring| ring.atoms.contains(&left) && ring.atoms.contains(&right))
}

fn aromatic_incident_bond_count(mol: &Molecule, atom_id: AtomId) -> usize {
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .filter(|(_, bond)| bond.aromatic)
        .count()
}

fn clear_ring_oxo_chalcogen_atoms(mol: &mut Molecule, rings: &[Ring]) {
    let atoms_to_clear = rings
        .iter()
        .flat_map(|ring| {
            ring.atoms
                .iter()
                .copied()
                .filter(|atom_id| is_ring_oxo_chalcogen(mol, ring, *atom_id))
        })
        .collect::<BTreeSet<_>>();

    for atom_id in &atoms_to_clear {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = false;
        }
    }
    for bond in mol.bonds.iter_mut().flatten() {
        if atoms_to_clear.contains(&bond.a()) || atoms_to_clear.contains(&bond.b()) {
            bond.aromatic = false;
        }
    }
}

fn clear_terminal_chalcogen_oxo_ring_atoms(mol: &mut Molecule, rings: &[Ring]) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        if ring.atoms.len() != 5
            || !ring_contains_element(mol, ring, "N")
            || !ring_has_terminal_chalcogen_exocyclic_pi_bond(mol, ring)
        {
            continue;
        }
        for atom_id in &ring.atoms {
            let retained_by_other_ring = rings.iter().enumerate().any(|(other_index, other)| {
                other_index != index
                    && other.atoms.contains(atom_id)
                    && !ring_has_terminal_chalcogen_exocyclic_pi_bond(mol, other)
            });
            if !retained_by_other_ring {
                atoms_to_clear.insert(*atom_id);
            }
        }
    }

    for atom_id in &atoms_to_clear {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = false;
        }
    }
    for bond in mol.bonds.iter_mut().flatten() {
        if atoms_to_clear.contains(&bond.a()) || atoms_to_clear.contains(&bond.b()) {
            bond.aromatic = false;
        }
    }
}

fn clear_fused_lactam_bridge_ring_atoms(mol: &mut Molecule, rings: &[Ring]) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        if ring.atoms.len() < 6
            || !ring_contains_element(mol, ring, "N")
            || !ring_contains_element(mol, ring, "O")
            || ring_terminal_exocyclic_pi_bond_count(mol, ring) == 0
        {
            continue;
        }
        for atom_id in &ring.atoms {
            let retained_by_other_ring = rings.iter().enumerate().any(|(other_index, other)| {
                other_index != index
                    && other.atoms.contains(atom_id)
                    && ring_terminal_exocyclic_pi_bond_count(mol, other) == 0
            });
            if !retained_by_other_ring {
                atoms_to_clear.insert(*atom_id);
            }
        }
    }

    for atom_id in &atoms_to_clear {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = false;
        }
    }
    for bond in mol.bonds.iter_mut().flatten() {
        if atoms_to_clear.contains(&bond.a()) || atoms_to_clear.contains(&bond.b()) {
            bond.aromatic = false;
        }
    }
}

fn clear_imide_carbonyl_ring_atoms(mol: &mut Molecule, rings: &[Ring]) {
    let atoms_to_clear = rings
        .iter()
        .enumerate()
        .filter(|(index, ring)| {
            ring_contains_element(mol, ring, "N")
                && (ring_has_or_is_fused_to_cationic_nitrogen(mol, rings, *index)
                    || ring.atoms.len() == 5 && ring_has_imide_nitrogen(mol, ring))
                && ring_terminal_exocyclic_pi_bond_count(mol, ring) >= 2
        })
        .flat_map(|(_, ring)| {
            ring.atoms.iter().copied().filter(|atom_id| {
                atom_has_terminal_exocyclic_pi_bond(mol, ring, *atom_id)
                    || is_saturated_ring_carbon(mol, ring, *atom_id)
            })
        })
        .collect::<BTreeSet<_>>();

    for atom_id in &atoms_to_clear {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = false;
        }
    }
    for bond in mol.bonds.iter_mut().flatten() {
        if atoms_to_clear.contains(&bond.a()) || atoms_to_clear.contains(&bond.b()) {
            bond.aromatic = false;
        }
    }
}

fn is_saturated_ring_carbon(mol: &Molecule, ring: &Ring, atom_id: AtomId) -> bool {
    mol.atom(atom_id)
        .is_ok_and(|atom| atom.element.symbol() == "C")
        && !ring_atom_has_pi_bond(mol, ring, atom_id)
        && !atom_has_exocyclic_pi_bond(mol, ring, atom_id)
}

fn ring_has_imide_nitrogen(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().any(|atom_id| {
        let Ok(atom) = mol.atom(*atom_id) else {
            return false;
        };
        atom.element.symbol() == "N"
            && atom.formal_charge == 0
            && !ring_atom_has_pi_bond(mol, ring, *atom_id)
            && !atom_has_exocyclic_pi_bond(mol, ring, *atom_id)
            && ring_carbonyl_neighbor_count(mol, ring, *atom_id) >= 2
    })
}

fn ring_carbonyl_neighbor_count(mol: &Molecule, ring: &Ring, atom_id: AtomId) -> usize {
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .filter(|(_, bond)| ring.atoms.contains(&bond.other_atom(atom_id)))
        .filter(|(_, bond)| {
            atom_has_terminal_exocyclic_pi_bond(mol, ring, bond.other_atom(atom_id))
        })
        .count()
}

fn ring_has_or_is_fused_to_cationic_nitrogen(mol: &Molecule, rings: &[Ring], index: usize) -> bool {
    ring_has_cationic_nitrogen(mol, &rings[index])
        || rings.iter().enumerate().any(|(other_index, other)| {
            other_index != index
                && rings_share_bond(&rings[index], other)
                && ring_has_cationic_nitrogen(mol, other)
        })
}

fn ring_has_cationic_nitrogen(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().any(|atom_id| {
        mol.atom(*atom_id)
            .is_ok_and(|atom| atom.element.symbol() == "N" && atom.formal_charge > 0)
    })
}

fn clear_fused_lactone_bridge_ring_atoms(mol: &mut Molecule, rings: &[Ring]) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        let fused = rings
            .iter()
            .enumerate()
            .any(|(other_index, other)| other_index != index && rings_share_bond(ring, other));
        if !fused
            || ring.atoms.len() < 6
            || !ring_has_chalcogen_donor(mol, ring)
            || ring_has_conjugated_atom_path(mol, ring)
            || !ring
                .atoms
                .iter()
                .any(|atom_id| is_chalcogen_bridge_without_pi(mol, ring, *atom_id))
            || ring_terminal_exocyclic_pi_bond_count(mol, ring) == 0
        {
            continue;
        }
        for atom_id in &ring.atoms {
            let retained_by_other_ring = rings.iter().enumerate().any(|(other_index, other)| {
                other_index != index
                    && other.atoms.contains(atom_id)
                    && ring_terminal_exocyclic_pi_bond_count(mol, other) == 0
                    && other
                        .atoms
                        .iter()
                        .any(|other_atom| mol.atom(*other_atom).is_ok_and(|atom| atom.aromatic))
            });
            if !retained_by_other_ring {
                atoms_to_clear.insert(*atom_id);
            }
        }
    }

    for atom_id in &atoms_to_clear {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = false;
        }
    }
    for bond in mol.bonds.iter_mut().flatten() {
        if atoms_to_clear.contains(&bond.a()) || atoms_to_clear.contains(&bond.b()) {
            bond.aromatic = false;
        }
    }
}

fn is_chalcogen_bridge_without_pi(mol: &Molecule, ring: &Ring, atom_id: AtomId) -> bool {
    mol.atom(atom_id)
        .is_ok_and(|atom| matches!(atom.element.symbol(), "O" | "S" | "Se" | "Te"))
        && !ring_atom_has_pi_bond(mol, ring, atom_id)
        && !atom_has_exocyclic_pi_bond(mol, ring, atom_id)
}

fn clear_saturated_tertiary_amine_ring_atoms(mol: &mut Molecule, rings: &[Ring]) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        if !ring_has_saturated_tertiary_amine_without_donor_chalcogen(mol, ring) {
            continue;
        }
        for atom_id in &ring.atoms {
            let retained_by_other_ring = rings.iter().enumerate().any(|(other_index, other)| {
                other_index != index
                    && other.atoms.contains(atom_id)
                    && !ring_has_saturated_tertiary_amine_without_donor_chalcogen(mol, other)
                    && other
                        .atoms
                        .iter()
                        .any(|other_atom| mol.atom(*other_atom).is_ok_and(|atom| atom.aromatic))
            });
            if !retained_by_other_ring {
                atoms_to_clear.insert(*atom_id);
            }
        }
    }

    for atom_id in &atoms_to_clear {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = false;
        }
    }
    for bond in mol.bonds.iter_mut().flatten() {
        if atoms_to_clear.contains(&bond.a()) || atoms_to_clear.contains(&bond.b()) {
            bond.aromatic = false;
        }
    }
}

fn clear_saturated_fused_carbon_ring_atoms(mol: &mut Molecule, rings: &[Ring]) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        let fused = rings
            .iter()
            .enumerate()
            .any(|(other_index, other)| other_index != index && rings_share_bond(ring, other));
        if !fused
            || !fused_component_is_all_carbon(mol, ring)
            || !ring_has_saturated_carbon_atom(mol, ring)
        {
            continue;
        }
        for atom_id in &ring.atoms {
            if !atom_is_retained_by_other_aromatic_ring(mol, rings, index, *atom_id) {
                atoms_to_clear.insert(*atom_id);
            }
        }
    }

    for atom_id in &atoms_to_clear {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = false;
        }
    }
    for bond in mol.bonds.iter_mut().flatten() {
        if atoms_to_clear.contains(&bond.a()) || atoms_to_clear.contains(&bond.b()) {
            bond.aromatic = false;
        }
    }
}

fn clear_fused_carbonyl_bridge_atoms(mol: &mut Molecule, rings: &[Ring]) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        let fused = rings
            .iter()
            .enumerate()
            .any(|(other_index, other)| other_index != index && rings_share_bond(ring, other));
        if !fused || !fused_component_is_all_carbon(mol, ring) {
            continue;
        }
        if ring.atoms.len() > 4 && ring_terminal_exocyclic_pi_bond_count(mol, ring) >= 2 {
            atoms_to_clear.extend(ring.atoms.iter().copied().filter(|atom_id| {
                !atom_is_retained_by_other_aromatic_ring(mol, rings, index, *atom_id)
            }));
            continue;
        }
        if ring.atoms.len() == 5 || ring_has_saturated_carbon_atom(mol, ring) {
            atoms_to_clear.extend(
                ring.atoms
                    .iter()
                    .copied()
                    .filter(|atom_id| atom_has_terminal_exocyclic_pi_bond(mol, ring, *atom_id)),
            );
        }
    }

    for atom_id in &atoms_to_clear {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = false;
        }
    }
    for bond in mol.bonds.iter_mut().flatten() {
        if atoms_to_clear.contains(&bond.a()) || atoms_to_clear.contains(&bond.b()) {
            bond.aromatic = false;
        }
    }
}

fn atom_is_retained_by_other_aromatic_ring(
    mol: &Molecule,
    rings: &[Ring],
    ring_index: usize,
    atom_id: AtomId,
) -> bool {
    rings.iter().enumerate().any(|(other_index, other)| {
        other_index != ring_index
            && other.atoms.contains(&atom_id)
            && other.bonds.iter().any(|bond_id| {
                mol.bond(*bond_id)
                    .is_ok_and(|bond| bond.aromatic && (bond.a() == atom_id || bond.b() == atom_id))
            })
    })
}

fn ring_has_saturated_carbon_atom(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().any(|atom_id| {
        mol.atom(*atom_id)
            .is_ok_and(|atom| atom.element.symbol() == "C")
            && !ring_atom_has_pi_bond(mol, ring, *atom_id)
            && !atom_has_exocyclic_pi_bond(mol, ring, *atom_id)
    })
}

fn perceive_fused_single_exocyclic_carbon_rings(
    mol: &mut Molecule,
    rings: &[Ring],
    protected_non_aromatic_bonds: &BTreeSet<BondId>,
) {
    let selected = rings
        .iter()
        .enumerate()
        .filter(|(index, ring)| {
            let fused = rings
                .iter()
                .enumerate()
                .any(|(other_index, other)| other_index != *index && rings_share_bond(ring, other));
            fused
                && ring_has_conjugated_atom_path(mol, ring)
                && ((ring.atoms.len() == 6
                    && fused_component_is_all_carbon(mol, ring)
                    && ring_exocyclic_pi_bond_count(mol, ring) > 0
                    && ring_terminal_exocyclic_pi_bond_count(mol, ring) <= 1
                    && ring_pi_bond_count(mol, ring) >= 1)
                    || (ring.atoms.len() == 4
                        && fused_component_is_all_carbon(mol, ring)
                        && ring_terminal_exocyclic_pi_bond_count(mol, ring) >= 2)
                    || (ring.atoms.len() == 5
                        && ring_contains_element(mol, ring, "N")
                        && ring_has_chalcogen_donor(mol, ring)
                        && ring_has_carbon_hetero_exocyclic_pi_bond(mol, ring)
                        && ring_exocyclic_pi_bond_count(mol, ring) > 0))
        })
        .map(|(_, ring)| ring.clone())
        .collect::<Vec<_>>();

    for ring in selected {
        for atom_id in &ring.atoms {
            if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
                atom.aromatic = true;
            }
        }
        for bond_id in &ring.bonds {
            if let Some(bond) = mol.bonds[bond_id.index()].as_mut() {
                bond.aromatic = !protected_non_aromatic_bonds.contains(bond_id);
            }
        }
    }
}

fn perceive_fused_aromatic_components(
    mol: &mut Molecule,
    rings: &[Ring],
    protected_non_aromatic_bonds: &BTreeSet<BondId>,
) -> std::result::Result<(), AromaticityError> {
    let candidates = rings
        .iter()
        .enumerate()
        .filter_map(|(index, ring)| aromatic_fused_candidate(mol, ring).then_some(index))
        .collect::<Vec<_>>();
    let mut components = (0..candidates.len()).collect::<Vec<_>>();
    for left in 0..candidates.len() {
        for right in (left + 1)..candidates.len() {
            if rings_share_bond(&rings[candidates[left]], &rings[candidates[right]]) {
                union_components(&mut components, left, right);
            }
        }
    }

    let mut component_rings = BTreeMap::<usize, Vec<usize>>::new();
    for (component_index, ring_index) in candidates.iter().copied().enumerate() {
        let root = find_component(&mut components, component_index);
        component_rings.entry(root).or_default().push(ring_index);
    }

    for indexes in component_rings.values() {
        if indexes.len() < 2 {
            continue;
        }
        let component = fused_component_ring(rings, indexes);
        if component.atoms.len() > 48 {
            continue;
        }
        let electrons = aromatic_fused_component_pi_electrons(mol, &component)?;
        let all_carbon_component = fused_component_is_all_carbon(mol, &component);
        let component_has_exocyclic_pi = ring_exocyclic_pi_bond_count(mol, &component) > 0;
        if electrons >= 6
            && ((electrons - 2) % 4 == 0 || all_carbon_component && !component_has_exocyclic_pi)
        {
            for atom_id in &component.atoms {
                if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
                    atom.aromatic = true;
                }
            }
            for bond_id in &component.bonds {
                if let Some(bond) = mol.bonds[bond_id.index()].as_mut() {
                    bond.aromatic = !protected_non_aromatic_bonds.contains(bond_id);
                }
            }
        } else if all_carbon_component && component_has_exocyclic_pi {
            let aromatic_atoms =
                atoms_in_limited_terminal_exocyclic_pi_rings(mol, rings, indexes, 1);
            if aromatic_atoms.len() >= 6 {
                for atom_id in &aromatic_atoms {
                    if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
                        atom.aromatic = true;
                    }
                }
                for bond_id in &component.bonds {
                    if let Some(bond) = mol.bonds[bond_id.index()].as_mut() {
                        bond.aromatic = aromatic_atoms.contains(&bond.a())
                            && aromatic_atoms.contains(&bond.b())
                            && !protected_non_aromatic_bonds.contains(bond_id);
                    }
                }
            }
        } else if fused_component_is_carbon_nitrogen(mol, &component) {
            let terminal_atoms_retained =
                terminal_exocyclic_atoms_in_nitrogen_rings(mol, rings, indexes);
            let atoms_retained_by_ring_context =
                atoms_in_nitrogen_or_terminal_pi_free_rings(mol, rings, indexes);
            let aromatic_atoms = component
                .atoms
                .iter()
                .copied()
                .filter(|atom_id| {
                    atoms_retained_by_ring_context.contains(atom_id)
                        && (!atom_has_terminal_exocyclic_pi_bond(mol, &component, *atom_id)
                            || terminal_atoms_retained.contains(atom_id))
                })
                .collect::<BTreeSet<_>>();
            if aromatic_atoms.len() >= 6 {
                for atom_id in &aromatic_atoms {
                    if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
                        atom.aromatic = true;
                    }
                }
                for bond_id in &component.bonds {
                    if let Some(bond) = mol.bonds[bond_id.index()].as_mut() {
                        bond.aromatic = aromatic_atoms.contains(&bond.a())
                            && aromatic_atoms.contains(&bond.b())
                            && !protected_non_aromatic_bonds.contains(bond_id);
                    }
                }
            }
        }
    }
    Ok(())
}

fn fused_component_is_all_carbon(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().all(|atom_id| {
        mol.atom(*atom_id)
            .map(|atom| atom.element.symbol() == "C")
            .unwrap_or(false)
    })
}

fn fused_component_is_carbon_nitrogen(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().all(|atom_id| {
        mol.atom(*atom_id)
            .map(|atom| matches!(atom.element.symbol(), "C" | "N"))
            .unwrap_or(false)
    })
}

fn terminal_exocyclic_atoms_in_nitrogen_rings(
    mol: &Molecule,
    rings: &[Ring],
    indexes: &[usize],
) -> BTreeSet<AtomId> {
    let mut atoms = BTreeSet::new();
    for index in indexes {
        let ring = &rings[*index];
        if !ring_contains_element(mol, ring, "N") {
            continue;
        }
        for atom_id in &ring.atoms {
            if atom_has_terminal_exocyclic_pi_bond(mol, ring, *atom_id) {
                atoms.insert(*atom_id);
            }
        }
    }
    atoms
}

fn atoms_in_limited_terminal_exocyclic_pi_rings(
    mol: &Molecule,
    rings: &[Ring],
    indexes: &[usize],
    max_terminal_exocyclic_pi: usize,
) -> BTreeSet<AtomId> {
    let mut atoms = BTreeSet::new();
    for index in indexes {
        let ring = &rings[*index];
        if ring_terminal_exocyclic_pi_bond_count(mol, ring) <= max_terminal_exocyclic_pi {
            atoms.extend(ring.atoms.iter().copied());
        }
    }
    atoms
}

fn atoms_in_nitrogen_or_terminal_pi_free_rings(
    mol: &Molecule,
    rings: &[Ring],
    indexes: &[usize],
) -> BTreeSet<AtomId> {
    let has_exocyclic_pi_ring = indexes
        .iter()
        .any(|index| ring_exocyclic_pi_bond_count(mol, &rings[*index]) > 0);
    if !has_exocyclic_pi_ring {
        return indexes
            .iter()
            .filter(|index| {
                !ring_has_saturated_tertiary_amine_without_donor_chalcogen(mol, &rings[**index])
            })
            .flat_map(|index| rings[*index].atoms.iter().copied())
            .collect();
    }

    let mut atoms = BTreeSet::new();
    for index in indexes {
        let ring = &rings[*index];
        if ring_has_saturated_tertiary_amine_without_donor_chalcogen(mol, ring) {
            continue;
        }
        let exocyclic_pi_count = ring_exocyclic_pi_bond_count(mol, ring);
        let contains_nitrogen = ring_contains_element(mol, ring, "N");
        if exocyclic_pi_count == 0
            || !contains_nitrogen && exocyclic_pi_count <= 1
            || contains_nitrogen && ring_hetero_donor_count(mol, ring) >= 2
            || ring.atoms.len() >= 6
                && contains_nitrogen
                && !ring_has_saturated_carbon_atom(mol, ring)
        {
            atoms.extend(ring.atoms.iter().copied());
        }
    }
    atoms
}

fn ring_has_saturated_tertiary_amine_without_donor_chalcogen(mol: &Molecule, ring: &Ring) -> bool {
    !ring_has_saturated_chalcogen_donor(mol, ring)
        && !ring_has_conjugated_atom_path(mol, ring)
        && ring_active_hetero_donor_count(mol, ring) == 1
        && ring
            .atoms
            .iter()
            .any(|atom_id| is_saturated_tertiary_amine(mol, ring, *atom_id))
}

fn ring_has_saturated_chalcogen_donor(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().any(|atom_id| {
        mol.atom(*atom_id)
            .is_ok_and(|atom| matches!(atom.element.symbol(), "O" | "S" | "Se" | "Te"))
            && !ring_atom_has_pi_bond(mol, ring, *atom_id)
            && !atom_has_exocyclic_pi_bond(mol, ring, *atom_id)
    })
}

fn ring_active_hetero_donor_count(mol: &Molecule, ring: &Ring) -> usize {
    ring.atoms
        .iter()
        .filter(|atom_id| {
            let Ok(atom) = mol.atom(**atom_id) else {
                return false;
            };
            match atom.element.symbol() {
                "N" | "P" => true,
                "O" | "S" | "Se" | "Te" => !atom_has_exocyclic_pi_bond(mol, ring, **atom_id),
                _ => false,
            }
        })
        .count()
}

fn aromatic_fused_component_pi_electrons(
    mol: &Molecule,
    ring: &Ring,
) -> std::result::Result<u8, AromaticityError> {
    let mut electrons = 0u8;
    let mut has_pi_bond = false;
    for bond_id in &ring.bonds {
        let bond = mol.bond(*bond_id).expect("ring bond should be live");
        if matches!(bond.order, BondOrder::Double | BondOrder::Aromatic) {
            has_pi_bond = true;
            electrons += 2;
        }
    }

    for atom_id in &ring.atoms {
        let atom = mol.atom(*atom_id).expect("ring atom should be live");
        match atom.element.symbol() {
            "C" => {}
            "N" => {
                if !(ring_atom_has_pi_bond(mol, ring, *atom_id)
                    || atom.formal_charge > 0 && atom_has_exocyclic_pi_bond(mol, ring, *atom_id))
                {
                    electrons += 2;
                }
            }
            "O" | "S" | "Se" | "Te" | "P" => {
                if !ring_atom_has_pi_bond(mol, ring, *atom_id) {
                    electrons += 2;
                }
            }
            _ => {
                if ring_has_aromatic_order(mol, ring) {
                    return Err(AromaticityError::UnsupportedElement(*atom_id));
                }
                return Ok(0);
            }
        }
    }

    if ring.atoms.len() == 6
        && fused_component_is_all_carbon(mol, ring)
        && ring_exocyclic_pi_bond_count(mol, ring) == 1
        && ring_pi_bond_count(mol, ring) >= 1
    {
        return Ok(6);
    }

    if has_pi_bond {
        Ok(electrons)
    } else {
        Ok(0)
    }
}

fn ring_pi_bond_count(mol: &Molecule, ring: &Ring) -> usize {
    ring.bonds
        .iter()
        .filter(|bond_id| {
            mol.bond(**bond_id)
                .map(|bond| matches!(bond.order, BondOrder::Double | BondOrder::Aromatic))
                .unwrap_or(false)
        })
        .count()
}

fn aromatic_fused_candidate(mol: &Molecule, ring: &Ring) -> bool {
    if ring.atoms.len() == 7 && ring_has_chalcogen_donor(mol, ring) {
        return false;
    }
    if ring.atoms.len() > 7 {
        return ring.atoms.len() <= 18
            && fused_component_is_carbon_nitrogen(mol, ring)
            && ring_has_anionic_nitrogen(mol, ring)
            && ring_has_conjugated_atom_path(mol, ring)
            && ring_pi_bond_count(mol, ring) >= 4;
    }
    let pi_bonds = ring_pi_bond_count(mol, ring);
    ring_has_conjugated_atom_path(mol, ring)
        && (pi_bonds >= 2
            || ring_hetero_donor_count(mol, ring) > 0
                && !ring_has_low_unsaturation_chalcogen_bridge_for_fused(mol, ring)
            || ring.atoms.len() == 5 && pi_bonds >= 1 && ring_contains_element(mol, ring, "N")
            || ring.atoms.len() == 6 && pi_bonds >= 1 && fused_component_is_all_carbon(mol, ring))
}

fn ring_has_conjugated_atom_path(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().all(|atom_id| {
        let Ok(atom) = mol.atom(*atom_id) else {
            return false;
        };
        match atom.element.symbol() {
            "C" => {
                ring_atom_has_pi_bond(mol, ring, *atom_id)
                    || atom_has_exocyclic_pi_bond(mol, ring, *atom_id)
            }
            "N" | "O" | "S" | "Se" | "Te" | "P" => true,
            _ => false,
        }
    })
}

fn is_saturated_tertiary_amine(mol: &Molecule, ring: &Ring, atom_id: AtomId) -> bool {
    let Ok(atom) = mol.atom(atom_id) else {
        return false;
    };
    atom.element.symbol() == "N"
        && atom.formal_charge == 0
        && atom.explicit_hydrogens == 0
        && !ring_atom_has_pi_bond(mol, ring, atom_id)
        && !atom_has_exocyclic_pi_bond(mol, ring, atom_id)
        && mol
            .incident_bonds(atom_id)
            .map(|bonds| {
                let mut degree = 0usize;
                for (_, bond) in bonds {
                    degree += 1;
                    let other = bond.other_atom(atom_id);
                    if !ring.atoms.contains(&other)
                        && !is_saturated_tertiary_amine_substituent(mol, ring, other)
                    {
                        return false;
                    }
                    if !matches!(bond.order, BondOrder::Single) {
                        return false;
                    }
                }
                degree >= 3
            })
            .unwrap_or(false)
}

fn is_saturated_tertiary_amine_substituent(mol: &Molecule, ring: &Ring, atom_id: AtomId) -> bool {
    let Ok(atom) = mol.atom(atom_id) else {
        return false;
    };
    atom.element.symbol() == "C"
        || matches!(atom.element.symbol(), "S" | "Se" | "Te")
            && atom_has_terminal_exocyclic_pi_bond(mol, ring, atom_id)
}

fn is_saturated_chalcogen_bridge(mol: &Molecule, ring: &Ring, atom_id: AtomId) -> bool {
    let Ok(atom) = mol.atom(atom_id) else {
        return false;
    };
    matches!(atom.element.symbol(), "S" | "Se" | "Te")
        && atom.formal_charge == 0
        && atom.explicit_hydrogens == 0
        && !ring_atom_has_pi_bond(mol, ring, atom_id)
        && !atom_has_exocyclic_pi_bond(mol, ring, atom_id)
        && mol
            .incident_bonds(atom_id)
            .map(|bonds| {
                let mut degree = 0usize;
                for (_, bond) in bonds {
                    degree += 1;
                    if !matches!(bond.order, BondOrder::Single) {
                        return false;
                    }
                }
                degree == 2
            })
            .unwrap_or(false)
}

fn is_ring_oxo_chalcogen(mol: &Molecule, ring: &Ring, atom_id: AtomId) -> bool {
    let Ok(atom) = mol.atom(atom_id) else {
        return false;
    };
    matches!(atom.element.symbol(), "S" | "Se" | "Te")
        && atom.formal_charge == 0
        && ring.atoms.contains(&atom_id)
        && atom_has_terminal_exocyclic_pi_bond(mol, ring, atom_id)
}

fn ring_hetero_donor_count(mol: &Molecule, ring: &Ring) -> usize {
    ring.atoms
        .iter()
        .filter(|atom_id| {
            mol.atom(**atom_id)
                .map(|atom| matches!(atom.element.symbol(), "N" | "O" | "S" | "Se" | "Te" | "P"))
                .unwrap_or(false)
        })
        .count()
}

fn ring_has_anionic_nitrogen(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().any(|atom_id| {
        mol.atom(*atom_id)
            .map(|atom| atom.element.symbol() == "N" && atom.formal_charge < 0)
            .unwrap_or(false)
    })
}

fn ring_has_nitrogen_lone_pair_donor(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().any(|atom_id| {
        mol.atom(*atom_id)
            .map(|atom| {
                atom.element.symbol() == "N"
                    && (atom.explicit_hydrogens > 0
                        || atom.formal_charge < 0
                        || atom.formal_charge == 0 && !ring_atom_has_pi_bond(mol, ring, *atom_id)
                        || mol.neighbors(*atom_id).is_ok_and(|mut neighbors| {
                            neighbors.any(|neighbor| {
                                mol.atom(neighbor).is_ok_and(|neighbor_atom| {
                                    neighbor_atom.element.symbol() == "H"
                                })
                            })
                        }))
            })
            .unwrap_or(false)
    })
}

fn ring_has_chalcogen_donor(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().any(|atom_id| {
        mol.atom(*atom_id)
            .map(|atom| matches!(atom.element.symbol(), "O" | "S" | "Se" | "Te"))
            .unwrap_or(false)
    })
}

fn ring_has_low_unsaturation_chalcogen_bridge(mol: &Molecule, ring: &Ring) -> bool {
    ring_pi_bond_count(mol, ring) + ring_terminal_exocyclic_pi_bond_count(mol, ring) < 2
        && ring_hetero_donor_count(mol, ring) > 1
        && ring_has_chalcogen_donor(mol, ring)
}

fn ring_has_low_unsaturation_chalcogen_bridge_for_fused(mol: &Molecule, ring: &Ring) -> bool {
    ring_pi_bond_count(mol, ring) < 2
        && ring_hetero_donor_count(mol, ring) > 1
        && ring_has_chalcogen_donor(mol, ring)
}

fn ring_exocyclic_pi_bond_count(mol: &Molecule, ring: &Ring) -> usize {
    ring.atoms
        .iter()
        .filter(|atom_id| atom_has_exocyclic_pi_bond(mol, ring, **atom_id))
        .count()
}

fn ring_has_carbon_hetero_exocyclic_pi_bond(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().any(|atom_id| {
        mol.atom(*atom_id)
            .map(|atom| atom.element.symbol() == "C")
            .unwrap_or(false)
            && mol
                .incident_bonds(*atom_id)
                .ok()
                .into_iter()
                .flatten()
                .any(|(bond_id, bond)| {
                    if ring.bonds.contains(&bond_id)
                        || ring.atoms.contains(&bond.other_atom(*atom_id))
                        || !matches!(bond.order, BondOrder::Double | BondOrder::Aromatic)
                    {
                        return false;
                    }
                    mol.atom(bond.other_atom(*atom_id))
                        .map(|atom| atom.element.symbol() != "C")
                        .unwrap_or(false)
                })
    })
}

fn ring_has_terminal_chalcogen_exocyclic_pi_bond(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().any(|atom_id| {
        mol.atom(*atom_id)
            .map(|atom| matches!(atom.element.symbol(), "S" | "Se" | "Te"))
            .unwrap_or(false)
            && atom_has_terminal_exocyclic_pi_bond(mol, ring, *atom_id)
    })
}

fn ring_terminal_exocyclic_pi_bond_count(mol: &Molecule, ring: &Ring) -> usize {
    ring.atoms
        .iter()
        .filter(|atom_id| atom_has_terminal_exocyclic_pi_bond(mol, ring, **atom_id))
        .count()
}

fn rings_share_bond(left: &Ring, right: &Ring) -> bool {
    left.bonds.iter().any(|bond| right.bonds.contains(bond))
}

fn fused_component_ring(rings: &[Ring], indexes: &[usize]) -> Ring {
    let mut atoms = BTreeSet::new();
    let mut bonds = BTreeSet::new();
    for left in indexes {
        atoms.extend(rings[*left].atoms.iter().copied());
        bonds.extend(rings[*left].bonds.iter().copied());
    }
    Ring {
        atoms: atoms.into_iter().collect(),
        bonds: bonds.into_iter().collect(),
    }
}

fn aromatic_ring_pi_electrons(
    mol: &Molecule,
    ring: &Ring,
) -> std::result::Result<u8, AromaticityError> {
    if ring_has_aromatic_order(mol, ring) {
        return aromatic_order_ring_pi_electrons(mol, ring);
    }
    if ring.atoms.len() == 7 && ring_has_chalcogen_donor(mol, ring) {
        return Ok(0);
    }
    if ring.atoms.len() == 5
        && ring_pi_bond_count(mol, ring) < 2
        && !ring_contains_element(mol, ring, "N")
    {
        return Ok(0);
    }
    if ring.atoms.len() == 5
        && ring_contains_element(mol, ring, "N")
        && ring_exocyclic_pi_bond_count(mol, ring) > 0
        && ((ring_hetero_donor_count(mol, ring) < 2 && !ring_has_chalcogen_donor(mol, ring))
            || (ring_has_chalcogen_donor(mol, ring)
                && !ring_has_carbon_hetero_exocyclic_pi_bond(mol, ring)))
    {
        return Ok(0);
    }
    if ring.atoms.len() == 5
        && ring_contains_element(mol, ring, "N")
        && ring_has_terminal_chalcogen_exocyclic_pi_bond(mol, ring)
    {
        return Ok(0);
    }
    if ring.atoms.len() > 5 && ring_has_low_unsaturation_chalcogen_bridge(mol, ring) {
        return Ok(0);
    }
    if ring.atoms.len() > 5
        && ring_exocyclic_pi_bond_count(mol, ring) == 0
        && ring_has_chalcogen_donor(mol, ring)
        && ring_contains_element(mol, ring, "N")
        && ring_hetero_donor_count(mol, ring) > 1
    {
        return Ok(0);
    }

    let mut electrons = 0u8;
    let mut has_pi_bond = false;
    for bond_id in &ring.bonds {
        let bond = mol.bond(*bond_id).expect("ring bond should be live");
        if matches!(bond.order, BondOrder::Double | BondOrder::Aromatic) {
            has_pi_bond = true;
            electrons += 2;
        }
    }

    for atom_id in &ring.atoms {
        let atom = mol.atom(*atom_id).expect("ring atom should be live");
        match atom.element.symbol() {
            "C" => {
                if !ring_atom_has_pi_bond(mol, ring, *atom_id)
                    && !atom_has_exocyclic_pi_bond(mol, ring, *atom_id)
                {
                    if atom.formal_charge < 0 {
                        electrons += 2;
                    } else {
                        return Ok(0);
                    }
                }
            }
            "N" => {
                if !(ring_atom_has_pi_bond(mol, ring, *atom_id)
                    || atom.formal_charge > 0 && atom_has_exocyclic_pi_bond(mol, ring, *atom_id))
                {
                    electrons += 2;
                }
            }
            "O" | "S" | "Se" | "Te" | "P" => {
                if !ring_atom_has_pi_bond(mol, ring, *atom_id) {
                    electrons += 2;
                }
            }
            _ => {
                if ring_has_aromatic_order(mol, ring) {
                    return Err(AromaticityError::UnsupportedElement(*atom_id));
                }
                return Ok(0);
            }
        }
    }

    if has_pi_bond {
        Ok(electrons)
    } else {
        Ok(0)
    }
}

fn aromatic_order_ring_pi_electrons(
    mol: &Molecule,
    ring: &Ring,
) -> std::result::Result<u8, AromaticityError> {
    if ring.atoms.len() == 6
        && ring_has_chalcogen_donor(mol, ring)
        && ring_has_carbon_hetero_exocyclic_pi_bond(mol, ring)
        && ring_terminal_exocyclic_pi_bond_count(mol, ring) == 1
        && !ring_has_external_aromatic_bond(mol, ring)
    {
        return Ok(6);
    }

    let mut min_electrons = 0u8;
    let mut max_electrons = 0u8;
    for atom_id in &ring.atoms {
        let atom = mol.atom(*atom_id).expect("ring atom should be live");
        let (min_contribution, max_contribution) = match atom.element.symbol() {
            "B" | "C" => (1, 1),
            "N" => {
                if atom.explicit_hydrogens > 0
                    || atom.formal_charge == 0
                        && aromatic_order_nitrogen_is_pyrrole_like(mol, *atom_id)
                {
                    (1, 2)
                } else {
                    (1, 1)
                }
            }
            "O" | "S" | "Se" | "Te" | "P" => (2, 2),
            _ => return Err(AromaticityError::UnsupportedElement(*atom_id)),
        };
        min_electrons += min_contribution;
        max_electrons += max_contribution;
    }
    if let Some(electrons) = huckel_electron_count_in_range(min_electrons, max_electrons) {
        Ok(electrons)
    } else {
        Ok(max_electrons)
    }
}

fn huckel_electron_count_in_range(min_electrons: u8, max_electrons: u8) -> Option<u8> {
    if max_electrons == 2 {
        return Some(2);
    }
    if max_electrons < 6 {
        return None;
    }
    (min_electrons..=max_electrons).find(|electrons| (electrons - 2) % 4 == 0)
}

fn aromatic_order_nitrogen_is_pyrrole_like(mol: &Molecule, atom_id: AtomId) -> bool {
    mol.incident_bonds(atom_id)
        .map(|bonds| bonds.count() >= 3)
        .unwrap_or(false)
}

fn ring_contains_element(mol: &Molecule, ring: &Ring, symbol: &str) -> bool {
    ring.atoms.iter().any(|atom_id| {
        mol.atom(*atom_id)
            .map(|atom| atom.element.symbol() == symbol)
            .unwrap_or(false)
    })
}

fn ring_has_aromatic_order(mol: &Molecule, ring: &Ring) -> bool {
    ring.bonds.iter().any(|bond_id| {
        mol.bond(*bond_id)
            .map(|bond| bond.order == BondOrder::Aromatic)
            .unwrap_or(false)
    })
}

fn ring_has_external_aromatic_bond(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().any(|atom_id| {
        mol.incident_bonds(*atom_id)
            .ok()
            .into_iter()
            .flatten()
            .any(|(bond_id, bond)| {
                !ring.bonds.contains(&bond_id)
                    && matches!(bond.order, BondOrder::Aromatic)
                    && mol
                        .atom(bond.other_atom(*atom_id))
                        .is_ok_and(|atom| atom.aromatic)
            })
    })
}

fn ring_atom_has_pi_bond(mol: &Molecule, ring: &Ring, atom_id: AtomId) -> bool {
    ring.bonds.iter().any(|bond_id| {
        mol.bond(*bond_id)
            .map(|bond| {
                (bond.a == atom_id || bond.b == atom_id)
                    && matches!(bond.order, BondOrder::Double | BondOrder::Aromatic)
            })
            .unwrap_or(false)
    })
}

fn atom_has_exocyclic_pi_bond(mol: &Molecule, ring: &Ring, atom_id: AtomId) -> bool {
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .any(|(bond_id, bond)| {
            !ring.bonds.contains(&bond_id)
                && !ring.atoms.contains(&bond.other_atom(atom_id))
                && matches!(bond.order, BondOrder::Double | BondOrder::Aromatic)
        })
}

fn atom_has_terminal_exocyclic_pi_bond(mol: &Molecule, ring: &Ring, atom_id: AtomId) -> bool {
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .any(|(bond_id, bond)| {
            if ring.bonds.contains(&bond_id) || !matches!(bond.order, BondOrder::Double) {
                return false;
            }
            let other = bond.other_atom(atom_id);
            if ring.atoms.contains(&other) {
                return false;
            }
            mol.incident_bonds(bond.other_atom(atom_id))
                .map(|bonds| bonds.count() == 1)
                .unwrap_or(false)
        })
}
