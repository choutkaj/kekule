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
    clear_fused_lactam_bridge_ring_atoms(mol, ring_set.rings());
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
    }

    mol.perception.aromaticity = ComputedState::Fresh;
    Ok(())
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
            .flat_map(|index| rings[*index].atoms.iter().copied())
            .collect();
    }

    let mut atoms = BTreeSet::new();
    for index in indexes {
        let ring = &rings[*index];
        let exocyclic_pi_count = ring_exocyclic_pi_bond_count(mol, ring);
        let contains_nitrogen = ring_contains_element(mol, ring, "N");
        if exocyclic_pi_count == 0
            || !contains_nitrogen && exocyclic_pi_count <= 1
            || contains_nitrogen && ring_hetero_donor_count(mol, ring) >= 2
            || ring.atoms.len() >= 6 && contains_nitrogen
        {
            atoms.extend(ring.atoms.iter().copied());
        }
    }
    atoms
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
    let mut electrons = 0u8;
    for atom_id in &ring.atoms {
        let atom = mol.atom(*atom_id).expect("ring atom should be live");
        electrons += match atom.element.symbol() {
            "B" | "C" => 1,
            "N" => {
                if atom.explicit_hydrogens > 0 {
                    2
                } else {
                    1
                }
            }
            "O" | "S" | "Se" | "Te" | "P" => 2,
            _ => return Err(AromaticityError::UnsupportedElement(*atom_id)),
        };
    }
    Ok(electrons)
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
