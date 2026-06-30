use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use super::*;
use crate::core::*;

const MAX_FUSED_AROMATIC_COMBINATION_RINGS: usize = 6;
const MAX_FUSED_AROMATIC_RING_SIZE: usize = 24;
const LARGE_FUSED_RING_SYSTEM_SEARCH_LIMIT: usize = 300;

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

#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AromaticElectronDonorType {
    Vacant,
    One,
    Two,
    OneOrTwo,
    Any,
    None,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct AromaticRingAtomDonor {
    atom: AtomId,
    donor: AromaticElectronDonorType,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct AromaticRingDonorAnalysis {
    atoms: Vec<AromaticRingAtomDonor>,
    fixed_electron_count: Option<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RingAromaticityAnalysis {
    aromatic: AromaticRingDonorAnalysis,
    localized: Option<AromaticRingDonorAnalysis>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RdkitAromaticCandidateOptions {
    allow_third_row: bool,
    allow_triple_bonds: bool,
    allow_higher_exceptions: bool,
    only_carbon_or_nitrogen: bool,
    allow_exocyclic_multiple_bonds: bool,
}

impl Default for RdkitAromaticCandidateOptions {
    fn default() -> Self {
        Self {
            allow_third_row: true,
            allow_triple_bonds: true,
            allow_higher_exceptions: true,
            only_carbon_or_nitrogen: false,
            allow_exocyclic_multiple_bonds: true,
        }
    }
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

    let ring_analyses = ring_set
        .rings()
        .iter()
        .map(|ring| RingAromaticityAnalysis::new(mol, ring))
        .collect::<std::result::Result<Vec<_>, AromaticityError>>()?;
    let ring_aromatic = ring_analyses
        .iter()
        .map(RingAromaticityAnalysis::is_huckel_aromatic)
        .collect::<Vec<_>>();
    let non_aromatic_fusion_singles = ring_set
        .rings()
        .iter()
        .flat_map(|ring| ring.bonds.iter().copied())
        .filter(|bond_id| {
            mol.bond(*bond_id)
                .is_ok_and(|bond| matches!(bond.order, BondOrder::Single))
        })
        .filter(|bond_id| {
            let containing_rings = ring_set
                .rings()
                .iter()
                .enumerate()
                .filter(|(_, ring)| ring.bonds.contains(bond_id))
                .collect::<Vec<_>>();
            if !bond_has_candidate_carbon_endpoints_in_containing_ring(
                mol,
                *bond_id,
                &containing_rings,
                &ring_analyses,
            ) {
                return false;
            }
            if containing_rings.iter().any(|(index, ring)| {
                ring_is_all_candidate_carbon(mol, ring, &ring_analyses[*index])
            }) {
                return false;
            }
            containing_rings.iter().any(|(index, ring)| {
                ring_protects_non_aromatic_fusion_single(
                    mol,
                    ring,
                    &ring_analyses[*index],
                    ring_aromatic[*index],
                    containing_rings.len(),
                )
            })
        })
        .collect::<BTreeSet<_>>();
    let protected_non_aromatic_bonds = imported_explicit_aromatic_singles
        .union(&non_aromatic_fusion_singles)
        .copied()
        .collect::<BTreeSet<_>>();

    for (ring, aromatic) in ring_set.rings().iter().zip(ring_aromatic) {
        if aromatic {
            mark_aromatic_atoms_and_bonds(
                mol,
                ring.atoms.iter().copied(),
                ring.bonds.iter().copied(),
                &protected_non_aromatic_bonds,
            );
        }
    }
    perceive_fused_aromatic_components(
        mol,
        ring_set.rings(),
        &ring_analyses,
        &protected_non_aromatic_bonds,
    )?;
    perceive_fused_single_exocyclic_carbon_rings(
        mol,
        ring_set.rings(),
        &ring_analyses,
        &protected_non_aromatic_bonds,
    );
    clear_terminal_chalcogen_oxo_ring_atoms(mol, ring_set.rings(), &ring_analyses);
    clear_ring_oxo_chalcogen_atoms(mol, ring_set.rings(), &ring_analyses);
    clear_fused_lactam_bridge_ring_atoms(mol, ring_set.rings(), &ring_analyses);
    clear_imide_carbonyl_ring_atoms(mol, ring_set.rings(), &ring_analyses);
    clear_fused_lactam_enone_atoms(mol, ring_set.rings(), &ring_analyses);
    clear_saturated_fused_lactam_carbonyl_ring_atoms(mol, ring_set.rings(), &ring_analyses);
    clear_fused_lactone_bridge_ring_atoms(mol, ring_set.rings(), &ring_analyses);
    clear_saturated_fused_ether_bridge_atoms(mol, ring_set.rings(), &ring_analyses);
    clear_saturated_tertiary_amine_ring_atoms(mol, ring_set.rings(), &ring_analyses);
    clear_exocyclic_alkene_chalcogen_ring_atoms(mol, ring_set.rings(), &ring_analyses);
    clear_saturated_chalcogen_bridge_atoms(mol, ring_set.rings(), &ring_analyses);
    clear_fused_carbonyl_bridge_atoms(mol, ring_set.rings(), &ring_analyses);
    clear_saturated_fused_carbon_ring_atoms(mol, ring_set.rings(), &ring_analyses);
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

fn mark_aromatic_atoms_and_bonds(
    mol: &mut Molecule,
    atoms: impl IntoIterator<Item = AtomId>,
    bonds: impl IntoIterator<Item = BondId>,
    protected_non_aromatic_bonds: &BTreeSet<BondId>,
) {
    for atom_id in atoms {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = true;
        }
    }
    for bond_id in bonds {
        if let Some(bond) = mol.bonds[bond_id.index()].as_mut() {
            bond.aromatic = !protected_non_aromatic_bonds.contains(&bond_id);
        }
    }
}

fn mark_aromatic_atom_set_with_internal_bonds(
    mol: &mut Molecule,
    atoms: &BTreeSet<AtomId>,
    bonds: impl IntoIterator<Item = BondId>,
    protected_non_aromatic_bonds: &BTreeSet<BondId>,
) {
    for atom_id in atoms {
        if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
            atom.aromatic = true;
        }
    }
    for bond_id in bonds {
        if let Some(bond) = mol.bonds[bond_id.index()].as_mut() {
            bond.aromatic = atoms.contains(&bond.a())
                && atoms.contains(&bond.b())
                && !protected_non_aromatic_bonds.contains(&bond_id);
        }
    }
}

fn mark_aromatic_fused_ring_system(
    mol: &mut Molecule,
    rings: &[Ring],
    indexes: &[usize],
    protected_non_aromatic_bonds: &BTreeSet<BondId>,
) {
    for bond_id in fused_perimeter_bonds(rings, indexes) {
        let Ok(bond) = mol.bond(bond_id) else {
            continue;
        };
        let [left, right] = [bond.a(), bond.b()];
        for atom_id in [left, right] {
            if let Some(atom) = mol.atoms[atom_id.index()].as_mut() {
                atom.aromatic = true;
            }
        }
        if let Some(bond) = mol.bonds[bond_id.index()].as_mut() {
            bond.aromatic = !protected_non_aromatic_bonds.contains(&bond_id);
        }
    }
}

fn fused_perimeter_bonds(rings: &[Ring], indexes: &[usize]) -> BTreeSet<BondId> {
    let mut bond_counts = BTreeMap::<BondId, usize>::new();
    for index in indexes {
        for bond_id in &rings[*index].bonds {
            *bond_counts.entry(*bond_id).or_default() += 1;
        }
    }
    bond_counts
        .into_iter()
        .filter_map(|(bond_id, count)| (count == 1).then_some(bond_id))
        .collect()
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

fn bond_has_candidate_carbon_endpoints_in_containing_ring(
    mol: &Molecule,
    bond_id: BondId,
    containing_rings: &[(usize, &Ring)],
    ring_analyses: &[RingAromaticityAnalysis],
) -> bool {
    let Ok(bond) = mol.bond(bond_id) else {
        return false;
    };
    containing_rings.iter().any(|(index, ring)| {
        ring.atoms.contains(&bond.a())
            && ring.atoms.contains(&bond.b())
            && ring_analyses[*index].localized_atom_is_element_candidate(mol, bond.a(), "C")
            && ring_analyses[*index].localized_atom_is_element_candidate(mol, bond.b(), "C")
    })
}

fn ring_protects_non_aromatic_fusion_single(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
    ring_is_aromatic: bool,
    containing_ring_count: usize,
) -> bool {
    if ring_is_aromatic {
        return false;
    }
    let active_hetero_donors = analysis.localized_active_hetero_donor_count(mol);
    let non_donor_five_ring = ring.atoms.len() == 5
        && !ring_is_all_candidate_carbon(mol, ring, analysis)
        && !analysis.localized_has_active_chalcogen_donor(mol)
        && active_hetero_donors < 2
        && !analysis.has_nitrogen_lone_pair_donor(mol);
    let multi_hetero_dione_ring = ring.atoms.len() == 6
        && containing_ring_count > 1
        && active_hetero_donors >= 2
        && ring_terminal_exocyclic_pi_bond_count(mol, ring) >= 2;
    non_donor_five_ring || multi_hetero_dione_ring
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

fn clear_saturated_chalcogen_bridge_atoms(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        if ring.atoms.len() < 6
            || !rings
                .iter()
                .enumerate()
                .any(|(other_index, other)| other_index != index && rings_share_bond(ring, other))
        {
            continue;
        }
        for atom_id in &ring.atoms {
            if is_saturated_chalcogen_bridge_cleanup_candidate(
                mol,
                ring,
                &ring_analyses[index],
                *atom_id,
            ) {
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
            if !atom_is_rdkit_aromatic_element_candidate(mol, nitrogen_id, "N")
                || !atom_is_rdkit_aromatic_element_candidate(mol, carbon_id, "C")
            {
                return None;
            }
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
            if atom.element.symbol() != "C" || !atom.aromatic || atom.formal_charge < 0 {
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
            if atom.element.symbol() != "C"
                || !atom.aromatic
                || !atom_is_rdkit_aromatic_candidate(mol, atom_id, atom)
            {
                return None;
            }
            let bonds = mol
                .incident_bonds(atom_id)
                .ok()?
                .map(|(_, bond)| bond)
                .collect::<Vec<_>>();
            let has_imine_nitrogen = bonds.iter().any(|bond| {
                let neighbor_id = bond.other_atom(atom_id);
                matches!(bond.order, BondOrder::Double)
                    && atom_is_rdkit_aromatic_element_candidate(mol, neighbor_id, "N")
            });
            let has_single_nitrogen = bonds.iter().any(|bond| {
                let neighbor_id = bond.other_atom(atom_id);
                matches!(bond.order, BondOrder::Single)
                    && atom_is_rdkit_aromatic_element_candidate(mol, neighbor_id, "N")
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
                    atom_is_rdkit_aromatic_element_candidate(mol, neighbor_id, "N")
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

fn clear_ring_oxo_chalcogen_atoms(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        for atom_id in &ring.atoms {
            if is_ring_oxo_chalcogen_cleanup_candidate(mol, ring, &ring_analyses[index], *atom_id) {
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

fn is_ring_oxo_chalcogen_cleanup_candidate(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
    atom_id: AtomId,
) -> bool {
    is_ring_oxo_chalcogen(mol, ring, atom_id)
        && mol
            .atom(atom_id)
            .is_ok_and(atom_passes_rdkit_aromatic_radical_eligibility)
        && !analysis.localized_atom_has_active_chalcogen_donor(mol, atom_id)
}

fn clear_terminal_chalcogen_oxo_ring_atoms(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        if ring.atoms.len() != 5
            || !ring_analyses[index].localized_has_active_element_donor(mol, "N")
            || !ring_has_terminal_inactive_chalcogen_exocyclic_pi_bond(
                mol,
                ring,
                &ring_analyses[index],
            )
        {
            continue;
        }
        for atom_id in &ring.atoms {
            let retained_by_other_ring = rings.iter().enumerate().any(|(other_index, other)| {
                other_index != index
                    && other.atoms.contains(atom_id)
                    && !ring_has_terminal_inactive_chalcogen_exocyclic_pi_bond(
                        mol,
                        other,
                        &ring_analyses[other_index],
                    )
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

fn clear_fused_lactam_bridge_ring_atoms(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        if !ring_is_fused_lactam_bridge_cleanup_candidate(mol, ring, &ring_analyses[index]) {
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

fn ring_is_fused_lactam_bridge_cleanup_candidate(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
) -> bool {
    ring.atoms.len() >= 6
        && analysis.localized_has_active_element_donor(mol, "N")
        && analysis.localized_has_active_element_donor(mol, "O")
        && ring_terminal_exocyclic_pi_bond_count(mol, ring) > 0
}

fn clear_imide_carbonyl_ring_atoms(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
) {
    let atoms_to_clear = rings
        .iter()
        .enumerate()
        .filter(|(index, ring)| {
            ring_is_imide_carbonyl_cleanup_candidate(
                mol,
                rings,
                *index,
                ring,
                &ring_analyses[*index],
                ring_analyses,
            )
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

fn ring_is_imide_carbonyl_cleanup_candidate(
    mol: &Molecule,
    rings: &[Ring],
    index: usize,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
    ring_analyses: &[RingAromaticityAnalysis],
) -> bool {
    analysis.localized_has_element_candidate(mol, "N")
        && (ring_has_or_is_fused_to_cationic_nitrogen(mol, rings, index, ring_analyses)
            || ring.atoms.len() == 5 && ring_has_imide_nitrogen(mol, ring, analysis))
        && ring_terminal_exocyclic_pi_bond_count(mol, ring) >= 2
}

fn clear_fused_lactam_enone_atoms(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
) {
    let atoms_to_clear = rings
        .iter()
        .enumerate()
        .filter(|(index, ring)| {
            ring_is_fused_lactam_enone_cleanup_candidate(
                mol,
                rings,
                *index,
                ring,
                &ring_analyses[*index],
            )
        })
        .flat_map(|(index, ring)| {
            let analysis = &ring_analyses[index];
            ring.atoms
                .iter()
                .copied()
                .filter(|atom_id| is_aromatic_lactam_enone_carbon(mol, ring, analysis, *atom_id))
                .collect::<Vec<_>>()
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

fn ring_is_fused_lactam_enone_cleanup_candidate(
    mol: &Molecule,
    rings: &[Ring],
    index: usize,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
) -> bool {
    analysis.localized_has_active_element_donor(mol, "N")
        && !analysis.localized_has_active_chalcogen_donor(mol)
        && ring_terminal_exocyclic_pi_bond_count(mol, ring) > 0
        && rings
            .iter()
            .enumerate()
            .any(|(other_index, other)| other_index != index && rings_share_bond(ring, other))
}

fn is_aromatic_lactam_enone_carbon(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
    atom_id: AtomId,
) -> bool {
    mol.atom(atom_id)
        .is_ok_and(|atom| atom.element.symbol() == "C" && atom.aromatic)
        && mol
            .incident_bonds(atom_id)
            .ok()
            .into_iter()
            .flatten()
            .any(|(bond_id, bond)| {
                ring.bonds.contains(&bond_id)
                    && matches!(bond.order, BondOrder::Double)
                    && mol
                        .atom(bond.other_atom(atom_id))
                        .is_ok_and(|other| other.element.symbol() == "C" && !other.aromatic)
            })
        && mol
            .incident_bonds(atom_id)
            .ok()
            .into_iter()
            .flatten()
            .any(|(bond_id, bond)| {
                ring.bonds.contains(&bond_id)
                    && matches!(bond.order, BondOrder::Single)
                    && analysis.localized_atom_has_active_element_donor(
                        mol,
                        bond.other_atom(atom_id),
                        "N",
                    )
            })
}

fn clear_saturated_fused_lactam_carbonyl_ring_atoms(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        if !ring_is_saturated_fused_lactam_carbonyl_cleanup_candidate(
            mol,
            rings,
            index,
            ring,
            &ring_analyses[index],
        ) {
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

fn ring_is_saturated_fused_lactam_carbonyl_cleanup_candidate(
    mol: &Molecule,
    rings: &[Ring],
    index: usize,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
) -> bool {
    rings
        .iter()
        .enumerate()
        .any(|(other_index, other)| other_index != index && rings_share_bond(ring, other))
        && analysis.localized_has_active_element_donor(mol, "N")
        && ring_terminal_exocyclic_pi_bond_count(mol, ring) > 0
        && ring_has_saturated_carbon_atom(mol, ring)
}

fn is_saturated_ring_carbon(mol: &Molecule, ring: &Ring, atom_id: AtomId) -> bool {
    mol.atom(atom_id).is_ok_and(|atom| {
        atom.element.symbol() == "C" && atom_passes_rdkit_aromatic_radical_eligibility(atom)
    }) && !ring_atom_has_pi_bond(mol, ring, atom_id)
        && !atom_has_exocyclic_pi_bond(mol, ring, atom_id)
}

fn ring_has_imide_nitrogen(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
) -> bool {
    ring.atoms.iter().any(|atom_id| {
        let Ok(atom) = mol.atom(*atom_id) else {
            return false;
        };
        analysis.localized_atom_is_element_candidate(mol, *atom_id, "N")
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

fn ring_has_or_is_fused_to_cationic_nitrogen(
    mol: &Molecule,
    rings: &[Ring],
    index: usize,
    ring_analyses: &[RingAromaticityAnalysis],
) -> bool {
    ring_has_cationic_nitrogen(mol, &rings[index], &ring_analyses[index])
        || rings.iter().enumerate().any(|(other_index, other)| {
            other_index != index
                && rings_share_bond(&rings[index], other)
                && ring_has_cationic_nitrogen(mol, other, &ring_analyses[other_index])
        })
}

fn ring_has_cationic_nitrogen(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
) -> bool {
    ring.atoms.iter().any(|atom_id| {
        mol.atom(*atom_id).is_ok_and(|atom| atom.formal_charge > 0)
            && analysis.localized_atom_is_element_candidate(mol, *atom_id, "N")
    })
}

fn clear_fused_lactone_bridge_ring_atoms(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        let fused = rings
            .iter()
            .enumerate()
            .any(|(other_index, other)| other_index != index && rings_share_bond(ring, other));
        if !ring_is_fused_lactone_bridge_cleanup_candidate(mol, ring, &ring_analyses[index], fused)
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

fn ring_is_fused_lactone_bridge_cleanup_candidate(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
    fused: bool,
) -> bool {
    fused
        && ring.atoms.len() >= 6
        && analysis.localized_has_active_chalcogen_donor(mol)
        && !analysis.localized_all_atoms_are_candidates()
        && ring
            .atoms
            .iter()
            .any(|atom_id| is_active_chalcogen_bridge_without_pi(mol, ring, analysis, *atom_id))
        && ring_terminal_exocyclic_pi_bond_count(mol, ring) > 0
}

fn clear_saturated_fused_ether_bridge_atoms(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
) {
    let atoms_to_clear = rings
        .iter()
        .enumerate()
        .filter(|(index, ring)| {
            ring_is_saturated_fused_ether_bridge_cleanup_candidate(
                mol,
                rings,
                *index,
                ring,
                &ring_analyses[*index],
            )
        })
        .flat_map(|(index, ring)| {
            let analysis = &ring_analyses[index];
            ring.atoms.iter().copied().filter(|atom_id| {
                mol.atom(*atom_id)
                    .is_ok_and(|atom| atom.element.symbol() == "O" && atom.aromatic)
                    && is_active_element_bridge_without_pi(mol, ring, analysis, *atom_id, "O")
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

fn ring_is_saturated_fused_ether_bridge_cleanup_candidate(
    mol: &Molecule,
    rings: &[Ring],
    index: usize,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
) -> bool {
    ring.atoms.len() >= 5
        && analysis.localized_has_active_element_donor(mol, "O")
        && analysis.localized_active_hetero_donor_count(mol) == 1
        && !analysis.localized_all_atoms_are_candidates()
        && ring_terminal_exocyclic_pi_bond_count(mol, ring) == 0
        && rings
            .iter()
            .enumerate()
            .any(|(other_index, other)| other_index != index && rings_share_bond(ring, other))
}

#[cfg(test)]
fn is_chalcogen_bridge_without_pi(mol: &Molecule, ring: &Ring, atom_id: AtomId) -> bool {
    mol.atom(atom_id).is_ok_and(|atom| {
        matches!(atom.element.symbol(), "O" | "S" | "Se" | "Te")
            && atom_passes_rdkit_aromatic_radical_eligibility(atom)
    }) && !ring_atom_has_pi_bond(mol, ring, atom_id)
        && !atom_has_exocyclic_pi_bond(mol, ring, atom_id)
}

fn is_active_chalcogen_bridge_without_pi(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
    atom_id: AtomId,
) -> bool {
    !ring_atom_has_pi_bond(mol, ring, atom_id)
        && !atom_has_exocyclic_pi_bond(mol, ring, atom_id)
        && analysis.localized_atom_has_active_chalcogen_donor(mol, atom_id)
}

fn is_active_element_bridge_without_pi(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
    atom_id: AtomId,
    symbol: &str,
) -> bool {
    !ring_atom_has_pi_bond(mol, ring, atom_id)
        && !atom_has_exocyclic_pi_bond(mol, ring, atom_id)
        && analysis.localized_atom_has_active_element_donor(mol, atom_id, symbol)
}

fn clear_saturated_tertiary_amine_ring_atoms(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        if !ring_has_saturated_tertiary_amine_without_donor_chalcogen(
            mol,
            ring,
            &ring_analyses[index],
        ) {
            continue;
        }
        for atom_id in &ring.atoms {
            let retained_by_other_ring = rings.iter().enumerate().any(|(other_index, other)| {
                other_index != index
                    && other.atoms.contains(atom_id)
                    && !ring_has_saturated_tertiary_amine_without_donor_chalcogen(
                        mol,
                        other,
                        &ring_analyses[other_index],
                    )
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

fn clear_exocyclic_alkene_chalcogen_ring_atoms(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        let analysis = &ring_analyses[index];
        if !analysis.localized_has_active_element_donor(mol, "N")
            || !analysis.localized_has_active_chalcogen_donor(mol)
        {
            continue;
        }
        for atom_id in ring.atoms.iter().copied().filter(|atom_id| {
            is_aromatic_ring_carbon_with_exocyclic_carbon_pi_bond(mol, ring, analysis, *atom_id)
        }) {
            atoms_to_clear.insert(atom_id);
            atoms_to_clear.extend(active_hetero_neighbors_for_exocyclic_alkene_cleanup(
                mol, ring, analysis, atom_id,
            ));
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

fn active_hetero_neighbors_for_exocyclic_alkene_cleanup(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
    atom_id: AtomId,
) -> BTreeSet<AtomId> {
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .filter_map(|(_, bond)| {
            let neighbor_id = bond.other_atom(atom_id);
            (ring.atoms.contains(&neighbor_id)
                && (analysis.localized_atom_has_active_element_donor(mol, neighbor_id, "N")
                    || analysis.localized_atom_has_active_chalcogen_donor(mol, neighbor_id)))
            .then_some(neighbor_id)
        })
        .collect()
}

fn is_aromatic_ring_carbon_with_exocyclic_carbon_pi_bond(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
    atom_id: AtomId,
) -> bool {
    mol.atom(atom_id)
        .is_ok_and(|atom| atom.element.symbol() == "C" && atom.aromatic)
        && ring_atom_has_active_nitrogen_and_chalcogen_neighbors(mol, ring, analysis, atom_id)
        && mol
            .incident_bonds(atom_id)
            .ok()
            .into_iter()
            .flatten()
            .any(|(bond_id, bond)| {
                !ring.bonds.contains(&bond_id)
                    && !ring.atoms.contains(&bond.other_atom(atom_id))
                    && matches!(bond.order, BondOrder::Double | BondOrder::Aromatic)
                    && mol
                        .atom(bond.other_atom(atom_id))
                        .is_ok_and(|other| other.element.symbol() == "C" && !other.aromatic)
            })
}

fn ring_atom_has_active_nitrogen_and_chalcogen_neighbors(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
    atom_id: AtomId,
) -> bool {
    let ring_neighbors = mol
        .incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .filter_map(|(bond_id, bond)| ring.bonds.contains(&bond_id).then_some(bond))
        .map(|bond| bond.other_atom(atom_id))
        .collect::<Vec<_>>();
    ring_neighbors
        .iter()
        .any(|neighbor| analysis.localized_atom_has_active_element_donor(mol, *neighbor, "N"))
        && ring_neighbors
            .iter()
            .any(|neighbor| analysis.localized_atom_has_active_chalcogen_donor(mol, *neighbor))
}

fn clear_saturated_fused_carbon_ring_atoms(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        let analysis = &ring_analyses[index];
        let fused = rings
            .iter()
            .enumerate()
            .any(|(other_index, other)| other_index != index && rings_share_bond(ring, other));
        if !fused
            || !ring_is_all_candidate_carbon(mol, ring, analysis)
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

fn clear_fused_carbonyl_bridge_atoms(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
) {
    let mut atoms_to_clear = BTreeSet::new();
    for (index, ring) in rings.iter().enumerate() {
        let analysis = &ring_analyses[index];
        let fused = rings
            .iter()
            .enumerate()
            .any(|(other_index, other)| other_index != index && rings_share_bond(ring, other));
        if !fused || !ring_is_all_candidate_carbon(mol, ring, analysis) {
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
    ring.atoms
        .iter()
        .any(|atom_id| is_saturated_ring_carbon(mol, ring, *atom_id))
}

fn perceive_fused_single_exocyclic_carbon_rings(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
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
                && ring_analyses[*index].localized_all_atoms_are_candidates()
                && ((ring.atoms.len() == 6
                    && ring_is_all_candidate_carbon(mol, ring, &ring_analyses[*index])
                    && ring_exocyclic_pi_bond_count(mol, ring) > 0
                    && ring_terminal_exocyclic_pi_bond_count(mol, ring) <= 1
                    && ring_pi_bond_count(mol, ring) >= 1)
                    || (ring.atoms.len() == 4
                        && ring_is_all_candidate_carbon(mol, ring, &ring_analyses[*index])
                        && ring_terminal_exocyclic_pi_bond_count(mol, ring) >= 2)
                    || (ring.atoms.len() == 5
                        && ring_analyses[*index].localized_has_active_element_donor(mol, "N")
                        && ring_analyses[*index].localized_has_active_chalcogen_donor(mol)
                        && ring_analyses[*index]
                            .localized_has_candidate_carbon_electronegative_exocyclic_pi_bond(
                                mol, ring,
                            )
                        && ring_exocyclic_pi_bond_count(mol, ring) > 0))
        })
        .map(|(_, ring)| ring.clone())
        .collect::<Vec<_>>();

    for ring in selected {
        mark_aromatic_atoms_and_bonds(
            mol,
            ring.atoms.iter().copied(),
            ring.bonds.iter().copied(),
            protected_non_aromatic_bonds,
        );
    }
}

fn perceive_fused_aromatic_components(
    mol: &mut Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
    protected_non_aromatic_bonds: &BTreeSet<BondId>,
) -> std::result::Result<(), AromaticityError> {
    let candidates = rings
        .iter()
        .enumerate()
        .filter_map(|(index, ring)| {
            aromatic_fused_candidate_from_analysis(mol, ring, &ring_analyses[index])
                .then_some(index)
        })
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
        let all_candidate_carbon_component =
            fused_component_is_all_candidate_carbon(mol, &component, rings, ring_analyses, indexes);
        let component_has_exocyclic_pi = ring_exocyclic_pi_bond_count(mol, &component) > 0;
        if indexes.len() <= MAX_FUSED_AROMATIC_COMBINATION_RINGS {
            let aromaticity_ring = fused_component_aromaticity_ring(rings, indexes);
            let analysis = localized_ring_donor_analysis(mol, &aromaticity_ring)?;
            if analysis.is_fused_huckel_aromatic() {
                if all_candidate_carbon_component {
                    mark_aromatic_fused_ring_system(
                        mol,
                        rings,
                        indexes,
                        protected_non_aromatic_bonds,
                    );
                } else {
                    mark_aromatic_atoms_and_bonds(
                        mol,
                        component.atoms.iter().copied(),
                        component.bonds.iter().copied(),
                        protected_non_aromatic_bonds,
                    );
                }
                continue;
            }
        }
        for subset in aromatic_fused_ring_subsets(mol, rings, indexes)? {
            let subset_ring = fused_component_ring(rings, &subset);
            if fused_component_is_all_candidate_carbon(
                mol,
                &subset_ring,
                rings,
                ring_analyses,
                &subset,
            ) {
                mark_aromatic_fused_ring_system(mol, rings, &subset, protected_non_aromatic_bonds);
            } else {
                mark_aromatic_atoms_and_bonds(
                    mol,
                    subset_ring.atoms.iter().copied(),
                    subset_ring.bonds.iter().copied(),
                    protected_non_aromatic_bonds,
                );
            }
        }
        if all_candidate_carbon_component && component_has_exocyclic_pi {
            let aromatic_atoms =
                atoms_in_limited_terminal_exocyclic_pi_rings(mol, rings, indexes, 1);
            if aromatic_atoms.len() >= 6 {
                mark_aromatic_atom_set_with_internal_bonds(
                    mol,
                    &aromatic_atoms,
                    component.bonds.iter().copied(),
                    protected_non_aromatic_bonds,
                );
            }
        } else if fused_component_is_carbon_with_candidate_nitrogen(
            mol,
            &component,
            rings,
            ring_analyses,
            indexes,
        ) {
            let terminal_atoms_retained =
                terminal_exocyclic_atoms_in_nitrogen_rings(mol, rings, ring_analyses, indexes);
            let atoms_retained_by_ring_context =
                atoms_in_nitrogen_or_terminal_pi_free_rings(mol, rings, ring_analyses, indexes);
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
                mark_aromatic_atom_set_with_internal_bonds(
                    mol,
                    &aromatic_atoms,
                    component.bonds.iter().copied(),
                    protected_non_aromatic_bonds,
                );
            }
        }
    }
    Ok(())
}

fn aromatic_fused_ring_subsets(
    mol: &Molecule,
    rings: &[Ring],
    indexes: &[usize],
) -> std::result::Result<Vec<Vec<usize>>, AromaticityError> {
    if indexes.len() < 3 {
        return Ok(Vec::new());
    }
    let mut accepted = Vec::new();
    let mut done_bonds = BTreeSet::new();
    let all_bonds = indexes
        .iter()
        .flat_map(|index| rings[*index].bonds.iter().copied())
        .collect::<BTreeSet<_>>();
    let max_subset_size = indexes.len().min(MAX_FUSED_AROMATIC_COMBINATION_RINGS);
    for subset_size in 2..=max_subset_size {
        if subset_size == indexes.len() {
            continue;
        }
        if subset_size > 2 && indexes.len() > LARGE_FUSED_RING_SYSTEM_SEARCH_LIMIT {
            break;
        }
        for subset in connected_ring_subsets(rings, indexes, subset_size) {
            let aromaticity_ring = fused_component_aromaticity_ring(rings, &subset);
            let analysis = localized_ring_donor_analysis(mol, &aromaticity_ring)?;
            if analysis.is_fused_huckel_aromatic() {
                done_bonds.extend(fused_perimeter_bonds(rings, &subset));
                accepted.push(subset);
                if done_bonds.len() >= all_bonds.len() {
                    return Ok(accepted);
                }
            }
        }
    }
    Ok(accepted)
}

fn connected_ring_subsets(
    rings: &[Ring],
    indexes: &[usize],
    subset_size: usize,
) -> Vec<Vec<usize>> {
    let mut subsets = Vec::new();
    let mut current = Vec::with_capacity(subset_size);
    collect_connected_ring_subsets(rings, indexes, subset_size, 0, &mut current, &mut subsets);
    subsets
}

fn collect_connected_ring_subsets(
    rings: &[Ring],
    indexes: &[usize],
    subset_size: usize,
    start: usize,
    current: &mut Vec<usize>,
    subsets: &mut Vec<Vec<usize>>,
) {
    if current.len() == subset_size {
        if ring_subset_is_connected(rings, current) {
            subsets.push(current.clone());
        }
        return;
    }
    for position in start..indexes.len() {
        current.push(indexes[position]);
        collect_connected_ring_subsets(rings, indexes, subset_size, position + 1, current, subsets);
        current.pop();
    }
}

fn ring_subset_is_connected(rings: &[Ring], indexes: &[usize]) -> bool {
    let mut visited = BTreeSet::new();
    let mut stack = vec![indexes[0]];
    while let Some(index) = stack.pop() {
        if !visited.insert(index) {
            continue;
        }
        for other in indexes {
            if !visited.contains(other) && rings_share_bond(&rings[index], &rings[*other]) {
                stack.push(*other);
            }
        }
    }
    visited.len() == indexes.len()
}

#[cfg(test)]
fn fused_component_is_all_carbon(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().all(|atom_id| {
        mol.atom(*atom_id)
            .map(|atom| atom.element.symbol() == "C")
            .unwrap_or(false)
    })
}

fn ring_is_all_candidate_carbon(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
) -> bool {
    ring.atoms
        .iter()
        .all(|atom_id| analysis.localized_atom_is_element_candidate(mol, *atom_id, "C"))
}

fn fused_component_is_carbon_or_candidate_nitrogen(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
) -> bool {
    ring.atoms.iter().all(|atom_id| {
        analysis.localized_atom_is_element_candidate(mol, *atom_id, "C")
            || analysis.localized_atom_is_element_candidate(mol, *atom_id, "N")
    })
}

fn fused_component_is_carbon_with_candidate_nitrogen(
    mol: &Molecule,
    component: &Ring,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
    indexes: &[usize],
) -> bool {
    component.atoms.iter().all(|atom_id| {
        atom_is_localized_element_candidate_in_any_ring(
            mol,
            rings,
            ring_analyses,
            indexes,
            *atom_id,
            "C",
        ) || atom_is_localized_element_candidate_in_any_ring(
            mol,
            rings,
            ring_analyses,
            indexes,
            *atom_id,
            "N",
        )
    })
}

fn fused_component_is_all_candidate_carbon(
    mol: &Molecule,
    component: &Ring,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
    indexes: &[usize],
) -> bool {
    component.atoms.iter().all(|atom_id| {
        atom_is_localized_element_candidate_in_any_ring(
            mol,
            rings,
            ring_analyses,
            indexes,
            *atom_id,
            "C",
        )
    })
}

fn atom_is_localized_element_candidate_in_any_ring(
    mol: &Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
    indexes: &[usize],
    atom_id: AtomId,
    symbol: &str,
) -> bool {
    indexes.iter().any(|index| {
        rings[*index].atoms.contains(&atom_id)
            && ring_analyses[*index].localized_atom_is_element_candidate(mol, atom_id, symbol)
    })
}

fn terminal_exocyclic_atoms_in_nitrogen_rings(
    mol: &Molecule,
    rings: &[Ring],
    ring_analyses: &[RingAromaticityAnalysis],
    indexes: &[usize],
) -> BTreeSet<AtomId> {
    let mut atoms = BTreeSet::new();
    for index in indexes {
        let ring = &rings[*index];
        if !ring_analyses[*index].localized_has_active_element_donor(mol, "N") {
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
    ring_analyses: &[RingAromaticityAnalysis],
    indexes: &[usize],
) -> BTreeSet<AtomId> {
    let has_exocyclic_pi_ring = indexes
        .iter()
        .any(|index| ring_exocyclic_pi_bond_count(mol, &rings[*index]) > 0);
    if !has_exocyclic_pi_ring {
        return indexes
            .iter()
            .filter(|index| {
                !ring_has_saturated_tertiary_amine_without_donor_chalcogen(
                    mol,
                    &rings[**index],
                    &ring_analyses[**index],
                )
            })
            .flat_map(|index| rings[*index].atoms.iter().copied())
            .collect();
    }

    let mut atoms = BTreeSet::new();
    for index in indexes {
        let ring = &rings[*index];
        if ring_has_saturated_tertiary_amine_without_donor_chalcogen(
            mol,
            ring,
            &ring_analyses[*index],
        ) {
            continue;
        }
        let exocyclic_pi_count = ring_exocyclic_pi_bond_count(mol, ring);
        let has_active_nitrogen_donor =
            ring_analyses[*index].localized_has_active_element_donor(mol, "N");
        let active_hetero_donors = ring_analyses[*index].localized_active_hetero_donor_count(mol);
        if exocyclic_pi_count == 0
            || !has_active_nitrogen_donor && exocyclic_pi_count <= 1
            || has_active_nitrogen_donor && active_hetero_donors >= 2
            || ring.atoms.len() >= 6
                && has_active_nitrogen_donor
                && !ring_has_saturated_carbon_atom(mol, ring)
        {
            atoms.extend(ring.atoms.iter().copied());
        }
    }
    atoms
}

fn ring_has_saturated_tertiary_amine_without_donor_chalcogen(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
) -> bool {
    if ring_has_saturated_active_chalcogen_donor(mol, ring, analysis) {
        return false;
    }
    !analysis.localized_all_atoms_are_candidates()
        && analysis.localized_active_hetero_donor_count(mol) == 1
        && ring
            .atoms
            .iter()
            .any(|atom_id| is_saturated_tertiary_amine(mol, ring, *atom_id))
}

fn ring_has_saturated_active_chalcogen_donor(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
) -> bool {
    ring.atoms.iter().any(|atom_id| {
        !ring_atom_has_pi_bond(mol, ring, *atom_id)
            && !atom_has_exocyclic_pi_bond(mol, ring, *atom_id)
            && analysis.localized_atom_has_active_chalcogen_donor(mol, *atom_id)
    })
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

fn aromatic_fused_candidate_from_analysis(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
) -> bool {
    if !analysis.localized_all_atoms_are_candidates() {
        return false;
    }
    if ring.atoms.len() > 7 {
        return ring.atoms.len() <= MAX_FUSED_AROMATIC_RING_SIZE
            && fused_component_is_carbon_or_candidate_nitrogen(mol, ring, analysis)
            && analysis.localized_has_active_anionic_nitrogen_donor(mol)
            && ring_pi_bond_count(mol, ring) >= 4;
    }
    let pi_bonds = ring_pi_bond_count(mol, ring);
    pi_bonds >= 2
        || analysis.localized_active_hetero_donor_count(mol) > 0
            && !ring_has_low_unsaturation_active_chalcogen_bridge_for_fused(mol, ring, analysis)
        || ring.atoms.len() == 5
            && pi_bonds >= 1
            && analysis.localized_has_active_element_donor(mol, "N")
        || ring.atoms.len() == 6
            && pi_bonds >= 1
            && ring_is_all_candidate_carbon(mol, ring, analysis)
}

fn is_saturated_tertiary_amine(mol: &Molecule, ring: &Ring, atom_id: AtomId) -> bool {
    let Ok(atom) = mol.atom(atom_id) else {
        return false;
    };
    atom.element.symbol() == "N"
        && atom.formal_charge == 0
        && atom_passes_rdkit_aromatic_radical_eligibility(atom)
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
    atom.element.symbol() == "C" && atom_passes_rdkit_aromatic_radical_eligibility(atom)
        || matches!(atom.element.symbol(), "S" | "Se" | "Te")
            && atom_passes_rdkit_aromatic_radical_eligibility(atom)
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

fn is_saturated_chalcogen_bridge_cleanup_candidate(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
    atom_id: AtomId,
) -> bool {
    is_saturated_chalcogen_bridge(mol, ring, atom_id)
        && analysis.localized_atom_is_chalcogen_candidate(mol, atom_id)
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

#[cfg(test)]
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

#[cfg(test)]
fn ring_has_chalcogen_donor(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().any(|atom_id| {
        mol.atom(*atom_id)
            .map(|atom| matches!(atom.element.symbol(), "O" | "S" | "Se" | "Te"))
            .unwrap_or(false)
    })
}

fn ring_has_low_unsaturation_active_chalcogen_bridge_for_fused(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
) -> bool {
    ring_pi_bond_count(mol, ring) < 2
        && analysis.localized_active_hetero_donor_count(mol) > 1
        && analysis.localized_has_active_chalcogen_donor(mol)
}

fn ring_exocyclic_pi_bond_count(mol: &Molecule, ring: &Ring) -> usize {
    ring.atoms
        .iter()
        .filter(|atom_id| atom_has_exocyclic_pi_bond(mol, ring, **atom_id))
        .count()
}

#[cfg(test)]
fn ring_has_carbon_electronegative_exocyclic_pi_bond(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().any(|atom_id| {
        let Ok(atom) = mol.atom(*atom_id) else {
            return false;
        };
        atom.element.symbol() == "C"
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
                    atom_is_more_electronegative_than(mol, bond.other_atom(*atom_id), atom)
                })
    })
}

fn ring_has_candidate_carbon_electronegative_exocyclic_pi_bond(
    mol: &Molecule,
    ring: &Ring,
    analysis: &AromaticRingDonorAnalysis,
) -> bool {
    ring.atoms.iter().any(|atom_id| {
        analysis.atom_is_element_candidate(mol, *atom_id, "C")
            && atom_has_electronegative_exocyclic_pi_bond_for_ring_gate(mol, ring, *atom_id)
    })
}

fn atom_has_electronegative_exocyclic_pi_bond_for_ring_gate(
    mol: &Molecule,
    ring: &Ring,
    atom_id: AtomId,
) -> bool {
    let Ok(atom) = mol.atom(atom_id) else {
        return false;
    };
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .any(|(bond_id, bond)| {
            if ring.bonds.contains(&bond_id)
                || ring.atoms.contains(&bond.other_atom(atom_id))
                || !matches!(bond.order, BondOrder::Double | BondOrder::Aromatic)
            {
                return false;
            }
            atom_is_more_electronegative_than(mol, bond.other_atom(atom_id), atom)
        })
}

#[cfg(test)]
fn ring_has_terminal_chalcogen_exocyclic_pi_bond(mol: &Molecule, ring: &Ring) -> bool {
    ring.atoms.iter().any(|atom_id| {
        mol.atom(*atom_id)
            .map(|atom| matches!(atom.element.symbol(), "S" | "Se" | "Te"))
            .unwrap_or(false)
            && atom_has_terminal_exocyclic_pi_bond(mol, ring, *atom_id)
    })
}

fn ring_has_terminal_candidate_chalcogen_exocyclic_pi_bond(
    mol: &Molecule,
    ring: &Ring,
    analysis: &AromaticRingDonorAnalysis,
) -> bool {
    ring.atoms.iter().any(|atom_id| {
        atom_has_terminal_exocyclic_pi_bond(mol, ring, *atom_id)
            && analysis.atom_is_heavy_chalcogen_candidate(mol, *atom_id)
    })
}

fn ring_has_terminal_inactive_chalcogen_exocyclic_pi_bond(
    mol: &Molecule,
    ring: &Ring,
    analysis: &RingAromaticityAnalysis,
) -> bool {
    ring.atoms.iter().any(|atom_id| {
        atom_has_terminal_exocyclic_pi_bond(mol, ring, *atom_id)
            && mol.atom(*atom_id).is_ok_and(|atom| {
                matches!(atom.element.symbol(), "S" | "Se" | "Te")
                    && atom_passes_rdkit_aromatic_radical_eligibility(atom)
            })
            && !analysis.localized_atom_has_active_chalcogen_donor(mol, *atom_id)
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

fn fused_component_aromaticity_ring(rings: &[Ring], indexes: &[usize]) -> Ring {
    let mut atom_counts = BTreeMap::<AtomId, usize>::new();
    let mut bonds = BTreeSet::new();
    for index in indexes {
        for atom_id in &rings[*index].atoms {
            *atom_counts.entry(*atom_id).or_default() += 1;
        }
        bonds.extend(rings[*index].bonds.iter().copied());
    }
    Ring {
        atoms: atom_counts
            .into_iter()
            .filter_map(|(atom_id, count)| (count <= 2).then_some(atom_id))
            .collect(),
        bonds: bonds.into_iter().collect(),
    }
}

fn aromatic_ring_donor_analysis(
    mol: &Molecule,
    ring: &Ring,
    localized: &AromaticRingDonorAnalysis,
) -> std::result::Result<AromaticRingDonorAnalysis, AromaticityError> {
    if ring_has_aromatic_order(mol, ring) {
        return aromatic_order_ring_donor_analysis(mol, ring, localized);
    }
    let active_hetero_donors = localized.active_hetero_donor_count(mol);
    let has_active_chalcogen_donor = localized.has_active_chalcogen_donor(mol);
    let has_active_nitrogen_donor = localized.has_active_element_donor(mol, "N");
    if ring.atoms.len() == 7 && has_active_chalcogen_donor {
        return Ok(AromaticRingDonorAnalysis::non_aromatic());
    }
    if ring.atoms.len() == 5 && ring_pi_bond_count(mol, ring) < 2 && !has_active_nitrogen_donor {
        return Ok(AromaticRingDonorAnalysis::non_aromatic());
    }
    if ring.atoms.len() == 5
        && has_active_nitrogen_donor
        && ring_exocyclic_pi_bond_count(mol, ring) > 0
        && ((active_hetero_donors < 2 && !has_active_chalcogen_donor)
            || (has_active_chalcogen_donor
                && !ring_has_candidate_carbon_electronegative_exocyclic_pi_bond(
                    mol, ring, localized,
                )))
    {
        return Ok(AromaticRingDonorAnalysis::non_aromatic());
    }
    if ring.atoms.len() == 5
        && has_active_nitrogen_donor
        && ring_has_terminal_candidate_chalcogen_exocyclic_pi_bond(mol, ring, localized)
    {
        return Ok(AromaticRingDonorAnalysis::non_aromatic());
    }
    if ring.atoms.len() > 5
        && ring_pi_bond_count(mol, ring) + ring_terminal_exocyclic_pi_bond_count(mol, ring) < 2
        && active_hetero_donors > 1
        && has_active_chalcogen_donor
    {
        return Ok(AromaticRingDonorAnalysis::non_aromatic());
    }
    if ring.atoms.len() > 5
        && ring_exocyclic_pi_bond_count(mol, ring) == 0
        && has_active_chalcogen_donor
        && has_active_nitrogen_donor
        && active_hetero_donors > 1
    {
        return Ok(AromaticRingDonorAnalysis::non_aromatic());
    }

    if ring_pi_bond_count(mol, ring) == 0 {
        return Ok(AromaticRingDonorAnalysis::non_aromatic());
    }

    Ok(localized.clone())
}

fn localized_ring_donor_analysis(
    mol: &Molecule,
    ring: &Ring,
) -> std::result::Result<AromaticRingDonorAnalysis, AromaticityError> {
    let atoms = ring
        .atoms
        .iter()
        .map(|atom_id| {
            localized_ring_atom_donor_type(mol, ring, *atom_id).map(|donor| AromaticRingAtomDonor {
                atom: *atom_id,
                donor,
            })
        })
        .collect::<std::result::Result<Vec<_>, _>>()?;
    Ok(AromaticRingDonorAnalysis {
        atoms,
        fixed_electron_count: None,
    })
}

impl RingAromaticityAnalysis {
    fn new(mol: &Molecule, ring: &Ring) -> std::result::Result<Self, AromaticityError> {
        let localized = localized_ring_donor_analysis(mol, ring)?;
        Ok(Self {
            aromatic: aromatic_ring_donor_analysis(mol, ring, &localized)?,
            localized: Some(localized),
        })
    }

    fn is_huckel_aromatic(&self) -> bool {
        self.aromatic.is_huckel_aromatic()
    }

    fn has_nitrogen_lone_pair_donor(&self, mol: &Molecule) -> bool {
        self.localized
            .as_ref()
            .is_some_and(|analysis| analysis.has_nitrogen_lone_pair_donor(mol))
    }

    fn localized_all_atoms_are_candidates(&self) -> bool {
        self.localized
            .as_ref()
            .is_some_and(AromaticRingDonorAnalysis::all_atoms_are_candidates)
    }

    fn localized_active_hetero_donor_count(&self, mol: &Molecule) -> usize {
        self.localized
            .as_ref()
            .map(|analysis| analysis.active_hetero_donor_count(mol))
            .unwrap_or(0)
    }

    fn localized_has_active_chalcogen_donor(&self, mol: &Molecule) -> bool {
        self.localized
            .as_ref()
            .is_some_and(|analysis| analysis.has_active_chalcogen_donor(mol))
    }

    fn localized_has_active_element_donor(&self, mol: &Molecule, symbol: &str) -> bool {
        self.localized
            .as_ref()
            .is_some_and(|analysis| analysis.has_active_element_donor(mol, symbol))
    }

    fn localized_has_element_candidate(&self, mol: &Molecule, symbol: &str) -> bool {
        self.localized
            .as_ref()
            .is_some_and(|analysis| analysis.has_element_candidate(mol, symbol))
    }

    fn localized_has_candidate_carbon_electronegative_exocyclic_pi_bond(
        &self,
        mol: &Molecule,
        ring: &Ring,
    ) -> bool {
        self.localized.as_ref().is_some_and(|analysis| {
            ring_has_candidate_carbon_electronegative_exocyclic_pi_bond(mol, ring, analysis)
        })
    }

    fn localized_has_active_anionic_nitrogen_donor(&self, mol: &Molecule) -> bool {
        self.localized
            .as_ref()
            .is_some_and(|analysis| analysis.has_active_anionic_nitrogen_donor(mol))
    }

    fn localized_atom_has_active_element_donor(
        &self,
        mol: &Molecule,
        atom_id: AtomId,
        symbol: &str,
    ) -> bool {
        self.localized
            .as_ref()
            .is_some_and(|analysis| analysis.atom_has_active_element_donor(mol, atom_id, symbol))
    }

    fn localized_atom_has_active_chalcogen_donor(&self, mol: &Molecule, atom_id: AtomId) -> bool {
        self.localized
            .as_ref()
            .is_some_and(|analysis| analysis.atom_has_active_chalcogen_donor(mol, atom_id))
    }

    fn localized_atom_is_chalcogen_candidate(&self, mol: &Molecule, atom_id: AtomId) -> bool {
        self.localized
            .as_ref()
            .is_some_and(|analysis| analysis.atom_is_chalcogen_candidate(mol, atom_id))
    }

    fn localized_atom_is_element_candidate(
        &self,
        mol: &Molecule,
        atom_id: AtomId,
        symbol: &str,
    ) -> bool {
        self.localized
            .as_ref()
            .is_some_and(|analysis| analysis.atom_is_element_candidate(mol, atom_id, symbol))
    }
}

impl AromaticRingDonorAnalysis {
    fn non_aromatic() -> Self {
        Self {
            atoms: Vec::new(),
            fixed_electron_count: Some(0),
        }
    }

    fn all_atoms_are_candidates(&self) -> bool {
        self.atoms
            .iter()
            .all(|atom| !matches!(atom.donor, AromaticElectronDonorType::None))
    }

    fn is_huckel_aromatic(&self) -> bool {
        self.electron_count()
            .is_some_and(|electrons| electrons >= 2 && (electrons - 2) % 4 == 0)
    }

    fn is_fused_huckel_aromatic(&self) -> bool {
        self.electron_count()
            .is_some_and(|electrons| electrons >= 6 && (electrons - 2) % 4 == 0)
    }

    fn electron_count(&self) -> Option<u8> {
        if let Some(electrons) = self.fixed_electron_count {
            return Some(electrons);
        }
        if !self.all_atoms_are_candidates() {
            return None;
        }
        let donors = self.atoms.iter().map(|atom| atom.donor).collect::<Vec<_>>();
        huckel_electron_count_for_donors(&donors).or_else(|| Some(self.max_electron_count()))
    }

    fn max_electron_count(&self) -> u8 {
        self.atoms
            .iter()
            .map(|atom| aromatic_donor_electron_range(atom.donor).1)
            .sum()
    }

    fn active_hetero_donor_count(&self, mol: &Molecule) -> usize {
        self.atoms
            .iter()
            .filter(|atom_donor| {
                mol.atom(atom_donor.atom).is_ok_and(|atom| {
                    matches!(atom.element.symbol(), "N" | "O" | "P" | "S" | "Se" | "Te")
                }) && aromatic_donor_electron_range(atom_donor.donor).1 > 0
            })
            .count()
    }

    fn has_nitrogen_lone_pair_donor(&self, mol: &Molecule) -> bool {
        self.atoms.iter().any(|atom_donor| {
            mol.atom(atom_donor.atom)
                .is_ok_and(|atom| atom.element.symbol() == "N")
                && aromatic_donor_electron_range(atom_donor.donor).1 >= 2
        })
    }

    fn has_active_chalcogen_donor(&self, mol: &Molecule) -> bool {
        self.atoms.iter().any(|atom_donor| {
            mol.atom(atom_donor.atom)
                .is_ok_and(|atom| matches!(atom.element.symbol(), "O" | "S" | "Se" | "Te"))
                && aromatic_donor_electron_range(atom_donor.donor).1 > 0
        })
    }

    fn has_active_element_donor(&self, mol: &Molecule, symbol: &str) -> bool {
        self.atoms.iter().any(|atom_donor| {
            mol.atom(atom_donor.atom)
                .is_ok_and(|atom| atom.element.symbol() == symbol)
                && aromatic_donor_electron_range(atom_donor.donor).1 > 0
        })
    }

    fn has_element_candidate(&self, mol: &Molecule, symbol: &str) -> bool {
        self.atoms.iter().any(|atom_donor| {
            mol.atom(atom_donor.atom)
                .is_ok_and(|atom| atom.element.symbol() == symbol)
                && !matches!(atom_donor.donor, AromaticElectronDonorType::None)
        })
    }

    fn has_active_anionic_nitrogen_donor(&self, mol: &Molecule) -> bool {
        self.atoms.iter().any(|atom_donor| {
            mol.atom(atom_donor.atom)
                .is_ok_and(|atom| atom.element.symbol() == "N" && atom.formal_charge < 0)
                && aromatic_donor_electron_range(atom_donor.donor).1 > 0
        })
    }

    fn atom_has_active_element_donor(&self, mol: &Molecule, atom_id: AtomId, symbol: &str) -> bool {
        self.atoms.iter().any(|atom_donor| {
            atom_donor.atom == atom_id
                && mol
                    .atom(atom_id)
                    .is_ok_and(|atom| atom.element.symbol() == symbol)
                && aromatic_donor_electron_range(atom_donor.donor).1 > 0
        })
    }

    fn atom_has_active_chalcogen_donor(&self, mol: &Molecule, atom_id: AtomId) -> bool {
        self.atoms.iter().any(|atom_donor| {
            atom_donor.atom == atom_id
                && mol
                    .atom(atom_id)
                    .is_ok_and(|atom| matches!(atom.element.symbol(), "O" | "S" | "Se" | "Te"))
                && aromatic_donor_electron_range(atom_donor.donor).1 > 0
        })
    }

    fn atom_is_chalcogen_candidate(&self, mol: &Molecule, atom_id: AtomId) -> bool {
        self.atoms.iter().any(|atom_donor| {
            atom_donor.atom == atom_id
                && mol
                    .atom(atom_id)
                    .is_ok_and(|atom| matches!(atom.element.symbol(), "O" | "S" | "Se" | "Te"))
                && !matches!(atom_donor.donor, AromaticElectronDonorType::None)
        })
    }

    fn atom_is_heavy_chalcogen_candidate(&self, mol: &Molecule, atom_id: AtomId) -> bool {
        self.atoms.iter().any(|atom_donor| {
            atom_donor.atom == atom_id
                && mol
                    .atom(atom_id)
                    .is_ok_and(|atom| matches!(atom.element.symbol(), "S" | "Se" | "Te"))
                && !matches!(atom_donor.donor, AromaticElectronDonorType::None)
        })
    }

    fn atom_is_element_candidate(&self, mol: &Molecule, atom_id: AtomId, symbol: &str) -> bool {
        self.atoms.iter().any(|atom_donor| {
            atom_donor.atom == atom_id
                && mol
                    .atom(atom_id)
                    .is_ok_and(|atom| atom.element.symbol() == symbol)
                && !matches!(atom_donor.donor, AromaticElectronDonorType::None)
        })
    }
}

fn localized_ring_atom_donor_type(
    mol: &Molecule,
    ring: &Ring,
    atom_id: AtomId,
) -> std::result::Result<AromaticElectronDonorType, AromaticityError> {
    let atom = mol.atom(atom_id).expect("ring atom should be live");
    if !aromaticity_supported_element(atom) {
        return if ring_has_aromatic_order(mol, ring) {
            Err(AromaticityError::UnsupportedElement(atom_id))
        } else {
            Ok(AromaticElectronDonorType::None)
        };
    }
    let donor = if ring_atom_has_pi_bond(mol, ring, atom_id) {
        AromaticElectronDonorType::One
    } else {
        rdkit_like_atom_donor_type(mol, ring, atom_id, atom, true)
    };
    if !atom_is_rdkit_aromatic_candidate_for_donor(
        mol,
        atom_id,
        atom,
        donor,
        RdkitAromaticCandidateOptions::default(),
    ) {
        return Ok(AromaticElectronDonorType::None);
    }

    Ok(donor)
}

fn aromaticity_supported_element(atom: &Atom) -> bool {
    matches!(
        atom.element.symbol(),
        "B" | "C" | "N" | "O" | "P" | "S" | "Se" | "Te"
    )
}

fn atom_is_rdkit_aromatic_candidate(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> bool {
    atom_is_rdkit_aromatic_candidate_with_options(
        mol,
        atom_id,
        atom,
        RdkitAromaticCandidateOptions::default(),
    )
}

fn atom_is_rdkit_aromatic_candidate_with_options(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    options: RdkitAromaticCandidateOptions,
) -> bool {
    atom_is_rdkit_aromatic_candidate_for_donor(
        mol,
        atom_id,
        atom,
        AromaticElectronDonorType::Any,
        options,
    )
}

fn atom_is_rdkit_aromatic_candidate_for_donor(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    donor: AromaticElectronDonorType,
    options: RdkitAromaticCandidateOptions,
) -> bool {
    if matches!(donor, AromaticElectronDonorType::None) {
        return false;
    }
    let atomic_number = atom.element.atomic_number();
    if options.only_carbon_or_nitrogen && !matches!(atomic_number, 6 | 7) {
        return false;
    }
    if !options.allow_third_row && atomic_number > 10 {
        return false;
    }
    if atomic_number > 18 && (!options.allow_higher_exceptions || !matches!(atomic_number, 34 | 52))
    {
        return false;
    }
    if atom_aromatic_candidate_degree(mol, atom_id, atom) > 3 {
        return false;
    }
    let Some(default_valence) = rdkit_default_valence(atom) else {
        return false;
    };
    let Some(charge_adjusted_default_valence) = rdkit_charge_adjusted_default_valence(atom) else {
        return false;
    };
    if default_valence > 0
        && atom_rdkit_aromatic_total_valence(mol, atom_id, atom) > charge_adjusted_default_valence
    {
        return false;
    }
    if atom_explicit_pi_bond_count(mol, atom_id) > 1 {
        return false;
    }
    if !options.allow_triple_bonds && atom_has_explicit_triple_bond(mol, atom_id) {
        return false;
    }
    if !options.allow_exocyclic_multiple_bonds && atom_has_non_ring_multiple_bond(mol, atom_id) {
        return false;
    }
    atom_passes_rdkit_aromatic_radical_eligibility(atom)
}

fn atom_passes_rdkit_aromatic_radical_eligibility(atom: &Atom) -> bool {
    let radical_electrons = atom.radical.map_or(0, AtomRadical::unpaired_electron_count);
    radical_electrons == 0 || atom.element.symbol() == "C" && atom.formal_charge == 0
}

fn atom_is_rdkit_aromatic_element_candidate(mol: &Molecule, atom_id: AtomId, symbol: &str) -> bool {
    mol.atom(atom_id).is_ok_and(|atom| {
        atom.element.symbol() == symbol && atom_is_rdkit_aromatic_candidate(mol, atom_id, atom)
    })
}

fn atom_aromatic_candidate_degree(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> u8 {
    let bonded_degree = mol
        .incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .filter(|(_, bond)| !matches!(bond.order, BondOrder::Zero | BondOrder::Dative))
        .count()
        .min(u8::MAX as usize) as u8;
    bonded_degree
        .saturating_add(atom.explicit_hydrogens)
        .saturating_add(aromaticity_implicit_hydrogen_count(mol, atom_id, atom))
}

fn aromaticity_implicit_hydrogen_count(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> u8 {
    if atom.no_implicit_hydrogens {
        return atom.implicit_hydrogens.unwrap_or(0);
    }
    if let Some(hydrogens) = atom.implicit_hydrogens {
        return hydrogens;
    }
    let Some(target) = aromaticity_valence_target(atom) else {
        return 0;
    };
    target.saturating_sub(explicit_valence(mol, atom_id).saturating_add(atom.explicit_hydrogens))
}

fn atom_rdkit_aromatic_total_valence(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> u8 {
    explicit_valence(mol, atom_id)
        .saturating_add(atom.explicit_hydrogens)
        .saturating_add(aromaticity_implicit_hydrogen_count(mol, atom_id, atom))
}

fn aromaticity_valence_target(atom: &Atom) -> Option<u8> {
    if atom.aromatic {
        return match atom.element.symbol() {
            "B" | "C" => Some(3),
            "N" => {
                if atom.explicit_hydrogens > 0 {
                    Some(3)
                } else {
                    Some(2)
                }
            }
            "O" | "S" | "Se" | "Te" => Some(2),
            "P" => Some(3),
            _ => None,
        };
    }

    match (atom.element.symbol(), atom.formal_charge) {
        ("B", -1) => Some(4),
        ("B", _) => Some(3),
        ("C", 1 | -1) => Some(3),
        ("C", _) => Some(4),
        ("N", 1) => Some(4),
        ("N", -1) => Some(2),
        ("N", _) => Some(3),
        ("O", -1) => Some(1),
        ("O", 1) => Some(3),
        ("O", _) => Some(2),
        ("P", 1) => Some(4),
        ("P", _) => Some(3),
        ("S" | "Se" | "Te", -1 | 1) => Some(1),
        ("S" | "Se" | "Te", _) => Some(2),
        _ => None,
    }
}

fn aromatic_order_ring_donor_analysis(
    mol: &Molecule,
    ring: &Ring,
    localized: &AromaticRingDonorAnalysis,
) -> std::result::Result<AromaticRingDonorAnalysis, AromaticityError> {
    let exocyclic_pi_steals_electrons =
        aromatic_order_ring_allows_exocyclic_carbon_zero_contribution(mol, ring, localized);
    let mut donors = Vec::with_capacity(ring.atoms.len());
    for atom_id in &ring.atoms {
        let atom = mol.atom(*atom_id).expect("ring atom should be live");
        let donor = aromatic_order_atom_donor_type(
            mol,
            ring,
            *atom_id,
            atom,
            exocyclic_pi_steals_electrons,
        )?;
        donors.push(donor);
    }
    Ok(AromaticRingDonorAnalysis {
        atoms: ring
            .atoms
            .iter()
            .copied()
            .zip(donors)
            .map(|(atom, donor)| AromaticRingAtomDonor { atom, donor })
            .collect(),
        fixed_electron_count: None,
    })
}

fn aromatic_order_atom_donor_type(
    mol: &Molecule,
    ring: &Ring,
    atom_id: AtomId,
    atom: &Atom,
    exocyclic_pi_steals_electrons: bool,
) -> std::result::Result<AromaticElectronDonorType, AromaticityError> {
    if !aromaticity_supported_element(atom) {
        return Err(AromaticityError::UnsupportedElement(atom_id));
    }
    let donor = if exocyclic_pi_steals_electrons
        && atom_has_electronegative_exocyclic_pi_bond(mol, ring, atom_id, atom)
    {
        AromaticElectronDonorType::Vacant
    } else {
        match atom.element.symbol() {
            "B" | "C" => AromaticElectronDonorType::One,
            "N" => {
                if atom.explicit_hydrogens > 0
                    || atom.formal_charge == 0
                        && aromatic_order_nitrogen_is_pyrrole_like(mol, atom_id)
                {
                    AromaticElectronDonorType::OneOrTwo
                } else {
                    AromaticElectronDonorType::One
                }
            }
            "O" | "S" | "Se" | "Te" | "P" => AromaticElectronDonorType::Two,
            _ => return Err(AromaticityError::UnsupportedElement(atom_id)),
        }
    };
    if !atom_is_rdkit_aromatic_candidate_for_donor(
        mol,
        atom_id,
        atom,
        donor,
        RdkitAromaticCandidateOptions::default(),
    ) {
        return Ok(AromaticElectronDonorType::None);
    }
    Ok(donor)
}

fn rdkit_like_atom_donor_type(
    mol: &Molecule,
    ring: &Ring,
    atom_id: AtomId,
    atom: &Atom,
    exocyclic_bonds_steal_electrons: bool,
) -> AromaticElectronDonorType {
    let Some(mut electrons) = count_rdkit_like_atom_pi_electrons(mol, atom_id, atom) else {
        return AromaticElectronDonorType::None;
    };
    let exocyclic_pi_neighbor = atom_exocyclic_pi_neighbor(mol, ring, atom_id);
    let has_exocyclic_pi = exocyclic_pi_neighbor.is_some();
    let has_incident_pi_bond = atom_explicit_pi_bond_count(mol, atom_id) > 0
        || mol
            .incident_bonds(atom_id)
            .ok()
            .into_iter()
            .flatten()
            .any(|(_, bond)| matches!(bond.order, BondOrder::Aromatic));

    if electrons == 0 {
        if has_exocyclic_pi {
            AromaticElectronDonorType::Vacant
        } else if ring_atom_has_pi_bond(mol, ring, atom_id) {
            AromaticElectronDonorType::One
        } else {
            AromaticElectronDonorType::None
        }
    } else if electrons == 1 {
        if exocyclic_pi_neighbor.is_some_and(|neighbor| {
            exocyclic_bonds_steal_electrons
                && atom_is_more_electronegative_than(mol, neighbor, atom)
        }) {
            AromaticElectronDonorType::Vacant
        } else if has_exocyclic_pi || has_incident_pi_bond {
            AromaticElectronDonorType::One
        } else if atom.formal_charge == 1 {
            AromaticElectronDonorType::Vacant
        } else {
            AromaticElectronDonorType::None
        }
    } else {
        if exocyclic_pi_neighbor.is_some_and(|neighbor| {
            exocyclic_bonds_steal_electrons
                && atom_is_more_electronegative_than(mol, neighbor, atom)
        }) {
            electrons -= 1;
        }
        if electrons % 2 == 1 {
            AromaticElectronDonorType::One
        } else {
            AromaticElectronDonorType::Two
        }
    }
}

fn count_rdkit_like_atom_pi_electrons(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> Option<u8> {
    let default_valence = rdkit_default_valence(atom)?;
    let degree = atom_aromatic_candidate_degree(mol, atom_id, atom);
    if default_valence <= 1 || degree > 3 {
        return None;
    }

    let lone_pair_electrons = (i16::from(rdkit_outer_electrons(atom)?)
        - i16::from(default_valence)
        - i16::from(atom.formal_charge))
    .max(0);
    let radical_electrons = i16::from(atom.radical.map_or(0, AtomRadical::unpaired_electron_count));
    let mut electrons =
        i16::from(default_valence) - i16::from(degree) + lone_pair_electrons - radical_electrons;
    if electrons < 0 {
        return None;
    }
    if electrons > 1 && atom_explicit_unsaturation(mol, atom_id) > 1 {
        electrons = 1;
    }
    u8::try_from(electrons).ok()
}

fn rdkit_default_valence(atom: &Atom) -> Option<u8> {
    rdkit_default_valence_for_element(atom.element)
}

fn rdkit_charge_adjusted_default_valence(atom: &Atom) -> Option<u8> {
    let atomic_number = i16::from(atom.element.atomic_number()) - i16::from(atom.formal_charge);
    let atomic_number = u8::try_from(atomic_number).ok()?;
    rdkit_default_valence_for_element(Element::from_atomic_number(atomic_number)?)
}

fn rdkit_default_valence_for_element(element: Element) -> Option<u8> {
    match element.symbol() {
        "B" => Some(3),
        "C" => Some(4),
        "N" | "P" => Some(3),
        "O" | "S" | "Se" | "Te" => Some(2),
        _ => None,
    }
}

fn rdkit_outer_electrons(atom: &Atom) -> Option<u8> {
    match atom.element.symbol() {
        "B" => Some(3),
        "C" => Some(4),
        "N" | "P" => Some(5),
        "O" | "S" | "Se" | "Te" => Some(6),
        _ => None,
    }
}

fn aromatic_donor_electron_range(donor: AromaticElectronDonorType) -> (u8, u8) {
    match donor {
        AromaticElectronDonorType::Vacant | AromaticElectronDonorType::None => (0, 0),
        AromaticElectronDonorType::One => (1, 1),
        AromaticElectronDonorType::Two => (2, 2),
        AromaticElectronDonorType::OneOrTwo => (1, 2),
        AromaticElectronDonorType::Any => (0, 2),
    }
}

fn huckel_electron_count_for_donors(donors: &[AromaticElectronDonorType]) -> Option<u8> {
    if donors
        .iter()
        .filter(|donor| matches!(donor, AromaticElectronDonorType::Any))
        .count()
        > 1
    {
        return None;
    }

    let min_electrons = donors
        .iter()
        .map(|donor| aromatic_donor_electron_range(*donor).0)
        .sum();
    let max_electrons = donors
        .iter()
        .map(|donor| aromatic_donor_electron_range(*donor).1)
        .sum();
    huckel_electron_count_in_range(min_electrons, max_electrons)
}

fn aromatic_order_ring_allows_exocyclic_carbon_zero_contribution(
    mol: &Molecule,
    ring: &Ring,
    localized: &AromaticRingDonorAnalysis,
) -> bool {
    if !localized.has_active_chalcogen_donor(mol)
        || !ring_has_candidate_carbon_electronegative_exocyclic_pi_bond(mol, ring, localized)
    {
        return false;
    }

    let terminal_exocyclic_pi_count = ring_terminal_exocyclic_pi_bond_count(mol, ring);
    ring.atoms.len() == 5
        || ring.atoms.len() == 6
            && (terminal_exocyclic_pi_count == 1
                || localized.has_active_element_donor(mol, "N") && terminal_exocyclic_pi_count >= 2)
}

fn atom_has_electronegative_exocyclic_pi_bond(
    mol: &Molecule,
    ring: &Ring,
    atom_id: AtomId,
    atom: &Atom,
) -> bool {
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .any(|(bond_id, bond)| {
            if ring.bonds.contains(&bond_id) || !matches!(bond.order, BondOrder::Double) {
                return false;
            }
            let other_id = bond.other_atom(atom_id);
            if ring.atoms.contains(&other_id) {
                return false;
            }
            atom_is_more_electronegative_than(mol, other_id, atom)
        })
}

fn huckel_electron_count_in_range(min_electrons: u8, max_electrons: u8) -> Option<u8> {
    if max_electrons == 2 {
        return Some(2);
    }
    if max_electrons < 6 {
        return None;
    }
    (min_electrons..=max_electrons)
        .filter(|electrons| *electrons >= 6)
        .find(|electrons| (electrons - 2) % 4 == 0)
}

fn aromatic_order_nitrogen_is_pyrrole_like(mol: &Molecule, atom_id: AtomId) -> bool {
    mol.incident_bonds(atom_id)
        .map(|bonds| bonds.count() >= 3)
        .unwrap_or(false)
}

#[cfg(test)]
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

fn atom_explicit_pi_bond_count(mol: &Molecule, atom_id: AtomId) -> usize {
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .filter(|(_, bond)| matches!(bond.order, BondOrder::Double | BondOrder::Triple))
        .count()
}

fn atom_has_explicit_triple_bond(mol: &Molecule, atom_id: AtomId) -> bool {
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .any(|(_, bond)| matches!(bond.order, BondOrder::Triple))
}

fn atom_has_non_ring_multiple_bond(mol: &Molecule, atom_id: AtomId) -> bool {
    let computed_membership;
    let membership = if let Some(membership) = mol.ring_membership() {
        membership
    } else {
        computed_membership = super::rings::compute_ring_membership(mol).0;
        &computed_membership
    };
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .any(|(bond_id, bond)| {
            matches!(bond.order, BondOrder::Double | BondOrder::Triple)
                && !membership.bond_in_ring(bond_id)
        })
}

fn atom_explicit_unsaturation(mol: &Molecule, atom_id: AtomId) -> u8 {
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .map(|(_, bond)| match bond.order {
            BondOrder::Double => 1,
            BondOrder::Triple => 2,
            BondOrder::Quadruple => 3,
            _ => 0,
        })
        .sum()
}

fn atom_exocyclic_pi_neighbor(mol: &Molecule, ring: &Ring, atom_id: AtomId) -> Option<AtomId> {
    mol.incident_bonds(atom_id)
        .ok()?
        .filter_map(|(bond_id, bond)| {
            if ring.bonds.contains(&bond_id)
                || ring.atoms.contains(&bond.other_atom(atom_id))
                || !matches!(bond.order, BondOrder::Double | BondOrder::Triple)
            {
                return None;
            }
            Some(bond.other_atom(atom_id))
        })
        .next()
}

fn atom_is_more_electronegative_than(mol: &Molecule, left: AtomId, right: &Atom) -> bool {
    mol.atom(left).is_ok_and(|left| {
        atom_electronegativity(left)
            .zip(atom_electronegativity(right))
            .is_some_and(|(left, right)| left > right)
    })
}

fn atom_electronegativity(atom: &Atom) -> Option<u16> {
    match atom.element.symbol() {
        "B" => Some(204),
        "C" | "Se" => Some(255),
        "N" => Some(304),
        "O" => Some(344),
        "P" | "Te" => Some(219),
        "S" => Some(258),
        _ => None,
    }
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn localized_ring_analysis_feeds_simple_huckel_path() {
        let mut mol = Molecule::new();
        let atoms = (0..6)
            .map(|_| mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element"))))
            .collect::<Vec<_>>();
        let bonds = (0..6)
            .map(|index| {
                let order = if index % 2 == 0 {
                    BondOrder::Double
                } else {
                    BondOrder::Single
                };
                mol.add_bond(atoms[index], atoms[(index + 1) % atoms.len()], order)
                    .expect("ring bond")
            })
            .collect::<Vec<_>>();
        let ring = Ring { atoms, bonds };

        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");
        let localized = analysis.localized.as_ref().expect("localized analysis");
        assert_eq!(&analysis.aromatic, localized);
        assert!(analysis.is_huckel_aromatic());
    }

    #[test]
    fn candidate_options_can_disallow_exocyclic_multiple_bonds() {
        let mut mol = Molecule::new();
        let carbonyl_carbon =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_e = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let oxygen = mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(carbonyl_carbon, carbon_b, BondOrder::Single)
            .expect("ring bond");
        mol.add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        mol.add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        mol.add_bond(carbon_d, carbon_e, BondOrder::Single)
            .expect("ring bond");
        mol.add_bond(carbon_e, carbonyl_carbon, BondOrder::Single)
            .expect("ring bond");
        let carbonyl_bond = mol
            .add_bond(carbonyl_carbon, oxygen, BondOrder::Double)
            .expect("carbonyl bond");
        let membership = perceive_ring_membership(&mut mol);
        let atom = mol.atom(carbonyl_carbon).expect("carbonyl carbon");

        assert!(!membership.bond_in_ring(carbonyl_bond));
        assert!(atom_is_rdkit_aromatic_candidate(
            &mol,
            carbonyl_carbon,
            atom
        ));
        assert!(!atom_is_rdkit_aromatic_candidate_for_donor(
            &mol,
            carbonyl_carbon,
            atom,
            AromaticElectronDonorType::None,
            RdkitAromaticCandidateOptions::default()
        ));
        assert!(atom_is_rdkit_aromatic_candidate_for_donor(
            &mol,
            carbonyl_carbon,
            atom,
            AromaticElectronDonorType::One,
            RdkitAromaticCandidateOptions::default()
        ));
        assert!(!atom_is_rdkit_aromatic_candidate_for_donor(
            &mol,
            carbonyl_carbon,
            atom,
            AromaticElectronDonorType::One,
            RdkitAromaticCandidateOptions {
                allow_exocyclic_multiple_bonds: false,
                ..RdkitAromaticCandidateOptions::default()
            }
        ));
    }

    #[test]
    fn fused_candidate_does_not_reject_seven_member_chalcogen_by_shape() {
        let mut mol = Molecule::new();
        let sulfur = mol.add_atom(aromatic_atom("S"));
        let carbon_a = mol.add_atom(aromatic_carbon());
        let carbon_b = mol.add_atom(aromatic_carbon());
        let carbon_c = mol.add_atom(aromatic_carbon());
        let carbon_d = mol.add_atom(aromatic_carbon());
        let carbon_e = mol.add_atom(aromatic_carbon());
        let carbon_f = mol.add_atom(aromatic_carbon());
        let bond_a = mol
            .add_bond(sulfur, carbon_a, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_d, carbon_e, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_e, carbon_f, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_g = mol
            .add_bond(carbon_f, sulfur, BondOrder::Aromatic)
            .expect("ring bond");
        let ring = Ring {
            atoms: vec![
                sulfur, carbon_a, carbon_b, carbon_c, carbon_d, carbon_e, carbon_f,
            ],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f, bond_g],
        };

        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("donor analysis");
        assert!(analysis.localized_all_atoms_are_candidates());
        assert!(analysis.localized_has_active_chalcogen_donor(&mol));
        assert_eq!(ring.atoms.len(), 7);
        assert!(aromatic_fused_candidate_from_analysis(
            &mol, &ring, &analysis
        ));
    }

    #[test]
    fn non_aromatic_fusion_protection_uses_active_donor_state() {
        let mut mol = Molecule::new();
        let nitrogen_a = mol.add_atom(cation_with_one_implicit_hydrogen("N"));
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let nitrogen_b = mol.add_atom(cation_with_one_implicit_hydrogen("N"));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(nitrogen_a, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_b, nitrogen_b, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(nitrogen_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_d, nitrogen_a, BondOrder::Single)
            .expect("ring bond");
        let oxygen_a = mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        let oxygen_b = mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(carbon_a, oxygen_a, BondOrder::Double)
            .expect("terminal carbonyl");
        mol.add_bond(carbon_c, oxygen_b, BondOrder::Double)
            .expect("terminal carbonyl");
        let ring = Ring {
            atoms: vec![
                nitrogen_a, carbon_a, carbon_b, nitrogen_b, carbon_c, carbon_d,
            ],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };

        assert_eq!(ring_hetero_donor_count(&mol, &ring), 2);
        assert_eq!(ring_terminal_exocyclic_pi_bond_count(&mol, &ring), 2);
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");
        assert_eq!(analysis.localized_active_hetero_donor_count(&mol), 0);
        assert!(!ring_protects_non_aromatic_fusion_single(
            &mol, &ring, &analysis, false, 2
        ));
        let rings = vec![ring];
        let analyses = vec![analysis];
        assert!(
            atoms_in_nitrogen_or_terminal_pi_free_rings(&mol, &rings, &analyses, &[0]).is_empty()
        );
    }

    #[test]
    fn non_aromatic_fusion_protection_requires_candidate_carbon_state() {
        let mut mol = Molecule::new();
        let mut disqualified_carbon = aromatic_carbon();
        disqualified_carbon.formal_charge = 1;
        disqualified_carbon.radical = Some(AtomRadical::Doublet);
        let carbon_a = mol.add_atom(disqualified_carbon);
        let carbon_b = mol.add_atom(aromatic_carbon());
        let carbon_c = mol.add_atom(aromatic_carbon());
        let carbon_d = mol.add_atom(aromatic_carbon());
        let carbon_e = mol.add_atom(aromatic_carbon());
        let bond_a = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_d, carbon_e, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_e, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let ring = Ring {
            atoms: vec![carbon_a, carbon_b, carbon_c, carbon_d, carbon_e],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");

        assert!(fused_component_is_all_carbon(&mol, &ring));
        assert!(!ring_is_all_candidate_carbon(&mol, &ring, &analysis));
        let containing_rings = vec![(0usize, &ring)];
        let analyses = vec![analysis.clone()];
        assert!(!bond_has_candidate_carbon_endpoints_in_containing_ring(
            &mol,
            bond_a,
            &containing_rings,
            &analyses,
        ));
        assert!(ring_protects_non_aromatic_fusion_single(
            &mol, &ring, &analysis, false, 2
        ));
    }

    #[test]
    fn fused_carbon_nitrogen_fallback_requires_candidate_nitrogen_state() {
        let mut mol = Molecule::new();
        let mut nitrogen_atom = Atom::new(Element::from_symbol("N").expect("test element"));
        nitrogen_atom.explicit_hydrogens = 2;
        nitrogen_atom.no_implicit_hydrogens = true;
        let nitrogen = mol.add_atom(nitrogen_atom);
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_e = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_d, carbon_e, BondOrder::Single)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_e, nitrogen, BondOrder::Single)
            .expect("ring bond");
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, carbon_b, carbon_c, carbon_d, carbon_e],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");
        let rings = vec![ring.clone()];
        let analyses = vec![analysis];

        assert!(!fused_component_is_carbon_or_candidate_nitrogen(
            &mol,
            &ring,
            &analyses[0]
        ));
        assert!(!analyses[0].localized_atom_is_element_candidate(&mol, nitrogen, "N"));
        assert!(!fused_component_is_carbon_with_candidate_nitrogen(
            &mol,
            &ring,
            &rings,
            &analyses,
            &[0]
        ));
    }

    #[test]
    fn fused_carbon_nitrogen_fallback_requires_candidate_carbon_state() {
        let mut mol = Molecule::new();
        let mut carbon_atom = aromatic_carbon();
        carbon_atom.formal_charge = 1;
        carbon_atom.radical = Some(AtomRadical::Doublet);
        let carbon = mol.add_atom(carbon_atom);
        let nitrogen = mol.add_atom(aromatic_atom("N"));
        let carbon_a = mol.add_atom(aromatic_carbon());
        let carbon_b = mol.add_atom(aromatic_carbon());
        let carbon_c = mol.add_atom(aromatic_carbon());
        let carbon_d = mol.add_atom(aromatic_carbon());
        let bond_a = mol
            .add_bond(carbon, nitrogen, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_d, carbon, BondOrder::Aromatic)
            .expect("ring bond");
        let ring = Ring {
            atoms: vec![carbon, nitrogen, carbon_a, carbon_b, carbon_c, carbon_d],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");
        let rings = vec![ring.clone()];
        let analyses = vec![analysis];

        assert!(mol
            .atom(carbon)
            .is_ok_and(|atom| atom.element.symbol() == "C"));
        assert!(!analyses[0].localized_atom_is_element_candidate(&mol, carbon, "C"));
        assert!(analyses[0].localized_atom_is_element_candidate(&mol, nitrogen, "N"));
        assert!(!fused_component_is_carbon_or_candidate_nitrogen(
            &mol,
            &ring,
            &analyses[0]
        ));
        assert!(!fused_component_is_carbon_with_candidate_nitrogen(
            &mol,
            &ring,
            &rings,
            &analyses,
            &[0]
        ));
    }

    #[test]
    fn fused_exocyclic_carbon_component_uses_donor_analysis() {
        let mut mol = Molecule::new();
        let carbon_a = mol.add_atom(aromatic_carbon());
        let carbon_b = mol.add_atom(aromatic_carbon());
        let carbon_c = mol.add_atom(aromatic_carbon());
        let carbon_d = mol.add_atom(aromatic_carbon());
        let carbon_e = mol.add_atom(aromatic_carbon());
        let carbon_f = mol.add_atom(aromatic_carbon());
        let bond_a = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_d, carbon_e, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_e, carbon_f, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_f, carbon_a, BondOrder::Aromatic)
            .expect("ring bond");
        let oxygen = mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(carbon_c, oxygen, BondOrder::Double)
            .expect("terminal carbonyl");
        let ring = Ring {
            atoms: vec![carbon_a, carbon_b, carbon_c, carbon_d, carbon_e, carbon_f],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };

        let analysis = localized_ring_donor_analysis(&mol, &ring).expect("fused donor analysis");

        assert!(fused_component_is_all_carbon(&mol, &ring));
        assert_eq!(ring_exocyclic_pi_bond_count(&mol, &ring), 1);
        assert_eq!(ring_pi_bond_count(&mol, &ring), 6);
        assert_eq!(analysis.fixed_electron_count, None);
        assert_eq!(analysis.atoms.len(), ring.atoms.len());
        assert!(analysis.is_fused_huckel_aromatic());
    }

    #[test]
    fn fused_component_search_uses_member_ring_size_not_component_atom_cap() {
        let mut mol = Molecule::new();
        let mut rings = Vec::new();
        let mut shared = None;
        for _ in 0..12 {
            let (ring, next_left, next_right, next_bond) = fused_carbon_ring(&mut mol, shared, 6);
            rings.push(ring);
            shared = Some((next_left, next_right, next_bond));
        }
        let analyses = rings
            .iter()
            .map(|ring| RingAromaticityAnalysis::new(&mol, ring).expect("ring analysis"))
            .collect::<Vec<_>>();
        let indexes = (0..rings.len()).collect::<Vec<_>>();
        let component = fused_component_ring(&rings, &indexes);

        assert!(component.atoms.len() > 48);
        assert!(rings
            .iter()
            .all(|ring| ring.atoms.len() <= MAX_FUSED_AROMATIC_RING_SIZE));
        assert!(rings.iter().enumerate().all(|(index, ring)| {
            aromatic_fused_candidate_from_analysis(&mol, ring, &analyses[index])
        }));
        let subsets = aromatic_fused_ring_subsets(&mol, &rings, &indexes).expect("fused subsets");
        assert!(!subsets.is_empty(), "expected at least one accepted subset");
        for bond_id in &component.bonds {
            mol.bonds[bond_id.index()]
                .as_mut()
                .expect("component bond")
                .aromatic = false;
        }
        assert_eq!(
            component
                .bonds
                .iter()
                .filter(|bond_id| mol.bond(**bond_id).expect("component bond").aromatic)
                .count(),
            0
        );
        perceive_fused_aromatic_components(&mut mol, &rings, &analyses, &BTreeSet::new())
            .expect("fused aromaticity");

        let aromatic_bond_count = component
            .bonds
            .iter()
            .filter(|bond_id| mol.bond(**bond_id).expect("component bond").aromatic)
            .count();
        assert!(
            aromatic_bond_count > 0,
            "aromatic_bond_count={aromatic_bond_count}"
        );
    }

    #[test]
    fn all_carbon_exocyclic_fallback_requires_candidate_carbon_state() {
        let mut mol = Molecule::new();
        let mut carbon_atom = aromatic_carbon();
        carbon_atom.formal_charge = 1;
        carbon_atom.radical = Some(AtomRadical::Doublet);
        let carbon_a = mol.add_atom(carbon_atom);
        let carbon_b = mol.add_atom(aromatic_carbon());
        let carbon_c = mol.add_atom(aromatic_carbon());
        let carbon_d = mol.add_atom(aromatic_carbon());
        let carbon_e = mol.add_atom(aromatic_carbon());
        let carbon_f = mol.add_atom(aromatic_carbon());
        let bond_a = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_d, carbon_e, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_e, carbon_f, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_f, carbon_a, BondOrder::Aromatic)
            .expect("ring bond");
        let oxygen = mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(carbon_c, oxygen, BondOrder::Double)
            .expect("terminal carbonyl");
        let ring = Ring {
            atoms: vec![carbon_a, carbon_b, carbon_c, carbon_d, carbon_e, carbon_f],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let rings = vec![ring.clone()];
        let analyses = vec![RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis")];

        assert!(fused_component_is_all_carbon(&mol, &ring));
        assert_eq!(ring_exocyclic_pi_bond_count(&mol, &ring), 1);
        assert!(!analyses[0].localized_atom_is_element_candidate(&mol, carbon_a, "C"));
        assert!(!fused_component_is_all_candidate_carbon(
            &mol,
            &ring,
            &rings,
            &analyses,
            &[0]
        ));
    }

    #[test]
    fn fused_carbonyl_bridge_cleanup_requires_candidate_carbon_state() {
        let mut mol = Molecule::new();
        let mut disqualified_carbon = aromatic_carbon();
        disqualified_carbon.formal_charge = 1;
        disqualified_carbon.radical = Some(AtomRadical::Doublet);
        let carbon_a = mol.add_atom(disqualified_carbon);
        let carbon_b = mol.add_atom(aromatic_carbon());
        let carbon_c = mol.add_atom(aromatic_carbon());
        let carbon_d = mol.add_atom(aromatic_carbon());
        let carbon_e = mol.add_atom(aromatic_carbon());
        let carbon_f = mol.add_atom(aromatic_carbon());
        let bond_a = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_d, carbon_e, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_e, carbon_f, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_f, carbon_a, BondOrder::Aromatic)
            .expect("ring bond");
        let oxygen_a = mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        let oxygen_b = mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(carbon_c, oxygen_a, BondOrder::Double)
            .expect("terminal carbonyl");
        mol.add_bond(carbon_e, oxygen_b, BondOrder::Double)
            .expect("terminal carbonyl");

        let carbon_g = mol.add_atom(aromatic_carbon());
        let carbon_h = mol.add_atom(aromatic_carbon());
        let bond_g = mol
            .add_bond(carbon_b, carbon_g, BondOrder::Aromatic)
            .expect("fused ring bond");
        let bond_h = mol
            .add_bond(carbon_g, carbon_h, BondOrder::Aromatic)
            .expect("fused ring bond");
        let bond_i = mol
            .add_bond(carbon_h, carbon_a, BondOrder::Aromatic)
            .expect("fused ring bond");

        let carbonyl_ring = Ring {
            atoms: vec![carbon_a, carbon_b, carbon_c, carbon_d, carbon_e, carbon_f],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let fused_ring = Ring {
            atoms: vec![carbon_a, carbon_b, carbon_g, carbon_h],
            bonds: vec![bond_a, bond_g, bond_h, bond_i],
        };
        let rings = vec![carbonyl_ring.clone(), fused_ring];
        let analyses = rings
            .iter()
            .map(|ring| RingAromaticityAnalysis::new(&mol, ring).expect("ring analysis"))
            .collect::<Vec<_>>();

        assert!(fused_component_is_all_carbon(&mol, &carbonyl_ring));
        assert!(!ring_is_all_candidate_carbon(
            &mol,
            &carbonyl_ring,
            &analyses[0]
        ));
        assert_eq!(
            ring_terminal_exocyclic_pi_bond_count(&mol, &carbonyl_ring),
            2
        );

        clear_fused_carbonyl_bridge_atoms(&mut mol, &rings, &analyses);

        assert!(mol.atom(carbon_c).expect("carbonyl carbon").aromatic);
        assert!(mol.atom(carbon_d).expect("bridge carbon").aromatic);
        assert!(mol.atom(carbon_e).expect("carbonyl carbon").aromatic);
    }

    #[test]
    fn saturated_fused_ether_cleanup_uses_active_oxygen_donor_state() {
        let mut mol = Molecule::new();
        let mut oxygen = Atom::new(Element::from_symbol("O").expect("test element"));
        oxygen.explicit_hydrogens = 2;
        oxygen.no_implicit_hydrogens = true;
        oxygen.aromatic = true;
        let oxygen = mol.add_atom(oxygen);
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(oxygen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_d, oxygen, BondOrder::Single)
            .expect("ring bond");
        let carbon_e = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_f = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_g = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_f = mol
            .add_bond(carbon_a, carbon_e, BondOrder::Single)
            .expect("fused ring bond");
        let bond_g = mol
            .add_bond(carbon_e, carbon_f, BondOrder::Single)
            .expect("fused ring bond");
        let bond_h = mol
            .add_bond(carbon_f, carbon_g, BondOrder::Single)
            .expect("fused ring bond");
        let bond_i = mol
            .add_bond(carbon_g, carbon_b, BondOrder::Single)
            .expect("fused ring bond");
        let ring = Ring {
            atoms: vec![oxygen, carbon_a, carbon_b, carbon_c, carbon_d],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e],
        };
        let fused_neighbor = Ring {
            atoms: vec![carbon_a, carbon_e, carbon_f, carbon_g, carbon_b],
            bonds: vec![bond_f, bond_g, bond_h, bond_i, bond_b],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");
        let rings = vec![ring.clone(), fused_neighbor];

        assert_eq!(ring_hetero_donor_count(&mol, &ring), 1);
        assert!(ring_has_chalcogen_donor(&mol, &ring));
        assert!(!analysis.localized_has_active_element_donor(&mol, "O"));
        assert_eq!(analysis.localized_active_hetero_donor_count(&mol), 0);
        assert!(!ring_is_saturated_fused_ether_bridge_cleanup_candidate(
            &mol, &rings, 0, &ring, &analysis
        ));
    }

    #[test]
    fn saturated_fused_ether_cleanup_is_not_blocked_by_unrelated_heavy_chalcogen() {
        let mut mol = Molecule::new();
        let mut oxygen = Atom::new(Element::from_symbol("O").expect("test element"));
        oxygen.aromatic = true;
        let oxygen = mol.add_atom(oxygen);
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(oxygen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_d, oxygen, BondOrder::Single)
            .expect("ring bond");
        let fused_carbon_a =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let fused_carbon_b =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let fused_carbon_c =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_f = mol
            .add_bond(carbon_a, fused_carbon_a, BondOrder::Single)
            .expect("fused ring bond");
        let bond_g = mol
            .add_bond(fused_carbon_a, fused_carbon_b, BondOrder::Single)
            .expect("fused ring bond");
        let bond_h = mol
            .add_bond(fused_carbon_b, fused_carbon_c, BondOrder::Single)
            .expect("fused ring bond");
        let bond_i = mol
            .add_bond(fused_carbon_c, carbon_b, BondOrder::Single)
            .expect("fused ring bond");
        let sulfur = mol.add_atom(Atom::new(Element::from_symbol("S").expect("test element")));
        let ring = Ring {
            atoms: vec![oxygen, carbon_a, carbon_b, carbon_c, carbon_d],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e],
        };
        let fused_neighbor = Ring {
            atoms: vec![
                carbon_a,
                fused_carbon_a,
                fused_carbon_b,
                fused_carbon_c,
                carbon_b,
            ],
            bonds: vec![bond_f, bond_g, bond_h, bond_i, bond_b],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");
        let rings = vec![ring.clone(), fused_neighbor];
        let analyses = vec![
            analysis,
            RingAromaticityAnalysis::new(&mol, &rings[1]).expect("fused analysis"),
        ];

        assert!(mol
            .atom(sulfur)
            .is_ok_and(|atom| atom.element.symbol() == "S"));
        assert!(analyses[0].localized_has_active_element_donor(&mol, "O"));
        assert!(is_active_element_bridge_without_pi(
            &mol,
            &ring,
            &analyses[0],
            oxygen,
            "O"
        ));
        assert_eq!(analyses[0].localized_active_hetero_donor_count(&mol), 1);
        assert!(!analyses[0].localized_all_atoms_are_candidates());
        assert!(ring_is_saturated_fused_ether_bridge_cleanup_candidate(
            &mol,
            &rings,
            0,
            &ring,
            &analyses[0]
        ));

        clear_saturated_fused_ether_bridge_atoms(&mut mol, &rings, &analyses);

        assert!(!mol.atom(oxygen).expect("ring oxygen").aromatic);
        assert!(mol.atom(sulfur).is_ok());
    }

    #[test]
    fn exocyclic_pi_stealing_requires_electronegative_neighbor() {
        let mut mol = Molecule::new();
        let carbon = mol.add_atom(aromatic_carbon());
        let nitrogen = mol.add_atom(aromatic_atom("N"));
        let sulfur = mol.add_atom(aromatic_atom("S"));
        let carbon_b = mol.add_atom(aromatic_carbon());
        let carbon_c = mol.add_atom(aromatic_carbon());
        let phosphorus = mol.add_atom(Atom::new(Element::from_symbol("P").expect("test element")));
        let bond_a = mol
            .add_bond(carbon, nitrogen, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(nitrogen, sulfur, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(sulfur, carbon_b, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_c, carbon, BondOrder::Aromatic)
            .expect("ring bond");
        mol.add_bond(carbon, phosphorus, BondOrder::Double)
            .expect("exocyclic phosphorus");
        let ring = Ring {
            atoms: vec![carbon, nitrogen, sulfur, carbon_b, carbon_c],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e],
        };
        let localized = localized_ring_donor_analysis(&mol, &ring).expect("localized analysis");

        assert!(localized.has_active_element_donor(&mol, "N"));
        assert!(localized.has_active_chalcogen_donor(&mol));
        assert!(!ring_has_carbon_electronegative_exocyclic_pi_bond(
            &mol, &ring
        ));
        assert!(
            !aromatic_order_ring_allows_exocyclic_carbon_zero_contribution(&mol, &ring, &localized)
        );
        assert_eq!(
            aromatic_order_atom_donor_type(&mol, &ring, carbon, mol.atom(carbon).unwrap(), true)
                .expect("donor type"),
            AromaticElectronDonorType::One
        );
    }

    #[test]
    fn exocyclic_carbon_gate_requires_candidate_carbon_state() {
        let mut mol = Molecule::new();
        let mut carbon_atom = aromatic_carbon();
        carbon_atom.formal_charge = 1;
        carbon_atom.radical = Some(AtomRadical::Doublet);
        let carbon = mol.add_atom(carbon_atom);
        let nitrogen = mol.add_atom(aromatic_atom("N"));
        let sulfur = mol.add_atom(aromatic_atom("S"));
        let carbon_b = mol.add_atom(aromatic_carbon());
        let carbon_c = mol.add_atom(aromatic_carbon());
        let bond_a = mol
            .add_bond(carbon, nitrogen, BondOrder::Double)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(nitrogen, sulfur, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(sulfur, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Double)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_c, carbon, BondOrder::Single)
            .expect("ring bond");
        let oxygen = mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(carbon, oxygen, BondOrder::Double)
            .expect("terminal carbonyl");
        let ring = Ring {
            atoms: vec![carbon, nitrogen, sulfur, carbon_b, carbon_c],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e],
        };
        let localized = localized_ring_donor_analysis(&mol, &ring).expect("localized analysis");

        assert!(ring_has_carbon_electronegative_exocyclic_pi_bond(
            &mol, &ring
        ));
        assert!(!localized.atom_is_element_candidate(&mol, carbon, "C"));
        assert!(
            !ring_has_candidate_carbon_electronegative_exocyclic_pi_bond(&mol, &ring, &localized)
        );
        assert!(
            !aromatic_order_ring_allows_exocyclic_carbon_zero_contribution(&mol, &ring, &localized)
        );

        let aromatic =
            aromatic_ring_donor_analysis(&mol, &ring, &localized).expect("aromatic analysis");
        assert_eq!(aromatic.fixed_electron_count, Some(0));
    }

    #[test]
    fn ring_oxo_chalcogen_cleanup_uses_inactive_donor_state() {
        let mut mol = Molecule::new();
        let mut sulfur = Atom::new(Element::from_symbol("S").expect("test element"));
        sulfur.aromatic = true;
        let sulfur = mol.add_atom(sulfur);
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(sulfur, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_d, sulfur, BondOrder::Single)
            .expect("ring bond");
        let exocyclic_oxygen =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(sulfur, exocyclic_oxygen, BondOrder::Double)
            .expect("terminal sulfoxide");
        mol.bonds[bond_a.index()]
            .as_mut()
            .expect("ring bond")
            .aromatic = true;

        let ring = Ring {
            atoms: vec![sulfur, carbon_a, carbon_b, carbon_c, carbon_d],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e],
        };
        let analyses = vec![RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis")];

        assert!(ring_has_chalcogen_donor(&mol, &ring));
        assert!(is_ring_oxo_chalcogen(&mol, &ring, sulfur));
        assert!(!analyses[0].localized_atom_has_active_chalcogen_donor(&mol, sulfur));
        assert!(is_ring_oxo_chalcogen_cleanup_candidate(
            &mol,
            &ring,
            &analyses[0],
            sulfur
        ));

        clear_ring_oxo_chalcogen_atoms(&mut mol, &[ring], &analyses);

        assert!(!mol.atom(sulfur).expect("ring sulfur").aromatic);
        assert!(!mol.bond(bond_a).expect("ring bond").aromatic);
    }

    #[test]
    fn saturated_chalcogen_bridge_cleanup_requires_candidate_state() {
        let mut mol = Molecule::new();
        let mut sulfur = aromatic_atom("S");
        sulfur.radical = Some(AtomRadical::Doublet);
        let sulfur = mol.add_atom(sulfur);
        let carbon_a = mol.add_atom(aromatic_carbon());
        let carbon_b = mol.add_atom(aromatic_carbon());
        let carbon_c = mol.add_atom(aromatic_carbon());
        let carbon_d = mol.add_atom(aromatic_carbon());
        let carbon_e = mol.add_atom(aromatic_carbon());
        let bond_a = mol
            .add_bond(sulfur, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_d, carbon_e, BondOrder::Single)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_e, sulfur, BondOrder::Single)
            .expect("ring bond");
        for bond_id in [bond_a, bond_b, bond_c, bond_d, bond_e, bond_f] {
            mol.bonds[bond_id.index()]
                .as_mut()
                .expect("live bond")
                .aromatic = true;
        }

        let fused_carbon = mol.add_atom(aromatic_carbon());
        let bond_g = mol
            .add_bond(carbon_a, fused_carbon, BondOrder::Single)
            .expect("fused bond");
        let bond_h = mol
            .add_bond(fused_carbon, carbon_b, BondOrder::Single)
            .expect("fused bond");
        let ring = Ring {
            atoms: vec![sulfur, carbon_a, carbon_b, carbon_c, carbon_d, carbon_e],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let fused_neighbor = Ring {
            atoms: vec![carbon_a, fused_carbon, carbon_b],
            bonds: vec![bond_g, bond_h, bond_b],
        };
        let rings = vec![ring.clone(), fused_neighbor];
        let analyses = rings
            .iter()
            .map(|ring| RingAromaticityAnalysis::new(&mol, ring).expect("ring analysis"))
            .collect::<Vec<_>>();

        assert!(is_saturated_chalcogen_bridge(&mol, &ring, sulfur));
        assert!(!analyses[0].localized_atom_is_chalcogen_candidate(&mol, sulfur));
        assert!(!is_saturated_chalcogen_bridge_cleanup_candidate(
            &mol,
            &ring,
            &analyses[0],
            sulfur
        ));

        clear_saturated_chalcogen_bridge_atoms(&mut mol, &rings, &analyses);

        assert!(mol.atom(sulfur).expect("ring sulfur").aromatic);
        assert!(mol.bond(bond_a).expect("ring bond").aromatic);
        assert!(mol.bond(bond_f).expect("ring bond").aromatic);
    }

    #[test]
    fn aromatic_amidine_cleanup_requires_candidate_nitrogen_state() {
        let mut mol = Molecule::new();
        let mut carbon_atom = aromatic_carbon();
        carbon_atom.formal_charge = 1;
        let carbon = mol.add_atom(carbon_atom);
        let mut imine_nitrogen_atom = aromatic_atom("N");
        imine_nitrogen_atom.radical = Some(AtomRadical::Doublet);
        let imine_nitrogen = mol.add_atom(imine_nitrogen_atom);
        let single_nitrogen = mol.add_atom(aromatic_atom("N"));
        let saturated_carbon =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(carbon, imine_nitrogen, BondOrder::Double)
            .expect("imine bond");
        let bond_b = mol
            .add_bond(carbon, single_nitrogen, BondOrder::Single)
            .expect("single nitrogen bond");
        let bond_c = mol
            .add_bond(carbon, saturated_carbon, BondOrder::Single)
            .expect("saturated carbon bond");
        let bond_d = mol
            .add_bond(saturated_carbon, single_nitrogen, BondOrder::Single)
            .expect("ring closure bond");
        for bond_id in [bond_a, bond_b, bond_d] {
            mol.bonds[bond_id.index()]
                .as_mut()
                .expect("live bond")
                .aromatic = true;
        }
        let ring = Ring {
            atoms: vec![carbon, saturated_carbon, single_nitrogen],
            bonds: vec![bond_b, bond_c, bond_d],
        };

        assert!(mol.atom(carbon).expect("amidine carbon").aromatic);
        assert!(
            mol.atom(imine_nitrogen)
                .expect("imine nitrogen")
                .element
                .symbol()
                == "N"
        );
        assert!(!atom_is_rdkit_aromatic_element_candidate(
            &mol,
            imine_nitrogen,
            "N"
        ));
        assert!(atom_is_rdkit_aromatic_element_candidate(
            &mol,
            single_nitrogen,
            "N"
        ));

        clear_aromatic_amidine_carbon_atoms(&mut mol, &[ring]);

        assert!(mol.atom(carbon).expect("amidine carbon").aromatic);
        assert!(
            mol.atom(imine_nitrogen)
                .expect("radical imine nitrogen")
                .aromatic
        );
        assert!(mol.atom(single_nitrogen).expect("single nitrogen").aromatic);
        assert!(mol.bond(bond_a).expect("imine bond").aromatic);
    }

    #[test]
    fn aromatic_amidine_cleanup_requires_candidate_carbon_state() {
        let mut mol = Molecule::new();
        let mut carbon_atom = aromatic_carbon();
        carbon_atom.formal_charge = 1;
        carbon_atom.radical = Some(AtomRadical::Doublet);
        let carbon = mol.add_atom(carbon_atom);
        let imine_nitrogen = mol.add_atom(aromatic_atom("N"));
        let single_nitrogen = mol.add_atom(aromatic_atom("N"));
        let saturated_carbon =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(carbon, imine_nitrogen, BondOrder::Double)
            .expect("imine bond");
        let bond_b = mol
            .add_bond(carbon, single_nitrogen, BondOrder::Single)
            .expect("single nitrogen bond");
        let bond_c = mol
            .add_bond(carbon, saturated_carbon, BondOrder::Single)
            .expect("saturated carbon bond");
        let bond_d = mol
            .add_bond(saturated_carbon, single_nitrogen, BondOrder::Single)
            .expect("ring closure bond");
        for bond_id in [bond_a, bond_b, bond_d] {
            mol.bonds[bond_id.index()]
                .as_mut()
                .expect("live bond")
                .aromatic = true;
        }
        let ring = Ring {
            atoms: vec![carbon, saturated_carbon, single_nitrogen],
            bonds: vec![bond_b, bond_c, bond_d],
        };

        assert!(!atom_is_rdkit_aromatic_element_candidate(&mol, carbon, "C"));
        assert!(atom_is_rdkit_aromatic_element_candidate(
            &mol,
            imine_nitrogen,
            "N"
        ));
        assert!(atom_is_rdkit_aromatic_element_candidate(
            &mol,
            single_nitrogen,
            "N"
        ));

        clear_aromatic_amidine_carbon_atoms(&mut mol, &[ring]);

        assert!(mol.atom(carbon).expect("amidine carbon").aromatic);
        assert!(mol.atom(imine_nitrogen).expect("imine nitrogen").aromatic);
        assert!(mol.atom(single_nitrogen).expect("single nitrogen").aromatic);
        assert!(mol.bond(bond_a).expect("imine bond").aromatic);
    }

    #[test]
    fn fused_lactone_cleanup_uses_active_chalcogen_donor_state() {
        let mut mol = Molecule::new();
        let mut oxygen = Atom::new(Element::from_symbol("O").expect("test element"));
        oxygen.explicit_hydrogens = 2;
        oxygen.no_implicit_hydrogens = true;
        let oxygen = mol.add_atom(oxygen);
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_e = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(oxygen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_d, carbon_e, BondOrder::Single)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_e, oxygen, BondOrder::Single)
            .expect("ring bond");
        let exocyclic_oxygen =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(carbon_c, exocyclic_oxygen, BondOrder::Double)
            .expect("terminal carbonyl");
        let ring = Ring {
            atoms: vec![oxygen, carbon_a, carbon_b, carbon_c, carbon_d, carbon_e],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");

        assert!(ring_has_chalcogen_donor(&mol, &ring));
        assert!(is_chalcogen_bridge_without_pi(&mol, &ring, oxygen));
        assert_eq!(ring_terminal_exocyclic_pi_bond_count(&mol, &ring), 1);
        assert!(!analysis.localized_has_active_chalcogen_donor(&mol));
        assert!(!is_active_chalcogen_bridge_without_pi(
            &mol, &ring, &analysis, oxygen
        ));
        assert!(!ring_is_fused_lactone_bridge_cleanup_candidate(
            &mol, &ring, &analysis, true
        ));
    }

    #[test]
    fn chalcogen_bridge_without_pi_uses_rdkit_radical_eligibility() {
        let mut mol = Molecule::new();
        let mut oxygen_atom = Atom::new(Element::from_symbol("O").expect("test element"));
        oxygen_atom.radical = Some(AtomRadical::Doublet);
        let oxygen = mol.add_atom(oxygen_atom);
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_e = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(oxygen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_d, carbon_e, BondOrder::Single)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_e, oxygen, BondOrder::Single)
            .expect("ring bond");
        let ring = Ring {
            atoms: vec![oxygen, carbon_a, carbon_b, carbon_c, carbon_d, carbon_e],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };

        assert!(!atom_passes_rdkit_aromatic_radical_eligibility(
            mol.atom(oxygen).expect("ring oxygen")
        ));
        assert!(!is_chalcogen_bridge_without_pi(&mol, &ring, oxygen));
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");
        assert!(!is_active_chalcogen_bridge_without_pi(
            &mol, &ring, &analysis, oxygen
        ));
    }

    #[test]
    fn fused_lactam_enone_cleanup_uses_active_nitrogen_donor_state() {
        let mut mol = Molecule::new();
        let nitrogen = mol.add_atom(cation_with_one_implicit_hydrogen("N"));
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_e = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_d, carbon_e, BondOrder::Single)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_e, nitrogen, BondOrder::Single)
            .expect("ring bond");
        let carbonyl_oxygen =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(carbon_c, carbonyl_oxygen, BondOrder::Double)
            .expect("terminal carbonyl");
        let fused_carbon_a =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let fused_carbon_b =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let fused_carbon_c =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_g = mol
            .add_bond(carbon_a, fused_carbon_a, BondOrder::Single)
            .expect("fused ring bond");
        let bond_h = mol
            .add_bond(fused_carbon_a, fused_carbon_b, BondOrder::Single)
            .expect("fused ring bond");
        let bond_i = mol
            .add_bond(fused_carbon_b, fused_carbon_c, BondOrder::Single)
            .expect("fused ring bond");
        let bond_j = mol
            .add_bond(fused_carbon_c, carbon_b, BondOrder::Single)
            .expect("fused ring bond");
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, carbon_b, carbon_c, carbon_d, carbon_e],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let fused_neighbor = Ring {
            atoms: vec![
                carbon_a,
                fused_carbon_a,
                fused_carbon_b,
                fused_carbon_c,
                carbon_b,
            ],
            bonds: vec![bond_g, bond_h, bond_i, bond_j, bond_b],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");
        let rings = vec![ring.clone(), fused_neighbor];

        assert!(ring_contains_element(&mol, &ring, "N"));
        assert_eq!(ring_terminal_exocyclic_pi_bond_count(&mol, &ring), 1);
        assert!(!analysis.localized_has_active_element_donor(&mol, "N"));
        assert!(!ring_is_fused_lactam_enone_cleanup_candidate(
            &mol, &rings, 0, &ring, &analysis
        ));
    }

    #[test]
    fn lactam_enone_carbon_uses_active_neighbor_nitrogen_state() {
        let mut mol = Molecule::new();
        let mut inactive_nitrogen_atom =
            Atom::new(Element::from_symbol("N").expect("test element"));
        inactive_nitrogen_atom.radical = Some(AtomRadical::Doublet);
        let inactive_nitrogen = mol.add_atom(inactive_nitrogen_atom);
        let target_carbon = mol.add_atom(aromatic_carbon());
        let saturated_carbon =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon = mol.add_atom(aromatic_carbon());
        let active_nitrogen = mol.add_atom(aromatic_atom("N"));
        let bond_a = mol
            .add_bond(inactive_nitrogen, target_carbon, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(target_carbon, saturated_carbon, BondOrder::Double)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(saturated_carbon, carbon, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon, active_nitrogen, BondOrder::Double)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(active_nitrogen, inactive_nitrogen, BondOrder::Single)
            .expect("ring bond");
        let ring = Ring {
            atoms: vec![
                inactive_nitrogen,
                target_carbon,
                saturated_carbon,
                carbon,
                active_nitrogen,
            ],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");

        assert!(analysis.localized_has_active_element_donor(&mol, "N"));
        assert!(!analysis.localized_atom_has_active_element_donor(&mol, inactive_nitrogen, "N"));
        assert!(!is_aromatic_lactam_enone_carbon(
            &mol,
            &ring,
            &analysis,
            target_carbon
        ));
    }

    #[test]
    fn terminal_aromatic_imine_cleanup_uses_rdkit_candidate_state() {
        let mut mol = Molecule::new();
        let mut nitrogen_atom = aromatic_atom("N");
        nitrogen_atom.radical = Some(AtomRadical::Doublet);
        let nitrogen = mol.add_atom(nitrogen_atom);
        let carbon = mol.add_atom(aromatic_carbon());
        let bond = mol
            .add_bond(nitrogen, carbon, BondOrder::Double)
            .expect("imine bond");
        mol.bonds[bond.index()]
            .as_mut()
            .expect("live bond")
            .aromatic = true;

        assert!(!atom_is_rdkit_aromatic_element_candidate(
            &mol, nitrogen, "N"
        ));
        clear_terminal_aromatic_imine_fragments(&mut mol);
        assert!(mol.atom(nitrogen).expect("nitrogen").aromatic);
        assert!(mol.atom(carbon).expect("carbon").aromatic);
        assert!(mol.bond(bond).expect("imine bond").aromatic);
    }

    #[test]
    fn terminal_chalcogen_oxo_cleanup_uses_active_nitrogen_donor_state() {
        let mut mol = Molecule::new();
        let mut nitrogen_atom = cation_with_one_implicit_hydrogen("N");
        nitrogen_atom.aromatic = true;
        let nitrogen = mol.add_atom(nitrogen_atom);
        let carbon_a = mol.add_atom(aromatic_carbon());
        let sulfur = mol.add_atom(aromatic_atom("S"));
        let carbon_b = mol.add_atom(aromatic_carbon());
        let carbon_c = mol.add_atom(aromatic_carbon());
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, sulfur, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(sulfur, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_c, nitrogen, BondOrder::Single)
            .expect("ring bond");
        let exocyclic_oxygen =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(sulfur, exocyclic_oxygen, BondOrder::Double)
            .expect("terminal sulfoxide");
        for bond_id in [bond_a, bond_b, bond_c, bond_d, bond_e] {
            mol.bonds[bond_id.index()]
                .as_mut()
                .expect("live bond")
                .aromatic = true;
        }
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, sulfur, carbon_b, carbon_c],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");
        let rings = vec![ring.clone()];
        let analyses = vec![analysis];

        assert!(ring_contains_element(&mol, &ring, "N"));
        assert!(ring_has_terminal_chalcogen_exocyclic_pi_bond(&mol, &ring));
        assert!(!analyses[0].localized_has_active_element_donor(&mol, "N"));
        clear_terminal_chalcogen_oxo_ring_atoms(&mut mol, &rings, &analyses);
        assert!(ring
            .atoms
            .iter()
            .all(|atom_id| mol.atom(*atom_id).is_ok_and(|atom| atom.aromatic)));
        assert!(ring
            .bonds
            .iter()
            .all(|bond_id| mol.bond(*bond_id).is_ok_and(|bond| bond.aromatic)));
    }

    #[test]
    fn terminal_chalcogen_oxo_cleanup_uses_inactive_chalcogen_state() {
        let mut mol = Molecule::new();
        let nitrogen = mol.add_atom(aromatic_atom("N"));
        let carbon_a = mol.add_atom(aromatic_carbon());
        let sulfur = mol.add_atom(aromatic_atom("S"));
        let carbon_b = mol.add_atom(aromatic_carbon());
        let carbon_c = mol.add_atom(aromatic_carbon());
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, sulfur, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(sulfur, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_c, nitrogen, BondOrder::Single)
            .expect("ring bond");
        let exocyclic_oxygen =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(sulfur, exocyclic_oxygen, BondOrder::Double)
            .expect("terminal sulfoxide");
        for bond_id in [bond_a, bond_b, bond_c, bond_d, bond_e] {
            mol.bonds[bond_id.index()]
                .as_mut()
                .expect("live bond")
                .aromatic = true;
        }
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, sulfur, carbon_b, carbon_c],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e],
        };
        let analyses = vec![RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis")];
        let rings = vec![ring.clone()];

        assert!(analyses[0].localized_has_active_element_donor(&mol, "N"));
        assert!(ring_has_terminal_chalcogen_exocyclic_pi_bond(&mol, &ring));
        assert!(!analyses[0].localized_atom_has_active_chalcogen_donor(&mol, sulfur));
        assert!(ring_has_terminal_inactive_chalcogen_exocyclic_pi_bond(
            &mol,
            &ring,
            &analyses[0]
        ));

        clear_terminal_chalcogen_oxo_ring_atoms(&mut mol, &rings, &analyses);

        assert!(ring
            .atoms
            .iter()
            .all(|atom_id| mol.atom(*atom_id).is_ok_and(|atom| !atom.aromatic)));
        assert!(ring
            .bonds
            .iter()
            .all(|bond_id| mol.bond(*bond_id).is_ok_and(|bond| !bond.aromatic)));
    }

    #[test]
    fn inactive_terminal_chalcogen_cleanup_uses_rdkit_radical_eligibility() {
        let mut mol = Molecule::new();
        let nitrogen = mol.add_atom(aromatic_atom("N"));
        let carbon_a = mol.add_atom(aromatic_carbon());
        let mut sulfur_atom = aromatic_atom("S");
        sulfur_atom.radical = Some(AtomRadical::Doublet);
        let sulfur = mol.add_atom(sulfur_atom);
        let carbon_b = mol.add_atom(aromatic_carbon());
        let carbon_c = mol.add_atom(aromatic_carbon());
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, sulfur, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(sulfur, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_c, nitrogen, BondOrder::Single)
            .expect("ring bond");
        let exocyclic_oxygen =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(sulfur, exocyclic_oxygen, BondOrder::Double)
            .expect("terminal sulfoxide");
        for bond_id in [bond_a, bond_b, bond_c, bond_d, bond_e] {
            mol.bonds[bond_id.index()]
                .as_mut()
                .expect("live bond")
                .aromatic = true;
        }
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, sulfur, carbon_b, carbon_c],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e],
        };
        let analyses = vec![RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis")];
        let rings = vec![ring.clone()];

        assert!(analyses[0].localized_has_active_element_donor(&mol, "N"));
        assert!(ring_has_terminal_chalcogen_exocyclic_pi_bond(&mol, &ring));
        assert!(!analyses[0].localized_atom_is_chalcogen_candidate(&mol, sulfur));
        assert!(!ring_has_terminal_inactive_chalcogen_exocyclic_pi_bond(
            &mol,
            &ring,
            &analyses[0]
        ));
        assert!(!is_ring_oxo_chalcogen_cleanup_candidate(
            &mol,
            &ring,
            &analyses[0],
            sulfur
        ));

        clear_terminal_chalcogen_oxo_ring_atoms(&mut mol, &rings, &analyses);
        clear_ring_oxo_chalcogen_atoms(&mut mol, &rings, &analyses);

        assert!(ring
            .atoms
            .iter()
            .all(|atom_id| mol.atom(*atom_id).is_ok_and(|atom| atom.aromatic)));
        assert!(ring
            .bonds
            .iter()
            .all(|bond_id| mol.bond(*bond_id).is_ok_and(|bond| bond.aromatic)));
    }

    #[test]
    fn saturated_ring_carbon_uses_rdkit_radical_eligibility() {
        let mut mol = Molecule::new();
        let mut carbon_atom = Atom::new(Element::from_symbol("C").expect("test element"));
        carbon_atom.formal_charge = 1;
        carbon_atom.radical = Some(AtomRadical::Doublet);
        let carbon = mol.add_atom(carbon_atom);
        let nitrogen = mol.add_atom(Atom::new(Element::from_symbol("N").expect("test element")));
        let oxygen = mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        let sulfur = mol.add_atom(Atom::new(Element::from_symbol("S").expect("test element")));
        let phosphorus = mol.add_atom(Atom::new(Element::from_symbol("P").expect("test element")));
        let bond_a = mol
            .add_bond(carbon, nitrogen, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(nitrogen, oxygen, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(oxygen, sulfur, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(sulfur, phosphorus, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(phosphorus, carbon, BondOrder::Single)
            .expect("ring bond");
        let ring = Ring {
            atoms: vec![carbon, nitrogen, oxygen, sulfur, phosphorus],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e],
        };

        assert!(!atom_passes_rdkit_aromatic_radical_eligibility(
            mol.atom(carbon).expect("charged carbon radical")
        ));
        assert!(!is_saturated_ring_carbon(&mol, &ring, carbon));
        assert!(!ring_has_saturated_carbon_atom(&mol, &ring));
    }

    #[test]
    fn terminal_chalcogen_gate_uses_candidate_state() {
        let mut mol = Molecule::new();
        let nitrogen = mol.add_atom(aromatic_atom("N"));
        let carbon_a = mol.add_atom(aromatic_carbon());
        let mut sulfur_atom = aromatic_atom("S");
        sulfur_atom.radical = Some(AtomRadical::Doublet);
        let sulfur = mol.add_atom(sulfur_atom);
        let oxygen = mol.add_atom(aromatic_atom("O"));
        let carbon_b = mol.add_atom(aromatic_carbon());
        let exocyclic_sulfur_oxygen =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        let exocyclic_carbonyl_oxygen =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Double)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, sulfur, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(sulfur, oxygen, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(oxygen, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_b, nitrogen, BondOrder::Single)
            .expect("ring bond");
        mol.add_bond(sulfur, exocyclic_sulfur_oxygen, BondOrder::Double)
            .expect("terminal sulfur oxo");
        mol.add_bond(carbon_b, exocyclic_carbonyl_oxygen, BondOrder::Double)
            .expect("terminal carbonyl");
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, sulfur, oxygen, carbon_b],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e],
        };
        let localized = localized_ring_donor_analysis(&mol, &ring).expect("localized analysis");

        assert!(localized.has_active_element_donor(&mol, "N"));
        assert!(localized.has_active_chalcogen_donor(&mol));
        assert!(ring_has_terminal_chalcogen_exocyclic_pi_bond(&mol, &ring));
        assert!(ring_has_carbon_electronegative_exocyclic_pi_bond(
            &mol, &ring
        ));
        assert!(!localized.atom_is_heavy_chalcogen_candidate(&mol, sulfur));
        assert!(!ring_has_terminal_candidate_chalcogen_exocyclic_pi_bond(
            &mol, &ring, &localized
        ));
        let aromatic =
            aromatic_ring_donor_analysis(&mol, &ring, &localized).expect("aromatic analysis");
        assert_eq!(aromatic.fixed_electron_count, None);
    }

    #[test]
    fn saturated_lactam_carbonyl_cleanup_uses_active_nitrogen_donor_state() {
        let mut mol = Molecule::new();
        let nitrogen = mol.add_atom(cation_with_one_implicit_hydrogen("N"));
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_e = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_d, carbon_e, BondOrder::Single)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_e, nitrogen, BondOrder::Single)
            .expect("ring bond");
        let carbonyl_oxygen =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(carbon_c, carbonyl_oxygen, BondOrder::Double)
            .expect("terminal carbonyl");
        let fused_carbon_a =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let fused_carbon_b =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let fused_carbon_c =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_g = mol
            .add_bond(carbon_a, fused_carbon_a, BondOrder::Single)
            .expect("fused ring bond");
        let bond_h = mol
            .add_bond(fused_carbon_a, fused_carbon_b, BondOrder::Single)
            .expect("fused ring bond");
        let bond_i = mol
            .add_bond(fused_carbon_b, fused_carbon_c, BondOrder::Single)
            .expect("fused ring bond");
        let bond_j = mol
            .add_bond(fused_carbon_c, carbon_b, BondOrder::Single)
            .expect("fused ring bond");
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, carbon_b, carbon_c, carbon_d, carbon_e],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let fused_neighbor = Ring {
            atoms: vec![
                carbon_a,
                fused_carbon_a,
                fused_carbon_b,
                fused_carbon_c,
                carbon_b,
            ],
            bonds: vec![bond_g, bond_h, bond_i, bond_j, bond_b],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");
        let rings = vec![ring.clone(), fused_neighbor];

        assert!(ring_contains_element(&mol, &ring, "N"));
        assert_eq!(ring_terminal_exocyclic_pi_bond_count(&mol, &ring), 1);
        assert!(ring_has_saturated_carbon_atom(&mol, &ring));
        assert!(!analysis.localized_has_active_element_donor(&mol, "N"));
        assert!(!ring_is_saturated_fused_lactam_carbonyl_cleanup_candidate(
            &mol, &rings, 0, &ring, &analysis
        ));
    }

    #[test]
    fn imide_carbonyl_cleanup_uses_candidate_nitrogen_state() {
        let mut mol = Molecule::new();
        let mut nitrogen_atom = Atom::new(Element::from_symbol("N").expect("test element"));
        nitrogen_atom.explicit_hydrogens = 2;
        nitrogen_atom.no_implicit_hydrogens = true;
        let nitrogen = mol.add_atom(nitrogen_atom);
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_d, nitrogen, BondOrder::Single)
            .expect("ring bond");
        let carbonyl_oxygen_a =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(carbon_a, carbonyl_oxygen_a, BondOrder::Double)
            .expect("terminal carbonyl");
        let carbonyl_oxygen_b =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(carbon_d, carbonyl_oxygen_b, BondOrder::Double)
            .expect("terminal carbonyl");
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, carbon_b, carbon_c, carbon_d],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");
        let rings = vec![ring.clone()];
        let analyses = vec![analysis];

        assert!(ring_contains_element(&mol, &ring, "N"));
        assert_eq!(ring_terminal_exocyclic_pi_bond_count(&mol, &ring), 2);
        assert_eq!(ring_carbonyl_neighbor_count(&mol, &ring, nitrogen), 2);
        assert!(!analyses[0].localized_has_element_candidate(&mol, "N"));
        assert!(!ring_has_imide_nitrogen(&mol, &ring, &analyses[0]));
        assert!(!ring_is_imide_carbonyl_cleanup_candidate(
            &mol,
            &rings,
            0,
            &ring,
            &analyses[0],
            &analyses
        ));
    }

    #[test]
    fn fused_lactam_bridge_cleanup_uses_active_nitrogen_and_oxygen_donor_state() {
        let mut mol = Molecule::new();
        let nitrogen = mol.add_atom(cation_with_one_implicit_hydrogen("N"));
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let mut oxygen_atom = Atom::new(Element::from_symbol("O").expect("test element"));
        oxygen_atom.explicit_hydrogens = 2;
        oxygen_atom.no_implicit_hydrogens = true;
        let oxygen = mol.add_atom(oxygen_atom);
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, oxygen, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(oxygen, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_d, nitrogen, BondOrder::Single)
            .expect("ring bond");
        let carbonyl_oxygen =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(carbon_c, carbonyl_oxygen, BondOrder::Double)
            .expect("terminal carbonyl");
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, oxygen, carbon_b, carbon_c, carbon_d],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");

        assert!(ring_contains_element(&mol, &ring, "N"));
        assert!(ring_contains_element(&mol, &ring, "O"));
        assert_eq!(ring_terminal_exocyclic_pi_bond_count(&mol, &ring), 1);
        assert!(!analysis.localized_has_active_element_donor(&mol, "N"));
        assert!(!analysis.localized_has_active_element_donor(&mol, "O"));
        assert!(!ring_is_fused_lactam_bridge_cleanup_candidate(
            &mol, &ring, &analysis
        ));
    }

    #[test]
    fn aromatic_order_exocyclic_carbon_uses_active_donor_state() {
        let mut mol = Molecule::new();
        let mut nitrogen = Atom::new(Element::from_symbol("N").expect("test element"));
        nitrogen.explicit_hydrogens = 2;
        nitrogen.no_implicit_hydrogens = true;
        let nitrogen = mol.add_atom(nitrogen);
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let mut oxygen = Atom::new(Element::from_symbol("O").expect("test element"));
        oxygen.explicit_hydrogens = 2;
        oxygen.no_implicit_hydrogens = true;
        let oxygen = mol.add_atom(oxygen);
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, carbon_b, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(carbon_b, oxygen, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(oxygen, carbon_c, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_d, nitrogen, BondOrder::Aromatic)
            .expect("ring bond");
        let exocyclic_oxygen_a =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        let exocyclic_oxygen_b =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(carbon_a, exocyclic_oxygen_a, BondOrder::Double)
            .expect("terminal carbonyl");
        mol.add_bond(carbon_c, exocyclic_oxygen_b, BondOrder::Double)
            .expect("terminal carbonyl");
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, carbon_b, oxygen, carbon_c, carbon_d],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let localized = localized_ring_donor_analysis(&mol, &ring).expect("localized analysis");

        assert!(ring_contains_element(&mol, &ring, "N"));
        assert!(ring_has_chalcogen_donor(&mol, &ring));
        assert_eq!(ring_terminal_exocyclic_pi_bond_count(&mol, &ring), 2);
        assert!(!localized.has_active_element_donor(&mol, "N"));
        assert!(!localized.has_active_chalcogen_donor(&mol));
        assert!(
            !aromatic_order_ring_allows_exocyclic_carbon_zero_contribution(&mol, &ring, &localized)
        );
    }

    #[test]
    fn aromatic_order_single_exocyclic_pi_uses_variable_donor_analysis() {
        let mut mol = Molecule::new();
        let nitrogen = mol.add_atom(aromatic_atom("N"));
        let carbon_a = mol.add_atom(aromatic_carbon());
        let sulfur = mol.add_atom(aromatic_atom("S"));
        let carbon_b = mol.add_atom(aromatic_carbon());
        let carbon_c = mol.add_atom(aromatic_carbon());
        let carbon_d = mol.add_atom(aromatic_carbon());
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, sulfur, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(sulfur, carbon_b, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Aromatic)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_d, nitrogen, BondOrder::Aromatic)
            .expect("ring bond");
        let exocyclic_oxygen =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        mol.add_bond(carbon_c, exocyclic_oxygen, BondOrder::Double)
            .expect("terminal carbonyl");
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, sulfur, carbon_b, carbon_c, carbon_d],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let localized = localized_ring_donor_analysis(&mol, &ring).expect("localized analysis");

        assert!(localized.has_active_element_donor(&mol, "N"));
        assert!(localized.has_active_chalcogen_donor(&mol));
        assert!(ring_has_carbon_electronegative_exocyclic_pi_bond(
            &mol, &ring
        ));
        assert_eq!(ring_terminal_exocyclic_pi_bond_count(&mol, &ring), 1);
        assert!(
            aromatic_order_ring_allows_exocyclic_carbon_zero_contribution(&mol, &ring, &localized)
        );

        let aromatic =
            aromatic_order_ring_donor_analysis(&mol, &ring, &localized).expect("aromatic analysis");

        assert_eq!(aromatic.fixed_electron_count, None);
        assert!(aromatic.is_huckel_aromatic());
    }

    #[test]
    fn saturated_tertiary_amine_guard_uses_active_chalcogen_donor_state() {
        let mut mol = Molecule::new();
        let nitrogen = mol.add_atom(Atom::new(Element::from_symbol("N").expect("test element")));
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let mut oxygen = Atom::new(Element::from_symbol("O").expect("test element"));
        oxygen.explicit_hydrogens = 2;
        oxygen.no_implicit_hydrogens = true;
        let oxygen = mol.add_atom(oxygen);
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let methyl = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, oxygen, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(oxygen, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_d, nitrogen, BondOrder::Single)
            .expect("ring bond");
        mol.add_bond(nitrogen, methyl, BondOrder::Single)
            .expect("tertiary amine substituent");
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, oxygen, carbon_b, carbon_c, carbon_d],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");

        assert!(ring_has_chalcogen_donor(&mol, &ring));
        assert!(!analysis.localized_atom_has_active_chalcogen_donor(&mol, oxygen));
        assert_eq!(analysis.localized_active_hetero_donor_count(&mol), 1);
        assert!(!ring_has_saturated_active_chalcogen_donor(
            &mol, &ring, &analysis
        ));
        assert!(ring_has_saturated_tertiary_amine_without_donor_chalcogen(
            &mol, &ring, &analysis
        ));
    }

    #[test]
    fn saturated_tertiary_amine_guard_uses_rdkit_radical_eligibility() {
        let mut mol = Molecule::new();
        let mut nitrogen_atom = Atom::new(Element::from_symbol("N").expect("test element"));
        nitrogen_atom.radical = Some(AtomRadical::Doublet);
        let nitrogen = mol.add_atom(nitrogen_atom);
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let oxygen = mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let methyl = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, oxygen, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(oxygen, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_d, nitrogen, BondOrder::Single)
            .expect("ring bond");
        mol.add_bond(nitrogen, methyl, BondOrder::Single)
            .expect("tertiary amine substituent");
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, oxygen, carbon_b, carbon_c, carbon_d],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");

        assert!(!atom_passes_rdkit_aromatic_radical_eligibility(
            mol.atom(nitrogen).expect("ring nitrogen")
        ));
        assert!(!is_saturated_tertiary_amine(&mol, &ring, nitrogen));
        assert!(!ring_has_saturated_tertiary_amine_without_donor_chalcogen(
            &mol, &ring, &analysis
        ));
    }

    #[test]
    fn tertiary_amine_chalcogen_substituent_uses_rdkit_radical_eligibility() {
        let mut mol = Molecule::new();
        let nitrogen = mol.add_atom(Atom::new(Element::from_symbol("N").expect("test element")));
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let oxygen = mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let mut sulfur_atom = Atom::new(Element::from_symbol("S").expect("test element"));
        sulfur_atom.radical = Some(AtomRadical::Doublet);
        let sulfur = mol.add_atom(sulfur_atom);
        let exocyclic_oxygen =
            mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, oxygen, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(oxygen, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_d, nitrogen, BondOrder::Single)
            .expect("ring bond");
        mol.add_bond(nitrogen, sulfur, BondOrder::Single)
            .expect("tertiary amine substituent");
        mol.add_bond(sulfur, exocyclic_oxygen, BondOrder::Double)
            .expect("terminal sulfoxide");
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, oxygen, carbon_b, carbon_c, carbon_d],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");

        assert!(atom_has_terminal_exocyclic_pi_bond(&mol, &ring, sulfur));
        assert!(!atom_passes_rdkit_aromatic_radical_eligibility(
            mol.atom(sulfur).expect("sulfur substituent")
        ));
        assert!(!is_saturated_tertiary_amine_substituent(
            &mol, &ring, sulfur
        ));
        assert!(!is_saturated_tertiary_amine(&mol, &ring, nitrogen));
        assert!(!ring_has_saturated_tertiary_amine_without_donor_chalcogen(
            &mol, &ring, &analysis
        ));
    }

    #[test]
    fn tertiary_amine_carbon_substituent_uses_rdkit_radical_eligibility() {
        let mut mol = Molecule::new();
        let nitrogen = mol.add_atom(Atom::new(Element::from_symbol("N").expect("test element")));
        let carbon_a = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let oxygen = mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        let carbon_b = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_c = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let carbon_d = mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let mut substituent_atom = Atom::new(Element::from_symbol("C").expect("test element"));
        substituent_atom.formal_charge = 1;
        substituent_atom.radical = Some(AtomRadical::Doublet);
        let substituent = mol.add_atom(substituent_atom);
        let bond_a = mol
            .add_bond(nitrogen, carbon_a, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon_a, oxygen, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(oxygen, carbon_b, BondOrder::Single)
            .expect("ring bond");
        let bond_d = mol
            .add_bond(carbon_b, carbon_c, BondOrder::Single)
            .expect("ring bond");
        let bond_e = mol
            .add_bond(carbon_c, carbon_d, BondOrder::Single)
            .expect("ring bond");
        let bond_f = mol
            .add_bond(carbon_d, nitrogen, BondOrder::Single)
            .expect("ring bond");
        mol.add_bond(nitrogen, substituent, BondOrder::Single)
            .expect("tertiary amine substituent");
        let ring = Ring {
            atoms: vec![nitrogen, carbon_a, oxygen, carbon_b, carbon_c, carbon_d],
            bonds: vec![bond_a, bond_b, bond_c, bond_d, bond_e, bond_f],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");

        assert!(!atom_passes_rdkit_aromatic_radical_eligibility(
            mol.atom(substituent).expect("carbon substituent")
        ));
        assert!(!is_saturated_tertiary_amine_substituent(
            &mol,
            &ring,
            substituent
        ));
        assert!(!is_saturated_tertiary_amine(&mol, &ring, nitrogen));
        assert!(!ring_has_saturated_tertiary_amine_without_donor_chalcogen(
            &mol, &ring, &analysis
        ));
    }

    #[test]
    fn exocyclic_alkene_chalcogen_cleanup_uses_active_neighbor_donors() {
        let mut mol = Molecule::new();
        let nitrogen = mol.add_atom(cation_with_one_implicit_hydrogen("N"));
        let mut carbon = Atom::new(Element::from_symbol("C").expect("test element"));
        carbon.aromatic = true;
        let carbon = mol.add_atom(carbon);
        let oxygen = mol.add_atom(Atom::new(Element::from_symbol("O").expect("test element")));
        let exocyclic_carbon =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(nitrogen, carbon, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon, oxygen, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(oxygen, nitrogen, BondOrder::Single)
            .expect("ring bond");
        mol.add_bond(carbon, exocyclic_carbon, BondOrder::Double)
            .expect("exocyclic alkene");
        let ring = Ring {
            atoms: vec![nitrogen, carbon, oxygen],
            bonds: vec![bond_a, bond_b, bond_c],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");

        assert!(analysis.localized_has_active_chalcogen_donor(&mol));
        assert!(!analysis.localized_has_active_element_donor(&mol, "N"));
        assert!(!ring_atom_has_active_nitrogen_and_chalcogen_neighbors(
            &mol, &ring, &analysis, carbon
        ));
        assert!(!is_aromatic_ring_carbon_with_exocyclic_carbon_pi_bond(
            &mol, &ring, &analysis, carbon
        ));
    }

    #[test]
    fn exocyclic_alkene_cleanup_clears_only_active_ring_hetero_neighbors() {
        let mut mol = Molecule::new();
        let mut carbon = Atom::new(Element::from_symbol("C").expect("test element"));
        carbon.aromatic = true;
        let carbon = mol.add_atom(carbon);
        let mut nitrogen_atom = Atom::new(Element::from_symbol("N").expect("test element"));
        nitrogen_atom.aromatic = true;
        let nitrogen = mol.add_atom(nitrogen_atom);
        let mut oxygen_atom = Atom::new(Element::from_symbol("O").expect("test element"));
        oxygen_atom.aromatic = true;
        let oxygen = mol.add_atom(oxygen_atom);
        let mut exocyclic_nitrogen_atom =
            Atom::new(Element::from_symbol("N").expect("test element"));
        exocyclic_nitrogen_atom.aromatic = true;
        let exocyclic_nitrogen = mol.add_atom(exocyclic_nitrogen_atom);
        let exocyclic_carbon =
            mol.add_atom(Atom::new(Element::from_symbol("C").expect("test element")));
        let bond_a = mol
            .add_bond(nitrogen, carbon, BondOrder::Single)
            .expect("ring bond");
        let bond_b = mol
            .add_bond(carbon, oxygen, BondOrder::Single)
            .expect("ring bond");
        let bond_c = mol
            .add_bond(oxygen, nitrogen, BondOrder::Single)
            .expect("ring bond");
        mol.add_bond(carbon, exocyclic_carbon, BondOrder::Double)
            .expect("exocyclic alkene");
        mol.add_bond(carbon, exocyclic_nitrogen, BondOrder::Single)
            .expect("exocyclic neutral hetero substituent");
        let ring = Ring {
            atoms: vec![nitrogen, carbon, oxygen],
            bonds: vec![bond_a, bond_b, bond_c],
        };
        let analysis = RingAromaticityAnalysis::new(&mol, &ring).expect("ring analysis");
        let rings = vec![ring.clone()];
        let analyses = vec![analysis];

        assert!(analyses[0].localized_atom_has_active_element_donor(&mol, nitrogen, "N"));
        assert!(analyses[0].localized_atom_has_active_chalcogen_donor(&mol, oxygen));
        assert!(ring_atom_has_active_nitrogen_and_chalcogen_neighbors(
            &mol,
            &ring,
            &analyses[0],
            carbon
        ));
        let active_neighbors =
            active_hetero_neighbors_for_exocyclic_alkene_cleanup(&mol, &ring, &analyses[0], carbon);
        assert!(active_neighbors.contains(&nitrogen));
        assert!(active_neighbors.contains(&oxygen));
        assert!(!active_neighbors.contains(&exocyclic_nitrogen));

        clear_exocyclic_alkene_chalcogen_ring_atoms(&mut mol, &rings, &analyses);

        assert!(!mol.atom(carbon).expect("seed carbon").aromatic);
        assert!(!mol.atom(nitrogen).expect("ring nitrogen").aromatic);
        assert!(!mol.atom(oxygen).expect("ring oxygen").aromatic);
        assert!(
            mol.atom(exocyclic_nitrogen)
                .expect("exocyclic nitrogen")
                .aromatic
        );
    }

    fn cation_with_one_implicit_hydrogen(symbol: &str) -> Atom {
        let mut atom = Atom::new(Element::from_symbol(symbol).expect("test element"));
        atom.formal_charge = 1;
        atom.implicit_hydrogens = Some(1);
        atom.no_implicit_hydrogens = true;
        atom
    }

    fn aromatic_atom(symbol: &str) -> Atom {
        let mut atom = Atom::new(Element::from_symbol(symbol).expect("test element"));
        atom.aromatic = true;
        atom
    }

    fn aromatic_carbon() -> Atom {
        aromatic_atom("C")
    }

    fn fused_carbon_ring(
        mol: &mut Molecule,
        shared: Option<(AtomId, AtomId, BondId)>,
        size: usize,
    ) -> (Ring, AtomId, AtomId, BondId) {
        assert!(size >= 6);
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        if let Some((left, right, bond)) = shared {
            atoms.push(left);
            atoms.push(right);
            bonds.push(bond);
        } else {
            atoms.push(mol.add_atom(aromatic_carbon()));
            atoms.push(mol.add_atom(aromatic_carbon()));
            bonds.push(
                mol.add_bond(atoms[0], atoms[1], BondOrder::Aromatic)
                    .expect("initial fused bond"),
            );
        }

        let first_new_index = atoms.len();
        for _ in 0..(size - atoms.len()) {
            atoms.push(mol.add_atom(aromatic_carbon()));
        }

        for index in 1..(atoms.len() - 1) {
            bonds.push(
                mol.add_bond(atoms[index], atoms[index + 1], BondOrder::Aromatic)
                    .expect("ring bond"),
            );
        }
        bonds.push(
            mol.add_bond(
                *atoms.last().expect("ring atom"),
                atoms[0],
                BondOrder::Aromatic,
            )
            .expect("ring closure"),
        );

        let next_left_index = first_new_index + 1;
        let next_left = atoms[next_left_index];
        let next_right = atoms[next_left_index + 1];
        let next_bond = bonds[next_left_index];
        (Ring { atoms, bonds }, next_left, next_right, next_bond)
    }
}
