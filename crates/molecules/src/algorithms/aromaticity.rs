use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use super::*;
use crate::core::*;

const MAX_FUSED_AROMATIC_COMBINATION_RINGS: usize = 6;
const MAX_FUSED_AROMATIC_RING_SIZE: usize = 24;
const LARGE_FUSED_RING_SYSTEM_SEARCH_LIMIT: usize = 300;
const MAX_IMPORTED_AROMATIC_MATCHING_STATES: usize = 100_000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AromaticityModel {
    RdkitLike,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AromaticityError {
    UnsupportedElement(AtomId),
    InvalidAromaticRepresentation(AtomId),
    ImportedAromaticKekulizationLimit {
        atom: AtomId,
        examined_states: usize,
        limit: usize,
    },
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
            Self::ImportedAromaticKekulizationLimit {
                atom,
                examined_states,
                limit,
            } => write!(
                f,
                "imported aromatic kekulization limit exceeded at atom {atom}: examined {examined_states} matching states, limit {limit}"
            ),
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
    let mut staged = mol.clone();
    match model {
        AromaticityModel::RdkitLike => perceive_rdkit_like_aromaticity(&mut staged, ring_options),
    }?;
    *mol = staged;
    Ok(())
}

fn perceive_rdkit_like_aromaticity(
    mol: &mut Molecule,
    ring_options: RingPerceptionOptions,
) -> std::result::Result<(), AromaticityError> {
    let ring_set = perceive_ring_set_with_options(mol, ring_options)
        .map_err(AromaticityError::RingPerception)?;
    let imported_aromatic_components = imported_aromatic_bond_components(mol);
    for component in imported_aromatic_components {
        if !try_kekulize_imported_component(mol, &component)? {
            return Err(AromaticityError::InvalidAromaticRepresentation(
                component[0],
            ));
        }
    }
    assign_rdkit_like_localized_aromaticity(mol, &ring_set)
}

fn try_kekulize_imported_component(
    mol: &mut Molecule,
    component: &[AtomId],
) -> std::result::Result<bool, AromaticityError> {
    try_kekulize_imported_component_with_limit(
        mol,
        component,
        MAX_IMPORTED_AROMATIC_MATCHING_STATES,
    )
}

fn try_kekulize_imported_component_with_limit(
    mol: &mut Molecule,
    component: &[AtomId],
    max_matching_states: usize,
) -> std::result::Result<bool, AromaticityError> {
    let component_atoms = component.iter().copied().collect::<BTreeSet<_>>();
    let mut demand = BTreeSet::new();
    for atom_id in component {
        let Ok(atom) = mol.atom(*atom_id) else {
            return Ok(false);
        };
        let Some(default_valence) = rdkit_charge_adjusted_default_valence(atom) else {
            return Ok(false);
        };
        let target_valence = default_valence
            .saturating_sub(atom.radical.map_or(0, AtomRadical::unpaired_electron_count));
        let bond_valence = mol
            .incident_bonds(*atom_id)
            .ok()
            .into_iter()
            .flatten()
            .map(|(_, bond)| match bond.order {
                BondOrder::Zero | BondOrder::Dative => 0,
                BondOrder::Single | BondOrder::Aromatic => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Quadruple => 4,
            })
            .sum::<usize>();
        let occupied_valence = bond_valence
            .saturating_add(usize::from(atom.explicit_hydrogens))
            .saturating_add(usize::from(
                mol.implicit_hydrogens(*atom_id).ok().flatten().unwrap_or(0),
            ));
        let required_double_bonds = usize::from(target_valence).checked_sub(occupied_valence);
        match required_double_bonds {
            Some(0) => {}
            Some(1) => {
                demand.insert(*atom_id);
            }
            _ => return Ok(false),
        }
    }
    if demand.len() % 2 != 0 {
        return Ok(false);
    }

    let mut adjacency = BTreeMap::<AtomId, Vec<(AtomId, BondId)>>::new();
    for (bond_id, bond) in mol.bonds().filter(|(_, bond)| {
        bond.order == BondOrder::Aromatic
            && component_atoms.contains(&bond.a())
            && component_atoms.contains(&bond.b())
    }) {
        if demand.contains(&bond.a()) && demand.contains(&bond.b()) {
            adjacency
                .entry(bond.a())
                .or_default()
                .push((bond.b(), bond_id));
            adjacency
                .entry(bond.b())
                .or_default()
                .push((bond.a(), bond_id));
        }
    }
    for neighbors in adjacency.values_mut() {
        neighbors.sort();
    }

    let mut stack = vec![(demand, Vec::<BondId>::new())];
    let mut examined_states = 0usize;
    let selected_double_bonds = loop {
        let Some((unmatched, selected)) = stack.pop() else {
            return Ok(false);
        };
        examined_states += 1;
        if examined_states > max_matching_states {
            return Err(AromaticityError::ImportedAromaticKekulizationLimit {
                atom: component[0],
                examined_states,
                limit: max_matching_states,
            });
        }
        if unmatched.is_empty() {
            break selected;
        }
        let Some(atom_id) = unmatched.iter().copied().min_by_key(|atom_id| {
            adjacency
                .get(atom_id)
                .map(|neighbors| {
                    neighbors
                        .iter()
                        .filter(|(neighbor, _)| unmatched.contains(neighbor))
                        .count()
                })
                .unwrap_or(0)
        }) else {
            return Ok(false);
        };
        let candidates = adjacency
            .get(&atom_id)
            .into_iter()
            .flatten()
            .filter(|(neighbor, _)| unmatched.contains(neighbor))
            .copied()
            .collect::<Vec<_>>();
        if candidates.is_empty() {
            continue;
        }
        for (neighbor, bond_id) in candidates.into_iter().rev() {
            let mut next_unmatched = unmatched.clone();
            next_unmatched.remove(&atom_id);
            next_unmatched.remove(&neighbor);
            let mut next_selected = selected.clone();
            next_selected.push(bond_id);
            stack.push((next_unmatched, next_selected));
        }
    };

    let selected_double_bonds = selected_double_bonds.into_iter().collect::<BTreeSet<_>>();
    for atom_id in component {
        mol.set_atom_aromatic(*atom_id, false);
    }
    for (bond_id, bond) in mol
        .bonds
        .iter_mut()
        .enumerate()
        .filter_map(|(index, bond)| bond.as_mut().map(|bond| (BondId::new(index as u32), bond)))
    {
        if bond.order == BondOrder::Aromatic
            && component_atoms.contains(&bond.a())
            && component_atoms.contains(&bond.b())
        {
            bond.order = if selected_double_bonds.contains(&bond_id) {
                BondOrder::Double
            } else {
                BondOrder::Single
            };
        }
    }
    Ok(true)
}

fn assign_rdkit_like_localized_aromaticity(
    mol: &mut Molecule,
    ring_set: &RingSet,
) -> std::result::Result<(), AromaticityError> {
    mol.begin_aromaticity(AromaticityProvenance::Perceived(
        AromaticityModel::RdkitLike,
    ));

    let mut donors = vec![AromaticElectronDonorType::None; mol.atoms.len()];
    let mut atom_candidates = vec![false; mol.atoms.len()];
    for (atom_id, atom) in mol.atoms() {
        let donor = rdkit_localized_atom_donor_type(mol, atom_id, atom);
        donors[atom_id.index()] = donor;
        atom_candidates[atom_id.index()] = atom_is_rdkit_aromatic_candidate_for_donor(
            mol,
            atom_id,
            atom,
            donor,
            RdkitAromaticCandidateOptions::default(),
        );
    }

    let candidates = ring_set
        .rings()
        .iter()
        .enumerate()
        .filter_map(|(index, ring)| {
            ring.atoms
                .iter()
                .all(|atom| atom_candidates[atom.index()])
                .then_some(index)
        })
        .collect::<Vec<_>>();
    let components = rdkit_fused_ring_components(ring_set.rings(), &candidates);
    for component in components {
        apply_rdkit_huckel_to_fused_component(mol, ring_set.rings(), &component, &donors);
    }

    Ok(())
}

fn rdkit_localized_atom_donor_type(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
) -> AromaticElectronDonorType {
    let Some(mut electrons) = count_rdkit_like_atom_pi_electrons(mol, atom_id, atom) else {
        return AromaticElectronDonorType::None;
    };
    let noncyclic_pi_neighbor = atom_noncyclic_pi_neighbor(mol, atom_id);
    let has_cyclic_pi_bond = atom_has_cyclic_pi_bond(mol, atom_id);
    let has_multiple_bond = atom_explicit_pi_bond_count(mol, atom_id) > 0;

    if electrons == 0 {
        if noncyclic_pi_neighbor.is_some() {
            AromaticElectronDonorType::None
        } else if has_cyclic_pi_bond {
            AromaticElectronDonorType::One
        } else {
            AromaticElectronDonorType::None
        }
    } else if electrons == 1 {
        if let Some(neighbor) = noncyclic_pi_neighbor {
            if atom_is_more_electronegative_than(mol, neighbor, atom) {
                AromaticElectronDonorType::Vacant
            } else {
                AromaticElectronDonorType::One
            }
        } else if has_multiple_bond {
            AromaticElectronDonorType::One
        } else if atom.formal_charge == 1 {
            AromaticElectronDonorType::Vacant
        } else {
            AromaticElectronDonorType::None
        }
    } else {
        if noncyclic_pi_neighbor
            .is_some_and(|neighbor| atom_is_more_electronegative_than(mol, neighbor, atom))
        {
            electrons -= 1;
        }
        if electrons % 2 == 1 {
            AromaticElectronDonorType::One
        } else {
            AromaticElectronDonorType::Two
        }
    }
}

fn atom_noncyclic_pi_neighbor(mol: &Molecule, atom_id: AtomId) -> Option<AtomId> {
    let membership = mol
        .ring_membership()
        .expect("ring membership is computed before aromatic donor assignment");
    mol.incident_bonds(atom_id)
        .ok()?
        .find_map(|(bond_id, bond)| {
            (!membership.bond_in_ring(bond_id)
                && matches!(
                    bond.order,
                    BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple
                ))
            .then_some(bond.other_atom(atom_id))
        })
}

fn atom_has_cyclic_pi_bond(mol: &Molecule, atom_id: AtomId) -> bool {
    let membership = mol
        .ring_membership()
        .expect("ring membership is computed before aromatic donor assignment");
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .any(|(bond_id, bond)| {
            membership.bond_in_ring(bond_id)
                && matches!(
                    bond.order,
                    BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple
                )
        })
}

fn rdkit_rings_are_fused(left: &Ring, right: &Ring) -> bool {
    if left.bonds.len() > MAX_FUSED_AROMATIC_RING_SIZE
        || right.bonds.len() > MAX_FUSED_AROMATIC_RING_SIZE
    {
        return false;
    }
    left.bonds
        .iter()
        .filter(|bond| right.bonds.contains(bond))
        .count()
        == 1
}

fn rdkit_fused_ring_components(rings: &[Ring], candidates: &[usize]) -> Vec<Vec<usize>> {
    let mut components = (0..candidates.len()).collect::<Vec<_>>();
    for left in 0..candidates.len() {
        for right in (left + 1)..candidates.len() {
            if rdkit_rings_are_fused(&rings[candidates[left]], &rings[candidates[right]]) {
                union_components(&mut components, left, right);
            }
        }
    }
    let mut grouped = BTreeMap::<usize, Vec<usize>>::new();
    for (position, ring_index) in candidates.iter().copied().enumerate() {
        let root = find_component(&mut components, position);
        grouped.entry(root).or_default().push(ring_index);
    }
    grouped.into_values().collect()
}

fn apply_rdkit_huckel_to_fused_component(
    mol: &mut Molecule,
    rings: &[Ring],
    component: &[usize],
    donors: &[AromaticElectronDonorType],
) {
    let component_bonds = component
        .iter()
        .flat_map(|index| rings[*index].bonds.iter().copied())
        .collect::<BTreeSet<_>>();
    let mut done_bonds = BTreeSet::new();
    let max_subset_size = component.len().min(MAX_FUSED_AROMATIC_COMBINATION_RINGS);
    for subset_size in 1..=max_subset_size {
        if subset_size > 2 && component.len() > LARGE_FUSED_RING_SYSTEM_SEARCH_LIMIT {
            break;
        }
        for subset in connected_ring_subsets(rings, component, subset_size) {
            if !rdkit_ring_subset_is_connected(rings, &subset) {
                continue;
            }
            let mut atom_counts = BTreeMap::<AtomId, usize>::new();
            for ring_index in &subset {
                for atom in &rings[*ring_index].atoms {
                    *atom_counts.entry(*atom).or_default() += 1;
                }
            }
            let subset_donors = atom_counts
                .into_iter()
                .filter_map(|(atom, count)| (count <= 2).then_some(donors[atom.index()]))
                .collect::<Vec<_>>();
            if huckel_electron_count_for_donors(&subset_donors).is_none() {
                continue;
            }
            mark_rdkit_aromatic_subset(mol, rings, &subset, &mut done_bonds);
            if done_bonds.len() >= component_bonds.len() {
                return;
            }
        }
    }
}

fn rdkit_ring_subset_is_connected(rings: &[Ring], indexes: &[usize]) -> bool {
    let mut visited = BTreeSet::new();
    let mut stack = vec![indexes[0]];
    while let Some(index) = stack.pop() {
        if !visited.insert(index) {
            continue;
        }
        for other in indexes {
            if !visited.contains(other) && rdkit_rings_are_fused(&rings[index], &rings[*other]) {
                stack.push(*other);
            }
        }
    }
    visited.len() == indexes.len()
}

fn mark_rdkit_aromatic_subset(
    mol: &mut Molecule,
    rings: &[Ring],
    indexes: &[usize],
    done_bonds: &mut BTreeSet<BondId>,
) {
    let mut bond_counts = BTreeMap::<BondId, usize>::new();
    for index in indexes {
        for bond in &rings[*index].bonds {
            *bond_counts.entry(*bond).or_default() += 1;
        }
    }
    for (bond_id, count) in bond_counts {
        if count != 1 {
            continue;
        }
        done_bonds.insert(bond_id);
        let Some(bond) = mol.bonds[bond_id.index()].as_ref() else {
            continue;
        };
        let order = bond.order;
        let (left, right) = bond.endpoints();
        mol.set_bond_aromatic(bond_id, true);
        if matches!(order, BondOrder::Single | BondOrder::Double) {
            mol.set_atom_aromatic(left, true);
            mol.set_atom_aromatic(right, true);
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

fn rings_share_bond(left: &Ring, right: &Ring) -> bool {
    left.bonds.iter().any(|bond| right.bonds.contains(bond))
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
        && atom_rdkit_aromatic_total_valence(mol, atom_id, atom)
            > usize::from(charge_adjusted_default_valence)
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

fn atom_aromatic_candidate_degree(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> usize {
    let bonded_degree = mol
        .incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .filter(|(_, bond)| !matches!(bond.order, BondOrder::Zero | BondOrder::Dative))
        .count();
    bonded_degree
        .saturating_add(usize::from(atom.explicit_hydrogens))
        .saturating_add(aromaticity_implicit_hydrogen_count(mol, atom_id, atom))
}

fn aromaticity_implicit_hydrogen_count(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> usize {
    if let Some(hydrogens) = mol.implicit_hydrogens(atom_id).ok().flatten() {
        return usize::from(hydrogens);
    }
    if atom.no_implicit_hydrogens {
        return 0;
    }
    let Some(target) = aromaticity_valence_target(mol, atom_id, atom) else {
        return 0;
    };
    usize::from(target).saturating_sub(
        explicit_valence(mol, atom_id).saturating_add(usize::from(atom.explicit_hydrogens)),
    )
}

fn atom_rdkit_aromatic_total_valence(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> usize {
    explicit_valence(mol, atom_id)
        .saturating_add(usize::from(atom.explicit_hydrogens))
        .saturating_add(aromaticity_implicit_hydrogen_count(mol, atom_id, atom))
}

fn aromaticity_valence_target(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> Option<u8> {
    if mol.atom_is_aromatic(atom_id).ok().flatten() == Some(true) {
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
    let mut electrons = i16::from(default_valence)
        - i16::try_from(degree).expect("candidate degree is at most three")
        + lone_pair_electrons
        - radical_electrons;
    if electrons < 0 {
        return None;
    }
    if electrons > 1 && atom_explicit_unsaturation(mol, atom_id) > 1 {
        electrons = 1;
    }
    u8::try_from(electrons).ok()
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

fn huckel_electron_count_for_donors(donors: &[AromaticElectronDonorType]) -> Option<usize> {
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
        .map(|donor| usize::from(aromatic_donor_electron_range(*donor).0))
        .sum();
    let max_electrons = donors
        .iter()
        .map(|donor| usize::from(aromatic_donor_electron_range(*donor).1))
        .sum();
    huckel_electron_count_in_range(min_electrons, max_electrons)
}

fn huckel_electron_count_in_range(min_electrons: usize, max_electrons: usize) -> Option<usize> {
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

fn atom_explicit_unsaturation(mol: &Molecule, atom_id: AtomId) -> usize {
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

fn atom_is_more_electronegative_than(mol: &Molecule, left: AtomId, right: &Atom) -> bool {
    mol.atom(left).is_ok_and(|left| {
        rdkit_outer_electrons(left)
            .zip(rdkit_outer_electrons(right))
            .is_some_and(|(left_electrons, right_electrons)| {
                left_electrons > right_electrons
                    || left_electrons == right_electrons
                        && left.element.atomic_number() < right.element.atomic_number()
            })
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn huckel_electron_count_does_not_narrow_large_candidate_systems() {
        let donors = vec![AromaticElectronDonorType::One; 130];

        assert_eq!(huckel_electron_count_for_donors(&donors), Some(130));
    }

    #[test]
    fn fused_ten_electron_perimeter_preserves_explicit_aromatic_fusion_single() {
        let mut molecule =
            crate::small::SmallMolecule::from_smiles("On2c1-c(ccc2)ccn1").expect("parses");
        let valence = perceive_valence(molecule.graph_mut(), ValenceModel::RdkitLike);
        assert!(valence.is_ok(), "{:#?}", valence.issues);
        perceive_ring_set(molecule.graph_mut()).expect("rings");
        let protected_single = molecule
            .graph()
            .bonds()
            .find_map(|(bond_id, bond)| {
                (bond.order == BondOrder::Single
                    && molecule
                        .graph()
                        .atom(bond.a())
                        .is_ok_and(|atom| atom.aromatic)
                    && molecule
                        .graph()
                        .atom(bond.b())
                        .is_ok_and(|atom| atom.aromatic))
                .then_some(bond_id)
            })
            .expect("explicit aromatic fusion single");

        perceive_aromaticity(molecule.graph_mut(), AromaticityModel::RdkitLike)
            .expect("fused aromaticity");

        let protected = molecule
            .graph()
            .bond(protected_single)
            .expect("protected fusion bond");
        assert_eq!(protected.order, BondOrder::Single);
        assert!(!protected.aromatic);
        assert_eq!(
            molecule
                .graph()
                .atoms()
                .filter(|(_, atom)| atom.element.symbol() != "O" && atom.aromatic)
                .count(),
            9
        );
    }

    #[test]
    fn imported_aromatic_bonds_keep_implicit_hydrogen_nitrogen_pyrrole_like() {
        let input = "N2c1c(Nc3c2c6c(OS(=O)(=O)[O-])c7c(cccc7)c(OS(=O)(=O)[O-])c6cc3Cl)c4c(OS(=O)(=O)[O-])c5c(cccc5)c(OS(=O)(=O)[O-])c4cc1Cl";
        let mut molecule = crate::small::SmallMolecule::from_smiles(input).expect("dye parses");
        let valence = perceive_valence(molecule.graph_mut(), ValenceModel::RdkitLike);
        assert!(valence.is_ok(), "{:#?}", valence.issues);
        perceive_aromaticity(molecule.graph_mut(), AromaticityModel::RdkitLike)
            .expect("aromaticity");

        let nitrogens = molecule
            .graph()
            .atoms()
            .filter(|(_, atom)| atom.element.symbol() == "N")
            .map(|(_, atom)| (atom.aromatic, atom.implicit_hydrogens))
            .collect::<Vec<_>>();
        assert_eq!(nitrogens, vec![(false, Some(1)), (false, Some(1))]);
    }

    #[test]
    fn imported_aromatic_matching_limit_is_structured_and_transactional() {
        let mut molecule =
            crate::small::SmallMolecule::from_smiles("c1ccccc1").expect("aromatic benzene parses");
        let valence = perceive_valence(molecule.graph_mut(), ValenceModel::RdkitLike);
        assert!(valence.is_ok(), "{:#?}", valence.issues);
        let component = imported_aromatic_bond_components(molecule.graph())
            .into_iter()
            .next()
            .expect("imported aromatic component");
        let before = molecule.graph().clone();

        let error = try_kekulize_imported_component_with_limit(molecule.graph_mut(), &component, 0)
            .expect_err("zero matching budget should fail structurally");

        assert!(matches!(
            error,
            AromaticityError::ImportedAromaticKekulizationLimit {
                examined_states: 1,
                limit: 0,
                ..
            }
        ));
        assert_eq!(molecule.graph(), &before);
    }

    #[test]
    fn neutral_carbon_radical_can_complete_an_aromatic_sextet() {
        let mut molecule = crate::small::SmallMolecule::from_smiles("C1=CC(=CC=[C]1)N")
            .expect("aminophenyl radical parses");
        let valence = perceive_valence(molecule.graph_mut(), ValenceModel::RdkitLike);
        assert!(valence.is_ok(), "{:#?}", valence.issues);
        perceive_ring_set(molecule.graph_mut()).expect("ring perception");
        let radical = molecule
            .graph()
            .atoms()
            .find_map(|(id, atom)| atom.radical.is_some().then_some(id))
            .expect("radical carbon");
        let radical_atom = molecule.graph().atom(radical).expect("radical atom");
        let donor = rdkit_localized_atom_donor_type(molecule.graph(), radical, radical_atom);
        assert_eq!(donor, AromaticElectronDonorType::One);
        assert!(atom_is_rdkit_aromatic_candidate_for_donor(
            molecule.graph(),
            radical,
            radical_atom,
            donor,
            RdkitAromaticCandidateOptions::default(),
        ));

        perceive_aromaticity(molecule.graph_mut(), AromaticityModel::RdkitLike)
            .expect("aromaticity perception");

        assert_eq!(
            molecule
                .graph()
                .atoms()
                .filter(|(_, atom)| atom.element.symbol() == "C" && atom.aromatic)
                .count(),
            6
        );
    }

    #[test]
    fn charge_adjusted_candidate_valence_does_not_change_carbocation_electron_count() {
        for smiles in ["C1=C[C+]=CC(=C1)N", "C1=C[C+]=CC(=C1)C=O"] {
            let mut molecule =
                crate::small::SmallMolecule::from_smiles(smiles).expect("carbocation parses");
            crate::perception::sanitize(&mut molecule).expect("carbocation sanitizes");
            assert!(
                molecule.graph().atoms().all(|(_, atom)| !atom.aromatic),
                "{smiles}"
            );
        }
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
        assert!(atom_is_rdkit_aromatic_candidate_for_donor(
            &mol,
            carbonyl_carbon,
            atom,
            AromaticElectronDonorType::Any,
            RdkitAromaticCandidateOptions::default()
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
}
