use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::ops::ControlFlow;

use super::*;
use crate::core::*;

// RDKit's default model searches combinations of up to six candidate rings.
// Rings larger than 24 are evaluated individually, and components above 300
// rings are searched only as singles and pairs (RDKit 2026.03.3 Aromaticity.cpp).
const MAX_FUSED_AROMATIC_COMBINATION_RINGS: usize = 6;
const MAX_FUSED_AROMATIC_RING_SIZE: usize = 24;
const LARGE_FUSED_RING_SYSTEM_SEARCH_LIMIT: usize = 300;

#[cfg(test)]
mod fused_tests;

#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AromaticityError {
    UnsupportedElement(AtomId),
    RingPerception(RingPerceptionError),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AromaticElectronDonorType {
    Vacant,
    One,
    Two,
    None,
}

impl fmt::Display for AromaticityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedElement(id) => {
                write!(f, "unsupported aromaticity element at atom {id}")
            }
            Self::RingPerception(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for AromaticityError {}

/// Perceives aromatic atom and bond membership without changing localized bonds.
///
/// Reuses installed rings and implicit hydrogen counts. Missing rings are
/// perceived with default limits; missing hydrogen counts use the default
/// valence calculation without installing valence state. On failure, all
/// previously installed perception is preserved.
pub fn perceive_aromaticity(
    mol: &mut Molecule,
    model: AromaticityModel,
) -> std::result::Result<(), AromaticityError> {
    perceive_aromaticity_with_ring_options(mol, model, RingPerceptionOptions::default())
}

/// Perceives aromaticity, using `ring_options` if a ring basis is not installed.
///
/// An installed ring basis is reused, so these limits do not revalidate it.
/// Like [`perceive_aromaticity`], this operation is transactional and leaves
/// represented atom and bond chemistry unchanged.
pub fn perceive_aromaticity_with_ring_options(
    mol: &mut Molecule,
    model: AromaticityModel,
    ring_options: RingPerceptionOptions,
) -> std::result::Result<(), AromaticityError> {
    let previous = mol.perception().clone();
    if let Err(error) = perceive_aromaticity_with_ring_options_in_place(mol, model, ring_options) {
        mol.install_perception(previous)
            .expect("previous perception state must remain valid");
        return Err(error);
    }
    Ok(())
}

pub(crate) fn perceive_aromaticity_in_place(
    mol: &mut Molecule,
    model: AromaticityModel,
) -> std::result::Result<(), AromaticityError> {
    perceive_aromaticity_with_ring_options_in_place(mol, model, RingPerceptionOptions::default())
}

fn perceive_aromaticity_with_ring_options_in_place(
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
    let ring_set = match mol.ring_set() {
        Some(ring_set) => ring_set.clone(),
        None => perceive_ring_set_with_options(mol, ring_options)
            .map_err(AromaticityError::RingPerception)?,
    };
    assign_rdkit_like_localized_aromaticity(mol, &ring_set);
    Ok(())
}

fn assign_rdkit_like_localized_aromaticity(mol: &mut Molecule, ring_set: &RingSet) {
    mol.begin_aromaticity(AromaticityModel::RdkitLike);

    let mut donors = vec![AromaticElectronDonorType::None; mol.graph.atom_slot_count()];
    let mut atom_candidates = vec![false; mol.graph.atom_slot_count()];
    for (atom_id, atom) in mol.atoms() {
        let donor = rdkit_localized_atom_donor_type(mol, atom_id, atom);
        donors[atom_id.index()] = donor;
        atom_candidates[atom_id.index()] =
            atom_is_rdkit_aromatic_candidate_for_donor(mol, atom_id, atom, donor);
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
    let neighbors = rdkit_fused_ring_neighbors(ring_set.rings(), &candidates);
    let components = rdkit_fused_ring_components(&neighbors, &candidates);
    for component in components {
        apply_rdkit_huckel_to_fused_component(
            mol,
            ring_set.rings(),
            &neighbors,
            &component,
            &donors,
        );
    }
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
    let has_multiple_bond = mol.incident_bonds(atom_id).is_ok_and(|mut bonds| {
        bonds.any(|(_, bond)| {
            matches!(
                bond.order,
                BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple
            )
        })
    });

    if electrons == 0 {
        if noncyclic_pi_neighbor.is_some() {
            AromaticElectronDonorType::Vacant
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

fn rdkit_fused_ring_neighbors(rings: &[Ring], candidates: &[usize]) -> Vec<Vec<usize>> {
    let mut neighbors = vec![Vec::new(); rings.len()];
    for (position, &left) in candidates.iter().enumerate() {
        for &right in &candidates[position + 1..] {
            if rdkit_rings_are_fused(&rings[left], &rings[right]) {
                neighbors[left].push(right);
                neighbors[right].push(left);
            }
        }
    }
    neighbors
}

fn rdkit_fused_ring_components(neighbors: &[Vec<usize>], candidates: &[usize]) -> Vec<Vec<usize>> {
    let mut visited = vec![false; neighbors.len()];
    let mut components = Vec::new();
    for &root in candidates {
        if visited[root] {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![root];
        visited[root] = true;
        while let Some(ring) = stack.pop() {
            component.push(ring);
            for &neighbor in &neighbors[ring] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        components.push(component);
    }
    components
}

fn apply_rdkit_huckel_to_fused_component(
    mol: &mut Molecule,
    rings: &[Ring],
    neighbors: &[Vec<usize>],
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
        let result =
            visit_connected_ring_subsets(neighbors, component, subset_size, &mut |subset| {
                let mut atom_counts = BTreeMap::<AtomId, usize>::new();
                for ring_index in subset {
                    for atom in &rings[*ring_index].atoms {
                        *atom_counts.entry(*atom).or_default() += 1;
                    }
                }
                let subset_donors = atom_counts
                    .into_iter()
                    .filter_map(|(atom, count)| (count <= 2).then_some(donors[atom.index()]))
                    .collect::<Vec<_>>();
                if huckel_electron_count_for_donors(&subset_donors).is_none() {
                    return ControlFlow::Continue(());
                }
                mark_rdkit_aromatic_subset(mol, rings, subset, &mut done_bonds);
                if done_bonds.len() >= component_bonds.len() {
                    return ControlFlow::Break(());
                }
                ControlFlow::Continue(())
            });
        if result.is_break() {
            return;
        }
    }
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
        let Some(bond) = mol.graph.bonds[bond_id.index()].as_ref() else {
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

// Enumerate connected subsets directly, without materializing all combinations.
// Each subset belongs to its smallest ring index. At each depth, excluded rings
// record earlier sibling choices, so a subset is visited through only one path.
// The callback can stop as soon as all component bonds have been marked.
fn visit_connected_ring_subsets(
    neighbors: &[Vec<usize>],
    indexes: &[usize],
    subset_size: usize,
    visit: &mut impl FnMut(&[usize]) -> ControlFlow<()>,
) -> ControlFlow<()> {
    let mut excluded = vec![false; neighbors.len()];
    let mut current = Vec::with_capacity(subset_size);
    for &root in indexes {
        current.push(root);
        let frontier = neighbors[root]
            .iter()
            .copied()
            .filter(|&ring| ring > root)
            .collect();
        extend_connected_ring_subset(
            neighbors,
            subset_size,
            frontier,
            &mut current,
            &mut excluded,
            visit,
        )?;
        current.pop();
    }
    ControlFlow::Continue(())
}

fn extend_connected_ring_subset(
    neighbors: &[Vec<usize>],
    subset_size: usize,
    mut frontier: Vec<usize>,
    current: &mut Vec<usize>,
    excluded: &mut [bool],
    visit: &mut impl FnMut(&[usize]) -> ControlFlow<()>,
) -> ControlFlow<()> {
    if current.len() == subset_size {
        return visit(current);
    }
    let mut blocked = Vec::new();
    while let Some(ring) = frontier.pop() {
        excluded[ring] = true;
        blocked.push(ring);
        current.push(ring);
        let mut next_frontier = frontier.clone();
        for &neighbor in &neighbors[ring] {
            if neighbor > current[0] && !excluded[neighbor] && !next_frontier.contains(&neighbor) {
                next_frontier.push(neighbor);
            }
        }
        extend_connected_ring_subset(
            neighbors,
            subset_size,
            next_frontier,
            current,
            excluded,
            visit,
        )?;
        current.pop();
    }
    for ring in blocked {
        excluded[ring] = false;
    }
    ControlFlow::Continue(())
}

fn atom_is_rdkit_aromatic_candidate_for_donor(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    donor: AromaticElectronDonorType,
) -> bool {
    if matches!(donor, AromaticElectronDonorType::None) {
        return false;
    }
    let atomic_number = atom.element.atomic_number();
    if atomic_number > 18 && !matches!(atomic_number, 34 | 52) {
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
    if atom_explicit_unsaturation(mol, atom_id, atom) > 1
        && atom_explicit_pi_bond_count(mol, atom_id) > 1
    {
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
        .saturating_add(usize::from(atom.hydrogens.explicit_count()))
        .saturating_add(aromaticity_implicit_hydrogen_count(mol, atom_id, atom))
}

fn aromaticity_implicit_hydrogen_count(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> usize {
    if let Some(hydrogens) = mol.implicit_hydrogens(atom_id).ok().flatten() {
        return usize::from(hydrogens);
    }
    usize::from(super::valence::rdkit_implicit_hydrogen_count(
        mol, atom_id, atom,
    ))
}

fn atom_rdkit_aromatic_total_valence(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> usize {
    explicit_valence(mol, atom_id)
        .saturating_add(usize::from(atom.hydrogens.explicit_count()))
        .saturating_add(aromaticity_implicit_hydrogen_count(mol, atom_id, atom))
}

fn count_rdkit_like_atom_pi_electrons(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> Option<u8> {
    let default_valence = rdkit_default_valence(atom)?;
    let degree = atom_aromatic_candidate_degree(mol, atom_id, atom);
    if default_valence <= 1 || degree > 3 {
        return None;
    }

    let lone_pair_electrons = (i16::from(rdkit_outer_electrons(atom))
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
    if electrons > 1 && atom_explicit_unsaturation(mol, atom_id, atom) > 1 {
        electrons = 1;
    }
    u8::try_from(electrons).ok()
}

fn rdkit_outer_electrons(atom: &Atom) -> u8 {
    // RDKit 2026.03.3 atomic_data.cpp, indexed by atomic number. All elements
    // are needed here: an exocyclic neighbor can withdraw electrons even when
    // that neighbor is not itself eligible for aromaticity.
    const OUTER_ELECTRONS: [u8; 119] = [
        0, 1, 2, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
        2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 3,
        4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 5, 6, 7, 8, 1,
        2, 3, 4, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
        2, 2, 2,
    ];
    OUTER_ELECTRONS[usize::from(atom.element.atomic_number())]
}

fn aromatic_donor_electron_count(donor: AromaticElectronDonorType) -> usize {
    match donor {
        AromaticElectronDonorType::Vacant | AromaticElectronDonorType::None => 0,
        AromaticElectronDonorType::One => 1,
        AromaticElectronDonorType::Two => 2,
    }
}

fn huckel_electron_count_for_donors(donors: &[AromaticElectronDonorType]) -> Option<usize> {
    let electrons = donors
        .iter()
        .copied()
        .map(aromatic_donor_electron_count)
        .sum();
    (electrons == 2 || (electrons >= 6 && (electrons - 2) % 4 == 0)).then_some(electrons)
}

fn atom_explicit_pi_bond_count(mol: &Molecule, atom_id: AtomId) -> usize {
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .filter(|(_, bond)| matches!(bond.order, BondOrder::Double | BondOrder::Triple))
        .count()
}

fn atom_explicit_unsaturation(mol: &Molecule, atom_id: AtomId, atom: &Atom) -> usize {
    // RDKit uses explicit valence minus the unadjusted graph degree here.
    // Declared H contributes to valence; zero/dative bonds contribute to degree.
    // Negative values behave like zero for both comparisons that use this.
    let degree = mol.incident_bonds(atom_id).map_or(0, Iterator::count);
    explicit_valence(mol, atom_id)
        .saturating_add(usize::from(atom.hydrogens.explicit_count()))
        .saturating_sub(degree)
}

fn atom_is_more_electronegative_than(mol: &Molecule, left: AtomId, right: &Atom) -> bool {
    mol.atom(left).is_ok_and(|left| {
        let left_electrons = rdkit_outer_electrons(left);
        let right_electrons = rdkit_outer_electrons(right);
        left_electrons > right_electrons
            || left_electrons == right_electrons
                && left.element.atomic_number() < right.element.atomic_number()
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
        let source = "On2c1-c(ccc2)ccn1";
        let document = crate::io::parse_smiles_document(source).expect("parses");
        let (mut molecule, report) = document
            .interpret()
            .expect("interprets")
            .into_parts()
            .expect("one component");
        let explicit_single_offset = source.find('-').expect("explicit single marker") + 1;
        let protected_single = report
            .bond_mappings()
            .iter()
            .find(|mapping| mapping.source_offset() == explicit_single_offset)
            .map(|mapping| mapping.bond())
            .expect("explicit aromatic fusion single mapping");
        assert_eq!(
            molecule.bond(protected_single).unwrap().order,
            BondOrder::Single
        );
        let valence = perceive_valence(&mut molecule, ValenceModel::RdkitLike);
        assert!(valence.is_ok(), "{valence:#?}");
        perceive_ring_set(&mut molecule).expect("rings");

        perceive_aromaticity(&mut molecule, AromaticityModel::RdkitLike)
            .expect("fused aromaticity");

        let protected = molecule
            .bond(protected_single)
            .expect("protected fusion bond");
        assert_eq!(protected.order, BondOrder::Single);
        assert_eq!(molecule.bond_is_aromatic(protected_single), Ok(Some(false)));
        assert_eq!(
            molecule
                .atoms()
                .filter(|(atom_id, atom)| {
                    atom.element.symbol() != "O"
                        && molecule.atom_is_aromatic(*atom_id) == Ok(Some(true))
                })
                .count(),
            9
        );
    }

    #[test]
    fn canonical_localized_dye_assigns_nitrogen_hydrogens_before_aromaticity() {
        let input = "N2c1c(Nc3c2c6c(OS(=O)(=O)[O-])c7c(cccc7)c(OS(=O)(=O)[O-])c6cc3Cl)c4c(OS(=O)(=O)[O-])c5c(cccc5)c(OS(=O)(=O)[O-])c4cc1Cl";
        let mut molecule = crate::tests::read_smiles(input).expect("dye parses");
        assert!(!molecule.perception().has_aromaticity());
        let valence = perceive_valence(&mut molecule, ValenceModel::RdkitLike);
        assert!(valence.is_ok(), "{valence:#?}");
        let valence_nitrogens = molecule
            .atoms()
            .filter(|(_, atom)| atom.element.symbol() == "N")
            .map(|(atom_id, _)| molecule.implicit_hydrogens(atom_id))
            .collect::<Vec<_>>();
        assert_eq!(valence_nitrogens, vec![Ok(Some(1)), Ok(Some(1))]);
        assert!(!molecule.perception().has_aromaticity());
        perceive_aromaticity(&mut molecule, AromaticityModel::RdkitLike).expect("aromaticity");

        let nitrogens = molecule
            .atoms()
            .filter(|(_, atom)| atom.element.symbol() == "N")
            .map(|(atom_id, _)| {
                (
                    molecule.atom_is_aromatic(atom_id),
                    molecule.implicit_hydrogens(atom_id),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            nitrogens,
            vec![
                (Ok(Some(false)), Ok(Some(1))),
                (Ok(Some(false)), Ok(Some(1)))
            ]
        );
    }

    #[test]
    fn neutral_carbon_radical_can_complete_an_aromatic_sextet() {
        let mut molecule =
            crate::tests::read_smiles("C1=CC(=CC=[C]1)N").expect("aminophenyl radical parses");
        molecule
            .atom_mut(AtomId::new(5))
            .expect("source-selected radical carbon")
            .radical = Some(AtomRadical::Doublet);
        let valence = perceive_valence(&mut molecule, ValenceModel::RdkitLike);
        assert!(valence.is_ok(), "{valence:#?}");
        perceive_ring_set(&mut molecule).expect("ring perception");
        let radical = molecule
            .atoms()
            .find_map(|(id, atom)| atom.radical.is_some().then_some(id))
            .expect("radical carbon");
        let radical_atom = molecule.atom(radical).expect("radical atom");
        let donor = rdkit_localized_atom_donor_type(&molecule, radical, radical_atom);
        assert_eq!(donor, AromaticElectronDonorType::One);
        assert!(atom_is_rdkit_aromatic_candidate_for_donor(
            &molecule,
            radical,
            radical_atom,
            donor,
        ));

        perceive_aromaticity(&mut molecule, AromaticityModel::RdkitLike)
            .expect("aromaticity perception");

        assert_eq!(
            molecule
                .atoms()
                .filter(|(atom_id, atom)| {
                    atom.element.symbol() == "C"
                        && molecule.atom_is_aromatic(*atom_id) == Ok(Some(true))
                })
                .count(),
            6
        );
    }

    #[test]
    fn charge_adjusted_candidate_valence_does_not_change_carbocation_electron_count() {
        for smiles in ["C1=C[C+]=CC(=C1)N", "C1=C[C+]=CC(=C1)C=O"] {
            let mut molecule = crate::tests::read_smiles(smiles).expect("carbocation parses");
            molecule.perceive().expect("carbocation perceives");
            assert!(
                molecule
                    .atom_ids()
                    .all(|atom_id| { molecule.atom_is_aromatic(atom_id) == Ok(Some(false)) }),
                "{smiles}"
            );
        }
    }

    #[test]
    fn rdkit_candidate_elements_include_light_main_group_and_heavy_chalcogens() {
        // Full atom/bond masks checked against RDKit 2026.03.3 MolFromSmiles.
        for (source, aromatic) in [
            ("[Be-]1=CC=CC=C1", true),
            ("[Mg-]1=CC=CC=C1", true),
            ("[Al]1=CC=CC=C1", true),
            ("[SiH]1=CC=CC=C1", true),
            ("[SiH-]1C=CC=C1", true),
            ("[PH]1C=CC=C1", true),
            ("[Se]1C=CC=C1", true),
            ("[Te]1C=CC=C1", true),
            ("[GeH]1=CC=CC=C1", false),
            ("[AsH]1C=CC=C1", false),
        ] {
            let mut molecule = crate::tests::read_smiles(source).expect(source);
            molecule.perceive().expect(source);
            assert!(
                molecule
                    .atom_ids()
                    .all(|atom| { molecule.atom_is_aromatic(atom) == Ok(Some(aromatic)) }),
                "atom mask for {source}"
            );
            assert!(
                molecule
                    .bond_ids()
                    .all(|bond| { molecule.bond_is_aromatic(bond) == Ok(Some(aromatic)) }),
                "bond mask for {source}"
            );
        }
    }

    #[test]
    fn exocyclic_electron_withdrawal_uses_all_neighbor_elements() {
        // Exocyclic halogens and transition metals participate in RDKit's
        // outer-electron electronegativity comparison despite not being ring
        // candidates. An equally or less electronegative neighbor does not
        // remove the carbon's electron.
        for (source, aromatic) in [
            ("C1(=[Cl+])C=CC=CC=C1", true),
            ("C1(=[I+])C=CC=CC=C1", true),
            ("C1(=[Fe])C=CC=CC=C1", true),
            ("C1(=[Se])C=CC=CC=C1", true),
            ("C1(=O)C=CC=CC=C1", true),
            ("C1(=C)C=CC=CC=C1", false),
            ("C1(=[SiH2])C=CC=CC=C1", false),
        ] {
            let mut molecule = crate::tests::read_smiles(source).expect(source);
            molecule.perceive().expect(source);
            for atom in molecule.atom_ids() {
                assert_eq!(
                    molecule.atom_is_aromatic(atom),
                    Ok(Some(aromatic && atom.index() != 1)),
                    "{source}: {atom}"
                );
            }
            for (bond_id, bond) in molecule.bonds() {
                let (left, right) = bond.endpoints();
                assert_eq!(
                    molecule.bond_is_aromatic(bond_id),
                    Ok(Some(aromatic && left.index() != 1 && right.index() != 1)),
                    "{source}: {bond_id}"
                );
            }
        }
    }

    #[test]
    fn zero_electron_carbon_with_exocyclic_multiple_bond_is_vacant() {
        let mut molecule = crate::tests::read_smiles("C1(=O)C=CC=CC=C1").expect("tropone");
        molecule.atom_mut(AtomId::new(0)).unwrap().radical = Some(AtomRadical::Doublet);
        perceive_ring_set(&mut molecule).expect("rings");
        let atom = molecule.atom(AtomId::new(0)).unwrap();
        assert_eq!(
            count_rdkit_like_atom_pi_electrons(&molecule, AtomId::new(0), atom),
            Some(0)
        );
        assert_eq!(
            rdkit_localized_atom_donor_type(&molecule, AtomId::new(0), atom),
            AromaticElectronDonorType::Vacant
        );

        // Matches Kekulize(clearAromaticFlags=True), setting atom zero's
        // radical count to one, then SetAromaticity in RDKit 2026.03.3.
        perceive_aromaticity(&mut molecule, AromaticityModel::RdkitLike).expect("aromaticity");
        for atom in molecule.atom_ids() {
            assert_eq!(molecule.atom_is_aromatic(atom), Ok(Some(atom.index() != 1)));
        }
        for (bond_id, bond) in molecule.bonds() {
            let (left, right) = bond.endpoints();
            assert_eq!(
                molecule.bond_is_aromatic(bond_id),
                Ok(Some(left.index() != 1 && right.index() != 1))
            );
        }
    }

    #[test]
    fn zero_bond_contributes_to_rdkit_unsaturation_degree() {
        let mut molecule = crate::tests::read_smiles("C1=C=CC=C1").expect("cyclic allene");
        let helium = molecule
            .add_atom(Atom::new(Element::from_symbol("He").unwrap()))
            .unwrap();
        let zero = molecule
            .add_bond(AtomId::new(1), helium, BondOrder::Zero)
            .unwrap();
        molecule.perceive().expect("zero bond model");

        // RDKit's raw graph degree includes the zero-order bond. The resulting
        // explicit-valence-minus-degree is one, so the two double bonds do not
        // trigger the multiple-unsaturation exclusion and donate two electrons.
        for atom in molecule.atom_ids() {
            assert_eq!(molecule.atom_is_aromatic(atom), Ok(Some(atom != helium)));
        }
        for bond in molecule.bond_ids() {
            assert_eq!(molecule.bond_is_aromatic(bond), Ok(Some(bond != zero)));
        }
    }
}
