use std::collections::{BTreeMap, BTreeSet};

use crate::algorithms::{
    allowed_valences, canonical_atom_ranking, ordered_atom_pair, rdkit_default_valence,
    CanonicalAtomRanking,
};
use crate::core::*;
use crate::io::MolWriteError;
use crate::small::model::SmallMolecule;

use super::write::{
    smiles_atom, smiles_bond_between, smiles_connected_components, smiles_incident_bonds_for_style,
    smiles_ring_number, validate_smiles_writeable, CanonicalAtomStyle, SmilesBondOrder,
    SmilesRingClosure, SmilesWritePlan, StereoWriteMode,
};

pub fn write_canonical_smiles(
    molecule: &SmallMolecule,
) -> std::result::Result<String, MolWriteError> {
    validate_smiles_writeable(molecule.graph(), StereoWriteMode::Ignore)?;
    let normalized = canonical_nonisomeric_graph(molecule.graph())?;
    let mol = &normalized;
    let ranking = canonical_atom_ranking(mol);
    let mut components = Vec::new();
    for component in smiles_connected_components(mol)? {
        let atom_style = canonical_component_atom_style(mol, &component)?;
        let mut candidates = Vec::new();
        for preference in [
            CanonicalBondTraversal::HighOrderFirst,
            CanonicalBondTraversal::LowOrderFirst,
        ] {
            candidates.extend(
                component
                    .iter()
                    .map(|root| {
                        write_canonical_smiles_component(
                            mol, *root, &ranking, preference, atom_style,
                        )
                    })
                    .collect::<std::result::Result<Vec<_>, _>>()?,
            );
        }
        candidates.sort_by_key(|candidate| canonical_smiles_candidate_key(candidate));
        candidates.dedup();
        if let Some(candidate) = candidates.into_iter().next() {
            components.push(candidate);
        }
    }
    components.sort();
    Ok(components.join("."))
}

fn canonical_nonisomeric_graph(mol: &Molecule) -> std::result::Result<Molecule, MolWriteError> {
    let mut normalized = mol.clone();
    let perception = normalized.perception.clone();
    let collapsible_hydrogens = mol
        .atoms()
        .filter_map(|(atom_id, atom)| {
            if atom.element.symbol() != "H" || atom.isotope.is_some() {
                return None;
            }
            let bonds = mol.incident_bonds(atom_id).ok()?.collect::<Vec<_>>();
            if bonds.len() != 1 || !matches!(bonds[0].1.order, BondOrder::Single) {
                return None;
            }
            let parent = bonds[0].1.other_atom(atom_id);
            mol.atom(parent)
                .is_ok_and(|parent_atom| parent_atom.element.symbol() != "H")
                .then_some((atom_id, parent))
        })
        .collect::<Vec<_>>();

    let mut implicit_by_parent = BTreeMap::new();
    for (hydrogen, parent) in collapsible_hydrogens {
        let implicit = mol
            .implicit_hydrogens(parent)
            .map_err(|error| MolWriteError::new(error.to_string()))?
            .unwrap_or(0);
        let count = implicit_by_parent.entry(parent).or_insert(implicit);
        *count = count.saturating_add(1);
        let parent_atom = normalized
            .atoms
            .get_mut(parent.index())
            .and_then(Option::as_mut)
            .ok_or_else(|| MolWriteError::new(format!("invalid hydrogen parent atom {parent}")))?;
        parent_atom.hydrogens = HydrogenDeclaration::Infer {
            explicit: parent_atom.hydrogens.explicit_count(),
        };
        normalized
            .delete_atom(hydrogen)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
    }
    normalized.perception = perception;
    for (parent, implicit) in implicit_by_parent {
        normalized.set_implicit_hydrogens(parent, implicit);
    }
    Ok(normalized)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CanonicalBondTraversal {
    HighOrderFirst,
    LowOrderFirst,
}

impl CanonicalBondTraversal {
    fn order_key(self, order: SmilesBondOrder) -> u8 {
        match self {
            Self::HighOrderFirst => reverse_bond_order_code(order),
            Self::LowOrderFirst => bond_order_code(order),
        }
    }
}

fn canonical_component_atom_style(
    mol: &Molecule,
    atom_ids: &[AtomId],
) -> std::result::Result<CanonicalAtomStyle, MolWriteError> {
    if canonical_component_has_aromatic_shorthand_sensitive_atom(mol, atom_ids)? {
        Ok(CanonicalAtomStyle::StoredKekule)
    } else {
        Ok(CanonicalAtomStyle::Aromatic)
    }
}

fn canonical_component_has_aromatic_shorthand_sensitive_atom(
    mol: &Molecule,
    atom_ids: &[AtomId],
) -> std::result::Result<bool, MolWriteError> {
    let atom_set = atom_ids.iter().copied().collect::<BTreeSet<_>>();
    let component_has_aromatic_atom = atom_ids
        .iter()
        .any(|atom_id| mol.atom_is_aromatic(*atom_id).ok().flatten() == Some(true));
    if !component_has_aromatic_atom {
        return Ok(false);
    }
    for atom_id in atom_ids {
        let atom = mol
            .atom(*atom_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        let aromatic = mol.atom_is_aromatic(*atom_id).ok().flatten() == Some(true);
        if aromatic && atom.formal_charge != 0 && matches!(atom.element.symbol(), "B" | "C") {
            return Ok(true);
        }
        if aromatic && atom_has_exocyclic_hetero_multiple_bond(mol, *atom_id, &atom_set)? {
            return Ok(true);
        }
        if aromatic {
            continue;
        }
        let mut aromatic_neighbors = 0usize;
        let mut pi_framework_neighbors = 0usize;
        let mut multiple_bond_to_non_aromatic_neighbor = false;
        for (_, bond) in mol
            .incident_bonds(*atom_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?
        {
            let neighbor_id = bond.other_atom(*atom_id);
            let neighbor = mol
                .atom(neighbor_id)
                .map_err(|error| MolWriteError::new(error.to_string()))?;
            let neighbor_aromatic = mol.atom_is_aromatic(neighbor_id).ok().flatten() == Some(true);
            if atom_set.contains(&neighbor_id) && neighbor_aromatic {
                aromatic_neighbors += 1;
            }
            if atom_set.contains(&neighbor_id)
                && matches!(neighbor.element.symbol(), "B" | "C" | "N" | "P" | "S")
            {
                pi_framework_neighbors += 1;
            }
            if matches!(bond.order, BondOrder::Double | BondOrder::Triple) && !neighbor_aromatic {
                multiple_bond_to_non_aromatic_neighbor = true;
            }
        }
        let unsupported_aromatic_ring_element = aromatic_neighbors > 0
            && mol
                .ring_membership()
                .is_some_and(|membership| membership.atom_in_ring(*atom_id))
            && !matches!(
                atom.element.symbol(),
                "B" | "C" | "N" | "O" | "P" | "S" | "Se" | "Te"
            );
        if unsupported_aromatic_ring_element {
            return Ok(true);
        }
        if atom.formal_charge == 0
            && (aromatic_neighbors > 0 || pi_framework_neighbors >= 3)
            && pi_framework_neighbors >= 2
            && multiple_bond_to_non_aromatic_neighbor
        {
            return Ok(true);
        }
    }
    Ok(false)
}

fn atom_has_exocyclic_hetero_multiple_bond(
    mol: &Molecule,
    atom_id: AtomId,
    atom_set: &BTreeSet<AtomId>,
) -> std::result::Result<bool, MolWriteError> {
    for (_, bond) in mol
        .incident_bonds(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?
    {
        if !matches!(bond.order, BondOrder::Double | BondOrder::Triple) {
            continue;
        }
        let neighbor_id = bond.other_atom(atom_id);
        let neighbor = mol
            .atom(neighbor_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        if !atom_set.contains(&neighbor_id)
            || mol.atom_is_aromatic(neighbor_id).ok().flatten() != Some(true)
                && !matches!(neighbor.element.symbol(), "B" | "C")
        {
            return Ok(true);
        }
    }
    Ok(false)
}

fn canonical_smiles_candidate_key(candidate: &str) -> (usize, usize, usize, String) {
    (
        candidate.matches('(').count(),
        explicit_ring_bond_marker_count(candidate),
        leading_ring_label_count(candidate),
        candidate.to_owned(),
    )
}

fn leading_ring_label_count(candidate: &str) -> usize {
    let bytes = candidate.as_bytes();
    let mut index = smiles_atom_token_end(candidate);
    let mut count = 0usize;
    while let Some(byte) = bytes.get(index) {
        if byte.is_ascii_digit() {
            count += 1;
            index += 1;
        } else if *byte == b'%' && bytes.get(index + 1).is_some_and(u8::is_ascii_digit) {
            count += 1;
            index += 3;
        } else {
            break;
        }
    }
    count
}

fn smiles_atom_token_end(candidate: &str) -> usize {
    let bytes = candidate.as_bytes();
    if bytes.first() == Some(&b'[') {
        return bytes
            .iter()
            .position(|byte| *byte == b']')
            .map(|index| index + 1)
            .unwrap_or(candidate.len());
    }
    if matches!(bytes.first(), Some(b'B' | b'C')) && matches!(bytes.get(1), Some(b'l' | b'r')) {
        2
    } else {
        bytes.first().map(|_| 1).unwrap_or(0)
    }
}

fn explicit_ring_bond_marker_count(candidate: &str) -> usize {
    let bytes = candidate.as_bytes();
    bytes
        .windows(2)
        .filter(|pair| matches!(pair[0], b'-' | b'=' | b'#' | b':') && pair[1].is_ascii_digit())
        .count()
        + bytes
            .windows(2)
            .filter(|pair| matches!(pair[0], b'-' | b'=' | b'#' | b':') && pair[1] == b'%')
            .count()
}

fn write_canonical_smiles_component(
    mol: &Molecule,
    root: AtomId,
    ranking: &CanonicalAtomRanking,
    preference: CanonicalBondTraversal,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<String, MolWriteError> {
    let plan = plan_canonical_smiles_component(mol, root, ranking, preference, atom_style)?;
    write_canonical_smiles_component_with_plan(mol, root, &plan, ranking, preference, atom_style)
}

fn plan_canonical_smiles_component(
    mol: &Molecule,
    root: AtomId,
    ranking: &CanonicalAtomRanking,
    preference: CanonicalBondTraversal,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<SmilesWritePlan, MolWriteError> {
    struct Frame {
        parent_bond: Option<BondId>,
        incident: Vec<(BondId, SmilesBondOrder, AtomId)>,
        next_edge: usize,
    }

    let mut visited = BTreeSet::<AtomId>::new();
    let mut tree_bonds = BTreeSet::<BondId>::new();
    let mut ring_bonds = BTreeMap::<BondId, (AtomId, AtomId, SmilesBondOrder)>::new();
    visited.insert(root);
    let mut stack = vec![Frame {
        parent_bond: None,
        incident: canonical_smiles_incident_bonds(mol, root, ranking, preference, atom_style)?,
        next_edge: 0,
    }];
    while let Some(frame) = stack.last_mut() {
        if frame.next_edge >= frame.incident.len() {
            stack.pop();
            continue;
        }
        let (bond_id, order, neighbor) = frame.incident[frame.next_edge];
        frame.next_edge += 1;
        if Some(bond_id) == frame.parent_bond {
            continue;
        }
        if visited.contains(&neighbor) {
            if !tree_bonds.contains(&bond_id) {
                let bond = mol
                    .bond(bond_id)
                    .map_err(|error| MolWriteError::new(error.to_string()))?;
                ring_bonds
                    .entry(bond_id)
                    .or_insert((bond.a(), bond.b(), order));
            }
            continue;
        }
        tree_bonds.insert(bond_id);
        visited.insert(neighbor);
        stack.push(Frame {
            parent_bond: Some(bond_id),
            incident: canonical_smiles_incident_bonds(
                mol, neighbor, ranking, preference, atom_style,
            )?,
            next_edge: 0,
        });
    }

    let mut ring_bonds = ring_bonds
        .into_iter()
        .map(|(bond_id, (a, b, order))| {
            let (first, second) = ordered_atom_pair(a, b);
            (bond_id, first, second, order)
        })
        .collect::<Vec<_>>();
    ring_bonds.sort_by_key(|(bond_id, first, second, order)| {
        (
            canonical_rank(ranking, *first),
            canonical_rank(ranking, *second),
            bond_order_code(*order),
            *first,
            *second,
            *bond_id,
        )
    });
    if ring_bonds.len() > 99 {
        return Err(MolWriteError::new(
            "SMILES writer supports at most 99 simultaneous ring closures",
        ));
    }

    let mut closures = BTreeMap::<AtomId, Vec<SmilesRingClosure>>::new();
    for (number, (bond_id, first, second, order)) in (1u64..).zip(ring_bonds) {
        closures.entry(first).or_default().push(SmilesRingClosure {
            bond: bond_id,
            number,
            order,
            other: second,
        });
        closures.entry(second).or_default().push(SmilesRingClosure {
            bond: bond_id,
            number,
            order,
            other: first,
        });
    }
    for (atom, closures) in &mut closures {
        closures.sort_by_key(|closure| {
            (
                canonical_rank(ranking, closure.other),
                bond_order_code(closure.order),
                closure.other,
                *atom,
            )
        });
    }

    Ok(SmilesWritePlan {
        roots: vec![root],
        tree_bonds,
        closures,
        subtree_sizes: BTreeMap::new(),
    })
}

fn write_canonical_smiles_component_with_plan(
    mol: &Molecule,
    root: AtomId,
    plan: &SmilesWritePlan,
    ranking: &CanonicalAtomRanking,
    preference: CanonicalBondTraversal,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<String, MolWriteError> {
    enum Action {
        Node {
            atom: AtomId,
            parent: Option<AtomId>,
        },
        Bond {
            order: SmilesBondOrder,
            left: AtomId,
            right: AtomId,
        },
        OpenBranch,
        CloseBranch,
    }

    let mut out = String::new();
    let mut actions = vec![Action::Node {
        atom: root,
        parent: None,
    }];
    while let Some(action) = actions.pop() {
        match action {
            Action::OpenBranch => out.push('('),
            Action::CloseBranch => out.push(')'),
            Action::Bond { order, left, right } => {
                out.push_str(smiles_bond_between(mol, order, left, right)?);
            }
            Action::Node { atom, parent } => {
                let atom_record = mol
                    .atom(atom)
                    .map_err(|error| MolWriteError::new(error.to_string()))?;
                out.push_str(&canonical_smiles_atom(mol, atom, atom_record, atom_style)?);
                if let Some(closures) = plan.closures.get(&atom) {
                    for closure in closures {
                        out.push_str(smiles_bond_between(
                            mol,
                            closure.order,
                            atom,
                            closure.other,
                        )?);
                        out.push_str(&smiles_ring_number(closure.number));
                    }
                }

                let mut children =
                    canonical_smiles_incident_bonds(mol, atom, ranking, preference, atom_style)?
                        .into_iter()
                        .filter(|(bond_id, _, neighbor)| {
                            plan.tree_bonds.contains(bond_id) && Some(*neighbor) != parent
                        })
                        .collect::<Vec<_>>();
                children.sort_by_key(|(bond_id, order, child)| {
                    (
                        !canonical_smiles_aromatic_continuation(mol, atom, *child, *order),
                        canonical_rank(ranking, *child),
                        canonical_smiles_atom_for_sort(mol, *child, atom_style),
                        preference.order_key(*order),
                        *child,
                        *bond_id,
                    )
                });
                let main_child = children.first().copied();
                if let Some((_, order, child)) = main_child {
                    actions.push(Action::Node {
                        atom: child,
                        parent: Some(atom),
                    });
                    actions.push(Action::Bond {
                        order,
                        left: atom,
                        right: child,
                    });
                }
                for (index, (_, order, child)) in children.into_iter().enumerate().rev() {
                    if index == 0 {
                        continue;
                    }
                    actions.push(Action::CloseBranch);
                    actions.push(Action::Node {
                        atom: child,
                        parent: Some(atom),
                    });
                    actions.push(Action::Bond {
                        order,
                        left: atom,
                        right: child,
                    });
                    actions.push(Action::OpenBranch);
                }
            }
        }
    }
    Ok(out)
}

fn canonical_smiles_aromatic_continuation(
    mol: &Molecule,
    left: AtomId,
    right: AtomId,
    order: SmilesBondOrder,
) -> bool {
    order == SmilesBondOrder::Aromatic
        && mol.atom_is_aromatic(left).ok().flatten() == Some(true)
        && mol.atom_is_aromatic(right).ok().flatten() == Some(true)
}

fn canonical_smiles_incident_bonds(
    mol: &Molecule,
    atom_id: AtomId,
    ranking: &CanonicalAtomRanking,
    preference: CanonicalBondTraversal,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<Vec<(BondId, SmilesBondOrder, AtomId)>, MolWriteError> {
    let mut incident = smiles_incident_bonds_for_style(mol, atom_id, atom_style)?;
    incident.sort_by_key(|(bond_id, order, atom)| {
        (
            canonical_rank(ranking, *atom),
            canonical_smiles_atom_for_sort(mol, *atom, atom_style),
            preference.order_key(*order),
            *atom,
            *bond_id,
        )
    });
    Ok(incident)
}

fn canonical_rank(ranking: &CanonicalAtomRanking, atom: AtomId) -> u32 {
    ranking
        .rank_of(atom)
        .expect("canonical ranking should cover every live atom")
}

fn bond_order_code(order: SmilesBondOrder) -> u8 {
    match order {
        SmilesBondOrder::Single => 1,
        SmilesBondOrder::Double => 2,
        SmilesBondOrder::Triple => 3,
        SmilesBondOrder::Aromatic => 5,
    }
}

fn reverse_bond_order_code(order: SmilesBondOrder) -> u8 {
    u8::MAX - bond_order_code(order)
}

fn canonical_smiles_atom(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<String, MolWriteError> {
    let mut normalized = atom.clone();
    let aromatic = mol.atom_is_aromatic(atom_id).ok().flatten() == Some(true);
    let mut implicit_hydrogens = mol
        .implicit_hydrogens(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?
        .unwrap_or(0);
    normalized.isotope = None;
    let represented_hydrogens = atom.hydrogens.explicit_count();
    if atom.isotope.is_some() && represented_hydrogens > 0 {
        implicit_hydrogens = represented_hydrogens.saturating_add(implicit_hydrogens);
        normalized.hydrogens = HydrogenDeclaration::Infer { explicit: 0 };
    }
    canonical_smiles_atom_normalized(
        mol,
        atom_id,
        &normalized,
        aromatic && !matches!(atom_style, CanonicalAtomStyle::StoredKekule),
        implicit_hydrogens,
        matches!(atom_style, CanonicalAtomStyle::StoredKekule),
    )
}

fn canonical_smiles_atom_normalized(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    aromatic: bool,
    implicit_hydrogens: u8,
    stored_kekule: bool,
) -> std::result::Result<String, MolWriteError> {
    if canonical_smiles_should_bracket_metal_bound_hydrogens(
        mol,
        atom_id,
        atom,
        aromatic,
        implicit_hydrogens,
    )? {
        let mut normalized = atom.clone();
        normalized.hydrogens = HydrogenDeclaration::Fixed(
            atom.hydrogens
                .explicit_count()
                .saturating_add(implicit_hydrogens),
        );
        return Ok(smiles_atom(&normalized, aromatic, 0));
    }
    if canonical_smiles_should_bracket_metal_bound_zero_hydrogens(
        mol,
        atom_id,
        atom,
        implicit_hydrogens,
    )? {
        let mut normalized = atom.clone();
        normalized.hydrogens = HydrogenDeclaration::Fixed(atom.hydrogens.explicit_count());
        return Ok(smiles_atom(&normalized, aromatic, 0));
    }
    if canonical_smiles_can_use_organic_form(
        mol,
        atom_id,
        atom,
        aromatic,
        implicit_hydrogens,
        stored_kekule,
    )? {
        let mut normalized = atom.clone();
        normalized.hydrogens = HydrogenDeclaration::Infer { explicit: 0 };
        return Ok(smiles_atom(&normalized, aromatic, implicit_hydrogens));
    }
    let mut normalized = atom.clone();
    if implicit_hydrogens > 0 {
        normalized.hydrogens = HydrogenDeclaration::Fixed(
            atom.hydrogens
                .explicit_count()
                .saturating_add(implicit_hydrogens),
        );
    }
    Ok(smiles_atom(&normalized, aromatic, 0))
}

fn canonical_smiles_should_bracket_metal_bound_hydrogens(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    aromatic: bool,
    implicit_hydrogens: u8,
) -> std::result::Result<bool, MolWriteError> {
    Ok(atom.formal_charge == 0
        && atom.radical.is_none()
        && atom.atom_map.is_none()
        && !aromatic
        && atom.hydrogens.allows_implicit()
        && atom.hydrogens.explicit_count() == 0
        && implicit_hydrogens > 0
        && matches!(atom.element.symbol(), "B" | "C" | "N" | "O" | "P" | "S")
        && atom_has_metal_neighbor(mol, atom_id)?)
}

fn canonical_smiles_should_bracket_metal_bound_zero_hydrogens(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    implicit_hydrogens: u8,
) -> std::result::Result<bool, MolWriteError> {
    Ok(atom.formal_charge == 0
        && atom.radical.is_none()
        && atom.atom_map.is_none()
        && atom.hydrogens.explicit_count() == 0
        && implicit_hydrogens == 0
        && matches!(
            atom.element.symbol(),
            "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
        )
        && atom_has_metal_neighbor(mol, atom_id)?)
}

fn canonical_smiles_atom_for_sort(
    mol: &Molecule,
    atom_id: AtomId,
    atom_style: CanonicalAtomStyle,
) -> String {
    let atom = mol
        .atom(atom_id)
        .expect("canonical atom sort should only use live atoms");
    canonical_smiles_atom(mol, atom_id, atom, atom_style)
        .expect("canonical atom sort should be encodable")
}

fn canonical_smiles_can_use_organic_form(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    aromatic: bool,
    implicit_hydrogens: u8,
    stored_kekule: bool,
) -> std::result::Result<bool, MolWriteError> {
    if atom.formal_charge != 0
        || atom.radical.is_some()
        || atom.atom_map.is_some()
        || (aromatic && atom.hydrogens.explicit_count() > 0)
    {
        return Ok(false);
    }
    if !matches!(
        atom.element.symbol(),
        "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
    ) {
        return Ok(false);
    }
    if (!atom.hydrogens.allows_implicit() || implicit_hydrogens == 0)
        && atom_has_metal_neighbor(mol, atom_id)?
    {
        return Ok(false);
    }
    let bond_valence = smiles_bond_valence_sum(mol, atom_id, stored_kekule)?;
    if aromatic {
        let Some(target) = canonical_organic_valence_target(atom, true) else {
            return Ok(false);
        };
        let total_hydrogens = atom
            .hydrogens
            .explicit_count()
            .saturating_add(implicit_hydrogens);
        return Ok(bond_valence.saturating_add(total_hydrogens) == target);
    }
    let total_hydrogens = atom
        .hydrogens
        .explicit_count()
        .saturating_add(implicit_hydrogens);
    let occupied_valence = bond_valence.saturating_add(total_hydrogens);
    Ok(
        allowed_valences(atom).is_some_and(|allowed| allowed.contains(&occupied_valence))
            && (total_hydrogens == 0 || rdkit_default_valence(atom) == Some(occupied_valence)),
    )
}

fn atom_has_metal_neighbor(
    mol: &Molecule,
    atom_id: AtomId,
) -> std::result::Result<bool, MolWriteError> {
    for (_, bond) in mol
        .incident_bonds(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?
    {
        let neighbor_id = bond.other_atom(atom_id);
        let neighbor = mol
            .atom(neighbor_id)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        if is_smiles_metal_like(neighbor.element.symbol()) {
            return Ok(true);
        }
    }
    Ok(false)
}

fn is_smiles_metal_like(symbol: &str) -> bool {
    matches!(
        symbol,
        "Li" | "Na"
            | "K"
            | "Rb"
            | "Cs"
            | "Fr"
            | "Be"
            | "Mg"
            | "Ca"
            | "Sr"
            | "Ba"
            | "Ra"
            | "Al"
            | "Ge"
            | "Ga"
            | "In"
            | "Tl"
            | "Sn"
            | "Pb"
            | "Sb"
            | "Bi"
            | "Po"
            | "Sc"
            | "Ti"
            | "V"
            | "Cr"
            | "Mn"
            | "Fe"
            | "Co"
            | "Ni"
            | "Cu"
            | "Zn"
            | "Y"
            | "Zr"
            | "Nb"
            | "Mo"
            | "Tc"
            | "Ru"
            | "Rh"
            | "Pd"
            | "Ag"
            | "Cd"
            | "La"
            | "Ce"
            | "Pr"
            | "Nd"
            | "Sm"
            | "Eu"
            | "Gd"
            | "Tb"
            | "Dy"
            | "Ho"
            | "Er"
            | "Tm"
            | "Yb"
            | "Lu"
            | "Ac"
            | "Th"
            | "Pa"
            | "U"
            | "Np"
            | "Pu"
            | "Am"
            | "Cm"
            | "Bk"
            | "Cf"
            | "Es"
            | "Fm"
            | "Md"
            | "No"
            | "Lr"
            | "Hf"
            | "Ta"
            | "W"
            | "Re"
            | "Os"
            | "Ir"
            | "Pt"
            | "Au"
            | "Hg"
    )
}

fn canonical_organic_valence_target(atom: &Atom, aromatic: bool) -> Option<u8> {
    match (atom.element.symbol(), aromatic) {
        ("B", false) => Some(3),
        ("C", false) => Some(4),
        ("N", false) | ("P", false) => Some(3),
        ("O", false) | ("S", false) => Some(2),
        ("F" | "Cl" | "Br" | "I", false) => Some(1),
        ("B" | "C", true) => Some(3),
        ("N" | "O" | "S" | "P", true) => Some(2),
        _ => None,
    }
}

fn smiles_bond_valence_sum(
    mol: &Molecule,
    atom_id: AtomId,
    stored_kekule: bool,
) -> std::result::Result<u8, MolWriteError> {
    mol.incident_bonds(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?
        .map(|(bond_id, bond)| {
            if mol.bond_is_aromatic(bond_id).ok().flatten() == Some(true) && !stored_kekule {
                return Ok(1);
            }
            Ok(match bond.order {
                BondOrder::Zero | BondOrder::Dative => 0,
                BondOrder::Single => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
                BondOrder::Quadruple => 4,
            })
        })
        .try_fold(0u8, |sum, value: std::result::Result<u8, MolWriteError>| {
            Ok(sum.saturating_add(value?))
        })
}
