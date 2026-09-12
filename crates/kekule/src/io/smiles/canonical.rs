use std::collections::{BTreeMap, BTreeSet};

use crate::algorithms::{
    allowed_valences, canonical_atom_ranking, rdkit_default_valence, CanonicalAtomRanking,
};
use crate::core::Molecule;
use crate::core::*;
use crate::io::MolWriteError;

use super::write::{
    collect_smiles_tree, smiles_atom, smiles_atom_requires_brackets, smiles_connected_components,
    smiles_incident_bonds_for_style, smiles_ring_closures, validate_smiles_writeable,
    write_smiles_component, CanonicalAtomStyle, SmilesBondOrder, SmilesStereoWriteContext,
    SmilesWritePlan, StereoWriteMode,
};

mod labeling;
use labeling::CanonicalOrder;

const MAX_CANDIDATE_VISITS: usize = 50_000_000;
const MAX_GRAPH_SLOTS: usize = 2_000_000;

pub fn write_canonical_smiles(molecule: &Molecule) -> std::result::Result<String, MolWriteError> {
    write_canonical_smiles_with_limits(molecule, MAX_CANDIDATE_VISITS, MAX_GRAPH_SLOTS)
}

fn write_canonical_smiles_with_limits(
    molecule: &Molecule,
    max_candidate_visits: usize,
    max_graph_slots: usize,
) -> std::result::Result<String, MolWriteError> {
    // Check before cloning, ranking or constructing any candidate. A sparse
    // graph can have many deleted slots, so live atom counts alone do not bound
    // the dense scratch arrays used by ranking and labeling.
    let slots = molecule
        .graph
        .atom_slot_count()
        .checked_add(molecule.graph.bond_slot_count());
    if slots.is_none_or(|slots| slots > max_graph_slots) {
        return Err(MolWriteError::resource_limit(format!(
            "canonical SMILES exceeds the graph storage limit ({max_graph_slots} atom/bond slots)"
        )));
    }
    let atoms = molecule.atom_count();
    let candidate_visits = molecule
        .bond_count()
        .checked_mul(2)
        .and_then(|edges| edges.checked_add(atoms))
        .and_then(|visits| visits.checked_mul(atoms))
        .and_then(|visits| visits.checked_mul(2));
    if candidate_visits.is_none_or(|visits| visits > max_candidate_visits) {
        return Err(MolWriteError::resource_limit(format!(
            "canonical SMILES exceeds the candidate traversal limit ({max_candidate_visits} atom/edge visits)"
        )));
    }
    validate_smiles_writeable(molecule, StereoWriteMode::Encode)?;
    let normalized = canonical_hydrogen_graph(molecule)?;
    let mol = &normalized;
    let mut components = Vec::new();
    for component in smiles_connected_components(mol)? {
        let atom_style = canonical_atom_style(mol);
        let projected = canonical_projection_graph(mol, atom_style)?;
        let mol = &projected;
        let ranking = canonical_atom_ranking(mol);
        let order = CanonicalOrder::new(mol, &ranking, atom_style)?;
        let stereo = SmilesStereoWriteContext::new(mol, |atom| canonical_label(&order, atom))?;
        let mut best = None;
        for preference in [
            CanonicalBondTraversal::HighOrderFirst,
            CanonicalBondTraversal::LowOrderFirst,
        ] {
            for root in &component {
                let candidate = write_canonical_smiles_component(
                    mol, *root, &order, preference, atom_style, &stereo,
                )?;
                let key = canonical_smiles_candidate_key(candidate);
                if best.as_ref().is_none_or(|best| key < *best) {
                    best = Some(key);
                }
            }
        }
        if let Some((_, _, _, candidate)) = best {
            components.push(candidate);
        }
    }
    components.sort();
    Ok(components.join("."))
}

fn canonical_hydrogen_graph(mol: &Molecule) -> std::result::Result<Molecule, MolWriteError> {
    let mut normalized = mol.clone();
    // Metadata does not affect a molecular identifier. Use the common hydrogen
    // transform so carrier remapping and count reconstruction have one owner.
    normalized.clear_properties();
    normalized.remove_hydrogens().map_err(|error| {
        MolWriteError::new(format!(
            "canonical SMILES hydrogen normalization requires known hydrogen perception: {error}"
        ))
    })?;
    restore_projection_aromaticity(&mut normalized, mol.perception());
    Ok(normalized)
}

fn canonical_projection_graph(
    mol: &Molecule,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<Molecule, MolWriteError> {
    // Ranking must see the same isotope and hydrogen projection as the output.
    // Keep the general atom-ranking API sensitive to authoritative chemistry;
    // only this private copy adopts the exported atom representation.
    let mut projected = mol.clone();
    for (atom_id, atom) in mol.atoms() {
        let (payload, _, implicit_hydrogens) =
            canonical_smiles_atom_representation(mol, atom_id, atom, atom_style)?;
        projected.graph.atoms[atom_id.index()] = Some(payload);
        projected.set_implicit_hydrogens(atom_id, implicit_hydrogens);
    }
    restore_projection_aromaticity(&mut projected, mol.perception());
    Ok(projected)
}

fn restore_projection_aromaticity(projected: &mut Molecule, original: &Perception) {
    // Changing the storage location of an unchanged hydrogen count invalidates
    // aromaticity through the ordinary edit helper. This private export copy
    // retains the original perceived aromatic system: isotope/H projection
    // neither changes its valence nor selects another aromaticity model.
    if let Some(aromaticity) = original.aromaticity_state() {
        projected.begin_aromaticity(aromaticity.model());
        for atom in aromaticity.atoms() {
            if projected.atom(atom).is_ok() {
                projected.set_atom_aromatic(atom, true);
            }
        }
        for bond in aromaticity.bonds() {
            if projected.bond(bond).is_ok() {
                projected.set_bond_aromatic(bond, true);
            }
        }
    }
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

fn canonical_atom_style(mol: &Molecule) -> CanonicalAtomStyle {
    if mol
        .stereo_elements()
        .any(|(_, element)| match &element.kind {
            StereoElementKind::DoubleBond(value) => [value.left, value.right]
                .iter()
                .any(|atom| mol.atom_is_aromatic(*atom).ok().flatten() == Some(true)),
            StereoElementKind::Tetrahedral(value) => {
                mol.atom_is_aromatic(value.center).ok().flatten() == Some(true)
            }
            StereoElementKind::Axis(_) => false,
        })
    {
        CanonicalAtomStyle::StoredKekule
    } else {
        CanonicalAtomStyle::Aromatic
    }
}

fn canonical_smiles_candidate_key(candidate: String) -> (usize, usize, usize, String) {
    (
        candidate.matches('(').count(),
        explicit_ring_bond_marker_count(&candidate),
        leading_ring_label_count(&candidate),
        candidate,
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
    ranking: &CanonicalOrder,
    preference: CanonicalBondTraversal,
    atom_style: CanonicalAtomStyle,
    stereo: &SmilesStereoWriteContext,
) -> std::result::Result<String, MolWriteError> {
    let plan = plan_canonical_smiles_component(mol, root, ranking, preference, atom_style)?;
    write_canonical_smiles_component_with_plan(
        mol, root, &plan, ranking, preference, atom_style, stereo,
    )
}

fn plan_canonical_smiles_component(
    mol: &Molecule,
    root: AtomId,
    ranking: &CanonicalOrder,
    preference: CanonicalBondTraversal,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<SmilesWritePlan, MolWriteError> {
    let mut visited = BTreeSet::<AtomId>::new();
    let mut tree_bonds = BTreeSet::<BondId>::new();
    let mut ring_bonds = BTreeMap::<BondId, (AtomId, AtomId, SmilesBondOrder)>::new();
    collect_smiles_tree(
        mol,
        root,
        None,
        &mut visited,
        &mut tree_bonds,
        &mut ring_bonds,
        |atom| canonical_smiles_incident_bonds(mol, atom, ranking, preference, atom_style),
    )?;

    let mut ring_bonds = ring_bonds
        .into_iter()
        .map(|(bond_id, (a, b, order))| {
            let (first, second) = if ranking.rank(a) < ranking.rank(b) {
                (a, b)
            } else {
                (b, a)
            };
            (bond_id, first, second, order)
        })
        .collect::<Vec<_>>();
    ring_bonds.sort_by_key(|(_, first, second, order)| {
        (
            canonical_rank(ranking, *first),
            canonical_rank(ranking, *second),
            bond_order_code(*order),
            canonical_label(ranking, *first),
            canonical_label(ranking, *second),
        )
    });
    let mut closures = smiles_ring_closures(ring_bonds);
    for closures in closures.values_mut() {
        closures.sort_by_key(|closure| {
            (
                canonical_rank(ranking, closure.other),
                bond_order_code(closure.order),
                canonical_label(ranking, closure.other),
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
    ranking: &CanonicalOrder,
    preference: CanonicalBondTraversal,
    atom_style: CanonicalAtomStyle,
    stereo: &SmilesStereoWriteContext,
) -> std::result::Result<String, MolWriteError> {
    write_smiles_component(
        mol,
        root,
        plan,
        Some(stereo),
        atom_style,
        |atom, children| {
            children.sort_by_key(|(_, order, child)| {
                (
                    !canonical_smiles_aromatic_continuation(mol, atom, *child, *order),
                    canonical_rank(ranking, *child),
                    canonical_smiles_atom_for_sort(mol, *child, atom_style),
                    preference.order_key(*order),
                    canonical_label(ranking, *child),
                )
            });
            (!children.is_empty()).then_some(0)
        },
    )
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
    ranking: &CanonicalOrder,
    preference: CanonicalBondTraversal,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<Vec<(BondId, SmilesBondOrder, AtomId)>, MolWriteError> {
    let mut incident = smiles_incident_bonds_for_style(mol, atom_id, atom_style)?;
    incident.sort_by_key(|(_, order, atom)| {
        (
            canonical_rank(ranking, *atom),
            canonical_smiles_atom_for_sort(mol, *atom, atom_style),
            preference.order_key(*order),
            canonical_label(ranking, *atom),
        )
    });
    Ok(incident)
}

fn canonical_rank(ranking: &CanonicalOrder, atom: AtomId) -> u32 {
    ranking.rank(atom).0
}

fn canonical_label(ranking: &CanonicalOrder, atom: AtomId) -> usize {
    ranking.rank(atom).1
}

fn bond_order_code(order: SmilesBondOrder) -> u8 {
    match order {
        SmilesBondOrder::Single => 1,
        SmilesBondOrder::Double => 2,
        SmilesBondOrder::Triple => 3,
        SmilesBondOrder::Quadruple => 4,
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
    let (atom, aromatic, implicit_hydrogens) =
        canonical_smiles_atom_representation(mol, atom_id, atom, atom_style)?;
    Ok(smiles_atom(&atom, aromatic, implicit_hydrogens))
}

fn canonical_smiles_atom_representation(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<(Atom, bool, u8), MolWriteError> {
    let normalized = atom.clone();
    let aromatic = mol.atom_is_aromatic(atom_id).ok().flatten() == Some(true);
    let perceived_hydrogens = mol
        .implicit_hydrogens(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?;
    let implicit_hydrogens = perceived_hydrogens.unwrap_or(0);
    atom.hydrogens
        .explicit_count()
        .checked_add(implicit_hydrogens)
        .ok_or_else(|| {
            MolWriteError::new("hydrogen count exceeds the SMILES representation limit")
        })?;
    let aromatic = aromatic && !matches!(atom_style, CanonicalAtomStyle::StoredKekule);
    let (mut payload, mut implicit_hydrogens) = canonical_smiles_atom_normalized(
        mol,
        atom_id,
        &normalized,
        aromatic,
        implicit_hydrogens,
        matches!(atom_style, CanonicalAtomStyle::StoredKekule),
    )?;
    if smiles_atom_requires_brackets(&payload, aromatic, implicit_hydrogens) {
        if atom.hydrogens.allows_implicit() && perceived_hydrogens.is_none() {
            return Err(MolWriteError::new(format!(
                "canonical SMILES bracket atom {atom_id} requires installed hydrogen perception; perceive the molecule before writing"
            )));
        }
        payload.hydrogens = HydrogenDeclaration::Fixed(
            payload
                .hydrogens
                .explicit_count()
                .saturating_add(implicit_hydrogens),
        );
        implicit_hydrogens = 0;
    }
    Ok((payload, aromatic, implicit_hydrogens))
}

fn canonical_smiles_atom_normalized(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    aromatic: bool,
    implicit_hydrogens: u8,
    stored_kekule: bool,
) -> std::result::Result<(Atom, u8), MolWriteError> {
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
        return Ok((normalized, 0));
    }
    if canonical_smiles_should_bracket_metal_bound_zero_hydrogens(
        mol,
        atom_id,
        atom,
        implicit_hydrogens,
    )? {
        let mut normalized = atom.clone();
        normalized.hydrogens = HydrogenDeclaration::Fixed(atom.hydrogens.explicit_count());
        return Ok((normalized, 0));
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
        return Ok((
            normalized,
            atom.hydrogens
                .explicit_count()
                .saturating_add(implicit_hydrogens),
        ));
    }
    let mut normalized = atom.clone();
    if implicit_hydrogens > 0 {
        normalized.hydrogens = HydrogenDeclaration::Fixed(
            atom.hydrogens
                .explicit_count()
                .saturating_add(implicit_hydrogens),
        );
    }
    Ok((normalized, 0))
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

#[cfg(test)]
mod resource_tests {
    use super::*;
    use crate::io::MolWriteErrorKind;

    #[test]
    fn canonical_export_checks_deleted_slots_before_dense_scratch_allocation() {
        let mut editor = MoleculeEditor::new();
        let carbon = Atom::new(Element::from_symbol("C").unwrap());
        editor.add_atom(carbon.clone()).unwrap();
        for _ in 0..8 {
            let deleted = editor.add_atom(carbon.clone()).unwrap();
            editor.delete_atom(deleted).unwrap();
        }
        let molecule = editor.finish().unwrap();
        assert_eq!(molecule.atom_count(), 1);
        let error = write_canonical_smiles_with_limits(&molecule, usize::MAX, 4).unwrap_err();
        assert_eq!(error.kind(), MolWriteErrorKind::ResourceLimit);
        assert!(error.to_string().contains("graph storage"));
    }
}
