use std::collections::{BTreeMap, BTreeSet};

use crate::algorithms::{allowed_valences, rdkit_default_valence};
use crate::core::Molecule;
use crate::core::*;
use crate::io::MolWriteError;

use super::write::{
    smiles_atom, smiles_atom_requires_brackets, smiles_incident_bonds_for_style,
    smiles_ring_closures, validate_smiles_bracket_radical, validate_smiles_writeable,
    write_smiles_component, CanonicalAtomStyle, SmilesBondOrder, SmilesStereoWriteContext,
    SmilesWritePlan, StereoWriteMode,
};

mod labeling;
use super::emit::Emission;
use labeling::CanonicalOrder;

const MAX_INPUT_COMPLEXITY: usize = 50_000_000;
const MAX_GRAPH_SLOTS: usize = 2_000_000;

pub fn write_canonical_smiles(molecule: &Molecule) -> std::result::Result<String, MolWriteError> {
    write_canonical_smiles_with_limits(molecule, MAX_INPUT_COMPLEXITY, MAX_GRAPH_SLOTS)
}

fn write_canonical_smiles_with_limits(
    molecule: &Molecule,
    max_input_complexity: usize,
    max_graph_slots: usize,
) -> std::result::Result<String, MolWriteError> {
    write_canonical_emission_with_limits(molecule, max_input_complexity, max_graph_slots)?.render()
}

pub(super) fn write_canonical_emission(
    molecule: &Molecule,
) -> std::result::Result<Emission, MolWriteError> {
    write_canonical_emission_with_limits(molecule, MAX_INPUT_COMPLEXITY, MAX_GRAPH_SLOTS)
}

fn write_canonical_emission_with_limits(
    molecule: &Molecule,
    max_input_complexity: usize,
    max_graph_slots: usize,
) -> std::result::Result<Emission, MolWriteError> {
    // Check before cloning or ranking. A sparse graph can have many deleted
    // slots, so live atom counts alone do not bound
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
    // This quadratic input-size guard bounds admission to canonicalization,
    // not its actual search work. Labeling separately meters refinement and
    // search states. Keep the established 2*n*(n+2*m) admission threshold.
    let input_complexity = molecule
        .bond_count()
        .checked_mul(2)
        .and_then(|edges| edges.checked_add(atoms))
        .and_then(|visits| visits.checked_mul(atoms))
        .and_then(|visits| visits.checked_mul(2));
    if input_complexity.is_none_or(|score| score > max_input_complexity) {
        return Err(MolWriteError::resource_limit(format!(
            "canonical SMILES exceeds the input complexity limit ({max_input_complexity})"
        )));
    }
    validate_smiles_writeable(molecule, StereoWriteMode::Encode)?;
    let normalized = canonical_hydrogen_graph(molecule)?;
    let atom_style = canonical_atom_style(&normalized);
    let projected = canonical_projection_graph(&normalized, atom_style)?;
    let mol = &projected;
    let order = CanonicalOrder::new(mol, atom_style)?;
    let stereo = SmilesStereoWriteContext::new(mol, |atom| canonical_label(&order, atom))?;
    // A Molecule is nonempty and connected; topology writing owns component
    // ordering. Like RDKit, start this traversal at the least-ranked atom.
    let root = mol
        .atom_ids()
        .min_by_key(|atom| order.rank(*atom))
        .expect("canonical projection preserves a nonempty molecule");
    write_canonical_smiles_component(mol, root, &order, atom_style, &stereo)
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
    // Ranking and labeling must see the same isotope and hydrogen projection
    // as the output; only this private copy adopts the exported representation.
    let mut projected = mol.clone();
    for (atom_id, atom) in mol.atoms() {
        let (payload, _, inferred_hydrogens) =
            canonical_smiles_atom_representation(mol, atom_id, atom, atom_style)?;
        projected.graph.atoms[atom_id.index()] = Some(payload);
        projected.set_inferred_hydrogens(atom_id, inferred_hydrogens);
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

fn write_canonical_smiles_component(
    mol: &Molecule,
    root: AtomId,
    ranking: &CanonicalOrder,
    atom_style: CanonicalAtomStyle,
    stereo: &SmilesStereoWriteContext,
) -> std::result::Result<Emission, MolWriteError> {
    let rings = crate::algorithms::compute_ring_membership(mol);
    let plan = plan_canonical_smiles_component(mol, root, ranking, atom_style, &rings)?;
    write_smiles_component(
        mol,
        root,
        &plan,
        Some(stereo),
        atom_style,
        true,
        |_, children| {
            sort_canonical_smiles_neighbors(children, ranking, &rings);
            children.len().checked_sub(1)
        },
    )
}

fn plan_canonical_smiles_component(
    mol: &Molecule,
    root: AtomId,
    ranking: &CanonicalOrder,
    atom_style: CanonicalAtomStyle,
    rings: &RingMembership,
) -> std::result::Result<SmilesWritePlan, MolWriteError> {
    // Record back edges in DFS discovery order. In particular, closures
    // sharing their opening atom must follow when their far endpoints are
    // encountered, not the ranks of those endpoints.
    enum Action {
        Enter(AtomId, Option<BondId>),
        Edge(AtomId, BondId, SmilesBondOrder, AtomId),
        Exit(AtomId),
    }
    let mut visited = BTreeSet::new();
    let mut active = BTreeSet::new();
    let mut tree_bonds = BTreeSet::new();
    let mut ring_bonds = Vec::new();
    let mut actions = vec![Action::Enter(root, None)];
    while let Some(action) = actions.pop() {
        match action {
            Action::Enter(atom, parent) => {
                visited.insert(atom);
                active.insert(atom);
                let mut incident = smiles_incident_bonds_for_style(mol, atom, atom_style)?;
                sort_canonical_smiles_neighbors(&mut incident, ranking, rings);
                // Existing ancestors precede new continuations. Stable sorting
                // preserves the bond-order/rank ordering within each class.
                incident.sort_by_key(|(_, _, other)| !active.contains(other));
                actions.push(Action::Exit(atom));
                for (bond, order, other) in incident.into_iter().rev() {
                    if Some(bond) != parent {
                        actions.push(Action::Edge(atom, bond, order, other));
                    }
                }
            }
            Action::Edge(atom, bond, order, other) => {
                if !visited.contains(&other) {
                    tree_bonds.insert(bond);
                    actions.push(Action::Enter(other, Some(bond)));
                } else if active.contains(&other) {
                    ring_bonds.push((bond, other, atom, order));
                }
            }
            Action::Exit(atom) => {
                active.remove(&atom);
            }
        }
    }
    let closures = smiles_ring_closures(ring_bonds);

    Ok(SmilesWritePlan {
        roots: vec![root],
        tree_bonds,
        closures,
        subtree_sizes: BTreeMap::new(),
    })
}

fn sort_canonical_smiles_neighbors(
    neighbors: &mut [(BondId, SmilesBondOrder, AtomId)],
    ranking: &CanonicalOrder,
    rings: &RingMembership,
) {
    // Both tree planning and emission must use this order: acyclic branches
    // precede ring continuations, whose higher bond orders precede lower ones.
    // The last child is emitted as the main continuation.
    neighbors.sort_by_key(|(bond, order, atom)| {
        (
            rings.bond_in_ring(*bond),
            if rings.bond_in_ring(*bond) {
                reverse_bond_order_code(*order)
            } else {
                0
            },
            ranking.rank(*atom),
        )
    });
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
    let (atom, aromatic, inferred_hydrogens) =
        canonical_smiles_atom_representation(mol, atom_id, atom, atom_style)?;
    Ok(smiles_atom(&atom, aromatic, inferred_hydrogens))
}

fn canonical_smiles_atom_representation(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<(Atom, bool, u8), MolWriteError> {
    let aromatic = mol.atom_is_aromatic(atom_id).ok().flatten() == Some(true);
    let perceived_hydrogens = mol
        .inferred_hydrogens(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?;
    let inferred_hydrogens = perceived_hydrogens.unwrap_or(0);
    atom.hydrogens
        .specified_count()
        .checked_add(inferred_hydrogens)
        .ok_or_else(|| {
            MolWriteError::new("hydrogen count exceeds the SMILES representation limit")
        })?;
    let aromatic = aromatic && !matches!(atom_style, CanonicalAtomStyle::StoredKekule);
    let (mut payload, mut inferred_hydrogens) = canonical_smiles_atom_normalized(
        mol,
        atom_id,
        atom,
        aromatic,
        inferred_hydrogens,
        matches!(atom_style, CanonicalAtomStyle::StoredKekule),
    )?;
    if smiles_atom_requires_brackets(&payload, aromatic, inferred_hydrogens) {
        if atom.hydrogens.allows_inference() && perceived_hydrogens.is_none() {
            return Err(MolWriteError::new(format!(
                "canonical SMILES bracket atom {atom_id} requires installed hydrogen perception; perceive the molecule before writing"
            )));
        }
        payload.hydrogens = HydrogenDeclaration::Fixed(
            payload
                .hydrogens
                .specified_count()
                .saturating_add(inferred_hydrogens),
        );
        inferred_hydrogens = 0;
        validate_smiles_bracket_radical(mol, atom_id, &payload, inferred_hydrogens)?;
    }
    Ok((payload, aromatic, inferred_hydrogens))
}

fn canonical_smiles_atom_normalized(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    aromatic: bool,
    inferred_hydrogens: u8,
    stored_kekule: bool,
) -> std::result::Result<(Atom, u8), MolWriteError> {
    if canonical_smiles_should_bracket_metal_bound_hydrogens(
        mol,
        atom_id,
        atom,
        aromatic,
        inferred_hydrogens,
    )? {
        let mut normalized = atom.clone();
        normalized.hydrogens = HydrogenDeclaration::Fixed(
            atom.hydrogens
                .specified_count()
                .saturating_add(inferred_hydrogens),
        );
        return Ok((normalized, 0));
    }
    if canonical_smiles_should_bracket_metal_bound_zero_hydrogens(
        mol,
        atom_id,
        atom,
        inferred_hydrogens,
    )? {
        let mut normalized = atom.clone();
        normalized.hydrogens = HydrogenDeclaration::Fixed(atom.hydrogens.specified_count());
        return Ok((normalized, 0));
    }
    if canonical_smiles_can_use_organic_form(
        mol,
        atom_id,
        atom,
        aromatic,
        inferred_hydrogens,
        stored_kekule,
    )? {
        let mut normalized = atom.clone();
        normalized.hydrogens = HydrogenDeclaration::Infer { specified: 0 };
        return Ok((
            normalized,
            atom.hydrogens
                .specified_count()
                .saturating_add(inferred_hydrogens),
        ));
    }
    let mut normalized = atom.clone();
    if inferred_hydrogens > 0 {
        normalized.hydrogens = HydrogenDeclaration::Fixed(
            atom.hydrogens
                .specified_count()
                .saturating_add(inferred_hydrogens),
        );
    }
    Ok((normalized, 0))
}

fn canonical_smiles_should_bracket_metal_bound_hydrogens(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    aromatic: bool,
    inferred_hydrogens: u8,
) -> std::result::Result<bool, MolWriteError> {
    Ok(atom.formal_charge == 0
        && atom.radical.is_none()
        && atom.atom_map.is_none()
        && !aromatic
        && atom.hydrogens.allows_inference()
        && atom.hydrogens.specified_count() == 0
        && inferred_hydrogens > 0
        && matches!(atom.element.symbol(), "B" | "C" | "N" | "O" | "P" | "S")
        && atom_has_metal_neighbor(mol, atom_id)?)
}

fn canonical_smiles_should_bracket_metal_bound_zero_hydrogens(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    inferred_hydrogens: u8,
) -> std::result::Result<bool, MolWriteError> {
    Ok(atom.formal_charge == 0
        && atom.radical.is_none()
        && atom.atom_map.is_none()
        && atom.hydrogens.specified_count() == 0
        && inferred_hydrogens == 0
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
    inferred_hydrogens: u8,
    stored_kekule: bool,
) -> std::result::Result<bool, MolWriteError> {
    if atom.formal_charge != 0
        || atom.radical.is_some()
        || atom.atom_map.is_some()
        || (aromatic && atom.hydrogens.specified_count() > 0)
    {
        return Ok(false);
    }
    if !matches!(
        atom.element.symbol(),
        "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
    ) {
        return Ok(false);
    }
    if (!atom.hydrogens.allows_inference() || inferred_hydrogens == 0)
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
            .specified_count()
            .saturating_add(inferred_hydrogens);
        return Ok(bond_valence.saturating_add(total_hydrogens) == target);
    }
    let total_hydrogens = atom
        .hydrogens
        .specified_count()
        .saturating_add(inferred_hydrogens);
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
    fn canonical_input_complexity_accepts_the_boundary() {
        let mut molecule = crate::smiles::to_molecules("CC").unwrap().pop().unwrap();
        molecule.perceive().unwrap();
        assert_eq!(
            write_canonical_smiles_with_limits(&molecule, 16, usize::MAX).unwrap(),
            "CC"
        );
        let error = write_canonical_smiles_with_limits(&molecule, 15, usize::MAX).unwrap_err();
        assert_eq!(error.kind(), MolWriteErrorKind::ResourceLimit);
        assert!(error.to_string().contains("input complexity"));
    }

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
