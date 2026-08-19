use std::collections::{BTreeMap, BTreeSet};

use crate::algorithms::ordered_atom_pair;
use crate::core::*;
use crate::io::MolWriteError;
use crate::small::model::SmallMolecule;

use super::parse::SmilesDirectionToken;

pub fn write_smiles(molecule: &SmallMolecule) -> std::result::Result<String, MolWriteError> {
    let mol = molecule.graph();
    let plan = plan_smiles_write(mol, StereoWriteMode::Reject)?;
    let mut parts = Vec::new();
    for start in &plan.roots {
        parts.push(write_smiles_component(
            mol,
            *start,
            None,
            &plan,
            None,
            CanonicalAtomStyle::Aromatic,
        )?);
    }
    Ok(parts.join("."))
}

pub fn write_isomeric_smiles(
    molecule: &SmallMolecule,
) -> std::result::Result<String, MolWriteError> {
    let mol = molecule.graph();
    let plan = plan_smiles_write(mol, StereoWriteMode::Encode)?;
    let stereo = SmilesStereoWriteContext::new(mol)?;
    let component_styles = smiles_connected_components(mol)?
        .into_iter()
        .map(|component| {
            (
                component.into_iter().collect::<BTreeSet<_>>(),
                CanonicalAtomStyle::StoredKekule,
            )
        })
        .collect::<Vec<_>>();
    let mut parts = Vec::new();
    for start in &plan.roots {
        let atom_style = component_styles
            .iter()
            .find_map(|(component, style)| component.contains(start).then_some(*style))
            .unwrap_or(CanonicalAtomStyle::Aromatic);
        parts.push(write_smiles_component(
            mol,
            *start,
            None,
            &plan,
            Some(&stereo),
            atom_style,
        )?);
    }
    Ok(parts.join("."))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum CanonicalAtomStyle {
    Aromatic,
    StoredKekule,
}

#[derive(Debug, Clone)]
pub(super) struct SmilesWritePlan {
    pub(super) roots: Vec<AtomId>,
    pub(super) tree_bonds: BTreeSet<BondId>,
    pub(super) closures: BTreeMap<AtomId, Vec<SmilesRingClosure>>,
    pub(super) subtree_sizes: BTreeMap<AtomId, usize>,
}

#[derive(Debug, Clone, Copy)]
pub(super) struct SmilesRingClosure {
    pub(super) bond: BondId,
    pub(super) number: u64,
    pub(super) order: SmilesBondOrder,
    pub(super) other: AtomId,
}

/// A format-local bond representation used while emitting SMILES.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SmilesBondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

fn plan_smiles_write(
    mol: &Molecule,
    stereo: StereoWriteMode,
) -> std::result::Result<SmilesWritePlan, MolWriteError> {
    validate_smiles_writeable(mol, stereo)?;
    let mut roots = Vec::new();
    let mut visited = BTreeSet::<AtomId>::new();
    let mut tree_bonds = BTreeSet::<BondId>::new();
    let mut ring_bonds = BTreeMap::<BondId, (AtomId, AtomId, SmilesBondOrder)>::new();

    for start in mol.atom_ids() {
        if visited.contains(&start) {
            continue;
        }
        roots.push(start);
        collect_smiles_tree(
            mol,
            start,
            None,
            &mut visited,
            &mut tree_bonds,
            &mut ring_bonds,
        )?;
    }

    let mut ring_bonds = ring_bonds
        .into_iter()
        .map(|(bond_id, (a, b, order))| {
            let (first, second) = ordered_atom_pair(a, b);
            (bond_id, first, second, order)
        })
        .collect::<Vec<_>>();
    ring_bonds.sort_by_key(|(bond_id, first, second, _)| (*first, *second, *bond_id));
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

    let mut subtree_sizes = BTreeMap::new();
    for root in &roots {
        compute_smiles_subtree_sizes(mol, *root, None, &tree_bonds, &mut subtree_sizes)?;
    }

    Ok(SmilesWritePlan {
        roots,
        tree_bonds,
        closures,
        subtree_sizes,
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum StereoWriteMode {
    Reject,
    Ignore,
    Encode,
}

pub(super) fn validate_smiles_writeable(
    mol: &Molecule,
    stereo: StereoWriteMode,
) -> std::result::Result<(), MolWriteError> {
    match stereo {
        StereoWriteMode::Reject if mol.stereo_elements().next().is_some() => {
            return Err(MolWriteError::new(
                "SMILES writer cannot encode stereochemistry",
            ));
        }
        StereoWriteMode::Encode => validate_isomeric_smiles_stereo(mol)?,
        StereoWriteMode::Reject | StereoWriteMode::Ignore => {}
    }
    for (_, atom) in mol.atoms() {
        if atom.radical.is_some() {
            return Err(MolWriteError::new(
                "SMILES writer cannot encode radicals without an explicit radical token",
            ));
        }
    }
    for (_, bond) in mol.bonds() {
        match bond.order {
            BondOrder::Single | BondOrder::Double | BondOrder::Triple => {}
            BondOrder::Zero | BondOrder::Dative | BondOrder::Quadruple => {
                return Err(MolWriteError::new(
                    "SMILES writer cannot encode zero, dative, or quadruple bonds",
                ));
            }
        }
    }
    Ok(())
}

fn validate_isomeric_smiles_stereo(mol: &Molecule) -> std::result::Result<(), MolWriteError> {
    if mol.stereo_groups().next().is_some() {
        return Err(MolWriteError::new(
            "isomeric SMILES writer cannot encode enhanced stereo groups",
        ));
    }
    for (_, element) in mol.stereo_elements() {
        if !element.is_specified() {
            return Err(MolWriteError::new(
                "isomeric SMILES writer cannot encode explicitly unknown stereo",
            ));
        }
        match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => {
                if stereo.carriers.len() != 4 {
                    return Err(MolWriteError::new(
                        "isomeric SMILES writer cannot encode invalid tetrahedral stereo",
                    ));
                }
                let hydrogen_count = stereo
                    .carriers
                    .iter()
                    .filter(|carrier| matches!(carrier, StereoCarrier::ImplicitHydrogen))
                    .count();
                if hydrogen_count > 1 {
                    return Err(MolWriteError::new(
                        "isomeric SMILES writer cannot encode tetrahedral stereo with repeated implicit hydrogens",
                    ));
                }
            }
            StereoElementKind::DoubleBond(stereo) => {
                validate_isomeric_double_bond_endpoint(
                    mol,
                    stereo.left,
                    stereo.right,
                    stereo.bond,
                    stereo.left_carrier,
                )?;
                validate_isomeric_double_bond_endpoint(
                    mol,
                    stereo.right,
                    stereo.left,
                    stereo.bond,
                    stereo.right_carrier,
                )?;
            }
            StereoElementKind::Axis(_) => {
                return Err(MolWriteError::new(
                    "isomeric SMILES writer cannot encode axial stereochemistry yet",
                ));
            }
        }
    }
    Ok(())
}

fn validate_isomeric_double_bond_endpoint(
    mol: &Molecule,
    endpoint: AtomId,
    other_endpoint: AtomId,
    focus_bond: BondId,
    carrier: StereoCarrier,
) -> std::result::Result<(), MolWriteError> {
    match carrier {
        StereoCarrier::Atom(atom) => {
            let bond = mol
                .bond_between(endpoint, atom)
                .map_err(|error| MolWriteError::new(error.to_string()))?
                .ok_or_else(|| MolWriteError::new("double-bond stereo carrier is not bonded"))?;
            let order = mol
                .bond(bond)
                .map_err(|error| MolWriteError::new(error.to_string()))?
                .order;
            if order != BondOrder::Single || atom == other_endpoint {
                return Err(MolWriteError::new(
                    "isomeric SMILES writer cannot encode invalid double-bond stereo carrier",
                ));
            }
        }
        StereoCarrier::ImplicitHydrogen => {
            let atom = mol
                .atom(endpoint)
                .map_err(|error| MolWriteError::new(error.to_string()))?;
            let hydrogens = atom
                .explicit_hydrogens
                .saturating_add(mol.implicit_hydrogens(endpoint).ok().flatten().unwrap_or(0));
            if hydrogens == 0 {
                return Err(MolWriteError::new(
                    "isomeric SMILES writer cannot encode unavailable implicit double-bond hydrogen carrier",
                ));
            }
            if implicit_double_bond_printable_carrier_bond(
                mol,
                endpoint,
                other_endpoint,
                focus_bond,
            )?
            .is_none()
            {
                return Err(MolWriteError::new(
                    "isomeric SMILES writer cannot encode implicit double-bond carrier without a unique explicit substituent bond",
                ));
            }
        }
        StereoCarrier::ImplicitLonePair => {
            return Err(MolWriteError::new(
                "isomeric SMILES writer cannot encode lone-pair double-bond carrier",
            ));
        }
    }
    Ok(())
}

pub(super) fn smiles_connected_components(
    mol: &Molecule,
) -> std::result::Result<Vec<Vec<AtomId>>, MolWriteError> {
    let mut components = Vec::new();
    let mut visited = BTreeSet::new();
    for start in mol.atom_ids() {
        if !visited.insert(start) {
            continue;
        }
        let mut component = Vec::new();
        let mut stack = vec![start];
        while let Some(atom) = stack.pop() {
            component.push(atom);
            for (_, _, neighbor) in smiles_incident_bonds(mol, atom)? {
                if visited.insert(neighbor) {
                    stack.push(neighbor);
                }
            }
        }
        component.sort();
        components.push(component);
    }
    Ok(components)
}

fn collect_smiles_tree(
    mol: &Molecule,
    atom_id: AtomId,
    parent_bond: Option<BondId>,
    visited: &mut BTreeSet<AtomId>,
    tree_bonds: &mut BTreeSet<BondId>,
    ring_bonds: &mut BTreeMap<BondId, (AtomId, AtomId, SmilesBondOrder)>,
) -> std::result::Result<(), MolWriteError> {
    struct Frame {
        parent_bond: Option<BondId>,
        incident: Vec<(BondId, SmilesBondOrder, AtomId)>,
        next_edge: usize,
    }

    visited.insert(atom_id);
    let mut stack = vec![Frame {
        parent_bond,
        incident: smiles_incident_bonds(mol, atom_id)?,
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
            incident: smiles_incident_bonds(mol, neighbor)?,
            next_edge: 0,
        });
    }
    Ok(())
}

fn compute_smiles_subtree_sizes(
    mol: &Molecule,
    atom_id: AtomId,
    parent: Option<AtomId>,
    tree_bonds: &BTreeSet<BondId>,
    subtree_sizes: &mut BTreeMap<AtomId, usize>,
) -> std::result::Result<usize, MolWriteError> {
    let mut stack = vec![(atom_id, parent, false)];
    while let Some((current, parent, expanded)) = stack.pop() {
        if expanded {
            let mut size = 1usize;
            for (bond_id, _, neighbor) in smiles_incident_bonds(mol, current)? {
                if tree_bonds.contains(&bond_id) && Some(neighbor) != parent {
                    size = size
                        .saturating_add(subtree_sizes.get(&neighbor).copied().unwrap_or_default());
                }
            }
            subtree_sizes.insert(current, size);
            continue;
        }
        stack.push((current, parent, true));
        let mut children = smiles_incident_bonds(mol, current)?
            .into_iter()
            .filter(|(bond_id, _, neighbor)| {
                tree_bonds.contains(bond_id) && Some(*neighbor) != parent
            })
            .map(|(_, _, neighbor)| neighbor)
            .collect::<Vec<_>>();
        children.sort();
        for child in children.into_iter().rev() {
            stack.push((child, Some(current), false));
        }
    }
    Ok(subtree_sizes.get(&atom_id).copied().unwrap_or_default())
}

#[derive(Debug, Clone)]
struct SmilesStereoWriteContext {
    tetrahedral: BTreeMap<AtomId, TetrahedralSmilesState>,
    directional: BTreeMap<BondId, Vec<DirectionalSmilesConstraint>>,
}

#[derive(Debug, Clone)]
struct TetrahedralSmilesState {
    carriers: Vec<StereoCarrier>,
    orientation: TetrahedralOrientation,
}

#[derive(Debug, Clone, Copy)]
struct ChiralAtomWriteState {
    orientation: TetrahedralOrientation,
    force_hydrogen: bool,
}

#[derive(Debug, Clone, Copy)]
struct DirectionalSmilesConstraint {
    endpoint: AtomId,
    direction_at_endpoint: SmilesDirectionToken,
}

impl SmilesStereoWriteContext {
    fn new(mol: &Molecule) -> std::result::Result<Self, MolWriteError> {
        let mut tetrahedral = BTreeMap::new();
        let mut directional = BTreeMap::<BondId, Vec<DirectionalSmilesConstraint>>::new();
        for (_, element) in mol.stereo_elements() {
            match &element.kind {
                StereoElementKind::Tetrahedral(stereo) => {
                    let Some(orientation) = stereo.orientation else {
                        return Err(MolWriteError::new(
                            "isomeric SMILES writer cannot encode explicitly unknown tetrahedral stereo",
                        ));
                    };
                    if tetrahedral
                        .insert(
                            stereo.center,
                            TetrahedralSmilesState {
                                carriers: stereo.carriers.clone(),
                                orientation,
                            },
                        )
                        .is_some()
                    {
                        return Err(MolWriteError::new(
                            "isomeric SMILES writer cannot encode multiple tetrahedral elements on one atom",
                        ));
                    }
                }
                StereoElementKind::DoubleBond(stereo) => {
                    add_double_bond_directional_constraints(mol, stereo, &mut directional)?;
                }
                StereoElementKind::Axis(_) => {}
            }
        }
        Ok(Self {
            tetrahedral,
            directional,
        })
    }

    fn atom_chirality(
        &self,
        atom: AtomId,
        parent: Option<AtomId>,
        closures: Option<&[SmilesRingClosure]>,
        children: &[(BondId, SmilesBondOrder, AtomId)],
        main_child_index: Option<usize>,
    ) -> Option<std::result::Result<ChiralAtomWriteState, MolWriteError>> {
        let stereo = self.tetrahedral.get(&atom)?;
        Some(tetrahedral_chirality_for_smiles_order(
            stereo,
            parent,
            closures,
            children,
            main_child_index,
        ))
    }

    fn directional_bond(
        &self,
        bond: BondId,
        left: AtomId,
        right: AtomId,
    ) -> std::result::Result<Option<SmilesDirectionToken>, MolWriteError> {
        let Some(constraints) = self.directional.get(&bond) else {
            return Ok(None);
        };
        let mut concrete = None;
        for constraint in constraints {
            let mark = directional_mark_for_emitted_bond(
                constraint.direction_at_endpoint,
                constraint.endpoint,
                left,
                right,
            )?;
            if let Some(previous) = concrete {
                if previous != mark {
                    return Err(MolWriteError::new(
                        "isomeric SMILES writer cannot encode conflicting double-bond stereo constraints",
                    ));
                }
            } else {
                concrete = Some(mark);
            }
        }
        Ok(concrete)
    }
}

fn add_double_bond_directional_constraints(
    mol: &Molecule,
    stereo: &DoubleBondStereo,
    directional: &mut BTreeMap<BondId, Vec<DirectionalSmilesConstraint>>,
) -> std::result::Result<(), MolWriteError> {
    let Some(orientation) = stereo.orientation else {
        return Err(MolWriteError::new(
            "isomeric SMILES writer cannot encode explicitly unknown double-bond stereo",
        ));
    };
    let left_carrier_bond = double_bond_printable_carrier_bond(
        mol,
        stereo.left,
        stereo.right,
        stereo.bond,
        stereo.left_carrier,
    )?;
    let right_carrier_bond = double_bond_printable_carrier_bond(
        mol,
        stereo.right,
        stereo.left,
        stereo.bond,
        stereo.right_carrier,
    )?;
    let left_direction = SmilesDirectionToken::Up;
    let right_direction = match orientation {
        DoubleBondOrientation::Together => left_direction,
        DoubleBondOrientation::Opposite => invert_directional_mark(left_direction),
    };
    directional
        .entry(left_carrier_bond.bond)
        .or_default()
        .push(DirectionalSmilesConstraint {
            endpoint: stereo.left,
            direction_at_endpoint: if left_carrier_bond.invert_direction {
                invert_directional_mark(left_direction)
            } else {
                left_direction
            },
        });
    directional
        .entry(right_carrier_bond.bond)
        .or_default()
        .push(DirectionalSmilesConstraint {
            endpoint: stereo.right,
            direction_at_endpoint: if right_carrier_bond.invert_direction {
                invert_directional_mark(right_direction)
            } else {
                right_direction
            },
        });
    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct DoubleBondPrintableCarrierBond {
    bond: BondId,
    invert_direction: bool,
}

fn double_bond_printable_carrier_bond(
    mol: &Molecule,
    endpoint: AtomId,
    other_endpoint: AtomId,
    focus_bond: BondId,
    carrier: StereoCarrier,
) -> std::result::Result<DoubleBondPrintableCarrierBond, MolWriteError> {
    match carrier {
        StereoCarrier::Atom(atom) => {
            let bond = mol
                .bond_between(endpoint, atom)
                .map_err(|error| MolWriteError::new(error.to_string()))?
                .ok_or_else(|| MolWriteError::new("double-bond stereo carrier is not bonded"))?;
            Ok(DoubleBondPrintableCarrierBond {
                bond,
                invert_direction: false,
            })
        }
        StereoCarrier::ImplicitHydrogen => {
            let Some(bond) = implicit_double_bond_printable_carrier_bond(
                mol,
                endpoint,
                other_endpoint,
                focus_bond,
            )?
            else {
                return Err(MolWriteError::new(
                    "isomeric SMILES writer cannot encode implicit double-bond carrier without a unique explicit substituent bond",
                ));
            };
            Ok(DoubleBondPrintableCarrierBond {
                bond,
                invert_direction: true,
            })
        }
        StereoCarrier::ImplicitLonePair => Err(MolWriteError::new(
            "isomeric SMILES writer cannot encode lone-pair double-bond carrier",
        )),
    }
}

fn implicit_double_bond_printable_carrier_bond(
    mol: &Molecule,
    endpoint: AtomId,
    other_endpoint: AtomId,
    focus_bond: BondId,
) -> std::result::Result<Option<BondId>, MolWriteError> {
    let mut candidates = Vec::new();
    for (bond_id, bond) in mol
        .incident_bonds(endpoint)
        .map_err(|error| MolWriteError::new(error.to_string()))?
    {
        if bond_id == focus_bond || bond.order != BondOrder::Single {
            continue;
        }
        let other = bond.other_atom(endpoint);
        if other != other_endpoint {
            candidates.push(bond_id);
        }
    }
    match candidates.as_slice() {
        [bond] => Ok(Some(*bond)),
        [] => Ok(None),
        _ => Err(MolWriteError::new(
            "isomeric SMILES writer cannot encode implicit double-bond carrier with multiple explicit substituent bonds",
        )),
    }
}

fn directional_mark_for_emitted_bond(
    direction_at_endpoint: SmilesDirectionToken,
    endpoint: AtomId,
    left: AtomId,
    right: AtomId,
) -> std::result::Result<SmilesDirectionToken, MolWriteError> {
    if endpoint == left {
        Ok(direction_at_endpoint)
    } else if endpoint == right {
        Ok(invert_directional_mark(direction_at_endpoint))
    } else {
        Err(MolWriteError::new(
            "double-bond stereo endpoint is not on emitted directional bond",
        ))
    }
}

fn invert_directional_mark(kind: SmilesDirectionToken) -> SmilesDirectionToken {
    match kind {
        SmilesDirectionToken::Up => SmilesDirectionToken::Down,
        SmilesDirectionToken::Down => SmilesDirectionToken::Up,
    }
}

fn tetrahedral_chirality_for_smiles_order(
    stereo: &TetrahedralSmilesState,
    parent: Option<AtomId>,
    closures: Option<&[SmilesRingClosure]>,
    children: &[(BondId, SmilesBondOrder, AtomId)],
    main_child_index: Option<usize>,
) -> std::result::Result<ChiralAtomWriteState, MolWriteError> {
    let force_hydrogen = stereo
        .carriers
        .iter()
        .any(|carrier| matches!(carrier, StereoCarrier::ImplicitHydrogen));
    let mut emitted = Vec::with_capacity(stereo.carriers.len());
    if let Some(parent) = parent {
        emitted.push(StereoCarrier::Atom(parent));
    }
    if force_hydrogen {
        emitted.push(StereoCarrier::ImplicitHydrogen);
    }
    if let Some(closures) = closures {
        emitted.extend(
            closures
                .iter()
                .map(|closure| StereoCarrier::Atom(closure.other)),
        );
    }
    emitted.extend(
        children
            .iter()
            .enumerate()
            .filter(|(index, _)| Some(*index) != main_child_index)
            .map(|(_, (_, _, child))| StereoCarrier::Atom(*child)),
    );
    if let Some(index) = main_child_index {
        emitted.push(StereoCarrier::Atom(children[index].2));
    }
    if stereo
        .carriers
        .iter()
        .any(|carrier| matches!(carrier, StereoCarrier::ImplicitLonePair))
    {
        emitted.push(StereoCarrier::ImplicitLonePair);
    }
    if emitted != stereo.carriers {
        let Some(odd) = carrier_permutation_is_odd(&stereo.carriers, &emitted) else {
            return Err(MolWriteError::new(
                "isomeric SMILES writer cannot encode tetrahedral carrier order",
            ));
        };
        Ok(ChiralAtomWriteState {
            orientation: if odd {
                flip_tetrahedral_orientation(stereo.orientation)
            } else {
                stereo.orientation
            },
            force_hydrogen,
        })
    } else {
        Ok(ChiralAtomWriteState {
            orientation: stereo.orientation,
            force_hydrogen,
        })
    }
}

fn carrier_permutation_is_odd(from: &[StereoCarrier], to: &[StereoCarrier]) -> Option<bool> {
    if from.len() != to.len() {
        return None;
    }
    let mut positions = Vec::with_capacity(to.len());
    let mut used = vec![false; to.len()];
    for carrier in from {
        let position = to
            .iter()
            .enumerate()
            .find(|(index, candidate)| !used[*index] && *candidate == carrier)
            .map(|(index, _)| index)?;
        used[position] = true;
        positions.push(position);
    }
    let mut odd = false;
    for left in 0..positions.len() {
        for right in (left + 1)..positions.len() {
            if positions[left] > positions[right] {
                odd = !odd;
            }
        }
    }
    Some(odd)
}

fn flip_tetrahedral_orientation(orientation: TetrahedralOrientation) -> TetrahedralOrientation {
    match orientation {
        TetrahedralOrientation::Clockwise => TetrahedralOrientation::CounterClockwise,
        TetrahedralOrientation::CounterClockwise => TetrahedralOrientation::Clockwise,
    }
}

fn write_smiles_component(
    mol: &Molecule,
    atom_id: AtomId,
    parent: Option<AtomId>,
    plan: &SmilesWritePlan,
    stereo: Option<&SmilesStereoWriteContext>,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<String, MolWriteError> {
    enum Action {
        Node {
            atom: AtomId,
            parent: Option<AtomId>,
        },
        Bond {
            bond: BondId,
            order: SmilesBondOrder,
            left: AtomId,
            right: AtomId,
        },
        OpenBranch,
        CloseBranch,
    }

    let mut out = String::new();
    let mut actions = vec![Action::Node {
        atom: atom_id,
        parent,
    }];
    while let Some(action) = actions.pop() {
        match action {
            Action::OpenBranch => out.push('('),
            Action::CloseBranch => out.push(')'),
            Action::Bond {
                bond,
                order,
                left,
                right,
            } => {
                let directional = stereo
                    .map(|context| context.directional_bond(bond, left, right))
                    .transpose()?
                    .flatten();
                out.push_str(smiles_bond_between_with_direction(
                    mol,
                    order,
                    left,
                    right,
                    directional,
                )?);
            }
            Action::Node { atom, parent } => {
                let atom_record = mol
                    .atom(atom)
                    .map_err(|error| MolWriteError::new(error.to_string()))?;
                let closures = plan.closures.get(&atom).map(Vec::as_slice);
                let mut children = smiles_incident_bonds_for_style(mol, atom, atom_style)?
                    .into_iter()
                    .filter(|(bond_id, _, neighbor)| {
                        plan.tree_bonds.contains(bond_id) && Some(*neighbor) != parent
                    })
                    .collect::<Vec<_>>();
                children.sort_by_key(|(bond_id, _, child)| (*child, *bond_id));
                let main_child_index = children
                    .iter()
                    .enumerate()
                    .max_by_key(|(_, child_entry)| {
                        let child = child_entry.2;
                        (plan.subtree_sizes.get(&child).copied().unwrap_or(0), child)
                    })
                    .map(|(index, _)| index);
                let chirality = stereo
                    .and_then(|context| {
                        context.atom_chirality(atom, parent, closures, &children, main_child_index)
                    })
                    .transpose()?;
                out.push_str(&smiles_atom_with_style_and_chirality(
                    mol,
                    atom,
                    atom_record,
                    atom_style,
                    chirality.map(|state| state.orientation),
                    chirality.is_some_and(|state| state.force_hydrogen),
                )?);
                if let Some(closures) = closures {
                    for closure in closures {
                        let closure_order = match atom_style {
                            CanonicalAtomStyle::Aromatic => closure.order,
                            CanonicalAtomStyle::StoredKekule => smiles_bond_order(
                                mol.bond(closure.bond)
                                    .map_err(|error| MolWriteError::new(error.to_string()))?
                                    .order,
                            )?,
                        };
                        let directional = stereo
                            .map(|context| {
                                context.directional_bond(closure.bond, atom, closure.other)
                            })
                            .transpose()?
                            .flatten();
                        out.push_str(smiles_bond_between_with_direction(
                            mol,
                            closure_order,
                            atom,
                            closure.other,
                            directional,
                        )?);
                        out.push_str(&smiles_ring_number(closure.number));
                    }
                }

                if let Some(index) = main_child_index {
                    let (bond, order, child) = children[index];
                    actions.push(Action::Node {
                        atom: child,
                        parent: Some(atom),
                    });
                    actions.push(Action::Bond {
                        bond,
                        order,
                        left: atom,
                        right: child,
                    });
                }
                for (index, (bond, order, child)) in children.into_iter().enumerate().rev() {
                    if Some(index) == main_child_index {
                        continue;
                    }
                    actions.push(Action::CloseBranch);
                    actions.push(Action::Node {
                        atom: child,
                        parent: Some(atom),
                    });
                    actions.push(Action::Bond {
                        bond,
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

fn smiles_incident_bonds(
    mol: &Molecule,
    atom_id: AtomId,
) -> std::result::Result<Vec<(BondId, SmilesBondOrder, AtomId)>, MolWriteError> {
    smiles_incident_bonds_for_style(mol, atom_id, CanonicalAtomStyle::Aromatic)
}

pub(super) fn smiles_incident_bonds_for_style(
    mol: &Molecule,
    atom_id: AtomId,
    atom_style: CanonicalAtomStyle,
) -> std::result::Result<Vec<(BondId, SmilesBondOrder, AtomId)>, MolWriteError> {
    let mut incident = Vec::new();
    for (bond_id, bond) in mol
        .incident_bonds(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?
    {
        let order = match atom_style {
            CanonicalAtomStyle::Aromatic
                if mol.bond_is_aromatic(bond_id).ok().flatten() == Some(true)
                    && !matches!(bond.order, BondOrder::Triple | BondOrder::Quadruple) =>
            {
                SmilesBondOrder::Aromatic
            }
            CanonicalAtomStyle::Aromatic | CanonicalAtomStyle::StoredKekule => {
                smiles_bond_order(bond.order)?
            }
        };
        incident.push((bond_id, order, bond.other_atom(atom_id)));
    }
    incident.sort_by_key(|(bond_id, _, atom)| (*atom, *bond_id));
    Ok(incident)
}

pub(super) fn smiles_ring_number(number: u64) -> String {
    if number < 10 {
        number.to_string()
    } else {
        format!("%{number}")
    }
}

fn smiles_bond_order(order: BondOrder) -> std::result::Result<SmilesBondOrder, MolWriteError> {
    match order {
        BondOrder::Single => Ok(SmilesBondOrder::Single),
        BondOrder::Double => Ok(SmilesBondOrder::Double),
        BondOrder::Triple => Ok(SmilesBondOrder::Triple),
        BondOrder::Zero | BondOrder::Dative | BondOrder::Quadruple => Err(MolWriteError::new(
            "SMILES writer cannot encode zero, dative, or quadruple bonds",
        )),
    }
}

fn smiles_bond(order: SmilesBondOrder) -> &'static str {
    match order {
        SmilesBondOrder::Single => "",
        SmilesBondOrder::Double => "=",
        SmilesBondOrder::Triple => "#",
        SmilesBondOrder::Aromatic => ":",
    }
}

pub(super) fn smiles_bond_between(
    mol: &Molecule,
    order: SmilesBondOrder,
    left: AtomId,
    right: AtomId,
) -> std::result::Result<&'static str, MolWriteError> {
    if matches!(order, SmilesBondOrder::Single | SmilesBondOrder::Aromatic) {
        mol.atom(left)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        mol.atom(right)
            .map_err(|error| MolWriteError::new(error.to_string()))?;
        if mol.atom_is_aromatic(left).ok().flatten() == Some(true)
            && mol.atom_is_aromatic(right).ok().flatten() == Some(true)
        {
            return Ok(if order == SmilesBondOrder::Single {
                "-"
            } else {
                ""
            });
        }
    }
    Ok(smiles_bond(order))
}

fn smiles_bond_between_with_direction(
    mol: &Molecule,
    order: SmilesBondOrder,
    left: AtomId,
    right: AtomId,
    directional: Option<SmilesDirectionToken>,
) -> std::result::Result<&'static str, MolWriteError> {
    if let Some(directional) = directional {
        if order != SmilesBondOrder::Single {
            return Err(MolWriteError::new(
                "isomeric SMILES writer cannot place directional stereo on a non-single bond",
            ));
        }
        return match directional {
            SmilesDirectionToken::Up => Ok("/"),
            SmilesDirectionToken::Down => Ok("\\"),
        };
    }
    smiles_bond_between(mol, order, left, right)
}

pub(super) fn smiles_atom(atom: &Atom, aromatic: bool, implicit_hydrogens: u8) -> String {
    smiles_atom_with_chirality(atom, aromatic, implicit_hydrogens, None, false)
}

fn smiles_atom_with_chirality(
    atom: &Atom,
    aromatic: bool,
    implicit_hydrogens: u8,
    chirality: Option<TetrahedralOrientation>,
    force_hydrogen: bool,
) -> String {
    let explicit_hydrogens = if force_hydrogen {
        smiles_atom_explicit_hydrogens(atom, aromatic, implicit_hydrogens).max(1)
    } else {
        smiles_atom_explicit_hydrogens(atom, aromatic, implicit_hydrogens)
    };
    let organic = atom.isotope.is_none()
        && atom.formal_charge == 0
        && explicit_hydrogens == 0
        && !atom.no_implicit_hydrogens
        && atom.atom_map.is_none()
        && chirality.is_none()
        && matches!(
            atom.element.symbol(),
            "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
        );
    if organic {
        if aromatic {
            atom.element.symbol().to_ascii_lowercase()
        } else {
            atom.element.symbol().to_owned()
        }
    } else {
        let mut out = String::from("[");
        if let Some(isotope) = atom.isotope {
            out.push_str(&isotope.to_string());
        }
        if aromatic {
            out.push_str(&atom.element.symbol().to_ascii_lowercase());
        } else {
            out.push_str(atom.element.symbol());
        }
        if let Some(chirality) = chirality {
            out.push('@');
            if chirality == TetrahedralOrientation::CounterClockwise {
                out.push('@');
            }
        }
        if explicit_hydrogens > 0 {
            out.push('H');
            if explicit_hydrogens > 1 {
                out.push_str(&explicit_hydrogens.to_string());
            }
        }
        if atom.formal_charge > 0 {
            out.push('+');
            if atom.formal_charge > 1 {
                out.push_str(&atom.formal_charge.to_string());
            }
        } else if atom.formal_charge < 0 {
            out.push('-');
            if atom.formal_charge < -1 {
                out.push_str(&(-atom.formal_charge).to_string());
            }
        }
        if let Some(map) = atom.atom_map {
            out.push(':');
            out.push_str(&map.to_string());
        }
        out.push(']');
        out
    }
}

fn smiles_atom_with_style_and_chirality(
    mol: &Molecule,
    atom_id: AtomId,
    atom: &Atom,
    atom_style: CanonicalAtomStyle,
    chirality: Option<TetrahedralOrientation>,
    force_hydrogen: bool,
) -> std::result::Result<String, MolWriteError> {
    let aromatic = mol.atom_is_aromatic(atom_id).ok().flatten() == Some(true);
    let implicit_hydrogens = mol
        .implicit_hydrogens(atom_id)
        .map_err(|error| MolWriteError::new(error.to_string()))?
        .unwrap_or(0);
    if matches!(atom_style, CanonicalAtomStyle::StoredKekule) && aromatic {
        let mut normalized = atom.clone();
        normalized.isotope = None;
        let mut normalized_implicit = implicit_hydrogens;
        if !matches!(atom.element.symbol(), "B" | "C") && implicit_hydrogens > 0 {
            normalized.explicit_hydrogens =
                atom.explicit_hydrogens.saturating_add(implicit_hydrogens);
            normalized_implicit = 0;
            normalized.no_implicit_hydrogens = true;
        }
        return Ok(smiles_atom_with_chirality(
            &normalized,
            false,
            normalized_implicit,
            chirality,
            force_hydrogen,
        ));
    }
    Ok(smiles_atom_with_chirality(
        atom,
        aromatic,
        implicit_hydrogens,
        chirality,
        force_hydrogen,
    ))
}

fn smiles_atom_explicit_hydrogens(atom: &Atom, aromatic: bool, implicit_hydrogens: u8) -> u8 {
    if atom.element.symbol() == "N"
        && aromatic
        && atom.explicit_hydrogens == 0
        && implicit_hydrogens == 1
    {
        1
    } else {
        atom.explicit_hydrogens
    }
}
