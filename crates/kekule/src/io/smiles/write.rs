use std::collections::{BTreeMap, BTreeSet};

use crate::algorithms::ordered_atom_pair;
use crate::core::Molecule;
use crate::core::*;
use crate::io::MolWriteError;

use super::parse::SmilesDirectionToken;

pub fn write_smiles(molecule: &Molecule) -> std::result::Result<String, MolWriteError> {
    write_source_order_smiles(
        molecule,
        StereoWriteMode::Reject,
        CanonicalAtomStyle::Aromatic,
    )
}

pub fn write_isomeric_smiles(molecule: &Molecule) -> std::result::Result<String, MolWriteError> {
    write_source_order_smiles(
        molecule,
        StereoWriteMode::Encode,
        CanonicalAtomStyle::StoredKekule,
    )
}

fn write_source_order_smiles(
    mol: &Molecule,
    mode: StereoWriteMode,
    style: CanonicalAtomStyle,
) -> std::result::Result<String, MolWriteError> {
    let plan = plan_smiles_write(mol, mode)?;
    let stereo = (mode == StereoWriteMode::Encode)
        .then(|| SmilesStereoWriteContext::new(mol, AtomId::index))
        .transpose()?;
    let mut parts = Vec::new();
    for start in &plan.roots {
        parts.push(write_smiles_component(
            mol,
            *start,
            &plan,
            stereo.as_ref(),
            style,
            |_, children| {
                children.sort_by_key(|(bond, _, atom)| (*atom, *bond));
                children
                    .iter()
                    .enumerate()
                    .max_by_key(|(_, (_, _, child))| {
                        (plan.subtree_sizes.get(child).copied().unwrap_or(0), *child)
                    })
                    .map(|(index, _)| index)
            },
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
    pub(super) order: SmilesBondOrder,
    pub(super) other: AtomId,
}

/// A format-local bond representation used while emitting SMILES.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum SmilesBondOrder {
    Single,
    Double,
    Triple,
    Quadruple,
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
            |atom| smiles_incident_bonds(mol, atom),
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
    let closures = smiles_ring_closures(ring_bonds);

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
        StereoWriteMode::Reject => {}
    }
    for (_, atom) in mol.atoms() {
        if atom.radical.is_some() {
            return Err(MolWriteError::new(
                "SMILES writer cannot encode radicals without an explicit radical token",
            ));
        }
        if matches!(atom.hydrogens, HydrogenDeclaration::Infer { explicit } if explicit > 0) {
            return Err(MolWriteError::new(
                "SMILES cannot encode represented hydrogens while leaving implicit-H inference enabled",
            ));
        }
    }
    for (_, bond) in mol.bonds() {
        match bond.order {
            BondOrder::Single | BondOrder::Double | BondOrder::Triple | BondOrder::Quadruple => {}
            BondOrder::Zero | BondOrder::Dative => {
                return Err(MolWriteError::new(
                    "SMILES writer cannot encode zero or dative bonds",
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
                let atom = mol
                    .atom(stereo.center)
                    .map_err(|error| MolWriteError::new(error.to_string()))?;
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
                let implicit = mol
                    .implicit_hydrogens(stereo.center)
                    .map_err(|error| MolWriteError::new(error.to_string()))?;
                if atom.hydrogens.allows_implicit() && implicit.is_none() {
                    return Err(MolWriteError::new(
                        "SMILES stereo atom requires installed hydrogen perception",
                    ));
                }
                if usize::from(atom.hydrogens.explicit_count()) + usize::from(implicit.unwrap_or(0))
                    != hydrogen_count
                {
                    return Err(MolWriteError::new(
                        "SMILES tetrahedral hydrogen carriers disagree with the atom hydrogen count",
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
                .hydrogens
                .explicit_count()
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

pub(super) fn collect_smiles_tree<F>(
    mol: &Molecule,
    atom_id: AtomId,
    parent_bond: Option<BondId>,
    visited: &mut BTreeSet<AtomId>,
    tree_bonds: &mut BTreeSet<BondId>,
    ring_bonds: &mut BTreeMap<BondId, (AtomId, AtomId, SmilesBondOrder)>,
    mut incident_bonds: F,
) -> std::result::Result<(), MolWriteError>
where
    F: FnMut(AtomId) -> std::result::Result<Vec<(BondId, SmilesBondOrder, AtomId)>, MolWriteError>,
{
    struct Frame {
        parent_bond: Option<BondId>,
        incident: Vec<(BondId, SmilesBondOrder, AtomId)>,
        next_edge: usize,
    }

    visited.insert(atom_id);
    let mut stack = vec![Frame {
        parent_bond,
        incident: incident_bonds(atom_id)?,
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
            incident: incident_bonds(neighbor)?,
            next_edge: 0,
        });
    }
    Ok(())
}

pub(super) fn smiles_ring_closures(
    ring_bonds: Vec<(BondId, AtomId, AtomId, SmilesBondOrder)>,
) -> BTreeMap<AtomId, Vec<SmilesRingClosure>> {
    let mut closures = BTreeMap::<AtomId, Vec<SmilesRingClosure>>::new();
    for (bond, first, second, order) in ring_bonds {
        closures.entry(first).or_default().push(SmilesRingClosure {
            bond,
            order,
            other: second,
        });
        closures.entry(second).or_default().push(SmilesRingClosure {
            bond,
            order,
            other: first,
        });
    }
    closures
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
pub(super) struct SmilesStereoWriteContext {
    tetrahedral: BTreeMap<AtomId, TetrahedralSmilesState>,
    directional: BTreeMap<BondId, DirectionalSmilesConstraint>,
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
    component: BondId,
}

struct DirectionalBondConstraints {
    preferred: DirectionalSmilesConstraint,
    // The Boolean records whether the two canonically directed bonds must
    // carry opposite marks. Flipping both marks preserves the stereo assertion.
    neighbors: Vec<(BondId, bool)>,
}

impl SmilesStereoWriteContext {
    pub(super) fn new(
        mol: &Molecule,
        rank: impl Fn(AtomId) -> usize,
    ) -> std::result::Result<Self, MolWriteError> {
        let mut tetrahedral = BTreeMap::new();
        let mut double_bonds = Vec::new();
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
                    double_bonds.push(stereo.clone());
                }
                StereoElementKind::Axis(_) => {}
            }
        }
        let directional = choose_directional_bonds(mol, double_bonds, &rank)?;
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
        phases: &mut BTreeMap<BondId, bool>,
    ) -> std::result::Result<Option<SmilesDirectionToken>, MolWriteError> {
        let Some(constraint) = self.directional.get(&bond) else {
            return Ok(None);
        };
        let direction = directional_mark_for_emitted_bond(
            constraint.direction_at_endpoint,
            constraint.endpoint,
            left,
            right,
        )?;
        let invert = *phases
            .entry(constraint.component)
            .or_insert(direction == SmilesDirectionToken::Down);
        Ok(Some(if invert {
            invert_directional_mark(direction)
        } else {
            direction
        }))
    }
}

fn validate_directional_projection(
    mol: &Molecule,
    directional: &BTreeMap<BondId, DirectionalSmilesConstraint>,
) -> std::result::Result<(), MolWriteError> {
    use crate::chemistry::{
        normalize_source_stereo, SourceStereoBondMark, SourceStereoBondMarkKind,
    };
    let expected = mol
        .stereo_elements()
        .filter_map(|(_, element)| match &element.kind {
            StereoElementKind::DoubleBond(value) => Some((value.bond, value)),
            _ => None,
        })
        .collect::<BTreeMap<_, _>>();
    if expected.is_empty() {
        return Ok(());
    }
    let marks = directional
        .iter()
        .map(|(&bond, value)| SourceStereoBondMark {
            bond,
            from: value.endpoint,
            kind: match value.direction_at_endpoint {
                SmilesDirectionToken::Up => SourceStereoBondMarkKind::DirectionalUp,
                SmilesDirectionToken::Down => SourceStereoBondMarkKind::DirectionalDown,
            },
        })
        .collect::<Vec<_>>();
    let mut decoded = mol.clone();
    decoded.graph.stereo_elements.clear();
    decoded.graph.stereo_groups.clear();
    normalize_source_stereo(&mut decoded, None, &marks).map_err(|error| {
        MolWriteError::new(format!(
            "SMILES directional stereo cannot be represented: {error}"
        ))
    })?;
    if decoded.stereo_elements().count() != expected.len() {
        return Err(MolWriteError::new(
            "SMILES directional bonds would add or lose a double-bond stereo assertion",
        ));
    }
    for (_, element) in decoded.stereo_elements() {
        let StereoElementKind::DoubleBond(actual) = &element.kind else {
            unreachable!()
        };
        let Some(source) = expected.get(&actual.bond) else {
            return Err(MolWriteError::new(
                "SMILES directional bonds would specify an unasserted double bond",
            ));
        };
        let inverted = (source.left_carrier != actual.left_carrier)
            != (source.right_carrier != actual.right_carrier);
        let orientation = source
            .orientation
            .map(|value| if inverted { value.inverted() } else { value });
        if actual.orientation != orientation {
            return Err(MolWriteError::new(
                "SMILES directional bonds would change double-bond configuration",
            ));
        }
    }
    Ok(())
}

fn add_double_bond_directional_constraints(
    mol: &Molecule,
    stereo: &DoubleBondStereo,
    [left_carrier_bond, right_carrier_bond]: [DoubleBondPrintableCarrierBond; 2],
    constraints: &mut BTreeMap<BondId, DirectionalBondConstraints>,
) -> std::result::Result<(), MolWriteError> {
    let orientation = stereo.orientation.expect("validated specified double bond");
    let left_direction = SmilesDirectionToken::Up;
    let right_direction = match orientation {
        DoubleBondOrientation::Together => left_direction,
        DoubleBondOrientation::Opposite => invert_directional_mark(left_direction),
    };
    let left =
        canonical_directional_constraint(mol, left_carrier_bond, stereo.left, left_direction)?;
    let right =
        canonical_directional_constraint(mol, right_carrier_bond, stereo.right, right_direction)?;
    let opposite = left.direction_at_endpoint != right.direction_at_endpoint;
    for (bond, preferred, neighbor) in [
        (left_carrier_bond.bond, left, right_carrier_bond.bond),
        (right_carrier_bond.bond, right, left_carrier_bond.bond),
    ] {
        constraints
            .entry(bond)
            .or_insert_with(|| DirectionalBondConstraints {
                preferred,
                neighbors: Vec::new(),
            })
            .neighbors
            .push((neighbor, opposite));
    }
    Ok(())
}

fn canonical_directional_constraint(
    mol: &Molecule,
    carrier: DoubleBondPrintableCarrierBond,
    endpoint: AtomId,
    direction: SmilesDirectionToken,
) -> std::result::Result<DirectionalSmilesConstraint, MolWriteError> {
    let bond = mol
        .bond(carrier.bond)
        .map_err(|error| MolWriteError::new(error.to_string()))?;
    let (left, right) = ordered_atom_pair(bond.a(), bond.b());
    let direction = if carrier.invert_direction {
        invert_directional_mark(direction)
    } else {
        direction
    };
    Ok(DirectionalSmilesConstraint {
        endpoint: left,
        direction_at_endpoint: directional_mark_for_emitted_bond(direction, endpoint, left, right)?,
        component: carrier.bond,
    })
}

fn solve_directional_constraints(
    constraints: BTreeMap<BondId, DirectionalBondConstraints>,
) -> std::result::Result<BTreeMap<BondId, DirectionalSmilesConstraint>, MolWriteError> {
    let mut assigned = BTreeMap::new();
    for (&seed, initial) in &constraints {
        if assigned.contains_key(&seed) {
            continue;
        }
        assigned.insert(seed, initial.preferred);
        let mut pending = vec![seed];
        while let Some(bond) = pending.pop() {
            let direction = assigned[&bond].direction_at_endpoint;
            for &(neighbor, opposite) in &constraints[&bond].neighbors {
                let expected = if opposite {
                    invert_directional_mark(direction)
                } else {
                    direction
                };
                if let Some(previous) = assigned.get(&neighbor) {
                    if previous.direction_at_endpoint != expected {
                        return Err(MolWriteError::new(
                            "isomeric SMILES writer cannot encode conflicting double-bond stereo constraints",
                        ));
                    }
                } else {
                    assigned.insert(
                        neighbor,
                        DirectionalSmilesConstraint {
                            endpoint: constraints[&neighbor].preferred.endpoint,
                            direction_at_endpoint: expected,
                            component: seed,
                        },
                    );
                    pending.push(neighbor);
                }
            }
        }
    }
    Ok(assigned)
}

#[derive(Debug, Clone, Copy)]
struct DoubleBondPrintableCarrierBond {
    bond: BondId,
    invert_direction: bool,
}

fn choose_directional_bonds(
    mol: &Molecule,
    mut stereo: Vec<DoubleBondStereo>,
    rank: &impl Fn(AtomId) -> usize,
) -> std::result::Result<BTreeMap<BondId, DirectionalSmilesConstraint>, MolWriteError> {
    if stereo.is_empty() {
        return Ok(BTreeMap::new());
    }
    // Normalize endpoint and assertion order before searching: the first valid
    // encoding must depend on canonical labels, not the source's carrier choice.
    for value in &mut stereo {
        if rank(value.left) > rank(value.right) {
            std::mem::swap(&mut value.left, &mut value.right);
            std::mem::swap(&mut value.left_carrier, &mut value.right_carrier);
        }
    }
    stereo.sort_by_key(|value| (rank(value.left), rank(value.right)));
    let mut candidates = Vec::new();
    for value in &stereo {
        for (endpoint, carrier) in [
            (value.left, value.left_carrier),
            (value.right, value.right_carrier),
        ] {
            let mut bonds = mol
                .incident_bonds(endpoint)
                .map_err(|error| MolWriteError::new(error.to_string()))?
                .filter(|(id, bond)| *id != value.bond && bond.order == BondOrder::Single)
                .map(|(id, bond)| {
                    (
                        rank(bond.other_atom(endpoint)),
                        DoubleBondPrintableCarrierBond {
                            bond: id,
                            invert_direction: carrier
                                != StereoCarrier::Atom(bond.other_atom(endpoint)),
                        },
                    )
                })
                .collect::<Vec<_>>();
            bonds.sort_by_key(|(label, _)| *label);
            if bonds.is_empty() {
                return Err(MolWriteError::new(
                    "SMILES double-bond stereo requires an explicit single-bond substituent",
                ));
            }
            candidates.push(bonds.into_iter().map(|(_, bond)| bond).collect::<Vec<_>>());
        }
    }
    let mut choices = vec![0; candidates.len()];
    for attempt in 0usize..4096 {
        if attempt
            .saturating_add(1)
            .saturating_mul(mol.atom_count().saturating_add(mol.bond_count()))
            > 50_000_000
        {
            return Err(MolWriteError::resource_limit(
                "SMILES directional assignment exceeds 50,000,000 graph visits",
            ));
        }
        let mut constraints = BTreeMap::new();
        for (index, value) in stereo.iter().enumerate() {
            add_double_bond_directional_constraints(
                mol,
                value,
                [
                    candidates[2 * index][choices[2 * index]],
                    candidates[2 * index + 1][choices[2 * index + 1]],
                ],
                &mut constraints,
            )?;
        }
        let result = solve_directional_constraints(constraints).and_then(|directional| {
            validate_directional_projection(mol, &directional)?;
            Ok(directional)
        });
        let error = match result {
            Ok(directional) => return Ok(directional),
            Err(error) => error,
        };
        let mut next = false;
        for index in (0..choices.len()).rev() {
            choices[index] += 1;
            if choices[index] < candidates[index].len() {
                next = true;
                break;
            }
            choices[index] = 0;
        }
        if !next {
            return Err(error);
        }
    }
    Err(MolWriteError::resource_limit(
        "SMILES directional assignment exceeds 4,096 carrier combinations",
    ))
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
    let lone_pair = stereo.carriers.contains(&StereoCarrier::ImplicitLonePair);
    let mut emitted = Vec::with_capacity(stereo.carriers.len());
    if let Some(parent) = parent {
        emitted.push(StereoCarrier::Atom(parent));
    }
    if force_hydrogen && !lone_pair {
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
    if lone_pair {
        if force_hydrogen {
            emitted.push(StereoCarrier::ImplicitHydrogen);
        }
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
                stereo.orientation.inverted()
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

pub(super) fn write_smiles_component(
    mol: &Molecule,
    atom_id: AtomId,
    plan: &SmilesWritePlan,
    stereo: Option<&SmilesStereoWriteContext>,
    atom_style: CanonicalAtomStyle,
    order_children: impl Fn(AtomId, &mut Vec<(BondId, SmilesBondOrder, AtomId)>) -> Option<usize>,
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
    let mut phases = BTreeMap::new();
    let mut open_rings = BTreeMap::new();
    let mut available_rings = (0..=99u64).collect::<BTreeSet<_>>();
    let mut actions = vec![Action::Node {
        atom: atom_id,
        parent: None,
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
                    .map(|context| context.directional_bond(bond, left, right, &mut phases))
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
                let main_child_index = order_children(atom, &mut children);
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
                                context.directional_bond(
                                    closure.bond,
                                    atom,
                                    closure.other,
                                    &mut phases,
                                )
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
                        let number = if let Some(number) = open_rings.remove(&closure.bond) {
                            available_rings.insert(number);
                            number
                        } else {
                            let number = available_rings
                                .range(1..)
                                .next()
                                .copied()
                                .or_else(|| available_rings.first().copied())
                                .ok_or_else(|| {
                                    MolWriteError::resource_limit(
                                        "SMILES requires more than 100 simultaneous ring labels",
                                    )
                                })?;
                            available_rings.remove(&number);
                            open_rings.insert(closure.bond, number);
                            number
                        };
                        out.push_str(&smiles_ring_number(number));
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
        BondOrder::Quadruple => Ok(SmilesBondOrder::Quadruple),
        BondOrder::Zero | BondOrder::Dative => Err(MolWriteError::new(
            "SMILES writer cannot encode zero or dative bonds",
        )),
    }
}

fn smiles_bond(order: SmilesBondOrder) -> &'static str {
    match order {
        SmilesBondOrder::Single => "",
        SmilesBondOrder::Double => "=",
        SmilesBondOrder::Triple => "#",
        SmilesBondOrder::Quadruple => "$",
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
    let organic =
        explicit_hydrogens == 0 && chirality.is_none() && smiles_atom_is_organic_subset(atom);
    if organic {
        if aromatic {
            atom.element.symbol().to_ascii_lowercase()
        } else {
            atom.element.symbol().to_owned()
        }
    } else {
        // Bracket atoms do not infer hydrogens in SMILES. Metadata (isotopes,
        // maps, charge or stereo) can require brackets even when the stored
        // atom permits inference, so materialize the installed count here.
        let explicit_hydrogens = atom
            .hydrogens
            .explicit_count()
            .saturating_add(implicit_hydrogens)
            .max(explicit_hydrogens);
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
                out.push_str(&atom.formal_charge.unsigned_abs().to_string());
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
    let written = if matches!(atom_style, CanonicalAtomStyle::StoredKekule) && aromatic {
        let mut normalized = atom.clone();
        let mut normalized_implicit = implicit_hydrogens;
        if !matches!(atom.element.symbol(), "B" | "C") && implicit_hydrogens > 0 {
            normalized.hydrogens = HydrogenDeclaration::Fixed(
                atom.hydrogens
                    .explicit_count()
                    .saturating_add(implicit_hydrogens),
            );
            normalized_implicit = 0;
        }
        smiles_atom_with_chirality(
            &normalized,
            false,
            normalized_implicit,
            chirality,
            force_hydrogen,
        )
    } else {
        smiles_atom_with_chirality(
            atom,
            aromatic,
            implicit_hydrogens,
            chirality,
            force_hydrogen,
        )
    };
    if atom.hydrogens.allows_implicit() && perceived_hydrogens.is_none() && written.starts_with('[')
    {
        return Err(MolWriteError::new(format!(
            "SMILES bracket atom {atom_id} requires installed hydrogen perception; perceive the molecule before writing"
        )));
    }
    Ok(written)
}

fn smiles_atom_explicit_hydrogens(atom: &Atom, aromatic: bool, implicit_hydrogens: u8) -> u8 {
    if atom.element.symbol() == "N"
        && aromatic
        && atom.hydrogens.explicit_count() == 0
        && implicit_hydrogens == 1
    {
        1
    } else {
        atom.hydrogens.explicit_count()
    }
}

fn smiles_atom_is_organic_subset(atom: &Atom) -> bool {
    atom.isotope.is_none()
        && atom.formal_charge == 0
        && atom.hydrogens.allows_implicit()
        && atom.atom_map.is_none()
        && matches!(
            atom.element.symbol(),
            "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
        )
}

pub(super) fn smiles_atom_requires_brackets(
    atom: &Atom,
    aromatic: bool,
    implicit_hydrogens: u8,
) -> bool {
    !smiles_atom_is_organic_subset(atom)
        || smiles_atom_explicit_hydrogens(atom, aromatic, implicit_hydrogens) > 0
}
