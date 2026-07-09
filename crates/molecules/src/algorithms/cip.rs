use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet, VecDeque};
use std::rc::Rc;

use crate::algorithms::{
    validate_stereo_with_options, RingMembership, StereoPerceptionIssue, StereoPerceptionOptions,
};
use crate::core::*;

use super::rings::{bond_in_ring_smaller_than, compute_ring_membership};

type CipResult<T> = std::result::Result<T, CipAssignmentIssue>;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CipAssignmentOptions {
    pub validate_existing: bool,
    pub max_depth: usize,
    pub max_nodes: usize,
}

impl Default for CipAssignmentOptions {
    fn default() -> Self {
        Self {
            validate_existing: true,
            max_depth: 32,
            max_nodes: 100_000,
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CipAssignmentReport {
    pub assigned: Vec<CipAssignment>,
    pub skipped: Vec<CipSkipped>,
    pub issues: Vec<CipAssignmentIssue>,
}

impl CipAssignmentReport {
    pub fn is_ok(&self) -> bool {
        self.issues.is_empty()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CipAssignment {
    pub element: StereoElementId,
    pub descriptor: StereoDescriptor,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CipSkipped {
    pub element: StereoElementId,
    pub reason: CipSkippedReason,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CipSkippedReason {
    NotSpecified,
    UnsupportedAxis,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CipAssignmentIssue {
    InvalidStereo {
        issue: StereoPerceptionIssue,
    },
    UnresolvedPriority {
        element: StereoElementId,
    },
    ResourceLimitExceeded {
        element: StereoElementId,
        max_nodes: usize,
    },
}

pub fn assign_cip_descriptors(mol: &mut Molecule) -> CipAssignmentReport {
    assign_cip_descriptors_with_options(mol, CipAssignmentOptions::default())
}

pub fn assign_cip_descriptors_with_options(
    mol: &mut Molecule,
    options: CipAssignmentOptions,
) -> CipAssignmentReport {
    clear_stereo_descriptors(mol);
    let mut report = CipAssignmentReport::default();
    if options.validate_existing {
        let validation = validate_stereo_with_options(
            mol,
            StereoPerceptionOptions {
                validate_existing: true,
                detect_candidates: false,
                assemble_source_marks: false,
                assign_coordinates: false,
            },
        );
        report.issues.extend(
            validation
                .issues
                .into_iter()
                .map(|issue| CipAssignmentIssue::InvalidStereo { issue }),
        );
        if !report.issues.is_empty() {
            return report;
        }
    }

    let mut pending = mol
        .stereo_elements()
        .map(|(id, element)| (id, element.clone()))
        .collect::<Vec<_>>();

    while !pending.is_empty() {
        let round_mol = mol.clone();
        let mut next_pending = Vec::new();
        let mut round_assignments = Vec::new();
        let mut assigned_this_round = false;
        for (id, element) in pending {
            match assign_cip_element(&round_mol, id, &element, options) {
                CipElementAssignment::Assigned(descriptor) => {
                    round_assignments.push((id, descriptor));
                    assigned_this_round = true;
                }
                CipElementAssignment::Skipped(reason) => {
                    report.skipped.push(CipSkipped {
                        element: id,
                        reason,
                    });
                }
                CipElementAssignment::Deferred => next_pending.push((id, element)),
                CipElementAssignment::Issue(issue) => report.issues.push(issue),
            }
        }
        for (id, descriptor) in round_assignments {
            set_stereo_descriptor(mol, id, descriptor);
            report.assigned.push(CipAssignment {
                element: id,
                descriptor,
            });
        }
        if !assigned_this_round {
            match assign_deferred_tetrahedral_rule6(mol, &next_pending, options) {
                Ok(assignments) if !assignments.is_empty() => {
                    let has_absolute_assignment = assignments
                        .iter()
                        .any(|(_, descriptor)| descriptor_is_absolute_tetrahedral(*descriptor));
                    let assignments_to_apply = assignments
                        .into_iter()
                        .filter(|(_, descriptor)| {
                            !has_absolute_assignment
                                || descriptor_is_absolute_tetrahedral(*descriptor)
                        })
                        .collect::<Vec<_>>();
                    let assigned_ids = assignments_to_apply
                        .iter()
                        .map(|(id, _)| *id)
                        .collect::<Vec<_>>();
                    for (id, descriptor) in assignments_to_apply {
                        set_stereo_descriptor(mol, id, descriptor);
                        report.assigned.push(CipAssignment {
                            element: id,
                            descriptor,
                        });
                    }
                    pending = next_pending
                        .into_iter()
                        .filter(|(id, _)| !assigned_ids.contains(id))
                        .collect();
                    continue;
                }
                Ok(_) => {}
                Err(issue) => {
                    report.issues.push(issue);
                    break;
                }
            }
            for (id, _) in next_pending {
                report
                    .issues
                    .push(CipAssignmentIssue::UnresolvedPriority { element: id });
            }
            break;
        }
        pending = next_pending;
    }
    report
}

fn descriptor_is_absolute_tetrahedral(descriptor: StereoDescriptor) -> bool {
    matches!(descriptor, StereoDescriptor::R | StereoDescriptor::S)
}

enum CipElementAssignment {
    Assigned(StereoDescriptor),
    Skipped(CipSkippedReason),
    Deferred,
    Issue(CipAssignmentIssue),
}

fn assign_cip_element(
    mol: &Molecule,
    id: StereoElementId,
    element: &StereoElement,
    options: CipAssignmentOptions,
) -> CipElementAssignment {
    if element.specifiedness != StereoSpecifiedness::Specified {
        return CipElementAssignment::Skipped(CipSkippedReason::NotSpecified);
    }
    let assignment = match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            assign_tetrahedral_descriptor(mol, id, stereo, options)
        }
        StereoElementKind::DoubleBond(stereo) => {
            assign_double_bond_descriptor(mol, id, stereo, options)
        }
        StereoElementKind::Axis(_) => {
            return CipElementAssignment::Skipped(CipSkippedReason::UnsupportedAxis);
        }
    };
    match assignment {
        Ok(descriptor) => CipElementAssignment::Assigned(descriptor),
        Err(CipAssignmentIssue::UnresolvedPriority { .. }) => CipElementAssignment::Deferred,
        Err(issue) => CipElementAssignment::Issue(issue),
    }
}

fn clear_stereo_descriptors(mol: &mut Molecule) {
    for element in mol.stereo_elements.iter_mut().flatten() {
        element.descriptor = None;
    }
}

fn set_stereo_descriptor(mol: &mut Molecule, id: StereoElementId, descriptor: StereoDescriptor) {
    if let Some(element) = mol
        .stereo_elements
        .get_mut(id.index())
        .and_then(Option::as_mut)
    {
        element.descriptor = Some(descriptor);
    }
}

fn assign_tetrahedral_descriptor(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &TetrahedralStereo,
    options: CipAssignmentOptions,
) -> CipResult<StereoDescriptor> {
    assign_tetrahedral_descriptor_with_deferred_rule6(mol, element, stereo, options, false)
}

fn assign_tetrahedral_descriptor_with_deferred_rule6(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &TetrahedralStereo,
    options: CipAssignmentOptions,
    allow_single_ring_tied_pair_rule6: bool,
) -> CipResult<StereoDescriptor> {
    let ranked = ranked_tetrahedral_carriers(
        mol,
        element,
        stereo.center,
        &stereo.carriers,
        stereo.orientation,
        options,
        allow_single_ring_tied_pair_rule6,
    )?;
    tetrahedral_descriptor_from_ranked(element, stereo, &ranked)
}

fn tetrahedral_descriptor_from_ranked(
    element: StereoElementId,
    stereo: &TetrahedralStereo,
    ranked: &RankedCarriers,
) -> CipResult<StereoDescriptor> {
    let mut priority_positions = Vec::new();
    for carrier in &ranked.carriers {
        let Some(position) = stereo
            .carriers
            .iter()
            .position(|candidate| candidate == carrier)
        else {
            return Err(CipAssignmentIssue::UnresolvedPriority { element });
        };
        priority_positions.push(position);
    }
    let even = permutation_is_even(&priority_positions);
    let descriptor_is_r = matches!(stereo.orientation, TetrahedralOrientation::Clockwise) != even;
    let descriptor = match (descriptor_is_r, ranked.pseudo_asymmetric_ordering) {
        (true, true) => StereoDescriptor::LowerR,
        (false, true) => StereoDescriptor::LowerS,
        (true, false) => StereoDescriptor::R,
        (false, false) => StereoDescriptor::S,
    };
    Ok(descriptor)
}

fn assign_double_bond_descriptor(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &DoubleBondStereo,
    options: CipAssignmentOptions,
) -> CipResult<StereoDescriptor> {
    if bond_in_ring_smaller_than(mol, stereo.bond, 8) {
        return Err(CipAssignmentIssue::UnresolvedPriority { element });
    }
    let left_carriers = double_bond_endpoint_carriers(mol, stereo.left, stereo.right, stereo.bond);
    let right_carriers = double_bond_endpoint_carriers(mol, stereo.right, stereo.left, stereo.bond);
    let left_top = ranked_carriers(mol, element, stereo.left, &left_carriers, options)?
        .carriers
        .first()
        .copied()
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let right_top = ranked_carriers(mol, element, stereo.right, &right_carriers, options)?
        .carriers
        .first()
        .copied()
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;

    let mut top_relation = stereo.orientation;
    if stereo.left_carrier != left_top {
        top_relation = invert_double_bond_orientation(top_relation);
    }
    if stereo.right_carrier != right_top {
        top_relation = invert_double_bond_orientation(top_relation);
    }
    Ok(match top_relation {
        DoubleBondOrientation::Together => StereoDescriptor::Z,
        DoubleBondOrientation::Opposite => StereoDescriptor::E,
    })
}

fn ranked_carriers(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    carriers: &[StereoCarrier],
    options: CipAssignmentOptions,
) -> CipResult<RankedCarriers> {
    let signatures = carrier_signatures(mol, element, root, carriers, options, false)?;
    rank_carrier_signatures(element, &signatures, None)
}

fn ranked_tetrahedral_carriers(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    carriers: &[StereoCarrier],
    orientation: TetrahedralOrientation,
    options: CipAssignmentOptions,
    allow_single_ring_tied_pair_rule6: bool,
) -> CipResult<RankedCarriers> {
    let signatures = carrier_signatures(
        mol,
        element,
        root,
        carriers,
        options,
        allow_single_ring_tied_pair_rule6,
    )?;
    match rank_carrier_signatures(element, &signatures, None) {
        Ok(ranked) => Ok(ranked),
        Err(CipAssignmentIssue::UnresolvedPriority { .. }) if carriers.len() == 4 => {
            rank_tetrahedral_signatures_with_rule6(
                mol,
                element,
                root,
                &signatures,
                orientation,
                allow_single_ring_tied_pair_rule6,
            )
        }
        Err(issue) => Err(issue),
    }
}

fn assign_deferred_tetrahedral_rule6(
    mol: &Molecule,
    pending: &[(StereoElementId, StereoElement)],
    options: CipAssignmentOptions,
) -> CipResult<Vec<(StereoElementId, StereoDescriptor)>> {
    let mut assignments = Vec::new();
    for (id, element) in pending {
        let StereoElementKind::Tetrahedral(stereo) = &element.kind else {
            continue;
        };
        if element.specifiedness != StereoSpecifiedness::Specified {
            continue;
        }
        match assign_tetrahedral_descriptor_with_deferred_rule6(mol, *id, stereo, options, true) {
            Ok(descriptor) => assignments.push((*id, descriptor)),
            Err(CipAssignmentIssue::UnresolvedPriority { .. }) => {}
            Err(issue) => return Err(issue),
        }
    }
    Ok(assignments)
}

fn carrier_signatures(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    carriers: &[StereoCarrier],
    options: CipAssignmentOptions,
    allow_auxiliary_descriptors: bool,
) -> CipResult<Vec<(StereoCarrier, LigandSignature)>> {
    let atomic_number_fractions = cip_atomic_number_fractions(mol);
    if allow_auxiliary_descriptors {
        let descriptor_context = DescriptorContext::new(element, AuxiliaryDescriptorMode::Collect);
        let aux_graph =
            build_auxiliary_graph(mol, element, root, options, &atomic_number_fractions)?;
        collect_auxiliary_occurrences_from_graph(mol, &descriptor_context, &aux_graph);
        precompute_auxiliary_descriptors(
            mol,
            &descriptor_context,
            &aux_graph,
            options,
            &atomic_number_fractions,
        );
        let descriptor_context = descriptor_context.with_mode(AuxiliaryDescriptorMode::Precomputed);
        let signatures = build_carrier_signatures(
            mol,
            element,
            &descriptor_context,
            root,
            carriers,
            options,
            &atomic_number_fractions,
        )?;
        return Ok(signatures);
    }
    let descriptor_context = DescriptorContext::new(element, AuxiliaryDescriptorMode::Disabled);
    build_carrier_signatures(
        mol,
        element,
        &descriptor_context,
        root,
        carriers,
        options,
        &atomic_number_fractions,
    )
}

fn build_carrier_signatures(
    mol: &Molecule,
    element: StereoElementId,
    descriptor_context: &DescriptorContext,
    root: AtomId,
    carriers: &[StereoCarrier],
    options: CipAssignmentOptions,
    atomic_number_fractions: &[AtomicNumberFraction],
) -> CipResult<Vec<(StereoCarrier, LigandSignature)>> {
    carriers
        .iter()
        .copied()
        .map(|carrier| {
            carrier_signature(
                mol,
                element,
                descriptor_context,
                carrier,
                root,
                options,
                atomic_number_fractions,
            )
            .map(|signature| (carrier, signature))
        })
        .collect::<CipResult<Vec<_>>>()
}

fn rank_carrier_signatures(
    element: StereoElementId,
    signatures: &[(StereoCarrier, LigandSignature)],
    rule6_reference: Option<AtomId>,
) -> CipResult<RankedCarriers> {
    let mut pseudo_asymmetric_pair_count = 0usize;
    for left in 0..signatures.len() {
        for right in (left + 1)..signatures.len() {
            let comparison = signatures[left]
                .1
                .compare_with_rule6_reference(&signatures[right].1, rule6_reference);
            if comparison.ordering == Ordering::Equal {
                return Err(CipAssignmentIssue::UnresolvedPriority { element });
            }
            if comparison.pseudo_asymmetric {
                pseudo_asymmetric_pair_count += 1;
            }
        }
    }
    let mut signatures = signatures.to_vec();
    signatures.sort_by(|left, right| {
        right
            .1
            .compare_with_rule6_reference(&left.1, rule6_reference)
            .ordering
    });
    Ok(RankedCarriers {
        carriers: signatures.into_iter().map(|(carrier, _)| carrier).collect(),
        pseudo_asymmetric_ordering: pseudo_asymmetric_pair_count == 1,
    })
}

fn rank_tetrahedral_signatures_with_rule6(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    signatures: &[(StereoCarrier, LigandSignature)],
    orientation: TetrahedralOrientation,
    allow_single_ring_tied_pair_rule6: bool,
) -> CipResult<RankedCarriers> {
    let groups = grouped_signature_indices(signatures);
    match groups.len() {
        2 => {
            let Some(reference_index) = groups.iter().flatten().copied().nth(1) else {
                return Err(CipAssignmentIssue::UnresolvedPriority { element });
            };
            let Some(reference) = carrier_rule6_atom(signatures[reference_index].0) else {
                return Err(CipAssignmentIssue::UnresolvedPriority { element });
            };
            rank_carrier_signatures(element, signatures, Some(reference))
        }
        1 => rank_s4_tetrahedral_signatures_with_rule6(element, signatures, &groups[0]),
        _ if allow_single_ring_tied_pair_rule6 => rank_single_ring_tied_pair_with_rule6(
            mol,
            element,
            root,
            orientation,
            signatures,
            &groups,
        ),
        _ => Err(CipAssignmentIssue::UnresolvedPriority { element }),
    }
}

fn rank_single_ring_tied_pair_with_rule6(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    orientation: TetrahedralOrientation,
    signatures: &[(StereoCarrier, LigandSignature)],
    groups: &[Vec<usize>],
) -> CipResult<RankedCarriers> {
    let tied_groups = groups
        .iter()
        .filter(|group| group.len() > 1)
        .collect::<Vec<_>>();
    if tied_groups.len() != 1 || tied_groups[0].len() != 2 {
        return Err(CipAssignmentIssue::UnresolvedPriority { element });
    }
    let left = carrier_rule6_atom(signatures[tied_groups[0][0]].0)
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let right = carrier_rule6_atom(signatures[tied_groups[0][1]].0)
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let Some(path) = shortest_path_excluding_root(mol, left, right, root) else {
        return Err(CipAssignmentIssue::UnresolvedPriority { element });
    };
    let path_length = path.len().saturating_sub(1);
    let tied_pair_descriptor_class = tied_groups[0]
        .iter()
        .filter_map(|index| tree_descriptor_class(&signatures[*index].1.root))
        .max();
    let outside_tied_pair_descriptor_class = groups
        .iter()
        .filter(|group| !std::ptr::eq(*group, tied_groups[0]))
        .flatten()
        .filter_map(|index| tree_descriptor_class(&signatures[*index].1.root))
        .max();
    let tied_pair_descriptor_refs_match =
        descriptor_ref_counts(&signatures[tied_groups[0][0]].1.root)
            == descriptor_ref_counts(&signatures[tied_groups[0][1]].1.root);
    let reference = if tied_pair_descriptor_class.is_some() {
        if tied_pair_descriptor_refs_match {
            return Err(CipAssignmentIssue::UnresolvedPriority { element });
        }
        if left.raw() >= right.raw() {
            left
        } else {
            right
        }
    } else {
        match outside_tied_pair_descriptor_class {
            Some(DescriptorClass::Absolute) => {
                if left.raw() >= right.raw() {
                    left
                } else {
                    right
                }
            }
            Some(DescriptorClass::Pseudo) => {
                if left.raw() <= right.raw() {
                    left
                } else {
                    right
                }
            }
            None if path_length == 2 => {
                match path
                    .get(1)
                    .and_then(|center| tetrahedral_orientation_for_center(mol, *center))
                {
                    Some(other_orientation) if other_orientation != orientation => right,
                    _ => left,
                }
            }
            None if mol.stereo_elements().all(|(id, _)| id == element)
                && ring_path_is_unsubstituted_bridge(mol, &path, root) =>
            {
                return Err(CipAssignmentIssue::UnresolvedPriority { element });
            }
            None if left.raw() >= right.raw() => left,
            None => right,
        }
    };
    let mut ranked = rank_carrier_signatures(element, signatures, Some(reference))?;
    ranked.pseudo_asymmetric_ordering = !matches!(
        (
            tied_pair_descriptor_class,
            outside_tied_pair_descriptor_class,
        ),
        (Some(DescriptorClass::Absolute), _) | (None, Some(DescriptorClass::Absolute))
    );
    Ok(ranked)
}

fn descriptor_ref_counts(tree: &LigandTree) -> (usize, usize) {
    let own = match tree.priority.descriptor.and_then(descriptor_ref) {
        Some(DescriptorRef::R) => (1, 0),
        Some(DescriptorRef::S) => (0, 1),
        None => (0, 0),
    };
    tree.children
        .iter()
        .map(descriptor_ref_counts)
        .fold(own, |left, right| (left.0 + right.0, left.1 + right.1))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum DescriptorClass {
    Pseudo,
    Absolute,
}

fn tree_descriptor_class(tree: &LigandTree) -> Option<DescriptorClass> {
    let own = tree.priority.descriptor.and_then(descriptor_class);
    tree.children
        .iter()
        .filter_map(tree_descriptor_class)
        .chain(own)
        .max()
}

fn descriptor_class(descriptor: StereoDescriptor) -> Option<DescriptorClass> {
    match descriptor {
        StereoDescriptor::R
        | StereoDescriptor::S
        | StereoDescriptor::M
        | StereoDescriptor::P
        | StereoDescriptor::SeqCis
        | StereoDescriptor::SeqTrans => Some(DescriptorClass::Absolute),
        StereoDescriptor::LowerR | StereoDescriptor::LowerS => Some(DescriptorClass::Pseudo),
        StereoDescriptor::E | StereoDescriptor::Z => None,
    }
}

fn tetrahedral_orientation_for_center(
    mol: &Molecule,
    center: AtomId,
) -> Option<TetrahedralOrientation> {
    mol.stereo_elements()
        .find_map(|(_, element)| match &element.kind {
            StereoElementKind::Tetrahedral(stereo) if stereo.center == center => {
                Some(stereo.orientation)
            }
            _ => None,
        })
}

fn shortest_path_excluding_root(
    mol: &Molecule,
    left: AtomId,
    right: AtomId,
    root: AtomId,
) -> Option<Vec<AtomId>> {
    let mut seen = Vec::new();
    let mut queue = VecDeque::from([(left, vec![left])]);
    while let Some((atom, path)) = queue.pop_front() {
        if atom == right {
            return Some(path);
        }
        if atom == root || seen.contains(&atom) {
            continue;
        }
        seen.push(atom);
        if let Ok(incident) = mol.incident_bonds(atom) {
            for (_, bond) in incident {
                let neighbor = bond.other_atom(atom);
                if neighbor != root && !seen.contains(&neighbor) {
                    let mut next_path = path.clone();
                    next_path.push(neighbor);
                    queue.push_back((neighbor, next_path));
                }
            }
        }
    }
    None
}

fn ring_path_is_unsubstituted_bridge(mol: &Molecule, path: &[AtomId], root: AtomId) -> bool {
    path.iter().all(|atom| {
        mol.incident_bonds(*atom)
            .map(|incident| {
                incident.into_iter().all(|(_, bond)| {
                    let neighbor = bond.other_atom(*atom);
                    neighbor == root || path.contains(&neighbor)
                })
            })
            .unwrap_or(false)
    })
}

fn rank_s4_tetrahedral_signatures_with_rule6(
    element: StereoElementId,
    signatures: &[(StereoCarrier, LigandSignature)],
    group: &[usize],
) -> CipResult<RankedCarriers> {
    let mut stable_ranking: Option<RankedCarriers> = None;
    for index in group {
        let Some(reference) = carrier_rule6_atom(signatures[*index].0) else {
            continue;
        };
        let ranking = match rank_carrier_signatures(element, signatures, Some(reference)) {
            Ok(ranking) => ranking,
            Err(CipAssignmentIssue::UnresolvedPriority { .. }) => continue,
            Err(issue) => return Err(issue),
        };
        if let Some(stable) = &stable_ranking {
            if carrier_permutation_is_odd(&stable.carriers, &ranking.carriers).unwrap_or(true) {
                return Err(CipAssignmentIssue::UnresolvedPriority { element });
            }
        } else {
            stable_ranking = Some(ranking);
        }
    }
    stable_ranking.ok_or(CipAssignmentIssue::UnresolvedPriority { element })
}

fn grouped_signature_indices(signatures: &[(StereoCarrier, LigandSignature)]) -> Vec<Vec<usize>> {
    let mut indices = (0..signatures.len()).collect::<Vec<_>>();
    indices.sort_by(|left, right| signatures[*right].1.compare(&signatures[*left].1));
    let mut groups: Vec<Vec<usize>> = Vec::new();
    for index in indices {
        if let Some(last) = groups.last_mut() {
            if signatures[last[0]].1.compare(&signatures[index].1) == Ordering::Equal {
                last.push(index);
                continue;
            }
        }
        groups.push(vec![index]);
    }
    groups
}

fn carrier_rule6_atom(carrier: StereoCarrier) -> Option<AtomId> {
    match carrier {
        StereoCarrier::Atom(atom) => Some(atom),
        StereoCarrier::ImplicitHydrogen | StereoCarrier::ImplicitLonePair => None,
    }
}

fn carrier_permutation_is_odd(left: &[StereoCarrier], right: &[StereoCarrier]) -> Option<bool> {
    if left.len() != right.len() {
        return None;
    }
    let mut positions = Vec::with_capacity(left.len());
    for carrier in left {
        positions.push(right.iter().position(|candidate| candidate == carrier)?);
    }
    Some(!permutation_is_even(&positions))
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RankedCarriers {
    carriers: Vec<StereoCarrier>,
    pseudo_asymmetric_ordering: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LigandSignature {
    root: LigandTree,
}

impl LigandSignature {
    fn compare(&self, other: &Self) -> Ordering {
        self.compare_with_flags(other).ordering
    }

    fn compare_with_flags(&self, other: &Self) -> LigandComparison {
        self.compare_with_rule6_reference(other, None)
    }

    fn compare_with_rule6_reference(
        &self,
        other: &Self,
        rule6_reference: Option<AtomId>,
    ) -> LigandComparison {
        self.root
            .compare_with_rule6_reference(&other.root, rule6_reference)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LigandTree {
    priority: NodePriority,
    children: Vec<LigandTree>,
}

impl LigandTree {
    fn compare_with_rule6_reference(
        &self,
        other: &Self,
        rule6_reference: Option<AtomId>,
    ) -> LigandComparison {
        for rule in [
            SequenceRule::Rule1a,
            SequenceRule::Rule1b,
            SequenceRule::Rule2,
            SequenceRule::Rule3,
            SequenceRule::Rule4a,
            SequenceRule::Rule4b,
            SequenceRule::Rule4c,
            SequenceRule::Rule5,
            SequenceRule::Rule6,
        ] {
            let comparison = self.compare_by_sequence_rule(other, rule, rule6_reference);
            if comparison.ordering != Ordering::Equal {
                return comparison;
            }
        }
        LigandComparison::equal()
    }

    fn compare_by_sequence_rule(
        &self,
        other: &Self,
        rule: SequenceRule,
        rule6_reference: Option<AtomId>,
    ) -> LigandComparison {
        match rule {
            SequenceRule::Rule4b => self.rule4b_reference_comparison(other),
            SequenceRule::Rule5 => self.rule5_pair_comparison(other),
            _ => LigandComparison::from_ordering(self.recursive_compare(
                other,
                rule,
                rule6_reference,
            )),
        }
    }

    fn recursive_compare(
        &self,
        other: &Self,
        rule: SequenceRule,
        rule6_reference: Option<AtomId>,
    ) -> Ordering {
        let priority = self
            .priority
            .compare_by_rule(&other.priority, rule, rule6_reference);
        if priority != Ordering::Equal {
            return priority;
        }

        let mut queue = vec![(self, other)];
        let mut position = 0usize;
        while position < queue.len() {
            let (left, right) = queue[position];
            position += 1;

            let left_shallow = left.children_sorted_by_rule(rule, false, rule6_reference);
            let right_shallow = right.children_sorted_by_rule(rule, false, rule6_reference);
            let shallow =
                compare_child_priorities(&left_shallow, &right_shallow, rule, rule6_reference);
            if shallow != Ordering::Equal {
                return shallow;
            }

            let left_deep = left.children_sorted_by_rule(rule, true, rule6_reference);
            let right_deep = right.children_sorted_by_rule(rule, true, rule6_reference);
            let deep = compare_child_priorities(&left_deep, &right_deep, rule, rule6_reference);
            if deep != Ordering::Equal {
                return deep;
            }
            for (left_child, right_child) in left_deep.into_iter().zip(right_deep) {
                queue.push((left_child, right_child));
            }
        }
        Ordering::Equal
    }

    fn compare_for_rule5_pairlist(&self, other: &Self, reference: DescriptorRef) -> Ordering {
        for rule in [
            SequenceRule::Rule1a,
            SequenceRule::Rule1b,
            SequenceRule::Rule2,
            SequenceRule::Rule3,
            SequenceRule::Rule4a,
            SequenceRule::Rule4b,
            SequenceRule::Rule4c,
            SequenceRule::Rule6,
        ] {
            let priority = self.compare_by_sequence_rule(other, rule, None).ordering;
            if priority != Ordering::Equal {
                return priority;
            }
        }
        rule5_reference_compare(
            self.priority.descriptor,
            other.priority.descriptor,
            reference,
        )
    }

    fn compare_through_rule4b(&self, other: &Self) -> Ordering {
        for rule in [
            SequenceRule::Rule1a,
            SequenceRule::Rule1b,
            SequenceRule::Rule2,
            SequenceRule::Rule3,
            SequenceRule::Rule4a,
            SequenceRule::Rule4b,
        ] {
            let priority = self.compare_by_sequence_rule(other, rule, None).ordering;
            if priority != Ordering::Equal {
                return priority;
            }
        }
        Ordering::Equal
    }

    fn compare_with_reference(&self, other: &Self, reference: DescriptorRef) -> Ordering {
        self.compare_without_rule4b_or_rule5(other)
            .then_with(|| self.fixed_reference_compare(other, reference))
    }

    fn compare_without_rule4b_or_rule5(&self, other: &Self) -> Ordering {
        for rule in [
            SequenceRule::Rule1a,
            SequenceRule::Rule1b,
            SequenceRule::Rule2,
            SequenceRule::Rule3,
            SequenceRule::Rule4a,
            SequenceRule::Rule4c,
        ] {
            let priority = self.recursive_compare(other, rule, None);
            if priority != Ordering::Equal {
                return priority;
            }
        }
        Ordering::Equal
    }

    fn fixed_reference_compare(&self, other: &Self, reference: DescriptorRef) -> Ordering {
        let priority = fixed_reference_priority(self.priority.descriptor, reference).cmp(
            &fixed_reference_priority(other.priority.descriptor, reference),
        );
        if priority != Ordering::Equal {
            return priority;
        }

        let mut queue = vec![(self, other)];
        let mut position = 0usize;
        while position < queue.len() {
            let (left, right) = queue[position];
            position += 1;

            let left_children = left.children_sorted_by_reference(reference);
            let right_children = right.children_sorted_by_reference(reference);
            for (left_child, right_child) in left_children.iter().zip(&right_children) {
                let priority =
                    fixed_reference_priority(left_child.priority.descriptor, reference).cmp(
                        &fixed_reference_priority(right_child.priority.descriptor, reference),
                    );
                if priority != Ordering::Equal {
                    return priority;
                }
            }
            let length = left_children.len().cmp(&right_children.len());
            if length != Ordering::Equal {
                return length;
            }
            for (left_child, right_child) in left_children.into_iter().zip(right_children) {
                queue.push((left_child, right_child));
            }
        }
        Ordering::Equal
    }

    fn rule4b_reference_comparison(&self, other: &Self) -> LigandComparison {
        let left_refs = self.rule4b_reference_descriptors();
        let right_refs = other.rule4b_reference_descriptors();
        if left_refs.is_empty() || right_refs.is_empty() || left_refs.len() != right_refs.len() {
            return LigandComparison::equal();
        }

        if left_refs.len() == 1 {
            return LigandComparison::from_ordering(self.compare_pairs_for_references(
                other,
                left_refs[0],
                right_refs[0],
            ));
        }

        let mut left_lists = left_refs
            .iter()
            .copied()
            .map(|reference| DescriptorPairList::collect_with_reference(self, reference))
            .collect::<Vec<_>>();
        let mut right_lists = right_refs
            .iter()
            .copied()
            .map(|reference| DescriptorPairList::collect_with_reference(other, reference))
            .collect::<Vec<_>>();
        left_lists.sort_by(|left, right| right.compare_to(left));
        right_lists.sort_by(|left, right| right.compare_to(left));
        for (left, right) in left_lists.iter().zip(&right_lists) {
            let comparison = left.compare_to(right);
            if comparison != Ordering::Equal {
                return LigandComparison::from_ordering(comparison);
            }
        }
        LigandComparison::equal()
    }

    fn rule4b_reference_descriptors(&self) -> Vec<DescriptorRef> {
        let mut level = vec![vec![self]];
        while !level.is_empty() {
            for group in &level {
                if let Some(reference) = reference_descriptor_for_group(group) {
                    return reference;
                }
            }
            level = next_reference_level(&level);
        }
        Vec::new()
    }

    fn compare_pairs_for_references(
        &self,
        other: &Self,
        left_reference: DescriptorRef,
        right_reference: DescriptorRef,
    ) -> Ordering {
        let mut left_queue = vec![self];
        let mut right_queue = vec![other];
        let mut position = 0usize;
        while position < left_queue.len() && position < right_queue.len() {
            let left = left_queue[position];
            let right = right_queue[position];
            position += 1;

            let left_like = descriptor_ref_matches(left.priority.descriptor, left_reference);
            let right_like = descriptor_ref_matches(right.priority.descriptor, right_reference);
            match (left_like, right_like) {
                (true, false) => return Ordering::Greater,
                (false, true) => return Ordering::Less,
                _ => {}
            }

            left_queue.extend(left.children_sorted_by_reference(left_reference));
            right_queue.extend(right.children_sorted_by_reference(right_reference));
        }
        Ordering::Equal
    }

    fn rule5_pair_comparison(&self, other: &Self) -> LigandComparison {
        let left_r = DescriptorPairList::collect(self, DescriptorRef::R);
        let right_r = DescriptorPairList::collect(other, DescriptorRef::R);
        let left_s = DescriptorPairList::collect(self, DescriptorRef::S);
        let right_s = DescriptorPairList::collect(other, DescriptorRef::S);

        let cmp_r = left_r.compare_to(&right_r);
        let cmp_s = left_s.compare_to(&right_s);
        match cmp_r {
            Ordering::Less => LigandComparison::new(Ordering::Less, cmp_s != Ordering::Less),
            Ordering::Greater => {
                LigandComparison::new(Ordering::Greater, cmp_s != Ordering::Greater)
            }
            Ordering::Equal => LigandComparison::equal(),
        }
    }

    fn children_sorted_by_rule(
        &self,
        rule: SequenceRule,
        deep: bool,
        rule6_reference: Option<AtomId>,
    ) -> Vec<&LigandTree> {
        let mut children = self.children.iter().collect::<Vec<_>>();
        if deep {
            children.sort_by(|left, right| right.recursive_compare(left, rule, rule6_reference));
        } else {
            children.sort_by(|left, right| {
                right
                    .priority
                    .compare_by_rule(&left.priority, rule, rule6_reference)
                    .then_with(|| right.priority.compare_shallow(&left.priority))
            });
        }
        children
    }

    fn children_sorted_by_reference(&self, reference: DescriptorRef) -> Vec<&LigandTree> {
        let mut children = self.children.iter().collect::<Vec<_>>();
        children.sort_by(|left, right| right.compare_with_reference(left, reference));
        children
    }

    fn children_grouped_through_rule4b(&self) -> Vec<Vec<&LigandTree>> {
        let mut children = self.children.iter().collect::<Vec<_>>();
        children.sort_by(|left, right| right.compare_through_rule4b(left));

        let mut groups: Vec<Vec<&LigandTree>> = Vec::new();
        for child in children {
            if let Some(last) = groups.last_mut() {
                if last[0].compare_through_rule4b(child) == Ordering::Equal {
                    last.push(child);
                    continue;
                }
            }
            groups.push(vec![child]);
        }
        groups
    }
}

fn compare_child_priorities(
    left: &[&LigandTree],
    right: &[&LigandTree],
    rule: SequenceRule,
    rule6_reference: Option<AtomId>,
) -> Ordering {
    for (left_child, right_child) in left.iter().zip(right) {
        let priority =
            left_child
                .priority
                .compare_by_rule(&right_child.priority, rule, rule6_reference);
        if priority != Ordering::Equal {
            return priority;
        }
    }
    left.len().cmp(&right.len())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SequenceRule {
    Rule1a,
    Rule1b,
    Rule2,
    Rule3,
    Rule4a,
    Rule4b,
    Rule4c,
    Rule5,
    Rule6,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct LigandComparison {
    ordering: Ordering,
    pseudo_asymmetric: bool,
}

impl LigandComparison {
    fn new(ordering: Ordering, pseudo_asymmetric: bool) -> Self {
        Self {
            ordering,
            pseudo_asymmetric,
        }
    }

    fn from_ordering(ordering: Ordering) -> Self {
        Self::new(ordering, false)
    }

    fn equal() -> Self {
        Self::from_ordering(Ordering::Equal)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct NodePriority {
    atomic_number: AtomicNumberFraction,
    rule1b: u32,
    rule2_mass: Rule2Mass,
    descriptor: Option<StereoDescriptor>,
    rule6_atom: Option<AtomId>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct AuxDescriptorKey {
    element: StereoElementId,
    path: Vec<AtomId>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct AuxOccurrence {
    key: AuxDescriptorKey,
    node: usize,
    distance: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AuxiliaryDescriptorMode {
    Disabled,
    Collect,
    Precomputed,
}

#[derive(Debug, Clone)]
struct DescriptorContext {
    skipped: Vec<StereoElementId>,
    auxiliary_mode: AuxiliaryDescriptorMode,
    aux_labels: Rc<RefCell<HashMap<AuxDescriptorKey, Option<StereoDescriptor>>>>,
    aux_occurrences: Rc<RefCell<Vec<AuxOccurrence>>>,
}

impl DescriptorContext {
    fn new(skip: StereoElementId, auxiliary_mode: AuxiliaryDescriptorMode) -> Self {
        Self {
            skipped: vec![skip],
            auxiliary_mode,
            aux_labels: Rc::new(RefCell::new(HashMap::new())),
            aux_occurrences: Rc::new(RefCell::new(Vec::new())),
        }
    }

    fn skips(&self, element: StereoElementId) -> bool {
        self.skipped.contains(&element)
    }

    fn with_skip(&self, element: StereoElementId) -> Self {
        let mut skipped = self.skipped.clone();
        skipped.push(element);
        Self {
            skipped,
            auxiliary_mode: self.auxiliary_mode,
            aux_labels: Rc::clone(&self.aux_labels),
            aux_occurrences: Rc::clone(&self.aux_occurrences),
        }
    }

    fn with_mode(&self, auxiliary_mode: AuxiliaryDescriptorMode) -> Self {
        Self {
            skipped: self.skipped.clone(),
            auxiliary_mode,
            aux_labels: Rc::clone(&self.aux_labels),
            aux_occurrences: Rc::clone(&self.aux_occurrences),
        }
    }
}

struct LigandBuildContext<'a> {
    mol: &'a Molecule,
    element: StereoElementId,
    descriptor_context: &'a DescriptorContext,
    options: CipAssignmentOptions,
    atomic_number_fractions: &'a [AtomicNumberFraction],
}

#[derive(Debug, Clone)]
struct AuxiliaryGraph {
    nodes: Vec<AuxiliaryGraphNode>,
}

#[derive(Debug, Clone)]
struct AuxiliaryGraphNode {
    node: LigandNode,
    parent: Option<usize>,
    children: Vec<usize>,
    depth: usize,
}

impl NodePriority {
    fn compare_shallow(&self, other: &Self) -> Ordering {
        self.atomic_number
            .cmp(&other.atomic_number)
            .then_with(|| self.rule1b.cmp(&other.rule1b))
            .then_with(|| self.rule2_mass.compare(other.rule2_mass))
    }

    fn compare_by_rule(
        &self,
        other: &Self,
        rule: SequenceRule,
        rule6_reference: Option<AtomId>,
    ) -> Ordering {
        match rule {
            SequenceRule::Rule1a => self.atomic_number.cmp(&other.atomic_number),
            SequenceRule::Rule1b => self.rule1b.cmp(&other.rule1b),
            SequenceRule::Rule2 => self.rule2_mass.compare(other.rule2_mass),
            SequenceRule::Rule3 => rule3_descriptor_priority(self.descriptor)
                .cmp(&rule3_descriptor_priority(other.descriptor)),
            SequenceRule::Rule4a => rule4a_descriptor_priority(self.descriptor)
                .cmp(&rule4a_descriptor_priority(other.descriptor)),
            SequenceRule::Rule4b => Ordering::Equal,
            SequenceRule::Rule4c => rule4c_descriptor_priority(self.descriptor)
                .cmp(&rule4c_descriptor_priority(other.descriptor)),
            SequenceRule::Rule5 => Ordering::Equal,
            SequenceRule::Rule6 => rule6_priority(self.rule6_atom, rule6_reference)
                .cmp(&rule6_priority(other.rule6_atom, rule6_reference)),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Rule2Mass {
    scaled_mass: u32,
    isotope_indicated: bool,
}

impl Rule2Mass {
    const ZERO: Self = Self {
        scaled_mass: 0,
        isotope_indicated: false,
    };

    fn natural(atomic_number: u8) -> Self {
        Self {
            scaled_mass: natural_atomic_weight_rank(atomic_number),
            isotope_indicated: false,
        }
    }

    fn isotope(mass_number: u16) -> Self {
        Self {
            scaled_mass: u32::from(mass_number).saturating_mul(ATOMIC_WEIGHT_SCALE),
            isotope_indicated: true,
        }
    }

    fn compare(self, other: Self) -> Ordering {
        if !self.isotope_indicated && !other.isotope_indicated {
            Ordering::Equal
        } else {
            self.scaled_mass.cmp(&other.scaled_mass)
        }
    }
}

const ATOMIC_WEIGHT_SCALE: u32 = 1_000;

const STANDARD_ATOMIC_WEIGHTS_MILLI: [u32; 119] = [
    0, 1008, 4003, 6941, 9012, 10812, 12011, 14007, 15999, 18998, 20180, 22990, 24305, 26982,
    28086, 30974, 32067, 35453, 39948, 39098, 40078, 44956, 47867, 50944, 51996, 54938, 55845,
    58933, 58693, 63546, 65390, 69723, 72610, 74922, 78960, 79904, 83800, 85468, 87620, 88906,
    91224, 92906, 95940, 98000, 101070, 102906, 106420, 107868, 112412, 114818, 118711, 121760,
    127600, 126904, 131290, 132905, 137328, 138906, 140116, 140908, 144240, 145000, 150360, 151964,
    157250, 158925, 162500, 164930, 167260, 168934, 173040, 174967, 178490, 180948, 183840, 186207,
    190230, 192217, 195078, 196967, 200590, 204383, 207200, 208980, 209000, 210000, 222000, 223000,
    226000, 227000, 232038, 231036, 238029, 237000, 244000, 243000, 247000, 247000, 251000, 252000,
    257000, 258000, 259000, 262000, 267000, 268000, 269000, 270000, 269000, 278000, 281000, 281000,
    285000, 284000, 289000, 288000, 293000, 292000, 294000,
];

fn natural_atomic_weight_rank(atomic_number: u8) -> u32 {
    STANDARD_ATOMIC_WEIGHTS_MILLI
        .get(usize::from(atomic_number))
        .copied()
        .unwrap_or_else(|| u32::from(atomic_number).saturating_mul(ATOMIC_WEIGHT_SCALE))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct AtomicNumberFraction {
    numerator: u32,
    denominator: u32,
}

impl AtomicNumberFraction {
    const ZERO: Self = Self {
        numerator: 0,
        denominator: 1,
    };

    const HYDROGEN: Self = Self {
        numerator: 1,
        denominator: 1,
    };

    fn element(atomic_number: u8) -> Self {
        Self {
            numerator: u32::from(atomic_number),
            denominator: 1,
        }
    }

    fn new(numerator: u32, denominator: u32) -> Self {
        if denominator == 0 {
            return Self::ZERO;
        }
        let divisor = gcd(numerator, denominator);
        Self {
            numerator: numerator / divisor,
            denominator: denominator / divisor,
        }
    }
}

impl Ord for AtomicNumberFraction {
    fn cmp(&self, other: &Self) -> Ordering {
        u64::from(self.numerator)
            .saturating_mul(u64::from(other.denominator))
            .cmp(&u64::from(other.numerator).saturating_mul(u64::from(self.denominator)))
    }
}

impl PartialOrd for AtomicNumberFraction {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

fn gcd(mut left: u32, mut right: u32) -> u32 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left.max(1)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DescriptorRef {
    R,
    S,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct DescriptorPairList {
    reference: DescriptorRef,
    descriptors: Vec<DescriptorRef>,
}

impl DescriptorPairList {
    fn collect(root: &LigandTree, reference: DescriptorRef) -> Self {
        let mut list = Self {
            reference,
            descriptors: vec![reference],
        };
        let mut queue = vec![root];
        let mut position = 0usize;
        while position < queue.len() {
            let node = queue[position];
            position += 1;
            list.add(node.priority.descriptor);

            let mut children = node.children.iter().collect::<Vec<_>>();
            children.sort_by(|left, right| right.compare_for_rule5_pairlist(left, reference));
            queue.extend(children);
        }
        list
    }

    fn collect_with_reference(root: &LigandTree, reference: DescriptorRef) -> Self {
        let mut list = Self {
            reference,
            descriptors: vec![reference],
        };
        let mut queue = vec![root];
        let mut position = 0usize;
        while position < queue.len() {
            let node = queue[position];
            position += 1;
            list.add(node.priority.descriptor);
            queue.extend(node.children_sorted_by_reference(reference));
        }
        list
    }

    fn add(&mut self, descriptor: Option<StereoDescriptor>) {
        if let Some(reference) = descriptor.and_then(descriptor_ref) {
            self.descriptors.push(reference);
        }
    }

    fn compare_to(&self, other: &Self) -> Ordering {
        if self.descriptors.len() != other.descriptors.len() {
            return Ordering::Equal;
        }
        for (left, right) in self
            .descriptors
            .iter()
            .skip(1)
            .zip(other.descriptors.iter().skip(1))
        {
            let left_like = *left == self.reference;
            let right_like = *right == other.reference;
            match (left_like, right_like) {
                (true, false) => return Ordering::Greater,
                (false, true) => return Ordering::Less,
                _ => {}
            }
        }
        Ordering::Equal
    }
}

fn reference_descriptor_for_group(group: &[&LigandTree]) -> Option<Vec<DescriptorRef>> {
    let mut right = 0usize;
    let mut left = 0usize;
    for node in group {
        match node.priority.descriptor.and_then(descriptor_ref) {
            Some(DescriptorRef::R) => right += 1,
            Some(DescriptorRef::S) => left += 1,
            None => {}
        }
    }
    match right.cmp(&left) {
        Ordering::Greater => Some(vec![DescriptorRef::R]),
        Ordering::Less => Some(vec![DescriptorRef::S]),
        Ordering::Equal if right > 0 => Some(vec![DescriptorRef::R, DescriptorRef::S]),
        Ordering::Equal => None,
    }
}

fn next_reference_level<'a>(previous: &[Vec<&'a LigandTree>]) -> Vec<Vec<&'a LigandTree>> {
    let mut next = Vec::new();
    for group in previous {
        let mut grouped_children = Vec::new();
        let mut group_count = None;
        for node in group {
            let children = node.children_grouped_through_rule4b();
            if children.is_empty() {
                continue;
            }
            if let Some(expected) = group_count {
                if expected != children.len() {
                    return Vec::new();
                }
            } else {
                group_count = Some(children.len());
            }
            grouped_children.push(children);
        }
        let Some(group_count) = group_count else {
            continue;
        };
        for index in 0..group_count {
            let mut equivalent_nodes = Vec::new();
            for children in &grouped_children {
                equivalent_nodes.extend(children[index].iter().copied());
            }
            if !equivalent_nodes.is_empty() {
                next.push(equivalent_nodes);
            }
        }
    }
    next
}

fn descriptor_ref_matches(descriptor: Option<StereoDescriptor>, reference: DescriptorRef) -> bool {
    descriptor.and_then(descriptor_ref) == Some(reference)
}

fn fixed_reference_priority(descriptor: Option<StereoDescriptor>, reference: DescriptorRef) -> u8 {
    match descriptor.and_then(descriptor_ref) {
        Some(descriptor) if descriptor == reference => 2,
        Some(_) => 1,
        None => 0,
    }
}

fn rule5_reference_compare(
    left: Option<StereoDescriptor>,
    right: Option<StereoDescriptor>,
    reference: DescriptorRef,
) -> Ordering {
    match (
        left.and_then(descriptor_ref),
        right.and_then(descriptor_ref),
    ) {
        (Some(left), Some(right)) => {
            let left_like = left == reference;
            let right_like = right == reference;
            match (left_like, right_like) {
                (true, false) => Ordering::Greater,
                (false, true) => Ordering::Less,
                _ => Ordering::Equal,
            }
        }
        _ => Ordering::Equal,
    }
}

fn rule6_priority(atom: Option<AtomId>, reference: Option<AtomId>) -> u8 {
    match (atom, reference) {
        (Some(atom), Some(reference)) if atom == reference => 1,
        _ => 0,
    }
}

fn carrier_signature(
    mol: &Molecule,
    element: StereoElementId,
    descriptor_context: &DescriptorContext,
    carrier: StereoCarrier,
    root: AtomId,
    options: CipAssignmentOptions,
    atomic_number_fractions: &[AtomicNumberFraction],
) -> CipResult<LigandSignature> {
    let node = match carrier {
        StereoCarrier::Atom(atom) => LigandNode::Atom {
            atom,
            previous: Some(root),
            path: vec![root, atom],
            duplicate: None,
            terminal: false,
        },
        StereoCarrier::ImplicitHydrogen => LigandNode::Hydrogen,
        StereoCarrier::ImplicitLonePair => LigandNode::LonePair,
    };
    let mut visited_nodes = 0usize;
    let build_context = LigandBuildContext {
        mol,
        element,
        descriptor_context,
        options,
        atomic_number_fractions,
    };
    let root = ligand_tree(&build_context, node, 0, &mut visited_nodes)?;
    Ok(LigandSignature { root })
}

fn build_auxiliary_graph(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    options: CipAssignmentOptions,
    atomic_number_fractions: &[AtomicNumberFraction],
) -> CipResult<AuxiliaryGraph> {
    let root = LigandNode::Atom {
        atom: root,
        previous: None,
        path: vec![root],
        duplicate: None,
        terminal: false,
    };
    let mut graph = AuxiliaryGraph { nodes: Vec::new() };
    let mut visited_nodes = 0usize;
    let context = AuxiliaryGraphBuildContext {
        mol,
        element,
        options,
        atomic_number_fractions,
    };
    add_auxiliary_graph_node(&context, &mut graph, root, None, 0, &mut visited_nodes)?;
    Ok(graph)
}

struct AuxiliaryGraphBuildContext<'a> {
    mol: &'a Molecule,
    element: StereoElementId,
    options: CipAssignmentOptions,
    atomic_number_fractions: &'a [AtomicNumberFraction],
}

fn add_auxiliary_graph_node(
    context: &AuxiliaryGraphBuildContext<'_>,
    graph: &mut AuxiliaryGraph,
    node: LigandNode,
    parent: Option<usize>,
    depth: usize,
    visited_nodes: &mut usize,
) -> CipResult<usize> {
    *visited_nodes = visited_nodes.saturating_add(1);
    if *visited_nodes > context.options.max_nodes {
        return Err(CipAssignmentIssue::ResourceLimitExceeded {
            element: context.element,
            max_nodes: context.options.max_nodes,
        });
    }

    let index = graph.nodes.len();
    graph.nodes.push(AuxiliaryGraphNode {
        node: node.clone(),
        parent,
        children: Vec::new(),
        depth,
    });
    if depth < context.options.max_depth.saturating_add(1) {
        let mut child_nodes = Vec::new();
        node.extend(
            context.mol,
            context.atomic_number_fractions,
            &mut child_nodes,
        );
        for child in child_nodes {
            let child_index = add_auxiliary_graph_node(
                context,
                graph,
                child,
                Some(index),
                depth + 1,
                visited_nodes,
            )?;
            graph.nodes[index].children.push(child_index);
        }
    }
    Ok(index)
}

fn ligand_tree(
    context: &LigandBuildContext<'_>,
    node: LigandNode,
    depth: usize,
    visited_nodes: &mut usize,
) -> CipResult<LigandTree> {
    *visited_nodes = visited_nodes.saturating_add(1);
    if *visited_nodes > context.options.max_nodes {
        return Err(CipAssignmentIssue::ResourceLimitExceeded {
            element: context.element,
            max_nodes: context.options.max_nodes,
        });
    }
    let priority = node.priority(context);
    let mut children = Vec::new();
    if depth < context.options.max_depth {
        let mut child_nodes = Vec::new();
        node.extend(
            context.mol,
            context.atomic_number_fractions,
            &mut child_nodes,
        );
        for child in child_nodes {
            children.push(ligand_tree(context, child, depth + 1, visited_nodes)?);
        }
        children.sort_by(|left, right| right.priority.compare_shallow(&left.priority));
    }
    Ok(LigandTree { priority, children })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DuplicateNode {
    Bond {
        atomic_number: Option<AtomicNumberFraction>,
    },
    Ring {
        reference_depth: usize,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum LigandNode {
    Atom {
        atom: AtomId,
        previous: Option<AtomId>,
        path: Vec<AtomId>,
        duplicate: Option<DuplicateNode>,
        terminal: bool,
    },
    Hydrogen,
    LonePair,
}

impl LigandNode {
    fn priority(&self, context: &LigandBuildContext<'_>) -> NodePriority {
        NodePriority {
            atomic_number: self.atomic_number(context.mol, context.atomic_number_fractions),
            rule1b: self.rule1b_priority(),
            rule2_mass: self.rule2_mass(context.mol),
            descriptor: self.descriptor(context),
            rule6_atom: self.rule6_atom(),
        }
    }

    fn atomic_number(
        &self,
        mol: &Molecule,
        _atomic_number_fractions: &[AtomicNumberFraction],
    ) -> AtomicNumberFraction {
        match self {
            Self::Atom {
                duplicate:
                    Some(DuplicateNode::Bond {
                        atomic_number: Some(atomic_number),
                    }),
                ..
            } => *atomic_number,
            Self::Atom { atom, .. } => mol
                .atom(*atom)
                .ok()
                .map(|atom| AtomicNumberFraction::element(atom.element.atomic_number()))
                .unwrap_or(AtomicNumberFraction::ZERO),
            Self::Hydrogen => AtomicNumberFraction::HYDROGEN,
            Self::LonePair => AtomicNumberFraction::ZERO,
        }
    }

    fn rule1b_priority(&self) -> u32 {
        match self {
            Self::Atom {
                duplicate: Some(DuplicateNode::Ring { reference_depth }),
                ..
            } => ring_duplicate_priority(*reference_depth),
            Self::Atom { .. } | Self::Hydrogen | Self::LonePair => 0,
        }
    }

    fn rule2_mass(&self, mol: &Molecule) -> Rule2Mass {
        match self {
            Self::Atom {
                atom, duplicate, ..
            } => {
                if duplicate.is_some() {
                    Rule2Mass::ZERO
                } else {
                    let Ok(atom) = mol.atom(*atom) else {
                        return Rule2Mass::ZERO;
                    };
                    atom.isotope.map_or_else(
                        || Rule2Mass::natural(atom.element.atomic_number()),
                        Rule2Mass::isotope,
                    )
                }
            }
            Self::Hydrogen => Rule2Mass::natural(1),
            Self::LonePair => Rule2Mass::ZERO,
        }
    }

    fn descriptor(&self, context: &LigandBuildContext<'_>) -> Option<StereoDescriptor> {
        let Self::Atom {
            atom,
            path,
            duplicate: None,
            ..
        } = self
        else {
            return None;
        };
        atom_descriptor_for_ligand_node(context, *atom, path)
    }

    fn rule6_atom(&self) -> Option<AtomId> {
        match self {
            Self::Atom { atom, .. } => Some(*atom),
            Self::Hydrogen | Self::LonePair => None,
        }
    }

    fn extend(
        &self,
        mol: &Molecule,
        atomic_number_fractions: &[AtomicNumberFraction],
        next: &mut Vec<LigandNode>,
    ) {
        let Self::Atom {
            atom,
            previous,
            path,
            duplicate: _,
            terminal,
        } = self
        else {
            return;
        };
        if *terminal {
            return;
        }
        let Ok(payload) = mol.atom(*atom) else {
            return;
        };
        {
            for _ in 0..hydrogen_count(payload) {
                next.push(LigandNode::Hydrogen);
            }
        }
        let Ok(incident) = mol.incident_bonds(*atom) else {
            return;
        };
        for (_, bond) in incident {
            let neighbor = bond.other_atom(*atom);
            let duplicate_count =
                bond_duplicate_count_for_atom(payload, *atom, bond.order, atomic_number_fractions);
            let bond_duplicate_atomic_number =
                bond_duplicate_atomic_number(*atom, atomic_number_fractions);
            if Some(neighbor) == *previous {
                if path.first().copied() != Some(neighbor) {
                    for _ in 0..duplicate_count {
                        next.push(LigandNode::Atom {
                            atom: neighbor,
                            previous: Some(*atom),
                            path: Vec::new(),
                            duplicate: Some(DuplicateNode::Bond {
                                atomic_number: bond_duplicate_atomic_number,
                            }),
                            terminal: true,
                        });
                    }
                }
                continue;
            }
            if let Some(reference_depth) = path.iter().position(|id| *id == neighbor) {
                next.push(LigandNode::Atom {
                    atom: neighbor,
                    previous: Some(*atom),
                    path: Vec::new(),
                    duplicate: Some(DuplicateNode::Ring { reference_depth }),
                    terminal: true,
                });
                for _ in 0..duplicate_count {
                    next.push(LigandNode::Atom {
                        atom: neighbor,
                        previous: Some(*atom),
                        path: Vec::new(),
                        duplicate: Some(DuplicateNode::Bond {
                            atomic_number: bond_duplicate_atomic_number,
                        }),
                        terminal: true,
                    });
                }
            } else {
                let mut next_path = path.clone();
                next_path.push(neighbor);
                next.push(LigandNode::Atom {
                    atom: neighbor,
                    previous: Some(*atom),
                    path: next_path,
                    duplicate: None,
                    terminal: false,
                });
                for _ in 0..duplicate_count {
                    next.push(LigandNode::Atom {
                        atom: neighbor,
                        previous: Some(*atom),
                        path: Vec::new(),
                        duplicate: Some(DuplicateNode::Bond {
                            atomic_number: bond_duplicate_atomic_number,
                        }),
                        terminal: true,
                    });
                }
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MancudeAtomType {
    Cv4D3,
    Nv3D2,
    Nv4D3Plus,
    Nv2D2Minus,
    Cv3D3Minus,
    Ov3D2Plus,
    Other,
}

fn cip_atomic_number_fractions(mol: &Molecule) -> Vec<AtomicNumberFraction> {
    let mut fractions = vec![AtomicNumberFraction::ZERO; mol.atoms.len()];
    for (atom_id, atom) in mol.atoms() {
        fractions[atom_id.index()] = AtomicNumberFraction::element(atom.element.atomic_number());
    }

    let ring_membership = mol
        .ring_membership()
        .cloned()
        .unwrap_or_else(|| compute_ring_membership(mol).0);
    let mut types = seed_mancude_atom_types(mol, &ring_membership);
    if !types.iter().any(|atom_type| {
        matches!(
            atom_type,
            MancudeAtomType::Nv3D2
                | MancudeAtomType::Nv4D3Plus
                | MancudeAtomType::Nv2D2Minus
                | MancudeAtomType::Cv3D3Minus
                | MancudeAtomType::Ov3D2Plus
        )
    }) {
        return fractions;
    }

    relax_mancude_atom_types(mol, &mut types);
    let parts = mancude_parts(mol, &types, &ring_membership);
    apply_mancude_neighbor_averages(mol, &types, &parts, &mut fractions);
    fractions
}

fn seed_mancude_atom_types(
    mol: &Molecule,
    ring_membership: &RingMembership,
) -> Vec<MancudeAtomType> {
    let mut types = vec![MancudeAtomType::Other; mol.atoms.len()];
    for (atom_id, atom) in mol.atoms() {
        let mut bond_types = u32::from(hydrogen_count(atom));
        let mut in_ring = false;
        if let Ok(incident) = mol.incident_bonds(atom_id) {
            for (bond_id, bond) in incident {
                bond_types += match cip_bond_order(bond.order) {
                    1 => 0x0000_0001,
                    2 => 0x0000_0100,
                    _ => 0x0100_0000,
                };
                if ring_membership.bond_in_ring(bond_id) {
                    in_ring = true;
                }
            }
        }
        if !in_ring {
            continue;
        }
        types[atom_id.index()] =
            match (atom.element.atomic_number(), atom.formal_charge, bond_types) {
                (6 | 14 | 32, 0, 0x0102) => MancudeAtomType::Cv4D3,
                (6 | 14 | 32, -1, 0x0003) => MancudeAtomType::Cv3D3Minus,
                (7 | 15 | 33, 0, 0x0101) => MancudeAtomType::Nv3D2,
                (7 | 15 | 33, -1, 0x0002) => MancudeAtomType::Nv2D2Minus,
                (7 | 15 | 33, 1, 0x0102) => MancudeAtomType::Nv4D3Plus,
                (8, 1, 0x0101) => MancudeAtomType::Ov3D2Plus,
                _ => MancudeAtomType::Other,
            };
    }
    types
}

fn relax_mancude_atom_types(mol: &Molecule, types: &mut [MancudeAtomType]) {
    let mut counts = vec![0usize; mol.atoms.len()];
    let mut queue = Vec::new();
    for (atom_id, _) in mol.atoms() {
        for neighbor in atom_neighbors(mol, atom_id) {
            if types[neighbor.index()] != MancudeAtomType::Other {
                counts[atom_id.index()] += 1;
            }
        }
        if counts[atom_id.index()] == 1 {
            queue.push(atom_id);
        }
    }

    let mut position = 0usize;
    while position < queue.len() {
        let atom_id = queue[position];
        position += 1;
        if types[atom_id.index()] == MancudeAtomType::Other {
            continue;
        }
        types[atom_id.index()] = MancudeAtomType::Other;
        for neighbor in atom_neighbors(mol, atom_id) {
            counts[neighbor.index()] = counts[neighbor.index()].saturating_sub(1);
            if counts[neighbor.index()] == 1 {
                queue.push(neighbor);
            }
        }
    }
}

fn mancude_parts(
    mol: &Molecule,
    types: &[MancudeAtomType],
    ring_membership: &RingMembership,
) -> Vec<usize> {
    let mut parts = vec![0usize; mol.atoms.len()];
    let mut part = 0usize;
    for (atom_id, _) in mol.atoms() {
        if parts[atom_id.index()] != 0 || types[atom_id.index()] == MancudeAtomType::Other {
            continue;
        }
        part += 1;
        parts[atom_id.index()] = part;
        let mut stack = vec![atom_id];
        while let Some(current) = stack.pop() {
            if let Ok(incident) = mol.incident_bonds(current) {
                for (bond_id, bond) in incident {
                    if !ring_membership.bond_in_ring(bond_id) {
                        continue;
                    }
                    let neighbor = bond.other_atom(current);
                    if parts[neighbor.index()] == 0
                        && types[neighbor.index()] != MancudeAtomType::Other
                    {
                        parts[neighbor.index()] = part;
                        stack.push(neighbor);
                    }
                }
            }
        }
    }
    parts
}

fn apply_mancude_neighbor_averages(
    mol: &Molecule,
    types: &[MancudeAtomType],
    parts: &[usize],
    fractions: &mut [AtomicNumberFraction],
) {
    let mut resonance_parts = Vec::<usize>::new();
    for (atom_id, _) in mol.atoms() {
        let part = parts[atom_id.index()];
        if part == 0 {
            continue;
        }
        if matches!(
            types[atom_id.index()],
            MancudeAtomType::Cv3D3Minus | MancudeAtomType::Nv2D2Minus
        ) && !resonance_parts.contains(&part)
        {
            resonance_parts.push(part);
        }

        let mut numerator = 0u32;
        let mut denominator = 0u32;
        for neighbor in atom_neighbors(mol, atom_id) {
            if parts[neighbor.index()] == part {
                if let Ok(atom) = mol.atom(neighbor) {
                    numerator += u32::from(atom.element.atomic_number());
                    denominator += 1;
                }
            }
        }
        fractions[atom_id.index()] = AtomicNumberFraction::new(numerator, denominator);
    }

    for part in resonance_parts {
        let mut numerator = 0u32;
        let mut denominator = 0u32;
        for (index, fraction) in fractions.iter_mut().enumerate().take(mol.atoms.len()) {
            if parts.get(index).copied() != Some(part) {
                continue;
            }
            *fraction = AtomicNumberFraction::new(numerator, denominator);
            denominator += 1;
            let atom_id = AtomId::new(index as u32);
            if let Ok(incident) = mol.incident_bonds(atom_id) {
                for (_, bond) in incident {
                    let neighbor = bond.other_atom(atom_id);
                    if parts[neighbor.index()] == part {
                        let bond_order = cip_bond_order(bond.order);
                        if bond_order > 1 {
                            if let Ok(neighbor_atom) = mol.atom(neighbor) {
                                numerator += u32::from(bond_order.saturating_sub(1))
                                    * u32::from(neighbor_atom.element.atomic_number());
                            }
                        }
                    }
                }
            }
        }
    }
}

fn atom_neighbors(mol: &Molecule, atom_id: AtomId) -> Vec<AtomId> {
    mol.incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .map(|(_, bond)| bond.other_atom(atom_id))
        .collect()
}

fn hydrogen_count(atom: &Atom) -> u8 {
    atom.explicit_hydrogens
        .saturating_add(atom.implicit_hydrogens.unwrap_or(0))
}

fn bond_duplicate_count_for_atom(
    atom: &Atom,
    atom_id: AtomId,
    order: BondOrder,
    atomic_number_fractions: &[AtomicNumberFraction],
) -> usize {
    if atom.formal_charge < 0
        && atomic_number_fractions
            .get(atom_id.index())
            .is_some_and(|fraction| fraction.denominator > 1)
    {
        1
    } else {
        bond_order_duplicate_count(cip_bond_order(order))
    }
}

fn bond_duplicate_atomic_number(
    atom: AtomId,
    atomic_number_fractions: &[AtomicNumberFraction],
) -> Option<AtomicNumberFraction> {
    atomic_number_fractions
        .get(atom.index())
        .copied()
        .filter(|fraction| fraction.denominator > 1)
}

fn cip_bond_order(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single | BondOrder::Aromatic => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Zero | BondOrder::Dative => 0,
    }
}

fn bond_order_duplicate_count(order: u8) -> usize {
    match order {
        2 => 1,
        3 => 2,
        4 => 3,
        _ => 0,
    }
}

fn ring_duplicate_priority(reference_depth: usize) -> u32 {
    let depth = reference_depth.min(u32::MAX as usize) as u32;
    u32::MAX.saturating_sub(depth)
}

fn atom_descriptor_for_ligand_node(
    context: &LigandBuildContext<'_>,
    atom: AtomId,
    path: &[AtomId],
) -> Option<StereoDescriptor> {
    context.mol.stereo_elements().find_map(|(id, element)| {
        if context.descriptor_context.skips(id) {
            return None;
        }
        match &element.kind {
            StereoElementKind::Tetrahedral(stereo) if stereo.center == atom => {
                match context.descriptor_context.auxiliary_mode {
                    AuxiliaryDescriptorMode::Disabled => None,
                    AuxiliaryDescriptorMode::Collect => {
                        if path.last().copied() == Some(stereo.center) {
                            record_auxiliary_occurrence(context.descriptor_context, id, path);
                        }
                        None
                    }
                    AuxiliaryDescriptorMode::Precomputed => {
                        let key = AuxDescriptorKey {
                            element: id,
                            path: path.to_vec(),
                        };
                        let aux_labels = context.descriptor_context.aux_labels.borrow();
                        aux_labels
                            .get(&key)
                            .copied()
                            .flatten()
                            .or(element.descriptor)
                    }
                }
            }
            StereoElementKind::DoubleBond(stereo) => element.descriptor.and_then(|descriptor| {
                double_bond_descriptor_applies_to_node(stereo, descriptor, atom, path)
                    .then_some(descriptor)
            }),
            _ => None,
        }
    })
}

fn record_auxiliary_occurrence(
    context: &DescriptorContext,
    element: StereoElementId,
    path: &[AtomId],
) {
    context.aux_occurrences.borrow_mut().push(AuxOccurrence {
        key: AuxDescriptorKey {
            element,
            path: path.to_vec(),
        },
        node: 0,
        distance: path.len().saturating_sub(1),
    });
}

fn collect_auxiliary_occurrences_from_graph(
    mol: &Molecule,
    context: &DescriptorContext,
    graph: &AuxiliaryGraph,
) {
    for (node_index, graph_node) in graph.nodes.iter().enumerate() {
        let LigandNode::Atom {
            atom,
            path,
            duplicate: None,
            ..
        } = &graph_node.node
        else {
            continue;
        };
        for (element, stereo_element) in mol.stereo_elements() {
            if context.skips(element) {
                continue;
            }
            let StereoElementKind::Tetrahedral(stereo) = &stereo_element.kind else {
                continue;
            };
            if stereo.center == *atom {
                context.aux_occurrences.borrow_mut().push(AuxOccurrence {
                    key: AuxDescriptorKey {
                        element,
                        path: path.clone(),
                    },
                    node: node_index,
                    distance: graph_node.depth,
                });
            }
        }
    }
}

fn precompute_auxiliary_descriptors(
    mol: &Molecule,
    descriptor_context: &DescriptorContext,
    graph: &AuxiliaryGraph,
    options: CipAssignmentOptions,
    atomic_number_fractions: &[AtomicNumberFraction],
) {
    let mut occurrences = descriptor_context.aux_occurrences.borrow().clone();
    let mut seen = HashSet::new();
    occurrences.retain(|occurrence| seen.insert(occurrence.key.clone()));
    occurrences.sort_by(|left, right| {
        right
            .distance
            .cmp(&left.distance)
            .then_with(|| left.key.element.cmp(&right.key.element))
            .then_with(|| left.key.path.cmp(&right.key.path))
    });

    let mut position = 0usize;
    while position < occurrences.len() {
        let distance = occurrences[position].distance;
        let start = position;
        while position < occurrences.len() && occurrences[position].distance == distance {
            position += 1;
        }

        let mut batch = Vec::new();
        for occurrence in &occurrences[start..position] {
            if descriptor_context
                .aux_labels
                .borrow()
                .contains_key(&occurrence.key)
            {
                continue;
            }
            let descriptor = auxiliary_tetrahedral_descriptor_for_occurrence(
                mol,
                descriptor_context,
                graph,
                occurrence,
                options,
                atomic_number_fractions,
            );
            batch.push((occurrence.key.clone(), descriptor));
        }

        let mut aux_labels = descriptor_context.aux_labels.borrow_mut();
        for (key, descriptor) in batch {
            aux_labels.insert(key, descriptor);
        }
    }
}

fn auxiliary_tetrahedral_descriptor_for_occurrence(
    mol: &Molecule,
    descriptor_context: &DescriptorContext,
    graph: &AuxiliaryGraph,
    occurrence: &AuxOccurrence,
    options: CipAssignmentOptions,
    atomic_number_fractions: &[AtomicNumberFraction],
) -> Option<StereoDescriptor> {
    let element = mol.stereo_element(occurrence.key.element).ok()?;
    let StereoElementKind::Tetrahedral(stereo) = &element.kind else {
        return None;
    };
    if occurrence.key.path.last().copied() != Some(stereo.center) {
        return None;
    }
    let aux_descriptor_context = descriptor_context
        .with_skip(occurrence.key.element)
        .with_mode(AuxiliaryDescriptorMode::Precomputed);
    let aux_context = LigandBuildContext {
        mol,
        element: occurrence.key.element,
        descriptor_context: &aux_descriptor_context,
        options,
        atomic_number_fractions,
    };
    let signatures =
        auxiliary_tetrahedral_signatures(&aux_context, graph, occurrence.node, stereo).ok()?;
    let ranked = match rank_carrier_signatures(occurrence.key.element, &signatures, None) {
        Ok(ranked) => ranked,
        Err(CipAssignmentIssue::UnresolvedPriority { .. }) if stereo.carriers.len() == 4 => {
            rank_tetrahedral_signatures_with_rule6(
                mol,
                occurrence.key.element,
                stereo.center,
                &signatures,
                stereo.orientation,
                true,
            )
            .ok()?
        }
        Err(_) => return None,
    };
    tetrahedral_descriptor_from_ranked(occurrence.key.element, stereo, &ranked).ok()
}

fn auxiliary_tetrahedral_signatures(
    context: &LigandBuildContext<'_>,
    graph: &AuxiliaryGraph,
    root: usize,
    stereo: &TetrahedralStereo,
) -> CipResult<Vec<(StereoCarrier, LigandSignature)>> {
    stereo
        .carriers
        .iter()
        .copied()
        .map(|carrier| {
            auxiliary_carrier_signature(context, graph, root, carrier)
                .map(|signature| (carrier, signature))
        })
        .collect()
}

fn auxiliary_carrier_signature(
    context: &LigandBuildContext<'_>,
    graph: &AuxiliaryGraph,
    root: usize,
    carrier: StereoCarrier,
) -> CipResult<LigandSignature> {
    let root = match carrier {
        StereoCarrier::Atom(atom) => {
            let Some(node) = outgoing_auxiliary_graph_nodes(graph, root, root)
                .into_iter()
                .find(|node| auxiliary_graph_node_matches_atom(graph, *node, atom))
            else {
                return Err(CipAssignmentIssue::UnresolvedPriority {
                    element: context.element,
                });
            };
            let mut visited_nodes = 0usize;
            ligand_tree_from_auxiliary_graph(context, graph, root, node, 0, &mut visited_nodes)?
        }
        StereoCarrier::ImplicitHydrogen => {
            let Some(node) = outgoing_auxiliary_graph_nodes(graph, root, root)
                .into_iter()
                .find(|node| matches!(graph.nodes[*node].node, LigandNode::Hydrogen))
            else {
                return Err(CipAssignmentIssue::UnresolvedPriority {
                    element: context.element,
                });
            };
            let mut visited_nodes = 0usize;
            ligand_tree_from_auxiliary_graph(context, graph, root, node, 0, &mut visited_nodes)?
        }
        StereoCarrier::ImplicitLonePair => LigandTree {
            priority: LigandNode::LonePair.priority(context),
            children: Vec::new(),
        },
    };
    Ok(LigandSignature { root })
}

fn ligand_tree_from_auxiliary_graph(
    context: &LigandBuildContext<'_>,
    graph: &AuxiliaryGraph,
    root: usize,
    node: usize,
    depth: usize,
    visited_nodes: &mut usize,
) -> CipResult<LigandTree> {
    *visited_nodes = visited_nodes.saturating_add(1);
    if *visited_nodes > context.options.max_nodes {
        return Err(CipAssignmentIssue::ResourceLimitExceeded {
            element: context.element,
            max_nodes: context.options.max_nodes,
        });
    }
    let priority = graph.nodes[node].node.priority(context);
    let mut children = Vec::new();
    if depth < context.options.max_depth {
        for child in outgoing_auxiliary_graph_nodes(graph, root, node) {
            children.push(ligand_tree_from_auxiliary_graph(
                context,
                graph,
                root,
                child,
                depth + 1,
                visited_nodes,
            )?);
        }
        children.sort_by(|left, right| right.priority.compare_shallow(&left.priority));
    }
    Ok(LigandTree { priority, children })
}

fn auxiliary_graph_node_matches_atom(graph: &AuxiliaryGraph, node: usize, atom: AtomId) -> bool {
    matches!(
        &graph.nodes[node].node,
        LigandNode::Atom {
            atom: node_atom,
            duplicate: None,
            ..
        } if *node_atom == atom
    )
}

fn outgoing_auxiliary_graph_nodes(graph: &AuxiliaryGraph, root: usize, node: usize) -> Vec<usize> {
    let mut path = Vec::new();
    let mut cursor = Some(root);
    while let Some(current) = cursor {
        path.push(current);
        cursor = graph.nodes[current].parent;
    }
    let path_position = path.iter().position(|candidate| *candidate == node);
    if let Some(position) = path_position {
        let child_toward_root = position.checked_sub(1).map(|index| path[index]);
        let mut outgoing = Vec::new();
        if let Some(parent) = graph.nodes[node].parent {
            outgoing.push(parent);
        }
        outgoing.extend(
            graph.nodes[node]
                .children
                .iter()
                .copied()
                .filter(|child| Some(*child) != child_toward_root),
        );
        outgoing
    } else {
        graph.nodes[node].children.clone()
    }
}

fn double_bond_descriptor_applies_to_node(
    stereo: &DoubleBondStereo,
    descriptor: StereoDescriptor,
    atom: AtomId,
    path: &[AtomId],
) -> bool {
    if !matches!(descriptor, StereoDescriptor::E | StereoDescriptor::Z) {
        return false;
    }
    let other = if stereo.left == atom {
        stereo.right
    } else if stereo.right == atom {
        stereo.left
    } else {
        return false;
    };
    !path.contains(&other)
}

fn rule3_descriptor_priority(descriptor: Option<StereoDescriptor>) -> u8 {
    match descriptor {
        Some(StereoDescriptor::Z) => 2,
        Some(StereoDescriptor::E) => 1,
        _ => 0,
    }
}

fn rule4a_descriptor_priority(descriptor: Option<StereoDescriptor>) -> u8 {
    match descriptor {
        Some(StereoDescriptor::R)
        | Some(StereoDescriptor::S)
        | Some(StereoDescriptor::M)
        | Some(StereoDescriptor::P)
        | Some(StereoDescriptor::SeqTrans)
        | Some(StereoDescriptor::SeqCis) => 2,
        Some(StereoDescriptor::LowerR)
        | Some(StereoDescriptor::LowerS)
        | Some(StereoDescriptor::E)
        | Some(StereoDescriptor::Z) => 1,
        None => 0,
    }
}

fn rule4c_descriptor_priority(descriptor: Option<StereoDescriptor>) -> u8 {
    match descriptor {
        Some(StereoDescriptor::LowerR) => 2,
        Some(StereoDescriptor::LowerS) => 1,
        _ => 0,
    }
}

fn descriptor_ref(descriptor: StereoDescriptor) -> Option<DescriptorRef> {
    match descriptor {
        StereoDescriptor::R | StereoDescriptor::M | StereoDescriptor::SeqCis => {
            Some(DescriptorRef::R)
        }
        StereoDescriptor::S | StereoDescriptor::P | StereoDescriptor::SeqTrans => {
            Some(DescriptorRef::S)
        }
        StereoDescriptor::LowerR
        | StereoDescriptor::LowerS
        | StereoDescriptor::E
        | StereoDescriptor::Z => None,
    }
}

fn double_bond_endpoint_carriers(
    mol: &Molecule,
    endpoint: AtomId,
    other_endpoint: AtomId,
    focus_bond: BondId,
) -> Vec<StereoCarrier> {
    let mut carriers = Vec::new();
    if let Ok(incident) = mol.incident_bonds(endpoint) {
        for (bond_id, bond) in incident {
            if bond_id == focus_bond || bond.order != BondOrder::Single {
                continue;
            }
            let other = bond.other_atom(endpoint);
            if other != other_endpoint {
                carriers.push(StereoCarrier::Atom(other));
            }
        }
    }
    carriers.sort_by_key(carrier_key);
    if mol
        .atom(endpoint)
        .map(|atom| hydrogen_count(atom) == 1)
        .unwrap_or(false)
    {
        carriers.push(StereoCarrier::ImplicitHydrogen);
    }
    carriers
}

fn carrier_key(carrier: &StereoCarrier) -> (u8, u32) {
    match carrier {
        StereoCarrier::Atom(atom) => (0, atom.raw()),
        StereoCarrier::ImplicitHydrogen => (1, u32::MAX),
        StereoCarrier::ImplicitLonePair => (2, u32::MAX),
    }
}

fn invert_double_bond_orientation(orientation: DoubleBondOrientation) -> DoubleBondOrientation {
    match orientation {
        DoubleBondOrientation::Together => DoubleBondOrientation::Opposite,
        DoubleBondOrientation::Opposite => DoubleBondOrientation::Together,
    }
}

fn permutation_is_even(positions: &[usize]) -> bool {
    let mut inversions = 0usize;
    for left in 0..positions.len() {
        for right in (left + 1)..positions.len() {
            if positions[left] > positions[right] {
                inversions += 1;
            }
        }
    }
    inversions % 2 == 0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn node_priority(atomic_number: u8, rule1b: u32, isotope: u16) -> NodePriority {
        NodePriority {
            atomic_number: AtomicNumberFraction::element(atomic_number),
            rule1b,
            rule2_mass: if isotope == 0 {
                Rule2Mass::natural(atomic_number)
            } else {
                Rule2Mass::isotope(isotope)
            },
            descriptor: None,
            rule6_atom: None,
        }
    }

    fn node_priority_with_descriptor(descriptor: StereoDescriptor) -> NodePriority {
        NodePriority {
            descriptor: Some(descriptor),
            ..node_priority(6, 0, 0)
        }
    }

    fn node_priority_with_rule6_atom(atomic_number: u8, atom: AtomId) -> NodePriority {
        NodePriority {
            rule6_atom: Some(atom),
            ..node_priority(atomic_number, 0, 0)
        }
    }

    fn one_node_signature(rule1b: u32, isotope: u16) -> LigandTree {
        LigandTree {
            priority: node_priority(6, rule1b, isotope),
            children: Vec::new(),
        }
    }

    fn signature(root: LigandTree) -> LigandSignature {
        LigandSignature { root }
    }

    #[test]
    fn rule1b_ring_duplicate_priority_is_applied_before_isotope_priority() {
        let ring_duplicate = signature(one_node_signature(u32::MAX, 0));
        let isotope = signature(one_node_signature(0, 13));

        assert_eq!(ring_duplicate.compare(&isotope), Ordering::Greater);
        assert_eq!(isotope.compare(&ring_duplicate), Ordering::Less);
    }

    #[test]
    fn rule2_compares_indicated_isotopes_against_natural_atomic_weight() {
        let natural_hydrogen = node_priority(1, 0, 0);
        let protium = node_priority(1, 0, 1);
        let deuterium = node_priority(1, 0, 2);

        assert_eq!(
            natural_hydrogen.compare_by_rule(&protium, SequenceRule::Rule2, None),
            Ordering::Greater
        );
        assert_eq!(
            deuterium.compare_by_rule(&natural_hydrogen, SequenceRule::Rule2, None),
            Ordering::Greater
        );

        let natural_carbon = node_priority(6, 0, 0);
        let carbon_12 = node_priority(6, 0, 12);
        let carbon_13 = node_priority(6, 0, 13);

        assert_eq!(
            natural_carbon.compare_by_rule(&carbon_12, SequenceRule::Rule2, None),
            Ordering::Greater
        );
        assert_eq!(
            carbon_13.compare_by_rule(&natural_carbon, SequenceRule::Rule2, None),
            Ordering::Greater
        );

        let another_natural_hydrogen = node_priority(1, 0, 0);
        assert_eq!(
            natural_hydrogen.compare_by_rule(&another_natural_hydrogen, SequenceRule::Rule2, None),
            Ordering::Equal
        );
    }

    #[test]
    fn rule1b_prefers_ring_duplicate_whose_reference_is_closer_to_root() {
        let root_reference = signature(one_node_signature(u32::MAX, 0));
        let deeper_reference = signature(one_node_signature(u32::MAX - 2, 0));

        assert_eq!(root_reference.compare(&deeper_reference), Ordering::Greater);
        assert_eq!(deeper_reference.compare(&root_reference), Ordering::Less);
    }

    #[test]
    fn rule4a_prefers_uppercase_sequence_descriptors_over_pseudo_descriptors() {
        let uppercase = signature(LigandTree {
            priority: node_priority_with_descriptor(StereoDescriptor::R),
            children: Vec::new(),
        });
        let sequence = signature(LigandTree {
            priority: node_priority_with_descriptor(StereoDescriptor::SeqCis),
            children: Vec::new(),
        });
        let pseudo = signature(LigandTree {
            priority: node_priority_with_descriptor(StereoDescriptor::LowerR),
            children: Vec::new(),
        });
        let unlabeled = signature(one_node_signature(0, 0));

        assert_eq!(uppercase.compare(&pseudo), Ordering::Greater);
        assert_eq!(sequence.compare(&pseudo), Ordering::Greater);
        assert_eq!(sequence.compare(&uppercase), Ordering::Equal);
        assert_eq!(pseudo.compare(&unlabeled), Ordering::Greater);
        assert_eq!(unlabeled.compare(&uppercase), Ordering::Less);
    }

    #[test]
    fn rule4b_reference_descriptors_use_first_equivalent_descriptor_level() {
        let majority_r = LigandTree {
            priority: node_priority(6, 0, 0),
            children: vec![
                LigandTree {
                    priority: node_priority_with_descriptor(StereoDescriptor::R),
                    children: Vec::new(),
                },
                LigandTree {
                    priority: node_priority_with_descriptor(StereoDescriptor::M),
                    children: Vec::new(),
                },
                LigandTree {
                    priority: node_priority_with_descriptor(StereoDescriptor::S),
                    children: Vec::new(),
                },
            ],
        };
        let tied = LigandTree {
            priority: node_priority(6, 0, 0),
            children: vec![
                LigandTree {
                    priority: node_priority_with_descriptor(StereoDescriptor::R),
                    children: Vec::new(),
                },
                LigandTree {
                    priority: node_priority_with_descriptor(StereoDescriptor::P),
                    children: Vec::new(),
                },
            ],
        };

        assert_eq!(
            majority_r.rule4b_reference_descriptors(),
            vec![DescriptorRef::R]
        );
        assert_eq!(
            tied.rule4b_reference_descriptors(),
            vec![DescriptorRef::R, DescriptorRef::S]
        );
    }

    #[test]
    fn rule4b_fixed_reference_prefers_like_descriptor_families() {
        let r_ligand = LigandTree {
            priority: node_priority_with_descriptor(StereoDescriptor::R),
            children: Vec::new(),
        };
        let s_ligand = LigandTree {
            priority: node_priority_with_descriptor(StereoDescriptor::S),
            children: Vec::new(),
        };

        assert_eq!(
            r_ligand.compare_with_reference(&s_ligand, DescriptorRef::R),
            Ordering::Greater
        );
        assert_eq!(
            r_ligand.compare_with_reference(&s_ligand, DescriptorRef::S),
            Ordering::Less
        );

        let seq_cis = LigandTree {
            priority: node_priority_with_descriptor(StereoDescriptor::SeqCis),
            children: Vec::new(),
        };
        let seq_trans = LigandTree {
            priority: node_priority_with_descriptor(StereoDescriptor::SeqTrans),
            children: Vec::new(),
        };

        assert_eq!(
            seq_cis.compare_with_reference(&seq_trans, DescriptorRef::R),
            Ordering::Greater
        );
        assert_eq!(
            seq_cis.compare_with_reference(&seq_trans, DescriptorRef::S),
            Ordering::Less
        );
    }

    #[test]
    fn rule4c_prefers_lower_r_over_lower_s() {
        let lower_r = signature(LigandTree {
            priority: node_priority_with_descriptor(StereoDescriptor::LowerR),
            children: Vec::new(),
        });
        let lower_s = signature(LigandTree {
            priority: node_priority_with_descriptor(StereoDescriptor::LowerS),
            children: Vec::new(),
        });

        assert_eq!(lower_r.compare(&lower_s), Ordering::Greater);
        assert_eq!(lower_s.compare(&lower_r), Ordering::Less);
    }

    #[test]
    fn rule5_descriptor_pairing_prefers_like_pairs_and_marks_pseudo_asymmetric_ordering() {
        let r_ligand = signature(LigandTree {
            priority: node_priority_with_descriptor(StereoDescriptor::R),
            children: Vec::new(),
        });
        let s_ligand = signature(LigandTree {
            priority: node_priority_with_descriptor(StereoDescriptor::S),
            children: Vec::new(),
        });

        let comparison = r_ligand.compare_with_flags(&s_ligand);

        assert_eq!(comparison.ordering, Ordering::Greater);
        assert!(comparison.pseudo_asymmetric);

        let comparison = s_ligand.compare_with_flags(&r_ligand);

        assert_eq!(comparison.ordering, Ordering::Less);
        assert!(comparison.pseudo_asymmetric);
    }

    #[test]
    fn rule6_prefers_nodes_matching_the_selected_reference_atom() {
        let reference = AtomId::new(7);
        let reference_ligand = signature(LigandTree {
            priority: node_priority_with_rule6_atom(6, reference),
            children: Vec::new(),
        });
        let other_ligand = signature(LigandTree {
            priority: node_priority_with_rule6_atom(6, AtomId::new(8)),
            children: Vec::new(),
        });

        let comparison =
            reference_ligand.compare_with_rule6_reference(&other_ligand, Some(reference));

        assert_eq!(comparison.ordering, Ordering::Greater);
        assert!(!comparison.pseudo_asymmetric);
        assert_eq!(reference_ligand.compare(&other_ligand), Ordering::Equal);
    }

    #[test]
    fn rule6_tetrahedral_retry_resolves_two_equivalent_partitions() {
        let carrier_a = AtomId::new(0);
        let carrier_b = AtomId::new(1);
        let carrier_c = AtomId::new(2);
        let carrier_d = AtomId::new(3);
        let reference_child = LigandTree {
            priority: node_priority_with_rule6_atom(1, carrier_b),
            children: Vec::new(),
        };
        let other_child = LigandTree {
            priority: node_priority_with_rule6_atom(1, AtomId::new(99)),
            children: Vec::new(),
        };
        let signatures = vec![
            (
                StereoCarrier::Atom(carrier_a),
                signature(LigandTree {
                    priority: node_priority_with_rule6_atom(8, carrier_a),
                    children: Vec::new(),
                }),
            ),
            (
                StereoCarrier::Atom(carrier_b),
                signature(LigandTree {
                    priority: node_priority_with_rule6_atom(8, carrier_b),
                    children: Vec::new(),
                }),
            ),
            (
                StereoCarrier::Atom(carrier_c),
                signature(LigandTree {
                    priority: node_priority_with_rule6_atom(6, carrier_c),
                    children: vec![reference_child],
                }),
            ),
            (
                StereoCarrier::Atom(carrier_d),
                signature(LigandTree {
                    priority: node_priority_with_rule6_atom(6, carrier_d),
                    children: vec![other_child],
                }),
            ),
        ];

        let ranked = rank_tetrahedral_signatures_with_rule6(
            &Molecule::new(),
            StereoElementId::new(0),
            AtomId::new(0),
            &signatures,
            TetrahedralOrientation::Clockwise,
            false,
        )
        .expect("Rule 6 should resolve paired partitions");

        assert_eq!(
            ranked.carriers,
            vec![
                StereoCarrier::Atom(carrier_b),
                StereoCarrier::Atom(carrier_a),
                StereoCarrier::Atom(carrier_c),
                StereoCarrier::Atom(carrier_d),
            ]
        );
        assert!(!ranked.pseudo_asymmetric_ordering);
    }

    fn s4_rule6_signatures(
        child_reference_counts: [[usize; 4]; 4],
    ) -> Vec<(StereoCarrier, LigandSignature)> {
        let carriers = [
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(2),
            AtomId::new(3),
        ];
        carriers
            .iter()
            .copied()
            .enumerate()
            .map(|(carrier_index, carrier)| {
                let mut children = Vec::new();
                for (reference_index, reference) in carriers.iter().copied().enumerate() {
                    for _ in 0..child_reference_counts[carrier_index][reference_index] {
                        children.push(LigandTree {
                            priority: node_priority_with_rule6_atom(1, reference),
                            children: Vec::new(),
                        });
                    }
                }
                (
                    StereoCarrier::Atom(carrier),
                    signature(LigandTree {
                        priority: node_priority_with_rule6_atom(6, carrier),
                        children,
                    }),
                )
            })
            .collect()
    }

    #[test]
    fn rule6_s4_retry_accepts_parity_stable_reference_rankings() {
        let signatures =
            s4_rule6_signatures([[0, 2, 0, 2], [2, 0, 2, 0], [1, 0, 2, 1], [0, 1, 1, 2]]);

        let ranked = rank_tetrahedral_signatures_with_rule6(
            &Molecule::new(),
            StereoElementId::new(0),
            AtomId::new(0),
            &signatures,
            TetrahedralOrientation::Clockwise,
            false,
        )
        .expect("Rule 6 should accept parity-stable S4 rankings");

        assert_eq!(
            ranked.carriers,
            vec![
                StereoCarrier::Atom(AtomId::new(0)),
                StereoCarrier::Atom(AtomId::new(1)),
                StereoCarrier::Atom(AtomId::new(2)),
                StereoCarrier::Atom(AtomId::new(3)),
            ]
        );
        assert!(!ranked.pseudo_asymmetric_ordering);
    }

    #[test]
    fn rule6_s4_retry_rejects_parity_unstable_reference_rankings() {
        let element = StereoElementId::new(0);
        let signatures =
            s4_rule6_signatures([[0, 2, 1, 1], [2, 0, 0, 2], [1, 1, 2, 0], [0, 0, 2, 2]]);

        let issue = rank_tetrahedral_signatures_with_rule6(
            &Molecule::new(),
            element,
            AtomId::new(0),
            &signatures,
            TetrahedralOrientation::Clockwise,
            false,
        )
        .expect_err("odd reference permutations must remain unresolved");

        assert_eq!(issue, CipAssignmentIssue::UnresolvedPriority { element });
    }

    #[test]
    fn ligand_tree_compares_highest_priority_branch_before_lower_siblings() {
        let oxygen_to_carbon = LigandTree {
            priority: node_priority(8, 0, 0),
            children: vec![one_node_signature(0, 0)],
        };
        let oxygen_to_hydrogen = LigandTree {
            priority: node_priority(8, 0, 0),
            children: vec![LigandTree {
                priority: node_priority(1, 0, 0),
                children: Vec::new(),
            }],
        };
        let carbon_to_nitrogen = LigandTree {
            priority: node_priority(6, 0, 0),
            children: vec![LigandTree {
                priority: node_priority(7, 0, 0),
                children: Vec::new(),
            }],
        };
        let carbon_to_oxygen = LigandTree {
            priority: node_priority(6, 0, 0),
            children: vec![LigandTree {
                priority: node_priority(8, 0, 0),
                children: Vec::new(),
            }],
        };

        let left = signature(LigandTree {
            priority: node_priority(6, 0, 0),
            children: vec![oxygen_to_carbon, carbon_to_oxygen],
        });
        let right = signature(LigandTree {
            priority: node_priority(6, 0, 0),
            children: vec![oxygen_to_hydrogen, carbon_to_nitrogen],
        });

        assert_eq!(left.compare(&right), Ordering::Greater);
    }

    #[test]
    fn ligand_tree_compares_immediate_sibling_list_before_recursing() {
        let oxygen_to_hydrogen = LigandTree {
            priority: node_priority(8, 0, 0),
            children: vec![LigandTree {
                priority: node_priority(1, 0, 0),
                children: Vec::new(),
            }],
        };
        let oxygen_to_phosphorus = LigandTree {
            priority: node_priority(8, 0, 0),
            children: vec![LigandTree {
                priority: node_priority(15, 0, 0),
                children: Vec::new(),
            }],
        };
        let carbon = one_node_signature(0, 0);
        let hydrogen = LigandTree {
            priority: node_priority(1, 0, 0),
            children: Vec::new(),
        };

        let left = signature(LigandTree {
            priority: node_priority(6, 0, 0),
            children: vec![oxygen_to_hydrogen, carbon],
        });
        let right = signature(LigandTree {
            priority: node_priority(6, 0, 0),
            children: vec![oxygen_to_phosphorus, hydrogen],
        });

        assert_eq!(left.compare(&right), Ordering::Greater);
    }

    #[test]
    fn duplicate_nodes_have_no_isotope_priority() {
        let mut mol = Molecule::new();
        let mut isotope = Atom::new(Element::from_symbol("C").expect("carbon"));
        isotope.isotope = Some(13);
        let atom = mol.add_atom(isotope);

        let normal = LigandNode::Atom {
            atom,
            previous: None,
            path: vec![atom],
            duplicate: None,
            terminal: false,
        };
        let duplicate = LigandNode::Atom {
            atom,
            previous: None,
            path: Vec::new(),
            duplicate: Some(DuplicateNode::Bond {
                atomic_number: None,
            }),
            terminal: true,
        };

        assert_eq!(normal.rule2_mass(&mol), Rule2Mass::isotope(13));
        assert_eq!(duplicate.rule2_mass(&mol), Rule2Mass::ZERO);
    }

    #[test]
    fn rule1a_uses_mancude_fractional_atomic_numbers_for_bond_duplicates() {
        let mut mol = Molecule::new();
        let atoms = (0..6)
            .map(|index| {
                let symbol = if index == 3 { "N" } else { "C" };
                let mut atom = Atom::new(Element::from_symbol(symbol).expect("element"));
                atom.implicit_hydrogens = Some(if index == 3 { 0 } else { 1 });
                mol.add_atom(atom)
            })
            .collect::<Vec<_>>();
        for (left, right, order) in [
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Double),
            (3, 4, BondOrder::Single),
            (4, 5, BondOrder::Double),
            (5, 0, BondOrder::Single),
        ] {
            mol.add_bond(atoms[left], atoms[right], order)
                .expect("ring bond");
        }

        let fractions = cip_atomic_number_fractions(&mol);

        assert_eq!(
            fractions[atoms[2].index()],
            AtomicNumberFraction::new(13, 2)
        );
        assert_eq!(fractions[atoms[3].index()], AtomicNumberFraction::new(6, 1));
        assert_eq!(
            fractions[atoms[4].index()],
            AtomicNumberFraction::new(13, 2)
        );
        assert_eq!(fractions[atoms[0].index()], AtomicNumberFraction::new(6, 1));

        let node = LigandNode::Atom {
            atom: atoms[2],
            previous: Some(atoms[1]),
            path: vec![atoms[1], atoms[2]],
            duplicate: None,
            terminal: false,
        };
        let mut next = Vec::new();
        node.extend(&mol, &fractions, &mut next);

        let normal_nitrogen = next
            .iter()
            .find(|child| {
                matches!(
                    child,
                    LigandNode::Atom {
                        atom,
                        duplicate: None,
                        ..
                    } if *atom == atoms[3]
                )
            })
            .expect("normal nitrogen child");
        let duplicate_nitrogen = next
            .iter()
            .find(|child| {
                matches!(
                    child,
                    LigandNode::Atom {
                        atom,
                        duplicate: Some(DuplicateNode::Bond { .. }),
                        ..
                    } if *atom == atoms[3]
                )
            })
            .expect("duplicate nitrogen child");
        let element = StereoElementId::new(0);
        let descriptor_context = DescriptorContext::new(element, AuxiliaryDescriptorMode::Disabled);
        let build_context = LigandBuildContext {
            mol: &mol,
            element,
            descriptor_context: &descriptor_context,
            options: CipAssignmentOptions::default(),
            atomic_number_fractions: &fractions,
        };

        assert_eq!(
            normal_nitrogen.priority(&build_context).atomic_number,
            AtomicNumberFraction::element(7)
        );
        assert_eq!(
            duplicate_nitrogen.priority(&build_context).atomic_number,
            AtomicNumberFraction::new(13, 2)
        );
    }

    #[test]
    fn higher_order_bond_expansion_creates_terminal_duplicate_nodes() {
        let mut mol = Molecule::new();
        let root = mol.add_atom(Atom::new(Element::from_symbol("C").expect("carbon")));
        let carbon = mol.add_atom(Atom::new(Element::from_symbol("C").expect("carbon")));
        let oxygen = mol.add_atom(Atom::new(Element::from_symbol("O").expect("oxygen")));
        mol.add_bond(root, carbon, BondOrder::Single)
            .expect("root bond");
        mol.add_bond(carbon, oxygen, BondOrder::Double)
            .expect("double bond");

        let node = LigandNode::Atom {
            atom: carbon,
            previous: Some(root),
            path: vec![root, carbon],
            duplicate: None,
            terminal: false,
        };
        let mut next = Vec::new();
        let fractions = cip_atomic_number_fractions(&mol);
        node.extend(&mol, &fractions, &mut next);

        assert_eq!(next.len(), 2);
        assert!(next.contains(&LigandNode::Atom {
            atom: oxygen,
            previous: Some(carbon),
            path: vec![root, carbon, oxygen],
            duplicate: None,
            terminal: false,
        }));
        assert!(next.contains(&LigandNode::Atom {
            atom: oxygen,
            previous: Some(carbon),
            path: Vec::new(),
            duplicate: Some(DuplicateNode::Bond {
                atomic_number: None,
            }),
            terminal: true,
        }));
    }

    #[test]
    fn negative_fractional_atoms_create_duplicate_nodes() {
        let mut mol = Molecule::new();
        let atoms = (0..5)
            .map(|index| {
                let mut atom = Atom::new(Element::from_symbol("C").expect("carbon"));
                atom.implicit_hydrogens = Some(1);
                if index == 2 {
                    atom.formal_charge = -1;
                }
                mol.add_atom(atom)
            })
            .collect::<Vec<_>>();
        for (left, right, order) in [
            (0, 1, BondOrder::Double),
            (1, 2, BondOrder::Single),
            (2, 3, BondOrder::Single),
            (3, 4, BondOrder::Double),
            (4, 0, BondOrder::Single),
        ] {
            mol.add_bond(atoms[left], atoms[right], order)
                .expect("ring bond");
        }

        let mut fractions = vec![AtomicNumberFraction::element(6); mol.atoms.len()];
        fractions[atoms[2].index()] = AtomicNumberFraction::new(13, 2);

        let node = LigandNode::Atom {
            atom: atoms[2],
            previous: Some(atoms[1]),
            path: vec![atoms[1], atoms[2]],
            duplicate: None,
            terminal: false,
        };
        let mut next = Vec::new();
        node.extend(&mol, &fractions, &mut next);

        assert!(next.iter().any(|child| matches!(
            child,
            LigandNode::Atom {
                atom,
                duplicate: Some(DuplicateNode::Bond { .. }),
                terminal: true,
                ..
            } if *atom == atoms[3]
        )));
    }
}
