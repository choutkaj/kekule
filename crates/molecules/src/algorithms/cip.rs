use std::cmp::Ordering;

use crate::algorithms::{
    validate_stereo_with_options, StereoPerceptionIssue, StereoPerceptionOptions,
};
use crate::core::*;

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
        let mut next_pending = Vec::new();
        let mut assigned_this_round = false;
        for (id, element) in pending {
            match assign_cip_element(mol, id, &element, options) {
                CipElementAssignment::Assigned(descriptor) => {
                    set_stereo_descriptor(mol, id, descriptor);
                    report.assigned.push(CipAssignment {
                        element: id,
                        descriptor,
                    });
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
        if !assigned_this_round {
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
    let ranked =
        ranked_tetrahedral_carriers(mol, element, stereo.center, &stereo.carriers, options)?;
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
    Ok(match (descriptor_is_r, ranked.pseudo_asymmetric_ordering) {
        (true, true) => StereoDescriptor::LowerR,
        (false, true) => StereoDescriptor::LowerS,
        (true, false) => StereoDescriptor::R,
        (false, false) => StereoDescriptor::S,
    })
}

fn assign_double_bond_descriptor(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &DoubleBondStereo,
    options: CipAssignmentOptions,
) -> CipResult<StereoDescriptor> {
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
    let signatures = carrier_signatures(mol, element, root, carriers, options)?;
    rank_carrier_signatures(element, &signatures, None)
}

fn ranked_tetrahedral_carriers(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    carriers: &[StereoCarrier],
    options: CipAssignmentOptions,
) -> CipResult<RankedCarriers> {
    let signatures = carrier_signatures(mol, element, root, carriers, options)?;
    match rank_carrier_signatures(element, &signatures, None) {
        Ok(ranked) => Ok(ranked),
        Err(CipAssignmentIssue::UnresolvedPriority { .. }) if carriers.len() == 4 => {
            rank_tetrahedral_signatures_with_rule6(element, &signatures)
        }
        Err(issue) => Err(issue),
    }
}

fn carrier_signatures(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    carriers: &[StereoCarrier],
    options: CipAssignmentOptions,
) -> CipResult<Vec<(StereoCarrier, LigandSignature)>> {
    carriers
        .iter()
        .copied()
        .map(|carrier| {
            carrier_signature(mol, element, carrier, root, options)
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
    element: StereoElementId,
    signatures: &[(StereoCarrier, LigandSignature)],
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
        _ => Err(CipAssignmentIssue::UnresolvedPriority { element }),
    }
}

fn rank_s4_tetrahedral_signatures_with_rule6(
    element: StereoElementId,
    signatures: &[(StereoCarrier, LigandSignature)],
    group: &[usize],
) -> CipResult<RankedCarriers> {
    let Some(first_reference) = group
        .first()
        .and_then(|index| carrier_rule6_atom(signatures[*index].0))
    else {
        return Err(CipAssignmentIssue::UnresolvedPriority { element });
    };
    let Some(second_reference) = group
        .get(1)
        .and_then(|index| carrier_rule6_atom(signatures[*index].0))
    else {
        return Err(CipAssignmentIssue::UnresolvedPriority { element });
    };
    let first = rank_carrier_signatures(element, signatures, Some(first_reference))?;
    let second = rank_carrier_signatures(element, signatures, Some(second_reference))?;
    if carrier_permutation_is_odd(&first.carriers, &second.carriers).unwrap_or(true) {
        Err(CipAssignmentIssue::UnresolvedPriority { element })
    } else {
        Ok(second)
    }
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

    fn compare_without_rule5(&self, other: &Self) -> Ordering {
        for rule in [
            SequenceRule::Rule1a,
            SequenceRule::Rule1b,
            SequenceRule::Rule2,
            SequenceRule::Rule3,
            SequenceRule::Rule4a,
            SequenceRule::Rule4b,
            SequenceRule::Rule4c,
        ] {
            let priority = self.compare_by_sequence_rule(other, rule, None).ordering;
            if priority != Ordering::Equal {
                return priority;
            }
        }
        Ordering::Equal
    }

    fn compare_before_rule4b(&self, other: &Self) -> Ordering {
        for rule in [
            SequenceRule::Rule1a,
            SequenceRule::Rule1b,
            SequenceRule::Rule2,
            SequenceRule::Rule3,
            SequenceRule::Rule4a,
        ] {
            let priority = self.recursive_compare(other, rule, None);
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

    fn children_grouped_before_rule4b(&self) -> Vec<Vec<&LigandTree>> {
        let mut children = self.children.iter().collect::<Vec<_>>();
        children.sort_by(|left, right| right.compare_before_rule4b(left));

        let mut groups: Vec<Vec<&LigandTree>> = Vec::new();
        for child in children {
            if let Some(last) = groups.last_mut() {
                if last[0].compare_before_rule4b(child) == Ordering::Equal {
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
    atomic_number: u8,
    rule1b: u32,
    isotope: u16,
    descriptor: Option<StereoDescriptor>,
    rule6_atom: Option<AtomId>,
}

impl NodePriority {
    fn compare_shallow(&self, other: &Self) -> Ordering {
        self.atomic_number
            .cmp(&other.atomic_number)
            .then_with(|| self.rule1b.cmp(&other.rule1b))
            .then_with(|| self.isotope.cmp(&other.isotope))
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
            SequenceRule::Rule2 => self.isotope.cmp(&other.isotope),
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
            children.sort_by(|left, right| right.compare_without_rule5(left));
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
            let children = node.children_grouped_before_rule4b();
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

fn rule6_priority(atom: Option<AtomId>, reference: Option<AtomId>) -> u8 {
    match (atom, reference) {
        (Some(atom), Some(reference)) if atom == reference => 1,
        _ => 0,
    }
}

fn carrier_signature(
    mol: &Molecule,
    element: StereoElementId,
    carrier: StereoCarrier,
    root: AtomId,
    options: CipAssignmentOptions,
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
    let root = ligand_tree(mol, element, node, options, 0, &mut visited_nodes)?;
    Ok(LigandSignature { root })
}

fn ligand_tree(
    mol: &Molecule,
    element: StereoElementId,
    node: LigandNode,
    options: CipAssignmentOptions,
    depth: usize,
    visited_nodes: &mut usize,
) -> CipResult<LigandTree> {
    *visited_nodes = visited_nodes.saturating_add(1);
    if *visited_nodes > options.max_nodes {
        return Err(CipAssignmentIssue::ResourceLimitExceeded {
            element,
            max_nodes: options.max_nodes,
        });
    }
    let priority = node.priority(mol);
    let mut children = Vec::new();
    if depth < options.max_depth {
        let mut child_nodes = Vec::new();
        node.extend(mol, &mut child_nodes);
        for child in child_nodes {
            children.push(ligand_tree(
                mol,
                element,
                child,
                options,
                depth + 1,
                visited_nodes,
            )?);
        }
        children.sort_by(|left, right| right.priority.compare_shallow(&left.priority));
    }
    Ok(LigandTree { priority, children })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DuplicateNode {
    Bond,
    Ring { reference_depth: usize },
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
    fn priority(&self, mol: &Molecule) -> NodePriority {
        NodePriority {
            atomic_number: self.atomic_number(mol),
            rule1b: self.rule1b_priority(),
            isotope: self.isotope(mol),
            descriptor: self.descriptor(mol),
            rule6_atom: self.rule6_atom(),
        }
    }

    fn atomic_number(&self, mol: &Molecule) -> u8 {
        match self {
            Self::Atom { atom, .. } => mol
                .atom(*atom)
                .map(|atom| atom.element.atomic_number())
                .unwrap_or(0),
            Self::Hydrogen => 1,
            Self::LonePair => 0,
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

    fn isotope(&self, mol: &Molecule) -> u16 {
        match self {
            Self::Atom {
                atom, duplicate, ..
            } => {
                if duplicate.is_some() {
                    0
                } else {
                    mol.atom(*atom)
                        .ok()
                        .and_then(|atom| atom.isotope)
                        .unwrap_or(0)
                }
            }
            Self::Hydrogen | Self::LonePair => 0,
        }
    }

    fn descriptor(&self, mol: &Molecule) -> Option<StereoDescriptor> {
        let Self::Atom {
            atom,
            path,
            duplicate: None,
            ..
        } = self
        else {
            return None;
        };
        atom_descriptor_for_ligand_node(mol, *atom, path)
    }

    fn rule6_atom(&self) -> Option<AtomId> {
        match self {
            Self::Atom { atom, .. } => Some(*atom),
            Self::Hydrogen | Self::LonePair => None,
        }
    }

    fn extend(&self, mol: &Molecule, next: &mut Vec<LigandNode>) {
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
        if let Ok(payload) = mol.atom(*atom) {
            for _ in 0..hydrogen_count(payload) {
                next.push(LigandNode::Hydrogen);
            }
        }
        let Ok(incident) = mol.incident_bonds(*atom) else {
            return;
        };
        for (_, bond) in incident {
            let neighbor = bond.other_atom(*atom);
            let duplicate_count = bond_order_duplicate_count(bond.order);
            if Some(neighbor) == *previous {
                if path.first().copied() != Some(neighbor) {
                    for _ in 0..duplicate_count {
                        next.push(LigandNode::Atom {
                            atom: neighbor,
                            previous: Some(*atom),
                            path: Vec::new(),
                            duplicate: Some(DuplicateNode::Bond),
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
                        duplicate: Some(DuplicateNode::Bond),
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
                        duplicate: Some(DuplicateNode::Bond),
                        terminal: true,
                    });
                }
            }
        }
    }
}

fn hydrogen_count(atom: &Atom) -> u8 {
    atom.explicit_hydrogens
        .saturating_add(atom.implicit_hydrogens.unwrap_or(0))
}

fn bond_order_duplicate_count(order: BondOrder) -> usize {
    match order {
        BondOrder::Double => 1,
        BondOrder::Triple => 2,
        BondOrder::Quadruple => 3,
        _ => 0,
    }
}

fn ring_duplicate_priority(reference_depth: usize) -> u32 {
    let depth = reference_depth.min(u32::MAX as usize) as u32;
    u32::MAX.saturating_sub(depth)
}

fn atom_descriptor_for_ligand_node(
    mol: &Molecule,
    atom: AtomId,
    path: &[AtomId],
) -> Option<StereoDescriptor> {
    mol.stereo_elements().find_map(|(_, element)| {
        let descriptor = element.descriptor?;
        match &element.kind {
            StereoElementKind::Tetrahedral(stereo) if stereo.center == atom => Some(descriptor),
            StereoElementKind::DoubleBond(stereo)
                if double_bond_descriptor_applies_to_node(stereo, descriptor, atom, path) =>
            {
                Some(descriptor)
            }
            _ => None,
        }
    })
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
        | Some(StereoDescriptor::P) => 2,
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
        StereoDescriptor::R | StereoDescriptor::M => Some(DescriptorRef::R),
        StereoDescriptor::S | StereoDescriptor::P => Some(DescriptorRef::S),
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
            atomic_number,
            rule1b,
            isotope,
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
        let pseudo = signature(LigandTree {
            priority: node_priority_with_descriptor(StereoDescriptor::LowerR),
            children: Vec::new(),
        });
        let unlabeled = signature(one_node_signature(0, 0));

        assert_eq!(uppercase.compare(&pseudo), Ordering::Greater);
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

        let ranked = rank_tetrahedral_signatures_with_rule6(StereoElementId::new(0), &signatures)
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
            duplicate: Some(DuplicateNode::Bond),
            terminal: true,
        };

        assert_eq!(normal.isotope(&mol), 13);
        assert_eq!(duplicate.isotope(&mol), 0);
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
        node.extend(&mol, &mut next);

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
            duplicate: Some(DuplicateNode::Bond),
            terminal: true,
        }));
    }
}
