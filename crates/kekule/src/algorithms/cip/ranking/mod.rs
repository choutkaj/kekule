use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

use crate::algorithms::{validate_stereo, RingMembership};
use crate::core::*;

use super::super::rings::compute_ring_membership;
use super::super::{atom_hydrogen_count, double_bond_endpoint_carriers};
use super::{
    CipAssignment, CipAssignmentError, CipAssignmentIssue, CipAssignmentOptions,
    CipAssignmentReport, CipResult, CipSkipped, CipSkippedReason,
};

mod assignment;
mod isotope_masses;

use assignment::{
    assign_cip_element, axis_descriptor_from_ranked, axis_endpoint_carriers,
    axis_reference_carriers, double_bond_descriptor_from_ranked, element_is_finally_nonstereogenic,
    rank_carrier_signatures, rank_tetrahedral_signatures_with_rule6, set_stereo_descriptor,
    tetrahedral_descriptor_from_ranked, CipElementAssignment,
};

pub(super) fn assign_cip_descriptors_with_options(
    mol: &mut Molecule,
    options: CipAssignmentOptions,
) -> std::result::Result<CipAssignmentReport, CipAssignmentError> {
    if let Err(error) = validate_stereo(mol) {
        return Err(CipAssignmentError {
            issues: error
                .issues
                .into_iter()
                .map(|issue| CipAssignmentIssue::InvalidStereo { issue })
                .collect(),
        });
    }

    let mut report = CipAssignmentReport::default();
    let mut issues = Vec::new();

    for (id, element) in mol.stereo_elements() {
        match assign_cip_element(mol, id, element, options) {
            CipElementAssignment::Assigned(descriptor) => report.assigned.push(CipAssignment {
                element: id,
                descriptor,
            }),
            CipElementAssignment::Skipped(reason) => report.skipped.push(CipSkipped {
                element: id,
                reason,
            }),
            CipElementAssignment::Unresolved => {
                match element_is_finally_nonstereogenic(mol, id, element, options) {
                    Ok(true) => report.skipped.push(CipSkipped {
                        element: id,
                        reason: CipSkippedReason::NotStereogenic,
                    }),
                    Ok(false) => {
                        issues.push(CipAssignmentIssue::UnresolvedPriority { element: id })
                    }
                    Err(issue) => issues.push(issue),
                }
            }
            CipElementAssignment::Issue(issue) => issues.push(issue),
        }
    }
    if issues.is_empty() {
        mol.replace_stereo_perception(Some(StereoPerception::default()));
        for assignment in &report.assigned {
            set_stereo_descriptor(mol, assignment.element, assignment.descriptor);
        }
        Ok(report)
    } else {
        Err(CipAssignmentError { issues })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RankedCarriers {
    carriers: Vec<StereoCarrier>,
    pseudo_asymmetric_ordering: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LigandSignature {
    root: LigandTree,
    truncated: bool,
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
        if self.truncated || other.truncated {
            return LigandComparison::from_ordering(self.root.recursive_compare(
                &other.root,
                SequenceRule::Rule1a,
                None,
            ));
        }
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
        for rule in SEQUENCE_RULES {
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
            SequenceRule::Rule6 => {
                let ordering = self.recursive_compare(other, rule, rule6_reference);
                // Like Rule 5, a Rule 6 preference is enantiomorphic. Preserve
                // that distinction when combining the ligand comparisons.
                LigandComparison::new(ordering, ordering != Ordering::Equal)
            }
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
        // Rules 4b and 5 may select reference descriptors only at the
        // comparison root. Descendants use the reference already selected.
        for rule in [
            SequenceRule::Rule1a,
            SequenceRule::Rule1b,
            SequenceRule::Rule2,
            SequenceRule::Rule3,
            SequenceRule::Rule4a,
            SequenceRule::Rule4c,
        ] {
            let priority = self.compare_by_sequence_rule(other, rule, None).ordering;
            if priority != Ordering::Equal {
                return priority;
            }
        }
        self.fixed_reference_compare(other, reference)
    }

    fn compare_through_rule4a(&self, other: &Self) -> Ordering {
        for rule in [
            SequenceRule::Rule1a,
            SequenceRule::Rule1b,
            SequenceRule::Rule2,
            SequenceRule::Rule3,
            SequenceRule::Rule4a,
        ] {
            let priority = self.compare_by_sequence_rule(other, rule, None).ordering;
            if priority != Ordering::Equal {
                return priority;
            }
        }
        Ordering::Equal
    }

    fn compare_with_reference(&self, other: &Self, reference: DescriptorRef) -> Ordering {
        self.compare_through_rule4a(other)
            .then_with(|| self.fixed_reference_compare(other, reference))
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
        if left_refs.is_empty() || right_refs.is_empty() {
            return LigandComparison::equal();
        }

        if left_refs.len() == 1 && right_refs.len() == 1 {
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
        children.sort_by(|left, right| {
            for preceding_rule in SEQUENCE_RULES {
                // Unreferenced Rules 4b and 5 do not rank descendants: a
                // nested reference choice would change the comparison root.
                if matches!(preceding_rule, SequenceRule::Rule4b | SequenceRule::Rule5) {
                    if preceding_rule == rule {
                        break;
                    }
                    continue;
                }
                let comparison = if deep {
                    right
                        .compare_by_sequence_rule(left, preceding_rule, rule6_reference)
                        .ordering
                } else {
                    right
                        .priority
                        .compare_by_rule(&left.priority, preceding_rule, rule6_reference)
                };
                if comparison != Ordering::Equal || preceding_rule == rule {
                    return comparison;
                }
            }
            Ordering::Equal
        });
        children
    }

    fn children_sorted_by_reference(&self, reference: DescriptorRef) -> Vec<&LigandTree> {
        let mut children = self.children.iter().collect::<Vec<_>>();
        children.sort_by(|left, right| right.compare_with_reference(left, reference));
        children
    }

    fn children_grouped_through_rule4b(&self) -> Vec<Vec<&LigandTree>> {
        let mut children = self.children.iter().collect::<Vec<_>>();
        children.sort_by(|left, right| right.compare_through_rule4a(left));

        let mut groups: Vec<Vec<&LigandTree>> = Vec::new();
        for child in children {
            if let Some(last) = groups.last_mut() {
                if last[0].compare_through_rule4a(child) == Ordering::Equal {
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

const SEQUENCE_RULES: [SequenceRule; 9] = [
    SequenceRule::Rule1a,
    SequenceRule::Rule1b,
    SequenceRule::Rule2,
    SequenceRule::Rule3,
    SequenceRule::Rule4a,
    SequenceRule::Rule4b,
    SequenceRule::Rule4c,
    SequenceRule::Rule5,
    SequenceRule::Rule6,
];

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

#[derive(Debug)]
struct DescriptorContext {
    skipped: StereoElementId,
    aux_labels: HashMap<AuxDescriptorKey, StereoDescriptor>,
}

impl DescriptorContext {
    fn new(skipped: StereoElementId) -> Self {
        Self {
            skipped,
            aux_labels: HashMap::new(),
        }
    }

    fn skips(&self, element: StereoElementId) -> bool {
        self.skipped == element
    }
}

struct LigandBuildContext<'a> {
    mol: &'a Molecule,
    element: StereoElementId,
    descriptor_context: &'a DescriptorContext,
    options: CipAssignmentOptions,
    atomic_number_fractions: &'a [AtomicNumberFraction],
    atropisomer_mode: bool,
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
    scaled_mass: u64,
    isotope_indicated: bool,
}

impl Rule2Mass {
    const ZERO: Self = Self {
        scaled_mass: 0,
        isotope_indicated: false,
    };

    fn natural(atomic_number: u8) -> Self {
        Self {
            scaled_mass: u64::from(natural_atomic_weight_rank(atomic_number)) * 1_000,
            isotope_indicated: false,
        }
    }

    fn isotope(atomic_number: u8, mass_number: u16) -> Self {
        if mass_number == 0 {
            return Self::natural(atomic_number);
        }
        let (nearest, mass) = isotope_masses::NEAREST_ISOTOPE_MASSES[usize::from(atomic_number)];
        Self {
            scaled_mass: if mass_number == nearest {
                u64::from(mass)
            } else {
                u64::from(mass_number) * 1_000_000
            },
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

fn rule6_priority(atom: Option<AtomId>, reference: Option<AtomId>) -> u8 {
    match (atom, reference) {
        (Some(atom), Some(reference)) if atom == reference => 1,
        _ => 0,
    }
}

fn carrier_signature(
    context: &LigandBuildContext<'_>,
    carrier: StereoCarrier,
    root: AtomId,
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
    let mut truncated = false;
    let root = ligand_tree(context, node, 0, &mut visited_nodes, &mut truncated)?;
    Ok(LigandSignature { root, truncated })
}

fn build_auxiliary_graph(
    mol: &Molecule,
    element: StereoElementId,
    root: AtomId,
    options: CipAssignmentOptions,
    atomic_number_fractions: &[AtomicNumberFraction],
    atropisomer_mode: bool,
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
        atropisomer_mode,
    };
    add_auxiliary_graph_node(&context, &mut graph, root, None, 0, &mut visited_nodes)?;
    Ok(graph)
}

struct AuxiliaryGraphBuildContext<'a> {
    mol: &'a Molecule,
    element: StereoElementId,
    options: CipAssignmentOptions,
    atomic_number_fractions: &'a [AtomicNumberFraction],
    atropisomer_mode: bool,
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
    let mut child_nodes = Vec::new();
    node.extend(
        context.mol,
        context.atomic_number_fractions,
        context.atropisomer_mode,
        &mut child_nodes,
    );
    if depth >= context.options.max_depth.saturating_add(1) && !child_nodes.is_empty() {
        return Err(CipAssignmentIssue::DepthLimitExceeded {
            element: context.element,
            max_depth: context.options.max_depth,
        });
    }
    if depth < context.options.max_depth.saturating_add(1) {
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
    truncated: &mut bool,
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
    let mut child_nodes = Vec::new();
    node.extend(
        context.mol,
        context.atomic_number_fractions,
        context.atropisomer_mode,
        &mut child_nodes,
    );
    if depth >= context.options.max_depth && !child_nodes.is_empty() {
        *truncated = true;
    }
    if depth < context.options.max_depth {
        for child in child_nodes {
            children.push(ligand_tree(
                context,
                child,
                depth + 1,
                visited_nodes,
                truncated,
            )?);
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
                        |isotope| Rule2Mass::isotope(atom.element.atomic_number(), isotope),
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
        atropisomer_mode: bool,
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
        for _ in 0..atom_hydrogen_count(mol, *atom) {
            next.push(LigandNode::Hydrogen);
        }
        let Ok(incident) = mol.incident_bonds(*atom) else {
            return;
        };
        for (_, bond) in incident {
            let neighbor = bond.other_atom(*atom);
            let duplicate_count =
                bond_duplicate_count_for_atom(payload, *atom, bond, atomic_number_fractions);
            let bond_duplicate_atomic_number =
                bond_duplicate_atomic_number(*atom, atomic_number_fractions);
            if Some(neighbor) == *previous {
                if path.first().copied() != Some(neighbor) || atropisomer_mode {
                    for _ in 0..bond_order_duplicate_count(cip_bond_order(bond.order)) {
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
            }
            let duplicate_count = if previous.is_none() && !atropisomer_mode {
                0
            } else {
                duplicate_count
            };
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
    let mut fractions = vec![AtomicNumberFraction::ZERO; mol.graph.atom_slot_count()];
    for (atom_id, atom) in mol.atoms() {
        fractions[atom_id.index()] = AtomicNumberFraction::element(atom.element.atomic_number());
    }

    let ring_membership = mol
        .ring_membership()
        .cloned()
        .unwrap_or_else(|| compute_ring_membership(mol));
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

    relax_mancude_atom_types(mol, &ring_membership, &mut types);
    let parts = mancude_parts(mol, &types, &ring_membership);
    apply_mancude_neighbor_averages(mol, &types, &parts, &mut fractions);
    fractions
}

fn seed_mancude_atom_types(
    mol: &Molecule,
    ring_membership: &RingMembership,
) -> Vec<MancudeAtomType> {
    let mut types = vec![MancudeAtomType::Other; mol.graph.atom_slot_count()];
    for (atom_id, atom) in mol.atoms() {
        let mut bond_types = u32::from(atom_hydrogen_count(mol, atom_id));
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

fn relax_mancude_atom_types(
    mol: &Molecule,
    ring_membership: &RingMembership,
    types: &mut [MancudeAtomType],
) {
    let mut counts = vec![0usize; mol.graph.atom_slot_count()];
    let mut queue = Vec::new();
    for (atom_id, _) in mol.atoms() {
        if types[atom_id.index()] == MancudeAtomType::Other {
            continue;
        }
        for (bond_id, bond) in mol.incident_bonds(atom_id).into_iter().flatten() {
            if ring_membership.bond_in_ring(bond_id)
                && types[bond.other_atom(atom_id).index()] != MancudeAtomType::Other
            {
                counts[atom_id.index()] += 1;
            }
        }
        if counts[atom_id.index()] <= 1 {
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
        for (bond_id, bond) in mol.incident_bonds(atom_id).into_iter().flatten() {
            let neighbor = bond.other_atom(atom_id);
            if !ring_membership.bond_in_ring(bond_id)
                || types[neighbor.index()] == MancudeAtomType::Other
            {
                continue;
            }
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
    let mut parts = vec![0usize; mol.graph.atom_slot_count()];
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
        for (atom_id, _) in mol.atoms() {
            if parts.get(atom_id.index()).copied() != Some(part) {
                continue;
            }
            denominator += 1;
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
        let average = AtomicNumberFraction::new(numerator, denominator);
        for (atom_id, _) in mol.atoms() {
            if parts[atom_id.index()] == part {
                fractions[atom_id.index()] = average;
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

fn bond_duplicate_count_for_atom(
    atom: &Atom,
    atom_id: AtomId,
    bond: &Bond,
    atomic_number_fractions: &[AtomicNumberFraction],
) -> usize {
    let negative_fractional_atom = atom.formal_charge < 0
        && atomic_number_fractions
            .get(atom_id.index())
            .is_some_and(|fraction| fraction.denominator > 1);
    if negative_fractional_atom {
        1
    } else {
        bond_order_duplicate_count(cip_bond_order(bond.order))
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
        BondOrder::Single => 1,
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
    let depth = u32::try_from(reference_depth).unwrap_or(u32::MAX);
    u32::MAX.saturating_sub(depth)
}

fn atom_descriptor_for_ligand_node(
    context: &LigandBuildContext<'_>,
    atom: AtomId,
    path: &[AtomId],
) -> Option<StereoDescriptor> {
    context.mol.stereo_elements().find_map(|(id, _)| {
        if context.descriptor_context.skips(id) {
            return None;
        }
        debug_assert_eq!(path.last(), Some(&atom));
        context
            .descriptor_context
            .aux_labels
            .get(&AuxDescriptorKey {
                element: id,
                path: path.to_vec(),
            })
            .copied()
    })
}

fn collect_auxiliary_occurrences_from_molecule(
    mol: &Molecule,
    context: &DescriptorContext,
    graph: &AuxiliaryGraph,
) -> Vec<AuxOccurrence> {
    let mut occurrences = Vec::new();
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
            let is_focus = match &stereo_element.kind {
                StereoElementKind::Tetrahedral(stereo) => stereo.center == *atom,
                StereoElementKind::DoubleBond(stereo) => {
                    auxiliary_bond_occurrence(graph, node_index, stereo.left, stereo.right)
                }
                StereoElementKind::Axis(stereo) => {
                    let Ok(bond) = mol.bond(stereo.axis) else {
                        continue;
                    };
                    auxiliary_bond_occurrence(graph, node_index, bond.a(), bond.b())
                }
            };
            if is_focus {
                occurrences.push(AuxOccurrence {
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
    occurrences
}

fn auxiliary_bond_occurrence(
    graph: &AuxiliaryGraph,
    node: usize,
    left: AtomId,
    right: AtomId,
) -> bool {
    let Some((left_node, right_node)) = auxiliary_bond_nodes(graph, node, left, right) else {
        return false;
    };
    graph.nodes[node].depth
        == graph.nodes[left_node]
            .depth
            .min(graph.nodes[right_node].depth)
}

fn auxiliary_bond_nodes(
    graph: &AuxiliaryGraph,
    node: usize,
    left: AtomId,
    right: AtomId,
) -> Option<(usize, usize)> {
    let LigandNode::Atom {
        atom,
        duplicate: None,
        ..
    } = graph.nodes[node].node
    else {
        return None;
    };
    let other_atom = if atom == left {
        right
    } else if atom == right {
        left
    } else {
        return None;
    };
    let other_node = graph.nodes[node].parent.into_iter().chain(graph.nodes[node].children.iter().copied()).find(|other| {
        matches!(graph.nodes[*other].node, LigandNode::Atom { atom, duplicate: None, .. } if atom == other_atom)
    })?;
    Some(if atom == left {
        (node, other_node)
    } else {
        (other_node, node)
    })
}

fn precompute_auxiliary_descriptors(
    mol: &Molecule,
    descriptor_context: &mut DescriptorContext,
    graph: &AuxiliaryGraph,
    options: CipAssignmentOptions,
    atomic_number_fractions: &[AtomicNumberFraction],
    atropisomer_mode: bool,
) {
    let mut occurrences =
        collect_auxiliary_occurrences_from_molecule(mol, descriptor_context, graph);
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
            let descriptor = auxiliary_descriptor_for_occurrence(
                mol,
                descriptor_context,
                graph,
                occurrence,
                options,
                atomic_number_fractions,
                atropisomer_mode,
            );
            if let Some(descriptor) = descriptor {
                batch.push((occurrence.key.clone(), descriptor));
            }
        }

        for (key, descriptor) in batch {
            descriptor_context.aux_labels.insert(key, descriptor);
        }
    }
}

fn auxiliary_descriptor_for_occurrence(
    mol: &Molecule,
    descriptor_context: &DescriptorContext,
    graph: &AuxiliaryGraph,
    occurrence: &AuxOccurrence,
    options: CipAssignmentOptions,
    atomic_number_fractions: &[AtomicNumberFraction],
    atropisomer_mode: bool,
) -> Option<StereoDescriptor> {
    let element = mol.stereo_element(occurrence.key.element).ok()?;
    let aux_context = LigandBuildContext {
        mol,
        element: occurrence.key.element,
        descriptor_context,
        options: CipAssignmentOptions {
            max_depth: graph.nodes.len(),
            ..options
        },
        atomic_number_fractions,
        atropisomer_mode,
    };
    let StereoElementKind::Tetrahedral(stereo) = &element.kind else {
        return auxiliary_bond_descriptor(&aux_context, graph, occurrence.node, &element.kind).ok();
    };
    let signatures =
        auxiliary_tetrahedral_signatures(&aux_context, graph, occurrence.node, stereo).ok()?;
    let ranked = match rank_carrier_signatures(occurrence.key.element, &signatures, None) {
        Ok(ranked) => ranked,
        Err(CipAssignmentIssue::UnresolvedPriority { .. }) if stereo.carriers.len() == 4 => {
            rank_tetrahedral_signatures_with_rule6(occurrence.key.element, &signatures).ok()?
        }
        Err(_) => return None,
    };
    tetrahedral_descriptor_from_ranked(occurrence.key.element, stereo, &ranked).ok()
}

fn auxiliary_bond_descriptor(
    context: &LigandBuildContext<'_>,
    graph: &AuxiliaryGraph,
    node: usize,
    kind: &StereoElementKind,
) -> CipResult<StereoDescriptor> {
    let element = context.element;
    let (left, right, left_reference, right_reference, left_carriers, right_carriers) = match kind {
        StereoElementKind::DoubleBond(stereo) => (
            stereo.left,
            stereo.right,
            stereo.left_carrier,
            stereo.right_carrier,
            double_bond_endpoint_carriers(context.mol, stereo.left, stereo.right, stereo.bond),
            double_bond_endpoint_carriers(context.mol, stereo.right, stereo.left, stereo.bond),
        ),
        StereoElementKind::Axis(stereo) => {
            let bond = context
                .mol
                .bond(stereo.axis)
                .map_err(|_| CipAssignmentIssue::UnresolvedPriority { element })?;
            let (left, right) = bond.endpoints();
            let (left_reference, right_reference) =
                axis_reference_carriers(context.mol, element, stereo, left, right)?;
            (
                left,
                right,
                left_reference,
                right_reference,
                axis_endpoint_carriers(context.mol, left, right, stereo.axis),
                axis_endpoint_carriers(context.mol, right, left, stereo.axis),
            )
        }
        StereoElementKind::Tetrahedral(_) => {
            return Err(CipAssignmentIssue::UnresolvedPriority { element })
        }
    };
    let (left_node, right_node) = auxiliary_bond_nodes(graph, node, left, right)
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let rank_endpoint = |root, carriers: Vec<StereoCarrier>| {
        let signatures = carriers
            .into_iter()
            .map(|carrier| {
                auxiliary_carrier_signature(context, graph, root, carrier)
                    .map(|signature| (carrier, signature))
            })
            .collect::<CipResult<Vec<_>>>()?;
        rank_carrier_signatures(element, &signatures, None)
    };
    let left_ranked = rank_endpoint(left_node, left_carriers)?;
    let right_ranked = rank_endpoint(right_node, right_carriers)?;
    match kind {
        StereoElementKind::DoubleBond(stereo) => double_bond_descriptor_from_ranked(
            element,
            stereo
                .orientation
                .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?,
            (left_reference, right_reference),
            &left_ranked,
            &right_ranked,
        ),
        StereoElementKind::Axis(stereo) => axis_descriptor_from_ranked(
            element,
            stereo
                .orientation
                .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?,
            (left_reference, right_reference),
            &left_ranked,
            &right_ranked,
        ),
        StereoElementKind::Tetrahedral(_) => unreachable!(),
    }
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
        StereoCarrier::ImplicitHydrogen => LigandTree {
            priority: LigandNode::Hydrogen.priority(context),
            children: Vec::new(),
        },
        StereoCarrier::ImplicitLonePair => LigandTree {
            priority: LigandNode::LonePair.priority(context),
            children: Vec::new(),
        },
    };
    Ok(LigandSignature {
        root,
        truncated: false,
    })
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
            duplicate: None | Some(DuplicateNode::Ring { .. }),
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
        | Some(StereoDescriptor::LowerM)
        | Some(StereoDescriptor::LowerP)
        | Some(StereoDescriptor::E)
        | Some(StereoDescriptor::Z) => 1,
        None => 0,
    }
}

fn rule4c_descriptor_priority(descriptor: Option<StereoDescriptor>) -> u8 {
    match descriptor {
        Some(StereoDescriptor::LowerR) | Some(StereoDescriptor::LowerM) => 2,
        Some(StereoDescriptor::LowerS) | Some(StereoDescriptor::LowerP) => 1,
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
        | StereoDescriptor::LowerM
        | StereoDescriptor::LowerP
        | StereoDescriptor::E
        | StereoDescriptor::Z => None,
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
    inversions.is_multiple_of(2)
}

#[cfg(test)]
mod tests;
