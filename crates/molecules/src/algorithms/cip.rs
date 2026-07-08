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
    let ranked = ranked_carriers(mol, element, stereo.center, &stereo.carriers, options)?;
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
    let mut signatures = carriers
        .iter()
        .copied()
        .map(|carrier| {
            carrier_signature(mol, element, carrier, root, options)
                .map(|signature| (carrier, signature))
        })
        .collect::<CipResult<Vec<_>>>()?;
    let mut pseudo_asymmetric_pair_count = 0usize;
    for left in 0..signatures.len() {
        for right in (left + 1)..signatures.len() {
            let comparison = signatures[left].1.compare_with_flags(&signatures[right].1);
            if comparison.ordering == Ordering::Equal {
                return Err(CipAssignmentIssue::UnresolvedPriority { element });
            }
            if comparison.pseudo_asymmetric {
                pseudo_asymmetric_pair_count += 1;
            }
        }
    }
    signatures.sort_by(|left, right| right.1.compare(&left.1));
    Ok(RankedCarriers {
        carriers: signatures.into_iter().map(|(carrier, _)| carrier).collect(),
        pseudo_asymmetric_ordering: pseudo_asymmetric_pair_count == 1,
    })
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
        self.root.compare_with_flags(&other.root)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LigandTree {
    priority: NodePriority,
    children: Vec<LigandTree>,
}

impl LigandTree {
    fn compare_with_flags(&self, other: &Self) -> LigandComparison {
        for rule in [
            SequenceRule::Rule1a,
            SequenceRule::Rule1b,
            SequenceRule::Rule2,
            SequenceRule::Rule3,
            SequenceRule::Rule4a,
            SequenceRule::Rule4c,
            SequenceRule::Rule5,
        ] {
            let comparison = self.compare_by_sequence_rule(other, rule);
            if comparison.ordering != Ordering::Equal {
                return comparison;
            }
        }
        LigandComparison::equal()
    }

    fn compare_by_sequence_rule(&self, other: &Self, rule: SequenceRule) -> LigandComparison {
        match rule {
            SequenceRule::Rule5 => self.rule5_pair_comparison(other),
            _ => LigandComparison::from_ordering(self.recursive_compare(other, rule)),
        }
    }

    fn recursive_compare(&self, other: &Self, rule: SequenceRule) -> Ordering {
        let priority = self.priority.compare_by_rule(&other.priority, rule);
        if priority != Ordering::Equal {
            return priority;
        }

        let mut queue = vec![(self, other)];
        let mut position = 0usize;
        while position < queue.len() {
            let (left, right) = queue[position];
            position += 1;

            let left_shallow = left.children_sorted_by_rule(rule, false);
            let right_shallow = right.children_sorted_by_rule(rule, false);
            let shallow = compare_child_priorities(&left_shallow, &right_shallow, rule);
            if shallow != Ordering::Equal {
                return shallow;
            }

            let left_deep = left.children_sorted_by_rule(rule, true);
            let right_deep = right.children_sorted_by_rule(rule, true);
            let deep = compare_child_priorities(&left_deep, &right_deep, rule);
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
            SequenceRule::Rule4c,
        ] {
            let priority = self.recursive_compare(other, rule);
            if priority != Ordering::Equal {
                return priority;
            }
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

    fn children_sorted_by_rule(&self, rule: SequenceRule, deep: bool) -> Vec<&LigandTree> {
        let mut children = self.children.iter().collect::<Vec<_>>();
        if deep {
            children.sort_by(|left, right| right.recursive_compare(left, rule));
        } else {
            children.sort_by(|left, right| {
                right
                    .priority
                    .compare_by_rule(&left.priority, rule)
                    .then_with(|| right.priority.compare_shallow(&left.priority))
            });
        }
        children
    }
}

fn compare_child_priorities(
    left: &[&LigandTree],
    right: &[&LigandTree],
    rule: SequenceRule,
) -> Ordering {
    for (left_child, right_child) in left.iter().zip(right) {
        let priority = left_child
            .priority
            .compare_by_rule(&right_child.priority, rule);
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
    Rule4c,
    Rule5,
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
}

impl NodePriority {
    fn compare_shallow(&self, other: &Self) -> Ordering {
        self.atomic_number
            .cmp(&other.atomic_number)
            .then_with(|| self.rule1b.cmp(&other.rule1b))
            .then_with(|| self.isotope.cmp(&other.isotope))
    }

    fn compare_by_rule(&self, other: &Self, rule: SequenceRule) -> Ordering {
        match rule {
            SequenceRule::Rule1a => self.atomic_number.cmp(&other.atomic_number),
            SequenceRule::Rule1b => self.rule1b.cmp(&other.rule1b),
            SequenceRule::Rule2 => self.isotope.cmp(&other.isotope),
            SequenceRule::Rule3 => rule3_descriptor_priority(self.descriptor)
                .cmp(&rule3_descriptor_priority(other.descriptor)),
            SequenceRule::Rule4a => rule4a_descriptor_priority(self.descriptor)
                .cmp(&rule4a_descriptor_priority(other.descriptor)),
            SequenceRule::Rule4c => rule4c_descriptor_priority(self.descriptor)
                .cmp(&rule4c_descriptor_priority(other.descriptor)),
            SequenceRule::Rule5 => Ordering::Equal,
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
        }
    }

    fn node_priority_with_descriptor(descriptor: StereoDescriptor) -> NodePriority {
        NodePriority {
            descriptor: Some(descriptor),
            ..node_priority(6, 0, 0)
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
