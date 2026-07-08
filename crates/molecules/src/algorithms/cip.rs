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

    let elements = mol
        .stereo_elements()
        .map(|(id, element)| (id, element.clone()))
        .collect::<Vec<_>>();
    for (id, element) in elements {
        if element.specifiedness != StereoSpecifiedness::Specified {
            report.skipped.push(CipSkipped {
                element: id,
                reason: CipSkippedReason::NotSpecified,
            });
            continue;
        }
        let assignment = match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => {
                assign_tetrahedral_descriptor(mol, id, stereo, options)
            }
            StereoElementKind::DoubleBond(stereo) => {
                assign_double_bond_descriptor(mol, id, stereo, options)
            }
            StereoElementKind::Axis(_) => {
                report.skipped.push(CipSkipped {
                    element: id,
                    reason: CipSkippedReason::UnsupportedAxis,
                });
                continue;
            }
        };
        match assignment {
            Ok(descriptor) => {
                set_stereo_descriptor(mol, id, descriptor);
                report.assigned.push(CipAssignment {
                    element: id,
                    descriptor,
                });
            }
            Err(issue) => report.issues.push(issue),
        }
    }
    report
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
    for carrier in ranked {
        let Some(position) = stereo
            .carriers
            .iter()
            .position(|candidate| *candidate == carrier)
        else {
            return Err(CipAssignmentIssue::UnresolvedPriority { element });
        };
        priority_positions.push(position);
    }
    let even = permutation_is_even(&priority_positions);
    let descriptor_is_r = matches!(stereo.orientation, TetrahedralOrientation::Clockwise) != even;
    Ok(if descriptor_is_r {
        StereoDescriptor::R
    } else {
        StereoDescriptor::S
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
        .first()
        .copied()
        .ok_or(CipAssignmentIssue::UnresolvedPriority { element })?;
    let right_top = ranked_carriers(mol, element, stereo.right, &right_carriers, options)?
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
) -> CipResult<Vec<StereoCarrier>> {
    let mut signatures = carriers
        .iter()
        .copied()
        .map(|carrier| {
            carrier_signature(mol, element, carrier, root, options)
                .map(|signature| (carrier, signature))
        })
        .collect::<CipResult<Vec<_>>>()?;
    for left in 0..signatures.len() {
        for right in (left + 1)..signatures.len() {
            if signatures[left].1.compare(&signatures[right].1) == Ordering::Equal {
                return Err(CipAssignmentIssue::UnresolvedPriority { element });
            }
        }
    }
    signatures.sort_by(|left, right| right.1.compare(&left.1));
    Ok(signatures.into_iter().map(|(carrier, _)| carrier).collect())
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LigandSignature {
    root: LigandTree,
}

impl LigandSignature {
    fn compare(&self, other: &Self) -> Ordering {
        self.root.compare(&other.root)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct LigandTree {
    priority: NodePriority,
    children: Vec<LigandTree>,
}

impl LigandTree {
    fn compare(&self, other: &Self) -> Ordering {
        let priority = self.priority.compare(&other.priority);
        if priority != Ordering::Equal {
            return priority;
        }
        let mut queue = vec![(self, other)];
        let mut position = 0usize;
        while position < queue.len() {
            let (left, right) = queue[position];
            position += 1;

            let left_shallow = left.children_sorted_by_priority();
            let right_shallow = right.children_sorted_by_priority();
            let shallow = compare_child_priorities(&left_shallow, &right_shallow);
            if shallow != Ordering::Equal {
                return shallow;
            }

            let left_deep = left.children_sorted_deep();
            let right_deep = right.children_sorted_deep();
            let deep = compare_child_priorities(&left_deep, &right_deep);
            if deep != Ordering::Equal {
                return deep;
            }
            for (left_child, right_child) in left_deep.into_iter().zip(right_deep) {
                queue.push((left_child, right_child));
            }
        }
        Ordering::Equal
    }

    fn children_sorted_by_priority(&self) -> Vec<&LigandTree> {
        let mut children = self.children.iter().collect::<Vec<_>>();
        children.sort_by(|left, right| right.priority.compare(&left.priority));
        children
    }

    fn children_sorted_deep(&self) -> Vec<&LigandTree> {
        let mut children = self.children.iter().collect::<Vec<_>>();
        children.sort_by(|left, right| right.compare(left));
        children
    }
}

fn compare_child_priorities(left: &[&LigandTree], right: &[&LigandTree]) -> Ordering {
    for (left_child, right_child) in left.iter().zip(right) {
        let priority = left_child.priority.compare(&right_child.priority);
        if priority != Ordering::Equal {
            return priority;
        }
    }
    left.len().cmp(&right.len())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct NodePriority {
    atomic_number: u8,
    rule1b: u32,
    isotope: u16,
}

impl NodePriority {
    fn compare(&self, other: &Self) -> Ordering {
        self.atomic_number
            .cmp(&other.atomic_number)
            .then_with(|| self.rule1b.cmp(&other.rule1b))
            .then_with(|| self.isotope.cmp(&other.isotope))
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
        children.sort_by(|left, right| right.priority.compare(&left.priority));
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
}

impl LigandNode {
    fn priority(&self, mol: &Molecule) -> NodePriority {
        NodePriority {
            atomic_number: self.atomic_number(mol),
            rule1b: self.rule1b_priority(),
            isotope: self.isotope(mol),
        }
    }

    fn atomic_number(&self, mol: &Molecule) -> u8 {
        match self {
            Self::Atom { atom, .. } => mol
                .atom(*atom)
                .map(|atom| atom.element.atomic_number())
                .unwrap_or(0),
            Self::Hydrogen => 1,
        }
    }

    fn rule1b_priority(&self) -> u32 {
        match self {
            Self::Atom {
                duplicate: Some(DuplicateNode::Ring { reference_depth }),
                ..
            } => ring_duplicate_priority(*reference_depth),
            Self::Atom { .. } | Self::Hydrogen => 0,
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
            Self::Hydrogen => 0,
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
                for _ in 0..duplicate_count {
                    next.push(LigandNode::Atom {
                        atom: neighbor,
                        previous: Some(*atom),
                        path: Vec::new(),
                        duplicate: Some(DuplicateNode::Bond),
                        terminal: true,
                    });
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

    fn one_node_signature(rule1b: u32, isotope: u16) -> LigandTree {
        LigandTree {
            priority: NodePriority {
                atomic_number: 6,
                rule1b,
                isotope,
            },
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
    fn ligand_tree_compares_highest_priority_branch_before_lower_siblings() {
        let oxygen_to_carbon = LigandTree {
            priority: NodePriority {
                atomic_number: 8,
                rule1b: 0,
                isotope: 0,
            },
            children: vec![one_node_signature(0, 0)],
        };
        let oxygen_to_hydrogen = LigandTree {
            priority: NodePriority {
                atomic_number: 8,
                rule1b: 0,
                isotope: 0,
            },
            children: vec![LigandTree {
                priority: NodePriority {
                    atomic_number: 1,
                    rule1b: 0,
                    isotope: 0,
                },
                children: Vec::new(),
            }],
        };
        let carbon_to_nitrogen = LigandTree {
            priority: NodePriority {
                atomic_number: 6,
                rule1b: 0,
                isotope: 0,
            },
            children: vec![LigandTree {
                priority: NodePriority {
                    atomic_number: 7,
                    rule1b: 0,
                    isotope: 0,
                },
                children: Vec::new(),
            }],
        };
        let carbon_to_oxygen = LigandTree {
            priority: NodePriority {
                atomic_number: 6,
                rule1b: 0,
                isotope: 0,
            },
            children: vec![LigandTree {
                priority: NodePriority {
                    atomic_number: 8,
                    rule1b: 0,
                    isotope: 0,
                },
                children: Vec::new(),
            }],
        };

        let left = signature(LigandTree {
            priority: NodePriority {
                atomic_number: 6,
                rule1b: 0,
                isotope: 0,
            },
            children: vec![oxygen_to_carbon, carbon_to_oxygen],
        });
        let right = signature(LigandTree {
            priority: NodePriority {
                atomic_number: 6,
                rule1b: 0,
                isotope: 0,
            },
            children: vec![oxygen_to_hydrogen, carbon_to_nitrogen],
        });

        assert_eq!(left.compare(&right), Ordering::Greater);
    }

    #[test]
    fn ligand_tree_compares_immediate_sibling_list_before_recursing() {
        let oxygen_to_hydrogen = LigandTree {
            priority: NodePriority {
                atomic_number: 8,
                rule1b: 0,
                isotope: 0,
            },
            children: vec![LigandTree {
                priority: NodePriority {
                    atomic_number: 1,
                    rule1b: 0,
                    isotope: 0,
                },
                children: Vec::new(),
            }],
        };
        let oxygen_to_phosphorus = LigandTree {
            priority: NodePriority {
                atomic_number: 8,
                rule1b: 0,
                isotope: 0,
            },
            children: vec![LigandTree {
                priority: NodePriority {
                    atomic_number: 15,
                    rule1b: 0,
                    isotope: 0,
                },
                children: Vec::new(),
            }],
        };
        let carbon = one_node_signature(0, 0);
        let hydrogen = LigandTree {
            priority: NodePriority {
                atomic_number: 1,
                rule1b: 0,
                isotope: 0,
            },
            children: Vec::new(),
        };

        let left = signature(LigandTree {
            priority: NodePriority {
                atomic_number: 6,
                rule1b: 0,
                isotope: 0,
            },
            children: vec![oxygen_to_hydrogen, carbon],
        });
        let right = signature(LigandTree {
            priority: NodePriority {
                atomic_number: 6,
                rule1b: 0,
                isotope: 0,
            },
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
