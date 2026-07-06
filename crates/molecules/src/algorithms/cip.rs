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
            max_nodes: 4096,
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
    atomic_spheres: Vec<Vec<u8>>,
    isotope_spheres: Vec<Vec<u16>>,
}

impl LigandSignature {
    fn compare(&self, other: &Self) -> Ordering {
        let len = self
            .atomic_spheres
            .len()
            .max(other.atomic_spheres.len())
            .max(self.isotope_spheres.len())
            .max(other.isotope_spheres.len());
        for index in 0..len {
            let atomic = compare_sphere(
                self.atomic_spheres.get(index).map(Vec::as_slice),
                other.atomic_spheres.get(index).map(Vec::as_slice),
            );
            if atomic != Ordering::Equal {
                return atomic;
            }
            let isotope = compare_sphere(
                self.isotope_spheres.get(index).map(Vec::as_slice),
                other.isotope_spheres.get(index).map(Vec::as_slice),
            );
            if isotope != Ordering::Equal {
                return isotope;
            }
        }
        Ordering::Equal
    }
}

fn compare_sphere<T: Ord>(left: Option<&[T]>, right: Option<&[T]>) -> Ordering {
    match (left, right) {
        (Some(left), Some(right)) => left.cmp(right),
        (Some(left), None) => {
            if left.is_empty() {
                Ordering::Equal
            } else {
                Ordering::Greater
            }
        }
        (None, Some(right)) => {
            if right.is_empty() {
                Ordering::Equal
            } else {
                Ordering::Less
            }
        }
        (None, None) => Ordering::Equal,
    }
}

fn carrier_signature(
    mol: &Molecule,
    element: StereoElementId,
    carrier: StereoCarrier,
    root: AtomId,
    options: CipAssignmentOptions,
) -> CipResult<LigandSignature> {
    let mut nodes = vec![match carrier {
        StereoCarrier::Atom(atom) => LigandNode::Atom {
            atom,
            previous: Some(root),
            path: vec![root, atom],
            terminal: false,
        },
        StereoCarrier::ImplicitHydrogen => LigandNode::Hydrogen,
    }];
    let mut atomic_spheres = Vec::new();
    let mut isotope_spheres = Vec::new();
    let mut visited_nodes = nodes.len();
    for depth in 0..=options.max_depth {
        let mut atomic = nodes
            .iter()
            .map(|node| node.atomic_number(mol))
            .collect::<Vec<_>>();
        atomic.sort_unstable_by(|left, right| right.cmp(left));
        atomic_spheres.push(atomic);

        let mut isotopes = nodes
            .iter()
            .map(|node| node.isotope(mol))
            .collect::<Vec<_>>();
        isotopes.sort_unstable_by(|left, right| right.cmp(left));
        isotope_spheres.push(isotopes);

        if depth == options.max_depth {
            break;
        }

        let mut next = Vec::new();
        for node in &nodes {
            node.extend(mol, &mut next);
        }
        if next.is_empty() {
            break;
        }
        visited_nodes = visited_nodes.saturating_add(next.len());
        if visited_nodes > options.max_nodes {
            return Err(CipAssignmentIssue::ResourceLimitExceeded {
                element,
                max_nodes: options.max_nodes,
            });
        }
        nodes = next;
    }
    Ok(LigandSignature {
        atomic_spheres,
        isotope_spheres,
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum LigandNode {
    Atom {
        atom: AtomId,
        previous: Option<AtomId>,
        path: Vec<AtomId>,
        terminal: bool,
    },
    Hydrogen,
}

impl LigandNode {
    fn atomic_number(&self, mol: &Molecule) -> u8 {
        match self {
            Self::Atom { atom, .. } => mol
                .atom(*atom)
                .map(|atom| atom.element.atomic_number())
                .unwrap_or(0),
            Self::Hydrogen => 1,
        }
    }

    fn isotope(&self, mol: &Molecule) -> u16 {
        match self {
            Self::Atom { atom, .. } => mol
                .atom(*atom)
                .ok()
                .and_then(|atom| atom.isotope)
                .unwrap_or(0),
            Self::Hydrogen => 0,
        }
    }

    fn extend(&self, mol: &Molecule, next: &mut Vec<LigandNode>) {
        let Self::Atom {
            atom,
            previous,
            path,
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
            if Some(neighbor) == *previous {
                continue;
            }
            let count = bond_order_duplicate_count(bond.order);
            for _ in 0..count {
                let closes_cycle = path.contains(&neighbor);
                let mut next_path = path.clone();
                if !closes_cycle {
                    next_path.push(neighbor);
                }
                next.push(LigandNode::Atom {
                    atom: neighbor,
                    previous: Some(*atom),
                    path: next_path,
                    terminal: closes_cycle,
                });
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
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        _ => 1,
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
