use crate::core::*;
use crate::geometry::Point3;
use std::fmt;

use super::RingMembership;

/// Options for read-only coordinate-derived stereo inference.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CoordinateStereoOptions {
    /// Infer the conservative coordinate-axis subset in addition to
    /// tetrahedral and double-bond stereo.
    pub infer_axes: bool,
}

/// Detached coordinate-derived stereo proposed for represented materialization.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CoordinateStereoResult {
    pub elements: Vec<StereoElement>,
}

/// Stereo elements created by an explicit coordinate-stereo materialization.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CoordinateStereoMaterializationReport {
    pub created_elements: Vec<StereoElementId>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StereoCandidate {
    Tetrahedral {
        center: AtomId,
        carriers: Vec<StereoCarrier>,
    },
    DoubleBond {
        bond: BondId,
        left: AtomId,
        right: AtomId,
        left_carriers: Vec<StereoCarrier>,
        right_carriers: Vec<StereoCarrier>,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StereoValidationIssue {
    MissingStereoAtom {
        element: StereoElementId,
        atom: AtomId,
    },
    MissingStereoBond {
        element: StereoElementId,
        bond: BondId,
    },
    InvalidTetrahedralCarrierCount {
        element: StereoElementId,
        center: AtomId,
        carrier_count: usize,
    },
    DuplicateTetrahedralCarrier {
        element: StereoElementId,
        center: AtomId,
        carrier: StereoCarrier,
    },
    TetrahedralCarrierNotAdjacent {
        element: StereoElementId,
        center: AtomId,
        carrier: StereoCarrier,
    },
    InvalidDoubleBondOrder {
        element: StereoElementId,
        bond: BondId,
        order: BondOrder,
    },
    DoubleBondFocusMismatch {
        element: StereoElementId,
        bond: BondId,
        left: AtomId,
        right: AtomId,
    },
    DoubleBondCarrierIsFocusAtom {
        element: StereoElementId,
        endpoint: AtomId,
        carrier: AtomId,
    },
    DoubleBondCarrierNotAdjacent {
        element: StereoElementId,
        endpoint: AtomId,
        carrier: StereoCarrier,
    },
    UnsupportedDoubleBondCarrier {
        element: StereoElementId,
        endpoint: AtomId,
        carrier: StereoCarrier,
    },
    InvalidAxisCarrierCount {
        element: StereoElementId,
        axis: BondId,
        carrier_count: usize,
    },
    AxisCarrierIsFocusAtom {
        element: StereoElementId,
        axis: BondId,
        carrier: AtomId,
    },
    AxisCarrierNotAdjacent {
        element: StereoElementId,
        axis: BondId,
        carrier: StereoCarrier,
    },
    UnsupportedAxisCarrier {
        element: StereoElementId,
        axis: BondId,
        carrier: StereoCarrier,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StereoValidationError {
    pub issues: Vec<StereoValidationIssue>,
}

impl fmt::Display for StereoValidationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "stereo validation reported {} issue(s)",
            self.issues.len()
        )
    }
}

impl std::error::Error for StereoValidationError {}

/// Failure while inferring or explicitly materializing coordinate stereo.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CoordinateStereoError {
    InvalidStereo(StereoValidationError),
    CouldNotCreateElement(MoleculeError),
}

impl fmt::Display for CoordinateStereoError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidStereo(error) => write!(formatter, "{error}"),
            Self::CouldNotCreateElement(error) => {
                write!(
                    formatter,
                    "could not materialize coordinate stereo: {error}"
                )
            }
        }
    }
}

impl std::error::Error for CoordinateStereoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::InvalidStereo(error) => Some(error),
            Self::CouldNotCreateElement(error) => Some(error),
        }
    }
}

impl From<StereoValidationError> for CoordinateStereoError {
    fn from(error: StereoValidationError) -> Self {
        Self::InvalidStereo(error)
    }
}

pub fn validate_stereo(mol: &Molecule) -> std::result::Result<(), StereoValidationError> {
    let mut issues = Vec::new();
    validate_existing_elements(mol, &mut issues);
    if issues.is_empty() {
        Ok(())
    } else {
        Err(StereoValidationError { issues })
    }
}

pub fn detect_stereo_candidates(mol: &Molecule) -> Vec<StereoCandidate> {
    let mut candidates = tetrahedral_candidates(mol);
    candidates.extend(double_bond_candidates(mol));
    candidates
}

/// Infer coordinate-derived stereo without changing the molecule.
pub fn infer_coordinate_stereo(
    mol: &Molecule,
) -> std::result::Result<CoordinateStereoResult, CoordinateStereoError> {
    infer_coordinate_stereo_with_options(mol, CoordinateStereoOptions::default())
}

/// Infer coordinate-derived stereo with an explicit axis policy.
pub fn infer_coordinate_stereo_with_options(
    mol: &Molecule,
    options: CoordinateStereoOptions,
) -> std::result::Result<CoordinateStereoResult, CoordinateStereoError> {
    validate_stereo(mol)?;
    Ok(CoordinateStereoResult {
        elements: infer_coordinate_stereo_elements(mol, options.infer_axes),
    })
}

/// Materialize inferred coordinate stereo as represented chemistry.
pub fn materialize_coordinate_stereo(
    mol: &mut Molecule,
) -> std::result::Result<CoordinateStereoMaterializationReport, CoordinateStereoError> {
    materialize_coordinate_stereo_with_options(mol, CoordinateStereoOptions::default())
}

/// Materialize inferred coordinate stereo with an explicit axis policy.
pub fn materialize_coordinate_stereo_with_options(
    mol: &mut Molecule,
    options: CoordinateStereoOptions,
) -> std::result::Result<CoordinateStereoMaterializationReport, CoordinateStereoError> {
    let inferred = infer_coordinate_stereo_with_options(mol, options)?;
    let mut staged = mol.clone();
    let mut created_elements = Vec::with_capacity(inferred.elements.len());
    for element in inferred.elements {
        let id = staged
            .add_stereo_element(element)
            .map_err(CoordinateStereoError::CouldNotCreateElement)?;
        created_elements.push(id);
    }
    validate_stereo(&staged)?;
    *mol = staged;
    Ok(CoordinateStereoMaterializationReport { created_elements })
}

fn validate_existing_elements(mol: &Molecule, issues: &mut Vec<StereoValidationIssue>) {
    for (id, element) in mol.stereo_elements() {
        match &element.kind {
            StereoElementKind::Tetrahedral(stereo) => validate_tetrahedral(mol, id, stereo, issues),
            StereoElementKind::DoubleBond(stereo) => validate_double_bond(mol, id, stereo, issues),
            StereoElementKind::Axis(stereo) => validate_axis(mol, id, stereo, issues),
        }
    }
}

fn validate_tetrahedral(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &TetrahedralStereo,
    issues: &mut Vec<StereoValidationIssue>,
) {
    if mol.atom(stereo.center).is_err() {
        issues.push(StereoValidationIssue::MissingStereoAtom {
            element,
            atom: stereo.center,
        });
        return;
    }
    if stereo.carriers.len() != 4 {
        issues.push(StereoValidationIssue::InvalidTetrahedralCarrierCount {
            element,
            center: stereo.center,
            carrier_count: stereo.carriers.len(),
        });
    }
    let mut seen = Vec::<StereoCarrier>::new();
    for carrier in &stereo.carriers {
        if seen.contains(carrier) {
            issues.push(StereoValidationIssue::DuplicateTetrahedralCarrier {
                element,
                center: stereo.center,
                carrier: *carrier,
            });
        } else {
            seen.push(*carrier);
        }
        match carrier {
            StereoCarrier::Atom(atom) => {
                if mol.atom(*atom).is_err() {
                    issues.push(StereoValidationIssue::MissingStereoAtom {
                        element,
                        atom: *atom,
                    });
                } else if mol
                    .bond_between(stereo.center, *atom)
                    .ok()
                    .flatten()
                    .is_none()
                {
                    issues.push(StereoValidationIssue::TetrahedralCarrierNotAdjacent {
                        element,
                        center: stereo.center,
                        carrier: *carrier,
                    });
                }
            }
            StereoCarrier::ImplicitHydrogen | StereoCarrier::ImplicitLonePair => {}
        }
    }
}

fn validate_double_bond(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &DoubleBondStereo,
    issues: &mut Vec<StereoValidationIssue>,
) {
    let Ok(bond) = mol.bond(stereo.bond) else {
        issues.push(StereoValidationIssue::MissingStereoBond {
            element,
            bond: stereo.bond,
        });
        return;
    };
    if bond.order != BondOrder::Double {
        issues.push(StereoValidationIssue::InvalidDoubleBondOrder {
            element,
            bond: stereo.bond,
            order: bond.order,
        });
    }
    if !bond_connects(bond, stereo.left, stereo.right) {
        issues.push(StereoValidationIssue::DoubleBondFocusMismatch {
            element,
            bond: stereo.bond,
            left: stereo.left,
            right: stereo.right,
        });
    }
    validate_double_bond_carrier(
        mol,
        element,
        stereo.left,
        stereo.right,
        stereo.left_carrier,
        issues,
    );
    validate_double_bond_carrier(
        mol,
        element,
        stereo.right,
        stereo.left,
        stereo.right_carrier,
        issues,
    );
}

fn validate_double_bond_carrier(
    mol: &Molecule,
    element: StereoElementId,
    endpoint: AtomId,
    other_endpoint: AtomId,
    carrier: StereoCarrier,
    issues: &mut Vec<StereoValidationIssue>,
) {
    match carrier {
        StereoCarrier::Atom(atom) => {
            if atom == endpoint || atom == other_endpoint {
                issues.push(StereoValidationIssue::DoubleBondCarrierIsFocusAtom {
                    element,
                    endpoint,
                    carrier: atom,
                });
            } else if mol.atom(atom).is_err() {
                issues.push(StereoValidationIssue::MissingStereoAtom { element, atom });
            } else if mol.bond_between(endpoint, atom).ok().flatten().is_none() {
                issues.push(StereoValidationIssue::DoubleBondCarrierNotAdjacent {
                    element,
                    endpoint,
                    carrier,
                });
            }
        }
        StereoCarrier::ImplicitHydrogen => {}
        StereoCarrier::ImplicitLonePair => {
            issues.push(StereoValidationIssue::UnsupportedDoubleBondCarrier {
                element,
                endpoint,
                carrier,
            });
        }
    }
}

fn validate_axis(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &AxisStereo,
    issues: &mut Vec<StereoValidationIssue>,
) {
    let Ok(bond) = mol.bond(stereo.axis) else {
        issues.push(StereoValidationIssue::MissingStereoBond {
            element,
            bond: stereo.axis,
        });
        return;
    };
    if stereo.carriers.len() != 2 {
        issues.push(StereoValidationIssue::InvalidAxisCarrierCount {
            element,
            axis: stereo.axis,
            carrier_count: stereo.carriers.len(),
        });
    }
    let (left, right) = bond.endpoints();
    for carrier in &stereo.carriers {
        validate_axis_carrier(mol, element, stereo.axis, left, right, *carrier, issues);
    }
}

fn validate_axis_carrier(
    mol: &Molecule,
    element: StereoElementId,
    axis: BondId,
    left: AtomId,
    right: AtomId,
    carrier: StereoCarrier,
    issues: &mut Vec<StereoValidationIssue>,
) {
    match carrier {
        StereoCarrier::Atom(atom) => {
            if atom == left || atom == right {
                issues.push(StereoValidationIssue::AxisCarrierIsFocusAtom {
                    element,
                    axis,
                    carrier: atom,
                });
            } else if mol.atom(atom).is_err() {
                issues.push(StereoValidationIssue::MissingStereoAtom { element, atom });
            } else {
                let adjacent_left = mol.bond_between(left, atom).ok().flatten().is_some();
                let adjacent_right = mol.bond_between(right, atom).ok().flatten().is_some();
                if adjacent_left == adjacent_right {
                    issues.push(StereoValidationIssue::AxisCarrierNotAdjacent {
                        element,
                        axis,
                        carrier,
                    });
                }
            }
        }
        StereoCarrier::ImplicitHydrogen | StereoCarrier::ImplicitLonePair => {
            issues.push(StereoValidationIssue::UnsupportedAxisCarrier {
                element,
                axis,
                carrier,
            });
        }
    }
}

fn tetrahedral_candidates(mol: &Molecule) -> Vec<StereoCandidate> {
    let mut candidates = Vec::new();
    for (center, atom) in mol.atoms() {
        if atom.element.symbol() == "H" {
            continue;
        }
        let Ok(incident) = mol.incident_bonds(center) else {
            continue;
        };
        let mut atom_carriers = Vec::new();
        let mut single_bonded = true;
        for (_, bond) in incident {
            single_bonded &= bond.order == BondOrder::Single;
            atom_carriers.push(StereoCarrier::Atom(bond.other_atom(center)));
        }
        atom_carriers.sort_by_key(|carrier| carrier.canonical_order_key());
        let hydrogens = atom_hydrogen_count(mol, center);
        if single_bonded && hydrogens <= 1 && atom_carriers.len() + usize::from(hydrogens) == 4 {
            if hydrogens == 1 {
                atom_carriers.push(StereoCarrier::ImplicitHydrogen);
            }
            candidates.push(StereoCandidate::Tetrahedral {
                center,
                carriers: atom_carriers,
            });
        }
    }
    candidates
}

fn double_bond_candidates(mol: &Molecule) -> Vec<StereoCandidate> {
    let mut candidates = Vec::new();
    for (bond_id, bond) in mol.bonds() {
        if double_bond_stereo_is_unsupported(mol, bond_id, bond) {
            continue;
        }
        let left = bond.a();
        let right = bond.b();
        let left_carriers = double_bond_endpoint_carriers(mol, left, right, bond_id);
        let right_carriers = double_bond_endpoint_carriers(mol, right, left, bond_id);
        if !left_carriers.is_empty() && !right_carriers.is_empty() {
            candidates.push(StereoCandidate::DoubleBond {
                bond: bond_id,
                left,
                right,
                left_carriers,
                right_carriers,
            });
        }
    }
    candidates
}

fn double_bond_stereo_is_unsupported(mol: &Molecule, bond_id: BondId, bond: &Bond) -> bool {
    bond.order != BondOrder::Double
        || mol.bond_is_aromatic(bond_id).ok().flatten() == Some(true)
        || double_bond_between_aromatic_atoms(mol, bond)
        || super::rings::bond_in_ring_smaller_than(mol, bond_id, 8)
        || (double_bond_is_in_ring(mol, bond_id) && double_bond_has_noncarbon_endpoint(mol, bond))
}

pub(crate) fn double_bond_between_aromatic_atoms(mol: &Molecule, bond: &Bond) -> bool {
    mol.atom_is_aromatic(bond.a()).ok().flatten() == Some(true)
        && mol.atom_is_aromatic(bond.b()).ok().flatten() == Some(true)
}

pub(crate) fn double_bond_is_in_ring(mol: &Molecule, bond: BondId) -> bool {
    mol.ring_membership()
        .map(|membership| membership.bond_in_ring(bond))
        .unwrap_or(false)
}

pub(crate) fn double_bond_has_noncarbon_endpoint(mol: &Molecule, bond: &Bond) -> bool {
    [bond.a(), bond.b()].into_iter().any(|atom_id| {
        mol.atom(atom_id)
            .map(|atom| atom.element.symbol() != "C")
            .unwrap_or(true)
    })
}

pub(crate) fn double_bond_endpoint_carriers(
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
    carriers.sort_by_key(|carrier| carrier.canonical_order_key());
    if atom_hydrogen_count(mol, endpoint) == 1 {
        carriers.push(StereoCarrier::ImplicitHydrogen);
    }
    carriers
}

pub(crate) fn atom_axis_carriers(
    mol: &Molecule,
    endpoint: AtomId,
    axis: BondId,
) -> Option<Vec<AtomId>> {
    let mut carriers = Vec::new();
    for (bond_id, bond) in mol.incident_bonds(endpoint).ok()? {
        if bond_id != axis {
            carriers.push(bond.other_atom(endpoint));
        }
    }
    carriers.sort();
    Some(carriers)
}

fn atom_is_atropisomeric_sp2_endpoint(
    mol: &Molecule,
    ring_membership: &RingMembership,
    atom_id: AtomId,
) -> bool {
    if mol.atom(atom_id).is_err() {
        return false;
    }
    let incident = mol
        .incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    let total_degree = incident
        .len()
        .saturating_add(usize::from(atom_hydrogen_count(mol, atom_id)));
    if !(2..=3).contains(&total_degree) {
        return false;
    }
    ring_membership.atom_in_ring(atom_id)
        || mol.atom_is_aromatic(atom_id).ok().flatten() == Some(true)
        || incident.iter().any(|(bond_id, bond)| {
            mol.bond_is_aromatic(*bond_id).ok().flatten() == Some(true)
                || bond.order == BondOrder::Double
        })
}

fn infer_coordinate_stereo_elements(mol: &Molecule, infer_axes: bool) -> Vec<StereoElement> {
    let Some((_, conformer)) = mol.first_conformer() else {
        return Vec::new();
    };
    let mut assigned = Vec::new();
    assigned.extend(infer_coordinate_tetrahedral(mol, conformer));
    assigned.extend(infer_coordinate_double_bonds(mol, conformer));
    if infer_axes {
        assigned.extend(infer_coordinate_axes(mol, conformer));
    }
    assigned
}

fn infer_coordinate_tetrahedral(mol: &Molecule, conformer: &Conformer) -> Vec<StereoElement> {
    let mut assigned = Vec::new();
    for candidate in tetrahedral_candidates(mol) {
        let StereoCandidate::Tetrahedral { center, carriers } = candidate else {
            continue;
        };
        if has_tetrahedral_stereo(mol, center) {
            continue;
        }
        let atom_carriers = carriers
            .iter()
            .map(|carrier| match carrier {
                StereoCarrier::Atom(atom) => Some(*atom),
                StereoCarrier::ImplicitHydrogen | StereoCarrier::ImplicitLonePair => None,
            })
            .collect::<Option<Vec<_>>>();
        let Some(atom_carriers) = atom_carriers else {
            continue;
        };
        let Some(points) = tetrahedral_points(conformer, center, &atom_carriers) else {
            continue;
        };
        let Some(orientation) = tetrahedral_orientation_from_points(points) else {
            continue;
        };
        assigned.push(StereoElement::new(StereoElementKind::Tetrahedral(
            TetrahedralStereo {
                center,
                carriers,
                orientation: Some(orientation),
            },
        )));
    }
    assigned
}

fn infer_coordinate_double_bonds(mol: &Molecule, conformer: &Conformer) -> Vec<StereoElement> {
    let mut assigned = Vec::new();
    for candidate in double_bond_candidates(mol) {
        let StereoCandidate::DoubleBond {
            bond,
            left,
            right,
            left_carriers,
            right_carriers,
        } = candidate
        else {
            continue;
        };
        if has_double_bond_stereo(mol, bond) {
            continue;
        }
        let Some(left_carrier) = only_atom_carrier(&left_carriers) else {
            continue;
        };
        let Some(right_carrier) = only_atom_carrier(&right_carriers) else {
            continue;
        };
        let Some(points) = double_bond_points(conformer, left, right, left_carrier, right_carrier)
        else {
            continue;
        };
        let Some(orientation) = double_bond_orientation_from_points(points) else {
            continue;
        };
        assigned.push(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond,
                left,
                right,
                left_carrier: StereoCarrier::Atom(left_carrier),
                right_carrier: StereoCarrier::Atom(right_carrier),
                orientation: Some(orientation),
            },
        )));
    }
    assigned
}

fn infer_coordinate_axes(mol: &Molecule, conformer: &Conformer) -> Vec<StereoElement> {
    let ring_membership = mol
        .ring_membership()
        .cloned()
        .unwrap_or_else(|| super::rings::compute_ring_membership(mol));
    let mut assigned = Vec::new();
    for (axis, bond) in mol.bonds() {
        if bond.order != BondOrder::Single || has_axis_stereo(mol, axis) {
            continue;
        }
        let (left, right) = bond.endpoints();
        if !atom_is_atropisomeric_sp2_endpoint(mol, &ring_membership, left)
            || !atom_is_atropisomeric_sp2_endpoint(mol, &ring_membership, right)
        {
            continue;
        }
        let Some(left_carriers) = atom_axis_carriers(mol, left, axis) else {
            continue;
        };
        let Some(right_carriers) = atom_axis_carriers(mol, right, axis) else {
            continue;
        };
        if left_carriers.len() != 2 || right_carriers.len() != 2 {
            continue;
        }
        let left_reference = left_carriers[0];
        let right_reference = right_carriers[0];
        let Some(orientation) =
            axis_orientation_from_3d_coordinates(conformer, bond, left_reference, right_reference)
        else {
            continue;
        };
        assigned.push(StereoElement::new(StereoElementKind::Axis(AxisStereo {
            axis,
            carriers: vec![
                StereoCarrier::Atom(left_reference),
                StereoCarrier::Atom(right_reference),
            ],
            orientation: Some(orientation),
        })));
    }
    assigned
}

fn only_atom_carrier(carriers: &[StereoCarrier]) -> Option<AtomId> {
    let mut atoms = carriers.iter().filter_map(|carrier| match carrier {
        StereoCarrier::Atom(atom) => Some(*atom),
        StereoCarrier::ImplicitHydrogen | StereoCarrier::ImplicitLonePair => None,
    });
    let atom = atoms.next()?;
    atoms.next().is_none().then_some(atom)
}

pub(crate) fn tetrahedral_points(
    conformer: &Conformer,
    center: AtomId,
    carriers: &[AtomId],
) -> Option<[Point3; 5]> {
    (carriers.len() == 4).then_some(())?;
    Some([
        conformer.position_value(center)?,
        conformer.position_value(carriers[0])?,
        conformer.position_value(carriers[1])?,
        conformer.position_value(carriers[2])?,
        conformer.position_value(carriers[3])?,
    ])
}

fn double_bond_points(
    conformer: &Conformer,
    left: AtomId,
    right: AtomId,
    left_carrier: AtomId,
    right_carrier: AtomId,
) -> Option<[Point3; 4]> {
    Some([
        conformer.position_value(left)?,
        conformer.position_value(right)?,
        conformer.position_value(left_carrier)?,
        conformer.position_value(right_carrier)?,
    ])
}

pub(crate) fn tetrahedral_orientation_from_points(
    points: [Point3; 5],
) -> Option<TetrahedralOrientation> {
    let a = points[1] - points[4];
    let b = points[2] - points[4];
    let c = points[3] - points[4];
    let volume = a.cross(b).dot(c);
    if volume.abs() <= COORDINATE_EPSILON {
        return None;
    }
    Some(if volume > 0.0 {
        TetrahedralOrientation::Clockwise
    } else {
        TetrahedralOrientation::CounterClockwise
    })
}

fn double_bond_orientation_from_points(points: [Point3; 4]) -> Option<DoubleBondOrientation> {
    let axis = points[1] - points[0];
    let left_vector = points[2] - points[0];
    let right_vector = points[3] - points[1];
    let sidedness = axis.cross(left_vector).dot(axis.cross(right_vector));
    if sidedness.abs() <= COORDINATE_EPSILON {
        return None;
    }
    Some(if sidedness > 0.0 {
        DoubleBondOrientation::Together
    } else {
        DoubleBondOrientation::Opposite
    })
}

fn axis_orientation_from_3d_coordinates(
    conformer: &Conformer,
    axis_bond: &Bond,
    left_reference: AtomId,
    right_reference: AtomId,
) -> Option<AxisOrientation> {
    let (left, right) = axis_bond.endpoints();
    let left_point = conformer.position_value(left)?;
    let right_point = conformer.position_value(right)?;
    let left_reference_point = conformer.position_value(left_reference)?;
    let right_reference_point = conformer.position_value(right_reference)?;
    let points = [
        left_point,
        right_point,
        left_reference_point,
        right_reference_point,
    ];
    if coordinates_are_planar(&points) {
        return None;
    }
    let axis = right_point - left_point;
    let left_vector = left_reference_point - left_point;
    let right_vector = right_reference_point - right_point;
    let handedness = axis.dot(left_vector.cross(right_vector));
    if handedness.abs() <= COORDINATE_EPSILON {
        return None;
    }
    Some(if handedness > 0.0 {
        AxisOrientation::Clockwise
    } else {
        AxisOrientation::CounterClockwise
    })
}

pub(crate) fn coordinates_are_planar(points: &[Point3]) -> bool {
    let Some(first) = points.first() else {
        return false;
    };
    points
        .iter()
        .all(|point| (point.z - first.z).abs() <= COORDINATE_EPSILON)
}

const COORDINATE_EPSILON: f64 = 1.0e-8;

pub(crate) fn atom_hydrogen_count(mol: &Molecule, atom: AtomId) -> u8 {
    let Ok(payload) = mol.atom(atom) else {
        return 0;
    };
    payload
        .hydrogens
        .explicit_count()
        .saturating_add(mol.implicit_hydrogens(atom).ok().flatten().unwrap_or(0))
}

fn bond_connects(bond: &Bond, a: AtomId, b: AtomId) -> bool {
    (bond.a() == a && bond.b() == b) || (bond.a() == b && bond.b() == a)
}
