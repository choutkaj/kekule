use crate::chemistry::AtomPositionSource;
use crate::geometry::{Point3, Vector3};
use crate::molecule::*;
use crate::structure::Positions;
use std::fmt;

mod axis_geometry;
pub(crate) use axis_geometry::atom_is_atropisomeric_sp2_endpoint;
mod eligibility;
mod nitrogen;
mod perception;
mod symmetry;
use eligibility::{
    has_repeated_terminal_ligands, tetrahedral_carriers, unclassified_tetrahedral_geometry,
};
use perception::prepared;
pub use perception::{
    cleanup_stereo, detect_stereo_candidates, StereoCleanupReport, StereoSearchError,
    StereoSearchOptions,
};

struct MoleculePositions<'a> {
    atom_to_dense: Vec<Option<usize>>,
    positions: &'a Positions,
}

impl<'a> MoleculePositions<'a> {
    fn new(molecule: &Molecule, positions: &'a Positions) -> Self {
        debug_assert_eq!(molecule.atom_count(), positions.len());
        let capacity = molecule
            .atom_ids()
            .map(AtomId::index)
            .max()
            .map_or(0, |index| index + 1);
        let mut atom_to_dense = vec![None; capacity];
        for (dense_index, atom) in molecule.atom_ids().enumerate() {
            atom_to_dense[atom.index()] = Some(dense_index);
        }
        Self {
            atom_to_dense,
            positions,
        }
    }
}

impl AtomPositionSource for MoleculePositions<'_> {
    fn position_value(&self, atom: AtomId) -> Option<Point3> {
        let index = self.atom_to_dense.get(atom.index()).copied().flatten()?;
        self.positions
            .position_at(index)
            .ok()
            .map(|point| point.into_value())
    }
}

/// Options for read-only coordinate-derived stereo inference.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CoordinateStereoOptions {
    /// Infer the conservative coordinate-axis subset in addition to
    /// tetrahedral and double-bond stereo.
    pub infer_axes: bool,
    /// Bounds for chemical candidate perception before coordinate inference.
    pub perception: StereoSearchOptions,
}

/// Coordinate-derived stereo for the foci the molecule leaves unasserted, and
/// the asserted stereo the coordinates contradict.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CoordinateStereoResult {
    /// New elements for stereocenters, double bonds, and (when enabled) axes
    /// without a stereo element, ready to materialize.
    pub elements: Vec<StereoElement>,
    /// Existing specified elements whose configuration is the opposite of the
    /// one their carriers have in the coordinates, in element order. Unknown
    /// configurations and degenerate or flat geometry are never contradictions.
    pub disagreements: Vec<StereoElementId>,
}

/// Stereo elements created by an explicit coordinate-stereo materialization.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct CoordinateStereoMaterializationReport {
    pub created_elements: Vec<StereoElementId>,
    /// Existing elements the coordinates contradict; materialization leaves
    /// them unchanged. See [`CoordinateStereoResult::disagreements`].
    pub disagreements: Vec<StereoElementId>,
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
    DuplicateStereoFocus {
        element: StereoElementId,
        previous: StereoElementId,
    },
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
    UnrepresentedTetrahedralNeighbor {
        element: StereoElementId,
        center: AtomId,
        neighbor: AtomId,
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
    DoubleBondEndpointOvercoordinated {
        element: StereoElementId,
        endpoint: AtomId,
        substituent_count: usize,
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
    AxisCarrierEndpointMismatch {
        element: StereoElementId,
        axis: BondId,
    },
    AxisEndpointOvercoordinated {
        element: StereoElementId,
        endpoint: AtomId,
        substituent_count: usize,
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
    Perception(StereoSearchError),
    PositionCountMismatch { expected: usize, actual: usize },
    InvalidStereo(StereoValidationError),
    CouldNotCreateElement(MoleculeError),
}

impl fmt::Display for CoordinateStereoError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Perception(error) => error.fmt(formatter),
            Self::PositionCountMismatch { expected, actual } => write!(
                formatter,
                "molecule requires {expected} positions, but received {actual}"
            ),
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
            Self::Perception(error) => Some(error),
            Self::PositionCountMismatch { .. } => None,
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

/// Infers coordinate-derived stereo without changing the molecule.
///
/// Only foci without a stereo element are inferred: an asserted center,
/// double bond, or axis keeps its element, including an explicitly unknown
/// configuration, so the result can be materialized as is. Asserted
/// configurations the coordinates contradict are reported in
/// [`CoordinateStereoResult::disagreements`] rather than overridden.
///
/// `options` selects whether axes are inferred. Asserted axes are compared
/// with the coordinates whether or not axis inference is enabled.
///
/// `positions` must be dense and follow the molecule's atom iteration order.
pub fn infer_coordinate_stereo(
    mol: &Molecule,
    positions: &Positions,
    options: CoordinateStereoOptions,
) -> std::result::Result<CoordinateStereoResult, CoordinateStereoError> {
    if positions.len() != mol.atom_count() {
        return Err(CoordinateStereoError::PositionCountMismatch {
            expected: mol.atom_count(),
            actual: positions.len(),
        });
    }
    validate_stereo(mol)?;
    let perceived = perception::prepared(mol).map_err(CoordinateStereoError::Perception)?;
    let mol = perceived.as_ref();
    let candidates = detect_stereo_candidates(mol, options.perception)
        .map_err(CoordinateStereoError::Perception)?;
    let source = MoleculePositions::new(mol, positions);
    let disagreements = mol
        .stereo_elements()
        .filter(|(_, element)| coordinates_contradict(mol, &source, &element.kind))
        .map(|(id, _)| id)
        .collect();
    Ok(CoordinateStereoResult {
        elements: infer_coordinate_stereo_elements(mol, &source, options.infer_axes, &candidates),
        disagreements,
    })
}

/// Materializes inferred coordinate stereo as represented chemistry.
///
/// Adds elements only for unasserted foci; asserted stereo, including
/// configurations the coordinates contradict, is left unchanged and the
/// contradictions are reported. See [`infer_coordinate_stereo`].
///
/// `options` selects whether axes are inferred, as for
/// [`infer_coordinate_stereo`].
///
/// `positions` must be dense and follow the edited molecule's atom iteration
/// order.
pub fn materialize_coordinate_stereo(
    editor: &mut MoleculeEditor,
    positions: &Positions,
    options: CoordinateStereoOptions,
) -> std::result::Result<CoordinateStereoMaterializationReport, CoordinateStereoError> {
    let inferred = infer_coordinate_stereo(editor.working(), positions, options)?;
    let mut staged = editor.clone();
    let mut created_elements = Vec::with_capacity(inferred.elements.len());
    for element in inferred.elements {
        let id = staged
            .add_stereo_element(element)
            .map_err(CoordinateStereoError::CouldNotCreateElement)?;
        created_elements.push(id);
    }
    validate_stereo(staged.working())?;
    *editor = staged;
    Ok(CoordinateStereoMaterializationReport {
        created_elements,
        disagreements: inferred.disagreements,
    })
}

fn validate_existing_elements(mol: &Molecule, issues: &mut Vec<StereoValidationIssue>) {
    let mut foci = std::collections::BTreeMap::new();
    for (id, element) in mol.stereo_elements() {
        if let Some(previous) = foci.insert(element.kind.focus(), id) {
            issues.push(StereoValidationIssue::DuplicateStereoFocus {
                element: id,
                previous,
            });
        }
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
    for neighbor in mol.neighbors(stereo.center).into_iter().flatten() {
        if !stereo.carriers.contains(&StereoCarrier::Atom(neighbor)) {
            issues.push(StereoValidationIssue::UnrepresentedTetrahedralNeighbor {
                element,
                center: stereo.center,
                neighbor,
            });
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
    for (endpoint, other) in [(stereo.left, stereo.right), (stereo.right, stereo.left)] {
        let substituent_count = mol
            .neighbors(endpoint)
            .into_iter()
            .flatten()
            .filter(|neighbor| *neighbor != other)
            .count();
        if substituent_count > 2 {
            issues.push(StereoValidationIssue::DoubleBondEndpointOvercoordinated {
                element,
                endpoint,
                substituent_count,
            });
        }
    }
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
    for (endpoint, other) in [(left, right), (right, left)] {
        let substituent_count = mol
            .neighbors(endpoint)
            .into_iter()
            .flatten()
            .filter(|neighbor| *neighbor != other)
            .count();
        if substituent_count > 2 {
            issues.push(StereoValidationIssue::AxisEndpointOvercoordinated {
                element,
                endpoint,
                substituent_count,
            });
        }
    }
    for carrier in &stereo.carriers {
        validate_axis_carrier(mol, element, stereo.axis, left, right, *carrier, issues);
    }
    if stereo.carriers.len() == 2
        && stereo
            .carriers
            .iter()
            .all(|carrier| matches!(carrier, StereoCarrier::Atom(_)))
    {
        let endpoint_counts = [left, right].map(|endpoint| {
            stereo
                .carriers
                .iter()
                .filter(|carrier| {
                    matches!(carrier, StereoCarrier::Atom(atom)
                    if mol.bond_between(endpoint, *atom).ok().flatten().is_some())
                })
                .count()
        });
        if endpoint_counts != [1, 1] {
            issues.push(StereoValidationIssue::AxisCarrierEndpointMismatch {
                element,
                axis: stereo.axis,
            });
        }
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
    for center in mol.atom_ids() {
        let Some(atom_carriers) = tetrahedral_carriers(mol, center) else {
            continue;
        };
        if has_repeated_terminal_ligands(mol, center, &atom_carriers) {
            continue;
        }
        candidates.push(StereoCandidate::Tetrahedral {
            center,
            carriers: atom_carriers,
        });
    }
    candidates
}

fn double_bond_candidates(mol: &Molecule) -> Vec<StereoCandidate> {
    let mut candidates = Vec::new();
    for (bond_id, bond) in mol.bonds() {
        if !double_bond_geometry_is_supported(mol, bond_id, bond)
            || mol.bond_is_aromatic(bond_id).ok().flatten() == Some(true)
        {
            continue;
        }
        let left = bond.a();
        let right = bond.b();
        let left_carriers = double_bond_endpoint_carriers(mol, left, right, bond_id);
        let right_carriers = double_bond_endpoint_carriers(mol, right, left, bond_id);
        if (1..=2).contains(&left_carriers.len())
            && (1..=2).contains(&right_carriers.len())
            && !has_repeated_terminal_ligands(mol, left, &left_carriers)
            && !has_repeated_terminal_ligands(mol, right, &right_carriers)
        {
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

/// Double bonds that are stereo candidates before graph-symmetry analysis:
/// supported geometry, not aromatic, and one or two carriers at each endpoint,
/// neither endpoint bearing two identical terminal ligands. Unlike the full
/// candidate set, this does not depend on the molecule's stereo elements.
pub(crate) fn locally_eligible_double_bonds(
    mol: &Molecule,
) -> std::result::Result<std::collections::BTreeSet<BondId>, StereoSearchError> {
    let perceived = prepared(mol)?;
    Ok(double_bond_candidates(&perceived)
        .into_iter()
        .filter_map(|candidate| match candidate {
            StereoCandidate::DoubleBond { bond, .. } => Some(bond),
            StereoCandidate::Tetrahedral { .. } => None,
        })
        .collect())
}

/// Shared represented-graph geometry bounds for source decoding and perception.
/// Ring size is graph-theoretic and independent of the selected ring basis.
/// A heteroatom endpoint does not by itself exclude cis/trans stereochemistry.
pub(crate) fn double_bond_geometry_is_supported(
    mol: &Molecule,
    bond_id: BondId,
    bond: &Bond,
) -> bool {
    bond.order == BondOrder::Double
        && [bond.a(), bond.b()].into_iter().all(|endpoint| {
            mol.neighbors(endpoint)
                .map(|neighbors| neighbors.count() <= 3)
                .unwrap_or(false)
        })
        && !super::rings::bond_in_ring_smaller_than(mol, bond_id, 8)
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

fn infer_coordinate_stereo_elements(
    mol: &Molecule,
    coordinates: &dyn AtomPositionSource,
    infer_axes: bool,
    candidates: &[StereoCandidate],
) -> Vec<StereoElement> {
    let mut assigned = Vec::new();
    assigned.extend(infer_coordinate_tetrahedral(mol, coordinates, candidates));
    assigned.extend(infer_coordinate_double_bonds(mol, coordinates, candidates));
    if infer_axes {
        assigned.extend(infer_coordinate_axes(mol, coordinates));
    }
    assigned
}

fn infer_coordinate_tetrahedral(
    mol: &Molecule,
    coordinates: &dyn AtomPositionSource,
    candidates: &[StereoCandidate],
) -> Vec<StereoElement> {
    let mut assigned = Vec::new();
    for candidate in candidates.iter().cloned() {
        let StereoCandidate::Tetrahedral { center, carriers } = candidate else {
            continue;
        };
        if has_tetrahedral_stereo(mol, center) {
            continue;
        }
        let Some(orientation) = tetrahedral_orientation_at(coordinates, center, &carriers) else {
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

fn infer_coordinate_double_bonds(
    mol: &Molecule,
    coordinates: &dyn AtomPositionSource,
    candidates: &[StereoCandidate],
) -> Vec<StereoElement> {
    let mut assigned = Vec::new();
    for candidate in candidates.iter().cloned() {
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
        let Some(left_carrier) = first_atom_carrier(&left_carriers) else {
            continue;
        };
        let Some(right_carrier) = first_atom_carrier(&right_carriers) else {
            continue;
        };
        let Some(points) =
            double_bond_points(coordinates, left, right, left_carrier, right_carrier)
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

fn infer_coordinate_axes(
    mol: &Molecule,
    coordinates: &dyn AtomPositionSource,
) -> Vec<StereoElement> {
    let mut assigned = Vec::new();
    for (axis, bond) in mol.bonds() {
        if bond.order != BondOrder::Single || has_axis_stereo(mol, axis) {
            continue;
        }
        let (left, right) = bond.endpoints();
        if !atom_is_atropisomeric_sp2_endpoint(mol, left, atom_hydrogen_count(mol, left))
            || !atom_is_atropisomeric_sp2_endpoint(mol, right, atom_hydrogen_count(mol, right))
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
        let Some(orientation) = axis_orientation_from_3d_coordinates(
            coordinates,
            bond,
            left_reference,
            right_reference,
        ) else {
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

/// Whether the coordinates give an asserted configuration's own carriers the
/// opposite configuration. Unknown configurations and degenerate geometry are
/// not contradictions.
pub(crate) fn coordinates_contradict(
    mol: &Molecule,
    coordinates: &dyn AtomPositionSource,
    kind: &StereoElementKind,
) -> bool {
    match kind {
        StereoElementKind::Tetrahedral(stereo) => stereo
            .orientation
            .zip(tetrahedral_orientation_at(
                coordinates,
                stereo.center,
                &stereo.carriers,
            ))
            .is_some_and(|(asserted, observed)| asserted != observed),
        StereoElementKind::DoubleBond(stereo) => stereo
            .orientation
            .zip(double_bond_orientation_at(mol, coordinates, stereo))
            .is_some_and(|(asserted, observed)| asserted != observed),
        StereoElementKind::Axis(stereo) => stereo
            .orientation
            .zip(axis_orientation_at(mol, coordinates, stereo))
            .is_some_and(|(asserted, observed)| asserted != observed),
    }
}

/// The configuration of four tetrahedral carriers in `coordinates`.
pub(crate) fn tetrahedral_orientation_at(
    coordinates: &dyn AtomPositionSource,
    center: AtomId,
    carriers: &[StereoCarrier],
) -> Option<TetrahedralOrientation> {
    let center_point = coordinates.position_value(center)?;
    // For three explicit ligands, using the center as the fourth point
    // has the same determinant sign as a virtual ligand opposite their
    // vector sum. No hydrogen coordinate needs to be guessed or stored.
    let point = |carrier: &StereoCarrier| match carrier {
        StereoCarrier::Atom(atom) => coordinates.position_value(*atom),
        StereoCarrier::ImplicitHydrogen | StereoCarrier::ImplicitLonePair => Some(center_point),
    };
    let [first, second, third, fourth] = carriers else {
        return None;
    };
    tetrahedral_orientation_from_points([
        center_point,
        point(first)?,
        point(second)?,
        point(third)?,
        point(fourth)?,
    ])
}

/// The configuration of a double bond's carriers in `coordinates`. An
/// implicit-hydrogen carrier is located through the endpoint's one explicit
/// substituent, which lies on the opposite side.
fn double_bond_orientation_at(
    mol: &Molecule,
    coordinates: &dyn AtomPositionSource,
    stereo: &DoubleBondStereo,
) -> Option<DoubleBondOrientation> {
    let reference = |endpoint: AtomId, other: AtomId, carrier: StereoCarrier| match carrier {
        StereoCarrier::Atom(atom) => Some((atom, false)),
        StereoCarrier::ImplicitHydrogen => {
            let mut substituents = mol
                .neighbors(endpoint)
                .ok()?
                .filter(|neighbor| *neighbor != other);
            let substituent = substituents.next()?;
            substituents.next().is_none().then_some((substituent, true))
        }
        StereoCarrier::ImplicitLonePair => None,
    };
    let (left, left_inverted) = reference(stereo.left, stereo.right, stereo.left_carrier)?;
    let (right, right_inverted) = reference(stereo.right, stereo.left, stereo.right_carrier)?;
    let orientation = double_bond_orientation_from_points(double_bond_points(
        coordinates,
        stereo.left,
        stereo.right,
        left,
        right,
    )?)?;
    Some(if left_inverted != right_inverted {
        orientation.inverted()
    } else {
        orientation
    })
}

/// The configuration of an axis's two atom carriers in `coordinates`. The
/// handedness is unchanged when both endpoint roles swap, so it does not depend
/// on which endpoint the first carrier is attached to.
fn axis_orientation_at(
    mol: &Molecule,
    coordinates: &dyn AtomPositionSource,
    stereo: &AxisStereo,
) -> Option<AxisOrientation> {
    let [StereoCarrier::Atom(first), StereoCarrier::Atom(second)] = stereo.carriers[..] else {
        return None;
    };
    let bond = mol.bond(stereo.axis).ok()?;
    let (start, _) = bond.endpoints();
    let first_on_start = mol.bond_between(start, first).ok()?.is_some();
    let (left_reference, right_reference) = if first_on_start {
        (first, second)
    } else {
        (second, first)
    };
    axis_orientation_from_3d_coordinates(coordinates, bond, left_reference, right_reference)
}

fn first_atom_carrier(carriers: &[StereoCarrier]) -> Option<AtomId> {
    carriers.iter().find_map(|carrier| match carrier {
        StereoCarrier::Atom(atom) => Some(*atom),
        StereoCarrier::ImplicitHydrogen | StereoCarrier::ImplicitLonePair => None,
    })
}

fn double_bond_points(
    coordinates: &dyn AtomPositionSource,
    left: AtomId,
    right: AtomId,
    left_carrier: AtomId,
    right_carrier: AtomId,
) -> Option<[Point3; 4]> {
    Some([
        coordinates.position_value(left)?,
        coordinates.position_value(right)?,
        coordinates.position_value(left_carrier)?,
        coordinates.position_value(right_carrier)?,
    ])
}

pub(crate) fn tetrahedral_orientation_from_points(
    points: [Point3; 5],
) -> Option<TetrahedralOrientation> {
    let a = normalized_direction(points[1] - points[4])?;
    let b = normalized_direction(points[2] - points[4])?;
    let c = normalized_direction(points[3] - points[4])?;
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

pub(crate) fn double_bond_orientation_from_points(
    points: [Point3; 4],
) -> Option<DoubleBondOrientation> {
    let axis = normalized_direction(points[1] - points[0])?;
    let left_vector = normalized_direction(points[2] - points[0])?;
    let right_vector = normalized_direction(points[3] - points[1])?;
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
    coordinates: &dyn AtomPositionSource,
    axis_bond: &Bond,
    left_reference: AtomId,
    right_reference: AtomId,
) -> Option<AxisOrientation> {
    let (left, right) = axis_bond.endpoints();
    let left_point = coordinates.position_value(left)?;
    let right_point = coordinates.position_value(right)?;
    let left_reference_point = coordinates.position_value(left_reference)?;
    let right_reference_point = coordinates.position_value(right_reference)?;
    let axis = normalized_direction(right_point - left_point)?;
    let left_vector = normalized_direction(left_reference_point - left_point)?;
    let right_vector = normalized_direction(right_reference_point - right_point)?;
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

// Test angular degeneracy instead of comparing a dimensional volume with a
// fixed threshold. The orientation must not depend on coordinate scale.
fn normalized_direction(vector: Vector3) -> Option<Vector3> {
    let norm = vector.norm();
    (norm.is_finite() && norm > 0.0).then(|| vector / norm)
}

/// Whether a nitrogen's lone pair is a tetrahedral stereo carrier under the
/// shared eligibility model (a three-membered ring or a bridgehead), so source
/// annotations decode the same nitrogen centers that coordinate inference
/// finds. Without installed perception it is read from a detached perceived
/// copy; a graph outside the valence model has no such nitrogen.
pub(crate) fn nitrogen_lone_pair_is_stereogenic(mol: &Molecule, center: AtomId) -> bool {
    mol.atom(center)
        .is_ok_and(|atom| atom.element.symbol() == "N")
        && prepared(mol).is_ok_and(|mol| nitrogen::lone_pair_is_stereogenic(&mol, center))
}

pub(crate) fn atom_hydrogen_count(mol: &Molecule, atom: AtomId) -> u8 {
    let count = mol.implicit_hydrogens(atom).ok().flatten().unwrap_or(0);
    // Stereo eligibility only needs small ligand counts. Preserve the existing
    // saturation for oversized expert-supplied perception assignments.
    u8::try_from(count).unwrap_or(u8::MAX)
}

fn bond_connects(bond: &Bond, a: AtomId, b: AtomId) -> bool {
    (bond.a() == a && bond.b() == b) || (bond.a() == b && bond.b() == a)
}
