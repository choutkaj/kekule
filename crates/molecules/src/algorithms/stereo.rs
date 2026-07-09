use crate::core::*;

use super::RingMembership;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StereoPerceptionOptions {
    pub validate_existing: bool,
    pub detect_candidates: bool,
    pub assemble_source_marks: bool,
    pub assign_coordinates: bool,
}

impl Default for StereoPerceptionOptions {
    fn default() -> Self {
        Self {
            validate_existing: true,
            detect_candidates: true,
            assemble_source_marks: true,
            assign_coordinates: true,
        }
    }
}

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct StereoPerceptionReport {
    pub candidates: Vec<StereoCandidate>,
    pub issues: Vec<StereoPerceptionIssue>,
    pub assembled_elements: Vec<StereoElement>,
    pub created_elements: Vec<StereoElementId>,
}

impl StereoPerceptionReport {
    pub fn is_ok(&self) -> bool {
        self.issues.is_empty()
    }
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
pub enum StereoPerceptionIssue {
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
    TetrahedralHydrogenCarrierUnavailable {
        element: StereoElementId,
        center: AtomId,
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
    DoubleBondHydrogenCarrierUnavailable {
        element: StereoElementId,
        endpoint: AtomId,
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
    AmbiguousTetrahedralWedgeMarks {
        center: AtomId,
        mark_count: usize,
    },
    UnassembledTetrahedralBondMark {
        bond: BondId,
        kind: StereoBondMarkKind,
    },
    AmbiguousDirectionalBondMarks {
        double_bond: BondId,
        endpoint: AtomId,
        mark_count: usize,
    },
    UnpairedDirectionalBondMark {
        bond: BondId,
    },
    UnsupportedSourceBondMark {
        bond: BondId,
        kind: StereoBondMarkKind,
    },
    CouldNotCreateElement {
        message: String,
    },
}

pub fn validate_stereo(mol: &Molecule) -> StereoPerceptionReport {
    validate_stereo_with_options(mol, StereoPerceptionOptions::default())
}

pub fn validate_stereo_with_options(
    mol: &Molecule,
    options: StereoPerceptionOptions,
) -> StereoPerceptionReport {
    stereo_report(mol, options)
}

pub fn perceive_stereo(mol: &mut Molecule) -> StereoPerceptionReport {
    perceive_stereo_with_options(mol, StereoPerceptionOptions::default())
}

pub fn perceive_stereo_with_options(
    mol: &mut Molecule,
    options: StereoPerceptionOptions,
) -> StereoPerceptionReport {
    let mut report = stereo_report(mol, options);
    for element in report.assembled_elements.clone() {
        match mol.add_stereo_element(element) {
            Ok(id) => report.created_elements.push(id),
            Err(error) => report
                .issues
                .push(StereoPerceptionIssue::CouldNotCreateElement {
                    message: error.to_string(),
                }),
        }
    }
    mol.perception.stereo = ComputedState::Fresh;
    report
}

fn stereo_report(mol: &Molecule, options: StereoPerceptionOptions) -> StereoPerceptionReport {
    let mut report = StereoPerceptionReport::default();
    if options.validate_existing {
        validate_existing_elements(mol, &mut report.issues);
    }
    if options.detect_candidates {
        report.candidates.extend(tetrahedral_candidates(mol));
        report.candidates.extend(double_bond_candidates(mol));
    }
    if options.assemble_source_marks {
        let mut used_marks = Vec::<BondId>::new();
        let axis_elements = assemble_atropisomeric_axes(mol, &mut used_marks);
        report
            .assembled_elements
            .extend(assemble_tetrahedral_wedges(
                mol,
                &mut report.issues,
                &mut used_marks,
            ));
        report.assembled_elements.extend(axis_elements);
        report
            .assembled_elements
            .extend(assemble_directional_double_bonds(
                mol,
                &mut report.issues,
                &mut used_marks,
            ));
        report_unassembled_source_marks(mol, &used_marks, &mut report.issues);
    }
    if options.assign_coordinates {
        report
            .assembled_elements
            .extend(assign_coordinate_stereo(mol, &report.assembled_elements));
    }
    report
}

fn validate_existing_elements(mol: &Molecule, issues: &mut Vec<StereoPerceptionIssue>) {
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
    issues: &mut Vec<StereoPerceptionIssue>,
) {
    if mol.atom(stereo.center).is_err() {
        issues.push(StereoPerceptionIssue::MissingStereoAtom {
            element,
            atom: stereo.center,
        });
        return;
    }
    if stereo.carriers.len() != 4 {
        issues.push(StereoPerceptionIssue::InvalidTetrahedralCarrierCount {
            element,
            center: stereo.center,
            carrier_count: stereo.carriers.len(),
        });
    }
    let mut seen = Vec::<StereoCarrier>::new();
    for carrier in &stereo.carriers {
        if seen.contains(carrier) {
            issues.push(StereoPerceptionIssue::DuplicateTetrahedralCarrier {
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
                    issues.push(StereoPerceptionIssue::MissingStereoAtom {
                        element,
                        atom: *atom,
                    });
                } else if mol
                    .bond_between(stereo.center, *atom)
                    .ok()
                    .flatten()
                    .is_none()
                {
                    issues.push(StereoPerceptionIssue::TetrahedralCarrierNotAdjacent {
                        element,
                        center: stereo.center,
                        carrier: *carrier,
                    });
                }
            }
            StereoCarrier::ImplicitHydrogen => {
                if hydrogen_count(mol, stereo.center) == 0 {
                    issues.push(
                        StereoPerceptionIssue::TetrahedralHydrogenCarrierUnavailable {
                            element,
                            center: stereo.center,
                        },
                    );
                }
            }
            StereoCarrier::ImplicitLonePair => {
                if !implicit_lone_pair_available(mol, stereo.center) {
                    issues.push(
                        StereoPerceptionIssue::TetrahedralHydrogenCarrierUnavailable {
                            element,
                            center: stereo.center,
                        },
                    );
                }
            }
        }
    }
}

fn validate_double_bond(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &DoubleBondStereo,
    issues: &mut Vec<StereoPerceptionIssue>,
) {
    let Ok(bond) = mol.bond(stereo.bond) else {
        issues.push(StereoPerceptionIssue::MissingStereoBond {
            element,
            bond: stereo.bond,
        });
        return;
    };
    if bond.order != BondOrder::Double {
        issues.push(StereoPerceptionIssue::InvalidDoubleBondOrder {
            element,
            bond: stereo.bond,
            order: bond.order,
        });
    }
    if !bond_connects(bond, stereo.left, stereo.right) {
        issues.push(StereoPerceptionIssue::DoubleBondFocusMismatch {
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
    issues: &mut Vec<StereoPerceptionIssue>,
) {
    match carrier {
        StereoCarrier::Atom(atom) => {
            if atom == endpoint || atom == other_endpoint {
                issues.push(StereoPerceptionIssue::DoubleBondCarrierIsFocusAtom {
                    element,
                    endpoint,
                    carrier: atom,
                });
            } else if mol.atom(atom).is_err() {
                issues.push(StereoPerceptionIssue::MissingStereoAtom { element, atom });
            } else if mol.bond_between(endpoint, atom).ok().flatten().is_none() {
                issues.push(StereoPerceptionIssue::DoubleBondCarrierNotAdjacent {
                    element,
                    endpoint,
                    carrier,
                });
            }
        }
        StereoCarrier::ImplicitHydrogen => {
            if hydrogen_count(mol, endpoint) == 0 {
                issues.push(
                    StereoPerceptionIssue::DoubleBondHydrogenCarrierUnavailable {
                        element,
                        endpoint,
                    },
                );
            }
        }
        StereoCarrier::ImplicitLonePair => {
            issues.push(
                StereoPerceptionIssue::DoubleBondHydrogenCarrierUnavailable { element, endpoint },
            );
        }
    }
}

fn validate_axis(
    mol: &Molecule,
    element: StereoElementId,
    stereo: &AxisStereo,
    issues: &mut Vec<StereoPerceptionIssue>,
) {
    let Ok(bond) = mol.bond(stereo.axis) else {
        issues.push(StereoPerceptionIssue::MissingStereoBond {
            element,
            bond: stereo.axis,
        });
        return;
    };
    if stereo.carriers.len() != 2 {
        issues.push(StereoPerceptionIssue::InvalidAxisCarrierCount {
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
    issues: &mut Vec<StereoPerceptionIssue>,
) {
    match carrier {
        StereoCarrier::Atom(atom) => {
            if atom == left || atom == right {
                issues.push(StereoPerceptionIssue::AxisCarrierIsFocusAtom {
                    element,
                    axis,
                    carrier: atom,
                });
            } else if mol.atom(atom).is_err() {
                issues.push(StereoPerceptionIssue::MissingStereoAtom { element, atom });
            } else {
                let adjacent_left = mol.bond_between(left, atom).ok().flatten().is_some();
                let adjacent_right = mol.bond_between(right, atom).ok().flatten().is_some();
                if adjacent_left == adjacent_right {
                    issues.push(StereoPerceptionIssue::AxisCarrierNotAdjacent {
                        element,
                        axis,
                        carrier,
                    });
                }
            }
        }
        StereoCarrier::ImplicitHydrogen | StereoCarrier::ImplicitLonePair => {
            issues.push(StereoPerceptionIssue::AxisCarrierNotAdjacent {
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
        atom_carriers.sort_by_key(carrier_key);
        let hydrogens = hydrogen_count(mol, center);
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

fn tetrahedral_carriers(mol: &Molecule, center: AtomId) -> Option<Vec<StereoCarrier>> {
    let atom = mol.atom(center).ok()?;
    if atom.element.symbol() == "H" {
        return None;
    }
    let mut carriers = Vec::new();
    for (_, bond) in mol.incident_bonds(center).ok()? {
        if bond.order != BondOrder::Single {
            return None;
        }
        carriers.push(StereoCarrier::Atom(bond.other_atom(center)));
    }
    carriers.sort_by_key(carrier_key);
    let hydrogens = hydrogen_count(mol, center);
    if hydrogens > 1 || carriers.len() + usize::from(hydrogens) != 4 {
        return None;
    }
    if hydrogens == 1 {
        carriers.push(StereoCarrier::ImplicitHydrogen);
    }
    Some(carriers)
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
        || bond.aromatic
        || double_bond_between_aromatic_atoms(mol, bond)
        || super::rings::bond_in_ring_smaller_than(mol, bond_id, 8)
        || (double_bond_is_in_ring(mol, bond_id) && double_bond_has_noncarbon_endpoint(mol, bond))
}

fn double_bond_between_aromatic_atoms(mol: &Molecule, bond: &Bond) -> bool {
    mol.atom(bond.a())
        .map(|atom| atom.aromatic)
        .unwrap_or(false)
        && mol
            .atom(bond.b())
            .map(|atom| atom.aromatic)
            .unwrap_or(false)
}

fn double_bond_is_in_ring(mol: &Molecule, bond: BondId) -> bool {
    mol.ring_membership()
        .map(|membership| membership.bond_in_ring(bond))
        .unwrap_or(false)
}

fn double_bond_has_noncarbon_endpoint(mol: &Molecule, bond: &Bond) -> bool {
    [bond.a(), bond.b()].into_iter().any(|atom_id| {
        mol.atom(atom_id)
            .map(|atom| atom.element.symbol() != "C")
            .unwrap_or(true)
    })
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
    if hydrogen_count(mol, endpoint) == 1 {
        carriers.push(StereoCarrier::ImplicitHydrogen);
    }
    carriers
}

fn assemble_tetrahedral_wedges(
    mol: &Molecule,
    issues: &mut Vec<StereoPerceptionIssue>,
    used_marks: &mut Vec<BondId>,
) -> Vec<StereoElement> {
    let mut marks = Vec::<TetrahedralWedgeMark<'_>>::new();
    for mark in mol.stereo_bond_marks() {
        if !matches!(
            mark.kind,
            StereoBondMarkKind::WedgeUp
                | StereoBondMarkKind::WedgeDown
                | StereoBondMarkKind::WedgeEither
        ) {
            continue;
        }
        if used_marks.contains(&mark.bond) {
            continue;
        }
        let Ok(bond) = mol.bond(mark.bond) else {
            continue;
        };
        if bond.order != BondOrder::Single {
            continue;
        }
        marks.push(TetrahedralWedgeMark {
            center: bond.a(),
            carrier: bond.b(),
            mark,
        });
    }
    marks.sort_by_key(|mark| (mark.center, mark.mark.bond));

    let mut assembled = Vec::new();
    let mut start = 0;
    while start < marks.len() {
        let center = marks[start].center;
        let end = marks[start..]
            .iter()
            .position(|mark| mark.center != center)
            .map_or(marks.len(), |offset| start + offset);
        let center_marks = &marks[start..end];
        if has_tetrahedral_element(mol, center) {
            used_marks.extend(center_marks.iter().map(|mark| mark.mark.bond));
            start = end;
            continue;
        }
        if center_marks.len() > 1 {
            issues.push(StereoPerceptionIssue::AmbiguousTetrahedralWedgeMarks {
                center,
                mark_count: center_marks.len(),
            });
            start = end;
            continue;
        }
        let mark = center_marks[0];
        if let Some(carriers) = tetrahedral_carriers_from_wedge(mol, center, mark.carrier) {
            used_marks.push(mark.mark.bond);
            assembled.push(tetrahedral_element_from_wedge(
                mol, mark.mark, center, carriers,
            ));
        }
        start = end;
    }
    assembled
}

#[derive(Clone, Copy)]
struct TetrahedralWedgeMark<'a> {
    center: AtomId,
    carrier: AtomId,
    mark: &'a StereoBondMark,
}

fn tetrahedral_carriers_from_wedge(
    mol: &Molecule,
    center: AtomId,
    marked_carrier: AtomId,
) -> Option<Vec<StereoCarrier>> {
    let marked = StereoCarrier::Atom(marked_carrier);
    let mut carriers = tetrahedral_carriers(mol, center)?;
    if !carriers.contains(&marked) {
        return None;
    }
    carriers.retain(|carrier| *carrier != marked);
    carriers.insert(0, marked);
    Some(carriers)
}

fn tetrahedral_element_from_wedge(
    mol: &Molecule,
    mark: &StereoBondMark,
    center: AtomId,
    carriers: Vec<StereoCarrier>,
) -> StereoElement {
    let (specifiedness, orientation) = match mark.kind {
        StereoBondMarkKind::WedgeUp | StereoBondMarkKind::WedgeDown => {
            let orientation = tetrahedral_wedge_orientation(mol, center, &carriers, mark.kind)
                .unwrap_or_else(|| match mark.kind {
                    StereoBondMarkKind::WedgeUp => TetrahedralOrientation::CounterClockwise,
                    StereoBondMarkKind::WedgeDown => TetrahedralOrientation::Clockwise,
                    _ => unreachable!("wedge orientation branch received non-wedge mark"),
                });
            (StereoSpecifiedness::Specified, orientation)
        }
        StereoBondMarkKind::WedgeEither => (
            StereoSpecifiedness::Unknown,
            TetrahedralOrientation::Clockwise,
        ),
        _ => unreachable!("non-wedge mark passed to tetrahedral wedge assembly"),
    };
    StereoElement {
        kind: StereoElementKind::Tetrahedral(TetrahedralStereo {
            center,
            carriers,
            orientation,
        }),
        specifiedness,
        source: mark.source,
        group: None,
        descriptor: None,
    }
}

fn tetrahedral_wedge_orientation(
    mol: &Molecule,
    center: AtomId,
    carriers: &[StereoCarrier],
    kind: StereoBondMarkKind,
) -> Option<TetrahedralOrientation> {
    let (_, conformer) = mol.first_conformer()?;
    let out_of_plane = match kind {
        StereoBondMarkKind::WedgeUp => 1.0,
        StereoBondMarkKind::WedgeDown => -1.0,
        _ => return None,
    };
    if let Some(atom_carriers) = carriers
        .iter()
        .map(|carrier| match carrier {
            StereoCarrier::Atom(atom) => Some(*atom),
            StereoCarrier::ImplicitHydrogen | StereoCarrier::ImplicitLonePair => None,
        })
        .collect::<Option<Vec<_>>>()
    {
        let mut points = tetrahedral_points(conformer, center, &atom_carriers)?;
        if matches!(coordinate_source(&points), StereoSource::Coordinates3D) {
            return tetrahedral_orientation_from_points(points);
        }
        points[1].z += out_of_plane;
        return tetrahedral_orientation_from_points(points);
    }
    let points =
        tetrahedral_points_with_virtual_implicit_h(conformer, center, carriers, out_of_plane)?;
    tetrahedral_orientation_from_points(points)
}

fn tetrahedral_points_with_virtual_implicit_h(
    conformer: &Conformer,
    center: AtomId,
    carriers: &[StereoCarrier],
    out_of_plane: f64,
) -> Option<[Point3; 5]> {
    (carriers.len() == 4).then_some(())?;
    let missing_hydrogen = carriers
        .iter()
        .enumerate()
        .filter_map(|(index, carrier)| {
            matches!(carrier, StereoCarrier::ImplicitHydrogen).then_some(index)
        })
        .collect::<Vec<_>>();
    if missing_hydrogen.len() != 1
        || carriers
            .iter()
            .any(|carrier| matches!(carrier, StereoCarrier::ImplicitLonePair))
    {
        return None;
    }

    let center_point = conformer.position(center)?;
    let mut carrier_points = [None; 4];
    for (index, carrier) in carriers.iter().enumerate() {
        if let StereoCarrier::Atom(atom) = carrier {
            carrier_points[index] = Some(conformer.position(*atom)?);
        }
    }

    let mut explicit_points = vec![center_point];
    explicit_points.extend(carrier_points.iter().filter_map(|point| *point));
    if matches!(
        coordinate_source(&explicit_points),
        StereoSource::Coordinates2D
    ) {
        carrier_points[0].as_mut()?.z += out_of_plane;
    }

    let mut vector_sum = Point3::new(0.0, 0.0, 0.0);
    for point in carrier_points.iter().filter_map(|point| *point) {
        let vector = vector_between(center_point, point);
        vector_sum.x += vector.x;
        vector_sum.y += vector.y;
        vector_sum.z += vector.z;
    }
    carrier_points[missing_hydrogen[0]] = Some(Point3::new(
        center_point.x - vector_sum.x,
        center_point.y - vector_sum.y,
        center_point.z - vector_sum.z,
    ));

    Some([
        center_point,
        carrier_points[0]?,
        carrier_points[1]?,
        carrier_points[2]?,
        carrier_points[3]?,
    ])
}

fn assemble_atropisomeric_axes(mol: &Molecule, used_marks: &mut Vec<BondId>) -> Vec<StereoElement> {
    let ring_membership = mol
        .ring_membership()
        .cloned()
        .unwrap_or_else(|| super::rings::compute_ring_membership(mol).0);
    let mut assembled = Vec::new();
    let mut assembled_axes = Vec::<BondId>::new();
    for mark in mol.stereo_bond_marks() {
        if used_marks.contains(&mark.bond)
            || !matches!(
                mark.kind,
                StereoBondMarkKind::WedgeUp | StereoBondMarkKind::WedgeDown
            )
        {
            continue;
        }

        let candidates = atropisomeric_axis_candidates(mol, &ring_membership, mark);
        if candidates.len() != 1 {
            continue;
        };
        let (axis, element) = candidates
            .into_iter()
            .next()
            .expect("one atrop axis candidate");
        used_marks.push(mark.bond);
        if assembled_axes.contains(&axis) {
            continue;
        }
        if let StereoElementKind::Axis(stereo) = &element.kind {
            assembled_axes.push(stereo.axis);
        }
        assembled.push(element);
    }
    assembled
}

fn atropisomeric_axis_candidates(
    mol: &Molecule,
    ring_membership: &RingMembership,
    mark: &StereoBondMark,
) -> Vec<(BondId, StereoElement)> {
    let Ok(marked_bond) = mol.bond(mark.bond) else {
        return Vec::new();
    };
    if marked_bond.order != BondOrder::Single {
        return Vec::new();
    }
    let near = marked_bond.a();
    let marked_carrier = marked_bond.b();
    mol.incident_bonds(near)
        .ok()
        .into_iter()
        .flatten()
        .filter_map(|(axis, bond)| {
            if axis == mark.bond {
                return None;
            }
            atropisomeric_axis_candidate(
                mol,
                ring_membership,
                mark,
                axis,
                bond,
                near,
                marked_carrier,
            )
            .map(|element| (axis, element))
        })
        .collect()
}

fn atropisomeric_axis_candidate(
    mol: &Molecule,
    ring_membership: &RingMembership,
    mark: &StereoBondMark,
    axis: BondId,
    axis_bond: &Bond,
    near: AtomId,
    marked_carrier: AtomId,
) -> Option<StereoElement> {
    if axis_bond.order != BondOrder::Single
        || ring_membership.bond_in_ring(axis)
        || has_axis_element(mol, axis)
    {
        return None;
    }
    let other = axis_bond.other_atom(near);
    if !atom_is_atropisomeric_sp2_endpoint(mol, ring_membership, near)
        || !atom_is_atropisomeric_sp2_endpoint(mol, ring_membership, other)
    {
        return None;
    }
    let left = axis_bond.a();
    let right = axis_bond.b();
    let left_carriers = atom_axis_carriers(mol, left, axis)?;
    let right_carriers = atom_axis_carriers(mol, right, axis)?;
    let marked_endpoint_carriers = if near == left {
        &left_carriers
    } else {
        &right_carriers
    };
    if left_carriers.len() != 2
        || right_carriers.len() != 2
        || !marked_endpoint_carriers.contains(&marked_carrier)
    {
        return None;
    }
    let left_reference = left_carriers[0];
    let right_reference = right_carriers[0];
    let orientation = axis_orientation_from_wedge(
        mol,
        axis_bond,
        left_reference,
        right_reference,
        near,
        marked_carrier,
        mark.kind,
    )?;
    Some(StereoElement {
        kind: StereoElementKind::Axis(AxisStereo {
            axis,
            carriers: vec![
                StereoCarrier::Atom(left_reference),
                StereoCarrier::Atom(right_reference),
            ],
            orientation,
        }),
        specifiedness: StereoSpecifiedness::Specified,
        source: mark.source,
        group: None,
        descriptor: None,
    })
}

fn atom_axis_carriers(mol: &Molecule, endpoint: AtomId, axis: BondId) -> Option<Vec<AtomId>> {
    let mut carriers = Vec::new();
    for (bond_id, bond) in mol.incident_bonds(endpoint).ok()? {
        if bond_id != axis {
            carriers.push(bond.other_atom(endpoint));
        }
    }
    carriers.sort();
    Some(carriers)
}

fn axis_orientation_from_wedge(
    mol: &Molecule,
    axis_bond: &Bond,
    left_reference: AtomId,
    right_reference: AtomId,
    marked_endpoint: AtomId,
    marked_carrier: AtomId,
    kind: StereoBondMarkKind,
) -> Option<AxisOrientation> {
    let (_, conformer) = mol.first_conformer()?;
    let (left, right) = axis_bond.endpoints();
    let left_point = conformer.position(left)?;
    let right_point = conformer.position(right)?;
    let mut left_reference_point = conformer.position(left_reference)?;
    let mut right_reference_point = conformer.position(right_reference)?;
    let marked_endpoint_point = conformer.position(marked_endpoint)?;
    let marked_point = conformer.position(marked_carrier)?;
    let coordinate_points = [
        left_point,
        right_point,
        left_reference_point,
        right_reference_point,
        marked_endpoint_point,
        marked_point,
    ];
    let axis = vector_between(left_point, right_point);
    if matches!(
        coordinate_source(&coordinate_points),
        StereoSource::Coordinates2D
    ) {
        let z_sign = match kind {
            StereoBondMarkKind::WedgeUp => 1.0,
            StereoBondMarkKind::WedgeDown => -1.0,
            _ => return None,
        };
        let marked_side = planar_cross(axis, vector_between(marked_endpoint_point, marked_point));
        if marked_side.abs() <= COORDINATE_EPSILON {
            return None;
        }
        left_reference_point.z += axis_reference_z_offset(
            axis,
            left_point,
            left_reference_point,
            left == marked_endpoint,
            marked_side,
            z_sign,
        )?;
        right_reference_point.z += axis_reference_z_offset(
            axis,
            right_point,
            right_reference_point,
            right == marked_endpoint,
            marked_side,
            z_sign,
        )?;
    } else if !matches!(
        kind,
        StereoBondMarkKind::WedgeUp | StereoBondMarkKind::WedgeDown
    ) {
        return None;
    }

    let left_vector = vector_between(left_point, left_reference_point);
    let right_vector = vector_between(right_point, right_reference_point);
    let handedness = dot(axis, cross(left_vector, right_vector));
    if handedness.abs() <= COORDINATE_EPSILON {
        return None;
    }
    Some(if handedness > 0.0 {
        AxisOrientation::Clockwise
    } else {
        AxisOrientation::CounterClockwise
    })
}

fn axis_reference_z_offset(
    axis: Point3,
    endpoint_point: Point3,
    reference_point: Point3,
    same_endpoint_as_mark: bool,
    marked_side: f64,
    marked_z: f64,
) -> Option<f64> {
    let side = planar_cross(axis, vector_between(endpoint_point, reference_point));
    if side.abs() <= COORDINATE_EPSILON {
        return None;
    }
    let same_side_as_mark = side.signum() == marked_side.signum();
    let side_factor = if same_side_as_mark { 1.0 } else { -1.0 };
    let endpoint_factor = if same_endpoint_as_mark { 1.0 } else { -1.0 };
    Some(marked_z * side_factor * endpoint_factor)
}

fn has_axis_element(mol: &Molecule, axis: BondId) -> bool {
    mol.stereo_elements().any(|(_, element)| {
        matches!(
            &element.kind,
            StereoElementKind::Axis(stereo) if stereo.axis == axis
        )
    })
}

fn atom_is_atropisomeric_sp2_endpoint(
    mol: &Molecule,
    ring_membership: &RingMembership,
    atom_id: AtomId,
) -> bool {
    let Ok(atom) = mol.atom(atom_id) else {
        return false;
    };
    let incident = mol
        .incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    let total_degree = incident
        .len()
        .saturating_add(usize::from(hydrogen_count(mol, atom_id)));
    if !(2..=3).contains(&total_degree) {
        return false;
    }
    ring_membership.atom_in_ring(atom_id)
        || atom.aromatic
        || incident.iter().any(|(_, bond)| {
            bond.aromatic || matches!(bond.order, BondOrder::Double | BondOrder::Aromatic)
        })
}

fn assemble_directional_double_bonds(
    mol: &Molecule,
    issues: &mut Vec<StereoPerceptionIssue>,
    used_marks: &mut Vec<BondId>,
) -> Vec<StereoElement> {
    let mut assembled = Vec::new();
    for (bond_id, bond) in mol.bonds() {
        if double_bond_stereo_is_unsupported(mol, bond_id, bond) {
            continue;
        }
        let left = bond.a();
        let right = bond.b();
        let left_marks = directional_marks_for_endpoint(mol, left, bond_id);
        let right_marks = directional_marks_for_endpoint(mol, right, bond_id);
        if has_double_bond_element(mol, bond_id) {
            used_marks.extend(left_marks.iter().map(|mark| mark.bond));
            used_marks.extend(right_marks.iter().map(|mark| mark.bond));
            continue;
        }
        let Some(left_mark) =
            select_directional_mark(mol, left, right, bond_id, &left_marks, issues)
        else {
            continue;
        };
        let Some(right_mark) =
            select_directional_mark(mol, right, left, bond_id, &right_marks, issues)
        else {
            continue;
        };
        let orientation = if left_mark.direction == right_mark.direction {
            DoubleBondOrientation::Together
        } else {
            DoubleBondOrientation::Opposite
        };
        used_marks.extend(left_marks.iter().map(|mark| mark.bond));
        used_marks.extend(right_marks.iter().map(|mark| mark.bond));
        assembled.push(StereoElement::specified(
            StereoElementKind::DoubleBond(DoubleBondStereo {
                bond: bond_id,
                left,
                right,
                left_carrier: StereoCarrier::Atom(left_mark.carrier),
                right_carrier: StereoCarrier::Atom(right_mark.carrier),
                orientation,
            }),
            common_source(left_mark.mark.source, right_mark.mark.source),
        ));
    }
    assembled
}

fn select_directional_mark<'a>(
    mol: &Molecule,
    endpoint: AtomId,
    other_endpoint: AtomId,
    focus_bond: BondId,
    marks: &'a [EndpointMark<'a>],
    issues: &mut Vec<StereoPerceptionIssue>,
) -> Option<EndpointMark<'a>> {
    match marks {
        [] => None,
        [mark] => Some(*mark),
        [first, second]
            if redundant_endpoint_directional_marks(
                mol,
                endpoint,
                other_endpoint,
                focus_bond,
                first,
                second,
            ) =>
        {
            Some(*first)
        }
        _ => {
            issues.push(StereoPerceptionIssue::AmbiguousDirectionalBondMarks {
                double_bond: focus_bond,
                endpoint,
                mark_count: marks.len(),
            });
            None
        }
    }
}

fn redundant_endpoint_directional_marks(
    mol: &Molecule,
    endpoint: AtomId,
    other_endpoint: AtomId,
    focus_bond: BondId,
    first: &EndpointMark<'_>,
    second: &EndpointMark<'_>,
) -> bool {
    if first.direction == second.direction {
        return false;
    }
    let mut marked_carriers = [first.carrier, second.carrier];
    marked_carriers.sort_unstable();
    let mut atom_carriers =
        double_bond_endpoint_carriers(mol, endpoint, other_endpoint, focus_bond)
            .into_iter()
            .filter_map(|carrier| match carrier {
                StereoCarrier::Atom(atom) => Some(atom),
                StereoCarrier::ImplicitHydrogen | StereoCarrier::ImplicitLonePair => None,
            })
            .collect::<Vec<_>>();
    atom_carriers.sort_unstable();
    marked_carriers.as_slice() == atom_carriers.as_slice()
}

fn report_unassembled_source_marks(
    mol: &Molecule,
    used_marks: &[BondId],
    issues: &mut Vec<StereoPerceptionIssue>,
) {
    for mark in mol.stereo_bond_marks() {
        match mark.kind {
            StereoBondMarkKind::DirectionalUp | StereoBondMarkKind::DirectionalDown => {
                if !used_marks.contains(&mark.bond) {
                    issues.push(StereoPerceptionIssue::UnpairedDirectionalBondMark {
                        bond: mark.bond,
                    });
                }
            }
            StereoBondMarkKind::WedgeUp
            | StereoBondMarkKind::WedgeDown
            | StereoBondMarkKind::WedgeEither => {
                if !used_marks.contains(&mark.bond) {
                    issues.push(StereoPerceptionIssue::UnassembledTetrahedralBondMark {
                        bond: mark.bond,
                        kind: mark.kind,
                    });
                }
            }
            StereoBondMarkKind::DoubleBondEither => {
                issues.push(StereoPerceptionIssue::UnsupportedSourceBondMark {
                    bond: mark.bond,
                    kind: mark.kind,
                });
            }
        }
    }
}

fn assign_coordinate_stereo(
    mol: &Molecule,
    planned_elements: &[StereoElement],
) -> Vec<StereoElement> {
    let Some((_, conformer)) = mol.first_conformer() else {
        return Vec::new();
    };
    let mut assigned = Vec::new();
    assigned.extend(assign_coordinate_tetrahedral(
        mol,
        conformer,
        planned_elements,
    ));
    assigned.extend(assign_coordinate_double_bonds(
        mol,
        conformer,
        planned_elements,
    ));
    assigned
}

fn assign_coordinate_tetrahedral(
    mol: &Molecule,
    conformer: &Conformer,
    planned_elements: &[StereoElement],
) -> Vec<StereoElement> {
    let mut assigned = Vec::new();
    for candidate in tetrahedral_candidates(mol) {
        let StereoCandidate::Tetrahedral { center, carriers } = candidate else {
            continue;
        };
        if has_tetrahedral_element(mol, center)
            || planned_elements
                .iter()
                .any(|element| planned_tetrahedral_center(element) == Some(center))
        {
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
        assigned.push(StereoElement::specified(
            StereoElementKind::Tetrahedral(TetrahedralStereo {
                center,
                carriers,
                orientation,
            }),
            StereoSource::Coordinates3D,
        ));
    }
    assigned
}

fn assign_coordinate_double_bonds(
    mol: &Molecule,
    conformer: &Conformer,
    planned_elements: &[StereoElement],
) -> Vec<StereoElement> {
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
        if has_double_bond_element(mol, bond)
            || planned_elements
                .iter()
                .any(|element| planned_double_bond(element) == Some(bond))
        {
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
        assigned.push(StereoElement::specified(
            StereoElementKind::DoubleBond(DoubleBondStereo {
                bond,
                left,
                right,
                left_carrier: StereoCarrier::Atom(left_carrier),
                right_carrier: StereoCarrier::Atom(right_carrier),
                orientation,
            }),
            coordinate_source(&points),
        ));
    }
    assigned
}

fn planned_tetrahedral_center(element: &StereoElement) -> Option<AtomId> {
    match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => Some(stereo.center),
        _ => None,
    }
}

fn planned_double_bond(element: &StereoElement) -> Option<BondId> {
    match &element.kind {
        StereoElementKind::DoubleBond(stereo) => Some(stereo.bond),
        _ => None,
    }
}

fn only_atom_carrier(carriers: &[StereoCarrier]) -> Option<AtomId> {
    let mut atoms = carriers.iter().filter_map(|carrier| match carrier {
        StereoCarrier::Atom(atom) => Some(*atom),
        StereoCarrier::ImplicitHydrogen | StereoCarrier::ImplicitLonePair => None,
    });
    let atom = atoms.next()?;
    atoms.next().is_none().then_some(atom)
}

fn tetrahedral_points(
    conformer: &Conformer,
    center: AtomId,
    carriers: &[AtomId],
) -> Option<[Point3; 5]> {
    (carriers.len() == 4).then_some(())?;
    Some([
        conformer.position(center)?,
        conformer.position(carriers[0])?,
        conformer.position(carriers[1])?,
        conformer.position(carriers[2])?,
        conformer.position(carriers[3])?,
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
        conformer.position(left)?,
        conformer.position(right)?,
        conformer.position(left_carrier)?,
        conformer.position(right_carrier)?,
    ])
}

fn tetrahedral_orientation_from_points(points: [Point3; 5]) -> Option<TetrahedralOrientation> {
    let a = vector_between(points[4], points[1]);
    let b = vector_between(points[4], points[2]);
    let c = vector_between(points[4], points[3]);
    let volume = dot(cross(a, b), c);
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
    let axis = vector_between(points[0], points[1]);
    let left_vector = vector_between(points[0], points[2]);
    let right_vector = vector_between(points[1], points[3]);
    let sidedness = dot(cross(axis, left_vector), cross(axis, right_vector));
    if sidedness.abs() <= COORDINATE_EPSILON {
        return None;
    }
    Some(if sidedness > 0.0 {
        DoubleBondOrientation::Together
    } else {
        DoubleBondOrientation::Opposite
    })
}

fn coordinate_source(points: &[Point3]) -> StereoSource {
    let Some(first) = points.first() else {
        return StereoSource::Coordinates3D;
    };
    if points
        .iter()
        .all(|point| (point.z - first.z).abs() <= COORDINATE_EPSILON)
    {
        StereoSource::Coordinates2D
    } else {
        StereoSource::Coordinates3D
    }
}

fn vector_between(origin: Point3, point: Point3) -> Point3 {
    Point3::new(point.x - origin.x, point.y - origin.y, point.z - origin.z)
}

fn planar_cross(a: Point3, b: Point3) -> f64 {
    a.x * b.y - a.y * b.x
}

fn cross(a: Point3, b: Point3) -> Point3 {
    Point3::new(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    )
}

fn dot(a: Point3, b: Point3) -> f64 {
    a.x * b.x + a.y * b.y + a.z * b.z
}

const COORDINATE_EPSILON: f64 = 1.0e-8;

#[derive(Clone, Copy)]
struct EndpointMark<'a> {
    bond: BondId,
    carrier: AtomId,
    mark: &'a StereoBondMark,
    direction: StereoBondMarkKind,
}

fn directional_marks_for_endpoint(
    mol: &Molecule,
    endpoint: AtomId,
    focus_bond: BondId,
) -> Vec<EndpointMark<'_>> {
    let mut marks = Vec::new();
    let Ok(incident) = mol.incident_bonds(endpoint) else {
        return marks;
    };
    for (bond_id, bond) in incident {
        if bond_id == focus_bond || bond.order != BondOrder::Single {
            continue;
        }
        let Some(mark) = mol.stereo_bond_mark(bond_id) else {
            continue;
        };
        if matches!(
            mark.kind,
            StereoBondMarkKind::DirectionalUp | StereoBondMarkKind::DirectionalDown
        ) {
            marks.push(EndpointMark {
                bond: bond_id,
                carrier: bond.other_atom(endpoint),
                mark,
                direction: directional_mark_at_endpoint(mark.kind, bond, endpoint),
            });
        }
    }
    marks.sort_by_key(|mark| (mark.bond, mark.carrier));
    marks
}

fn directional_mark_at_endpoint(
    kind: StereoBondMarkKind,
    bond: &Bond,
    endpoint: AtomId,
) -> StereoBondMarkKind {
    if bond.b() == endpoint {
        invert_directional_mark(kind)
    } else {
        kind
    }
}

fn invert_directional_mark(kind: StereoBondMarkKind) -> StereoBondMarkKind {
    match kind {
        StereoBondMarkKind::DirectionalUp => StereoBondMarkKind::DirectionalDown,
        StereoBondMarkKind::DirectionalDown => StereoBondMarkKind::DirectionalUp,
        _ => kind,
    }
}

fn has_double_bond_element(mol: &Molecule, bond: BondId) -> bool {
    mol.stereo_elements().any(|(_, element)| {
        matches!(
            &element.kind,
            StereoElementKind::DoubleBond(stereo) if stereo.bond == bond
        )
    })
}

fn has_tetrahedral_element(mol: &Molecule, center: AtomId) -> bool {
    mol.stereo_elements().any(|(_, element)| {
        matches!(
            &element.kind,
            StereoElementKind::Tetrahedral(stereo) if stereo.center == center
        )
    })
}

fn hydrogen_count(mol: &Molecule, atom: AtomId) -> u8 {
    let Ok(atom) = mol.atom(atom) else {
        return 0;
    };
    atom.explicit_hydrogens
        .saturating_add(atom.implicit_hydrogens.unwrap_or(0))
}

fn implicit_lone_pair_available(mol: &Molecule, atom: AtomId) -> bool {
    mol.atom(atom)
        .map(|atom_payload| {
            matches!(
                atom_payload.element.symbol(),
                "N" | "P" | "As" | "Sb" | "O" | "S" | "Se" | "Te"
            ) && hydrogen_count(mol, atom) == 0
        })
        .unwrap_or(false)
}

fn carrier_key(carrier: &StereoCarrier) -> (u8, u32) {
    match carrier {
        StereoCarrier::Atom(atom) => (0, atom.raw()),
        StereoCarrier::ImplicitHydrogen => (1, u32::MAX),
        StereoCarrier::ImplicitLonePair => (2, u32::MAX),
    }
}

fn bond_connects(bond: &Bond, a: AtomId, b: AtomId) -> bool {
    (bond.a() == a && bond.b() == b) || (bond.a() == b && bond.b() == a)
}

fn common_source(left: StereoSource, right: StereoSource) -> StereoSource {
    if left == right {
        left
    } else {
        StereoSource::User
    }
}
