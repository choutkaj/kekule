use std::collections::{BTreeMap, BTreeSet};

use crate::algorithms::{
    compute_graph_ring_membership, graph_bond_in_ring_smaller_than, validate_stereo,
    RingMembership, StereoValidationError,
};
use crate::core::*;
use crate::geometry::Point3;

use super::normalization::{
    NormalizationReport, NormalizationWarning, SourceStereoNormalizationError,
    SourceStereoNormalizationIssue,
};

/// One detached source-format bond stereo assertion awaiting interpretation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct SourceStereoBondMark {
    pub(crate) bond: BondId,
    pub(crate) from: AtomId,
    pub(crate) kind: SourceStereoBondMarkKind,
}

/// A format-local Molfile bond mark and the endpoint that must be emitted first.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct MolfileStereoBondProjection {
    pub(crate) from: AtomId,
    pub(crate) kind: SourceStereoBondMarkKind,
}

/// Source-format bond stereo syntax interpreted before molecule publication.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub(crate) enum SourceStereoBondMarkKind {
    DirectionalUp,
    DirectionalDown,
    WedgeUp,
    WedgeDown,
    WedgeEither,
    DoubleBondEither,
}

/// Canonicalize source-declared stereo using represented state only.
///
/// This kernel is deliberately separate from coordinate perception. It reads
/// coordinates only when a source wedge/hash mark needs its format-local
/// drawing geometry decoded; it never infers stereo from unmarked coordinates.
pub(super) fn normalize_source_stereo(
    molecule: &mut Molecule,
    source_marks: &[SourceStereoBondMark],
) -> std::result::Result<NormalizationReport, SourceStereoNormalizationError> {
    validate_stereo(molecule).map_err(source_validation_error)?;
    validate_source_stereo_marks(molecule, source_marks)?;

    let ring_membership = compute_graph_ring_membership(molecule);
    let mut warnings = Vec::new();
    let mut issues = Vec::new();
    let mut used_marks = Vec::<BondId>::new();
    let axis_elements =
        assemble_atropisomeric_axes(molecule, &ring_membership, source_marks, &mut used_marks);
    let mut planned_elements =
        assemble_tetrahedral_wedges(molecule, source_marks, &mut warnings, &mut used_marks);
    planned_elements.extend(axis_elements);
    planned_elements.extend(assemble_unknown_double_bonds(
        molecule,
        source_marks,
        &mut used_marks,
    ));
    planned_elements.extend(assemble_directional_double_bonds(
        molecule,
        &ring_membership,
        source_marks,
        &mut issues,
        &mut used_marks,
    ));
    report_unassembled_source_marks(source_marks, &used_marks, &mut issues);
    if !issues.is_empty() {
        return Err(SourceStereoNormalizationError { issues });
    }

    let mut created_stereo_elements = Vec::with_capacity(planned_elements.len());
    for element in planned_elements {
        match molecule.add_stereo_element(element) {
            Ok(id) => created_stereo_elements.push(id),
            Err(error) => {
                return Err(SourceStereoNormalizationError {
                    issues: vec![SourceStereoNormalizationIssue::CouldNotCreateStereoElement(
                        error,
                    )],
                });
            }
        }
    }
    validate_stereo(molecule).map_err(source_validation_error)?;

    Ok(NormalizationReport {
        created_stereo_elements,
        warnings,
    })
}

/// Project canonical stereo elements into format-local Molfile bond marks.
///
/// The marks are detached writer state and are never stored back on the
/// canonical molecule. Candidate projections are verified by decoding them
/// with the same source-stereo kernel used during interpretation.
pub(crate) fn project_molfile_stereo_bond_marks(
    molecule: &Molecule,
) -> std::result::Result<BTreeMap<BondId, MolfileStereoBondProjection>, String> {
    if molecule.stereo_elements().next().is_none() {
        return Ok(BTreeMap::new());
    }
    if molecule.stereo_groups().next().is_some() {
        return Err("Molfile writer does not support stereo groups".to_owned());
    }

    let mut projected = BTreeMap::new();
    let mut occupied = BTreeSet::new();
    for (_, target) in molecule.stereo_elements() {
        let candidates = molfile_projection_candidates(molecule, target)?;
        let mut selected = None;
        for (bond, projection) in candidates {
            if occupied.contains(&bond) {
                continue;
            }
            let mut staged = molecule.clone();
            staged.stereo_elements.clear();
            staged.stereo_groups.clear();
            staged.perception = PerceptionState::default();
            let source_marks = [SourceStereoBondMark {
                bond,
                from: projection.from,
                kind: projection.kind,
            }];
            let Ok(report) = normalize_source_stereo(&mut staged, &source_marks) else {
                continue;
            };
            if !report.warnings.is_empty() || report.created_stereo_elements.len() != 1 {
                continue;
            }
            let Ok(decoded) = staged.stereo_element(report.created_stereo_elements[0]) else {
                continue;
            };
            if decoded.kind == target.kind {
                selected = Some((bond, projection));
                break;
            }
        }
        let Some((bond, projection)) = selected else {
            return Err(format!(
                "Molfile writer cannot encode canonical stereo element {:?}",
                target.kind
            ));
        };
        occupied.insert(bond);
        projected.insert(bond, projection);
    }
    Ok(projected)
}

fn molfile_projection_candidates(
    molecule: &Molecule,
    element: &StereoElement,
) -> std::result::Result<Vec<(BondId, MolfileStereoBondProjection)>, String> {
    match &element.kind {
        StereoElementKind::Tetrahedral(stereo) => {
            let kinds = match stereo.orientation {
                Some(_) => vec![
                    SourceStereoBondMarkKind::WedgeUp,
                    SourceStereoBondMarkKind::WedgeDown,
                ],
                None => vec![SourceStereoBondMarkKind::WedgeEither],
            };
            Ok(molecule
                .incident_bonds(stereo.center)
                .map_err(|error| error.to_string())?
                .filter(|(_, bond)| bond.order == BondOrder::Single)
                .flat_map(|(bond, _)| {
                    kinds.iter().copied().map(move |kind| {
                        (
                            bond,
                            MolfileStereoBondProjection {
                                from: stereo.center,
                                kind,
                            },
                        )
                    })
                })
                .collect())
        }
        StereoElementKind::Axis(stereo) => {
            if stereo.orientation.is_none() {
                return Err("Molfile writer cannot encode unknown axis stereo".to_owned());
            }
            let kinds = [
                SourceStereoBondMarkKind::WedgeUp,
                SourceStereoBondMarkKind::WedgeDown,
            ];
            let axis = molecule
                .bond(stereo.axis)
                .map_err(|error| error.to_string())?;
            Ok([axis.a(), axis.b()]
                .into_iter()
                .flat_map(|endpoint| {
                    molecule
                        .incident_bonds(endpoint)
                        .ok()
                        .into_iter()
                        .flatten()
                        .map(move |(bond, value)| (endpoint, bond, value))
                })
                .filter(|(_, bond, value)| *bond != stereo.axis && value.order == BondOrder::Single)
                .flat_map(|(from, bond, _)| {
                    kinds
                        .iter()
                        .copied()
                        .map(move |kind| (bond, MolfileStereoBondProjection { from, kind }))
                })
                .collect())
        }
        StereoElementKind::DoubleBond(stereo) if stereo.orientation.is_none() => Ok(vec![(
            stereo.bond,
            MolfileStereoBondProjection {
                from: stereo.left,
                kind: SourceStereoBondMarkKind::DoubleBondEither,
            },
        )]),
        StereoElementKind::DoubleBond(_) => {
            Err("Molfile writer cannot encode specified double-bond stereo".to_owned())
        }
    }
}

fn validate_source_stereo_marks(
    molecule: &Molecule,
    source_marks: &[SourceStereoBondMark],
) -> std::result::Result<(), SourceStereoNormalizationError> {
    let issues = source_marks
        .iter()
        .filter_map(|mark| match molecule.bond(mark.bond) {
            Err(_) => {
                Some(SourceStereoNormalizationIssue::MissingSourceBondMark { bond: mark.bond })
            }
            Ok(bond) if ![bond.a(), bond.b()].contains(&mark.from) => Some(
                SourceStereoNormalizationIssue::InvalidSourceBondMarkEndpoint {
                    bond: mark.bond,
                    from: mark.from,
                },
            ),
            Ok(_) => None,
        })
        .collect::<Vec<_>>();
    if issues.is_empty() {
        Ok(())
    } else {
        Err(SourceStereoNormalizationError { issues })
    }
}

fn source_validation_error(error: StereoValidationError) -> SourceStereoNormalizationError {
    SourceStereoNormalizationError {
        issues: error
            .issues
            .into_iter()
            .map(SourceStereoNormalizationIssue::InvalidStereo)
            .collect(),
    }
}

fn assemble_tetrahedral_wedges(
    molecule: &Molecule,
    source_marks: &[SourceStereoBondMark],
    warnings: &mut Vec<NormalizationWarning>,
    used_marks: &mut Vec<BondId>,
) -> Vec<StereoElement> {
    let mut marks = Vec::<TetrahedralWedgeMark<'_>>::new();
    for mark in source_marks_in_bond_order(source_marks) {
        if !matches!(
            mark.kind,
            SourceStereoBondMarkKind::WedgeUp
                | SourceStereoBondMarkKind::WedgeDown
                | SourceStereoBondMarkKind::WedgeEither
        ) || used_marks.contains(&mark.bond)
        {
            continue;
        }
        let Ok(bond) = molecule.bond(mark.bond) else {
            continue;
        };
        if bond.order != BondOrder::Single {
            continue;
        }
        marks.push(TetrahedralWedgeMark {
            center: mark.from,
            carrier: bond.other_atom(mark.from),
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
        if has_tetrahedral_element(molecule, center) {
            used_marks.extend(center_marks.iter().map(|mark| mark.mark.bond));
            start = end;
            continue;
        }
        if center_marks.len() > 1 {
            used_marks.extend(center_marks.iter().map(|mark| mark.mark.bond));
            warnings.push(NormalizationWarning::AmbiguousTetrahedralWedgeMarks {
                center,
                mark_count: center_marks.len(),
            });
            start = end;
            continue;
        }
        let mark = center_marks[0];
        if let Some(carriers) = tetrahedral_carriers_from_wedge(molecule, center, mark.carrier) {
            used_marks.push(mark.mark.bond);
            assembled.push(tetrahedral_element_from_wedge(
                molecule, mark.mark, center, carriers,
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
    mark: &'a SourceStereoBondMark,
}

fn tetrahedral_carriers_from_wedge(
    molecule: &Molecule,
    center: AtomId,
    marked_carrier: AtomId,
) -> Option<Vec<StereoCarrier>> {
    let marked = StereoCarrier::Atom(marked_carrier);
    let mut carriers = source_tetrahedral_carriers(molecule, center)?;
    if !carriers.contains(&marked) {
        return None;
    }
    carriers.retain(|carrier| *carrier != marked);
    carriers.insert(0, marked);
    Some(carriers)
}

fn source_tetrahedral_carriers(molecule: &Molecule, center: AtomId) -> Option<Vec<StereoCarrier>> {
    let atom = molecule.atom(center).ok()?;
    if atom.element.symbol() == "H" {
        return None;
    }
    let mut carriers = Vec::new();
    for (_, bond) in molecule.incident_bonds(center).ok()? {
        if matches!(bond.order, BondOrder::Zero | BondOrder::Dative) {
            return None;
        }
        carriers.push(StereoCarrier::Atom(bond.other_atom(center)));
    }
    carriers.sort_by_key(carrier_key);

    // Only a primary explicit-H declaration can assert a virtual carrier here.
    // Installed implicit-H assignments are deliberately invisible to this
    // normalization kernel.
    let declared_hydrogens = atom.hydrogens.explicit_count();
    if declared_hydrogens == 0
        && carriers.len() == 3
        && stable_tetrahedral_lone_pair_center(atom.element.symbol())
    {
        carriers.push(StereoCarrier::ImplicitLonePair);
        return Some(carriers);
    }
    if declared_hydrogens > 1 || carriers.len() + usize::from(declared_hydrogens) != 4 {
        return None;
    }
    if declared_hydrogens == 1 {
        carriers.push(StereoCarrier::ImplicitHydrogen);
    }
    Some(carriers)
}

fn stable_tetrahedral_lone_pair_center(symbol: &str) -> bool {
    matches!(symbol, "P" | "As" | "Sb" | "S" | "Se" | "Te")
}

fn tetrahedral_element_from_wedge(
    molecule: &Molecule,
    mark: &SourceStereoBondMark,
    center: AtomId,
    carriers: Vec<StereoCarrier>,
) -> StereoElement {
    let orientation = match mark.kind {
        SourceStereoBondMarkKind::WedgeUp | SourceStereoBondMarkKind::WedgeDown => {
            let orientation = tetrahedral_wedge_orientation(molecule, center, &carriers, mark.kind)
                .unwrap_or_else(|| match mark.kind {
                    SourceStereoBondMarkKind::WedgeUp => TetrahedralOrientation::CounterClockwise,
                    SourceStereoBondMarkKind::WedgeDown => TetrahedralOrientation::Clockwise,
                    _ => unreachable!("wedge orientation branch received non-wedge mark"),
                });
            Some(orientation)
        }
        SourceStereoBondMarkKind::WedgeEither => None,
        _ => unreachable!("non-wedge mark passed to tetrahedral wedge assembly"),
    };
    StereoElement::new(StereoElementKind::Tetrahedral(TetrahedralStereo {
        center,
        carriers,
        orientation,
    }))
}

fn tetrahedral_wedge_orientation(
    molecule: &Molecule,
    center: AtomId,
    carriers: &[StereoCarrier],
    kind: SourceStereoBondMarkKind,
) -> Option<TetrahedralOrientation> {
    let (_, conformer) = molecule.first_conformer()?;
    let out_of_plane = match kind {
        SourceStereoBondMarkKind::WedgeUp => 1.0,
        SourceStereoBondMarkKind::WedgeDown => -1.0,
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
        if !coordinates_are_planar(&points) {
            return tetrahedral_orientation_from_points(points);
        }
        points[1].z += out_of_plane;
        return tetrahedral_orientation_from_points(points);
    }
    let points = tetrahedral_points_with_virtual_declared_hydrogen(
        conformer,
        center,
        carriers,
        out_of_plane,
    )?;
    tetrahedral_orientation_from_points(points)
}

fn tetrahedral_points_with_virtual_declared_hydrogen(
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

    let center_point = conformer.position_value(center)?;
    let mut carrier_points = [None; 4];
    for (index, carrier) in carriers.iter().enumerate() {
        if let StereoCarrier::Atom(atom) = carrier {
            carrier_points[index] = Some(conformer.position_value(*atom)?);
        }
    }

    let mut explicit_points = vec![center_point];
    explicit_points.extend(carrier_points.iter().filter_map(|point| *point));
    if coordinates_are_planar(&explicit_points) {
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

fn assemble_atropisomeric_axes(
    molecule: &Molecule,
    ring_membership: &RingMembership,
    source_marks: &[SourceStereoBondMark],
    used_marks: &mut Vec<BondId>,
) -> Vec<StereoElement> {
    let mut assembled = Vec::new();
    let mut assembled_axes = Vec::<BondId>::new();
    for mark in source_marks_in_bond_order(source_marks) {
        if used_marks.contains(&mark.bond)
            || !matches!(
                mark.kind,
                SourceStereoBondMarkKind::WedgeUp | SourceStereoBondMarkKind::WedgeDown
            )
        {
            continue;
        }

        let candidates = atropisomeric_axis_candidates(molecule, ring_membership, mark);
        if candidates.len() != 1 {
            continue;
        }
        let (axis, element) = candidates
            .into_iter()
            .next()
            .expect("one atrop axis candidate");
        used_marks.push(mark.bond);
        if assembled_axes.contains(&axis) {
            continue;
        }
        assembled_axes.push(axis);
        assembled.push(element);
    }
    assembled
}

fn atropisomeric_axis_candidates(
    molecule: &Molecule,
    ring_membership: &RingMembership,
    mark: &SourceStereoBondMark,
) -> Vec<(BondId, StereoElement)> {
    let Ok(marked_bond) = molecule.bond(mark.bond) else {
        return Vec::new();
    };
    if marked_bond.order != BondOrder::Single {
        return Vec::new();
    }
    let near = mark.from;
    let marked_carrier = marked_bond.other_atom(mark.from);
    let mut candidates = molecule
        .incident_bonds(near)
        .ok()
        .into_iter()
        .flatten()
        .filter_map(|(axis, bond)| {
            if axis == mark.bond {
                return None;
            }
            atropisomeric_axis_candidate(
                molecule,
                ring_membership,
                mark,
                axis,
                bond,
                near,
                marked_carrier,
            )
            .map(|element| (axis, element))
        })
        .collect::<Vec<_>>();
    candidates.sort_by_key(|(axis, _)| *axis);
    let has_non_ring_axis = candidates
        .iter()
        .any(|(axis, _)| !ring_membership.bond_in_ring(*axis));
    candidates
        .into_iter()
        .filter(|(axis, _)| !has_non_ring_axis || !ring_membership.bond_in_ring(*axis))
        .collect()
}

fn atropisomeric_axis_candidate(
    molecule: &Molecule,
    ring_membership: &RingMembership,
    mark: &SourceStereoBondMark,
    axis: BondId,
    axis_bond: &Bond,
    near: AtomId,
    marked_carrier: AtomId,
) -> Option<StereoElement> {
    if axis_bond.order != BondOrder::Single || has_axis_element(molecule, axis) {
        return None;
    }
    let other = axis_bond.other_atom(near);
    if !source_atom_is_atropisomeric_sp2_endpoint(molecule, ring_membership, near)
        || !source_atom_is_atropisomeric_sp2_endpoint(molecule, ring_membership, other)
    {
        return None;
    }
    let left = axis_bond.a();
    let right = axis_bond.b();
    let left_carriers = atom_axis_carriers(molecule, left, axis)?;
    let right_carriers = atom_axis_carriers(molecule, right, axis)?;
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
        molecule,
        axis_bond,
        left_reference,
        right_reference,
        near,
        marked_carrier,
        mark.kind,
    )?;
    Some(StereoElement::new(StereoElementKind::Axis(AxisStereo {
        axis,
        carriers: vec![
            StereoCarrier::Atom(left_reference),
            StereoCarrier::Atom(right_reference),
        ],
        orientation: Some(orientation),
    })))
}

fn source_atom_is_atropisomeric_sp2_endpoint(
    molecule: &Molecule,
    ring_membership: &RingMembership,
    atom_id: AtomId,
) -> bool {
    let Ok(atom) = molecule.atom(atom_id) else {
        return false;
    };
    let incident = molecule
        .incident_bonds(atom_id)
        .ok()
        .into_iter()
        .flatten()
        .collect::<Vec<_>>();
    // Source-declared atom hydrogens participate in represented degree;
    // perceived implicit hydrogens do not.
    let total_degree = incident
        .len()
        .saturating_add(usize::from(atom.hydrogens.explicit_count()));
    if !(2..=3).contains(&total_degree) {
        return false;
    }
    ring_membership.atom_in_ring(atom_id)
        || incident
            .iter()
            .any(|(_, bond)| bond.order == BondOrder::Double)
}

fn atom_axis_carriers(molecule: &Molecule, endpoint: AtomId, axis: BondId) -> Option<Vec<AtomId>> {
    let mut carriers = Vec::new();
    for (bond_id, bond) in molecule.incident_bonds(endpoint).ok()? {
        if bond_id != axis {
            carriers.push(bond.other_atom(endpoint));
        }
    }
    carriers.sort_unstable();
    Some(carriers)
}

fn axis_orientation_from_wedge(
    molecule: &Molecule,
    axis_bond: &Bond,
    left_reference: AtomId,
    right_reference: AtomId,
    marked_endpoint: AtomId,
    marked_carrier: AtomId,
    kind: SourceStereoBondMarkKind,
) -> Option<AxisOrientation> {
    let (_, conformer) = molecule.first_conformer()?;
    let (left, right) = axis_bond.endpoints();
    let left_point = conformer.position_value(left)?;
    let right_point = conformer.position_value(right)?;
    let mut left_reference_point = conformer.position_value(left_reference)?;
    let mut right_reference_point = conformer.position_value(right_reference)?;
    let marked_endpoint_point = conformer.position_value(marked_endpoint)?;
    let marked_point = conformer.position_value(marked_carrier)?;
    let coordinate_points = [
        left_point,
        right_point,
        left_reference_point,
        right_reference_point,
        marked_endpoint_point,
        marked_point,
    ];
    let axis = vector_between(left_point, right_point);
    if coordinates_are_planar(&coordinate_points) {
        let z_sign = match kind {
            SourceStereoBondMarkKind::WedgeUp => 1.0,
            SourceStereoBondMarkKind::WedgeDown => -1.0,
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
    let side_factor = if side.signum() == marked_side.signum() {
        1.0
    } else {
        -1.0
    };
    let endpoint_factor = if same_endpoint_as_mark { 1.0 } else { -1.0 };
    Some(marked_z * side_factor * endpoint_factor)
}

fn assemble_directional_double_bonds(
    molecule: &Molecule,
    ring_membership: &RingMembership,
    source_marks: &[SourceStereoBondMark],
    issues: &mut Vec<SourceStereoNormalizationIssue>,
    used_marks: &mut Vec<BondId>,
) -> Vec<StereoElement> {
    let mut assembled = Vec::new();
    for (bond_id, bond) in molecule.bonds() {
        if source_double_bond_stereo_is_unsupported(molecule, ring_membership, bond_id, bond) {
            continue;
        }
        let left = bond.a();
        let right = bond.b();
        let left_marks = directional_marks_for_endpoint(molecule, source_marks, left, bond_id);
        let right_marks = directional_marks_for_endpoint(molecule, source_marks, right, bond_id);
        if has_double_bond_element(molecule, bond_id) {
            used_marks.extend(left_marks.iter().map(|mark| mark.bond));
            used_marks.extend(right_marks.iter().map(|mark| mark.bond));
            continue;
        }
        let Some(left_mark) =
            select_directional_mark(molecule, left, right, bond_id, &left_marks, issues)
        else {
            continue;
        };
        let Some(right_mark) =
            select_directional_mark(molecule, right, left, bond_id, &right_marks, issues)
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
        assembled.push(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond: bond_id,
                left,
                right,
                left_carrier: StereoCarrier::Atom(left_mark.carrier),
                right_carrier: StereoCarrier::Atom(right_mark.carrier),
                orientation: Some(orientation),
            },
        )));
    }
    assembled
}

fn assemble_unknown_double_bonds(
    molecule: &Molecule,
    source_marks: &[SourceStereoBondMark],
    used_marks: &mut Vec<BondId>,
) -> Vec<StereoElement> {
    let mut assembled = Vec::new();
    for mark in source_marks_in_bond_order(source_marks) {
        if mark.kind != SourceStereoBondMarkKind::DoubleBondEither {
            continue;
        }
        let Ok(bond) = molecule.bond(mark.bond) else {
            continue;
        };
        if bond.order != BondOrder::Double {
            continue;
        }
        if has_double_bond_element(molecule, mark.bond) {
            used_marks.push(mark.bond);
            continue;
        }

        let left = bond.a();
        let right = bond.b();
        let Some(left_carrier) =
            source_double_bond_endpoint_carriers(molecule, left, right, mark.bond)
                .into_iter()
                .next()
        else {
            continue;
        };
        let Some(right_carrier) =
            source_double_bond_endpoint_carriers(molecule, right, left, mark.bond)
                .into_iter()
                .next()
        else {
            continue;
        };

        used_marks.push(mark.bond);
        assembled.push(StereoElement::new(StereoElementKind::DoubleBond(
            DoubleBondStereo {
                bond: mark.bond,
                left,
                right,
                left_carrier,
                right_carrier,
                orientation: None,
            },
        )));
    }
    assembled
}

fn source_double_bond_stereo_is_unsupported(
    molecule: &Molecule,
    ring_membership: &RingMembership,
    bond_id: BondId,
    bond: &Bond,
) -> bool {
    bond.order != BondOrder::Double
        || graph_bond_in_ring_smaller_than(molecule, bond_id, 8)
        || (ring_membership.bond_in_ring(bond_id)
            && [bond.a(), bond.b()].into_iter().any(|atom_id| {
                molecule
                    .atom(atom_id)
                    .map(|atom| atom.element.symbol() != "C")
                    .unwrap_or(true)
            }))
}

fn source_double_bond_endpoint_carriers(
    molecule: &Molecule,
    endpoint: AtomId,
    other_endpoint: AtomId,
    focus_bond: BondId,
) -> Vec<StereoCarrier> {
    let mut carriers = Vec::new();
    if let Ok(incident) = molecule.incident_bonds(endpoint) {
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
    // A virtual endpoint carrier is accepted only when the represented atom
    // explicitly declares exactly one hydrogen.
    if molecule
        .atom(endpoint)
        .is_ok_and(|atom| atom.hydrogens.explicit_count() == 1)
    {
        carriers.push(StereoCarrier::ImplicitHydrogen);
    }
    carriers
}

fn select_directional_mark(
    molecule: &Molecule,
    endpoint: AtomId,
    other_endpoint: AtomId,
    focus_bond: BondId,
    marks: &[EndpointMark],
    issues: &mut Vec<SourceStereoNormalizationIssue>,
) -> Option<EndpointMark> {
    match marks {
        [] => None,
        [mark] => Some(*mark),
        [first, second]
            if redundant_endpoint_directional_marks(
                molecule,
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
            issues.push(
                SourceStereoNormalizationIssue::AmbiguousDirectionalBondMarks {
                    double_bond: focus_bond,
                    endpoint,
                    mark_count: marks.len(),
                },
            );
            None
        }
    }
}

fn redundant_endpoint_directional_marks(
    molecule: &Molecule,
    endpoint: AtomId,
    other_endpoint: AtomId,
    focus_bond: BondId,
    first: &EndpointMark,
    second: &EndpointMark,
) -> bool {
    if first.direction == second.direction {
        return false;
    }
    let mut marked_carriers = [first.carrier, second.carrier];
    marked_carriers.sort_unstable();
    let mut atom_carriers =
        source_double_bond_endpoint_carriers(molecule, endpoint, other_endpoint, focus_bond)
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
    source_marks: &[SourceStereoBondMark],
    used_marks: &[BondId],
    issues: &mut Vec<SourceStereoNormalizationIssue>,
) {
    for mark in source_marks_in_bond_order(source_marks) {
        match mark.kind {
            SourceStereoBondMarkKind::DirectionalUp | SourceStereoBondMarkKind::DirectionalDown => {
                if !used_marks.contains(&mark.bond) {
                    issues.push(
                        SourceStereoNormalizationIssue::UnpairedDirectionalBondMark {
                            bond: mark.bond,
                        },
                    );
                }
            }
            SourceStereoBondMarkKind::WedgeUp
            | SourceStereoBondMarkKind::WedgeDown
            | SourceStereoBondMarkKind::WedgeEither => {
                if !used_marks.contains(&mark.bond) {
                    issues.push(
                        SourceStereoNormalizationIssue::UnassembledTetrahedralBondMark {
                            bond: mark.bond,
                            kind: mark.kind,
                        },
                    );
                }
            }
            SourceStereoBondMarkKind::DoubleBondEither => {
                if !used_marks.contains(&mark.bond) {
                    issues.push(SourceStereoNormalizationIssue::UnsupportedSourceBondMark {
                        bond: mark.bond,
                        kind: mark.kind,
                    });
                }
            }
        }
    }
}

#[derive(Clone, Copy)]
struct EndpointMark {
    bond: BondId,
    carrier: AtomId,
    direction: SourceStereoBondMarkKind,
}

fn directional_marks_for_endpoint(
    molecule: &Molecule,
    source_marks: &[SourceStereoBondMark],
    endpoint: AtomId,
    focus_bond: BondId,
) -> Vec<EndpointMark> {
    let mut marks = Vec::new();
    let Ok(incident) = molecule.incident_bonds(endpoint) else {
        return marks;
    };
    for (bond_id, bond) in incident {
        if bond_id == focus_bond || bond.order != BondOrder::Single {
            continue;
        }
        let Some(mark) = source_marks.iter().find(|mark| mark.bond == bond_id) else {
            continue;
        };
        if matches!(
            mark.kind,
            SourceStereoBondMarkKind::DirectionalUp | SourceStereoBondMarkKind::DirectionalDown
        ) {
            marks.push(EndpointMark {
                bond: bond_id,
                carrier: bond.other_atom(endpoint),
                direction: directional_mark_at_endpoint(mark.kind, mark.from, endpoint),
            });
        }
    }
    marks.sort_by_key(|mark| (mark.bond, mark.carrier));
    marks
}

fn directional_mark_at_endpoint(
    kind: SourceStereoBondMarkKind,
    from: AtomId,
    endpoint: AtomId,
) -> SourceStereoBondMarkKind {
    if from == endpoint {
        kind
    } else {
        invert_directional_mark(kind)
    }
}

fn invert_directional_mark(kind: SourceStereoBondMarkKind) -> SourceStereoBondMarkKind {
    match kind {
        SourceStereoBondMarkKind::DirectionalUp => SourceStereoBondMarkKind::DirectionalDown,
        SourceStereoBondMarkKind::DirectionalDown => SourceStereoBondMarkKind::DirectionalUp,
        _ => kind,
    }
}

fn has_double_bond_element(molecule: &Molecule, bond: BondId) -> bool {
    molecule.stereo_elements().any(|(_, element)| {
        matches!(
            &element.kind,
            StereoElementKind::DoubleBond(stereo) if stereo.bond == bond
        )
    })
}

fn source_marks_in_bond_order(source_marks: &[SourceStereoBondMark]) -> Vec<&SourceStereoBondMark> {
    let mut marks = source_marks.iter().collect::<Vec<_>>();
    marks.sort_by_key(|mark| mark.bond);
    marks
}

fn has_tetrahedral_element(molecule: &Molecule, center: AtomId) -> bool {
    molecule.stereo_elements().any(|(_, element)| {
        matches!(
            &element.kind,
            StereoElementKind::Tetrahedral(stereo) if stereo.center == center
        )
    })
}

fn has_axis_element(molecule: &Molecule, axis: BondId) -> bool {
    molecule.stereo_elements().any(|(_, element)| {
        matches!(
            &element.kind,
            StereoElementKind::Axis(stereo) if stereo.axis == axis
        )
    })
}

fn tetrahedral_points(
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

fn coordinates_are_planar(points: &[Point3]) -> bool {
    let Some(first) = points.first() else {
        return false;
    };
    points
        .iter()
        .all(|point| (point.z - first.z).abs() <= COORDINATE_EPSILON)
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

fn carrier_key(carrier: &StereoCarrier) -> (u8, u32) {
    match carrier {
        StereoCarrier::Atom(atom) => (0, atom.raw()),
        StereoCarrier::ImplicitHydrogen => (1, u32::MAX),
        StereoCarrier::ImplicitLonePair => (2, u32::MAX),
    }
}

const COORDINATE_EPSILON: f64 = 1.0e-8;
