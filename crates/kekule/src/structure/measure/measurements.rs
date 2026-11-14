//! Distance, angle, dihedral, and contact measurements over structure views,
//! with the vector kernels that geometry edits share.

use crate::diagnostics::TopologyMismatch;
use crate::math::{Point3, Vector3};
use crate::topology::{AtomSelection, InstanceAtomId};
use crate::units::{Quantity, CANONICAL_ANGLE_UNIT, CANONICAL_LENGTH_UNIT};

use crate::structure::{CoordinateError, StructureView};

/// Optional topology validation for Cartesian measurements.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum ConnectivityCheck {
    /// Allow arbitrary atoms, without requiring bonds. This is the default.
    #[default]
    Unrestricted,
    /// Require A-B for distance, A-B and B-C for angle, and A-B, B-C, C-D
    /// for dihedral. Bond order is irrelevant. Different molecule instances
    /// and repeated consecutive atoms cannot satisfy this requirement.
    /// Nonconsecutive atoms need not be bonded or distinct.
    ///
    /// All atom IDs are validated before connectivity. The first missing pair
    /// in argument order is reported before any geometric calculation.
    ConsecutiveBonds,
}

fn point(view: StructureView<'_>, atom: InstanceAtomId) -> Result<Point3, CoordinateError> {
    let index = view
        .topology()
        .atom_index(atom)
        .ok_or(CoordinateError::InvalidAtomId(atom))?;
    Ok(view.positions().values().value()[index.index()])
}

fn check_connectivity(
    view: StructureView<'_>,
    atoms: &[InstanceAtomId],
    connectivity: ConnectivityCheck,
) -> Result<(), CoordinateError> {
    if connectivity == ConnectivityCheck::Unrestricted {
        return Ok(());
    }
    let topology = view.topology();
    for &atom in atoms {
        if topology.atom_index(atom).is_none() {
            return Err(CoordinateError::InvalidAtomId(atom));
        }
    }
    for pair in atoms.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        let bonded = topology
            .atom(a)
            .ok_or(CoordinateError::InvalidAtomId(a))?
            .neighbors()
            .any(|neighbor| neighbor.id() == b);
        if !bonded {
            return Err(CoordinateError::MissingBond { a, b });
        }
    }
    Ok(())
}

pub(in crate::structure) fn norm(vector: Vector3) -> Result<f64, CoordinateError> {
    let value = vector.x.hypot(vector.y).hypot(vector.z);
    if !value.is_finite() {
        return Err(CoordinateError::NumericalFailure);
    }
    Ok(value)
}

pub(in crate::structure) fn normalized(vector: Vector3) -> Result<Vector3, CoordinateError> {
    let length = norm(vector)?;
    if length == 0.0 {
        return Err(CoordinateError::DegenerateGeometry);
    }
    // Divide components directly: forming 1 / length first can overflow for
    // finite subnormal vectors, even though the normalized vector is finite.
    Ok(Vector3::new(
        vector.x / length,
        vector.y / length,
        vector.z / length,
    ))
}

/// Measures the stored Cartesian distance. Coincident atoms have distance zero.
///
/// `connectivity` selects whether a bond A-B is required. Connectivity
/// validation precedes geometric calculation; cells are ignored.
pub fn distance(
    view: StructureView<'_>,
    a: InstanceAtomId,
    b: InstanceAtomId,
    connectivity: ConnectivityCheck,
) -> Result<Quantity<f64>, CoordinateError> {
    check_connectivity(view, &[a, b], connectivity)?;
    Ok(Quantity::new(
        norm(point(view, b)? - point(view, a)?)?,
        CANONICAL_LENGTH_UNIT,
    ))
}

/// Measures the angle A-B-C, with B as vertex, in [0, pi] radians.
/// A zero-length arm is an error; a straight angle is valid.
///
/// `connectivity` selects whether bonds A-B and B-C are required.
/// Connectivity validation precedes geometric calculation; cells are ignored.
pub fn angle(
    view: StructureView<'_>,
    a: InstanceAtomId,
    b: InstanceAtomId,
    c: InstanceAtomId,
    connectivity: ConnectivityCheck,
) -> Result<Quantity<f64>, CoordinateError> {
    check_connectivity(view, &[a, b, c], connectivity)?;
    Ok(Quantity::new(
        angle_points(point(view, a)?, point(view, b)?, point(view, c)?)?,
        CANONICAL_ANGLE_UNIT,
    ))
}

pub(in crate::structure) fn angle_points(
    a: Point3,
    b: Point3,
    c: Point3,
) -> Result<f64, CoordinateError> {
    let u = normalized(a - b)?;
    let v = normalized(c - b)?;
    Ok(norm(u.cross(v))?.atan2(u.dot(v)))
}

/// Measures the signed A-B-C-D dihedral in [-pi, pi] radians.
/// With consecutive bond vectors u, v, w, the sign follows
/// `atan2(((u x v) x (v x w)) . v_hat, (u x v) . (v x w))`.
/// Zero-length bonds and collinear defining triples are errors.
///
/// `connectivity` selects whether bonds A-B, B-C, and C-D are required.
/// Connectivity validation precedes geometric calculation; cells are ignored.
pub fn dihedral(
    view: StructureView<'_>,
    a: InstanceAtomId,
    b: InstanceAtomId,
    c: InstanceAtomId,
    d: InstanceAtomId,
    connectivity: ConnectivityCheck,
) -> Result<Quantity<f64>, CoordinateError> {
    check_connectivity(view, &[a, b, c, d], connectivity)?;
    let [a, b, c, d] = [
        point(view, a)?,
        point(view, b)?,
        point(view, c)?,
        point(view, d)?,
    ];
    Ok(Quantity::new(
        dihedral_points(a, b, c, d)?,
        CANONICAL_ANGLE_UNIT,
    ))
}

pub(in crate::structure) fn dihedral_points(
    a: Point3,
    b: Point3,
    c: Point3,
    d: Point3,
) -> Result<f64, CoordinateError> {
    let u = normalized(b - a)?;
    let v = normalized(c - b)?;
    let w = normalized(d - c)?;
    let n0 = normalized(u.cross(v))?;
    let n1 = normalized(v.cross(w))?;
    Ok(n0.cross(n1).dot(v).atan2(n0.dot(n1)))
}

/// Selects candidate atoms at Cartesian distance **<= cutoff** from any
/// reference atom. Both selections must share the view's topology layout.
/// Overlapping candidate/reference atoms match even at zero cutoff. Empty
/// references produce an empty selection. The result is a static atom set for
/// this view; call again for each frame when membership should change.
/// Comparisons use canonical floating-point values without an added tolerance.
///
/// Whole-residue expansion is separate: call [`AtomSelection::expand_to_residues`]
/// on the result. Candidate restrictions apply before expansion.
pub fn within(
    view: StructureView<'_>,
    candidates: &AtomSelection,
    reference: &AtomSelection,
    cutoff: Quantity<f64>,
) -> Result<AtomSelection, CoordinateError> {
    let topology = view.shared_topology();
    for selection in [candidates, reference] {
        selection
            .ensure_compatible(&topology)
            .map_err(|_| CoordinateError::TopologyMismatch(TopologyMismatch::Selection))?;
    }
    let cutoff = cutoff
        .value_in(CANONICAL_LENGTH_UNIT)
        .map_err(CoordinateError::Unit)?;
    if !cutoff.is_finite() || cutoff < 0.0 {
        return Err(CoordinateError::InvalidCutoff);
    }
    let points = reference
        .atom_ids()
        .map(|atom| point(view, atom))
        .collect::<Result<Vec<_>, _>>()?;
    let mut selected = Vec::new();
    if !points.is_empty() {
        // Candidates are measured up to the first missing one, which then
        // fails.
        let mut measured = Vec::new();
        let mut missing = None;
        for atom in candidates.atom_ids() {
            match point(view, atom) {
                Ok(p) => measured.push((atom, p)),
                Err(error) => {
                    missing = Some(error);
                    break;
                }
            }
        }
        // Coordinates this small keep every distance finite, so measuring
        // only the pairs a grid finds near selects the same atoms and fails
        // alike; few pairs are measured directly.
        let bounded = |p: &Point3| p.x.abs() <= 1e300 && p.y.abs() <= 1e300 && p.z.abs() <= 1e300;
        if measured.len().saturating_mul(points.len()) > 1024
            && points
                .iter()
                .chain(measured.iter().map(|(_, p)| p))
                .all(bounded)
        {
            let mut grid = crate::math::PointGrid::default();
            grid.build(
                points.iter().map(|&q| (q, ())),
                crate::math::bounding_box(points.iter().copied()),
                cutoff,
                measured.len(),
            );
            for &(atom, p) in &measured {
                if grid.any_near(&p, |q, ()| norm(p - *q).is_ok_and(|d| d <= cutoff)) {
                    selected.push(atom);
                }
            }
        } else {
            for &(atom, p) in &measured {
                for q in &points {
                    if norm(p - *q)? <= cutoff {
                        selected.push(atom);
                        break;
                    }
                }
            }
        }
        if let Some(error) = missing {
            return Err(error);
        }
    }
    Ok(AtomSelection::from_atoms(topology.clone(), selected)
        .expect("selected atoms are atoms of the view's topology"))
}
