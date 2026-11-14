//! Distance, angle, dihedral, and contact measurements over structure views,
//! with the vector kernels that geometry edits share.

use std::fmt;

use crate::diagnostics::Capitalized;
use crate::math::{Point3, Vector3};
use crate::topology::{AtomSelection, InstanceAtomId, SelectionError};
use crate::units::{Quantity, UnitError, CANONICAL_ANGLE_UNIT, CANONICAL_LENGTH_UNIT};

use crate::structure::StructureView;

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

/// A failed Cartesian measurement or spatial selection.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum MeasurementError {
    /// A measurement atom is not in the view's topology.
    InvalidAtomId(InstanceAtomId),
    /// Consecutive measurement atoms lack a direct bond, in argument order.
    MissingBond {
        /// The earlier atom of the pair.
        a: InstanceAtomId,
        /// The later atom of the pair.
        b: InstanceAtomId,
    },
    /// A selection does not match the view's topology; see the source.
    Selection(SelectionError),
    /// The cutoff could not be converted to the canonical length unit; see
    /// the source.
    Unit(UnitError),
    /// The cutoff is negative or not finite.
    InvalidCutoff,
    /// A direction needed by the measurement has zero length.
    DegenerateGeometry,
    /// An intermediate length is not finite.
    NumericalFailure,
}

impl fmt::Display for MeasurementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(atom) => write!(f, "{} is not in the topology.", Capitalized(atom)),
            Self::MissingBond { a, b } => write!(f, "{} and {b} are not bonded.", Capitalized(a)),
            Self::Selection(error) => error.fmt(f),
            Self::Unit(error) => error.fmt(f),
            Self::InvalidCutoff => f.write_str("The distance cutoff is negative or not finite."),
            Self::DegenerateGeometry => {
                f.write_str("The atoms define no angle or dihedral: a direction has zero length.")
            }
            Self::NumericalFailure => {
                f.write_str("The measurement exceeds the finite numerical range.")
            }
        }
    }
}

impl std::error::Error for MeasurementError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Selection(error) => error.source(),
            Self::Unit(error) => error.source(),
            _ => None,
        }
    }
}

fn point(view: StructureView<'_>, atom: InstanceAtomId) -> Result<Point3, MeasurementError> {
    let index = view
        .topology()
        .atom_index(atom)
        .ok_or(MeasurementError::InvalidAtomId(atom))?;
    Ok(view.positions().values().value()[index.index()])
}

fn check_connectivity(
    view: StructureView<'_>,
    atoms: &[InstanceAtomId],
    connectivity: ConnectivityCheck,
) -> Result<(), MeasurementError> {
    if connectivity == ConnectivityCheck::Unrestricted {
        return Ok(());
    }
    let topology = view.topology();
    for &atom in atoms {
        if topology.atom_index(atom).is_none() {
            return Err(MeasurementError::InvalidAtomId(atom));
        }
    }
    for pair in atoms.windows(2) {
        let (a, b) = (pair[0], pair[1]);
        let bonded = topology
            .atom(a)
            .ok_or(MeasurementError::InvalidAtomId(a))?
            .neighbors()
            .any(|neighbor| neighbor.id() == b);
        if !bonded {
            return Err(MeasurementError::MissingBond { a, b });
        }
    }
    Ok(())
}

pub(in crate::structure) fn norm(vector: Vector3) -> Result<f64, MeasurementError> {
    let value = vector.x.hypot(vector.y).hypot(vector.z);
    if !value.is_finite() {
        return Err(MeasurementError::NumericalFailure);
    }
    Ok(value)
}

pub(in crate::structure) fn normalized(vector: Vector3) -> Result<Vector3, MeasurementError> {
    let length = norm(vector)?;
    if length == 0.0 {
        return Err(MeasurementError::DegenerateGeometry);
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
) -> Result<Quantity<f64>, MeasurementError> {
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
) -> Result<Quantity<f64>, MeasurementError> {
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
) -> Result<f64, MeasurementError> {
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
) -> Result<Quantity<f64>, MeasurementError> {
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
) -> Result<f64, MeasurementError> {
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
) -> Result<AtomSelection, MeasurementError> {
    let topology = view.shared_topology();
    candidates
        .ensure_compatible(&topology)
        .map_err(MeasurementError::Selection)?;
    reference
        .ensure_compatible(&topology)
        .map_err(MeasurementError::Selection)?;
    let cutoff = cutoff
        .value_in(CANONICAL_LENGTH_UNIT)
        .map_err(MeasurementError::Unit)?;
    if !cutoff.is_finite() || cutoff < 0.0 {
        return Err(MeasurementError::InvalidCutoff);
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
    AtomSelection::from_atoms(topology.clone(), selected).map_err(MeasurementError::Selection)
}
