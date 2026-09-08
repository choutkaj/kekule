//! Checked measurements of stored Cartesian coordinates through [`ModelView`].
//!
//! Cells are ignored: these operations never wrap, image, unwrap, or align
//! coordinates. Apply the desired preprocessing before measuring. Angles use
//! radians and distances use the canonical length unit; results support ordinary
//! [`Quantity`] conversions. No structural ownership is materialized.
//!
//! A whole-residue pocket around one source-identified ligand:
//! ```
//! use kekule::{
//!     core::Element,
//!     structure::{measure, Model},
//!     topology::AtomSelection,
//!     units::{Quantity, ANGSTROM},
//! };
//! # fn pocket(model: &Model) -> Result<Model, Box<dyn std::error::Error>> {
//! let topology = model.shared_topology();
//! let ligand = topology.residue_by_author("A", "215", None)?;
//! let reference = AtomSelection::for_residues(&topology, [ligand.id()])?;
//! let hydrogen = AtomSelection::for_elements(
//!     &topology, [Element::from_symbol("H").unwrap()],
//! )?;
//! let candidates = AtomSelection::all(&topology).difference(&hydrogen)?;
//! let atoms = measure::within(
//!     model.view(), &candidates, &reference, Quantity::new(5.0, ANGSTROM),
//! )?.expand_to_residues();
//! let pocket = model.slice(&atoms)?;
//! # Ok(pocket)
//! # }
//! ```

use std::fmt;

use crate::geometry::{Point3, Vector3};
use crate::topology::{AtomSelection, InstanceAtomId, SelectionError};
use crate::units::{Quantity, UnitError, CANONICAL_ANGLE_UNIT, CANONICAL_LENGTH_UNIT};

use super::ModelView;

/// A failed Cartesian measurement or spatial selection.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum MeasurementError {
    InvalidAtomId(InstanceAtomId),
    Selection(SelectionError),
    Unit(UnitError),
    InvalidCutoff,
    DegenerateGeometry,
    NumericalFailure,
}

impl fmt::Display for MeasurementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidAtomId(atom) => write!(f, "invalid measurement atom: {atom}"),
            Self::Selection(error) => write!(f, "measurement selection: {error}"),
            Self::Unit(error) => write!(f, "measurement unit: {error}"),
            Self::InvalidCutoff => f.write_str("distance cutoff must be finite and nonnegative"),
            Self::DegenerateGeometry => {
                f.write_str("angle or dihedral is undefined for this geometry")
            }
            Self::NumericalFailure => {
                f.write_str("Cartesian measurement exceeds finite numerical range")
            }
        }
    }
}

impl std::error::Error for MeasurementError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Selection(error) => Some(error),
            Self::Unit(error) => Some(error),
            _ => None,
        }
    }
}

fn point(view: ModelView<'_>, atom: InstanceAtomId) -> Result<Point3, MeasurementError> {
    let index = view
        .topology()
        .atom_index(atom)
        .ok_or(MeasurementError::InvalidAtomId(atom))?;
    Ok(view.positions().values().value()[index.index()])
}

fn norm(vector: Vector3) -> Result<f64, MeasurementError> {
    let value = vector.x.hypot(vector.y).hypot(vector.z);
    if !value.is_finite() {
        return Err(MeasurementError::NumericalFailure);
    }
    Ok(value)
}

fn normalized(vector: Vector3) -> Result<Vector3, MeasurementError> {
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
pub fn distance(
    view: ModelView<'_>,
    a: InstanceAtomId,
    b: InstanceAtomId,
) -> Result<Quantity<f64>, MeasurementError> {
    Ok(Quantity::new(
        norm(point(view, b)? - point(view, a)?)?,
        CANONICAL_LENGTH_UNIT,
    ))
}

/// Measures the angle A-B-C, with B as vertex, in [0, pi] radians.
/// A zero-length arm is an error; a straight angle is valid.
pub fn angle(
    view: ModelView<'_>,
    a: InstanceAtomId,
    b: InstanceAtomId,
    c: InstanceAtomId,
) -> Result<Quantity<f64>, MeasurementError> {
    let b = point(view, b)?;
    let u = normalized(point(view, a)? - b)?;
    let v = normalized(point(view, c)? - b)?;
    Ok(Quantity::new(
        norm(u.cross(v))?.atan2(u.dot(v)),
        CANONICAL_ANGLE_UNIT,
    ))
}

/// Measures the signed A-B-C-D dihedral in [-pi, pi] radians.
/// With consecutive bond vectors u, v, w, the sign follows
/// `atan2(((u x v) x (v x w)) . v_hat, (u x v) . (v x w))`.
/// Zero-length bonds and collinear defining triples are errors.
pub fn dihedral(
    view: ModelView<'_>,
    a: InstanceAtomId,
    b: InstanceAtomId,
    c: InstanceAtomId,
    d: InstanceAtomId,
) -> Result<Quantity<f64>, MeasurementError> {
    let [a, b, c, d] = [
        point(view, a)?,
        point(view, b)?,
        point(view, c)?,
        point(view, d)?,
    ];
    let u = normalized(b - a)?;
    let v = normalized(c - b)?;
    let w = normalized(d - c)?;
    let n0 = normalized(u.cross(v))?;
    let n1 = normalized(v.cross(w))?;
    Ok(Quantity::new(
        n0.cross(n1).dot(v).atan2(n0.dot(n1)),
        CANONICAL_ANGLE_UNIT,
    ))
}

/// Selects candidate atoms at Cartesian distance **<= cutoff** from any
/// reference atom. Both selections must belong to the view's exact topology.
/// Overlapping candidate/reference atoms match even at zero cutoff. Empty
/// references produce an empty selection. The result is a static atom set for
/// this view; call again for each frame when membership should change.
/// Comparisons use canonical floating-point values without an added tolerance.
///
/// Whole-residue expansion is separate: call [`AtomSelection::expand_to_residues`]
/// on the result. Candidate restrictions apply before expansion.
pub fn within(
    view: ModelView<'_>,
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
        for atom in candidates.atom_ids() {
            let p = point(view, atom)?;
            for q in &points {
                if norm(p - *q)? <= cutoff {
                    selected.push(atom);
                    break;
                }
            }
        }
    }
    AtomSelection::from_atoms(&topology, selected).map_err(MeasurementError::Selection)
}
