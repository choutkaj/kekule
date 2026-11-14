use super::{center, SolvationError};
use crate::structure::Structure;
use crate::{
    math::{Lattice, PeriodicCell, Vector3},
    units::{Quantity, CANONICAL_LENGTH_UNIT},
};

/// The cell shape for [`PeriodicBoxOptions::Padding`]. Every shape has the same
/// minimum image distance, the computed width.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum BoxShape {
    /// A cubic cell, the default.
    #[default]
    Cube,
    /// A rhombic dodecahedron with a square xy face, about 71% of the cube volume.
    Dodecahedron,
    /// A truncated octahedron, about 77% of the cube volume.
    Octahedron,
}
/// Exactly one sizing mode, so conflicting sizing inputs are unrepresentable.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum PeriodicBoxOptions {
    /// An orthorhombic box with these edge lengths, in any length unit.
    Dimensions(Quantity<Vector3>),
    /// An explicit cell, which must be periodic along all three axes.
    Cell(PeriodicCell),
    /// A box of `shape` sized from the solute; see [`add_periodic_box`].
    Padding {
        /// Minimum separation between periodic solute copies; finite, positive, any
        /// length unit.
        padding: Quantity<f64>,
        /// The cell shape.
        shape: BoxShape,
    },
}
/// Sets or replaces a fully periodic box without moving any atom or replacing topology.
///
/// Padding uses a sphere centered at the Cartesian bounding-box midpoint. For
/// diameter D and padding p, width is max(D + p, 2p). The padding is the minimum
/// separation guaranteed between solute copies, not a margin on each box face.
/// Wrapped solutes should be made whole before padding-based sizing. Explicit
/// dimensions/vectors do not assert a solute clearance. Failure changes nothing.
///
/// # Errors
///
/// Returns [`SolvationError::InvalidPadding`] for a padding that is not finite
/// and positive, [`SolvationError::Unit`] when it is not a length,
/// [`SolvationError::Cell`] when the box vectors form no valid cell or
/// centering or the cell geometry fails numerically, and
/// [`SolvationError::NotFullyPeriodic`] for a [`PeriodicBoxOptions::Cell`]
/// that is not periodic along all three axes.
pub fn add_periodic_box(
    structure: &mut Structure,
    options: &PeriodicBoxOptions,
) -> Result<(), SolvationError> {
    let cell = match options {
        PeriodicBoxOptions::Dimensions(lengths) => PeriodicCell::orthorhombic(*lengths, [true; 3])?,
        PeriodicBoxOptions::Cell(cell) => *cell,
        PeriodicBoxOptions::Padding { padding, shape } => {
            let padding = padding.into_unit(CANONICAL_LENGTH_UNIT)?.into_value();
            if !padding.is_finite() || padding <= 0.0 {
                return Err(SolvationError::InvalidPadding);
            }
            let points = structure.positions().values();
            let center = center(points.value())?;
            let radius = points
                .value()
                .iter()
                .map(|p| (*p - center).norm())
                .fold(0.0_f64, f64::max);
            let width = (2.0 * radius + padding).max(2.0 * padding);
            let vectors = match shape {
                BoxShape::Cube => [
                    Vector3::new(1.0, 0.0, 0.0),
                    Vector3::new(0.0, 1.0, 0.0),
                    Vector3::new(0.0, 0.0, 1.0),
                ],
                BoxShape::Dodecahedron => [
                    Vector3::new(1.0, 0.0, 0.0),
                    Vector3::new(0.0, 1.0, 0.0),
                    Vector3::new(0.5, 0.5, 0.5 * 2.0_f64.sqrt()),
                ],
                BoxShape::Octahedron => [
                    Vector3::new(1.0, 0.0, 0.0),
                    Vector3::new(1.0 / 3.0, 2.0 * 2.0_f64.sqrt() / 3.0, 0.0),
                    Vector3::new(-1.0 / 3.0, 2.0_f64.sqrt() / 3.0, 6.0_f64.sqrt() / 3.0),
                ],
            }
            .map(|v| v * width);
            PeriodicCell::new(Quantity::new(vectors, CANONICAL_LENGTH_UNIT), [true; 3])?
        }
    };
    if cell.periodic_axes() != [true; 3] {
        return Err(SolvationError::NotFullyPeriodic);
    }
    Lattice::new(cell)?;
    structure.geometry_mut().set_cell(Some(cell));
    Ok(())
}
