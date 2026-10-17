use super::center;
use crate::structure::Model;
use crate::{
    geometry::{PeriodicCell, PeriodicCellError, PeriodicGeometry, PeriodicGeometryError, Vector3},
    units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT},
};
use std::fmt;

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
/// A periodic box could not be built for the model.
#[derive(Debug)]
#[non_exhaustive]
pub enum PeriodicBoxError {
    /// The padding is not finite and positive.
    InvalidPadding,
    /// The supplied cell is not periodic along all three axes.
    NotFullyPeriodic,
    /// The padding could not be converted to the canonical length unit; see
    /// the source.
    Unit(UnitError),
    /// The box vectors do not form a valid cell; see the source.
    Cell(PeriodicCellError),
    /// Centering the solute or validating the cell geometry failed; see the
    /// source.
    Geometry(PeriodicGeometryError),
}
impl fmt::Display for PeriodicBoxError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPadding => f.write_str("padding must be finite and positive"),
            Self::NotFullyPeriodic => f.write_str("box construction requires three periodic axes"),
            Self::Unit(_) => f.write_str("invalid box unit"),
            Self::Cell(_) => f.write_str("invalid periodic cell"),
            Self::Geometry(_) => f.write_str("periodic geometry failed"),
        }
    }
}
impl std::error::Error for PeriodicBoxError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit(e) => Some(e),
            Self::Cell(e) => Some(e),
            Self::Geometry(e) => Some(e),
            _ => None,
        }
    }
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
/// Returns [`PeriodicBoxError::InvalidPadding`] for a padding that is not finite
/// and positive, [`PeriodicBoxError::Unit`] when it is not a length,
/// [`PeriodicBoxError::Cell`] when the box vectors form no valid cell,
/// [`PeriodicBoxError::NotFullyPeriodic`] for a [`PeriodicBoxOptions::Cell`]
/// that is not periodic along all three axes, and [`PeriodicBoxError::Geometry`]
/// when centering or the cell geometry fails numerically.
pub fn add_periodic_box(
    model: &mut Model,
    options: &PeriodicBoxOptions,
) -> Result<(), PeriodicBoxError> {
    let cell = match options {
        PeriodicBoxOptions::Dimensions(lengths) => {
            PeriodicCell::orthorhombic(*lengths, [true; 3]).map_err(PeriodicBoxError::Cell)?
        }
        PeriodicBoxOptions::Cell(cell) => *cell,
        PeriodicBoxOptions::Padding { padding, shape } => {
            let padding = padding
                .into_unit(CANONICAL_LENGTH_UNIT)
                .map_err(PeriodicBoxError::Unit)?
                .into_value();
            if !padding.is_finite() || padding <= 0.0 {
                return Err(PeriodicBoxError::InvalidPadding);
            }
            let points = model.positions().values();
            let center = center(points.value()).map_err(PeriodicBoxError::Geometry)?;
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
            PeriodicCell::new(Quantity::new(vectors, CANONICAL_LENGTH_UNIT), [true; 3])
                .map_err(PeriodicBoxError::Cell)?
        }
    };
    if cell.periodic_axes() != [true; 3] {
        return Err(PeriodicBoxError::NotFullyPeriodic);
    }
    PeriodicGeometry::new(cell).map_err(PeriodicBoxError::Geometry)?;
    model.set_cell(Some(cell));
    Ok(())
}
