use super::{center, Model};
use crate::{
    geometry::{PeriodicCell, PeriodicCellError, PeriodicGeometry, PeriodicGeometryError, Vector3},
    units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT},
};
use std::fmt;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum BoxShape {
    #[default]
    Cube,
    Dodecahedron,
    Octahedron,
}
/// Exactly one sizing mode, so conflicting sizing inputs are unrepresentable.
#[derive(Debug, Clone)]
pub enum PeriodicBoxOptions {
    Dimensions(Quantity<Vector3>),
    Cell(PeriodicCell),
    Padding {
        padding: Quantity<f64>,
        shape: BoxShape,
    },
}
#[derive(Debug)]
#[non_exhaustive]
pub enum PeriodicBoxError {
    InvalidPadding,
    NotFullyPeriodic,
    Unit(UnitError),
    Cell(PeriodicCellError),
    Geometry(PeriodicGeometryError),
}
impl fmt::Display for PeriodicBoxError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPadding => f.write_str("padding must be finite and positive"),
            Self::NotFullyPeriodic => f.write_str("box construction requires three periodic axes"),
            Self::Unit(e) => e.fmt(f),
            Self::Cell(e) => e.fmt(f),
            Self::Geometry(e) => e.fmt(f),
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
impl Model {
    /// Sets or replaces a fully periodic box without moving any atom or replacing topology.
    ///
    /// Padding uses a sphere centered at the Cartesian bounding-box midpoint. For
    /// diameter D and padding p, width is max(D + p, 2p). The padding is the minimum
    /// separation guaranteed between solute copies, not a margin on each box face.
    /// Wrapped solutes should be made whole before padding-based sizing. Explicit
    /// dimensions/vectors do not assert a solute clearance. Failure changes nothing.
    pub fn add_periodic_box(
        &mut self,
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
                let points = self.positions().values();
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
        self.set_cell(Some(cell));
        Ok(())
    }
}
