//! [`GeometryError`]: invalid geometry state.

use crate::diagnostics::LengthMismatch;
use std::fmt;

use crate::math::PeriodicCellError;
use crate::properties::PropertyError;
use crate::structure::DenseArrayError;
use crate::units::UnitError;

/// Invalid geometry state.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum GeometryError {
    /// Dense per-atom state is invalid.
    Array(DenseArrayError),
    /// Bound bond property rows have the wrong length.
    BondCountMismatch(LengthMismatch),
    /// The geometry exceeds addressable capacity.
    CapacityOverflow,
    /// A geometry property edit failed; see the source.
    Property(Box<PropertyError>),
    /// The periodic cell is invalid; see the source.
    Cell(Box<PeriodicCellError>),
    /// A quantity could not be converted to the required unit; see the source.
    Unit(Box<UnitError>),
}

impl fmt::Display for GeometryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Array(error) => error.fmt(formatter),
            Self::BondCountMismatch(error) => error.fmt(formatter),
            Self::CapacityOverflow => {
                formatter.write_str("The geometry exceeds the addressable capacity.")
            }
            Self::Property(error) => error.fmt(formatter),
            Self::Cell(error) => error.fmt(formatter),
            Self::Unit(error) => error.fmt(formatter),
        }
    }
}

impl std::error::Error for GeometryError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Array(error) => error.source(),
            Self::Property(error) => error.source(),
            Self::Cell(error) => error.source(),
            Self::Unit(error) => error.source(),
            _ => None,
        }
    }
}

impl From<DenseArrayError> for GeometryError {
    fn from(error: DenseArrayError) -> Self {
        Self::Array(error)
    }
}

impl From<PropertyError> for GeometryError {
    fn from(error: PropertyError) -> Self {
        Self::Property(Box::new(error))
    }
}

impl From<PeriodicCellError> for GeometryError {
    fn from(error: PeriodicCellError) -> Self {
        Self::Cell(Box::new(error))
    }
}

impl From<UnitError> for GeometryError {
    fn from(error: UnitError) -> Self {
        Self::Unit(Box::new(error))
    }
}
