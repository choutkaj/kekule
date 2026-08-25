use std::fmt;

use crate::geometry::Point3;
use crate::units::{Quantity, UnitError, MODEL_LENGTH_UNIT};

/// A dense numerical coordinate array in canonical model length units.
///
/// `Positions` has no topology context. Structural owners such as [`super::Model`]
/// validate its length and translate semantic atom identifiers to dense indices.
#[derive(Debug, Clone, PartialEq)]
pub struct Positions {
    values: Vec<Point3>,
}

impl Positions {
    pub(super) fn from_model_values(values: Vec<Point3>) -> Self {
        debug_assert!(values.iter().all(|point| point.is_finite()));
        Self { values }
    }

    /// Constructs positions from numerical coordinates alone.
    pub fn new<T>(positions: Quantity<T>) -> Result<Self, PositionError>
    where
        T: AsRef<[Point3]>,
    {
        let factor = positions.unit().conversion_factor_to(MODEL_LENGTH_UNIT)?;
        let source = positions.value().as_ref();
        let values = source
            .iter()
            .copied()
            .enumerate()
            .map(|(index, point)| {
                let point = Point3::new(point.x * factor, point.y * factor, point.z * factor);
                if !point.is_finite() {
                    return Err(PositionError::NonFinitePosition { index });
                }
                Ok(point)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self { values })
    }

    /// Constructs a zero-filled coordinate array with `len` entries.
    pub fn zeros(len: usize) -> Self {
        Self {
            values: vec![Point3::origin(); len],
        }
    }

    pub fn len(&self) -> usize {
        self.values.len()
    }

    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    pub fn values(&self) -> Quantity<&[Point3]> {
        Quantity::new(self.values.as_slice(), MODEL_LENGTH_UNIT)
    }

    /// Copies a deterministic dense projection in the requested index order.
    pub fn select_indices(&self, indices: &[usize]) -> Result<Self, PositionError> {
        let values = indices
            .iter()
            .map(|index| {
                self.values
                    .get(*index)
                    .copied()
                    .ok_or(PositionError::InvalidIndex { index: *index })
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self { values })
    }

    pub fn position_at(&self, index: usize) -> Result<Quantity<Point3>, PositionError> {
        self.values
            .get(index)
            .copied()
            .map(|point| Quantity::new(point, MODEL_LENGTH_UNIT))
            .ok_or(PositionError::InvalidIndex { index })
    }

    pub fn set_position_at(
        &mut self,
        index: usize,
        position: Quantity<Point3>,
    ) -> Result<(), PositionError> {
        let point = position.to_unit(MODEL_LENGTH_UNIT)?.to_value();
        if !point.is_finite() {
            return Err(PositionError::NonFinitePosition { index });
        }
        let destination = self
            .values
            .get_mut(index)
            .ok_or(PositionError::InvalidIndex { index })?;
        *destination = point;
        Ok(())
    }

    /// Replaces all positions transactionally while reusing the current
    /// allocation when the capacity permits.
    pub fn set_all<T>(&mut self, positions: Quantity<T>) -> Result<(), PositionError>
    where
        T: AsRef<[Point3]>,
    {
        let factor = self.validate_replacement(&positions)?;
        self.copy_from_validated(positions.value().as_ref(), factor);
        Ok(())
    }

    /// Validates a complete replacement without changing this array.
    pub fn validate_all<T>(&self, positions: &Quantity<T>) -> Result<(), PositionError>
    where
        T: AsRef<[Point3]>,
    {
        self.validate_replacement(positions).map(drop)
    }

    fn validate_replacement<T>(&self, positions: &Quantity<T>) -> Result<f64, PositionError>
    where
        T: AsRef<[Point3]>,
    {
        let factor = positions.unit().conversion_factor_to(MODEL_LENGTH_UNIT)?;
        let source = positions.value().as_ref();
        if source.len() != self.len() {
            return Err(PositionError::PositionCountMismatch {
                expected: self.len(),
                actual: source.len(),
            });
        }
        for (index, point) in source.iter().copied().enumerate() {
            let converted = Point3::new(point.x * factor, point.y * factor, point.z * factor);
            if !converted.is_finite() {
                return Err(PositionError::NonFinitePosition { index });
            }
        }
        Ok(factor)
    }

    pub(crate) fn copy_from_validated(&mut self, source: &[Point3], factor: f64) {
        for (destination, source) in self.values.iter_mut().zip(source.iter().copied()) {
            *destination = Point3::new(source.x * factor, source.y * factor, source.z * factor);
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PositionError {
    InvalidIndex { index: usize },
    PositionCountMismatch { expected: usize, actual: usize },
    NonFinitePosition { index: usize },
    Unit(UnitError),
}

impl fmt::Display for PositionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidIndex { index } => write!(formatter, "invalid position index {index}"),
            Self::PositionCountMismatch { expected, actual } => write!(
                formatter,
                "positions require {expected} coordinates, but received {actual}"
            ),
            Self::NonFinitePosition { index } => {
                write!(formatter, "position at {index} is not finite")
            }
            Self::Unit(error) => write!(formatter, "invalid position unit: {error}"),
        }
    }
}

impl std::error::Error for PositionError {}

impl From<UnitError> for PositionError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}
