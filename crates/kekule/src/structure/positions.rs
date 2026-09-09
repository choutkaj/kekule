use std::fmt;

use crate::geometry::Point3;
use crate::units::{Quantity, UnitError, CANONICAL_LENGTH_UNIT};

/// A dense numerical coordinate array in Kekule's canonical length unit.
///
/// `Positions` has no topology context. Structural owners such as [`super::Model`]
/// validate its length and translate semantic atom identifiers to dense indices.
/// Construction accepts any compatible length [`crate::units::Unit`] and
/// converts values once to the library-wide canonical unit. Non-finite values
/// are rejected.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Positions {
    values: Vec<Point3>,
}

impl Positions {
    pub(crate) fn into_canonical_values(self) -> Vec<Point3> {
        self.values
    }
    pub(crate) fn try_reserve(
        &mut self,
        additional: usize,
    ) -> Result<(), std::collections::TryReserveError> {
        self.values.try_reserve(additional)
    }
    pub(crate) fn extend_canonical(&mut self, positions: &Positions) {
        self.values.extend_from_slice(&positions.values);
    }
    pub(super) fn from_canonical_values(values: Vec<Point3>) -> Self {
        debug_assert!(values.iter().all(|point| point.is_finite()));
        Self { values }
    }

    /// Copies numerical coordinates into canonical storage.
    ///
    /// Use [`Self::from_vec`] to transfer an owned vector without allocating a
    /// second coordinate array.
    pub fn new<T>(positions: Quantity<T>) -> Result<Self, PositionError>
    where
        T: AsRef<[Point3]>,
    {
        let factor = positions
            .unit()
            .conversion_factor_to(CANONICAL_LENGTH_UNIT)?;
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

    /// Consumes a coordinate vector, converting and validating it in place.
    ///
    /// The vector's allocation and capacity are retained. An incompatible unit
    /// or non-finite converted coordinate returns an error.
    pub fn from_vec(positions: Quantity<Vec<Point3>>) -> Result<Self, PositionError> {
        let factor = positions
            .unit()
            .conversion_factor_to(CANONICAL_LENGTH_UNIT)?;
        let mut values = positions.into_value();
        for (index, point) in values.iter_mut().enumerate() {
            *point = Point3::new(point.x * factor, point.y * factor, point.z * factor);
            if !point.is_finite() {
                return Err(PositionError::NonFinitePosition { index });
            }
        }
        Ok(Self { values })
    }

    /// Transfers the coordinate vector in its canonical length unit.
    pub fn into_values(self) -> Quantity<Vec<Point3>> {
        Quantity::new(self.values, CANONICAL_LENGTH_UNIT)
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
        Quantity::new(self.values.as_slice(), CANONICAL_LENGTH_UNIT)
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
            .map(|point| Quantity::new(point, CANONICAL_LENGTH_UNIT))
            .ok_or(PositionError::InvalidIndex { index })
    }

    pub fn set_position_at(
        &mut self,
        index: usize,
        position: Quantity<Point3>,
    ) -> Result<(), PositionError> {
        let point = position.into_unit(CANONICAL_LENGTH_UNIT)?.into_value();
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
        // AsRef is user-implementable and can return a different slice on each
        // call. Validate and copy the same borrow, not two separate conversions.
        let source = positions.value().as_ref();
        let factor = self.validate_replacement(&Quantity::new(source, positions.unit()))?;
        self.copy_from_validated(source, factor);
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
        let factor = positions
            .unit()
            .conversion_factor_to(CANONICAL_LENGTH_UNIT)?;
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

impl std::error::Error for PositionError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Unit(error) => Some(error),
            _ => None,
        }
    }
}

impl From<UnitError> for PositionError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}
