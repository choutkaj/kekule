use crate::geometry::Vector3;
use crate::units::{Quantity, Unit, CANONICAL_FORCE_UNIT, CANONICAL_VELOCITY_UNIT};

use super::ConformationError;

#[derive(Debug, Clone, PartialEq)]
struct DenseVectors {
    values: Vec<Vector3>,
    unit: Unit,
}

impl DenseVectors {
    fn from_vec(values: Quantity<Vec<Vector3>>, unit: Unit) -> Result<Self, ConformationError> {
        let factor = values.unit().conversion_factor_to(unit)?;
        let mut values = values.into_value();
        for (index, vector) in values.iter_mut().enumerate() {
            let converted = *vector * factor;
            if !converted.is_finite() {
                return Err(ConformationError::NonFiniteVector { index });
            }
            *vector = converted;
        }
        Ok(Self { values, unit })
    }

    fn new<T>(values: Quantity<T>, unit: Unit) -> Result<Self, ConformationError>
    where
        T: AsRef<[Vector3]>,
    {
        let values = Quantity::new(values.value().as_ref().to_vec(), values.unit());
        Self::from_vec(values, unit)
    }

    fn select(&self, rows: &[usize]) -> Self {
        Self {
            values: rows.iter().map(|row| self.values[*row]).collect(),
            unit: self.unit,
        }
    }

    fn set_all<T>(&mut self, values: Quantity<T>) -> Result<(), ConformationError>
    where
        T: AsRef<[Vector3]>,
    {
        // AsRef may use interior state. Validate and copy the same borrowed slice.
        let factor = values.unit().conversion_factor_to(self.unit)?;
        let source = values.value().as_ref();
        if source.len() != self.values.len() {
            return Err(ConformationError::AtomCountMismatch {
                expected: self.values.len(),
                actual: source.len(),
            });
        }
        if let Some(index) = source
            .iter()
            .position(|vector| !(*vector * factor).is_finite())
        {
            return Err(ConformationError::NonFiniteVector { index });
        }
        for (destination, source) in self.values.iter_mut().zip(source) {
            *destination = *source * factor;
        }
        Ok(())
    }

    fn rotated(&self, rotate: impl Fn(Vector3) -> Vector3) -> Result<Self, ConformationError> {
        let values = self.values.iter().map(|vector| rotate(*vector)).collect();
        Self::from_vec(Quantity::new(values, self.unit), self.unit)
    }
}

macro_rules! vector_array {
    ($name:ident, $unit:expr, $doc:literal) => {
        #[doc = $doc]
        ///
        /// Values are stored in the canonical unit in dense atom order and carry
        /// no topology.
        #[derive(Debug, Clone, PartialEq)]
        pub struct $name(DenseVectors);

        impl $name {
            /// Copies dense vectors and converts them to canonical units.
            pub fn new<T>(values: Quantity<T>) -> Result<Self, ConformationError>
            where
                T: AsRef<[Vector3]>,
            {
                Ok(Self(DenseVectors::new(values, $unit)?))
            }

            /// Takes ownership of dense vectors, converting and validating in place.
            pub fn from_vec(values: Quantity<Vec<Vector3>>) -> Result<Self, ConformationError> {
                Ok(Self(DenseVectors::from_vec(values, $unit)?))
            }

            pub fn zeros(len: usize) -> Self {
                Self(DenseVectors {
                    values: vec![Vector3::zero(); len],
                    unit: $unit,
                })
            }

            pub fn len(&self) -> usize {
                self.0.values.len()
            }

            pub fn is_empty(&self) -> bool {
                self.0.values.is_empty()
            }

            pub fn values(&self) -> Quantity<&[Vector3]> {
                Quantity::new(self.0.values.as_slice(), self.0.unit)
            }

            /// Consumes the container, returning its canonical-unit vector without copying.
            pub fn into_values(self) -> Quantity<Vec<Vector3>> {
                Quantity::new(self.0.values, self.0.unit)
            }

            /// Replaces every vector without changing the length.
            pub fn set_all<T>(&mut self, values: Quantity<T>) -> Result<(), ConformationError>
            where
                T: AsRef<[Vector3]>,
            {
                self.0.set_all(values)
            }

            pub(crate) fn select(&self, rows: &[usize]) -> Self {
                Self(self.0.select(rows))
            }

            pub(crate) fn rotated(
                &self,
                rotate: impl Fn(Vector3) -> Vector3,
            ) -> Result<Self, ConformationError> {
                Ok(Self(self.0.rotated(rotate)?))
            }
        }
    };
}

vector_array!(
    Velocities,
    CANONICAL_VELOCITY_UNIT,
    "Dense per-atom velocities."
);
vector_array!(Forces, CANONICAL_FORCE_UNIT, "Dense per-atom forces.");
