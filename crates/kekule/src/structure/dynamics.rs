use crate::diagnostics::{Countable, LengthMismatch};
use std::sync::Arc;

use crate::math::Vector3;
use crate::structure::DenseArrayError;
use crate::units::{Quantity, Unit, CANONICAL_FORCE_UNIT, CANONICAL_VELOCITY_UNIT};

/// Dense vectors in one canonical unit, shared by clones until one is
/// written.
#[derive(Debug, Clone, PartialEq)]
struct DenseVectors {
    values: Arc<Vec<Vector3>>,
}

impl DenseVectors {
    fn from_vec(values: Quantity<Vec<Vector3>>, unit: Unit) -> Result<Self, DenseArrayError> {
        let factor = values.unit().conversion_factor_to(unit)?;
        let mut values = values.into_value();
        if factor == 1.0 {
            if let Some(index) = values.iter().position(|vector| !vector.is_finite()) {
                return Err(DenseArrayError::NonFinite { index });
            }
        } else {
            for (index, vector) in values.iter_mut().enumerate() {
                let converted = *vector * factor;
                if !converted.is_finite() {
                    return Err(DenseArrayError::NonFinite { index });
                }
                *vector = converted;
            }
        }
        Ok(Self {
            values: Arc::new(values),
        })
    }

    fn new<T>(values: Quantity<T>, unit: Unit) -> Result<Self, DenseArrayError>
    where
        T: AsRef<[Vector3]>,
    {
        let values = Quantity::new(values.value().as_ref().to_vec(), values.unit());
        Self::from_vec(values, unit)
    }

    fn select(&self, rows: &[usize]) -> Self {
        if crate::structure::positions::is_identity(rows, self.values.len()) {
            return self.clone();
        }
        Self {
            values: Arc::new(rows.iter().map(|row| self.values[*row]).collect()),
        }
    }

    fn set_all<T>(&mut self, values: Quantity<T>, unit: Unit) -> Result<(), DenseArrayError>
    where
        T: AsRef<[Vector3]>,
    {
        // AsRef may use interior state. Validate and copy the same borrowed slice.
        let factor = values.unit().conversion_factor_to(unit)?;
        let source = values.value().as_ref();
        if source.len() != self.values.len() {
            return Err(DenseArrayError::CountMismatch(LengthMismatch::new(
                Countable::Row,
                self.values.len(),
                source.len(),
            )));
        }
        if let Some(index) = source
            .iter()
            .position(|vector| !(*vector * factor).is_finite())
        {
            return Err(DenseArrayError::NonFinite { index });
        }
        let values = Arc::make_mut(&mut self.values);
        for (destination, source) in values.iter_mut().zip(source) {
            *destination = *source * factor;
        }
        Ok(())
    }

    fn rotated(&self, rotate: impl Fn(Vector3) -> Vector3) -> Result<Self, DenseArrayError> {
        let mut rotated = Self {
            values: Arc::new(self.values.as_ref().clone()),
        };
        rotated.rotate_in_place(rotate)?;
        Ok(rotated)
    }

    fn rotate_in_place(
        &mut self,
        rotate: impl Fn(Vector3) -> Vector3,
    ) -> Result<(), DenseArrayError> {
        for (index, vector) in Arc::make_mut(&mut self.values).iter_mut().enumerate() {
            *vector = rotate(*vector);
            if !vector.is_finite() {
                return Err(DenseArrayError::NonFinite { index });
            }
        }
        Ok(())
    }
}

macro_rules! vector_array {
    ($name:ident, $unit:expr, $doc:literal) => {
        #[doc = $doc]
        ///
        /// Values are stored in the canonical unit in dense atom order and carry
        /// no topology. Clones share the values until one of them is written,
        /// as for [`Positions`](super::Positions).
        #[derive(Debug, Clone, PartialEq)]
        pub struct $name(DenseVectors);

        impl $name {
            /// Copies dense vectors and converts them to the canonical unit.
            ///
            /// # Errors
            ///
            /// Returns [`DenseArrayError::Unit`] when the unit has the wrong
            /// dimension and [`DenseArrayError::NonFinite`] for the first
            /// non-finite converted vector.
            pub fn new<T>(values: Quantity<T>) -> Result<Self, DenseArrayError>
            where
                T: AsRef<[Vector3]>,
            {
                Ok(Self(DenseVectors::new(values, $unit)?))
            }

            /// Takes ownership of dense vectors, converting and validating in place.
            ///
            /// # Errors
            ///
            /// As for [`Self::new`].
            pub fn from_vec(values: Quantity<Vec<Vector3>>) -> Result<Self, DenseArrayError> {
                Ok(Self(DenseVectors::from_vec(values, $unit)?))
            }

            /// Creates `len` zero vectors.
            pub fn zeros(len: usize) -> Self {
                Self(DenseVectors {
                    values: Arc::new(vec![Vector3::zero(); len]),
                })
            }

            /// Number of per-atom vectors.
            pub fn len(&self) -> usize {
                self.0.values.len()
            }

            /// Whether the array holds no vectors.
            pub fn is_empty(&self) -> bool {
                self.0.values.is_empty()
            }

            /// Borrows the dense vectors in the canonical unit.
            pub fn values(&self) -> Quantity<&[Vector3]> {
                Quantity::new(self.0.values.as_slice(), $unit)
            }

            /// Consumes the container, returning its canonical-unit vector,
            /// without copying unless a clone shares it.
            pub fn into_values(self) -> Quantity<Vec<Vector3>> {
                Quantity::new(Arc::unwrap_or_clone(self.0.values), $unit)
            }

            /// Replaces every vector without changing the length.
            ///
            /// The array is unchanged when an error is returned.
            ///
            /// # Errors
            ///
            /// Returns [`DenseArrayError::Unit`] for a unit of the wrong dimension,
            /// [`DenseArrayError::CountMismatch`] when the length differs, and
            /// [`DenseArrayError::NonFinite`] for the first non-finite converted vector.
            pub fn set_all<T>(&mut self, values: Quantity<T>) -> Result<(), DenseArrayError>
            where
                T: AsRef<[Vector3]>,
            {
                self.0.set_all(values, $unit)
            }

            pub(crate) fn select(&self, rows: &[usize]) -> Self {
                Self(self.0.select(rows))
            }

            pub(crate) fn rotated(
                &self,
                rotate: impl Fn(Vector3) -> Vector3,
            ) -> Result<Self, DenseArrayError> {
                Ok(Self(self.0.rotated(rotate)?))
            }

            /// Rotates every vector in place. On failure some vectors are
            /// rotated, so the caller drops the array.
            pub(crate) fn rotate_in_place(
                &mut self,
                rotate: impl Fn(Vector3) -> Vector3,
            ) -> Result<(), DenseArrayError> {
                self.0.rotate_in_place(rotate)
            }
        }
    };
}

vector_array!(
    Velocities,
    CANONICAL_VELOCITY_UNIT,
    "Dense per-atom velocities; the canonical unit is nm/ps."
);
vector_array!(
    Forces,
    CANONICAL_FORCE_UNIT,
    "Dense per-atom forces; the canonical unit is kJ/(mol·nm)."
);
