use crate::topology::TopologyAtomIndex;
use crate::units::{Quantity, SQUARE_NANOMETER};

use super::state::atom_row;
use super::StructureError;

/// One optional per-atom scalar column in dense atom order, such as the
/// occupancies or B-factors of a model or ensemble member.
///
/// The column is allocated on its first recorded value; a column whose values
/// are all missing equals one that was never recorded. Atom indices are checked
/// against the owner's atom count `len`.
#[derive(Debug, Clone, Default)]
pub(crate) struct AtomScalars(Option<Vec<Option<f64>>>);

/// The column of an owner that never records one, such as a trajectory frame.
pub(crate) static ABSENT: AtomScalars = AtomScalars(None);

impl PartialEq for AtomScalars {
    fn eq(&self, other: &Self) -> bool {
        match (&self.0, &other.0) {
            (Some(left), Some(right)) => left == right,
            (Some(values), None) | (None, Some(values)) => values.iter().all(Option::is_none),
            (None, None) => true,
        }
    }
}

impl AtomScalars {
    /// Every value in dense atom order, or `None` when none was recorded.
    pub(crate) fn values(&self) -> Option<&[Option<f64>]> {
        self.0.as_deref()
    }

    pub(crate) fn get(
        &self,
        atom: TopologyAtomIndex,
        len: usize,
    ) -> Result<Option<f64>, StructureError> {
        let row = atom_row(len, atom)?;
        Ok(self.0.as_ref().and_then(|values| values[row]))
    }

    /// Sets or clears one value; values must be finite.
    pub(crate) fn set(
        &mut self,
        atom: TopologyAtomIndex,
        len: usize,
        value: Option<f64>,
    ) -> Result<(), StructureError> {
        let row = atom_row(len, atom)?;
        self.set_row(row, len, value)
    }

    /// Replaces every value in dense atom order, or clears them all.
    pub(crate) fn replace(
        &mut self,
        len: usize,
        values: Option<Vec<Option<f64>>>,
    ) -> Result<(), StructureError> {
        if let Some(values) = &values {
            if values.len() != len {
                return Err(StructureError::AtomCountMismatch {
                    expected: len,
                    actual: values.len(),
                });
            }
            if let Some(index) = values
                .iter()
                .position(|value| value.is_some_and(|value| !value.is_finite()))
            {
                return Err(StructureError::NonFiniteValue { index });
            }
        }
        self.0 = values.filter(|values| values.iter().any(Option::is_some));
        Ok(())
    }

    pub(crate) fn select(&self, rows: &[usize]) -> Self {
        Self(
            self.0
                .as_ref()
                .map(|values| rows.iter().map(|row| values[*row]).collect()),
        )
    }

    pub(crate) fn resize(&mut self, len: usize) {
        if let Some(values) = &mut self.0 {
            values.resize(len, None);
        }
    }

    pub(crate) fn clear(&mut self, row: usize) {
        if let Some(values) = &mut self.0 {
            values[row] = None;
        }
    }

    /// Copies every source row into the given destination rows of an owner
    /// with `len` atoms.
    pub(crate) fn copy_rows_from(&mut self, source: &Self, rows: &[usize], len: usize) {
        if let Some(values) = source.values() {
            for (&row, value) in rows.iter().zip(values) {
                self.set_row(row, len, *value)
                    .expect("source values are finite");
            }
        }
    }

    fn set_row(
        &mut self,
        row: usize,
        len: usize,
        value: Option<f64>,
    ) -> Result<(), StructureError> {
        if value.is_some_and(|value| !value.is_finite()) {
            return Err(StructureError::NonFiniteValue { index: row });
        }
        match (&mut self.0, value) {
            (Some(values), value) => values[row] = value,
            (None, Some(value)) => {
                let mut values = vec![None; len];
                values[row] = Some(value);
                self.0 = Some(values);
            }
            (None, None) => {}
        }
        Ok(())
    }
}

/// Converts one B-factor to the stored square nanometres.
pub(crate) fn canonical_b_factor(
    value: Option<Quantity<f64>>,
) -> Result<Option<f64>, StructureError> {
    Ok(value
        .map(|value| value.into_unit(SQUARE_NANOMETER).map(Quantity::into_value))
        .transpose()?)
}

/// Converts a B-factor column to the stored square nanometres.
pub(crate) fn canonical_b_factors(
    values: Option<Quantity<Vec<Option<f64>>>>,
) -> Result<Option<Vec<Option<f64>>>, StructureError> {
    values
        .map(|values| {
            let factor = values.unit().conversion_factor_to(SQUARE_NANOMETER)?;
            Ok(values
                .into_value()
                .into_iter()
                .map(|value| value.map(|value| value * factor))
                .collect())
        })
        .transpose()
}

/// Reads one stored B-factor with its unit.
pub(crate) fn b_factor_quantity(value: Option<f64>) -> Option<Quantity<f64>> {
    value.map(|value| Quantity::new(value, SQUARE_NANOMETER))
}

/// Reads a stored B-factor column with its unit.
pub(crate) fn b_factor_quantities(values: &AtomScalars) -> Option<Quantity<&[Option<f64>]>> {
    values
        .values()
        .map(|values| Quantity::new(values, SQUARE_NANOMETER))
}
