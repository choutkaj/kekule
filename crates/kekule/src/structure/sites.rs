use crate::topology::TopologyAtomIndex;
use crate::units::{Quantity, SQUARE_NANOMETER};

use super::state::atom_row;
use super::ConformationError;

/// Crystallographic site parameters of one model or ensemble member:
/// per-atom occupancies and isotropic B-factors in dense atom order.
///
/// Trajectory frames never carry site parameters. Each column is allocated on
/// its first recorded value; B-factors are stored in square nanometres. Atom
/// indices are checked against the owner's atom count `len`.
#[derive(Debug, Clone, Default, PartialEq)]
pub(crate) struct SiteParameters {
    occupancies: AtomScalars,
    b_factors: AtomScalars,
}

impl SiteParameters {
    pub(crate) fn occupancy(
        &self,
        atom: TopologyAtomIndex,
        len: usize,
    ) -> Result<Option<f64>, ConformationError> {
        Ok(self.occupancies.get(atom_row(len, atom)?))
    }

    pub(crate) fn occupancies(&self) -> Option<&[Option<f64>]> {
        self.occupancies.values()
    }

    pub(crate) fn b_factor(
        &self,
        atom: TopologyAtomIndex,
        len: usize,
    ) -> Result<Option<Quantity<f64>>, ConformationError> {
        Ok(self
            .b_factors
            .get(atom_row(len, atom)?)
            .map(|value| Quantity::new(value, SQUARE_NANOMETER)))
    }

    pub(crate) fn b_factors(&self) -> Option<Quantity<&[Option<f64>]>> {
        self.b_factors
            .values()
            .map(|values| Quantity::new(values, SQUARE_NANOMETER))
    }

    pub(crate) fn set_occupancy(
        &mut self,
        atom: TopologyAtomIndex,
        len: usize,
        value: Option<f64>,
    ) -> Result<(), ConformationError> {
        self.occupancies.set(atom_row(len, atom)?, len, value)
    }

    pub(crate) fn set_occupancies(
        &mut self,
        len: usize,
        values: Option<Vec<Option<f64>>>,
    ) -> Result<(), ConformationError> {
        self.occupancies.replace(len, values)
    }

    /// Converts the value to square nanometres before storing it.
    pub(crate) fn set_b_factor(
        &mut self,
        atom: TopologyAtomIndex,
        len: usize,
        value: Option<Quantity<f64>>,
    ) -> Result<(), ConformationError> {
        let value = value
            .map(|value| value.into_unit(SQUARE_NANOMETER).map(Quantity::into_value))
            .transpose()?;
        self.b_factors.set(atom_row(len, atom)?, len, value)
    }

    /// Converts every value to square nanometres before storing them.
    pub(crate) fn set_b_factors(
        &mut self,
        len: usize,
        values: Option<Quantity<Vec<Option<f64>>>>,
    ) -> Result<(), ConformationError> {
        let values = values
            .map(|values| {
                let factor = values.unit().conversion_factor_to(SQUARE_NANOMETER)?;
                Ok::<_, ConformationError>(
                    values
                        .into_value()
                        .into_iter()
                        .map(|value| value.map(|value| value * factor))
                        .collect::<Vec<_>>(),
                )
            })
            .transpose()?;
        self.b_factors.replace(len, values)
    }

    pub(crate) fn select(&self, rows: &[usize]) -> Self {
        Self {
            occupancies: self.occupancies.select(rows),
            b_factors: self.b_factors.select(rows),
        }
    }

    pub(crate) fn resize(&mut self, len: usize) {
        self.occupancies.resize(len);
        self.b_factors.resize(len);
    }

    pub(crate) fn clear_row(&mut self, row: usize) {
        self.occupancies.clear(row);
        self.b_factors.clear(row);
    }

    /// Copies every source row into the given destination rows of an owner
    /// with `len` atoms.
    pub(crate) fn copy_rows_from(&mut self, source: &Self, rows: &[usize], len: usize) {
        for (scalars, source) in [
            (&mut self.occupancies, &source.occupancies),
            (&mut self.b_factors, &source.b_factors),
        ] {
            if let Some(values) = source.values() {
                for (&row, value) in rows.iter().zip(values) {
                    scalars
                        .set(row, len, *value)
                        .expect("source values are finite");
                }
            }
        }
    }
}

/// Optional dense per-atom scalars, allocated on the first recorded value.
#[derive(Debug, Clone, Default)]
struct AtomScalars(Option<Vec<Option<f64>>>);

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
    fn values(&self) -> Option<&[Option<f64>]> {
        self.0.as_deref()
    }

    fn get(&self, row: usize) -> Option<f64> {
        self.0.as_ref().and_then(|values| values[row])
    }

    fn set(&mut self, row: usize, len: usize, value: Option<f64>) -> Result<(), ConformationError> {
        if value.is_some_and(|value| !value.is_finite()) {
            return Err(ConformationError::NonFiniteValue { index: row });
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

    fn replace(
        &mut self,
        len: usize,
        values: Option<Vec<Option<f64>>>,
    ) -> Result<(), ConformationError> {
        if let Some(values) = &values {
            if values.len() != len {
                return Err(ConformationError::AtomCountMismatch {
                    expected: len,
                    actual: values.len(),
                });
            }
            if let Some(index) = values
                .iter()
                .position(|value| value.is_some_and(|value| !value.is_finite()))
            {
                return Err(ConformationError::NonFiniteValue { index });
            }
        }
        self.0 = values.filter(|values| values.iter().any(Option::is_some));
        Ok(())
    }

    fn select(&self, rows: &[usize]) -> Self {
        Self(
            self.0
                .as_ref()
                .map(|values| rows.iter().map(|row| values[*row]).collect()),
        )
    }

    fn resize(&mut self, len: usize) {
        if let Some(values) = &mut self.0 {
            values.resize(len, None);
        }
    }

    fn clear(&mut self, row: usize) {
        if let Some(values) = &mut self.0 {
            values[row] = None;
        }
    }
}
