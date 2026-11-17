//! The per-geometry columns of a set: [`GeometryIndex`], [`SetColumnsMut`],
//! and [`ColumnRef`].

use crate::properties::{
    LiveRows, PropertyColumn, PropertyError, PropertyKey, PropertyTable, PropertyTableMut,
    PropertyValue,
};
use crate::structure::sampling::is_energy;
use crate::structure::Sampling;
use crate::units::{Quantity, CANONICAL_ENERGY_UNIT, CANONICAL_TIME_UNIT, DIMENSIONLESS};

use super::{GeometryStorage, Order, Requirement, StructureSet, StructureSetError};

impl From<PropertyError> for StructureSetError {
    fn from(error: PropertyError) -> Self {
        Self::Column(Box::new(error))
    }
}

/// The position of one geometry in its set: the row of the set's columns
/// ([`StructureSet::columns`]).
///
/// It is the geometry's index in the set. It shifts when an earlier geometry
/// is removed ([`StructureSet::remove`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GeometryIndex(usize);

impl GeometryIndex {
    /// The row of geometry `index`.
    pub const fn new(index: usize) -> Self {
        Self(index)
    }

    /// The geometry's index in its set.
    pub const fn index(self) -> usize {
        self.0
    }
}

/// A value per geometry that [`StructureSet::column_values`] reads: a value
/// of the order, or a set column.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ColumnRef {
    /// The times of a time order, in ps.
    Time,
    /// The steps of a step order or of a time order with steps.
    Step,
    /// The values of a coordinate order, in its unit.
    Coordinate,
    /// A real or integer set column.
    Key(PropertyKey),
}

/// Length-preserving mutable access to a set's columns: one row per geometry.
///
/// Reads go through `Deref` to the [`PropertyTable`]. Every mutation is
/// transactional, as [`PropertyTableMut`]'s: on error the columns are
/// unchanged. A failure of the table's own checks is
/// [`StructureSetError::Column`].
/// A column the set's [`Sampling`] record names stays a real column of
/// energies: an edit that would put other values in it fails with
/// [`StructureSetError::SamplingColumn`], and removing it fails with
/// [`StructureSetError::ColumnInUse`]. Clearing all its cells is allowed; the
/// column is then absent, which the record reads as missing values.
#[derive(Debug)]
pub struct SetColumnsMut<'a> {
    table: PropertyTableMut<'a, GeometryIndex>,
    sampling: Option<&'a Sampling>,
}

impl std::ops::Deref for SetColumnsMut<'_> {
    type Target = PropertyTable<GeometryIndex>;

    fn deref(&self) -> &Self::Target {
        &self.table
    }
}

impl SetColumnsMut<'_> {
    /// Inserts or replaces a complete column of one cell per geometry, as
    /// [`PropertyTableMut::insert`] does, returning the replaced one.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::SamplingColumn`] for values other than real
    /// energies in a column the sampling record names, and
    /// [`StructureSetError::Column`] with the errors of
    /// [`PropertyTableMut::insert`].
    pub fn insert(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, StructureSetError> {
        if let Some(sampling) = self.sampling {
            let energies = matches!(column, PropertyColumn::Real { unit, .. } if is_energy(unit));
            if !energies && !column.is_all_missing() && sampling.names(&key) {
                return Err(StructureSetError::SamplingColumn { key });
            }
        }
        Ok(self.table.insert(key, column)?)
    }

    /// Removes and returns the column under `key`, if present.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::ColumnInUse`] for a present column the sampling
    /// record names.
    pub fn remove(
        &mut self,
        key: &PropertyKey,
    ) -> Result<Option<PropertyColumn>, StructureSetError> {
        if self.table.get(key).is_some()
            && self.sampling.is_some_and(|sampling| sampling.names(key))
        {
            return Err(StructureSetError::ColumnInUse { key: key.clone() });
        }
        Ok(self.table.remove(key))
    }

    /// Removes every column.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::ColumnInUse`] when a column the sampling record
    /// names is present.
    pub fn clear(&mut self) -> Result<(), StructureSetError> {
        if let Some(sampling) = self.sampling {
            if let Some(key) = sampling.columns().find(|key| self.table.get(key).is_some()) {
                return Err(StructureSetError::ColumnInUse { key: key.clone() });
            }
        }
        self.table.clear();
        Ok(())
    }

    /// Sets or clears the cell of one geometry, as
    /// [`PropertyTableMut::set_value`] does.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::SamplingColumn`] for a value other than a real
    /// energy in a column the sampling record names, and
    /// [`StructureSetError::Column`] with the errors of
    /// [`PropertyTableMut::set_value`].
    pub fn set_value(
        &mut self,
        key: PropertyKey,
        geometry: GeometryIndex,
        value: Option<PropertyValue>,
    ) -> Result<(), StructureSetError> {
        if let Some(sampling) = self.sampling {
            sampling.check_value(&key, value.as_ref())?;
        }
        Ok(self.table.set_value(key, geometry, value)?)
    }

    /// Clears the cell of one geometry, as [`PropertyTableMut::clear_value`]
    /// does.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::Column`] with the errors of
    /// [`PropertyTableMut::clear_value`].
    pub fn clear_value(
        &mut self,
        key: PropertyKey,
        geometry: GeometryIndex,
    ) -> Result<(), StructureSetError> {
        Ok(self.table.clear_value(key, geometry)?)
    }

    /// Applies several cell updates to one column transactionally, as
    /// [`PropertyTableMut::set_values`] does.
    ///
    /// # Errors
    ///
    /// As [`Self::set_value`] for the first update that fails.
    pub fn set_values(
        &mut self,
        key: PropertyKey,
        values: impl IntoIterator<Item = (GeometryIndex, Option<PropertyValue>)>,
    ) -> Result<(), StructureSetError> {
        let values = values.into_iter().collect::<Vec<_>>();
        if let Some(sampling) = self.sampling {
            for (_, value) in &values {
                sampling.check_value(&key, value.as_ref())?;
            }
        }
        Ok(self.table.set_values(key, values)?)
    }
}

impl<S: GeometryStorage> StructureSet<S> {
    /// The set's columns: per-geometry scalars, one row per geometry, each
    /// cell possibly missing.
    ///
    /// A column belongs to the set, as the order does: selecting keeps each
    /// selected geometry's row, and a removed geometry takes its row with
    /// it. Geometries added later have missing cells.
    ///
    /// ```
    /// use kekule::properties::{PropertyColumn, PropertyKey};
    /// use kekule::structure::{GeometryIndex, Positions, Geometry, StructureSet};
    /// use kekule::units::KILOJOULE_PER_MOLE;
    ///
    /// let topology = kekule::smiles::to_topology("O")?;
    /// let geometry = || Geometry::new(Positions::zeros(topology.atom_count()));
    /// let mut set = StructureSet::from_geometries(topology.clone(), [geometry(), geometry()])?;
    /// let energy = PropertyKey::new("energy")?;
    /// set.columns_mut().insert(
    ///     energy.clone(),
    ///     PropertyColumn::Real { unit: KILOJOULE_PER_MOLE, values: vec![Some(-1.5), None] },
    /// )?;
    /// set.push(geometry(), None)?;
    /// assert_eq!(set.columns().len(), 3);
    /// assert_eq!(set.columns().value(&energy, GeometryIndex::new(2))?, None);
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn columns(&self) -> &PropertyTable<GeometryIndex> {
        &self.columns
    }

    /// Length-preserving mutable access to the set's columns.
    pub fn columns_mut(&mut self) -> SetColumnsMut<'_> {
        SetColumnsMut {
            table: PropertyTableMut::new(&mut self.columns, LiveRows::All, None),
            sampling: self.sampling.as_ref(),
        }
    }

    /// The values of `column` for every geometry, as numbers in their unit:
    /// times in ps, steps and integer cells dimensionless, coordinates and
    /// real cells in their stored unit. A missing cell is `None`; order
    /// values are never missing. Integers beyond 2^53 lose precision.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::Missing`] for a set without the order values
    /// asked for: `TimeOrder` or `Times` for [`ColumnRef::Time`], `Steps` for
    /// [`ColumnRef::Step`], and `CoordinateOrder` for
    /// [`ColumnRef::Coordinate`]; [`StructureSetError::MissingColumn`] for a
    /// key the set has no column of and its sampling record does not name (a
    /// named column without values reads as missing energies, in kJ/mol);
    /// and
    /// [`StructureSetError::NonNumericColumn`] for a boolean or string column.
    pub fn column_values(
        &self,
        column: &ColumnRef,
    ) -> Result<Quantity<Vec<Option<f64>>>, StructureSetError> {
        match column {
            ColumnRef::Time => {
                let Some(Order::Time(axis)) = self.order() else {
                    return Err(Requirement::TimeOrder.into());
                };
                let times = axis.times().ok_or(Requirement::Times)?;
                Ok(Quantity::new(
                    times.value().iter().copied().map(Some).collect(),
                    CANONICAL_TIME_UNIT,
                ))
            }
            ColumnRef::Step => {
                let steps = self.steps().ok_or(Requirement::Steps)?;
                Ok(Quantity::new(
                    steps.iter().map(|&step| Some(step as f64)).collect(),
                    DIMENSIONLESS,
                ))
            }
            ColumnRef::Coordinate => {
                let Some(Order::Coordinate { unit, values }) = self.order() else {
                    return Err(Requirement::CoordinateOrder.into());
                };
                Ok(Quantity::new(
                    values.iter().copied().map(Some).collect(),
                    *unit,
                ))
            }
            ColumnRef::Key(key) => match self.columns.get(key) {
                None if self
                    .sampling
                    .as_ref()
                    .is_some_and(|sampling| sampling.names(key)) =>
                {
                    Ok(Quantity::new(vec![None; self.len()], CANONICAL_ENERGY_UNIT))
                }
                None => Err(StructureSetError::MissingColumn { key: key.clone() }),
                Some(PropertyColumn::Real { unit, values }) => {
                    Ok(Quantity::new(values.clone(), *unit))
                }
                Some(PropertyColumn::Int(values)) => Ok(Quantity::new(
                    values
                        .iter()
                        .map(|value| value.map(|value| value as f64))
                        .collect(),
                    DIMENSIONLESS,
                )),
                Some(_) => Err(StructureSetError::NonNumericColumn { key: key.clone() }),
            },
        }
    }

    /// The MD step or step number of every geometry: the steps of a time
    /// order with steps, or of a step order. `None` without them.
    pub fn steps(&self) -> Option<&[u64]> {
        match self.order()? {
            Order::Time(axis) => axis.steps(),
            Order::Step(steps) => Some(steps),
            Order::Coordinate { .. } => None,
        }
    }
}
