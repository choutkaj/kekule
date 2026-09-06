use kekule::geometry::{PeriodicCell, Point3};
use kekule::properties::{Properties, PropertyColumn, PropertyError, PropertyKey, PropertyValue};
use kekule::units::Quantity;

use super::{Forces, FrameError, TrajectoryFrame, Velocities};

/// Dimension-preserving mutable access to a stored trajectory frame.
///
/// Read access dereferences to [`TrajectoryFrame`]. Mutations take effect only
/// after validation; dropping or forgetting this editor cannot bypass checks.
/// Whole-payload replacement goes through [`super::Trajectory::replace_frame`].
#[derive(Debug)]
pub struct TrajectoryFrameMut<'a> {
    frame: &'a mut TrajectoryFrame,
}

impl std::ops::Deref for TrajectoryFrameMut<'_> {
    type Target = TrajectoryFrame;
    fn deref(&self) -> &Self::Target {
        self.frame
    }
}

impl<'a> TrajectoryFrameMut<'a> {
    pub(super) fn new(frame: &'a mut TrajectoryFrame) -> Self {
        Self { frame }
    }

    pub fn set_positions<T: AsRef<[Point3]>>(
        &mut self,
        positions: Quantity<T>,
    ) -> Result<(), FrameError> {
        self.frame.set_positions(positions)
    }

    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.frame.set_cell(cell);
    }
    pub fn set_velocities(&mut self, values: Option<Velocities>) -> Result<(), FrameError> {
        self.frame.set_velocities(values)
    }
    pub fn set_forces(&mut self, values: Option<Forces>) -> Result<(), FrameError> {
        self.frame.set_forces(values)
    }
    pub fn set_time(&mut self, time: Option<Quantity<f64>>) -> Result<(), FrameError> {
        self.frame.set_time(time)
    }
    pub fn set_step(&mut self, step: Option<u64>) {
        self.frame.set_step(step);
    }

    pub fn insert_property(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>, PropertyError> {
        self.frame.insert_property(key, value)
    }
    pub fn remove_property(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        self.frame.remove_property(key)
    }
    pub fn clear_properties(&mut self) {
        self.frame.clear_properties();
    }

    /// Replaces properties while preserving the stored atom and bond domains,
    /// including their dimensions when no columns are populated.
    pub fn set_properties(&mut self, properties: Properties) -> Result<(), FrameError> {
        super::validate_atom_count(
            self.frame.positions.len(),
            properties.realization_atom_properties().len(),
        )?;
        let expected = self.frame.bond_properties().len();
        let actual = properties.realization_bond_properties().len();
        if expected != actual {
            return Err(FrameError::BondCountMismatch { expected, actual });
        }
        properties.validate_realization_properties()?;
        self.frame.properties = properties;
        Ok(())
    }

    pub fn set_atom_property(
        &mut self,
        index: usize,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), FrameError> {
        self.frame.set_atom_property(index, key, value)
    }
    pub fn insert_atom_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, FrameError> {
        self.frame.insert_atom_property_column(key, column)
    }
    pub fn remove_atom_property_column(
        &mut self,
        key: &PropertyKey,
    ) -> Result<Option<PropertyColumn>, FrameError> {
        self.frame.remove_atom_property_column(key)
    }
    pub fn set_bond_property(
        &mut self,
        index: usize,
        key: PropertyKey,
        value: Option<PropertyValue>,
    ) -> Result<(), FrameError> {
        self.frame.set_bond_property(index, key, value)
    }
    pub fn insert_bond_property_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, FrameError> {
        // Detached-frame insertion may establish a bond dimension. Stored frames
        // already have an authoritative domain and must use strict insertion.
        Ok(self
            .frame
            .properties
            .insert_realization_bond_column(key, column)?)
    }
    pub fn remove_bond_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
        self.frame.remove_bond_property_column(key)
    }
    pub fn set_occupancy_at(&mut self, index: usize, value: Option<f64>) -> Result<(), FrameError> {
        self.frame.set_occupancy_at(index, value)
    }
    pub fn set_b_factor_at(
        &mut self,
        index: usize,
        value: Option<Quantity<f64>>,
    ) -> Result<(), FrameError> {
        self.frame.set_b_factor_at(index, value)
    }
}
