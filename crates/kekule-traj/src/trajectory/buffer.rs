use std::sync::Arc;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{
    Properties, PropertyColumn, PropertyError, PropertyKey, PropertyTable, PropertyValue,
};
use kekule::structure::{ModelView, Positions};
use kekule::topology::Topology;
use kekule::units::{Quantity, CANONICAL_TIME_UNIT};

use super::frame::{Forces, TrajectoryFrameView, Velocities};
use super::{validate_atom_count, FrameError};

/// Complete borrowed frame state ready for transactional publication.
///
/// A decoder builds this value only after it has read one complete frame into
/// reusable scratch. [`FrameBuffer::replace_from_data`] validates every field,
/// including atom and bond property-table dimensions, before changing the destination,
/// converts units once, reuses dense-array allocations, and clears optional
/// fields omitted from this value.
#[derive(Debug, Clone, Copy)]
pub struct FrameBufferData<'a> {
    positions: Quantity<&'a [Point3]>,
    cell: Option<PeriodicCell>,
    velocities: Option<Quantity<&'a [Vector3]>>,
    forces: Option<Quantity<&'a [Vector3]>>,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
    properties: Option<&'a Properties>,
}

impl<'a> FrameBufferData<'a> {
    /// Starts complete frame data with required dense positions.
    pub const fn new(positions: Quantity<&'a [Point3]>) -> Self {
        Self {
            positions,
            cell: None,
            velocities: None,
            forces: None,
            time: None,
            step: None,
            properties: None,
        }
    }

    /// Borrows every field from an already validated frame view.
    pub fn from_frame_view(frame: TrajectoryFrameView<'a>) -> Self {
        Self {
            positions: frame.positions.values(),
            cell: frame.cell.copied(),
            velocities: frame.velocities,
            forces: frame.forces,
            time: frame.time,
            step: frame.step,
            properties: Some(frame.properties),
        }
    }

    pub const fn with_cell(mut self, cell: PeriodicCell) -> Self {
        self.cell = Some(cell);
        self
    }

    pub const fn with_velocities(mut self, velocities: Quantity<&'a [Vector3]>) -> Self {
        self.velocities = Some(velocities);
        self
    }

    pub const fn with_forces(mut self, forces: Quantity<&'a [Vector3]>) -> Self {
        self.forces = Some(forces);
        self
    }

    pub const fn with_time(mut self, time: Quantity<f64>) -> Self {
        self.time = Some(time);
        self
    }

    pub const fn with_step(mut self, step: u64) -> Self {
        self.step = Some(step);
        self
    }

    pub const fn with_properties(mut self, properties: &'a Properties) -> Self {
        self.properties = Some(properties);
        self
    }
}

/// Reusable caller-owned frame storage owning one exact topology context.
#[derive(Debug, Clone)]
pub struct FrameBuffer {
    pub(super) topology: Arc<Topology>,
    positions: Positions,
    cell: Option<PeriodicCell>,
    properties: Properties,
    pub(super) velocities: Velocities,
    has_velocities: bool,
    pub(super) forces: Forces,
    has_forces: bool,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
}

impl FrameBuffer {
    pub fn new(topology: Arc<Topology>) -> Self {
        Self {
            positions: Positions::zeros(topology.atom_count()),
            cell: None,
            properties: Properties::realization(topology.atom_count(), topology.bond_count()),
            velocities: Velocities::zeros(topology.atom_count()),
            has_velocities: false,
            forces: Forces::zeros(topology.atom_count()),
            has_forces: false,
            time: None,
            step: None,
            topology,
        }
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    pub fn positions(&self) -> &Positions {
        &self.positions
    }

    pub const fn cell(&self) -> Option<&PeriodicCell> {
        self.cell.as_ref()
    }

    pub const fn properties(&self) -> &Properties {
        &self.properties
    }

    pub fn insert_property(
        &mut self,
        key: PropertyKey,
        value: PropertyValue,
    ) -> Result<Option<PropertyValue>, PropertyError> {
        self.properties.insert(key, value)
    }

    pub fn remove_property(&mut self, key: &PropertyKey) -> Option<PropertyValue> {
        self.properties.remove(key)
    }

    pub fn clear_properties(&mut self) {
        self.properties.clear_owner();
    }

    realization_property_api!();

    fn insert_bond_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, FrameError> {
        Ok(self
            .properties
            .insert_realization_bond_column(key, column)?)
    }

    pub fn set_positions<T>(&mut self, positions: Quantity<T>) -> Result<(), FrameError>
    where
        T: AsRef<[Point3]>,
    {
        self.positions.set_all(positions)?;
        Ok(())
    }

    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.cell = cell;
    }

    pub fn set_velocities<T>(&mut self, velocities: Option<Quantity<T>>) -> Result<(), FrameError>
    where
        T: AsRef<[Vector3]>,
    {
        match velocities {
            Some(values) => {
                self.velocities.set_all(values)?;
                self.has_velocities = true;
            }
            None => self.has_velocities = false,
        }
        Ok(())
    }

    pub fn set_forces<T>(&mut self, forces: Option<Quantity<T>>) -> Result<(), FrameError>
    where
        T: AsRef<[Vector3]>,
    {
        match forces {
            Some(values) => {
                self.forces.set_all(values)?;
                self.has_forces = true;
            }
            None => self.has_forces = false,
        }
        Ok(())
    }

    pub fn set_time(&mut self, time: Option<Quantity<f64>>) -> Result<(), FrameError> {
        self.time = match time {
            Some(time) => {
                let time = time.to_unit(CANONICAL_TIME_UNIT)?;
                if !time.value().is_finite() {
                    return Err(FrameError::NonFiniteTime);
                }
                Some(time)
            }
            None => None,
        };
        Ok(())
    }

    pub fn set_step(&mut self, step: Option<u64>) {
        self.step = step;
    }

    pub fn set_properties(&mut self, properties: Properties) -> Result<(), FrameError> {
        if properties.realization_atom_properties().len() != self.topology.atom_count() {
            return Err(FrameError::AtomCountMismatch {
                expected: self.topology.atom_count(),
                actual: properties.realization_atom_properties().len(),
            });
        }
        if properties.realization_bond_properties().len() != self.topology.bond_count() {
            return Err(FrameError::BondCountMismatch {
                expected: self.topology.bond_count(),
                actual: properties.realization_bond_properties().len(),
            });
        }
        properties.validate_realization_canonical_properties()?;
        self.properties = properties;
        Ok(())
    }

    /// Clears all per-frame state except positions while retaining reusable
    /// array allocations and the bound topology.
    pub fn reset_dynamic_state(&mut self) {
        self.cell = None;
        self.properties =
            Properties::realization(self.topology.atom_count(), self.topology.bond_count());
        self.has_velocities = false;
        self.has_forces = false;
        self.time = None;
        self.step = None;
    }

    pub fn model_view(&self) -> ModelView<'_> {
        ModelView::new(
            &self.topology,
            &self.positions,
            self.cell.as_ref(),
            &self.properties,
        )
        .expect("frame buffer state is bound to its topology")
    }

    pub fn frame_view(&self) -> TrajectoryFrameView<'_> {
        TrajectoryFrameView {
            topology: &self.topology,
            positions: &self.positions,
            cell: self.cell.as_ref(),
            properties: &self.properties,
            velocities: self.has_velocities.then(|| self.velocities.values()),
            forces: self.has_forces.then(|| self.forces.values()),
            time: self.time,
            step: self.step,
        }
    }

    /// Replaces the complete visible frame transactionally.
    ///
    /// All count, unit, finite-value, property-table, and optional-array
    /// validation completes before any destination field changes.
    /// Existing position, velocity, and force allocations are reused. Optional
    /// fields absent from `data`, including properties, are cleared.
    pub fn replace_from_data(&mut self, data: FrameBufferData<'_>) -> Result<(), FrameError> {
        self.positions.validate_all(&data.positions)?;
        let velocities = data
            .velocities
            .map(|values| {
                self.velocities
                    .0
                    .validate_replacement(&values)
                    .map(|factor| (values, factor))
            })
            .transpose()?;
        let forces = data
            .forces
            .map(|values| {
                self.forces
                    .0
                    .validate_replacement(&values)
                    .map(|factor| (values, factor))
            })
            .transpose()?;
        let time = data
            .time
            .map(|time| {
                let time = time.to_unit(CANONICAL_TIME_UNIT)?;
                if !time.value().is_finite() {
                    return Err(FrameError::NonFiniteTime);
                }
                Ok(time)
            })
            .transpose()?;
        if let Some(properties) = data.properties {
            validate_atom_count(
                self.topology.atom_count(),
                properties.realization_atom_properties().len(),
            )?;
            if properties.realization_bond_properties().len() != self.topology.bond_count() {
                return Err(FrameError::BondCountMismatch {
                    expected: self.topology.bond_count(),
                    actual: properties.realization_bond_properties().len(),
                });
            }
            properties.validate_realization_canonical_properties()?;
        }
        let properties = data.properties.cloned().unwrap_or_else(|| {
            Properties::realization(self.topology.atom_count(), self.topology.bond_count())
        });

        self.positions.set_all(data.positions)?;
        self.cell = data.cell;
        match velocities {
            Some((values, factor)) => {
                self.velocities
                    .0
                    .copy_from_validated(values.value(), factor);
                self.has_velocities = true;
            }
            None => self.has_velocities = false,
        }
        match forces {
            Some((values, factor)) => {
                self.forces.0.copy_from_validated(values.value(), factor);
                self.has_forces = true;
            }
            None => self.has_forces = false,
        }
        self.time = time;
        self.step = data.step;
        self.properties = properties;
        Ok(())
    }

    pub fn copy_from(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), FrameError> {
        if !Arc::ptr_eq(&self.topology, frame.topology) {
            return Err(FrameError::TopologyMismatch);
        }
        self.replace_from_data(FrameBufferData::from_frame_view(frame))
    }
}
