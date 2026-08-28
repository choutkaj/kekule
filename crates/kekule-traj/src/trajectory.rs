//! Fixed-topology trajectory frames, reusable buffers, in-memory storage, and
//! streaming reader/writer contracts.

use std::{fmt, io, sync::Arc};

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{
    Properties, PropertyColumn, PropertyError, PropertyKey, PropertyTable, PropertyValue,
};
use kekule::structure::{ModelError, ModelView, PositionError, Positions};
use kekule::topology::transform::TopologySubsetError;
use kekule::topology::{AtomSelection, InstanceAtomId, Topology};
use kekule::units::{
    Quantity, Unit, UnitError, CANONICAL_FORCE_UNIT, CANONICAL_TIME_UNIT, CANONICAL_VELOCITY_UNIT,
};

#[derive(Debug, Clone, PartialEq)]
struct DenseVectors {
    values: Vec<Vector3>,
    unit: Unit,
}

impl DenseVectors {
    fn new<T>(values: Quantity<T>, unit: Unit) -> Result<Self, FrameError>
    where
        T: AsRef<[Vector3]>,
    {
        let factor = values.unit().conversion_factor_to(unit)?;
        let values = values
            .value()
            .as_ref()
            .iter()
            .copied()
            .enumerate()
            .map(|(index, vector)| {
                let vector = Vector3::new(vector.x * factor, vector.y * factor, vector.z * factor);
                if !vector.is_finite() {
                    return Err(FrameError::NonFiniteVector { index });
                }
                Ok(vector)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self { values, unit })
    }

    fn zeros(len: usize, unit: Unit) -> Self {
        Self {
            values: vec![Vector3::zero(); len],
            unit,
        }
    }

    fn len(&self) -> usize {
        self.values.len()
    }

    fn values(&self) -> Quantity<&[Vector3]> {
        Quantity::new(self.values.as_slice(), self.unit)
    }

    fn select_indices(&self, indices: &[usize]) -> Result<Self, FrameError> {
        let values = indices
            .iter()
            .map(|index| {
                self.values
                    .get(*index)
                    .copied()
                    .ok_or(FrameError::InvalidIndex { index: *index })
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            values,
            unit: self.unit,
        })
    }

    fn set_all<T>(&mut self, values: Quantity<T>) -> Result<(), FrameError>
    where
        T: AsRef<[Vector3]>,
    {
        let factor = self.validate_replacement(&values)?;
        self.copy_from_validated(values.value().as_ref(), factor);
        Ok(())
    }

    fn validate_replacement<T>(&self, values: &Quantity<T>) -> Result<f64, FrameError>
    where
        T: AsRef<[Vector3]>,
    {
        let factor = values.unit().conversion_factor_to(self.unit)?;
        let source = values.value().as_ref();
        if source.len() != self.values.len() {
            return Err(FrameError::AtomCountMismatch {
                expected: self.values.len(),
                actual: source.len(),
            });
        }
        for (index, vector) in source.iter().copied().enumerate() {
            let converted = Vector3::new(vector.x * factor, vector.y * factor, vector.z * factor);
            if !converted.is_finite() {
                return Err(FrameError::NonFiniteVector { index });
            }
        }
        Ok(factor)
    }

    fn copy_from_validated(&mut self, source: &[Vector3], factor: f64) {
        for (destination, source) in self.values.iter_mut().zip(source.iter().copied()) {
            *destination = Vector3::new(source.x * factor, source.y * factor, source.z * factor);
        }
    }
}

macro_rules! vector_array {
    ($name:ident, $unit:expr) => {
        #[doc = concat!("Dense numerical ", stringify!($name), "; contains no topology.")]
        #[derive(Debug, Clone, PartialEq)]
        pub struct $name(DenseVectors);

        impl $name {
            pub fn new<T>(values: Quantity<T>) -> Result<Self, FrameError>
            where
                T: AsRef<[Vector3]>,
            {
                Ok(Self(DenseVectors::new(values, $unit)?))
            }

            pub fn zeros(len: usize) -> Self {
                Self(DenseVectors::zeros(len, $unit))
            }

            pub fn len(&self) -> usize {
                self.0.len()
            }

            pub fn is_empty(&self) -> bool {
                self.len() == 0
            }

            pub fn values(&self) -> Quantity<&[Vector3]> {
                self.0.values()
            }

            /// Copies a dense projection in the requested index order.
            pub fn select_indices(&self, indices: &[usize]) -> Result<Self, FrameError> {
                Ok(Self(self.0.select_indices(indices)?))
            }

            pub fn set_all<T>(&mut self, values: Quantity<T>) -> Result<(), FrameError>
            where
                T: AsRef<[Vector3]>,
            {
                self.0.set_all(values)
            }
        }
    };
}

vector_array!(Velocities, CANONICAL_VELOCITY_UNIT);
vector_array!(Forces, CANONICAL_FORCE_UNIT);

macro_rules! realization_property_api {
    () => {
        pub const fn atom_properties(&self) -> &PropertyTable {
            self.properties.realization_atom_properties()
        }

        pub const fn bond_properties(&self) -> &PropertyTable {
            self.properties.realization_bond_properties()
        }

        pub fn atom_property(
            &self,
            key: &PropertyKey,
            index: usize,
        ) -> Result<Option<PropertyValue>, FrameError> {
            Ok(self.atom_properties().value(key, index)?)
        }

        pub fn set_atom_property(
            &mut self,
            key: PropertyKey,
            index: usize,
            value: Option<PropertyValue>,
        ) -> Result<(), FrameError> {
            Ok(self
                .properties
                .set_realization_atom_value(key, index, value)?)
        }

        pub fn insert_atom_property_column(
            &mut self,
            key: PropertyKey,
            column: PropertyColumn,
        ) -> Result<Option<PropertyColumn>, FrameError> {
            Ok(self
                .properties
                .insert_realization_atom_column(key, column)?)
        }

        pub fn remove_atom_property_column(
            &mut self,
            key: &PropertyKey,
        ) -> Result<Option<PropertyColumn>, FrameError> {
            Ok(self.properties.remove_realization_atom_column(key)?)
        }

        pub fn bond_property(
            &self,
            key: &PropertyKey,
            index: usize,
        ) -> Result<Option<PropertyValue>, FrameError> {
            Ok(self.bond_properties().value(key, index)?)
        }

        pub fn set_bond_property(
            &mut self,
            key: PropertyKey,
            index: usize,
            value: Option<PropertyValue>,
        ) -> Result<(), FrameError> {
            Ok(self
                .properties
                .set_realization_bond_value(key, index, value)?)
        }

        pub fn insert_bond_property_column(
            &mut self,
            key: PropertyKey,
            column: PropertyColumn,
        ) -> Result<Option<PropertyColumn>, FrameError> {
            Ok(self
                .properties
                .insert_realization_bond_column(key, column)?)
        }

        pub fn remove_bond_property_column(&mut self, key: &PropertyKey) -> Option<PropertyColumn> {
            self.properties.remove_realization_bond_column(key)
        }

        pub fn occupancy_at(&self, index: usize) -> Result<Option<f64>, FrameError> {
            Ok(self.properties.occupancy_at(index)?)
        }

        pub fn set_occupancy_at(
            &mut self,
            index: usize,
            value: Option<f64>,
        ) -> Result<(), FrameError> {
            Ok(self.properties.set_occupancy_at(index, value)?)
        }

        pub fn b_factor_at(&self, index: usize) -> Result<Option<Quantity<f64>>, FrameError> {
            Ok(self.properties.b_factor_at(index)?)
        }

        pub fn set_b_factor_at(
            &mut self,
            index: usize,
            value: Option<Quantity<f64>>,
        ) -> Result<(), FrameError> {
            Ok(self.properties.set_b_factor_at(index, value)?)
        }
    };
}

#[derive(Debug, Clone, PartialEq)]
pub struct TrajectoryFrame {
    positions: Positions,
    cell: Option<PeriodicCell>,
    properties: Properties,
    velocities: Option<Velocities>,
    forces: Option<Forces>,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
}

impl TrajectoryFrame {
    pub fn new(positions: Positions, bond_count: usize) -> Self {
        let properties = Properties::realization(positions.len(), bond_count);
        Self {
            positions,
            cell: None,
            properties,
            velocities: None,
            forces: None,
            time: None,
            step: None,
        }
    }

    pub fn positions(&self) -> &Positions {
        &self.positions
    }

    pub const fn cell(&self) -> Option<&PeriodicCell> {
        self.cell.as_ref()
    }

    pub fn set_cell(&mut self, cell: Option<PeriodicCell>) {
        self.cell = cell;
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

    pub fn set_properties(&mut self, properties: Properties) -> Result<(), FrameError> {
        if properties.realization_atom_properties().len() != self.positions.len() {
            return Err(FrameError::AtomCountMismatch {
                expected: self.positions.len(),
                actual: properties.realization_atom_properties().len(),
            });
        }
        if properties.realization_bond_properties().len() != self.bond_properties().len() {
            return Err(FrameError::BondCountMismatch {
                expected: self.bond_properties().len(),
                actual: properties.realization_bond_properties().len(),
            });
        }
        properties.validate_realization_canonical_properties()?;
        self.properties = properties;
        Ok(())
    }

    pub fn velocities(&self) -> Option<&Velocities> {
        self.velocities.as_ref()
    }

    pub fn forces(&self) -> Option<&Forces> {
        self.forces.as_ref()
    }

    pub const fn time(&self) -> Option<Quantity<f64>> {
        self.time
    }

    pub const fn step(&self) -> Option<u64> {
        self.step
    }

    pub fn set_velocities(&mut self, velocities: Option<Velocities>) -> Result<(), FrameError> {
        if let Some(values) = velocities.as_ref() {
            validate_atom_count(self.positions.len(), values.len())?;
        }
        self.velocities = velocities;
        Ok(())
    }

    pub fn set_forces(&mut self, forces: Option<Forces>) -> Result<(), FrameError> {
        if let Some(values) = forces.as_ref() {
            validate_atom_count(self.positions.len(), values.len())?;
        }
        self.forces = forces;
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

    pub fn validate(&self, topology: &Arc<Topology>) -> Result<(), FrameError> {
        validate_atom_count(topology.atom_count(), self.positions.len())?;
        validate_atom_count(topology.atom_count(), self.atom_properties().len())?;
        if self.bond_properties().len() != topology.bond_count() {
            return Err(FrameError::BondCountMismatch {
                expected: topology.bond_count(),
                actual: self.bond_properties().len(),
            });
        }
        if let Some(values) = &self.velocities {
            validate_atom_count(topology.atom_count(), values.len())?;
        }
        if let Some(values) = &self.forces {
            validate_atom_count(topology.atom_count(), values.len())?;
        }
        self.properties
            .validate_realization_canonical_properties()?;
        Ok(())
    }

    pub fn view<'a>(
        &'a self,
        topology: &'a Arc<Topology>,
    ) -> Result<TrajectoryFrameView<'a>, FrameError> {
        self.validate(topology)?;
        Ok(TrajectoryFrameView {
            topology,
            positions: &self.positions,
            cell: self.cell.as_ref(),
            properties: &self.properties,
            velocities: self.velocities.as_ref().map(Velocities::values),
            forces: self.forces.as_ref().map(Forces::values),
            time: self.time,
            step: self.step,
        })
    }
}

/// Borrowed trajectory frame state.
#[derive(Debug, Clone, Copy)]
pub struct TrajectoryFrameView<'a> {
    topology: &'a Arc<Topology>,
    positions: &'a Positions,
    cell: Option<&'a PeriodicCell>,
    properties: &'a Properties,
    velocities: Option<Quantity<&'a [Vector3]>>,
    forces: Option<Quantity<&'a [Vector3]>>,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
}

impl<'a> TrajectoryFrameView<'a> {
    pub fn topology(self) -> &'a Topology {
        self.topology
    }

    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.topology)
    }

    pub(crate) const fn topology_arc(self) -> &'a Arc<Topology> {
        self.topology
    }

    pub const fn positions(self) -> &'a Positions {
        self.positions
    }

    pub const fn cell(self) -> Option<&'a PeriodicCell> {
        self.cell
    }

    pub const fn properties(self) -> &'a Properties {
        self.properties
    }

    pub const fn atom_properties(self) -> &'a PropertyTable {
        self.properties.realization_atom_properties()
    }

    pub const fn bond_properties(self) -> &'a PropertyTable {
        self.properties.realization_bond_properties()
    }

    pub fn atom_property(
        self,
        key: &PropertyKey,
        index: usize,
    ) -> Result<Option<PropertyValue>, FrameError> {
        Ok(self.atom_properties().value(key, index)?)
    }

    pub fn bond_property(
        self,
        key: &PropertyKey,
        index: usize,
    ) -> Result<Option<PropertyValue>, FrameError> {
        Ok(self.bond_properties().value(key, index)?)
    }

    pub fn occupancy_at(self, index: usize) -> Result<Option<f64>, FrameError> {
        Ok(self.properties.occupancy_at(index)?)
    }

    pub fn b_factor_at(self, index: usize) -> Result<Option<Quantity<f64>>, FrameError> {
        Ok(self.properties.b_factor_at(index)?)
    }

    pub const fn velocities(self) -> Option<Quantity<&'a [Vector3]>> {
        self.velocities
    }

    pub const fn forces(self) -> Option<Quantity<&'a [Vector3]>> {
        self.forces
    }

    pub const fn time(self) -> Option<Quantity<f64>> {
        self.time
    }

    pub const fn step(self) -> Option<u64> {
        self.step
    }

    pub fn model_view(self) -> ModelView<'a> {
        ModelView::new(self.topology, self.positions, self.cell, self.properties)
            .expect("trajectory frame view has validated topology")
    }
}

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
    topology: Arc<Topology>,
    positions: Positions,
    cell: Option<PeriodicCell>,
    properties: Properties,
    velocities: Velocities,
    has_velocities: bool,
    forces: Forces,
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

/// Deliberately loaded finite in-memory trajectory.
#[derive(Debug, Clone)]
pub struct Trajectory {
    topology: Arc<Topology>,
    properties: Properties,
    frames: Vec<TrajectoryFrame>,
}

impl Trajectory {
    pub fn new(topology: Arc<Topology>) -> Self {
        Self {
            topology,
            properties: Properties::new(),
            frames: Vec::new(),
        }
    }

    pub fn from_frames(
        topology: Arc<Topology>,
        frames: impl IntoIterator<Item = TrajectoryFrame>,
    ) -> Result<Self, TrajectoryError> {
        let mut trajectory = Self::new(topology);
        for frame in frames {
            trajectory.push(frame)?;
        }
        Ok(trajectory)
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
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

    /// Constructs one topology subset and applies its dense mapping to every frame.
    pub fn slice(&self, selection: &AtomSelection) -> Result<Self, TrajectorySliceError> {
        let subset = self.topology.subset(selection)?;
        let atom_indices = subset
            .correspondence()
            .source_atom_indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>();
        let bond_indices = subset
            .correspondence()
            .source_bond_indices()
            .iter()
            .map(|index| index.index())
            .collect::<Vec<_>>();
        let mut target = Self::new(subset.shared_topology());
        for frame in &self.frames {
            target.push(TrajectoryFrame {
                positions: frame.positions.select_indices(&atom_indices)?,
                cell: frame.cell,
                properties: frame
                    .properties
                    .project_realization(&atom_indices, &bond_indices)?,
                velocities: frame
                    .velocities
                    .as_ref()
                    .map(|values| values.select_indices(&atom_indices))
                    .transpose()?,
                forces: frame
                    .forces
                    .as_ref()
                    .map(|values| values.select_indices(&atom_indices))
                    .transpose()?,
                time: frame.time,
                step: frame.step,
            })?;
        }
        Ok(target)
    }

    pub fn len(&self) -> usize {
        self.frames.len()
    }

    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    pub fn frame(&self, index: usize) -> Option<&TrajectoryFrame> {
        self.frames.get(index)
    }

    pub fn frames(&self) -> impl ExactSizeIterator<Item = TrajectoryFrameView<'_>> {
        self.frames.iter().map(|frame| {
            frame
                .view(&self.topology)
                .expect("trajectory validates frame topology on insertion")
        })
    }

    pub fn push(&mut self, frame: TrajectoryFrame) -> Result<(), TrajectoryError> {
        frame.validate(&self.topology)?;
        self.frames.push(frame);
        Ok(())
    }

    /// Explicitly validates that all present times are non-decreasing.
    pub fn validate_monotonic_time(&self, require_all: bool) -> Result<(), TrajectoryError> {
        let mut previous = None;
        for (index, frame) in self.frames.iter().enumerate() {
            let Some(time) = frame.time else {
                if require_all {
                    return Err(TrajectoryError::MissingTime { frame: index });
                }
                continue;
            };
            let value = time.to_value();
            if previous.is_some_and(|previous| value < previous) {
                return Err(TrajectoryError::NonMonotonicTime { frame: index });
            }
            previous = Some(value);
        }
        Ok(())
    }
}

pub trait TrajectoryReader {
    fn topology(&self) -> &Topology;

    fn shared_topology(&self) -> Arc<Topology>;

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError>;
}

pub trait SeekableTrajectoryReader: TrajectoryReader {
    fn frame_count(&self) -> Option<u64>;

    fn read_frame(
        &mut self,
        index: u64,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryError>;
}

pub trait TrajectoryWriter {
    fn topology(&self) -> &Topology;

    fn shared_topology(&self) -> Arc<Topology>;

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError>;
}

/// Sequential and seekable reader over an in-memory trajectory.
pub struct MemoryTrajectoryReader<'a> {
    trajectory: &'a Trajectory,
    cursor: usize,
}

impl<'a> MemoryTrajectoryReader<'a> {
    pub fn new(trajectory: &'a Trajectory) -> Self {
        Self {
            trajectory,
            cursor: 0,
        }
    }
}

impl TrajectoryReader for MemoryTrajectoryReader<'_> {
    fn topology(&self) -> &Topology {
        self.trajectory.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.trajectory.shared_topology()
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        let Some(frame) = self.trajectory.frames.get(self.cursor) else {
            return Ok(false);
        };
        destination.copy_from(frame.view(&self.trajectory.topology)?)?;
        self.cursor += 1;
        Ok(true)
    }
}

impl SeekableTrajectoryReader for MemoryTrajectoryReader<'_> {
    fn frame_count(&self) -> Option<u64> {
        u64::try_from(self.trajectory.len()).ok()
    }

    fn read_frame(
        &mut self,
        index: u64,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryError> {
        let index =
            usize::try_from(index).map_err(|_| TrajectoryError::FrameIndexOutOfRange(index))?;
        let frame = self
            .trajectory
            .frames
            .get(index)
            .ok_or(TrajectoryError::FrameIndexOutOfRange(index as u64))?;
        destination.copy_from(frame.view(&self.trajectory.topology)?)?;
        self.cursor = index.saturating_add(1);
        Ok(())
    }
}

/// In-memory writer that preserves every core frame field.
pub struct MemoryTrajectoryWriter {
    trajectory: Trajectory,
}

impl MemoryTrajectoryWriter {
    pub fn new(topology: Arc<Topology>) -> Self {
        Self {
            trajectory: Trajectory::new(topology),
        }
    }

    pub fn to_trajectory(self) -> Trajectory {
        self.trajectory
    }
}

impl TrajectoryWriter for MemoryTrajectoryWriter {
    fn topology(&self) -> &Topology {
        self.trajectory.topology()
    }

    fn shared_topology(&self) -> Arc<Topology> {
        self.trajectory.shared_topology()
    }

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError> {
        if !Arc::ptr_eq(&self.trajectory.topology, frame.topology) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        let positions = Positions::new(frame.positions.values())?;
        let mut owned = TrajectoryFrame::new(positions, self.trajectory.topology.bond_count());
        owned.cell = frame.cell.copied();
        owned.properties = frame.properties.clone();
        owned.velocities = frame.velocities.map(Velocities::new).transpose()?;
        owned.forces = frame.forces.map(Forces::new).transpose()?;
        owned.time = frame.time;
        owned.step = frame.step;
        self.trajectory.push(owned)
    }
}

/// Proof that a coordinate-only source order exactly matches one topology.
#[derive(Debug, Clone)]
pub struct AtomOrderAssertion {
    topology: Arc<Topology>,
    kind: AtomOrderAssertionKind,
}

impl PartialEq for AtomOrderAssertion {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology) && self.kind == other.kind
    }
}

impl Eq for AtomOrderAssertion {}

/// Evidence represented by an [`AtomOrderAssertion`].
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum AtomOrderAssertionKind {
    SemanticOrder,
    DeclaredTopologyOrder,
}

impl AtomOrderAssertion {
    /// Proves that an explicit semantic atom sequence is the topology's exact
    /// authoritative dense order.
    pub fn from_semantic_order(
        topology: &Arc<Topology>,
        atom_order: &[InstanceAtomId],
    ) -> Result<Self, TrajectoryError> {
        if topology.atom_ids() != atom_order {
            return Err(TrajectoryError::AtomOrderMismatch);
        }
        Ok(Self {
            topology: Arc::clone(topology),
            kind: AtomOrderAssertionKind::SemanticOrder,
        })
    }

    /// Records the caller's explicit assertion that a topology-free file uses
    /// this topology's authoritative dense atom order.
    ///
    /// This is evidence supplied by the caller, not an inference from atom
    /// count. Format readers must still validate all stronger file metadata.
    pub fn assert_file_uses_topology_order(topology: &Arc<Topology>) -> Self {
        Self {
            topology: Arc::clone(topology),
            kind: AtomOrderAssertionKind::DeclaredTopologyOrder,
        }
    }

    pub fn is_compatible(&self, topology: &Arc<Topology>) -> bool {
        Arc::ptr_eq(&self.topology, topology)
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub const fn kind(&self) -> AtomOrderAssertionKind {
        self.kind
    }

    /// Backward-compatible spelling for [`Self::from_semantic_order`].
    pub fn new(
        topology: &Arc<Topology>,
        atom_order: &[InstanceAtomId],
    ) -> Result<Self, TrajectoryError> {
        Self::from_semantic_order(topology, atom_order)
    }
}

/// Reference reader for a topology-free coordinate source.
pub struct CoordinateFrameReader {
    topology: Arc<Topology>,
    frames: Vec<Vec<Point3>>,
    cursor: usize,
}

impl CoordinateFrameReader {
    pub fn new(
        topology: Arc<Topology>,
        assertion: AtomOrderAssertion,
        frames: impl IntoIterator<Item = Quantity<Vec<Point3>>>,
    ) -> Result<Self, TrajectoryError> {
        if !assertion.is_compatible(&topology) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        let frames = frames
            .into_iter()
            .map(|frame| {
                let positions = Positions::new(frame).map_err(TrajectoryError::Position)?;
                validate_atom_count(topology.atom_count(), positions.len())?;
                Ok::<_, TrajectoryError>(positions.values().value().to_vec())
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            topology,
            frames,
            cursor: 0,
        })
    }
}

impl TrajectoryReader for CoordinateFrameReader {
    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        if !Arc::ptr_eq(&self.topology, &destination.topology) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        let Some(frame) = self.frames.get(self.cursor) else {
            return Ok(false);
        };
        destination.set_positions(Quantity::new(
            frame.as_slice(),
            kekule::units::CANONICAL_LENGTH_UNIT,
        ))?;
        destination.reset_dynamic_state();
        self.cursor += 1;
        Ok(true)
    }
}

fn validate_atom_count(expected: usize, actual: usize) -> Result<(), FrameError> {
    if actual != expected {
        return Err(FrameError::AtomCountMismatch { expected, actual });
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum FrameError {
    TopologyMismatch,
    InvalidIndex { index: usize },
    AtomCountMismatch { expected: usize, actual: usize },
    BondCountMismatch { expected: usize, actual: usize },
    NonFiniteVector { index: usize },
    NonFiniteTime,
    Position(PositionError),
    Property(PropertyError),
    Unit(UnitError),
}

impl fmt::Display for FrameError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => {
                formatter.write_str("trajectory frame belongs to a different topology")
            }
            Self::InvalidIndex { index } => write!(formatter, "invalid dense frame index {index}"),
            Self::AtomCountMismatch { expected, actual } => write!(
                formatter,
                "trajectory array requires {expected} atoms, but received {actual}"
            ),
            Self::BondCountMismatch { expected, actual } => write!(
                formatter,
                "trajectory array requires {expected} bonds, but received {actual}"
            ),
            Self::NonFiniteVector { index } => {
                write!(formatter, "trajectory vector at {index} is not finite")
            }
            Self::NonFiniteTime => formatter.write_str("trajectory time must be finite"),
            Self::Position(error) => write!(formatter, "invalid frame positions: {error}"),
            Self::Property(error) => write!(formatter, "invalid frame property: {error}"),
            Self::Unit(error) => write!(formatter, "invalid frame quantity unit: {error}"),
        }
    }
}

impl std::error::Error for FrameError {}

impl From<PositionError> for FrameError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}

impl From<PropertyError> for FrameError {
    fn from(error: PropertyError) -> Self {
        Self::Property(error)
    }
}

impl From<UnitError> for FrameError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

/// Stable identity for a trajectory file format.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum TrajectoryFormat {
    Xyz,
    Dcd,
    Xtc,
    Trr,
}

impl TrajectoryFormat {
    pub const fn name(self) -> &'static str {
        match self {
            Self::Xyz => "XYZ",
            Self::Dcd => "DCD",
            Self::Xtc => "XTC",
            Self::Trr => "TRR",
        }
    }
}

impl fmt::Display for TrajectoryFormat {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.name())
    }
}

/// File or stream operation active when trajectory I/O failed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum TrajectoryIoOperation {
    Detect,
    Open,
    Index,
    ReadHeader,
    ReadFrame,
    WriteHeader,
    WriteFrame,
    Finish,
}

impl fmt::Display for TrajectoryIoOperation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let name = match self {
            Self::Detect => "detect",
            Self::Open => "open",
            Self::Index => "index",
            Self::ReadHeader => "read header",
            Self::ReadFrame => "read frame",
            Self::WriteHeader => "write header",
            Self::WriteFrame => "write frame",
            Self::Finish => "finish",
        };
        formatter.write_str(name)
    }
}

/// Typed classification for malformed, unsupported, or unsafe codec input.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum TrajectoryCodecErrorKind {
    UnknownFormat,
    FormatMismatch,
    InvalidHeader,
    UnsupportedVariant,
    TruncatedRecord,
    InvalidRecordLength,
    RecordMarkerMismatch,
    InvalidFrame,
    InconsistentAtomCount,
    InconsistentMetadata,
    InvalidPrecision,
    ResourceLimitExceeded,
    UnsupportedField,
    NegativeOrUnrepresentableStep,
    CorruptCompressedData,
}

impl fmt::Display for TrajectoryCodecErrorKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let description = match self {
            Self::UnknownFormat => "unknown trajectory format",
            Self::FormatMismatch => "trajectory format mismatch",
            Self::InvalidHeader => "invalid trajectory header",
            Self::UnsupportedVariant => "unsupported trajectory variant",
            Self::TruncatedRecord => "truncated trajectory record",
            Self::InvalidRecordLength => "invalid trajectory record length",
            Self::RecordMarkerMismatch => "trajectory record markers do not match",
            Self::InvalidFrame => "invalid trajectory frame",
            Self::InconsistentAtomCount => "inconsistent trajectory atom count",
            Self::InconsistentMetadata => "inconsistent trajectory metadata",
            Self::InvalidPrecision => "invalid trajectory precision",
            Self::ResourceLimitExceeded => "trajectory resource limit exceeded",
            Self::UnsupportedField => "unsupported trajectory field",
            Self::NegativeOrUnrepresentableStep => "negative or unrepresentable trajectory step",
            Self::CorruptCompressedData => "corrupt compressed trajectory data",
        };
        formatter.write_str(description)
    }
}

/// Cloneable typed context for an underlying file or stream error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrajectoryIoErrorContext {
    operation: TrajectoryIoOperation,
    format: Option<TrajectoryFormat>,
    source_label: Option<String>,
    frame: Option<u64>,
    byte_offset: Option<u64>,
    error_kind: io::ErrorKind,
    message: String,
}

impl TrajectoryIoErrorContext {
    pub fn new(
        operation: TrajectoryIoOperation,
        error_kind: io::ErrorKind,
        message: impl Into<String>,
    ) -> Self {
        Self {
            operation,
            format: None,
            source_label: None,
            frame: None,
            byte_offset: None,
            error_kind,
            message: message.into(),
        }
    }

    pub fn with_format(mut self, format: TrajectoryFormat) -> Self {
        self.format = Some(format);
        self
    }

    pub fn with_source_label(mut self, source_label: impl Into<String>) -> Self {
        self.source_label = Some(source_label.into());
        self
    }

    pub const fn with_frame(mut self, frame: u64) -> Self {
        self.frame = Some(frame);
        self
    }

    pub const fn with_byte_offset(mut self, byte_offset: u64) -> Self {
        self.byte_offset = Some(byte_offset);
        self
    }

    pub const fn operation(&self) -> TrajectoryIoOperation {
        self.operation
    }

    pub const fn format(&self) -> Option<TrajectoryFormat> {
        self.format
    }

    pub fn source_label(&self) -> Option<&str> {
        self.source_label.as_deref()
    }

    pub const fn frame(&self) -> Option<u64> {
        self.frame
    }

    pub const fn byte_offset(&self) -> Option<u64> {
        self.byte_offset
    }

    pub const fn error_kind(&self) -> io::ErrorKind {
        self.error_kind
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

/// Typed context for a codec validation or capability error.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrajectoryCodecErrorContext {
    kind: TrajectoryCodecErrorKind,
    operation: TrajectoryIoOperation,
    format: Option<TrajectoryFormat>,
    source_label: Option<String>,
    frame: Option<u64>,
    byte_offset: Option<u64>,
    expected: Option<u64>,
    actual: Option<u64>,
    detail: Option<String>,
}

impl TrajectoryCodecErrorContext {
    pub const fn new(
        kind: TrajectoryCodecErrorKind,
        operation: TrajectoryIoOperation,
        format: Option<TrajectoryFormat>,
    ) -> Self {
        Self {
            kind,
            operation,
            format,
            source_label: None,
            frame: None,
            byte_offset: None,
            expected: None,
            actual: None,
            detail: None,
        }
    }

    pub fn with_source_label(mut self, source_label: impl Into<String>) -> Self {
        self.source_label = Some(source_label.into());
        self
    }

    pub const fn with_frame(mut self, frame: u64) -> Self {
        self.frame = Some(frame);
        self
    }

    pub const fn with_byte_offset(mut self, byte_offset: u64) -> Self {
        self.byte_offset = Some(byte_offset);
        self
    }

    pub const fn with_counts(mut self, expected: u64, actual: u64) -> Self {
        self.expected = Some(expected);
        self.actual = Some(actual);
        self
    }

    pub fn with_detail(mut self, detail: impl Into<String>) -> Self {
        self.detail = Some(detail.into());
        self
    }

    pub const fn kind(&self) -> TrajectoryCodecErrorKind {
        self.kind
    }

    pub const fn operation(&self) -> TrajectoryIoOperation {
        self.operation
    }

    pub const fn format(&self) -> Option<TrajectoryFormat> {
        self.format
    }

    pub fn source_label(&self) -> Option<&str> {
        self.source_label.as_deref()
    }

    pub const fn frame(&self) -> Option<u64> {
        self.frame
    }

    pub const fn byte_offset(&self) -> Option<u64> {
        self.byte_offset
    }

    pub const fn expected(&self) -> Option<u64> {
        self.expected
    }

    pub const fn actual(&self) -> Option<u64> {
        self.actual
    }

    pub fn detail(&self) -> Option<&str> {
        self.detail.as_deref()
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TrajectoryError {
    TopologyMismatch,
    AtomOrderMismatch,
    FrameIndexOutOfRange(u64),
    UnsupportedRandomAccess,
    MissingRequiredTopology,
    MissingTime { frame: usize },
    NonMonotonicTime { frame: usize },
    UnsupportedField(&'static str),
    Frame(Box<FrameError>),
    Position(PositionError),
    Io(Box<TrajectoryIoErrorContext>),
    Codec(Box<TrajectoryCodecErrorContext>),
}

/// Failure to subset a trajectory topology or transfer frame state.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TrajectorySliceError {
    Topology(TopologySubsetError),
    Position(PositionError),
    Property(PropertyError),
    Frame(Box<FrameError>),
    Trajectory(TrajectoryError),
}

impl fmt::Display for TrajectorySliceError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "cannot slice trajectory: {self:?}")
    }
}

impl std::error::Error for TrajectorySliceError {}

impl From<TopologySubsetError> for TrajectorySliceError {
    fn from(error: TopologySubsetError) -> Self {
        Self::Topology(error)
    }
}
impl From<PositionError> for TrajectorySliceError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}
impl From<PropertyError> for TrajectorySliceError {
    fn from(error: PropertyError) -> Self {
        Self::Property(error)
    }
}
impl From<FrameError> for TrajectorySliceError {
    fn from(error: FrameError) -> Self {
        Self::Frame(Box::new(error))
    }
}
impl From<TrajectoryError> for TrajectorySliceError {
    fn from(error: TrajectoryError) -> Self {
        Self::Trajectory(error)
    }
}

impl fmt::Display for TrajectoryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => {
                formatter.write_str("trajectory object belongs to a different topology")
            }
            Self::AtomOrderMismatch => {
                formatter.write_str("coordinate-source atom order does not match topology order")
            }
            Self::FrameIndexOutOfRange(index) => {
                write!(formatter, "trajectory frame index {index} is out of range")
            }
            Self::UnsupportedRandomAccess => {
                formatter.write_str("trajectory source does not support random access")
            }
            Self::MissingRequiredTopology => {
                formatter.write_str("coordinate-only trajectory source requires a topology")
            }
            Self::MissingTime { frame } => {
                write!(formatter, "trajectory frame {frame} has no time")
            }
            Self::NonMonotonicTime { frame } => {
                write!(formatter, "trajectory time decreases at frame {frame}")
            }
            Self::UnsupportedField(field) => {
                write!(formatter, "trajectory writer does not support {field}")
            }
            Self::Frame(error) => write!(formatter, "invalid trajectory frame: {error}"),
            Self::Position(error) => write!(formatter, "invalid trajectory positions: {error}"),
            Self::Io(context) => {
                write!(formatter, "trajectory {} I/O failed", context.operation)?;
                if let Some(format) = context.format {
                    write!(formatter, " for {format}")?;
                }
                if let Some(source) = &context.source_label {
                    write!(formatter, " at {source}")?;
                }
                if let Some(frame) = context.frame {
                    write!(formatter, " in frame {frame}")?;
                }
                if let Some(offset) = context.byte_offset {
                    write!(formatter, " at byte {offset}")?;
                }
                write!(
                    formatter,
                    ": {} ({:?})",
                    context.message, context.error_kind
                )
            }
            Self::Codec(context) => {
                write!(formatter, "{}", context.kind)?;
                if let Some(format) = context.format {
                    write!(formatter, " for {format}")?;
                }
                write!(formatter, " while attempting to {}", context.operation)?;
                if let Some(source) = &context.source_label {
                    write!(formatter, " at {source}")?;
                }
                if let Some(frame) = context.frame {
                    write!(formatter, " in frame {frame}")?;
                }
                if let Some(offset) = context.byte_offset {
                    write!(formatter, " at byte {offset}")?;
                }
                if let (Some(expected), Some(actual)) = (context.expected, context.actual) {
                    write!(formatter, " (expected {expected}, actual {actual})")?;
                }
                if let Some(detail) = &context.detail {
                    write!(formatter, ": {detail}")?;
                }
                Ok(())
            }
        }
    }
}

impl std::error::Error for TrajectoryError {}

impl From<FrameError> for TrajectoryError {
    fn from(error: FrameError) -> Self {
        match error {
            FrameError::TopologyMismatch => Self::TopologyMismatch,
            error => Self::Frame(Box::new(error)),
        }
    }
}

impl From<PositionError> for TrajectoryError {
    fn from(error: PositionError) -> Self {
        Self::Position(error)
    }
}

impl From<ModelError> for TrajectoryError {
    fn from(_: ModelError) -> Self {
        Self::TopologyMismatch
    }
}

impl From<TrajectoryIoErrorContext> for TrajectoryError {
    fn from(context: TrajectoryIoErrorContext) -> Self {
        Self::Io(Box::new(context))
    }
}

impl From<TrajectoryCodecErrorContext> for TrajectoryError {
    fn from(context: TrajectoryCodecErrorContext) -> Self {
        Self::Codec(Box::new(context))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use kekule::core::{Atom, BondOrder, Element, MoleculeEditor};
    use kekule::geometry::Point3;
    use kekule::properties::{Properties, PropertyKey, PropertyValue};
    use kekule::topology::{AtomSelection, TopologyBuilder};
    use kekule::units::{
        ANGSTROM, DIMENSIONLESS, KELVIN, KILOJOULE_PER_MOLE, NANOMETER, PICOSECOND,
        SQUARE_ANGSTROM, SQUARE_NANOMETER,
    };
    use std::sync::Arc;

    fn make_topology(with_bond: bool) -> Arc<Topology> {
        let mut editor = MoleculeEditor::new();
        let first = editor
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .unwrap();
        if with_bond {
            let second = editor
                .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
                .unwrap();
            editor.add_bond(first, second, BondOrder::Single).unwrap();
        }
        let molecule = editor.finish().unwrap();
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_molecule_definition(&molecule).unwrap();
        builder.add_instance(definition).unwrap();
        Arc::new(builder.build().unwrap())
    }

    fn positions(values: &[Point3]) -> Positions {
        Positions::new(Quantity::new(values, NANOMETER)).unwrap()
    }

    fn key(value: &str) -> PropertyKey {
        PropertyKey::new(value).unwrap()
    }

    fn real(value: f64) -> PropertyValue {
        PropertyValue::Real {
            value,
            unit: DIMENSIONLESS,
        }
    }

    #[test]
    fn vector_arrays_are_topology_free_unit_aware_and_equal_by_values() {
        let vectors = [Vector3::new(1.0, 2.0, 3.0)];
        let velocities = Velocities::new(Quantity::new(vectors, ANGSTROM / PICOSECOND)).unwrap();
        let same = Velocities::new(Quantity::new(vectors, CANONICAL_VELOCITY_UNIT)).unwrap();
        let converted_velocity = velocities.values().value()[0];
        assert!((converted_velocity.x - 0.1).abs() < 1.0e-12);
        assert!((converted_velocity.y - 0.2).abs() < 1.0e-12);
        assert!((converted_velocity.z - 0.3).abs() < 1.0e-12);
        assert_ne!(velocities, same);
        assert_eq!(velocities.len(), 1);
        assert!(!velocities.is_empty());

        let forces = Forces::new(Quantity::new(vectors, KILOJOULE_PER_MOLE / ANGSTROM)).unwrap();
        let converted_force = forces.values().value()[0];
        assert!((converted_force.x - 10.0).abs() < 1.0e-12);
        assert!((converted_force.y - 20.0).abs() < 1.0e-12);
        assert!((converted_force.z - 30.0).abs() < 1.0e-12);
        assert!(matches!(
            Velocities::new(Quantity::new(
                [Vector3::new(f64::NAN, 0.0, 0.0)],
                CANONICAL_VELOCITY_UNIT
            )),
            Err(FrameError::NonFiniteVector { index: 0 })
        ));
        assert!(matches!(
            Forces::new(Quantity::new(vectors, KELVIN)),
            Err(FrameError::Unit(UnitError::IncompatibleUnits { .. }))
        ));
    }

    #[test]
    fn frames_validate_all_dense_dimensions_at_the_owner_boundary() {
        let topology = make_topology(true);
        let valid_positions = positions(&[Point3::origin(), Point3::origin()]);
        let mut frame = TrajectoryFrame::new(valid_positions, topology.bond_count());
        frame
            .set_velocities(Some(Velocities::zeros(topology.atom_count())))
            .unwrap();
        frame
            .set_forces(Some(Forces::zeros(topology.atom_count())))
            .unwrap();
        frame.validate(&topology).unwrap();

        assert!(matches!(
            TrajectoryFrame::new(positions(&[Point3::origin()]), topology.bond_count())
                .validate(&topology),
            Err(FrameError::AtomCountMismatch {
                expected: 2,
                actual: 1
            })
        ));
        assert!(matches!(
            frame.set_velocities(Some(Velocities::zeros(1))),
            Err(FrameError::AtomCountMismatch {
                expected: 2,
                actual: 1
            })
        ));
        assert!(matches!(
            frame.set_forces(Some(Forces::zeros(1))),
            Err(FrameError::AtomCountMismatch {
                expected: 2,
                actual: 1
            })
        ));
        assert!(matches!(
            frame.set_properties(Properties::realization(1, topology.bond_count())),
            Err(FrameError::AtomCountMismatch { .. })
        ));
        assert!(matches!(
            frame.set_properties(Properties::realization(topology.atom_count(), 0)),
            Err(FrameError::BondCountMismatch { .. })
        ));
    }

    #[test]
    fn frame_view_borrows_model_state_and_preserves_time_and_step() {
        let topology = make_topology(false);
        let atom = topology.atom_ids()[0];
        let mut frame = TrajectoryFrame::new(
            positions(&[Point3::new(3.0, 0.0, 0.0)]),
            topology.bond_count(),
        );
        frame
            .set_atom_property(key("score"), 0, Some(real(0.8)))
            .unwrap();
        frame
            .set_time(Some(Quantity::new(2.5, PICOSECOND)))
            .unwrap();
        frame.set_step(Some(7));

        let view = frame.view(&topology).unwrap();
        assert_eq!(view.model_view().position(atom).unwrap().value().x, 3.0);
        assert_eq!(
            view.atom_property(&key("score"), 0).unwrap(),
            Some(PropertyValue::Real {
                value: 0.8,
                unit: DIMENSIONLESS,
            })
        );
        assert_eq!(view.time(), Some(Quantity::new(2.5, PICOSECOND)));
        assert_eq!(view.step(), Some(7));
        assert_eq!(
            view.positions().values().value().as_ptr(),
            frame.positions().values().value().as_ptr()
        );
    }

    #[test]
    fn trajectory_rejects_bad_frame_dimensions_and_keeps_one_topology() {
        let topology = make_topology(true);
        let mut trajectory = Trajectory::new(Arc::clone(&topology));
        let mut valid = TrajectoryFrame::new(
            positions(&[Point3::origin(), Point3::origin()]),
            topology.bond_count(),
        );
        valid.set_step(Some(1));
        trajectory.push(valid).unwrap();
        assert!(Arc::ptr_eq(&trajectory.shared_topology(), &topology));
        assert_eq!(trajectory.frames().next().unwrap().step(), Some(1));

        let wrong = TrajectoryFrame::new(positions(&[Point3::origin()]), topology.bond_count());
        assert!(matches!(
            trajectory.push(wrong),
            Err(TrajectoryError::Frame(error))
                if matches!(*error, FrameError::AtomCountMismatch { .. })
        ));
    }

    #[test]
    fn frame_and_buffer_scope_canonical_atom_properties_to_semantic_apis() {
        let topology = make_topology(true);
        let mut frame = TrajectoryFrame::new(
            positions(&[Point3::origin(), Point3::origin()]),
            topology.bond_count(),
        );
        assert!(matches!(
            frame.set_atom_property(key("occupancy"), 0, Some(PropertyValue::Int(1))),
            Err(FrameError::Property(PropertyError::ReservedKey(_)))
        ));
        frame.set_occupancy_at(0, Some(0.8)).unwrap();
        frame
            .set_b_factor_at(0, Some(Quantity::new(25.0, SQUARE_ANGSTROM)))
            .unwrap();
        assert_eq!(frame.occupancy_at(0).unwrap(), Some(0.8));
        let b_factor = frame.b_factor_at(0).unwrap().unwrap();
        assert_eq!(b_factor.unit(), SQUARE_NANOMETER);
        assert!((*b_factor.value() - 0.25).abs() < 1.0e-12);
        assert!(frame.set_occupancy_at(0, Some(f64::NAN)).is_err());
        assert!(frame
            .set_b_factor_at(0, Some(Quantity::new(1.0, KELVIN)))
            .is_err());

        let mut buffer = FrameBuffer::new(Arc::clone(&topology));
        assert!(matches!(
            buffer.insert_atom_property_column(
                key("b_factor"),
                PropertyColumn::String(vec![Some("bad".into()), None]),
            ),
            Err(FrameError::Property(PropertyError::ReservedKey(_)))
        ));
        buffer.set_properties(frame.properties().clone()).unwrap();
        assert_eq!(buffer.occupancy_at(0).unwrap(), Some(0.8));
        assert_eq!(buffer.b_factor_at(0).unwrap(), Some(b_factor));
        assert!(matches!(
            buffer.remove_atom_property_column(&key("occupancy")),
            Err(FrameError::Property(PropertyError::ReservedKey(_)))
        ));
    }

    #[test]
    fn trajectory_slice_transfers_every_per_atom_frame_field() {
        let topology = make_topology(true);
        let atoms = topology.atom_ids();
        let selection = AtomSelection::from_atoms(&topology, [atoms[1]]).unwrap();
        let mut frame = TrajectoryFrame::new(
            positions(&[Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)]),
            topology.bond_count(),
        );
        frame
            .set_atom_property(key("score"), 0, Some(real(0.25)))
            .unwrap();
        frame
            .set_atom_property(key("score"), 1, Some(real(0.75)))
            .unwrap();
        frame
            .set_bond_property(key("score"), 0, Some(real(9.0)))
            .unwrap();
        frame
            .insert_property(key("frame_energy"), real(12.0))
            .unwrap();
        frame
            .set_velocities(Some(
                Velocities::new(Quantity::new(
                    [Vector3::new(3.0, 0.0, 0.0), Vector3::new(4.0, 0.0, 0.0)],
                    CANONICAL_VELOCITY_UNIT,
                ))
                .unwrap(),
            ))
            .unwrap();
        frame
            .set_forces(Some(
                Forces::new(Quantity::new(
                    [Vector3::new(5.0, 0.0, 0.0), Vector3::new(6.0, 0.0, 0.0)],
                    CANONICAL_FORCE_UNIT,
                ))
                .unwrap(),
            ))
            .unwrap();
        frame
            .set_time(Some(Quantity::new(2.0, PICOSECOND)))
            .unwrap();
        frame.set_step(Some(8));
        let mut trajectory = Trajectory::from_frames(Arc::clone(&topology), [frame]).unwrap();
        trajectory
            .insert_property(
                key("collection_source"),
                PropertyValue::String("test".into()),
            )
            .unwrap();

        let sliced = trajectory.slice(&selection).unwrap();
        assert!(sliced.properties().owner_is_empty());
        assert_eq!(sliced.topology().atom_count(), 1);
        assert_eq!(sliced.topology().bond_count(), 0);
        let frame = sliced.frame(0).unwrap();
        assert!(frame.properties().owner_is_empty());
        assert_eq!(
            frame.positions().values().value(),
            &[Point3::new(2.0, 0.0, 0.0)]
        );
        assert_eq!(
            frame.atom_properties().value(&key("score"), 0).unwrap(),
            Some(real(0.75))
        );
        assert!(!frame.bond_properties().has_data());
        assert_eq!(frame.velocities().unwrap().values().value()[0].x, 4.0);
        assert_eq!(frame.forces().unwrap().values().value()[0].x, 6.0);
        assert_eq!(frame.time(), Some(Quantity::new(2.0, PICOSECOND)));
        assert_eq!(frame.step(), Some(8));
    }

    #[test]
    fn frame_buffer_publication_is_transactional_reuses_arrays_and_clears_optionals() {
        let topology = make_topology(true);
        let mut buffer = FrameBuffer::new(Arc::clone(&topology));
        let position_ptr = buffer.positions().values().value().as_ptr();
        let velocity_ptr = buffer.velocities.values().value().as_ptr();
        let force_ptr = buffer.forces.values().value().as_ptr();
        let points = [Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)];
        let vectors = [Vector3::new(1.0, 0.0, 0.0); 2];
        let mut properties = Properties::realization(topology.atom_count(), topology.bond_count());
        properties
            .set_realization_atom_value(key("atom_score"), 0, Some(real(0.75)))
            .unwrap();
        properties
            .set_realization_bond_value(key("bond_score"), 0, Some(real(2.0)))
            .unwrap();
        properties
            .insert(key("codec_value"), PropertyValue::Int(7))
            .unwrap();

        buffer
            .replace_from_data(
                FrameBufferData::new(Quantity::new(points.as_slice(), ANGSTROM))
                    .with_velocities(Quantity::new(vectors.as_slice(), CANONICAL_VELOCITY_UNIT))
                    .with_forces(Quantity::new(vectors.as_slice(), CANONICAL_FORCE_UNIT))
                    .with_time(Quantity::new(1.0, PICOSECOND))
                    .with_step(4)
                    .with_properties(&properties),
            )
            .unwrap();
        assert!(buffer.frame_view().velocities().is_some());
        assert!(buffer.frame_view().forces().is_some());
        assert_eq!(buffer.positions().values().value().as_ptr(), position_ptr);
        assert_eq!(buffer.velocities.values().value().as_ptr(), velocity_ptr);
        assert_eq!(buffer.forces.values().value().as_ptr(), force_ptr);
        assert!(buffer.atom_properties().has_data());
        assert!(buffer.bond_properties().has_data());
        assert_eq!(
            buffer.properties().get(&key("codec_value")),
            Some(&PropertyValue::Int(7))
        );

        let before = buffer.frame_view().positions().values().value().to_vec();
        assert!(matches!(
            buffer.replace_from_data(FrameBufferData::new(Quantity::new(
                [Point3::new(f64::NAN, 0.0, 0.0)].as_slice(),
                ANGSTROM
            ))),
            Err(FrameError::Position(_))
        ));
        assert_eq!(*buffer.positions().values().value(), before.as_slice());

        buffer
            .replace_from_data(FrameBufferData::new(Quantity::new(
                points.as_slice(),
                ANGSTROM,
            )))
            .unwrap();
        assert!(buffer.frame_view().velocities().is_none());
        assert!(buffer.frame_view().forces().is_none());
        assert!(buffer.frame_view().time().is_none());
        assert!(buffer.frame_view().step().is_none());
        assert!(buffer.properties().is_empty());
    }

    #[test]
    fn frame_buffer_copy_requires_the_exact_frame_topology() {
        let topology = make_topology(false);
        let independent = make_topology(false);
        let frame = TrajectoryFrame::new(positions(&[Point3::origin()]), topology.bond_count());
        let view = frame.view(&topology).unwrap();
        let mut correct = FrameBuffer::new(Arc::clone(&topology));
        correct.copy_from(view).unwrap();

        let mut wrong = FrameBuffer::new(independent);
        assert_eq!(wrong.copy_from(view), Err(FrameError::TopologyMismatch));
    }

    #[test]
    fn memory_reader_and_writer_round_trip_validated_frames() {
        let topology = make_topology(true);
        let mut trajectory = Trajectory::new(Arc::clone(&topology));
        let mut frame = TrajectoryFrame::new(
            positions(&[Point3::new(5.0, 0.0, 0.0), Point3::new(6.0, 0.0, 0.0)]),
            topology.bond_count(),
        );
        frame.set_step(Some(9));
        frame
            .set_atom_property(key("atom_score"), 0, Some(real(0.6)))
            .unwrap();
        frame
            .set_bond_property(key("bond_score"), 0, Some(real(4.0)))
            .unwrap();
        assert_eq!(
            frame.atom_property(&key("atom_score"), 0).unwrap(),
            Some(real(0.6))
        );
        assert_eq!(
            frame.bond_property(&key("bond_score"), 0).unwrap(),
            Some(real(4.0))
        );
        trajectory.push(frame).unwrap();

        let mut reader = MemoryTrajectoryReader::new(&trajectory);
        let mut buffer = FrameBuffer::new(Arc::clone(&topology));
        assert!(reader.read_next(&mut buffer).unwrap());
        assert_eq!(buffer.positions().values().value()[0].x, 5.0);
        assert_eq!(buffer.frame_view().step(), Some(9));
        assert_eq!(
            buffer.atom_property(&key("atom_score"), 0).unwrap(),
            Some(real(0.6))
        );
        assert_eq!(
            buffer.bond_property(&key("bond_score"), 0).unwrap(),
            Some(real(4.0))
        );
        assert!(!reader.read_next(&mut buffer).unwrap());

        let mut writer = MemoryTrajectoryWriter::new(Arc::clone(&topology));
        writer.write_frame(buffer.frame_view()).unwrap();
        let written = writer.to_trajectory();
        assert_eq!(written.len(), 1);
        assert_eq!(written.frames().next().unwrap().step(), Some(9));
        assert_eq!(
            written
                .frame(0)
                .unwrap()
                .atom_properties()
                .value(&key("atom_score"), 0)
                .unwrap(),
            Some(real(0.6))
        );
        assert_eq!(
            written
                .frame(0)
                .unwrap()
                .bond_properties()
                .value(&key("bond_score"), 0)
                .unwrap(),
            Some(real(4.0))
        );

        let independent = make_topology(true);
        let other = TrajectoryFrame::new(
            positions(&[Point3::origin(), Point3::origin()]),
            independent.bond_count(),
        );
        let mut writer = MemoryTrajectoryWriter::new(Arc::clone(&topology));
        assert!(matches!(
            writer.write_frame(other.view(&independent).unwrap()),
            Err(TrajectoryError::TopologyMismatch)
        ));
    }

    #[test]
    fn coordinate_reader_requires_explicit_order_and_matching_buffer_topology() {
        let topology = make_topology(false);
        let assertion = AtomOrderAssertion::assert_file_uses_topology_order(&topology);
        let mut reader = CoordinateFrameReader::new(
            Arc::clone(&topology),
            assertion,
            [Quantity::new(vec![Point3::new(6.0, 0.0, 0.0)], ANGSTROM)],
        )
        .unwrap();
        let mut buffer = FrameBuffer::new(Arc::clone(&topology));
        assert!(reader.read_next(&mut buffer).unwrap());
        assert!((buffer.positions().values().value()[0].x - 0.6).abs() < 1.0e-15);

        let assertion = AtomOrderAssertion::assert_file_uses_topology_order(&topology);
        let mut reader = CoordinateFrameReader::new(
            Arc::clone(&topology),
            assertion,
            [Quantity::new(vec![Point3::origin()], ANGSTROM)],
        )
        .unwrap();
        let mut wrong_buffer = FrameBuffer::new(make_topology(false));
        assert!(matches!(
            reader.read_next(&mut wrong_buffer),
            Err(TrajectoryError::TopologyMismatch)
        ));

        let independent = make_topology(false);
        let assertion = AtomOrderAssertion::assert_file_uses_topology_order(&topology);
        assert!(matches!(
            CoordinateFrameReader::new(
                independent,
                assertion,
                [Quantity::new(vec![Point3::origin()], ANGSTROM)]
            ),
            Err(TrajectoryError::TopologyMismatch)
        ));
    }
}
