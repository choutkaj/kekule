use std::sync::Arc;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{
    Properties, PropertyColumn, PropertyError, PropertyKey, PropertyTable, PropertyValue,
};
use kekule::structure::{Model, ModelView, Positions};
use kekule::topology::Topology;
use kekule::units::{
    Quantity, Unit, CANONICAL_FORCE_UNIT, CANONICAL_TIME_UNIT, CANONICAL_VELOCITY_UNIT,
};

use super::{validate_atom_count, FrameError};

#[derive(Debug, Clone, PartialEq)]
pub(super) struct DenseVectors {
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

    pub(super) fn validate_replacement<T>(&self, values: &Quantity<T>) -> Result<f64, FrameError>
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

    pub(super) fn copy_from_validated(&mut self, source: &[Vector3], factor: f64) {
        for (destination, source) in self.values.iter_mut().zip(source.iter().copied()) {
            *destination = Vector3::new(source.x * factor, source.y * factor, source.z * factor);
        }
    }
}

macro_rules! vector_array {
    ($name:ident, $unit:expr) => {
        #[doc = concat!("Dense numerical ", stringify!($name), "; contains no topology.")]
        #[derive(Debug, Clone, PartialEq)]
        pub struct $name(pub(super) DenseVectors);

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

/// One topology-free trajectory realization payload.
///
/// A frame always contains positions and may additionally contain a periodic
/// cell, velocities, forces, time, step, and realization-scoped properties.
/// Atom and bond dimensions are checked when the frame is viewed with a
/// topology or inserted into a [`super::Trajectory`].
///
/// Use [`Self::view`] for a topology-bound [`TrajectoryFrameView`], or insert
/// the frame into a trajectory and use [`super::Trajectory::frame`].
#[derive(Debug, Clone, PartialEq)]
pub struct TrajectoryFrame {
    pub(super) positions: Positions,
    pub(super) cell: Option<PeriodicCell>,
    pub(super) properties: Properties,
    pub(super) velocities: Option<Velocities>,
    pub(super) forces: Option<Forces>,
    pub(super) time: Option<Quantity<f64>>,
    pub(super) step: Option<u64>,
}

impl TrajectoryFrame {
    /// Constructs a detached frame. Empty property domains are sized on insertion.
    pub fn new(positions: Positions) -> Self {
        let properties = Properties::realization(positions.len(), 0);
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

    /// Replaces coordinates without changing the frame's atom count.
    /// Units, dimensions, and finite values are checked before mutation.
    pub fn set_positions<T: AsRef<[Point3]>>(
        &mut self,
        positions: Quantity<T>,
    ) -> Result<(), FrameError> {
        self.positions.set_all(positions)?;
        Ok(())
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

    /// Replaces detached properties. Populated atom columns must match positions;
    /// bond columns are checked against topology when inserted into a trajectory.
    pub fn set_properties(&mut self, mut properties: Properties) -> Result<(), FrameError> {
        if properties.realization_atom_properties().has_data()
            && properties.realization_atom_properties().len() != self.positions.len()
        {
            return Err(FrameError::AtomCountMismatch {
                expected: self.positions.len(),
                actual: properties.realization_atom_properties().len(),
            });
        }
        properties.normalize_realization_dimensions(
            self.positions.len(),
            properties.realization_bond_properties().len(),
        )?;
        self.properties = properties;
        Ok(())
    }

    fn insert_bond_column(
        &mut self,
        key: PropertyKey,
        column: PropertyColumn,
    ) -> Result<Option<PropertyColumn>, FrameError> {
        if !self.bond_properties().has_data() {
            let mut properties = self.properties.clone();
            properties.normalize_realization_dimensions(self.positions.len(), column.len())?;
            let previous = properties.insert_realization_bond_column(key, column)?;
            self.properties = properties;
            return Ok(previous);
        }
        Ok(self
            .properties
            .insert_realization_bond_column(key, column)?)
    }

    pub(super) fn prepare(&mut self, topology: &Arc<Topology>) -> Result<(), FrameError> {
        validate_atom_count(topology.atom_count(), self.positions.len())?;
        if self.atom_properties().has_data() {
            validate_atom_count(topology.atom_count(), self.atom_properties().len())?;
        }
        if self.bond_properties().has_data()
            && self.bond_properties().len() != topology.bond_count()
        {
            return Err(FrameError::BondCountMismatch {
                expected: topology.bond_count(),
                actual: self.bond_properties().len(),
            });
        }
        self.properties
            .normalize_realization_dimensions(topology.atom_count(), topology.bond_count())?;
        self.validate(topology)
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

    /// Validates already dimensioned state without changing empty domains.
    /// [`super::Trajectory::push`] establishes those domains before validation.
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
        self.properties.validate_realization_properties()?;
        Ok(())
    }

    /// Borrows already dimensioned state. Prefer [`super::Trajectory::frame`]
    /// to obtain a view after insertion has established empty property domains.
    pub fn view<'a>(
        &'a self,
        topology: &'a Arc<Topology>,
    ) -> Result<TrajectoryFrameView<'a>, FrameError> {
        self.validate(topology)?;
        Ok(self.validated_view(topology))
    }

    /// Only owners that already enforce all publication invariants may use this.
    pub(super) fn validated_view<'a>(
        &'a self,
        topology: &'a Arc<Topology>,
    ) -> TrajectoryFrameView<'a> {
        TrajectoryFrameView {
            topology,
            positions: &self.positions,
            cell: self.cell.as_ref(),
            properties: &self.properties,
            velocities: self.velocities.as_ref().map(Velocities::values),
            forces: self.forces.as_ref().map(Forces::values),
            time: self.time,
            step: self.step,
        }
    }
}

/// Borrowed topology-bound trajectory frame state.
///
/// The view retains exact topology identity and can be projected to
/// [`kekule::structure::ModelView`] with [`Self::as_model`] without allocating
/// or copying coordinate arrays.
#[derive(Debug, Clone, Copy)]
pub struct TrajectoryFrameView<'a> {
    pub(super) topology: &'a Arc<Topology>,
    pub(super) positions: &'a Positions,
    pub(super) cell: Option<&'a PeriodicCell>,
    pub(super) properties: &'a Properties,
    pub(super) velocities: Option<Quantity<&'a [Vector3]>>,
    pub(super) forces: Option<Quantity<&'a [Vector3]>>,
    pub(super) time: Option<Quantity<f64>>,
    pub(super) step: Option<u64>,
}

impl<'a> TrajectoryFrameView<'a> {
    /// Copies the complete realization, including velocities, forces, time, step,
    /// and every property. The detached payload carries no topology.
    pub fn to_frame(self) -> TrajectoryFrame {
        self.with_positions(self.positions.clone())
    }

    /// Internal publication of validated coordinate-only transformations.
    pub(crate) fn with_positions(self, positions: Positions) -> TrajectoryFrame {
        assert_eq!(
            positions.len(),
            self.positions.len(),
            "coordinate transformation preserves atom count"
        );
        TrajectoryFrame {
            positions,
            cell: self.cell.copied(),
            properties: self.properties.clone(),
            velocities: self.velocities.map(|values| {
                Velocities(DenseVectors {
                    values: values.value().to_vec(),
                    unit: values.unit(),
                })
            }),
            forces: self.forces.map(|values| {
                Forces(DenseVectors {
                    values: values.value().to_vec(),
                    unit: values.unit(),
                })
            }),
            time: self.time,
            step: self.step,
        }
    }

    pub fn topology(self) -> &'a Topology {
        self.topology
    }

    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.topology)
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
        index: usize,
        key: &PropertyKey,
    ) -> Result<Option<PropertyValue>, FrameError> {
        Ok(self.atom_properties().value(key, index)?)
    }

    pub fn bond_property(
        self,
        index: usize,
        key: &PropertyKey,
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

    /// Projects this frame into zero-copy borrowed model semantics.
    pub fn as_model(self) -> ModelView<'a> {
        ModelView::new(self.topology, self.positions, self.cell, self.properties)
            .expect("trajectory frame view has validated topology")
    }

    /// Materializes model-relevant frame state as an owned model.
    ///
    /// Velocities, forces, time, and step are trajectory-specific and are not
    /// represented by the returned model.
    pub fn to_model(self) -> Model {
        self.as_model().to_model()
    }
}
