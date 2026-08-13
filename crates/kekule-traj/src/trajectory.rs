//! Fixed-topology trajectory frames, reusable buffers, in-memory storage, and
//! streaming reader/writer contracts.

use std::{fmt, io, sync::Arc};

use kekule::core::PropMap;
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::structure::{
    remap::{dense_atom_values, validate_complete_atom_mapping},
    AtomData, AtomDataError, BondData, ModelError, ModelView, PositionError, Positions,
    TopologyRemapError,
};
use kekule::topology::{InstanceAtomId, Topology, TopologyMapping};
use kekule::units::{
    Quantity, Unit, UnitError, MODEL_FORCE_UNIT, MODEL_TIME_UNIT, MODEL_VELOCITY_UNIT,
};

#[derive(Debug, Clone)]
struct TopologyVectors {
    topology: Arc<Topology>,
    values: Vec<Vector3>,
    unit: Unit,
}

impl PartialEq for TopologyVectors {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.topology, &other.topology)
            && self.values == other.values
            && self.unit == other.unit
    }
}

impl TopologyVectors {
    fn new<T>(topology: &Arc<Topology>, values: Quantity<T>, unit: Unit) -> Result<Self, FrameError>
    where
        T: AsRef<[Vector3]>,
    {
        let factor = values.unit().conversion_factor_to(unit)?;
        let source = values.value().as_ref();
        if source.len() != topology.atom_count() {
            return Err(FrameError::AtomCountMismatch {
                expected: topology.atom_count(),
                actual: source.len(),
            });
        }
        let values = source
            .iter()
            .copied()
            .enumerate()
            .map(|(index, vector)| {
                let vector = Vector3::new(vector.x * factor, vector.y * factor, vector.z * factor);
                if !vector.is_finite() {
                    return Err(FrameError::NonFiniteVector {
                        atom: topology.atom_ids()[index],
                    });
                }
                Ok(vector)
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            topology: Arc::clone(topology),
            values,
            unit,
        })
    }

    fn zeros(topology: &Arc<Topology>, unit: Unit) -> Self {
        Self {
            topology: Arc::clone(topology),
            values: vec![Vector3::zero(); topology.atom_count()],
            unit,
        }
    }

    fn is_compatible(&self, topology: &Arc<Topology>) -> bool {
        Arc::ptr_eq(&self.topology, topology)
    }

    fn topology(&self) -> &Topology {
        &self.topology
    }

    fn values(&self) -> Quantity<&[Vector3]> {
        Quantity::new(self.values.as_slice(), self.unit)
    }

    fn set_all<T>(
        &mut self,
        topology: &Arc<Topology>,
        values: Quantity<T>,
    ) -> Result<(), FrameError>
    where
        T: AsRef<[Vector3]>,
    {
        let factor = self.validate_replacement(topology, &values)?;
        self.copy_from_validated(values.value().as_ref(), factor);
        Ok(())
    }

    fn validate_replacement<T>(
        &self,
        topology: &Arc<Topology>,
        values: &Quantity<T>,
    ) -> Result<f64, FrameError>
    where
        T: AsRef<[Vector3]>,
    {
        if !self.is_compatible(topology) {
            return Err(FrameError::TopologyMismatch);
        }
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
                return Err(FrameError::NonFiniteVector {
                    atom: topology.atom_ids()[index],
                });
            }
        }
        Ok(factor)
    }

    fn copy_from_validated(&mut self, source: &[Vector3], factor: f64) {
        for (destination, source) in self.values.iter_mut().zip(source.iter().copied()) {
            *destination = Vector3::new(source.x * factor, source.y * factor, source.z * factor);
        }
    }

    fn remap_to(
        &self,
        source: &Arc<Topology>,
        target: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<Self, TrajectoryRemapError> {
        if !self.is_compatible(source) {
            return Err(TrajectoryRemapError::SourceFrameTopologyMismatch);
        }
        Ok(Self {
            topology: Arc::clone(target),
            values: dense_atom_values(&self.values, source, target, mapping)?,
            unit: self.unit,
        })
    }

    fn copy_remapped_from_validated(
        &mut self,
        source: Quantity<&[Vector3]>,
        mapping: &TopologyMapping,
        factor: f64,
    ) {
        for (source_index, target_index) in mapping.atom_index_pairs() {
            let vector = source.value()[source_index.index()];
            self.values[target_index.index()] =
                Vector3::new(vector.x * factor, vector.y * factor, vector.z * factor);
        }
    }
}

macro_rules! vector_array {
    ($name:ident, $unit:expr) => {
        #[derive(Debug, Clone, PartialEq)]
        pub struct $name(TopologyVectors);

        impl $name {
            pub fn new<T>(topology: &Arc<Topology>, values: Quantity<T>) -> Result<Self, FrameError>
            where
                T: AsRef<[Vector3]>,
            {
                Ok(Self(TopologyVectors::new(topology, values, $unit)?))
            }

            pub fn zeros(topology: &Arc<Topology>) -> Self {
                Self(TopologyVectors::zeros(topology, $unit))
            }

            pub fn is_compatible(&self, topology: &Arc<Topology>) -> bool {
                self.0.is_compatible(topology)
            }

            pub(crate) fn topology(&self) -> &Topology {
                self.0.topology()
            }

            pub fn values(&self) -> Quantity<&[Vector3]> {
                self.0.values()
            }

            pub fn set_all<T>(
                &mut self,
                topology: &Arc<Topology>,
                values: Quantity<T>,
            ) -> Result<(), FrameError>
            where
                T: AsRef<[Vector3]>,
            {
                self.0.set_all(topology, values)
            }
        }
    };
}

vector_array!(Velocities, MODEL_VELOCITY_UNIT);
vector_array!(Forces, MODEL_FORCE_UNIT);

/// One owned frame over one exact topology.
#[derive(Debug, Clone, PartialEq)]
pub struct TrajectoryFrame {
    positions: Positions,
    cell: Option<PeriodicCell>,
    atom_data: AtomData,
    bond_data: BondData,
    velocities: Option<Velocities>,
    forces: Option<Forces>,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
    props: PropMap,
}

impl TrajectoryFrame {
    pub fn new(positions: Positions) -> Self {
        let topology = positions.shared_topology();
        let atom_data = AtomData::new(&topology);
        let bond_data = BondData::new(&topology);
        Self {
            positions,
            cell: None,
            atom_data,
            bond_data,
            velocities: None,
            forces: None,
            time: None,
            step: None,
            props: PropMap::new(),
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

    pub fn atom_data(&self) -> &AtomData {
        &self.atom_data
    }

    pub fn atom_data_mut(&mut self) -> &mut AtomData {
        &mut self.atom_data
    }

    pub fn set_atom_data(&mut self, atom_data: AtomData) -> Result<(), FrameError> {
        if !std::ptr::eq(self.positions.topology(), atom_data.topology()) {
            return Err(FrameError::TopologyMismatch);
        }
        self.atom_data = atom_data;
        Ok(())
    }

    pub fn bond_data(&self) -> &BondData {
        &self.bond_data
    }

    pub fn bond_data_mut(&mut self) -> &mut BondData {
        &mut self.bond_data
    }

    pub fn set_bond_data(&mut self, bond_data: BondData) -> Result<(), FrameError> {
        if !std::ptr::eq(self.positions.topology(), bond_data.topology()) {
            return Err(FrameError::TopologyMismatch);
        }
        self.bond_data = bond_data;
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

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }

    pub fn set_velocities(&mut self, velocities: Option<Velocities>) -> Result<(), FrameError> {
        if velocities
            .as_ref()
            .is_some_and(|values| !std::ptr::eq(self.positions.topology(), values.topology()))
        {
            return Err(FrameError::TopologyMismatch);
        }
        self.velocities = velocities;
        Ok(())
    }

    pub fn set_forces(&mut self, forces: Option<Forces>) -> Result<(), FrameError> {
        if forces
            .as_ref()
            .is_some_and(|values| !std::ptr::eq(self.positions.topology(), values.topology()))
        {
            return Err(FrameError::TopologyMismatch);
        }
        self.forces = forces;
        Ok(())
    }

    pub fn set_time(&mut self, time: Option<Quantity<f64>>) -> Result<(), FrameError> {
        self.time = match time {
            Some(time) => {
                let time = time.into_unit(MODEL_TIME_UNIT)?;
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
        if !self.positions.is_compatible(topology)
            || !self.atom_data.is_compatible(topology)
            || !self.bond_data.is_compatible(topology)
            || self
                .velocities
                .as_ref()
                .is_some_and(|values| !values.is_compatible(topology))
            || self
                .forces
                .as_ref()
                .is_some_and(|values| !values.is_compatible(topology))
        {
            return Err(FrameError::TopologyMismatch);
        }
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
            atom_data: &self.atom_data,
            bond_data: &self.bond_data,
            velocities: self.velocities.as_ref().map(Velocities::values),
            forces: self.forces.as_ref().map(Forces::values),
            time: self.time,
            step: self.step,
            props: &self.props,
        })
    }

    /// Remaps this owned frame to an explicitly related target topology.
    pub fn remap_to(
        &self,
        source: &Arc<Topology>,
        target: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<Self, TrajectoryRemapError> {
        if self.validate(source).is_err() {
            return Err(TrajectoryRemapError::SourceFrameTopologyMismatch);
        }
        let positions = self.positions.remap_to(source, target, mapping)?;
        let atom_data = self.atom_data.remap_to(source, target, mapping)?;
        let bond_data = self.bond_data.remap_to(source, target, mapping)?;
        let velocities = self
            .velocities
            .as_ref()
            .map(|values| values.0.remap_to(source, target, mapping).map(Velocities))
            .transpose()?;
        let forces = self
            .forces
            .as_ref()
            .map(|values| values.0.remap_to(source, target, mapping).map(Forces))
            .transpose()?;
        Ok(Self {
            positions,
            cell: self.cell,
            atom_data,
            bond_data,
            velocities,
            forces,
            time: self.time,
            step: self.step,
            props: self.props.clone(),
        })
    }
}

/// Borrowed trajectory frame state.
#[derive(Debug, Clone, Copy)]
pub struct TrajectoryFrameView<'a> {
    topology: &'a Arc<Topology>,
    positions: &'a Positions,
    cell: Option<&'a PeriodicCell>,
    atom_data: &'a AtomData,
    bond_data: &'a BondData,
    velocities: Option<Quantity<&'a [Vector3]>>,
    forces: Option<Quantity<&'a [Vector3]>>,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
    props: &'a PropMap,
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

    pub const fn atom_data(self) -> &'a AtomData {
        self.atom_data
    }

    pub const fn bond_data(self) -> &'a BondData {
        self.bond_data
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

    pub const fn props(self) -> &'a PropMap {
        self.props
    }

    pub fn model_view(self) -> ModelView<'a> {
        ModelView::new(
            self.topology,
            self.positions,
            self.cell,
            self.atom_data,
            self.bond_data,
        )
        .expect("trajectory frame view has validated topology")
    }
}

/// Complete borrowed frame state ready for transactional publication.
///
/// A decoder builds this value only after it has read one complete frame into
/// reusable scratch. [`FrameBuffer::replace_from_data`] validates every field,
/// including topology-bound atom and bond data, before changing the destination,
/// converts units once, reuses dense-array allocations, and clears optional
/// fields omitted from this value.
#[derive(Debug, Clone, Copy)]
pub struct FrameBufferData<'a> {
    topology: &'a Arc<Topology>,
    positions: Quantity<&'a [Point3]>,
    cell: Option<PeriodicCell>,
    velocities: Option<Quantity<&'a [Vector3]>>,
    forces: Option<Quantity<&'a [Vector3]>>,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
    atom_data: Option<&'a AtomData>,
    bond_data: Option<&'a BondData>,
    props: Option<&'a PropMap>,
}

impl<'a> FrameBufferData<'a> {
    /// Starts complete frame data with required topology-bound positions.
    pub const fn new(topology: &'a Arc<Topology>, positions: Quantity<&'a [Point3]>) -> Self {
        Self {
            topology,
            positions,
            cell: None,
            velocities: None,
            forces: None,
            time: None,
            step: None,
            atom_data: None,
            bond_data: None,
            props: None,
        }
    }

    /// Borrows every field from an already validated frame view.
    pub fn from_frame_view(frame: TrajectoryFrameView<'a>) -> Self {
        Self {
            topology: frame.topology,
            positions: frame.positions.values(),
            cell: frame.cell.copied(),
            velocities: frame.velocities,
            forces: frame.forces,
            time: frame.time,
            step: frame.step,
            atom_data: Some(frame.atom_data),
            bond_data: Some(frame.bond_data),
            props: Some(frame.props),
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

    pub const fn with_atom_data(mut self, atom_data: &'a AtomData) -> Self {
        self.atom_data = Some(atom_data);
        self
    }

    pub const fn with_bond_data(mut self, bond_data: &'a BondData) -> Self {
        self.bond_data = Some(bond_data);
        self
    }

    pub const fn with_props(mut self, props: &'a PropMap) -> Self {
        self.props = Some(props);
        self
    }
}

/// Reusable caller-owned frame storage bound to one exact topology.
#[derive(Debug, Clone)]
pub struct FrameBuffer {
    topology: Arc<Topology>,
    positions: Positions,
    cell: Option<PeriodicCell>,
    atom_data: AtomData,
    bond_data: BondData,
    velocities: Velocities,
    has_velocities: bool,
    forces: Forces,
    has_forces: bool,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
    props: PropMap,
}

impl FrameBuffer {
    pub fn new(topology: Arc<Topology>) -> Self {
        Self {
            positions: Positions::zeros(&topology),
            cell: None,
            atom_data: AtomData::new(&topology),
            bond_data: BondData::new(&topology),
            velocities: Velocities::zeros(&topology),
            has_velocities: false,
            forces: Forces::zeros(&topology),
            has_forces: false,
            time: None,
            step: None,
            props: PropMap::new(),
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

    pub fn atom_data(&self) -> &AtomData {
        &self.atom_data
    }

    pub fn atom_data_mut(&mut self) -> &mut AtomData {
        &mut self.atom_data
    }

    pub fn bond_data(&self) -> &BondData {
        &self.bond_data
    }

    pub fn bond_data_mut(&mut self) -> &mut BondData {
        &mut self.bond_data
    }

    pub fn set_positions<T>(&mut self, positions: Quantity<T>) -> Result<(), FrameError>
    where
        T: AsRef<[Point3]>,
    {
        self.positions.set_all(&self.topology, positions)?;
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
                self.velocities.set_all(&self.topology, values)?;
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
                self.forces.set_all(&self.topology, values)?;
                self.has_forces = true;
            }
            None => self.has_forces = false,
        }
        Ok(())
    }

    pub fn set_time(&mut self, time: Option<Quantity<f64>>) -> Result<(), FrameError> {
        self.time = match time {
            Some(time) => {
                let time = time.into_unit(MODEL_TIME_UNIT)?;
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

    pub fn set_atom_data(&mut self, atom_data: AtomData) -> Result<(), FrameError> {
        if !atom_data.is_compatible(&self.topology) {
            return Err(FrameError::TopologyMismatch);
        }
        self.atom_data = atom_data;
        Ok(())
    }

    pub fn set_bond_data(&mut self, bond_data: BondData) -> Result<(), FrameError> {
        if !bond_data.is_compatible(&self.topology) {
            return Err(FrameError::TopologyMismatch);
        }
        self.bond_data = bond_data;
        Ok(())
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }

    /// Clears all per-frame state except positions while retaining reusable
    /// array allocations and the bound topology.
    pub fn reset_dynamic_state(&mut self) {
        self.cell = None;
        self.atom_data = AtomData::new(&self.topology);
        self.bond_data = BondData::new(&self.topology);
        self.has_velocities = false;
        self.has_forces = false;
        self.time = None;
        self.step = None;
        self.props.clear();
    }

    pub fn model_view(&self) -> ModelView<'_> {
        ModelView::new(
            &self.topology,
            &self.positions,
            self.cell.as_ref(),
            &self.atom_data,
            &self.bond_data,
        )
        .expect("frame buffer state is bound to its topology")
    }

    pub fn frame_view(&self) -> TrajectoryFrameView<'_> {
        TrajectoryFrameView {
            topology: &self.topology,
            positions: &self.positions,
            cell: self.cell.as_ref(),
            atom_data: &self.atom_data,
            bond_data: &self.bond_data,
            velocities: self.has_velocities.then(|| self.velocities.values()),
            forces: self.has_forces.then(|| self.forces.values()),
            time: self.time,
            step: self.step,
            props: &self.props,
        }
    }

    /// Replaces the complete visible frame transactionally.
    ///
    /// All topology, count, unit, finite-value, atom-data, bond-data, and
    /// optional-array validation completes before any destination field changes.
    /// Existing position, velocity, and force allocations are reused. Optional
    /// fields absent from `data`, including properties, are cleared.
    pub fn replace_from_data(&mut self, data: FrameBufferData<'_>) -> Result<(), FrameError> {
        if !Arc::ptr_eq(&self.topology, data.topology) {
            return Err(FrameError::TopologyMismatch);
        }

        self.positions
            .validate_all(&self.topology, &data.positions)?;
        let velocities = data
            .velocities
            .map(|values| {
                self.velocities
                    .0
                    .validate_replacement(&self.topology, &values)
                    .map(|factor| (values, factor))
            })
            .transpose()?;
        let forces = data
            .forces
            .map(|values| {
                self.forces
                    .0
                    .validate_replacement(&self.topology, &values)
                    .map(|factor| (values, factor))
            })
            .transpose()?;
        let time = data
            .time
            .map(|time| {
                let time = time.into_unit(MODEL_TIME_UNIT)?;
                if !time.value().is_finite() {
                    return Err(FrameError::NonFiniteTime);
                }
                Ok(time)
            })
            .transpose()?;
        if data
            .atom_data
            .is_some_and(|atom_data| !atom_data.is_compatible(&self.topology))
        {
            return Err(FrameError::TopologyMismatch);
        }
        if data
            .bond_data
            .is_some_and(|bond_data| !bond_data.is_compatible(&self.topology))
        {
            return Err(FrameError::TopologyMismatch);
        }
        let atom_data = data
            .atom_data
            .cloned()
            .unwrap_or_else(|| AtomData::new(&self.topology));
        let bond_data = data
            .bond_data
            .cloned()
            .unwrap_or_else(|| BondData::new(&self.topology));
        let props = data.props.cloned().unwrap_or_default();

        self.positions.set_all(&self.topology, data.positions)?;
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
        self.atom_data = atom_data;
        self.bond_data = bond_data;
        self.props = props;
        Ok(())
    }

    pub fn copy_from(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), FrameError> {
        self.replace_from_data(FrameBufferData::from_frame_view(frame))
    }

    /// Copies a borrowed source frame through explicit topology lineage.
    ///
    /// All fallible validation and metadata staging completes before
    /// destination-visible state changes. Position, velocity, and force
    /// allocations owned by this buffer are reused.
    pub fn copy_remapped_from(
        &mut self,
        frame: TrajectoryFrameView<'_>,
        mapping: &TopologyMapping,
    ) -> Result<(), TrajectoryRemapError> {
        if !mapping.is_target(&self.topology) {
            return Err(TrajectoryRemapError::IncompatibleDestinationBuffer);
        }
        if !frame.positions.is_compatible(frame.topology)
            || !frame.atom_data.is_compatible(frame.topology)
            || !frame.bond_data.is_compatible(frame.topology)
        {
            return Err(TrajectoryRemapError::SourceFrameTopologyMismatch);
        }
        validate_complete_atom_mapping(frame.topology, &self.topology, mapping)?;
        validate_borrowed_array(frame.velocities, frame.topology)?;
        validate_borrowed_array(frame.forces, frame.topology)?;

        let velocity_factor = frame
            .velocities
            .map(|values| values.unit().conversion_factor_to(MODEL_VELOCITY_UNIT))
            .transpose()?;
        let force_factor = frame
            .forces
            .map(|values| values.unit().conversion_factor_to(MODEL_FORCE_UNIT))
            .transpose()?;
        let atom_data = frame
            .atom_data
            .remap_to(frame.topology, &self.topology, mapping)?;
        let bond_data = frame
            .bond_data
            .remap_to(frame.topology, &self.topology, mapping)?;
        let props = frame.props.clone();

        self.positions.copy_remapped_from(
            frame.positions,
            frame.topology,
            &self.topology,
            mapping,
        )?;
        self.cell = frame.cell.copied();
        match (frame.velocities, velocity_factor) {
            (Some(values), Some(factor)) => {
                self.velocities
                    .0
                    .copy_remapped_from_validated(values, mapping, factor);
                self.has_velocities = true;
            }
            (None, None) => self.has_velocities = false,
            _ => unreachable!("velocity factor is staged with its source array"),
        }
        match (frame.forces, force_factor) {
            (Some(values), Some(factor)) => {
                self.forces
                    .0
                    .copy_remapped_from_validated(values, mapping, factor);
                self.has_forces = true;
            }
            (None, None) => self.has_forces = false,
            _ => unreachable!("force factor is staged with its source array"),
        }
        self.time = frame.time;
        self.step = frame.step;
        self.atom_data = atom_data;
        self.bond_data = bond_data;
        self.props = props;
        Ok(())
    }
}

fn validate_borrowed_array(
    values: Option<Quantity<&[Vector3]>>,
    source: &Topology,
) -> Result<(), TrajectoryRemapError> {
    if let Some(values) = values {
        if values.value().len() != source.atom_count() {
            return Err(TopologyRemapError::SourceAtomCountMismatch {
                expected: source.atom_count(),
                actual: values.value().len(),
            }
            .into());
        }
    }
    Ok(())
}

/// Deliberately loaded finite in-memory trajectory.
#[derive(Debug, Clone)]
pub struct Trajectory {
    topology: Arc<Topology>,
    frames: Vec<TrajectoryFrame>,
}

impl Trajectory {
    pub fn new(topology: Arc<Topology>) -> Self {
        Self {
            topology,
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
            let value = time.into_value();
            if previous.is_some_and(|previous| value < previous) {
                return Err(TrajectoryError::NonMonotonicTime { frame: index });
            }
            previous = Some(value);
        }
        Ok(())
    }

    /// Remaps every frame to one exact target topology while preserving order
    /// and complete frame state.
    pub fn remap_to(
        &self,
        target: &Arc<Topology>,
        mapping: &TopologyMapping,
    ) -> Result<Self, TrajectoryRemapError> {
        if !mapping.is_source(&self.topology) {
            return Err(TopologyRemapError::MappingSourceMismatch.into());
        }
        if !mapping.is_target(target) {
            return Err(TopologyRemapError::MappingTargetMismatch.into());
        }
        let mut frames = Vec::with_capacity(self.frames.len());
        for (frame_index, frame) in self.frames.iter().enumerate() {
            frames.push(
                frame
                    .remap_to(&self.topology, target, mapping)
                    .map_err(|error| TrajectoryRemapError::Frame {
                        frame: frame_index,
                        error: Box::new(error),
                    })?,
            );
        }
        Ok(Self {
            topology: Arc::clone(target),
            frames,
        })
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

    pub fn into_trajectory(self) -> Trajectory {
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
        let positions = Positions::new(&self.trajectory.topology, frame.positions.values())?;
        let mut owned = TrajectoryFrame::new(positions);
        owned.cell = frame.cell.copied();
        owned.atom_data = frame.atom_data.clone();
        owned.bond_data = frame.bond_data.clone();
        owned.velocities = frame
            .velocities
            .map(|values| Velocities::new(&self.trajectory.topology, values))
            .transpose()?;
        owned.forces = frame
            .forces
            .map(|values| Forces::new(&self.trajectory.topology, values))
            .transpose()?;
        owned.time = frame.time;
        owned.step = frame.step;
        owned.props = frame.props.clone();
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
                Positions::new(&topology, frame)
                    .map(|positions| positions.values().value().to_vec())
                    .map_err(TrajectoryError::Position)
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
            kekule::units::MODEL_LENGTH_UNIT,
        ))?;
        destination.reset_dynamic_state();
        self.cursor += 1;
        Ok(true)
    }
}

/// Failure to remap owned or borrowed fixed-topology trajectory state.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TrajectoryRemapError {
    /// An owned or borrowed source frame is not bound to the supplied source.
    SourceFrameTopologyMismatch,
    /// A caller-owned destination buffer is not bound to the mapping target.
    IncompatibleDestinationBuffer,
    /// Complete topology-bound state could not be remapped.
    Topology(TopologyRemapError),
    /// One trajectory frame could not be remapped.
    Frame {
        frame: usize,
        error: Box<TrajectoryRemapError>,
    },
    /// A borrowed vector array has an incompatible unit.
    Unit(UnitError),
}

impl fmt::Display for TrajectoryRemapError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::SourceFrameTopologyMismatch => {
                formatter.write_str("source frame does not belong to the supplied topology")
            }
            Self::IncompatibleDestinationBuffer => {
                formatter.write_str("destination buffer does not match the mapping target")
            }
            Self::Topology(error) => {
                write!(formatter, "cannot remap topology-bound state: {error}")
            }
            Self::Frame { frame, error } => {
                write!(formatter, "cannot remap trajectory frame {frame}: {error}")
            }
            Self::Unit(error) => write!(formatter, "cannot remap trajectory units: {error}"),
        }
    }
}

impl std::error::Error for TrajectoryRemapError {}

impl From<TopologyRemapError> for TrajectoryRemapError {
    fn from(error: TopologyRemapError) -> Self {
        Self::Topology(error)
    }
}

impl From<UnitError> for TrajectoryRemapError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum FrameError {
    TopologyMismatch,
    AtomCountMismatch { expected: usize, actual: usize },
    NonFiniteVector { atom: InstanceAtomId },
    NonFiniteTime,
    Position(PositionError),
    AtomData(AtomDataError),
    Unit(UnitError),
}

impl fmt::Display for FrameError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyMismatch => {
                formatter.write_str("trajectory frame belongs to a different topology")
            }
            Self::AtomCountMismatch { expected, actual } => write!(
                formatter,
                "trajectory array requires {expected} atoms, but received {actual}"
            ),
            Self::NonFiniteVector { atom } => {
                write!(formatter, "trajectory vector for atom {atom} is not finite")
            }
            Self::NonFiniteTime => formatter.write_str("trajectory time must be finite"),
            Self::Position(error) => write!(formatter, "invalid frame positions: {error}"),
            Self::AtomData(error) => {
                write!(formatter, "invalid frame atom data: {error}")
            }
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

impl From<AtomDataError> for FrameError {
    fn from(error: AtomDataError) -> Self {
        Self::AtomData(error)
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
    use kekule::core::{Atom, BondOrder, Element, Molecule, PropValue};
    use kekule::small::SmallMolecule;
    use kekule::topology::{
        transform::retain_instances, MoleculeInstanceMetadata, TopologyBuilder,
    };
    use kekule::units::{ANGSTROM, KELVIN, KILOJOULE_PER_MOLE, PICOSECOND};

    fn one_atom_topology() -> Arc<Topology> {
        let mut graph = Molecule::builder();
        graph
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .expect("atom identifier capacity");
        let molecule = SmallMolecule::from_graph(graph.build().unwrap());
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_small_molecule_definition(&molecule).unwrap();
        builder
            .add_instance(definition, MoleculeInstanceMetadata::default())
            .unwrap();
        Arc::new(builder.build().unwrap())
    }

    fn one_bond_topology() -> Arc<Topology> {
        let mut graph = Molecule::builder();
        let carbon = graph
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .unwrap();
        let oxygen = graph
            .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
            .unwrap();
        graph.add_bond(carbon, oxygen, BondOrder::Single).unwrap();
        let molecule = SmallMolecule::from_graph(graph.build().unwrap());
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_small_molecule_definition(&molecule).unwrap();
        builder
            .add_instance(definition, MoleculeInstanceMetadata::default())
            .unwrap();
        Arc::new(builder.build().unwrap())
    }

    fn positions(topology: &Arc<Topology>, x: f64) -> Positions {
        Positions::new(
            topology,
            Quantity::new(vec![Point3::new(x, 0.0, 0.0)], ANGSTROM),
        )
        .unwrap()
    }

    fn frame(topology: &Arc<Topology>, x: f64, time: f64) -> TrajectoryFrame {
        let mut frame = TrajectoryFrame::new(positions(topology, x));
        frame
            .set_velocities(Some(
                Velocities::new(
                    topology,
                    Quantity::new(vec![Vector3::new(x + 1.0, 0.0, 0.0)], MODEL_VELOCITY_UNIT),
                )
                .unwrap(),
            ))
            .unwrap();
        frame
            .set_time(Some(Quantity::new(time, PICOSECOND)))
            .unwrap();
        frame.set_step(Some(time as u64));
        frame
    }

    fn assert_same_buffer_state(actual: &FrameBuffer, expected: &FrameBuffer) {
        assert!(Arc::ptr_eq(&actual.topology, &expected.topology));
        assert_eq!(actual.positions, expected.positions);
        assert_eq!(actual.cell, expected.cell);
        assert_eq!(actual.atom_data, expected.atom_data);
        assert_eq!(actual.bond_data, expected.bond_data);
        assert_eq!(actual.velocities, expected.velocities);
        assert_eq!(actual.has_velocities, expected.has_velocities);
        assert_eq!(actual.forces, expected.forces);
        assert_eq!(actual.has_forces, expected.has_forces);
        assert_eq!(actual.time, expected.time);
        assert_eq!(actual.step, expected.step);
        assert_eq!(actual.props, expected.props);
    }

    #[test]
    fn frames_validate_optional_arrays_and_exact_topology() {
        let topology = one_atom_topology();
        assert_eq!(
            Velocities::new(
                &topology,
                Quantity::new(Vec::<Vector3>::new(), MODEL_VELOCITY_UNIT)
            ),
            Err(FrameError::AtomCountMismatch {
                expected: 1,
                actual: 0
            })
        );
        assert!(matches!(
            Forces::new(
                &topology,
                Quantity::new(vec![Vector3::new(f64::NAN, 0.0, 0.0)], MODEL_FORCE_UNIT)
            ),
            Err(FrameError::NonFiniteVector { .. })
        ));

        let independent = one_atom_topology();
        let mut frame = TrajectoryFrame::new(positions(&topology, 0.0));
        assert_eq!(
            frame.set_velocities(Some(Velocities::zeros(&independent))),
            Err(FrameError::TopologyMismatch)
        );
        assert_eq!(
            frame.set_atom_data(AtomData::new(&independent)),
            Err(FrameError::TopologyMismatch)
        );
        assert_eq!(
            frame.set_bond_data(BondData::new(&independent)),
            Err(FrameError::TopologyMismatch)
        );
        assert_eq!(
            frame.set_time(Some(Quantity::new(f64::INFINITY, PICOSECOND))),
            Err(FrameError::NonFiniteTime)
        );
    }

    #[test]
    fn bond_data_survives_frames_views_memory_streaming_and_buffer_reset() {
        let topology = one_bond_topology();
        let positions = Positions::new(
            &topology,
            Quantity::new(
                vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
                ANGSTROM,
            ),
        )
        .unwrap();
        let entropy_unit = KILOJOULE_PER_MOLE / KELVIN;
        let mut frame = TrajectoryFrame::new(positions);
        frame
            .bond_data_mut()
            .set_property(
                "conformational_entropy",
                Quantity::new(vec![Some(0.025)], entropy_unit),
            )
            .unwrap();
        let frame_view = frame.view(&topology).unwrap();
        assert!(std::ptr::eq(frame.bond_data(), frame_view.bond_data()));
        assert!(std::ptr::eq(
            frame_view.bond_data(),
            frame_view.model_view().bond_data()
        ));

        let trajectory = Trajectory::from_frames(Arc::clone(&topology), [frame.clone()]).unwrap();
        let mut writer = MemoryTrajectoryWriter::new(Arc::clone(&topology));
        writer
            .write_frame(trajectory.frames().next().unwrap())
            .unwrap();
        let written = writer.into_trajectory();
        assert_eq!(
            written
                .frame(0)
                .unwrap()
                .bond_data()
                .property("conformational_entropy")
                .unwrap()
                .unwrap()
                .value(),
            &[Some(0.025)]
        );

        let target = one_bond_topology();
        let mapping = TopologyMapping::between_identical_layouts(&topology, &target).unwrap();
        let remapped = frame.remap_to(&topology, &target, &mapping).unwrap();
        assert_eq!(
            remapped
                .bond_data()
                .property("conformational_entropy")
                .unwrap()
                .unwrap()
                .value(),
            &[Some(0.025)]
        );
        let mut remapped_buffer = FrameBuffer::new(Arc::clone(&target));
        remapped_buffer
            .copy_remapped_from(frame_view, &mapping)
            .unwrap();
        assert_eq!(
            remapped_buffer
                .bond_data()
                .property("conformational_entropy")
                .unwrap()
                .unwrap()
                .value(),
            &[Some(0.025)]
        );

        let mut buffer = FrameBuffer::new(Arc::clone(&topology));
        buffer.copy_from(frame_view).unwrap();
        assert_eq!(
            buffer
                .bond_data()
                .property("conformational_entropy")
                .unwrap()
                .unwrap()
                .value(),
            &[Some(0.025)]
        );
        assert!(std::ptr::eq(
            buffer.bond_data(),
            buffer.model_view().bond_data()
        ));

        let before = buffer.clone();
        let independent = one_bond_topology();
        let wrong_bond_data = BondData::new(&independent);
        assert_eq!(
            buffer.replace_from_data(
                FrameBufferData::new(&topology, frame.positions().values())
                    .with_bond_data(&wrong_bond_data)
            ),
            Err(FrameError::TopologyMismatch)
        );
        assert_same_buffer_state(&buffer, &before);

        buffer.reset_dynamic_state();
        assert!(buffer.bond_data().is_empty());
        assert!(buffer.frame_view().bond_data().is_empty());
    }

    #[test]
    fn reusable_buffer_streams_and_seeks_without_reallocating_positions() {
        let topology = one_atom_topology();
        let mut first = frame(&topology, 0.0, 0.0);
        first
            .atom_data
            .set_occupancy(&topology, topology.atom_ids()[0], Some(0.5))
            .unwrap();
        let cell = kekule::geometry::PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(10.0, 10.0, 10.0), ANGSTROM),
            [true; 3],
        )
        .unwrap();
        let second_positions = Positions::new(
            &topology,
            Quantity::new(vec![Point3::new(2.0, 0.0, 0.0)], ANGSTROM),
        )
        .unwrap();
        let mut second = TrajectoryFrame::new(second_positions);
        second.set_cell(Some(cell));
        second
            .set_time(Some(Quantity::new(1.0, PICOSECOND)))
            .unwrap();
        second.set_step(Some(1));
        let trajectory = Trajectory::from_frames(Arc::clone(&topology), [first, second]).unwrap();
        assert_eq!(trajectory.len(), 2);
        trajectory.validate_monotonic_time(true).unwrap();

        let mut reader = MemoryTrajectoryReader::new(&trajectory);
        let mut buffer = FrameBuffer::new(Arc::clone(&topology));
        let pointer = buffer.positions().values().value().as_ptr();
        assert!(reader.read_next(&mut buffer).unwrap());
        assert_eq!(buffer.positions().values().value().as_ptr(), pointer);
        assert!(buffer.frame_view().velocities().is_some());
        assert!(!buffer.frame_view().atom_data().is_empty());
        assert!(reader.read_next(&mut buffer).unwrap());
        assert_eq!(buffer.positions().values().value().as_ptr(), pointer);
        assert_eq!(buffer.model_view().positions().values().value()[0].x, 2.0);
        assert_eq!(buffer.frame_view().cell(), Some(&cell));
        assert!(!reader.read_next(&mut buffer).unwrap());

        reader.read_frame(0, &mut buffer).unwrap();
        assert_eq!(buffer.model_view().positions().values().value()[0].x, 0.0);
        assert_eq!(reader.frame_count(), Some(2));
    }

    #[test]
    fn complete_buffer_publication_is_transactional_and_reuses_allocations() {
        let topology = one_atom_topology();
        let cell = PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(10.0, 11.0, 12.0), ANGSTROM),
            [true; 3],
        )
        .unwrap();
        let positions = [Point3::new(1.0, 2.0, 3.0)];
        let velocities = [Vector3::new(4.0, 5.0, 6.0)];
        let forces = [Vector3::new(7.0, 8.0, 9.0)];
        let mut props = PropMap::new();
        props.insert("codec:value".into(), PropValue::Int(7));

        let mut buffer = FrameBuffer::new(Arc::clone(&topology));
        buffer
            .replace_from_data(
                FrameBufferData::new(&topology, Quantity::new(&positions, ANGSTROM))
                    .with_cell(cell)
                    .with_velocities(Quantity::new(&velocities, MODEL_VELOCITY_UNIT))
                    .with_forces(Quantity::new(&forces, MODEL_FORCE_UNIT))
                    .with_time(Quantity::new(2.0, PICOSECOND))
                    .with_step(8)
                    .with_props(&props),
            )
            .unwrap();

        let position_pointer = buffer.positions.values().value().as_ptr();
        let velocity_pointer = buffer.velocities.0.values.as_ptr();
        let velocity_capacity = buffer.velocities.0.values.capacity();
        let force_pointer = buffer.forces.0.values.as_ptr();
        let force_capacity = buffer.forces.0.values.capacity();

        let before_failure = buffer.clone();
        let replacement_positions = [Point3::new(20.0, 21.0, 22.0)];
        let invalid_forces = [Vector3::new(f64::NAN, 0.0, 0.0)];
        assert!(matches!(
            buffer.replace_from_data(
                FrameBufferData::new(&topology, Quantity::new(&replacement_positions, ANGSTROM))
                    .with_velocities(Quantity::new(&velocities, MODEL_VELOCITY_UNIT))
                    .with_forces(Quantity::new(&invalid_forces, MODEL_FORCE_UNIT))
                    .with_time(Quantity::new(3.0, PICOSECOND))
            ),
            Err(FrameError::NonFiniteVector { .. })
        ));
        assert_same_buffer_state(&buffer, &before_failure);

        buffer
            .replace_from_data(FrameBufferData::new(
                &topology,
                Quantity::new(&replacement_positions, ANGSTROM),
            ))
            .unwrap();
        assert_eq!(buffer.positions.values().value().as_ptr(), position_pointer);
        assert_eq!(buffer.velocities.0.values.as_ptr(), velocity_pointer);
        assert_eq!(buffer.velocities.0.values.capacity(), velocity_capacity);
        assert_eq!(buffer.forces.0.values.as_ptr(), force_pointer);
        assert_eq!(buffer.forces.0.values.capacity(), force_capacity);
        assert_eq!(buffer.positions.values().value(), &replacement_positions);
        assert!(buffer.cell.is_none());
        assert!(buffer.frame_view().velocities().is_none());
        assert!(buffer.frame_view().forces().is_none());
        assert!(buffer.frame_view().time().is_none());
        assert!(buffer.frame_view().step().is_none());
        assert!(buffer.props().is_empty());
    }

    #[test]
    fn copy_from_uses_complete_transactional_publication() {
        let topology = one_atom_topology();
        let mut buffer = FrameBuffer::new(Arc::clone(&topology));
        buffer
            .set_positions(positions(&topology, 5.0).values())
            .unwrap();
        buffer
            .props_mut()
            .insert("old".into(), PropValue::Bool(true));
        let before = buffer.clone();

        let source = positions(&topology, 10.0);
        let atom_data = AtomData::new(&topology);
        let bond_data = BondData::new(&topology);
        let props = PropMap::new();
        let invalid = TrajectoryFrameView {
            topology: &topology,
            positions: &source,
            cell: None,
            atom_data: &atom_data,
            bond_data: &bond_data,
            velocities: None,
            forces: None,
            time: Some(Quantity::new(f64::INFINITY, PICOSECOND)),
            step: Some(9),
            props: &props,
        };
        assert_eq!(buffer.copy_from(invalid), Err(FrameError::NonFiniteTime));
        assert_same_buffer_state(&buffer, &before);
    }

    #[test]
    fn atom_order_helpers_bind_exact_shared_topology() {
        let topology = one_atom_topology();
        let independent = one_atom_topology();
        let semantic =
            AtomOrderAssertion::from_semantic_order(&topology, topology.atom_ids()).unwrap();
        assert!(semantic.is_compatible(&topology));
        assert!(!semantic.is_compatible(&independent));
        assert!(Arc::ptr_eq(&semantic.topology, &topology));
        assert_eq!(semantic.kind(), AtomOrderAssertionKind::SemanticOrder);

        let asserted = AtomOrderAssertion::assert_file_uses_topology_order(&topology);
        assert!(asserted.is_compatible(&topology));
        assert!(!asserted.is_compatible(&independent));
        assert_eq!(
            asserted.kind(),
            AtomOrderAssertionKind::DeclaredTopologyOrder
        );
    }

    #[test]
    fn file_and_codec_error_context_is_typed_and_preserved() {
        let io_context = TrajectoryIoErrorContext::new(
            TrajectoryIoOperation::ReadFrame,
            io::ErrorKind::Other,
            "disk",
        )
        .with_format(TrajectoryFormat::Dcd)
        .with_source_label("sample.dcd")
        .with_frame(4)
        .with_byte_offset(128);
        assert_eq!(io_context.operation(), TrajectoryIoOperation::ReadFrame);
        assert_eq!(io_context.format(), Some(TrajectoryFormat::Dcd));
        assert_eq!(io_context.source_label(), Some("sample.dcd"));
        assert_eq!(io_context.frame(), Some(4));
        assert_eq!(io_context.byte_offset(), Some(128));
        assert_eq!(io_context.error_kind(), io::ErrorKind::Other);

        let codec_context = TrajectoryCodecErrorContext::new(
            TrajectoryCodecErrorKind::InconsistentAtomCount,
            TrajectoryIoOperation::Index,
            Some(TrajectoryFormat::Xtc),
        )
        .with_source_label("sample.xtc")
        .with_frame(2)
        .with_byte_offset(64)
        .with_counts(10, 9)
        .with_detail("repeated atom count differs");
        assert_eq!(
            codec_context.kind(),
            TrajectoryCodecErrorKind::InconsistentAtomCount
        );
        assert_eq!(codec_context.expected(), Some(10));
        assert_eq!(codec_context.actual(), Some(9));
        assert!(TrajectoryError::from(codec_context)
            .to_string()
            .contains("expected 10, actual 9"));
    }

    #[test]
    fn memory_writer_round_trips_frames_and_rejects_other_topologies() {
        let topology = one_atom_topology();
        let trajectory =
            Trajectory::from_frames(Arc::clone(&topology), [frame(&topology, 3.0, 2.0)]).unwrap();
        let mut writer = MemoryTrajectoryWriter::new(Arc::clone(&topology));
        writer
            .write_frame(trajectory.frames().next().unwrap())
            .unwrap();
        let written = writer.into_trajectory();
        assert_eq!(written.len(), 1);
        assert_eq!(
            written
                .frames()
                .next()
                .unwrap()
                .model_view()
                .positions()
                .values()
                .value()[0]
                .x,
            3.0
        );

        let independent = one_atom_topology();
        let foreign =
            Trajectory::from_frames(Arc::clone(&independent), [frame(&independent, 4.0, 3.0)])
                .unwrap();
        let mut writer = MemoryTrajectoryWriter::new(topology);
        assert_eq!(
            writer.write_frame(foreign.frames().next().unwrap()),
            Err(TrajectoryError::TopologyMismatch)
        );
    }

    #[test]
    fn coordinate_only_reader_requires_atom_order_and_matching_buffer_topology() {
        let topology = one_atom_topology();
        assert_eq!(
            AtomOrderAssertion::new(&topology, &[]),
            Err(TrajectoryError::AtomOrderMismatch)
        );
        let assertion = AtomOrderAssertion::new(&topology, topology.atom_ids()).unwrap();
        let mut reader = CoordinateFrameReader::new(
            Arc::clone(&topology),
            assertion,
            [Quantity::new(vec![Point3::new(5.0, 0.0, 0.0)], ANGSTROM)],
        )
        .unwrap();
        let independent = one_atom_topology();
        let mut wrong_buffer = FrameBuffer::new(independent);
        assert_eq!(
            reader.read_next(&mut wrong_buffer),
            Err(TrajectoryError::TopologyMismatch)
        );

        let mut buffer = FrameBuffer::new(topology);
        assert!(reader.read_next(&mut buffer).unwrap());
        assert_eq!(buffer.model_view().positions().values().value()[0].x, 5.0);
        assert!(!reader.read_next(&mut buffer).unwrap());
    }

    #[test]
    fn coordinate_only_reader_clears_all_dynamic_state_without_reallocating_positions() {
        let topology = one_atom_topology();
        let assertion = AtomOrderAssertion::new(&topology, topology.atom_ids()).unwrap();
        let mut reader = CoordinateFrameReader::new(
            Arc::clone(&topology),
            assertion,
            [Quantity::new(vec![Point3::new(7.0, 0.0, 0.0)], ANGSTROM)],
        )
        .unwrap();
        let mut buffer = FrameBuffer::new(Arc::clone(&topology));
        let cell = kekule::geometry::PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(10.0, 10.0, 10.0), ANGSTROM),
            [true; 3],
        )
        .unwrap();
        buffer.set_cell(Some(cell));
        buffer
            .set_velocities(Some(Quantity::new(
                vec![Vector3::new(1.0, 2.0, 3.0)],
                MODEL_VELOCITY_UNIT,
            )))
            .unwrap();
        buffer
            .set_forces(Some(Quantity::new(
                vec![Vector3::new(4.0, 5.0, 6.0)],
                MODEL_FORCE_UNIT,
            )))
            .unwrap();
        buffer
            .set_time(Some(Quantity::new(2.0, PICOSECOND)))
            .unwrap();
        buffer.set_step(Some(3));
        let mut atom_data = AtomData::new(&topology);
        atom_data
            .set_b_factor(
                &topology,
                topology.atom_ids()[0],
                Some(Quantity::new(2.0, kekule::units::SQUARE_ANGSTROM)),
            )
            .unwrap();
        buffer.set_atom_data(atom_data).unwrap();
        buffer
            .props_mut()
            .insert("source".to_owned(), PropValue::String("stale".to_owned()));
        let positions_pointer = buffer.positions().values().value().as_ptr();

        assert!(reader.read_next(&mut buffer).unwrap());

        let frame = buffer.frame_view();
        assert_eq!(frame.cell(), None);
        assert_eq!(frame.velocities(), None);
        assert_eq!(frame.forces(), None);
        assert_eq!(frame.time(), None);
        assert_eq!(frame.step(), None);
        assert!(frame.atom_data().is_empty());
        assert!(frame.props().is_empty());
        assert_eq!(
            buffer.positions().values().value().as_ptr(),
            positions_pointer
        );
        assert_eq!(buffer.model_view().positions().values().value()[0].x, 7.0);
    }

    #[test]
    fn repeated_borrowed_remaps_reuse_all_dense_buffer_allocations() {
        let mut graph = Molecule::builder();
        graph
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .expect("atom identifier capacity");
        let molecule = SmallMolecule::from_graph(graph.build().unwrap());
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_small_molecule_definition(&molecule).unwrap();
        builder
            .add_instance(definition, MoleculeInstanceMetadata::default())
            .unwrap();
        let retained = builder
            .add_instance(definition, MoleculeInstanceMetadata::default())
            .unwrap();
        let source = Arc::new(builder.build().unwrap());
        let edit = retain_instances(&source, [retained]).unwrap();
        let positions = Positions::new(
            &source,
            Quantity::new(
                vec![Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)],
                ANGSTROM,
            ),
        )
        .unwrap();
        let mut frame = TrajectoryFrame::new(positions);
        frame
            .set_velocities(Some(
                Velocities::new(
                    &source,
                    Quantity::new(
                        vec![Vector3::new(3.0, 0.0, 0.0), Vector3::new(4.0, 0.0, 0.0)],
                        MODEL_VELOCITY_UNIT,
                    ),
                )
                .unwrap(),
            ))
            .unwrap();
        frame
            .set_forces(Some(
                Forces::new(
                    &source,
                    Quantity::new(
                        vec![Vector3::new(5.0, 0.0, 0.0), Vector3::new(6.0, 0.0, 0.0)],
                        MODEL_FORCE_UNIT,
                    ),
                )
                .unwrap(),
            ))
            .unwrap();

        let mut buffer = FrameBuffer::new(edit.shared_topology());
        let view = frame.view(&source).unwrap();
        buffer
            .copy_remapped_from(view, edit.mapping())
            .expect("warm-up remap");
        let position_pointer = buffer.positions.values().value().as_ptr();
        let velocity_pointer = buffer.velocities.0.values.as_ptr();
        let force_pointer = buffer.forces.0.values.as_ptr();
        let velocity_capacity = buffer.velocities.0.values.capacity();
        let force_capacity = buffer.forces.0.values.capacity();

        for _ in 0..32 {
            buffer
                .copy_remapped_from(view, edit.mapping())
                .expect("repeated remap");
            assert_eq!(buffer.positions.values().value().as_ptr(), position_pointer);
            assert_eq!(buffer.velocities.0.values.as_ptr(), velocity_pointer);
            assert_eq!(buffer.forces.0.values.as_ptr(), force_pointer);
            assert_eq!(buffer.velocities.0.values.capacity(), velocity_capacity);
            assert_eq!(buffer.forces.0.values.capacity(), force_capacity);
        }
    }
}
