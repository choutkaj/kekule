//! Fixed-topology trajectory frames, reusable buffers, in-memory storage, and
//! streaming reader/writer contracts.

use std::{fmt, io, sync::Arc};

use kekule::core::PropMap;
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::structure::{
    AtomData, AtomDataError, BondData, BondDataError, ModelError, ModelView, PositionError,
    Positions,
};
use kekule::topology::{InstanceAtomId, Topology};
use kekule::units::{
    Quantity, Unit, UnitError, MODEL_FORCE_UNIT, MODEL_TIME_UNIT, MODEL_VELOCITY_UNIT,
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

            pub fn set_all<T>(&mut self, values: Quantity<T>) -> Result<(), FrameError>
            where
                T: AsRef<[Vector3]>,
            {
                self.0.set_all(values)
            }
        }
    };
}

vector_array!(Velocities, MODEL_VELOCITY_UNIT);
vector_array!(Forces, MODEL_FORCE_UNIT);

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
    pub fn new(positions: Positions, bond_count: usize) -> Self {
        let atom_data = AtomData::new(positions.len());
        let bond_data = BondData::new(bond_count);
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
        if atom_data.len() != self.positions.len() {
            return Err(FrameError::AtomCountMismatch {
                expected: self.positions.len(),
                actual: atom_data.len(),
            });
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
        if bond_data.len() != self.bond_data.len() {
            return Err(FrameError::BondCountMismatch {
                expected: self.bond_data.len(),
                actual: bond_data.len(),
            });
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
                let time = time.to_unit(MODEL_TIME_UNIT)?;
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
        validate_atom_count(topology.atom_count(), self.atom_data.len())?;
        if self.bond_data.len() != topology.bond_count() {
            return Err(FrameError::BondCountMismatch {
                expected: topology.bond_count(),
                actual: self.bond_data.len(),
            });
        }
        if let Some(values) = &self.velocities {
            validate_atom_count(topology.atom_count(), values.len())?;
        }
        if let Some(values) = &self.forces {
            validate_atom_count(topology.atom_count(), values.len())?;
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
/// including atom and bond data dimensions, before changing the destination,
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
    atom_data: Option<&'a AtomData>,
    bond_data: Option<&'a BondData>,
    props: Option<&'a PropMap>,
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
            atom_data: None,
            bond_data: None,
            props: None,
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

/// Reusable caller-owned frame storage owning one exact topology context.
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
            positions: Positions::zeros(topology.atom_count()),
            cell: None,
            atom_data: AtomData::new(topology.atom_count()),
            bond_data: BondData::new(topology.bond_count()),
            velocities: Velocities::zeros(topology.atom_count()),
            has_velocities: false,
            forces: Forces::zeros(topology.atom_count()),
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
                let time = time.to_unit(MODEL_TIME_UNIT)?;
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
        if atom_data.len() != self.topology.atom_count() {
            return Err(FrameError::AtomCountMismatch {
                expected: self.topology.atom_count(),
                actual: atom_data.len(),
            });
        }
        self.atom_data = atom_data;
        Ok(())
    }

    pub fn set_bond_data(&mut self, bond_data: BondData) -> Result<(), FrameError> {
        if bond_data.len() != self.topology.bond_count() {
            return Err(FrameError::BondCountMismatch {
                expected: self.topology.bond_count(),
                actual: bond_data.len(),
            });
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
        self.atom_data = AtomData::new(self.topology.atom_count());
        self.bond_data = BondData::new(self.topology.bond_count());
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
    /// All count, unit, finite-value, atom-data, bond-data, and optional-array
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
                let time = time.to_unit(MODEL_TIME_UNIT)?;
                if !time.value().is_finite() {
                    return Err(FrameError::NonFiniteTime);
                }
                Ok(time)
            })
            .transpose()?;
        if let Some(atom_data) = data.atom_data {
            validate_atom_count(self.topology.atom_count(), atom_data.len())?;
        }
        if let Some(bond_data) = data.bond_data {
            if bond_data.len() != self.topology.bond_count() {
                return Err(FrameError::BondCountMismatch {
                    expected: self.topology.bond_count(),
                    actual: bond_data.len(),
                });
            }
        }
        let atom_data = data
            .atom_data
            .cloned()
            .unwrap_or_else(|| AtomData::new(self.topology.atom_count()));
        let bond_data = data
            .bond_data
            .cloned()
            .unwrap_or_else(|| BondData::new(self.topology.bond_count()));
        let props = data.props.cloned().unwrap_or_default();

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
        self.atom_data = atom_data;
        self.bond_data = bond_data;
        self.props = props;
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
        owned.atom_data = frame.atom_data.clone();
        owned.bond_data = frame.bond_data.clone();
        owned.velocities = frame.velocities.map(Velocities::new).transpose()?;
        owned.forces = frame.forces.map(Forces::new).transpose()?;
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
            kekule::units::MODEL_LENGTH_UNIT,
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
    AtomCountMismatch { expected: usize, actual: usize },
    BondCountMismatch { expected: usize, actual: usize },
    NonFiniteVector { index: usize },
    NonFiniteTime,
    Position(PositionError),
    AtomData(AtomDataError),
    BondData(BondDataError),
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
            Self::BondCountMismatch { expected, actual } => write!(
                formatter,
                "trajectory array requires {expected} bonds, but received {actual}"
            ),
            Self::NonFiniteVector { index } => {
                write!(formatter, "trajectory vector at {index} is not finite")
            }
            Self::NonFiniteTime => formatter.write_str("trajectory time must be finite"),
            Self::Position(error) => write!(formatter, "invalid frame positions: {error}"),
            Self::AtomData(error) => {
                write!(formatter, "invalid frame atom data: {error}")
            }
            Self::BondData(error) => {
                write!(formatter, "invalid frame bond data: {error}")
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

impl From<BondDataError> for FrameError {
    fn from(error: BondDataError) -> Self {
        Self::BondData(error)
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
    use kekule::core::{Atom, BondOrder, Element, MoleculeEditor, PropValue};
    use kekule::geometry::Point3;
    use kekule::topology::TopologyBuilder;
    use kekule::units::{
        ANGSTROM, DIMENSIONLESS, KELVIN, KILOJOULE_PER_MOLE, NANOMETER, PICOSECOND,
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
        Positions::new(Quantity::new(values, ANGSTROM)).unwrap()
    }

    #[test]
    fn vector_arrays_are_topology_free_unit_aware_and_equal_by_values() {
        let vectors = [Vector3::new(1.0, 2.0, 3.0)];
        let velocities = Velocities::new(Quantity::new(vectors, NANOMETER / PICOSECOND)).unwrap();
        let same = Velocities::new(Quantity::new(vectors, MODEL_VELOCITY_UNIT)).unwrap();
        let converted_velocity = velocities.values().value()[0];
        assert!((converted_velocity.x - 10.0).abs() < 1.0e-12);
        assert!((converted_velocity.y - 20.0).abs() < 1.0e-12);
        assert!((converted_velocity.z - 30.0).abs() < 1.0e-12);
        assert_ne!(velocities, same);
        assert_eq!(velocities.len(), 1);
        assert!(!velocities.is_empty());

        let forces = Forces::new(Quantity::new(vectors, KILOJOULE_PER_MOLE / NANOMETER)).unwrap();
        let converted_force = forces.values().value()[0];
        assert!((converted_force.x - 0.1).abs() < 1.0e-12);
        assert!((converted_force.y - 0.2).abs() < 1.0e-12);
        assert!((converted_force.z - 0.3).abs() < 1.0e-12);
        assert!(matches!(
            Velocities::new(Quantity::new(
                [Vector3::new(f64::NAN, 0.0, 0.0)],
                MODEL_VELOCITY_UNIT
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
            frame.set_atom_data(AtomData::new(1)),
            Err(FrameError::AtomCountMismatch { .. })
        ));
        assert!(matches!(
            frame.set_bond_data(BondData::new(0)),
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
            .atom_data_mut()
            .set_occupancy_at(0, Some(0.8))
            .unwrap();
        frame
            .set_time(Some(Quantity::new(2.5, PICOSECOND)))
            .unwrap();
        frame.set_step(Some(7));

        let view = frame.view(&topology).unwrap();
        assert_eq!(view.model_view().position(atom).unwrap().value().x, 3.0);
        assert_eq!(view.model_view().occupancy(atom).unwrap(), Some(0.8));
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
    fn frame_buffer_publication_is_transactional_reuses_arrays_and_clears_optionals() {
        let topology = make_topology(true);
        let mut buffer = FrameBuffer::new(Arc::clone(&topology));
        let position_ptr = buffer.positions().values().value().as_ptr();
        let velocity_ptr = buffer.velocities.values().value().as_ptr();
        let force_ptr = buffer.forces.values().value().as_ptr();
        let points = [Point3::new(1.0, 0.0, 0.0), Point3::new(2.0, 0.0, 0.0)];
        let vectors = [Vector3::new(1.0, 0.0, 0.0); 2];
        let mut atom_data = AtomData::new(topology.atom_count());
        atom_data.set_occupancy_at(0, Some(0.75)).unwrap();
        let mut bond_data = BondData::new(topology.bond_count());
        bond_data
            .set_property("score", Quantity::new([Some(2.0)], DIMENSIONLESS))
            .unwrap();
        let mut props = PropMap::new();
        props.insert("codec:value".into(), PropValue::Int(7));

        buffer
            .replace_from_data(
                FrameBufferData::new(Quantity::new(points.as_slice(), ANGSTROM))
                    .with_velocities(Quantity::new(vectors.as_slice(), MODEL_VELOCITY_UNIT))
                    .with_forces(Quantity::new(vectors.as_slice(), MODEL_FORCE_UNIT))
                    .with_time(Quantity::new(1.0, PICOSECOND))
                    .with_step(4)
                    .with_atom_data(&atom_data)
                    .with_bond_data(&bond_data)
                    .with_props(&props),
            )
            .unwrap();
        assert!(buffer.frame_view().velocities().is_some());
        assert!(buffer.frame_view().forces().is_some());
        assert_eq!(buffer.positions().values().value().as_ptr(), position_ptr);
        assert_eq!(buffer.velocities.values().value().as_ptr(), velocity_ptr);
        assert_eq!(buffer.forces.values().value().as_ptr(), force_ptr);
        assert!(buffer.atom_data().has_data());
        assert!(buffer.bond_data().has_data());
        assert_eq!(buffer.props().get("codec:value"), Some(&PropValue::Int(7)));

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
        assert!(!buffer.atom_data().has_data());
        assert!(!buffer.bond_data().has_data());
        assert!(buffer.props().is_empty());
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
            .atom_data_mut()
            .set_occupancy_at(0, Some(0.6))
            .unwrap();
        frame
            .bond_data_mut()
            .set_property("score", Quantity::new([Some(4.0)], DIMENSIONLESS))
            .unwrap();
        trajectory.push(frame).unwrap();

        let mut reader = MemoryTrajectoryReader::new(&trajectory);
        let mut buffer = FrameBuffer::new(Arc::clone(&topology));
        assert!(reader.read_next(&mut buffer).unwrap());
        assert_eq!(buffer.positions().values().value()[0].x, 5.0);
        assert_eq!(buffer.frame_view().step(), Some(9));
        assert_eq!(buffer.atom_data().occupancy_at(0).unwrap(), Some(0.6));
        assert_eq!(
            buffer.bond_data().property_value_at("score", 0).unwrap(),
            Some(Quantity::new(4.0, DIMENSIONLESS))
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
                .atom_data()
                .occupancy_at(0)
                .unwrap(),
            Some(0.6)
        );
        assert_eq!(
            written
                .frame(0)
                .unwrap()
                .bond_data()
                .property_value_at("score", 0)
                .unwrap(),
            Some(Quantity::new(4.0, DIMENSIONLESS))
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
        assert_eq!(buffer.positions().values().value()[0].x, 6.0);

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
