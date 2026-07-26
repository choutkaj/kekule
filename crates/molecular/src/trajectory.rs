//! Fixed-topology trajectory frames, reusable buffers, in-memory storage, and
//! streaming reader/writer contracts.

use std::fmt;

use crate::core::PropMap;
use crate::geometry::{Point3, Vector3};
use crate::structure::{
    Configuration, ConfigurationView, ModelError, ModelView, ObservationError, PositionError,
    Positions, StructureObservation,
};
use crate::topology::{InstanceAtomId, Topology, TopologyIdentity};
use crate::units::{
    Quantity, Unit, UnitError, MODEL_FORCE_UNIT, MODEL_TIME_UNIT, MODEL_VELOCITY_UNIT,
};

#[derive(Debug, Clone, PartialEq)]
struct TopologyVectors {
    topology: TopologyIdentity,
    values: Vec<Vector3>,
    unit: Unit,
}

impl TopologyVectors {
    fn new<T>(topology: &Topology, values: Quantity<T>, unit: Unit) -> Result<Self, FrameError>
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
            topology: topology.identity(),
            values,
            unit,
        })
    }

    fn zeros(topology: &Topology, unit: Unit) -> Self {
        Self {
            topology: topology.identity(),
            values: vec![Vector3::zero(); topology.atom_count()],
            unit,
        }
    }

    fn is_compatible(&self, topology: &Topology) -> bool {
        self.topology == topology.identity()
    }

    fn values(&self) -> Quantity<&[Vector3]> {
        Quantity::new(self.values.as_slice(), self.unit)
    }

    fn set_all<T>(&mut self, topology: &Topology, values: Quantity<T>) -> Result<(), FrameError>
    where
        T: AsRef<[Vector3]>,
    {
        if !self.is_compatible(topology) {
            return Err(FrameError::TopologyIdentityMismatch);
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
        for (destination, source) in self.values.iter_mut().zip(source.iter().copied()) {
            *destination = Vector3::new(source.x * factor, source.y * factor, source.z * factor);
        }
        Ok(())
    }
}

macro_rules! vector_array {
    ($name:ident, $unit:expr) => {
        #[derive(Debug, Clone, PartialEq)]
        pub struct $name(TopologyVectors);

        impl $name {
            pub fn new<T>(topology: &Topology, values: Quantity<T>) -> Result<Self, FrameError>
            where
                T: AsRef<[Vector3]>,
            {
                Ok(Self(TopologyVectors::new(topology, values, $unit)?))
            }

            pub fn zeros(topology: &Topology) -> Self {
                Self(TopologyVectors::zeros(topology, $unit))
            }

            pub fn is_compatible(&self, topology: &Topology) -> bool {
                self.0.is_compatible(topology)
            }

            pub fn values(&self) -> Quantity<&[Vector3]> {
                self.0.values()
            }

            pub fn set_all<T>(
                &mut self,
                topology: &Topology,
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
    configuration: Configuration,
    velocities: Option<Velocities>,
    forces: Option<Forces>,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
    observation: Option<StructureObservation>,
    props: PropMap,
}

impl TrajectoryFrame {
    pub fn new(configuration: Configuration) -> Self {
        Self {
            configuration,
            velocities: None,
            forces: None,
            time: None,
            step: None,
            observation: None,
            props: PropMap::new(),
        }
    }

    pub fn configuration(&self) -> &Configuration {
        &self.configuration
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

    pub fn observation(&self) -> Option<&StructureObservation> {
        self.observation.as_ref()
    }

    pub fn props(&self) -> &PropMap {
        &self.props
    }

    pub fn props_mut(&mut self) -> &mut PropMap {
        &mut self.props
    }

    pub fn set_velocities(&mut self, velocities: Option<Velocities>) -> Result<(), FrameError> {
        if velocities.as_ref().is_some_and(|values| {
            values.0.topology != *self.configuration.positions().topology_identity()
        }) {
            return Err(FrameError::TopologyIdentityMismatch);
        }
        self.velocities = velocities;
        Ok(())
    }

    pub fn set_forces(&mut self, forces: Option<Forces>) -> Result<(), FrameError> {
        if forces.as_ref().is_some_and(|values| {
            values.0.topology != *self.configuration.positions().topology_identity()
        }) {
            return Err(FrameError::TopologyIdentityMismatch);
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

    pub fn set_observation(
        &mut self,
        observation: Option<StructureObservation>,
    ) -> Result<(), FrameError> {
        if observation.as_ref().is_some_and(|observation| {
            observation.topology_identity() != self.configuration.positions().topology_identity()
        }) {
            return Err(FrameError::TopologyIdentityMismatch);
        }
        self.observation = observation;
        Ok(())
    }

    pub fn validate(&self, topology: &Topology) -> Result<(), FrameError> {
        if !self.configuration.positions().is_compatible(topology)
            || self
                .velocities
                .as_ref()
                .is_some_and(|values| !values.is_compatible(topology))
            || self
                .forces
                .as_ref()
                .is_some_and(|values| !values.is_compatible(topology))
            || self
                .observation
                .as_ref()
                .is_some_and(|observation| !observation.is_compatible(topology))
        {
            return Err(FrameError::TopologyIdentityMismatch);
        }
        Ok(())
    }

    pub fn view<'a>(
        &'a self,
        topology: &'a Topology,
    ) -> Result<TrajectoryFrameView<'a>, FrameError> {
        self.validate(topology)?;
        Ok(TrajectoryFrameView {
            topology,
            configuration: self.configuration.view(),
            velocities: self.velocities.as_ref().map(Velocities::values),
            forces: self.forces.as_ref().map(Forces::values),
            time: self.time,
            step: self.step,
            observation: self.observation.as_ref(),
            props: &self.props,
        })
    }
}

/// Borrowed trajectory frame state.
#[derive(Debug, Clone, Copy)]
pub struct TrajectoryFrameView<'a> {
    topology: &'a Topology,
    configuration: ConfigurationView<'a>,
    velocities: Option<Quantity<&'a [Vector3]>>,
    forces: Option<Quantity<&'a [Vector3]>>,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
    observation: Option<&'a StructureObservation>,
    props: &'a PropMap,
}

impl<'a> TrajectoryFrameView<'a> {
    pub const fn topology(self) -> &'a Topology {
        self.topology
    }

    pub const fn configuration(self) -> ConfigurationView<'a> {
        self.configuration
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

    pub const fn observation(self) -> Option<&'a StructureObservation> {
        self.observation
    }

    pub const fn props(self) -> &'a PropMap {
        self.props
    }

    pub fn model_view(self) -> ModelView<'a> {
        ModelView::new(self.topology, self.configuration)
            .expect("trajectory frame view has validated topology")
    }
}

/// Reusable caller-owned frame storage bound to one exact topology.
#[derive(Debug, Clone)]
pub struct FrameBuffer {
    topology: Topology,
    configuration: Configuration,
    velocities: Velocities,
    has_velocities: bool,
    forces: Forces,
    has_forces: bool,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
    observation: Option<StructureObservation>,
    props: PropMap,
}

impl FrameBuffer {
    pub fn new(topology: Topology) -> Self {
        Self {
            configuration: Configuration::new(Positions::zeros(&topology)),
            velocities: Velocities::zeros(&topology),
            has_velocities: false,
            forces: Forces::zeros(&topology),
            has_forces: false,
            time: None,
            step: None,
            observation: None,
            props: PropMap::new(),
            topology,
        }
    }

    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    pub fn configuration(&self) -> &Configuration {
        &self.configuration
    }

    pub fn set_positions<T>(&mut self, positions: Quantity<T>) -> Result<(), FrameError>
    where
        T: AsRef<[Point3]>,
    {
        self.configuration
            .positions_mut()
            .set_all(&self.topology, positions)?;
        Ok(())
    }

    pub fn set_cell(&mut self, cell: Option<crate::geometry::PeriodicCell>) {
        self.configuration.set_cell(cell);
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

    pub fn set_observation(
        &mut self,
        observation: Option<StructureObservation>,
    ) -> Result<(), FrameError> {
        if observation
            .as_ref()
            .is_some_and(|observation| !observation.is_compatible(&self.topology))
        {
            return Err(FrameError::TopologyIdentityMismatch);
        }
        self.observation = observation;
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
        self.configuration.set_cell(None);
        self.has_velocities = false;
        self.has_forces = false;
        self.time = None;
        self.step = None;
        self.observation = None;
        self.props.clear();
    }

    pub fn model_view(&self) -> ModelView<'_> {
        ModelView::new(&self.topology, self.configuration.view())
            .expect("frame buffer configuration is bound to its topology")
    }

    pub fn frame_view(&self) -> TrajectoryFrameView<'_> {
        TrajectoryFrameView {
            topology: &self.topology,
            configuration: self.configuration.view(),
            velocities: self.has_velocities.then(|| self.velocities.values()),
            forces: self.has_forces.then(|| self.forces.values()),
            time: self.time,
            step: self.step,
            observation: self.observation.as_ref(),
            props: &self.props,
        }
    }

    pub fn copy_from(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), FrameError> {
        if !self.topology.same_identity(frame.topology) {
            return Err(FrameError::TopologyIdentityMismatch);
        }
        self.set_positions(frame.configuration.positions().values())?;
        self.set_cell(frame.configuration.cell().copied());
        self.set_velocities(frame.velocities)?;
        self.set_forces(frame.forces)?;
        self.set_time(frame.time)?;
        self.step = frame.step;
        self.observation = frame.observation.cloned();
        self.props = frame.props.clone();
        Ok(())
    }
}

/// Deliberately loaded finite in-memory trajectory.
#[derive(Debug, Clone)]
pub struct Trajectory {
    topology: Topology,
    frames: Vec<TrajectoryFrame>,
}

impl Trajectory {
    pub fn new(topology: Topology) -> Self {
        Self {
            topology,
            frames: Vec::new(),
        }
    }

    pub fn from_frames(
        topology: Topology,
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
}

pub trait TrajectoryReader {
    fn topology(&self) -> &Topology;

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

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        let Some(frame) = self.trajectory.frames.get(self.cursor) else {
            return Ok(false);
        };
        destination.copy_from(frame.view(self.trajectory.topology())?)?;
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
        destination.copy_from(frame.view(self.trajectory.topology())?)?;
        self.cursor = index.saturating_add(1);
        Ok(())
    }
}

/// In-memory writer that preserves every core frame field.
pub struct MemoryTrajectoryWriter {
    trajectory: Trajectory,
}

impl MemoryTrajectoryWriter {
    pub fn new(topology: Topology) -> Self {
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

    fn write_frame(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError> {
        if !self.trajectory.topology.same_identity(frame.topology) {
            return Err(TrajectoryError::TopologyIdentityMismatch);
        }
        let positions = Positions::new(
            &self.trajectory.topology,
            frame.configuration.positions().values(),
        )?;
        let configuration = match frame.configuration.cell().copied() {
            Some(cell) => Configuration::with_cell(positions, cell),
            None => Configuration::new(positions),
        };
        let mut owned = TrajectoryFrame::new(configuration);
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
        owned.observation = frame.observation.cloned();
        owned.props = frame.props.clone();
        self.trajectory.push(owned)
    }
}

/// Proof that a coordinate-only source order exactly matches one topology.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomOrderAssertion {
    topology: TopologyIdentity,
}

impl AtomOrderAssertion {
    pub fn new(
        topology: &Topology,
        atom_order: &[InstanceAtomId],
    ) -> Result<Self, TrajectoryError> {
        if topology.atom_ids() != atom_order {
            return Err(TrajectoryError::AtomOrderMismatch);
        }
        Ok(Self {
            topology: topology.identity(),
        })
    }
}

/// Reference reader for a topology-free coordinate source.
pub struct CoordinateFrameReader {
    topology: Topology,
    frames: Vec<Vec<Point3>>,
    cursor: usize,
}

impl CoordinateFrameReader {
    pub fn new(
        topology: Topology,
        assertion: AtomOrderAssertion,
        frames: impl IntoIterator<Item = Quantity<Vec<Point3>>>,
    ) -> Result<Self, TrajectoryError> {
        if assertion.topology != topology.identity() {
            return Err(TrajectoryError::TopologyIdentityMismatch);
        }
        let frames = frames
            .into_iter()
            .map(|frame| {
                Positions::new(&topology, frame)
                    .map(|positions| positions.values_raw().to_vec())
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

    fn read_next(&mut self, destination: &mut FrameBuffer) -> Result<bool, TrajectoryError> {
        if !self.topology.same_identity(destination.topology()) {
            return Err(TrajectoryError::TopologyIdentityMismatch);
        }
        let Some(frame) = self.frames.get(self.cursor) else {
            return Ok(false);
        };
        destination.set_positions(Quantity::new(
            frame.as_slice(),
            crate::units::MODEL_LENGTH_UNIT,
        ))?;
        destination.reset_dynamic_state();
        self.cursor += 1;
        Ok(true)
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum FrameError {
    TopologyIdentityMismatch,
    AtomCountMismatch { expected: usize, actual: usize },
    NonFiniteVector { atom: InstanceAtomId },
    NonFiniteTime,
    Position(PositionError),
    Observation(ObservationError),
    Unit(UnitError),
}

impl fmt::Display for FrameError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyIdentityMismatch => {
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
            Self::Observation(error) => {
                write!(formatter, "invalid frame observation: {error}")
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

impl From<ObservationError> for FrameError {
    fn from(error: ObservationError) -> Self {
        Self::Observation(error)
    }
}

impl From<UnitError> for FrameError {
    fn from(error: UnitError) -> Self {
        Self::Unit(error)
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TrajectoryError {
    TopologyIdentityMismatch,
    AtomOrderMismatch,
    FrameIndexOutOfRange(u64),
    UnsupportedRandomAccess,
    MissingRequiredTopology,
    MissingTime { frame: usize },
    NonMonotonicTime { frame: usize },
    UnsupportedField(&'static str),
    Frame(Box<FrameError>),
    Position(PositionError),
}

impl fmt::Display for TrajectoryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TopologyIdentityMismatch => {
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
        }
    }
}

impl std::error::Error for TrajectoryError {}

impl From<FrameError> for TrajectoryError {
    fn from(error: FrameError) -> Self {
        match error {
            FrameError::TopologyIdentityMismatch => Self::TopologyIdentityMismatch,
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
        Self::TopologyIdentityMismatch
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::{Atom, Element, Molecule, PropValue};
    use crate::small::SmallMolecule;
    use crate::topology::{MoleculeInstanceMetadata, TopologyBuilder};
    use crate::units::{ANGSTROM, PICOSECOND};

    fn one_atom_topology() -> Topology {
        let mut graph = Molecule::new();
        graph.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
        let molecule = SmallMolecule::from_graph(graph);
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_small_molecule_definition(&molecule).unwrap();
        builder
            .add_instance(definition, MoleculeInstanceMetadata::default())
            .unwrap();
        builder.build().unwrap()
    }

    fn configuration(topology: &Topology, x: f64) -> Configuration {
        Configuration::new(
            Positions::new(
                topology,
                Quantity::new(vec![Point3::new(x, 0.0, 0.0)], ANGSTROM),
            )
            .unwrap(),
        )
    }

    fn frame(topology: &Topology, x: f64, time: f64) -> TrajectoryFrame {
        let mut frame = TrajectoryFrame::new(configuration(topology, x));
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
        let mut frame = TrajectoryFrame::new(configuration(&topology, 0.0));
        assert_eq!(
            frame.set_velocities(Some(Velocities::zeros(&independent))),
            Err(FrameError::TopologyIdentityMismatch)
        );
        assert_eq!(
            frame.set_observation(Some(StructureObservation::empty(&independent))),
            Err(FrameError::TopologyIdentityMismatch)
        );
        assert_eq!(
            frame.set_time(Some(Quantity::new(f64::INFINITY, PICOSECOND))),
            Err(FrameError::NonFiniteTime)
        );
    }

    #[test]
    fn reusable_buffer_streams_and_seeks_without_reallocating_positions() {
        let topology = one_atom_topology();
        let mut first = frame(&topology, 0.0, 0.0);
        let observation = StructureObservation::empty(&topology);
        first.set_observation(Some(observation)).unwrap();
        let cell = crate::geometry::PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(10.0, 10.0, 10.0), ANGSTROM),
            [true; 3],
        )
        .unwrap();
        let second_positions = Positions::new(
            &topology,
            Quantity::new(vec![Point3::new(2.0, 0.0, 0.0)], ANGSTROM),
        )
        .unwrap();
        let mut second = TrajectoryFrame::new(Configuration::with_cell(second_positions, cell));
        second
            .set_time(Some(Quantity::new(1.0, PICOSECOND)))
            .unwrap();
        second.set_step(Some(1));
        let trajectory = Trajectory::from_frames(topology.clone(), [first, second]).unwrap();
        assert_eq!(trajectory.len(), 2);
        trajectory.validate_monotonic_time(true).unwrap();

        let mut reader = MemoryTrajectoryReader::new(&trajectory);
        let mut buffer = FrameBuffer::new(topology.clone());
        let pointer = buffer.configuration().positions().values().value().as_ptr();
        assert!(reader.read_next(&mut buffer).unwrap());
        assert_eq!(
            buffer.configuration().positions().values().value().as_ptr(),
            pointer
        );
        assert!(buffer.frame_view().velocities().is_some());
        assert!(buffer.frame_view().observation().is_some());
        assert!(reader.read_next(&mut buffer).unwrap());
        assert_eq!(
            buffer.configuration().positions().values().value().as_ptr(),
            pointer
        );
        assert_eq!(buffer.model_view().positions().value()[0].x, 2.0);
        assert_eq!(buffer.frame_view().configuration().cell(), Some(&cell));
        assert!(!reader.read_next(&mut buffer).unwrap());

        reader.read_frame(0, &mut buffer).unwrap();
        assert_eq!(buffer.model_view().positions().value()[0].x, 0.0);
        assert_eq!(reader.frame_count(), Some(2));
    }

    #[test]
    fn memory_writer_round_trips_frames_and_rejects_other_topologies() {
        let topology = one_atom_topology();
        let trajectory =
            Trajectory::from_frames(topology.clone(), [frame(&topology, 3.0, 2.0)]).unwrap();
        let mut writer = MemoryTrajectoryWriter::new(topology.clone());
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
                .value()[0]
                .x,
            3.0
        );

        let independent = one_atom_topology();
        let foreign =
            Trajectory::from_frames(independent.clone(), [frame(&independent, 4.0, 3.0)]).unwrap();
        let mut writer = MemoryTrajectoryWriter::new(topology);
        assert_eq!(
            writer.write_frame(foreign.frames().next().unwrap()),
            Err(TrajectoryError::TopologyIdentityMismatch)
        );
    }

    #[test]
    fn coordinate_only_reader_requires_atom_order_and_matching_buffer_identity() {
        let topology = one_atom_topology();
        assert_eq!(
            AtomOrderAssertion::new(&topology, &[]),
            Err(TrajectoryError::AtomOrderMismatch)
        );
        let assertion = AtomOrderAssertion::new(&topology, topology.atom_ids()).unwrap();
        let mut reader = CoordinateFrameReader::new(
            topology.clone(),
            assertion,
            [Quantity::new(vec![Point3::new(5.0, 0.0, 0.0)], ANGSTROM)],
        )
        .unwrap();
        let independent = one_atom_topology();
        let mut wrong_buffer = FrameBuffer::new(independent);
        assert_eq!(
            reader.read_next(&mut wrong_buffer),
            Err(TrajectoryError::TopologyIdentityMismatch)
        );

        let mut buffer = FrameBuffer::new(topology);
        assert!(reader.read_next(&mut buffer).unwrap());
        assert_eq!(buffer.model_view().positions().value()[0].x, 5.0);
        assert!(!reader.read_next(&mut buffer).unwrap());
    }

    #[test]
    fn coordinate_only_reader_clears_all_dynamic_state_without_reallocating_positions() {
        let topology = one_atom_topology();
        let assertion = AtomOrderAssertion::new(&topology, topology.atom_ids()).unwrap();
        let mut reader = CoordinateFrameReader::new(
            topology.clone(),
            assertion,
            [Quantity::new(vec![Point3::new(7.0, 0.0, 0.0)], ANGSTROM)],
        )
        .unwrap();
        let mut buffer = FrameBuffer::new(topology.clone());
        let cell = crate::geometry::PeriodicCell::orthorhombic(
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
        buffer
            .set_observation(Some(StructureObservation::empty(&topology)))
            .unwrap();
        buffer
            .props_mut()
            .insert("source".to_owned(), PropValue::String("stale".to_owned()));
        let positions_pointer = buffer.configuration().positions().values().value().as_ptr();

        assert!(reader.read_next(&mut buffer).unwrap());

        let frame = buffer.frame_view();
        assert_eq!(frame.configuration().cell(), None);
        assert_eq!(frame.velocities(), None);
        assert_eq!(frame.forces(), None);
        assert_eq!(frame.time(), None);
        assert_eq!(frame.step(), None);
        assert_eq!(frame.observation(), None);
        assert!(frame.props().is_empty());
        assert_eq!(
            buffer.configuration().positions().values().value().as_ptr(),
            positions_pointer
        );
        assert_eq!(buffer.model_view().positions().value()[0].x, 7.0);
    }
}
