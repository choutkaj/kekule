use std::fmt;
use std::sync::Arc;

use crate::geometry::{PeriodicCell, RigidTransform};
use crate::properties::{OwnerProperties, StructureProperties, TopologyPropertiesMut};
use crate::topology::TopologySubsetError;
use crate::topology::{AtomSelection, Topology, TopologyPerceptionError};
use crate::units::{Quantity, CANONICAL_TIME_UNIT};

use super::scalars::ABSENT;
use super::state::{
    is_bound, structure_mutators, structure_readers, subset_rows, transform_cell,
    transform_positions, validate_binding,
};
use super::{
    AsModelView, Ensemble, EnsembleMember, Forces, Model, ModelView, Positions, StructureError,
    Velocities,
};

/// A time-ordered sequence of frames of one shared topology, as produced by
/// dynamics.
///
/// Frame order is temporal: frame `i + 1` follows frame `i`. Frames may carry
/// time, step, velocities, and forces, and selection keeps each frame's stored
/// time and step. Use [`super::Ensemble`] for a weighted, unordered sample;
/// [`Self::into_ensemble`] is the explicit, lossy projection from one to the
/// other.
///
/// Frames are validated against the shared topology on insertion.
/// Superposition and RMSD are inherent methods; see [`crate::alignment`].
/// In-memory storage lives here; file codecs, streaming readers and writers,
/// periodic preprocessing, and trajectory reductions live in the
/// `kekule-traj` companion crate.
#[derive(Debug, Clone)]
pub struct Trajectory {
    topology: Arc<Topology>,
    properties: OwnerProperties,
    frames: Vec<TrajectoryFrame>,
}

impl Trajectory {
    /// An empty trajectory of one shared topology.
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        Self {
            topology: topology.into(),
            properties: OwnerProperties::new(),
            frames: Vec::new(),
        }
    }

    /// Validates and collects each frame in order, binding its bond property
    /// rows. Fails on the first invalid frame.
    pub fn from_items(
        topology: impl Into<Arc<Topology>>,
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

    /// Trajectory-level annotations.
    pub fn properties(&self) -> &OwnerProperties {
        &self.properties
    }

    pub fn properties_mut(&mut self) -> &mut OwnerProperties {
        &mut self.properties
    }

    pub fn len(&self) -> usize {
        self.frames.len()
    }

    pub fn is_empty(&self) -> bool {
        self.frames.is_empty()
    }

    /// One topology-bound frame by stable index.
    pub fn get(&self, index: usize) -> Option<TrajectoryFrameView<'_>> {
        self.frames.get(index).map(|frame| TrajectoryFrameView {
            topology: &self.topology,
            frame,
        })
    }

    /// A dimension-preserving editor for one frame. Whole-frame replacement
    /// goes through [`Self::replace`].
    pub fn get_mut(&mut self, index: usize) -> Option<TrajectoryFrameMut<'_>> {
        self.frames
            .get_mut(index)
            .map(|frame| TrajectoryFrameMut { frame })
    }

    /// Topology-bound frames in temporal order.
    pub fn iter(&self) -> TrajectoryIter<'_> {
        TrajectoryIter {
            topology: &self.topology,
            frames: self.frames.iter(),
        }
    }

    /// Validates and appends one frame, binding its bond property rows.
    pub fn push(&mut self, mut frame: TrajectoryFrame) -> Result<(), TrajectoryError> {
        self.bind(&mut frame)?;
        self.frames.push(frame);
        Ok(())
    }

    /// Validates and replaces one frame, returning the previous one. The
    /// index is checked first; failure leaves the trajectory unchanged.
    pub fn replace(
        &mut self,
        index: usize,
        mut frame: TrajectoryFrame,
    ) -> Result<TrajectoryFrame, TrajectoryError> {
        self.check_index(index)?;
        self.bind(&mut frame)?;
        Ok(std::mem::replace(&mut self.frames[index], frame))
    }

    /// Removes and returns one frame, shifting later ones down.
    pub fn remove(&mut self, index: usize) -> Result<TrajectoryFrame, TrajectoryError> {
        self.check_index(index)?;
        Ok(self.frames.remove(index))
    }

    /// Replaces the positions of every frame transactionally, in order.
    /// Other state is kept.
    pub fn replace_positions(
        &mut self,
        positions: impl IntoIterator<Item = Positions>,
    ) -> Result<(), TrajectoryError> {
        let positions = positions.into_iter().collect::<Vec<_>>();
        if positions.len() != self.frames.len() {
            return Err(TrajectoryError::FrameCountMismatch {
                expected: self.frames.len(),
                actual: positions.len(),
            });
        }
        for values in &positions {
            if values.len() != self.topology.atom_count() {
                return Err(StructureError::AtomCountMismatch {
                    expected: self.topology.atom_count(),
                    actual: values.len(),
                }
                .into());
            }
        }
        for (frame, positions) in self.frames.iter_mut().zip(positions) {
            frame.positions = positions;
        }
        Ok(())
    }

    /// Copies frames in the requested order, sharing the topology. Ranges and
    /// strides can be written as `(start..end).step_by(stride)`; duplicates
    /// and reordering are allowed, and each frame keeps its stored time and
    /// step. Trajectory properties are kept.
    pub fn select(
        &self,
        indices: impl IntoIterator<Item = usize>,
    ) -> Result<Self, TrajectoryError> {
        let frames = indices
            .into_iter()
            .map(|index| {
                self.check_index(index)?;
                Ok(self.frames[index].clone())
            })
            .collect::<Result<Vec<_>, TrajectoryError>>()?;
        Ok(Self {
            topology: self.shared_topology(),
            properties: self.properties.clone(),
            frames,
        })
    }

    /// Publishes one induced topology subset and projects every frame onto
    /// it. Trajectory and per-frame owner properties are dropped.
    pub fn subset(&self, selection: &AtomSelection) -> Result<Self, TrajectoryError> {
        let subset = self.topology.subset(selection)?;
        let (atoms, bonds) = subset_rows(&subset);
        let frames = self
            .frames
            .iter()
            .map(|frame| frame.project(&atoms, &bonds))
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            topology: subset.shared_topology(),
            properties: OwnerProperties::new(),
            frames,
        })
    }

    /// Installs default perception through one new shared topology snapshot.
    ///
    /// Delegates to [`Topology::perceived`] once and keeps every frame and the
    /// trajectory properties without copying. Failure leaves the trajectory
    /// unchanged. The new snapshot shares the original layout, so selections,
    /// buffers, and prepared potentials stay usable.
    pub fn perceive(&mut self) -> Result<(), TopologyPerceptionError> {
        self.topology = Arc::new(self.topology.perceived()?);
        Ok(())
    }

    /// Length-preserving mutable access to the shared topology's static
    /// annotations.
    ///
    /// A topology shared with other owners is first replaced by a snapshot of
    /// it ([`Arc::make_mut`]), so those owners are unaffected. Either way the
    /// layout is kept, so every frame, selection, buffer, and prepared
    /// potential remains usable. See [`Topology::properties_mut`].
    pub fn topology_properties_mut(&mut self) -> TopologyPropertiesMut<'_> {
        Arc::make_mut(&mut self.topology).properties_mut()
    }

    /// Consumes the trajectory, transferring topology, trajectory properties,
    /// and every frame without copying.
    pub fn into_parts(self) -> (Arc<Topology>, OwnerProperties, Vec<TrajectoryFrame>) {
        (self.topology, self.properties, self.frames)
    }

    pub fn into_items(self) -> Vec<TrajectoryFrame> {
        self.frames
    }

    /// Checks that every present time is non-decreasing; with `require_all`,
    /// a frame without time is an error.
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

    /// Reinterprets the frames as an equally weighted ensemble (every member
    /// has weight `1.0`), dropping temporal order, velocities, forces, time,
    /// and step. Positions, cells, and properties move without copying.
    ///
    /// Equal weights are right for frames sampled from equilibrium at a
    /// constant interval; reweight explicitly otherwise.
    pub fn into_ensemble(self) -> Ensemble {
        let members = self
            .frames
            .into_iter()
            .map(|frame| {
                EnsembleMember::from_frame_state(frame.positions, frame.cell, frame.properties)
            })
            .collect();
        Ensemble::from_bound_members(self.topology, self.properties, members)
    }

    /// Applies one rigid transform per frame transactionally.
    pub(crate) fn transform_frames(
        &mut self,
        transforms: &[RigidTransform],
    ) -> Result<(), StructureError> {
        debug_assert_eq!(transforms.len(), self.frames.len());
        let frames = self
            .frames
            .iter()
            .zip(transforms)
            .map(|(frame, transform)| frame.transformed(*transform))
            .collect::<Result<Vec<_>, _>>()?;
        self.frames = frames;
        Ok(())
    }

    fn check_index(&self, index: usize) -> Result<(), TrajectoryError> {
        if index >= self.frames.len() {
            return Err(TrajectoryError::IndexOutOfRange {
                index,
                len: self.frames.len(),
            });
        }
        Ok(())
    }

    fn bind(&self, frame: &mut TrajectoryFrame) -> Result<(), TrajectoryError> {
        validate_binding(frame.atom_count(), &frame.properties, &self.topology)?;
        frame.properties.bind_bond_rows(self.topology.bond_count());
        Ok(())
    }
}

impl<'a> IntoIterator for &'a Trajectory {
    type Item = TrajectoryFrameView<'a>;
    type IntoIter = TrajectoryIter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// One trajectory frame: positions in dense topology atom order, an optional
/// periodic cell, structure properties, and optional velocities, forces,
/// time, and step.
///
/// Frames carry no crystallographic occupancies or B-factors. The frame
/// carries no topology; it is bound when inserted into a [`Trajectory`].
///
/// ```compile_fail,E0599
/// use kekule::structure::{Positions, TrajectoryFrame};
/// let mut frame = TrajectoryFrame::new(Positions::zeros(1));
/// frame.set_occupancies(None);
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct TrajectoryFrame {
    positions: Positions,
    cell: Option<PeriodicCell>,
    properties: StructureProperties,
    velocities: Option<Velocities>,
    forces: Option<Forces>,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
}

impl TrajectoryFrame {
    /// Creates a non-periodic frame without properties or dynamics.
    pub fn new(positions: Positions) -> Self {
        let properties = StructureProperties::new(positions.len());
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

    /// Returns this frame with a periodic cell.
    #[must_use]
    pub fn with_cell(mut self, cell: PeriodicCell) -> Self {
        self.cell = Some(cell);
        self
    }

    pub fn atom_count(&self) -> usize {
        self.positions.len()
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

    /// Sets or clears velocities; they must have one vector per atom.
    pub fn set_velocities(&mut self, velocities: Option<Velocities>) -> Result<(), StructureError> {
        if let Some(values) = &velocities {
            self.check_len(values.len())?;
        }
        self.velocities = velocities;
        Ok(())
    }

    /// Sets or clears forces; they must have one vector per atom.
    pub fn set_forces(&mut self, forces: Option<Forces>) -> Result<(), StructureError> {
        if let Some(values) = &forces {
            self.check_len(values.len())?;
        }
        self.forces = forces;
        Ok(())
    }

    /// Removes velocities, returning their storage for reuse.
    pub fn take_velocities(&mut self) -> Option<Velocities> {
        self.velocities.take()
    }

    /// Removes forces, returning their storage for reuse.
    pub fn take_forces(&mut self) -> Option<Forces> {
        self.forces.take()
    }

    /// Sets or clears a finite time, converted to the canonical time unit.
    pub fn set_time(&mut self, time: Option<Quantity<f64>>) -> Result<(), StructureError> {
        self.time = time
            .map(|time| {
                let time = time.into_unit(CANONICAL_TIME_UNIT)?;
                if !time.value().is_finite() {
                    return Err(StructureError::NonFiniteTime);
                }
                Ok(time)
            })
            .transpose()?;
        Ok(())
    }

    pub fn set_step(&mut self, step: Option<u64>) {
        self.step = step;
    }

    /// Applies a rigid transform transactionally: positions move; cells,
    /// velocities, and forces rotate.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), StructureError> {
        *self = self.transformed(*transform)?;
        Ok(())
    }

    fn check_len(&self, actual: usize) -> Result<(), StructureError> {
        if actual != self.positions.len() {
            return Err(StructureError::AtomCountMismatch {
                expected: self.positions.len(),
                actual,
            });
        }
        Ok(())
    }

    fn model_view<'a>(&'a self, topology: &'a Arc<Topology>) -> ModelView<'a> {
        ModelView::bound(
            topology,
            &self.positions,
            self.cell.as_ref(),
            &ABSENT,
            &ABSENT,
            &self.properties,
        )
    }

    /// Projects every dense row through a subset correspondence.
    fn project(&self, atoms: &[usize], bonds: &[usize]) -> Result<Self, StructureError> {
        Ok(Self {
            positions: self.positions.select_indices(atoms)?,
            cell: self.cell,
            properties: self.properties.project(atoms, bonds)?,
            velocities: self.velocities.as_ref().map(|values| values.select(atoms)),
            forces: self.forces.as_ref().map(|values| values.select(atoms)),
            time: self.time,
            step: self.step,
        })
    }

    /// Positions move; cells, velocities, and forces rotate.
    fn transformed(&self, transform: RigidTransform) -> Result<Self, StructureError> {
        let rotate = |vector| transform.transform_vector(vector);
        Ok(Self {
            cell: transform_cell(self.cell, transform)?,
            positions: transform_positions(&self.positions, transform)?,
            properties: self.properties.clone(),
            velocities: self
                .velocities
                .as_ref()
                .map(|values| values.rotated(rotate))
                .transpose()?,
            forces: self
                .forces
                .as_ref()
                .map(|values| values.rotated(rotate))
                .transpose()?,
            time: self.time,
            step: self.step,
        })
    }

    fn storage_ref(&self) -> &Self {
        self
    }

    fn storage(&mut self) -> &mut Self {
        self
    }

    structure_readers!();
    structure_mutators!();
}

/// One topology-bound [`TrajectoryFrame`] borrowed from a [`Trajectory`].
///
/// Reads frame state through `Deref` and converts to a [`ModelView`] through
/// [`AsModelView`] without copying.
#[derive(Debug, Clone, Copy)]
pub struct TrajectoryFrameView<'a> {
    topology: &'a Arc<Topology>,
    frame: &'a TrajectoryFrame,
}

impl std::ops::Deref for TrajectoryFrameView<'_> {
    type Target = TrajectoryFrame;

    fn deref(&self) -> &Self::Target {
        self.frame
    }
}

impl<'a> TrajectoryFrameView<'a> {
    /// Binds a frame that already has every row `topology` requires.
    pub fn new(
        topology: &'a Arc<Topology>,
        frame: &'a TrajectoryFrame,
    ) -> Result<Self, StructureError> {
        validate_binding(frame.atom_count(), &frame.properties, topology)?;
        if !is_bound(frame.atom_count(), &frame.properties, topology) {
            return Err(StructureError::BondCountMismatch {
                expected: topology.bond_count(),
                actual: frame.properties.bonds().len(),
            });
        }
        Ok(Self { topology, frame })
    }

    pub fn topology(self) -> &'a Topology {
        self.topology
    }

    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.topology)
    }

    /// The borrowed frame with the view's full lifetime.
    pub fn payload(self) -> &'a TrajectoryFrame {
        self.frame
    }

    /// Borrows the frame as the common kernel input, with the view's full
    /// lifetime.
    pub fn as_model_view(self) -> ModelView<'a> {
        self.frame.model_view(self.topology)
    }

    /// Materializes the frame as an owned model, dropping velocities, forces,
    /// time, and step.
    pub fn to_model(self) -> Model {
        self.as_model_view().to_model()
    }
}

impl AsModelView for TrajectoryFrameView<'_> {
    fn as_model_view(&self) -> ModelView<'_> {
        TrajectoryFrameView::as_model_view(*self)
    }
}

/// Dimension-preserving mutable access to one frame of a [`Trajectory`].
///
/// Reads use `Deref`. There is deliberately no `DerefMut`; replace a whole
/// frame through [`Trajectory::replace`].
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

impl TrajectoryFrameMut<'_> {
    pub fn set_velocities(&mut self, velocities: Option<Velocities>) -> Result<(), StructureError> {
        self.frame.set_velocities(velocities)
    }

    pub fn set_forces(&mut self, forces: Option<Forces>) -> Result<(), StructureError> {
        self.frame.set_forces(forces)
    }

    pub fn set_time(&mut self, time: Option<Quantity<f64>>) -> Result<(), StructureError> {
        self.frame.set_time(time)
    }

    pub fn set_step(&mut self, step: Option<u64>) {
        self.frame.set_step(step);
    }

    /// Removes velocities, returning their storage for reuse.
    pub fn take_velocities(&mut self) -> Option<Velocities> {
        self.frame.take_velocities()
    }

    /// Removes forces, returning their storage for reuse.
    pub fn take_forces(&mut self) -> Option<Forces> {
        self.frame.take_forces()
    }

    /// Applies a rigid transform transactionally: positions move; cells,
    /// velocities, and forces rotate.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), StructureError> {
        self.frame.apply_rigid_transform(transform)
    }

    fn storage(&mut self) -> &mut TrajectoryFrame {
        self.frame
    }

    structure_mutators!();
}

/// Iterator over the topology-bound frames of a [`Trajectory`], in temporal
/// order.
#[derive(Debug, Clone)]
pub struct TrajectoryIter<'a> {
    topology: &'a Arc<Topology>,
    frames: std::slice::Iter<'a, TrajectoryFrame>,
}

impl<'a> Iterator for TrajectoryIter<'a> {
    type Item = TrajectoryFrameView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let topology = self.topology;
        self.frames
            .next()
            .map(|frame| TrajectoryFrameView { topology, frame })
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.frames.size_hint()
    }
}

impl ExactSizeIterator for TrajectoryIter<'_> {}

/// Failure of a [`Trajectory`] operation.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum TrajectoryError {
    IndexOutOfRange {
        index: usize,
        len: usize,
    },
    /// A trajectory-wide replacement has the wrong number of frames.
    FrameCountMismatch {
        expected: usize,
        actual: usize,
    },
    MissingTime {
        frame: usize,
    },
    NonMonotonicTime {
        frame: usize,
    },
    Structure(StructureError),
    Subset(TopologySubsetError),
}

impl fmt::Display for TrajectoryError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IndexOutOfRange { index, len } => {
                write!(formatter, "index {index} is out of range for {len} frames")
            }
            Self::FrameCountMismatch { expected, actual } => write!(
                formatter,
                "trajectory requires {expected} replacement frames, but received {actual}"
            ),
            Self::MissingTime { frame } => write!(formatter, "frame {frame} has no time"),
            Self::NonMonotonicTime { frame } => {
                write!(formatter, "time decreases at frame {frame}")
            }
            Self::Structure(_) => formatter.write_str("invalid trajectory frame"),
            Self::Subset(_) => formatter.write_str("cannot subset trajectory"),
        }
    }
}

impl std::error::Error for TrajectoryError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Structure(error) => Some(error),
            Self::Subset(error) => Some(error),
            _ => None,
        }
    }
}

impl From<StructureError> for TrajectoryError {
    fn from(error: StructureError) -> Self {
        Self::Structure(error)
    }
}

impl From<TopologySubsetError> for TrajectoryError {
    fn from(error: TopologySubsetError) -> Self {
        Self::Subset(error)
    }
}
