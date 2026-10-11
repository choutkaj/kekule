use std::sync::Arc;

use crate::geometry::{PeriodicCell, RigidTransform};
use crate::properties::RealizationProperties;
use crate::topology::Topology;
use crate::units::{Quantity, CANONICAL_TIME_UNIT};

use super::realizations::sealed::Payload;
use super::realizations::{realization_collection, RealizationStore};
use super::scalars::ABSENT;
use super::state::{
    realization_mutators, realization_readers, transform_cell, transform_positions,
};
use super::{
    ConformationError, Ensemble, EnsembleMember, Forces, ModelView, Positions, Realization,
    RealizationError, RealizationMut, RealizationView, Velocities,
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
    store: RealizationStore<TrajectoryFrame>,
}

realization_collection!(
    Trajectory,
    TrajectoryFrame,
    TrajectoryFrameView<'_>,
    TrajectoryFrameMut<'_>,
    "frame"
);

/// One topology-bound [`TrajectoryFrame`] borrowed from a [`Trajectory`].
pub type TrajectoryFrameView<'a> = RealizationView<'a, TrajectoryFrame>;

/// Dimension-preserving editor for one [`TrajectoryFrame`] of a [`Trajectory`].
pub type TrajectoryFrameMut<'a> = RealizationMut<'a, TrajectoryFrame>;

/// One trajectory frame: positions in dense topology atom order, an optional
/// periodic cell, realization properties, and optional velocities, forces,
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
    properties: RealizationProperties,
    velocities: Option<Velocities>,
    forces: Option<Forces>,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
}

impl TrajectoryFrame {
    /// Creates a non-periodic frame without properties or dynamics.
    pub fn new(positions: Positions) -> Self {
        let properties = RealizationProperties::new(positions.len());
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
    pub fn set_velocities(
        &mut self,
        velocities: Option<Velocities>,
    ) -> Result<(), ConformationError> {
        if let Some(values) = &velocities {
            self.check_len(values.len())?;
        }
        self.velocities = velocities;
        Ok(())
    }

    /// Sets or clears forces; they must have one vector per atom.
    pub fn set_forces(&mut self, forces: Option<Forces>) -> Result<(), ConformationError> {
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
    pub fn set_time(&mut self, time: Option<Quantity<f64>>) -> Result<(), ConformationError> {
        self.time = time
            .map(|time| {
                let time = time.into_unit(CANONICAL_TIME_UNIT)?;
                if !time.value().is_finite() {
                    return Err(ConformationError::NonFiniteTime);
                }
                Ok(time)
            })
            .transpose()?;
        Ok(())
    }

    pub fn set_step(&mut self, step: Option<u64>) {
        self.step = step;
    }

    fn check_len(&self, actual: usize) -> Result<(), ConformationError> {
        if actual != self.positions.len() {
            return Err(ConformationError::AtomCountMismatch {
                expected: self.positions.len(),
                actual,
            });
        }
        Ok(())
    }

    fn storage_ref(&self) -> &Self {
        self
    }

    fn storage(&mut self) -> &mut Self {
        self
    }

    realization_readers!();
    realization_mutators!();
}

impl Payload for TrajectoryFrame {
    fn row_count(&self) -> usize {
        self.positions.len()
    }

    fn realization_properties(&self) -> &RealizationProperties {
        &self.properties
    }

    fn positions_storage_mut(&mut self) -> &mut Positions {
        &mut self.positions
    }

    fn properties_storage_mut(&mut self) -> &mut RealizationProperties {
        &mut self.properties
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

    fn project(&self, atoms: &[usize], bonds: &[usize]) -> Result<Self, ConformationError> {
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
    fn transformed(&self, transform: RigidTransform) -> Result<Self, ConformationError> {
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
}

impl Realization for TrajectoryFrame {}

impl TrajectoryFrame {
    /// Applies a rigid transform transactionally: positions move; cells,
    /// velocities, and forces rotate.
    pub fn apply_rigid_transform(
        &mut self,
        transform: &RigidTransform,
    ) -> Result<(), ConformationError> {
        *self = self.transformed(*transform)?;
        Ok(())
    }
}

impl RealizationMut<'_, TrajectoryFrame> {
    fn storage(&mut self) -> &mut TrajectoryFrame {
        self.item
    }

    realization_mutators!();

    pub fn set_velocities(
        &mut self,
        velocities: Option<Velocities>,
    ) -> Result<(), ConformationError> {
        self.item.set_velocities(velocities)
    }

    pub fn set_forces(&mut self, forces: Option<Forces>) -> Result<(), ConformationError> {
        self.item.set_forces(forces)
    }

    pub fn set_time(&mut self, time: Option<Quantity<f64>>) -> Result<(), ConformationError> {
        self.item.set_time(time)
    }

    pub fn set_step(&mut self, step: Option<u64>) {
        self.item.set_step(step);
    }

    /// Removes velocities, returning their storage for reuse.
    pub fn take_velocities(&mut self) -> Option<Velocities> {
        self.item.take_velocities()
    }

    /// Removes forces, returning their storage for reuse.
    pub fn take_forces(&mut self) -> Option<Forces> {
        self.item.take_forces()
    }
}

impl Trajectory {
    /// Checks that every present time is non-decreasing; with `require_all`,
    /// a frame without time is an error.
    pub fn validate_monotonic_time(&self, require_all: bool) -> Result<(), RealizationError> {
        let mut previous = None;
        for (index, frame) in self.iter().enumerate() {
            let Some(time) = frame.time() else {
                if require_all {
                    return Err(RealizationError::MissingTime { frame: index });
                }
                continue;
            };
            let value = time.into_value();
            if previous.is_some_and(|previous| value < previous) {
                return Err(RealizationError::NonMonotonicTime { frame: index });
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
        Ensemble::from_store(self.store.map_items(|frame| {
            EnsembleMember::from_frame_state(frame.positions, frame.cell, frame.properties)
        }))
    }

    pub(crate) fn store(&self) -> &RealizationStore<TrajectoryFrame> {
        &self.store
    }

    pub(crate) fn store_mut(&mut self) -> &mut RealizationStore<TrajectoryFrame> {
        &mut self.store
    }
}
