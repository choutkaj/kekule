use crate::geometry::RigidTransform;
use crate::units::{Quantity, CANONICAL_TIME_UNIT};

use super::realizations::sealed::Payload;
use super::{
    Conformation, ConformationError, ConformationMut, Ensemble, EnsembleMember, Forces,
    Realization, RealizationError, RealizationMut, RealizationView, Realizations, Velocities,
};

/// An ordered temporal sequence of realizations of one shared topology.
///
/// In-memory storage and transforms live here; streaming readers, writers,
/// codecs, and periodic preprocessing live in the `kekule-traj` companion
/// crate.
pub type Trajectory = Realizations<TrajectoryFrame>;

/// One topology-bound [`TrajectoryFrame`] borrowed from a [`Trajectory`].
pub type TrajectoryFrameView<'a> = RealizationView<'a, TrajectoryFrame>;

/// Dimension-preserving editor for one [`TrajectoryFrame`] of a [`Trajectory`].
pub type TrajectoryFrameMut<'a> = RealizationMut<'a, TrajectoryFrame>;

/// One trajectory frame: a [`Conformation`] plus optional velocities, forces,
/// time, and step.
///
/// The frame carries no topology; it is bound when inserted into a
/// [`Trajectory`]. Conformation state reads through `Deref`.
#[derive(Debug, Clone, PartialEq)]
pub struct TrajectoryFrame {
    conformation: Conformation,
    velocities: Option<Velocities>,
    forces: Option<Forces>,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
}

impl From<Conformation> for TrajectoryFrame {
    fn from(conformation: Conformation) -> Self {
        Self {
            conformation,
            velocities: None,
            forces: None,
            time: None,
            step: None,
        }
    }
}

impl std::ops::Deref for TrajectoryFrame {
    type Target = Conformation;

    fn deref(&self) -> &Self::Target {
        &self.conformation
    }
}

impl TrajectoryFrame {
    /// Creates a frame without dynamics; accepts [`Conformation`] or
    /// [`super::Positions`].
    pub fn new(conformation: impl Into<Conformation>) -> Self {
        Self::from(conformation.into())
    }

    pub fn conformation(&self) -> &Conformation {
        &self.conformation
    }

    pub fn conformation_mut(&mut self) -> ConformationMut<'_> {
        ConformationMut::new(&mut self.conformation)
    }

    pub fn into_conformation(self) -> Conformation {
        self.conformation
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
        if actual != self.conformation.atom_count() {
            return Err(ConformationError::AtomCountMismatch {
                expected: self.conformation.atom_count(),
                actual,
            });
        }
        Ok(())
    }
}

impl Payload for TrajectoryFrame {
    fn conformation_storage(&mut self) -> &mut Conformation {
        &mut self.conformation
    }

    fn project(&self, atoms: &[usize], bonds: &[usize]) -> Result<Self, ConformationError> {
        Ok(Self {
            conformation: self.conformation.project(atoms, bonds)?,
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
            conformation: self.conformation.transformed(transform)?,
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

    /// Reinterprets the frames as an unweighted ensemble, dropping velocities,
    /// forces, time, and step. Conformations move without copying.
    pub fn into_ensemble(self) -> Ensemble {
        self.map_items(|frame| EnsembleMember::from(frame.into_conformation()))
    }
}
