use std::sync::Arc;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::OwnerProperties;
use kekule::structure::{
    AsModelView, Conformation, ConformationError, Forces, ModelView, Positions, Trajectory,
    TrajectoryFrame, TrajectoryFrameMut, TrajectoryFrameView, Velocities,
};
use kekule::topology::Topology;
use kekule::units::Quantity;

use super::TrajectoryError;

/// Complete borrowed frame state ready for transactional publication.
///
/// A decoder builds this value only after it has read one complete frame into
/// reusable scratch. [`FrameBuffer::replace_from_data`] validates every field
/// before changing the destination, converts units once, reuses dense-array
/// allocations, and clears optional fields omitted from this value.
#[derive(Debug, Clone, Copy)]
pub struct FrameBufferData<'a> {
    positions: Quantity<&'a [Point3]>,
    cell: Option<PeriodicCell>,
    velocities: Option<Quantity<&'a [Vector3]>>,
    forces: Option<Quantity<&'a [Vector3]>>,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
    owner: Option<&'a OwnerProperties>,
    annotations: Option<&'a Conformation>,
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
            owner: None,
            annotations: None,
        }
    }

    /// Borrows every field, including occupancies, B-factors, and all
    /// realization properties, from a topology-bound frame.
    pub fn from_frame_view(frame: TrajectoryFrameView<'a>) -> Self {
        let frame = frame.payload();
        Self {
            positions: frame.positions().values(),
            cell: frame.cell().copied(),
            velocities: frame.velocities().map(Velocities::values),
            forces: frame.forces().map(Forces::values),
            time: frame.time(),
            step: frame.step(),
            owner: None,
            annotations: Some(frame.conformation()),
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

    /// Frame-level owner annotations, such as a decoded format scalar.
    pub const fn with_owner_properties(mut self, owner: &'a OwnerProperties) -> Self {
        self.owner = Some(owner);
        self
    }
}

/// Reusable caller-owned frame storage bound to one topology snapshot.
///
/// It accepts frames and readers of any snapshot sharing that layout. Frame
/// state reads through `Deref` (`buffer.positions()`, `buffer.velocities()`,
/// ...); edit it through [`Self::frame_mut`] or publish a complete frame with
/// [`Self::replace_from_data`]. Velocity and force allocations are reused
/// across frames.
#[derive(Debug, Clone)]
pub struct FrameBuffer {
    /// Exactly one frame; the collection binds and validates it.
    storage: Trajectory,
    spare_velocities: Option<Velocities>,
    spare_forces: Option<Forces>,
}

impl std::ops::Deref for FrameBuffer {
    type Target = TrajectoryFrame;

    fn deref(&self) -> &Self::Target {
        self.frame_view().payload()
    }
}

impl AsModelView for FrameBuffer {
    fn as_model_view(&self) -> ModelView<'_> {
        FrameBuffer::as_model_view(self)
    }
}

impl FrameBuffer {
    pub fn new(topology: Arc<Topology>) -> Self {
        let frame = TrajectoryFrame::new(Positions::zeros(topology.atom_count()));
        let storage = Trajectory::from_items(topology, [frame])
            .expect("zero positions match the topology atom count");
        Self {
            storage,
            spare_velocities: None,
            spare_forces: None,
        }
    }

    pub fn topology(&self) -> &Topology {
        self.storage.topology()
    }

    pub fn shared_topology(&self) -> Arc<Topology> {
        self.storage.shared_topology()
    }

    /// The buffered frame bound to its topology.
    pub fn frame_view(&self) -> TrajectoryFrameView<'_> {
        self.storage.get(0).expect("frame buffer holds one frame")
    }

    pub fn as_model_view(&self) -> ModelView<'_> {
        self.frame_view().as_model_view()
    }

    /// Dimension-preserving editor for the buffered frame.
    pub fn frame_mut(&mut self) -> TrajectoryFrameMut<'_> {
        self.storage
            .get_mut(0)
            .expect("frame buffer holds one frame")
    }

    /// Copies velocities into reusable storage, or clears them with `None`.
    pub fn set_velocities<T>(
        &mut self,
        velocities: Option<Quantity<T>>,
    ) -> Result<(), ConformationError>
    where
        T: AsRef<[Vector3]>,
    {
        let staged = velocities
            .map(|values| self.stage_velocities(values))
            .transpose()?;
        self.publish_velocities(staged);
        Ok(())
    }

    /// Copies forces into reusable storage, or clears them with `None`.
    pub fn set_forces<T>(&mut self, forces: Option<Quantity<T>>) -> Result<(), ConformationError>
    where
        T: AsRef<[Vector3]>,
    {
        let staged = forces.map(|values| self.stage_forces(values)).transpose()?;
        self.publish_forces(staged);
        Ok(())
    }

    /// Clears velocities while retaining the reusable backing allocation.
    pub fn clear_velocities(&mut self) {
        self.publish_velocities(None);
    }

    /// Clears forces while retaining the reusable backing allocation.
    pub fn clear_forces(&mut self) {
        self.publish_forces(None);
    }

    /// Clears all per-frame state except positions, keeping allocations and
    /// the bound topology.
    pub fn reset_dynamic_state(&mut self) {
        self.publish_velocities(None);
        self.publish_forces(None);
        let mut frame = self.frame_mut();
        frame.set_step(None);
        frame.set_time(None).expect("clearing time is always valid");
        let mut conformation = frame.conformation_mut();
        conformation.set_cell(None);
        clear_annotations(&mut conformation);
    }

    /// Replaces the complete visible frame transactionally.
    ///
    /// All count, unit, finite-value, and annotation validation completes
    /// before any field changes. Optional fields absent from `data`, including
    /// annotations, are cleared.
    pub fn replace_from_data(
        &mut self,
        data: FrameBufferData<'_>,
    ) -> Result<(), ConformationError> {
        self.positions().validate_all(&data.positions)?;
        let time = validated_time(data.time)?;
        if let Some(source) = data.annotations {
            if source.atom_count() != self.atom_count() {
                return Err(ConformationError::AtomCountMismatch {
                    expected: self.atom_count(),
                    actual: source.atom_count(),
                });
            }
            if source.properties().bonds().len() != self.properties().bonds().len() {
                return Err(ConformationError::BondCountMismatch {
                    expected: self.properties().bonds().len(),
                    actual: source.properties().bonds().len(),
                });
            }
        }
        let velocities = data
            .velocities
            .map(|values| self.stage_velocities(values))
            .transpose()?;
        let forces = match data
            .forces
            .map(|values| self.stage_forces(values))
            .transpose()
        {
            Ok(forces) => forces,
            Err(error) => {
                self.spare_velocities = velocities.or(self.spare_velocities.take());
                return Err(error);
            }
        };

        // Every check passed; publication cannot fail from here on.
        self.publish_velocities(velocities);
        self.publish_forces(forces);
        let mut frame = self.frame_mut();
        frame.set_step(data.step);
        frame.set_time(time).expect("time was validated");
        let mut conformation = frame.conformation_mut();
        conformation
            .set_positions(data.positions)
            .expect("positions were validated");
        conformation.set_cell(data.cell);
        match data.annotations {
            Some(source) => {
                conformation
                    .set_occupancies(source.occupancies().map(<[_]>::to_vec))
                    .expect("source occupancies are valid");
                conformation
                    .set_b_factors(
                        source
                            .b_factors()
                            .map(|values| Quantity::new(values.value().to_vec(), values.unit())),
                    )
                    .expect("source B-factors are valid");
                conformation
                    .set_properties(source.properties().clone())
                    .expect("annotation dimensions were validated");
            }
            None => {
                clear_annotations(&mut conformation);
                if let Some(owner) = data.owner {
                    *conformation.properties_mut().owner_mut() = owner.clone();
                }
            }
        }
        Ok(())
    }

    /// Copies a complete frame of any snapshot sharing this buffer's layout.
    pub fn copy_from(&mut self, frame: TrajectoryFrameView<'_>) -> Result<(), TrajectoryError> {
        if !self.topology().shares_layout(frame.topology()) {
            return Err(TrajectoryError::TopologyMismatch);
        }
        Ok(self.replace_from_data(FrameBufferData::from_frame_view(frame))?)
    }

    fn stage_velocities<T: AsRef<[Vector3]>>(
        &mut self,
        values: Quantity<T>,
    ) -> Result<Velocities, ConformationError> {
        let atoms = self.atom_count();
        let mut storage = self
            .spare_velocities
            .take()
            .unwrap_or_else(|| Velocities::zeros(atoms));
        match storage.set_all(values) {
            Ok(()) => Ok(storage),
            Err(error) => {
                self.spare_velocities = Some(storage);
                Err(error)
            }
        }
    }

    fn stage_forces<T: AsRef<[Vector3]>>(
        &mut self,
        values: Quantity<T>,
    ) -> Result<Forces, ConformationError> {
        let atoms = self.atom_count();
        let mut storage = self
            .spare_forces
            .take()
            .unwrap_or_else(|| Forces::zeros(atoms));
        match storage.set_all(values) {
            Ok(()) => Ok(storage),
            Err(error) => {
                self.spare_forces = Some(storage);
                Err(error)
            }
        }
    }

    /// Swaps staged storage in; the previous allocation becomes the spare.
    fn publish_velocities(&mut self, velocities: Option<Velocities>) {
        let mut frame = self.frame_mut();
        let previous = frame.take_velocities();
        if velocities.is_some() {
            frame
                .set_velocities(velocities)
                .expect("staged velocities have one vector per atom");
        }
        if previous.is_some() {
            self.spare_velocities = previous;
        }
    }

    fn publish_forces(&mut self, forces: Option<Forces>) {
        let mut frame = self.frame_mut();
        let previous = frame.take_forces();
        if forces.is_some() {
            frame
                .set_forces(forces)
                .expect("staged forces have one vector per atom");
        }
        if previous.is_some() {
            self.spare_forces = previous;
        }
    }
}

fn validated_time(time: Option<Quantity<f64>>) -> Result<Option<Quantity<f64>>, ConformationError> {
    let mut probe = TrajectoryFrame::new(Positions::zeros(0));
    probe.set_time(time)?;
    Ok(probe.time())
}

fn clear_annotations(conformation: &mut kekule::structure::ConformationMut<'_>) {
    conformation
        .set_occupancies(None)
        .expect("clearing occupancies is always valid");
    conformation
        .set_b_factors(None)
        .expect("clearing B-factors is always valid");
    conformation.properties_mut().clear();
}
