use std::sync::Arc;

use crate::math::{PeriodicCell, Point3, Vector3};
use crate::properties::{GeometryProperties, GeometryPropertiesMut, OwnerProperties};
use crate::structure::scalars::AtomScalars;
use crate::structure::{
    DenseArrayError, Forces, Geometry, GeometryError, GeometryMut, Place, Positions, StructureView,
    Velocities,
};
use crate::topology::Topology;
use crate::traj::TrajectoryIoErrorKind;
use crate::units::Quantity;

use super::frame::canonical_time;
use super::{FrameView, TrajectoryIoError};

/// Complete borrowed frame state ready for transactional publication.
///
/// A decoder builds this value only after it has read one complete frame into
/// reusable scratch. [`FrameBuffer::replace_from_data`] validates every field
/// before changing the destination, converts units once, reuses dense-array
/// allocations, and clears optional fields omitted from this value. Only
/// [`Self::from_frame_view`] carries occupancies and B-factors, so a decoded
/// frame has none.
#[derive(Debug, Clone, Copy)]
pub struct FrameData<'a> {
    positions: Quantity<&'a [Point3]>,
    cell: Option<PeriodicCell>,
    velocities: Option<Quantity<&'a [Vector3]>>,
    forces: Option<Quantity<&'a [Vector3]>>,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
    owner: Option<&'a OwnerProperties>,
    annotations: Option<&'a GeometryProperties>,
    /// Occupancies and B-factors, shared rather than copied.
    sites: Option<(&'a AtomScalars, &'a AtomScalars)>,
}

impl<'a> FrameData<'a> {
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
            sites: None,
        }
    }

    /// Borrows every field, including occupancies, B-factors, and all
    /// geometry properties, from a topology-bound frame.
    pub fn from_frame_view(frame: FrameView<'a>) -> Self {
        let geometry = frame.payload();
        let state = geometry.state();
        Self {
            positions: geometry.positions().values(),
            cell: geometry.cell().copied(),
            velocities: geometry.velocities().map(Velocities::values),
            forces: geometry.forces().map(Forces::values),
            time: frame.time(),
            step: frame.step(),
            owner: None,
            annotations: Some(geometry.properties()),
            sites: Some((&state.occupancies, &state.b_factors)),
        }
    }

    /// Sets the periodic cell.
    pub const fn with_cell(mut self, cell: PeriodicCell) -> Self {
        self.cell = Some(cell);
        self
    }

    /// Sets per-atom velocities, one per atom in dense order; any velocity
    /// unit is accepted and converted on publication.
    pub const fn with_velocities(mut self, velocities: Quantity<&'a [Vector3]>) -> Self {
        self.velocities = Some(velocities);
        self
    }

    /// Sets per-atom forces, one per atom in dense order; any force unit is
    /// accepted and converted on publication.
    pub const fn with_forces(mut self, forces: Quantity<&'a [Vector3]>) -> Self {
        self.forces = Some(forces);
        self
    }

    /// Sets the simulation time; any time unit is accepted and converted on
    /// publication.
    pub const fn with_time(mut self, time: Quantity<f64>) -> Self {
        self.time = Some(time);
        self
    }

    /// Sets the integration step.
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

/// Reusable caller-owned frame storage bound to one topology snapshot: one
/// [`Geometry`] with its time and MD step.
///
/// It accepts frames and readers of any snapshot sharing that layout. Frame
/// state reads through `Deref` (`buffer.positions()`, `buffer.velocities()`,
/// ...) and [`Self::time`] and [`Self::step`]; edit it through
/// [`Self::geometry_mut`], [`Self::set_time`], and [`Self::set_step`], or
/// publish a complete frame with [`Self::replace_from_data`]. Velocity and
/// force allocations are reused across frames. [`Self::frame_view`] lends
/// the frame to a [`TrajectoryWriter`](super::TrajectoryWriter).
///
/// To keep a decoded frame, for example in a structure set or a frame cache,
/// take it out with [`Self::into_geometry`] rather than copying it: the geometry
/// is bound to the buffer's topology, and its place goes with it into a
/// time-ordered set ([`StructureSet::push`](crate::structure::StructureSet::push)).
/// A geometry that is no longer needed becomes a buffer again with
/// [`Self::from_geometry`], so the next frame decodes into it without
/// allocating. A buffer whose geometry is shared, such as one cloned out of
/// it, copies the geometry once on the next edit.
#[derive(Debug, Clone)]
pub struct FrameBuffer {
    topology: Arc<Topology>,
    /// Bound to `topology`.
    geometry: Geometry,
    /// In the canonical time unit.
    time: Option<Quantity<f64>>,
    step: Option<u64>,
    spare_velocities: Option<Velocities>,
    spare_forces: Option<Forces>,
}

impl std::ops::Deref for FrameBuffer {
    type Target = Geometry;

    fn deref(&self) -> &Self::Target {
        &self.geometry
    }
}

impl FrameBuffer {
    /// A buffer bound to `topology` holding one frame with every position at
    /// the origin and no optional state.
    pub fn new(topology: impl Into<Arc<Topology>>) -> Self {
        let topology = topology.into();
        let geometry = Geometry::new(Positions::zeros(topology.atom_count()));
        Self::from_geometry(topology, geometry)
            .expect("zero positions match the topology atom count")
    }

    /// A buffer holding `geometry` without copying it, bound to `topology`,
    /// without a time or step.
    ///
    /// The geometry is validated and its bond property rows bound once, as
    /// [`StructureSet::push`](crate::structure::StructureSet::push) does; it is
    /// the buffered frame until the next read replaces it, reusing its
    /// allocations.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidGeometry`] when the geometry's
    /// atom count or property rows do not match `topology`.
    pub fn from_geometry(
        topology: impl Into<Arc<Topology>>,
        mut geometry: Geometry,
    ) -> Result<Self, TrajectoryIoError> {
        let topology = topology.into();
        geometry.bind(&topology)?;
        Ok(Self {
            topology,
            geometry,
            time: None,
            step: None,
            spare_velocities: None,
            spare_forces: None,
        })
    }

    /// The buffered frame, without copying it: its geometry, bound to the
    /// buffer's topology, and its place on a time order.
    ///
    /// The geometry keeps every field of the buffered frame. Spare velocity
    /// and force allocations are dropped; clear velocities or forces first
    /// to drop those of the frame as well.
    pub fn into_geometry(self) -> (Geometry, Place) {
        let place = Place::Time {
            time: self.time,
            step: self.step,
        };
        (self.geometry, place)
    }

    /// The topology snapshot the buffer is bound to.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// A shared handle to [`Self::topology`].
    pub fn shared_topology(&self) -> Arc<Topology> {
        Arc::clone(&self.topology)
    }

    /// The buffered frame bound to its topology.
    pub fn frame_view(&self) -> FrameView<'_> {
        FrameView::bound(&self.topology, &self.geometry, self.time, self.step)
    }

    /// The buffered frame as a [`StructureView`], so structure analyses apply to it.
    pub fn as_structure_view(&self) -> StructureView<'_> {
        self.frame_view().as_structure_view()
    }

    /// The simulation time in the canonical time unit, when known.
    pub const fn time(&self) -> Option<Quantity<f64>> {
        self.time
    }

    /// The MD step, when known.
    pub const fn step(&self) -> Option<u64> {
        self.step
    }

    /// Dimension-preserving editor for the buffered geometry.
    pub fn geometry_mut(&mut self) -> GeometryMut<'_> {
        GeometryMut::new(&mut self.geometry)
    }

    /// Sets or clears a finite time, converted to the canonical time unit.
    ///
    /// # Errors
    ///
    /// [`GeometryError::NonFiniteTime`] or [`GeometryError::Unit`] for an
    /// invalid time; the buffer is then unchanged.
    pub fn set_time(&mut self, time: Option<Quantity<f64>>) -> Result<(), GeometryError> {
        self.time = canonical_time(time)?;
        Ok(())
    }

    /// Sets or clears the MD step.
    pub fn set_step(&mut self, step: Option<u64>) {
        self.step = step;
    }

    /// Copies velocities into reusable storage, or clears them with `None`.
    ///
    /// # Errors
    ///
    /// Returns [`GeometryError::Array`] when the count differs from the atom
    /// count, a value is not finite, or the unit is not a velocity; the
    /// buffer is then unchanged.
    pub fn set_velocities<T>(
        &mut self,
        velocities: Option<Quantity<T>>,
    ) -> Result<(), GeometryError>
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
    ///
    /// # Errors
    ///
    /// Returns [`GeometryError::Array`] when the count differs from the atom
    /// count, a value is not finite, or the unit is not a force; the buffer is
    /// then unchanged.
    pub fn set_forces<T>(&mut self, forces: Option<Quantity<T>>) -> Result<(), GeometryError>
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
    /// the bound topology: velocities, forces, time, step, cell, occupancies,
    /// B-factors, and properties.
    pub fn reset_dynamic_state(&mut self) {
        self.publish_velocities(None);
        self.publish_forces(None);
        self.step = None;
        self.time = None;
        let state = self.geometry.state_in_place();
        state.cell = None;
        state.occupancies = AtomScalars::default();
        state.b_factors = AtomScalars::default();
        GeometryPropertiesMut::new(&mut state.properties).clear();
    }

    /// Replaces the complete visible frame transactionally.
    ///
    /// All count, unit, finite-value, and annotation validation completes
    /// before any field changes. Optional fields absent from `data`, including
    /// annotations, occupancies, and B-factors, are cleared.
    ///
    /// # Errors
    ///
    /// Returns [`GeometryError::Array`] when positions, velocities, forces or
    /// atom annotations have the wrong count, a non-finite value or a unit of
    /// the wrong dimension, [`GeometryError::BondCountMismatch`] when bond
    /// annotations have the wrong row count, and
    /// [`GeometryError::NonFiniteTime`] or [`GeometryError::Unit`] for an
    /// invalid time. On error the buffer is unchanged.
    pub fn replace_from_data(&mut self, data: FrameData<'_>) -> Result<(), GeometryError> {
        self.positions().validate_all(&data.positions)?;
        let time = canonical_time(data.time)?;
        if let Some(source) = data.annotations {
            if source.atoms().len() != self.atom_count() {
                return Err(GeometryError::Array(DenseArrayError::CountMismatch {
                    expected: self.atom_count(),
                    actual: source.atoms().len(),
                }));
            }
            if source.bonds().len() != self.properties().bonds().len() {
                return Err(GeometryError::BondCountMismatch {
                    expected: self.properties().bonds().len(),
                    actual: source.bonds().len(),
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
        self.step = data.step;
        self.time = time;
        let state = self.geometry.state_in_place();
        state
            .positions
            .set_all(data.positions)
            .expect("positions were validated");
        state.cell = data.cell;
        match data.sites {
            Some((occupancies, b_factors)) => {
                state.occupancies.clone_from(occupancies);
                state.b_factors.clone_from(b_factors);
            }
            None => {
                state.occupancies = AtomScalars::default();
                state.b_factors = AtomScalars::default();
            }
        }
        match data.annotations {
            Some(source) => state.properties.clone_from(source),
            None => {
                let mut properties = GeometryPropertiesMut::new(&mut state.properties);
                properties.clear();
                if let Some(owner) = data.owner {
                    *properties.owner_mut() = owner.clone();
                }
            }
        }
        Ok(())
    }

    /// Copies a complete frame of any snapshot sharing this buffer's layout,
    /// every field of it: what the frame lacks, the buffer then lacks too.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::TopologyMismatch`] when the frame's
    /// topology does not share this buffer's layout; the buffer is then
    /// unchanged.
    pub fn copy_from(&mut self, frame: FrameView<'_>) -> Result<(), TrajectoryIoError> {
        if !self.topology().shares_layout(frame.topology()) {
            return Err(TrajectoryIoError::new(
                TrajectoryIoErrorKind::TopologyMismatch,
            ));
        }
        Ok(self.replace_from_data(FrameData::from_frame_view(frame))?)
    }

    fn stage_velocities<T: AsRef<[Vector3]>>(
        &mut self,
        values: Quantity<T>,
    ) -> Result<Velocities, GeometryError> {
        let atoms = self.atom_count();
        let mut storage = self
            .spare_velocities
            .take()
            .unwrap_or_else(|| Velocities::zeros(atoms));
        match storage.set_all(values) {
            Ok(()) => Ok(storage),
            Err(error) => {
                self.spare_velocities = Some(storage);
                Err(error.into())
            }
        }
    }

    fn stage_forces<T: AsRef<[Vector3]>>(
        &mut self,
        values: Quantity<T>,
    ) -> Result<Forces, GeometryError> {
        let atoms = self.atom_count();
        let mut storage = self
            .spare_forces
            .take()
            .unwrap_or_else(|| Forces::zeros(atoms));
        match storage.set_all(values) {
            Ok(()) => Ok(storage),
            Err(error) => {
                self.spare_forces = Some(storage);
                Err(error.into())
            }
        }
    }

    /// Swaps staged storage in; the previous allocation becomes the spare.
    fn publish_velocities(&mut self, velocities: Option<Velocities>) {
        if velocities.is_none() && self.geometry.velocities().is_none() {
            return;
        }
        let previous =
            std::mem::replace(&mut self.geometry.state_in_place().velocities, velocities);
        if previous.is_some() {
            self.spare_velocities = previous;
        }
    }

    fn publish_forces(&mut self, forces: Option<Forces>) {
        if forces.is_none() && self.geometry.forces().is_none() {
            return;
        }
        let previous = std::mem::replace(&mut self.geometry.state_in_place().forces, forces);
        if previous.is_some() {
            self.spare_forces = previous;
        }
    }
}
