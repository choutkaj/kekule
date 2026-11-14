use std::sync::Arc;

use crate::structure::{
    check_bound, AsStructureView, Geometry, GeometryView, Place, SetContext, Structure,
    StructureKind, StructureSetError, StructureView, StructureViewSealed,
};
use crate::topology::Topology;
use crate::units::{Quantity, CANONICAL_TIME_UNIT};

use super::TrajectoryIoError;
use crate::diagnostics::GeometryField;

/// One trajectory frame as codecs read and write it: a topology-bound
/// [`Geometry`] with its time and MD step, borrowed without copying.
///
/// It is the input of [`TrajectoryWriter::write_frame`](super::TrajectoryWriter::write_frame).
/// A [`FrameBuffer`](super::FrameBuffer) lends its frame
/// ([`FrameBuffer::frame_view`](super::FrameBuffer::frame_view)); a geometry
/// of a time-ordered [`StructureSet`](crate::structure::StructureSet) becomes one
/// with [`Self::from_geometry`], and any bound geometry with [`Self::new`]
/// ([`Geometry::bind`] binds one made outside a set).
/// Frame state reads through `Deref` (`frame.positions()`,
/// `frame.velocities()`, ...), and [`AsStructureView`] lends it to kernels as a
/// geometry at a time place.
///
/// ```
/// use kekule::structure::{Order, Place, Positions, Geometry, StructureSet, TimeAxis};
/// use kekule::traj::{FrameView, MemoryTrajectoryWriter, TrajectoryWriter};
/// use kekule::units::{Quantity, PICOSECOND};
///
/// let topology = kekule::smiles::to_topology("O")?;
/// let mut set = StructureSet::new(topology.clone());
/// let times = Quantity::new(Vec::new(), PICOSECOND);
/// set.set_order(Some(Order::Time(TimeAxis::new(Some(times), Some(Vec::new()))?)))?;
/// let place = Place::Time { time: Some(Quantity::new(2.0, PICOSECOND)), step: Some(1000) };
/// set.push(Geometry::new(Positions::zeros(topology.atom_count())), Some(place))?;
///
/// let mut writer = MemoryTrajectoryWriter::new(topology);
/// for geometry in &set {
///     writer.write_frame(FrameView::from_geometry(geometry)?)?;
/// }
/// assert_eq!(writer.into_structure_set().order(), set.order());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy)]
pub struct FrameView<'a> {
    topology: &'a Arc<Topology>,
    geometry: &'a Geometry,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
}

impl std::ops::Deref for FrameView<'_> {
    type Target = Geometry;

    fn deref(&self) -> &Self::Target {
        self.geometry
    }
}

impl<'a> FrameView<'a> {
    /// A frame of `geometry`, bound to `topology`, at `time` and `step`.
    ///
    /// The geometry must be bound to `topology`, as a geometry of a set of it
    /// is. Bind a geometry made with [`Geometry::new`] first with
    /// [`Geometry::bind`]; the frame borrows the geometry and copies nothing.
    ///
    /// # Errors
    ///
    /// [`StructureSetError::Geometry`] when the geometry does not have a
    /// position per atom and a property row per atom and bond of `topology`,
    /// such as an unbound geometry of a topology with bonds
    /// ([`GeometryError::BondCountMismatch`](crate::structure::GeometryError::BondCountMismatch)),
    /// and [`StructureSetError::NonFiniteTime`] or [`StructureSetError::Unit`]
    /// for an invalid time.
    pub fn new(
        topology: &'a Arc<Topology>,
        geometry: &'a Geometry,
        time: Option<Quantity<f64>>,
        step: Option<u64>,
    ) -> Result<Self, StructureSetError> {
        check_bound(topology, geometry)?;
        Ok(Self::bound(topology, geometry, canonical_time(time)?, step))
    }

    /// A frame of a bound geometry at a canonical time.
    pub(crate) const fn bound(
        topology: &'a Arc<Topology>,
        geometry: &'a Geometry,
        time: Option<Quantity<f64>>,
        step: Option<u64>,
    ) -> Self {
        Self {
            topology,
            geometry,
            time,
            step,
        }
    }

    /// The frame of a geometry of a set or another reader: its time and step
    /// on a time order, or neither in an unordered set.
    ///
    /// # Errors
    ///
    /// [`TrajectoryIoErrorKind::UnsupportedField`](super::TrajectoryIoErrorKind::UnsupportedField)
    /// for a geometry that has a place on a step or coordinate order ([`GeometryField::OrderStep`],
    /// [`GeometryField::OrderCoordinate`]): a frame cannot hold them.
    pub fn from_geometry<R: ?Sized>(view: GeometryView<'a, R>) -> Result<Self, TrajectoryIoError> {
        let view = view.as_dyn();
        let (time, step) = match view.place() {
            None => (None, None),
            Some(Place::Time { time, step }) => (time, step),
            Some(Place::Step(_)) => {
                return Err(TrajectoryIoError::unsupported_field(
                    GeometryField::OrderStep,
                ))
            }
            Some(Place::Coordinate(_)) => {
                return Err(TrajectoryIoError::unsupported_field(
                    GeometryField::OrderCoordinate,
                ))
            }
        };
        Ok(Self::bound(
            view.topology_handle(),
            view.payload(),
            time,
            step,
        ))
    }

    /// The topology the frame is bound to, with the view's full lifetime.
    pub fn topology(self) -> &'a Topology {
        self.topology
    }

    /// A shared handle to the bound topology.
    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.topology)
    }

    /// The borrowed geometry with the view's full lifetime.
    pub const fn payload(self) -> &'a Geometry {
        self.geometry
    }

    /// The simulation time in [`CANONICAL_TIME_UNIT`], when known.
    pub const fn time(self) -> Option<Quantity<f64>> {
        self.time
    }

    /// The MD step, when known.
    pub const fn step(self) -> Option<u64> {
        self.step
    }

    /// The frame's place on a time order: its time and step.
    pub const fn place(self) -> Place {
        Place::Time {
            time: self.time,
            step: self.step,
        }
    }

    /// Borrows the frame as the common kernel input, with the view's full
    /// lifetime: a geometry at its time place, without a membership.
    pub fn as_structure_view(self) -> StructureView<'a> {
        self.geometry.structure_view(
            self.topology,
            StructureKind::InSet(SetContext::of_frame(self.place())),
        )
    }

    /// Copies the frame out as an owned structure on the shared topology.
    pub fn to_structure(self) -> Structure {
        self.as_structure_view().to_structure()
    }
}

/// `time` in [`CANONICAL_TIME_UNIT`], checked finite.
pub(super) fn canonical_time(
    time: Option<Quantity<f64>>,
) -> Result<Option<Quantity<f64>>, StructureSetError> {
    time.map(|time| {
        let time = time.into_unit(CANONICAL_TIME_UNIT)?;
        if !time.value().is_finite() {
            return Err(StructureSetError::NonFiniteTime);
        }
        Ok(time)
    })
    .transpose()
}

impl StructureViewSealed for FrameView<'_> {}

impl AsStructureView for FrameView<'_> {
    fn as_structure_view(&self) -> StructureView<'_> {
        FrameView::as_structure_view(*self)
    }
}
