use std::sync::Arc;

use crate::structure::{
    check_bound, AsModelView, Model, ModelView, ModelViewSealed, Place, Sample, SampleContext,
    SampleView, StructureError, StructureKind,
};
use crate::topology::Topology;
use crate::units::{Quantity, CANONICAL_TIME_UNIT};

use super::{TrajectoryField, TrajectoryIoError};

/// One trajectory frame as codecs read and write it: a topology-bound
/// [`Sample`] with its time and MD step, borrowed without copying.
///
/// It is the input of [`TrajectoryWriter::write_frame`](super::TrajectoryWriter::write_frame).
/// A [`FrameBuffer`](super::FrameBuffer) lends its frame
/// ([`FrameBuffer::frame_view`](super::FrameBuffer::frame_view)); a sample
/// of a time-ordered [`SampleSet`](crate::structure::SampleSet) becomes one
/// with [`Self::from_sample`], and any bound sample with [`Self::new`].
/// Frame state reads through `Deref` (`frame.positions()`,
/// `frame.velocities()`, ...), and [`AsModelView`] lends it to kernels as a
/// sample at a time place.
///
/// ```
/// use kekule::structure::{Order, Place, Positions, Sample, SampleSet, TimeAxis};
/// use kekule::traj::{FrameView, MemoryTrajectoryWriter, TrajectoryWriter};
/// use kekule::units::{Quantity, PICOSECOND};
///
/// let topology = kekule::smiles::to_topology("O")?;
/// let mut set = SampleSet::new(topology.clone());
/// let times = Quantity::new(Vec::new(), PICOSECOND);
/// set.set_order(Some(Order::Time(TimeAxis::new(Some(times), Some(Vec::new()))?)))?;
/// let place = Place::Time { time: Some(Quantity::new(2.0, PICOSECOND)), step: Some(1000) };
/// set.push(Sample::new(Positions::zeros(topology.atom_count())), Some(place))?;
///
/// let mut writer = MemoryTrajectoryWriter::new(topology);
/// for sample in &set {
///     writer.write_frame(FrameView::from_sample(sample)?)?;
/// }
/// assert_eq!(writer.into_sample_set().order(), set.order());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, Copy)]
pub struct FrameView<'a> {
    topology: &'a Arc<Topology>,
    sample: &'a Sample,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
}

impl std::ops::Deref for FrameView<'_> {
    type Target = Sample;

    fn deref(&self) -> &Self::Target {
        self.sample
    }
}

impl<'a> FrameView<'a> {
    /// A frame of `sample`, bound to `topology`, at `time` and `step`.
    ///
    /// # Errors
    ///
    /// [`StructureError`] when the sample does not have a position per atom
    /// and a property row per atom and bond of `topology`, as a sample of a
    /// set does, and [`StructureError::NonFiniteTime`] or
    /// [`StructureError::Unit`] for an invalid time.
    pub fn new(
        topology: &'a Arc<Topology>,
        sample: &'a Sample,
        time: Option<Quantity<f64>>,
        step: Option<u64>,
    ) -> Result<Self, StructureError> {
        check_bound(topology, sample)?;
        Ok(Self::bound(topology, sample, canonical_time(time)?, step))
    }

    /// A frame of a bound sample at a canonical time.
    pub(crate) const fn bound(
        topology: &'a Arc<Topology>,
        sample: &'a Sample,
        time: Option<Quantity<f64>>,
        step: Option<u64>,
    ) -> Self {
        Self {
            topology,
            sample,
            time,
            step,
        }
    }

    /// The frame of a sample of a set or another reader: its time and step
    /// on a time order, or neither in an unordered set.
    ///
    /// # Errors
    ///
    /// [`TrajectoryIoErrorKind::UnsupportedField`](super::TrajectoryIoErrorKind::UnsupportedField)
    /// for a sample that is a member of its set's distribution
    /// ([`TrajectoryField::Membership`]), or that has a place on a step or
    /// coordinate order ([`TrajectoryField::OrderStep`],
    /// [`TrajectoryField::OrderCoordinate`]): a frame cannot hold them.
    pub fn from_sample<R: ?Sized>(view: SampleView<'a, R>) -> Result<Self, TrajectoryIoError> {
        let view = view.as_dyn();
        if view.weight().is_some() {
            return Err(TrajectoryIoError::unsupported_field(
                TrajectoryField::Membership,
            ));
        }
        let (time, step) = match view.place() {
            None => (None, None),
            Some(Place::Time { time, step }) => (time, step),
            Some(Place::Step(_)) => {
                return Err(TrajectoryIoError::unsupported_field(
                    TrajectoryField::OrderStep,
                ))
            }
            Some(Place::Coordinate(_)) => {
                return Err(TrajectoryIoError::unsupported_field(
                    TrajectoryField::OrderCoordinate,
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

    /// The borrowed sample with the view's full lifetime.
    pub const fn payload(self) -> &'a Sample {
        self.sample
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
    /// lifetime: a sample at its time place, without a membership.
    pub fn as_model_view(self) -> ModelView<'a> {
        self.sample.structure().model_view(
            self.topology,
            StructureKind::Sample(SampleContext::of_frame(self.sample, self.place())),
        )
    }

    /// Copies the frame out as an owned model on the shared topology.
    pub fn to_model(self) -> Model {
        self.as_model_view().to_model()
    }
}

/// `time` in [`CANONICAL_TIME_UNIT`], checked finite.
pub(super) fn canonical_time(
    time: Option<Quantity<f64>>,
) -> Result<Option<Quantity<f64>>, StructureError> {
    time.map(|time| {
        let time = time.into_unit(CANONICAL_TIME_UNIT)?;
        if !time.value().is_finite() {
            return Err(StructureError::NonFiniteTime);
        }
        Ok(time)
    })
    .transpose()
}

impl ModelViewSealed for FrameView<'_> {}

impl AsModelView for FrameView<'_> {
    fn as_model_view(&self) -> ModelView<'_> {
        FrameView::as_model_view(*self)
    }
}
