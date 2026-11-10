//! What a view of one sample reports about it: [`SampleContext`].

use std::fmt;

use crate::structure::series::ItemPlace;
use crate::structure::{Ensemble, Forces, SeriesItem, TrajectoryFrame, Velocities};
use crate::units::Quantity;

use super::{Place, SampleSet};

/// What a borrowed sample carries besides its structure: its place in its
/// set's order, its membership in the set's distribution, and its
/// velocities and forces when present.
///
/// [`StructureKind::Sample`](crate::structure::StructureKind::Sample) holds
/// it. Writers and reductions read it so that this state is written,
/// rejected, or deliberately omitted, never silently dropped. A sample of a
/// set without an order has no place, and one outside the set's
/// distribution has no membership. Views of an ensemble sample, a trajectory
/// frame, and a series item give the same context: membership, a time place
/// with velocities and forces, and a step or coordinate place.
///
/// The context borrows the sample's container and reads these on request,
/// so a view stays small.
#[derive(Clone, Copy)]
pub struct SampleContext<'a> {
    source: Source<'a>,
}

/// Where a context reads the sample's state.
#[derive(Clone, Copy)]
enum Source<'a> {
    Set(&'a SampleSet, usize),
    Ensemble(&'a Ensemble, usize),
    Frame(&'a TrajectoryFrame),
    Item(&'a SeriesItem),
}

impl<'a> SampleContext<'a> {
    /// The context of sample `index` of `set`, which must exist.
    pub(crate) const fn of_set(set: &'a SampleSet, index: usize) -> Self {
        Self {
            source: Source::Set(set, index),
        }
    }

    /// The context of sample `index` of `ensemble`, which must exist.
    pub(crate) const fn of_ensemble(ensemble: &'a Ensemble, index: usize) -> Self {
        Self {
            source: Source::Ensemble(ensemble, index),
        }
    }

    /// The context of a trajectory frame.
    pub(crate) const fn of_frame(frame: &'a TrajectoryFrame) -> Self {
        Self {
            source: Source::Frame(frame),
        }
    }

    /// The context of a series item.
    pub(crate) const fn of_item(item: &'a SeriesItem) -> Self {
        Self {
            source: Source::Item(item),
        }
    }

    /// The sample's place in its set's order: a time and step, a step, or a
    /// coordinate in the order's unit. `None` in an unordered set.
    pub fn place(self) -> Option<Place> {
        match self.source {
            Source::Set(set, index) => set.order().map(|order| order.place(index)),
            Source::Ensemble(..) => None,
            Source::Frame(frame) => Some(Place::Time {
                time: frame.time,
                step: frame.step,
            }),
            Source::Item(item) => Some(match item.place {
                ItemPlace::Step(step) => Place::Step(step),
                ItemPlace::Coordinate(coordinate) => Place::Coordinate(coordinate),
            }),
        }
    }

    /// The sample's membership in its set's distribution, or `None` when
    /// the set has none or the sample is not a member.
    pub fn membership(self) -> Option<Membership> {
        match self.source {
            Source::Set(set, index) => set
                .distribution()
                .and_then(|distribution| distribution.membership(index)),
            Source::Ensemble(ensemble, index) => Some(ensemble.membership(index)),
            Source::Frame(_) | Source::Item(_) => None,
        }
    }

    /// The sample's velocities, when present.
    pub fn velocities(self) -> Option<&'a Velocities> {
        match self.source {
            Source::Set(set, index) => set.samples[index].velocities(),
            Source::Frame(frame) => frame.velocities.as_ref(),
            Source::Ensemble(..) | Source::Item(_) => None,
        }
    }

    /// The sample's forces, when present.
    pub fn forces(self) -> Option<&'a Forces> {
        match self.source {
            Source::Set(set, index) => set.samples[index].forces(),
            Source::Frame(frame) => frame.forces.as_ref(),
            Source::Ensemble(..) | Source::Item(_) => None,
        }
    }
}

/// Two contexts are equal when they report the same place, membership,
/// velocities, and forces, whatever container they borrow.
impl PartialEq for SampleContext<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.place() == other.place()
            && self.membership() == other.membership()
            && self.velocities() == other.velocities()
            && self.forces() == other.forces()
    }
}

impl fmt::Debug for SampleContext<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SampleContext")
            .field("place", &self.place())
            .field("membership", &self.membership())
            .field("velocities", &self.velocities())
            .field("forces", &self.forces())
            .finish()
    }
}

/// A sample's membership in a distribution: its weight and state, and what
/// the distribution records about that state.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Membership {
    weight: f64,
    state: usize,
    representative: bool,
    energy: Option<Quantity<f64>>,
    temperature: Option<Quantity<f64>>,
}

impl Membership {
    pub(crate) const fn new(
        weight: f64,
        state: usize,
        representative: bool,
        energy: Option<Quantity<f64>>,
        temperature: Option<Quantity<f64>>,
    ) -> Self {
        Self {
            weight,
            state,
            representative,
            energy,
            temperature,
        }
    }

    /// The sample's finite, positive relative weight; `1.0` when the
    /// distribution's weights are equal.
    pub const fn weight(self) -> f64 {
        self.weight
    }

    /// The index of the sample's state.
    pub const fn state(self) -> usize {
        self.state
    }

    /// Whether the sample represents its state.
    pub const fn is_representative(self) -> bool {
        self.representative
    }

    /// The energy of the sample's state, when recorded.
    pub const fn energy(self) -> Option<Quantity<f64>> {
        self.energy
    }

    /// The temperature the distribution's populations describe, when
    /// recorded.
    pub const fn temperature(self) -> Option<Quantity<f64>> {
        self.temperature
    }
}
