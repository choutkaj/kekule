//! What a view of one sample reports about it: [`SampleContext`].

use std::fmt;

use crate::properties::OwnerProperties;
use crate::structure::{Distribution, Forces, Velocities};
use crate::units::Quantity;

use super::{Place, Sample, SampleRead, SampleView};

/// What a borrowed sample carries besides its structure: its place in its
/// set's order, its membership in the set's distribution, and its
/// velocities and forces when present.
///
/// [`StructureKind::Sample`](crate::structure::StructureKind::Sample) holds
/// it. Writers and reductions read it so that this state is written,
/// rejected, or deliberately omitted, never silently dropped. A sample of a
/// set without an order has no place, and one outside the set's
/// distribution has no membership. A trajectory frame of a codec
/// ([`FrameView`](crate::traj::FrameView)) has a time place and no
/// membership.
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
    Set(SampleView<'a, dyn SampleRead>),
    Frame(&'a Sample, Place),
}

impl<'a> SampleContext<'a> {
    /// The context of a sample of a set or another reader.
    pub(crate) const fn of_view(view: SampleView<'a, dyn SampleRead>) -> Self {
        Self {
            source: Source::Set(view),
        }
    }

    /// The context of a codec's frame: `sample` at a time `place`.
    pub(crate) const fn of_frame(sample: &'a Sample, place: Place) -> Self {
        Self {
            source: Source::Frame(sample, place),
        }
    }

    /// Whether the sample belongs to a set whose distribution leaves it out.
    pub(crate) fn is_outside_distribution(self) -> bool {
        match self.source {
            Source::Set(view) => view
                .reader
                .distribution()
                .is_some_and(|distribution| distribution.membership(view.index).is_none()),
            Source::Frame(..) => false,
        }
    }

    /// The distribution of the sample's set, whether or not the sample is a
    /// member.
    pub(crate) fn distribution(self) -> Option<&'a Distribution> {
        match self.source {
            Source::Set(view) => view.reader.distribution(),
            Source::Frame(..) => None,
        }
    }

    /// The properties of the sample's set, when its reader has them.
    pub(crate) fn set_properties(self) -> Option<&'a OwnerProperties> {
        match self.source {
            Source::Set(view) => view.reader.properties(),
            Source::Frame(..) => None,
        }
    }

    /// The sample's place in its set's order: a time and step, a step, or a
    /// coordinate in the order's unit. `None` in an unordered set.
    pub fn place(self) -> Option<Place> {
        match self.source {
            Source::Set(view) => view.place(),
            Source::Frame(_, place) => Some(place),
        }
    }

    /// The sample's membership in its set's distribution, or `None` when
    /// the set has none or the sample is not a member.
    pub fn membership(self) -> Option<Membership> {
        match self.source {
            Source::Set(view) => view
                .reader
                .distribution()
                .and_then(|distribution| distribution.membership(view.index)),
            Source::Frame(..) => None,
        }
    }

    /// The sample's velocities, when present.
    pub fn velocities(self) -> Option<&'a Velocities> {
        match self.source {
            Source::Set(view) => view.payload().velocities(),
            Source::Frame(sample, _) => sample.velocities(),
        }
    }

    /// The sample's forces, when present.
    pub fn forces(self) -> Option<&'a Forces> {
        match self.source {
            Source::Set(view) => view.payload().forces(),
            Source::Frame(sample, _) => sample.forces(),
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

    /// The energy of the sample's state, when recorded. It is the state's
    /// energy, never the sample's own.
    pub const fn state_energy(self) -> Option<Quantity<f64>> {
        self.energy
    }

    /// The temperature the distribution's populations describe, when
    /// recorded.
    pub const fn temperature(self) -> Option<Quantity<f64>> {
        self.temperature
    }
}
