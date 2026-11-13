//! What a view of one geometry reports about it: [`SetContext`].

use std::fmt;

use crate::properties::OwnerProperties;
use crate::structure::{Distribution, Forces, Velocities};
use crate::units::Quantity;

use super::{Geometry, GeometryView, Place, StructureSetRead};

/// What a borrowed geometry carries besides its structure: its place in its
/// set's order, its membership in the set's distribution, and its
/// velocities and forces when present.
///
/// [`StructureKind::InSet`](crate::structure::StructureKind::InSet) holds
/// it. Writers and reductions read it so that this state is written,
/// rejected, or deliberately omitted, never silently dropped. A geometry of a
/// set without an order has no place, and one outside the set's
/// distribution has no membership. A trajectory frame of a codec
/// ([`FrameView`](crate::traj::FrameView)) has a time place and no
/// membership.
///
/// The context borrows the geometry's container and reads these on request,
/// so a view stays small.
#[derive(Clone, Copy)]
pub struct SetContext<'a> {
    source: Source<'a>,
}

/// Where a context reads the geometry's state.
#[derive(Clone, Copy)]
enum Source<'a> {
    Set(GeometryView<'a, dyn StructureSetRead>),
    Frame(&'a Geometry, Place),
}

impl<'a> SetContext<'a> {
    /// The context of a geometry of a set or another reader.
    pub(crate) const fn of_view(view: GeometryView<'a, dyn StructureSetRead>) -> Self {
        Self {
            source: Source::Set(view),
        }
    }

    /// The context of a codec's frame: `geometry` at a time `place`.
    pub(crate) const fn of_frame(geometry: &'a Geometry, place: Place) -> Self {
        Self {
            source: Source::Frame(geometry, place),
        }
    }

    /// Whether the geometry belongs to a set whose distribution leaves it out.
    pub(crate) fn is_outside_distribution(self) -> bool {
        match self.source {
            Source::Set(view) => view
                .reader
                .distribution()
                .is_some_and(|distribution| distribution.membership(view.index).is_none()),
            Source::Frame(..) => false,
        }
    }

    /// The distribution of the geometry's set, whether or not the geometry is a
    /// member.
    pub(crate) fn distribution(self) -> Option<&'a Distribution> {
        match self.source {
            Source::Set(view) => view.reader.distribution(),
            Source::Frame(..) => None,
        }
    }

    /// The properties of the geometry's set, when its reader has them.
    pub(crate) fn set_properties(self) -> Option<&'a OwnerProperties> {
        match self.source {
            Source::Set(view) => view.reader.properties(),
            Source::Frame(..) => None,
        }
    }

    /// The geometry's place in its set's order: a time and step, a step, or a
    /// coordinate in the order's unit. `None` in an unordered set.
    pub fn place(self) -> Option<Place> {
        match self.source {
            Source::Set(view) => view.place(),
            Source::Frame(_, place) => Some(place),
        }
    }

    /// The geometry's membership in its set's distribution, or `None` when
    /// the set has none or the geometry is not a member.
    pub fn membership(self) -> Option<Membership> {
        match self.source {
            Source::Set(view) => view
                .reader
                .distribution()
                .and_then(|distribution| distribution.membership(view.index)),
            Source::Frame(..) => None,
        }
    }

    /// The geometry's velocities, when present.
    pub fn velocities(self) -> Option<&'a Velocities> {
        match self.source {
            Source::Set(view) => view.payload().velocities(),
            Source::Frame(geometry, _) => geometry.velocities(),
        }
    }

    /// The geometry's forces, when present.
    pub fn forces(self) -> Option<&'a Forces> {
        match self.source {
            Source::Set(view) => view.payload().forces(),
            Source::Frame(geometry, _) => geometry.forces(),
        }
    }
}

/// Two contexts are equal when they report the same place, membership,
/// velocities, and forces, whatever container they borrow.
impl PartialEq for SetContext<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.place() == other.place()
            && self.membership() == other.membership()
            && self.velocities() == other.velocities()
            && self.forces() == other.forces()
    }
}

impl fmt::Debug for SetContext<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SetContext")
            .field("place", &self.place())
            .field("membership", &self.membership())
            .field("velocities", &self.velocities())
            .field("forces", &self.forces())
            .finish()
    }
}

/// A geometry's membership in a distribution: its weight and state, and what
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

    /// The geometry's finite, positive relative weight; `1.0` when the
    /// distribution's weights are equal.
    pub const fn weight(self) -> f64 {
        self.weight
    }

    /// The index of the geometry's state.
    pub const fn state(self) -> usize {
        self.state
    }

    /// Whether the geometry represents its state.
    pub const fn is_representative(self) -> bool {
        self.representative
    }

    /// The energy of the geometry's state, when recorded. It is the state's
    /// energy, never the geometry's own.
    pub const fn state_energy(self) -> Option<Quantity<f64>> {
        self.energy
    }

    /// The temperature the distribution's populations describe, when
    /// recorded.
    pub const fn temperature(self) -> Option<Quantity<f64>> {
        self.temperature
    }
}
