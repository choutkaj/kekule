//! What a view of one geometry reports about it: [`SetContext`].

use std::fmt;

use crate::properties::OwnerProperties;
use crate::structure::{Forces, Velocities};

use super::{Geometry, GeometryView, Place, StructureSetRead};

/// What a borrowed geometry carries besides its state: its place in its
/// set's order, and its velocities and forces when present.
///
/// [`StructureKind::InSet`](crate::structure::StructureKind::InSet) holds
/// it. Writers and reductions read it so that this state is written,
/// rejected, or deliberately omitted, never silently dropped. A geometry of a
/// set without an order has no place. A trajectory frame of a codec
/// ([`FrameView`](crate::traj::FrameView)) has a time place.
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

/// Two contexts are equal when they report the same place, velocities, and
/// forces, whatever container they borrow.
impl PartialEq for SetContext<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.place() == other.place()
            && self.velocities() == other.velocities()
            && self.forces() == other.forces()
    }
}

impl fmt::Debug for SetContext<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SetContext")
            .field("place", &self.place())
            .field("velocities", &self.velocities())
            .field("forces", &self.forces())
            .finish()
    }
}
