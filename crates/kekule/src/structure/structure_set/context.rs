//! What a view of one geometry reports about its set: [`SetContext`].

use std::fmt;

use crate::properties::OwnerProperties;

use super::{GeometryView, Place, StructureSetRead};

/// What a borrowed geometry's set records about it: its place in the set's
/// order, and the set's properties.
///
/// [`StructureKind::InSet`](crate::structure::StructureKind::InSet) holds
/// it. Writers and reductions read it so that this state is written,
/// rejected, or deliberately omitted, never silently dropped. A geometry of a
/// set without an order has no place. A trajectory frame of a codec
/// ([`FrameView`](crate::traj::FrameView)) has a time place. The geometry's
/// own state, velocities and forces included, reads through the
/// [`StructureView`](crate::structure::StructureView).
///
/// The context borrows the geometry's container and reads these on request,
/// so a view stays small.
#[derive(Clone, Copy)]
pub struct SetContext<'a> {
    source: Source<'a>,
}

/// Where a context reads the geometry's place.
#[derive(Clone, Copy)]
enum Source<'a> {
    Set(GeometryView<'a, dyn StructureSetRead>),
    Frame(Place),
}

impl<'a> SetContext<'a> {
    /// The context of a geometry of a set or another reader.
    pub(crate) const fn of_view(view: GeometryView<'a, dyn StructureSetRead>) -> Self {
        Self {
            source: Source::Set(view),
        }
    }

    /// The context of a codec's frame at a time `place`.
    pub(crate) const fn of_frame(place: Place) -> Self {
        Self {
            source: Source::Frame(place),
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
            Source::Frame(place) => Some(place),
        }
    }
}

/// Two contexts are equal when they report the same place, whatever
/// container they borrow.
impl PartialEq for SetContext<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.place() == other.place()
    }
}

impl fmt::Debug for SetContext<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SetContext")
            .field("place", &self.place())
            .finish_non_exhaustive()
    }
}
