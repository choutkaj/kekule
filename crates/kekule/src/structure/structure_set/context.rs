//! What a view of one geometry reports about its set: [`SetContext`].

use std::fmt;

use crate::properties::{OwnerProperties, PropertyTable};

use super::{GeometryIndex, GeometryView, Place, Sampling, StructureSetRead};

/// What a borrowed geometry's set records about it: its place in the set's
/// order, and the set's properties, columns, and sampling record.
///
/// [`StructureKind::InSet`](crate::structure::StructureKind::InSet) holds
/// it. Writers and reductions read it so that this state is written,
/// rejected, or deliberately omitted, never silently dropped. A geometry of a
/// set without an order has no place. The geometry's own state, velocities
/// and forces included, reads through the
/// [`StructureView`](crate::structure::StructureView).
#[derive(Clone, Copy)]
pub struct SetContext<'a> {
    place: Option<Place>,
    properties: Option<&'a OwnerProperties>,
    columns: Option<&'a PropertyTable<GeometryIndex>>,
    sampling: Option<&'a Sampling>,
}

impl<'a> SetContext<'a> {
    /// The context of a geometry at `place` in a set without properties,
    /// columns, or a sampling record.
    pub(crate) const fn new(place: Option<Place>) -> Self {
        Self {
            place,
            properties: None,
            columns: None,
            sampling: None,
        }
    }

    /// The context of a geometry of a set or another reader.
    pub(crate) fn of_view(view: GeometryView<'a, dyn StructureSetRead>) -> Self {
        Self {
            place: view.place(),
            properties: view.reader.properties(),
            columns: view.reader.columns(),
            sampling: view.reader.sampling(),
        }
    }

    /// The properties of the geometry's set, when its reader has them.
    pub(crate) const fn set_properties(self) -> Option<&'a OwnerProperties> {
        self.properties
    }

    /// The columns of the geometry's set, when its reader has them.
    pub(crate) const fn set_columns(self) -> Option<&'a PropertyTable<GeometryIndex>> {
        self.columns
    }

    /// The sampling record of the geometry's set, when its reader has one.
    pub(crate) const fn set_sampling(self) -> Option<&'a Sampling> {
        self.sampling
    }

    /// The geometry's place in its set's order: a time and step, a step, or a
    /// coordinate in the order's unit. `None` in an unordered set.
    pub const fn place(self) -> Option<Place> {
        self.place
    }
}

/// Two contexts are equal when they report the same place, whatever set
/// properties, columns, and sampling record they borrow.
impl PartialEq for SetContext<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.place == other.place
    }
}

impl fmt::Debug for SetContext<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SetContext")
            .field("place", &self.place)
            .finish_non_exhaustive()
    }
}
