//! Borrowed geometries and states of a [`StructureSet`] or another
//! [`StructureSetRead`].

use std::fmt;
use std::marker::PhantomData;
use std::sync::Arc;

use crate::structure::model::sealed::Sealed;
use crate::structure::state::{is_bound, validate_binding};
use crate::structure::{AsModelView, GeometryError, Model, ModelView, State, StructureKind};
use crate::topology::Topology;

use super::{Geometry, Place, SetContext, StructureSet, StructureSetRead, StructureSetReadError};

/// One topology-bound [`Geometry`] borrowed from a [`StructureSet`] or another
/// [`StructureSetRead`], with its place in the reader's order and its membership
/// in the reader's distribution.
///
/// Reads geometry state through `Deref` and converts to a [`ModelView`]
/// through [`AsModelView`] without copying; the view's
/// [`StructureKind::InSet`] carries the geometry's [`SetContext`].
///
/// The type parameter says where the view comes from. `GeometryView<'_>` is a
/// geometry of a set in memory ([`StructureSet::get`]), whose state lends its
/// representative geometry ([`StateView::representative_sample`]).
/// `GeometryView<'_, dyn StructureSetRead>` is a geometry of any reader, as
/// [`StructureSetRead::peek`] lends it; [`Self::as_dyn`] turns the first into the
/// second.
pub struct GeometryView<'a, R: ?Sized = StructureSet> {
    pub(super) reader: &'a dyn StructureSetRead,
    pub(super) geometry: &'a Geometry,
    pub(super) index: usize,
    origin: PhantomData<&'a R>,
}

impl<R: ?Sized> Clone for GeometryView<'_, R> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<R: ?Sized> Copy for GeometryView<'_, R> {}

impl<R: ?Sized> fmt::Debug for GeometryView<'_, R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("GeometryView")
            .field("index", &self.index)
            .field("geometry", self.geometry)
            .finish_non_exhaustive()
    }
}

impl<R: ?Sized> std::ops::Deref for GeometryView<'_, R> {
    type Target = Geometry;

    fn deref(&self) -> &Self::Target {
        self.geometry
    }
}

impl<'a> GeometryView<'a, dyn StructureSetRead> {
    /// A view of `geometry` as geometry `index` of `reader`, for a
    /// [`StructureSetRead::peek`] that lends a geometry the reader holds, such as a
    /// host's cache of geometries read from another reader.
    ///
    /// # Errors
    ///
    /// [`StructureSetReadError::IndexOutOfRange`] past the reader's last geometry,
    /// and [`StructureSetReadError::Geometry`] for a geometry not bound to the
    /// reader's topology: one without a position per atom and a property
    /// row per atom and bond.
    pub fn new(
        reader: &'a dyn StructureSetRead,
        index: usize,
        geometry: &'a Geometry,
    ) -> Result<Self, StructureSetReadError> {
        let len = reader.len();
        if index >= len {
            return Err(StructureSetReadError::IndexOutOfRange { index, len });
        }
        check_bound(reader.topology(), geometry).map_err(|source| {
            StructureSetReadError::Geometry {
                geometry: index,
                source: Box::new(source),
            }
        })?;
        Ok(Self::of_reader(reader, index, geometry))
    }

    /// A view of `geometry`, geometry `index` of `reader` and bound to its
    /// topology.
    pub(crate) fn of_reader(
        reader: &'a dyn StructureSetRead,
        index: usize,
        geometry: &'a Geometry,
    ) -> Self {
        Self {
            reader,
            geometry,
            index,
            origin: PhantomData,
        }
    }
}

impl<'a> GeometryView<'a> {
    /// A view of geometry `index` of an in-memory `set`, which holds `geometry`.
    pub(crate) fn of_set(set: &'a StructureSet, index: usize, geometry: &'a Geometry) -> Self {
        Self {
            reader: set,
            geometry,
            index,
            origin: PhantomData,
        }
    }
}

impl<'a, R: ?Sized> GeometryView<'a, R> {
    /// The same view as one of any reader, as [`StructureSetRead::peek`] returns
    /// it.
    pub fn as_dyn(self) -> GeometryView<'a, dyn StructureSetRead> {
        GeometryView::of_reader(self.reader, self.index, self.geometry)
    }

    /// The geometry's index in its set.
    pub const fn index(self) -> usize {
        self.index
    }

    /// The geometry's place in the set's order, or `None` in an unordered set.
    pub fn place(self) -> Option<Place> {
        self.reader.order().map(|order| order.place(self.index))
    }

    /// The geometry's relative weight in the set's distribution (`1.0` when
    /// every member weighs the same), or `None` when the set has no
    /// distribution or the geometry is not a member. Weighted statistics walk
    /// [`StructureSet::weighted_samples`] instead, which counts every geometry of
    /// a set without a distribution.
    pub fn weight(self) -> Option<f64> {
        self.reader.distribution()?.weight(self.index)
    }

    /// The state the geometry belongs to, or `None` when it is not a member
    /// of the set's distribution.
    pub fn state(self) -> Option<StateView<'a, R>> {
        let state = self.reader.distribution()?.state_of(self.index)?;
        Some(StateView::new(self.reader, state))
    }

    /// Whether the geometry represents its state.
    pub fn is_representative(self) -> bool {
        self.state()
            .is_some_and(|state| state.representative() == self.index)
    }

    /// The topology the geometry is bound to, with the view's full lifetime.
    pub fn topology(self) -> &'a Topology {
        self.reader.topology()
    }

    /// A shared handle to the bound topology.
    pub fn shared_topology(self) -> Arc<Topology> {
        Arc::clone(self.reader.topology())
    }

    /// The bound topology's handle, with the view's full lifetime.
    pub(crate) fn topology_handle(self) -> &'a Arc<Topology> {
        self.reader.topology()
    }

    /// The borrowed geometry with the view's full lifetime.
    pub fn payload(self) -> &'a Geometry {
        self.geometry
    }

    /// Borrows the geometry as the common kernel input, with the view's full
    /// lifetime.
    pub fn as_model_view(self) -> ModelView<'a> {
        self.geometry.structure().model_view(
            self.reader.topology(),
            StructureKind::InSet(SetContext::of_view(self.as_dyn())),
        )
    }

    /// Materializes the geometry as an owned model, dropping its place,
    /// membership, velocities, and forces.
    pub fn to_model(self) -> Model {
        self.as_model_view().to_model()
    }
}

/// Fails unless `geometry` is bound to `topology`.
pub(crate) fn check_bound(topology: &Topology, geometry: &Geometry) -> Result<(), GeometryError> {
    let structure = geometry.structure();
    let atoms = structure.positions.len();
    validate_binding(atoms, &structure.properties, topology)?;
    if !is_bound(atoms, &structure.properties, topology) {
        return Err(GeometryError::BondCountMismatch {
            expected: topology.bond_count(),
            actual: structure.properties.bonds().len(),
        });
    }
    Ok(())
}

impl<R: ?Sized> Sealed for GeometryView<'_, R> {}

impl<R: ?Sized> AsModelView for GeometryView<'_, R> {
    fn as_model_view(&self) -> ModelView<'_> {
        GeometryView::as_model_view(*self)
    }
}

/// One state of the distribution of a [`StructureSet`] or another
/// [`StructureSetRead`], with its members, representative, and population.
///
/// Reads the state's energy and properties through `Deref`. The type
/// parameter says where the state comes from, as for [`GeometryView`]: a
/// state of a set in memory, `StateView<'_>`, lends its representative
/// geometry ([`Self::representative_sample`]); the state of a set in a file,
/// or of another reader, gives its index ([`Self::representative`]), which
/// the host reads.
pub struct StateView<'a, R: ?Sized = StructureSet> {
    reader: &'a dyn StructureSetRead,
    index: usize,
    origin: PhantomData<&'a R>,
}

impl<R: ?Sized> Clone for StateView<'_, R> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<R: ?Sized> Copy for StateView<'_, R> {}

impl<R: ?Sized> fmt::Debug for StateView<'_, R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("StateView")
            .field("index", &self.index)
            .field("state", self.payload())
            .finish_non_exhaustive()
    }
}

impl<R: ?Sized> std::ops::Deref for StateView<'_, R> {
    type Target = State;

    fn deref(&self) -> &Self::Target {
        self.payload()
    }
}

impl<'a, R: ?Sized> StateView<'a, R> {
    /// State `index` of `reader`'s distribution, which `R` describes.
    pub(crate) fn new(reader: &'a dyn StructureSetRead, index: usize) -> Self {
        Self {
            reader,
            index,
            origin: PhantomData,
        }
    }

    fn distribution(self) -> &'a super::Distribution {
        self.reader
            .distribution()
            .expect("a state view's reader has a distribution")
    }

    /// The state's index in the distribution.
    pub const fn index(self) -> usize {
        self.index
    }

    /// The geometry index of the state's representative.
    pub fn representative(self) -> usize {
        let distribution = self.distribution();
        distribution.member(distribution.partition().representatives()[self.index])
    }

    /// The geometry indices of the state's members, in geometry order.
    pub fn geometries(self) -> impl Iterator<Item = usize> + 'a {
        let distribution = self.distribution();
        distribution
            .partition()
            .samples(self.index)
            .map(move |position| distribution.member(position))
    }

    /// The population: the total weight of the state's members, divided by
    /// the largest member weight so that it cannot overflow. Only its ratio
    /// to other populations carries meaning;
    /// [`super::Distribution::state_probabilities`] gives every state's share
    /// of the total.
    pub fn population(self) -> f64 {
        self.distribution().population(self.index)
    }

    /// The borrowed state with the view's full lifetime.
    pub fn payload(self) -> &'a State {
        self.distribution().state(self.index)
    }
}

impl<'a> StateView<'a> {
    /// The state's representative geometry.
    pub fn representative_sample(self) -> GeometryView<'a> {
        let index = self.representative();
        let view = self
            .reader
            .peek(index)
            .expect("a set in memory lends every geometry");
        GeometryView {
            reader: view.reader,
            geometry: view.geometry,
            index,
            origin: PhantomData,
        }
    }
}

/// Iterator over the topology-bound geometries of a [`StructureSet`], in set
/// order.
#[derive(Debug, Clone)]
pub struct GeometryIter<'a> {
    set: &'a StructureSet,
    indices: std::ops::Range<usize>,
}

impl<'a> GeometryIter<'a> {
    pub(crate) fn new(set: &'a StructureSet) -> Self {
        Self {
            set,
            indices: 0..set.len(),
        }
    }
}

impl<'a> Iterator for GeometryIter<'a> {
    type Item = GeometryView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let set = self.set;
        self.indices.next().and_then(|index| set.get(index))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.indices.size_hint()
    }
}

impl ExactSizeIterator for GeometryIter<'_> {}

/// Iterator over the states of a set's distribution, in state order; empty
/// without a distribution. Its type parameter is the set's, as for
/// [`StateView`].
pub struct StateIter<'a, R: ?Sized = StructureSet> {
    reader: &'a dyn StructureSetRead,
    indices: std::ops::Range<usize>,
    origin: PhantomData<&'a R>,
}

impl<R: ?Sized> Clone for StateIter<'_, R> {
    fn clone(&self) -> Self {
        Self {
            reader: self.reader,
            indices: self.indices.clone(),
            origin: PhantomData,
        }
    }
}

impl<R: ?Sized> fmt::Debug for StateIter<'_, R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("StateIter")
            .field("indices", &self.indices)
            .finish_non_exhaustive()
    }
}

impl<'a, R: ?Sized> StateIter<'a, R> {
    /// The states of `reader`'s distribution, which `R` describes.
    pub(crate) fn new(reader: &'a dyn StructureSetRead) -> Self {
        let states = reader
            .distribution()
            .map_or(0, super::Distribution::state_count);
        Self {
            reader,
            indices: 0..states,
            origin: PhantomData,
        }
    }
}

impl<'a, R: ?Sized> Iterator for StateIter<'a, R> {
    type Item = StateView<'a, R>;

    fn next(&mut self) -> Option<Self::Item> {
        let reader = self.reader;
        self.indices
            .next()
            .map(|index| StateView::new(reader, index))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.indices.size_hint()
    }
}

impl<R: ?Sized> ExactSizeIterator for StateIter<'_, R> {}
