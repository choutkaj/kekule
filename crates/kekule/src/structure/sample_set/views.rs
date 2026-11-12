//! Borrowed samples and states of a [`SampleSet`] or another
//! [`SampleRead`].

use std::fmt;
use std::marker::PhantomData;
use std::sync::Arc;

use crate::structure::model::sealed::Sealed;
use crate::structure::state::{is_bound, validate_binding};
use crate::structure::{AsModelView, Model, ModelView, State, StructureError, StructureKind};
use crate::topology::Topology;

use super::{Place, Sample, SampleContext, SampleRead, SampleReadError, SampleSet};

/// One topology-bound [`Sample`] borrowed from a [`SampleSet`] or another
/// [`SampleRead`], with its place in the reader's order and its membership
/// in the reader's distribution.
///
/// Reads sample state through `Deref` and converts to a [`ModelView`]
/// through [`AsModelView`] without copying; the view's
/// [`StructureKind::Sample`] carries the sample's [`SampleContext`].
///
/// The type parameter says where the view comes from. `SampleView<'_>` is a
/// sample of a set in memory ([`SampleSet::get`]), whose state lends its
/// representative sample ([`StateView::representative_sample`]).
/// `SampleView<'_, dyn SampleRead>` is a sample of any reader, as
/// [`SampleRead::peek`] lends it; [`Self::as_dyn`] turns the first into the
/// second.
pub struct SampleView<'a, R: ?Sized = SampleSet> {
    pub(super) reader: &'a dyn SampleRead,
    pub(super) sample: &'a Sample,
    pub(super) index: usize,
    origin: PhantomData<&'a R>,
}

impl<R: ?Sized> Clone for SampleView<'_, R> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<R: ?Sized> Copy for SampleView<'_, R> {}

impl<R: ?Sized> fmt::Debug for SampleView<'_, R> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SampleView")
            .field("index", &self.index)
            .field("sample", self.sample)
            .finish_non_exhaustive()
    }
}

impl<R: ?Sized> std::ops::Deref for SampleView<'_, R> {
    type Target = Sample;

    fn deref(&self) -> &Self::Target {
        self.sample
    }
}

impl<'a> SampleView<'a, dyn SampleRead> {
    /// A view of `sample` as sample `index` of `reader`, for a
    /// [`SampleRead::peek`] that lends a sample the reader holds, such as a
    /// host's cache of samples read from another reader.
    ///
    /// # Errors
    ///
    /// [`SampleReadError::IndexOutOfRange`] past the reader's last sample,
    /// and [`SampleReadError::Structure`] for a sample not bound to the
    /// reader's topology: one without a position per atom and a property
    /// row per atom and bond.
    pub fn new(
        reader: &'a dyn SampleRead,
        index: usize,
        sample: &'a Sample,
    ) -> Result<Self, SampleReadError> {
        let len = reader.len();
        if index >= len {
            return Err(SampleReadError::IndexOutOfRange { index, len });
        }
        check_bound(reader.topology(), sample).map_err(|source| SampleReadError::Structure {
            sample: index,
            source: Box::new(source),
        })?;
        Ok(Self::of_reader(reader, index, sample))
    }

    /// A view of `sample`, sample `index` of `reader` and bound to its
    /// topology.
    pub(crate) fn of_reader(reader: &'a dyn SampleRead, index: usize, sample: &'a Sample) -> Self {
        Self {
            reader,
            sample,
            index,
            origin: PhantomData,
        }
    }
}

impl<'a> SampleView<'a> {
    /// A view of sample `index` of an in-memory `set`, which holds `sample`.
    pub(crate) fn of_set(set: &'a SampleSet, index: usize, sample: &'a Sample) -> Self {
        Self {
            reader: set,
            sample,
            index,
            origin: PhantomData,
        }
    }
}

impl<'a, R: ?Sized> SampleView<'a, R> {
    /// The same view as one of any reader, as [`SampleRead::peek`] returns
    /// it.
    pub fn as_dyn(self) -> SampleView<'a, dyn SampleRead> {
        SampleView::of_reader(self.reader, self.index, self.sample)
    }

    /// The sample's index in its set.
    pub const fn index(self) -> usize {
        self.index
    }

    /// The sample's place in the set's order, or `None` in an unordered set.
    pub fn place(self) -> Option<Place> {
        self.reader.order().map(|order| order.place(self.index))
    }

    /// The sample's relative weight in the set's distribution (`1.0` when
    /// every member weighs the same), or `None` when the set has no
    /// distribution or the sample is not a member. Weighted statistics walk
    /// [`SampleSet::weighted_samples`] instead, which counts every sample of
    /// a set without a distribution.
    pub fn weight(self) -> Option<f64> {
        self.reader.distribution()?.weight(self.index)
    }

    /// The state the sample belongs to, or `None` when it is not a member
    /// of the set's distribution.
    pub fn state(self) -> Option<StateView<'a, R>> {
        let state = self.reader.distribution()?.state_of(self.index)?;
        Some(StateView::new(self.reader, state))
    }

    /// Whether the sample represents its state.
    pub fn is_representative(self) -> bool {
        self.state()
            .is_some_and(|state| state.representative() == self.index)
    }

    /// The topology the sample is bound to, with the view's full lifetime.
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

    /// The borrowed sample with the view's full lifetime.
    pub fn payload(self) -> &'a Sample {
        self.sample
    }

    /// Borrows the sample as the common kernel input, with the view's full
    /// lifetime.
    pub fn as_model_view(self) -> ModelView<'a> {
        self.sample.structure().model_view(
            self.reader.topology(),
            StructureKind::Sample(SampleContext::of_view(self.as_dyn())),
        )
    }

    /// Materializes the sample as an owned model, dropping its place,
    /// membership, velocities, and forces.
    pub fn to_model(self) -> Model {
        self.as_model_view().to_model()
    }
}

/// Fails unless `sample` is bound to `topology`.
pub(crate) fn check_bound(topology: &Topology, sample: &Sample) -> Result<(), StructureError> {
    let structure = sample.structure();
    let atoms = structure.positions.len();
    validate_binding(atoms, &structure.properties, topology)?;
    if !is_bound(atoms, &structure.properties, topology) {
        return Err(StructureError::BondCountMismatch {
            expected: topology.bond_count(),
            actual: structure.properties.bonds().len(),
        });
    }
    Ok(())
}

impl<R: ?Sized> Sealed for SampleView<'_, R> {}

impl<R: ?Sized> AsModelView for SampleView<'_, R> {
    fn as_model_view(&self) -> ModelView<'_> {
        SampleView::as_model_view(*self)
    }
}

/// One state of the distribution of a [`SampleSet`] or another
/// [`SampleRead`], with its members, representative, and population.
///
/// Reads the state's energy and properties through `Deref`. The type
/// parameter says where the state comes from, as for [`SampleView`]: a
/// state of a set in memory, `StateView<'_>`, lends its representative
/// sample ([`Self::representative_sample`]); the state of a set in a file,
/// or of another reader, gives its index ([`Self::representative`]), which
/// the host reads.
pub struct StateView<'a, R: ?Sized = SampleSet> {
    reader: &'a dyn SampleRead,
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
    pub(crate) fn new(reader: &'a dyn SampleRead, index: usize) -> Self {
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

    /// The sample index of the state's representative.
    pub fn representative(self) -> usize {
        let distribution = self.distribution();
        distribution.member(distribution.partition().representatives()[self.index])
    }

    /// The sample indices of the state's members, in sample order.
    pub fn samples(self) -> impl Iterator<Item = usize> + 'a {
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
    /// The state's representative sample.
    pub fn representative_sample(self) -> SampleView<'a> {
        let index = self.representative();
        let view = self
            .reader
            .peek(index)
            .expect("a set in memory lends every sample");
        SampleView {
            reader: view.reader,
            sample: view.sample,
            index,
            origin: PhantomData,
        }
    }
}

/// Iterator over the topology-bound samples of a [`SampleSet`], in set
/// order.
#[derive(Debug, Clone)]
pub struct SampleIter<'a> {
    set: &'a SampleSet,
    indices: std::ops::Range<usize>,
}

impl<'a> SampleIter<'a> {
    pub(crate) fn new(set: &'a SampleSet) -> Self {
        Self {
            set,
            indices: 0..set.len(),
        }
    }
}

impl<'a> Iterator for SampleIter<'a> {
    type Item = SampleView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let set = self.set;
        self.indices.next().and_then(|index| set.get(index))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.indices.size_hint()
    }
}

impl ExactSizeIterator for SampleIter<'_> {}

/// Iterator over the states of a set's distribution, in state order; empty
/// without a distribution. Its type parameter is the set's, as for
/// [`StateView`].
pub struct StateIter<'a, R: ?Sized = SampleSet> {
    reader: &'a dyn SampleRead,
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
    pub(crate) fn new(reader: &'a dyn SampleRead) -> Self {
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
