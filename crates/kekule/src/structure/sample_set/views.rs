//! Borrowed samples and states of a [`SampleSet`] or another
//! [`SampleRead`].

use std::fmt;
use std::sync::Arc;

use crate::structure::model::sealed::Sealed;
use crate::structure::state::{is_bound, validate_binding};
use crate::structure::{
    AsModelView, EnsembleState, Model, ModelView, StructureError, StructureKind,
};
use crate::topology::Topology;

use super::{Place, Sample, SampleContext, SampleRead, SampleReadError, SampleSet};

/// One topology-bound [`Sample`] borrowed from a [`SampleSet`] or another
/// [`SampleRead`], with its place in the reader's order and its membership
/// in the reader's distribution.
///
/// Reads sample state through `Deref` and converts to a [`ModelView`]
/// through [`AsModelView`] without copying; the view's
/// [`StructureKind::Sample`] carries the sample's [`SampleContext`].
#[derive(Clone, Copy)]
pub struct SampleView<'a> {
    pub(super) reader: &'a dyn SampleRead,
    pub(super) sample: &'a Sample,
    pub(super) index: usize,
}

impl fmt::Debug for SampleView<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SampleView")
            .field("index", &self.index)
            .field("sample", self.sample)
            .finish_non_exhaustive()
    }
}

impl std::ops::Deref for SampleView<'_> {
    type Target = Sample;

    fn deref(&self) -> &Self::Target {
        self.sample
    }
}

impl<'a> SampleView<'a> {
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
        }
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
    pub fn state(self) -> Option<StateView<'a>> {
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

    /// The borrowed sample with the view's full lifetime.
    pub fn payload(self) -> &'a Sample {
        self.sample
    }

    /// Borrows the sample as the common kernel input, with the view's full
    /// lifetime.
    pub fn as_model_view(self) -> ModelView<'a> {
        self.sample.structure().model_view(
            self.reader.topology(),
            StructureKind::Sample(SampleContext::of_view(self)),
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

impl Sealed for SampleView<'_> {}

impl AsModelView for SampleView<'_> {
    fn as_model_view(&self) -> ModelView<'_> {
        SampleView::as_model_view(*self)
    }
}

/// One state of the distribution of a [`SampleSet`] or another
/// [`SampleRead`], with its members, representative, and population.
///
/// Reads the state's energy and properties through `Deref`.
#[derive(Clone, Copy)]
pub struct StateView<'a> {
    reader: &'a dyn SampleRead,
    index: usize,
}

impl fmt::Debug for StateView<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("StateView")
            .field("index", &self.index)
            .field("state", self.payload())
            .finish_non_exhaustive()
    }
}

impl std::ops::Deref for StateView<'_> {
    type Target = EnsembleState;

    fn deref(&self) -> &Self::Target {
        self.payload()
    }
}

impl<'a> StateView<'a> {
    pub(crate) fn new(reader: &'a dyn SampleRead, index: usize) -> Self {
        Self { reader, index }
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

    /// The state's representative sample, lent by the reader
    /// ([`SampleRead::peek`]): always by a set in memory, and `None` when
    /// the reader does not hold it.
    pub fn representative_sample(self) -> Option<SampleView<'a>> {
        self.reader.peek(self.representative())
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
    pub fn payload(self) -> &'a EnsembleState {
        self.distribution().state(self.index)
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

/// Iterator over the states of a [`SampleSet`]'s distribution, in state
/// order; empty without a distribution.
#[derive(Debug, Clone)]
pub struct StateIter<'a> {
    set: &'a SampleSet,
    indices: std::ops::Range<usize>,
}

impl<'a> StateIter<'a> {
    pub(crate) fn new(set: &'a SampleSet) -> Self {
        let states = set
            .distribution()
            .map_or(0, super::Distribution::state_count);
        Self {
            set,
            indices: 0..states,
        }
    }
}

impl<'a> Iterator for StateIter<'a> {
    type Item = StateView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let set = self.set;
        self.indices.next().map(|index| StateView::new(set, index))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.indices.size_hint()
    }
}

impl ExactSizeIterator for StateIter<'_> {}
