//! Borrowed samples and states of a [`SampleSet`].

use std::sync::Arc;

use crate::structure::model::sealed::Sealed;
use crate::structure::{AsModelView, EnsembleState, Model, ModelView, StructureKind};
use crate::topology::Topology;

use super::{Place, Sample, SampleContext, SampleSet};

/// One topology-bound [`Sample`] borrowed from a [`SampleSet`], with its
/// place in the set's order and its membership in the set's distribution.
///
/// Reads sample state through `Deref` and converts to a [`ModelView`]
/// through [`AsModelView`] without copying; the view's
/// [`StructureKind::Sample`] carries the sample's [`SampleContext`].
#[derive(Debug, Clone, Copy)]
pub struct SampleView<'a> {
    set: &'a SampleSet,
    index: usize,
}

impl std::ops::Deref for SampleView<'_> {
    type Target = Sample;

    fn deref(&self) -> &Self::Target {
        &self.set.samples[self.index]
    }
}

impl<'a> SampleView<'a> {
    pub(crate) const fn new(set: &'a SampleSet, index: usize) -> Self {
        Self { set, index }
    }

    /// The sample's index in its set.
    pub const fn index(self) -> usize {
        self.index
    }

    /// The sample's place in the set's order, or `None` in an unordered set.
    pub fn place(self) -> Option<Place> {
        self.set.order.as_ref().map(|order| order.place(self.index))
    }

    /// The sample's relative weight in the set's distribution (`1.0` when
    /// every member weighs the same), or `None` when it is not a member.
    pub fn weight(self) -> Option<f64> {
        self.set.distribution.as_ref()?.weight(self.index)
    }

    /// The state the sample belongs to, or `None` when it is not a member
    /// of the set's distribution.
    pub fn state(self) -> Option<StateView<'a>> {
        let state = self.set.distribution.as_ref()?.state_of(self.index)?;
        Some(StateView::new(self.set, state))
    }

    /// Whether the sample represents its state.
    pub fn is_representative(self) -> bool {
        self.state()
            .is_some_and(|state| state.representative() == self.index)
    }

    /// The topology the sample is bound to, with the view's full lifetime.
    pub fn topology(self) -> &'a Topology {
        &self.set.topology
    }

    /// A shared handle to the bound topology.
    pub fn shared_topology(self) -> Arc<Topology> {
        self.set.shared_topology()
    }

    /// The borrowed sample with the view's full lifetime.
    pub fn payload(self) -> &'a Sample {
        &self.set.samples[self.index]
    }

    /// Borrows the sample as the common kernel input, with the view's full
    /// lifetime.
    pub fn as_model_view(self) -> ModelView<'a> {
        self.payload().structure().model_view(
            &self.set.topology,
            StructureKind::Sample(SampleContext::of_set(self.set, self.index)),
        )
    }

    /// Materializes the sample as an owned model, dropping its place,
    /// membership, velocities, and forces.
    pub fn to_model(self) -> Model {
        self.as_model_view().to_model()
    }
}

impl Sealed for SampleView<'_> {}

impl AsModelView for SampleView<'_> {
    fn as_model_view(&self) -> ModelView<'_> {
        SampleView::as_model_view(*self)
    }
}

/// One state of a [`SampleSet`]'s distribution, with its members,
/// representative, and population.
///
/// Reads the state's energy and properties through `Deref`.
#[derive(Debug, Clone, Copy)]
pub struct StateView<'a> {
    set: &'a SampleSet,
    index: usize,
}

impl std::ops::Deref for StateView<'_> {
    type Target = EnsembleState;

    fn deref(&self) -> &Self::Target {
        self.payload()
    }
}

impl<'a> StateView<'a> {
    pub(crate) const fn new(set: &'a SampleSet, index: usize) -> Self {
        Self { set, index }
    }

    fn distribution(self) -> &'a super::Distribution {
        self.set
            .distribution
            .as_ref()
            .expect("a state view's set has a distribution")
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

    /// The state's representative sample.
    pub fn representative_sample(self) -> SampleView<'a> {
        SampleView::new(self.set, self.representative())
    }

    /// The sample indices of the state's members, in sample order.
    pub fn samples(self) -> impl Iterator<Item = usize> + 'a {
        let distribution = self.distribution();
        distribution
            .partition()
            .samples(self.index)
            .map(move |position| distribution.member(position))
    }

    /// The population: the total weight of the state's members, relative
    /// like the weights. [`super::Distribution::state_probabilities`] gives
    /// every state's share of the total.
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
            indices: 0..set.samples.len(),
        }
    }
}

impl<'a> Iterator for SampleIter<'a> {
    type Item = SampleView<'a>;

    fn next(&mut self) -> Option<Self::Item> {
        let set = self.set;
        self.indices.next().map(|index| SampleView::new(set, index))
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
            .distribution
            .as_ref()
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
