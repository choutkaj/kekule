//! Borrowed states of an [`Ensemble`].

use std::fmt;

use crate::structure::{GeometryView, InMemory, State};

use super::Ensemble;

/// One state of an [`Ensemble`]'s distribution, with its samples,
/// representative, and population.
///
/// Reads the state's energy and properties through `Deref`. The type
/// parameter is the ensemble's storage: a state of an ensemble in memory,
/// `StateView<'_>`, lends its representative geometry
/// ([`Self::representative_sample`]); the state of an ensemble over a file
/// gives its index ([`Self::representative`]), which the host reads.
pub struct StateView<'a, S = InMemory> {
    ensemble: &'a Ensemble<S>,
    index: usize,
}

impl<S> Clone for StateView<'_, S> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<S> Copy for StateView<'_, S> {}

impl<S> fmt::Debug for StateView<'_, S> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("StateView")
            .field("index", &self.index)
            .field("state", self.payload())
            .finish_non_exhaustive()
    }
}

impl<S> std::ops::Deref for StateView<'_, S> {
    type Target = State;

    fn deref(&self) -> &Self::Target {
        self.payload()
    }
}

impl<'a, S> StateView<'a, S> {
    /// State `index`, which exists, of `ensemble`.
    pub(super) fn new(ensemble: &'a Ensemble<S>, index: usize) -> Self {
        Self { ensemble, index }
    }

    /// The state's index in the distribution.
    pub const fn index(self) -> usize {
        self.index
    }

    /// The geometry index of the state's representative.
    pub fn representative(self) -> usize {
        let distribution = &self.ensemble.distribution;
        distribution.sample(distribution.partition().representatives()[self.index])
    }

    /// The geometry indices of the state's samples, in geometry order.
    pub fn geometries(self) -> impl Iterator<Item = usize> + 'a {
        let distribution = &self.ensemble.distribution;
        distribution
            .partition()
            .samples(self.index)
            .map(move |position| distribution.sample(position))
    }

    /// The natural logarithm of the state's population, the total weight of
    /// its samples, on the scale of the log weights. Only its differences
    /// to other states' carry meaning;
    /// [`Distribution::state_probabilities`](crate::structure::Distribution::state_probabilities)
    /// gives every state's share of the total.
    pub fn log_population(self) -> f64 {
        self.ensemble.distribution.state_log_population(self.index)
    }

    /// The borrowed state with the view's full lifetime.
    pub fn payload(self) -> &'a State {
        self.ensemble.distribution.state(self.index)
    }
}

impl<'a> StateView<'a> {
    /// The state's representative geometry.
    pub fn representative_sample(self) -> GeometryView<'a> {
        self.ensemble
            .set
            .get(self.representative())
            .expect("a representative is a geometry of the set")
    }
}

/// Iterator over the states of an [`Ensemble`]'s distribution, in state
/// order. Its type parameter is the ensemble's storage, as for
/// [`StateView`].
pub struct StateIter<'a, S = InMemory> {
    ensemble: &'a Ensemble<S>,
    indices: std::ops::Range<usize>,
}

impl<S> Clone for StateIter<'_, S> {
    fn clone(&self) -> Self {
        Self {
            ensemble: self.ensemble,
            indices: self.indices.clone(),
        }
    }
}

impl<S> fmt::Debug for StateIter<'_, S> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("StateIter")
            .field("indices", &self.indices)
            .finish_non_exhaustive()
    }
}

impl<'a, S> StateIter<'a, S> {
    /// Every state of `ensemble`.
    pub(super) fn new(ensemble: &'a Ensemble<S>) -> Self {
        Self {
            ensemble,
            indices: 0..ensemble.distribution.state_count(),
        }
    }
}

impl<'a, S> Iterator for StateIter<'a, S> {
    type Item = StateView<'a, S>;

    fn next(&mut self) -> Option<Self::Item> {
        let ensemble = self.ensemble;
        self.indices
            .next()
            .map(|index| StateView::new(ensemble, index))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.indices.size_hint()
    }
}

impl<S> ExactSizeIterator for StateIter<'_, S> {}
