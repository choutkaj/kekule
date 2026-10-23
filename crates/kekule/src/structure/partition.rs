use std::fmt;

/// A grouping of every sample of one ensemble into states, with one
/// representative sample per state.
///
/// States are numbered from zero, and every sample is in exactly one of
/// them. A state always contains its representative, so none is empty. The
/// partition records sample indices only: it applies to any ensemble with
/// the same sample count, and [`Ensemble::regroup`](super::Ensemble::regroup)
/// checks that count.
///
/// Every [`Ensemble`](super::Ensemble) owns one partition, its states.
/// Clustering returns a new one without changing the ensemble, and
/// [`one_state`](Self::one_state) and [`state_per_sample`](Self::state_per_sample)
/// say what a source's samples are: draws of one distribution, such as the
/// models of an NMR bundle or the frames of a trajectory, or distinct
/// conformers.
///
/// ```
/// use kekule::structure::EnsemblePartition;
/// let partition = EnsemblePartition::new(vec![0, 1, 0, 1], vec![2, 1])?;
/// assert_eq!(partition.state_count(), 2);
/// assert_eq!(partition.samples(0).collect::<Vec<_>>(), [0, 2]);
/// assert_eq!(EnsemblePartition::one_state(3).assignments(), [0, 0, 0]);
/// assert_eq!(EnsemblePartition::state_per_sample(3).representatives(), [0, 1, 2]);
/// # Ok::<(), kekule::structure::PartitionError>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnsemblePartition {
    assignments: Vec<usize>,
    representatives: Vec<usize>,
}

impl EnsemblePartition {
    /// Validates a partition built elsewhere: `assignments[s]` is the state
    /// of sample `s`, and `representatives[t]` is a sample assigned to state
    /// `t`.
    ///
    /// # Errors
    ///
    /// [`PartitionError::StateOutOfRange`] when a sample's state has no
    /// representative, [`PartitionError::RepresentativeOutOfRange`] when a
    /// representative is not a sample, and
    /// [`PartitionError::RepresentativeNotInState`] when a representative is
    /// assigned to another state.
    pub fn new(
        assignments: Vec<usize>,
        representatives: Vec<usize>,
    ) -> Result<Self, PartitionError> {
        for (sample, state) in assignments.iter().enumerate() {
            if *state >= representatives.len() {
                return Err(PartitionError::StateOutOfRange {
                    sample,
                    state: *state,
                    states: representatives.len(),
                });
            }
        }
        for (state, sample) in representatives.iter().enumerate() {
            match assignments.get(*sample) {
                None => {
                    return Err(PartitionError::RepresentativeOutOfRange {
                        state,
                        sample: *sample,
                        samples: assignments.len(),
                    })
                }
                Some(assigned) if *assigned != state => {
                    return Err(PartitionError::RepresentativeNotInState {
                        state,
                        sample: *sample,
                    })
                }
                Some(_) => {}
            }
        }
        Ok(Self {
            assignments,
            representatives,
        })
    }

    /// One state holding all `samples` samples, represented by the first:
    /// samples drawn from one distribution, such as the models of an NMR
    /// bundle or the frames of a trajectory. No samples give no states.
    pub fn one_state(samples: usize) -> Self {
        Self {
            assignments: vec![0; samples],
            representatives: if samples == 0 { Vec::new() } else { vec![0] },
        }
    }

    /// One state per sample, each represented by its only sample: distinct
    /// conformations, such as the conformers of a conformer search.
    pub fn state_per_sample(samples: usize) -> Self {
        Self {
            assignments: (0..samples).collect(),
            representatives: (0..samples).collect(),
        }
    }

    /// A partition whose invariants the caller has established.
    pub(crate) fn from_parts_unchecked(
        assignments: Vec<usize>,
        representatives: Vec<usize>,
    ) -> Self {
        debug_assert!(Self::new(assignments.clone(), representatives.clone()).is_ok());
        Self {
            assignments,
            representatives,
        }
    }

    /// The assignments and representatives, to be rebuilt.
    pub(crate) fn into_parts(self) -> (Vec<usize>, Vec<usize>) {
        (self.assignments, self.representatives)
    }

    /// The number of partitioned samples.
    pub fn sample_count(&self) -> usize {
        self.assignments.len()
    }

    /// The number of states, one per representative.
    pub fn state_count(&self) -> usize {
        self.representatives.len()
    }

    /// The state of every sample, in sample order.
    pub fn assignments(&self) -> &[usize] {
        &self.assignments
    }

    /// The representative sample of every state, in state order.
    pub fn representatives(&self) -> &[usize] {
        &self.representatives
    }

    /// The samples of one state in sample order; empty for a state that
    /// does not exist.
    pub fn samples(&self, state: usize) -> impl Iterator<Item = usize> + '_ {
        self.assignments
            .iter()
            .enumerate()
            .filter(move |(_, assigned)| **assigned == state)
            .map(|(sample, _)| sample)
    }

    /// The number of samples of every state, in state order.
    pub(crate) fn sizes(&self) -> Vec<usize> {
        let mut sizes = vec![0; self.state_count()];
        for state in &self.assignments {
            sizes[*state] += 1;
        }
        sizes
    }
}

/// An invalid [`EnsemblePartition`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum PartitionError {
    /// A sample is assigned to a state without a representative.
    StateOutOfRange {
        /// The assigned sample.
        sample: usize,
        /// The state it is assigned to.
        state: usize,
        /// The number of states, one per representative.
        states: usize,
    },
    /// A representative is not a partitioned sample.
    RepresentativeOutOfRange {
        /// The represented state.
        state: usize,
        /// The representative sample index.
        sample: usize,
        /// The number of partitioned samples.
        samples: usize,
    },
    /// A representative is assigned to another state.
    RepresentativeNotInState {
        /// The represented state.
        state: usize,
        /// The representative sample index.
        sample: usize,
    },
}

impl fmt::Display for PartitionError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::StateOutOfRange {
                sample,
                state,
                states,
            } => write!(
                formatter,
                "sample {sample} is assigned to state {state}, but there are {states} states"
            ),
            Self::RepresentativeOutOfRange {
                state,
                sample,
                samples,
            } => write!(
                formatter,
                "state {state} is represented by sample {sample}, but there are {samples} samples"
            ),
            Self::RepresentativeNotInState { state, sample } => write!(
                formatter,
                "sample {sample} represents state {state} but is assigned to another"
            ),
        }
    }
}

impl std::error::Error for PartitionError {}
