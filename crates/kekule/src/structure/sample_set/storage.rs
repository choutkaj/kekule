//! Where a [`super::SampleSet`] keeps its samples: [`InMemory`], or a file.

use std::fmt;

use super::Sample;

pub(crate) mod sealed {
    /// Restricts [`super::SampleStorage`] to Kekule's storages.
    pub trait Storage {
        /// The number of samples stored.
        fn sample_count(&self) -> usize;
    }
}

/// Where a [`SampleSet`](super::SampleSet) keeps its samples, its type
/// parameter: [`InMemory`], the default, or an indexed trajectory file
/// ([`IndexedFile`](crate::traj::io::IndexedFile)).
///
/// Every storage shares the topology, the set properties, the order, the
/// distribution, and the operations on them. Only an in-memory set lends
/// its samples by reference; any set is read through
/// [`SampleRead`](super::SampleRead). The trait is sealed.
pub trait SampleStorage: sealed::Storage + fmt::Debug + Send + Sync {}

/// Samples held in memory, the default storage of a
/// [`SampleSet`](super::SampleSet): they are lent by reference, and reading
/// them cannot fail.
#[derive(Debug, Clone, Default)]
pub struct InMemory {
    pub(super) samples: Vec<Sample>,
}

impl InMemory {
    pub(super) const fn new(samples: Vec<Sample>) -> Self {
        Self { samples }
    }
}

impl sealed::Storage for InMemory {
    fn sample_count(&self) -> usize {
        self.samples.len()
    }
}

impl SampleStorage for InMemory {}
