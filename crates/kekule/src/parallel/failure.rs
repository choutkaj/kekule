//! The first failure among a stage's work items, which stops the items past
//! it.

use std::sync::atomic::{AtomicUsize, Ordering};

/// The earliest failed work index of one stage, so work items past it stop:
/// their results cannot matter once an earlier index has failed.
#[derive(Debug)]
pub(crate) struct FirstFailure(AtomicUsize);

impl FirstFailure {
    pub(crate) const fn new() -> Self {
        Self(AtomicUsize::new(usize::MAX))
    }

    pub(crate) fn record(&self, index: usize) {
        self.0.fetch_min(index, Ordering::Relaxed);
    }

    /// Whether an index before `index` has failed.
    pub(crate) fn precedes(&self, index: usize) -> bool {
        self.0.load(Ordering::Relaxed) < index
    }
}

/// Why a work item stopped before its last index.
#[derive(Debug)]
pub(crate) enum Stop<E> {
    /// An index failed with this error, or cancellation was requested.
    Error(E),
    /// An earlier work item's index failed, so this item's results cannot
    /// matter.
    Superseded,
}

impl<E> From<E> for Stop<E> {
    fn from(error: E) -> Self {
        Self::Error(error)
    }
}

impl<E> Stop<E> {
    /// The error of the first stopped work item in order. A superseded item
    /// always follows the failed one that superseded it.
    pub(crate) fn into_error(self) -> E {
        match self {
            Self::Error(error) => error,
            Self::Superseded => unreachable!("a failed work item precedes every superseded one"),
        }
    }
}
