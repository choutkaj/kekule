//! Executor, cancellation and progress of a clustering run, and their
//! bookkeeping across work items.

use std::fmt;
use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Mutex, PoisonError};

use crate::diagnostics::Cancellation;
use crate::parallel::{Executor, ScopedThreads};

use super::ClusteringError;

/// The calling thread alone, the executor of a default run.
static CALLING_THREAD: ScopedThreads = ScopedThreads::new(NonZeroUsize::MIN);

/// Atom pairs fitted between two checkpoints: under a millisecond of work.
const CHECKPOINT_WORK: usize = 1 << 14;

/// Atom pairs fitted by one work item at least, a few checkpoints' worth, so
/// splitting pays for starting threads.
const WORK_ITEM_WORK: usize = 1 << 16;

/// Where clustering fits run, and how a host follows and stops them, for
/// [`Ensemble::pairwise_fitted_rmsd_with`](crate::structure::Ensemble::pairwise_fitted_rmsd_with)
/// and
/// [`Ensemble::cluster_by_radius_with`](crate::structure::Ensemble::cluster_by_radius_with).
///
/// By default fits run on the calling thread, without cancellation or
/// progress. The executor decides only where and when fits run: the result,
/// including any error, is the same for every executor and thread count.
///
/// ```
/// use std::num::NonZeroUsize;
/// use std::sync::atomic::{AtomicUsize, Ordering};
/// use kekule::clustering::{ClusteringProgress, ClusteringRun};
/// use kekule::diagnostics::Cancellation;
/// use kekule::parallel::ScopedThreads;
///
/// let threads = ScopedThreads::new(NonZeroUsize::new(4).unwrap());
/// let cancellation = Cancellation::new(|| false);
/// let completed = AtomicUsize::new(0);
/// let progress = |progress: ClusteringProgress| {
///     completed.store(progress.completed(), Ordering::Relaxed);
/// };
/// let run = ClusteringRun::default()
///     .with_executor(&threads)
///     .with_cancellation(&cancellation)
///     .with_progress(&progress);
/// assert_eq!(run.cancellation(), Some(&cancellation));
/// ```
#[derive(Clone, Copy)]
pub struct ClusteringRun<'a> {
    executor: &'a dyn Executor,
    cancellation: Option<&'a Cancellation>,
    progress: Option<&'a (dyn Fn(ClusteringProgress) + Sync)>,
}

impl Default for ClusteringRun<'_> {
    fn default() -> Self {
        Self {
            executor: &CALLING_THREAD,
            cancellation: None,
            progress: None,
        }
    }
}

impl fmt::Debug for ClusteringRun<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("ClusteringRun")
            .field("parallelism", &self.executor.parallelism())
            .field("cancellation", &self.cancellation)
            .field("progress", &self.progress.map(|_| "callback"))
            .finish()
    }
}

impl<'a> ClusteringRun<'a> {
    /// The executor that runs the fits; by default the calling thread
    /// alone.
    pub fn executor(&self) -> &'a dyn Executor {
        self.executor
    }

    /// Sets [`Self::executor`].
    ///
    /// Fits split into at most [`Executor::parallelism`] contiguous ranges,
    /// each one work item, whose results combine in order. Too few fits for
    /// two ranges of a few milliseconds each run on the calling thread
    /// without the executor.
    #[must_use]
    pub fn with_executor(mut self, executor: &'a dyn Executor) -> Self {
        self.executor = executor;
        self
    }

    /// The cancellation polled before the first fit and about once per
    /// millisecond of fits on every thread; once requested, the run fails
    /// with [`ClusteringError::Cancelled`] and returns no partial result.
    /// None by default.
    pub fn cancellation(&self) -> Option<&'a Cancellation> {
        self.cancellation
    }

    /// Sets [`Self::cancellation`].
    #[must_use]
    pub fn with_cancellation(mut self, cancellation: &'a Cancellation) -> Self {
        self.cancellation = Some(cancellation);
        self
    }

    /// Sets the callback called with a [`ClusteringProgress`] before the
    /// first fit, about once per millisecond of fits, and once all are done.
    ///
    /// It runs on the threads doing the fits, one call at a time, with
    /// counts that never decrease. It should be cheap, such as storing
    /// atomics a progress bar reads. None by default.
    #[must_use]
    pub fn with_progress(mut self, progress: &'a (dyn Fn(ClusteringProgress) + Sync)) -> Self {
        self.progress = Some(progress);
        self
    }
}

/// Progress of a clustering run: pairs fitted for
/// [`Ensemble::pairwise_fitted_rmsd_with`](crate::structure::Ensemble::pairwise_fitted_rmsd_with),
/// members placed for
/// [`Ensemble::cluster_by_radius_with`](crate::structure::Ensemble::cluster_by_radius_with).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ClusteringProgress {
    completed: usize,
    total: usize,
}

impl ClusteringProgress {
    /// Pairs fitted or members placed so far.
    pub const fn completed(&self) -> usize {
        self.completed
    }

    /// Pairs or members in all, which [`Self::completed`] reaches when the
    /// run succeeds.
    pub const fn total(&self) -> usize {
        self.total
    }
}

/// Progress and cancellation shared by the work items of one run.
pub(super) struct Tracker<'r> {
    run: ClusteringRun<'r>,
    total: usize,
    completed: Mutex<usize>,
    fits_per_checkpoint: usize,
}

impl<'r> Tracker<'r> {
    /// Reports that nothing of `total` is done yet and polls cancellation.
    /// Each fit pairs `pair_count` atoms.
    pub(super) fn start(
        run: ClusteringRun<'r>,
        total: usize,
        pair_count: usize,
    ) -> Result<Self, ClusteringError> {
        let tracker = Self {
            run,
            total,
            completed: Mutex::new(0),
            fits_per_checkpoint: (CHECKPOINT_WORK / pair_count.max(1)).max(1),
        };
        tracker.report(0);
        tracker.poll()?;
        Ok(tracker)
    }

    pub(super) fn executor(&self) -> &'r dyn Executor {
        self.run.executor
    }

    /// Fits between two checkpoints.
    pub(super) fn fits_per_checkpoint(&self) -> usize {
        self.fits_per_checkpoint
    }

    /// Fits a work item makes at least.
    pub(super) fn fits_per_work_item(&self) -> usize {
        self.fits_per_checkpoint * (WORK_ITEM_WORK / CHECKPOINT_WORK)
    }

    /// Fails once cancellation has been requested.
    pub(super) fn poll(&self) -> Result<(), ClusteringError> {
        Cancellation::poll(self.run.cancellation).map_err(ClusteringError::Cancelled)
    }

    /// Counts `completed` more units done, reports the new count, and polls
    /// cancellation.
    pub(super) fn checkpoint(&self, completed: usize) -> Result<(), ClusteringError> {
        self.finish(completed);
        self.poll()
    }

    /// Counts `completed` more units done and reports the new count; for the
    /// end of a work item, where polling would only discard finished work.
    pub(super) fn finish(&self, completed: usize) {
        if completed > 0 {
            self.report(completed);
        }
    }

    fn report(&self, completed: usize) {
        let Some(progress) = self.run.progress else {
            return;
        };
        // Counting and reporting under one lock keeps calls one at a time and
        // their counts nondecreasing.
        let mut done = self
            .completed
            .lock()
            .unwrap_or_else(PoisonError::into_inner);
        *done += completed;
        progress(ClusteringProgress {
            completed: *done,
            total: self.total,
        });
    }
}

/// The earliest failed work index of one stage, so work items past it stop:
/// their results cannot matter once an earlier index has failed.
pub(super) struct FirstFailure(AtomicUsize);

impl FirstFailure {
    pub(super) fn new() -> Self {
        Self(AtomicUsize::new(usize::MAX))
    }

    pub(super) fn record(&self, index: usize) {
        self.0.fetch_min(index, Ordering::Relaxed);
    }

    /// Whether an index before `index` has failed.
    pub(super) fn precedes(&self, index: usize) -> bool {
        self.0.load(Ordering::Relaxed) < index
    }
}

/// Why a work item stopped before its last fit.
pub(super) enum Stop {
    /// A fit failed or cancellation was requested.
    Error(ClusteringError),
    /// An earlier work item's fit failed, so this item's results cannot
    /// matter.
    Superseded,
}

impl From<ClusteringError> for Stop {
    fn from(error: ClusteringError) -> Self {
        Self::Error(error)
    }
}

impl Stop {
    /// The error of the first stopped work item in order. A superseded item
    /// always follows the failed one that superseded it.
    pub(super) fn into_error(self) -> ClusteringError {
        match self {
            Self::Error(error) => error,
            Self::Superseded => unreachable!("a failed work item precedes every superseded one"),
        }
    }
}
