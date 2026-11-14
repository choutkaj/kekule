//! The run control of a long operation's `*_with` form, and its progress
//! counts.

use std::fmt;
use std::num::NonZeroUsize;

use crate::diagnostics::{Cancellation, Cancelled};

use super::{Executor, ScopedThreads};

/// The calling thread alone, the executor of a run without one.
static CALLING_THREAD: ScopedThreads = ScopedThreads::new(NonZeroUsize::MIN);

/// Where a long operation runs, and how a host stops and follows it; see
/// [run control](crate::parallel#run-control).
///
/// [`Run::default`] and [`Run::new`] run on the calling thread, without
/// cancellation or progress, as the operation's plain form does.
///
/// ```
/// use std::num::NonZeroUsize;
/// use std::sync::atomic::{AtomicU64, Ordering};
/// use kekule::diagnostics::Cancellation;
/// use kekule::parallel::{Progress, Run, ScopedThreads};
///
/// let threads = ScopedThreads::new(NonZeroUsize::new(4).unwrap());
/// let cancellation = Cancellation::new(|| false);
/// let completed = AtomicU64::new(0);
/// let progress = |progress: Progress| completed.store(progress.completed(), Ordering::Relaxed);
/// let run = Run::new()
///     .with_executor(&threads)
///     .with_cancellation(&cancellation)
///     .with_progress(&progress);
/// assert_eq!(run.cancellation(), Some(&cancellation));
/// assert_eq!(run.executor().parallelism().get(), 4);
/// ```
#[derive(Clone, Copy, Default)]
pub struct Run<'a> {
    executor: Option<&'a dyn Executor>,
    cancellation: Option<&'a Cancellation>,
    progress: Option<&'a (dyn Fn(Progress) + Sync)>,
}

impl fmt::Debug for Run<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("Run")
            .field("parallelism", &self.executor().parallelism())
            .field("cancellation", &self.cancellation)
            .field("progress", &self.progress.map(|_| "callback"))
            .finish()
    }
}

impl<'a> Run<'a> {
    /// A run on the calling thread, without cancellation or progress.
    pub const fn new() -> Self {
        Self {
            executor: None,
            cancellation: None,
            progress: None,
        }
    }

    /// Sets [`Self::executor`].
    #[must_use]
    pub const fn with_executor(mut self, executor: &'a dyn Executor) -> Self {
        self.executor = Some(executor);
        self
    }

    /// Sets [`Self::cancellation`].
    #[must_use]
    pub const fn with_cancellation(mut self, cancellation: &'a Cancellation) -> Self {
        self.cancellation = Some(cancellation);
        self
    }

    /// Sets the callback that receives the operation's [`Progress`]. None by
    /// default.
    #[must_use]
    pub const fn with_progress(mut self, progress: &'a (dyn Fn(Progress) + Sync)) -> Self {
        self.progress = Some(progress);
        self
    }

    /// The executor that runs the operation's work items; the calling thread
    /// alone by default.
    pub fn executor(&self) -> &'a dyn Executor {
        match self.executor {
            Some(executor) => executor,
            None => &CALLING_THREAD,
        }
    }

    /// The cancellation the operation polls, if any.
    pub const fn cancellation(&self) -> Option<&'a Cancellation> {
        self.cancellation
    }

    /// Fails once [`Self::cancellation`] has been requested, running its
    /// check once; for long operations outside this crate.
    ///
    /// # Errors
    ///
    /// Returns [`Cancelled`] when the cancellation's check returns `true`.
    pub fn check(&self) -> Result<(), Cancelled> {
        Cancellation::poll(self.cancellation)
    }

    /// Passes `progress` to the progress callback, if any; for long
    /// operations outside this crate, which keep the calls one at a time and
    /// their counts nondecreasing.
    pub fn report(&self, progress: Progress) {
        if let Some(callback) = self.progress {
            callback(progress);
        }
    }

    /// Whether a progress callback is set, so counting progress is worth
    /// its cost.
    pub(crate) const fn reports_progress(&self) -> bool {
        self.progress.is_some()
    }
}

/// How far one long operation has come: [`Self::completed`] units of
/// [`Self::total`], in the unit its `*_with` form names, such as fitted
/// pairs or scanned bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Progress {
    completed: u64,
    total: u64,
}

impl Progress {
    /// `completed` units of `total`; a `total` below `completed` is raised
    /// to it.
    pub const fn new(completed: u64, total: u64) -> Self {
        Self {
            completed,
            total: if total < completed { completed } else { total },
        }
    }

    /// The progress of `completed` units of `total`, counted as `usize`.
    pub(crate) fn of(completed: usize, total: usize) -> Self {
        let wide = |count: usize| u64::try_from(count).unwrap_or(u64::MAX);
        Self::new(wide(completed), wide(total))
    }

    /// Units done so far.
    pub const fn completed(&self) -> u64 {
        self.completed
    }

    /// Units in all, never below [`Self::completed`]. An operation whose
    /// work is not known in advance, such as an iteration that may converge
    /// early, states a bound it may end below.
    pub const fn total(&self) -> u64 {
        self.total
    }
}
