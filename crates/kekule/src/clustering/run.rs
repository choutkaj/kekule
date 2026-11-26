//! The bookkeeping of a clustering run's cancellation and progress across
//! work items.

use std::sync::{Mutex, PoisonError};

use crate::parallel::{Executor, Progress, Run};

use super::ClusteringError;

/// Atom pairs fitted between two checkpoints: under a millisecond of work.
const CHECKPOINT_WORK: usize = 1 << 14;

/// Atom pairs fitted by one work item at least, a few checkpoints' worth, so
/// splitting pays for starting threads.
const WORK_ITEM_WORK: usize = 1 << 16;

/// Progress and cancellation shared by the work items of one run.
pub(super) struct Tracker<'r> {
    run: Run<'r>,
    total: usize,
    completed: Mutex<usize>,
    fits_per_checkpoint: usize,
}

impl<'r> Tracker<'r> {
    /// Reports that nothing of `total` is done yet and polls cancellation.
    /// Each fit pairs `pair_count` atoms.
    pub(super) fn start(
        run: Run<'r>,
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
        self.run.executor()
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
        self.run.check().map_err(ClusteringError::Cancelled)
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
        if !self.run.reports_progress() {
            return;
        }
        // Counting and reporting under one lock keeps calls one at a time and
        // their counts nondecreasing.
        let mut done = self
            .completed
            .lock()
            .unwrap_or_else(PoisonError::into_inner);
        *done += completed;
        self.run.report(Progress::of(*done, self.total));
    }
}
