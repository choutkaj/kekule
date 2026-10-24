//! Host-supplied parallel-for, and ordered fork-join over index ranges on it.

use std::num::NonZeroUsize;
use std::ops::Range;
use std::sync::{Mutex, PoisonError};
use std::thread;

/// Runs the independent work items of one DSSP assignment, for
/// [`assign_with`](super::assign_with).
///
/// A host with its own thread pool implements this to run assignments on
/// it. Assignment splits each stage into at most [`Self::parallelism`] work
/// items and combines their results in order, so the result never depends
/// on where, when, or in which order the items run.
///
/// # Contract
///
/// [`Self::run`] calls `body` exactly once for each index in `0..count`, on
/// any threads and in any order, concurrently or not, and returns only after
/// every call has returned. When a call of `body` panics, `run` panics too,
/// after or instead of the remaining calls; it never returns normally. An
/// executor that skips or repeats an index makes assignment panic.
pub trait DsspExecutor: Sync {
    /// How many work items are worth running at once; assignment splits a
    /// stage into at most this many.
    fn parallelism(&self) -> NonZeroUsize;

    /// Runs `body(index)` for every index in `0..count`; see the
    /// [contract](Self#contract).
    fn run(&self, count: usize, body: &(dyn Fn(usize) + Sync));
}

/// The default [`DsspExecutor`]: scoped standard-library threads, started
/// for each stage and joined before it ends.
///
/// The calling thread takes part, so one thread runs everything on the
/// calling thread and spawns none. An item whose thread cannot be spawned,
/// as on targets without threads, runs on the calling thread. A panic in a
/// work item is resumed on the calling thread once the other threads have
/// finished.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScopedThreads {
    threads: NonZeroUsize,
}

impl ScopedThreads {
    /// An executor of at most `threads` threads, the calling thread
    /// included.
    pub const fn new(threads: NonZeroUsize) -> Self {
        Self { threads }
    }

    /// Most threads, the calling thread included.
    pub const fn threads(&self) -> NonZeroUsize {
        self.threads
    }
}

impl DsspExecutor for ScopedThreads {
    fn parallelism(&self) -> NonZeroUsize {
        self.threads
    }

    fn run(&self, count: usize, body: &(dyn Fn(usize) + Sync)) {
        let threads = self.threads.get().min(count);
        // Thread `worker` runs items `worker`, `worker + threads`, ...
        let items = move |worker: usize| {
            for index in (worker..count).step_by(threads.max(1)) {
                body(index);
            }
        };
        if threads <= 1 {
            return items(0);
        }
        thread::scope(|scope| {
            let spawned = (1..threads)
                .map(|worker| {
                    let handle = thread::Builder::new()
                        .spawn_scoped(scope, move || items(worker))
                        .ok();
                    (worker, handle)
                })
                .collect::<Vec<_>>();
            items(0);
            for (worker, handle) in spawned {
                match handle {
                    Some(handle) => {
                        if let Err(panic) = handle.join() {
                            std::panic::resume_unwind(panic);
                        }
                    }
                    None => items(worker),
                }
            }
        });
    }
}

/// Splits `0..len` into at most `executor.parallelism()` contiguous ranges
/// of at least `min_chunk` items, runs `work` on each as one work item, and
/// returns the results in range order.
///
/// Too little work for two ranges runs `work(0..len)` on the calling thread
/// without the executor.
pub(super) fn map_ranges<R: Send>(
    executor: &dyn DsspExecutor,
    len: usize,
    min_chunk: usize,
    work: impl Fn(Range<usize>) -> R + Sync,
) -> Vec<R> {
    let chunks = executor
        .parallelism()
        .get()
        .min(len / min_chunk.max(1))
        .max(1);
    if chunks == 1 {
        return vec![work(0..len)];
    }
    let slots = (0..chunks).map(|_| Mutex::new(None)).collect::<Vec<_>>();
    executor.run(chunks, &|chunk| {
        let result = work(chunk * len / chunks..(chunk + 1) * len / chunks);
        let previous = slots[chunk]
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .replace(result);
        assert!(
            previous.is_none(),
            "DSSP executor ran work item {chunk} twice"
        );
    });
    slots
        .into_iter()
        .enumerate()
        .map(|(chunk, slot)| {
            slot.into_inner()
                .unwrap_or_else(PoisonError::into_inner)
                .unwrap_or_else(|| panic!("DSSP executor skipped work item {chunk}"))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ranges_cover_every_index_once_in_order() {
        for threads in 1..=9 {
            let executor = ScopedThreads::new(NonZeroUsize::new(threads).unwrap());
            for len in [0, 1, 7, 64, 1000] {
                let ranges = map_ranges(&executor, len, 4, |range| range);
                assert!(ranges.len() <= threads);
                assert!(ranges
                    .iter()
                    .all(|range| range.len() >= 4 || ranges.len() == 1));
                let indices = ranges.into_iter().flatten().collect::<Vec<_>>();
                assert_eq!(indices, (0..len).collect::<Vec<_>>());
            }
        }
    }
}
