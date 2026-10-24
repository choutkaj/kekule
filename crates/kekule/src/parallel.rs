//! Host-supplied parallel-for for Kekule's parallel analyses.
//!
//! An analysis that can spread over threads splits its work into indexed
//! items and runs them through an [`Executor`]. It combines the items'
//! results in a fixed order, so the result never depends on the executor, on
//! how many threads it has, or on where, when, or in which order the items
//! run. [`ScopedThreads`], scoped standard-library threads, is the default;
//! a host with its own thread pool implements [`Executor`] to run analyses
//! on it.
//!
//! [`dssp::assign_with`](crate::dssp::assign_with),
//! [`Ensemble::pairwise_fitted_rmsd_with`](crate::structure::Ensemble::pairwise_fitted_rmsd_with)
//! and
//! [`Ensemble::cluster_by_radius_with`](crate::structure::Ensemble::cluster_by_radius_with)
//! take an executor, and selections measure `within` on one given by
//! [`SelectEvaluateOptions::with_executor`](crate::select::SelectEvaluateOptions::with_executor).
//!
//! ```
//! use std::num::NonZeroUsize;
//! use std::sync::atomic::{AtomicUsize, Ordering};
//! use kekule::parallel::{Executor, ScopedThreads};
//!
//! let executor = ScopedThreads::new(NonZeroUsize::new(4).unwrap());
//! let sum = AtomicUsize::new(0);
//! executor.run(10, &|index| {
//!     sum.fetch_add(index, Ordering::Relaxed);
//! });
//! assert_eq!(sum.into_inner(), 45);
//! ```

use std::num::NonZeroUsize;
use std::ops::Range;
use std::sync::{Mutex, PoisonError};
use std::thread;

/// Runs the independent work items of one parallel analysis.
///
/// An analysis splits each of its stages into at most [`Self::parallelism`]
/// work items and combines their results in order, so the result never
/// depends on where, when, or in which order the items run.
///
/// # Contract
///
/// [`Self::run`] calls `body` exactly once for each index in `0..count`, on
/// any threads and in any order, concurrently or not, and returns only after
/// every call has returned. When a call of `body` panics, `run` panics too,
/// after or instead of the remaining calls; it never returns normally. An
/// executor that skips or repeats an index makes the analysis panic.
pub trait Executor: Sync {
    /// How many work items are worth running at once; an analysis splits a
    /// stage into at most this many.
    fn parallelism(&self) -> NonZeroUsize;

    /// Runs `body(index)` for every index in `0..count`; see the
    /// [contract](Self#contract).
    fn run(&self, count: usize, body: &(dyn Fn(usize) + Sync));
}

/// The default [`Executor`]: scoped standard-library threads, started for
/// each stage and joined before it ends.
///
/// The calling thread takes part, so one thread runs everything on the
/// calling thread and spawns none. An item whose thread cannot be spawned,
/// as on targets without threads, runs on the calling thread. A panic in a
/// work item is resumed on the calling thread once the other threads have
/// finished. Starting threads costs a fraction of a millisecond per stage,
/// which a host's own pool avoids.
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

impl Executor for ScopedThreads {
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
///
/// # Panics
///
/// When a work item panics, or when `executor` skips or repeats one.
pub(crate) fn map_ranges<R: Send>(
    executor: &dyn Executor,
    len: usize,
    min_chunk: usize,
    work: impl Fn(Range<usize>) -> R + Sync,
) -> Vec<R> {
    let chunks = chunk_count(executor, len, min_chunk);
    run_chunks(executor, chunks, |chunk| {
        work(chunk_range(chunk, chunks, len))
    })
}

/// [`map_ranges`] over `0..weights.len()`, with ranges of near-equal total
/// weight instead of near-equal length, and no least length.
pub(crate) fn map_weighted_ranges<R: Send>(
    executor: &dyn Executor,
    weights: &[usize],
    work: impl Fn(Range<usize>) -> R + Sync,
) -> Vec<R> {
    let chunks = chunk_count(executor, weights.len(), 1);
    let mut ends = Vec::with_capacity(chunks + 1);
    ends.push(0);
    let total = weights.iter().sum::<usize>();
    let mut sum = 0;
    for (index, weight) in weights.iter().enumerate() {
        // Range `chunk` ends once the weight so far reaches its share.
        while ends.len() < chunks && sum * chunks >= total * ends.len() && sum > 0 {
            ends.push(index);
        }
        sum += weight;
    }
    ends.resize(chunks + 1, weights.len());
    let ranges = ends
        .windows(2)
        .map(|end| end[0]..end[1])
        .collect::<Vec<_>>();
    run_chunks(executor, ranges.len(), |chunk| work(ranges[chunk].clone()))
}

/// How many ranges [`map_ranges`] splits `0..len` into.
fn chunk_count(executor: &dyn Executor, len: usize, min_chunk: usize) -> usize {
    executor
        .parallelism()
        .get()
        .min(len / min_chunk.max(1))
        .max(1)
}

/// Range `chunk` of `chunks` near-equal contiguous ranges of `0..len`.
fn chunk_range(chunk: usize, chunks: usize, len: usize) -> Range<usize> {
    chunk * len / chunks..(chunk + 1) * len / chunks
}

/// Runs `work` on each of `chunks` work items and returns the results in
/// item order; a single item runs on the calling thread.
fn run_chunks<R: Send>(
    executor: &dyn Executor,
    chunks: usize,
    work: impl Fn(usize) -> R + Sync,
) -> Vec<R> {
    if chunks == 1 {
        return vec![work(0)];
    }
    let slots = (0..chunks).map(|_| Mutex::new(None)).collect::<Vec<_>>();
    executor.run(chunks, &|chunk| {
        let result = work(chunk);
        let previous = slots[chunk]
            .lock()
            .unwrap_or_else(PoisonError::into_inner)
            .replace(result);
        assert!(previous.is_none(), "executor ran work item {chunk} twice");
    });
    slots
        .into_iter()
        .enumerate()
        .map(|(chunk, slot)| {
            slot.into_inner()
                .unwrap_or_else(PoisonError::into_inner)
                .unwrap_or_else(|| panic!("executor skipped work item {chunk}"))
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

    #[test]
    fn weighted_ranges_cover_every_index_once_and_share_the_weight() {
        let skewed = (0..60)
            .map(|index| if index < 10 { 50 } else { 1 })
            .collect::<Vec<_>>();
        for threads in 1..=9 {
            let executor = ScopedThreads::new(NonZeroUsize::new(threads).unwrap());
            for weights in [vec![], vec![0; 5], vec![7], vec![1; 100], skewed.clone()] {
                let ranges = map_weighted_ranges(&executor, &weights, |range| range);
                assert!(ranges.len() <= threads.max(1));
                let indices = ranges.iter().cloned().flatten().collect::<Vec<_>>();
                assert_eq!(indices, (0..weights.len()).collect::<Vec<_>>());
            }
            let ranges = map_weighted_ranges(&executor, &skewed, |range| range);
            let total = skewed.iter().sum::<usize>();
            let heaviest = skewed.iter().max().unwrap();
            for range in ranges {
                let weight = skewed[range].iter().sum::<usize>();
                assert!(
                    weight <= total / threads + heaviest,
                    "{threads} threads: {weight}"
                );
            }
        }
    }
}
