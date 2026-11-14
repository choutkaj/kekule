//! Run control for long operations, and the host-supplied parallel-for they
//! spread their work on.
//!
//! # Run control
//!
//! An operation that can take long has a plain form, which runs on the
//! calling thread, cannot be cancelled, and reports no progress, and one
//! `*_with` form that also takes a [`Run`]: an optional [`Executor`], an
//! optional [`Cancellation`](crate::diagnostics::Cancellation), and an
//! optional [`Progress`] callback. `Run::default()` makes the `*_with` form
//! behave as the plain one. Options never hold executors or thread counts.
//!
//! - **Executor.** An operation that can spread over threads splits its work
//!   into indexed items, runs them on the executor, and combines their
//!   results in a fixed order, so its result, including which error it
//!   reports, never depends on the executor, on how many threads it has, or
//!   on where, when, or in which order the items run. Work too small to
//!   split runs on the calling thread, and an operation that does not split
//!   its work ignores the executor. [`ScopedThreads`], scoped
//!   standard-library threads, is ready to use; a host with its own thread
//!   pool implements [`Executor`]. A panicking work item, or an executor
//!   that breaks the [contract](Executor#contract), makes the operation
//!   panic.
//! - **Cancellation.** The operation polls the cancellation before its first
//!   stage, between stages, and periodically inside long loops, on every
//!   thread doing its work. Once requested, it fails with its error for
//!   [`Cancelled`](crate::diagnostics::Cancelled), returns no partial result,
//!   and leaves its receiver unchanged. A completed run returns exactly what
//!   the plain form returns.
//! - **Progress.** The callback receives a [`Progress`] before the first
//!   unit of work and as units complete, on the threads doing the work, one
//!   call at a time, with counts that never decrease. It should be cheap,
//!   such as storing atomics a progress bar reads. Each `*_with` form names
//!   the unit it counts; one that names none reports no progress.
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

mod executor;
mod run;

pub(crate) use executor::{map_ranges, map_weighted_ranges};
pub use executor::{Executor, ScopedThreads};
pub use run::{Progress, Run};
