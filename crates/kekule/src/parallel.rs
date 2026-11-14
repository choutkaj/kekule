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
//! [`StructureSet::pairwise_fitted_rmsd_with`](crate::structure::StructureSet::pairwise_fitted_rmsd_with)
//! and
//! [`StructureSet::cluster_by_radius_with`](crate::structure::StructureSet::cluster_by_radius_with)
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

mod executor;

pub(crate) use executor::{map_ranges, map_weighted_ranges};
pub use executor::{Executor, ScopedThreads};
