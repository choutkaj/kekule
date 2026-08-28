//! Dense coordinate/data storage and topology-owning models and ensembles.
//!
//! Primitive dense containers are topology-agnostic numerical and data
//! storage. [`Model`] and [`Ensemble`] own topology context and validate dense
//! state against it. Semantic atom, bond, and hierarchy navigation occurs
//! through topology-owning aggregates rather than detached dense arrays.

mod ensemble;
mod model;
mod positions;

pub use ensemble::*;
pub use model::*;
pub use positions::*;

#[cfg(test)]
mod tests;
