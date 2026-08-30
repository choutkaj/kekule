//! Geometry-bearing realizations of coordinate-free topologies.
//!
//! [`Positions`] is a topology-agnostic dense coordinate array. [`Model`]
//! combines one immutable [`crate::topology::Topology`] with one complete set
//! of positions, an optional periodic cell, and realization-scoped properties.
//! [`Ensemble`] stores several non-temporal realizations of one shared topology.
//!
//! Dense arrays intentionally carry no atom identity. Topology-owning values
//! validate their dimensions and translate semantic atom or bond identifiers to
//! dense indexes. Coordinate-dependent algorithms accept [`ModelView`], which
//! lets models, ensemble members, and companion-crate trajectory frames share
//! kernels without copying coordinates.

mod ensemble;
mod model;
mod positions;

pub use ensemble::*;
pub use model::*;
pub use positions::*;

#[cfg(test)]
mod tests;
