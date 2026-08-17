//! Topology-bound coordinate state, models, borrowed views, atom data, and
//! finite non-temporal ensembles.

mod data;
mod ensemble;
mod model;
mod positions;
pub mod remap;

pub use data::*;
pub use ensemble::*;
pub use model::*;
pub use positions::*;
pub use remap::TopologyRemapError;

#[cfg(test)]
mod tests;
