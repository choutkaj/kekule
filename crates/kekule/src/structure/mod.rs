//! Dense coordinate/data storage and topology-owning models and ensembles.
//!
//! Primitive arrays deliberately expose no topology ownership or generic
//! remapping API:
//!
//! ```compile_fail
//! # use kekule::structure::Positions;
//! # let positions: Positions = todo!();
//! let _ = positions.shared_topology();
//! ```
//!
//! ```compile_fail
//! # use kekule::structure::AtomData;
//! # use kekule::topology::InstanceAtomId;
//! # let data = AtomData::new(1);
//! # let atom: InstanceAtomId = todo!();
//! let _ = data.occupancy(atom);
//! ```
//!
//! Removed generic remapping entry points are also deliberately unavailable:
//!
//! ```compile_fail
//! # use std::sync::Arc;
//! # use kekule::structure::Positions;
//! # use kekule::topology::Topology;
//! # let positions: Positions = todo!();
//! # let topology: Arc<Topology> = todo!();
//! let _ = positions.remap_to(&topology);
//! ```
//!
//! ```compile_fail
//! # use std::sync::Arc;
//! # use kekule::structure::AtomData;
//! # use kekule::topology::Topology;
//! # let data = AtomData::new(1);
//! # let topology: Arc<Topology> = todo!();
//! let _ = data.remap_to(&topology);
//! ```
//!
//! ```compile_fail
//! # use std::sync::Arc;
//! # use kekule::structure::BondData;
//! # use kekule::topology::Topology;
//! # let data = BondData::new(1);
//! # let topology: Arc<Topology> = todo!();
//! let _ = data.remap_to(&topology);
//! ```
//!
//! ```compile_fail
//! # use std::sync::Arc;
//! # use kekule::structure::Model;
//! # use kekule::topology::Topology;
//! # let model: Model = todo!();
//! # let topology: Arc<Topology> = todo!();
//! let _ = model.remap_to(&topology);
//! ```
//!
//! ```compile_fail
//! # use std::sync::Arc;
//! # use kekule::structure::Ensemble;
//! # use kekule::topology::Topology;
//! # let ensemble: Ensemble = todo!();
//! # let topology: Arc<Topology> = todo!();
//! let _ = ensemble.remap_to(&topology);
//! ```

mod data;
mod ensemble;
mod model;
mod positions;

pub use data::*;
pub use ensemble::*;
pub use model::*;
pub use positions::*;

#[cfg(test)]
mod tests;
