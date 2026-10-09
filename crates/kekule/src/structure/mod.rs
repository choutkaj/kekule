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
//!
//! # Geometry editing and scans
//!
//! [`Model::set_distance`], [`Model::set_angle`], and [`Model::set_dihedral`]
//! move connected fragments. Prepare a [`DihedralEdit`] (or [`DistanceEdit`] /
//! [`AngleEdit`]) once for repeated absolute targets on the same topology:
//! ```
//! use kekule::{smiles, geometry::Point3, structure::{Model, Positions, DihedralEdit},
//!     units::{Quantity, ANGSTROM, DEGREE}};
//! let topology = smiles::to_topology("CCCC")?;
//! let positions = Positions::new(Quantity::new([
//!     Point3::new(0.0, 1.0, 0.0), Point3::origin(),
//!     Point3::new(1.5, 0.0, 0.0), Point3::new(1.5, 1.0, 1.0),
//! ], ANGSTROM))?;
//! let mut model = Model::new(topology, positions)?;
//! let [a, b, c, d] = model.topology().atom_ids().try_into().unwrap();
//! model.set_distance(b, c, Quantity::new(1.6, ANGSTROM))?;
//! let edit = DihedralEdit::new(&model.shared_topology(), [a, b, c, d])?;
//! for degrees in [-180.0, -60.0, 60.0, 180.0] {
//!     let mut sample = model.clone();
//!     edit.apply(&mut sample, Quantity::new(degrees, DEGREE))?;
//!     // Evaluate or collect sample; model remains the scan baseline.
//! }
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! Explicit selections can intentionally deform rings or boundary bonds:
//! ```
//! use kekule::{structure::{AngleEdit, Model}, topology::AtomSelection,
//!     units::{Quantity, DEGREE}};
//! # fn example(model: &mut Model) -> Result<(), Box<dyn std::error::Error>> {
//! let [a, b, c] = model.topology().atom_ids()[..3].try_into().unwrap();
//! let topology = model.shared_topology();
//! let moving = AtomSelection::from_atoms(&topology, [c])?;
//! let edit = AngleEdit::with_moving_atoms(&topology, [a, b, c], &moving)?;
//! edit.apply(model, Quantity::new(120.0, DEGREE))?;
//! // Only C moves; bonds to other unselected neighbors can deform.
//! # Ok(())
//! # }
//! ```
//!
//! Geometry edits are atomic, use stored Cartesian coordinates, and preserve
//! topology, properties, and cell. Apply periodic preprocessing explicitly.
//! [`Model::translate`], [`Model::rotate`], and [`Model::apply_transform`] also
//! support rigid selection manipulation. No relaxation is performed.

mod conformation;
mod dynamics;
mod ensemble;
mod geometry_edit;
pub mod measure;
mod model;
mod model_editor;
mod positions;
mod realizations;
mod solvation;
mod trajectory;
pub use solvation::*;

pub use conformation::*;
pub use dynamics::*;
pub use ensemble::*;
pub use geometry_edit::*;
pub use model::*;
pub use model_editor::*;
pub use positions::*;
pub(crate) use realizations::RealizationStore;
pub use realizations::{
    Realization, RealizationError, RealizationIter, RealizationMut, RealizationView,
};
pub use trajectory::*;

#[cfg(test)]
mod model_tests;
#[cfg(test)]
mod tests;
