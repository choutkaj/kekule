//! Connected molecular graphs and their represented and perceived chemistry.
//!
//! [`Molecule`] is the foundational chemical object. Its [`Graph`] stores
//! authoritative atoms, bonds, connectivity, and represented stereochemistry;
//! [`Perception`] stores reconstructible derived chemistry. Use
//! [`MoleculeEditor`] for construction and structural mutation so the
//! non-empty, connected publication invariant is checked transactionally.
//!
//! Coordinates and residue/chain hierarchy deliberately do not live here. See
//! [`crate::structure`] and [`crate::topology`] for those layers.

mod atom_bond;
mod element;
mod element_reference;
mod graph;
mod ids;
mod molecule;
mod molecule_edit;
mod molecule_editor_ops;
mod molecule_workflows;
mod perception;
mod stereo;

pub use atom_bond::*;
pub use element::*;
pub use element_reference::*;
pub use graph::*;
pub use ids::*;
pub use molecule::*;
pub use molecule_edit::*;
pub use molecule_editor_ops::*;
pub use perception::*;
pub use stereo::*;
