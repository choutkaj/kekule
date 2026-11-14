//! Pure-Rust foundations for molecular graphs, chemical perception, structure
//! I/O, structural bioinformatics, and molecular modelling.
//!
//! # Object model
//!
//! Kekule keeps chemical identity, system organization, and geometry separate:
//!
//! - [`molecule::Molecule`] is one non-empty connected molecular graph. It owns
//!   represented chemistry, derived [`perception::Perception`], and
//!   geometry-independent properties.
//! - [`topology::Topology`] is one coordinate-free system containing one or
//!   more explicit molecule instances plus an optional biological hierarchy.
//! - [`structure::Structure`] is a topology plus one geometry of it: positions,
//!   an optional periodic cell, and geometry properties.
//! - [`structure::StructureSet`] stores [`structure::Geometry`]s of one shared
//!   topology, each with positions and optional cell, velocities, and forces.
//!   An optional [`structure::Order`] places them in time, as the frames of a
//!   trajectory, or at steps or coordinate values, as the points of an
//!   optimization or a scan. File codecs and streaming analysis for
//!   trajectories live in [`traj`].
//! - [`structure::Ensemble`] reads a structure set as draws from a
//!   [`structure::Distribution`]: which geometries are samples, what each
//!   weighs, and how they group into states that may carry energies.
//!
//! A disconnected salt, solvent box, or protein-ligand complex is therefore a
//! [`topology::Topology`] containing several connected molecules, not one
//! disconnected [`molecule::Molecule`]. Instances may explicitly share one
//! reusable molecule definition. Geometry and numeric arrays address atoms by
//! dense topology index in the topology's stored layout order; molecule atom
//! IDs and hierarchy IDs are separate semantic identities.
//!
//! # Modules
//!
//! | Domain | Modules |
//! | --- | --- |
//! | Foundations | [`units`], [`math`], [`properties`], [`diagnostics`], [`parallel`] |
//! | Chemistry | [`molecule`], [`stereo`] ([`stereo::cip`], [`stereo::coordinates`], [`stereo::cleanup`]), [`perception`] ([`perception::valence`], [`perception::rings`], [`perception::aromaticity`], [`perception::conjugation`], [`perception::resonance`]), [`hydrogens`], [`descriptors`], [`sequence`] |
//! | Queries | [`query`], [`smarts`], [`substructure`], [`select`] |
//! | Systems and structures | [`topology`] ([`topology::edit`]), [`structure`] ([`structure::measure`], [`structure::placement`]), [`solvation`] |
//! | Analysis | [`alignment`], [`clustering`], [`statistics`], [`dssp`], [`interactions`] |
//! | Formats | [`smiles`], [`molfile`], [`sdf`], [`mmcif`], [`pdb`], [`xyz`], [`macromolecular`], [`components`], [`native`] |
//! | Trajectories | [`traj`] ([`traj::io`], [`traj::periodic`], [`traj::analysis`]) |
//!
//! [`prelude`] re-exports the most used items.
//!
//! # Typical workflow
//!
//! Format APIs deliberately separate parsing from chemical interpretation.
//! Interpretation publishes represented chemistry but does not run perception
//! implicitly. Geometry is supplied only when constructing a structure.
//!
//! ```
//! use kekule::{
//!     math::Point3,
//!     smiles,
//!     structure::{Structure, Positions},
//!     units::{Quantity, ANGSTROM},
//! };
//!
//! let mut molecules = smiles::to_molecules("CCO")?;
//! let mut ethanol = molecules.pop().expect("one connected component");
//! ethanol.perceive()?;
//!
//! let positions = Positions::new(Quantity::new(
//!     vec![
//!         Point3::new(0.0, 0.0, 0.0),
//!         Point3::new(1.5, 0.0, 0.0),
//!         Point3::new(2.8, 0.0, 0.0),
//!     ],
//!     ANGSTROM,
//! ))?;
//! let structure = Structure::from_molecule(ethanol, &positions)?;
//!
//! assert_eq!(structure.topology().instance_count(), 1);
//! assert_eq!(structure.atom_count(), 3);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Explicit operations
//!
//! Kekule avoids hidden chemistry changes. Parsing, interpretation,
//! perception, hydrogen transforms, coordinate stereo materialization, and
//! force-field preparation are separate operations. Structural mutation uses
//! transactional builders or editors so published values retain their
//! invariants.
//!
//! | Task | Interface | Publication |
//! | --- | --- | --- |
//! | Construct or edit one connected molecule | [`molecule::MoleculeEditor`] | `finish()` |
//! | Assemble complete molecules and reusable instances | [`topology::TopologyBuilder`] | `build()` |
//! | Assemble a system with explicit geometry | [`structure::StructureBuilder`] | `build()` |
//! | Edit system atoms, bonds, components, and hierarchy | [`topology::edit::TopologyEditor`] | `finish()` |
//! | Edit a system together with its coordinates | [`structure::StructureEditor`] | `finish()` |
//! | Change geometry or structure annotations | [`structure::Structure`] setters | Immediate checked update |
//!
//! Builders and editors offer non-consuming `validate()` and recoverable
//! `try_build()` / `try_finish()`. Use `to_editor()` for a detached draft or
//! `into_editor()` to move an owner into one. System editors resolve source IDs
//! to stable draft-only editing handles and publish completed values through
//! `finish()`. A bond deletion can split a system molecule;
//! a bond addition can join two occurrences. New structure atoms require
//! coordinates.
//!
//! Coordinate-dependent algorithms consume [`structure::StructureView`]. A
//! structure, a geometry of a structure set, or a streamed trajectory frame can
//! therefore share kernels without copying coordinates. APIs that address atoms
//! by index accept any topology snapshot sharing their layout
//! ([`topology::Topology::shares_layout`]), so perception never detaches
//! selections, buffers, or prepared potentials. Use
//! [`topology::Topology::same_layout`] only when complete static layout equality
//! between independent publications is intended.
#![forbid(unsafe_code)]
#![warn(rustdoc::broken_intra_doc_links)]

#[macro_use]
mod fixed_u32_id;

mod algorithms;
pub mod alignment;
mod chemistry;
pub mod clustering;
pub mod components;
pub mod descriptors;
pub mod diagnostics;
mod draft;
pub mod dssp;
pub mod hydrogens;
pub mod interactions;
mod io;
pub mod macromolecular;
pub mod math;
pub mod mmcif;
pub mod molecule;
pub mod molfile;
pub mod native;
pub mod parallel;
pub mod pdb;
pub mod perception;
pub mod prelude;
pub mod properties;
pub mod query;
pub mod sdf;
pub mod select;
mod selection_language;
pub mod sequence;
pub mod smarts;
pub mod smiles;
pub mod solvation;
pub mod statistics;
pub mod stereo;
pub mod structure;
pub mod substructure;
pub mod topology;
pub mod traj;
pub mod units;
pub mod xyz;

pub use draft::DraftError;

#[cfg(test)]
mod tests;

// Compile and run the README examples so they cannot drift from the API. The
// README is the workspace root's in the repository and the package root's once
// packaged; `CARGO_PKG_README` names it relative to the manifest in both.
#[cfg(doctest)]
#[doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/", env!("CARGO_PKG_README")))]
struct ReadmeExamples;
