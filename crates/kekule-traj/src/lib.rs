//! Streaming trajectory I/O, reusable frame buffers, periodic preprocessing,
//! and trajectory reductions built on [`kekule`].
//!
//! # Data model
//!
//! In-memory trajectories live in [`kekule::structure`]:
//! [`kekule::structure::Trajectory`] is an ordered temporal collection of
//! [`kekule::structure::TrajectoryFrame`] payloads sharing one immutable
//! topology. A frame is a [`kekule::structure::Conformation`] plus optional
//! velocities, forces, time, and step; it carries no topology of its own.
//! Frame order is temporal. [`kekule::structure::Ensemble`] is the distinct
//! type for weighted, unordered samples; both offer superposition and RMSD.
//!
//! This crate adds what needs files or bounded memory: format codecs, the
//! [`TrajectoryReader`] and [`TrajectoryWriter`] contracts, and the reusable,
//! topology-bound [`FrameBuffer`]. A buffer reads frame state through `Deref`
//! and implements [`kekule::structure::AsModelView`], so analyses and prepared
//! potentials operate on models, ensemble members, trajectory frames, and
//! buffered frames through one contract.
//!
//! # In-memory trajectory
//!
//! ```
//! use std::sync::Arc;
//!
//! use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
//! use kekule::{smiles, topology::Topology};
//!
//! let molecule = smiles::to_molecules("CC")?.pop().unwrap();
//! let topology = Arc::new(Topology::from_molecule(molecule)?);
//! let frame = TrajectoryFrame::new(Positions::zeros(topology.atom_count()));
//!
//! let mut trajectory = Trajectory::new(Arc::clone(&topology));
//! trajectory.push(frame)?;
//! assert_eq!(trajectory.len(), 1);
//! assert_eq!(trajectory.get(0).unwrap().as_model_view().atom_count(), 2);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # File I/O and analysis
//!
//! Use [`io::read_trajectory`] and [`io::write_trajectory`] for loaded trajectories, or
//! [`io::open_trajectory`] and a reusable [`FrameBuffer`] to process a large
//! file one frame at a time.
//!
//! ```no_run
//! use kekule::{mmcif, topology::AtomSelection};
//! use kekule_traj::io::{read_trajectory, write_trajectory};
//!
//! let document = mmcif::parse_str(&std::fs::read_to_string("system.cif")?)?;
//! let topology = document.interpret()?.into_topology();
//! let mut trajectory = read_trajectory("trajectory.xyz", topology.clone())?;
//! println!("{} frames, {} atoms", trajectory.len(), topology.atom_count());
//!
//! let fit = AtomSelection::all(&topology);
//! // Requires a nonempty trajectory and a non-collinear fitting selection.
//! // Coordinates are fitted as stored; repair split molecules first if needed.
//! trajectory.superpose(0, &fit)?;
//! let sampled = trajectory.select((0..trajectory.len()).step_by(10))?;
//! write_trajectory("aligned.xyz", &sampled)?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! Format-agnostic path readers and writers plus pure-Rust XYZ, DCD, TRR, and
//! XTC codecs live in [`io`]. Readers take a topology directly and interpret file
//! coordinates in its dense atom order. They check counts and available format
//! metadata automatically; matching counts alone cannot establish atom identity.
//! Interpreted structure files keep their source atom-row order, so a trajectory
//! written in a structure file's atom order lines up with that file's topology.
//!
//! In-memory superposition, direct RMSD, and fit-then-measure RMSD are methods
//! of trajectories and ensembles (see [`kekule::alignment`]). They mutate in
//! place and transactionally; clone first to keep the original. Direct RMSD
//! never performs an implicit fit. An independently loaded reference pairs its
//! atoms with a [`kekule::alignment::AtomCorrespondence`]; correspondence does
//! not assert chemical equality.
//!
//! [`analysis`] adds the streaming [`analysis::FrameSuperposer`] and per-atom
//! reductions. Molecular reconstruction, imaging, and temporal unwrapping live
//! in [`periodic`] and are explicit preprocessing steps, independent of
//! alignment. [`analysis::FrameSuperposer`], [`periodic::MoleculeImager`], and
//! [`periodic::TrajectoryUnwrapper`] apply the in-memory operations to
//! streamed frames.
#![forbid(unsafe_code)]
#![warn(rustdoc::broken_intra_doc_links)]

mod trajectory;

pub mod analysis;
pub mod io;
pub mod periodic;

pub use trajectory::*;
