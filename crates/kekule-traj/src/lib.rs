//! Fixed-topology molecular trajectories, streaming I/O, and trajectory
//! analysis built on [`kekule`].
//!
//! # Data model
//!
//! [`Trajectory`] is an ordered temporal sequence sharing one immutable
//! [`kekule::topology::Topology`]. Each [`TrajectoryFrame`] stores dense
//! positions and optional cell, velocity, force, time, step, and property state;
//! it carries no independent topology. Insertion validates every dense
//! dimension against the trajectory topology.
//!
//! [`TrajectoryFrameView`] and [`FrameBuffer`] expose
//! [`kekule::structure::ModelView`] without copying coordinates, so structural
//! analyses and prepared potentials can operate on models, ensemble members,
//! and trajectory frames through one contract.
//!
//! # In-memory trajectory
//!
//! ```
//! use std::sync::Arc;
//!
//! use kekule::{smiles, structure::Positions, topology::Topology};
//! use kekule_traj::{Trajectory, TrajectoryFrame};
//!
//! let molecule = smiles::to_molecules("CC")?.pop().unwrap();
//! let topology = Arc::new(Topology::from_molecule(&molecule)?);
//! let frame = TrajectoryFrame::new(Positions::zeros(topology.atom_count()));
//!
//! let mut trajectory = Trajectory::new(Arc::clone(&topology));
//! trajectory.push(frame)?;
//! assert_eq!(trajectory.len(), 1);
//! assert_eq!(trajectory.frame(0).unwrap().as_model().atom_count(), 2);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # File I/O and analysis
//!
//! Use [`io::read_trajectory`] for a fully loaded trajectory, or
//! [`io::open_trajectory`] and a reusable [`FrameBuffer`] to process a large
//! file one frame at a time.
//!
//! ```no_run
//! use kekule::{mmcif, topology::AtomSelection};
//! use kekule_traj::io::read_trajectory;
//!
//! let document = mmcif::parse_str(&std::fs::read_to_string("system.cif")?)?;
//! let topology = document.interpret()?.to_topology();
//! let trajectory = read_trajectory("trajectory.xyz", topology.clone())?;
//! println!("{} frames, {} atoms", trajectory.len(), topology.atom_count());
//!
//! let fit = AtomSelection::all(&topology);
//! // Requires a nonempty trajectory and a non-collinear fitting selection.
//! // Coordinates are fitted as stored; repair split molecules first if needed.
//! let aligned = trajectory.superpose_to_frame(0, &fit)?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! Format-agnostic path readers and writers plus pure-Rust XYZ, DCD, TRR, and
//! XTC codecs live in [`io`]. Readers take a topology directly and interpret file
//! coordinates in its dense atom order. They check counts and available format
//! metadata automatically; matching counts alone cannot establish atom identity.
//! In-memory superposition and direct or
//! aligned RMSD workflows live in [`analysis`]. Direct RMSD never performs an
//! implicit fit. Coordinate transformations return a new trajectory by default;
//! explicit `_in_place` methods mutate transactionally. Superposition reports
//! are available through [`Trajectory::superpose_to_frame_with_report`].
//! Molecular reconstruction, imaging, and temporal unwrapping live in [`periodic`]
//! and are explicit preprocessing steps, independent of alignment.
#![forbid(unsafe_code)]
#![warn(rustdoc::broken_intra_doc_links)]
// Kekule consistently names owned conversions `to_*`, including consuming ones.
#![allow(clippy::wrong_self_convention)]

mod trajectory;

pub mod analysis;
pub mod io;
pub mod periodic;

pub use trajectory::*;
