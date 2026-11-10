//! Trajectory file codecs, streaming reader and writer contracts, reusable
//! frame buffers, periodic preprocessing, and streaming analysis.
//!
//! # Data model
//!
//! In-memory trajectories live in [`crate::structure`]:
//! [`crate::structure::Trajectory`] is an ordered temporal collection of
//! [`crate::structure::TrajectoryFrame`] payloads sharing one immutable
//! topology. A frame holds positions, an optional cell, structure
//! properties, and optional velocities, forces, time, and step; it carries no
//! topology of its own.
//! Frame order is temporal. [`crate::structure::Ensemble`] is the distinct
//! type for weighted, unordered samples; both offer superposition and RMSD.
//!
//! This module adds what needs files or bounded memory: pure-Rust XYZ, DCD,
//! TRR and XTC codecs in [`io`]; the [`TrajectoryReader`],
//! [`SeekableTrajectoryReader`] and [`TrajectoryWriter`] contracts, with
//! [`MemoryTrajectoryWriter`] collecting frames into a
//! [`crate::structure::Trajectory`]; the reusable, topology-bound
//! [`FrameBuffer`], filled transactionally from [`FrameData`]; and
//! [`TrajectoryIoError`], the error of every I/O operation, classified by
//! [`TrajectoryIoErrorKind`]. [`read_ensemble`] and [`read_ensemble_seekable`]
//! sample chosen frames of a reader into an ensemble of one state, as
//! [`crate::structure::Trajectory::to_ensemble`] does for a loaded
//! trajectory. A buffer reads frame state through `Deref`
//! and lends a [`crate::structure::ModelView`] through
//! [`FrameBuffer::as_model_view`], so analyses and prepared potentials operate
//! on models, ensemble samples, trajectory frames, and buffered frames through
//! one contract.
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
//! This example loads a topology from an mmCIF file, reads a trajectory, aligns
//! all frames to the first, and saves every tenth frame. The mmCIF must contain
//! one structural block and model and match the trajectory's atom order.
//!
//! ```no_run
//! use kekule::traj::{
//!     io::{read_trajectory, write_trajectory},
//!     periodic,
//! };
//! use kekule::{mmcif, topology::AtomSelection, units::ANGSTROM};
//!
//! let document = mmcif::parse_str(&std::fs::read_to_string("system.cif")?, Default::default())?;
//! let topology = document.interpret(Default::default())?.into_topology();
//! let mut trajectory = read_trajectory("trajectory.xtc", topology.clone(), Default::default())?;
//! println!("{} frames, {} atoms", trajectory.len(), topology.atom_count());
//!
//! // Reconstruct molecules split across periodic boundaries first.
//! periodic::make_molecules_whole(&mut trajectory)?;
//!
//! // Fit all atoms onto frame 0, or a protein or backbone selection for a
//! // solvated system. Superposition changes the trajectory in place.
//! let fit = AtomSelection::all(topology.clone());
//! trajectory.superpose(0, &fit, Default::default())?;
//! let rmsd = trajectory.rmsd(0, &fit, Default::default())?.value_in(ANGSTROM)?;
//! for (index, value) in rmsd.iter().enumerate() {
//!     println!("Frame {index}: fitted RMSD = {value:.3} A");
//! }
//!
//! // Frame selection keeps the original times and simulation steps.
//! let sampled = trajectory.select((0..trajectory.len()).step_by(10))?;
//! write_trajectory("aligned.xtc", &sampled, Default::default())?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! Alignment requires a nonempty trajectory and a non-collinear fitting
//! selection. Superposition and periodic operations change the trajectory in
//! place and leave it unchanged if they fail; clone it first to keep the
//! original. Imaging and temporal unwrapping are explicit preprocessing steps;
//! unwrap before discarding intermediate frames.
//!
//! Saving infers the format from the extension, protects existing files, and
//! rejects unsupported metadata. [`io::TrajectoryOpenOptions`] and
//! [`io::TrajectoryWriteOptions`] set explicit codec policies and precision.
//! The `trajectory_workflow` and `trajectory_streaming` examples of the
//! repository show a loaded and a streamed workflow.
//!
//! Format-agnostic path readers and writers plus pure-Rust XYZ, DCD, TRR, and
//! XTC codecs live in [`io`]. Readers take a topology directly and interpret file
//! coordinates in its dense atom order. They check counts and available format
//! metadata automatically; matching counts alone cannot establish atom identity.
//! Interpreted structure files keep their source atom-row order, so a trajectory
//! written in a structure file's atom order lines up with that file's topology.
//!
//! In-memory superposition, direct RMSD, and fit-then-measure RMSD are methods
//! of trajectories and ensembles (see [`crate::alignment`]). They mutate in
//! place and transactionally; clone first to keep the original. Direct RMSD
//! never performs an implicit fit. An independently loaded reference pairs its
//! atoms with a [`crate::alignment::AtomCorrespondence`]; correspondence does
//! not assert chemical equality.
//!
//! [`analysis`] adds the streaming [`analysis::FrameSuperposer`]. Weighted
//! RMSF and contact occupancy are [`crate::statistics`]; their accumulators
//! observe streamed frames through `buffer.as_model_view()` with explicit
//! weights. Molecular reconstruction, imaging around anchors or onto a
//! reference frame, and temporal unwrapping live in [`periodic`] and are
//! explicit preprocessing steps, independent of alignment.
//! [`analysis::FrameSuperposer`], [`periodic::MoleculeImager`], and
//! [`periodic::TrajectoryUnwrapper`] apply the in-memory operations to
//! streamed frames.
mod trajectory;

pub mod analysis;
pub mod io;
pub mod periodic;

pub use trajectory::{
    read_ensemble, read_ensemble_seekable, EnsembleReadOptions, EnsembleReadProgress, FrameBuffer,
    FrameData, MemoryTrajectoryWriter, SeekableTrajectoryReader, TrajectoryField, TrajectoryFormat,
    TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation, TrajectoryReader,
    TrajectoryWriter,
};
