//! Trajectory file codecs, streaming reader and writer contracts, reusable
//! frame buffers, periodic preprocessing, and streaming analysis.
//!
//! # Data structure
//!
//! A loaded trajectory is a [`crate::structure::StructureSet`] with a time
//! [`crate::structure::Order`]: its [`crate::structure::Geometry`]s share one
//! immutable topology, and each holds positions, an optional cell, structure
//! properties, and optional velocities and forces. Each frame's time and MD
//! step belong to the order's [`crate::structure::TimeAxis`]. A geometry
//! carries no topology of its own.
//!
//! This module adds what needs files or bounded memory: pure-Rust XYZ, DCD,
//! TRR and XTC codecs in [`io`]; the [`TrajectoryReader`],
//! [`SeekableTrajectoryReader`] and [`TrajectoryWriter`] contracts, with
//! [`MemoryTrajectoryWriter`] collecting frames into a time-ordered geometry
//! set; the reusable, topology-bound [`FrameBuffer`], filled
//! transactionally from [`FrameData`]; and [`TrajectoryIoError`], the error
//! of every I/O operation, classified by [`TrajectoryIoErrorKind`].
//!
//! A writer takes a [`FrameView`]: a borrowed geometry with its time and step.
//! [`FrameView::from_geometry`] lends one of a set and [`FrameView::new`] one
//! held apart from a set, neither copying arrays.
//! A buffer reads its geometry through `Deref`, lends a
//! [`crate::structure::StructureView`] through [`FrameBuffer::as_structure_view`],
//! and becomes a geometry with its place through [`FrameBuffer::into_geometry`].
//! Analyses and prepared potentials therefore operate on structures, geometries,
//! and buffered frames through one contract.
//!
//! # In-memory trajectory
//!
//! ```
//! use std::sync::Arc;
//!
//! use kekule::structure::{Order, Place, Positions, Geometry, StructureSet, TimeAxis};
//! use kekule::units::{Quantity, PICOSECOND};
//! use kekule::{smiles, topology::Topology};
//!
//! let molecule = smiles::to_molecules("CC")?.pop().unwrap();
//! let topology = Arc::new(Topology::from_molecule(molecule)?);
//!
//! let mut trajectory = StructureSet::new(Arc::clone(&topology));
//! let times = Quantity::new(Vec::new(), PICOSECOND);
//! trajectory.set_order(Some(Order::Time(TimeAxis::new(Some(times), None)?)))?;
//! let place = Place::Time { time: Some(Quantity::new(0.0, PICOSECOND)), step: None };
//! trajectory.push(Geometry::new(Positions::zeros(topology.atom_count())), Some(place))?;
//! assert_eq!(trajectory.len(), 1);
//! assert_eq!(trajectory.get(0).unwrap().as_structure_view().atom_count(), 2);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # File I/O and analysis
//!
//! Use [`io::read_structure_set`] and [`io::write_structure_set`] for loaded
//! trajectories, or [`io::open_trajectory`] and a reusable [`FrameBuffer`]
//! to process a large file one frame at a time.
//!
//! This example loads a topology from an mmCIF file, reads a trajectory, aligns
//! all frames to the first, and saves every tenth frame. The mmCIF must contain
//! one structural block and structure and match the trajectory's atom order.
//!
//! ```no_run
//! use kekule::traj::{
//!     io::{read_structure_set, write_structure_set},
//!     periodic,
//! };
//! use kekule::{mmcif, topology::AtomSelection, units::ANGSTROM};
//!
//! let document = mmcif::parse_str(&std::fs::read_to_string("system.cif")?, Default::default())?;
//! let topology = document.interpret(Default::default())?.into_topology();
//! let mut trajectory = read_structure_set("trajectory.xtc", topology.clone(), Default::default())?;
//! println!("{} frames, {} atoms", trajectory.len(), topology.atom_count());
//!
//! // Reconstruct molecules split across periodic boundaries first.
//! periodic::make_molecules_whole(&mut trajectory)?;
//!
//! // Fit all atoms onto frame 0, or a protein or backbone selection for a
//! // solvated system. Superposition changes the set in place.
//! let fit = AtomSelection::all(topology.clone());
//! trajectory.superpose(0, &fit, Default::default())?;
//! let rmsd = trajectory.rmsd(0, &fit, Default::default())?.value_in(ANGSTROM)?;
//! for (index, value) in rmsd.iter().enumerate() {
//!     println!("Frame {index}: fitted RMSD = {value:.3} A");
//! }
//!
//! // Selection keeps each frame's time and simulation step.
//! let sampled = trajectory.select((0..trajectory.len()).step_by(10))?;
//! write_structure_set("aligned.xtc", &sampled, Default::default())?;
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! Alignment requires a nonempty set and a non-collinear fitting selection.
//! Superposition and periodic operations change the set in place and leave
//! it unchanged if they fail; clone it first to keep the original. Imaging
//! and temporal unwrapping are explicit preprocessing steps; unwrap before
//! discarding intermediate frames.
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
//! of structure sets (see [`crate::alignment`]). They mutate in place and
//! transactionally; clone first to keep the original. Direct RMSD never
//! performs an implicit fit. An independently loaded reference pairs its
//! atoms with a [`crate::alignment::AtomCorrespondence`]; correspondence does
//! not assert chemical equality.
//!
//! [`analysis`] adds the streaming [`analysis::FrameSuperposer`]. Weighted
//! RMSF and contact occupancy are [`crate::statistics`]; their accumulators
//! observe streamed frames through `buffer.as_structure_view()` with explicit
//! weights. Molecular reconstruction, imaging around anchors or onto a
//! reference frame, and temporal unwrapping live in [`periodic`] and are
//! explicit preprocessing steps, independent of alignment.
//! [`analysis::FrameSuperposer`], [`periodic::MoleculeImager`], and
//! [`periodic::TrajectoryUnwrapper`] apply the in-memory operations to
//! streamed frames.
//!
//! # Structure sets of files
//!
//! [`io::open_indexed_structure_set`] opens a file as a
//! [`crate::structure::StructureSet`] of [`io::IndexedFile`] storage, which
//! decodes each geometry when it is read through
//! [`crate::structure::StructureSetRead`]. [`Processed`] runs
//! [`ProcessingSteps`] on the geometries of any reader as they are read, with
//! the kernels of [`periodic`] and [`analysis::FrameSuperposer`].
//! [`crate::structure::StructureSet::from_reader`] loads chosen geometries of any
//! reader, by range, stride, or explicit indices, into memory.
mod processed;
mod trajectory;

pub mod analysis;
pub mod io;
pub mod periodic;

pub use processed::{FitReference, Processed, ProcessingError, ProcessingStep, ProcessingSteps};
pub use trajectory::{
    FrameBuffer, FrameData, FrameView, MemoryTrajectoryWriter, SeekableTrajectoryReader,
    TrajectoryField, TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind,
    TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
