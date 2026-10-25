//! Bounded pure-Rust file codecs for fixed-topology trajectories.
//!
//! The codecs implement this crate's reusable frame-buffer streaming
//! contracts and depend on [`kekule`] only for topology, structural state,
//! geometry, and units.
//!
//! Start with [`read_trajectory`] to load an entire file into an in-memory
//! [`kekule::structure::Trajectory`], and [`write_trajectory`] to save it
//! atomically. Every operation takes its options, which set the codec and
//! field policies; the defaults detect the format.
//! Use [`open_trajectory`] for sequential processing with
//! a reusable frame buffer, [`open_indexed_trajectory`] for verified random access,
//! [`create_trajectory_writer`] for atomic path-backed writing, or
//! [`create_trajectory_stream_writer`] and [`write_trajectory_to`] for writing
//! into a host's stream, so that the host places the file. Each operation
//! takes a topology directly or from its trajectory. Path readers expose opening
//! diagnostics through `open_report()` and verified format metadata through
//! `metadata()`. The [`xyz`], [`dcd`], [`trr`] and [`xtc`] modules hold the
//! codec-level readers and writers over caller-supplied streams.
//!
//! # Examples
//!
//! ```
//! use std::sync::Arc;
//!
//! use kekule::geometry::Point3;
//! use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
//! use kekule::units::{Quantity, ANGSTROM};
//! use kekule::{smiles, topology::Topology};
//! use kekule_traj::io::{read_trajectory, write_trajectory};
//!
//! let molecule = smiles::to_molecules("[He]")?.pop().unwrap();
//! let topology = Arc::new(Topology::from_molecule(molecule)?);
//! let points = vec![Point3::new(1.5, 0.0, 0.0)];
//! let frame = TrajectoryFrame::new(Positions::from_vec(Quantity::new(points, ANGSTROM))?);
//! let trajectory = Trajectory::from_items(Arc::clone(&topology), [frame])?;
//!
//! let path = std::env::temp_dir().join(format!("kekule-traj-io-{}.xyz", std::process::id()));
//! # let _ = std::fs::remove_file(&path);
//! // The `.xyz` extension selects the codec; an existing file is never replaced.
//! write_trajectory(&path, &trajectory, Default::default())?;
//! let loaded = read_trajectory(&path, Arc::clone(&topology), Default::default());
//! std::fs::remove_file(&path)?;
//! let loaded = loaded?;
//!
//! assert_eq!(loaded.len(), 1);
//! let position = loaded.get(0).unwrap().positions().get(0).unwrap().value_in(ANGSTROM)?;
//! assert_eq!(position, Point3::new(1.5, 0.0, 0.0));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Supported profiles
//!
//! | Format | Reader | Writer |
//! |---|---|---|
//! | XYZ | strict constant-count multi-frame element/x/y/z text; configured length unit (angstrom by default) | deterministic strict text; optional frame state is rejected unless omitted, and omitted frame state and topology content are reported |
//! | DCD | common 32-bit-record `CORD` files in either byte order, common unit cells, fixed-atom reconstruction, and strict `NSET` | canonical all-atom `CORD` in either byte order with optional unit cells; velocities, forces, properties, and times without a `DELTA` policy rejected unless omitted |
//! | TRR | GROMACS XDR frames with f32 or f64 position, box, velocity, and force blocks | one explicit f32 or f64 precision with per-frame optional blocks; properties rejected unless omitted; frames without a time or step rejected unless numbered |
//! | XTC | GROMACS magic 1995/2023, signed nonnegative i32 counts/steps, small uncompressed and ordinary compressed coordinates | magic 1995/2023 at explicit lossy precision, bounding scaled coordinates and per-axis spans as GROMACS does; velocities, forces, and properties rejected unless omitted; frames without a time or step rejected unless numbered |
//!
//! Readers accept an owned topology or a shared `Arc<Topology>`. File coordinate
//! index `i` means topology dense atom index `i`. Counts and available metadata
//! (including XYZ element order) are checked automatically; matching counts
//! alone cannot establish atom identity. Native units are converted once
//! at the codec boundary: DCD and default XYZ use angstrom, while TRR/XTC use
//! GROMACS nanometre/picosecond conventions. XTC coordinate resolution is
//! nominally `1 / precision` nanometres.
//!
//! Sequential readers retain one file handle and avoid a whole-file scan.
//! Codec-level indexed readers, such as [`xtc::IndexedXtcReader::new`], take
//! the stream position at opening as frame 0. Indexed readers retain one handle,
//! fully verify every frame from there during an
//! O(file-size) index build, store bounded checked offsets with capped geometric
//! growth, and then decode one complete frame per random read. Decoding
//! validates into reusable private scratch and publishes transactionally into
//! the caller's [`crate::FrameBuffer`]. Random reads restore all sequential reader
//! state before publication. Clean EOF is accepted only between frames,
//! including through a bounded probe at exact frame/index limits.
//!
//! The path-backed scan reports [`TrajectoryIndexProgress`] and stops on a
//! cancellation set in [`TrajectoryOpenOptions`]. Its offsets export as a
//! portable, versioned [`TrajectoryFrameIndex`] bound to a
//! [`TrajectoryFileFingerprint`] of the file, and
//! [`open_indexed_trajectory_from_index`] reopens the file from it without
//! scanning. Such a reader still verifies every frame it decodes and checks
//! that each random read ends where the index places the next frame.
//!
//! Path writers stage a temporary sibling. Only consuming
//! [`FileTrajectoryWriter::finish`] flushes, synchronizes, finalizes format
//! metadata, and publishes a nonempty trajectory. Any failed frame write or an
//! empty finish prevents publication. Its single rename can fail on Windows
//! while another program briefly holds the destination; a host that retries
//! or places files itself writes through a [`TrajectoryStreamWriter`], whose
//! consuming [`TrajectoryStreamWriter::finish`] finalizes and flushes the
//! host's stream and returns it.
//! The saving helpers reject unsupported collection properties before opening a
//! writer. Format policies also reject unsupported frame fields unless their
//! `frame_metadata` or `frame_properties` options omit them; they never
//! silently discard them, and the writers' reports name what they omit.
//! XTC and TRR number frames without a time or step only under an explicit
//! [`FrameTimingPolicy::Number`], and DCD numbers frames whose steps are missing or off its
//! header's sequence only under [`dcd::DcdFrameTimingPolicy::Number`]. Defaults preserve an existing destination, use TRR f32
//! with its lambda frame property, DCD without cells, and XTC resolution of
//! 0.001 nm. Use explicit write options when those profiles do not fit the data.
//!
//! # Limits and unsupported formats
//!
//! [`TrajectoryIoLimits`] bounds attacker-controlled atoms, frames, records,
//! scratch, index storage, text, comments, and detection before allocation or
//! seeking. Amber NetCDF, PDB, GRO/G96, Amber ASCII, LAMMPS dump, reactive
//! trajectories, and compressed wrappers are outside this initial profile and
//! return structured unsupported-format or unsupported-variant errors.
pub mod dcd;
mod detect;
mod file_reader;
mod file_writer;
mod frame_index;
mod frame_state;
mod limits;
mod metadata;
pub mod trr;
mod util;
pub mod xtc;
pub mod xyz;

pub use file_reader::{
    open_indexed_trajectory, open_indexed_trajectory_from_index, open_trajectory, read_trajectory,
    IndexedFileTrajectoryReader, SequentialFileTrajectoryReader,
};
pub use file_writer::{
    create_trajectory_stream_writer, create_trajectory_writer, write_trajectory,
    write_trajectory_to, FileTrajectoryWriter, OverwritePolicy, TrajectoryStreamWriter,
    TrajectoryWriteOptions,
};
pub use frame_index::{TrajectoryFileFingerprint, TrajectoryFrameIndex, TrajectoryIndexProgress};
pub use frame_state::FrameTimingPolicy;
pub use limits::TrajectoryIoLimits;
pub use metadata::{
    detect_trajectory_format, CoordinateEncoding, FieldAvailability, FileTrajectoryMetadata,
    FormatDetectionEvidence, RandomAccessCapability, ScalarPrecision, TrajectoryFieldAvailability,
    TrajectoryFormatDetection, TrajectoryFormatHint, TrajectoryOpenOptions, TrajectoryOpenReport,
};

pub(crate) use frame_index::{index_mismatch, FrameOffsets, FrameStamp, IndexBuildControl};
pub(crate) use limits::{projected_index_limit, reserve_index_for_push};
pub(crate) use util::{
    codec_context, frame_offset_context, io_context, probe_seekable_eof, require_no_partial_frame,
    require_nonempty_writer, validation_detail,
};

#[cfg(test)]
mod tests;
