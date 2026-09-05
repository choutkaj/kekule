//! Bounded pure-Rust file codecs for fixed-topology trajectories.
//!
//! The codecs implement this crate's reusable frame-buffer streaming
//! contracts and depend on [`kekule`] only for topology, structural state,
//! geometry, and units.
//!
//! Start with [`open_trajectory`] for sequential reading,
//! [`open_indexed_trajectory`] for verified random access, or
//! [`create_trajectory_writer`] for atomic path-backed writing. Each operation
//! takes a topology. Readers expose opening diagnostics through `open_report()`
//! and verified format metadata through `metadata()`.
//!
//! # Supported profiles
//!
//! | Format | Reader | Writer |
//! |---|---|---|
//! | XYZ | strict constant-count multi-frame element/x/y/z text; configured length unit (angstrom by default) | deterministic strict text; optional frame state is rejected |
//! | DCD | common 32-bit-record `CORD` files in either byte order, common unit cells, fixed-atom reconstruction, and strict `NSET` | canonical all-atom `CORD` in either byte order with optional unit cells |
//! | TRR | GROMACS XDR frames with f32 or f64 position, box, velocity, and force blocks | one explicit f32 or f64 precision with per-frame optional blocks |
//! | XTC | GROMACS magic 1995/2023, signed nonnegative i32 counts/steps, small uncompressed and ordinary compressed coordinates | magic 1995/2023 through the private audited encoder adapter at explicit lossy precision |
//!
//! Readers accept an owned topology or a shared `Arc<Topology>`. File coordinate
//! index `i` means topology dense atom index `i`. Counts and available metadata
//! (including XYZ element order) are checked automatically; matching counts
//! alone cannot establish atom identity. [`crate::validate_atom_order`] can check
//! an independently supplied sequence of semantic atom IDs. Native units are converted once
//! at the codec boundary: DCD and default XYZ use angstrom, while TRR/XTC use
//! GROMACS nanometre/picosecond conventions. XTC coordinate resolution is
//! nominally `1 / precision` nanometres.
//!
//! Sequential readers retain one file handle and avoid a whole-file scan.
//! Indexed readers retain one handle, fully verify every frame during an
//! O(file-size) index build, store bounded checked offsets with capped geometric
//! growth, and then decode one complete frame per random read. Decoding
//! validates into reusable private scratch and publishes transactionally into
//! the caller's [`crate::FrameBuffer`]. Random reads restore all sequential reader
//! state before publication. Clean EOF is accepted only between frames,
//! including through a bounded probe at exact frame/index limits.
//!
//! Path writers stage a temporary sibling. Only consuming
//! [`FileTrajectoryWriter::finish`] flushes, synchronizes, finalizes format
//! metadata, and publishes a nonempty trajectory. Any failed frame write or an
//! empty finish prevents publication.
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
mod limits;
mod metadata;
pub mod trr;
mod util;
pub mod xtc;
pub mod xyz;

pub use file_reader::{
    open_indexed_trajectory, open_indexed_trajectory_with_options, open_trajectory,
    open_trajectory_with_options, IndexedFileTrajectoryReader, SequentialFileTrajectoryReader,
};
pub use file_writer::{
    create_trajectory_writer, FileTrajectoryWriter, OverwritePolicy, TrajectoryWriteOptions,
};
pub use limits::TrajectoryIoLimits;
pub use metadata::{
    detect_trajectory_format, CoordinateEncoding, FieldAvailability, FileTrajectoryMetadata,
    FormatDetectionEvidence, RandomAccessCapability, ScalarPrecision, TrajectoryFieldAvailability,
    TrajectoryFormatDetection, TrajectoryFormatHint, TrajectoryOpenOptions, TrajectoryOpenReport,
};

pub(crate) use limits::{projected_index_limit, reserve_index_for_push};
pub(crate) use util::{
    codec_context, frame_offset_context, io_context, probe_seekable_eof, require_nonempty_writer,
};

#[cfg(test)]
mod tests;
