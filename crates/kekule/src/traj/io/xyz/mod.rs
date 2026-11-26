//! Strict bounded multi-frame XYZ trajectory I/O.
//!
//! [`XyzReader`] and [`IndexedXyzReader`] read constant-count
//! element/x/y/z frames whose atom lines follow the topology's dense order and
//! elements; [`XyzWriter`] streams trajectory frames in the text
//! [`kekule::xyz::write`](crate::xyz::write) writes for a structure or a
//! structure set. XYZ stores positions only: the writer rejects a frame's
//! velocities, forces, time, step, occupancies, B-factors and properties,
//! and its cell, unless its options omit them, and reports everything it
//! leaves out, such as bonds, charges and omitted cells, in
//! [`XyzWriter::report`]. To
//! obtain a topology from an XYZ file itself, interpret one of its frames
//! with `kekule::xyz`, which keeps atom-line order as dense order and infers
//! bonds from distances. These readers read atom lines as `kekule::xyz`
//! parses them: an element symbol in any letter case or an atomic number,
//! three coordinates, and columns after them as
//! [`XyzReadOptions::extra_columns`] allows. Lines end with LF, CRLF, or CR,
//! and blank lines may trail the last frame.
//!
//! The scratch limit bounds the dense position array, one point per atom,
//! and the reusable line buffer, which grows to hold the longest line read so
//! far, terminator included.
//!
//! # Examples
//!
//! ```
//! use std::io::Cursor;
//! use std::sync::Arc;
//!
//! use kekule::math::Point3;
//! use kekule::units::{Quantity, ANGSTROM};
//! use kekule::{smiles, topology::Topology};
//! use kekule::traj::io::xyz::{XyzReadOptions, XyzReader, XyzWriter};
//! use kekule::xyz::XyzWriteOptions;
//! use kekule::traj::{FrameBuffer, TrajectoryReader, TrajectoryWriter};
//!
//! let molecule = smiles::to_molecules("[He]")?.pop().unwrap();
//! let topology = Arc::new(Topology::from_molecule(molecule)?);
//! let mut frame = FrameBuffer::new(Arc::clone(&topology));
//! frame.geometry_mut().set_positions(Quantity::new(vec![Point3::new(1.5, 0.0, 0.0)], ANGSTROM))?;
//!
//! let options = XyzWriteOptions::default().with_decimal_places(3)?.with_comment("helium")?;
//! let mut writer = XyzWriter::new(Vec::new(), Arc::clone(&topology), options)?;
//! writer.write_frame(frame.frame_view())?;
//! let bytes = writer.finish()?;
//! assert_eq!(String::from_utf8(bytes.clone())?, "1\nhelium\nHe 1.500 0.000 0.000\n");
//!
//! let mut reader = XyzReader::new(Cursor::new(bytes), topology, XyzReadOptions::default())?;
//! let mut frame = reader.frame_buffer();
//! assert!(reader.read_next(&mut frame)?);
//! let position = frame.positions().get(0).unwrap().value_in(ANGSTROM)?;
//! assert_eq!(position, Point3::new(1.5, 0.0, 0.0));
//! assert!(!reader.read_next(&mut frame)?);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::diagnostics::{LimitExceeded, Resource, UnsupportedContent, WriteOmission, WriteReport};
use crate::io::xyz::atom_line::{element_of_token, parse_atom_line, AtomLineError};
use crate::traj::trajectory::reason::Reason;
use crate::traj::FrameView;
use std::fmt;
use std::io::{self, BufRead, Seek, SeekFrom, Write};
use std::sync::Arc;

use crate::math::Point3;
use crate::structure::Positions;
use crate::topology::Topology;
use crate::traj::io::DecodedGeometry;
use crate::traj::{
    FrameBuffer, FrameData, TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind,
    TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
use crate::units::{Quantity, Unit, ANGSTROM, CANONICAL_LENGTH_UNIT};

use super::codec::frame_error;
use super::frame_state::{check_frame_state, record};
use super::indexed_reader::{
    indexed_reader_traits, IndexFacts, IndexedCodec, IndexedReader, IndexedSource,
};
use super::limits::{exceeded_bytes, filled_bytes, reserve_scratch};
use super::metadata::{FieldAvailability, TrajectoryFieldAvailability};
use super::options::ReadSource;
use super::TrajectoryWriteOptions;
use super::{
    codec_context, io_context, require_no_partial_frame, require_nonempty_writer, FrameOffsets,
    FrameStamp, IndexBuildControl, TrajectoryIoLimits, TrajectoryOpenOptions,
};

/// The format of this codec's streams.
const FORMAT: TrajectoryFormat = TrajectoryFormat::Xyz;

mod indexed;
mod read;
mod write;

pub use indexed::IndexedXyzReader;
use read::ScanPoll;
pub use read::XyzReader;
pub use write::XyzWriter;

pub(crate) use crate::io::XyzWriteOptions;

/// What an XYZ frame carries: positions alone.
pub(super) const FIELDS: TrajectoryFieldAvailability = TrajectoryFieldAvailability {
    positions: FieldAvailability::Required,
    cell: FieldAvailability::Absent,
    velocities: FieldAvailability::Absent,
    forces: FieldAvailability::Absent,
    time: FieldAvailability::Absent,
    step: FieldAvailability::Absent,
    properties: FieldAvailability::Absent,
};

/// Fails unless `unit` is a length.
fn validate_length_unit(unit: Unit, source_label: &str) -> Result<(), TrajectoryIoError> {
    unit.conversion_factor_to(CANONICAL_LENGTH_UNIT)
        .map(|_| ())
        .map_err(|error| {
            codec_context(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::Open,
                Some(FORMAT),
                source_label,
                Reason::XyzLengthUnitNotLength,
            )
            .with_source(error)
        })
}

/// XYZ length-unit interpretation, and what happens to columns after the
/// coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct XyzReadOptions {
    length_unit: Unit,
    extra_columns: UnsupportedContent,
}

impl Default for XyzReadOptions {
    fn default() -> Self {
        Self {
            length_unit: ANGSTROM,
            extra_columns: UnsupportedContent::Reject,
        }
    }
}

impl XyzReadOptions {
    /// Sets [`Self::length_unit`].
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InconsistentMetadata`] when
    /// `length_unit` is not a length.
    pub fn with_length_unit(mut self, length_unit: Unit) -> Result<Self, TrajectoryIoError> {
        validate_length_unit(length_unit, "XYZ reader options")?;
        self.length_unit = length_unit;
        Ok(self)
    }

    /// The length unit of file coordinates; angstrom by default. Positions are
    /// converted to the canonical length unit on publication.
    pub const fn length_unit(&self) -> Unit {
        self.length_unit
    }

    /// Sets [`Self::extra_columns`].
    #[must_use]
    pub const fn with_extra_columns(mut self, extra_columns: UnsupportedContent) -> Self {
        self.extra_columns = extra_columns;
        self
    }

    /// What happens to columns after the coordinates, such as extended XYZ
    /// per-atom properties, which Kekule does not read.
    /// [`UnsupportedContent::Reject`] by default: such a frame fails to read
    /// with [`TrajectoryIoErrorKind::UnsupportedVariant`].
    /// [`UnsupportedContent::Report`] omits them, and an indexed reader's
    /// open report counts the frames that have them
    /// ([`TrajectoryOpenNote::XyzExtraColumns`]).
    ///
    /// [`TrajectoryOpenNote::XyzExtraColumns`]: super::TrajectoryOpenNote::XyzExtraColumns
    pub const fn extra_columns(&self) -> UnsupportedContent {
        self.extra_columns
    }
}
