//! Portable verified frame indexes, file fingerprints, and index-build control.
//!
//! An indexed reader's frame offsets are exported as a [`TrajectoryFrameIndex`]
//! bound to a [`TrajectoryFileFingerprint`] of the file. Opening from an index
//! checks the fingerprint, the format, and the first frame's position; every
//! random read then verifies its frame as usual and checks that it ends exactly
//! where the index places the next frame, or the file's end.

use std::fmt;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};
use std::sync::Arc;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

use kekule::diagnostics::{Cancellation, LimitExceeded, Resource};
use kekule::units::{Quantity, CANONICAL_TIME_UNIT};

use crate::{TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation};

use super::limits::index_hard_capacity;
use super::{io_context, reserve_index_for_push, TrajectoryIoLimits};

/// Files up to this size are hashed whole.
const WHOLE_FILE_HASH_BYTES: u64 = 256 * 1024;
/// Bytes hashed at each end of a larger file.
const EDGE_HASH_BYTES: u64 = 64 * 1024;
/// Bytes of each evenly spaced interior sample of a larger file.
const SAMPLE_HASH_BYTES: u64 = 4 * 1024;
/// Interior samples of a larger file.
const SAMPLE_COUNT: u64 = 16;

const MAGIC: [u8; 8] = *b"KEKTRIDX";
const VERSION: u32 = 2;
const HEADER_BYTES: usize = 56;
const CHECKSUM_BYTES: usize = 8;
const FLAG_MIXED_PRECISION: u8 = 1;
const FLAG_STEPS: u8 = 1 << 1;
const FLAG_TIMES: u8 = 1 << 2;

/// The identity of a file's content when its frame index was built.
///
/// It records the file's size, its modification time where the platform
/// reports one, and a 64-bit FNV-1a hash of the size followed by sampled
/// content: the whole file up to 256 KiB; otherwise its first and last
/// 64 KiB and sixteen evenly spaced 4 KiB interior blocks. It detects a
/// replaced, appended, truncated, or rewritten file; it is not a
/// cryptographic digest, and the per-frame checks of an index-backed reader
/// remain the guard against frames that are not where the index says.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct TrajectoryFileFingerprint {
    byte_count: u64,
    /// Nanoseconds from the Unix epoch, negative before it.
    modified: Option<i128>,
    content_hash: u64,
}

impl TrajectoryFileFingerprint {
    /// The file's size in bytes.
    pub const fn byte_count(&self) -> u64 {
        self.byte_count
    }

    /// The file's modification time, when the platform reports one.
    pub fn modified(&self) -> Option<SystemTime> {
        let nanos = self.modified?;
        let magnitude = nanos.unsigned_abs();
        let duration = Duration::new(
            u64::try_from(magnitude / 1_000_000_000).ok()?,
            u32::try_from(magnitude % 1_000_000_000).ok()?,
        );
        if nanos >= 0 {
            UNIX_EPOCH.checked_add(duration)
        } else {
            UNIX_EPOCH.checked_sub(duration)
        }
    }

    /// The FNV-1a hash of the size and the sampled content.
    pub const fn content_hash(&self) -> u64 {
        self.content_hash
    }

    /// Fingerprints `file` and leaves it positioned at its start.
    pub(crate) fn of_file(file: &mut File, source_label: &str) -> Result<Self, TrajectoryIoError> {
        let io_error = |error| io_context(TrajectoryIoOperation::Open, None, source_label, error);
        let metadata = file.metadata().map_err(io_error)?;
        let byte_count = metadata.len();
        let modified = metadata.modified().ok().and_then(unix_nanos);
        let mut hash = Fnv1a::new();
        hash.write(&byte_count.to_le_bytes());
        let mut buffer = [0_u8; 8192];
        for (start, length) in sample_ranges(byte_count) {
            file.seek(SeekFrom::Start(start)).map_err(io_error)?;
            let mut remaining = length;
            while remaining > 0 {
                let chunk = usize::try_from(remaining.min(buffer.len() as u64))
                    .expect("a chunk is no larger than the buffer");
                file.read_exact(&mut buffer[..chunk]).map_err(io_error)?;
                hash.write(&buffer[..chunk]);
                remaining -= chunk as u64;
            }
        }
        file.seek(SeekFrom::Start(0)).map_err(io_error)?;
        Ok(Self {
            byte_count,
            modified,
            content_hash: hash.finish(),
        })
    }
}

fn unix_nanos(time: SystemTime) -> Option<i128> {
    match time.duration_since(UNIX_EPOCH) {
        Ok(after) => i128::try_from(after.as_nanos()).ok(),
        Err(before) => i128::try_from(before.duration().as_nanos())
            .ok()
            .map(|nanos| -nanos),
    }
}

/// The `(start, length)` byte ranges a fingerprint hashes, in order.
fn sample_ranges(byte_count: u64) -> impl Iterator<Item = (u64, u64)> {
    let whole = (byte_count <= WHOLE_FILE_HASH_BYTES).then_some((0, byte_count));
    let sampled = (byte_count > WHOLE_FILE_HASH_BYTES).then(move || {
        let interior = (1..=SAMPLE_COUNT).map(move |sample| {
            let start = u128::from(byte_count) * u128::from(sample) / u128::from(SAMPLE_COUNT + 1);
            (
                u64::try_from(start).expect("a sample starts inside the file"),
                SAMPLE_HASH_BYTES,
            )
        });
        std::iter::once((0, EDGE_HASH_BYTES))
            .chain(interior)
            .chain(std::iter::once((
                byte_count - EDGE_HASH_BYTES,
                EDGE_HASH_BYTES,
            )))
    });
    whole.into_iter().chain(sampled.into_iter().flatten())
}

/// 64-bit FNV-1a, the fingerprint and encoding checksum hash.
struct Fnv1a(u64);

impl Fnv1a {
    const fn new() -> Self {
        Self(0xcbf2_9ce4_8422_2325)
    }

    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.0 = (self.0 ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3);
        }
    }

    const fn finish(&self) -> u64 {
        self.0
    }
}

/// Progress of building a frame index, reported after each verified frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TrajectoryIndexProgress {
    frame_count: usize,
    byte_count: u64,
    total_byte_count: u64,
}

impl TrajectoryIndexProgress {
    /// Frames verified and indexed so far.
    pub const fn frame_count(&self) -> usize {
        self.frame_count
    }

    /// Bytes of the file verified so far, from its start through the last
    /// verified frame, headers included.
    pub const fn byte_count(&self) -> u64 {
        self.byte_count
    }

    /// The file's size in bytes, which [`Self::byte_count`] reaches when the
    /// index is complete.
    pub const fn total_byte_count(&self) -> u64 {
        self.total_byte_count
    }
}

/// A host's observer of index-building progress; clones share one callback.
#[derive(Clone)]
pub(crate) struct IndexProgressCallback(Arc<dyn Fn(TrajectoryIndexProgress) + Send + Sync>);

impl IndexProgressCallback {
    pub(crate) fn new(report: impl Fn(TrajectoryIndexProgress) + Send + Sync + 'static) -> Self {
        Self(Arc::new(report))
    }
}

impl fmt::Debug for IndexProgressCallback {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str("IndexProgressCallback")
    }
}

// As for `Cancellation`: Kekule only calls the callback, so a panic inside it
// cannot leave Kekule state broken, and options holding one stay unwind safe.
impl std::panic::UnwindSafe for IndexProgressCallback {}
impl std::panic::RefUnwindSafe for IndexProgressCallback {}

/// Cancellation and progress reporting for one index build.
#[derive(Clone, Copy, Default)]
pub(crate) struct IndexBuildControl<'a> {
    cancellation: Option<&'a Cancellation>,
    progress: Option<&'a IndexProgressCallback>,
    total_byte_count: u64,
}

impl<'a> IndexBuildControl<'a> {
    pub(crate) const fn new(
        cancellation: Option<&'a Cancellation>,
        progress: Option<&'a IndexProgressCallback>,
        total_byte_count: u64,
    ) -> Self {
        Self {
            cancellation,
            progress,
            total_byte_count,
        }
    }

    /// Reports `frame_count` frames verified through byte `byte_count`, then
    /// polls the cancellation.
    pub(crate) fn checkpoint(
        &self,
        frame_count: usize,
        byte_count: u64,
        format: TrajectoryFormat,
        source_label: &str,
    ) -> Result<(), TrajectoryIoError> {
        if let Some(progress) = self.progress {
            (progress.0)(TrajectoryIndexProgress {
                frame_count,
                byte_count,
                total_byte_count: self.total_byte_count.max(byte_count),
            });
        }
        if let Some(cancellation) = self.cancellation {
            cancellation.check().map_err(|cancelled| {
                TrajectoryIoError::codec(
                    TrajectoryIoErrorKind::Cancelled,
                    TrajectoryIoOperation::Index,
                    Some(format),
                )
                .with_source_label(source_label)
                .with_frame(frame_count)
                .with_source(cancelled)
            })?;
        }
        Ok(())
    }
}

/// The time and step a frame records, where its format has them; the time is
/// in the canonical time unit.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub(crate) struct FrameStamp {
    pub(crate) time: Option<f64>,
    pub(crate) step: Option<u64>,
}

/// Frame start offsets, the end of the last frame, and each frame's time and
/// step where every frame records them, of one indexed stream.
#[derive(Debug, Clone)]
pub(crate) struct FrameOffsets {
    starts: Vec<u64>,
    end: u64,
    /// Finite, in the canonical time unit; as long as `starts` when present.
    times: Option<Vec<f64>>,
    /// As long as `starts` when present.
    steps: Option<Vec<u64>>,
}

/// Times compare by their bits, which are finite.
impl PartialEq for FrameOffsets {
    fn eq(&self, other: &Self) -> bool {
        let bits = |times: &Option<Vec<f64>>| {
            times
                .as_ref()
                .map(|times| times.iter().map(|time| time.to_bits()).collect::<Vec<_>>())
        };
        self.starts == other.starts
            && self.end == other.end
            && self.steps == other.steps
            && bits(&self.times) == bits(&other.times)
    }
}

impl Eq for FrameOffsets {}

impl FrameOffsets {
    /// An empty index of a stream whose frames would start at `start`.
    pub(crate) const fn empty(start: u64) -> Self {
        Self {
            starts: Vec::new(),
            end: start,
            times: None,
            steps: None,
        }
    }

    pub(crate) fn len(&self) -> usize {
        self.starts.len()
    }

    pub(crate) fn is_empty(&self) -> bool {
        self.starts.is_empty()
    }

    pub(crate) fn starts(&self) -> &[u64] {
        &self.starts
    }

    pub(crate) const fn end(&self) -> u64 {
        self.end
    }

    pub(crate) fn times(&self) -> Option<&[f64]> {
        self.times.as_deref()
    }

    pub(crate) fn steps(&self) -> Option<&[u64]> {
        self.steps.as_deref()
    }

    /// Appends a verified frame occupying `start..end` and recording `stamp`
    /// within the limits. Frame 0 decides whether times and steps are
    /// indexed, and every later frame must record the same ones.
    pub(crate) fn push(
        &mut self,
        start: u64,
        end: u64,
        stamp: FrameStamp,
        limits: &TrajectoryIoLimits,
        format: TrajectoryFormat,
        source_label: &str,
    ) -> Result<(), TrajectoryIoError> {
        let frame = self.starts.len();
        if frame == 0 {
            self.times = stamp.time.map(|_| Vec::new());
            self.steps = stamp.step.map(|_| Vec::new());
        } else if stamp.time.is_some() != self.times.is_some()
            || stamp.step.is_some() != self.steps.is_some()
        {
            return Err(TrajectoryIoError::codec(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::Index,
                Some(format),
            )
            .with_source_label(source_label)
            .with_frame(frame)
            .with_detail(format!(
                "{format} frame records a different set of time and step than frame 0"
            )));
        }
        reserve_index_for_push(&mut self.starts, limits, format, source_label, frame)?;
        let capacity = self.starts.capacity();
        if let (Some(times), Some(time)) = (&mut self.times, stamp.time) {
            reserve_stamps(times, capacity, format, source_label, frame)?;
            times.push(time);
        }
        if let (Some(steps), Some(step)) = (&mut self.steps, stamp.step) {
            reserve_stamps(steps, capacity, format, source_label, frame)?;
            steps.push(step);
        }
        self.starts.push(start);
        self.end = end;
        Ok(())
    }

    /// Checks that frame `index`, just decoded, records the time and step
    /// the index holds for it.
    pub(crate) fn verify_stamp(
        &self,
        index: usize,
        stamp: FrameStamp,
        operation: TrajectoryIoOperation,
        format: TrajectoryFormat,
        source_label: &str,
    ) -> Result<(), TrajectoryIoError> {
        let indexed = FrameStamp {
            time: self
                .times
                .as_ref()
                .and_then(|times| times.get(index).copied()),
            step: self
                .steps
                .as_ref()
                .and_then(|steps| steps.get(index).copied()),
        };
        let same_time = match (indexed.time, stamp.time) {
            (Some(indexed), Some(decoded)) => indexed.to_bits() == decoded.to_bits(),
            (None, None) => true,
            _ => false,
        };
        if !same_time || indexed.step != stamp.step {
            return Err(index_mismatch(
                operation,
                format,
                source_label,
                format!("{format} frame time or step differs from the index's"),
            )
            .with_frame(index));
        }
        Ok(())
    }

    /// The start of frame `index` and the offset its verified frame must end
    /// at: the next frame's start, or the end of the last frame.
    pub(crate) fn frame(&self, index: usize) -> Option<(u64, u64)> {
        let start = *self.starts.get(index)?;
        let end = self.starts.get(index + 1).copied().unwrap_or(self.end);
        Some((start, end))
    }

    /// Checks that frame `index`, just decoded from its start, ended at the
    /// current position of `stream`.
    pub(crate) fn verify_stream_frame_end(
        &self,
        index: usize,
        stream: &mut impl Seek,
        format: TrajectoryFormat,
        source_label: &str,
    ) -> Result<(), TrajectoryIoError> {
        let end = stream.stream_position().map_err(|error| {
            io_context(
                TrajectoryIoOperation::ReadFrame,
                Some(format),
                source_label,
                error,
            )
        })?;
        self.verify_frame_end(index, end, format, source_label)
    }

    /// Checks that frame `index`, decoded from its start, ended at `actual_end`.
    fn verify_frame_end(
        &self,
        index: usize,
        actual_end: u64,
        format: TrajectoryFormat,
        source_label: &str,
    ) -> Result<(), TrajectoryIoError> {
        let Some((_, expected_end)) = self.frame(index) else {
            return Err(TrajectoryIoError::frame_index_out_of_range(index));
        };
        if actual_end != expected_end {
            return Err(TrajectoryIoError::codec(
                TrajectoryIoErrorKind::IndexMismatch,
                TrajectoryIoOperation::ReadFrame,
                Some(format),
            )
            .with_source_label(source_label)
            .with_frame(index)
            .with_byte_offset(actual_end)
            .with_counts(expected_end, actual_end)
            .with_detail(format!(
                "{format} frame does not end where the index places the next frame"
            )));
        }
        Ok(())
    }

    /// Checks an imported index against a freshly opened codec reader whose
    /// first frame starts at `stream_start`. Only DCD, whose header declares
    /// the frame count, may index no frame.
    pub(crate) fn verify_start(
        &self,
        stream_start: u64,
        format: TrajectoryFormat,
        source_label: &str,
    ) -> Result<(), TrajectoryIoError> {
        let consistent = match self.starts.first() {
            Some(&first) => first == stream_start,
            None => format == TrajectoryFormat::Dcd && self.end == stream_start,
        };
        if !consistent {
            return Err(index_mismatch(
                TrajectoryIoOperation::Open,
                format,
                source_label,
                format!("{format} index does not start at the file's first frame"),
            ));
        }
        Ok(())
    }
}

/// Grows a full table of times or steps to the offset table's `capacity`.
fn reserve_stamps<T>(
    stamps: &mut Vec<T>,
    capacity: usize,
    format: TrajectoryFormat,
    source_label: &str,
    frame: usize,
) -> Result<(), TrajectoryIoError> {
    if stamps.len() < stamps.capacity() {
        return Ok(());
    }
    stamps
        .try_reserve_exact(capacity.saturating_sub(stamps.len()))
        .map_err(|_| {
            TrajectoryIoError::codec(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::Index,
                Some(format),
            )
            .with_source_label(source_label)
            .with_frame(frame)
            .with_detail(format!(
                "could not grow {format} index times or steps toward its bounded capacity"
            ))
        })
}

pub(crate) fn index_mismatch(
    operation: TrajectoryIoOperation,
    format: TrajectoryFormat,
    source_label: &str,
    detail: impl Into<String>,
) -> TrajectoryIoError {
    TrajectoryIoError::codec(
        TrajectoryIoErrorKind::IndexMismatch,
        operation,
        Some(format),
    )
    .with_source_label(source_label)
    .with_detail(detail)
}

/// A portable, versioned frame index of one trajectory file.
///
/// It holds the format, the byte offset of every frame, the frame count, each
/// frame's time and step where the format records them, and the
/// [`TrajectoryFileFingerprint`] of the file it was built for. Export one
/// from a fully verified reader with
/// [`IndexedFileTrajectoryReader::to_frame_index`], store it with
/// [`Self::to_bytes`], and later open the same file without scanning it
/// again with [`open_indexed_trajectory_from_index`].
///
/// An imported index is untrusted input. [`Self::from_bytes`] bounds its
/// allocation by the [`TrajectoryIoLimits`] and rejects malformed or
/// corrupted encodings through a checksum; opening rejects an index whose
/// fingerprint or format disagree with the file; and every frame read through
/// it is still fully verified and must end exactly where the index places the
/// next frame and record the time and step the index holds for it. A stale
/// or corrupted index is therefore a typed
/// [`TrajectoryIoErrorKind::IndexMismatch`] or
/// [`TrajectoryIoErrorKind::InvalidIndex`] error, never a wrong frame. Only
/// the index records which frame is which, so an index deliberately forged
/// with a valid checksum can at most renumber whole, valid frames of the file;
/// keep indexes where only the host writes them.
///
/// # Encoding
///
/// [`Self::to_bytes`] writes version 2 of a little-endian layout:
///
/// | Bytes | Field |
/// |---|---|
/// | 8 | magic `KEKTRIDX` |
/// | 4 | version, `u32` 2 |
/// | 1 | format: 1 XYZ, 2 DCD, 3 TRR, 4 XTC |
/// | 1 | flags: bit 0 set when TRR frames mix f32 and f64 scalars, bit 1 when steps follow the offsets, bit 2 when times follow them; other bits 0 |
/// | 1 | 1 when the modification time is present, else 0 |
/// | 1 | reserved, 0 |
/// | 8 | file size in bytes, `u64` |
/// | 16 | modification time in nanoseconds from the Unix epoch, `i128`; 0 when absent |
/// | 8 | content hash, `u64` |
/// | 8 | frame count `n`, `u64` |
/// | 8 `n` | strictly increasing frame start offsets, `u64`, each below the file size |
/// | 8 `n` | with flag bit 1, frame steps, `u64` |
/// | 8 `n` | with flag bit 2, frame times in picoseconds, finite `f64` |
/// | 8 | FNV-1a 64 checksum of all preceding bytes, `u64` |
///
/// The last frame ends at the end of the file. Steps and times are present
/// exactly when [`Self::steps`] and [`Self::times`] are: never for XYZ or an
/// empty index; steps always, and times under [`DcdTimePolicy::HeaderDelta`],
/// for DCD; and both for TRR and XTC. Other versions, including version 1,
/// which recorded no times or steps, are rejected with
/// [`TrajectoryIoErrorKind::InvalidIndex`]; build such an index again.
///
/// [`IndexedFileTrajectoryReader::to_frame_index`]:
///     super::IndexedFileTrajectoryReader::to_frame_index
/// [`open_indexed_trajectory_from_index`]: super::open_indexed_trajectory_from_index
/// [`DcdTimePolicy::HeaderDelta`]: super::dcd::DcdTimePolicy::HeaderDelta
///
/// # Examples
///
/// ```no_run
/// use kekule::{smiles, topology::Topology};
/// use kekule_traj::io::{
///     open_indexed_trajectory, open_indexed_trajectory_from_index, TrajectoryFrameIndex,
///     TrajectoryIoLimits,
/// };
/// use kekule_traj::{SeekableTrajectoryReader, TrajectoryReader};
///
/// let molecule = smiles::to_molecules("CC")?.pop().unwrap();
/// let topology = std::sync::Arc::new(Topology::from_molecule(molecule)?);
/// let reader = open_indexed_trajectory("ethane.xtc", topology.clone(), Default::default())?;
/// std::fs::write("ethane.xtc.index", reader.to_frame_index().to_bytes())?;
///
/// // Later: reopen without scanning the file again.
/// let bytes = std::fs::read("ethane.xtc.index")?;
/// let index = TrajectoryFrameIndex::from_bytes(&bytes, &TrajectoryIoLimits::default())?;
/// let mut reader =
///     open_indexed_trajectory_from_index("ethane.xtc", topology, index, Default::default())?;
/// let mut frame = reader.frame_buffer();
/// reader.read_frame(0, &mut frame)?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone, PartialEq, Eq)]
pub struct TrajectoryFrameIndex {
    format: TrajectoryFormat,
    fingerprint: TrajectoryFileFingerprint,
    /// Ends at `fingerprint.byte_count`.
    offsets: FrameOffsets,
    precision_mixed: bool,
}

impl fmt::Debug for TrajectoryFrameIndex {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("TrajectoryFrameIndex")
            .field("format", &self.format)
            .field("frame_count", &self.frame_count())
            .field("fingerprint", &self.fingerprint)
            .finish_non_exhaustive()
    }
}

impl TrajectoryFrameIndex {
    pub(crate) fn new(
        format: TrajectoryFormat,
        fingerprint: TrajectoryFileFingerprint,
        offsets: FrameOffsets,
        precision_mixed: bool,
    ) -> Self {
        debug_assert_eq!(offsets.end(), fingerprint.byte_count);
        Self {
            format,
            fingerprint,
            offsets,
            precision_mixed: precision_mixed && format == TrajectoryFormat::Trr,
        }
    }

    pub(crate) fn into_parts(self) -> (FrameOffsets, bool) {
        (self.offsets, self.precision_mixed)
    }

    /// The format of the indexed file.
    pub const fn format(&self) -> TrajectoryFormat {
        self.format
    }

    /// The number of indexed frames.
    pub fn frame_count(&self) -> usize {
        self.offsets.len()
    }

    /// The byte offset at which each frame starts, in frame order.
    pub fn frame_offsets(&self) -> &[u64] {
        self.offsets.starts()
    }

    /// Each frame's time in the canonical time unit, in frame order, or
    /// `None` when the frames record none.
    ///
    /// TRR and XTC frames record their time. DCD times are the header
    /// `DELTA` times the step, present only when the index was built under
    /// [`DcdTimePolicy::HeaderDelta`], which a reader opened from the index
    /// must then use too. XYZ frames record no time. The times are those the
    /// frames report when read, and need be neither evenly spaced nor
    /// increasing.
    ///
    /// [`DcdTimePolicy::HeaderDelta`]: super::dcd::DcdTimePolicy::HeaderDelta
    pub fn times(&self) -> Option<Quantity<&[f64]>> {
        self.offsets
            .times()
            .map(|times| Quantity::new(times, CANONICAL_TIME_UNIT))
    }

    /// Each frame's integration step, in frame order, or `None` when the
    /// frames record none.
    ///
    /// TRR and XTC frames record their step, DCD steps follow the header's
    /// `ISTART` and `NSAVC`, and XYZ frames record no step.
    pub fn steps(&self) -> Option<&[u64]> {
        self.offsets.steps()
    }

    /// The fingerprint of the file the index was built for.
    pub const fn fingerprint(&self) -> &TrajectoryFileFingerprint {
        &self.fingerprint
    }

    /// Encodes the index in the versioned layout described on the type.
    pub fn to_bytes(&self) -> Vec<u8> {
        let times = self.offsets.times();
        let steps = self.offsets.steps();
        let tables = 1 + usize::from(steps.is_some()) + usize::from(times.is_some());
        let mut bytes = Vec::with_capacity(
            HEADER_BYTES
                + CHECKSUM_BYTES
                + tables * self.offsets.len() * std::mem::size_of::<u64>(),
        );
        bytes.extend_from_slice(&MAGIC);
        bytes.extend_from_slice(&VERSION.to_le_bytes());
        bytes.push(format_code(self.format));
        let mut flags = 0;
        if self.precision_mixed {
            flags |= FLAG_MIXED_PRECISION;
        }
        if steps.is_some() {
            flags |= FLAG_STEPS;
        }
        if times.is_some() {
            flags |= FLAG_TIMES;
        }
        bytes.push(flags);
        bytes.push(u8::from(self.fingerprint.modified.is_some()));
        bytes.push(0);
        bytes.extend_from_slice(&self.fingerprint.byte_count.to_le_bytes());
        bytes.extend_from_slice(&self.fingerprint.modified.unwrap_or(0).to_le_bytes());
        bytes.extend_from_slice(&self.fingerprint.content_hash.to_le_bytes());
        bytes.extend_from_slice(&(self.offsets.len() as u64).to_le_bytes());
        for start in self.offsets.starts() {
            bytes.extend_from_slice(&start.to_le_bytes());
        }
        for step in steps.unwrap_or_default() {
            bytes.extend_from_slice(&step.to_le_bytes());
        }
        for time in times.unwrap_or_default() {
            bytes.extend_from_slice(&time.to_le_bytes());
        }
        let mut checksum = Fnv1a::new();
        checksum.write(&bytes);
        bytes.extend_from_slice(&checksum.finish().to_le_bytes());
        bytes
    }

    /// Decodes an index written by [`Self::to_bytes`].
    ///
    /// The frame count is checked against the frame and index limits of
    /// `limits` before anything is allocated; the offsets must be strictly
    /// increasing and lie inside the fingerprinted file, the times must be
    /// finite, and the times and steps present must be those the format
    /// records.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidIndex`] for a truncated,
    /// corrupt, or malformed encoding or an unsupported version, and
    /// [`TrajectoryIoErrorKind::LimitExceeded`] when the frame count exceeds
    /// the limits.
    pub fn from_bytes(
        bytes: &[u8],
        limits: &TrajectoryIoLimits,
    ) -> Result<Self, TrajectoryIoError> {
        if bytes.len() < HEADER_BYTES + CHECKSUM_BYTES {
            return Err(invalid_index("frame index is truncated"));
        }
        if bytes[..8] != MAGIC {
            return Err(invalid_index("frame index has no KEKTRIDX magic"));
        }
        let version = u32::from_le_bytes(array(&bytes[8..12]));
        if version != VERSION {
            return Err(invalid_index(format!(
                "frame index version {version} is not supported"
            )));
        }
        let (body, checksum) = bytes.split_at(bytes.len() - CHECKSUM_BYTES);
        let mut expected = Fnv1a::new();
        expected.write(body);
        if expected.finish() != u64::from_le_bytes(array(checksum)) {
            return Err(invalid_index("frame index checksum does not match"));
        }
        let format = match bytes[12] {
            1 => TrajectoryFormat::Xyz,
            2 => TrajectoryFormat::Dcd,
            3 => TrajectoryFormat::Trr,
            4 => TrajectoryFormat::Xtc,
            code => {
                return Err(invalid_index(format!(
                    "frame index format {code} is unknown"
                )))
            }
        };
        let flags = bytes[13];
        if flags & !(FLAG_MIXED_PRECISION | FLAG_STEPS | FLAG_TIMES) != 0
            || (flags & FLAG_MIXED_PRECISION != 0 && format != TrajectoryFormat::Trr)
        {
            return Err(invalid_index(format!(
                "frame index flags {flags:#04x} are invalid for {format}"
            )));
        }
        let has_steps = flags & FLAG_STEPS != 0;
        let has_times = flags & FLAG_TIMES != 0;
        let modified = i128::from_le_bytes(array(&bytes[24..40]));
        let modified = match bytes[14] {
            1 => Some(modified),
            0 if modified == 0 => None,
            _ => return Err(invalid_index("frame index modification time is malformed")),
        };
        if bytes[15] != 0 {
            return Err(invalid_index("frame index reserved byte is nonzero"));
        }
        let fingerprint = TrajectoryFileFingerprint {
            byte_count: u64::from_le_bytes(array(&bytes[16..24])),
            modified,
            content_hash: u64::from_le_bytes(array(&bytes[40..48])),
        };
        let declared = u64::from_le_bytes(array(&bytes[48..56]));
        let capacity = index_hard_capacity(limits);
        let frame_count = usize::try_from(declared)
            .ok()
            .filter(|&count| count <= capacity)
            .ok_or_else(|| {
                TrajectoryIoError::codec(
                    TrajectoryIoErrorKind::LimitExceeded,
                    TrajectoryIoOperation::Index,
                    Some(format),
                )
                .with_detail("frame index frame count exceeds the configured index limits")
                .with_limit(LimitExceeded::new(
                    Resource::Frames,
                    capacity,
                    usize::try_from(declared).unwrap_or(usize::MAX),
                ))
            })?;
        let recorded = match format {
            _ if frame_count == 0 => (false, false),
            TrajectoryFormat::Xyz => (false, false),
            TrajectoryFormat::Dcd => (true, has_times),
            TrajectoryFormat::Trr | TrajectoryFormat::Xtc => (true, true),
        };
        if (has_steps, has_times) != recorded {
            return Err(invalid_index(format!(
                "frame index flags {flags:#04x} disagree with the steps and times \
                 {frame_count} {format} frames record"
            )));
        }
        let table_count = 1 + usize::from(has_steps) + usize::from(has_times);
        let tables = &body[HEADER_BYTES..];
        let table_bytes = frame_count * std::mem::size_of::<u64>();
        if table_bytes
            .checked_mul(table_count)
            .is_none_or(|length| length != tables.len())
        {
            return Err(invalid_index(
                "frame index length disagrees with its frame count",
            ));
        }
        let (offset_bytes, stamp_bytes) = tables.split_at(table_bytes);
        let (step_bytes, time_bytes) =
            stamp_bytes.split_at(if has_steps { table_bytes } else { 0 });
        let mut starts = reserve_decoded(frame_count, format)?;
        for chunk in offset_bytes.as_chunks::<8>().0 {
            let start = u64::from_le_bytes(*chunk);
            if starts.last().is_some_and(|&previous| start <= previous)
                || start >= fingerprint.byte_count
            {
                return Err(invalid_index(format!(
                    "frame index offset {} is out of order or outside the file",
                    starts.len()
                )));
            }
            starts.push(start);
        }
        let steps = if has_steps {
            let mut steps = reserve_decoded(frame_count, format)?;
            steps.extend(
                step_bytes
                    .as_chunks::<8>()
                    .0
                    .iter()
                    .map(|chunk| u64::from_le_bytes(*chunk)),
            );
            Some(steps)
        } else {
            None
        };
        let times = if has_times {
            let mut times = reserve_decoded(frame_count, format)?;
            for chunk in time_bytes.as_chunks::<8>().0 {
                let time = f64::from_le_bytes(*chunk);
                if !time.is_finite() {
                    return Err(invalid_index(format!(
                        "frame index time {} is not finite",
                        times.len()
                    )));
                }
                times.push(time);
            }
            Some(times)
        } else {
            None
        };
        Ok(Self {
            format,
            offsets: FrameOffsets {
                starts,
                end: fingerprint.byte_count,
                times,
                steps,
            },
            fingerprint,
            precision_mixed: flags & FLAG_MIXED_PRECISION != 0,
        })
    }
}

/// An empty table with room for `frame_count` decoded entries.
fn reserve_decoded<T>(
    frame_count: usize,
    format: TrajectoryFormat,
) -> Result<Vec<T>, TrajectoryIoError> {
    let mut table = Vec::new();
    table.try_reserve_exact(frame_count).map_err(|_| {
        TrajectoryIoError::codec(
            TrajectoryIoErrorKind::LimitExceeded,
            TrajectoryIoOperation::Index,
            Some(format),
        )
        .with_detail("could not reserve frame index tables")
    })?;
    Ok(table)
}

fn format_code(format: TrajectoryFormat) -> u8 {
    match format {
        TrajectoryFormat::Xyz => 1,
        TrajectoryFormat::Dcd => 2,
        TrajectoryFormat::Trr => 3,
        TrajectoryFormat::Xtc => 4,
    }
}

fn array<const N: usize>(bytes: &[u8]) -> [u8; N] {
    bytes
        .try_into()
        .expect("encoded index fields have their documented widths")
}

fn invalid_index(detail: impl Into<String>) -> TrajectoryIoError {
    TrajectoryIoError::codec(
        TrajectoryIoErrorKind::InvalidIndex,
        TrajectoryIoOperation::Index,
        None,
    )
    .with_detail(detail)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// An index of `format` whose frames record a time and step exactly when
    /// `stamped`.
    fn index_of(
        format: TrajectoryFormat,
        starts: Vec<u64>,
        byte_count: u64,
        stamped: bool,
    ) -> TrajectoryFrameIndex {
        TrajectoryFrameIndex::new(
            format,
            TrajectoryFileFingerprint {
                byte_count,
                modified: Some(-1_500_000_000_123),
                content_hash: 0x0123_4567_89ab_cdef,
            },
            FrameOffsets {
                times: stamped.then(|| {
                    starts
                        .iter()
                        .map(|&start| start as f64 * 0.25 - 1.0)
                        .collect()
                }),
                steps: stamped.then(|| starts.iter().map(|&start| start * 10).collect()),
                starts,
                end: byte_count,
            },
            true,
        )
    }

    fn index(starts: Vec<u64>, byte_count: u64) -> TrajectoryFrameIndex {
        index_of(TrajectoryFormat::Trr, starts, byte_count, true)
    }

    fn reencode(bytes: &mut [u8]) {
        let body = bytes.len() - CHECKSUM_BYTES;
        let mut checksum = Fnv1a::new();
        checksum.write(&bytes[..body]);
        bytes[body..].copy_from_slice(&checksum.finish().to_le_bytes());
    }

    fn kind(bytes: &[u8], limits: &TrajectoryIoLimits) -> TrajectoryIoErrorKind {
        TrajectoryFrameIndex::from_bytes(bytes, limits)
            .unwrap_err()
            .kind()
    }

    #[test]
    fn encoding_round_trips_exactly() {
        let original = index(vec![0, 100, 250], 400);
        let bytes = original.to_bytes();
        assert_eq!(bytes.len(), HEADER_BYTES + 3 * 3 * 8 + CHECKSUM_BYTES);
        let decoded =
            TrajectoryFrameIndex::from_bytes(&bytes, &TrajectoryIoLimits::default()).unwrap();
        assert_eq!(decoded, original);
        assert_eq!(decoded.to_bytes(), bytes);
        assert_eq!(decoded.frame_offsets(), [0, 100, 250]);
        assert_eq!(decoded.steps(), Some(&[0, 1000, 2500][..]));
        let times = decoded.times().unwrap();
        assert_eq!(times.unit(), CANONICAL_TIME_UNIT);
        assert_eq!(*times.value(), [-1.0, 24.0, 61.5]);
        assert!(decoded.fingerprint().modified().unwrap() < UNIX_EPOCH);

        // Formats without times or steps, and empty indexes, have neither.
        for (format, starts) in [
            (TrajectoryFormat::Xyz, vec![0, 100]),
            (TrajectoryFormat::Xtc, Vec::new()),
        ] {
            let original = index_of(format, starts, 400, false);
            let bytes = original.to_bytes();
            assert_eq!(bytes[13], 0);
            let decoded =
                TrajectoryFrameIndex::from_bytes(&bytes, &TrajectoryIoLimits::default()).unwrap();
            assert_eq!(decoded, original);
            assert!(decoded.times().is_none() && decoded.steps().is_none());
        }
        // DCD records steps, and times only under a header-delta policy.
        let mut stepped = index_of(TrajectoryFormat::Dcd, vec![0, 100], 400, true);
        stepped.offsets.times = None;
        for original in [
            stepped,
            index_of(TrajectoryFormat::Dcd, vec![0, 100], 400, true),
        ] {
            let decoded = TrajectoryFrameIndex::from_bytes(
                &original.to_bytes(),
                &TrajectoryIoLimits::default(),
            )
            .unwrap();
            assert_eq!(decoded, original);
        }
    }

    #[test]
    fn version_1_encodings_are_refused() {
        // Version 1 had this header and the offsets alone.
        let mut bytes = index_of(TrajectoryFormat::Xtc, vec![0, 100], 400, false).to_bytes();
        bytes[8..12].copy_from_slice(&1_u32.to_le_bytes());
        reencode(&mut bytes);
        let error =
            TrajectoryFrameIndex::from_bytes(&bytes, &TrajectoryIoLimits::default()).unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidIndex);
        assert!(error.detail().unwrap().contains("version 1"));
    }

    #[test]
    fn steps_and_times_must_match_the_format() {
        let limits = TrajectoryIoLimits::default();
        let stamped = |format, stamped| index_of(format, vec![0, 100], 400, stamped).to_bytes();
        let with_flags = |mut bytes: Vec<u8>, flags| {
            bytes[13] = flags;
            reencode(&mut bytes);
            bytes
        };
        // Tables the format does not record, or missing ones it does.
        let invalid = [
            stamped(TrajectoryFormat::Xyz, true),
            with_flags(stamped(TrajectoryFormat::Xtc, true), FLAG_STEPS),
            with_flags(stamped(TrajectoryFormat::Trr, true), FLAG_TIMES),
            with_flags(stamped(TrajectoryFormat::Dcd, true), FLAG_TIMES),
            stamped(TrajectoryFormat::Xtc, false),
            with_flags(
                index_of(TrajectoryFormat::Xtc, Vec::new(), 400, false).to_bytes(),
                FLAG_STEPS | FLAG_TIMES,
            ),
        ];
        for bytes in invalid {
            assert_eq!(kind(&bytes, &limits), TrajectoryIoErrorKind::InvalidIndex);
        }
        // Non-finite times.
        for time in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut bytes = stamped(TrajectoryFormat::Xtc, true);
            let last_time = bytes.len() - CHECKSUM_BYTES - 8;
            bytes[last_time..last_time + 8].copy_from_slice(&time.to_le_bytes());
            reencode(&mut bytes);
            let error = TrajectoryFrameIndex::from_bytes(&bytes, &limits).unwrap_err();
            assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidIndex);
            assert!(error.detail().unwrap().contains("time 1"));
        }
    }

    #[test]
    fn every_frame_must_record_the_stamps_of_frame_0() {
        let limits = TrajectoryIoLimits::default();
        let stamp = FrameStamp {
            time: Some(1.5),
            step: Some(3),
        };
        let mut offsets = FrameOffsets::empty(0);
        offsets
            .push(0, 10, stamp, &limits, TrajectoryFormat::Xtc, "x")
            .unwrap();
        for different in [
            FrameStamp::default(),
            FrameStamp {
                time: None,
                ..stamp
            },
        ] {
            let error = offsets
                .push(10, 20, different, &limits, TrajectoryFormat::Xtc, "x")
                .unwrap_err();
            assert_eq!(error.kind(), TrajectoryIoErrorKind::InconsistentMetadata);
            assert_eq!(error.frame(), Some(1));
        }
        offsets
            .push(10, 20, stamp, &limits, TrajectoryFormat::Xtc, "x")
            .unwrap();
        assert_eq!(
            (offsets.times(), offsets.steps()),
            (Some(&[1.5, 1.5][..]), Some(&[3, 3][..]))
        );

        let verify = |stamp| {
            offsets.verify_stamp(
                1,
                stamp,
                TrajectoryIoOperation::ReadFrame,
                TrajectoryFormat::Xtc,
                "x",
            )
        };
        verify(stamp).unwrap();
        for different in [
            FrameStamp {
                time: Some(1.25),
                ..stamp
            },
            FrameStamp {
                step: Some(4),
                ..stamp
            },
            FrameStamp::default(),
        ] {
            let error = verify(different).unwrap_err();
            assert_eq!(error.kind(), TrajectoryIoErrorKind::IndexMismatch);
            assert_eq!(error.frame(), Some(1));
        }
    }

    #[test]
    fn malformed_encodings_are_typed_errors() {
        let limits = TrajectoryIoLimits::default();
        let bytes = index(vec![0, 100], 400).to_bytes();
        for length in 0..bytes.len() {
            assert_eq!(
                kind(&bytes[..length], &limits),
                TrajectoryIoErrorKind::InvalidIndex
            );
        }
        let mut corrupt = bytes.clone();
        corrupt[60] ^= 1;
        assert_eq!(kind(&corrupt, &limits), TrajectoryIoErrorKind::InvalidIndex);

        let mut version = bytes.clone();
        version[8] = 3;
        reencode(&mut version);
        let error = TrajectoryFrameIndex::from_bytes(&version, &limits).unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidIndex);
        assert!(error.detail().unwrap().contains("version 3"));

        for (position, value) in [(12, 9), (13, 2), (14, 2), (15, 1)] {
            let mut field = bytes.clone();
            field[position] = value;
            reencode(&mut field);
            assert_eq!(kind(&field, &limits), TrajectoryIoErrorKind::InvalidIndex);
        }
        let mut flags_on_xtc = bytes.clone();
        flags_on_xtc[12] = 4;
        reencode(&mut flags_on_xtc);
        assert_eq!(
            kind(&flags_on_xtc, &limits),
            TrajectoryIoErrorKind::InvalidIndex
        );

        for starts in [vec![100, 100], vec![100, 0], vec![0, 400]] {
            let bytes = index(starts, 400).to_bytes();
            assert_eq!(kind(&bytes, &limits), TrajectoryIoErrorKind::InvalidIndex);
        }
    }

    #[test]
    fn declared_frame_counts_are_bounded_before_allocation() {
        let bytes = index(vec![0, 100, 200], 400).to_bytes();
        let error = TrajectoryFrameIndex::from_bytes(
            &bytes,
            &TrajectoryIoLimits::default().with_max_index_entries(2),
        )
        .unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::LimitExceeded);
        assert_eq!(error.limit_exceeded().unwrap().observed(), 3);

        let mut huge = bytes;
        huge[48..56].copy_from_slice(&u64::MAX.to_le_bytes());
        reencode(&mut huge);
        assert_eq!(
            kind(&huge, &TrajectoryIoLimits::default()),
            TrajectoryIoErrorKind::LimitExceeded
        );
        let mut short = index(vec![0, 100, 200], 400).to_bytes();
        short[48..56].copy_from_slice(&2_u64.to_le_bytes());
        reencode(&mut short);
        assert_eq!(
            kind(&short, &TrajectoryIoLimits::default()),
            TrajectoryIoErrorKind::InvalidIndex
        );
    }

    #[test]
    fn sample_ranges_cover_small_files_whole_and_stay_inside_large_ones() {
        assert_eq!(sample_ranges(0).collect::<Vec<_>>(), [(0, 0)]);
        assert_eq!(
            sample_ranges(WHOLE_FILE_HASH_BYTES).collect::<Vec<_>>(),
            [(0, WHOLE_FILE_HASH_BYTES)]
        );
        let large = WHOLE_FILE_HASH_BYTES + 1;
        let ranges = sample_ranges(large).collect::<Vec<_>>();
        assert_eq!(ranges.len(), 18);
        assert!(ranges
            .iter()
            .all(|&(start, length)| start + length <= large));
    }

    #[test]
    fn frame_ends_must_meet_the_next_start() {
        let offsets = FrameOffsets {
            starts: vec![10, 30],
            end: 50,
            times: None,
            steps: None,
        };
        assert_eq!(offsets.frame(0), Some((10, 30)));
        assert_eq!(offsets.frame(1), Some((30, 50)));
        assert_eq!(offsets.frame(2), None);
        offsets
            .verify_frame_end(1, 50, TrajectoryFormat::Xtc, "x")
            .unwrap();
        let error = offsets
            .verify_frame_end(0, 29, TrajectoryFormat::Xtc, "x")
            .unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::IndexMismatch);
        assert_eq!((error.expected(), error.actual()), (Some(30), Some(29)));
    }
}
