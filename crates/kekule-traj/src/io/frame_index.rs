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
use std::path::Path;
use std::str::FromStr;
use std::sync::Arc;
use std::time::{Duration, SystemTime, UNIX_EPOCH};

use kekule::diagnostics::{Cancellation, LimitExceeded, Resource};

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

const INDEX_MAGIC: [u8; 8] = *b"KEKTRIDX";
const FINGERPRINT_MAGIC: [u8; 8] = *b"KEKTRFPR";
/// The version of both the index and the fingerprint encodings.
const VERSION: u32 = 1;
const HEADER_BYTES: usize = 56;
const CHECKSUM_BYTES: usize = 8;
const FINGERPRINT_BYTES: usize = 56;
/// Where the fingerprint fields lie in both encodings.
const FINGERPRINT_FIELDS: std::ops::Range<usize> = 14..48;
const FLAG_MIXED_PRECISION: u8 = 1;

/// The identity of a trajectory file's content.
///
/// It records the file's size, its modification time where the platform
/// reports one, and a 64-bit FNV-1a hash of the size followed by sampled
/// content: the whole file up to 256 KiB; otherwise its first and last
/// 64 KiB and sixteen evenly spaced 4 KiB interior blocks. It detects a
/// replaced, appended, truncated, or rewritten file; it is not a
/// cryptographic digest, and the per-frame checks of an index-backed reader
/// remain the guard against frames that are not where the index says.
///
/// Every [`TrajectoryFrameIndex`] carries the fingerprint of the file it was
/// built for. A host can also compute one with [`Self::of_path`] or
/// [`Self::of_file`], check a file against one with [`Self::matches_path`] or
/// [`Self::matches_file`], and keep one, for example in a project file or a
/// cache key, through its stable encoding.
///
/// Equality compares all three fields, as opening from an index does, so a
/// file that is touched, or copied without its modification time, no longer
/// matches. [`Self::same_content`] ignores the modification time, for a host
/// asking whether a moved or copied file is the one it recorded.
///
/// # Encoding
///
/// [`Self::to_bytes`] writes version 1 of a little-endian layout. It is
/// versioned together with the [`TrajectoryFrameIndex`] encoding, whose bytes
/// 14 to 47 hold the same fields at the same offsets:
///
/// | Bytes | Field |
/// |---|---|
/// | 8 | magic `KEKTRFPR` |
/// | 4 | version, `u32` 1 |
/// | 2 | reserved, 0 |
/// | 1 | 1 when the modification time is present, else 0 |
/// | 1 | reserved, 0 |
/// | 8 | file size in bytes, `u64` |
/// | 16 | modification time in nanoseconds from the Unix epoch, `i128`; 0 when absent |
/// | 8 | content hash, `u64` |
/// | 8 | FNV-1a 64 checksum of all preceding bytes, `u64` |
///
/// [`Display`](fmt::Display) writes these 56 bytes as 112 lowercase
/// hexadecimal digits, a form that fits a text file, a file name, or a cache
/// key, and [`FromStr`] reads it back. Both decoders take untrusted input:
/// anything but a valid version 1 encoding is a
/// [`TrajectoryIoErrorKind::InvalidIndex`] error, never a panic.
///
/// # Examples
///
/// ```no_run
/// use kekule_traj::io::TrajectoryFileFingerprint;
///
/// let recorded = TrajectoryFileFingerprint::of_path("run.xtc")?.to_string();
///
/// // Later, after the user relinks the trajectory:
/// let recorded: TrajectoryFileFingerprint = recorded.parse()?;
/// let relinked = TrajectoryFileFingerprint::of_path("moved/run.xtc")?;
/// assert!(recorded.same_content(&relinked));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct TrajectoryFileFingerprint {
    byte_count: u64,
    /// Nanoseconds from the Unix epoch, negative before it.
    modified: Option<i128>,
    content_hash: u64,
}

impl TrajectoryFileFingerprint {
    /// Fingerprints the file at `path`.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::Io`] when the file cannot be opened
    /// or read.
    pub fn of_path(path: impl AsRef<Path>) -> Result<Self, TrajectoryIoError> {
        let path = path.as_ref();
        let label = path.display().to_string();
        let mut file = File::open(path)
            .map_err(|error| io_context(TrajectoryIoOperation::Open, None, &label, error))?;
        Self::of_labelled_file(&mut file, &label)
    }

    /// Fingerprints an open file from its start and leaves it positioned at
    /// its start.
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::Io`] when the file cannot be read or
    /// repositioned.
    pub fn of_file(file: &mut File) -> Result<Self, TrajectoryIoError> {
        Self::fingerprint(file)
            .map_err(|error| TrajectoryIoError::io(TrajectoryIoOperation::Open, error))
    }

    /// As [`Self::of_file`], labelling errors with `source_label`.
    pub(crate) fn of_labelled_file(
        file: &mut File,
        source_label: &str,
    ) -> Result<Self, TrajectoryIoError> {
        Self::fingerprint(file)
            .map_err(|error| io_context(TrajectoryIoOperation::Open, None, source_label, error))
    }

    fn fingerprint(file: &mut File) -> std::io::Result<Self> {
        let metadata = file.metadata()?;
        let byte_count = metadata.len();
        let modified = metadata.modified().ok().and_then(unix_nanos);
        let mut hash = Fnv1a::new();
        hash.write(&byte_count.to_le_bytes());
        let mut buffer = [0_u8; 8192];
        for (start, length) in sample_ranges(byte_count) {
            file.seek(SeekFrom::Start(start))?;
            let mut remaining = length;
            while remaining > 0 {
                let chunk = usize::try_from(remaining.min(buffer.len() as u64))
                    .expect("a chunk is no larger than the buffer");
                file.read_exact(&mut buffer[..chunk])?;
                hash.write(&buffer[..chunk]);
                remaining -= chunk as u64;
            }
        }
        file.seek(SeekFrom::Start(0))?;
        Ok(Self {
            byte_count,
            modified,
            content_hash: hash.finish(),
        })
    }

    /// Whether the file at `path` has this fingerprint.
    ///
    /// # Errors
    ///
    /// As [`Self::of_path`].
    pub fn matches_path(&self, path: impl AsRef<Path>) -> Result<bool, TrajectoryIoError> {
        Ok(Self::of_path(path)? == *self)
    }

    /// Whether an open file has this fingerprint; the file is left
    /// positioned at its start.
    ///
    /// # Errors
    ///
    /// As [`Self::of_file`].
    pub fn matches_file(&self, file: &mut File) -> Result<bool, TrajectoryIoError> {
        Ok(Self::of_file(file)? == *self)
    }

    /// Whether both fingerprints have the same size and content hash,
    /// whatever their modification times.
    pub const fn same_content(&self, other: &Self) -> bool {
        self.byte_count == other.byte_count && self.content_hash == other.content_hash
    }

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

    /// Encodes the fingerprint in the versioned layout described on the type.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(FINGERPRINT_BYTES);
        bytes.extend_from_slice(&FINGERPRINT_MAGIC);
        bytes.extend_from_slice(&VERSION.to_le_bytes());
        bytes.extend_from_slice(&[0, 0]);
        self.write_fields(&mut bytes);
        push_checksum(&mut bytes);
        bytes
    }

    /// Decodes a fingerprint written by [`Self::to_bytes`].
    ///
    /// # Errors
    ///
    /// Returns [`TrajectoryIoErrorKind::InvalidIndex`] for a truncated,
    /// corrupt, or malformed encoding or an unsupported version.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, TrajectoryIoError> {
        const WHAT: &str = "file fingerprint";
        checked_body(bytes, FINGERPRINT_MAGIC, FINGERPRINT_BYTES, WHAT)?;
        if bytes.len() != FINGERPRINT_BYTES {
            return Err(invalid_index("file fingerprint has trailing bytes"));
        }
        if bytes[12..14] != [0, 0] {
            return Err(invalid_index("file fingerprint reserved bytes are nonzero"));
        }
        Self::read_fields(&bytes[FINGERPRINT_FIELDS], WHAT)
    }

    /// Writes the fields at [`FINGERPRINT_FIELDS`] of both encodings.
    fn write_fields(&self, bytes: &mut Vec<u8>) {
        bytes.push(u8::from(self.modified.is_some()));
        bytes.push(0);
        bytes.extend_from_slice(&self.byte_count.to_le_bytes());
        bytes.extend_from_slice(&self.modified.unwrap_or(0).to_le_bytes());
        bytes.extend_from_slice(&self.content_hash.to_le_bytes());
    }

    /// Reads the fields written by [`Self::write_fields`] into an encoded
    /// `what`.
    fn read_fields(fields: &[u8], what: &str) -> Result<Self, TrajectoryIoError> {
        let modified = i128::from_le_bytes(array(&fields[10..26]));
        let modified = match fields[0] {
            1 => Some(modified),
            0 if modified == 0 => None,
            _ => {
                return Err(invalid_index(format!(
                    "{what} modification time is malformed"
                )))
            }
        };
        if fields[1] != 0 {
            return Err(invalid_index(format!("{what} reserved byte is nonzero")));
        }
        Ok(Self {
            byte_count: u64::from_le_bytes(array(&fields[2..10])),
            modified,
            content_hash: u64::from_le_bytes(array(&fields[26..34])),
        })
    }
}

impl fmt::Display for TrajectoryFileFingerprint {
    /// Writes the encoding of [`Self::to_bytes`] as lowercase hexadecimal.
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.to_bytes()
            .iter()
            .try_for_each(|byte| write!(formatter, "{byte:02x}"))
    }
}

impl FromStr for TrajectoryFileFingerprint {
    type Err = TrajectoryIoError;

    /// Reads the lowercase hexadecimal form written by
    /// [`Display`](fmt::Display).
    fn from_str(text: &str) -> Result<Self, Self::Err> {
        if text.len() != 2 * FINGERPRINT_BYTES {
            return Err(invalid_index(
                "file fingerprint text is not 112 hexadecimal digits",
            ));
        }
        let digit = |byte: u8| match byte {
            b'0'..=b'9' => Some(byte - b'0'),
            b'a'..=b'f' => Some(byte - b'a' + 10),
            _ => None,
        };
        let bytes = text
            .as_bytes()
            .as_chunks::<2>()
            .0
            .iter()
            .map(|&[high, low]| Some((digit(high)? << 4) | digit(low)?))
            .collect::<Option<Vec<u8>>>()
            .ok_or_else(|| invalid_index("file fingerprint text is not lowercase hexadecimal"))?;
        Self::from_bytes(&bytes)
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

/// Frame start offsets and the end of the last frame of one indexed stream.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct FrameOffsets {
    starts: Vec<u64>,
    end: u64,
}

impl FrameOffsets {
    /// An empty index of a stream whose frames would start at `start`.
    pub(crate) const fn empty(start: u64) -> Self {
        Self {
            starts: Vec::new(),
            end: start,
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

    /// Appends a verified frame occupying `start..end` within the limits.
    pub(crate) fn push(
        &mut self,
        start: u64,
        end: u64,
        limits: &TrajectoryIoLimits,
        format: TrajectoryFormat,
        source_label: &str,
    ) -> Result<(), TrajectoryIoError> {
        let frame = self.starts.len();
        reserve_index_for_push(&mut self.starts, limits, format, source_label, frame)?;
        self.starts.push(start);
        self.end = end;
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
/// It holds the format, the byte offset of every frame, the frame count, and
/// the [`TrajectoryFileFingerprint`] of the file it was built for. Export one
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
/// next frame. A stale or corrupted index is therefore a typed
/// [`TrajectoryIoErrorKind::IndexMismatch`] or
/// [`TrajectoryIoErrorKind::InvalidIndex`] error, never a wrong frame. Only
/// the index records which frame is which, so an index deliberately forged
/// with a valid checksum can at most renumber whole, valid frames of the file;
/// keep indexes where only the host writes them.
///
/// # Encoding
///
/// [`Self::to_bytes`] writes version 1 of a little-endian layout:
///
/// | Bytes | Field |
/// |---|---|
/// | 8 | magic `KEKTRIDX` |
/// | 4 | version, `u32` 1 |
/// | 1 | format: 1 XYZ, 2 DCD, 3 TRR, 4 XTC |
/// | 1 | flags: bit 0 set when TRR frames mix f32 and f64 scalars; other bits 0 |
/// | 1 | 1 when the modification time is present, else 0 |
/// | 1 | reserved, 0 |
/// | 8 | file size in bytes, `u64` |
/// | 16 | modification time in nanoseconds from the Unix epoch, `i128`; 0 when absent |
/// | 8 | content hash, `u64` |
/// | 8 | frame count `n`, `u64` |
/// | 8 `n` | strictly increasing frame start offsets, `u64`, each below the file size |
/// | 8 | FNV-1a 64 checksum of all preceding bytes, `u64` |
///
/// The last frame ends at the end of the file. Bytes 14 to 47 are the
/// fingerprint fields of the [`TrajectoryFileFingerprint`] encoding, which
/// shares this version. Other versions are rejected with
/// [`TrajectoryIoErrorKind::InvalidIndex`].
///
/// [`IndexedFileTrajectoryReader::to_frame_index`]:
///     super::IndexedFileTrajectoryReader::to_frame_index
/// [`open_indexed_trajectory_from_index`]: super::open_indexed_trajectory_from_index
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

    /// The fingerprint of the file the index was built for.
    pub const fn fingerprint(&self) -> &TrajectoryFileFingerprint {
        &self.fingerprint
    }

    /// Encodes the index in the versioned layout described on the type.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::with_capacity(
            HEADER_BYTES + CHECKSUM_BYTES + self.offsets.len() * std::mem::size_of::<u64>(),
        );
        bytes.extend_from_slice(&INDEX_MAGIC);
        bytes.extend_from_slice(&VERSION.to_le_bytes());
        bytes.push(format_code(self.format));
        bytes.push(if self.precision_mixed {
            FLAG_MIXED_PRECISION
        } else {
            0
        });
        self.fingerprint.write_fields(&mut bytes);
        bytes.extend_from_slice(&(self.offsets.len() as u64).to_le_bytes());
        for start in self.offsets.starts() {
            bytes.extend_from_slice(&start.to_le_bytes());
        }
        push_checksum(&mut bytes);
        bytes
    }

    /// Decodes an index written by [`Self::to_bytes`].
    ///
    /// The frame count is checked against the frame and index limits of
    /// `limits` before anything is allocated, and the offsets must be strictly
    /// increasing and lie inside the fingerprinted file.
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
        let body = checked_body(
            bytes,
            INDEX_MAGIC,
            HEADER_BYTES + CHECKSUM_BYTES,
            "frame index",
        )?;
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
        if flags & !FLAG_MIXED_PRECISION != 0 || (flags != 0 && format != TrajectoryFormat::Trr) {
            return Err(invalid_index(format!(
                "frame index flags {flags:#04x} are invalid for {format}"
            )));
        }
        let fingerprint =
            TrajectoryFileFingerprint::read_fields(&bytes[FINGERPRINT_FIELDS], "frame index")?;
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
        let offset_bytes = &body[HEADER_BYTES..];
        if offset_bytes.len() / std::mem::size_of::<u64>() != frame_count
            || offset_bytes.len() % std::mem::size_of::<u64>() != 0
        {
            return Err(invalid_index(
                "frame index length disagrees with its frame count",
            ));
        }
        let mut starts = Vec::new();
        starts.try_reserve_exact(frame_count).map_err(|_| {
            TrajectoryIoError::codec(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::Index,
                Some(format),
            )
            .with_detail("could not reserve frame index offsets")
        })?;
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
        Ok(Self {
            format,
            offsets: FrameOffsets {
                starts,
                end: fingerprint.byte_count,
            },
            fingerprint,
            precision_mixed: flags & FLAG_MIXED_PRECISION != 0,
        })
    }
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
        .expect("encoded fields have their documented widths")
}

/// Appends the FNV-1a 64 checksum of `bytes`.
fn push_checksum(bytes: &mut Vec<u8>) {
    let mut checksum = Fnv1a::new();
    checksum.write(bytes);
    bytes.extend_from_slice(&checksum.finish().to_le_bytes());
}

/// Checks the length, magic, version, and checksum of an encoded `what` of
/// at least `min_len` bytes and returns the bytes before its checksum.
fn checked_body<'a>(
    bytes: &'a [u8],
    magic: [u8; 8],
    min_len: usize,
    what: &str,
) -> Result<&'a [u8], TrajectoryIoError> {
    if bytes.len() < min_len {
        return Err(invalid_index(format!("{what} is truncated")));
    }
    if bytes[..8] != magic {
        return Err(invalid_index(format!(
            "{what} has no {} magic",
            String::from_utf8_lossy(&magic)
        )));
    }
    let version = u32::from_le_bytes(array(&bytes[8..12]));
    if version != VERSION {
        return Err(invalid_index(format!(
            "{what} version {version} is not supported"
        )));
    }
    let (body, checksum) = bytes.split_at(bytes.len() - CHECKSUM_BYTES);
    let mut expected = Fnv1a::new();
    expected.write(body);
    if expected.finish() != u64::from_le_bytes(array(checksum)) {
        return Err(invalid_index(format!("{what} checksum does not match")));
    }
    Ok(body)
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

    fn index(starts: Vec<u64>, byte_count: u64) -> TrajectoryFrameIndex {
        TrajectoryFrameIndex::new(
            TrajectoryFormat::Trr,
            TrajectoryFileFingerprint {
                byte_count,
                modified: Some(-1_500_000_000_123),
                content_hash: 0x0123_4567_89ab_cdef,
            },
            FrameOffsets {
                starts,
                end: byte_count,
            },
            true,
        )
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
        assert_eq!(bytes.len(), HEADER_BYTES + 3 * 8 + CHECKSUM_BYTES);
        let decoded =
            TrajectoryFrameIndex::from_bytes(&bytes, &TrajectoryIoLimits::default()).unwrap();
        assert_eq!(decoded, original);
        assert_eq!(decoded.to_bytes(), bytes);
        assert_eq!(decoded.frame_offsets(), [0, 100, 250]);
        assert!(decoded.fingerprint().modified().unwrap() < UNIX_EPOCH);
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
        version[8] = 2;
        reencode(&mut version);
        let error = TrajectoryFrameIndex::from_bytes(&version, &limits).unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidIndex);
        assert!(error.detail().unwrap().contains("version 2"));

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
    fn fingerprint_encodings_round_trip_and_match_the_index_fields() {
        let index = index(vec![0, 100], 400);
        let fingerprint = index.fingerprint();
        let bytes = fingerprint.to_bytes();
        assert_eq!(bytes.len(), FINGERPRINT_BYTES);
        assert_eq!(bytes[..8], FINGERPRINT_MAGIC);
        assert_eq!(
            bytes[FINGERPRINT_FIELDS],
            index.to_bytes()[FINGERPRINT_FIELDS]
        );
        assert_eq!(
            &TrajectoryFileFingerprint::from_bytes(&bytes).unwrap(),
            fingerprint
        );
        let text = fingerprint.to_string();
        assert_eq!(text.len(), 2 * FINGERPRINT_BYTES);
        assert_eq!(
            &text.parse::<TrajectoryFileFingerprint>().unwrap(),
            fingerprint
        );

        let absent = TrajectoryFileFingerprint {
            modified: None,
            ..fingerprint.clone()
        };
        let decoded = TrajectoryFileFingerprint::from_bytes(&absent.to_bytes()).unwrap();
        assert_eq!(decoded, absent);
        assert!(decoded.same_content(fingerprint));
        assert_ne!(&decoded, fingerprint);
    }

    #[test]
    fn malformed_fingerprints_are_typed_errors() {
        let fingerprint = index(vec![0], 400).fingerprint().clone();
        let bytes = fingerprint.to_bytes();
        let kind = |bytes: &[u8]| {
            TrajectoryFileFingerprint::from_bytes(bytes)
                .unwrap_err()
                .kind()
        };
        for length in 0..bytes.len() {
            assert_eq!(kind(&bytes[..length]), TrajectoryIoErrorKind::InvalidIndex);
        }
        let mut long = bytes.clone();
        long.extend_from_slice(&[0; 8]);
        reencode(&mut long);
        assert_eq!(kind(&long), TrajectoryIoErrorKind::InvalidIndex);
        let mut corrupt = bytes.clone();
        corrupt[20] ^= 1;
        assert_eq!(kind(&corrupt), TrajectoryIoErrorKind::InvalidIndex);
        assert_eq!(
            kind(&index(vec![0], 400).to_bytes()),
            TrajectoryIoErrorKind::InvalidIndex
        );
        for (position, value) in [(8, 2), (12, 1), (13, 1), (14, 2), (15, 1)] {
            let mut field = bytes.clone();
            field[position] = value;
            reencode(&mut field);
            assert_eq!(kind(&field), TrajectoryIoErrorKind::InvalidIndex);
        }

        let text = fingerprint.to_string();
        for bad in [
            String::new(),
            text[1..].to_owned(),
            text.to_uppercase(),
            format!("{}g", &text[1..]),
            format!("{}é", &text[2..]),
            format!(
                "{}{}",
                &text[..text.len() - 1],
                if text.ends_with('0') { '1' } else { '0' }
            ),
        ] {
            assert_eq!(
                bad.parse::<TrajectoryFileFingerprint>().unwrap_err().kind(),
                TrajectoryIoErrorKind::InvalidIndex,
                "{bad}"
            );
        }
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
