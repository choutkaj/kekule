//! Continuing a trajectory file after a restart: keep the frames a
//! checkpoint covers, remove the rest, and append after them.

use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Cursor, Seek, SeekFrom};
use std::path::Path;
use std::sync::Arc;

use kekule::topology::Topology;

use crate::{
    TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation,
    TrajectoryWriter,
};

use super::file_writer::{resolve_format, resumed_stream_writer};
use super::{
    codec_context, index_mismatch, io_context, open_indexed_trajectory, TrajectoryFileFingerprint,
    TrajectoryFormatHint, TrajectoryOpenOptions, TrajectoryStreamWriter, TrajectoryWriteOptions,
};

/// Which frames of a trajectory file [`append_trajectory`] keeps; it removes
/// every frame after them.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum KeepFrames {
    /// The first this many frames. The file must hold at least as many
    /// complete frames.
    First(usize),
    /// The frames before the first whose step exceeds this one: those a run
    /// had written when it checkpointed at this step, which a run restarted
    /// from the checkpoint does not write again. Only formats whose frames
    /// record a step keep frames by step: XTC, TRR, and DCD, whose frames
    /// follow its header's step sequence.
    ThroughStep(u64),
}

/// What [`append_trajectory`] kept of a trajectory file, and what it
/// removed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrajectoryAppendReport {
    kept_frames: usize,
    last_kept_step: Option<u64>,
    dropped_frames: usize,
    dropped_bytes: u64,
}

impl TrajectoryAppendReport {
    /// The frames the file still holds, which the writer appends after.
    pub const fn kept_frames(&self) -> usize {
        self.kept_frames
    }

    /// The step of the last kept frame, when a frame was kept and the format
    /// records steps; a run compares it with its checkpoint to see that no
    /// frame it counts is missing.
    pub const fn last_kept_step(&self) -> Option<u64> {
        self.last_kept_step
    }

    /// The complete frames removed after the kept ones.
    pub const fn dropped_frames(&self) -> usize {
        self.dropped_frames
    }

    /// The bytes removed: those of the dropped frames and of a frame the
    /// file's writer had not finished, or the whole file when it did not yet
    /// hold a complete frame, or for DCD a complete header.
    pub const fn dropped_bytes(&self) -> u64 {
        self.dropped_bytes
    }
}

/// The frames a file holds and where the kept ones end.
struct Kept {
    frames: usize,
    last_step: Option<u64>,
    dropped_frames: usize,
    /// The end of the kept frames, and of a DCD's header without them.
    end: u64,
    /// The fingerprint of the frames the file holds, which must not change
    /// before the file is cut; `None` for a file started again.
    fingerprint: Option<TrajectoryFileFingerprint>,
}

/// Opens the trajectory file at `path` to continue it, as a run restarted
/// from a checkpoint does: verifies and indexes the frames it holds, keeps
/// those `keep` selects, removes the rest, and returns a writer that appends
/// after the kept frames, with a report of what was removed.
///
/// A run checkpointed at step N and restarted from there passes
/// [`KeepFrames::ThroughStep`]`(N)`, or [`KeepFrames::First`] with its frame
/// count for XYZ, and writes the frames after the checkpoint as it computes
/// them; with the `options` the file was written with, the file then holds,
/// byte for byte, what one uninterrupted run would have written.
///
/// The frames are verified as [`open_indexed_trajectory`] verifies those of a
/// file still being written ([`TrajectoryOpenOptions::with_growing`]): a
/// final frame its writer did not finish is no error and goes with the
/// frames removed, a file without a complete frame, or for DCD a complete
/// header, is started again, and corrupt data anywhere is an error. Where the
/// file records its writer's settings they must be those of `options`: the
/// DCD header (endianness, step sequence, `DELTA`, unit-cell flag, title),
/// and the magic and precision of the first XTC frame or the scalar width of
/// the TRR frames. The returned writer numbers its frames after the kept
/// ones, as the writer that wrote them would have
/// ([`super::FrameTimingPolicy::Number`],
/// [`super::dcd::DcdFrameTimingPolicy`]), and its
/// [`TrajectoryStreamWriter::report`] covers the frames it appends. `options`
/// takes the format, or the file extension does; its overwrite policy does
/// not apply.
///
/// The file is cut after the kept frames before the writer is returned, and
/// a DCD header then counts them, so the file holds a complete trajectory of
/// them. Appending survives a process that dies at any point: before the
/// cut, the file is as it was; after it, it holds the kept frames and maybe
/// some appended ones, possibly with a partial frame or a stale DCD frame
/// count, which another call takes as this one took the first file. Flush
/// the writer after each frame ([`TrajectoryWriter::flush`]) to let a reader
/// follow the file. A reader following it under
/// [`TrajectoryOpenOptions::with_growing`] whose index covers removed frames
/// gets [`TrajectoryIoErrorKind::IndexMismatch`] from its next refresh; it
/// opens the file again and refreshes to pick up the appended frames. No
/// other process may write the file meanwhile. For the frames a checkpoint
/// counts to survive the failure of a node, not only of a process, flush the
/// writer and sync the file (`writer.writer().get_ref().sync_data()`) before
/// recording the checkpoint.
///
/// # Errors
///
/// Returns [`TrajectoryIoErrorKind::UnknownFormat`] when neither `options`
/// nor the extension names a format;
/// [`TrajectoryIoErrorKind::UnsupportedVariant`] for
/// [`KeepFrames::ThroughStep`] on XYZ;
/// [`TrajectoryIoErrorKind::InconsistentMetadata`] when the file holds fewer
/// complete frames than [`KeepFrames::First`] asks for, or records other
/// writer settings than `options`; the decoding error of the first corrupt
/// frame and the errors of [`open_indexed_trajectory`], as for a file that
/// does not fit the topology; [`TrajectoryIoErrorKind::IndexMismatch`] when
/// the file changes before it is cut; the option errors of the codec
/// writer's constructor; and [`TrajectoryIoErrorKind::Io`] when the file
/// does not exist or cannot be read, cut, or written. The file changes only
/// once every check has passed; an I/O error after that leaves it cut after
/// the kept frames, which another call completes.
///
/// ```no_run
/// use std::sync::Arc;
///
/// use kekule::{smiles, topology::Topology};
/// use kekule_traj::io::{append_trajectory, KeepFrames, TrajectoryWriteOptions};
/// use kekule_traj::{TrajectoryFormat, TrajectoryWriter};
///
/// let molecule = smiles::to_molecules("CC")?.pop().unwrap();
/// let topology = Arc::new(Topology::from_molecule(molecule)?);
/// // The run was checkpointed at step 5000 and restarts from there.
/// let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc);
/// let (mut writer, report) =
///     append_trajectory("run.xtc", topology, KeepFrames::ThroughStep(5000), options)?;
/// eprintln!(
///     "kept {} frames, dropped {} frames and {} bytes",
///     report.kept_frames(),
///     report.dropped_frames(),
///     report.dropped_bytes()
/// );
/// # let frames: Vec<kekule::structure::TrajectoryFrame> = Vec::new();
/// # let trajectory = kekule::structure::Trajectory::from_items(writer.shared_topology(), frames)?;
/// for frame in &trajectory {
///     writer.write_frame(frame)?;
///     writer.flush()?;
/// }
/// writer.finish()?;
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn append_trajectory(
    path: impl AsRef<Path>,
    topology: impl Into<Arc<Topology>>,
    keep: KeepFrames,
    options: TrajectoryWriteOptions,
) -> Result<
    (
        TrajectoryStreamWriter<BufWriter<File>>,
        TrajectoryAppendReport,
    ),
    TrajectoryIoError,
> {
    let path = path.as_ref();
    let topology = topology.into();
    let label = path.display().to_string();
    let format = resolve_format(path, options.format())?;
    if matches!(keep, KeepFrames::ThroughStep(_)) && format == TrajectoryFormat::Xyz {
        return Err(codec_context(
            TrajectoryIoErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::Open,
            Some(format),
            &label,
            "XYZ frames record no step; keep frames by count with KeepFrames::First",
        ));
    }
    let kept = kept_frames(path, Arc::clone(&topology), format, keep, &options, &label)?;
    // The writer's own option errors come before the file changes.
    resumed_stream_writer(
        format,
        Cursor::new(Vec::new()),
        Arc::clone(&topology),
        options.clone(),
        label.clone(),
        kept.frames,
        kept.fingerprint.is_some(),
    )?;
    let io = |error| io_context(TrajectoryIoOperation::Open, Some(format), &label, error);
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .open(path)
        .map_err(io)?;
    if let Some(fingerprint) = &kept.fingerprint {
        if !fingerprint.is_prefix_of_file(&mut file)? {
            return Err(index_mismatch(
                TrajectoryIoOperation::Open,
                format,
                &label,
                "the file changed before it was cut",
            ));
        }
    }
    let length = file.metadata().map_err(io)?.len();
    file.set_len(kept.end).map_err(io)?;
    let mut stream = BufWriter::new(file);
    stream.seek(SeekFrom::Start(kept.end)).map_err(io)?;
    let mut writer = resumed_stream_writer(
        format,
        stream,
        topology,
        options,
        label,
        kept.frames,
        kept.fingerprint.is_some(),
    )?;
    // A DCD header counts the kept frames from here on.
    writer.flush()?;
    let report = TrajectoryAppendReport {
        kept_frames: kept.frames,
        last_kept_step: kept.last_step,
        dropped_frames: kept.dropped_frames,
        dropped_bytes: length.saturating_sub(kept.end),
    };
    Ok((writer, report))
}

/// Scans the file and chooses the frames to keep, checking that `options`
/// continue them.
fn kept_frames(
    path: &Path,
    topology: Arc<Topology>,
    format: TrajectoryFormat,
    keep: KeepFrames,
    options: &TrajectoryWriteOptions,
    label: &str,
) -> Result<Kept, TrajectoryIoError> {
    let open = TrajectoryOpenOptions::default()
        .with_format_hint(TrajectoryFormatHint::Explicit(format))
        .with_growing(true);
    let mut reader = match open_indexed_trajectory(path, topology, open) {
        Ok(reader) => reader,
        // Nothing complete to keep: the file starts again.
        Err(error) if error.kind() == TrajectoryIoErrorKind::Incomplete => {
            if let KeepFrames::First(frames @ 1..) = keep {
                return Err(too_few_frames(frames, 0, format, label));
            }
            return Ok(Kept {
                frames: 0,
                last_step: None,
                dropped_frames: 0,
                end: 0,
                fingerprint: None,
            });
        }
        Err(error) => return Err(error),
    };
    let index = reader.to_frame_index();
    let count = index.frame_count();
    let steps = index.steps();
    let frames = match keep {
        KeepFrames::First(frames) if frames > count => {
            return Err(too_few_frames(frames, count, format, label))
        }
        KeepFrames::First(frames) => frames,
        KeepFrames::ThroughStep(last) => steps
            .unwrap_or_default()
            .iter()
            .take_while(|&&step| step <= last)
            .count(),
    };
    // The kept frames, and a DCD's header, stay: their writer's settings
    // must be those of `options`.
    if frames > 0 || format == TrajectoryFormat::Dcd {
        reader.check_append(options)?;
    }
    let fingerprint = index.fingerprint().clone();
    let end = index
        .frame_offsets()
        .get(frames)
        .copied()
        .unwrap_or_else(|| fingerprint.byte_count());
    Ok(Kept {
        frames,
        last_step: frames
            .checked_sub(1)
            .and_then(|last| steps.map(|steps| steps[last])),
        dropped_frames: count - frames,
        end,
        fingerprint: Some(fingerprint),
    })
}

fn too_few_frames(
    frames: usize,
    count: usize,
    format: TrajectoryFormat,
    label: &str,
) -> TrajectoryIoError {
    codec_context(
        TrajectoryIoErrorKind::InconsistentMetadata,
        TrajectoryIoOperation::Open,
        Some(format),
        label,
        format!("cannot keep {frames} frames of a file that holds {count} complete frames"),
    )
    .with_counts(frames as u64, count as u64)
}
