use crate::{TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation};
use std::io;

pub(crate) fn frame_offset_context(
    mut error: TrajectoryIoError,
    frame: u64,
    byte_offset: u64,
) -> TrajectoryIoError {
    if matches!(
        error.kind(),
        TrajectoryIoErrorKind::InvalidStructure
            | TrajectoryIoErrorKind::InvalidTrajectory
            | TrajectoryIoErrorKind::TopologyMismatch
            | TrajectoryIoErrorKind::FrameIndexOutOfRange
            | TrajectoryIoErrorKind::UnsupportedField
    ) {
        return error;
    }
    if error.frame().is_none() {
        error = error.with_frame(frame);
    }
    if error.byte_offset().is_none() {
        error = error.with_byte_offset(byte_offset);
    }
    error
}

pub(crate) fn probe_seekable_eof<R: io::Read + io::Seek>(
    reader: &mut R,
    operation: TrajectoryIoOperation,
    format: TrajectoryFormat,
    source_label: &str,
) -> Result<bool, TrajectoryIoError> {
    let start = reader
        .stream_position()
        .map_err(|error| io_context(operation, Some(format), source_label, error))?;
    let mut byte = [0_u8; 1];
    let read = reader.read(&mut byte);
    let restore = reader
        .seek(io::SeekFrom::Start(start))
        .map_err(|error| io_context(operation, Some(format), source_label, error));
    match (read, restore) {
        (_, Err(error)) => Err(error),
        (Err(error), Ok(_)) => Err(io_context(operation, Some(format), source_label, error)),
        (Ok(0), Ok(_)) => Ok(true),
        (Ok(_), Ok(_)) => Ok(false),
    }
}

pub(crate) fn require_nonempty_writer(
    frame_count: u64,
    format: TrajectoryFormat,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    if frame_count == 0 {
        return Err(codec_context(
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::Finish,
            Some(format),
            source_label,
            format!("{format} production trajectories must contain at least one frame"),
        ));
    }
    Ok(())
}

/// Refuses further use of a codec writer whose earlier write failed partway
/// through a frame: its stream may end in a partial frame, so neither another
/// frame nor finishing can produce a valid file.
pub(crate) fn require_no_partial_frame(
    partial: bool,
    operation: TrajectoryIoOperation,
    format: TrajectoryFormat,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    if partial {
        return Err(codec_context(
            TrajectoryIoErrorKind::InvalidFrame,
            operation,
            Some(format),
            source_label,
            format!("an earlier {format} write failed partway; the stream may be incomplete"),
        ));
    }
    Ok(())
}

/// The cause of a reused reader-side validation failure without the reader's
/// operation, frame, and source context, which would misdescribe a write.
pub(crate) fn validation_detail(error: &TrajectoryIoError) -> String {
    match error.operation() {
        Some(_) if error.kind() != TrajectoryIoErrorKind::Io => error
            .detail()
            .map_or_else(|| error.kind().to_string(), str::to_owned),
        _ => error.to_string(),
    }
}

pub(crate) fn io_context(
    operation: TrajectoryIoOperation,
    format: Option<TrajectoryFormat>,
    source_label: &str,
    error: io::Error,
) -> TrajectoryIoError {
    let mut error = TrajectoryIoError::io(operation, error).with_source_label(source_label);
    if let Some(format) = format {
        error = error.with_format(format);
    }
    error
}

pub(crate) fn codec_context(
    kind: TrajectoryIoErrorKind,
    operation: TrajectoryIoOperation,
    format: Option<TrajectoryFormat>,
    source_label: &str,
    detail: impl Into<String>,
) -> TrajectoryIoError {
    TrajectoryIoError::codec(kind, operation, format)
        .with_source_label(source_label)
        .with_detail(detail)
}
