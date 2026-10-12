use crate::{
    TrajectoryCodecErrorContext, TrajectoryCodecErrorKind, TrajectoryFormat, TrajectoryIoError,
    TrajectoryIoErrorContext, TrajectoryIoOperation,
};
use std::io;

pub(crate) fn frame_offset_context(
    error: TrajectoryIoError,
    frame: u64,
    byte_offset: u64,
) -> TrajectoryIoError {
    match error {
        TrajectoryIoError::Io(context) => {
            let mut context = *context;
            if context.frame().is_none() {
                context = context.with_frame(frame);
            }
            if context.byte_offset().is_none() {
                context = context.with_byte_offset(byte_offset);
            }
            TrajectoryIoError::Io(Box::new(context))
        }
        TrajectoryIoError::Codec(context) => {
            let mut context = *context;
            if context.frame().is_none() {
                context = context.with_frame(frame);
            }
            if context.byte_offset().is_none() {
                context = context.with_byte_offset(byte_offset);
            }
            TrajectoryIoError::Codec(Box::new(context))
        }
        error => error,
    }
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
            TrajectoryCodecErrorKind::InvalidFrame,
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
            TrajectoryCodecErrorKind::InvalidFrame,
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
    match error {
        TrajectoryIoError::Codec(context) => context
            .detail()
            .map_or_else(|| context.kind().to_string(), str::to_owned),
        other => other.to_string(),
    }
}

pub(crate) fn io_context(
    operation: TrajectoryIoOperation,
    format: Option<TrajectoryFormat>,
    source_label: &str,
    error: io::Error,
) -> TrajectoryIoError {
    let mut context = TrajectoryIoErrorContext::new(operation, error.kind(), error.to_string())
        .with_source_label(source_label);
    if let Some(format) = format {
        context = context.with_format(format);
    }
    context.into()
}

pub(crate) fn codec_context(
    kind: TrajectoryCodecErrorKind,
    operation: TrajectoryIoOperation,
    format: Option<TrajectoryFormat>,
    source_label: &str,
    detail: impl Into<String>,
) -> TrajectoryIoError {
    TrajectoryCodecErrorContext::new(kind, operation, format)
        .with_source_label(source_label)
        .with_detail(detail)
        .into()
}
