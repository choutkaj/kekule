use crate::{
    TrajectoryCodecErrorContext, TrajectoryCodecErrorKind, TrajectoryError, TrajectoryFormat,
    TrajectoryIoErrorContext, TrajectoryIoOperation,
};
use std::io;

pub(crate) fn frame_offset_context(
    error: TrajectoryError,
    frame: u64,
    byte_offset: u64,
) -> TrajectoryError {
    match error {
        TrajectoryError::Io(context) => {
            let mut context = *context;
            if context.frame().is_none() {
                context = context.with_frame(frame);
            }
            if context.byte_offset().is_none() {
                context = context.with_byte_offset(byte_offset);
            }
            TrajectoryError::Io(Box::new(context))
        }
        TrajectoryError::Codec(context) => {
            let mut context = *context;
            if context.frame().is_none() {
                context = context.with_frame(frame);
            }
            if context.byte_offset().is_none() {
                context = context.with_byte_offset(byte_offset);
            }
            TrajectoryError::Codec(Box::new(context))
        }
        error => error,
    }
}

pub(crate) fn probe_seekable_eof<R: io::Read + io::Seek>(
    reader: &mut R,
    operation: TrajectoryIoOperation,
    format: TrajectoryFormat,
    source_label: &str,
) -> Result<bool, TrajectoryError> {
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
) -> Result<(), TrajectoryError> {
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

pub(crate) fn io_context(
    operation: TrajectoryIoOperation,
    format: Option<TrajectoryFormat>,
    source_label: &str,
    error: io::Error,
) -> TrajectoryError {
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
) -> TrajectoryError {
    TrajectoryCodecErrorContext::new(kind, operation, format)
        .with_source_label(source_label)
        .with_detail(detail)
        .into()
}
