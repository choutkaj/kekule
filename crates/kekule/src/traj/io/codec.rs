//! Errors and checks every codec shares; each takes the codec's format.

use std::io::{self, Read};

use crate::diagnostics::{LimitExceeded, Resource};
use crate::traj::{
    TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation,
};

use super::{codec_context, io_context, TrajectoryIoLimits};

/// The operation of a read of `frame`, or of the header before any frame.
pub(super) const fn read_operation(frame: Option<usize>) -> TrajectoryIoOperation {
    match frame {
        Some(_) => TrajectoryIoOperation::ReadFrame,
        None => TrajectoryIoOperation::ReadHeader,
    }
}

/// A malformed header.
pub(super) fn header_error(
    format: TrajectoryFormat,
    source_label: &str,
    detail: impl Into<String>,
) -> TrajectoryIoError {
    codec_context(
        TrajectoryIoErrorKind::InvalidHeader,
        TrajectoryIoOperation::ReadHeader,
        Some(format),
        source_label,
        detail,
    )
}

/// A failure of `kind` in frame `frame`.
pub(super) fn frame_error(
    format: TrajectoryFormat,
    kind: TrajectoryIoErrorKind,
    operation: TrajectoryIoOperation,
    source_label: &str,
    frame: usize,
    detail: impl Into<String>,
) -> TrajectoryIoError {
    TrajectoryIoError::codec(kind, operation, Some(format))
        .with_source_label(source_label)
        .with_frame(frame)
        .with_detail(detail)
}

/// A failure of `kind` in frame `frame`, when there is one.
pub(super) fn stream_error(
    format: TrajectoryFormat,
    kind: TrajectoryIoErrorKind,
    operation: TrajectoryIoOperation,
    source_label: &str,
    frame: Option<usize>,
    detail: impl Into<String>,
) -> TrajectoryIoError {
    let error = codec_context(kind, operation, Some(format), source_label, detail);
    match frame {
        Some(frame) => error.with_frame(frame),
        None => error,
    }
}

/// A limit or range a read exceeds, in frame `frame` when there is one.
pub(super) fn limit_error(
    format: TrajectoryFormat,
    operation: TrajectoryIoOperation,
    source_label: &str,
    frame: Option<usize>,
    detail: impl Into<String>,
) -> TrajectoryIoError {
    stream_error(
        format,
        TrajectoryIoErrorKind::LimitExceeded,
        operation,
        source_label,
        frame,
        detail,
    )
}

/// A frame a writer cannot write within the format's or the platform's
/// ranges.
pub(super) fn writer_limit(
    format: TrajectoryFormat,
    source_label: &str,
    detail: impl Into<String>,
) -> TrajectoryIoError {
    codec_context(
        TrajectoryIoErrorKind::LimitExceeded,
        TrajectoryIoOperation::WriteFrame,
        Some(format),
        source_label,
        detail,
    )
}

/// Fails unless a header's `count` atoms are at least one and at most the
/// limit, or `ceiling` where the format holds fewer.
pub(super) fn validate_atom_count(
    format: TrajectoryFormat,
    count: usize,
    ceiling: usize,
    limits: &TrajectoryIoLimits,
    source_label: &str,
) -> Result<(), TrajectoryIoError> {
    let error = |detail| {
        limit_error(
            format,
            TrajectoryIoOperation::ReadHeader,
            source_label,
            None,
            detail,
        )
    };
    if count == 0 {
        return Err(error("The header declares no atoms."));
    }
    let limit = limits.max_atoms().min(ceiling);
    if count > limit {
        return Err(error("").exceeding(
            "The atom count",
            LimitExceeded::new(Resource::Atoms, limit, count),
        ));
    }
    Ok(())
}

/// Fills `bytes`; a stream that ends first is a truncated record described
/// by `detail`.
pub(super) fn read_exact(
    format: TrajectoryFormat,
    reader: &mut impl Read,
    bytes: &mut [u8],
    operation: TrajectoryIoOperation,
    source_label: &str,
    frame: Option<usize>,
    detail: &str,
) -> Result<(), TrajectoryIoError> {
    reader.read_exact(bytes).map_err(|error| {
        if error.kind() == io::ErrorKind::UnexpectedEof {
            stream_error(
                format,
                TrajectoryIoErrorKind::TruncatedRecord,
                operation,
                source_label,
                frame,
                detail,
            )
        } else {
            io_context(operation, Some(format), source_label, error)
        }
    })
}
