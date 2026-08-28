use crate::{
    TrajectoryCodecErrorContext, TrajectoryCodecErrorKind, TrajectoryError, TrajectoryFormat,
    TrajectoryIoOperation,
};

/// Limits applied before attacker-controlled allocation, scanning, or seeking.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrajectoryIoLimits {
    pub max_atoms: usize,
    pub max_frames: u64,
    pub max_frame_bytes: u64,
    pub max_record_bytes: u64,
    pub max_scratch_bytes: usize,
    pub max_index_entries: usize,
    pub max_index_bytes: usize,
    pub max_text_line_bytes: usize,
    pub max_comment_bytes: usize,
    pub max_detection_bytes: usize,
}

impl Default for TrajectoryIoLimits {
    fn default() -> Self {
        Self {
            max_atoms: 10_000_000,
            max_frames: 100_000_000,
            max_frame_bytes: 4 * 1024 * 1024 * 1024,
            max_record_bytes: 4 * 1024 * 1024 * 1024,
            max_scratch_bytes: usize::try_from(4_u64 * 1024 * 1024 * 1024).unwrap_or(usize::MAX),
            max_index_entries: 100_000_000,
            max_index_bytes: 800_000_000,
            max_text_line_bytes: 1_048_576,
            max_comment_bytes: 1_048_576,
            max_detection_bytes: 4096,
        }
    }
}

pub(crate) fn projected_index_limit(
    current_entries: usize,
    limits: &TrajectoryIoLimits,
) -> Option<&'static str> {
    if u64::try_from(current_entries).map_or(true, |entries| entries >= limits.max_frames) {
        return Some("frame count");
    }
    let Some(projected_entries) = current_entries.checked_add(1) else {
        return Some("entry count");
    };
    if projected_entries > limits.max_index_entries {
        return Some("entry count");
    }
    if projected_entries
        .checked_mul(std::mem::size_of::<u64>())
        .is_none_or(|bytes| bytes > limits.max_index_bytes)
    {
        return Some("byte count");
    }
    None
}

pub(super) fn index_hard_capacity(limits: &TrajectoryIoLimits) -> usize {
    usize::try_from(limits.max_frames)
        .unwrap_or(usize::MAX)
        .min(limits.max_index_entries)
        .min(limits.max_index_bytes / std::mem::size_of::<u64>())
}

pub(super) fn next_index_capacity(
    current_entries: usize,
    current_capacity: usize,
    hard_capacity: usize,
) -> Option<usize> {
    if current_entries < current_capacity {
        return None;
    }
    let minimum = current_entries.checked_add(1)?;
    if minimum > hard_capacity {
        return None;
    }
    let geometric = if current_capacity == 0 {
        8
    } else {
        current_capacity.saturating_mul(2)
    };
    Some(geometric.max(minimum).min(hard_capacity))
}

pub(crate) fn reserve_index_for_push(
    offsets: &mut Vec<u64>,
    limits: &TrajectoryIoLimits,
    format: TrajectoryFormat,
    source_label: &str,
    frame: u64,
) -> Result<(), TrajectoryError> {
    if offsets.len() < offsets.capacity() {
        return Ok(());
    }
    let hard_capacity = index_hard_capacity(limits);
    let Some(target_capacity) =
        next_index_capacity(offsets.len(), offsets.capacity(), hard_capacity)
    else {
        return Err(TrajectoryCodecErrorContext::new(
            TrajectoryCodecErrorKind::ResourceLimitExceeded,
            TrajectoryIoOperation::Index,
            Some(format),
        )
        .with_source_label(source_label)
        .with_frame(frame)
        .with_detail(format!(
            "{format} index reached its configured hard capacity"
        ))
        .into());
    };
    offsets
        .try_reserve_exact(target_capacity.saturating_sub(offsets.len()))
        .map_err(|_| {
            TrajectoryCodecErrorContext::new(
                TrajectoryCodecErrorKind::ResourceLimitExceeded,
                TrajectoryIoOperation::Index,
                Some(format),
            )
            .with_source_label(source_label)
            .with_frame(frame)
            .with_detail(format!(
                "could not grow {format} index toward its bounded capacity"
            ))
            .into()
        })
}
