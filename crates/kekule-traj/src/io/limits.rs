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
    /// Combined retained capacities of private decode and validation buffers.
    /// Sequential and indexed access share this budget; index offsets are
    /// bounded separately by `max_index_entries` and `max_index_bytes`.
    pub max_scratch_bytes: usize,
    pub max_index_entries: usize,
    pub max_index_bytes: usize,
    pub max_text_line_bytes: usize,
    pub max_comment_bytes: usize,
    pub max_detection_bytes: usize,
}

/// Retained vector storage, including spare capacity.
pub(super) fn scratch_bytes<T>(values: &Vec<T>) -> Option<usize> {
    values.capacity().checked_mul(std::mem::size_of::<T>())
}

/// Reserves only after accounting for every other simultaneously live buffer.
/// The error factory keeps format-specific context out of the allocation logic.
pub(super) fn reserve_scratch<T>(
    values: &mut Vec<T>,
    required_capacity: usize,
    other_bytes: Option<usize>,
    limit: usize,
    error: impl Fn(&str) -> TrajectoryError,
) -> Result<(), TrajectoryError> {
    let fits = |capacity: usize| {
        capacity
            .checked_mul(std::mem::size_of::<T>())
            .and_then(|bytes| other_bytes.and_then(|other| bytes.checked_add(other)))
            .is_some_and(|bytes| bytes <= limit)
    };
    if !fits(values.capacity().max(required_capacity)) {
        return Err(error(
            "aggregate decode scratch exceeds the configured limit",
        ));
    }
    if values.capacity() < required_capacity {
        values
            .try_reserve_exact(required_capacity - values.len())
            .map_err(|_| error("could not reserve decode scratch"))?;
        // Vec permits an allocator to provide more capacity than requested.
        if !fits(values.capacity()) {
            *values = Vec::new();
            return Err(error(
                "allocated decode scratch exceeds the configured limit",
            ));
        }
    }
    Ok(())
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

#[cfg(test)]
mod scratch_tests {
    use super::*;

    fn error(detail: &str) -> TrajectoryError {
        TrajectoryCodecErrorContext::new(
            TrajectoryCodecErrorKind::ResourceLimitExceeded,
            TrajectoryIoOperation::ReadFrame,
            Some(TrajectoryFormat::Xtc),
        )
        .with_detail(detail)
        .into()
    }

    #[test]
    fn spare_capacity_counts_and_over_budget_reservations_allocate_nothing() {
        let mut retained = vec![0_u64; 16];
        retained.truncate(1);
        let required = scratch_bytes(&retained).unwrap() + 64;
        assert!(reserve_scratch(&mut retained, 1, Some(64), required - 1, error).is_err());
        reserve_scratch(&mut retained, 1, Some(64), required, error).unwrap();
        assert_eq!(retained, [0]);

        let mut empty = Vec::<u64>::new();
        assert!(reserve_scratch(&mut empty, 32, Some(64), 319, error).is_err());
        assert_eq!(empty.capacity(), 0);
        assert!(reserve_scratch(&mut empty, usize::MAX, Some(0), usize::MAX, error).is_err());
        assert_eq!(empty.capacity(), 0);
        assert!(reserve_scratch(&mut empty, 1, None, usize::MAX, error).is_err());
        assert_eq!(empty.capacity(), 0);
    }
}
