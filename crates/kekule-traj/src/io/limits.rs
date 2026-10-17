use crate::{TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation};

/// Limits applied before attacker-controlled allocation, scanning, or seeking.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrajectoryIoLimits {
    max_atoms: usize,
    max_frames: u64,
    max_frame_bytes: u64,
    max_record_bytes: u64,
    max_scratch_bytes: usize,
    max_index_entries: usize,
    max_index_bytes: usize,
    max_text_line_bytes: usize,
    max_comment_bytes: usize,
    max_detection_bytes: usize,
}

impl TrajectoryIoLimits {
    pub const fn max_atoms(&self) -> usize {
        self.max_atoms
    }

    /// Sets [`Self::max_atoms`].
    #[must_use]
    pub const fn with_max_atoms(mut self, max_atoms: usize) -> Self {
        self.max_atoms = max_atoms;
        self
    }

    pub const fn max_frames(&self) -> u64 {
        self.max_frames
    }

    /// Sets [`Self::max_frames`].
    #[must_use]
    pub const fn with_max_frames(mut self, max_frames: u64) -> Self {
        self.max_frames = max_frames;
        self
    }

    pub const fn max_frame_bytes(&self) -> u64 {
        self.max_frame_bytes
    }

    /// Sets [`Self::max_frame_bytes`].
    #[must_use]
    pub const fn with_max_frame_bytes(mut self, max_frame_bytes: u64) -> Self {
        self.max_frame_bytes = max_frame_bytes;
        self
    }

    pub const fn max_record_bytes(&self) -> u64 {
        self.max_record_bytes
    }

    /// Sets [`Self::max_record_bytes`].
    #[must_use]
    pub const fn with_max_record_bytes(mut self, max_record_bytes: u64) -> Self {
        self.max_record_bytes = max_record_bytes;
        self
    }

    /// Combined retained capacities of private decode and validation buffers.
    /// Sequential and indexed access share this budget; index offsets are
    /// bounded separately by `max_index_entries` and `max_index_bytes`.
    pub const fn max_scratch_bytes(&self) -> usize {
        self.max_scratch_bytes
    }

    /// Sets [`Self::max_scratch_bytes`].
    #[must_use]
    pub const fn with_max_scratch_bytes(mut self, max_scratch_bytes: usize) -> Self {
        self.max_scratch_bytes = max_scratch_bytes;
        self
    }

    pub const fn max_index_entries(&self) -> usize {
        self.max_index_entries
    }

    /// Sets [`Self::max_index_entries`].
    #[must_use]
    pub const fn with_max_index_entries(mut self, max_index_entries: usize) -> Self {
        self.max_index_entries = max_index_entries;
        self
    }

    pub const fn max_index_bytes(&self) -> usize {
        self.max_index_bytes
    }

    /// Sets [`Self::max_index_bytes`].
    #[must_use]
    pub const fn with_max_index_bytes(mut self, max_index_bytes: usize) -> Self {
        self.max_index_bytes = max_index_bytes;
        self
    }

    pub const fn max_text_line_bytes(&self) -> usize {
        self.max_text_line_bytes
    }

    /// Sets [`Self::max_text_line_bytes`].
    #[must_use]
    pub const fn with_max_text_line_bytes(mut self, max_text_line_bytes: usize) -> Self {
        self.max_text_line_bytes = max_text_line_bytes;
        self
    }

    pub const fn max_comment_bytes(&self) -> usize {
        self.max_comment_bytes
    }

    /// Sets [`Self::max_comment_bytes`].
    #[must_use]
    pub const fn with_max_comment_bytes(mut self, max_comment_bytes: usize) -> Self {
        self.max_comment_bytes = max_comment_bytes;
        self
    }

    pub const fn max_detection_bytes(&self) -> usize {
        self.max_detection_bytes
    }

    /// Sets [`Self::max_detection_bytes`].
    #[must_use]
    pub const fn with_max_detection_bytes(mut self, max_detection_bytes: usize) -> Self {
        self.max_detection_bytes = max_detection_bytes;
        self
    }
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
    error: impl Fn(&str) -> TrajectoryIoError,
) -> Result<(), TrajectoryIoError> {
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
) -> Result<(), TrajectoryIoError> {
    if offsets.len() < offsets.capacity() {
        return Ok(());
    }
    let hard_capacity = index_hard_capacity(limits);
    let Some(target_capacity) =
        next_index_capacity(offsets.len(), offsets.capacity(), hard_capacity)
    else {
        return Err(TrajectoryIoError::codec(
            TrajectoryIoErrorKind::LimitExceeded,
            TrajectoryIoOperation::Index,
            Some(format),
        )
        .with_source_label(source_label)
        .with_frame(frame)
        .with_detail(format!(
            "{format} index reached its configured hard capacity"
        )));
    };
    offsets
        .try_reserve_exact(target_capacity.saturating_sub(offsets.len()))
        .map_err(|_| {
            TrajectoryIoError::codec(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::Index,
                Some(format),
            )
            .with_source_label(source_label)
            .with_frame(frame)
            .with_detail(format!(
                "could not grow {format} index toward its bounded capacity"
            ))
        })
}

#[cfg(test)]
mod scratch_tests {
    use super::*;

    fn error(detail: &str) -> TrajectoryIoError {
        TrajectoryIoError::codec(
            TrajectoryIoErrorKind::LimitExceeded,
            TrajectoryIoOperation::ReadFrame,
            Some(TrajectoryFormat::Xtc),
        )
        .with_detail(detail)
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
