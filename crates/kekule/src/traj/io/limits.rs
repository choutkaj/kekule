use crate::diagnostics::{LimitExceeded, Resource};
use crate::traj::{
    TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation,
};

/// Limits applied before attacker-controlled allocation, scanning, or seeking.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TrajectoryIoLimits {
    max_atoms: usize,
    max_frames: usize,
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
    /// Largest atom count a file header or frame may declare; default 10 000 000.
    pub const fn max_atoms(&self) -> usize {
        self.max_atoms
    }

    /// Sets [`Self::max_atoms`].
    #[must_use]
    pub const fn with_max_atoms(mut self, max_atoms: usize) -> Self {
        self.max_atoms = max_atoms;
        self
    }

    /// Largest number of frames read or indexed from one stream; default
    /// 100 000 000.
    pub const fn max_frames(&self) -> usize {
        self.max_frames
    }

    /// Sets [`Self::max_frames`].
    #[must_use]
    pub const fn with_max_frames(mut self, max_frames: usize) -> Self {
        self.max_frames = max_frames;
        self
    }

    /// Largest encoded size of one frame in bytes; default 4 GiB.
    pub const fn max_frame_bytes(&self) -> u64 {
        self.max_frame_bytes
    }

    /// Sets [`Self::max_frame_bytes`].
    #[must_use]
    pub const fn with_max_frame_bytes(mut self, max_frame_bytes: u64) -> Self {
        self.max_frame_bytes = max_frame_bytes;
        self
    }

    /// Largest declared size of one header or data record in bytes; default
    /// 4 GiB.
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
    /// A frame's dense arrays count as one value per atom even when a read
    /// allocates them anew, so a frame passes or fails the same way whatever
    /// the reader read before. Sequential and indexed access share this
    /// budget; index offsets are bounded separately by `max_index_entries`
    /// and `max_index_bytes`.
    pub const fn max_scratch_bytes(&self) -> usize {
        self.max_scratch_bytes
    }

    /// Sets [`Self::max_scratch_bytes`].
    #[must_use]
    pub const fn with_max_scratch_bytes(mut self, max_scratch_bytes: usize) -> Self {
        self.max_scratch_bytes = max_scratch_bytes;
        self
    }

    /// Largest number of frame offsets an indexed reader stores; default
    /// 100 000 000.
    pub const fn max_index_entries(&self) -> usize {
        self.max_index_entries
    }

    /// Sets [`Self::max_index_entries`].
    #[must_use]
    pub const fn with_max_index_entries(mut self, max_index_entries: usize) -> Self {
        self.max_index_entries = max_index_entries;
        self
    }

    /// Largest memory in bytes an indexed reader's offset table may occupy;
    /// default 800 000 000. Frame times and steps, where the format records
    /// them, are tables of the same length and size.
    pub const fn max_index_bytes(&self) -> usize {
        self.max_index_bytes
    }

    /// Sets [`Self::max_index_bytes`].
    #[must_use]
    pub const fn with_max_index_bytes(mut self, max_index_bytes: usize) -> Self {
        self.max_index_bytes = max_index_bytes;
        self
    }

    /// Longest text line, such as an XYZ atom line or TRR version string, in
    /// bytes; default 1 MiB.
    pub const fn max_text_line_bytes(&self) -> usize {
        self.max_text_line_bytes
    }

    /// Sets [`Self::max_text_line_bytes`].
    #[must_use]
    pub const fn with_max_text_line_bytes(mut self, max_text_line_bytes: usize) -> Self {
        self.max_text_line_bytes = max_text_line_bytes;
        self
    }

    /// Longest XYZ comment line in bytes; default 1 MiB.
    pub const fn max_comment_bytes(&self) -> usize {
        self.max_comment_bytes
    }

    /// Sets [`Self::max_comment_bytes`].
    #[must_use]
    pub const fn with_max_comment_bytes(mut self, max_comment_bytes: usize) -> Self {
        self.max_comment_bytes = max_comment_bytes;
        self
    }

    /// Bytes inspected when detecting a format from content; default 4096.
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

/// The storage of a dense array once a read fills it with `len` values: its
/// capacity, or `len` values when the array is smaller, as it is after the
/// reader gave its storage to a sample. A limit counted this way does not
/// depend on what the reader read before.
pub(super) fn filled_bytes<T>(values: &Vec<T>, len: usize) -> Option<usize> {
    values
        .capacity()
        .max(len)
        .checked_mul(std::mem::size_of::<T>())
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
        return Err(error(&format!(
            "The buffers would exceed their limit of {limit} bytes."
        )));
    }
    if values.capacity() < required_capacity {
        values
            .try_reserve_exact(required_capacity - values.len())
            .map_err(|_| error("Cannot allocate the buffers."))?;
        // Vec permits an allocator to provide more capacity than requested.
        if !fits(values.capacity()) {
            *values = Vec::new();
            return Err(error(&format!(
                "The allocated buffers exceed their limit of {limit} bytes."
            )));
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

/// The limit a frame index of `current_entries` frames exceeds by indexing
/// one more frame: the frame limit, the index entry limit, or the frames the
/// index byte limit holds.
pub(crate) fn projected_index_limit(
    current_entries: usize,
    limits: &TrajectoryIoLimits,
) -> Option<LimitExceeded> {
    let projected_entries = current_entries.saturating_add(1);
    if current_entries >= limits.max_frames {
        return Some(LimitExceeded::new(
            Resource::Frames,
            limits.max_frames,
            projected_entries,
        ));
    }
    if projected_entries > limits.max_index_entries || current_entries == usize::MAX {
        return Some(LimitExceeded::new(
            Resource::Frames,
            limits.max_index_entries,
            projected_entries,
        ));
    }
    let byte_capacity = limits.max_index_bytes / std::mem::size_of::<u64>();
    if projected_entries > byte_capacity {
        return Some(LimitExceeded::new(
            Resource::Frames,
            byte_capacity,
            projected_entries,
        ));
    }
    None
}

/// A byte limit of `limit` exceeded by `observed` bytes, saturating where
/// the platform's address range is smaller.
pub(crate) fn exceeded_bytes(limit: u64, observed: u64) -> LimitExceeded {
    let size = |bytes: u64| usize::try_from(bytes).unwrap_or(usize::MAX);
    LimitExceeded::new(Resource::RecordBytes, size(limit), size(observed))
}

pub(super) fn index_hard_capacity(limits: &TrajectoryIoLimits) -> usize {
    limits
        .max_frames
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
    frame: usize,
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
            "The frame index is full at {hard_capacity} frames."
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
            .with_detail("Cannot allocate the frame index.")
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
