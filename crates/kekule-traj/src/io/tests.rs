use super::limits::{index_hard_capacity, next_index_capacity};
use super::*;
use crate::TrajectoryFormat;

#[test]
fn bounded_index_reservation_grows_logarithmically() {
    let entry_count = 100_000;
    let limits = TrajectoryIoLimits {
        max_frames: entry_count as u64,
        max_index_entries: entry_count,
        max_index_bytes: entry_count * std::mem::size_of::<u64>(),
        ..TrajectoryIoLimits::default()
    };
    let mut offsets = Vec::new();
    let mut growth_events = 0;
    for frame in 0..entry_count {
        assert_eq!(projected_index_limit(offsets.len(), &limits), None);
        let previous_capacity = offsets.capacity();
        reserve_index_for_push(
            &mut offsets,
            &limits,
            TrajectoryFormat::Xyz,
            "capacity-test.xyz",
            frame as u64,
        )
        .unwrap();
        growth_events += usize::from(offsets.capacity() != previous_capacity);
        offsets.push(frame as u64);
    }
    assert_eq!(offsets.len(), entry_count);
    assert!(
        growth_events <= 16,
        "{growth_events} growth events are not logarithmic"
    );
}

#[test]
fn index_hard_capacity_uses_the_smallest_configured_bound() {
    let limits = TrajectoryIoLimits {
        max_frames: 200,
        max_index_entries: 150,
        max_index_bytes: 125 * std::mem::size_of::<u64>(),
        ..TrajectoryIoLimits::default()
    };
    assert_eq!(index_hard_capacity(&limits), 125);
    assert_eq!(next_index_capacity(64, 64, 125), Some(125));
    assert_eq!(next_index_capacity(125, 125, 125), None);
}
