use super::*;

fn error(reason: Reason) -> TrajectoryIoError {
    TrajectoryIoError::codec(
        TrajectoryIoErrorKind::LimitExceeded,
        TrajectoryIoOperation::ReadFrame,
        Some(TrajectoryFormat::Xtc),
    )
    .with_reason(reason)
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
