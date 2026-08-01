#![no_main]

use std::io::{Cursor, Seek};

use libfuzzer_sys::fuzz_target;
use kekule_traj::io::{detect_trajectory_format, TrajectoryFormatHint, TrajectoryIoLimits};

fuzz_target!(|data: &[u8]| {
    let mut limits = TrajectoryIoLimits::default();
    limits.max_detection_bytes = 512;
    limits.max_scratch_bytes = 512;
    let mut cursor = Cursor::new(data);
    let start = cursor.position();
    let _ = detect_trajectory_format(
        &mut cursor,
        "fuzz-input",
        TrajectoryFormatHint::Auto,
        &limits,
    );
    assert_eq!(cursor.stream_position().expect("cursor seek"), start);
});
