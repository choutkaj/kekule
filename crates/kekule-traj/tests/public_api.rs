use std::fs;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::{SystemTime, UNIX_EPOCH};

use kekule::topology::Topology;
use kekule_traj::io::{
    open_indexed_trajectory, open_trajectory, FieldAvailability, RandomAccessCapability,
    TrajectoryFormatHint, TrajectoryOpenOptions,
};
use kekule_traj::{FrameBuffer, SeekableTrajectoryReader, TrajectoryFormat, TrajectoryReader};

mod support;
use support::{binding, topology as build_topology};

fn topology() -> Arc<Topology> {
    build_topology(&["C"], &[])
}

fn temporary_xyz() -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    std::env::temp_dir().join(format!(
        "kekule-trajectory-public-api-{}-{nonce}.xyz",
        std::process::id()
    ))
}

#[test]
fn format_agnostic_public_api_opens_sequential_and_indexed_readers() {
    let path = temporary_xyz();
    fs::write(&path, b"1\npublic API\nC 1.0 2.0 3.0\n").unwrap();
    let topology = topology();
    let options = TrajectoryOpenOptions::default()
        .with_format_hint(TrajectoryFormatHint::Explicit(TrajectoryFormat::Xyz));

    let (mut sequential, report) =
        open_trajectory(&path, binding(&topology), options.clone()).unwrap();
    assert_eq!(report.selected_format(), TrajectoryFormat::Xyz);
    assert_eq!(
        sequential.metadata().fields().positions,
        FieldAvailability::Required
    );
    assert_eq!(
        sequential.metadata().random_access(),
        RandomAccessCapability::SequentialOnly
    );
    let mut destination = FrameBuffer::new(Arc::clone(&topology));
    assert!(sequential.read_next(&mut destination).unwrap());
    assert!(!sequential.read_next(&mut destination).unwrap());

    let (mut indexed, _) = open_indexed_trajectory(&path, binding(&topology), options).unwrap();
    assert_eq!(indexed.frame_count(), Some(1));
    assert_eq!(
        indexed.metadata().random_access(),
        RandomAccessCapability::Indexed
    );
    indexed.read_frame(0, &mut destination).unwrap();
    assert!((destination.positions().values().value()[0].x - 0.1).abs() < 1.0e-15);

    fs::remove_file(path).unwrap();
}
