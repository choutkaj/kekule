use std::fs;
use std::path::PathBuf;
use std::sync::Arc;
use std::time::{SystemTime, UNIX_EPOCH};

use kekule::topology::Topology;
use kekule_traj::io::{
    open_indexed_trajectory, open_trajectory, FieldAvailability, RandomAccessCapability,
    TrajectoryFormatHint, TrajectoryIoLimits, TrajectoryOpenOptions,
};
use kekule_traj::{
    SeekableTrajectoryReader, TrajectoryFormat, TrajectoryIoErrorKind, TrajectoryReader,
};

use crate::support;
use support::topology as build_topology;

fn topology() -> Arc<Topology> {
    build_topology(&["C"], &[])
}

#[test]
fn trajectory_payload_and_buffer_reject_mismatched_annotation_rows() {
    use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
    use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
    use kekule_traj::{FrameBuffer, FrameData};

    let topology = topology();
    let key = PropertyKey::new("tag").unwrap();
    let positions = Positions::zeros(1);
    // Annotations of a larger structure do not fit this topology's rows.
    let mut larger = TrajectoryFrame::new(Positions::zeros(2));
    larger
        .properties_mut()
        .atoms_mut()
        .insert(key.clone(), PropertyColumn::Int(vec![Some(7), Some(8)]))
        .unwrap();
    let foreign = larger.properties().clone();
    let mut frame = TrajectoryFrame::new(positions.clone());
    let original = frame.properties().clone();
    assert!(frame.set_properties(foreign.clone()).is_err());
    assert_eq!(frame.properties(), &original);
    let mut buffer = FrameBuffer::new(topology);
    assert!(buffer.frame_mut().set_properties(foreign).is_err());
    assert!(buffer
        .replace_from_data(FrameData::from_frame_view(
            Trajectory::from_items(kekule::smiles::to_topology("CC").unwrap(), [larger])
                .unwrap()
                .get(0)
                .unwrap()
        ))
        .is_err());
    assert_eq!(buffer.properties(), &original);
    assert_eq!(buffer.positions(), &positions);
    // Owner annotations publish with the decoded frame.
    let mut owner = kekule::properties::OwnerProperties::new();
    owner.insert(key.clone(), PropertyValue::Int(3)).unwrap();
    buffer
        .replace_from_data(FrameData::new(positions.values()).with_owner_properties(&owner))
        .unwrap();
    assert_eq!(buffer.properties().owner(), &owner);
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

    let mut sequential = open_trajectory(&path, Arc::clone(&topology), options.clone()).unwrap();
    let report = sequential.open_report();
    assert_eq!(report.selected_format(), TrajectoryFormat::Xyz);
    assert_eq!(
        sequential.metadata().fields().positions(),
        FieldAvailability::Required
    );
    assert_eq!(
        sequential.metadata().random_access(),
        RandomAccessCapability::SequentialOnly
    );
    let mut destination = sequential.frame_buffer();
    assert!(Arc::ptr_eq(&topology, &destination.shared_topology()));
    assert!(sequential.read_next(&mut destination).unwrap());
    assert!(!sequential.read_next(&mut destination).unwrap());

    let mut indexed = open_indexed_trajectory(&path, Arc::clone(&topology), options).unwrap();
    assert_eq!(indexed.frame_count(), Some(1));
    assert_eq!(
        indexed.metadata().random_access(),
        RandomAccessCapability::Indexed
    );
    indexed.read_frame(0, &mut destination).unwrap();
    assert!((destination.positions().values().value()[0].x - 0.1).abs() < 1.0e-15);

    let owned_topology = Arc::try_unwrap(build_topology(&["C"], &[])).unwrap();
    let mut default_sequential =
        open_trajectory(&path, owned_topology, Default::default()).unwrap();
    assert_eq!(
        default_sequential.open_report().selected_format(),
        TrajectoryFormat::Xyz
    );
    let mut owned_buffer = default_sequential.frame_buffer();
    assert!(Arc::ptr_eq(
        &default_sequential.shared_topology(),
        &owned_buffer.shared_topology()
    ));
    // A separately owned topology remains a different buffer context, even at
    // the same atom count. Rejection must not consume the pending frame.
    assert!(matches!(
        default_sequential.read_next(&mut destination),
        Err(error) if error.kind() == TrajectoryIoErrorKind::TopologyMismatch
    ));
    assert!(default_sequential.read_next(&mut owned_buffer).unwrap());
    assert_eq!(
        owned_buffer.positions().values(),
        destination.positions().values()
    );
    assert!(!default_sequential.read_next(&mut owned_buffer).unwrap());

    let mut default_indexed =
        open_indexed_trajectory(&path, Arc::clone(&topology), Default::default()).unwrap();
    assert_eq!(default_indexed.frame_count(), Some(1));
    assert_eq!(
        default_indexed.open_report().selected_format(),
        TrajectoryFormat::Xyz
    );
    default_indexed.read_frame(0, &mut destination).unwrap();

    fs::remove_file(path).unwrap();
}

#[test]
fn path_options_enforce_global_limits_and_path_diagnostics() {
    use kekule_traj::io::xyz::XyzReadOptions;

    let path = temporary_xyz();
    fs::write(&path, b"1\nlimits\nC 1 2 3\n").unwrap();
    let topology = topology();
    let options = TrajectoryOpenOptions::default()
        .with_limits(TrajectoryIoLimits::default().with_max_atoms(0))
        .with_xyz_options(XyzReadOptions::default().with_source_label("custom stream"));
    let errors = [
        open_trajectory(&path, Arc::clone(&topology), options.clone())
            .err()
            .unwrap(),
        open_indexed_trajectory(&path, topology, options)
            .err()
            .unwrap(),
    ];
    for error in errors {
        assert_eq!(
            support::codec_kind(&error),
            Some(TrajectoryIoErrorKind::LimitExceeded)
        );
        let context = error;
        assert_eq!(context.source_label(), Some(path.to_str().unwrap()));
    }
    fs::remove_file(path).unwrap();
}
