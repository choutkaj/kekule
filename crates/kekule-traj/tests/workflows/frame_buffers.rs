//! Frame buffers handing decoded frames out and taking storage back, and
//! error chains of streamed trajectory frames. Collection behavior itself is
//! covered by the `kekule` crate.

use std::error::Error;
use std::io::Cursor;
use std::sync::Arc;

use kekule::geometry::Point3;
use kekule::structure::{
    DenseArrayError, Positions, StructureError, StructureKind, TrajectoryError, TrajectoryFrame,
};
use kekule::units::{Quantity, NANOMETER};
use kekule_traj::io::xyz::{XyzReadOptions, XyzReader};
use kekule_traj::{FrameBuffer, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryReader};

use crate::support;

const TWO_FRAMES: &str = "2
first
C 0.0 1.0 2.0
H 3.0 4.0 5.0
2
second
C 1.0 2.0 3.0
H 4.0 5.0 6.0
";

#[test]
fn decoded_frames_leave_and_return_to_a_buffer_without_copying() {
    let topology = support::topology(&["C", "H"], &[(0, 1)]);
    let mut reader = XyzReader::new(
        Cursor::new(TWO_FRAMES.as_bytes()),
        Arc::clone(&topology),
        XyzReadOptions::default(),
    )
    .unwrap();
    let mut buffer = reader.frame_buffer();
    assert!(reader.read_next(&mut buffer).unwrap());
    let decoded = buffer.positions().values().value().as_ptr();

    // The frame leaves the buffer bound, with the storage it was decoded into.
    let kept = buffer.into_trajectory();
    assert_eq!(kept.len(), 1);
    assert!(Arc::ptr_eq(&kept.shared_topology(), &topology));
    let frame = kept.get(0).unwrap();
    assert_eq!(frame.positions().values().value().as_ptr(), decoded);
    assert_eq!(frame.properties().bonds().len(), topology.bond_count());
    let view = frame.as_model_view();
    assert!(matches!(view.kind(), StructureKind::Frame { .. }));
    assert!((view.positions().values().value()[1].x - 0.3).abs() < 1.0e-12);

    // Storage no longer needed becomes a buffer again; the next frame
    // decodes into the same allocation.
    let (topology, _, mut frames) = kept.into_parts();
    let mut buffer = FrameBuffer::from_frame(topology, frames.pop().unwrap()).unwrap();
    assert_eq!(buffer.positions().values().value().as_ptr(), decoded);
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_eq!(buffer.positions().values().value().as_ptr(), decoded);
    let xs = support::x_coordinates(&buffer);
    assert!((xs[0] - 0.1).abs() < 1.0e-12 && (xs[1] - 0.4).abs() < 1.0e-12);
}

#[test]
fn buffers_from_frames_validate_them_against_the_topology() {
    let topology = support::topology(&["C", "H"], &[(0, 1)]);
    // A detached frame acquires its bond rows.
    let buffer = FrameBuffer::from_frame(
        Arc::clone(&topology),
        TrajectoryFrame::new(Positions::zeros(2)),
    )
    .unwrap();
    assert_eq!(buffer.properties().bonds().len(), 1);
    assert_eq!(
        FrameBuffer::from_frame(topology, TrajectoryFrame::new(Positions::zeros(3))).unwrap_err(),
        TrajectoryIoError::from(StructureError::Array(DenseArrayError::CountMismatch {
            expected: 2,
            actual: 3,
        }))
    );
}

#[test]
fn errors_retain_nested_causes() {
    let position_error =
        Positions::new(Quantity::new([Point3::new(f64::NAN, 0.0, 0.0)], NANOMETER)).unwrap_err();
    let frame_error = StructureError::from(position_error);
    let source_text = frame_error.source().unwrap().to_string();
    let trajectory_error = TrajectoryIoError::from(frame_error.clone());
    assert_eq!(
        trajectory_error.kind(),
        TrajectoryIoErrorKind::InvalidStructure
    );
    assert_eq!(
        trajectory_error
            .source()
            .unwrap()
            .source()
            .unwrap()
            .to_string(),
        source_text
    );
    // A collection error wrapping a structure error flattens to it.
    let flattened = TrajectoryIoError::from(TrajectoryError::from(frame_error));
    assert_eq!(flattened, trajectory_error);
    let collection_error = TrajectoryIoError::from(TrajectoryError::NonMonotonicTime { frame: 1 });
    assert_eq!(
        collection_error.kind(),
        TrajectoryIoErrorKind::InvalidTrajectory
    );
    assert!(collection_error.source().unwrap().is::<TrajectoryError>());
}
