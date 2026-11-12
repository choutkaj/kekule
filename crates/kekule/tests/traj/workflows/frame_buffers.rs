//! Frame buffers handing decoded frames out and taking storage back, frames
//! reaching a writer without a copy, and error chains of streamed
//! trajectory frames. Collection behavior itself is covered by the `kekule`
//! crate.

use std::error::Error;
use std::io::Cursor;
use std::sync::Arc;

use kekule::geometry::Point3;
use kekule::structure::{DenseArrayError, Place, Positions, Sample, StructureError, StructureKind};
use kekule::traj::io::xyz::{XyzReadOptions, XyzReader};
use kekule::traj::{
    FrameBuffer, FrameView, MemoryTrajectoryWriter, TrajectoryIoError, TrajectoryIoErrorKind,
    TrajectoryReader, TrajectoryWriter,
};
use kekule::units::{Quantity, NANOMETER};

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
    let (sample, place) = buffer.into_sample();
    assert_eq!(
        place,
        Place::Time {
            time: None,
            step: None
        }
    );
    assert_eq!(sample.positions().values().value().as_ptr(), decoded);
    assert_eq!(sample.properties().bonds().len(), topology.bond_count());
    let frame = FrameView::new(&topology, &sample, None, None).unwrap();
    let view = frame.as_model_view();
    assert!(matches!(
        view.kind(),
        StructureKind::Sample(context)
            if matches!(context.place(), Some(Place::Time { .. }))
    ));
    assert!((view.positions().values().value()[1].x - 0.3).abs() < 1.0e-12);

    // Storage no longer needed becomes a buffer again; the next frame
    // decodes into the same allocation.
    let mut buffer = FrameBuffer::from_sample(Arc::clone(&topology), sample).unwrap();
    assert_eq!(buffer.positions().values().value().as_ptr(), decoded);
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_eq!(buffer.positions().values().value().as_ptr(), decoded);
    let xs = support::x_coordinates(&buffer);
    assert!((xs[0] - 0.1).abs() < 1.0e-12 && (xs[1] - 0.4).abs() < 1.0e-12);
}

/// A frame written to memory, from a buffer or a set, is collected as the
/// sample it shows: the writer path copies no frame.
#[test]
fn written_frames_reach_a_memory_writer_without_copying() {
    let topology = support::topology(&["C", "H"], &[(0, 1)]);
    let mut reader = XyzReader::new(
        Cursor::new(TWO_FRAMES.as_bytes()),
        Arc::clone(&topology),
        XyzReadOptions::default(),
    )
    .unwrap();
    let mut buffer = reader.frame_buffer();
    let mut writer = MemoryTrajectoryWriter::new(Arc::clone(&topology));
    let mut decoded = Vec::new();
    while reader.read_next(&mut buffer).unwrap() {
        writer.write_frame(buffer.frame_view()).unwrap();
        decoded.push(Sample::clone(&buffer));
    }
    let set = writer.into_sample_set();
    assert_eq!(set.len(), 2);
    for (sample, decoded) in set.iter().zip(&decoded) {
        assert!(sample.payload().shares_structure_with(decoded));
    }

    // A frame of a set shows the set's own sample, and writing the set
    // again shares every sample with it.
    let mut writer = MemoryTrajectoryWriter::new(Arc::clone(&topology));
    for sample in &set {
        let frame = FrameView::from_sample(sample).unwrap();
        assert!(std::ptr::eq(frame.payload(), sample.payload()));
        writer.write_frame(frame).unwrap();
    }
    let copy = writer.into_sample_set();
    for (copied, sample) in copy.iter().zip(&set) {
        assert!(copied.payload().shares_structure_with(sample.payload()));
    }
}

#[test]
fn buffers_from_samples_validate_them_against_the_topology() {
    let topology = support::topology(&["C", "H"], &[(0, 1)]);
    // A detached sample acquires its bond rows.
    let buffer =
        FrameBuffer::from_sample(Arc::clone(&topology), Sample::new(Positions::zeros(2))).unwrap();
    assert_eq!(buffer.properties().bonds().len(), 1);
    assert_eq!(
        FrameBuffer::from_sample(topology, Sample::new(Positions::zeros(3))).unwrap_err(),
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
    let trajectory_error = TrajectoryIoError::from(frame_error.clone());
    assert_eq!(
        trajectory_error.kind(),
        TrajectoryIoErrorKind::InvalidStructure
    );
    assert_eq!(
        trajectory_error
            .source()
            .unwrap()
            .downcast_ref::<StructureError>()
            .map(ToString::to_string),
        Some(frame_error.to_string())
    );
}
