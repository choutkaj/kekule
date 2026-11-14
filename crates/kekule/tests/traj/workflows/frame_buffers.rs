//! Frame buffers handing decoded frames out and taking storage back,
//! keeping no stale occupancies or B-factors, frames reaching a writer
//! without a copy, and error chains of streamed trajectory frames. Collection behavior itself is covered by the `kekule`
//! crate.

use kekule::diagnostics::Countable;
use kekule::diagnostics::LengthMismatch;
use std::error::Error;
use std::io::Cursor;
use std::sync::Arc;

use kekule::math::Point3;
use kekule::structure::{
    DenseArrayError, Geometry, GeometryError, Place, Positions, StructureKind, StructureSetError,
};
use kekule::traj::io::xyz::{XyzReadOptions, XyzReader};
use kekule::traj::{
    FrameBuffer, FrameView, MemoryTrajectoryWriter, TrajectoryIoError, TrajectoryIoErrorKind,
    TrajectoryReader, TrajectoryWriter,
};
use kekule::units::{Quantity, NANOMETER, PICOSECOND, SQUARE_NANOMETER};

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
    let (geometry, place) = buffer.into_geometry();
    assert_eq!(
        place,
        Place::Time {
            time: None,
            step: None
        }
    );
    assert_eq!(geometry.positions().values().value().as_ptr(), decoded);
    assert_eq!(geometry.properties().bonds().len(), topology.bond_count());
    let frame = FrameView::new(&topology, &geometry, None, None).unwrap();
    let view = frame.as_structure_view();
    assert!(matches!(
        view.kind(),
        StructureKind::InSet(context)
            if matches!(context.place(), Some(Place::Time { .. }))
    ));
    assert!((view.positions().values().value()[1].x - 0.3).abs() < 1.0e-12);

    // Storage no longer needed becomes a buffer again; the next frame
    // decodes into the same allocation.
    let mut buffer = FrameBuffer::from_geometry(Arc::clone(&topology), geometry).unwrap();
    assert_eq!(buffer.positions().values().value().as_ptr(), decoded);
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_eq!(buffer.positions().values().value().as_ptr(), decoded);
    let xs = support::x_coordinates(&buffer);
    assert!((xs[0] - 0.1).abs() < 1.0e-12 && (xs[1] - 0.4).abs() < 1.0e-12);
}

/// A frame written to memory, from a buffer or a set, is collected as the
/// geometry it shows: the writer path copies no frame.
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
        decoded.push(Geometry::clone(&buffer));
    }
    let set = writer.into_structure_set();
    assert_eq!(set.len(), 2);
    for (geometry, decoded) in set.iter().zip(&decoded) {
        assert!(geometry.payload().shares_storage_with(decoded));
    }

    // A frame of a set shows the set's own geometry, and writing the set
    // again shares every geometry with it.
    let mut writer = MemoryTrajectoryWriter::new(Arc::clone(&topology));
    for geometry in &set {
        let frame = FrameView::from_geometry(geometry).unwrap();
        assert!(std::ptr::eq(frame.payload(), geometry.payload()));
        writer.write_frame(frame).unwrap();
    }
    let copy = writer.into_structure_set();
    for (copied, geometry) in copy.iter().zip(&set) {
        assert!(copied.payload().shares_storage_with(geometry.payload()));
    }
}

/// A frame without a time or step that the first frame has, or with one it
/// lacks, is refused, naming the frame; the frames before it stay.
#[test]
fn a_memory_writer_refuses_frames_whose_time_or_step_columns_differ() {
    let topology = support::topology(&["C", "H"], &[(0, 1)]);
    let mut writer = MemoryTrajectoryWriter::new(Arc::clone(&topology));
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    buffer
        .set_time(Some(Quantity::new(1.0, PICOSECOND)))
        .unwrap();
    writer.write_frame(buffer.frame_view()).unwrap();

    buffer.set_time(None).unwrap();
    let error = writer.write_frame(buffer.frame_view()).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::InconsistentMetadata);
    assert_eq!(error.frame(), Some(1));
    assert_eq!(
        error
            .source()
            .and_then(|source| source.downcast_ref::<StructureSetError>()),
        Some(&StructureSetError::MissingTime { geometry: 1 })
    );
    assert_eq!(
        error.to_string(),
        "Cannot write frame 1 of the trajectory. The frame has a time or step the first frame \
         lacks, or lacks one it has."
    );

    buffer
        .set_time(Some(Quantity::new(2.0, PICOSECOND)))
        .unwrap();
    buffer.set_step(Some(2));
    let error = writer.write_frame(buffer.frame_view()).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::InconsistentMetadata);
    assert_eq!(error.frame(), Some(1));
    assert_eq!(
        error
            .source()
            .and_then(|source| source.downcast_ref::<StructureSetError>()),
        Some(&StructureSetError::UnexpectedStep { geometry: 1 })
    );
    assert_eq!(writer.into_structure_set().len(), 1);
}

/// Occupancies and B-factors are frame state like any other: a decoded
/// frame has none, a reset clears them, and a copy takes the source's.
#[test]
fn buffers_keep_no_stale_occupancies_or_b_factors() {
    let topology = support::topology(&["C", "H"], &[(0, 1)]);
    let mut sited = Geometry::new(Positions::zeros(2));
    sited.set_occupancies(Some(vec![Some(0.5), None])).unwrap();
    sited
        .set_b_factors(Some(Quantity::new(vec![None, Some(0.2)], SQUARE_NANOMETER)))
        .unwrap();
    let sites = |buffer: &FrameBuffer| {
        (
            buffer.occupancies().map(<[_]>::to_vec),
            buffer.b_factors().map(|values| values.value().to_vec()),
        )
    };
    let none = (None, None);
    let expected = (Some(vec![Some(0.5), None]), Some(vec![None, Some(0.2)]));

    // A decoded frame replaces the geometry's occupancies and B-factors.
    let mut reader = XyzReader::new(
        Cursor::new(TWO_FRAMES.as_bytes()),
        Arc::clone(&topology),
        XyzReadOptions::default(),
    )
    .unwrap();
    let mut buffer = FrameBuffer::from_geometry(Arc::clone(&topology), sited.clone()).unwrap();
    assert_eq!(sites(&buffer), expected);
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_eq!(sites(&buffer), none);

    // A reset clears those set through the geometry editor.
    let mut geometry = buffer.geometry_mut();
    geometry.set_occupancies(Some(vec![Some(1.0); 2])).unwrap();
    geometry
        .set_b_factors(Some(Quantity::new(vec![Some(0.1); 2], SQUARE_NANOMETER)))
        .unwrap();
    buffer.reset_dynamic_state();
    assert_eq!(sites(&buffer), none);

    // A copy takes the source's, and a source without them clears them.
    let source = FrameBuffer::from_geometry(Arc::clone(&topology), sited).unwrap();
    buffer.copy_from(source.frame_view()).unwrap();
    assert_eq!(sites(&buffer), expected);
    buffer
        .copy_from(FrameBuffer::new(Arc::clone(&topology)).frame_view())
        .unwrap();
    assert_eq!(sites(&buffer), none);
}

#[test]
fn buffers_from_geometries_validate_them_against_the_topology() {
    let topology = support::topology(&["C", "H"], &[(0, 1)]);
    // A detached geometry acquires its bond rows.
    let buffer =
        FrameBuffer::from_geometry(Arc::clone(&topology), Geometry::new(Positions::zeros(2)))
            .unwrap();
    assert_eq!(buffer.properties().bonds().len(), 1);
    assert_eq!(
        FrameBuffer::from_geometry(topology, Geometry::new(Positions::zeros(3))).unwrap_err(),
        TrajectoryIoError::from(GeometryError::Array(Box::new(
            DenseArrayError::CountMismatch(LengthMismatch::new(Countable::Row, 2, 3))
        )))
    );
}

#[test]
fn errors_retain_nested_causes() {
    let position_error =
        Positions::new(Quantity::new([Point3::new(f64::NAN, 0.0, 0.0)], NANOMETER)).unwrap_err();
    let frame_error = GeometryError::from(position_error);
    let trajectory_error = TrajectoryIoError::from(frame_error.clone());
    assert_eq!(
        trajectory_error.kind(),
        TrajectoryIoErrorKind::InvalidGeometry
    );
    assert_eq!(
        trajectory_error
            .source()
            .unwrap()
            .downcast_ref::<GeometryError>()
            .map(ToString::to_string),
        Some(frame_error.to_string())
    );
}
