use kekule::{
    geometry::{PeriodicCell, Point3, Vector3},
    properties::{PropertyKey, PropertyValue},
    structure::Positions,
    topology::AtomSelection,
    units::{Quantity, NANOMETER, PICOSECOND},
};
use kekule_traj::{
    analysis::FrameSuperposer,
    periodic::{MoleculeImager, PeriodicError, TrajectoryUnwrapper},
    Forces, FrameBuffer, MemoryTrajectoryReader, Trajectory, TrajectoryFrame, TrajectoryReader,
    Velocities,
};
use std::sync::Arc;

mod support;
use support::linear_carbon_topology;

fn source() -> Trajectory {
    let mut trajectory = Trajectory::new(linear_carbon_topology(3));
    for (index, x) in [0.9, 0.1, 0.4, 0.8, 0.2].into_iter().enumerate() {
        let mut frame = TrajectoryFrame::new(
            Positions::new(Quantity::new(
                [
                    Point3::new(x, 0.2, 0.2),
                    Point3::new((x + 0.2) % 1.0, 0.2, 0.2),
                    Point3::new(x, 0.4, 0.2),
                ],
                NANOMETER,
            ))
            .unwrap(),
        );
        frame.set_cell(Some(
            PeriodicCell::orthorhombic(
                Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER),
                [true; 3],
            )
            .unwrap(),
        ));
        frame
            .set_time(Some(Quantity::new(index as f64, PICOSECOND)))
            .unwrap();
        frame.set_step(Some(index as u64 * 10));
        frame.set_velocities(Some(Velocities::zeros(3))).unwrap();
        frame.set_forces(Some(Forces::zeros(3))).unwrap();
        frame
            .insert_property(
                PropertyKey::new("label").unwrap(),
                PropertyValue::Int(index as i64),
            )
            .unwrap();
        trajectory.push(frame).unwrap();
    }
    trajectory
}

#[test]
fn streaming_and_loaded_whole_image_and_unwrap_preserve_identical_complete_frames() {
    let trajectory = source();
    let topology = trajectory.shared_topology();
    let anchors = AtomSelection::all(&topology);
    let imager = MoleculeImager::new(topology.clone());
    let expected_whole = trajectory.make_molecules_whole().unwrap();
    let expected_image = trajectory.image_molecules(&anchors).unwrap();
    let expected_unwrap = expected_whole.unwrap().unwrap();
    let mut unwrapper = TrajectoryUnwrapper::new(topology.clone());
    let mut reader = MemoryTrajectoryReader::new(&trajectory);
    let mut buffer = reader.frame_buffer();
    for index in 0..trajectory.len() {
        assert!(reader.read_next(&mut buffer).unwrap());
        let image = imager.image(index, buffer.frame_view(), &anchors).unwrap();
        assert_eq!(image, expected_image.frame(index).unwrap().to_frame());
        let whole = imager.make_whole(index, buffer.frame_view()).unwrap();
        imager.make_whole_in_place(index, &mut buffer).unwrap();
        assert_eq!(whole, buffer.frame_view().to_frame());
        assert_eq!(whole, expected_whole.frame(index).unwrap().to_frame());
        unwrapper.unwrap_in_place(index, &mut buffer).unwrap();
        assert_eq!(
            buffer.frame_view().to_frame(),
            expected_unwrap.frame(index).unwrap().to_frame()
        );
        // The same unwrapper survives processing-chunk boundaries.
        assert_eq!(unwrapper.last_frame_index(), Some(index));
        buffer.copy_from(trajectory.frame(index).unwrap()).unwrap();
        imager.image_in_place(index, &mut buffer, &anchors).unwrap();
        assert_eq!(buffer.frame_view().to_frame(), image);
        assert!(Arc::ptr_eq(&buffer.shared_topology(), &topology));
    }
    assert!(!reader.read_next(&mut buffer).unwrap());
    let mut copying = TrajectoryUnwrapper::new(topology);
    for (index, frame) in expected_whole.frames().enumerate() {
        assert_eq!(
            copying.unwrap_frame(index, frame).unwrap(),
            expected_unwrap.frame(index).unwrap().to_frame()
        );
    }
}

#[test]
fn unwrapper_rejects_skips_reordering_and_failures_without_advancing_and_supports_reset() {
    let trajectory = source();
    let topology = trajectory.shared_topology();
    let mut unwrapper = TrajectoryUnwrapper::new(topology.clone());
    let mut buffer = FrameBuffer::new(topology.clone());
    buffer.copy_from(trajectory.frame(0).unwrap()).unwrap();
    unwrapper.unwrap_in_place(10, &mut buffer).unwrap();
    buffer.copy_from(trajectory.frame(1).unwrap()).unwrap();
    let before = format!("{buffer:?}");
    for index in [10, 9, 12] {
        assert_eq!(
            unwrapper.unwrap_in_place(index, &mut buffer),
            Err(PeriodicError::NonSequentialFrame {
                previous: 10,
                frame: index
            })
        );
        assert_eq!(format!("{buffer:?}"), before);
        assert_eq!(unwrapper.last_frame_index(), Some(10));
    }
    let cell = buffer.cell().copied();
    buffer.set_cell(None);
    let missing = format!("{buffer:?}");
    assert_eq!(
        unwrapper.unwrap_in_place(11, &mut buffer),
        Err(PeriodicError::MissingCell { frame: 11 })
    );
    assert_eq!(format!("{buffer:?}"), missing);
    assert_eq!(unwrapper.last_frame_index(), Some(10));
    buffer.set_cell(cell);
    unwrapper.unwrap_in_place(11, &mut buffer).unwrap();
    assert_eq!(unwrapper.last_frame_index(), Some(11));
    unwrapper.reset();
    buffer.copy_from(trajectory.frame(0).unwrap()).unwrap();
    let before = buffer.frame_view().to_frame();
    unwrapper.unwrap_in_place(100, &mut buffer).unwrap();
    assert_eq!(buffer.frame_view().to_frame(), before);
    let mut foreign = FrameBuffer::new(linear_carbon_topology(3));
    assert_eq!(
        unwrapper.unwrap_in_place(101, &mut foreign),
        Err(PeriodicError::TopologyMismatch { frame: 101 })
    );
    assert_eq!(unwrapper.last_frame_index(), Some(100));
    assert_eq!(
        MoleculeImager::new(topology).make_whole_in_place(8, &mut foreign),
        Err(PeriodicError::TopologyMismatch { frame: 8 })
    );
}

#[test]
fn unwrapping_checks_available_time_across_missing_times_transactionally() {
    let mut trajectory = source();
    trajectory.frame_mut(1).unwrap().set_time(None).unwrap();
    trajectory
        .frame_mut(2)
        .unwrap()
        .set_time(Some(Quantity::new(-1.0, PICOSECOND)))
        .unwrap();
    let before = format!("{trajectory:?}");
    assert_eq!(
        trajectory.unwrap_in_place(),
        Err(PeriodicError::NonMonotonicTime { frame: 2 })
    );
    assert_eq!(format!("{trajectory:?}"), before);
}

#[test]
fn ambiguous_crossings_and_changed_axes_do_not_advance_streaming_state() {
    let trajectory = source();
    let topology = trajectory.shared_topology();
    let mut unwrapper = TrajectoryUnwrapper::new(topology.clone());
    let mut buffer = FrameBuffer::new(topology);
    buffer.copy_from(trajectory.frame(0).unwrap()).unwrap();
    unwrapper.unwrap_in_place(0, &mut buffer).unwrap();
    buffer.copy_from(trajectory.frame(1).unwrap()).unwrap();
    let mut points = buffer.positions().values().value().to_vec();
    points[0].x = 1.4; // Half a cell from the previous x = 0.9.
    buffer
        .set_positions(Quantity::new(points, NANOMETER))
        .unwrap();
    let before = buffer.frame_view().to_frame();
    assert!(matches!(
        unwrapper.unwrap_in_place(1, &mut buffer),
        Err(PeriodicError::AmbiguousDisplacement {
            frame: 1,
            axis: 0,
            ..
        })
    ));
    assert_eq!(buffer.frame_view().to_frame(), before);
    assert_eq!(unwrapper.last_frame_index(), Some(0));
    buffer.copy_from(trajectory.frame(1).unwrap()).unwrap();
    buffer.set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER),
            [true, true, false],
        )
        .unwrap(),
    ));
    let before = buffer.frame_view().to_frame();
    assert_eq!(
        unwrapper.unwrap_in_place(1, &mut buffer),
        Err(PeriodicError::PeriodicAxesChanged { frame: 1 })
    );
    assert_eq!(buffer.frame_view().to_frame(), before);
    assert_eq!(unwrapper.last_frame_index(), Some(0));
    buffer.copy_from(trajectory.frame(1).unwrap()).unwrap();
    unwrapper.unwrap_in_place(1, &mut buffer).unwrap();
    assert_eq!(
        buffer.frame_view().to_frame(),
        trajectory.unwrap().unwrap().frame(1).unwrap().to_frame()
    );
}

#[test]
fn superposition_streams_the_same_transforms_and_metadata_and_reports_source_indices() {
    let trajectory = source().make_molecules_whole().unwrap();
    let selection = AtomSelection::all(&trajectory.shared_topology());
    let (expected, reports) = trajectory
        .superpose_to_frame_with_report(0, &selection)
        .unwrap();
    let reference = trajectory.frame(0).unwrap().to_frame();
    let topology = trajectory.shared_topology();
    let superposer = FrameSuperposer::new(reference.view(&topology).unwrap(), &selection);
    let mut buffer = FrameBuffer::new(topology.clone());
    for (index, frame) in trajectory.frames().enumerate() {
        assert_eq!(
            superposer.superpose(index, frame).unwrap(),
            expected.frame(index).unwrap().to_frame()
        );
        buffer.copy_from(frame).unwrap();
        let report = superposer
            .superpose_in_place_with_report(index, &mut buffer)
            .unwrap();
        assert_eq!(&report, reports.alignment(index).unwrap());
        assert_eq!(
            buffer.frame_view().to_frame(),
            expected.frame(index).unwrap().to_frame()
        );
    }
    buffer
        .set_positions(Quantity::new([Point3::origin(); 3], NANOMETER))
        .unwrap();
    let before = format!("{buffer:?}");
    assert!(matches!(
        superposer.superpose_in_place(42, &mut buffer),
        Err(kekule_traj::analysis::SuperpositionError::Alignment { frame: 42, .. })
    ));
    assert_eq!(format!("{buffer:?}"), before);
}
