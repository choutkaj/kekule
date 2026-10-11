//! Streaming reader contracts and error chains over in-memory trajectories.
//! Collection behavior itself is covered by the `kekule` crate.

use std::error::Error;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{
    Forces, Positions, StructureError, Trajectory, TrajectoryError, TrajectoryFrame, Velocities,
};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};
use kekule_traj::{
    MemoryTrajectoryReader, SeekableTrajectoryReader, TrajectoryIoError, TrajectoryReader,
};

use crate::support;
use support::linear_carbon_topology;

fn annotated() -> Trajectory {
    let mut trajectory = Trajectory::new(linear_carbon_topology(3));
    trajectory
        .properties_mut()
        .insert(PropertyKey::new("run").unwrap(), PropertyValue::Int(7))
        .unwrap();
    for index in 0..4 {
        let mut frame = TrajectoryFrame::new(
            Positions::new(Quantity::new(
                [
                    Point3::new(index as f64, 0.0, 0.0),
                    Point3::new(0.1, 0.0, 0.0),
                    Point3::new(0.0, 0.1, 0.0),
                ],
                NANOMETER,
            ))
            .unwrap(),
        );
        frame
            .set_time(Some(Quantity::new(index as f64 * 0.5, PICOSECOND)))
            .unwrap();
        frame.set_step(Some(index * 10));
        frame.set_velocities(Some(Velocities::zeros(3))).unwrap();
        frame.set_forces(Some(Forces::zeros(3))).unwrap();
        frame.set_cell(Some(
            PeriodicCell::orthorhombic(
                Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
                [true; 3],
            )
            .unwrap(),
        ));
        let mut properties = frame.properties_mut();
        properties
            .owner_mut()
            .insert(
                PropertyKey::new("frame").unwrap(),
                PropertyValue::Int(index as i64),
            )
            .unwrap();
        properties
            .atoms_mut()
            .insert(
                PropertyKey::new("atom").unwrap(),
                PropertyColumn::Int(vec![Some(1), None, Some(3)]),
            )
            .unwrap();
        trajectory.push(frame).unwrap();
        trajectory
            .get_mut(index as usize)
            .unwrap()
            .properties_mut()
            .bonds_mut()
            .insert(
                PropertyKey::new("bond").unwrap(),
                PropertyColumn::Int(vec![Some(4), Some(5)]),
            )
            .unwrap();
    }
    trajectory
}

#[test]
fn memory_random_reads_preserve_sequential_cursor_and_error_destinations() {
    let trajectory = annotated();
    let mut reader = MemoryTrajectoryReader::new(&trajectory);
    let mut destination = reader.frame_buffer();
    reader.read_frame(3, &mut destination).unwrap();
    assert_eq!(
        (*destination).clone(),
        trajectory.get(3).unwrap().payload().clone()
    );
    assert!(reader.read_next(&mut destination).unwrap());
    assert_eq!(
        (*destination).clone(),
        trajectory.get(0).unwrap().payload().clone()
    );

    reader.read_frame(2, &mut destination).unwrap();
    let before = (*destination).clone();
    assert_eq!(
        reader.read_frame(99, &mut destination),
        Err(TrajectoryIoError::FrameIndexOutOfRange(99))
    );
    assert_eq!((*destination).clone(), before);
    let unrelated = annotated();
    let mut unrelated_buffer = MemoryTrajectoryReader::new(&unrelated).frame_buffer();
    let unrelated_before = (*unrelated_buffer).clone();
    assert_eq!(
        reader.read_frame(0, &mut unrelated_buffer),
        Err(TrajectoryIoError::TopologyMismatch)
    );
    assert_eq!((*unrelated_buffer).clone(), unrelated_before);
    for index in 1..4 {
        assert!(reader.read_next(&mut destination).unwrap());
        assert_eq!(
            (*destination).clone(),
            trajectory.get(index).unwrap().payload().clone()
        );
    }
    reader.read_frame(0, &mut destination).unwrap();
    assert!(!reader.read_next(&mut destination).unwrap());
}

#[test]
fn errors_retain_nested_causes() {
    let position_error =
        Positions::new(Quantity::new([Point3::new(f64::NAN, 0.0, 0.0)], NANOMETER)).unwrap_err();
    let frame_error = StructureError::from(position_error);
    let source_text = frame_error.source().unwrap().to_string();
    let trajectory_error = TrajectoryIoError::from(frame_error.clone());
    assert!(matches!(trajectory_error, TrajectoryIoError::Structure(_)));
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
    assert!(matches!(collection_error, TrajectoryIoError::Trajectory(_)));
    assert!(collection_error.source().unwrap().is::<TrajectoryError>());
}
