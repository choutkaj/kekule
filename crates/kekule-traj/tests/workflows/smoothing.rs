use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{Forces, Positions, Trajectory, TrajectoryFrame, Velocities};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};
use kekule_traj::periodic::{self, FrameSmoother, PeriodicError};
use kekule_traj::FrameBuffer;

use crate::support;
use support::linear_carbon_topology;

fn cell() -> PeriodicCell {
    PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER),
        [true; 3],
    )
    .unwrap()
}

/// A whole C–C molecule along x whose first atom starts at each `x`.
fn trajectory(xs: &[f64], cell: Option<PeriodicCell>) -> Trajectory {
    let mut trajectory = Trajectory::new(linear_carbon_topology(2));
    for (index, &x) in xs.iter().enumerate() {
        let points = [Point3::new(x, 0.5, 0.5), Point3::new(x + 0.15, 0.5, 0.5)];
        let mut frame =
            TrajectoryFrame::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap());
        frame.set_cell(cell);
        frame
            .set_time(Some(Quantity::new(index as f64, PICOSECOND)))
            .unwrap();
        frame.set_step(Some(index as u64 * 10));
        frame
            .set_velocities(Some(
                Velocities::new(Quantity::new(
                    [Vector3::new(1.0, 2.0, 3.0); 2],
                    kekule::units::CANONICAL_VELOCITY_UNIT,
                ))
                .unwrap(),
            ))
            .unwrap();
        frame.set_forces(Some(Forces::zeros(2))).unwrap();
        frame
            .properties_mut()
            .owner_mut()
            .insert(
                PropertyKey::new("label").unwrap(),
                PropertyValue::Int(index as i64),
            )
            .unwrap();
        trajectory.push(frame).unwrap();
    }
    trajectory
}

fn xs(trajectory: &Trajectory, frame: usize) -> Vec<f64> {
    let frame = trajectory.get(frame).unwrap();
    let positions = frame.positions().values();
    positions.value().iter().map(|point| point.x).collect()
}

fn smooth_frame(
    smoother: &FrameSmoother,
    trajectory: &Trajectory,
    index: usize,
) -> Result<FrameBuffer, PeriodicError> {
    let mut buffer = FrameBuffer::new(trajectory.shared_topology());
    buffer.copy_from(trajectory.get(index).unwrap()).unwrap();
    let neighbors = smoother
        .window(index, trajectory.len())
        .filter(|&neighbor| neighbor != index)
        .map(|neighbor| (neighbor, trajectory.get(neighbor).unwrap()));
    smoother.smooth(index, &mut buffer, neighbors)?;
    Ok(buffer)
}

fn close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (a, e) in actual.iter().zip(expected) {
        assert!((a - e).abs() < 1.0e-12, "{actual:?} != {expected:?}");
    }
}

#[test]
fn windows_are_truncated_at_both_trajectory_ends() {
    let smoother = FrameSmoother::new(linear_carbon_topology(1), 2);
    assert_eq!(smoother.half_width(), 2);
    assert_eq!(smoother.window(0, 10), 0..3);
    assert_eq!(smoother.window(1, 10), 0..4);
    assert_eq!(smoother.window(5, 10), 3..8);
    assert_eq!(smoother.window(9, 10), 7..10);
    assert_eq!(smoother.window(0, 1), 0..1);
    assert!(smoother.window(10, 10).is_empty());
    let wide = FrameSmoother::new(linear_carbon_topology(1), usize::MAX);
    assert_eq!(wide.window(usize::MAX - 1, usize::MAX), 0..usize::MAX);
}

#[test]
fn a_molecule_crossing_a_periodic_boundary_is_averaged_through_its_images() {
    // The molecule moves +0.1 nm per frame and wraps across x = 1 after frame 1.
    let source = trajectory(&[0.85, 0.95, 0.05, 0.15], Some(cell()));
    let mut wrapped = source.clone();
    for frame in 0..wrapped.len() {
        let mut frame = wrapped.get_mut(frame).unwrap();
        let points: Vec<_> = frame
            .positions()
            .values()
            .value()
            .iter()
            .map(|point| Point3::new(point.x.rem_euclid(1.0), point.y, point.z))
            .collect();
        frame
            .set_positions(Quantity::new(points, NANOMETER))
            .unwrap();
    }
    let smoother = FrameSmoother::new(source.shared_topology(), 1);

    // Each atom of an interior frame is imaged onto that frame's own atom.
    let middle = smooth_frame(&smoother, &wrapped, 2).unwrap();
    let positions = middle.positions().values();
    let middle: Vec<_> = positions.value().iter().map(|point| point.x).collect();
    close(&middle, &[0.05, 0.2]);
    let middle = smooth_frame(&smoother, &wrapped, 1).unwrap();
    let positions = middle.positions().values();
    let middle: Vec<_> = positions.value().iter().map(|point| point.x).collect();
    close(&middle, &[0.95, 0.1]);

    // The loaded form matches the streaming form frame by frame, ends included.
    let mut smoothed = wrapped.clone();
    periodic::smooth_trajectory(&mut smoothed, 1).unwrap();
    close(&xs(&smoothed, 0), &[0.9, 0.05]);
    close(&xs(&smoothed, 1), &[0.95, 0.1]);
    close(&xs(&smoothed, 2), &[0.05, 0.2]);
    close(&xs(&smoothed, 3), &[0.1, 0.25]);
    for index in 0..wrapped.len() {
        let buffer = smooth_frame(&smoother, &wrapped, index).unwrap();
        assert_eq!(
            buffer.frame_view().payload(),
            smoothed.get(index).unwrap().payload()
        );
    }
}

#[test]
fn cell_free_frames_take_a_plain_average_and_keep_all_other_state() {
    let source = trajectory(&[0.0, 0.9, 0.3, 3.0], None);
    let mut smoothed = source.clone();
    periodic::smooth_trajectory(&mut smoothed, 1).unwrap();
    close(&xs(&smoothed, 0), &[0.45, 0.6]);
    close(&xs(&smoothed, 1), &[0.4, 0.55]);
    close(&xs(&smoothed, 3), &[1.65, 1.8]);
    for index in 0..source.len() {
        let original = source.get(index).unwrap();
        let frame = smoothed.get(index).unwrap();
        assert_eq!(frame.cell(), original.cell());
        assert_eq!(frame.time(), original.time());
        assert_eq!(frame.step(), original.step());
        assert_eq!(frame.velocities(), original.velocities());
        assert_eq!(frame.forces(), original.forces());
        assert_eq!(frame.properties(), original.properties());
    }
    assert_eq!(smoothed.properties(), source.properties());

    // A zero half width, or a lone frame, leaves positions unchanged.
    let mut unchanged = source.clone();
    periodic::smooth_trajectory(&mut unchanged, 0).unwrap();
    assert_eq!(format!("{unchanged:?}"), format!("{source:?}"));
    let mut empty = Trajectory::new(source.shared_topology());
    periodic::smooth_trajectory(&mut empty, 3).unwrap();
    assert!(empty.is_empty());
}

#[test]
fn smoothing_failures_leave_the_buffer_and_trajectory_unchanged() {
    let source = trajectory(&[0.1, 0.2, 0.3, 0.4, 0.5], Some(cell()));
    let smoother = FrameSmoother::new(source.shared_topology(), 1);
    let mut buffer = FrameBuffer::new(source.shared_topology());
    buffer.copy_from(source.get(2).unwrap()).unwrap();
    let before = format!("{buffer:?}");
    let frame = |index| (index, source.get(index).unwrap());

    for (neighbors, neighbor) in [
        (vec![frame(0)], 0),
        (vec![frame(2)], 2),
        (vec![frame(3), frame(1)], 1),
        (vec![frame(1), frame(1)], 1),
        (vec![frame(1), frame(4)], 4),
    ] {
        assert_eq!(
            smoother.smooth(2, &mut buffer, neighbors),
            Err(PeriodicError::InvalidNeighbor { frame: 2, neighbor })
        );
        assert_eq!(format!("{buffer:?}"), before);
    }

    let foreign = trajectory(&[0.1], Some(cell()));
    assert_eq!(
        smoother.smooth(2, &mut buffer, [(3, foreign.get(0).unwrap())]),
        Err(PeriodicError::TopologyMismatch { frame: 3 })
    );
    let mut foreign_buffer = FrameBuffer::new(foreign.shared_topology());
    assert_eq!(
        smoother.smooth(7, &mut foreign_buffer, []),
        Err(PeriodicError::TopologyMismatch { frame: 7 })
    );

    let mut mixed = source.clone();
    mixed.get_mut(3).unwrap().set_cell(None);
    assert_eq!(
        smoother.smooth(2, &mut buffer, [frame(1), (3, mixed.get(3).unwrap())]),
        Err(PeriodicError::MissingCell { frame: 3 })
    );
    assert_eq!(format!("{buffer:?}"), before);
    let mut cell_free = FrameBuffer::new(source.shared_topology());
    cell_free.copy_from(mixed.get(3).unwrap()).unwrap();
    let cell_free_before = format!("{cell_free:?}");
    assert_eq!(
        smoother.smooth(3, &mut cell_free, [frame(2)]),
        Err(PeriodicError::MissingCell { frame: 3 })
    );
    assert_eq!(format!("{cell_free:?}"), cell_free_before);

    let loaded_before = format!("{mixed:?}");
    assert_eq!(
        periodic::smooth_trajectory(&mut mixed, 1),
        Err(PeriodicError::MissingCell { frame: 3 })
    );
    assert_eq!(format!("{mixed:?}"), loaded_before);

    // The failed calls left the buffer usable for a valid one.
    smoother
        .smooth(2, &mut buffer, [frame(1), frame(3)])
        .unwrap();
    let positions = buffer.positions().values();
    let x: Vec<_> = positions.value().iter().map(|point| point.x).collect();
    close(&x, &[0.3, 0.45]);
}
