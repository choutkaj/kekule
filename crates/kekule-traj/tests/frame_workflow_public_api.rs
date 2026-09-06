use std::error::Error;
use std::sync::Arc;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{Properties, PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{Model, Positions};
use kekule::topology::AtomSelection;
use kekule::units::{Quantity, NANOMETER, PICOSECOND};
use kekule_traj::{
    Forces, FrameError, Trajectory, TrajectoryError, TrajectoryFrame, TrajectorySliceError,
    Velocities,
};

mod support;
use support::linear_carbon_topology;

fn annotated() -> Trajectory {
    let topology = linear_carbon_topology(3);
    let mut trajectory = Trajectory::new(topology);
    trajectory
        .insert_property(PropertyKey::new("run").unwrap(), PropertyValue::Int(7))
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
        frame.set_cell(Some(
            PeriodicCell::orthorhombic(
                Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
                [true; 3],
            )
            .unwrap(),
        ));
        frame
            .set_time(Some(Quantity::new(index as f64 * 0.5, PICOSECOND)))
            .unwrap();
        frame.set_step(Some(index * 10));
        frame.set_velocities(Some(Velocities::zeros(3))).unwrap();
        frame.set_forces(Some(Forces::zeros(3))).unwrap();
        frame
            .insert_property(
                PropertyKey::new("frame").unwrap(),
                PropertyValue::Int(index as i64),
            )
            .unwrap();
        frame
            .insert_atom_property_column(
                PropertyKey::new("atom").unwrap(),
                PropertyColumn::Int(vec![Some(1), None, Some(3)]),
            )
            .unwrap();
        frame
            .insert_bond_property_column(
                PropertyKey::new("bond").unwrap(),
                PropertyColumn::Int(vec![Some(4), Some(5)]),
            )
            .unwrap();
        trajectory.push(frame).unwrap();
    }
    trajectory
}

#[test]
fn frame_selection_retains_exact_state_order_duplicates_and_empty_topology_binding() {
    let original = annotated();
    let before = format!("{original:?}");
    let selected = original.select_frames([3, 1, 1, 0]).unwrap();
    assert!(Arc::ptr_eq(
        &selected.shared_topology(),
        &original.shared_topology()
    ));
    assert_eq!(selected.properties(), original.properties());
    for (actual, index) in selected.frames().zip([3, 1, 1, 0]) {
        assert_eq!(actual.to_frame(), original.frame(index).unwrap().to_frame());
    }
    assert_eq!(
        selected.validate_monotonic_time(true),
        Err(TrajectoryError::NonMonotonicTime { frame: 1 })
    );
    let strided = original.select_frames((0..4).step_by(2)).unwrap();
    assert_eq!(strided.frame(1).unwrap().step(), Some(20));
    assert_eq!(
        strided.frame(1).unwrap().time(),
        Some(Quantity::new(1.0, PICOSECOND))
    );
    let empty = original.select_frames([]).unwrap();
    assert!(empty.is_empty());
    assert!(Arc::ptr_eq(
        &empty.shared_topology(),
        &original.shared_topology()
    ));
    assert_eq!(empty.properties(), original.properties());
    assert!(matches!(
        original.select_frames([0, 4]),
        Err(TrajectoryError::FrameIndexOutOfRange(4))
    ));
    assert_eq!(format!("{original:?}"), before);
    let atoms = AtomSelection::all(&original.shared_topology());
    assert_eq!(original.slice(&atoms).unwrap().len(), 4);
}

#[test]
fn owned_frames_are_complete_and_independent_and_replacement_is_transactional() {
    let mut original = annotated();
    let old = original.frame(1).unwrap().to_frame();
    let mut owned = old.clone();
    owned
        .set_positions(Quantity::new([Point3::origin(); 3], NANOMETER))
        .unwrap();
    owned.set_step(Some(999));
    assert_eq!(original.frame(1).unwrap().to_frame(), old);
    assert_eq!(owned.properties(), old.properties());
    assert_eq!(owned.velocities(), old.velocities());
    assert_eq!(owned.forces(), old.forces());
    assert_eq!(owned.cell(), old.cell());
    assert_eq!(original.replace_frame(1, owned.clone()).unwrap(), old);
    assert_eq!(original.frame(1).unwrap().to_frame(), owned);
    let before = format!("{original:?}");
    let invalid = TrajectoryFrame::new(Positions::zeros(2));
    assert_eq!(
        original.replace_frame(4, invalid.clone()),
        Err(TrajectoryError::FrameIndexOutOfRange(4))
    );
    assert!(original.replace_frame(1, invalid).is_err());
    assert_eq!(format!("{original:?}"), before);
}

#[test]
#[allow(clippy::forget_non_drop)] // Regression: forgetting an editor must never bypass validation.
fn stored_editor_keeps_dimensions_even_after_columns_are_removed_or_editor_is_forgotten() {
    let mut trajectory = annotated();
    let key = PropertyKey::new("bond").unwrap();
    let before = trajectory.frame(0).unwrap().to_frame();
    {
        let mut frame = trajectory.frame_mut(0).unwrap();
        assert!(frame
            .set_positions(Quantity::new([Point3::origin(); 2], NANOMETER))
            .is_err());
        assert!(frame.set_velocities(Some(Velocities::zeros(4))).is_err());
        assert!(frame.set_forces(Some(Forces::zeros(4))).is_err());
        assert!(frame
            .set_time(Some(Quantity::new(f64::NAN, PICOSECOND)))
            .is_err());
        assert!(frame.set_occupancy_at(0, Some(f64::NAN)).is_err());
        assert!(frame.set_properties(Properties::realization(3, 1)).is_err());
    }
    assert_eq!(trajectory.frame(0).unwrap().to_frame(), before);
    {
        let mut frame = trajectory.frame_mut(0).unwrap();
        frame.remove_bond_property_column(&key);
        assert!(frame
            .insert_bond_property_column(key.clone(), PropertyColumn::Int(vec![None; 1]))
            .is_err());
        assert_eq!(frame.bond_properties().len(), 2);
        frame
            .insert_bond_property_column(key.clone(), PropertyColumn::Int(vec![None; 2]))
            .unwrap();
        frame
            .set_bond_property(1, key, Some(PropertyValue::Int(9)))
            .unwrap();
        frame
            .set_positions(Quantity::new([Point3::origin(); 3], NANOMETER))
            .unwrap();
        frame.set_step(Some(500));
        // Invariants do not depend on Drop.
        std::mem::forget(frame);
    }
    assert_eq!(trajectory.frame(0).unwrap().step(), Some(500));
    trajectory
        .frame(0)
        .unwrap()
        .to_frame()
        .validate(&trajectory.shared_topology())
        .unwrap();
    assert!(trajectory.frame_mut(10).is_none());
}

#[test]
fn errors_retain_nested_causes_and_model_validation_is_not_misreported_as_topology_mismatch() {
    let topology = linear_carbon_topology(3);
    let model_error = Model::new(topology, Positions::zeros(2)).unwrap_err();
    let error = TrajectoryError::from(model_error.clone());
    assert!(matches!(error, TrajectoryError::Model(_)));
    assert_eq!(error.source().unwrap().to_string(), model_error.to_string());
    let position_error =
        Positions::new(Quantity::new([Point3::new(f64::NAN, 0.0, 0.0)], NANOMETER)).unwrap_err();
    let frame_error = FrameError::from(position_error);
    let source_text = frame_error.source().unwrap().to_string();
    let trajectory_error = TrajectoryError::from(frame_error);
    assert_eq!(
        trajectory_error
            .source()
            .unwrap()
            .source()
            .unwrap()
            .to_string(),
        source_text
    );
    let slice_error = TrajectorySliceError::from(trajectory_error);
    assert!(slice_error
        .source()
        .unwrap()
        .source()
        .unwrap()
        .source()
        .is_some());
}
