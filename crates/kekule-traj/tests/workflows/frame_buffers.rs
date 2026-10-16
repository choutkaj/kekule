//! Error chains of streamed trajectory frames. Collection behavior itself is
//! covered by the `kekule` crate.

use std::error::Error;

use kekule::geometry::Point3;
use kekule::structure::{Positions, StructureError, TrajectoryError};
use kekule::units::{Quantity, NANOMETER};
use kekule_traj::TrajectoryIoError;

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
