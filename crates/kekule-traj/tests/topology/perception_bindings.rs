//! Bindings made before perception keep working with the perceived trajectory.

use std::fs;
use std::path::PathBuf;
use std::sync::Arc;

use kekule::geometry::Point3;
use kekule::smiles;
use kekule::statistics::RmsfAccumulator;
use kekule::structure::Positions;
use kekule::structure::{Trajectory, TrajectoryFrame};
use kekule::topology::AtomSelection;
use kekule::units::{Quantity, ANGSTROM};
use kekule_traj::io::{open_trajectory, write_trajectory};
use kekule_traj::{FrameBuffer, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryReader};

struct TemporaryPath(PathBuf);

impl Drop for TemporaryPath {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.0);
    }
}

fn ring_frame(offset: f64) -> TrajectoryFrame {
    let points = (0..6)
        .map(|i| {
            let angle = f64::from(i) * std::f64::consts::FRAC_PI_3;
            Point3::new(1.4 * angle.cos() + offset, 1.4 * angle.sin(), 0.1 * offset)
        })
        .collect::<Vec<_>>();
    TrajectoryFrame::new(Positions::new(Quantity::new(points, ANGSTROM)).unwrap())
}

#[test]
fn selections_buffers_and_accumulators_survive_trajectory_perception() {
    let topology = smiles::to_topology("c1ccccc1").unwrap();
    let mut trajectory =
        Trajectory::from_items(Arc::clone(&topology), [ring_frame(0.0), ring_frame(0.5)]).unwrap();
    let selection = AtomSelection::all(Arc::clone(&topology));
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    let mut rmsf = RmsfAccumulator::new(&selection).unwrap();

    trajectory.perceive().unwrap();
    assert!(!Arc::ptr_eq(&trajectory.shared_topology(), &topology));
    assert!(trajectory.topology().shares_layout(&topology));

    // In-memory analysis accepts the selection made before perception.
    let perceived = trajectory.shared_topology();
    trajectory
        .superpose(0, &selection, Default::default())
        .unwrap();
    assert!(Arc::ptr_eq(&perceived, &trajectory.shared_topology()));
    for (index, frame) in trajectory.iter().enumerate() {
        rmsf.observe(index, &frame, 1.0).unwrap();
    }
    assert_eq!(rmsf.finish().unwrap().observation_count(), 2);

    // Streaming through a reader bound to the perceived snapshot fills the
    // buffer bound to the original one, which keeps its own snapshot.
    let path = TemporaryPath(std::env::temp_dir().join(format!(
        "kekule-perception-bindings-{}.xyz",
        std::process::id()
    )));
    write_trajectory(&path.0, &trajectory, Default::default()).unwrap();
    let mut reader =
        open_trajectory(&path.0, trajectory.shared_topology(), Default::default()).unwrap();
    let mut frames = 0;
    while reader.read_next(&mut buffer).unwrap() {
        frames += 1;
    }
    assert_eq!(frames, 2);
    assert!(Arc::ptr_eq(&buffer.shared_topology(), &topology));

    // An independently published equal topology is still a different layout.
    let independent = smiles::to_topology("c1ccccc1").unwrap();
    assert!(independent.same_layout(&topology));
    assert_eq!(
        FrameBuffer::new(independent).copy_from(trajectory.get(0).unwrap()),
        Err(TrajectoryIoError::new(
            TrajectoryIoErrorKind::TopologyMismatch
        ))
    );
}
