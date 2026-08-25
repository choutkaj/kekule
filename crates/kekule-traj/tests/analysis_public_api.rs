use std::sync::Arc;

use kekule::alignment::PeriodicAlignmentPolicy;
use kekule::geometry::Point3;
use kekule::structure::Positions;
use kekule::topology::{AtomSelection, Topology};
use kekule::units::{Quantity, ANGSTROM};
use kekule_traj::analysis::{
    AlignedRmsdOptions, PeriodicRmsdPolicy, RmsdOptions, SuperpositionOptions,
};
use kekule_traj::{Trajectory, TrajectoryFrame};

mod support;
use support::linear_carbon_topology;

fn topology() -> Arc<Topology> {
    linear_carbon_topology(3)
}

fn frame(topology: &Arc<Topology>, points: [Point3; 3]) -> TrajectoryFrame {
    TrajectoryFrame::new(
        Positions::new(Quantity::new(points, ANGSTROM)).unwrap(),
        topology.bond_count(),
    )
}

#[test]
fn downstream_code_can_split_or_fuse_superposition_and_rmsd() {
    let topology = topology();
    let reference = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let moving = reference.map(|point| Point3::new(point.x + 4.0, point.y - 2.0, point.z + 1.0));
    let trajectory = Trajectory::from_frames(
        Arc::clone(&topology),
        [frame(&topology, reference), frame(&topology, moving)],
    )
    .unwrap();
    let selection =
        AtomSelection::from_atoms(&topology, topology.atom_ids().iter().copied()).unwrap();

    let direct = trajectory
        .rmsd_to_frame_with_options(
            0,
            &selection,
            RmsdOptions {
                periodic_policy: PeriodicRmsdPolicy::RejectPeriodic,
                ..RmsdOptions::default()
            },
        )
        .unwrap();
    assert!(direct.value()[1] > 0.4);

    let fused = trajectory
        .aligned_rmsd_to_frame_with_options(
            0,
            &selection,
            &selection,
            AlignedRmsdOptions {
                superposition: SuperpositionOptions {
                    periodic_policy: PeriodicAlignmentPolicy::RejectPeriodic,
                    ..SuperpositionOptions::default()
                },
                ..AlignedRmsdOptions::default()
            },
        )
        .unwrap();
    assert!(fused.value()[1] < 1.0e-12);

    let mut split = trajectory;
    let report = split.superpose_to_frame(0, &selection).unwrap();
    assert_eq!(report.alignments().len(), split.len());
    let measured = split.rmsd_to_frame(0, &selection).unwrap();
    assert!(measured.value()[1] < 1.0e-12);
}
