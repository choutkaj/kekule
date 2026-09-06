use std::sync::Arc;

use kekule::alignment::PeriodicAlignmentPolicy;
use kekule::geometry::Point3;
use kekule::properties::{PropertyKey, PropertyValue};
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

fn frame(points: [Point3; 3]) -> TrajectoryFrame {
    TrajectoryFrame::new(Positions::new(Quantity::new(points, ANGSTROM)).unwrap())
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
    let key = PropertyKey::new("simulation").unwrap();
    let mut annotated = frame(moving);
    annotated
        .insert_property(key.clone(), PropertyValue::Int(7))
        .unwrap();
    let mut trajectory =
        Trajectory::from_frames(Arc::clone(&topology), [frame(reference), annotated]).unwrap();
    trajectory
        .insert_property(key.clone(), PropertyValue::String("run_1".into()))
        .unwrap();
    let owner_properties = trajectory.properties().clone();
    let frame_properties = trajectory.frame(1).unwrap().properties().clone();
    let selection = AtomSelection::all(&topology);

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
    assert!(split.superpose_to_frame_in_place(99, &selection).is_err());
    assert_eq!(split.properties(), &owner_properties);
    assert_eq!(
        split.frame(1).unwrap().positions().values(),
        frame(moving).positions().values()
    );
    let (aligned, report) = split.superpose_to_frame_with_report(0, &selection).unwrap();
    assert_eq!(
        split.frame(1).unwrap().positions().values(),
        frame(moving).positions().values()
    );
    split.superpose_to_frame_in_place(0, &selection).unwrap();
    assert_eq!(
        split.frame(1).unwrap().positions(),
        aligned.frame(1).unwrap().positions()
    );
    assert_eq!(report.alignments().len(), split.len());
    assert!(Arc::ptr_eq(&topology, &split.shared_topology()));
    assert_eq!(split.properties(), &owner_properties);
    assert_eq!(split.frame(1).unwrap().properties(), &frame_properties);
    let measured = split.rmsd_to_frame(0, &selection).unwrap();
    assert!(measured.value()[1] < 1.0e-12);
}

#[test]
fn ordinary_superposition_returns_a_copy_and_accepts_periodic_coordinates() {
    use kekule::alignment::AlignmentWeighting;
    use kekule::geometry::{PeriodicCell, Vector3};
    let topology = topology();
    let points = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let mut reference = frame(points);
    let mut moving = frame(points.map(|p| Point3::new(p.x + 3.0, p.y - 2.0, p.z + 1.0)));
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(10.0, 10.0, 10.0), ANGSTROM),
        [true; 3],
    )
    .unwrap();
    reference.set_cell(Some(cell));
    moving.set_cell(Some(cell));
    let mut original = Trajectory::from_frames(topology.clone(), [reference, moving]).unwrap();
    original
        .insert_property(PropertyKey::new("run").unwrap(), PropertyValue::Int(42))
        .unwrap();
    let before = format!("{original:?}");
    let fit = AtomSelection::all(&topology);
    assert!(original.rmsd_to_frame(0, &fit).unwrap().value()[1] > 0.3);
    let aligned = original.superpose_to_frame(0, &fit).unwrap();
    assert_eq!(format!("{original:?}"), before);
    assert_eq!(aligned.properties(), original.properties());
    assert!(Arc::ptr_eq(&topology, &aligned.shared_topology()));
    assert!(aligned.rmsd_to_frame(0, &fit).unwrap().value()[1] < 1.0e-12);
    assert!(
        original
            .aligned_rmsd_to_frame(0, &fit, &fit)
            .unwrap()
            .value()[1]
            < 1.0e-12
    );
    let weighted = original
        .superpose_to_frame_with_options(
            0,
            &fit,
            SuperpositionOptions {
                weighting: AlignmentWeighting::Explicit(&[1.0, 2.0, 3.0]),
                ..Default::default()
            },
        )
        .unwrap();
    assert!(weighted.rmsd_to_frame(0, &fit).unwrap().value()[1] < 1.0e-12);
    assert!(original
        .superpose_to_frame_in_place_with_options(
            0,
            &fit,
            SuperpositionOptions {
                periodic_policy: PeriodicAlignmentPolicy::RejectPeriodic,
                ..Default::default()
            }
        )
        .is_err());
    assert_eq!(format!("{original:?}"), before);
    original.superpose_to_frame_in_place(0, &fit).unwrap();
    for (actual, expected) in original.frames().zip(aligned.frames()) {
        assert_eq!(actual.positions(), expected.positions());
        assert_eq!(actual.cell(), expected.cell());
        assert_eq!(actual.properties(), expected.properties());
    }
}
