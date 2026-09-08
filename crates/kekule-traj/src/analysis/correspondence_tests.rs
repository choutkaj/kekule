use super::*;
use crate::{MemoryTrajectoryReader, TrajectoryReader};
use kekule::alignment::CorrespondenceSide;
use kekule::core::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::geometry::{Point3, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::Model;
use kekule::topology::Topology;
use kekule::units::{CANONICAL_FORCE_UNIT, CANONICAL_VELOCITY_UNIT, NANOMETER, PICOSECOND};
use std::sync::Arc;

fn model(points: &[Point3]) -> Model {
    let mut editor = MoleculeEditor::new();
    let mut previous = None;
    for _ in points {
        let atom = editor
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .unwrap();
        if let Some(previous) = previous {
            editor.add_bond(previous, atom, BondOrder::Single).unwrap();
        }
        previous = Some(atom);
    }
    Model::new(
        Topology::from_molecule(&editor.finish().unwrap()).unwrap(),
        Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
    )
    .unwrap()
}

fn transform(p: Point3) -> Point3 {
    Point3::new(-p.y + 2., p.x - 3., p.z + 1.)
}

fn fixture() -> (Trajectory, Model, AtomCorrespondence, AtomCorrespondence) {
    let points = [
        Point3::new(0., 0., 0.),
        Point3::new(2., 0., 0.),
        Point3::new(0., 3., 0.),
        Point3::new(0., 0., 4.),
        Point3::new(2., 2., 2.),
    ];
    let moving = model(&points);
    let reference = model(&[
        transform(points[2]),
        transform(points[0]),
        transform(points[3]),
        transform(points[1]),
        transform(points[4]),
        Point3::new(99., 99., 99.),
    ]);
    let top = moving.shared_topology();
    let fit = AtomCorrespondence::from_pairs(
        &top,
        &reference.shared_topology(),
        [(3, 2), (0, 1), (2, 0), (1, 3)].map(|(a, b)| (top.atom_ids()[a], reference.atom_ids()[b])),
    )
    .unwrap();
    let measurement = AtomCorrespondence::from_pairs(
        &top,
        &reference.shared_topology(),
        [(top.atom_ids()[4], reference.atom_ids()[4])],
    )
    .unwrap();
    let frames = (0..2)
        .map(|index| {
            let mut points = points;
            points[4].x += index as f64;
            let mut frame =
                TrajectoryFrame::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap());
            frame
                .set_time(Some(Quantity::new(index as f64 + 2., PICOSECOND)))
                .unwrap();
            frame.set_step(Some(index as u64 + 20));
            frame.set_cell(Some(
                PeriodicCell::new(
                    Quantity::new(
                        [
                            Vector3::new(20., 0., 0.),
                            Vector3::new(0., 20., 0.),
                            Vector3::new(0., 0., 20.),
                        ],
                        NANOMETER,
                    ),
                    [true; 3],
                )
                .unwrap(),
            ));
            frame
                .set_velocities(Some(
                    Velocities::new(Quantity::new(
                        vec![Vector3::new(1., 0., 0.); 5],
                        CANONICAL_VELOCITY_UNIT,
                    ))
                    .unwrap(),
                ))
                .unwrap();
            frame
                .set_forces(Some(
                    Forces::new(Quantity::new(
                        vec![Vector3::new(0., 1., 0.); 5],
                        CANONICAL_FORCE_UNIT,
                    ))
                    .unwrap(),
                ))
                .unwrap();
            frame
                .insert_property(
                    PropertyKey::new("label").unwrap(),
                    PropertyValue::String(format!("frame {index}")),
                )
                .unwrap();
            frame
                .set_atom_property(
                    4,
                    PropertyKey::new("tag").unwrap(),
                    Some(PropertyValue::String("moving site".into())),
                )
                .unwrap();
            frame
        })
        .collect::<Vec<_>>();
    let mut trajectory = Trajectory::from_frames(top, frames).unwrap();
    trajectory
        .insert_property(
            PropertyKey::new("source").unwrap(),
            PropertyValue::String("moving system".into()),
        )
        .unwrap();
    (trajectory, reference, fit, measurement)
}

#[test]
fn correspondence_superposition_agrees_loaded_streaming_and_in_place_with_all_fields() {
    let (trajectory, reference, fit, _) = fixture();
    let before = trajectory
        .frames()
        .map(|f| f.to_frame())
        .collect::<Vec<_>>();
    let reference_before = reference.clone();
    let weights = [8., 1., 4., 2.];
    let options = SuperpositionOptions {
        weighting: kekule::alignment::AlignmentWeighting::Explicit(&weights),
        ..Default::default()
    };
    let loaded = trajectory
        .superpose_to_model_with_options(reference.view(), &fit, options)
        .unwrap();
    assert!(Arc::ptr_eq(
        &trajectory.shared_topology(),
        &loaded.shared_topology()
    ));
    assert_eq!(loaded.properties(), trajectory.properties());
    let superposer =
        FrameSuperposer::with_correspondence_and_options(reference.view(), &fit, options);
    let mut reader = MemoryTrajectoryReader::new(&trajectory);
    let mut buffer = reader.frame_buffer();
    let mut index = 0;
    while reader.read_next(&mut buffer).unwrap() {
        let source = buffer.frame_view().to_frame();
        superposer.superpose_in_place(index, &mut buffer).unwrap();
        let actual = buffer.frame_view();
        assert_eq!(actual.to_frame(), loaded.frame(index).unwrap().to_frame());
        assert_eq!(actual.properties(), source.properties());
        assert_eq!(actual.time(), source.time());
        assert_eq!(actual.step(), source.step());
        for (&p, &q) in actual
            .positions()
            .values()
            .value()
            .iter()
            .zip(source.positions().values().value().iter())
        {
            assert!((p - transform(q)).norm() < 1e-12);
        }
        for v in actual.velocities().unwrap().value().iter() {
            assert!((*v - Vector3::new(0., 1., 0.)).norm() < 1e-12);
        }
        for f in actual.forces().unwrap().value().iter() {
            assert!((*f - Vector3::new(-1., 0., 0.)).norm() < 1e-12);
        }
        let vectors = actual.cell().unwrap().vectors();
        assert!((vectors.value()[0] - Vector3::new(0., 20., 0.)).norm() < 1e-12);
        assert_eq!(actual.cell().unwrap().periodic_axes(), [true; 3]);
        index += 1;
    }
    assert_eq!(index, 2);
    let mut in_place = trajectory.clone();
    in_place
        .superpose_to_model_in_place_with_options(reference.view(), &fit, options)
        .unwrap();
    assert_eq!(
        in_place.frames().map(|f| f.to_frame()).collect::<Vec<_>>(),
        loaded.frames().map(|f| f.to_frame()).collect::<Vec<_>>()
    );
    assert_eq!(
        trajectory
            .frames()
            .map(|f| f.to_frame())
            .collect::<Vec<_>>(),
        before
    );
    assert_eq!(reference, reference_before);
}

#[test]
fn correspondence_fit_and_measurement_are_independent_and_weighted_in_pair_order() {
    let (trajectory, reference, fit, measurement) = fixture();
    let fused = trajectory
        .aligned_rmsd_to_model(reference.view(), &fit, &measurement)
        .unwrap();
    assert!(fused.value()[0].abs() < 1e-12);
    assert!((fused.value()[1] - 1.).abs() < 1e-12);
    let loaded = trajectory
        .superpose_to_model(reference.view(), &fit)
        .unwrap();
    let split = loaded
        .rmsd_to_model(reference.view(), &measurement)
        .unwrap();
    for (a, b) in fused.value().iter().zip(split.value()) {
        assert!((a - b).abs() < 1e-12);
    }
    let top = trajectory.shared_topology();
    let measured = AtomCorrespondence::from_pairs(
        &top,
        &reference.shared_topology(),
        [
            (top.atom_ids()[4], reference.atom_ids()[4]),
            (top.atom_ids()[0], reference.atom_ids()[1]),
        ],
    )
    .unwrap();
    let options = AlignedRmsdOptions {
        measurement_weighting: RmsdWeighting::Explicit(&[3., 1.]),
        ..Default::default()
    };
    let weighted = trajectory
        .aligned_rmsd_to_model_with_options(reference.view(), &fit, &measured, options)
        .unwrap();
    assert!((weighted.value()[1] - 0.75_f64.sqrt()).abs() < 1e-12);
    let direct = trajectory
        .rmsd_to_model(reference.view(), &measured)
        .unwrap();
    for (frame, actual) in trajectory.frames().zip(direct.value()) {
        let expected = measured
            .index_pairs()
            .iter()
            .map(|&(a, b)| {
                (frame.positions().values().value()[a.index()]
                    - reference.positions().values().value()[b.index()])
                .norm_squared()
            })
            .sum::<f64>()
            / 2.;
        assert!((actual - expected.sqrt()).abs() < 1e-12);
    }
    assert!(direct.value()[0] > 1.);
}

#[test]
fn correspondence_failures_preserve_buffers_and_whole_trajectories() {
    let (mut trajectory, reference, fit, measurement) = fixture();
    let top = trajectory.shared_topology();
    let invalid = TrajectoryFrame::new(Positions::zeros(top.atom_count()));
    trajectory.push(invalid).unwrap();
    let before = trajectory
        .frames()
        .map(|f| f.to_frame())
        .collect::<Vec<_>>();
    assert!(matches!(
        trajectory.superpose_to_model_in_place(reference.view(), &fit),
        Err(SuperpositionError::Alignment {
            frame: 2,
            source: AlignmentError::DegenerateGeometry { .. }
        })
    ));
    assert_eq!(
        trajectory
            .frames()
            .map(|f| f.to_frame())
            .collect::<Vec<_>>(),
        before
    );
    assert!(matches!(
        trajectory.aligned_rmsd_to_model(reference.view(), &fit, &measurement),
        Err(RmsdError::Alignment { frame: 2, .. })
    ));
    let mut buffer = FrameBuffer::new(top);
    buffer.copy_from(trajectory.frame(2).unwrap()).unwrap();
    let original = buffer.frame_view().to_frame();
    assert!(FrameSuperposer::with_correspondence(reference.view(), &fit)
        .superpose_in_place(17, &mut buffer)
        .is_err());
    assert_eq!(buffer.frame_view().to_frame(), original);
    let wrong = model(reference.positions().values().value());
    assert!(matches!(
        trajectory.rmsd_to_model(wrong.view(), &measurement),
        Err(RmsdError::Correspondence(
            AtomCorrespondenceError::TopologyMismatch {
                side: CorrespondenceSide::Reference
            }
        ))
    ));
    let empty = Trajectory::new(wrong.shared_topology());
    assert!(matches!(
        empty.superpose_to_model(reference.view(), &fit),
        Err(SuperpositionError::Correspondence(
            AtomCorrespondenceError::TopologyMismatch {
                side: CorrespondenceSide::Moving
            }
        ))
    ));
}

#[test]
fn correspondence_measurement_validates_weights_empty_pairs_and_periodic_policy() {
    let (trajectory, reference, fit, measurement) = fixture();
    let empty = AtomCorrespondence::from_pairs(
        &trajectory.shared_topology(),
        &reference.shared_topology(),
        [],
    )
    .unwrap();
    assert_eq!(
        trajectory.rmsd_to_model(reference.view(), &empty),
        Err(RmsdError::EmptySelection)
    );
    assert!(matches!(
        trajectory.rmsd_to_model_with_options(
            reference.view(),
            &measurement,
            RmsdOptions {
                weighting: RmsdWeighting::Explicit(&[]),
                ..Default::default()
            }
        ),
        Err(RmsdError::WeightCountMismatch { .. })
    ));
    assert!(matches!(
        trajectory.rmsd_to_model_with_options(
            reference.view(),
            &measurement,
            RmsdOptions {
                periodic_policy: PeriodicRmsdPolicy::RejectPeriodic,
                ..Default::default()
            }
        ),
        Err(RmsdError::PeriodicCoordinates { frame: 0, .. })
    ));
    assert!(matches!(
        trajectory.aligned_rmsd_to_model_with_options(
            reference.view(),
            &fit,
            &measurement,
            AlignedRmsdOptions {
                superposition: SuperpositionOptions {
                    periodic_policy: kekule::alignment::PeriodicAlignmentPolicy::RejectPeriodic,
                    ..Default::default()
                },
                ..Default::default()
            }
        ),
        Err(RmsdError::Alignment {
            frame: 0,
            source: AlignmentError::PeriodicCoordinates { .. }
        })
    ));
}
