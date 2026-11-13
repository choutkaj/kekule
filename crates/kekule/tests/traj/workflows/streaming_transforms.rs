use std::sync::Arc;

use kekule::alignment::{AlignmentError, AlignmentOptions, AtomCorrespondence, Weighting};
use kekule::math::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{Forces, Model, Place, Positions, SampleSet, Velocities};
use kekule::topology::AtomSelection;
use kekule::traj::analysis::FrameSuperposer;
use kekule::traj::periodic::{self, MoleculeImager, PeriodicError, TrajectoryUnwrapper};
use kekule::traj::{FrameBuffer, FrameView};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};

use crate::support::{self, dynamics, linear_carbon_topology, time_ordered, Frame, Recorded};

fn source() -> SampleSet {
    let mut frames = Vec::new();
    for (index, x) in [0.9, 0.1, 0.4, 0.8, 0.2].into_iter().enumerate() {
        let mut frame = Frame::new(
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
            .properties_mut()
            .owner_mut()
            .insert(
                PropertyKey::new("label").unwrap(),
                PropertyValue::Int(index as i64),
            )
            .unwrap();
        frames.push(frame);
    }
    time_ordered(linear_carbon_topology(3), frames)
}

/// `set` made whole, the next result `recorded`.
fn whole(set: &SampleSet, recorded: &mut Recorded) -> SampleSet {
    let mut whole = set.clone();
    periodic::make_molecules_whole(&mut whole).unwrap();
    recorded.check(&dynamics(&whole));
    whole
}

/// Sample `index` of `set` as a frame.
fn frame_of(set: &SampleSet, index: usize) -> FrameView<'_> {
    FrameView::from_sample(set.get(index).unwrap()).unwrap()
}

#[test]
fn streaming_and_loaded_whole_image_and_unwrap_preserve_identical_complete_frames() {
    let mut recorded = Recorded::new(&[
        0xf5ab_a216_67c8_6875, // whole
        0x7242_7eee_5264_a0e3, // imaged
        0x73a8_6f7c_f4ee_e247, // whole, then unwrapped
    ]);
    let set = source();
    let topology = set.shared_topology();
    let anchors = AtomSelection::all(Arc::clone(&topology));
    let imager = MoleculeImager::new(topology.clone());
    let expected_whole = whole(&set, &mut recorded);
    let mut expected_image = set.clone();
    periodic::image_molecules(&mut expected_image, &anchors).unwrap();
    recorded.check(&dynamics(&expected_image));
    let mut expected_unwrap = expected_whole.clone();
    periodic::unwrap_trajectory(&mut expected_unwrap).unwrap();
    recorded.check(&dynamics(&expected_unwrap));
    let mut unwrapper = TrajectoryUnwrapper::new(topology.clone());
    let mut buffer = FrameBuffer::new(set.shared_topology());
    for index in 0..set.len() {
        buffer.copy_from(frame_of(&set, index)).unwrap();
        let mut image = buffer.clone();
        imager.image(index, &mut image, &anchors).unwrap();
        assert_eq!(*image, *expected_image.get(index).unwrap().payload());
        imager.make_whole(index, &mut buffer).unwrap();
        assert_eq!(*buffer, *expected_whole.get(index).unwrap().payload());
        unwrapper.unwrap(index, &mut buffer).unwrap();
        assert_eq!(*buffer, *expected_unwrap.get(index).unwrap().payload());
        // The same unwrapper survives processing-chunk boundaries.
        assert_eq!(unwrapper.last_frame_index(), Some(index));
        buffer.copy_from(frame_of(&set, index)).unwrap();
        imager.image(index, &mut buffer, &anchors).unwrap();
        assert_eq!(*buffer, *image);
        assert!(Arc::ptr_eq(&buffer.shared_topology(), &topology));
    }
}

#[test]
fn unwrapper_rejects_skips_reordering_and_failures_without_advancing_and_supports_reset() {
    let set = source();
    let topology = set.shared_topology();
    let mut unwrapper = TrajectoryUnwrapper::new(topology.clone());
    let mut buffer = FrameBuffer::new(topology.clone());
    buffer.copy_from(frame_of(&set, 0)).unwrap();
    unwrapper.unwrap(10, &mut buffer).unwrap();
    buffer.copy_from(frame_of(&set, 1)).unwrap();
    let before = format!("{buffer:?}");
    for index in [10, 9, 12] {
        assert_eq!(
            unwrapper.unwrap(index, &mut buffer),
            Err(PeriodicError::NonSequentialFrame {
                previous: 10,
                frame: index
            })
        );
        assert_eq!(format!("{buffer:?}"), before);
        assert_eq!(unwrapper.last_frame_index(), Some(10));
    }
    let cell = buffer.cell().copied();
    buffer.sample_mut().set_cell(None);
    let missing = format!("{buffer:?}");
    assert_eq!(
        unwrapper.unwrap(11, &mut buffer),
        Err(PeriodicError::MissingCell { frame: 11 })
    );
    assert_eq!(format!("{buffer:?}"), missing);
    assert_eq!(unwrapper.last_frame_index(), Some(10));
    buffer.sample_mut().set_cell(cell);
    unwrapper.unwrap(11, &mut buffer).unwrap();
    assert_eq!(unwrapper.last_frame_index(), Some(11));
    unwrapper.reset();
    buffer.copy_from(frame_of(&set, 0)).unwrap();
    let before = (*buffer).clone();
    unwrapper.unwrap(100, &mut buffer).unwrap();
    assert_eq!((*buffer).clone(), before);
    let mut foreign = FrameBuffer::new(linear_carbon_topology(3));
    assert_eq!(
        unwrapper.unwrap(101, &mut foreign),
        Err(PeriodicError::TopologyMismatch { frame: 101 })
    );
    assert_eq!(unwrapper.last_frame_index(), Some(100));
    assert_eq!(
        MoleculeImager::new(topology).make_whole(8, &mut foreign),
        Err(PeriodicError::TopologyMismatch { frame: 8 })
    );
}

#[test]
fn unwrapping_checks_that_time_increases_transactionally() {
    let mut set = source();
    set.set_place(
        2,
        Some(Place::Time {
            time: Some(Quantity::new(-1.0, PICOSECOND)),
            step: Some(20),
        }),
    )
    .unwrap();
    let before = format!("{set:?}");
    assert_eq!(
        periodic::unwrap_trajectory(&mut set),
        Err(PeriodicError::NonMonotonicTime { frame: 2 })
    );
    assert_eq!(format!("{set:?}"), before);
}

#[test]
fn ambiguous_crossings_and_changed_axes_do_not_advance_streaming_state() {
    let mut recorded = Recorded::new(&[
        0xfc7d_d02b_e2ff_b45c, // unwrapped
    ]);
    let set = source();
    let topology = set.shared_topology();
    let mut unwrapper = TrajectoryUnwrapper::new(topology.clone());
    let mut buffer = FrameBuffer::new(topology);
    buffer.copy_from(frame_of(&set, 0)).unwrap();
    unwrapper.unwrap(0, &mut buffer).unwrap();
    buffer.copy_from(frame_of(&set, 1)).unwrap();
    let mut points = buffer.positions().values().value().to_vec();
    points[0].x = 1.4; // Half a cell from the previous x = 0.9.
    buffer
        .sample_mut()
        .set_positions(Quantity::new(points, NANOMETER))
        .unwrap();
    let before = (*buffer).clone();
    assert!(matches!(
        unwrapper.unwrap(1, &mut buffer),
        Err(PeriodicError::AmbiguousDisplacement {
            frame: 1,
            axis: 0,
            ..
        })
    ));
    assert_eq!((*buffer).clone(), before);
    assert_eq!(unwrapper.last_frame_index(), Some(0));
    buffer.copy_from(frame_of(&set, 1)).unwrap();
    buffer.sample_mut().set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER),
            [true, true, false],
        )
        .unwrap(),
    ));
    let before = (*buffer).clone();
    assert_eq!(
        unwrapper.unwrap(1, &mut buffer),
        Err(PeriodicError::PeriodicAxesChanged { frame: 1 })
    );
    assert_eq!((*buffer).clone(), before);
    assert_eq!(unwrapper.last_frame_index(), Some(0));
    buffer.copy_from(frame_of(&set, 1)).unwrap();
    unwrapper.unwrap(1, &mut buffer).unwrap();
    let mut unwrapped = set.clone();
    periodic::unwrap_trajectory(&mut unwrapped).unwrap();
    recorded.check(&dynamics(&unwrapped));
    assert_eq!(*buffer, *unwrapped.get(1).unwrap().payload());
}

#[test]
fn superposition_streams_the_same_transforms_and_metadata_and_reports_source_indices() {
    let mut recorded = Recorded::new(&[
        0xf5ab_a216_67c8_6875, // whole
        0x8b5f_d339_0a74_9dfc, // superposed
    ]);
    let set = whole(&source(), &mut recorded);
    let topology = set.shared_topology();
    let selection = AtomSelection::all(Arc::clone(&topology));
    let mut expected = set.clone();
    let reports = expected
        .superpose(0, &selection, Default::default())
        .unwrap();
    recorded.check(&dynamics(&expected));
    let reference = set.get(0).unwrap();
    let superposer = FrameSuperposer::new(&reference, &selection, Default::default());
    let mut buffer = FrameBuffer::new(topology.clone());
    for (index, frame) in set.iter().enumerate() {
        buffer
            .copy_from(FrameView::from_sample(frame).unwrap())
            .unwrap();
        let report = superposer.superpose(index, &mut buffer).unwrap();
        assert_eq!(&report, reports.alignment(index).unwrap());
        assert_eq!(*buffer, *expected.get(index).unwrap().payload());
    }
    buffer
        .sample_mut()
        .set_positions(Quantity::new([Point3::origin(); 3], NANOMETER))
        .unwrap();
    let before = format!("{buffer:?}");
    assert!(matches!(
        superposer.superpose(42, &mut buffer),
        Err(AlignmentError::Item { index: 42, .. })
    ));
    assert_eq!(format!("{buffer:?}"), before);
}

#[test]
fn correspondence_streaming_matches_loaded_superposition_and_keeps_failed_buffers() {
    let mut recorded = Recorded::new(&[
        0xf5ab_a216_67c8_6875, // whole
        0x84cc_499d_3ac3_1bba, // superposed
    ]);
    let set = whole(&source(), &mut recorded);
    let moving = set.shared_topology();
    let first = set.get(0).unwrap().positions().values().value().to_vec();
    let turn = |p: Point3| Point3::new(-p.y + 2.0, p.x - 3.0, p.z + 1.0);
    // An independent reference with permuted atoms and an unpaired extra atom.
    let reference = Model::new(
        support::linear_carbon_topology(4),
        Positions::new(Quantity::new(
            [
                turn(first[2]),
                turn(first[0]),
                Point3::new(9.0, 9.0, 9.0),
                turn(first[1]),
            ],
            NANOMETER,
        ))
        .unwrap(),
    )
    .unwrap();
    let reference_ids = reference.topology().atom_ids();
    let fit = AtomCorrespondence::from_pairs(
        Arc::clone(&moving),
        reference.shared_topology(),
        [(0, 1), (1, 3), (2, 0)].map(|(a, b)| (moving.atom_ids()[a], reference_ids[b])),
    )
    .unwrap();
    let options = AlignmentOptions::default().with_weighting(Weighting::Explicit(&[4.0, 1.0, 2.0]));
    let mut loaded = set.clone();
    loaded.superpose(&reference, &fit, options).unwrap();
    recorded.check(&dynamics(&loaded));
    let superposer = FrameSuperposer::new(&reference, &fit, options);
    let mut buffer = FrameBuffer::new(set.shared_topology());
    let mut index = 0;
    for frame in set.iter() {
        buffer
            .copy_from(FrameView::from_sample(frame).unwrap())
            .unwrap();
        let source = (*buffer).clone();
        let (time, step) = (buffer.time(), buffer.step());
        superposer.superpose(index, &mut buffer).unwrap();
        assert_eq!(*buffer, *loaded.get(index).unwrap().payload());
        assert_eq!(buffer.properties(), source.properties());
        assert_eq!((buffer.time(), buffer.step()), (time, step));
        assert!(time.is_some() && step.is_some());
        index += 1;
    }
    assert_eq!(index, set.len());
    for (actual, expected) in loaded.get(0).unwrap().positions().values().value()[..2]
        .iter()
        .zip([turn(first[0]), turn(first[1])])
    {
        assert!((*actual - expected).norm() < 1.0e-12);
    }

    let mut collapsed = FrameBuffer::new(moving);
    let before = format!("{collapsed:?}");
    assert!(matches!(
        superposer.superpose(7, &mut collapsed),
        Err(AlignmentError::Item { index: 7, source })
            if matches!(*source, AlignmentError::DegenerateGeometry { .. })
    ));
    assert_eq!(format!("{collapsed:?}"), before);
}
