use std::sync::Arc;

use kekule::alignment::{AlignmentError, AlignmentOptions, AtomCorrespondence, Weighting};
use kekule::math::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{Forces, Place, Positions, Structure, StructureSet, Velocities};
use kekule::topology::AtomSelection;
use kekule::traj::analysis::FrameSuperposer;
use kekule::traj::periodic::{self, MoleculeImager, PeriodicError, TrajectoryUnwrapper};
use kekule::traj::{FrameBuffer, FrameView};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};

use crate::support::{self, linear_carbon_topology, time_ordered, Frame};

/// The x coordinate of atom 0 in each frame of [`source`].
const XS: [f64; 5] = [0.9, 0.1, 0.4, 0.8, 0.2];

/// Five frames of a three-carbon chain in a 1 nm cubic cell: atom 0 at
/// `(x, 0.2, 0.2)`, atom 1 at `(x + 0.2, 0.2, 0.2)` wrapped into the cell,
/// and atom 2 at `(x, 0.4, 0.2)`, for each `x` of [`XS`].
fn source() -> StructureSet {
    let mut frames = Vec::new();
    for (index, x) in XS.into_iter().enumerate() {
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

/// Asserts that every geometry of `set` has the positions `expected` gives
/// for its index, and the zero velocities and forces of [`source`].
fn assert_positions(set: &StructureSet, expected: impl Fn(usize) -> [[f64; 3]; 3]) {
    assert_eq!(set.len(), XS.len());
    for (index, geometry) in set.iter().enumerate() {
        let positions = geometry.positions().values();
        for (atom, (a, e)) in positions.value().iter().zip(expected(index)).enumerate() {
            let e = Point3::new(e[0], e[1], e[2]);
            assert!(
                (*a - e).norm() < 1.0e-12,
                "geometry {index} atom {atom}: {a:?} != {e:?}"
            );
        }
        assert_eq!(geometry.velocities(), Some(&Velocities::zeros(3)));
        assert_eq!(geometry.forces(), Some(&Forces::zeros(3)));
    }
}

const X: Vector3 = Vector3::new(1.0, 0.0, 0.0);
const Y: Vector3 = Vector3::new(0.0, 1.0, 0.0);
const Z: Vector3 = Vector3::new(0.0, 0.0, 1.0);

/// Asserts that geometry `index` of `set` has a cell of these vectors.
fn assert_cell(set: &StructureSet, index: usize, expected: [Vector3; 3]) {
    let vectors = set.get(index).unwrap().cell().unwrap().vectors();
    for (actual, expected) in vectors.value().iter().zip(expected) {
        assert!(
            (*actual - expected).norm() < 1.0e-12,
            "geometry {index}: {actual:?} != {expected:?}"
        );
    }
}

/// `set` made whole. The first atom of each frame stays and atom 1 takes
/// its image 0.2 nm along x from it, so no bond crosses the cell.
fn whole(set: &StructureSet) -> StructureSet {
    let mut whole = set.clone();
    periodic::make_molecules_whole(&mut whole).unwrap();
    assert_positions(&whole, |index| {
        let x = XS[index];
        [[x, 0.2, 0.2], [x + 0.2, 0.2, 0.2], [x, 0.4, 0.2]]
    });
    whole
}

/// Geometry `index` of `set` as a frame.
fn frame_of(set: &StructureSet, index: usize) -> FrameView<'_> {
    FrameView::from_geometry(set.get(index).unwrap()).unwrap()
}

#[test]
fn streaming_and_loaded_whole_image_and_unwrap_preserve_identical_complete_frames() {
    let set = source();
    let topology = set.shared_topology();
    let anchors = AtomSelection::all(Arc::clone(&topology));
    let imager = MoleculeImager::new(topology.clone());
    let expected_whole = whole(&set);
    let mut expected_image = set.clone();
    periodic::image_molecules(&mut expected_image, &anchors).unwrap();
    // Imaging moves the whole molecule's geometric center, (x + 1/15, 4/15,
    // 0.2), to the cell center in every frame.
    let (dx, dy) = (1.0 / 15.0, 4.0 / 15.0);
    assert_positions(&expected_image, |_| {
        [
            [0.5 - dx, 0.7 - dy, 0.5],
            [0.7 - dx, 0.7 - dy, 0.5],
            [0.5 - dx, 0.9 - dy, 0.5],
        ]
    });
    let mut expected_unwrap = expected_whole.clone();
    periodic::unwrap_trajectory(&mut expected_unwrap).unwrap();
    // Each atom moves by under half a cell between frames: atom 0 by +0.2
    // (crossing x = 1), +0.3, +0.4, and +0.4 (crossing again).
    let unwrapped_xs = [0.9, 1.1, 1.4, 1.8, 2.2];
    assert_positions(&expected_unwrap, |index| {
        let x = unwrapped_xs[index];
        [[x, 0.2, 0.2], [x + 0.2, 0.2, 0.2], [x, 0.4, 0.2]]
    });
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
    buffer.geometry_mut().set_cell(None);
    let missing = format!("{buffer:?}");
    assert_eq!(
        unwrapper.unwrap(11, &mut buffer),
        Err(PeriodicError::MissingCell { frame: 11 })
    );
    assert_eq!(format!("{buffer:?}"), missing);
    assert_eq!(unwrapper.last_frame_index(), Some(10));
    buffer.geometry_mut().set_cell(cell);
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
        .geometry_mut()
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
    buffer.geometry_mut().set_cell(Some(
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
    // Unwrapped atom by atom: atom 0 and atom 2 move as in the whole
    // molecule; atom 1 starts wrapped at 0.1 and moves by +0.2, +0.3, +0.4
    // (from 0.6 to 0.0, crossing x = 1), and +0.4.
    let atom_0_xs = [0.9, 1.1, 1.4, 1.8, 2.2];
    let atom_1_xs = [0.1, 0.3, 0.6, 1.0, 1.4];
    assert_positions(&unwrapped, |index| {
        let (x, x1) = (atom_0_xs[index], atom_1_xs[index]);
        [[x, 0.2, 0.2], [x1, 0.2, 0.2], [x, 0.4, 0.2]]
    });
    assert_eq!(*buffer, *unwrapped.get(1).unwrap().payload());
}

#[test]
fn superposition_streams_the_same_transforms_and_metadata_and_reports_source_indices() {
    let set = whole(&source());
    let topology = set.shared_topology();
    let selection = AtomSelection::all(Arc::clone(&topology));
    let mut expected = set.clone();
    let reports = expected
        .superpose(0, &selection, Default::default())
        .unwrap();
    // Every whole frame is frame 0 translated along x, so each superposes
    // exactly onto it, with zero RMSD.
    assert_positions(&expected, |_| {
        let x = XS[0];
        [[x, 0.2, 0.2], [x + 0.2, 0.2, 0.2], [x, 0.4, 0.2]]
    });
    for index in 0..expected.len() {
        assert!(reports.alignment(index).unwrap().rmsd().into_value() < 1.0e-12);
        // The fitted rotation is the identity to rounding, so the cell is
        // the unit cube still.
        assert_cell(&expected, index, [X, Y, Z]);
    }
    let reference = set.get(0).unwrap();
    let superposer = FrameSuperposer::new(&reference, &selection, Default::default());
    let mut buffer = FrameBuffer::new(topology.clone());
    for (index, frame) in set.iter().enumerate() {
        buffer
            .copy_from(FrameView::from_geometry(frame).unwrap())
            .unwrap();
        let report = superposer.superpose(index, &mut buffer).unwrap();
        assert_eq!(&report, reports.alignment(index).unwrap());
        assert_eq!(*buffer, *expected.get(index).unwrap().payload());
    }
    buffer
        .geometry_mut()
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
    let set = whole(&source());
    let moving = set.shared_topology();
    let first = set.get(0).unwrap().positions().values().value().to_vec();
    let turn = |p: Point3| Point3::new(-p.y + 2.0, p.x - 3.0, p.z + 1.0);
    // An independent reference with permuted atoms and an unpaired extra atom.
    let reference = Structure::new(
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
    // Every whole frame is frame 0 translated, and the paired reference
    // atoms are frame 0 turned rigidly, so every frame lands on them.
    assert_positions(&loaded, |_| {
        [turn(first[0]), turn(first[1]), turn(first[2])].map(|p| [p.x, p.y, p.z])
    });
    // The cell turns with the frame: a to +y, b to -x, c unchanged.
    for index in 0..loaded.len() {
        assert_cell(&loaded, index, [Y, -X, Z]);
    }
    let superposer = FrameSuperposer::new(&reference, &fit, options);
    let mut buffer = FrameBuffer::new(set.shared_topology());
    let mut index = 0;
    for frame in set.iter() {
        buffer
            .copy_from(FrameView::from_geometry(frame).unwrap())
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
