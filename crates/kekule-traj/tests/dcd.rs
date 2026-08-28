use std::io::Cursor;
use std::sync::Arc;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::topology::Topology;
use kekule::units::{Quantity, ANGSTROM, PICOSECOND};
use kekule_traj::io::dcd::{
    DcdEndian, DcdReadOptions, DcdReader, DcdTimePolicy, DcdWriteOptions, DcdWriter,
};
use kekule_traj::io::TrajectoryIoLimits;
use kekule_traj::{
    FrameBuffer, SeekableTrajectoryReader, TrajectoryCodecErrorKind, TrajectoryError,
    TrajectoryReader, TrajectoryWriter,
};
use sha2::{Digest, Sha256};

mod support;
use support::{
    binding, buffer_snapshot, codec_kind, topology as build_topology, x_coordinates as xs,
    GuardedCursor, NoBackwardSeekCursor, RestoreSeekFailure,
};

fn topology() -> Arc<Topology> {
    build_topology(&["C", "H", "O"], &[(0, 1), (0, 2)])
}

fn assert_xs_close(buffer: &FrameBuffer, expected: &[f64]) {
    let actual = xs(buffer);
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.into_iter().zip(expected) {
        assert!((actual - expected).abs() < 1.0e-12);
    }
}

fn set_frame(
    buffer: &mut FrameBuffer,
    coordinates: [[f64; 3]; 3],
    step: u64,
    time: Option<f64>,
    cell: Option<PeriodicCell>,
) {
    buffer
        .set_positions(Quantity::new(
            coordinates.map(|[x, y, z]| Point3::new(x, y, z)),
            ANGSTROM,
        ))
        .unwrap();
    buffer.set_step(Some(step));
    buffer
        .set_time(time.map(|value| Quantity::new(value, PICOSECOND)))
        .unwrap();
    buffer.set_cell(cell);
}

#[test]
fn canonical_dcd_round_trips_both_endians_cells_steps_and_explicit_time() {
    for endian in [DcdEndian::Little, DcdEndian::Big] {
        let topology = topology();
        let options = DcdWriteOptions::default()
            .with_endian(endian)
            .with_cells(true)
            .with_step_sequence(10, 2)
            .with_header_delta(0.5, DcdTimePolicy::HeaderDelta { unit: PICOSECOND });
        let mut writer = DcdWriter::new(
            Cursor::new(Vec::new()),
            Arc::clone(&topology),
            options,
            "memory.dcd",
        )
        .unwrap();
        let cell = PeriodicCell::new(
            Quantity::new(
                [
                    Vector3::new(10.0, 0.0, 0.0),
                    Vector3::new(2.0, 11.0, 0.0),
                    Vector3::new(1.0, 3.0, 12.0),
                ],
                ANGSTROM,
            ),
            [true; 3],
        )
        .unwrap();
        let mut source = FrameBuffer::new(Arc::clone(&topology));
        set_frame(
            &mut source,
            [[0.0, 1.0, 2.0], [3.0, 4.0, 5.0], [6.0, 7.0, 8.0]],
            10,
            Some(5.0),
            Some(cell),
        );
        writer.write_frame(source.frame_view()).unwrap();
        set_frame(
            &mut source,
            [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
            12,
            Some(6.0),
            Some(cell),
        );
        writer.write_frame(source.frame_view()).unwrap();
        let bytes = writer.finish().unwrap().into_inner();

        let options = DcdReadOptions::default()
            .with_time_policy(DcdTimePolicy::HeaderDelta { unit: PICOSECOND });
        let mut reader = DcdReader::new(
            Cursor::new(bytes.clone()),
            binding(&topology),
            options,
            TrajectoryIoLimits::default(),
            "memory.dcd",
        )
        .unwrap();
        let mut destination = FrameBuffer::new(Arc::clone(&topology));
        assert!(reader.read_next(&mut destination).unwrap());
        assert_xs_close(&destination, &[0.0, 0.3, 0.6]);
        assert_eq!(destination.frame_view().step(), Some(10));
        assert_eq!(destination.frame_view().time().unwrap().value(), &5.0);
        assert!(destination.cell().is_some());
        assert!(reader.read_next(&mut destination).unwrap());
        assert_xs_close(&destination, &[0.1, 0.4, 0.7]);
        assert_eq!(destination.frame_view().step(), Some(12));
        assert!(!reader.read_next(&mut destination).unwrap());

        let mut indexed = DcdReader::new(
            Cursor::new(bytes),
            binding(&topology),
            options,
            TrajectoryIoLimits::default(),
            "memory.dcd",
        )
        .unwrap()
        .to_indexed()
        .unwrap();
        assert_eq!(indexed.frame_count(), Some(2));
        indexed.read_frame(1, &mut destination).unwrap();
        assert_xs_close(&destination, &[0.1, 0.4, 0.7]);
        assert!(indexed.read_next(&mut destination).unwrap());
        assert_xs_close(&destination, &[0.0, 0.3, 0.6]);
    }
}

#[test]
fn fixed_atom_dcd_reconstructs_complete_frames_and_random_access() {
    let topology = topology();
    let bytes = fixed_atom_fixture(DcdEndian::Little);
    let reader = DcdReader::new(
        Cursor::new(bytes),
        binding(&topology),
        DcdReadOptions::default(),
        TrajectoryIoLimits::default(),
        "fixed.dcd",
    )
    .unwrap();
    let mut reader = reader.to_indexed().unwrap();
    let mut buffer = FrameBuffer::new(topology);
    reader.read_frame(1, &mut buffer).unwrap();
    assert_xs_close(&buffer, &[0.0, 1.0, 0.2]);
    assert_eq!(buffer.frame_view().step(), Some(1));
    reader.read_frame(0, &mut buffer).unwrap();
    assert_xs_close(&buffer, &[0.0, 0.1, 0.2]);
}

#[test]
fn dcd_exact_frame_and_index_limits_still_allow_clean_eof() {
    let topology = topology();
    let bytes = fixed_atom_fixture(DcdEndian::Little);
    let limits = TrajectoryIoLimits {
        max_frames: 2,
        max_index_entries: 2,
        max_index_bytes: 2 * std::mem::size_of::<u64>(),
        ..TrajectoryIoLimits::default()
    };
    let mut reader = DcdReader::new(
        Cursor::new(bytes.clone()),
        binding(&topology),
        DcdReadOptions::default(),
        limits.clone(),
        "exact-limit.dcd",
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    assert!(reader.read_next(&mut buffer).unwrap());
    assert!(reader.read_next(&mut buffer).unwrap());
    assert!(!reader.read_next(&mut buffer).unwrap());

    let indexed = DcdReader::new(
        Cursor::new(bytes),
        binding(&topology),
        DcdReadOptions::default(),
        limits,
        "exact-index-limit.dcd",
    )
    .unwrap()
    .to_indexed()
    .unwrap();
    assert_eq!(indexed.frame_count(), Some(2));
}

#[test]
fn dcd_declared_frame_count_is_strict_in_sequential_and_indexed_modes() {
    let topology = topology();
    let valid = fixed_atom_fixture(DcdEndian::Little);
    for (declared, frames_before_error) in [(1_i32, 1_usize), (3_i32, 2_usize)] {
        let mut bytes = valid.clone();
        bytes[8..12].copy_from_slice(&declared.to_le_bytes());

        let mut sequential = DcdReader::new(
            Cursor::new(bytes.clone()),
            binding(&topology),
            DcdReadOptions::default(),
            TrajectoryIoLimits::default(),
            "declared-sequential.dcd",
        )
        .unwrap();
        let mut destination = FrameBuffer::new(Arc::clone(&topology));
        for _ in 0..frames_before_error {
            assert!(sequential.read_next(&mut destination).unwrap());
        }
        assert_eq!(
            codec_kind(&sequential.read_next(&mut destination).unwrap_err()),
            Some(TrajectoryCodecErrorKind::InconsistentMetadata)
        );

        let error = DcdReader::new(
            Cursor::new(bytes),
            binding(&topology),
            DcdReadOptions::default(),
            TrajectoryIoLimits::default(),
            "declared-indexed.dcd",
        )
        .unwrap()
        .to_indexed()
        .err()
        .unwrap();
        assert_eq!(
            codec_kind(&error),
            Some(TrajectoryCodecErrorKind::InconsistentMetadata)
        );
    }
}

#[test]
fn indexed_dcd_restoration_failure_does_not_publish_or_change_destination() {
    let topology = topology();
    let (stream, control) = RestoreSeekFailure::new(fixed_atom_fixture(DcdEndian::Little));
    let mut indexed = DcdReader::new(
        stream,
        binding(&topology),
        DcdReadOptions::default(),
        TrajectoryIoLimits::default(),
        "restore-failure.dcd",
    )
    .unwrap()
    .to_indexed()
    .unwrap();
    let mut destination = FrameBuffer::new(topology);
    set_frame(
        &mut destination,
        [[9.0, 8.0, 7.0], [6.0, 5.0, 4.0], [3.0, 2.0, 1.0]],
        99,
        Some(12.5),
        None,
    );
    destination
        .insert_property(
            PropertyKey::new("sentinel").unwrap(),
            PropertyValue::Bool(true),
        )
        .unwrap();
    let before = buffer_snapshot(&destination);
    control.arm_at_current_position();
    let error = indexed.read_frame(1, &mut destination).unwrap_err();
    assert!(matches!(error, TrajectoryError::Io(_)));
    assert_eq!(buffer_snapshot(&destination), before);
}

#[test]
fn dcd_limits_probe_but_do_not_decode_or_consume_frame_n_plus_one() {
    let topology = topology();
    let options = DcdWriteOptions::default().with_step_sequence(0, 1);
    let mut writer = DcdWriter::new(
        Cursor::new(Vec::new()),
        Arc::clone(&topology),
        options,
        "guarded.dcd",
    )
    .unwrap();
    let mut frame = FrameBuffer::new(Arc::clone(&topology));
    set_frame(
        &mut frame,
        [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]],
        0,
        None,
        None,
    );
    writer.write_frame(frame.frame_view()).unwrap();
    let second_offset = writer.writer().position();
    set_frame(
        &mut frame,
        [[3.0, 3.0, 3.0], [4.0, 4.0, 4.0], [5.0, 5.0, 5.0]],
        1,
        None,
        None,
    );
    writer.write_frame(frame.frame_view()).unwrap();
    let bytes = writer.finish().unwrap().into_inner();

    let sequential_limits = TrajectoryIoLimits {
        max_frames: 1,
        ..TrajectoryIoLimits::default()
    };
    let (stream, control) = GuardedCursor::new(bytes.clone(), second_offset);
    let mut reader = DcdReader::new(
        stream,
        binding(&topology),
        DcdReadOptions::default(),
        sequential_limits,
        "guarded-sequential.dcd",
    )
    .unwrap();
    let mut destination = FrameBuffer::new(Arc::clone(&topology));
    assert!(reader.read_next(&mut destination).unwrap());
    assert_eq!(
        codec_kind(&reader.read_next(&mut destination).unwrap_err()),
        Some(TrajectoryCodecErrorKind::ResourceLimitExceeded)
    );
    assert!(!control.violated());
    assert_eq!(control.probed_bytes(), 1);

    for limits in [
        TrajectoryIoLimits {
            max_frames: 1,
            ..TrajectoryIoLimits::default()
        },
        TrajectoryIoLimits {
            max_index_entries: 1,
            ..TrajectoryIoLimits::default()
        },
        TrajectoryIoLimits {
            max_index_bytes: std::mem::size_of::<u64>(),
            ..TrajectoryIoLimits::default()
        },
    ] {
        let (stream, control) = GuardedCursor::new(bytes.clone(), second_offset);
        let error = DcdReader::new(
            stream,
            binding(&topology),
            DcdReadOptions::default(),
            limits,
            "guarded-index.dcd",
        )
        .unwrap()
        .to_indexed()
        .err()
        .unwrap();
        assert_eq!(
            codec_kind(&error),
            Some(TrajectoryCodecErrorKind::ResourceLimitExceeded)
        );
        assert!(!control.violated());
        assert_eq!(control.probed_bytes(), 1);
    }
}

#[test]
fn dcd_ordinary_frames_do_not_require_backward_eof_probe_seeks() {
    let topology = topology();
    let mut writer = DcdWriter::new(
        Cursor::new(Vec::new()),
        Arc::clone(&topology),
        DcdWriteOptions::default(),
        "no-probe.dcd",
    )
    .unwrap();
    let mut frame = FrameBuffer::new(Arc::clone(&topology));
    set_frame(&mut frame, [[0.0; 3], [1.0; 3], [2.0; 3]], 0, None, None);
    writer.write_frame(frame.frame_view()).unwrap();
    set_frame(&mut frame, [[3.0; 3], [4.0; 3], [5.0; 3]], 1, None, None);
    writer.write_frame(frame.frame_view()).unwrap();
    let bytes = writer.finish().unwrap().into_inner();

    let (stream, control) = NoBackwardSeekCursor::new(bytes);
    let mut reader = DcdReader::new(
        stream,
        binding(&topology),
        DcdReadOptions::default(),
        TrajectoryIoLimits::default(),
        "no-probe.dcd",
    )
    .unwrap();
    control.arm();
    let mut destination = FrameBuffer::new(topology);
    assert!(reader.read_next(&mut destination).unwrap());
    assert!(reader.read_next(&mut destination).unwrap());
}

#[test]
fn empty_dcd_writer_is_rejected() {
    let topology = topology();
    let error = DcdWriter::new(
        Cursor::new(Vec::new()),
        topology,
        DcdWriteOptions::default(),
        "empty.dcd",
    )
    .unwrap()
    .finish()
    .unwrap_err();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::InvalidFrame)
    );
}

#[test]
fn dcd_truncation_marker_counts_limits_and_publication_are_strict() {
    let topology = topology();
    let valid = fixed_atom_fixture(DcdEndian::Big);

    let mut wrong_count = valid.clone();
    wrong_count[8..12].copy_from_slice(&3_i32.to_be_bytes());
    let error = DcdReader::new(
        Cursor::new(wrong_count),
        binding(&topology),
        DcdReadOptions::default(),
        TrajectoryIoLimits::default(),
        "count.dcd",
    )
    .unwrap()
    .to_indexed()
    .err()
    .unwrap();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::InconsistentMetadata)
    );

    let mut bad_marker = valid.clone();
    *bad_marker.last_mut().unwrap() ^= 1;
    let mut reader = DcdReader::new(
        Cursor::new(bad_marker),
        binding(&topology),
        DcdReadOptions::default(),
        TrajectoryIoLimits::default(),
        "marker.dcd",
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    assert!(reader.read_next(&mut buffer).unwrap());
    let before = xs(&buffer);
    let error = reader.read_next(&mut buffer).unwrap_err();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::RecordMarkerMismatch)
    );
    let TrajectoryError::Codec(context) = &error else {
        panic!("expected typed DCD codec context");
    };
    assert_eq!(context.frame(), Some(1));
    assert!(context.byte_offset().is_some());
    assert_eq!(xs(&buffer), before);

    let mut truncated = valid.clone();
    truncated.pop();
    let error = DcdReader::new(
        Cursor::new(truncated),
        binding(&topology),
        DcdReadOptions::default(),
        TrajectoryIoLimits::default(),
        "truncated.dcd",
    )
    .unwrap()
    .to_indexed()
    .err()
    .unwrap();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::TruncatedRecord)
    );

    let limits = TrajectoryIoLimits {
        max_frame_bytes: 8,
        ..TrajectoryIoLimits::default()
    };
    let mut reader = DcdReader::new(
        Cursor::new(valid),
        binding(&topology),
        DcdReadOptions::default(),
        limits,
        "limited.dcd",
    )
    .unwrap();
    let error = reader.read_next(&mut buffer).unwrap_err();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::ResourceLimitExceeded)
    );
}

#[test]
fn dcd_writer_validates_the_complete_frame_before_writing_any_record() {
    let topology = topology();
    let options = DcdWriteOptions::default()
        .with_cells(true)
        .with_step_sequence(0, 1);
    let mut writer = DcdWriter::new(
        Cursor::new(Vec::new()),
        Arc::clone(&topology),
        options,
        "late-invalid.dcd",
    )
    .unwrap();
    let header_bytes = writer.writer().get_ref().len();
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(2.0, 2.0, 2.0), ANGSTROM),
        [true; 3],
    )
    .unwrap();
    let mut frame = FrameBuffer::new(Arc::clone(&topology));
    set_frame(
        &mut frame,
        [[f64::MAX, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]],
        0,
        None,
        Some(cell),
    );
    assert_eq!(
        codec_kind(&writer.write_frame(frame.frame_view()).unwrap_err()),
        Some(TrajectoryCodecErrorKind::InvalidFrame)
    );
    assert_eq!(writer.writer().get_ref().len(), header_bytes);

    let mut bond_annotated = FrameBuffer::new(Arc::clone(&topology));
    set_frame(
        &mut bond_annotated,
        [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 2.0, 2.0]],
        0,
        None,
        Some(cell),
    );
    bond_annotated
        .insert_bond_property_column(
            PropertyKey::new("conformational_entropy").unwrap(),
            kekule::properties::PropertyColumn::Real {
                unit: ANGSTROM,
                values: vec![Some(1.0); topology.bond_count()],
            },
        )
        .unwrap();
    assert_eq!(
        codec_kind(&writer.write_frame(bond_annotated.frame_view()).unwrap_err()),
        Some(TrajectoryCodecErrorKind::UnsupportedField)
    );
    assert_eq!(writer.writer().get_ref().len(), header_bytes);
}

#[test]
fn independently_generated_mdanalysis_fixture_is_interoperable() {
    let topology = topology();
    let fixture = include_bytes!("fixtures/mdanalysis-2.9.0-three-atoms.dcd");
    let digest = Sha256::digest(fixture);
    let actual_digest = digest
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect::<String>();
    assert_eq!(
        actual_digest,
        "baabffaca935dff654d1466fe43ba469fb0b1b84d077b6afc3c10f5208cc7d0e"
    );
    let mut reader = DcdReader::new(
        Cursor::new(fixture),
        binding(&topology),
        DcdReadOptions::default(),
        TrajectoryIoLimits::default(),
        "mdanalysis-2.9.0-three-atoms.dcd",
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(topology);
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_xs_close(&buffer, &[0.0, 0.3, 0.6]);
    assert!(buffer.cell().is_some());
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_xs_close(&buffer, &[0.1, 0.4, 0.7]);
    assert!(!reader.read_next(&mut buffer).unwrap());
}

fn fixed_atom_fixture(endian: DcdEndian) -> Vec<u8> {
    let i32_bytes = |value: i32| match endian {
        DcdEndian::Little => value.to_le_bytes(),
        DcdEndian::Big => value.to_be_bytes(),
        _ => unreachable!("test covers current DCD endian variants"),
    };
    let f32_bytes = |value: f32| match endian {
        DcdEndian::Little => value.to_le_bytes(),
        DcdEndian::Big => value.to_be_bytes(),
        _ => unreachable!("test covers current DCD endian variants"),
    };
    let mut bytes = Vec::new();
    let mut record = |payload: &[u8]| {
        bytes.extend(i32_bytes(payload.len() as i32));
        bytes.extend(payload);
        bytes.extend(i32_bytes(payload.len() as i32));
    };

    let mut header = [0_u8; 84];
    header[..4].copy_from_slice(b"CORD");
    let mut controls = [0_i32; 20];
    controls[0] = 2;
    controls[2] = 1;
    controls[8] = 2;
    controls[19] = 24;
    for (chunk, value) in header[4..].as_chunks_mut::<4>().0.iter_mut().zip(controls) {
        chunk.copy_from_slice(&i32_bytes(value));
    }
    header[40..44].copy_from_slice(&f32_bytes(1.0));
    record(&header);
    let mut title = [0_u8; 84];
    title[..4].copy_from_slice(&i32_bytes(1));
    title[4..9].copy_from_slice(b"fixed");
    record(&title);
    record(&i32_bytes(3));
    record(&i32_bytes(2));
    for values in [
        vec![0.0_f32, 1.0, 2.0],
        vec![0.0, 1.0, 2.0],
        vec![0.0, 1.0, 2.0],
        vec![10.0],
        vec![11.0],
        vec![12.0],
    ] {
        let payload = values.into_iter().flat_map(f32_bytes).collect::<Vec<_>>();
        record(&payload);
    }
    bytes
}
