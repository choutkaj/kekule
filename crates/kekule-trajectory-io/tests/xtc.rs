use std::io::Cursor;

use kekule::core::{Atom, Element, Molecule};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::small::SmallMolecule;
use kekule::topology::{MoleculeInstanceMetadata, Topology, TopologyBuilder};
use kekule::trajectory::{
    AtomOrderAssertion, FrameBuffer, SeekableTrajectoryReader, TrajectoryCodecErrorKind,
    TrajectoryError, TrajectoryReader, TrajectoryWriter,
};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};
use kekule_trajectory_io::xtc::{XtcMagic, XtcReadOptions, XtcReader, XtcWriteOptions, XtcWriter};
use kekule_trajectory_io::{TrajectoryIoLimits, TrajectoryTopologyBinding};
use sha2::{Digest, Sha256};

mod support;
use support::{buffer_snapshot, GuardedCursor, RestoreSeekFailure};

fn topology(atom_count: usize) -> Topology {
    let mut graph = Molecule::new();
    for _ in 0..atom_count {
        graph
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .unwrap();
    }
    let molecule = SmallMolecule::from_graph(graph);
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_small_molecule_definition(&molecule).unwrap();
    builder
        .add_instance(definition, MoleculeInstanceMetadata::default())
        .unwrap();
    builder.build().unwrap()
}

fn binding(topology: &Topology) -> TrajectoryTopologyBinding {
    TrajectoryTopologyBinding::new(
        topology.clone(),
        AtomOrderAssertion::assert_file_uses_topology_order(topology),
    )
    .unwrap()
}

fn codec_kind(error: &TrajectoryError) -> Option<TrajectoryCodecErrorKind> {
    match error {
        TrajectoryError::Codec(context) => Some(context.kind()),
        _ => None,
    }
}

fn source_frame(topology: &Topology, shift: f64, step: u64) -> FrameBuffer {
    let mut frame = FrameBuffer::new(topology.clone());
    let positions = (0..topology.atom_count())
        .map(|index| {
            let index = index as f64;
            Point3::new(
                0.1 * index + shift,
                0.2 * index + shift,
                0.3 * index + shift,
            )
        })
        .collect::<Vec<_>>();
    frame
        .set_positions(Quantity::new(positions, NANOMETER))
        .unwrap();
    frame.set_cell(Some(
        PeriodicCell::new(
            Quantity::new(
                [
                    Vector3::new(2.0, 0.0, 0.0),
                    Vector3::new(0.1, 2.1, 0.0),
                    Vector3::new(0.2, 0.3, 2.2),
                ],
                NANOMETER,
            ),
            [true; 3],
        )
        .unwrap(),
    ));
    frame
        .set_time(Some(Quantity::new(step as f64 * 0.25, PICOSECOND)))
        .unwrap();
    frame.set_step(Some(step));
    frame
}

fn x_values(buffer: &FrameBuffer) -> Vec<f64> {
    buffer
        .configuration()
        .positions()
        .values()
        .value()
        .iter()
        .map(|point| point.x)
        .collect()
}

fn assert_x_close(buffer: &FrameBuffer, expected: &[f64], tolerance: f64) {
    let actual = x_values(buffer);
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.into_iter().zip(expected) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "{actual} differs from {expected}"
        );
    }
}

fn encoded(atom_count: usize, magic: XtcMagic) -> (Topology, Vec<u8>) {
    let topology = topology(atom_count);
    let options = XtcWriteOptions::default()
        .with_magic(magic)
        .with_precision(1000.0)
        .unwrap();
    let mut writer = XtcWriter::new(
        Cursor::new(Vec::new()),
        topology.clone(),
        options,
        "memory.xtc",
    )
    .unwrap();
    let first = source_frame(&topology, 0.0, 4);
    let second = source_frame(&topology, 0.01, 5);
    writer.write_frame(first.frame_view()).unwrap();
    writer.write_frame(second.frame_view()).unwrap();
    (topology, writer.finish().unwrap().into_inner())
}

fn encoded_frame(
    topology: &Topology,
    magic: XtcMagic,
    precision: f32,
    shift: f64,
    step: u64,
) -> Vec<u8> {
    let options = XtcWriteOptions::default()
        .with_magic(magic)
        .with_precision(precision)
        .unwrap();
    let mut writer = XtcWriter::new(
        Cursor::new(Vec::new()),
        topology.clone(),
        options,
        "one-frame.xtc",
    )
    .unwrap();
    writer
        .write_frame(source_frame(topology, shift, step).frame_view())
        .unwrap();
    writer.finish().unwrap().into_inner()
}

#[test]
fn xtc_round_trips_small_and_compressed_frames_with_both_magic_variants() {
    for (atom_count, magic) in [
        (3, XtcMagic::Xtc1995),
        (12, XtcMagic::Xtc1995),
        (12, XtcMagic::Xtc2023),
    ] {
        let (topology, bytes) = encoded(atom_count, magic);
        let mut reader = XtcReader::new(
            Cursor::new(bytes.clone()),
            binding(&topology),
            XtcReadOptions::default(),
            TrajectoryIoLimits::default(),
            "memory.xtc",
        )
        .unwrap();
        let mut destination = FrameBuffer::new(topology.clone());
        let pointer = destination
            .configuration()
            .positions()
            .values()
            .value()
            .as_ptr();
        assert!(reader.read_next(&mut destination).unwrap());
        let expected = (0..atom_count)
            .map(|index| index as f64)
            .collect::<Vec<_>>();
        assert_x_close(&destination, &expected, 0.011);
        assert_eq!(destination.frame_view().step(), Some(4));
        assert_eq!(destination.frame_view().time().unwrap().value(), &1.0);
        assert!(destination.configuration().cell().is_some());
        assert!(reader.read_next(&mut destination).unwrap());
        let expected = (0..atom_count)
            .map(|index| index as f64 + 0.1)
            .collect::<Vec<_>>();
        assert_x_close(&destination, &expected, 0.011);
        assert_eq!(
            destination
                .configuration()
                .positions()
                .values()
                .value()
                .as_ptr(),
            pointer
        );
        assert!(!reader.read_next(&mut destination).unwrap());

        let mut indexed = XtcReader::new(
            Cursor::new(bytes),
            binding(&topology),
            XtcReadOptions::default(),
            TrajectoryIoLimits::default(),
            "memory.xtc",
        )
        .unwrap()
        .into_indexed()
        .unwrap();
        assert_eq!(indexed.frame_count(), Some(2));
        indexed.read_frame(1, &mut destination).unwrap();
        assert_x_close(&destination, &expected, 0.011);
        assert!(indexed.read_next(&mut destination).unwrap());
        let expected = (0..atom_count)
            .map(|index| index as f64)
            .collect::<Vec<_>>();
        assert_x_close(&destination, &expected, 0.011);
    }
}

#[test]
fn xtc_exact_frame_and_index_limits_still_allow_clean_eof() {
    let (topology, bytes) = encoded(12, XtcMagic::Xtc1995);
    let limits = TrajectoryIoLimits {
        max_frames: 2,
        max_index_entries: 2,
        max_index_bytes: 2 * std::mem::size_of::<u64>(),
        ..TrajectoryIoLimits::default()
    };
    let mut reader = XtcReader::new(
        Cursor::new(bytes.clone()),
        binding(&topology),
        XtcReadOptions::default(),
        limits.clone(),
        "exact-limit.xtc",
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(topology.clone());
    assert!(reader.read_next(&mut buffer).unwrap());
    assert!(reader.read_next(&mut buffer).unwrap());
    assert!(!reader.read_next(&mut buffer).unwrap());

    let indexed = XtcReader::new(
        Cursor::new(bytes),
        binding(&topology),
        XtcReadOptions::default(),
        limits,
        "exact-index-limit.xtc",
    )
    .unwrap()
    .into_indexed()
    .unwrap();
    assert_eq!(indexed.frame_count(), Some(2));
}

#[test]
fn xtc_signed_xdr_counts_and_steps_are_validated_before_private_adaptation() {
    let topology = topology(3);
    let valid = encoded_frame(&topology, XtcMagic::Xtc1995, 1000.0, 0.0, 4);

    for (range, expected_offset, label) in [
        (4..8, 4, "negative-count.xtc"),
        (52..56, 52, "negative-repeat.xtc"),
    ] {
        let mut negative = valid.clone();
        negative[range].copy_from_slice(&(-1_i32).to_be_bytes());
        let error = XtcReader::new(
            Cursor::new(negative),
            binding(&topology),
            XtcReadOptions::default(),
            TrajectoryIoLimits::default(),
            label,
        )
        .err()
        .unwrap();
        assert_eq!(
            codec_kind(&error),
            Some(TrajectoryCodecErrorKind::InvalidHeader)
        );
        let TrajectoryError::Codec(context) = error else {
            panic!("expected typed XTC count error");
        };
        assert_eq!(context.frame(), Some(0));
        assert_eq!(context.byte_offset(), Some(expected_offset));
    }

    let mut negative_step = valid.clone();
    negative_step[8..12].copy_from_slice(&(-1_i32).to_be_bytes());
    let error = XtcReader::new(
        Cursor::new(negative_step),
        binding(&topology),
        XtcReadOptions::default(),
        TrajectoryIoLimits::default(),
        "negative-step.xtc",
    )
    .err()
    .unwrap();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::NegativeOrUnrepresentableStep)
    );
    let TrajectoryError::Codec(context) = error else {
        panic!("expected typed XTC step error");
    };
    assert_eq!(context.frame(), Some(0));
    assert_eq!(context.byte_offset(), Some(8));

    let mut maximum_step = valid.clone();
    maximum_step[8..12].copy_from_slice(&i32::MAX.to_be_bytes());
    let mut reader = XtcReader::new(
        Cursor::new(maximum_step),
        binding(&topology),
        XtcReadOptions::default(),
        TrajectoryIoLimits::default(),
        "maximum-step.xtc",
    )
    .unwrap();
    let mut destination = FrameBuffer::new(topology.clone());
    assert!(reader.read_next(&mut destination).unwrap());
    assert_eq!(destination.frame_view().step(), Some(i32::MAX as u64));

    let mut maximum_count = valid.clone();
    maximum_count[4..8].copy_from_slice(&i32::MAX.to_be_bytes());
    maximum_count[52..56].copy_from_slice(&i32::MAX.to_be_bytes());
    let error = XtcReader::new(
        Cursor::new(maximum_count),
        binding(&topology),
        XtcReadOptions::default(),
        TrajectoryIoLimits {
            max_atoms: i32::MAX as usize,
            ..TrajectoryIoLimits::default()
        },
        "maximum-count.xtc",
    )
    .err()
    .unwrap();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::InconsistentAtomCount)
    );

    let mut writer = XtcWriter::new(
        Cursor::new(Vec::new()),
        topology.clone(),
        XtcWriteOptions::default(),
        "signed-writer.xtc",
    )
    .unwrap();
    let maximum = source_frame(&topology, 0.0, i32::MAX as u64);
    writer.write_frame(maximum.frame_view()).unwrap();
    let bytes = writer.finish().unwrap().into_inner();
    assert_eq!(
        i32::from_be_bytes(bytes[8..12].try_into().unwrap()),
        i32::MAX
    );

    let mut writer = XtcWriter::new(
        Cursor::new(Vec::new()),
        topology.clone(),
        XtcWriteOptions::default(),
        "overflow-writer.xtc",
    )
    .unwrap();
    let overflow = source_frame(&topology, 0.0, i32::MAX as u64 + 1);
    assert_eq!(
        codec_kind(&writer.write_frame(overflow.frame_view()).unwrap_err()),
        Some(TrajectoryCodecErrorKind::NegativeOrUnrepresentableStep)
    );
    assert!(writer.writer().get_ref().is_empty());
}

#[test]
fn indexed_xtc_restoration_failure_does_not_publish_or_change_destination() {
    let (topology, bytes) = encoded(12, XtcMagic::Xtc1995);
    let (stream, control) = RestoreSeekFailure::new(bytes);
    let mut indexed = XtcReader::new(
        stream,
        binding(&topology),
        XtcReadOptions::default(),
        TrajectoryIoLimits::default(),
        "restore-failure.xtc",
    )
    .unwrap()
    .into_indexed()
    .unwrap();
    let mut destination = source_frame(&topology, 9.0, 99);
    destination
        .props_mut()
        .insert("sentinel".into(), kekule::core::PropValue::Bool(true));
    let before = buffer_snapshot(&destination);
    control.arm_at_current_position();
    let error = indexed.read_frame(1, &mut destination).unwrap_err();
    assert!(matches!(error, TrajectoryError::Io(_)));
    assert_eq!(buffer_snapshot(&destination), before);
}

#[test]
fn xtc_limits_probe_but_do_not_decode_or_consume_frame_n_plus_one() {
    let topology = topology(12);
    let options = XtcWriteOptions::default().with_precision(1000.0).unwrap();
    let mut writer = XtcWriter::new(
        Cursor::new(Vec::new()),
        topology.clone(),
        options,
        "guarded.xtc",
    )
    .unwrap();
    writer
        .write_frame(source_frame(&topology, 0.0, 0).frame_view())
        .unwrap();
    let second_offset = writer.writer().position();
    writer
        .write_frame(source_frame(&topology, 0.01, 1).frame_view())
        .unwrap();
    let bytes = writer.finish().unwrap().into_inner();

    let (stream, control) = GuardedCursor::new(bytes.clone(), second_offset);
    let mut reader = XtcReader::new(
        stream,
        binding(&topology),
        XtcReadOptions::default(),
        TrajectoryIoLimits {
            max_frames: 1,
            ..TrajectoryIoLimits::default()
        },
        "guarded-sequential.xtc",
    )
    .unwrap();
    let mut destination = FrameBuffer::new(topology.clone());
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
        let error = XtcReader::new(
            stream,
            binding(&topology),
            XtcReadOptions::default(),
            limits,
            "guarded-index.xtc",
        )
        .unwrap()
        .into_indexed()
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
fn xtc_fuzz_regression_rejects_compressed_bitstream_underflow_without_panic() {
    let topology = topology(12);
    let fuzz_artifact = [
        0x00, 0x00, 0x07, 0xcb, 0x00, 0x00, 0x00, 0x0c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x3f, 0x0b,
        0x24, 0x21, 0x40, 0x01, 0xd2, 0x08, 0x00, 0x00, 0x00, 0x00, 0x3e, 0x44, 0x58, 0x2f, 0x3e,
        0xb0, 0x31, 0x2b, 0x40, 0x0a, 0x86, 0x3b, 0x00, 0x00, 0x00, 0x0c, 0x44, 0x7a, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x6e, 0x00, 0x00, 0x00, 0xdc, 0x00, 0x00, 0x01, 0x4a, 0x00, 0x00, 0x00, 0x12, 0x00, 0x00,
        0x00, 0x24, 0x70, 0x43, 0x17, 0x20, 0x13, 0x1a, 0xa1, 0x94, 0x86, 0x50, 0x26, 0x34, 0xc1,
        0x45, 0xc4, 0x00, 0x00, 0x00, 0x12, 0x00, 0x00, 0x00, 0x24, 0x70, 0x43, 0x17, 0x20, 0x13,
        0x1a, 0xa1, 0x94, 0x86, 0x50, 0x26, 0x34, 0xc1, 0x45, 0xc4, 0xa9, 0x5e, 0x88, 0x43, 0x62,
        0x74, 0xba, 0x48, 0xad, 0xe0, 0xbd, 0x96, 0xb2, 0x29, 0x66, 0xbe, 0x1a, 0x55, 0x7a, 0xaa,
        0x40,
    ];
    let error = XtcReader::new(
        Cursor::new(fuzz_artifact),
        binding(&topology),
        XtcReadOptions::default(),
        TrajectoryIoLimits::default(),
        "fuzz-underflow.xtc",
    )
    .err()
    .expect("fuzz artifact must be rejected during preflight");
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::CorruptCompressedData)
    );
}

#[test]
fn xtc_preflight_rejects_header_precision_truncation_corruption_and_limits() {
    let (topology, valid) = encoded(12, XtcMagic::Xtc1995);

    let mut repeated = valid.clone();
    repeated[52..56].copy_from_slice(&11_u32.to_be_bytes());
    let error = XtcReader::new(
        Cursor::new(repeated),
        binding(&topology),
        XtcReadOptions::default(),
        TrajectoryIoLimits::default(),
        "count.xtc",
    )
    .err()
    .unwrap();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::InconsistentAtomCount)
    );

    let mut precision = valid.clone();
    precision[56..60].copy_from_slice(&0_f32.to_be_bytes());
    let error = XtcReader::new(
        Cursor::new(precision),
        binding(&topology),
        XtcReadOptions::default(),
        TrajectoryIoLimits::default(),
        "precision.xtc",
    )
    .err()
    .unwrap();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::InvalidPrecision)
    );

    let mut small_index = valid.clone();
    small_index[84..88].copy_from_slice(&73_u32.to_be_bytes());
    let error = XtcReader::new(
        Cursor::new(small_index),
        binding(&topology),
        XtcReadOptions::default(),
        TrajectoryIoLimits::default(),
        "corrupt.xtc",
    )
    .err()
    .unwrap();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::CorruptCompressedData)
    );

    let mut truncated = valid.clone();
    truncated.pop();
    let error = XtcReader::new(
        Cursor::new(truncated),
        binding(&topology),
        XtcReadOptions::default(),
        TrajectoryIoLimits::default(),
        "truncated.xtc",
    )
    .unwrap()
    .into_indexed()
    .err()
    .unwrap();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::TruncatedRecord)
    );

    let limits = TrajectoryIoLimits {
        max_frame_bytes: 100,
        ..TrajectoryIoLimits::default()
    };
    let error = XtcReader::new(
        Cursor::new(valid),
        binding(&topology),
        XtcReadOptions::default(),
        limits,
        "limited.xtc",
    )
    .err()
    .unwrap();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::ResourceLimitExceeded)
    );
}

#[test]
fn xtc_rejects_trailing_compressed_data_and_mixed_file_profiles() {
    let topology = topology(12);
    let mut trailing = encoded_frame(&topology, XtcMagic::Xtc1995, 1000.0, 0.0, 0);
    let payload_bytes =
        usize::try_from(u32::from_be_bytes(trailing[88..92].try_into().unwrap())).unwrap();
    let payload_start = 92;
    trailing.splice(
        payload_start + payload_bytes..payload_start + payload_bytes,
        [0_u8; 4],
    );
    trailing[88..92].copy_from_slice(&u32::try_from(payload_bytes + 4).unwrap().to_be_bytes());
    let error = XtcReader::new(
        Cursor::new(trailing),
        binding(&topology),
        XtcReadOptions::default(),
        TrajectoryIoLimits::default(),
        "trailing.xtc",
    )
    .err()
    .unwrap();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::CorruptCompressedData)
    );
    let TrajectoryError::Codec(context) = &error else {
        panic!("expected typed XTC codec context");
    };
    assert_eq!(context.frame(), Some(0));
    assert_eq!(context.byte_offset(), Some(0));

    for second in [
        encoded_frame(&topology, XtcMagic::Xtc1995, 100.0, 0.01, 1),
        encoded_frame(&topology, XtcMagic::Xtc2023, 1000.0, 0.01, 1),
    ] {
        let mut mixed = encoded_frame(&topology, XtcMagic::Xtc1995, 1000.0, 0.0, 0);
        mixed.extend(second);
        let mut reader = XtcReader::new(
            Cursor::new(mixed),
            binding(&topology),
            XtcReadOptions::default(),
            TrajectoryIoLimits::default(),
            "mixed-profile.xtc",
        )
        .unwrap();
        let mut destination = FrameBuffer::new(topology.clone());
        assert!(reader.read_next(&mut destination).unwrap());
        let before = x_values(&destination);
        let error = reader.read_next(&mut destination).unwrap_err();
        assert_eq!(
            codec_kind(&error),
            Some(TrajectoryCodecErrorKind::InconsistentMetadata)
        );
        assert_eq!(x_values(&destination), before);
    }
}

#[test]
fn xtc_writer_rejects_unrepresentable_or_unpreserved_state() {
    let topology = topology(12);
    assert!(XtcWriteOptions::default().with_precision(0.0).is_err());
    let mut writer = XtcWriter::new(
        Cursor::new(Vec::new()),
        topology.clone(),
        XtcWriteOptions::default(),
        "strict.xtc",
    )
    .unwrap();
    let mut frame = source_frame(&topology, 0.0, 0);
    frame.set_cell(None);
    assert_eq!(
        codec_kind(&writer.write_frame(frame.frame_view()).unwrap_err()),
        Some(TrajectoryCodecErrorKind::InconsistentMetadata)
    );
    frame.set_cell(
        source_frame(&topology, 0.0, 0)
            .configuration()
            .cell()
            .copied(),
    );
    frame
        .props_mut()
        .insert("unsupported".into(), kekule::core::PropValue::Bool(true));
    assert_eq!(
        codec_kind(&writer.write_frame(frame.frame_view()).unwrap_err()),
        Some(TrajectoryCodecErrorKind::UnsupportedField)
    );
}

#[test]
fn empty_xtc_writer_is_rejected() {
    let topology = topology(3);
    let error = XtcWriter::new(
        Cursor::new(Vec::new()),
        topology,
        XtcWriteOptions::default(),
        "empty.xtc",
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
fn independently_generated_mdanalysis_xtc_matches_lossy_profile() {
    let topology = topology(12);
    let fixture = include_bytes!("fixtures/mdanalysis-2.9.0-twelve-atoms.xtc");
    let digest = Sha256::digest(fixture);
    let actual_digest = digest
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect::<String>();
    assert_eq!(
        actual_digest,
        "4983aea35f003d170bca0933942559a6d1968d6ec7030f7def0cc0fdfc757fec"
    );
    let mut reader = XtcReader::new(
        Cursor::new(fixture),
        binding(&topology),
        XtcReadOptions::default(),
        TrajectoryIoLimits::default(),
        "mdanalysis-2.9.0-twelve-atoms.xtc",
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(topology);
    assert!(reader.read_next(&mut buffer).unwrap());
    let expected = (0..12).map(|index| index as f64 * 0.1).collect::<Vec<_>>();
    assert_x_close(&buffer, &expected, 0.011);
    assert!(buffer.configuration().cell().is_some());
    assert_eq!(buffer.frame_view().step(), Some(0));
    assert!(reader.read_next(&mut buffer).unwrap());
    let expected = (0..12)
        .map(|index| index as f64 * 0.1 + 0.01)
        .collect::<Vec<_>>();
    assert_x_close(&buffer, &expected, 0.011);
    assert_eq!(buffer.frame_view().step(), Some(1));
    assert!(!reader.read_next(&mut buffer).unwrap());
}
