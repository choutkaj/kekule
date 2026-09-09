use std::io::Cursor;
use std::sync::Arc;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::topology::Topology;
use kekule::units::{Quantity, NANOMETER, PICOSECOND};
use kekule_traj::io::xtc::{XtcMagic, XtcReadOptions, XtcReader, XtcWriteOptions, XtcWriter};
use kekule_traj::io::TrajectoryIoLimits;
use kekule_traj::{
    FrameBuffer, SeekableTrajectoryReader, TrajectoryCodecErrorKind, TrajectoryError,
    TrajectoryReader, TrajectoryWriter,
};
use sha2::{Digest, Sha256};

mod support;
use support::{
    buffer_snapshot, codec_kind, linear_carbon_topology as topology, x_coordinates as x_values,
    GuardedCursor, RestoreSeekFailure,
};

fn source_frame(topology: &Arc<Topology>, shift: f64, step: u64) -> FrameBuffer {
    let mut frame = FrameBuffer::new(Arc::clone(topology));
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

fn encoded(atom_count: usize, magic: XtcMagic) -> (Arc<Topology>, Vec<u8>) {
    let topology = topology(atom_count);
    let options = XtcWriteOptions::default()
        .with_magic(magic)
        .with_precision(1000.0)
        .unwrap();
    let mut writer = XtcWriter::new(
        Cursor::new(Vec::new()),
        Arc::clone(&topology),
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
    topology: &Arc<Topology>,
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
        Arc::clone(topology),
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
fn xtc_rejects_cells_degenerate_after_f32_rounding_before_appending_bytes() {
    let topology = topology(4);
    let mut writer = XtcWriter::new(
        Cursor::new(Vec::new()),
        Arc::clone(&topology),
        XtcWriteOptions::default(),
        "box-rounding.xtc",
    )
    .unwrap();
    let mut frame = source_frame(&topology, 0.0, 0);
    let valid_cell = frame.cell().copied();
    writer.write_frame(frame.frame_view()).unwrap();
    let before = writer.writer().clone();
    frame.set_cell(Some(
        PeriodicCell::new(
            Quantity::new(
                [
                    Vector3::new(1.0, 1.0, 0.0),
                    Vector3::new(1.0, 1.0 + 1e-8, 0.0),
                    Vector3::new(0.0, 0.0, 1.0),
                ],
                NANOMETER,
            ),
            [true; 3],
        )
        .unwrap(),
    ));
    assert_eq!(
        codec_kind(&writer.write_frame(frame.frame_view()).unwrap_err()),
        Some(TrajectoryCodecErrorKind::InvalidFrame)
    );
    assert_eq!(writer.writer(), &before);
    frame.set_cell(valid_cell);
    frame.set_step(Some(1));
    writer.write_frame(frame.frame_view()).unwrap();
    let mut reader = XtcReader::new(
        Cursor::new(writer.finish().unwrap().into_inner()),
        topology,
        XtcReadOptions::default(),
    )
    .unwrap();
    let mut destination = reader.frame_buffer();
    for step in [0, 1] {
        assert!(reader.read_next(&mut destination).unwrap());
        assert_eq!(destination.frame_view().step(), Some(step));
        assert!(destination.cell().is_some());
    }
    assert!(!reader.read_next(&mut destination).unwrap());
}

#[test]
fn xtc_aggregate_scratch_is_bounded_for_small_compressed_and_indexed_reads() {
    for atom_count in [4, 12] {
        for magic in [XtcMagic::Xtc1995, XtcMagic::Xtc2023] {
            let topology = topology(atom_count);
            let mut combined = Vec::new();
            let mut raw_sizes = Vec::new();
            for step in 0..3 {
                let mut frame = source_frame(&topology, step as f64 * 0.01, step);
                if step == 1 {
                    frame
                        .set_positions(Quantity::new(
                            (0..atom_count)
                                .map(|i| {
                                    Point3::new(
                                        (i * 37 % 101) as f64,
                                        (i * 53 % 97) as f64,
                                        (i * 71 % 89) as f64,
                                    )
                                })
                                .collect::<Vec<_>>(),
                            NANOMETER,
                        ))
                        .unwrap();
                }
                let mut writer = XtcWriter::new(
                    Cursor::new(Vec::new()),
                    Arc::clone(&topology),
                    XtcWriteOptions::default().with_magic(magic),
                    "scratch.xtc",
                )
                .unwrap();
                writer.write_frame(frame.frame_view()).unwrap();
                let bytes = writer.finish().unwrap().into_inner();
                let raw_size = if atom_count <= 9 {
                    atom_count * 12
                } else {
                    match magic {
                        XtcMagic::Xtc1995 => {
                            u32::from_be_bytes(bytes[88..92].try_into().unwrap()) as usize
                        }
                        XtcMagic::Xtc2023 => {
                            usize::try_from(u64::from_be_bytes(bytes[88..96].try_into().unwrap()))
                                .unwrap()
                        }
                        _ => unreachable!("covered XTC magic variants"),
                    }
                };
                raw_sizes.push(raw_size);
                combined.extend(bytes);
            }
            let dense_bytes =
                atom_count * (std::mem::size_of::<Point3>() + 3 * std::mem::size_of::<f32>());
            let first_total = dense_bytes + raw_sizes[0];
            let total = dense_bytes + raw_sizes.iter().max().unwrap();
            let options = |limit| {
                XtcReadOptions::default().with_limits(TrajectoryIoLimits {
                    max_scratch_bytes: limit,
                    ..TrajectoryIoLimits::default()
                })
            };
            let open = |limit| {
                XtcReader::new(
                    Cursor::new(combined.clone()),
                    Arc::clone(&topology),
                    options(limit),
                )
            };
            let payload_start = if atom_count <= 9 {
                56
            } else if magic == XtcMagic::Xtc1995 {
                92
            } else {
                96
            };
            let (stream, control) = GuardedCursor::new(combined.clone(), payload_start);
            let error = XtcReader::new(stream, Arc::clone(&topology), options(first_total - 1))
                .err()
                .unwrap();
            assert_eq!(
                codec_kind(&error),
                Some(TrajectoryCodecErrorKind::ResourceLimitExceeded)
            );
            assert!(!control.violated(), "over-budget payload must not be read");

            let mut sequential = open(total).unwrap();
            let mut destination = sequential.frame_buffer();
            let mut expected = Vec::new();
            while sequential.read_next(&mut destination).unwrap() {
                expected.push(destination.frame_view().to_frame());
            }
            if atom_count > 9 {
                assert!(
                    raw_sizes[1] > raw_sizes[0],
                    "the second compressed payload exercises capacity growth"
                );
                let mut bounded = open(total - 1).unwrap();
                assert!(bounded.read_next(&mut destination).unwrap());
                let before = buffer_snapshot(&destination);
                assert_eq!(
                    codec_kind(&bounded.read_next(&mut destination).unwrap_err()),
                    Some(TrajectoryCodecErrorKind::ResourceLimitExceeded)
                );
                assert_eq!(buffer_snapshot(&destination), before);
                assert_eq!(
                    codec_kind(&open(total - 1).unwrap().into_indexed().err().unwrap()),
                    Some(TrajectoryCodecErrorKind::ResourceLimitExceeded)
                );
            }
            let mut indexed = open(total).unwrap().into_indexed().unwrap();
            // A smaller frame follows the largest payload, so its retained raw
            // capacity still counts. Random reads must not duplicate any array
            // or reuse cached coordinates for the wrong pending frame.
            for (random, next) in [(1, 0), (0, 1), (1, 2)] {
                indexed.read_frame(random as u64, &mut destination).unwrap();
                assert_eq!(destination.frame_view().to_frame(), expected[random]);
                assert!(indexed.read_next(&mut destination).unwrap());
                assert_eq!(destination.frame_view().to_frame(), expected[next]);
            }
            indexed.read_frame(0, &mut destination).unwrap();
            assert!(!indexed.read_next(&mut destination).unwrap());
        }
    }
}

#[test]
fn indexed_xtc_rejects_changed_metadata_when_refreshing_a_pending_frame() {
    let (topology, bytes) = encoded(4, XtcMagic::Xtc1995);
    let stream = support::SharedCursor::new(bytes);
    let mut reader = XtcReader::new(stream.clone(), topology, XtcReadOptions::default())
        .unwrap()
        .into_indexed()
        .unwrap();
    let mut destination = reader.frame_buffer();
    // Indexing cached the first frame's original box. Random access replaces
    // its coordinate cache; the pending sequential frame must be revalidated.
    reader.read_frame(1, &mut destination).unwrap();
    stream.overwrite(16, &3.0_f32.to_be_bytes());
    let before = buffer_snapshot(&destination);
    assert!(reader.read_next(&mut destination).is_err());
    assert_eq!(buffer_snapshot(&destination), before);
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
            Arc::clone(&topology),
            XtcReadOptions::default()
                .with_limits(TrajectoryIoLimits::default())
                .with_source_label("memory.xtc"),
        )
        .unwrap();
        support::assert_rejects_unrelated_buffer(&mut reader);
        let mut destination = FrameBuffer::new(Arc::clone(&topology));
        let pointer = destination.positions().values().value().as_ptr();
        assert!(reader.read_next(&mut destination).unwrap());
        let expected = (0..atom_count)
            .map(|index| 0.1 * index as f64)
            .collect::<Vec<_>>();
        assert_x_close(&destination, &expected, 0.0011);
        assert_eq!(destination.frame_view().step(), Some(4));
        assert_eq!(destination.frame_view().time().unwrap().value(), &1.0);
        assert!(destination.cell().is_some());
        assert!(reader.read_next(&mut destination).unwrap());
        let expected = (0..atom_count)
            .map(|index| 0.1 * index as f64 + 0.01)
            .collect::<Vec<_>>();
        assert_x_close(&destination, &expected, 0.0011);
        assert_eq!(destination.positions().values().value().as_ptr(), pointer);
        assert!(!reader.read_next(&mut destination).unwrap());

        let mut indexed = XtcReader::new(
            Cursor::new(bytes),
            Arc::clone(&topology),
            XtcReadOptions::default()
                .with_limits(TrajectoryIoLimits::default())
                .with_source_label("memory.xtc"),
        )
        .unwrap()
        .into_indexed()
        .unwrap();
        assert_eq!(indexed.frame_count(), Some(2));
        indexed.read_frame(1, &mut destination).unwrap();
        assert_x_close(&destination, &expected, 0.0011);
        assert!(indexed.read_next(&mut destination).unwrap());
        let expected = (0..atom_count)
            .map(|index| 0.1 * index as f64)
            .collect::<Vec<_>>();
        assert_x_close(&destination, &expected, 0.0011);
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
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(limits.clone())
            .with_source_label("exact-limit.xtc"),
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    assert!(reader.read_next(&mut buffer).unwrap());
    assert!(reader.read_next(&mut buffer).unwrap());
    assert!(!reader.read_next(&mut buffer).unwrap());

    let indexed = XtcReader::new(
        Cursor::new(bytes),
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(limits)
            .with_source_label("exact-index-limit.xtc"),
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
            Arc::clone(&topology),
            XtcReadOptions::default()
                .with_limits(TrajectoryIoLimits::default())
                .with_source_label(label),
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
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("negative-step.xtc"),
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
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("maximum-step.xtc"),
    )
    .unwrap();
    let mut destination = FrameBuffer::new(Arc::clone(&topology));
    assert!(reader.read_next(&mut destination).unwrap());
    assert_eq!(destination.frame_view().step(), Some(i32::MAX as u64));

    let mut maximum_count = valid.clone();
    maximum_count[4..8].copy_from_slice(&i32::MAX.to_be_bytes());
    maximum_count[52..56].copy_from_slice(&i32::MAX.to_be_bytes());
    let error = XtcReader::new(
        Cursor::new(maximum_count),
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(TrajectoryIoLimits {
                max_atoms: i32::MAX as usize,
                ..TrajectoryIoLimits::default()
            })
            .with_source_label("maximum-count.xtc"),
    )
    .err()
    .unwrap();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::InconsistentAtomCount)
    );

    let mut writer = XtcWriter::new(
        Cursor::new(Vec::new()),
        Arc::clone(&topology),
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
        Arc::clone(&topology),
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
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("restore-failure.xtc"),
    )
    .unwrap()
    .into_indexed()
    .unwrap();
    let mut destination = source_frame(&topology, 9.0, 99);
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
fn xtc_limits_probe_but_do_not_decode_or_consume_frame_n_plus_one() {
    let topology = topology(12);
    let options = XtcWriteOptions::default().with_precision(1000.0).unwrap();
    let mut writer = XtcWriter::new(
        Cursor::new(Vec::new()),
        Arc::clone(&topology),
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
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(TrajectoryIoLimits {
                max_frames: 1,
                ..TrajectoryIoLimits::default()
            })
            .with_source_label("guarded-sequential.xtc"),
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
        let error = XtcReader::new(
            stream,
            Arc::clone(&topology),
            XtcReadOptions::default()
                .with_limits(limits)
                .with_source_label("guarded-index.xtc"),
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
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("fuzz-underflow.xtc"),
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
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("count.xtc"),
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
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("precision.xtc"),
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
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("corrupt.xtc"),
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
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("truncated.xtc"),
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
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(limits)
            .with_source_label("limited.xtc"),
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
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("trailing.xtc"),
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
            Arc::clone(&topology),
            XtcReadOptions::default()
                .with_limits(TrajectoryIoLimits::default())
                .with_source_label("mixed-profile.xtc"),
        )
        .unwrap();
        let mut destination = FrameBuffer::new(Arc::clone(&topology));
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
        Arc::clone(&topology),
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
    frame.set_cell(source_frame(&topology, 0.0, 0).cell().copied());
    frame
        .insert_property(
            PropertyKey::new("unsupported").unwrap(),
            PropertyValue::Bool(true),
        )
        .unwrap();
    assert_eq!(
        codec_kind(&writer.write_frame(frame.frame_view()).unwrap_err()),
        Some(TrajectoryCodecErrorKind::UnsupportedField)
    );
    frame.clear_properties();
    frame
        .insert_bond_property_column(
            PropertyKey::new("conformational_entropy").unwrap(),
            PropertyColumn::Real {
                unit: NANOMETER,
                values: vec![Some(1.0); topology.bond_count()],
            },
        )
        .unwrap();
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
        Arc::clone(&topology),
        XtcReadOptions::default()
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("mdanalysis-2.9.0-twelve-atoms.xtc"),
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(topology);
    assert!(reader.read_next(&mut buffer).unwrap());
    let expected = (0..12).map(|index| index as f64 * 0.01).collect::<Vec<_>>();
    assert_x_close(&buffer, &expected, 0.0011);
    assert!(buffer.cell().is_some());
    assert_eq!(buffer.frame_view().step(), Some(0));
    assert!(reader.read_next(&mut buffer).unwrap());
    let expected = (0..12)
        .map(|index| index as f64 * 0.01 + 0.001)
        .collect::<Vec<_>>();
    assert_x_close(&buffer, &expected, 0.0011);
    assert_eq!(buffer.frame_view().step(), Some(1));
    assert!(!reader.read_next(&mut buffer).unwrap());
}
