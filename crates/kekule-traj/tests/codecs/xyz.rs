use std::io::{self, BufRead, BufReader, Cursor, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::{SystemTime, UNIX_EPOCH};

use kekule::diagnostics::{LimitExceeded, Resource, WriteOmission, WriteReport};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::topology::Topology;
use kekule::units::{Quantity, ANGSTROM, CANONICAL_VELOCITY_UNIT, NANOMETER, PICOSECOND};
use kekule_traj::io::xyz::{
    IndexedXyzReader, XyzReadOptions, XyzReader, XyzWriteOptions, XyzWriter,
};
use kekule_traj::io::{
    create_trajectory_writer, detect_trajectory_format, open_indexed_trajectory, open_trajectory,
    FieldAvailability, FormatDetectionEvidence, RandomAccessCapability, TrajectoryFormatHint,
    TrajectoryIoLimits, TrajectoryOpenOptions, TrajectoryWriteOptions,
};
use kekule_traj::{
    FrameBuffer, FrameData, SeekableTrajectoryReader, TrajectoryFormat, TrajectoryIoError,
    TrajectoryIoErrorKind, TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};
use sha2::{Digest, Sha256};

use crate::support;
use support::{
    codec_error, codec_kind, topology as build_topology, x_coordinates as point_xs, GuardedCursor,
};

const TWO_FRAMES: &str = "2\r\nfirst\r\nC 0.0 1.0 2.0\r\nH 3.0 4.0 5.0\r\n\
2\nsecond\nC 1.0 2.0 3.0\nH 4.0 5.0 6.0";

fn topology() -> Arc<Topology> {
    build_topology(&["C", "H"], &[(0, 1)])
}

fn water_topology() -> Arc<Topology> {
    build_topology(&["O", "H", "H"], &[(0, 1), (0, 2)])
}

fn assert_xs_close(buffer: &FrameBuffer, expected: &[f64]) {
    let actual = point_xs(buffer);
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.into_iter().zip(expected) {
        assert!((actual - expected).abs() < 1.0e-12);
    }
}

fn temporary_path(extension: Option<&str>) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let mut name = format!("kekule-traj-{}-{nonce}", std::process::id());
    if let Some(extension) = extension {
        name.push('.');
        name.push_str(extension);
    }
    std::env::temp_dir().join(name)
}

struct FailingDetectionReader {
    cursor: Cursor<Vec<u8>>,
}

impl Read for FailingDetectionReader {
    fn read(&mut self, _buffer: &mut [u8]) -> io::Result<usize> {
        Err(io::Error::new(
            io::ErrorKind::ConnectionReset,
            "injected detection read failure",
        ))
    }
}

impl Seek for FailingDetectionReader {
    fn seek(&mut self, position: SeekFrom) -> io::Result<u64> {
        self.cursor.seek(position)
    }
}

#[test]
fn sequential_xyz_is_transactional_reuses_positions_and_clears_stale_state() {
    let topology = topology();
    let mut reader = XyzReader::new(
        Cursor::new(TWO_FRAMES.as_bytes()),
        Arc::clone(&topology),
        XyzReadOptions::default()
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("memory.xyz"),
    )
    .unwrap();
    support::assert_rejects_unrelated_buffer(&mut reader);
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    let pointer = buffer.positions().values().value().as_ptr();
    buffer.frame_mut().set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER),
            [true; 3],
        )
        .unwrap(),
    ));
    buffer
        .set_velocities(Some(Quantity::new(
            [Vector3::new(1.0, 0.0, 0.0), Vector3::new(2.0, 0.0, 0.0)],
            CANONICAL_VELOCITY_UNIT,
        )))
        .unwrap();
    buffer
        .frame_mut()
        .set_time(Some(Quantity::new(1.0, PICOSECOND)))
        .unwrap();
    buffer.frame_mut().set_step(Some(7));
    buffer
        .frame_mut()
        .properties_mut()
        .owner_mut()
        .insert(
            PropertyKey::new("stale").unwrap(),
            PropertyValue::Bool(true),
        )
        .unwrap();

    assert!(reader.read_next(&mut buffer).unwrap());
    assert_xs_close(&buffer, &[0.0, 0.3]);
    assert_eq!(buffer.positions().values().value().as_ptr(), pointer);
    assert!(buffer.cell().is_none());
    assert!(buffer.frame_view().velocities().is_none());
    assert!(buffer.frame_view().time().is_none());
    assert!(buffer.frame_view().step().is_none());
    assert!(buffer.properties().is_empty());

    assert!(reader.read_next(&mut buffer).unwrap());
    assert_xs_close(&buffer, &[0.1, 0.4]);
    let before_eof = point_xs(&buffer);
    assert!(!reader.read_next(&mut buffer).unwrap());
    assert_eq!(point_xs(&buffer), before_eof);
}

#[test]
fn xyz_units_elements_limits_and_late_failures_are_explicit() {
    let topology = topology();
    let input = "2\nnm\nC 0.1 0.2 0.3\nH 0.4 0.5 0.6\n";
    let mut reader = XyzReader::new(
        Cursor::new(input.as_bytes()),
        Arc::clone(&topology),
        XyzReadOptions::default()
            .with_length_unit(NANOMETER)
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("nanometers.xyz"),
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    reader.read_next(&mut buffer).unwrap();
    assert_xs_close(&buffer, &[0.1, 0.4]);

    for (input, expected) in [
        (
            "2\nbad\nC 0 0 0\nO 1 1 1\n",
            TrajectoryIoErrorKind::InconsistentMetadata,
        ),
        ("2\nbad\nC 0 0 0\n", TrajectoryIoErrorKind::TruncatedRecord),
        (
            "2\nbad\nC 0 0 0\nH NaN 1 1\n",
            TrajectoryIoErrorKind::InvalidFrame,
        ),
        (
            "3\nbad\nC 0 0 0\nH 1 1 1\nH 2 2 2\n",
            TrajectoryIoErrorKind::InconsistentAtomCount,
        ),
    ] {
        let mut reader = XyzReader::new(
            Cursor::new(input.as_bytes()),
            Arc::clone(&topology),
            XyzReadOptions::default()
                .with_limits(TrajectoryIoLimits::default())
                .with_source_label("bad.xyz"),
        )
        .unwrap();
        let mut destination = FrameBuffer::new(Arc::clone(&topology));
        destination
            .replace_from_data(FrameData::new(Quantity::new(
                &[Point3::new(9.0, 0.0, 0.0), Point3::new(8.0, 0.0, 0.0)],
                ANGSTROM,
            )))
            .unwrap();
        let before = point_xs(&destination);
        let error = reader.read_next(&mut destination).unwrap_err();
        assert_eq!(codec_kind(&error), Some(expected));
        assert_eq!(point_xs(&destination), before);
    }

    let limits = TrajectoryIoLimits::default().with_max_text_line_bytes(8);
    let mut reader = XyzReader::new(
        Cursor::new("2\ncomment\nC 000000000 0 0\nH 1 1 1\n".as_bytes()),
        Arc::clone(&topology),
        XyzReadOptions::default()
            .with_limits(limits)
            .with_source_label("limited.xyz"),
    )
    .unwrap();
    assert_eq!(
        codec_kind(
            &reader
                .read_next(&mut FrameBuffer::new(topology))
                .unwrap_err()
        ),
        Some(TrajectoryIoErrorKind::LimitExceeded)
    );
}

#[test]
fn indexed_xyz_matches_sequential_and_random_reads_preserve_cursor() {
    let topology = topology();
    let mut reader = IndexedXyzReader::new(
        Cursor::new(TWO_FRAMES.as_bytes()),
        Arc::clone(&topology),
        XyzReadOptions::default()
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("indexed.xyz"),
    )
    .unwrap();
    assert_eq!(reader.frame_count(), Some(2));
    let mut buffer = FrameBuffer::new(topology);

    reader.read_frame(1, &mut buffer).unwrap();
    assert_xs_close(&buffer, &[0.1, 0.4]);
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_xs_close(&buffer, &[0.0, 0.3]);
    assert!(matches!(
        reader.read_frame(2, &mut buffer),
        Err(error) if error.kind() == TrajectoryIoErrorKind::FrameIndexOutOfRange && error.frame() == Some(2)
    ));
    assert_xs_close(&buffer, &[0.0, 0.3]);
}

#[test]
fn indexed_xyz_addresses_frames_from_the_stream_position_it_was_opened_at() {
    // Nanometre input keeps the decoded coordinates exact.
    let topology = topology();
    let bytes = b"2
first
C 0.0 1.0 2.0
H 4.0 5.0 6.0
2
second
C 1.0 2.0 3.0
H 4.0 5.0 6.0
2
third
C 2.0 3.0 4.0
H 4.0 5.0 6.0
";
    let expected = support::read_all(
        &mut XyzReader::new(
            Cursor::new(bytes),
            Arc::clone(&topology),
            XyzReadOptions::default().with_length_unit(NANOMETER),
        )
        .unwrap(),
    );
    assert_eq!(
        expected
            .iter()
            .map(|frame| frame.positions().values().value()[0].x)
            .collect::<Vec<_>>(),
        [0.0, 1.0, 2.0]
    );
    let mut indexed = IndexedXyzReader::new(
        support::cursor_after_leading_bytes(bytes),
        topology,
        XyzReadOptions::default().with_length_unit(NANOMETER),
    )
    .unwrap();
    support::assert_indexes_every_frame_from_its_start(&mut indexed, &expected);
}

#[test]
fn xyz_writer_is_strict_and_round_trips_without_owned_frames() {
    let topology = topology();
    let points = [Point3::new(1.25, 2.5, 3.75), Point3::new(4.0, 5.0, 6.0)];
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    buffer
        .replace_from_data(FrameData::new(Quantity::new(&points, ANGSTROM)))
        .unwrap();

    let mut writer = XyzWriter::new(
        Vec::new(),
        Arc::clone(&topology),
        XyzWriteOptions::default()
            .with_decimal_places(4)
            .with_comment("round trip"),
        "memory.xyz",
    )
    .unwrap();
    writer.write_frame(buffer.frame_view()).unwrap();
    writer.write_frame(buffer.frame_view()).unwrap();
    let bytes = writer.finish().unwrap();

    let mut reader = XyzReader::new(
        Cursor::new(bytes),
        Arc::clone(&topology),
        XyzReadOptions::default()
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("round-trip.xyz"),
    )
    .unwrap();
    let mut decoded = FrameBuffer::new(Arc::clone(&topology));
    assert!(reader.read_next(&mut decoded).unwrap());
    assert_eq!(
        decoded.positions().values().value(),
        buffer.positions().values().value()
    );
    assert!(reader.read_next(&mut decoded).unwrap());
    assert!(!reader.read_next(&mut decoded).unwrap());

    buffer.frame_mut().set_step(Some(1));
    let mut strict = XyzWriter::new(
        Vec::new(),
        topology,
        XyzWriteOptions::default(),
        "strict.xyz",
    )
    .unwrap();
    assert_eq!(
        codec_kind(&strict.write_frame(buffer.frame_view()).unwrap_err()),
        Some(TrajectoryIoErrorKind::UnsupportedField)
    );
    buffer.frame_mut().set_step(None);
    buffer
        .frame_mut()
        .properties_mut()
        .bonds_mut()
        .insert(
            PropertyKey::new("conformational_entropy").unwrap(),
            PropertyColumn::Real {
                unit: ANGSTROM,
                values: vec![Some(1.0)],
            },
        )
        .unwrap();
    assert_eq!(
        codec_kind(&strict.write_frame(buffer.frame_view()).unwrap_err()),
        Some(TrajectoryIoErrorKind::UnsupportedField)
    );
}

#[test]
fn path_detection_metadata_indexing_and_atomic_finish_are_bounded() {
    let topology = topology();
    let path = temporary_path(None);
    std::fs::write(&path, TWO_FRAMES).unwrap();
    let mut sequential = open_trajectory(
        &path,
        Arc::clone(&topology),
        TrajectoryOpenOptions::default(),
    )
    .unwrap();
    let report = sequential.open_report();
    assert_eq!(report.selected_format(), TrajectoryFormat::Xyz);
    assert!(report
        .detection_evidence()
        .contains(&FormatDetectionEvidence::MissingExtension));
    assert_eq!(sequential.metadata().atom_count(), 2);
    assert_eq!(
        sequential.metadata().fields().positions(),
        FieldAvailability::Required
    );
    assert_eq!(
        sequential.metadata().random_access(),
        RandomAccessCapability::SequentialOnly
    );
    assert!(sequential
        .read_next(&mut FrameBuffer::new(Arc::clone(&topology)))
        .unwrap());

    let indexed = open_indexed_trajectory(
        &path,
        Arc::clone(&topology),
        TrajectoryOpenOptions::default(),
    )
    .unwrap();
    assert_eq!(indexed.frame_count(), Some(2));
    assert_eq!(indexed.metadata().indexed_frame_count(), Some(2));
    assert_eq!(
        indexed.metadata().random_access(),
        RandomAccessCapability::Indexed
    );

    let mismatch = path.with_extension("dcd");
    std::fs::write(&mismatch, TWO_FRAMES).unwrap();
    assert_eq!(
        codec_kind(
            &open_trajectory(
                &mismatch,
                Arc::clone(&topology),
                TrajectoryOpenOptions::default()
            )
            .err()
            .unwrap()
        ),
        Some(TrajectoryIoErrorKind::FormatMismatch)
    );
    open_trajectory(
        &mismatch,
        Arc::clone(&topology),
        TrajectoryOpenOptions::default()
            .with_format_hint(TrajectoryFormatHint::Explicit(TrajectoryFormat::Xyz)),
    )
    .unwrap();

    let output = temporary_path(Some("xyz"));
    let points = [Point3::new(1.0, 2.0, 3.0), Point3::new(4.0, 5.0, 6.0)];
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    buffer
        .replace_from_data(FrameData::new(Quantity::new(&points, ANGSTROM)))
        .unwrap();
    let mut writer = create_trajectory_writer(
        &output,
        Arc::clone(&topology),
        TrajectoryWriteOptions::new(TrajectoryFormat::Xyz),
    )
    .unwrap();
    writer.write_frame(buffer.frame_view()).unwrap();
    assert!(!output.exists());
    writer.finish().unwrap();
    assert!(output.exists());
    open_trajectory(
        &output,
        Arc::clone(&topology),
        TrajectoryOpenOptions::default(),
    )
    .unwrap();

    let unfinished = temporary_path(Some("xyz"));
    {
        let mut writer = create_trajectory_writer(
            &unfinished,
            topology,
            TrajectoryWriteOptions::new(TrajectoryFormat::Xyz),
        )
        .unwrap();
        writer.write_frame(buffer.frame_view()).unwrap();
    }
    assert!(!unfinished.exists());

    for file in [&path, &mismatch, &output, &unfinished] {
        let _ = std::fs::remove_file(file);
    }
}

#[test]
fn detection_restores_position_on_read_failure_and_honors_a_zero_byte_limit() {
    let mut reader = FailingDetectionReader {
        cursor: Cursor::new(b"prefix".to_vec()),
    };
    reader.seek(SeekFrom::Start(3)).unwrap();
    let error = detect_trajectory_format(
        &mut reader,
        "failing.xyz",
        TrajectoryFormatHint::Auto,
        &TrajectoryIoLimits::default(),
    )
    .unwrap_err();
    assert!(matches!(
        error,
        error if error.io_error_kind() == Some(io::ErrorKind::ConnectionReset)
    ));
    assert_eq!(reader.stream_position().unwrap(), 3);

    let mut cursor = Cursor::new(TWO_FRAMES.as_bytes());
    cursor.seek(SeekFrom::Start(2)).unwrap();
    let error = detect_trajectory_format(
        &mut cursor,
        "limited.xyz",
        TrajectoryFormatHint::Auto,
        &TrajectoryIoLimits::default().with_max_detection_bytes(0),
    )
    .unwrap_err();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryIoErrorKind::LimitExceeded)
    );
    assert_eq!(cursor.stream_position().unwrap(), 2);
}

#[test]
fn xyz_exact_frame_and_index_limits_still_allow_clean_eof() {
    let topology = topology();
    let limits = TrajectoryIoLimits::default()
        .with_max_frames(2)
        .with_max_index_entries(2)
        .with_max_index_bytes(2 * std::mem::size_of::<u64>());
    let mut reader = XyzReader::new(
        Cursor::new(TWO_FRAMES.as_bytes()),
        Arc::clone(&topology),
        XyzReadOptions::default()
            .with_limits(limits.clone())
            .with_source_label("exact-limit.xyz"),
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    assert!(reader.read_next(&mut buffer).unwrap());
    assert!(reader.read_next(&mut buffer).unwrap());
    assert!(!reader.read_next(&mut buffer).unwrap());

    let indexed = IndexedXyzReader::new(
        Cursor::new(TWO_FRAMES.as_bytes()),
        Arc::clone(&topology),
        XyzReadOptions::default()
            .with_limits(limits)
            .with_source_label("exact-index-limit.xyz"),
    )
    .unwrap();
    assert_eq!(indexed.frame_count(), Some(2));
}

#[test]
fn xyz_frame_byte_limit_charges_exact_source_bytes_including_terminators() {
    let topology = support::linear_carbon_topology(1);
    // Each row is (input, its exact byte length); the frame is accepted at
    // that limit and rejected one byte below it.
    for (input, frame_bytes) in [
        ("1\r\n\r\nC 0 0 0\r\n", 14),
        ("1\n\nC 0 0 0\n", 11),
        ("1\n\nC 0 0 0", 10),
    ] {
        assert_eq!(input.len() as u64, frame_bytes);
        let read = |max_frame_bytes| {
            let mut reader = XyzReader::new(
                Cursor::new(input.as_bytes()),
                Arc::clone(&topology),
                XyzReadOptions::default()
                    .with_limits(
                        TrajectoryIoLimits::default().with_max_frame_bytes(max_frame_bytes),
                    )
                    .with_source_label("frame-bytes.xyz"),
            )
            .unwrap();
            let mut buffer = FrameBuffer::new(Arc::clone(&topology));
            let first = reader.read_next(&mut buffer);
            (first, point_xs(&buffer), reader.read_next(&mut buffer))
        };

        let (first, xs, eof) = read(frame_bytes);
        assert!(first.unwrap(), "{input:?}");
        assert_eq!(xs, [0.0]);
        assert!(!eof.unwrap(), "{input:?}");

        let (first, _, _) = read(frame_bytes - 1);
        assert_eq!(
            first,
            Err(codec_error(
                TrajectoryIoErrorKind::LimitExceeded,
                TrajectoryIoOperation::ReadHeader,
                TrajectoryFormat::Xyz,
                "frame-bytes.xyz",
                Some(0),
                &format!(
                    "XYZ frame exceeds configured {}-byte limit",
                    frame_bytes - 1
                ),
            )
            .with_limit(LimitExceeded::new(
                Resource::RecordBytes,
                usize::try_from(frame_bytes - 1).unwrap(),
                usize::try_from(frame_bytes).unwrap(),
            ))),
            "{input:?}"
        );
    }
}

fn read_first_xyz_frame<R: BufRead>(
    stream: R,
    topology: &Arc<Topology>,
    max_scratch_bytes: usize,
) -> (Result<bool, TrajectoryIoError>, FrameBuffer) {
    let mut reader = XyzReader::new(
        stream,
        Arc::clone(topology),
        XyzReadOptions::default()
            .with_limits(TrajectoryIoLimits::default().with_max_scratch_bytes(max_scratch_bytes))
            .with_source_label("scratch.xyz"),
    )
    .unwrap();
    let mut destination = FrameBuffer::new(Arc::clone(topology));
    let result = reader.read_next(&mut destination);
    if result == Ok(true) {
        assert_eq!(reader.read_next(&mut destination), Ok(false));
    }
    (result, destination)
}

#[test]
fn xyz_line_buffer_shares_the_scratch_budget_with_positions() {
    let topology = support::linear_carbon_topology(1);
    let comment = "x".repeat(4096);
    let input = format!("1\n{comment}\nC 0 0 0\n").into_bytes();
    let positions = std::mem::size_of::<Point3>();
    // The retained line buffer must hold the longest line, here the comment
    // with its terminator; the count and atom lines alone would fit.
    let required = positions + comment.len() + 1;
    let exceeded = codec_error(
        TrajectoryIoErrorKind::LimitExceeded,
        TrajectoryIoOperation::ReadHeader,
        TrajectoryFormat::Xyz,
        "scratch.xyz",
        Some(0),
        "XYZ line buffer: aggregate decode scratch exceeds the configured limit",
    );
    let untouched = support::buffer_snapshot(&FrameBuffer::new(Arc::clone(&topology)));

    // A seven-byte buffer splits the comment across many chunks; the result
    // must not depend on how the stream is buffered.
    let streams = || -> [Box<dyn BufRead>; 2] {
        [
            Box::new(Cursor::new(input.clone())),
            Box::new(BufReader::with_capacity(7, Cursor::new(input.clone()))),
        ]
    };
    for max_scratch_bytes in [positions, required - 1] {
        for stream in streams() {
            let (result, destination) = read_first_xyz_frame(stream, &topology, max_scratch_bytes);
            assert_eq!(result, Err(exceeded.clone()), "{max_scratch_bytes}");
            assert_eq!(support::buffer_snapshot(&destination), untouched);
        }
        let indexed = IndexedXyzReader::new(
            Cursor::new(input.clone()),
            Arc::clone(&topology),
            XyzReadOptions::default()
                .with_limits(
                    TrajectoryIoLimits::default().with_max_scratch_bytes(max_scratch_bytes),
                )
                .with_source_label("scratch.xyz"),
        );
        assert_eq!(indexed.err(), Some(exceeded.clone()), "{max_scratch_bytes}");
    }
    for stream in streams() {
        let (result, destination) = read_first_xyz_frame(stream, &topology, required);
        assert_eq!(result, Ok(true));
        assert_eq!(point_xs(&destination), [0.0]);
    }
}

#[test]
fn xyz_limits_probe_but_do_not_parse_or_consume_frame_n_plus_one() {
    let topology = topology();
    let first = "2\nfirst\nC 0 1 2\nH 3 4 5\n";
    let second = "not-an-atom-count\nthis frame must never be parsed\n";
    let bytes = format!("{first}{second}").into_bytes();
    let second_offset = first.len() as u64;

    let (stream, control) = GuardedCursor::new(bytes.clone(), second_offset);
    let mut reader = XyzReader::new(
        stream,
        Arc::clone(&topology),
        XyzReadOptions::default()
            .with_limits(TrajectoryIoLimits::default().with_max_frames(1))
            .with_source_label("guarded-sequential.xyz"),
    )
    .unwrap();
    let mut destination = FrameBuffer::new(Arc::clone(&topology));
    assert!(reader.read_next(&mut destination).unwrap());
    assert_eq!(
        codec_kind(&reader.read_next(&mut destination).unwrap_err()),
        Some(TrajectoryIoErrorKind::LimitExceeded)
    );
    assert!(!control.violated());
    assert_eq!(control.probed_bytes(), 0);

    for limits in [
        TrajectoryIoLimits::default().with_max_frames(1),
        TrajectoryIoLimits::default().with_max_index_entries(1),
        TrajectoryIoLimits::default().with_max_index_bytes(std::mem::size_of::<u64>()),
    ] {
        let (stream, control) = GuardedCursor::new(bytes.clone(), second_offset);
        let error = IndexedXyzReader::new(
            stream,
            Arc::clone(&topology),
            XyzReadOptions::default()
                .with_limits(limits)
                .with_source_label("guarded-index.xyz"),
        )
        .err()
        .unwrap();
        assert_eq!(
            codec_kind(&error),
            Some(TrajectoryIoErrorKind::LimitExceeded)
        );
        assert!(!control.violated());
        assert_eq!(control.probed_bytes(), 0);
    }
}

#[test]
fn compressed_wrappers_and_insufficient_signatures_are_not_extension_dispatched() {
    let topology = topology();
    for (bytes, extension) in [(&b"\x1f\x8bgarbage"[..], "xyz"), (&b"not xyz"[..], "xyz")] {
        let path = temporary_path(Some(extension));
        std::fs::write(&path, bytes).unwrap();
        let error = open_trajectory(
            &path,
            Arc::clone(&topology),
            TrajectoryOpenOptions::default(),
        )
        .err()
        .unwrap();
        assert!(matches!(
            codec_kind(&error),
            Some(TrajectoryIoErrorKind::UnsupportedVariant | TrajectoryIoErrorKind::UnknownFormat)
        ));
        let _ = std::fs::remove_file(path);
    }
}

#[test]
fn independently_generated_ase_fixture_matches_expected_frames() {
    let topology = water_topology();
    let fixture = include_str!("../fixtures/ase-3.26.0-water.xyz");
    let digest = Sha256::digest(fixture.as_bytes());
    let actual_digest = digest
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect::<String>();
    assert_eq!(
        actual_digest,
        "91c41f2f0b02034c507bf564dc4084be991f4cc21418604fcb4864e0f2e805ea"
    );
    let mut reader = XyzReader::new(
        Cursor::new(fixture.as_bytes()),
        Arc::clone(&topology),
        XyzReadOptions::default()
            .with_limits(TrajectoryIoLimits::default())
            .with_source_label("ase-3.26.0-water.xyz"),
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(topology);
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_xs_close(&buffer, &[0.0, 0.09572, -0.0239987]);
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_xs_close(&buffer, &[0.01, 0.10572, -0.0139987]);
    assert!(!reader.read_next(&mut buffer).unwrap());
}

#[allow(dead_code)]
fn assert_path_is_inside_temp(path: &Path) {
    assert!(path.starts_with(std::env::temp_dir()));
}

#[test]
fn xyz_streams_debug_their_options_and_progress_without_their_stream() {
    let topology = topology();
    let options = XyzReadOptions::default().with_source_label("memory.xyz");
    let mut reader = XyzReader::new(
        Cursor::new(TWO_FRAMES.as_bytes()),
        Arc::clone(&topology),
        options.clone(),
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_eq!(
        format!("{reader:?}"),
        format!("XyzReader {{ atom_count: 2, options: {options:?}, frame_cursor: 1, .. }}")
    );

    let indexed = IndexedXyzReader::new(
        Cursor::new(TWO_FRAMES.as_bytes()),
        Arc::clone(&topology),
        options.clone(),
    )
    .unwrap();
    assert_eq!(
        format!("{indexed:?}"),
        format!(
            "IndexedXyzReader {{ reader: XyzReader {{ atom_count: 2, options: {options:?}, \
             frame_cursor: 0, .. }}, frame_count: 2, .. }}"
        )
    );

    let write_options = XyzWriteOptions::default();
    let writer = XyzWriter::new(
        Vec::new(),
        Arc::clone(&topology),
        write_options.clone(),
        "out.xyz",
    )
    .unwrap();
    assert_eq!(
        format!("{writer:?}"),
        format!(
            "XyzWriter {{ atom_count: 2, options: {write_options:?}, source_label: \"out.xyz\", \
             frame_count: 0, partial: false, .. }}"
        )
    );
}

#[test]
fn interpreted_xyz_topology_streams_every_frame_of_its_file() {
    // `kekule::xyz` interprets one frame into a topology in atom-line order;
    // the streaming reader then decodes all frames against it.
    let text = "3\nfirst\nO 0 0 0\nH 0.957 0 0\nH -0.24 0.927 0\n\
                3\nsecond\nO 0 0 0.1\nH 0.96 0 0.1\nH -0.25 0.93 0.1\n";
    let interpretation = kekule::xyz::parse_str(text, Default::default())
        .unwrap()
        .interpret(
            kekule::xyz::XyzInterpretOptions::default()
                .with_frame(kekule::xyz::XyzFrameSelection::First),
        )
        .unwrap();
    let topology = interpretation.model().shared_topology();
    assert_eq!(topology.bond_count(), 2);
    let mut reader = XyzReader::new(
        Cursor::new(text.as_bytes()),
        topology,
        XyzReadOptions::default(),
    )
    .unwrap();
    let mut frame = reader.frame_buffer();
    let mut z = Vec::new();
    while reader.read_next(&mut frame).unwrap() {
        z.push(
            frame
                .positions()
                .get(0)
                .unwrap()
                .value_in(ANGSTROM)
                .unwrap()
                .z,
        );
    }
    assert_eq!(z.len(), 2);
    assert!((z[1] - 0.1).abs() < 1e-12);
}

/// Writes one all-zero frame of `topology` and returns the text and report.
fn write_one_frame(topology: Arc<Topology>) -> (String, WriteReport) {
    let points = vec![Point3::new(0.0, 0.0, 0.0); topology.atom_count()];
    let mut buffer = FrameBuffer::new(Arc::clone(&topology));
    buffer
        .replace_from_data(FrameData::new(Quantity::new(&points, ANGSTROM)))
        .unwrap();
    let mut writer = XyzWriter::new(
        Vec::new(),
        topology,
        XyzWriteOptions::default()
            .with_decimal_places(1)
            .with_comment("report"),
        "report.xyz",
    )
    .unwrap();
    writer.write_frame(buffer.frame_view()).unwrap();
    let report = writer.report().clone();
    let text = String::from_utf8(writer.finish().unwrap()).unwrap();
    (text, report)
}

#[test]
fn xyz_writer_reports_the_charges_isotopes_and_hydrogens_it_leaves_out() {
    // Acetate labeled with carbon-13.
    let topology = kekule::smiles::to_topology("[13CH3]C(=O)[O-]").unwrap();
    let (text, report) = write_one_frame(Arc::clone(&topology));
    // Before perception the carbonyl carbon and oxygen have unknown inferred
    // counts, which XYZ readers would fix at zero, as well as the methyl's
    // three fixed hydrogens.
    assert_eq!(
        report.omissions(),
        [
            WriteOmission::Bonds { bonds: 3 },
            WriteOmission::ImplicitHydrogens { atoms: 3 },
            WriteOmission::FormalCharges { atoms: 1 },
            WriteOmission::Isotopes { atoms: 1 },
        ]
    );
    // The text is elements and coordinates only, whatever the report says.
    assert_eq!(
        text,
        "4
report
C 0.0 0.0 0.0
C 0.0 0.0 0.0
O 0.0 0.0 0.0
O 0.0 0.0 0.0
"
    );

    // Perceived counts are known: only the methyl carbon has hydrogens.
    let (perceived_text, perceived) = write_one_frame(Arc::new(topology.perceived().unwrap()));
    assert_eq!(
        perceived.omissions(),
        [
            WriteOmission::Bonds { bonds: 3 },
            WriteOmission::ImplicitHydrogens { atoms: 1 },
            WriteOmission::FormalCharges { atoms: 1 },
            WriteOmission::Isotopes { atoms: 1 },
        ]
    );
    assert_eq!(perceived_text, text);
}

const DIPEPTIDE: &str =
    "ATOM      1  N   ALA A   1      -0.677  -1.230  -0.491  1.00  0.00           N
ATOM      2  CA  ALA A   1      -0.001   0.064  -0.491  1.00  0.00           C
ATOM      3  C   ALA A   1       1.499  -0.110  -0.491  1.00  0.00           C
ATOM      4  O   ALA A   1       2.030  -1.227  -0.502  1.00  0.00           O
ATOM      5  CB  ALA A   1      -0.509   0.856   0.727  1.00  0.00           C
ATOM      6  N   GLY A   2       2.250   0.989  -0.479  1.00  0.00           N
ATOM      7  CA  GLY A   2       3.700   0.871  -0.479  1.00  0.00           C
ATOM      8  C   GLY A   2       4.250   2.290  -0.479  1.00  0.00           C
ATOM      9  O   GLY A   2       3.500   3.250  -0.479  1.00  0.00           O
ATOM     10  OXT GLY A   2       5.480   2.400  -0.479  1.00  0.00           O
END
";

#[test]
fn xyz_writer_reports_the_atom_sites_and_hydrogens_of_a_protein() {
    let interpretation = kekule::pdb::parse_str(DIPEPTIDE, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    let topology = interpretation.model().shared_topology();
    let (_, report) = write_one_frame(Arc::clone(&topology));
    // Heavy atoms only: every implicit count is inferred and unperceived.
    assert_eq!(
        report.omissions(),
        [
            WriteOmission::Bonds { bonds: 9 },
            WriteOmission::ImplicitHydrogens { atoms: 10 },
            WriteOmission::AtomSites { atoms: 10 },
        ]
    );
    // Perceived, the hydrogens of N, CA and CB of alanine, N and CA of
    // glycine, and the terminal OXT.
    let (_, perceived) = write_one_frame(Arc::new(topology.perceived().unwrap()));
    assert_eq!(
        perceived.omissions(),
        [
            WriteOmission::Bonds { bonds: 9 },
            WriteOmission::ImplicitHydrogens { atoms: 6 },
            WriteOmission::AtomSites { atoms: 10 },
        ]
    );
}

#[test]
fn an_interpreted_xyz_structure_loses_only_its_bonds() {
    let text = "3
water
O 0 0 0
H 0.957 0 0
H -0.24 0.927 0
";
    let interpretation = kekule::xyz::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    let topology = interpretation.model().shared_topology();
    // Interpretation fixes every count at zero and gives unnamed sites.
    let bonds = [WriteOmission::Bonds { bonds: 2 }];
    let (_, report) = write_one_frame(Arc::clone(&topology));
    assert_eq!(report.omissions(), bonds);

    // Path writers and `write_trajectory` return the codec's report; DCD
    // stores coordinates for a separately supplied topology, so it omits
    // nothing.
    for (extension, expected) in [("xyz", &bonds[..]), ("dcd", &[])] {
        let path = temporary_path(Some(extension));
        let writer =
            create_trajectory_writer(&path, Arc::clone(&topology), Default::default()).unwrap();
        assert_eq!(writer.report().omissions(), expected);
    }
    let trajectory = kekule::structure::Trajectory::from_items(
        Arc::clone(&topology),
        [kekule::structure::TrajectoryFrame::new(
            interpretation.model().positions().clone(),
        )],
    )
    .unwrap();
    let path = temporary_path(Some("xyz"));
    let saved = kekule_traj::io::write_trajectory(&path, &trajectory, Default::default());
    let _ = std::fs::remove_file(&path);
    assert_eq!(saved.unwrap().omissions(), bonds);
}
