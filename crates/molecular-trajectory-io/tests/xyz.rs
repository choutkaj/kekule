use std::io::{self, Cursor, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use molecular::core::{Atom, Element, Molecule, PropValue};
use molecular::geometry::{PeriodicCell, Point3, Vector3};
use molecular::small::SmallMolecule;
use molecular::topology::{MoleculeInstanceMetadata, Topology, TopologyBuilder};
use molecular::trajectory::{
    AtomOrderAssertion, FrameBuffer, FrameBufferData, SeekableTrajectoryReader,
    TrajectoryCodecErrorKind, TrajectoryError, TrajectoryFormat, TrajectoryReader,
    TrajectoryWriter,
};
use molecular::units::{Quantity, ANGSTROM, MODEL_VELOCITY_UNIT, NANOMETER, PICOSECOND};
use molecular_trajectory_io::xyz::{XyzReadOptions, XyzReader, XyzWriteOptions, XyzWriter};
use molecular_trajectory_io::{
    create_trajectory_writer, detect_trajectory_format, open_indexed_trajectory, open_trajectory,
    FieldAvailability, FormatDetectionEvidence, RandomAccessCapability, TrajectoryFormatHint,
    TrajectoryIoLimits, TrajectoryOpenOptions, TrajectoryTopologyBinding, TrajectoryWriteOptions,
};
use sha2::{Digest, Sha256};

const TWO_FRAMES: &str = "2\r\nfirst\r\nC 0.0 1.0 2.0\r\nH 3.0 4.0 5.0\r\n\
2\nsecond\nC 1.0 2.0 3.0\nH 4.0 5.0 6.0";

fn topology() -> Topology {
    let mut graph = Molecule::new();
    graph
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    graph
        .add_atom(Atom::new(Element::from_symbol("H").unwrap()))
        .unwrap();
    let molecule = SmallMolecule::from_graph(graph);
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_small_molecule_definition(&molecule).unwrap();
    builder
        .add_instance(definition, MoleculeInstanceMetadata::default())
        .unwrap();
    builder.build().unwrap()
}

fn water_topology() -> Topology {
    let mut graph = Molecule::new();
    for symbol in ["O", "H", "H"] {
        graph
            .add_atom(Atom::new(Element::from_symbol(symbol).unwrap()))
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

fn point_xs(buffer: &FrameBuffer) -> Vec<f64> {
    buffer
        .configuration()
        .positions()
        .values()
        .value()
        .iter()
        .map(|point| point.x)
        .collect()
}

fn temporary_path(extension: Option<&str>) -> PathBuf {
    let nonce = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let mut name = format!("molecular-trajectory-io-{}-{nonce}", std::process::id());
    if let Some(extension) = extension {
        name.push('.');
        name.push_str(extension);
    }
    std::env::temp_dir().join(name)
}

fn codec_kind(error: &TrajectoryError) -> Option<TrajectoryCodecErrorKind> {
    match error {
        TrajectoryError::Codec(context) => Some(context.kind()),
        _ => None,
    }
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
        binding(&topology),
        XyzReadOptions::default(),
        TrajectoryIoLimits::default(),
        "memory.xyz",
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(topology.clone());
    let pointer = buffer.configuration().positions().values().value().as_ptr();
    buffer.set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(1.0, 1.0, 1.0), NANOMETER),
            [true; 3],
        )
        .unwrap(),
    ));
    buffer
        .set_velocities(Some(Quantity::new(
            [Vector3::new(1.0, 0.0, 0.0), Vector3::new(2.0, 0.0, 0.0)],
            MODEL_VELOCITY_UNIT,
        )))
        .unwrap();
    buffer
        .set_time(Some(Quantity::new(1.0, PICOSECOND)))
        .unwrap();
    buffer.set_step(Some(7));
    buffer
        .props_mut()
        .insert("stale".into(), PropValue::Bool(true));

    assert!(reader.read_next(&mut buffer).unwrap());
    assert_eq!(point_xs(&buffer), vec![0.0, 3.0]);
    assert_eq!(
        buffer.configuration().positions().values().value().as_ptr(),
        pointer
    );
    assert!(buffer.configuration().cell().is_none());
    assert!(buffer.frame_view().velocities().is_none());
    assert!(buffer.frame_view().time().is_none());
    assert!(buffer.frame_view().step().is_none());
    assert!(buffer.props().is_empty());

    assert!(reader.read_next(&mut buffer).unwrap());
    assert_eq!(point_xs(&buffer), vec![1.0, 4.0]);
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
        binding(&topology),
        XyzReadOptions::default().with_length_unit(NANOMETER),
        TrajectoryIoLimits::default(),
        "nanometers.xyz",
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(topology.clone());
    reader.read_next(&mut buffer).unwrap();
    assert_eq!(point_xs(&buffer), vec![1.0, 4.0]);

    for (input, expected) in [
        (
            "2\nbad\nC 0 0 0\nO 1 1 1\n",
            TrajectoryCodecErrorKind::InconsistentMetadata,
        ),
        (
            "2\nbad\nC 0 0 0\n",
            TrajectoryCodecErrorKind::TruncatedRecord,
        ),
        (
            "2\nbad\nC 0 0 0\nH NaN 1 1\n",
            TrajectoryCodecErrorKind::InvalidFrame,
        ),
        (
            "3\nbad\nC 0 0 0\nH 1 1 1\nH 2 2 2\n",
            TrajectoryCodecErrorKind::InconsistentAtomCount,
        ),
    ] {
        let mut reader = XyzReader::new(
            Cursor::new(input.as_bytes()),
            binding(&topology),
            XyzReadOptions::default(),
            TrajectoryIoLimits::default(),
            "bad.xyz",
        )
        .unwrap();
        let mut destination = FrameBuffer::new(topology.clone());
        destination
            .replace_from_data(FrameBufferData::new(
                &topology,
                Quantity::new(
                    &[Point3::new(9.0, 0.0, 0.0), Point3::new(8.0, 0.0, 0.0)],
                    ANGSTROM,
                ),
            ))
            .unwrap();
        let before = point_xs(&destination);
        let error = reader.read_next(&mut destination).unwrap_err();
        assert_eq!(codec_kind(&error), Some(expected));
        assert_eq!(point_xs(&destination), before);
    }

    let limits = TrajectoryIoLimits {
        max_text_line_bytes: 8,
        ..TrajectoryIoLimits::default()
    };
    let mut reader = XyzReader::new(
        Cursor::new("2\ncomment\nC 000000000 0 0\nH 1 1 1\n".as_bytes()),
        binding(&topology),
        XyzReadOptions::default(),
        limits,
        "limited.xyz",
    )
    .unwrap();
    assert_eq!(
        codec_kind(
            &reader
                .read_next(&mut FrameBuffer::new(topology))
                .unwrap_err()
        ),
        Some(TrajectoryCodecErrorKind::ResourceLimitExceeded)
    );
}

#[test]
fn indexed_xyz_matches_sequential_and_random_reads_preserve_cursor() {
    let topology = topology();
    let reader = XyzReader::new(
        Cursor::new(TWO_FRAMES.as_bytes()),
        binding(&topology),
        XyzReadOptions::default(),
        TrajectoryIoLimits::default(),
        "indexed.xyz",
    )
    .unwrap();
    let mut reader = reader.into_indexed().unwrap();
    assert_eq!(reader.frame_count(), Some(2));
    let mut buffer = FrameBuffer::new(topology);

    reader.read_frame(1, &mut buffer).unwrap();
    assert_eq!(point_xs(&buffer), vec![1.0, 4.0]);
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_eq!(point_xs(&buffer), vec![0.0, 3.0]);
    assert!(matches!(
        reader.read_frame(2, &mut buffer),
        Err(TrajectoryError::FrameIndexOutOfRange(2))
    ));
    assert_eq!(point_xs(&buffer), vec![0.0, 3.0]);
}

#[test]
fn xyz_writer_is_strict_and_round_trips_without_owned_frames() {
    let topology = topology();
    let points = [Point3::new(1.25, 2.5, 3.75), Point3::new(4.0, 5.0, 6.0)];
    let mut buffer = FrameBuffer::new(topology.clone());
    buffer
        .replace_from_data(FrameBufferData::new(
            &topology,
            Quantity::new(&points, ANGSTROM),
        ))
        .unwrap();

    let mut writer = XyzWriter::new(
        Vec::new(),
        topology.clone(),
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
        binding(&topology),
        XyzReadOptions::default(),
        TrajectoryIoLimits::default(),
        "round-trip.xyz",
    )
    .unwrap();
    let mut decoded = FrameBuffer::new(topology.clone());
    assert!(reader.read_next(&mut decoded).unwrap());
    assert_eq!(
        decoded.configuration().positions().values().value(),
        &points
    );
    assert!(reader.read_next(&mut decoded).unwrap());
    assert!(!reader.read_next(&mut decoded).unwrap());

    buffer.set_step(Some(1));
    let mut strict = XyzWriter::new(
        Vec::new(),
        topology,
        XyzWriteOptions::default(),
        "strict.xyz",
    )
    .unwrap();
    assert_eq!(
        codec_kind(&strict.write_frame(buffer.frame_view()).unwrap_err()),
        Some(TrajectoryCodecErrorKind::UnsupportedField)
    );
}

#[test]
fn path_detection_metadata_indexing_and_atomic_finish_are_bounded() {
    let topology = topology();
    let path = temporary_path(None);
    std::fs::write(&path, TWO_FRAMES).unwrap();
    let (mut sequential, report) =
        open_trajectory(&path, binding(&topology), TrajectoryOpenOptions::default()).unwrap();
    assert_eq!(report.selected_format(), TrajectoryFormat::Xyz);
    assert!(report
        .detection_evidence()
        .contains(&FormatDetectionEvidence::MissingExtension));
    assert_eq!(sequential.metadata().atom_count(), 2);
    assert_eq!(
        sequential.metadata().fields().positions,
        FieldAvailability::Required
    );
    assert_eq!(
        sequential.metadata().random_access(),
        RandomAccessCapability::SequentialOnly
    );
    assert!(sequential
        .read_next(&mut FrameBuffer::new(topology.clone()))
        .unwrap());

    let (indexed, _) =
        open_indexed_trajectory(&path, binding(&topology), TrajectoryOpenOptions::default())
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
                binding(&topology),
                TrajectoryOpenOptions::default()
            )
            .err()
            .unwrap()
        ),
        Some(TrajectoryCodecErrorKind::FormatMismatch)
    );
    open_trajectory(
        &mismatch,
        binding(&topology),
        TrajectoryOpenOptions::default()
            .with_format_hint(TrajectoryFormatHint::Explicit(TrajectoryFormat::Xyz)),
    )
    .unwrap();

    let output = temporary_path(Some("xyz"));
    let points = [Point3::new(1.0, 2.0, 3.0), Point3::new(4.0, 5.0, 6.0)];
    let mut buffer = FrameBuffer::new(topology.clone());
    buffer
        .replace_from_data(FrameBufferData::new(
            &topology,
            Quantity::new(&points, ANGSTROM),
        ))
        .unwrap();
    let mut writer = create_trajectory_writer(
        &output,
        topology.clone(),
        TrajectoryWriteOptions::new(TrajectoryFormat::Xyz),
    )
    .unwrap();
    writer.write_frame(buffer.frame_view()).unwrap();
    assert!(!output.exists());
    writer.finish().unwrap();
    assert!(output.exists());
    open_trajectory(
        &output,
        binding(&topology),
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
        TrajectoryError::Io(context) if context.error_kind() == io::ErrorKind::ConnectionReset
    ));
    assert_eq!(reader.stream_position().unwrap(), 3);

    let mut cursor = Cursor::new(TWO_FRAMES.as_bytes());
    cursor.seek(SeekFrom::Start(2)).unwrap();
    let error = detect_trajectory_format(
        &mut cursor,
        "limited.xyz",
        TrajectoryFormatHint::Auto,
        &TrajectoryIoLimits {
            max_detection_bytes: 0,
            ..TrajectoryIoLimits::default()
        },
    )
    .unwrap_err();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::ResourceLimitExceeded)
    );
    assert_eq!(cursor.stream_position().unwrap(), 2);
}

#[test]
fn path_writer_failure_poisoning_prevents_partial_publication() {
    let topology = topology();
    let output = temporary_path(Some("xyz"));
    let mut frame = FrameBuffer::new(topology.clone());
    frame
        .set_positions(Quantity::new(
            [Point3::new(1.0, 2.0, 3.0), Point3::new(4.0, 5.0, 6.0)],
            ANGSTROM,
        ))
        .unwrap();
    let mut writer = create_trajectory_writer(
        &output,
        topology,
        TrajectoryWriteOptions::new(TrajectoryFormat::Xyz),
    )
    .unwrap();
    writer.write_frame(frame.frame_view()).unwrap();
    frame.set_step(Some(1));
    assert_eq!(
        codec_kind(&writer.write_frame(frame.frame_view()).unwrap_err()),
        Some(TrajectoryCodecErrorKind::UnsupportedField)
    );
    assert_eq!(
        codec_kind(&writer.finish().unwrap_err()),
        Some(TrajectoryCodecErrorKind::InvalidFrame)
    );
    assert!(!output.exists());
}

#[test]
fn xyz_exact_frame_and_index_limits_still_allow_clean_eof() {
    let topology = topology();
    let limits = TrajectoryIoLimits {
        max_frames: 2,
        max_index_entries: 2,
        max_index_bytes: 2 * std::mem::size_of::<u64>(),
        ..TrajectoryIoLimits::default()
    };
    let mut reader = XyzReader::new(
        Cursor::new(TWO_FRAMES.as_bytes()),
        binding(&topology),
        XyzReadOptions::default(),
        limits.clone(),
        "exact-limit.xyz",
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(topology.clone());
    assert!(reader.read_next(&mut buffer).unwrap());
    assert!(reader.read_next(&mut buffer).unwrap());
    assert!(!reader.read_next(&mut buffer).unwrap());

    let indexed = XyzReader::new(
        Cursor::new(TWO_FRAMES.as_bytes()),
        binding(&topology),
        XyzReadOptions::default(),
        limits,
        "exact-index-limit.xyz",
    )
    .unwrap()
    .into_indexed()
    .unwrap();
    assert_eq!(indexed.frame_count(), Some(2));
}

#[test]
fn compressed_wrappers_and_insufficient_signatures_are_not_extension_dispatched() {
    let topology = topology();
    for (bytes, extension) in [(&b"\x1f\x8bgarbage"[..], "xyz"), (&b"not xyz"[..], "xyz")] {
        let path = temporary_path(Some(extension));
        std::fs::write(&path, bytes).unwrap();
        let error = open_trajectory(&path, binding(&topology), TrajectoryOpenOptions::default())
            .err()
            .unwrap();
        assert!(matches!(
            codec_kind(&error),
            Some(
                TrajectoryCodecErrorKind::UnsupportedVariant
                    | TrajectoryCodecErrorKind::UnknownFormat
            )
        ));
        let _ = std::fs::remove_file(path);
    }
}

#[test]
fn independently_generated_ase_fixture_matches_expected_frames() {
    let topology = water_topology();
    let fixture = include_str!("fixtures/ase-3.26.0-water.xyz");
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
        binding(&topology),
        XyzReadOptions::default(),
        TrajectoryIoLimits::default(),
        "ase-3.26.0-water.xyz",
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(topology);
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_eq!(point_xs(&buffer), vec![0.0, 0.9572, -0.239987]);
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_eq!(point_xs(&buffer), vec![0.1, 1.0572, -0.139987]);
    assert!(!reader.read_next(&mut buffer).unwrap());
}

#[allow(dead_code)]
fn assert_path_is_inside_temp(path: &Path) {
    assert!(path.starts_with(std::env::temp_dir()));
}
