use std::io::Cursor;

use molecular::core::{Atom, Element, Molecule};
use molecular::geometry::{PeriodicCell, Point3, Vector3};
use molecular::small::SmallMolecule;
use molecular::topology::{MoleculeInstanceMetadata, Topology, TopologyBuilder};
use molecular::trajectory::{
    AtomOrderAssertion, FrameBuffer, SeekableTrajectoryReader, TrajectoryCodecErrorKind,
    TrajectoryError, TrajectoryReader, TrajectoryWriter,
};
use molecular::units::{Quantity, NANOMETER, PICOSECOND};
use molecular_trajectory_io::xtc::{
    XtcMagic, XtcReadOptions, XtcReader, XtcWriteOptions, XtcWriter,
};
use molecular_trajectory_io::{TrajectoryIoLimits, TrajectoryTopologyBinding};
use sha2::{Digest, Sha256};

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
        .insert("unsupported".into(), molecular::core::PropValue::Bool(true));
    assert_eq!(
        writer.write_frame(frame.frame_view()).unwrap_err(),
        TrajectoryError::UnsupportedField("properties")
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
