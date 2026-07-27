use std::io::Cursor;

use molecular::core::{Atom, Element, Molecule, PropValue};
use molecular::geometry::{PeriodicCell, Point3, Vector3};
use molecular::small::SmallMolecule;
use molecular::topology::{MoleculeInstanceMetadata, Topology, TopologyBuilder};
use molecular::trajectory::{
    AtomOrderAssertion, FrameBuffer, SeekableTrajectoryReader, TrajectoryCodecErrorKind,
    TrajectoryError, TrajectoryReader, TrajectoryWriter,
};
use molecular::units::{Quantity, MODEL_FORCE_UNIT, MODEL_VELOCITY_UNIT, NANOMETER, PICOSECOND};
use molecular_trajectory_io::trr::{
    TrrLambdaPolicy, TrrReadOptions, TrrReader, TrrScalarPrecision, TrrWriteOptions, TrrWriter,
};
use molecular_trajectory_io::{TrajectoryIoLimits, TrajectoryTopologyBinding};
use sha2::{Digest, Sha256};

fn topology() -> Topology {
    let mut graph = Molecule::new();
    for symbol in ["C", "H", "O"] {
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

fn codec_kind(error: &TrajectoryError) -> Option<TrajectoryCodecErrorKind> {
    match error {
        TrajectoryError::Codec(context) => Some(context.kind()),
        _ => None,
    }
}

fn xs(buffer: &FrameBuffer) -> Vec<f64> {
    buffer
        .configuration()
        .positions()
        .values()
        .value()
        .iter()
        .map(|point| point.x)
        .collect()
}

fn assert_xs_close(buffer: &FrameBuffer, expected: &[f64]) {
    let actual = xs(buffer);
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.into_iter().zip(expected) {
        assert!(
            (actual - expected).abs() <= 1.0e-5,
            "{actual} differs from {expected}"
        );
    }
}

fn populated_frame(topology: &Topology, shift: f64, step: u64) -> FrameBuffer {
    let mut frame = FrameBuffer::new(topology.clone());
    frame
        .set_positions(Quantity::new(
            [
                Point3::new(0.0 + shift, 1.0, 2.0),
                Point3::new(3.0 + shift, 4.0, 5.0),
                Point3::new(6.0 + shift, 7.0, 8.0),
            ],
            NANOMETER,
        ))
        .unwrap();
    frame
        .set_velocities(Some(Quantity::new(
            [
                Vector3::new(1.0, 2.0, 3.0),
                Vector3::new(4.0, 5.0, 6.0),
                Vector3::new(7.0, 8.0, 9.0),
            ],
            MODEL_VELOCITY_UNIT,
        )))
        .unwrap();
    frame
        .set_forces(Some(Quantity::new(
            [
                Vector3::new(10.0, 20.0, 30.0),
                Vector3::new(40.0, 50.0, 60.0),
                Vector3::new(70.0, 80.0, 90.0),
            ],
            MODEL_FORCE_UNIT,
        )))
        .unwrap();
    frame.set_cell(Some(
        PeriodicCell::new(
            Quantity::new(
                [
                    Vector3::new(2.0, 0.0, 0.0),
                    Vector3::new(0.2, 2.1, 0.0),
                    Vector3::new(0.3, 0.4, 2.2),
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
        .props_mut()
        .insert("gromacs.trr.lambda".into(), PropValue::Float(0.125));
    frame
}

#[test]
fn trr_f32_and_f64_round_trip_all_fields_and_clear_absent_state() {
    for precision in [TrrScalarPrecision::Float32, TrrScalarPrecision::Float64] {
        let topology = topology();
        let options = TrrWriteOptions::default().with_precision(precision);
        let mut writer = TrrWriter::new(
            Cursor::new(Vec::new()),
            topology.clone(),
            options,
            "memory.trr",
        )
        .unwrap();
        let first = populated_frame(&topology, 0.0, 4);
        writer.write_frame(first.frame_view()).unwrap();
        let mut second = populated_frame(&topology, 1.0, 5);
        second.set_cell(None);
        second.set_velocities::<&[Vector3]>(None).unwrap();
        second.set_forces::<&[Vector3]>(None).unwrap();
        second
            .props_mut()
            .insert("gromacs.trr.lambda".into(), PropValue::Float(0.25));
        writer.write_frame(second.frame_view()).unwrap();
        let bytes = writer.finish().unwrap().into_inner();

        let mut reader = TrrReader::new(
            Cursor::new(bytes.clone()),
            binding(&topology),
            TrrReadOptions::default(),
            TrajectoryIoLimits::default(),
            "memory.trr",
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
        assert_xs_close(&destination, &[0.0, 30.0, 60.0]);
        assert_eq!(destination.frame_view().step(), Some(4));
        assert_eq!(destination.frame_view().time().unwrap().value(), &1.0);
        assert!(destination.configuration().cell().is_some());
        assert!(destination.frame_view().velocities().is_some());
        assert!(destination.frame_view().forces().is_some());
        assert_eq!(
            destination.props().get("gromacs.trr.lambda"),
            Some(&PropValue::Float(0.125))
        );
        assert!(reader.read_next(&mut destination).unwrap());
        assert_xs_close(&destination, &[10.0, 40.0, 70.0]);
        assert!(destination.configuration().cell().is_none());
        assert!(destination.frame_view().velocities().is_none());
        assert!(destination.frame_view().forces().is_none());
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

        let mut indexed = TrrReader::new(
            Cursor::new(bytes),
            binding(&topology),
            TrrReadOptions::default(),
            TrajectoryIoLimits::default(),
            "memory.trr",
        )
        .unwrap()
        .into_indexed()
        .unwrap();
        assert_eq!(indexed.frame_count(), Some(2));
        indexed.read_frame(1, &mut destination).unwrap();
        assert_xs_close(&destination, &[10.0, 40.0, 70.0]);
        assert!(indexed.read_next(&mut destination).unwrap());
        assert_xs_close(&destination, &[0.0, 30.0, 60.0]);
    }
}

#[test]
fn trr_lambda_policy_and_writer_contract_are_explicit() {
    let topology = topology();
    let mut frame = populated_frame(&topology, 0.0, 0);
    let mut writer = TrrWriter::new(
        Cursor::new(Vec::new()),
        topology.clone(),
        TrrWriteOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero),
        "zero.trr",
    )
    .unwrap();
    assert_eq!(
        writer.write_frame(frame.frame_view()).unwrap_err(),
        TrajectoryError::UnsupportedField("properties")
    );
    frame.props_mut().clear();
    writer.write_frame(frame.frame_view()).unwrap();
    let bytes = writer.finish().unwrap().into_inner();
    let mut reader = TrrReader::new(
        Cursor::new(bytes),
        binding(&topology),
        TrrReadOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero),
        TrajectoryIoLimits::default(),
        "zero.trr",
    )
    .unwrap();
    let mut destination = FrameBuffer::new(topology);
    reader.read_next(&mut destination).unwrap();
    assert!(destination.props().is_empty());
}

#[test]
fn trr_malformed_sizes_truncation_limits_and_eof_are_transactional() {
    let topology = topology();
    let mut writer = TrrWriter::new(
        Cursor::new(Vec::new()),
        topology.clone(),
        TrrWriteOptions::default(),
        "memory.trr",
    )
    .unwrap();
    let first = populated_frame(&topology, 0.0, 0);
    let second = populated_frame(&topology, 1.0, 1);
    writer.write_frame(first.frame_view()).unwrap();
    writer.write_frame(second.frame_view()).unwrap();
    let valid = writer.finish().unwrap().into_inner();

    let mut invalid_size = valid.clone();
    invalid_size[52..56].copy_from_slice(&5_i32.to_be_bytes());
    let error = TrrReader::new(
        Cursor::new(invalid_size),
        binding(&topology),
        TrrReadOptions::default(),
        TrajectoryIoLimits::default(),
        "size.trr",
    )
    .err()
    .unwrap();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::InvalidRecordLength)
    );

    let mut truncated = valid.clone();
    truncated.pop();
    let mut reader = TrrReader::new(
        Cursor::new(truncated),
        binding(&topology),
        TrrReadOptions::default(),
        TrajectoryIoLimits::default(),
        "truncated.trr",
    )
    .unwrap();
    let mut destination = FrameBuffer::new(topology.clone());
    assert!(reader.read_next(&mut destination).unwrap());
    let before = xs(&destination);
    let error = reader.read_next(&mut destination).unwrap_err();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::TruncatedRecord)
    );
    assert_eq!(xs(&destination), before);

    let limits = TrajectoryIoLimits {
        max_frame_bytes: 128,
        ..TrajectoryIoLimits::default()
    };
    let mut reader = TrrReader::new(
        Cursor::new(valid),
        binding(&topology),
        TrrReadOptions::default(),
        limits,
        "limited.trr",
    )
    .unwrap();
    let error = reader.read_next(&mut destination).unwrap_err();
    assert_eq!(
        codec_kind(&error),
        Some(TrajectoryCodecErrorKind::ResourceLimitExceeded)
    );
}

#[test]
fn indexed_trr_accepts_per_frame_precision_and_verifies_both_payloads() {
    let topology = topology();
    let mut combined = Vec::new();
    for (precision, shift, step) in [
        (TrrScalarPrecision::Float32, 0.0, 0),
        (TrrScalarPrecision::Float64, 1.0, 1),
    ] {
        let mut writer = TrrWriter::new(
            Cursor::new(Vec::new()),
            topology.clone(),
            TrrWriteOptions::default().with_precision(precision),
            "mixed.trr",
        )
        .unwrap();
        let frame = populated_frame(&topology, shift, step);
        writer.write_frame(frame.frame_view()).unwrap();
        combined.extend(writer.finish().unwrap().into_inner());
    }
    let mut reader = TrrReader::new(
        Cursor::new(combined),
        binding(&topology),
        TrrReadOptions::default(),
        TrajectoryIoLimits::default(),
        "mixed.trr",
    )
    .unwrap()
    .into_indexed()
    .unwrap();
    let mut destination = FrameBuffer::new(topology);
    assert_eq!(reader.frame_count(), Some(2));
    reader.read_frame(1, &mut destination).unwrap();
    assert_xs_close(&destination, &[10.0, 40.0, 70.0]);
}

#[test]
fn independently_generated_mdanalysis_trr_preserves_all_supported_fields() {
    let topology = topology();
    let fixture = include_bytes!("fixtures/mdanalysis-2.9.0-three-atoms.trr");
    let digest = Sha256::digest(fixture);
    let actual_digest = digest
        .iter()
        .map(|byte| format!("{byte:02x}"))
        .collect::<String>();
    assert_eq!(
        actual_digest,
        "ff960f93b192e71f9b8962a7f8fe70c0aae23fd842e15aa65f3d3fd26a6b07c0"
    );
    let mut reader = TrrReader::new(
        Cursor::new(fixture),
        binding(&topology),
        TrrReadOptions::default(),
        TrajectoryIoLimits::default(),
        "mdanalysis-2.9.0-three-atoms.trr",
    )
    .unwrap();
    let mut buffer = FrameBuffer::new(topology);
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_xs_close(&buffer, &[0.0, 3.0, 6.0]);
    assert!(buffer.configuration().cell().is_some());
    assert!(buffer.frame_view().velocities().is_some());
    assert!(buffer.frame_view().forces().is_some());
    assert_eq!(buffer.frame_view().step(), Some(0));
    assert_eq!(
        buffer.props().get("gromacs.trr.lambda"),
        Some(&PropValue::Float(0.125))
    );
    assert!(reader.read_next(&mut buffer).unwrap());
    assert_xs_close(&buffer, &[1.0, 4.0, 7.0]);
    assert_eq!(buffer.frame_view().step(), Some(1));
    assert_eq!(
        buffer.props().get("gromacs.trr.lambda"),
        Some(&PropValue::Float(0.25))
    );
    assert!(!reader.read_next(&mut buffer).unwrap());
}
