//! Frame state a codec cannot carry, omitted and reported under the
//! `frame_metadata` policies, and frames without a time or step numbered for
//! XTC and TRR under `FrameTimingPolicy::Number`.

use std::io::Cursor;
use std::sync::{Arc, OnceLock};

use kekule::diagnostics::{StructureField, UnrepresentablePolicy, WriteOmission, WriteReport};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{Forces, Positions, Trajectory, TrajectoryFrame, Velocities};
use kekule::topology::Topology;
use kekule::units::{
    Quantity, ANGSTROM, CANONICAL_FORCE_UNIT, CANONICAL_VELOCITY_UNIT, DIMENSIONLESS, NANOMETER,
    PICOSECOND,
};
use kekule_traj::io::dcd::{DcdWriteOptions, DcdWriter};
use kekule_traj::io::trr::{
    TrrLambdaPolicy, TrrReadOptions, TrrReader, TrrWriteOptions, TrrWriter, TRR_LAMBDA_PROPERTY,
};
use kekule_traj::io::xtc::{XtcCellPolicy, XtcReadOptions, XtcReader, XtcWriteOptions, XtcWriter};
use kekule_traj::io::xyz::{XyzWriteOptions, XyzWriter};
use kekule_traj::io::FrameTimingPolicy;
use kekule_traj::{
    TrajectoryField, TrajectoryFormat, TrajectoryIoErrorKind, TrajectoryIoOperation,
    TrajectoryWriter,
};

use crate::support;
use support::{codec_error, read_all};

const OMIT: UnrepresentablePolicy = UnrepresentablePolicy::Omit;

/// One topology for every writer and frame, whose layouts must agree.
fn topology() -> Arc<Topology> {
    static TOPOLOGY: OnceLock<Arc<Topology>> = OnceLock::new();
    Arc::clone(TOPOLOGY.get_or_init(|| support::topology(&["C", "H", "O"], &[(0, 1), (0, 2)])))
}

/// Frame `index` with positions only.
fn bare_frame(index: usize) -> TrajectoryFrame {
    let shift = index as f64;
    TrajectoryFrame::new(
        Positions::new(Quantity::new(
            [
                Point3::new(shift, 1.0, 2.0),
                Point3::new(3.0, 4.0 + shift, 5.0),
                Point3::new(6.0, 7.0, 8.0 + shift),
            ],
            ANGSTROM,
        ))
        .unwrap(),
    )
}

fn set_velocities(frame: &mut TrajectoryFrame) {
    frame
        .set_velocities(Some(
            Velocities::new(Quantity::new(
                [Vector3::new(1.0, 2.0, 3.0); 3],
                CANONICAL_VELOCITY_UNIT,
            ))
            .unwrap(),
        ))
        .unwrap();
}

fn set_forces(frame: &mut TrajectoryFrame) {
    frame
        .set_forces(Some(
            Forces::new(Quantity::new(
                [Vector3::new(4.0, 5.0, 6.0); 3],
                CANONICAL_FORCE_UNIT,
            ))
            .unwrap(),
        ))
        .unwrap();
}

/// Gives every frame atom, bond, and frame properties.
fn with_properties(mut trajectory: Trajectory) -> Trajectory {
    for index in 0..trajectory.len() {
        let mut frame = trajectory.get_mut(index).unwrap();
        let mut properties = frame.properties_mut();
        properties
            .owner_mut()
            .insert(PropertyKey::new("frame").unwrap(), PropertyValue::Int(7))
            .unwrap();
        properties
            .atoms_mut()
            .insert(
                PropertyKey::new("atom").unwrap(),
                PropertyColumn::Int(vec![Some(1), Some(2), Some(3)]),
            )
            .unwrap();
        properties
            .bonds_mut()
            .insert(
                PropertyKey::new("bond").unwrap(),
                PropertyColumn::Int(vec![Some(4), None]),
            )
            .unwrap();
    }
    trajectory
}

fn trajectory(frames: impl IntoIterator<Item = TrajectoryFrame>) -> Trajectory {
    Trajectory::from_items(topology(), frames).unwrap()
}

fn time(picoseconds: f64) -> Option<Quantity<f64>> {
    Some(Quantity::new(picoseconds, PICOSECOND))
}

fn xtc_cell() -> PeriodicCell {
    PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(3.0, 3.0, 3.0), NANOMETER),
        [true; 3],
    )
    .unwrap()
}

fn omitted(field: StructureField, structures: usize) -> WriteOmission {
    WriteOmission::StructureState { field, structures }
}

/// Writes every frame of `trajectory` with `writer` and returns its stream
/// and report.
fn write_all<W: TrajectoryWriter>(
    mut writer: W,
    trajectory: &Trajectory,
    finish: impl FnOnce(W) -> (Vec<u8>, WriteReport),
) -> (Vec<u8>, WriteReport) {
    for frame in trajectory {
        writer.write_frame(frame).unwrap();
    }
    finish(writer)
}

fn xyz(options: XyzWriteOptions, trajectory: &Trajectory) -> (Vec<u8>, WriteReport) {
    let writer = XyzWriter::new(Vec::new(), topology(), options, "state.xyz").unwrap();
    write_all(writer, trajectory, |writer| {
        let report = writer.report().clone();
        (writer.finish().unwrap(), report)
    })
}

fn dcd(options: DcdWriteOptions, trajectory: &Trajectory) -> (Vec<u8>, WriteReport) {
    let writer = DcdWriter::new(Cursor::new(Vec::new()), topology(), options, "state.dcd").unwrap();
    write_all(writer, trajectory, |writer| {
        let report = writer.report().clone();
        (writer.finish().unwrap().into_inner(), report)
    })
}

fn trr(options: TrrWriteOptions, trajectory: &Trajectory) -> (Vec<u8>, WriteReport) {
    let writer = TrrWriter::new(Vec::new(), topology(), options, "state.trr").unwrap();
    write_all(writer, trajectory, |writer| {
        let report = writer.report().clone();
        (writer.finish().unwrap(), report)
    })
}

fn xtc(options: XtcWriteOptions, trajectory: &Trajectory) -> (Vec<u8>, WriteReport) {
    let writer = XtcWriter::new(Vec::new(), topology(), options, "state.xtc").unwrap();
    write_all(writer, trajectory, |writer| {
        let report = writer.report().clone();
        (writer.finish().unwrap(), report)
    })
}

/// Two frames as `base` makes them, and the same frames with `extra` state
/// and properties.
fn plain_and_rich(
    base: impl Fn(usize) -> TrajectoryFrame,
    extra: impl Fn(&mut TrajectoryFrame),
) -> (Trajectory, Trajectory) {
    let plain = trajectory((0..2).map(&base));
    let rich = with_properties(trajectory((0..2).map(|index| {
        let mut frame = base(index);
        extra(&mut frame);
        frame
    })));
    (plain, rich)
}

#[test]
fn xyz_omits_and_reports_every_frame_state_it_cannot_carry() {
    let (plain, rich) = plain_and_rich(bare_frame, |frame| {
        set_velocities(frame);
        set_forces(frame);
        frame.set_time(time(1.0)).unwrap();
        frame.set_step(Some(10));
    });
    let (expected, topology_report) = xyz(XyzWriteOptions::default(), &plain);
    let (bytes, report) = xyz(
        XyzWriteOptions::default()
            .with_frame_metadata(OMIT)
            .with_frame_properties(OMIT),
        &rich,
    );
    assert_eq!(bytes, expected);
    let state = [
        StructureField::Velocities,
        StructureField::Forces,
        StructureField::Time,
        StructureField::Step,
        StructureField::AtomProperties,
        StructureField::BondProperties,
        StructureField::Properties,
    ]
    .map(|field| omitted(field, 2));
    let expected_report: WriteReport = topology_report
        .omissions()
        .iter()
        .copied()
        .chain(state)
        .collect();
    assert_eq!(report, expected_report);
    assert_ne!(report, topology_report);

    let writer = XyzWriter::new(
        Vec::new(),
        topology(),
        XyzWriteOptions::default(),
        "state.xyz",
    );
    assert_eq!(
        writer.unwrap().write_frame(rich.get(0).unwrap()),
        Err(codec_error(
            TrajectoryIoErrorKind::UnsupportedField,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Xyz,
            "state.xyz",
            Some(0),
            "XYZ cannot preserve velocities",
        )
        .with_field(TrajectoryField::Velocities))
    );
}

#[test]
fn dcd_omits_and_reports_times_under_the_absent_policy_and_frame_state() {
    let stepped = |index: usize| {
        let mut frame = bare_frame(index);
        frame.set_step(Some(index as u64));
        frame
    };
    let (plain, rich) = plain_and_rich(stepped, |frame| {
        set_velocities(frame);
        set_forces(frame);
        frame.set_time(time(0.5)).unwrap();
    });
    let (expected, plain_report) = dcd(DcdWriteOptions::default(), &plain);
    assert!(plain_report.is_empty());
    let (bytes, report) = dcd(
        DcdWriteOptions::default()
            .with_frame_metadata(OMIT)
            .with_frame_properties(OMIT),
        &rich,
    );
    assert_eq!(bytes, expected);
    let expected_report: WriteReport = [
        StructureField::Velocities,
        StructureField::Forces,
        StructureField::AtomProperties,
        StructureField::BondProperties,
        StructureField::Properties,
        StructureField::Time,
    ]
    .map(|field| omitted(field, 2))
    .into_iter()
    .collect();
    assert_eq!(report, expected_report);

    // A frame time alone, as frames read from XTC have.
    let timed = trajectory((0..2).map(|index| {
        let mut frame = stepped(index);
        frame.set_time(time(0.5)).unwrap();
        frame
    }));
    let mut writer = DcdWriter::new(
        Cursor::new(Vec::new()),
        topology(),
        DcdWriteOptions::default(),
        "state.dcd",
    )
    .unwrap();
    assert_eq!(
        writer.write_frame(timed.get(0).unwrap()),
        Err(codec_error(
            TrajectoryIoErrorKind::UnsupportedField,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Dcd,
            "state.dcd",
            Some(0),
            "DCD cannot preserve time",
        )
        .with_field(TrajectoryField::Time))
    );
    let (bytes, report) = dcd(DcdWriteOptions::default().with_frame_metadata(OMIT), &timed);
    assert_eq!(bytes, expected);
    assert_eq!(report.omissions(), [omitted(StructureField::Time, 2)]);
}

#[test]
fn xtc_omits_and_reports_velocities_forces_and_properties() {
    let timed = |index: usize| {
        let mut frame = bare_frame(index);
        frame.set_step(Some(index as u64));
        frame.set_time(time(index as f64)).unwrap();
        frame.set_cell(Some(xtc_cell()));
        frame
    };
    let (plain, rich) = plain_and_rich(timed, |frame| {
        set_velocities(frame);
        set_forces(frame);
    });
    let (expected, plain_report) = xtc(XtcWriteOptions::default(), &plain);
    assert!(plain_report.is_empty());
    let (bytes, report) = xtc(
        XtcWriteOptions::default()
            .with_frame_metadata(OMIT)
            .with_frame_properties(OMIT),
        &rich,
    );
    assert_eq!(bytes, expected);
    let expected_report: WriteReport = [
        StructureField::Velocities,
        StructureField::Forces,
        StructureField::AtomProperties,
        StructureField::BondProperties,
        StructureField::Properties,
    ]
    .map(|field| omitted(field, 2))
    .into_iter()
    .collect();
    assert_eq!(report, expected_report);

    // Omitting one kind of frame state never omits the other.
    for (options, field) in [
        (
            XtcWriteOptions::default().with_frame_metadata(OMIT),
            TrajectoryField::AtomProperties,
        ),
        (
            XtcWriteOptions::default().with_frame_properties(OMIT),
            TrajectoryField::Velocities,
        ),
    ] {
        let mut writer = XtcWriter::new(Vec::new(), topology(), options, "state.xtc").unwrap();
        assert_eq!(
            writer.write_frame(rich.get(0).unwrap()),
            Err(codec_error(
                TrajectoryIoErrorKind::UnsupportedField,
                TrajectoryIoOperation::WriteFrame,
                TrajectoryFormat::Xtc,
                "state.xtc",
                Some(0),
                &format!("XTC cannot preserve {field}"),
            )
            .with_field(field))
        );
    }
}

#[test]
fn trr_omits_properties_but_keeps_its_lambda() {
    let timed = |index: usize| {
        let mut frame = bare_frame(index);
        frame.set_step(Some(index as u64));
        frame.set_time(time(index as f64)).unwrap();
        frame
    };
    let lambda = |frame: &mut TrajectoryFrame| {
        frame
            .properties_mut()
            .owner_mut()
            .insert(
                PropertyKey::new(TRR_LAMBDA_PROPERTY).unwrap(),
                PropertyValue::Real {
                    value: 0.25,
                    unit: DIMENSIONLESS,
                },
            )
            .unwrap();
    };
    let with_lambda = |index: usize| {
        let mut frame = timed(index);
        lambda(&mut frame);
        frame
    };
    let (plain, rich) = plain_and_rich(with_lambda, |_| {});

    let (expected, _) = trr(TrrWriteOptions::default(), &plain);
    let (bytes, report) = trr(
        TrrWriteOptions::default().with_frame_properties(OMIT),
        &rich,
    );
    assert_eq!(bytes, expected);
    let expected_report: WriteReport = [
        StructureField::AtomProperties,
        StructureField::BondProperties,
        StructureField::Properties,
    ]
    .map(|field| omitted(field, 2))
    .into_iter()
    .collect();
    assert_eq!(report, expected_report);
    let read = read_all(
        &mut TrrReader::new(Cursor::new(bytes), topology(), TrrReadOptions::default()).unwrap(),
    );
    assert_eq!(read[1].properties(), plain.get(1).unwrap().properties());

    // The lambda policy's own requirements are never omitted, and
    // RequireZero leaves out every frame property.
    let (_, without_lambda) = plain_and_rich(timed, |_| {});
    let mut writer = TrrWriter::new(
        Vec::new(),
        topology(),
        TrrWriteOptions::default().with_frame_properties(OMIT),
        "state.trr",
    )
    .unwrap();
    assert_eq!(
        writer.write_frame(without_lambda.get(0).unwrap()),
        Err(codec_error(
            TrajectoryIoErrorKind::InconsistentMetadata,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Trr,
            "state.trr",
            None,
            "TRR FrameProperty policy requires exactly gromacs.trr.lambda",
        ))
    );
    let zero = TrrWriteOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero);
    let (expected, _) = trr(zero, &trajectory((0..2).map(timed)));
    let (bytes, report) = trr(zero.with_frame_properties(OMIT), &rich);
    assert_eq!(bytes, expected);
    assert_eq!(report, expected_report);
}

#[test]
fn xtc_and_trr_number_frames_without_a_time_or_step() {
    let interval = Quantity::new(2.0, PICOSECOND);
    let numbered = FrameTimingPolicy::Number { interval };
    // Frames as read from XYZ, and one with a step of its own as read from
    // DCD, which keeps it.
    let unnumbered = trajectory((0..3).map(|index| {
        let mut frame = bare_frame(index);
        if index == 2 {
            frame.set_step(Some(500));
        }
        frame
    }));
    let explicit = trajectory((0..3).map(|index| {
        let mut frame = bare_frame(index);
        frame.set_step(Some(if index == 2 { 500 } else { index as u64 }));
        frame.set_time(time(2.0 * index as f64)).unwrap();
        frame
    }));

    let xtc_options =
        XtcWriteOptions::default().with_cell_policy(XtcCellPolicy::ZeroMatrixAsAbsent);
    let (expected, _) = xtc(xtc_options, &explicit);
    let (bytes, report) = xtc(xtc_options.with_frame_timing(numbered), &unnumbered);
    assert_eq!(bytes, expected);
    assert!(report.is_empty());
    let read = read_all(
        &mut XtcReader::new(
            Cursor::new(bytes),
            topology(),
            XtcReadOptions::default().with_cell_policy(XtcCellPolicy::ZeroMatrixAsAbsent),
        )
        .unwrap(),
    );
    let timing = read
        .iter()
        .map(|frame| {
            (
                frame.step(),
                frame.time().unwrap().value_in(PICOSECOND).unwrap(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(timing, [(Some(0), 0.0), (Some(1), 2.0), (Some(500), 4.0)]);

    let trr_options = TrrWriteOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero);
    let (expected, _) = trr(trr_options, &explicit);
    let (bytes, _) = trr(trr_options.with_frame_timing(numbered), &unnumbered);
    assert_eq!(bytes, expected);

    // By default such frames still fail.
    let mut writer = XtcWriter::new(Vec::new(), topology(), xtc_options, "state.xtc").unwrap();
    assert_eq!(
        writer.write_frame(unnumbered.get(0).unwrap()),
        Err(codec_error(
            TrajectoryIoErrorKind::InconsistentMetadata,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Xtc,
            "state.xtc",
            None,
            "XTC requires explicit frame step; FrameTimingPolicy::Number numbers frames without one",
        ))
    );
}

#[test]
fn frame_numbering_needs_a_finite_nonnegative_time_interval() {
    for interval in [
        Quantity::new(1.0, NANOMETER),
        Quantity::new(-1.0, PICOSECOND),
        Quantity::new(f64::NAN, PICOSECOND),
        Quantity::new(f64::INFINITY, PICOSECOND),
    ] {
        let timing = FrameTimingPolicy::Number { interval };
        let detail = |format| {
            codec_error(
                TrajectoryIoErrorKind::InconsistentMetadata,
                TrajectoryIoOperation::Open,
                format,
                "state",
                None,
                &format!("{format} frame numbering interval must be a finite, nonnegative time"),
            )
        };
        let xtc = XtcWriteOptions::default().with_frame_timing(timing);
        assert_eq!(
            XtcWriter::new(Vec::new(), topology(), xtc, "state").map(|_| ()),
            Err(detail(TrajectoryFormat::Xtc)),
            "{interval:?}"
        );
        let trr = TrrWriteOptions::default().with_frame_timing(timing);
        assert_eq!(
            TrrWriter::new(Vec::new(), topology(), trr, "state").map(|_| ()),
            Err(detail(TrajectoryFormat::Trr)),
            "{interval:?}"
        );
    }
    let zero = FrameTimingPolicy::Number {
        interval: Quantity::new(0.0, PICOSECOND),
    };
    let xtc = XtcWriteOptions::default().with_frame_timing(zero);
    assert!(XtcWriter::new(Vec::new(), topology(), xtc, "state").is_ok());
}
