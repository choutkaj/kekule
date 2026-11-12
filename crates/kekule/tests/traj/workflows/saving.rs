use std::{
    fs, io,
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
    sync::Arc,
};

use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{Forces, Positions, SampleSet, SampleView, Velocities};
use kekule::traj::io::dcd::{DcdCellPolicy, DcdWriteOptions};
use kekule::traj::io::trr::{TrrScalarPrecision, TrrWriteOptions};
use kekule::traj::io::{
    create_trajectory_writer, read_sample_set, write_sample_set, OverwritePolicy,
    TrajectoryFormatHint, TrajectoryOpenOptions, TrajectoryWriteOptions,
};
use kekule::traj::{
    FrameView, TrajectoryField, TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind,
    TrajectoryIoOperation, TrajectoryWriter,
};

use crate::support::{
    codec_error, contents, dynamics, fingerprint_bytes, frames_of, linear_carbon_topology,
    time_ordered, topology, Frame, Recorded,
};

struct Directory(PathBuf);
impl Directory {
    fn new() -> Self {
        static SEQUENCE: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "kekule-save-{}-{}",
            std::process::id(),
            SEQUENCE.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
    fn file(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
    /// Every entry, sorted, so a stray temporary file fails by name.
    fn names(&self) -> Vec<String> {
        let mut names = fs::read_dir(&self.0)
            .unwrap()
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .collect::<Vec<_>>();
        names.sort();
        names
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}
/// The source label a path writer reports for `path`.
fn label(path: &Path) -> String {
    path.display().to_string()
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/traj")
        .join(name)
}

/// Sample `index` of `set` as the frame a writer takes.
fn frame(set: &SampleSet, index: usize) -> FrameView<'_> {
    FrameView::from_sample(set.get(index).unwrap()).unwrap()
}

/// The MD steps of a set's samples.
fn steps(set: &SampleSet) -> Vec<Option<u64>> {
    frames_of(set).iter().map(Frame::step).collect()
}

fn same_frame(actual: SampleView<'_>, expected: SampleView<'_>, tolerance: f64) {
    // Codec precision and unit conversions permit numerical rounding. Field
    // presence, dimensions, units, steps, and every property remain exact.
    let actual = FrameView::from_sample(actual).unwrap();
    let expected = FrameView::from_sample(expected).unwrap();
    let close =
        |a: f64, b: f64| assert!((a - b).abs() <= tolerance * b.abs().max(1.0), "{a} != {b}");
    assert_eq!(actual.positions().len(), expected.positions().len());
    for (a, b) in actual
        .positions()
        .values()
        .value()
        .iter()
        .zip(expected.positions().values().value().iter())
    {
        close(a.x, b.x);
        close(a.y, b.y);
        close(a.z, b.z);
    }
    assert_eq!(actual.cell().is_some(), expected.cell().is_some());
    if let (Some(a), Some(b)) = (actual.cell(), expected.cell()) {
        assert_eq!(a.periodic_axes(), b.periodic_axes());
        for (a, b) in a
            .vectors()
            .into_value()
            .into_iter()
            .zip(b.vectors().into_value())
        {
            close(a.x, b.x);
            close(a.y, b.y);
            close(a.z, b.z);
        }
    }
    for (actual, expected) in [
        (
            actual.velocities().map(Velocities::values),
            expected.velocities().map(Velocities::values),
        ),
        (
            actual.forces().map(Forces::values),
            expected.forces().map(Forces::values),
        ),
    ] {
        assert_eq!(actual.is_some(), expected.is_some());
        if let (Some(a), Some(b)) = (actual, expected) {
            assert_eq!(a.unit(), b.unit());
            assert_eq!(a.value().len(), b.value().len());
            for (a, b) in a.value().iter().zip(b.value().iter()) {
                close(a.x, b.x);
                close(a.y, b.y);
                close(a.z, b.z);
            }
        }
    }
    assert_eq!(actual.time().is_some(), expected.time().is_some());
    if let (Some(a), Some(b)) = (actual.time(), expected.time()) {
        assert_eq!(a.unit(), b.unit());
        close(a.into_value(), b.into_value());
    }
    assert_eq!(actual.step(), expected.step());
    assert_eq!(actual.properties(), expected.properties());
}

#[test]
fn saving_round_trips_all_codecs_using_defaults_or_explicit_native_options() {
    let directory = Directory::new();
    for (source, extension, topology, options) in [
        (
            "ase-3.26.0-water.xyz",
            "XYZ",
            topology(&["O", "H", "H"], &[(0, 1), (0, 2)]),
            None,
        ),
        (
            "mdanalysis-2.9.0-three-atoms.trr",
            "trr",
            topology(&["C", "H", "O"], &[(0, 1), (0, 2)]),
            None,
        ),
        (
            "mdanalysis-2.9.0-twelve-atoms.xtc",
            "xtc",
            linear_carbon_topology(12),
            None,
        ),
        (
            "mdanalysis-2.9.0-three-atoms.dcd",
            "dcd",
            topology(&["C", "H", "O"], &[(0, 1), (0, 2)]),
            Some(
                TrajectoryWriteOptions::new(TrajectoryFormat::Dcd).with_dcd_options(
                    DcdWriteOptions::default().with_cell_policy(DcdCellPolicy::Write),
                ),
            ),
        ),
    ] {
        let trajectory =
            read_sample_set(fixture(source), topology.clone(), Default::default()).unwrap();
        let path = directory.file(&format!("output.{extension}"));
        match options {
            Some(options) => write_sample_set(&path, &trajectory, options).unwrap(),
            None => write_sample_set(&path, &trajectory, Default::default()).unwrap(),
        };
        let read = read_sample_set(&path, topology, Default::default()).unwrap();
        assert_eq!(read.len(), trajectory.len());
        for (actual, expected) in read.iter().zip(trajectory.iter()) {
            same_frame(actual, expected, 1.0e-6);
        }
    }
    assert_eq!(
        directory.names(),
        ["output.XYZ", "output.dcd", "output.trr", "output.xtc"]
    );
}

#[test]
fn explicit_format_overrides_extension_and_precision_is_preserved() {
    let directory = Directory::new();
    let topology = topology(&["C", "H", "O"], &[(0, 1), (0, 2)]);
    let trajectory = read_sample_set(
        fixture("mdanalysis-2.9.0-three-atoms.trr"),
        topology.clone(),
        Default::default(),
    )
    .unwrap();
    let path = directory.file("trajectory.data");
    let unknown: TrajectoryIoError = TrajectoryIoError::codec(
        TrajectoryIoErrorKind::UnknownFormat,
        TrajectoryIoOperation::Open,
        None,
    )
    .with_source_label(label(&path))
    .with_detail(
        "The destination's extension names no supported format. Choose the format explicitly.",
    );
    assert_eq!(
        write_sample_set(&path, &trajectory, Default::default()),
        Err(unknown.clone())
    );
    // Default options take the path writer's format from the extension too.
    assert_eq!(
        create_trajectory_writer(&path, topology.clone(), TrajectoryWriteOptions::default())
            .unwrap_err(),
        unknown
    );
    let saved = directory.file("saved.trr");
    write_sample_set(&saved, &trajectory, Default::default()).unwrap();
    let streamed = directory.file("streamed.trr");
    let mut writer = create_trajectory_writer(
        &streamed,
        topology.clone(),
        TrajectoryWriteOptions::default(),
    )
    .unwrap();
    for sample in &trajectory {
        writer
            .write_frame(FrameView::from_sample(sample).unwrap())
            .unwrap();
    }
    writer.finish().unwrap();
    assert_eq!(fs::read(&streamed).unwrap(), fs::read(&saved).unwrap());
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Trr)
        .with_trr_options(TrrWriteOptions::default().with_precision(TrrScalarPrecision::Float64));
    write_sample_set(&path, &trajectory, options).unwrap();
    let read = read_sample_set(
        path,
        topology,
        TrajectoryOpenOptions::default()
            .with_format_hint(TrajectoryFormatHint::Explicit(TrajectoryFormat::Trr)),
    )
    .unwrap();
    for (actual, expected) in read.iter().zip(trajectory.iter()) {
        same_frame(actual, expected, 1.0e-12);
    }
}

#[test]
fn every_empty_path_writer_rejects_publication_and_removes_temporary_files() {
    use kekule::traj::io::create_trajectory_writer;

    let directory = Directory::new();
    let topology = linear_carbon_topology(3);
    for (format, extension) in [
        (TrajectoryFormat::Xyz, "xyz"),
        (TrajectoryFormat::Dcd, "dcd"),
        (TrajectoryFormat::Trr, "trr"),
        (TrajectoryFormat::Xtc, "xtc"),
    ] {
        let path = directory.file(&format!("empty.{extension}"));
        let empty =
            create_trajectory_writer(&path, topology.clone(), TrajectoryWriteOptions::new(format))
                .unwrap();
        assert_eq!(
            empty.finish(),
            Err(codec_error(
                TrajectoryIoErrorKind::InvalidFrame,
                TrajectoryIoOperation::Finish,
                format,
                &label(&path),
                None,
                &format!("The {format} trajectory has no frame. It needs at least one."),
            ))
        );
        assert_eq!(directory.names(), Vec::<String>::new());
    }
}

#[test]
fn rejected_saves_leave_no_destination_or_temporary_file() {
    let directory = Directory::new();
    let path = directory.file("out.xyz");
    let empty = SampleSet::new(linear_carbon_topology(3));
    let mut moving = Frame::new(Positions::zeros(3));
    moving.set_velocities(Some(Velocities::zeros(3))).unwrap(); // XYZ cannot carry them.
    let moving = time_ordered(
        linear_carbon_topology(3),
        [Frame::new(Positions::zeros(3)), moving],
    );
    let mut annotated = time_ordered(linear_carbon_topology(3), [Frame::new(Positions::zeros(3))]);
    annotated
        .properties_mut()
        .insert(PropertyKey::new("run").unwrap(), PropertyValue::Int(2))
        .unwrap();
    for (trajectory, expected) in [
        (
            &empty,
            codec_error(
                TrajectoryIoErrorKind::InvalidFrame,
                TrajectoryIoOperation::Finish,
                TrajectoryFormat::Xyz,
                &label(&path),
                None,
                "The XYZ trajectory has no frame. It needs at least one.",
            ),
        ),
        (
            &moving,
            codec_error(
                TrajectoryIoErrorKind::UnsupportedField,
                TrajectoryIoOperation::WriteFrame,
                TrajectoryFormat::Xyz,
                &label(&path),
                Some(1),
                "XYZ cannot store the velocities.",
            )
            .with_field(TrajectoryField::Velocities),
        ),
        (
            &annotated,
            TrajectoryIoError::unsupported_field(TrajectoryField::CollectionProperties),
        ),
    ] {
        assert_eq!(
            write_sample_set(&path, trajectory, Default::default()),
            Err(expected)
        );
        assert_eq!(directory.names(), Vec::<String>::new());
    }
}

/// One frame of three atoms whose first x coordinate is `shift` angstroms.
fn single_frame(shift: f64) -> Frame {
    use kekule::geometry::Point3;
    use kekule::units::{Quantity, ANGSTROM};

    let positions = Positions::new(Quantity::new(
        [
            Point3::new(shift, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ],
        ANGSTROM,
    ))
    .unwrap();
    Frame::new(positions)
}

/// A set of [`single_frame`].
fn single(shift: f64) -> SampleSet {
    time_ordered(linear_carbon_topology(3), [single_frame(shift)])
}

#[test]
fn a_failed_frame_write_prevents_publication_of_the_frames_before_it() {
    use kekule::traj::io::create_trajectory_writer;
    use kekule::traj::TrajectoryWriter;

    let directory = Directory::new();
    let path = directory.file("poisoned.xyz");
    let topology = linear_carbon_topology(3);
    let mut frame = single_frame(0.0).bound(&topology);
    let mut writer = create_trajectory_writer(
        &path,
        Arc::clone(&topology),
        TrajectoryWriteOptions::new(TrajectoryFormat::Xyz),
    )
    .unwrap();
    writer.write_frame(frame.view(&topology)).unwrap();
    frame.set_step(Some(1));
    assert_eq!(
        writer.write_frame(frame.view(&topology)),
        Err(codec_error(
            TrajectoryIoErrorKind::UnsupportedField,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Xyz,
            &label(&path),
            Some(1),
            "XYZ cannot store the step.",
        )
        .with_field(TrajectoryField::Step))
    );
    assert_eq!(
        writer.finish(),
        Err(codec_error(
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::Finish,
            TrajectoryFormat::Xyz,
            &label(&path),
            None,
            "An earlier frame failed to write. The trajectory is incomplete.",
        ))
    );
    assert_eq!(directory.names(), Vec::<String>::new());
}

#[test]
fn forbid_policy_never_replaces_an_existing_destination() {
    use kekule::traj::io::create_trajectory_writer;
    use kekule::traj::TrajectoryWriter;

    let directory = Directory::new();
    let path = directory.file("out.xyz");
    write_sample_set(&path, &single(0.0), Default::default()).unwrap();
    let original = fs::read(&path).unwrap();

    // An existing destination is refused before any frame is written.
    assert_eq!(
        write_sample_set(&path, &single(5.0), Default::default()),
        Err(TrajectoryIoError::io(
            TrajectoryIoOperation::Open,
            io::Error::new(
                io::ErrorKind::AlreadyExists,
                "The destination already exists."
            ),
        )
        .with_format(TrajectoryFormat::Xyz)
        .with_source_label(label(&path)))
    );
    assert_eq!(fs::read(&path).unwrap(), original);

    // A destination that appears while the writer is open is not replaced
    // either: publication itself refuses to overwrite.
    let late = directory.file("late.xyz");
    let replacement = single(5.0);
    let mut writer = create_trajectory_writer(
        &late,
        replacement.shared_topology(),
        TrajectoryWriteOptions::new(TrajectoryFormat::Xyz),
    )
    .unwrap();
    writer.write_frame(frame(&replacement, 0)).unwrap();
    fs::write(&late, &original).unwrap();
    let Err(context) = writer.finish() else {
        panic!("publishing over a file that appeared later must fail with an I/O error");
    };
    // The message is the operating system's own text, so compare the rest.
    assert_eq!(
        (
            context.kind(),
            context.operation(),
            context.format(),
            context.source_label(),
            context.io_error_kind(),
        ),
        (
            TrajectoryIoErrorKind::Io,
            Some(TrajectoryIoOperation::Finish),
            Some(TrajectoryFormat::Xyz),
            Some(label(&late).as_str()),
            Some(io::ErrorKind::AlreadyExists),
        )
    );
    assert_eq!(fs::read(&late).unwrap(), original);
    assert_eq!(directory.names(), ["late.xyz", "out.xyz"]);
}

#[test]
fn replace_policy_atomically_replaces_an_existing_destination() {
    let directory = Directory::new();
    let path = directory.file("out.xyz");
    write_sample_set(&path, &single(0.0), Default::default()).unwrap();
    let replacement = single(5.0);
    write_sample_set(
        &path,
        &replacement,
        TrajectoryWriteOptions::new(TrajectoryFormat::Xyz)
            .with_overwrite_policy(OverwritePolicy::Replace),
    )
    .unwrap();
    let read = read_sample_set(&path, replacement.shared_topology(), Default::default()).unwrap();
    assert_eq!(read.len(), 1);
    // XYZ keeps six decimals of these exactly representable coordinates.
    same_frame(read.get(0).unwrap(), replacement.get(0).unwrap(), 0.0);
    assert_eq!(directory.names(), ["out.xyz"]);
}

#[test]
fn complete_loaded_and_streaming_workflows_save_equivalent_frames() {
    use kekule::traj::analysis::FrameSuperposer;
    use kekule::traj::io::{create_trajectory_writer, open_trajectory, trr::TRR_LAMBDA_PROPERTY};
    use kekule::traj::periodic::{self, MoleculeImager, TrajectoryUnwrapper};
    use kekule::traj::{TrajectoryReader, TrajectoryWriter};
    use kekule::{
        geometry::{PeriodicCell, Vector3},
        structure::{Order, TimeAxis},
        topology::AtomSelection,
        units::{Quantity, DIMENSIONLESS, NANOMETER, PICOSECOND},
    };
    let mut recorded = Recorded::new(&[
        0xe706_edac_04c6_dc3d, // loaded, whole, unwrapped, and superposed
    ]);
    let directory = Directory::new();
    let topology = topology(&["O", "H", "H"], &[(0, 1), (0, 2)]);
    let mut source = read_sample_set(
        fixture("ase-3.26.0-water.xyz"),
        topology.clone(),
        Default::default(),
    )
    .unwrap();
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    let count = source.len();
    source
        .set_order(Some(Order::Time(
            TimeAxis::new(
                Some(Quantity::new(
                    (0..count).map(|index| index as f64).collect(),
                    PICOSECOND,
                )),
                Some((0..count).map(|index| index as u64 * 5).collect()),
            )
            .unwrap(),
        )))
        .unwrap();
    for index in 0..count {
        let mut sample = source.get_mut(index).unwrap();
        sample.set_cell(Some(cell));
        sample.set_velocities(Some(Velocities::zeros(3))).unwrap();
        sample.set_forces(Some(Forces::zeros(3))).unwrap();
    }
    // Atom slicing produces a new topology, while frame selection retains it.
    let sliced = source
        .subset(&AtomSelection::all(topology.clone()))
        .unwrap();
    let mut selected = sliced.select([0, 1, 1]).unwrap();
    assert!(std::sync::Arc::ptr_eq(
        &sliced.shared_topology(),
        &selected.shared_topology()
    ));
    let topology = selected.shared_topology();
    for index in 0..selected.len() {
        selected
            .get_mut(index)
            .unwrap()
            .properties_mut()
            .owner_mut()
            .insert(
                PropertyKey::new(TRR_LAMBDA_PROPERTY).unwrap(),
                PropertyValue::Real {
                    value: 0.0,
                    unit: DIMENSIONLESS,
                },
            )
            .unwrap();
    }
    // Save the annotated input so both workflows exercise actual file readers.
    let input = directory.file("prepared.trr");
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Trr)
        .with_trr_options(TrrWriteOptions::default().with_precision(TrrScalarPrecision::Float64));
    write_sample_set(&input, &selected, options.clone()).unwrap();
    let loaded = read_sample_set(&input, topology.clone(), Default::default()).unwrap();
    let atoms = AtomSelection::all(topology.clone());
    let mut expected = loaded;
    periodic::make_molecules_whole(&mut expected).unwrap();
    periodic::unwrap_trajectory(&mut expected).unwrap();
    expected.superpose(0, &atoms, Default::default()).unwrap();
    recorded.check(&dynamics(&expected));
    let loaded_output = directory.file("loaded.trr");
    write_sample_set(&loaded_output, &expected, options.clone()).unwrap();

    let mut reader = open_trajectory(&input, topology.clone(), Default::default()).unwrap();
    let mut frame = reader.frame_buffer();
    let imager = MoleculeImager::new(topology.clone());
    let mut unwrapper = TrajectoryUnwrapper::new(topology.clone());
    assert!(reader.read_next(&mut frame).unwrap());
    imager.make_whole(0, &mut frame).unwrap();
    unwrapper.unwrap(0, &mut frame).unwrap();
    let reference = frame.as_model_view().to_model();
    let fitter = FrameSuperposer::new(&reference, &atoms, Default::default());
    let streamed_output = directory.file("streamed.trr");
    let mut writer = create_trajectory_writer(&streamed_output, topology.clone(), options).unwrap();
    let mut index = 0;
    loop {
        fitter.superpose(index, &mut frame).unwrap();
        writer.write_frame(frame.frame_view()).unwrap();
        index += 1;
        if !reader.read_next(&mut frame).unwrap() {
            break;
        }
        imager.make_whole(index, &mut frame).unwrap();
        unwrapper.unwrap(index, &mut frame).unwrap();
    }
    writer.finish().unwrap();
    assert_eq!(index, selected.len());
    for path in [loaded_output, streamed_output] {
        let reread = read_sample_set(path, topology.clone(), Default::default()).unwrap();
        assert_eq!(reread.len(), expected.len());
        for (actual, expected) in reread.iter().zip(expected.iter()) {
            same_frame(actual, expected, 1.0e-12);
        }
        assert_eq!(steps(&reread), [Some(0), Some(5), Some(5)]);
    }
}

/// Write options under which the frames of the DCD fixture, with steps and
/// cells but no times, save as `format`: XYZ leaving the step and cell out,
/// TRR and XTC numbering their times every 2 ps.
fn converting_options(format: TrajectoryFormat) -> TrajectoryWriteOptions {
    use kekule::diagnostics::UnrepresentablePolicy;
    use kekule::traj::io::trr::TrrLambdaPolicy;
    use kekule::traj::io::xyz::XyzWriteOptions;
    use kekule::traj::io::FrameTimingPolicy;
    use kekule::units::{Quantity, PICOSECOND};

    TrajectoryWriteOptions::new(format)
        .with_xyz_options(XyzWriteOptions::default().with_cell(UnrepresentablePolicy::Omit))
        .with_dcd_options(DcdWriteOptions::default().with_cell_policy(DcdCellPolicy::Write))
        .with_trr_options(
            TrrWriteOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero),
        )
        .with_frame_metadata(UnrepresentablePolicy::Omit)
        .with_frame_timing(FrameTimingPolicy::Number {
            interval: Quantity::new(2.0, PICOSECOND),
        })
}

fn dcd_fixture() -> SampleSet {
    read_sample_set(
        fixture("mdanalysis-2.9.0-three-atoms.dcd"),
        topology(&["C", "H", "O"], &[(0, 1), (0, 2)]),
        Default::default(),
    )
    .unwrap()
}

#[test]
fn frames_without_times_save_as_every_format_by_omitting_or_numbering_state() {
    use kekule::diagnostics::{StructureField, WriteOmission};
    use kekule::units::PICOSECOND;

    let directory = Directory::new();
    let source = dcd_fixture();
    let frames = source.len();
    assert!(frames_of(&source)
        .iter()
        .all(|frame| frame.time().is_none()));
    for format in [
        TrajectoryFormat::Xyz,
        TrajectoryFormat::Dcd,
        TrajectoryFormat::Trr,
        TrajectoryFormat::Xtc,
    ] {
        let path = directory.file(&format!("converted.{}", format.name()));
        let report = write_sample_set(&path, &source, converting_options(format)).unwrap();
        let read = read_sample_set(&path, source.shared_topology(), Default::default()).unwrap();
        assert_eq!(read.len(), frames, "{format}");
        let times = frames_of(&read)
            .iter()
            .map(|frame| frame.time().map(|time| time.value_in(PICOSECOND).unwrap()))
            .collect::<Vec<_>>();
        match format {
            TrajectoryFormat::Xyz => {
                assert!(report.omissions().contains(&WriteOmission::StructureState {
                    field: StructureField::Step,
                    structures: frames,
                }));
                assert!(report
                    .omissions()
                    .contains(&WriteOmission::Cell { structures: frames }));
                assert_eq!(steps(&read), vec![None; frames]);
            }
            TrajectoryFormat::Dcd => {
                assert!(report.is_empty());
                assert_eq!(steps(&read), steps(&source));
                assert_eq!(times, vec![None; frames]);
            }
            _ => {
                assert!(report.is_empty(), "{format}");
                assert_eq!(steps(&read), steps(&source), "{format}");
                let numbered = (0..frames)
                    .map(|frame| Some(2.0 * frame as f64))
                    .collect::<Vec<_>>();
                assert_eq!(times, numbered, "{format}");
            }
        }
    }
}

#[test]
fn stream_writers_write_from_the_stream_position_what_path_writers_publish() {
    use kekule::traj::io::write_sample_set_to;

    let directory = Directory::new();
    let source = dcd_fixture();
    for format in [
        TrajectoryFormat::Xyz,
        TrajectoryFormat::Dcd,
        TrajectoryFormat::Trr,
        TrajectoryFormat::Xtc,
    ] {
        let path = directory.file(&format!("published.{}", format.name()));
        let published = write_sample_set(&path, &source, converting_options(format)).unwrap();
        let prefix = b"host header";
        let mut stream = io::Cursor::new(prefix.to_vec());
        stream.set_position(prefix.len() as u64);
        let streamed =
            write_sample_set_to(&mut stream, &source, converting_options(format)).unwrap();
        assert_eq!(streamed, published, "{format}");
        assert_eq!(stream.position(), stream.get_ref().len() as u64, "{format}");
        let bytes = stream.into_inner();
        assert_eq!(&bytes[..prefix.len()], prefix, "{format}");
        assert_eq!(bytes[prefix.len()..], fs::read(&path).unwrap(), "{format}");
    }
}

/// A set saves the bytes the trajectory it replaces saved and reads back
/// as the set that trajectory read back as, both recorded before
/// trajectories became sample sets; a time order and no distribution.
#[test]
fn sample_sets_save_and_load_as_the_trajectories_they_replace_did() {
    use kekule::structure::{Distribution, Members, Order, Partition};

    let mut recorded = Recorded::new(&[
        0xded4_b250_42dd_4141, // the DCD fixture
        0xa9d0_5a35_2281_6863, // read back from XYZ
        0x1279_da89_45bc_8d39, // read back from DCD
        0x70b9_984f_7a16_b98a, // read back from TRR
        0x2b6f_fb22_6b5f_a5a8, // read back from XTC
    ]);
    let directory = Directory::new();
    let set = dcd_fixture();
    recorded.check(&contents(&set));
    for (format, bytes) in [
        (TrajectoryFormat::Xyz, 0xe04e_bbe3_5823_f5fa),
        (TrajectoryFormat::Dcd, 0xf5a0_5c2f_1292_84d9),
        (TrajectoryFormat::Trr, 0x6486_cf10_9168_e9eb),
        (TrajectoryFormat::Xtc, 0xd547_6671_a026_a2eb),
    ] {
        let path = directory.file(&format!("samples.{}", format.name()));
        write_sample_set(&path, &set, converting_options(format)).unwrap();
        assert_eq!(
            fingerprint_bytes(&fs::read(&path).unwrap()),
            bytes,
            "{format}"
        );
        let read = read_sample_set(&path, set.shared_topology(), Default::default()).unwrap();
        recorded.check(&contents(&read));
        assert!(matches!(read.order(), Some(Order::Time(_))), "{format}");
        assert_eq!(read.distribution(), None);
    }

    // A distribution and a step or coordinate order have no trajectory
    // representation; nothing is written.
    let mut weighted = set.clone();
    weighted
        .set_distribution(Some(
            Distribution::new(Members::All, Partition::one_state(set.len())).unwrap(),
        ))
        .unwrap();
    let mut stepped = set.clone();
    stepped
        .set_order(Some(Order::Step((0..set.len() as u64).collect())))
        .unwrap();
    let mut scanned = set.clone();
    scanned
        .set_order(Some(Order::Coordinate {
            unit: kekule::units::ANGSTROM,
            values: (0..set.len()).map(|index| index as f64).collect(),
        }))
        .unwrap();
    for (refused, field) in [
        (weighted, TrajectoryField::Distribution),
        (stepped, TrajectoryField::OrderStep),
        (scanned, TrajectoryField::OrderCoordinate),
    ] {
        let path = directory.file("refused.dcd");
        let error = write_sample_set(&path, &refused, Default::default()).unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::UnsupportedField);
        assert_eq!(error.field(), Some(field));
        assert!(!path.exists());
    }
}

#[test]
fn stream_writers_need_a_format_and_refuse_to_finish_after_a_failed_frame() {
    use kekule::traj::io::create_trajectory_stream_writer;

    let source = dcd_fixture();
    assert_eq!(
        create_trajectory_stream_writer(
            io::Cursor::new(Vec::new()),
            source.shared_topology(),
            TrajectoryWriteOptions::default(),
            "host stream",
        )
        .map(|_| ()),
        Err(TrajectoryIoError::codec(
            TrajectoryIoErrorKind::UnknownFormat,
            TrajectoryIoOperation::Open,
            None,
        )
        .with_source_label("host stream")
        .with_detail("A stream has no extension. Choose the format explicitly."))
    );

    // XYZ under its default options refuses the cell of every DCD frame.
    let mut writer = create_trajectory_stream_writer(
        io::Cursor::new(Vec::new()),
        source.shared_topology(),
        TrajectoryWriteOptions::new(TrajectoryFormat::Xyz),
        "host stream",
    )
    .unwrap();
    assert_eq!(writer.format(), TrajectoryFormat::Xyz);
    let error = writer.write_frame(frame(&source, 0)).unwrap_err();
    assert_eq!(error.field(), Some(TrajectoryField::PeriodicCell));
    let poisoned = codec_error(
        TrajectoryIoErrorKind::InvalidFrame,
        TrajectoryIoOperation::WriteFrame,
        TrajectoryFormat::Xyz,
        "host stream",
        None,
        "An earlier frame failed to write. The trajectory is incomplete.",
    );
    assert_eq!(writer.write_frame(frame(&source, 0)), Err(poisoned));
    assert!(writer.writer().get_ref().is_empty());
    assert_eq!(
        writer.finish().map(|_| ()),
        Err(codec_error(
            TrajectoryIoErrorKind::InvalidFrame,
            TrajectoryIoOperation::Finish,
            TrajectoryFormat::Xyz,
            "host stream",
            None,
            "An earlier frame failed to write. The trajectory is incomplete.",
        ))
    );
}
