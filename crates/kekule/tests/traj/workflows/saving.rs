use kekule::diagnostics::LengthMismatch;
use kekule::diagnostics::{Countable, GeometryField};
use std::{
    fs, io,
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
    sync::Arc,
};

use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{Forces, GeometryView, Positions, StructureSet, Velocities};
use kekule::traj::io::dcd::{DcdCellPolicy, DcdWriteOptions};
use kekule::traj::io::trr::{TrrScalarPrecision, TrrWriteOptions};
use kekule::traj::io::{
    create_trajectory_writer, read_structure_set, write_structure_set, OverwritePolicy,
    TrajectoryFormatHint, TrajectoryOpenOptions, TrajectoryWriteOptions,
};
use kekule::traj::{
    FrameView, TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation,
    TrajectoryWriter,
};

use crate::support::{
    codec_error, frames_of, linear_carbon_topology, time_ordered, topology, Frame,
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

/// Geometry `index` of `set` as the frame a writer takes.
fn frame(set: &StructureSet, index: usize) -> FrameView<'_> {
    FrameView::from_geometry(set.get(index).unwrap()).unwrap()
}

/// The MD steps of a set's geometries.
fn steps(set: &StructureSet) -> Vec<Option<u64>> {
    frames_of(set).iter().map(Frame::step).collect()
}

fn same_frame(actual: GeometryView<'_>, expected: GeometryView<'_>, tolerance: f64) {
    // Codec precision and unit conversions permit numerical rounding. Field
    // presence, dimensions, units, steps, and every property remain exact.
    let actual = FrameView::from_geometry(actual).unwrap();
    let expected = FrameView::from_geometry(expected).unwrap();
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
            read_structure_set(fixture(source), topology.clone(), Default::default()).unwrap();
        let path = directory.file(&format!("output.{extension}"));
        match options {
            Some(options) => write_structure_set(&path, &trajectory, options).unwrap(),
            None => write_structure_set(&path, &trajectory, Default::default()).unwrap(),
        };
        let read = read_structure_set(&path, topology, Default::default()).unwrap();
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
    let trajectory = read_structure_set(
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
        write_structure_set(&path, &trajectory, Default::default()),
        Err(unknown.clone())
    );
    // Default options take the path writer's format from the extension too.
    assert_eq!(
        create_trajectory_writer(&path, topology.clone(), TrajectoryWriteOptions::default())
            .unwrap_err(),
        unknown
    );
    let saved = directory.file("saved.trr");
    write_structure_set(&saved, &trajectory, Default::default()).unwrap();
    let streamed = directory.file("streamed.trr");
    let mut writer = create_trajectory_writer(
        &streamed,
        topology.clone(),
        TrajectoryWriteOptions::default(),
    )
    .unwrap();
    for geometry in &trajectory {
        writer
            .write_frame(FrameView::from_geometry(geometry).unwrap())
            .unwrap();
    }
    writer.finish().unwrap();
    assert_eq!(fs::read(&streamed).unwrap(), fs::read(&saved).unwrap());
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Trr)
        .with_trr_options(TrrWriteOptions::default().with_precision(TrrScalarPrecision::Float64));
    write_structure_set(&path, &trajectory, options).unwrap();
    let read = read_structure_set(
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
    let empty = StructureSet::new(linear_carbon_topology(3));
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
            .with_field(GeometryField::Velocities),
        ),
        (
            &annotated,
            TrajectoryIoError::unsupported_field(GeometryField::SetProperties),
        ),
    ] {
        assert_eq!(
            write_structure_set(&path, trajectory, Default::default()),
            Err(expected)
        );
        assert_eq!(directory.names(), Vec::<String>::new());
    }
}

/// One frame of three atoms whose first x coordinate is `shift` angstroms.
fn single_frame(shift: f64) -> Frame {
    use kekule::math::Point3;
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
fn single(shift: f64) -> StructureSet {
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
            "XYZ cannot store the MD step.",
        )
        .with_field(GeometryField::Step))
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
    write_structure_set(&path, &single(0.0), Default::default()).unwrap();
    let original = fs::read(&path).unwrap();

    // An existing destination is refused before any frame is written.
    assert_eq!(
        write_structure_set(&path, &single(5.0), Default::default()),
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
    write_structure_set(&path, &single(0.0), Default::default()).unwrap();
    let replacement = single(5.0);
    write_structure_set(
        &path,
        &replacement,
        TrajectoryWriteOptions::new(TrajectoryFormat::Xyz)
            .with_overwrite_policy(OverwritePolicy::Replace),
    )
    .unwrap();
    let read =
        read_structure_set(&path, replacement.shared_topology(), Default::default()).unwrap();
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
        math::{PeriodicCell, Vector3},
        structure::{Order, TimeAxis},
        topology::AtomSelection,
        units::{Quantity, DIMENSIONLESS, NANOMETER, PICOSECOND},
    };
    let directory = Directory::new();
    let topology = topology(&["O", "H", "H"], &[(0, 1), (0, 2)]);
    let mut source = read_structure_set(
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
        let mut geometry = source.get_mut(index).unwrap();
        geometry.set_cell(Some(cell));
        geometry.set_velocities(Some(Velocities::zeros(3))).unwrap();
        geometry.set_forces(Some(Forces::zeros(3))).unwrap();
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
    write_structure_set(&input, &selected, options.clone()).unwrap();
    let loaded = read_structure_set(&input, topology.clone(), Default::default()).unwrap();
    let atoms = AtomSelection::all(topology.clone());
    let mut expected = loaded;
    periodic::make_molecules_whole(&mut expected).unwrap();
    periodic::unwrap_trajectory(&mut expected).unwrap();
    expected.superpose(0, &atoms, Default::default()).unwrap();
    // The fixture's water is whole in the 2 nm cell and moves by less than
    // half a cell, so making it whole and unwrapping keep it; each of its
    // frames is the first translated, so every frame fitted onto the first
    // lies on it. The cell, the zero velocities and forces, and the steps
    // are kept.
    let first = selected.get(0).unwrap();
    let first = first.positions().values();
    assert_eq!(expected.len(), 3);
    for (index, geometry) in expected.iter().enumerate() {
        let points = geometry.positions().values();
        for (atom, (point, reference)) in
            points.value().iter().zip(first.value().iter()).enumerate()
        {
            assert!(
                (*point - *reference).norm() < 1.0e-12,
                "geometry {index}, atom {atom}: {point:?} against {reference:?}"
            );
        }
        // The fit rotates the cell with the frame, here by no more than
        // rounding.
        let vectors = geometry.cell().unwrap().vectors().into_value();
        for (vector, expected) in vectors.into_iter().zip(cell.vectors().into_value()) {
            assert!((vector - expected).norm() < 1.0e-12, "geometry {index}");
        }
        assert_eq!(
            geometry.velocities(),
            Some(&Velocities::zeros(3)),
            "geometry {index}"
        );
        assert_eq!(
            geometry.forces(),
            Some(&Forces::zeros(3)),
            "geometry {index}"
        );
    }
    assert_eq!(steps(&expected), [Some(0), Some(5), Some(5)]);
    let loaded_output = directory.file("loaded.trr");
    write_structure_set(&loaded_output, &expected, options.clone()).unwrap();

    let mut reader = open_trajectory(&input, topology.clone(), Default::default()).unwrap();
    let mut frame = reader.frame_buffer();
    let imager = MoleculeImager::new(topology.clone());
    let mut unwrapper = TrajectoryUnwrapper::new(topology.clone());
    assert!(reader.read_next(&mut frame).unwrap());
    imager.make_whole(0, &mut frame).unwrap();
    unwrapper.unwrap(0, &mut frame).unwrap();
    let reference = frame.as_structure_view().to_structure();
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
        let reread = read_structure_set(path, topology.clone(), Default::default()).unwrap();
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

fn dcd_fixture() -> StructureSet {
    read_structure_set(
        fixture("mdanalysis-2.9.0-three-atoms.dcd"),
        topology(&["C", "H", "O"], &[(0, 1), (0, 2)]),
        Default::default(),
    )
    .unwrap()
}

#[test]
fn frames_without_times_save_as_every_format_by_omitting_or_numbering_state() {
    use kekule::diagnostics::{GeometryField, WriteOmission};
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
        let report = write_structure_set(&path, &source, converting_options(format)).unwrap();
        let read = read_structure_set(&path, source.shared_topology(), Default::default()).unwrap();
        assert_eq!(read.len(), frames, "{format}");
        let times = frames_of(&read)
            .iter()
            .map(|frame| frame.time().map(|time| time.value_in(PICOSECOND).unwrap()))
            .collect::<Vec<_>>();
        match format {
            TrajectoryFormat::Xyz => {
                assert!(report.omissions().contains(&WriteOmission::GeometryState {
                    field: GeometryField::Step,
                    geometries: frames,
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
    use kekule::traj::io::write_structure_set_to;

    let directory = Directory::new();
    let source = dcd_fixture();
    for format in [
        TrajectoryFormat::Xyz,
        TrajectoryFormat::Dcd,
        TrajectoryFormat::Trr,
        TrajectoryFormat::Xtc,
    ] {
        let path = directory.file(&format!("published.{}", format.name()));
        let published = write_structure_set(&path, &source, converting_options(format)).unwrap();
        let prefix = b"host header";
        let mut stream = io::Cursor::new(prefix.to_vec());
        stream.set_position(prefix.len() as u64);
        let streamed =
            write_structure_set_to(&mut stream, &source, converting_options(format)).unwrap();
        assert_eq!(streamed, published, "{format}");
        assert_eq!(stream.position(), stream.get_ref().len() as u64, "{format}");
        let bytes = stream.into_inner();
        assert_eq!(&bytes[..prefix.len()], prefix, "{format}");
        assert_eq!(bytes[prefix.len()..], fs::read(&path).unwrap(), "{format}");
    }
}

/// A set saves as every format and reads back as the set it saved, to the
/// format's precision, with a time order.
#[test]
fn structure_sets_save_and_load_as_every_format() {
    use kekule::structure::Order;

    let directory = Directory::new();
    let set = dcd_fixture();
    assert_eq!(set.topology().atom_count(), 3);
    assert!(set.len() >= 2);
    // The largest difference in nanometres each format's precision allows:
    // XYZ writes eight decimals of ångströms, DCD single-precision ångströms
    // as the fixture stores them, TRR single-precision nanometres, and XTC a
    // thousandth of a nanometre.
    for (format, tolerance) in [
        (TrajectoryFormat::Xyz, 1.0e-9),
        (TrajectoryFormat::Dcd, 0.0),
        (TrajectoryFormat::Trr, 1.0e-6),
        (TrajectoryFormat::Xtc, 5.0e-4),
    ] {
        let path = directory.file(&format!("geometries.{}", format.name()));
        write_structure_set(&path, &set, converting_options(format)).unwrap();
        // Writing again gives the same bytes.
        let again = directory.file(&format!("again.{}", format.name()));
        write_structure_set(&again, &set, converting_options(format)).unwrap();
        assert_eq!(
            fs::read(&path).unwrap(),
            fs::read(&again).unwrap(),
            "{format}"
        );
        let read = read_structure_set(&path, set.shared_topology(), Default::default()).unwrap();
        assert!(matches!(read.order(), Some(Order::Time(_))), "{format}");
        assert_eq!(read.len(), set.len(), "{format}");
        assert!(read.properties().is_empty(), "{format}");
        for (index, (read, stored)) in read.iter().zip(&set).enumerate() {
            let values = read.positions().values();
            let expected = stored.positions().values();
            for (atom, (point, stored)) in values
                .value()
                .iter()
                .zip(expected.value().iter())
                .enumerate()
            {
                assert!(
                    (*point - *stored).norm() <= tolerance,
                    "{format}: geometry {index}, atom {atom}: {point:?} against {stored:?}"
                );
            }
            // XYZ leaves the cell out; the others store it.
            match format {
                TrajectoryFormat::Xyz => assert!(read.cell().is_none()),
                _ => assert!(read.cell().is_some(), "{format}: geometry {index}"),
            }
            assert!(read.velocities().is_none() && read.forces().is_none());
        }
    }

    // A step or coordinate order has no trajectory representation; nothing
    // is written.
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
        (stepped, GeometryField::OrderStep),
        (scanned, GeometryField::OrderCoordinate),
    ] {
        let path = directory.file("refused.dcd");
        let error = write_structure_set(&path, &refused, Default::default()).unwrap_err();
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
    assert_eq!(error.field(), Some(GeometryField::Cell));
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

#[test]
fn a_standalone_geometry_bound_to_a_solvated_topology_writes_as_a_frame() {
    use kekule::math::{PeriodicCell, Point3, Vector3};
    use kekule::structure::{Geometry, GeometryError, StructureSetError};
    use kekule::topology::TopologyBuilder;
    use kekule::units::{Quantity, NANOMETER, PICOSECOND};

    // Ten thousand waters: 30,000 atoms and 20,000 bonds.
    let waters = 10_000;
    let mut builder = TopologyBuilder::new();
    let water = builder
        .add_molecule_definition(kekule::smiles::to_molecules("[H]O[H]").unwrap().remove(0))
        .unwrap();
    for _ in 0..waters {
        builder.add_instance(water).unwrap();
    }
    let topology = Arc::new(builder.build().unwrap());
    assert_eq!(topology.bond_count(), 2 * waters);
    let points = (0..topology.atom_count())
        .map(|atom| {
            let atom = atom as f64;
            Point3::new(atom * 0.001, (atom * 0.37) % 3.0, (atom * 0.11) % 3.0)
        })
        .collect::<Vec<_>>();
    let positions = Positions::new(Quantity::new(points, NANOMETER)).unwrap();
    // A solvated system in a box, as XTC stores it.
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(30.0, 3.0, 3.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    let mut geometry = Geometry::new(positions.clone()).with_cell(cell);
    assert_eq!(
        FrameView::new(&topology, &geometry, None, Some(5)).unwrap_err(),
        StructureSetError::Geometry(GeometryError::BondCountMismatch(LengthMismatch::new(
            Countable::Row,
            2 * waters,
            0
        )))
    );
    let stored = geometry.positions().values().value().as_ptr();
    geometry.bind(&topology).unwrap();
    // Binding gives the geometry bond rows without copying its positions.
    assert_eq!(geometry.positions().values().value().as_ptr(), stored);

    let time = Quantity::new(10.0, PICOSECOND);
    let directory = Directory::new();
    let path = directory.file("water.xtc");
    let mut writer = create_trajectory_writer(
        &path,
        Arc::clone(&topology),
        TrajectoryWriteOptions::new(TrajectoryFormat::Xtc),
    )
    .unwrap();
    writer
        .write_frame(FrameView::new(&topology, &geometry, Some(time), Some(5)).unwrap())
        .unwrap();
    writer.finish().unwrap();
    let read = read_structure_set(
        &path,
        Arc::clone(&topology),
        TrajectoryOpenOptions::default(),
    )
    .unwrap();
    assert_eq!(read.len(), 1);
    let frame = read.get(0).unwrap();
    assert_eq!(frame.atom_count(), topology.atom_count());
    for (read, written) in frame
        .positions()
        .values()
        .value()
        .iter()
        .zip(positions.values().value().iter())
    {
        assert!((*read - *written).norm() < 1.0e-3);
    }
}

#[test]
fn a_failed_binding_changes_nothing_and_a_binding_copies_a_shared_geometry() {
    use kekule::structure::{Geometry, GeometryError, StructureSetError};
    use kekule::units::{Quantity, NANOMETER};

    let ethanol = kekule::smiles::to_topology("CCO").unwrap();
    let water = kekule::smiles::to_topology("O").unwrap();
    let points = (0..ethanol.atom_count())
        .map(|atom| kekule::math::Point3::new(atom as f64, 0.5, -0.25))
        .collect::<Vec<_>>();
    let mut geometry = Geometry::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap());
    let before = geometry.clone();
    let stored = geometry.positions().values().value().as_ptr();

    // Another atom count fails and leaves the geometry as it was.
    assert!(matches!(
        geometry.bind(&water).unwrap_err(),
        GeometryError::Array(_)
    ));
    assert_eq!(geometry, before);
    assert!(geometry.shares_storage_with(&before));
    assert_eq!(geometry.positions().values().value().as_ptr(), stored);

    // Bound, the geometry refuses the bond rows of another topology with as
    // many atoms and keeps its own.
    geometry.bind(&ethanol).unwrap();
    let bound = geometry.clone();
    let dimethyl_ether = kekule::smiles::to_topology("COC").unwrap();
    let ring = kekule::smiles::to_topology("C1CO1").unwrap();
    assert_eq!(ring.atom_count(), ethanol.atom_count());
    assert!(matches!(
        geometry.bind(&ring).unwrap_err(),
        GeometryError::BondCountMismatch { .. }
    ));
    assert_eq!(geometry, bound);
    assert!(geometry.shares_storage_with(&bound));
    // A topology with the same bond count binds as the geometry already is.
    geometry.bind(&dimethyl_ether).unwrap();
    assert!(geometry.shares_storage_with(&bound));

    // Binding a geometry that shares its storage copies it once and leaves the
    // other geometry unbound and unchanged.
    let mut copy = before.clone();
    assert!(copy.shares_storage_with(&before));
    copy.bind(&ethanol).unwrap();
    assert!(!copy.shares_storage_with(&before));
    assert_eq!(copy.positions(), before.positions());
    assert_eq!(
        FrameView::new(&ethanol, &before, None, None).unwrap_err(),
        StructureSetError::Geometry(GeometryError::BondCountMismatch(LengthMismatch::new(
            Countable::Row,
            ethanol.bond_count(),
            0
        )))
    );
    assert!(FrameView::new(&ethanol, &copy, None, None).is_ok());
    assert_eq!(before.positions().values().value().as_ptr(), stored);
}
