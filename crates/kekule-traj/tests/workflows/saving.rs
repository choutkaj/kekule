use std::{
    fs, io,
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
};

use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::Positions;
use kekule::structure::{Forces, Trajectory, TrajectoryFrame, TrajectoryFrameView, Velocities};
use kekule_traj::io::dcd::DcdWriteOptions;
use kekule_traj::io::trr::{TrrScalarPrecision, TrrWriteOptions};
use kekule_traj::io::{
    read_trajectory, read_trajectory_with_options, write_trajectory, write_trajectory_with_options,
    OverwritePolicy, TrajectoryFormatHint, TrajectoryOpenOptions, TrajectoryWriteOptions,
};
use kekule_traj::{
    TrajectoryCodecErrorContext, TrajectoryCodecErrorKind, TrajectoryFormat, TrajectoryIoError,
    TrajectoryIoErrorContext, TrajectoryIoOperation,
};

use crate::support;
use support::{codec_error, linear_carbon_topology, topology};

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
        .join("tests/fixtures")
        .join(name)
}

fn same_frame(actual: TrajectoryFrameView<'_>, expected: TrajectoryFrameView<'_>, tolerance: f64) {
    // Codec precision and unit conversions permit numerical rounding. Field
    // presence, dimensions, units, steps, and every property remain exact.
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
                TrajectoryWriteOptions::new(TrajectoryFormat::Dcd)
                    .with_dcd_options(DcdWriteOptions::default().with_cells(true)),
            ),
        ),
    ] {
        let trajectory = read_trajectory(fixture(source), topology.clone()).unwrap();
        let path = directory.file(&format!("output.{extension}"));
        match options {
            Some(options) => write_trajectory_with_options(&path, &trajectory, options).unwrap(),
            None => write_trajectory(&path, &trajectory).unwrap(),
        }
        let read = read_trajectory(&path, topology).unwrap();
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
    let trajectory = read_trajectory(
        fixture("mdanalysis-2.9.0-three-atoms.trr"),
        topology.clone(),
    )
    .unwrap();
    let path = directory.file("trajectory.data");
    assert_eq!(
        write_trajectory(&path, &trajectory),
        Err(TrajectoryCodecErrorContext::new(
            TrajectoryCodecErrorKind::UnknownFormat,
            TrajectoryIoOperation::Open,
            None,
        )
        .with_source_label(label(&path))
        .with_detail(
            "destination extension does not identify a supported format; supply explicit write options"
        )
        .into())
    );
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Trr)
        .with_trr_options(TrrWriteOptions::default().with_precision(TrrScalarPrecision::Float64));
    write_trajectory_with_options(&path, &trajectory, options).unwrap();
    let read = read_trajectory_with_options(
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
    use kekule_traj::io::create_trajectory_writer;

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
                TrajectoryCodecErrorKind::InvalidFrame,
                TrajectoryIoOperation::Finish,
                format,
                &label(&path),
                None,
                &format!("{format} production trajectories must contain at least one frame"),
            ))
        );
        assert_eq!(directory.names(), Vec::<String>::new());
    }
}

#[test]
fn rejected_saves_leave_no_destination_or_temporary_file() {
    let directory = Directory::new();
    let path = directory.file("out.xyz");
    let empty = Trajectory::new(linear_carbon_topology(3));
    let mut stepped = empty.clone();
    stepped
        .push(TrajectoryFrame::new(Positions::zeros(3)))
        .unwrap();
    let mut late = TrajectoryFrame::new(Positions::zeros(3));
    late.set_step(Some(1)); // XYZ cannot carry a step.
    stepped.push(late).unwrap();
    let mut annotated = empty.clone();
    annotated
        .push(TrajectoryFrame::new(Positions::zeros(3)))
        .unwrap();
    annotated
        .properties_mut()
        .insert(PropertyKey::new("run").unwrap(), PropertyValue::Int(2))
        .unwrap();
    for (trajectory, expected) in [
        (
            &empty,
            codec_error(
                TrajectoryCodecErrorKind::InvalidFrame,
                TrajectoryIoOperation::Finish,
                TrajectoryFormat::Xyz,
                &label(&path),
                None,
                "XYZ production trajectories must contain at least one frame",
            ),
        ),
        (
            &stepped,
            codec_error(
                TrajectoryCodecErrorKind::UnsupportedField,
                TrajectoryIoOperation::WriteFrame,
                TrajectoryFormat::Xyz,
                &label(&path),
                Some(1),
                "XYZ cannot preserve step",
            ),
        ),
        (
            &annotated,
            TrajectoryIoError::UnsupportedField("collection properties"),
        ),
    ] {
        assert_eq!(write_trajectory(&path, trajectory), Err(expected));
        assert_eq!(directory.names(), Vec::<String>::new());
    }
}

/// One frame of three atoms whose first x coordinate is `shift` angstroms.
fn single_frame(shift: f64) -> Trajectory {
    use kekule::geometry::Point3;
    use kekule::units::{Quantity, ANGSTROM};

    let mut trajectory = Trajectory::new(linear_carbon_topology(3));
    let positions = Positions::new(Quantity::new(
        [
            Point3::new(shift, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
        ],
        ANGSTROM,
    ))
    .unwrap();
    trajectory.push(TrajectoryFrame::new(positions)).unwrap();
    trajectory
}

#[test]
fn a_failed_frame_write_prevents_publication_of_the_frames_before_it() {
    use kekule_traj::io::create_trajectory_writer;
    use kekule_traj::TrajectoryWriter;

    let directory = Directory::new();
    let path = directory.file("poisoned.xyz");
    let mut trajectory = single_frame(0.0);
    let mut writer = create_trajectory_writer(
        &path,
        trajectory.shared_topology(),
        TrajectoryWriteOptions::new(TrajectoryFormat::Xyz),
    )
    .unwrap();
    writer.write_frame(trajectory.get(0).unwrap()).unwrap();
    trajectory.get_mut(0).unwrap().set_step(Some(1));
    assert_eq!(
        writer.write_frame(trajectory.get(0).unwrap()),
        Err(codec_error(
            TrajectoryCodecErrorKind::UnsupportedField,
            TrajectoryIoOperation::WriteFrame,
            TrajectoryFormat::Xyz,
            &label(&path),
            Some(1),
            "XYZ cannot preserve step",
        ))
    );
    assert_eq!(
        writer.finish(),
        Err(codec_error(
            TrajectoryCodecErrorKind::InvalidFrame,
            TrajectoryIoOperation::Finish,
            TrajectoryFormat::Xyz,
            &label(&path),
            None,
            "cannot publish a trajectory after an earlier frame-write failure",
        ))
    );
    assert_eq!(directory.names(), Vec::<String>::new());
}

#[test]
fn forbid_policy_never_replaces_an_existing_destination() {
    use kekule_traj::io::create_trajectory_writer;
    use kekule_traj::TrajectoryWriter;

    let directory = Directory::new();
    let path = directory.file("out.xyz");
    write_trajectory(&path, &single_frame(0.0)).unwrap();
    let original = fs::read(&path).unwrap();

    // An existing destination is refused before any frame is written.
    assert_eq!(
        write_trajectory(&path, &single_frame(5.0)),
        Err(TrajectoryIoErrorContext::new(
            TrajectoryIoOperation::Open,
            io::ErrorKind::AlreadyExists,
            "destination already exists",
        )
        .with_format(TrajectoryFormat::Xyz)
        .with_source_label(label(&path))
        .into())
    );
    assert_eq!(fs::read(&path).unwrap(), original);

    // A destination that appears while the writer is open is not replaced
    // either: publication itself refuses to overwrite.
    let late = directory.file("late.xyz");
    let replacement = single_frame(5.0);
    let mut writer = create_trajectory_writer(
        &late,
        replacement.shared_topology(),
        TrajectoryWriteOptions::new(TrajectoryFormat::Xyz),
    )
    .unwrap();
    writer.write_frame(replacement.get(0).unwrap()).unwrap();
    fs::write(&late, &original).unwrap();
    let Err(TrajectoryIoError::Io(context)) = writer.finish() else {
        panic!("publishing over a file that appeared later must fail with an I/O error");
    };
    // The message is the operating system's own text, so compare the rest.
    assert_eq!(
        (
            context.operation(),
            context.format(),
            context.source_label(),
            context.error_kind(),
        ),
        (
            TrajectoryIoOperation::Finish,
            Some(TrajectoryFormat::Xyz),
            Some(label(&late).as_str()),
            io::ErrorKind::AlreadyExists,
        )
    );
    assert_eq!(fs::read(&late).unwrap(), original);
    assert_eq!(directory.names(), ["late.xyz", "out.xyz"]);
}

#[test]
fn replace_policy_atomically_replaces_an_existing_destination() {
    let directory = Directory::new();
    let path = directory.file("out.xyz");
    write_trajectory(&path, &single_frame(0.0)).unwrap();
    let replacement = single_frame(5.0);
    write_trajectory_with_options(
        &path,
        &replacement,
        TrajectoryWriteOptions::new(TrajectoryFormat::Xyz)
            .with_overwrite_policy(OverwritePolicy::Replace),
    )
    .unwrap();
    let read = read_trajectory(&path, replacement.shared_topology()).unwrap();
    assert_eq!(read.len(), 1);
    // XYZ keeps six decimals of these exactly representable coordinates.
    same_frame(read.get(0).unwrap(), replacement.get(0).unwrap(), 0.0);
    assert_eq!(directory.names(), ["out.xyz"]);
}

#[test]
fn complete_loaded_and_streaming_workflows_save_equivalent_frames() {
    use kekule::{
        geometry::{PeriodicCell, Vector3},
        topology::AtomSelection,
        units::{Quantity, DIMENSIONLESS, NANOMETER, PICOSECOND},
    };
    use kekule_traj::analysis::FrameSuperposer;
    use kekule_traj::io::{create_trajectory_writer, open_trajectory, trr::TRR_LAMBDA_PROPERTY};
    use kekule_traj::periodic::{self, MoleculeImager, TrajectoryUnwrapper};
    use kekule_traj::{TrajectoryReader, TrajectoryWriter};
    let directory = Directory::new();
    let topology = topology(&["O", "H", "H"], &[(0, 1), (0, 2)]);
    let mut source = read_trajectory(fixture("ase-3.26.0-water.xyz"), topology.clone()).unwrap();
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    for index in 0..source.len() {
        let mut frame = source.get_mut(index).unwrap();
        frame.set_cell(Some(cell));
        frame
            .set_time(Some(Quantity::new(index as f64, PICOSECOND)))
            .unwrap();
        frame.set_step(Some(index as u64 * 5));
        frame.set_velocities(Some(Velocities::zeros(3))).unwrap();
        frame.set_forces(Some(Forces::zeros(3))).unwrap();
    }
    // Atom slicing produces a new topology, while frame selection retains it.
    let sliced = source.subset(&AtomSelection::all(&topology)).unwrap();
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
    write_trajectory_with_options(&input, &selected, options.clone()).unwrap();
    let loaded = read_trajectory(&input, topology.clone()).unwrap();
    let atoms = AtomSelection::all(&topology);
    let mut expected = loaded;
    periodic::make_molecules_whole(&mut expected).unwrap();
    periodic::unwrap(&mut expected).unwrap();
    expected.superpose(0, &atoms, Default::default()).unwrap();
    let loaded_output = directory.file("loaded.trr");
    write_trajectory_with_options(&loaded_output, &expected, options.clone()).unwrap();

    let mut reader = open_trajectory(&input, topology.clone()).unwrap();
    let mut frame = reader.frame_buffer();
    let imager = MoleculeImager::new(topology.clone());
    let mut unwrapper = TrajectoryUnwrapper::new(topology.clone());
    assert!(reader.read_next(&mut frame).unwrap());
    imager.make_whole(0, &mut frame).unwrap();
    unwrapper.unwrap(0, &mut frame).unwrap();
    let reference = frame.frame_view().to_model();
    let fitter = FrameSuperposer::new(&reference, &atoms);
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
        let reread = read_trajectory(path, topology.clone()).unwrap();
        assert_eq!(reread.len(), expected.len());
        for (actual, expected) in reread.iter().zip(expected.iter()) {
            same_frame(actual, expected, 1.0e-12);
        }
        assert_eq!(
            reread.iter().map(|frame| frame.step()).collect::<Vec<_>>(),
            [Some(0), Some(5), Some(5)]
        );
    }
}
