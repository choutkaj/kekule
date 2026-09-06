use std::{
    fs,
    path::{Path, PathBuf},
    sync::atomic::{AtomicU64, Ordering},
};

use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::Positions;
use kekule_traj::io::dcd::DcdWriteOptions;
use kekule_traj::io::trr::{TrrScalarPrecision, TrrWriteOptions};
use kekule_traj::io::{
    read_trajectory, read_trajectory_with_options, write_trajectory, write_trajectory_with_options,
    OverwritePolicy, TrajectoryFormatHint, TrajectoryOpenOptions, TrajectoryWriteOptions,
};
use kekule_traj::{
    Trajectory, TrajectoryError, TrajectoryFormat, TrajectoryFrame, TrajectoryFrameView,
};

mod support;
use support::{linear_carbon_topology, topology};

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
    fn count(&self) -> usize {
        fs::read_dir(&self.0).unwrap().count()
    }
}
impl Drop for Directory {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
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
            .to_value()
            .into_iter()
            .zip(b.vectors().to_value())
        {
            close(a.x, b.x);
            close(a.y, b.y);
            close(a.z, b.z);
        }
    }
    for (actual, expected) in [
        (actual.velocities(), expected.velocities()),
        (actual.forces(), expected.forces()),
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
        close(a.to_value(), b.to_value());
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
        for (actual, expected) in read.frames().zip(trajectory.frames()) {
            same_frame(actual, expected, 1.0e-6);
        }
    }
    assert_eq!(directory.count(), 4);
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
    assert!(write_trajectory(&path, &trajectory).is_err());
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
    for (actual, expected) in read.frames().zip(trajectory.frames()) {
        same_frame(actual, expected, 1.0e-12);
    }
}

#[test]
fn failed_saves_leave_no_output_or_temporary_files_and_preserve_existing_destinations() {
    let directory = Directory::new();
    let topology = linear_carbon_topology(3);
    let path = directory.file("out.xyz");
    let mut trajectory = Trajectory::new(topology);
    assert!(write_trajectory(&path, &trajectory).is_err());
    assert_eq!(directory.count(), 0);
    trajectory
        .push(TrajectoryFrame::new(Positions::zeros(3)))
        .unwrap();
    let mut late = TrajectoryFrame::new(Positions::zeros(3));
    late.set_step(Some(1)); // XYZ cannot carry this state.
    trajectory.push(late).unwrap();
    assert!(write_trajectory(&path, &trajectory).is_err());
    assert_eq!(directory.count(), 0);
    trajectory.frame_mut(1).unwrap().set_step(None);
    trajectory
        .insert_property(PropertyKey::new("run").unwrap(), PropertyValue::Int(2))
        .unwrap();
    assert_eq!(
        write_trajectory(&path, &trajectory),
        Err(TrajectoryError::UnsupportedField("collection properties"))
    );
    assert_eq!(directory.count(), 0);
    trajectory.clear_properties();
    write_trajectory(&path, &trajectory).unwrap();
    let before = fs::read(&path).unwrap();
    assert!(write_trajectory(&path, &trajectory).is_err());
    assert_eq!(fs::read(&path).unwrap(), before);
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xyz)
        .with_overwrite_policy(OverwritePolicy::Replace);
    let replaced = write_trajectory_with_options(&path, &trajectory, options);
    if cfg!(windows) {
        assert!(replaced.is_err());
    } else {
        replaced.unwrap();
    }
    assert_eq!(fs::read(&path).unwrap(), before);
    assert_eq!(directory.count(), 1);
}

#[test]
fn complete_loaded_and_streaming_workflows_save_equivalent_frames() {
    use kekule::{
        geometry::{PeriodicCell, Vector3},
        topology::AtomSelection,
        units::{Quantity, DIMENSIONLESS, NANOMETER, PICOSECOND},
    };
    use kekule_traj::{
        analysis::FrameSuperposer,
        io::{create_trajectory_writer, open_trajectory, trr::TRR_LAMBDA_PROPERTY},
        periodic::{MoleculeImager, TrajectoryUnwrapper},
        Forces, TrajectoryReader, TrajectoryWriter, Velocities,
    };
    let directory = Directory::new();
    let topology = topology(&["O", "H", "H"], &[(0, 1), (0, 2)]);
    let mut source = read_trajectory(fixture("ase-3.26.0-water.xyz"), topology.clone()).unwrap();
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    for index in 0..source.len() {
        let mut frame = source.frame_mut(index).unwrap();
        frame.set_cell(Some(cell));
        frame
            .set_time(Some(Quantity::new(index as f64, PICOSECOND)))
            .unwrap();
        frame.set_step(Some(index as u64 * 5));
        frame.set_velocities(Some(Velocities::zeros(3))).unwrap();
        frame.set_forces(Some(Forces::zeros(3))).unwrap();
    }
    // Atom slicing produces a new topology, while frame selection retains it.
    let sliced = source.slice(&AtomSelection::all(&topology)).unwrap();
    let mut selected = sliced.select_frames([0, 1, 1]).unwrap();
    assert!(std::sync::Arc::ptr_eq(
        &sliced.shared_topology(),
        &selected.shared_topology()
    ));
    let topology = selected.shared_topology();
    for index in 0..selected.len() {
        selected
            .frame_mut(index)
            .unwrap()
            .insert_property(
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
    let whole = loaded.make_molecules_whole().unwrap();
    let continuous = whole.unwrap().unwrap();
    let expected = continuous.superpose_to_frame(0, &atoms).unwrap();
    let loaded_output = directory.file("loaded.trr");
    write_trajectory_with_options(&loaded_output, &expected, options.clone()).unwrap();

    let mut reader = open_trajectory(&input, topology.clone()).unwrap();
    let mut frame = reader.frame_buffer();
    let imager = MoleculeImager::new(topology.clone());
    let mut unwrapper = TrajectoryUnwrapper::new(topology.clone());
    assert!(reader.read_next(&mut frame).unwrap());
    imager.make_whole_in_place(0, &mut frame).unwrap();
    unwrapper.unwrap_in_place(0, &mut frame).unwrap();
    let reference = frame.frame_view().to_frame();
    let fitter = FrameSuperposer::new(reference.view(&topology).unwrap(), &atoms);
    let streamed_output = directory.file("streamed.trr");
    let mut writer = create_trajectory_writer(&streamed_output, topology.clone(), options).unwrap();
    let mut index = 0;
    loop {
        fitter.superpose_in_place(index, &mut frame).unwrap();
        writer.write_frame(frame.frame_view()).unwrap();
        index += 1;
        if !reader.read_next(&mut frame).unwrap() {
            break;
        }
        imager.make_whole_in_place(index, &mut frame).unwrap();
        unwrapper.unwrap_in_place(index, &mut frame).unwrap();
    }
    writer.finish().unwrap();
    assert_eq!(index, selected.len());
    for path in [loaded_output, streamed_output] {
        let reread = read_trajectory(path, topology.clone()).unwrap();
        assert_eq!(reread.len(), expected.len());
        for (actual, expected) in reread.frames().zip(expected.frames()) {
            same_frame(actual, expected, 1.0e-12);
        }
        assert_eq!(
            reread
                .frames()
                .map(|frame| frame.step())
                .collect::<Vec<_>>(),
            [Some(0), Some(5), Some(5)]
        );
    }
}
