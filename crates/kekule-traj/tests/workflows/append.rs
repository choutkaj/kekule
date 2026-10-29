//! Continuing a trajectory file after a restart: the frames a checkpoint
//! covers are kept, the rest removed, and new frames appended, so the file
//! ends up as one uninterrupted run would have written it.

use std::fs::{self, File};
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, OnceLock};

use kekule::diagnostics::{Cancellation, UnrepresentablePolicy};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
use kekule::topology::Topology;
use kekule::units::{Quantity, DIMENSIONLESS, NANOMETER, PER_NANOMETER, PICOSECOND};
use kekule_traj::io::dcd::{DcdCellPolicy, DcdWriteOptions};
use kekule_traj::io::trr::{TrrLambdaPolicy, TrrScalarPrecision, TrrWriteOptions};
use kekule_traj::io::xtc::{XtcCellPolicy, XtcReadOptions, XtcWriteOptions};
use kekule_traj::io::xyz::XyzWriteOptions;
use kekule_traj::io::{
    append_trajectory, create_trajectory_stream_writer, open_indexed_trajectory, KeepFrames,
    TrajectoryAppendOptions, TrajectoryAppendReport, TrajectoryFileFingerprint,
    TrajectoryOpenOptions, TrajectoryWriteOptions,
};
use kekule_traj::{
    SeekableTrajectoryReader, TrajectoryFormat, TrajectoryIoErrorKind, TrajectoryReader,
    TrajectoryWriter,
};

use crate::support::linear_carbon_topology;

const FORMATS: [TrajectoryFormat; 4] = [
    TrajectoryFormat::Xtc,
    TrajectoryFormat::Trr,
    TrajectoryFormat::Dcd,
    TrajectoryFormat::Xyz,
];

/// Steps between saved frames.
const INTERVAL: u64 = 10;

/// A scratch directory removed on drop.
struct Directory(PathBuf);

impl Directory {
    fn new() -> Self {
        static SEQUENCE: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "kekule-append-{}-{}",
            std::process::id(),
            SEQUENCE.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }

    /// A path for a trajectory of `format`, named with its extension.
    fn file(&self, format: TrajectoryFormat) -> PathBuf {
        self.0
            .join(format!("run.{}", format.name().to_ascii_lowercase()))
    }
}

impl Drop for Directory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// One topology for every reader and writer, so their layouts agree.
fn topology() -> Arc<Topology> {
    static TOPOLOGY: OnceLock<Arc<Topology>> = OnceLock::new();
    Arc::clone(TOPOLOGY.get_or_init(|| linear_carbon_topology(12)))
}

/// Saved frame `index` of a run, at step `index × INTERVAL`, its atoms
/// displaced by `shift`.
fn frame_with(index: u64, shift: f64) -> TrajectoryFrame {
    let moved = index as f64 + shift;
    let points = (0..12)
        .map(|atom| {
            let atom = atom as f64;
            Point3::new(
                0.15 * atom + 0.011 * moved,
                0.2 + 0.007 * moved * atom.sin(),
                0.3 - 0.005 * moved,
            )
        })
        .collect::<Vec<_>>();
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(3.0, 3.5, 4.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    let mut frame =
        TrajectoryFrame::new(Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap())
            .with_cell(cell);
    let step = index * INTERVAL;
    frame.set_step(Some(step));
    frame
        .set_time(Some(Quantity::new(0.02 * step as f64, PICOSECOND)))
        .unwrap();
    frame
}

fn frame(index: u64) -> TrajectoryFrame {
    frame_with(index, 0.0)
}

/// Write options keeping what each format holds of [`frame`].
fn options(format: TrajectoryFormat) -> TrajectoryWriteOptions {
    TrajectoryWriteOptions::new(format)
        .with_dcd_options(
            DcdWriteOptions::default()
                .with_cell_policy(DcdCellPolicy::Write)
                .with_step_sequence(0, INTERVAL)
                .with_frame_metadata(UnrepresentablePolicy::Omit),
        )
        .with_trr_options(
            TrrWriteOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero),
        )
        .with_xyz_options(
            XyzWriteOptions::default()
                .with_cell(UnrepresentablePolicy::Omit)
                .with_frame_metadata(UnrepresentablePolicy::Omit),
        )
}

/// Writes `frames` to `writer`, flushing after each as a live run does.
fn write_frames(
    writer: &mut impl TrajectoryWriter,
    frames: impl IntoIterator<Item = TrajectoryFrame>,
) {
    let frames = frames.into_iter().collect::<Vec<_>>();
    let trajectory = Trajectory::from_items(topology(), frames).unwrap();
    for frame in &trajectory {
        writer.write_frame(frame).unwrap();
        writer.flush().unwrap();
    }
}

/// The file one uninterrupted run writes of frames `0..count`.
fn straight_through(format: TrajectoryFormat, path: &Path, count: u64) -> Vec<u8> {
    let file = File::create(path).unwrap();
    let mut writer =
        create_trajectory_stream_writer(BufWriter::new(file), topology(), options(format), "run")
            .unwrap();
    write_frames(&mut writer, (0..count).map(frame));
    drop(writer.finish().unwrap());
    fs::read(path).unwrap()
}

/// Where each frame of the complete file at `path` ends, and where the first
/// starts: after a DCD's header, at the start for the other formats.
fn frame_ends(path: &Path) -> (u64, Vec<u64>) {
    frame_ends_read(path, Default::default())
}

/// [`frame_ends`], read with `options`.
fn frame_ends_read(path: &Path, options: TrajectoryOpenOptions) -> (u64, Vec<u64>) {
    let index = open_indexed_trajectory(path, topology(), options)
        .unwrap()
        .to_frame_index();
    let starts = index.frame_offsets();
    let mut ends = starts[1..].to_vec();
    ends.push(index.fingerprint().byte_count());
    (starts[0], ends)
}

/// Appends to `path`, keeping what `keep` selects, then writes frames
/// `from..count` and finishes, returning the report.
fn continue_run(
    format: TrajectoryFormat,
    path: &Path,
    keep: KeepFrames,
    from: u64,
    count: u64,
) -> TrajectoryAppendReport {
    let (mut writer, report) = append_trajectory(path, topology(), keep, options(format))
        .unwrap_or_else(|error| panic!("{format} {keep:?}: {error}"));
    write_frames(&mut writer, (from..count).map(frame));
    drop(writer.finish().unwrap());
    report
}

#[test]
fn a_run_continued_after_a_cut_at_any_byte_is_the_run_written_straight_through() {
    const FRAMES: u64 = 3;
    for format in FORMATS {
        let directory = Directory::new();
        let path = directory.file(format);
        let straight = straight_through(format, &path, FRAMES);
        let (first, ends) = frame_ends(&path);
        for cut in 0..=straight.len() {
            let context = format!("{format} cut at {cut}");
            let complete = ends.iter().filter(|&&end| end <= cut as u64).count();
            // Where the frames a restart keeps end: a DCD keeps its header
            // once it is whole; a file without one starts again.
            let kept_end = |kept: usize| match kept {
                0 if format == TrajectoryFormat::Dcd && cut as u64 >= first => first,
                0 => 0,
                kept => ends[kept - 1],
            };

            // Keep every complete frame: by step where frames record one.
            fs::write(&path, &straight[..cut]).unwrap();
            let keep = match (format, complete) {
                (TrajectoryFormat::Xyz, _) => KeepFrames::First(complete),
                (_, 0) => KeepFrames::ThroughStep(0),
                _ => KeepFrames::ThroughStep((complete as u64 - 1) * INTERVAL + INTERVAL / 2),
            };
            let report = continue_run(format, &path, keep, complete as u64, FRAMES);
            assert_eq!(report.kept_frames(), complete, "{context}");
            assert_eq!(report.dropped_frames(), 0, "{context}");
            assert_eq!(
                report.dropped_bytes(),
                cut as u64 - kept_end(complete),
                "{context}"
            );
            assert!(fs::read(&path).unwrap() == straight, "{context}");

            // Keep one complete frame fewer, by count, and write it again.
            if complete > 0 {
                fs::write(&path, &straight[..cut]).unwrap();
                let report = continue_run(
                    format,
                    &path,
                    KeepFrames::First(complete - 1),
                    complete as u64 - 1,
                    FRAMES,
                );
                assert_eq!(report.dropped_frames(), 1, "{context}");
                assert_eq!(
                    report.dropped_bytes(),
                    cut as u64 - kept_end(complete - 1),
                    "{context}"
                );
                assert!(fs::read(&path).unwrap() == straight, "{context}");
            }
        }
    }
}

#[test]
fn a_run_restarted_from_checkpoints_again_and_again_writes_one_trajectory() {
    for format in FORMATS {
        let directory = Directory::new();
        let path = directory.file(format);
        let straight = straight_through(format, &path, 8);
        let (_, ends) = frame_ends(&path);

        // A first run writes five frames, the last of them half, and dies.
        straight_through(format, &path, 5);
        let partial = usize::try_from(ends[3]).unwrap() + 7;
        let bytes = fs::read(&path).unwrap();
        fs::write(&path, &bytes[..partial]).unwrap();
        // Restarted at the checkpoint of frame 2 (step 20); it writes two
        // frames and dies.
        let keep = |frames: usize| match format {
            TrajectoryFormat::Xyz => KeepFrames::First(frames),
            _ => KeepFrames::ThroughStep((frames as u64 - 1) * INTERVAL),
        };
        let (mut writer, report) =
            append_trajectory(&path, topology(), keep(3), options(format)).unwrap();
        assert_eq!(
            (report.kept_frames(), report.dropped_frames()),
            (3, 1),
            "{format}"
        );
        assert_eq!(report.dropped_bytes(), partial as u64 - ends[2], "{format}");
        let expected_step = (format != TrajectoryFormat::Xyz).then_some(20);
        assert_eq!(report.last_kept_step(), expected_step, "{format}");
        write_frames(&mut writer, (3..5).map(frame));
        drop(writer);
        // Restarted at frame 4 (step 40); it dies right after the cut.
        drop(append_trajectory(&path, topology(), keep(5), options(format)).unwrap());
        // And once more from there, to the end.
        let report = continue_run(format, &path, keep(5), 5, 8);
        assert_eq!(report.kept_frames(), 5, "{format}");
        assert_eq!(report.dropped_bytes(), 0, "{format}");
        assert!(fs::read(&path).unwrap() == straight, "{format}");
    }
}

#[test]
fn a_non_periodic_xtc_run_continues_under_the_writer_cell_policy() {
    let format = TrajectoryFormat::Xtc;
    let zero_boxes = options(format).with_xtc_options(
        XtcWriteOptions::default().with_cell_policy(XtcCellPolicy::ZeroMatrixAsAbsent),
    );
    let open_frame = |index: u64| {
        let mut frame = frame(index);
        frame.set_cell(None);
        frame
    };
    let directory = Directory::new();
    let path = directory.file(format);
    let write = |count: u64| {
        let file = File::create(&path).unwrap();
        let mut writer = create_trajectory_stream_writer(
            BufWriter::new(file),
            topology(),
            zero_boxes.clone(),
            "run",
        )
        .unwrap();
        write_frames(&mut writer, (0..count).map(open_frame));
        drop(writer.finish().unwrap());
        fs::read(&path).unwrap()
    };
    let straight = write(4);
    let zero_box_reads = TrajectoryOpenOptions::default().with_xtc_options(
        XtcReadOptions::default().with_cell_policy(XtcCellPolicy::ZeroMatrixAsAbsent),
    );
    let (_, ends) = frame_ends_read(&path, zero_box_reads);
    // Cut inside the third frame, as a run killed while writing it.
    let cut = usize::try_from(ends[1]).unwrap() + 3;
    fs::write(&path, &straight[..cut]).unwrap();
    let (mut writer, report) = append_trajectory(
        &path,
        topology(),
        KeepFrames::ThroughStep(INTERVAL),
        zero_boxes.clone(),
    )
    .unwrap();
    assert_eq!((report.kept_frames(), report.dropped_frames()), (2, 0));
    write_frames(&mut writer, (2..4).map(open_frame));
    drop(writer.finish().unwrap());
    assert!(fs::read(&path).unwrap() == straight);
    // Options that require a periodic cell do not continue zero boxes.
    assert_refused(
        &path,
        KeepFrames::First(2),
        options(format),
        TrajectoryIoErrorKind::UnsupportedVariant,
        "XTC zero boxes",
    );
}

#[test]
fn the_scan_of_the_file_to_continue_can_be_cancelled_and_reports_progress() {
    let format = TrajectoryFormat::Trr;
    let directory = Directory::new();
    let path = directory.file(format);
    let straight = straight_through(format, &path, 4);
    // A cancelled scan returns before the file changes, as it does for a
    // file without a complete frame or DCD header to keep.
    let cancel = Arc::new(AtomicBool::new(true));
    let cancelled = TrajectoryAppendOptions::new(options(format))
        .with_cancellation(Cancellation::from_flag(Arc::clone(&cancel)));
    assert_refused(
        &path,
        KeepFrames::First(2),
        cancelled.clone(),
        TrajectoryIoErrorKind::Cancelled,
        "cancelled",
    );
    for unfinished in [TrajectoryFormat::Dcd, TrajectoryFormat::Xtc] {
        let other = directory.file(unfinished);
        let bytes = straight_through(unfinished, &other, 1);
        fs::write(&other, &bytes[..40]).unwrap();
        let options = TrajectoryAppendOptions::new(options(unfinished))
            .with_cancellation(Cancellation::from_flag(Arc::clone(&cancel)));
        assert_refused(
            &other,
            KeepFrames::First(0),
            options,
            TrajectoryIoErrorKind::Cancelled,
            &format!("cancelled unfinished {unfinished}"),
        );
    }
    // Progress reaches every frame of the file, kept or not.
    cancel.store(false, Ordering::Relaxed);
    let frames = Arc::new(AtomicUsize::new(0));
    let seen = Arc::clone(&frames);
    let followed = cancelled.with_progress(move |step| {
        seen.store(step.frame_count(), Ordering::Relaxed);
    });
    let (mut writer, report) =
        append_trajectory(&path, topology(), KeepFrames::First(2), followed).unwrap();
    assert_eq!(frames.load(Ordering::Relaxed), 4);
    assert_eq!(report.dropped_frames(), 2);
    write_frames(&mut writer, (2..4).map(frame));
    drop(writer.finish().unwrap());
    assert!(fs::read(&path).unwrap() == straight);
}

#[test]
fn a_reader_following_the_file_sees_the_cut_and_then_the_appended_frames() {
    let growing = || TrajectoryOpenOptions::default().with_growing(true);
    for format in FORMATS {
        let directory = Directory::new();
        let path = directory.file(format);
        straight_through(format, &path, 5);
        let mut ahead = open_indexed_trajectory(&path, topology(), growing()).unwrap();
        assert_eq!(ahead.frame_count(), Some(5));

        // Restarted after frame 2; the frames it writes again differ.
        let (mut writer, _) =
            append_trajectory(&path, topology(), KeepFrames::First(2), options(format)).unwrap();
        let error = ahead.refresh().unwrap_err();
        assert_eq!(
            error.kind(),
            TrajectoryIoErrorKind::IndexMismatch,
            "{format}"
        );
        // A reader opened after the cut, as the first one opens again.
        let mut reader = open_indexed_trajectory(&path, topology(), growing()).unwrap();
        assert_eq!(reader.frame_count(), Some(2), "{format}");
        write_frames(&mut writer, (2..6).map(|index| frame_with(index, 0.5)));
        assert_eq!(reader.refresh().unwrap(), 4, "{format}");
        drop(writer.finish().unwrap());
        assert_eq!(reader.refresh().unwrap(), 0, "{format}");

        let mut buffer = reader.frame_buffer();
        let expected = Trajectory::from_items(
            topology(),
            (0..6).map(|index| frame_with(index, if index < 2 { 0.0 } else { 0.5 })),
        )
        .unwrap();
        for index in 0..6 {
            reader.read_frame(index, &mut buffer).unwrap();
            let read = buffer.frame_view().positions().values().value()[3];
            let written = expected.get(index).unwrap().positions().values().value()[3];
            assert!(
                (read - written).norm() < 2e-3,
                "{format} frame {index}: {read:?} != {written:?}"
            );
        }
        assert!(reader
            .to_frame_index()
            .fingerprint()
            .same_content(&TrajectoryFileFingerprint::of_path(&path).unwrap()));
    }
}

/// Asserts that appending fails with `kind` and leaves the file unchanged.
fn assert_refused(
    path: &Path,
    keep: KeepFrames,
    options: impl Into<TrajectoryAppendOptions>,
    kind: TrajectoryIoErrorKind,
    context: &str,
) {
    let before = fs::read(path).unwrap();
    let error = append_trajectory(path, topology(), keep, options).unwrap_err();
    assert_eq!(error.kind(), kind, "{context}: {error}");
    assert!(fs::read(path).unwrap() == before, "{context}: file changed");
}

#[test]
fn appending_with_settings_other_than_the_file_records_is_refused() {
    let directory = Directory::new();
    let keep = KeepFrames::First(2);
    let mismatch = TrajectoryIoErrorKind::InconsistentMetadata;

    let path = directory.file(TrajectoryFormat::Xtc);
    straight_through(TrajectoryFormat::Xtc, &path, 3);
    let precision = XtcWriteOptions::default()
        .with_precision(Quantity::new(100.0, PER_NANOMETER))
        .unwrap();
    let other = options(TrajectoryFormat::Xtc).with_xtc_options(precision);
    assert_refused(&path, keep, other, mismatch, "XTC precision");
    let magic = XtcWriteOptions::default().with_magic(kekule_traj::io::xtc::XtcMagic::Xtc2023);
    let other = options(TrajectoryFormat::Xtc).with_xtc_options(magic);
    assert_refused(&path, keep, other, mismatch, "XTC magic");
    // Dropping every frame drops their settings too.
    let other = options(TrajectoryFormat::Xtc).with_xtc_options(precision);
    let (writer, report) =
        append_trajectory(&path, topology(), KeepFrames::First(0), other).unwrap();
    assert_eq!(report.dropped_frames(), 3);
    drop(writer);
    assert!(fs::read(&path).unwrap().is_empty());

    let path = directory.file(TrajectoryFormat::Trr);
    straight_through(TrajectoryFormat::Trr, &path, 3);
    let double = TrrWriteOptions::default()
        .with_lambda_policy(TrrLambdaPolicy::RequireZero)
        .with_precision(TrrScalarPrecision::Float64);
    let other = options(TrajectoryFormat::Trr).with_trr_options(double);
    assert_refused(&path, keep, other, mismatch, "TRR precision");
    // Lambdas the writer would not have written.
    let lambdas = TrajectoryWriteOptions::new(TrajectoryFormat::Trr);
    let lambda = |index: u64| {
        let mut frame = frame(index);
        frame
            .properties_mut()
            .owner_mut()
            .insert(
                PropertyKey::new("gromacs.trr.lambda").unwrap(),
                PropertyValue::real(0.25, DIMENSIONLESS).unwrap(),
            )
            .unwrap();
        frame
    };
    let file = File::create(&path).unwrap();
    let mut writer =
        create_trajectory_stream_writer(BufWriter::new(file), topology(), lambdas.clone(), "run")
            .unwrap();
    write_frames(&mut writer, (0..3).map(lambda));
    drop(writer.finish().unwrap());
    assert_refused(
        &path,
        keep,
        options(TrajectoryFormat::Trr),
        mismatch,
        "TRR lambda",
    );
    let (writer, report) = append_trajectory(&path, topology(), keep, lambdas).unwrap();
    assert_eq!(report.dropped_frames(), 1);
    drop(writer);

    let path = directory.file(TrajectoryFormat::Dcd);
    straight_through(TrajectoryFormat::Dcd, &path, 3);
    let dcd = DcdWriteOptions::default()
        .with_cell_policy(DcdCellPolicy::Write)
        .with_frame_metadata(UnrepresentablePolicy::Omit);
    for (other, setting) in [
        (dcd.with_step_sequence(0, 5), "step interval"),
        (dcd.with_step_sequence(10, INTERVAL), "first step"),
        (
            dcd.with_step_sequence(0, INTERVAL)
                .with_cell_policy(DcdCellPolicy::Omit),
            "unit-cell flag",
        ),
    ] {
        let other = options(TrajectoryFormat::Dcd).with_dcd_options(other);
        assert_refused(&path, keep, other, mismatch, setting);
        // A DCD keeps its header with no frame kept.
        assert_refused(
            &path,
            KeepFrames::First(0),
            options(TrajectoryFormat::Dcd).with_dcd_options(dcd.with_step_sequence(0, 5)),
            mismatch,
            setting,
        );
    }
}

#[test]
fn appending_refuses_what_it_cannot_continue_and_leaves_the_file() {
    let directory = Directory::new();
    for format in FORMATS {
        let path = directory.file(format);
        straight_through(format, &path, 3);
        let context = format.to_string();
        // More frames than the file holds.
        assert_refused(
            &path,
            KeepFrames::First(4),
            options(format),
            TrajectoryIoErrorKind::InconsistentMetadata,
            &context,
        );
        // Another topology.
        let other = linear_carbon_topology(11);
        let before = fs::read(&path).unwrap();
        assert!(
            append_trajectory(&path, other, KeepFrames::First(1), options(format)).is_err(),
            "{context}"
        );
        assert!(fs::read(&path).unwrap() == before, "{context}");
        // Corrupt data, even after the kept frames.
        let (_, ends) = frame_ends(&path);
        let mut corrupt = before.clone();
        let inside_last = usize::try_from(ends[1]).unwrap() + 2;
        corrupt[inside_last] = match format {
            TrajectoryFormat::Xyz => b'x',
            _ => corrupt[inside_last] ^ 0x40,
        };
        fs::write(&path, &corrupt).unwrap();
        let error = append_trajectory(&path, topology(), KeepFrames::First(1), options(format));
        assert!(error.is_err(), "{context}");
        assert!(fs::read(&path).unwrap() == corrupt, "{context}");
    }
    // XYZ frames record no step.
    let path = directory.file(TrajectoryFormat::Xyz);
    straight_through(TrajectoryFormat::Xyz, &path, 2);
    assert_refused(
        &path,
        KeepFrames::ThroughStep(10),
        options(TrajectoryFormat::Xyz),
        TrajectoryIoErrorKind::UnsupportedVariant,
        "XYZ steps",
    );
    // A missing file is no trajectory to continue.
    let missing = directory.0.join("missing.xtc");
    let error = append_trajectory(
        &missing,
        topology(),
        KeepFrames::First(0),
        options(TrajectoryFormat::Xtc),
    )
    .unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::Io);
    assert!(!missing.exists());
    // Invalid options fail before an unfinished header is replaced.
    let path = directory.file(TrajectoryFormat::Dcd);
    fs::write(
        &path,
        &straight_through(TrajectoryFormat::Dcd, &path, 1)[..100],
    )
    .unwrap();
    let invalid = options(TrajectoryFormat::Dcd)
        .with_dcd_options(DcdWriteOptions::default().with_step_sequence(0, 0));
    assert_refused(
        &path,
        KeepFrames::First(0),
        invalid,
        TrajectoryIoErrorKind::InconsistentMetadata,
        "DCD step interval 0",
    );
}
