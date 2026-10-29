//! Following trajectory files that a writer is still appending to: indexed
//! readers pick up appended frames without scanning indexed ones again, take
//! a partial final frame as not yet written, and still reject corrupt data,
//! and writers leave whole frames after every flush.

use std::fs::{self, File, OpenOptions};
use std::io::{BufWriter, Cursor, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{mpsc, Arc, OnceLock};
use std::thread;
use std::time::{Duration, Instant};

use kekule::diagnostics::UnrepresentablePolicy;
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
use kekule::topology::Topology;
use kekule::units::{Quantity, NANOMETER, PICOSECOND};
use kekule_traj::io::dcd::{DcdCellPolicy, DcdWriteOptions};
use kekule_traj::io::trr::{TrrLambdaPolicy, TrrWriteOptions};
use kekule_traj::io::xyz::XyzWriteOptions;
use kekule_traj::io::{
    create_trajectory_stream_writer, open_indexed_trajectory, open_indexed_trajectory_from_index,
    IndexedFileTrajectoryReader, TrajectoryFileFingerprint, TrajectoryOpenOptions,
    TrajectoryWriteOptions,
};
use kekule_traj::{
    SeekableTrajectoryReader, TrajectoryFormat, TrajectoryIoErrorKind, TrajectoryReader,
    TrajectoryWriter,
};

use crate::support::{self, linear_carbon_topology};

const FORMATS: [TrajectoryFormat; 4] = [
    TrajectoryFormat::Xtc,
    TrajectoryFormat::Trr,
    TrajectoryFormat::Dcd,
    TrajectoryFormat::Xyz,
];

/// A scratch directory removed on drop.
struct Directory(PathBuf);

impl Directory {
    fn new() -> Self {
        static SEQUENCE: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "kekule-growing-{}-{}",
            std::process::id(),
            SEQUENCE.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }

    /// A path for a trajectory of `format`, named with its extension.
    fn file(&self, format: TrajectoryFormat) -> PathBuf {
        self.0
            .join(format!("live.{}", format.name().to_ascii_lowercase()))
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

/// Frame `step` of a moving chain, with the time, step, and cell XTC needs.
fn frame(step: u64) -> TrajectoryFrame {
    let shift = step as f64;
    let points = (0..12)
        .map(|atom| {
            let atom = atom as f64;
            Point3::new(
                0.15 * atom + 0.011 * shift,
                0.2 + 0.007 * shift * atom.sin(),
                0.3 - 0.005 * shift,
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
    frame.set_step(Some(step));
    frame
        .set_time(Some(Quantity::new(0.5 * shift, PICOSECOND)))
        .unwrap();
    frame
}

fn trajectory(steps: std::ops::Range<u64>) -> Trajectory {
    Trajectory::from_items(topology(), steps.map(frame)).unwrap()
}

/// Frame `step` with every atom at one point, which XTC and XYZ write in
/// fewer bytes than [`frame`].
fn squeezed(step: u64) -> TrajectoryFrame {
    let mut squeezed = frame(step);
    squeezed
        .set_positions(Quantity::new([Point3::new(0.1, 0.1, 0.1); 12], NANOMETER))
        .unwrap();
    squeezed
}

/// Write options keeping what each format holds of [`frame`].
fn options(format: TrajectoryFormat) -> TrajectoryWriteOptions {
    TrajectoryWriteOptions::new(format)
        .with_dcd_options(
            DcdWriteOptions::default()
                .with_cell_policy(DcdCellPolicy::Write)
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

fn growing() -> TrajectoryOpenOptions {
    TrajectoryOpenOptions::default().with_growing(true)
}

/// The complete encoding of `trajectory` in `format`.
fn encoded(format: TrajectoryFormat, trajectory: &Trajectory) -> Vec<u8> {
    let mut stream = Cursor::new(Vec::new());
    kekule_traj::io::write_trajectory_to(&mut stream, trajectory, options(format)).unwrap();
    stream.into_inner()
}

/// Every indexed frame, read at random.
fn frames(reader: &mut IndexedFileTrajectoryReader) -> Vec<TrajectoryFrame> {
    let mut buffer = reader.frame_buffer();
    (0..reader.frame_count().unwrap())
        .map(|index| {
            reader.read_frame(index, &mut buffer).unwrap();
            buffer.frame_view().payload().clone()
        })
        .collect()
}

/// The frames a strict reader finds in the complete file at `path`.
fn scanned(path: &Path) -> Vec<TrajectoryFrame> {
    frames(&mut open_indexed_trajectory(path, topology(), Default::default()).unwrap())
}

fn append(path: &Path, bytes: &[u8]) {
    let mut file = OpenOptions::new().append(true).open(path).unwrap();
    file.write_all(bytes).unwrap();
}

/// The x coordinates of frame `frame`, to the precision every format keeps.
fn assert_close(read: &TrajectoryFrame, step: u64, context: &str) {
    let expected = frame(step);
    let pairs = read
        .positions()
        .values()
        .value()
        .iter()
        .zip(expected.positions().values().value().iter());
    for (read, written) in pairs {
        assert!(
            (read.x - written.x).abs() < 2e-3 && (read.z - written.z).abs() < 2e-3,
            "{context}: frame {step} {read:?} != {written:?}"
        );
    }
}

#[test]
fn every_format_follows_a_file_its_writer_appends_to() {
    for format in FORMATS {
        let directory = Directory::new();
        let path = directory.file(format);
        let written = trajectory(0..6);
        let file = File::create(&path).unwrap();
        let mut writer = create_trajectory_stream_writer(
            BufWriter::new(file),
            topology(),
            options(format),
            "live",
        )
        .unwrap();
        for frame in written.iter().take(2) {
            writer.write_frame(frame).unwrap();
        }
        writer.flush().unwrap();

        let mut reader = open_indexed_trajectory(&path, topology(), growing()).unwrap();
        assert_eq!(reader.frame_count(), Some(2), "{format}");
        assert_eq!(reader.pending_byte_count(), 0, "{format}");
        // Sequential reads end with the indexed frames.
        let mut buffer = reader.frame_buffer();
        assert!(reader.read_next(&mut buffer).unwrap());
        assert!(reader.read_next(&mut buffer).unwrap());
        assert!(!reader.read_next(&mut buffer).unwrap(), "{format}");
        assert_eq!(reader.refresh().unwrap(), 0, "{format}");

        for frame in written.iter().skip(2) {
            writer.write_frame(frame).unwrap();
        }
        writer.flush().unwrap();
        assert_eq!(reader.refresh().unwrap(), 4, "{format}");
        assert_eq!(reader.frame_count(), Some(6));
        assert_eq!(reader.metadata().indexed_frame_count(), Some(6));
        // Sequential reads continue into the new frames.
        let mut sequential = 2;
        while reader.read_next(&mut buffer).unwrap() {
            assert_close(
                buffer.frame_view().payload(),
                sequential,
                &format.to_string(),
            );
            sequential += 1;
        }
        assert_eq!(sequential, 6, "{format}");
        drop(writer.finish().unwrap());
        // Finishing a DCD rewrites only the header's frame count, which is
        // no content: the indexed bytes still match.
        assert_eq!(reader.refresh().unwrap(), 0, "{format}");
        let read = frames(&mut reader);
        assert_eq!(read, scanned(&path), "{format}");
        for (step, frame) in read.iter().enumerate() {
            assert_close(frame, step as u64, &format.to_string());
        }

        // The exported index covers the whole file, as a fresh scan's does;
        // only the modification time can differ, after a DCD's count.
        let index = reader.to_frame_index();
        assert!(
            index
                .fingerprint()
                .same_content(&TrajectoryFileFingerprint::of_path(&path).unwrap()),
            "{format}"
        );
        let mut reopened =
            open_indexed_trajectory_from_index(&path, topology(), index, Default::default())
                .unwrap();
        assert_eq!(frames(&mut reopened), read, "{format}");
    }
}

/// The start of every frame and the end of the last in `bytes`.
fn frame_bounds(format: TrajectoryFormat, bytes: &[u8]) -> (u64, Vec<u64>) {
    let directory = Directory::new();
    let path = directory.file(format);
    fs::write(&path, bytes).unwrap();
    let index = open_indexed_trajectory(&path, topology(), Default::default())
        .unwrap()
        .to_frame_index();
    let starts = index.frame_offsets().to_vec();
    (starts[0], starts)
}

#[test]
fn a_partial_final_frame_is_pending_at_every_byte_and_complete_once_written() {
    for format in FORMATS {
        let bytes = encoded(format, &trajectory(0..3));
        let (first, mut bounds) = frame_bounds(format, &bytes);
        bounds.push(bytes.len() as u64);
        let directory = Directory::new();
        let path = directory.file(format);
        let complete = {
            fs::write(&path, &bytes).unwrap();
            scanned(&path)
        };
        for cut in 0..bytes.len() {
            let context = format!("{format} cut at {cut}");
            fs::write(&path, &bytes[..cut]).unwrap();
            let opened = open_indexed_trajectory(&path, topology(), growing());
            // Frames that end at or before the cut are complete.
            let whole = bounds[1..].iter().filter(|&&end| end <= cut as u64).count();
            let indexed_end = if whole == 0 { first } else { bounds[whole] };
            let header_complete = format == TrajectoryFormat::Dcd && cut as u64 >= first;
            if whole == 0 && !header_complete {
                let error = opened.unwrap_err();
                assert_eq!(error.kind(), TrajectoryIoErrorKind::Incomplete, "{context}");
                continue;
            }
            let mut reader = opened.unwrap_or_else(|error| panic!("{context}: {error}"));
            assert_eq!(reader.frame_count(), Some(whole), "{context}");
            assert_eq!(
                reader.pending_byte_count(),
                cut as u64 - indexed_end,
                "{context}"
            );
            assert_eq!(frames(&mut reader), complete[..whole], "{context}");
            // A strict reader takes the same bytes as corrupt, unless they
            // end between frames of a format without a header count.
            let strict = open_indexed_trajectory(&path, topology(), Default::default());
            let between_frames = cut as u64 == indexed_end && format != TrajectoryFormat::Dcd;
            if format != TrajectoryFormat::Xyz {
                assert_eq!(strict.is_ok(), between_frames, "{context}");
            }

            // Once the rest is written, refreshing indexes the frames from
            // the cut without reading the indexed ones again.
            append(&path, &bytes[cut..]);
            assert_eq!(reader.refresh().unwrap(), 3 - whole, "{context}");
            assert_eq!(reader.pending_byte_count(), 0, "{context}");
            assert_eq!(frames(&mut reader), complete, "{context}");
        }
    }
}

#[test]
fn corrupt_frames_of_a_growing_file_are_errors_and_refresh_changes_nothing_on_them() {
    for format in FORMATS {
        let bytes = encoded(format, &trajectory(0..3));
        let (_, starts) = frame_bounds(format, &bytes);
        let second = usize::try_from(starts[1]).unwrap();
        let third = usize::try_from(starts[2]).unwrap();
        // A first byte no valid frame starts with, whose frame is otherwise
        // complete: a magic, a record marker that disagrees with its
        // trailing marker, or an atom-count line that is not a number.
        let mut corrupt = bytes[second..third].to_vec();
        match format {
            TrajectoryFormat::Dcd => corrupt[0] += 1,
            TrajectoryFormat::Xyz => corrupt[0] = b'x',
            _ => corrupt[0] = 0xff,
        }
        let directory = Directory::new();
        let path = directory.file(format);

        // Inside a file, opening fails however the file is read.
        let mut inside = bytes[..second].to_vec();
        inside.extend_from_slice(&corrupt);
        inside.extend_from_slice(&bytes[third..]);
        fs::write(&path, &inside).unwrap();
        let error = open_indexed_trajectory(&path, topology(), growing()).unwrap_err();
        assert!(
            !matches!(
                error.kind(),
                TrajectoryIoErrorKind::TruncatedRecord | TrajectoryIoErrorKind::Incomplete
            ),
            "{format}: {error}"
        );

        // Appended, a refresh fails and leaves the reader as it was.
        fs::write(&path, &bytes[..second]).unwrap();
        let mut reader = open_indexed_trajectory(&path, topology(), growing()).unwrap();
        let before = reader.to_frame_index();
        append(&path, &corrupt);
        let error = reader.refresh().unwrap_err();
        assert!(
            error.kind() != TrajectoryIoErrorKind::TruncatedRecord,
            "{format}: {error}"
        );
        assert_eq!(reader.to_frame_index(), before, "{format}");
        assert_eq!(reader.frame_count(), Some(1));

        // Replacing the corrupt bytes lets the next refresh succeed.
        fs::write(&path, &bytes).unwrap();
        assert_eq!(reader.refresh().unwrap(), 2, "{format}");

        // A frame the file ends inside is rejected once what it holds is
        // wrong, such as an XTC or TRR magic.
        if matches!(format, TrajectoryFormat::Xtc | TrajectoryFormat::Trr) {
            fs::write(&path, &bytes[..second]).unwrap();
            let mut reader = open_indexed_trajectory(&path, topology(), growing()).unwrap();
            append(&path, &corrupt[..6]);
            let error = reader.refresh().unwrap_err();
            assert_eq!(
                error.kind(),
                TrajectoryIoErrorKind::InvalidHeader,
                "{format}"
            );
            append(&path, &corrupt[6..]);
            let error = open_indexed_trajectory(&path, topology(), growing()).unwrap_err();
            assert_eq!(
                error.kind(),
                TrajectoryIoErrorKind::InvalidHeader,
                "{format}"
            );
        }
    }
}

#[test]
fn refresh_refuses_a_file_whose_indexed_bytes_changed() {
    for format in FORMATS {
        let bytes = encoded(format, &trajectory(0..3));
        let (_, starts) = frame_bounds(format, &bytes);
        let directory = Directory::new();
        let path = directory.file(format);
        let second = usize::try_from(starts[1]).unwrap();
        fs::write(&path, &bytes[..second]).unwrap();
        let mut reader = open_indexed_trajectory(&path, topology(), growing()).unwrap();
        let before = reader.to_frame_index();
        // The same bytes and the rest of the trajectory, one indexed byte
        // changed: not the indexed file grown.
        let mut changed = bytes.clone();
        changed[second - 5] ^= 0x01;
        fs::write(&path, &changed).unwrap();
        let error = reader.refresh().unwrap_err();
        assert_eq!(
            error.kind(),
            TrajectoryIoErrorKind::IndexMismatch,
            "{format}"
        );
        assert_eq!(reader.to_frame_index(), before);
        let error =
            open_indexed_trajectory_from_index(&path, topology(), before.clone(), growing())
                .unwrap_err();
        assert_eq!(
            error.kind(),
            TrajectoryIoErrorKind::IndexMismatch,
            "{format}"
        );
        // Nor is a file truncated below the indexed bytes.
        fs::write(&path, &bytes[..second - 1]).unwrap();
        let error = reader.refresh().unwrap_err();
        assert_eq!(
            error.kind(),
            TrajectoryIoErrorKind::IndexMismatch,
            "{format}"
        );
        // The grown file is.
        fs::write(&path, &bytes).unwrap();
        assert_eq!(reader.refresh().unwrap(), 2, "{format}");
    }
}

#[test]
fn a_pending_frame_is_read_again_though_the_file_kept_its_length() {
    // The bytes a frame was pending in may hold a whole frame by the next
    // refresh without the file growing, as when the scan that found the
    // frame partial measured the file after the writer completed it.
    for format in [TrajectoryFormat::Xtc, TrajectoryFormat::Xyz] {
        let bytes = encoded(format, &trajectory(0..3));
        let shorter = encoded(
            format,
            &Trajectory::from_items(topology(), [frame(0), frame(1), squeezed(2)]).unwrap(),
        );
        let (_, starts) = frame_bounds(format, &bytes);
        let third = usize::try_from(starts[2]).unwrap();
        assert_eq!(bytes[..third], shorter[..third], "{format}");
        assert!(shorter.len() < bytes.len(), "{format}");
        let directory = Directory::new();
        let path = directory.file(format);
        fs::write(&path, &bytes[..shorter.len()]).unwrap();
        let mut reader = open_indexed_trajectory(&path, topology(), growing()).unwrap();
        assert_eq!(reader.frame_count(), Some(2), "{format}");
        assert_eq!(
            reader.pending_byte_count(),
            (shorter.len() - third) as u64,
            "{format}"
        );
        fs::write(&path, &shorter).unwrap();
        assert_eq!(reader.refresh().unwrap(), 1, "{format}");
        assert_eq!(reader.pending_byte_count(), 0, "{format}");
        assert_eq!(frames(&mut reader), scanned(&path), "{format}");
    }
}

#[test]
fn refresh_reads_a_file_moved_into_the_path_when_it_holds_the_indexed_bytes() {
    for format in FORMATS {
        let directory = Directory::new();
        let path = directory.file(format);
        let moved = directory.0.join("moved");
        fs::write(&path, encoded(format, &trajectory(0..2))).unwrap();
        let mut reader = open_indexed_trajectory(&path, topology(), growing()).unwrap();
        assert_eq!(reader.frame_count(), Some(2), "{format}");

        // A longer copy of the trajectory replaces the file.
        fs::write(&moved, encoded(format, &trajectory(0..5))).unwrap();
        fs::rename(&moved, &path).unwrap();
        assert_eq!(reader.refresh().unwrap(), 3, "{format}");
        let read = frames(&mut reader);
        assert_eq!(read, scanned(&path), "{format}");
        assert!(
            reader
                .to_frame_index()
                .fingerprint()
                .same_content(&TrajectoryFileFingerprint::of_path(&path).unwrap()),
            "{format}"
        );

        // Another trajectory does not, and the reader keeps reading the
        // file it had.
        let index = reader.to_frame_index();
        let other = Trajectory::from_items(
            topology(),
            [squeezed(0)].into_iter().chain((1..6).map(frame)),
        )
        .unwrap();
        fs::write(&moved, encoded(format, &other)).unwrap();
        fs::rename(&moved, &path).unwrap();
        let error = reader.refresh().unwrap_err();
        assert_eq!(
            error.kind(),
            TrajectoryIoErrorKind::IndexMismatch,
            "{format}"
        );
        assert_eq!(reader.to_frame_index(), index, "{format}");
        assert_eq!(frames(&mut reader), read, "{format}");
    }
}

#[test]
fn a_dcd_record_marker_of_another_size_is_corrupt_where_the_file_ends() {
    let format = TrajectoryFormat::Dcd;
    let bytes = encoded(format, &trajectory(0..3));
    let (_, starts) = frame_bounds(format, &bytes);
    let second = usize::try_from(starts[1]).unwrap();
    let directory = Directory::new();
    let path = directory.file(format);
    // The second frame's first record marker, intact, is a frame still
    // being written.
    fs::write(&path, &bytes[..second + 4]).unwrap();
    let mut reader = open_indexed_trajectory(&path, topology(), growing()).unwrap();
    assert_eq!(reader.frame_count(), Some(1));
    assert_eq!(reader.pending_byte_count(), 4);
    // Claiming more bytes than a unit-cell record holds, it is corrupt,
    // though the file ends before them.
    let mut corrupt = bytes[..second + 4].to_vec();
    corrupt[second + 2] ^= 0x40;
    fs::write(&path, &corrupt).unwrap();
    let error = reader.refresh().unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidFrame);
    let error = open_indexed_trajectory(&path, topology(), growing()).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidFrame);
}

#[test]
fn files_without_a_complete_frame_or_header_are_incomplete() {
    for format in FORMATS {
        let directory = Directory::new();
        let path = directory.file(format);
        File::create(&path).unwrap();
        let error = open_indexed_trajectory(&path, topology(), growing()).unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::Incomplete, "{format}");
        // A complete file opened without the option is a different error.
        let error = open_indexed_trajectory(&path, topology(), Default::default()).unwrap_err();
        assert_ne!(error.kind(), TrajectoryIoErrorKind::Incomplete, "{format}");
    }
    // A DCD whose header is complete opens without frames.
    let directory = Directory::new();
    let path = directory.file(TrajectoryFormat::Dcd);
    let file = File::create(&path).unwrap();
    let mut writer = create_trajectory_stream_writer(
        BufWriter::new(file),
        topology(),
        options(TrajectoryFormat::Dcd),
        "live",
    )
    .unwrap();
    writer.flush().unwrap();
    let mut reader = open_indexed_trajectory(&path, topology(), growing()).unwrap();
    assert_eq!(reader.frame_count(), Some(0));
    assert_eq!(reader.metadata().declared_frame_count(), Some(0));
    writer
        .write_frame(trajectory(0..1).get(0).unwrap())
        .unwrap();
    writer.flush().unwrap();
    assert_eq!(reader.refresh().unwrap(), 1);
    assert_eq!(reader.metadata().declared_frame_count(), Some(1));
}

#[test]
fn growing_dcd_frames_are_counted_in_the_file_whatever_the_header_says() {
    // Writers such as OpenMM's count a frame in the header before writing
    // it, kekule's after; neither header is the frame count while the
    // file grows.
    let bytes = encoded(TrajectoryFormat::Dcd, &trajectory(0..3));
    let (_, starts) = frame_bounds(TrajectoryFormat::Dcd, &bytes);
    let directory = Directory::new();
    let path = directory.file(TrajectoryFormat::Dcd);
    let third = usize::try_from(starts[2]).unwrap();
    let with_count = |count: u32, end: usize| {
        let mut prefix = bytes[..end].to_vec();
        prefix[8..12].copy_from_slice(&count.to_le_bytes());
        prefix
    };
    for count in [0, 1, 2, 3, 7] {
        fs::write(&path, with_count(count, third + 10)).unwrap();
        let mut reader = open_indexed_trajectory(&path, topology(), growing()).unwrap();
        assert_eq!(reader.frame_count(), Some(2), "NSET {count}");
        assert_eq!(reader.pending_byte_count(), 10);
        assert_eq!(
            reader.metadata().declared_frame_count(),
            Some(count as usize)
        );
        fs::write(&path, with_count(3, bytes.len())).unwrap();
        assert_eq!(reader.refresh().unwrap(), 1, "NSET {count}");
        assert_eq!(reader.metadata().declared_frame_count(), Some(3));
    }
    // A strict reader holds the file to its header.
    fs::write(&path, with_count(2, bytes.len())).unwrap();
    let error = open_indexed_trajectory(&path, topology(), Default::default()).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::InconsistentMetadata);
}

#[test]
fn fingerprints_take_a_file_grown_by_appending_for_the_same_file() {
    for format in FORMATS {
        let bytes = encoded(format, &trajectory(0..3));
        let (_, starts) = frame_bounds(format, &bytes);
        let directory = Directory::new();
        let path = directory.file(format);
        let second = usize::try_from(starts[1]).unwrap();
        fs::write(&path, &bytes[..second]).unwrap();
        let recorded = TrajectoryFileFingerprint::of_path(&path).unwrap();
        assert!(recorded.is_prefix_of_path(&path).unwrap());
        append(&path, &bytes[second..]);
        assert!(recorded.is_prefix_of_path(&path).unwrap(), "{format}");
        let mut file = File::open(&path).unwrap();
        assert!(recorded.is_prefix_of_file(&mut file).unwrap());
        assert_eq!(file.stream_position().unwrap(), 0);
        let grown = TrajectoryFileFingerprint::of_path(&path).unwrap();
        assert!(!grown.same_content(&recorded));
        assert!(!grown.is_prefix_of_path(directory.0.join("missing")).is_ok());
        // A changed byte of the recorded content, or a shorter file, is not.
        let mut changed = bytes.clone();
        changed[second / 2] ^= 0x10;
        fs::write(&path, &changed).unwrap();
        assert!(!recorded.is_prefix_of_path(&path).unwrap(), "{format}");
        fs::write(&path, &bytes[..second - 1]).unwrap();
        assert!(!recorded.is_prefix_of_path(&path).unwrap(), "{format}");
    }
    // A DCD header's frame counts are not content; its other fields are.
    let bytes = encoded(TrajectoryFormat::Dcd, &trajectory(0..2));
    let directory = Directory::new();
    let path = directory.file(TrajectoryFormat::Dcd);
    fs::write(&path, &bytes).unwrap();
    let recorded = TrajectoryFileFingerprint::of_path(&path).unwrap();
    for (offset, same) in [(8, true), (20, true), (12, false), (16, false), (24, false)] {
        let mut changed = bytes.clone();
        changed[offset] ^= 0x04;
        fs::write(&path, &changed).unwrap();
        let fingerprint = TrajectoryFileFingerprint::of_path(&path).unwrap();
        assert_eq!(fingerprint.same_content(&recorded), same, "byte {offset}");
        assert_eq!(
            recorded.is_prefix_of_path(&path).unwrap(),
            same,
            "byte {offset}"
        );
    }
}

#[test]
fn every_flush_leaves_a_complete_trajectory_a_strict_reader_opens() {
    for format in FORMATS {
        let directory = Directory::new();
        let path = directory.file(format);
        let written = trajectory(0..4);
        let file = File::create(&path).unwrap();
        // A buffer smaller than a frame writes frames in pieces.
        let mut writer = create_trajectory_stream_writer(
            BufWriter::with_capacity(16, file),
            topology(),
            options(format),
            "live",
        )
        .unwrap();
        for (count, frame) in written.iter().enumerate() {
            writer.write_frame(frame).unwrap();
            writer.flush().unwrap();
            let read = scanned(&path);
            assert_eq!(read.len(), count + 1, "{format}");
            assert_close(&read[count], count as u64, &format.to_string());
        }
        let mut file = writer.finish().unwrap().into_inner().unwrap();
        assert_eq!(
            file.stream_position().unwrap(),
            file.seek(SeekFrom::End(0)).unwrap()
        );
    }
}

#[test]
fn a_reader_refreshing_while_a_writer_appends_sees_every_frame_whole() {
    const FRAMES: u64 = 40;
    for format in FORMATS {
        let directory = Directory::new();
        let path = directory.file(format);
        let (started, first_frame) = mpsc::channel();
        let target = path.clone();
        let writer = thread::spawn(move || {
            let file = File::create(&target).unwrap();
            // A small buffer reaches the file partway through frames, so
            // the reader meets partial frames.
            let mut writer = create_trajectory_stream_writer(
                BufWriter::with_capacity(37, file),
                topology(),
                options(format),
                "live",
            )
            .unwrap();
            let written = trajectory(0..FRAMES);
            for (index, frame) in written.iter().enumerate() {
                writer.write_frame(frame).unwrap();
                if index % 3 != 1 {
                    writer.flush().unwrap();
                }
                if index == 0 {
                    started.send(()).unwrap();
                }
                thread::sleep(Duration::from_micros(200));
            }
            writer.finish().unwrap();
        });
        first_frame.recv().unwrap();
        let mut reader = open_indexed_trajectory(&path, topology(), growing())
            .unwrap_or_else(|error| panic!("{format}: {error}"));
        let mut buffer = reader.frame_buffer();
        let mut seen = 0;
        let deadline = Instant::now() + Duration::from_secs(60);
        loop {
            let finished = writer.is_finished();
            reader
                .refresh()
                .unwrap_or_else(|error| panic!("{format}: {error}"));
            let count = reader.frame_count().unwrap();
            for index in seen..count {
                reader.read_frame(index, &mut buffer).unwrap();
                assert_close(
                    buffer.frame_view().payload(),
                    index as u64,
                    &format.to_string(),
                );
            }
            seen = count;
            if finished && count == FRAMES as usize {
                break;
            }
            assert!(Instant::now() < deadline, "{format}: {seen} frames");
        }
        writer.join().unwrap();
        assert_eq!(reader.pending_byte_count(), 0, "{format}");
        assert_eq!(frames(&mut reader), scanned(&path), "{format}");
        // The index exported while following the file reopens it.
        let mut reopened = open_indexed_trajectory_from_index(
            &path,
            topology(),
            reader.to_frame_index(),
            Default::default(),
        )
        .unwrap();
        assert_eq!(support::read_all(&mut reopened).len(), FRAMES as usize);
    }
}
