use std::fs::{self, File, OpenOptions};
use std::io::{Cursor, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use kekule::diagnostics::{Cancellation, Cancelled};
use kekule::geometry::Point3;
use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
use kekule::topology::Topology;
use kekule::units::{Quantity, FEMTOSECOND, NANOMETER, PICOSECOND};
use kekule_traj::io::dcd::{DcdReadOptions, DcdTimePolicy};
use kekule_traj::io::trr::{TrrLambdaPolicy, TrrScalarPrecision, TrrWriteOptions, TrrWriter};
use kekule_traj::io::{
    open_indexed_trajectory, open_indexed_trajectory_from_index, CoordinateEncoding,
    IndexedFileTrajectoryReader, ScalarPrecision, TrajectoryFrameIndex, TrajectoryIndexProgress,
    TrajectoryIoLimits, TrajectoryOpenOptions,
};
use kekule_traj::{
    SeekableTrajectoryReader, TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind,
    TrajectoryIoOperation, TrajectoryReader, TrajectoryWriter,
};

use crate::support;
use support::{linear_carbon_topology, topology};

/// A scratch directory removed on drop.
struct Directory(PathBuf);

impl Directory {
    fn new() -> Self {
        static SEQUENCE: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "kekule-frame-index-{}-{}",
            std::process::id(),
            SEQUENCE.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }

    /// A copy of the named fixture in this directory.
    fn fixture(&self, name: &str) -> PathBuf {
        let path = self.0.join(name);
        fs::copy(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/fixtures")
                .join(name),
            &path,
        )
        .unwrap();
        path
    }
}

impl Drop for Directory {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}

/// The externally supplied fixture of every format with its topology.
fn fixtures() -> [(&'static str, Arc<Topology>); 4] {
    [
        (
            "ase-3.26.0-water.xyz",
            topology(&["O", "H", "H"], &[(0, 1), (0, 2)]),
        ),
        (
            "mdanalysis-2.9.0-three-atoms.dcd",
            topology(&["C", "H", "O"], &[(0, 1), (0, 2)]),
        ),
        (
            "mdanalysis-2.9.0-three-atoms.trr",
            topology(&["C", "H", "O"], &[(0, 1), (0, 2)]),
        ),
        (
            "mdanalysis-2.9.0-twelve-atoms.xtc",
            linear_carbon_topology(12),
        ),
    ]
}

fn frames(reader: &mut IndexedFileTrajectoryReader) -> Vec<TrajectoryFrame> {
    let mut buffer = reader.frame_buffer();
    (0..reader.frame_count().unwrap())
        .map(|index| {
            reader.read_frame(index, &mut buffer).unwrap();
            buffer.frame_view().payload().clone()
        })
        .collect()
}

fn reopen(
    path: &Path,
    topology: &Arc<Topology>,
    bytes: &[u8],
) -> Result<IndexedFileTrajectoryReader, TrajectoryIoError> {
    let index = TrajectoryFrameIndex::from_bytes(bytes, &TrajectoryIoLimits::default())?;
    open_indexed_trajectory_from_index(path, Arc::clone(topology), index, Default::default())
}

/// The documented FNV-1a 64 checksum, so tests can forge valid encodings.
fn fnv1a(bytes: &[u8]) -> u64 {
    bytes.iter().fold(0xcbf2_9ce4_8422_2325, |hash, &byte| {
        (hash ^ u64::from(byte)).wrapping_mul(0x0000_0100_0000_01b3)
    })
}

/// Re-encodes `index` with `starts` as its frame offsets, the steps and times
/// of its first `starts.len()` frames, and a valid checksum.
fn forge(index: &TrajectoryFrameIndex, starts: &[u64]) -> Vec<u8> {
    let mut forged = index.to_bytes()[..48].to_vec();
    forged.extend_from_slice(&(starts.len() as u64).to_le_bytes());
    for start in starts {
        forged.extend_from_slice(&start.to_le_bytes());
    }
    for step in index
        .steps()
        .map(|steps| &steps[..starts.len()])
        .unwrap_or_default()
    {
        forged.extend_from_slice(&step.to_le_bytes());
    }
    let times = index.times();
    for time in times
        .as_ref()
        .map(|times| &times.value()[..starts.len()])
        .unwrap_or_default()
    {
        forged.extend_from_slice(&time.to_le_bytes());
    }
    let checksum = fnv1a(&forged);
    forged.extend_from_slice(&checksum.to_le_bytes());
    forged
}

/// Recomputes the checksum of an edited encoding.
fn reencode(bytes: &mut [u8]) {
    let body = bytes.len() - 8;
    let checksum = fnv1a(&bytes[..body]);
    bytes[body..].copy_from_slice(&checksum.to_le_bytes());
}

/// Checks that `index` holds exactly the steps and times, in picoseconds, of
/// `frames`.
fn assert_stamps(index: &TrajectoryFrameIndex, frames: &[TrajectoryFrame]) {
    let (steps, times): (Vec<_>, Vec<_>) = frames
        .iter()
        .map(|frame| {
            (
                frame.step(),
                frame.time().map(|time| time.value_in(PICOSECOND).unwrap()),
            )
        })
        .unzip();
    let indexed_steps = index.steps().map_or_else(
        || vec![None; frames.len()],
        |steps| steps.iter().copied().map(Some).collect(),
    );
    let indexed_times = index.times().map_or_else(
        || vec![None; frames.len()],
        |times| {
            assert_eq!(times.unit(), PICOSECOND);
            times.value().iter().copied().map(Some).collect()
        },
    );
    assert_eq!(indexed_steps, steps);
    assert_eq!(indexed_times, times);
}

fn kind<T>(result: Result<T, TrajectoryIoError>) -> TrajectoryIoErrorKind {
    match result {
        Ok(_) => panic!("expected an error"),
        Err(error) => error.kind(),
    }
}

#[test]
fn exported_indexes_reopen_every_format_with_identical_frames_and_metadata() {
    let directory = Directory::new();
    for (name, topology) in fixtures() {
        let path = directory.fixture(name);
        let size = fs::metadata(&path).unwrap().len();
        let steps = Arc::new(Mutex::new(Vec::<TrajectoryIndexProgress>::new()));
        let recorded = Arc::clone(&steps);
        let options = TrajectoryOpenOptions::default()
            .with_index_progress(move |step| recorded.lock().unwrap().push(step));
        let mut built = open_indexed_trajectory(&path, Arc::clone(&topology), options).unwrap();
        let expected = frames(&mut built);
        assert!(expected.len() >= 2, "{name}");

        // Progress starts before the first frame, then advances one verified
        // frame at a time to the whole file.
        let steps = steps.lock().unwrap().clone();
        assert_eq!(steps.len(), expected.len() + 1, "{name}");
        for (count, step) in steps.iter().enumerate() {
            assert_eq!(step.frame_count(), count);
            assert_eq!(step.total_byte_count(), size);
        }
        assert!(steps
            .windows(2)
            .all(|pair| pair[0].byte_count() < pair[1].byte_count()));
        assert_eq!(steps.last().unwrap().byte_count(), size);

        let index = built.to_frame_index();
        assert_eq!(index.format(), built.metadata().format());
        assert_eq!(index.frame_count(), expected.len());
        assert_eq!(index.fingerprint().byte_count(), size);
        assert_eq!(
            index.fingerprint().modified(),
            fs::metadata(&path).unwrap().modified().ok()
        );
        // Steps and times are those the frames record, and survive encoding.
        let (has_steps, has_times) = match index.format() {
            TrajectoryFormat::Xyz => (false, false),
            TrajectoryFormat::Dcd => (true, false),
            _ => (true, true),
        };
        assert_eq!(index.steps().is_some(), has_steps, "{name}");
        assert_eq!(index.times().is_some(), has_times, "{name}");
        assert_stamps(&index, &expected);
        let bytes = index.to_bytes();
        let decoded =
            TrajectoryFrameIndex::from_bytes(&bytes, &TrajectoryIoLimits::default()).unwrap();
        assert_eq!(decoded, index);
        assert_stamps(&decoded, &expected);

        let mut reopened = reopen(&path, &topology, &bytes).unwrap();
        assert_eq!(reopened.metadata(), built.metadata(), "{name}");
        assert_eq!(reopened.to_frame_index(), index);
        assert!(reopened.open_report().notes()[0].contains("imported"));
        // Random reads in any order, then sequential reads, match the scan.
        let mut buffer = reopened.frame_buffer();
        for index in (0..expected.len()).rev() {
            reopened.read_frame(index, &mut buffer).unwrap();
            assert_eq!(buffer.frame_view().payload(), &expected[index], "{name}");
        }
        assert_eq!(support::read_all(&mut reopened), expected, "{name}");
        assert_eq!(
            reopened
                .read_frame(expected.len(), &mut buffer)
                .unwrap_err()
                .kind(),
            TrajectoryIoErrorKind::FrameIndexOutOfRange
        );
    }
}

#[test]
fn an_index_never_opens_a_changed_or_different_file() {
    let directory = Directory::new();
    let xtc = directory.fixture("mdanalysis-2.9.0-twelve-atoms.xtc");
    let trr = directory.fixture("mdanalysis-2.9.0-three-atoms.trr");
    let xtc_topology = linear_carbon_topology(12);
    let trr_topology = topology(&["C", "H", "O"], &[(0, 1), (0, 2)]);
    let bytes = open_indexed_trajectory(&xtc, Arc::clone(&xtc_topology), Default::default())
        .unwrap()
        .to_frame_index()
        .to_bytes();
    let mismatch = |path: &Path, topology: &Arc<Topology>| {
        let error = reopen(path, topology, &bytes).unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::IndexMismatch);
        assert_eq!(error.operation(), Some(TrajectoryIoOperation::Open));
        error
    };

    // Another format.
    let error = mismatch(&trr, &trr_topology);
    assert_eq!(error.format(), Some(TrajectoryFormat::Trr));

    // The same content rewritten in place keeps its size; restoring the
    // modification time leaves only the content hash to notice.
    let modified = fs::metadata(&xtc).unwrap().modified().unwrap();
    let original = fs::read(&xtc).unwrap();
    let mut changed = original.clone();
    let middle = changed.len() / 2;
    changed[middle] ^= 0x40;
    fs::write(&xtc, &changed).unwrap();
    File::options()
        .write(true)
        .open(&xtc)
        .unwrap()
        .set_modified(modified)
        .unwrap();
    mismatch(&xtc, &xtc_topology);

    // Restoring the bytes and the time makes the index valid again.
    fs::write(&xtc, &original).unwrap();
    File::options()
        .write(true)
        .open(&xtc)
        .unwrap()
        .set_modified(modified)
        .unwrap();
    reopen(&xtc, &xtc_topology, &bytes).unwrap();

    // An appended frame changes the size.
    let mut file = OpenOptions::new().append(true).open(&xtc).unwrap();
    file.write_all(&original[..original.len() / 2]).unwrap();
    drop(file);
    mismatch(&xtc, &xtc_topology);
}

#[test]
fn forged_offsets_fail_the_frames_they_misplace_and_no_other() {
    let directory = Directory::new();
    for (name, topology) in fixtures() {
        let path = directory.fixture(name);
        let mut built =
            open_indexed_trajectory(&path, Arc::clone(&topology), Default::default()).unwrap();
        let expected = frames(&mut built);
        let index = built.to_frame_index();
        let starts = index.frame_offsets().to_vec();

        // A first offset that is not the file's first frame fails to open.
        for delta in [1_u64, 4] {
            let mut forged = starts.clone();
            forged[0] += delta;
            let error = reopen(&path, &topology, &forge(&index, &forged)).unwrap_err();
            assert_eq!(error.kind(), TrajectoryIoErrorKind::IndexMismatch, "{name}");
        }

        // A later offset moved off its frame boundary fails the frame before
        // it, which no longer ends there, and the frame itself, which does
        // not decode; every other frame still reads correctly.
        for moved in 1..starts.len() {
            for delta in [-4_i64, -1, 1, 4] {
                let mut forged = starts.clone();
                forged[moved] = forged[moved].checked_add_signed(delta).unwrap();
                let mut reader = reopen(&path, &topology, &forge(&index, &forged)).unwrap();
                let mut buffer = reader.frame_buffer();
                for (frame, expected) in expected.iter().enumerate() {
                    let result = reader.read_frame(frame, &mut buffer);
                    if frame == moved - 1 {
                        let error = result.unwrap_err();
                        assert_eq!(error.kind(), TrajectoryIoErrorKind::IndexMismatch, "{name}");
                        assert_eq!(error.frame(), Some(frame));
                    } else if frame == moved {
                        assert!(result.is_err(), "{name} frame {frame} delta {delta}");
                    } else {
                        result.unwrap();
                        assert_eq!(buffer.frame_view().payload(), expected, "{name}");
                    }
                }
            }
        }

        // Out-of-order offsets never decode.
        let mut swapped = starts.clone();
        swapped.swap(0, 1);
        assert_eq!(
            kind(reopen(&path, &topology, &forge(&index, &swapped))),
            TrajectoryIoErrorKind::InvalidIndex
        );
    }
}

#[test]
fn a_dcd_index_must_count_the_declared_frames() {
    let directory = Directory::new();
    let path = directory.fixture("mdanalysis-2.9.0-three-atoms.dcd");
    let topology = topology(&["C", "H", "O"], &[(0, 1), (0, 2)]);
    let index = open_indexed_trajectory(&path, Arc::clone(&topology), Default::default())
        .unwrap()
        .to_frame_index();
    let bytes = forge(&index, &index.frame_offsets()[..1]);
    let error = reopen(&path, &topology, &bytes).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::IndexMismatch);
    assert_eq!((error.expected(), error.actual()), (Some(2), Some(1)));
}

#[test]
fn indexes_of_the_previous_version_are_refused() {
    let directory = Directory::new();
    for (name, topology) in fixtures() {
        let path = directory.fixture(name);
        let index = open_indexed_trajectory(&path, Arc::clone(&topology), Default::default())
            .unwrap()
            .to_frame_index();
        // Version 1 had the same header without step and time flags, then
        // the offsets alone.
        let mut previous = index.to_bytes()[..56 + 8 * index.frame_count()].to_vec();
        previous[8..12].copy_from_slice(&1_u32.to_le_bytes());
        previous[13] &= 1;
        previous.extend_from_slice(&[0; 8]);
        reencode(&mut previous);
        let error = reopen(&path, &topology, &previous).unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::InvalidIndex, "{name}");
        assert!(error.detail().unwrap().contains("version 1"), "{name}");
    }
}

#[test]
fn forged_times_and_steps_fail_their_frame_and_no_other() {
    let directory = Directory::new();
    for (name, topology) in &fixtures()[2..] {
        let path = directory.fixture(name);
        let mut built =
            open_indexed_trajectory(&path, Arc::clone(topology), Default::default()).unwrap();
        let expected = frames(&mut built);
        let bytes = built.to_frame_index().to_bytes();
        let count = expected.len();
        let steps_at = 56 + 8 * count;
        let times_at = steps_at + 8 * count;
        for (table, frame) in [(steps_at, 1), (times_at, count - 1)] {
            let mut forged = bytes.clone();
            forged[table + 8 * frame] ^= 1;
            reencode(&mut forged);
            let mut reader = reopen(&path, topology, &forged).unwrap();
            let mut buffer = reader.frame_buffer();
            for (index, expected) in expected.iter().enumerate() {
                let result = reader.read_frame(index, &mut buffer);
                if index == frame {
                    let error = result.unwrap_err();
                    assert_eq!(error.kind(), TrajectoryIoErrorKind::IndexMismatch, "{name}");
                    assert_eq!(error.operation(), Some(TrajectoryIoOperation::ReadFrame));
                    assert_eq!(error.frame(), Some(frame));
                } else {
                    result.unwrap();
                    assert_eq!(buffer.frame_view().payload(), expected, "{name}");
                }
            }
        }
    }
}

#[test]
fn dcd_indexes_record_times_under_the_header_delta_policy_they_were_built_with() {
    let directory = Directory::new();
    let path = directory.fixture("mdanalysis-2.9.0-three-atoms.dcd");
    let topology = topology(&["C", "H", "O"], &[(0, 1), (0, 2)]);
    let with_policy = |policy| {
        TrajectoryOpenOptions::default()
            .with_dcd_options(DcdReadOptions::default().with_time_policy(policy))
    };
    let timed = with_policy(DcdTimePolicy::HeaderDelta { unit: FEMTOSECOND });
    let mut built = open_indexed_trajectory(&path, Arc::clone(&topology), timed.clone()).unwrap();
    let expected = frames(&mut built);
    let index = built.to_frame_index();
    assert!(index.times().is_some());
    assert_stamps(&index, &expected);
    let bytes = index.to_bytes();
    assert_eq!(bytes[13], 0b110);
    let decoded = TrajectoryFrameIndex::from_bytes(&bytes, &TrajectoryIoLimits::default()).unwrap();
    assert_stamps(&decoded, &expected);

    // The same policy reopens it with the same frames.
    let mut reopened =
        open_indexed_trajectory_from_index(&path, Arc::clone(&topology), decoded.clone(), timed)
            .unwrap();
    assert_eq!(frames(&mut reopened), expected);
    assert_eq!(reopened.to_frame_index(), index);

    // Another policy would give the frames other times than the index.
    for policy in [
        DcdTimePolicy::Absent,
        DcdTimePolicy::HeaderDelta { unit: PICOSECOND },
    ] {
        let error = open_indexed_trajectory_from_index(
            &path,
            Arc::clone(&topology),
            decoded.clone(),
            with_policy(policy),
        )
        .unwrap_err();
        assert_eq!(error.kind(), TrajectoryIoErrorKind::IndexMismatch);
        assert_eq!(error.operation(), Some(TrajectoryIoOperation::Open));
        assert!(error.frame().is_some());
    }
}

#[test]
fn opening_from_an_index_respects_the_index_limits() {
    let directory = Directory::new();
    let path = directory.fixture("mdanalysis-2.9.0-twelve-atoms.xtc");
    let topology = linear_carbon_topology(12);
    let index = open_indexed_trajectory(&path, Arc::clone(&topology), Default::default())
        .unwrap()
        .to_frame_index();
    let limits = TrajectoryIoLimits::default().with_max_index_entries(index.frame_count() - 1);
    assert_eq!(
        kind(TrajectoryFrameIndex::from_bytes(&index.to_bytes(), &limits)),
        TrajectoryIoErrorKind::LimitExceeded
    );
    let error = open_indexed_trajectory_from_index(
        &path,
        topology,
        index,
        TrajectoryOpenOptions::default().with_limits(limits),
    )
    .unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::LimitExceeded);
    assert!(error.limit_exceeded().is_some());
}

#[test]
fn cancellation_stops_the_index_build_without_a_reader() {
    let directory = Directory::new();
    for (name, topology) in fixtures() {
        let path = directory.fixture(name);
        for after_frames in [0, 1] {
            let verified = Arc::new(AtomicUsize::new(0));
            let progress = Arc::clone(&verified);
            let observed = Arc::clone(&verified);
            let polls = Arc::new(AtomicUsize::new(0));
            let counted = Arc::clone(&polls);
            let options = TrajectoryOpenOptions::default()
                .with_index_progress(move |step| {
                    progress.store(step.frame_count(), Ordering::Relaxed);
                })
                .with_index_cancellation(Cancellation::new(move || {
                    counted.fetch_add(1, Ordering::Relaxed);
                    observed.load(Ordering::Relaxed) >= after_frames
                }));
            let error = open_indexed_trajectory(&path, Arc::clone(&topology), options).unwrap_err();
            assert_eq!(error.kind(), TrajectoryIoErrorKind::Cancelled, "{name}");
            assert_eq!(error.operation(), Some(TrajectoryIoOperation::Index));
            assert_eq!(error.frame(), Some(after_frames));
            assert!(std::error::Error::source(&error)
                .unwrap()
                .downcast_ref::<Cancelled>()
                .is_some());
            assert_eq!(polls.load(Ordering::Relaxed), after_frames + 1);
        }
        // A cancellation that is never requested lets the build finish.
        let options =
            TrajectoryOpenOptions::default().with_index_cancellation(Cancellation::new(|| false));
        open_indexed_trajectory(&path, Arc::clone(&topology), options).unwrap();
    }
}

fn mixed_precision_trr(path: &Path) -> Arc<Topology> {
    let topology = topology(&["C", "H", "O"], &[(0, 1), (0, 2)]);
    let mut combined = Vec::new();
    for (step, precision) in [
        TrrScalarPrecision::Float32,
        TrrScalarPrecision::Float64,
        TrrScalarPrecision::Float32,
    ]
    .into_iter()
    .enumerate()
    {
        let mut writer = TrrWriter::new(
            Cursor::new(Vec::new()),
            Arc::clone(&topology),
            TrrWriteOptions::default()
                .with_precision(precision)
                .with_lambda_policy(TrrLambdaPolicy::RequireZero),
            "mixed.trr",
        )
        .unwrap();
        let mut frame = TrajectoryFrame::new(
            Positions::from_vec(Quantity::new(
                vec![Point3::new(step as f64, 0.0, 0.0); 3],
                NANOMETER,
            ))
            .unwrap(),
        );
        frame.set_step(Some(step as u64));
        frame
            .set_time(Some(Quantity::new(step as f64, PICOSECOND)))
            .unwrap();
        let trajectory = Trajectory::from_items(Arc::clone(&topology), [frame]).unwrap();
        writer.write_frame(trajectory.get(0).unwrap()).unwrap();
        combined.extend(writer.finish().unwrap().into_inner());
    }
    fs::write(path, combined).unwrap();
    topology
}

#[test]
fn trr_indexes_keep_mixed_precision_and_verify_it_per_frame() {
    let directory = Directory::new();
    let path = directory.0.join("mixed.trr");
    let topology = mixed_precision_trr(&path);
    let mut built =
        open_indexed_trajectory(&path, Arc::clone(&topology), Default::default()).unwrap();
    let expected = frames(&mut built);
    let mixed = CoordinateEncoding::Lossless {
        precision: ScalarPrecision::Mixed,
    };
    assert_eq!(built.metadata().coordinate_encoding(), mixed);
    let bytes = built.to_frame_index().to_bytes();
    assert_eq!(bytes[13] & 1, 1);
    let mut reopened = reopen(&path, &topology, &bytes).unwrap();
    assert_eq!(reopened.metadata().coordinate_encoding(), mixed);
    assert_eq!(frames(&mut reopened), expected);

    // An index claiming uniform width rejects the f64 frame.
    let mut uniform = bytes.clone();
    uniform[13] &= !1;
    reencode(&mut uniform);
    let mut reader = reopen(&path, &topology, &uniform).unwrap();
    let mut buffer = reader.frame_buffer();
    reader.read_frame(0, &mut buffer).unwrap();
    reader.read_frame(2, &mut buffer).unwrap();
    let error = reader.read_frame(1, &mut buffer).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::IndexMismatch);
    assert_eq!(error.frame(), Some(1));
}

#[test]
fn a_file_that_grows_while_indexed_is_rejected() {
    let directory = Directory::new();
    let path = directory.fixture("mdanalysis-2.9.0-twelve-atoms.xtc");
    let topology = linear_carbon_topology(12);
    let original = fs::read(&path).unwrap();
    let appended = Arc::new(AtomicUsize::new(0));
    let target = path.clone();
    let once = Arc::clone(&appended);
    let options = TrajectoryOpenOptions::default().with_index_progress(move |step| {
        if step.frame_count() == 1 && once.fetch_add(1, Ordering::Relaxed) == 0 {
            // Append one more copy of the first frame's bytes after the scan
            // has passed the start, as a running simulation would.
            let mut file = OpenOptions::new().write(true).open(&target).unwrap();
            file.seek(SeekFrom::End(0)).unwrap();
            file.write_all(&original[..usize::try_from(step.byte_count()).unwrap()])
                .unwrap();
        }
    });
    let error = open_indexed_trajectory(&path, topology, options).unwrap_err();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::IndexMismatch);
    assert_eq!(error.operation(), Some(TrajectoryIoOperation::Index));
}
