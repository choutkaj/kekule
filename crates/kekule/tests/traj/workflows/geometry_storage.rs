//! Structure sets read from indexed trajectory files: every read is the indexed
//! reader's frame, reads run from many threads at once, a refresh takes the
//! frames a writer appended, and loading gives the set in memory.

use std::fs::{self, OpenOptions};
use std::io::{Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::thread;

use kekule::diagnostics::UnrepresentablePolicy;
use kekule::math::{Matrix3, PeriodicCell, Point3, RigidTransform, Vector3};
use kekule::native::StructureSetDescription;
use kekule::structure::{
    Distribution, Geometry, Members, Order, Partition, Positions, StructureSet, StructureSetRead,
    StructureSetReadError,
};
use kekule::topology::{AtomSelection, Topology};
use kekule::traj::io::dcd::{DcdCellPolicy, DcdEndian, DcdWriteOptions};
use kekule::traj::io::trr::{TrrLambdaPolicy, TrrWriteOptions};
use kekule::traj::io::xyz::XyzWriteOptions;
use kekule::traj::io::{
    append_trajectory, open_indexed_structure_set, open_indexed_structure_set_from_index,
    open_indexed_trajectory, read_structure_set, write_structure_set, IndexedFile, KeepFrames,
    TrajectoryIoLimits, TrajectoryOpenOptions, TrajectoryWriteOptions,
};
use kekule::traj::{
    FrameView, Processed, ProcessingSteps, SeekableTrajectoryReader, TrajectoryFormat,
    TrajectoryIoErrorKind, TrajectoryReader, TrajectoryWriter,
};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};

use crate::support::{
    fixed_atom_fixture, frames_of, linear_carbon_topology, time_ordered, topology, Frame,
};

/// The externally supplied fixture of every format with its topology.
fn fixtures() -> [(PathBuf, Arc<Topology>); 4] {
    let path = |name: &str| {
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests/fixtures/traj")
            .join(name)
    };
    [
        (
            path("ase-3.26.0-water.xyz"),
            topology(&["O", "H", "H"], &[(0, 1), (0, 2)]),
        ),
        (
            path("mdanalysis-2.9.0-three-atoms.dcd"),
            topology(&["C", "H", "O"], &[(0, 1), (0, 2)]),
        ),
        (
            path("mdanalysis-2.9.0-three-atoms.trr"),
            topology(&["C", "H", "O"], &[(0, 1), (0, 2)]),
        ),
        (
            path("mdanalysis-2.9.0-twelve-atoms.xtc"),
            linear_carbon_topology(12),
        ),
    ]
}

/// A scratch directory removed on drop.
struct Directory(PathBuf);

impl Directory {
    fn new() -> Self {
        static SEQUENCE: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "kekule-geometry-storage-{}-{}",
            std::process::id(),
            SEQUENCE.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }
}

impl Drop for Directory {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// `count` frames of a moving 12-carbon chain in a cell, with times and
/// steps.
fn chain(topology: &Arc<Topology>, count: u64) -> StructureSet {
    let frames = (0..count).map(|step| {
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
        let mut frame = Frame::new(Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap())
            .with_cell(cell);
        frame.set_step(Some(10 * step));
        frame
            .set_time(Some(Quantity::new(0.5 * shift, PICOSECOND)))
            .unwrap();
        frame
    });
    time_ordered(Arc::clone(topology), frames)
}

/// Every frame of the file at `path`, read through its indexed reader.
fn indexed_frames(path: &Path, topology: &Arc<Topology>) -> Vec<Frame> {
    let mut reader =
        open_indexed_trajectory(path, Arc::clone(topology), Default::default()).unwrap();
    let mut buffer = reader.frame_buffer();
    (0..reader.frame_count().unwrap())
        .map(|index| {
            reader.read_frame(index, &mut buffer).unwrap();
            Frame::of(&buffer)
        })
        .collect()
}

/// Checks that `geometry`, geometry `index` of `set`, is `frame` with its time
/// and step in the set's order.
fn assert_frame(set: &impl StructureSetRead, index: usize, geometry: &Geometry, frame: &Frame) {
    assert_eq!(geometry.positions(), frame.positions(), "geometry {index}");
    assert_eq!(geometry.cell(), frame.cell(), "geometry {index}");
    assert_eq!(
        geometry.velocities(),
        frame.velocities(),
        "geometry {index}"
    );
    assert_eq!(geometry.forces(), frame.forces(), "geometry {index}");
    assert_eq!(
        geometry.properties(),
        frame.properties(),
        "geometry {index}"
    );
    let Some(Order::Time(axis)) = set.order() else {
        panic!("a file's geometries are in time order");
    };
    assert_eq!(
        axis.time(index)
            .map(|time| time.value_in(PICOSECOND).unwrap()),
        frame.time().map(|time| time.value_in(PICOSECOND).unwrap()),
        "geometry {index}"
    );
    assert_eq!(axis.step(index), frame.step(), "geometry {index}");
}

#[test]
fn every_read_is_the_indexed_frame_of_every_format() {
    for (path, topology) in fixtures() {
        let frames = indexed_frames(&path, &topology);
        let set =
            open_indexed_structure_set(&path, Arc::clone(&topology), Default::default()).unwrap();
        assert_eq!(set.len(), frames.len(), "{path:?}");
        assert!(set.len() >= 2, "{path:?}");
        assert!(Arc::ptr_eq(&set.shared_topology(), &topology));
        assert!(set.distribution().is_none());
        assert_eq!(set.storage().path(), path);
        assert_eq!(set.storage().frame_index().frame_count(), frames.len());
        assert_eq!(
            set.storage().fingerprint(),
            set.storage().frame_index().fingerprint()
        );
        // Out of order, and twice.
        for index in (0..set.len()).rev().chain(0..set.len()) {
            assert!(set.peek(index).is_none());
            let geometry = set.read(index).unwrap();
            assert_frame(&set, index, &geometry, &frames[index]);
        }
        assert_eq!(
            set.read(set.len()),
            Err(StructureSetReadError::IndexOutOfRange {
                index: set.len(),
                len: set.len()
            })
        );

        // Opened again from the kept index, the set reads the same geometries.
        let index = set.storage().frame_index();
        let reopened = open_indexed_structure_set_from_index(
            &path,
            Arc::clone(&topology),
            index,
            Default::default(),
        )
        .unwrap();
        assert_eq!(reopened.order(), set.order());
        for (index, frame) in frames.iter().enumerate() {
            assert_frame(&reopened, index, &reopened.read(index).unwrap(), frame);
        }
    }
}

#[test]
fn loading_gives_the_set_in_memory_with_its_order_and_distribution() {
    for (path, topology) in fixtures() {
        let set =
            open_indexed_structure_set(&path, Arc::clone(&topology), Default::default()).unwrap();
        let loaded = set.load(0..set.len()).unwrap();
        let read = read_structure_set(&path, Arc::clone(&topology), Default::default()).unwrap();
        assert_eq!(loaded.len(), read.len(), "{path:?}");
        assert_eq!(loaded.order(), read.order(), "{path:?}");
        for index in 0..read.len() {
            assert_eq!(
                loaded.get(index).unwrap().payload(),
                read.get(index).unwrap().payload(),
                "{path:?} geometry {index}"
            );
        }
        let described = StructureSetDescription::of(&set);
        assert_eq!(described.geometry_count(), set.len());
        assert_eq!(described.order(), set.order());
    }

    // Listed geometries keep their order values and the distribution's
    // members among them, as selecting them in memory does.
    let topology = linear_carbon_topology(12);
    let directory = Directory::new();
    let path = directory.0.join("run.xtc");
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc);
    write_structure_set(&path, &chain(&topology, 8), options).unwrap();
    let mut set =
        open_indexed_structure_set(&path, Arc::clone(&topology), Default::default()).unwrap();
    let kept = (0..set.len())
        .filter(|index| index % 3 != 1)
        .collect::<Vec<_>>();
    let count = kept.len();
    let members = Members::Subset(kept);
    let weights = (0..count).map(|index| 1.0 + index as f64).collect();
    let distribution = Distribution::new(members, Partition::one_state(count))
        .unwrap()
        .with_weights(weights)
        .unwrap();
    set.set_distribution(Some(distribution)).unwrap();
    let mut memory = set.load(0..set.len()).unwrap();
    assert_eq!(memory.distribution(), set.distribution());
    let listed = [4, 1, 0, 4];
    let loaded = set.load(listed).unwrap();
    let selected = memory.select(listed).unwrap();
    assert_eq!(loaded.order(), selected.order());
    assert_eq!(loaded.distribution(), selected.distribution());
    for (position, index) in listed.into_iter().enumerate() {
        assert_eq!(
            loaded.get(position).unwrap().payload(),
            memory.get(index).unwrap().payload()
        );
    }
    assert!(matches!(
        set.load([0, set.len()]),
        Err(StructureSetReadError::IndexOutOfRange { .. })
    ));
    // The set in memory is the set's own: editing it leaves the file alone.
    memory.set_distribution(None).unwrap();
    assert!(set.distribution().is_some());
}

#[test]
fn reads_from_many_threads_at_once_equal_sequential_reads() {
    for (path, topology) in fixtures() {
        let set = open_indexed_structure_set(&path, topology, Default::default()).unwrap();
        let sequential = (0..set.len())
            .map(|index| set.read(index).unwrap())
            .collect::<Vec<_>>();
        let reader: &dyn StructureSetRead = &set;
        thread::scope(|scope| {
            let threads = (0..6)
                .map(|offset| {
                    let sequential = &sequential;
                    scope.spawn(move || {
                        for round in 0..20 {
                            let index = (offset * 7 + round * 3) % reader.len();
                            assert_eq!(reader.read(index).unwrap(), sequential[index]);
                        }
                    })
                })
                .collect::<Vec<_>>();
            for thread in threads {
                thread.join().unwrap();
            }
        });
    }
}

#[test]
fn refresh_takes_appended_frames_and_keeps_what_the_distribution_covered() {
    let topology = linear_carbon_topology(12);
    let all = chain(&topology, 8);
    let frames = frames_of(&all);
    let directory = Directory::new();
    let path = directory.0.join("run.xtc");
    let head = time_ordered(Arc::clone(&topology), frames[..3].iter().cloned());
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc);
    write_structure_set(&path, &head, options.clone()).unwrap();

    let growing = TrajectoryOpenOptions::default().with_growing(true);
    let mut set = open_indexed_structure_set(&path, Arc::clone(&topology), growing).unwrap();
    assert_eq!(set.len(), 3);
    assert_eq!(set.refresh().unwrap(), 0);
    let distribution = Distribution::new(Members::All, Partition::one_state(3)).unwrap();
    set.set_distribution(Some(distribution)).unwrap();
    let before = (0..3)
        .map(|index| set.read(index).unwrap())
        .collect::<Vec<_>>();

    // A restarted run continues the file after its three frames.
    let (mut writer, report) =
        append_trajectory(&path, Arc::clone(&topology), KeepFrames::First(3), options).unwrap();
    assert_eq!(report.kept_frames(), 3);
    for frame in all.iter().skip(3).take(2) {
        writer
            .write_frame(FrameView::from_geometry(frame).unwrap())
            .unwrap();
    }
    writer.flush().unwrap();
    assert_eq!(set.refresh().unwrap(), 2);
    assert_eq!(set.len(), 5);
    for frame in all.iter().skip(5) {
        writer
            .write_frame(FrameView::from_geometry(frame).unwrap())
            .unwrap();
    }
    writer.finish().unwrap();
    assert_eq!(set.refresh().unwrap(), frames.len() - 5);
    assert_eq!(set.len(), frames.len());

    let written = indexed_frames(&path, &topology);
    for (index, frame) in written.iter().enumerate() {
        assert_frame(&set, index, &set.read(index).unwrap(), frame);
    }
    for (index, geometry) in before.iter().enumerate() {
        assert_eq!(&set.read(index).unwrap(), geometry);
    }
    // The distribution over every geometry covers the geometries it covered.
    let distribution = set.distribution().unwrap();
    assert_eq!(distribution.members(), &Members::Subset(vec![0, 1, 2]));
    assert_eq!(set.storage().frame_index().frame_count(), frames.len());
    assert_eq!(
        open_indexed_structure_set(&path, topology, Default::default())
            .unwrap()
            .order(),
        set.order()
    );
}

/// Reads on many threads open readers while a writer appends; each refresh
/// takes the appended frame, and the readers kept from before read it.
#[test]
fn refresh_counts_from_the_set_after_concurrent_reads() {
    let topology = linear_carbon_topology(12);
    let all = chain(&topology, 24);
    let directory = Directory::new();
    let path = directory.0.join("run.xtc");
    let head = time_ordered(Arc::clone(&topology), frames_of(&all).into_iter().take(2));
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc);
    write_structure_set(&path, &head, options.clone()).unwrap();
    let growing = TrajectoryOpenOptions::default().with_growing(true);
    let mut set = open_indexed_structure_set(&path, Arc::clone(&topology), growing).unwrap();
    let (mut writer, _) =
        append_trajectory(&path, Arc::clone(&topology), KeepFrames::First(2), options).unwrap();
    for (count, frame) in all.iter().enumerate().skip(2) {
        writer
            .write_frame(FrameView::from_geometry(frame).unwrap())
            .unwrap();
        writer.flush().unwrap();
        let barrier = std::sync::Barrier::new(8);
        thread::scope(|scope| {
            for _ in 0..8 {
                let (set, barrier) = (&set, &barrier);
                scope.spawn(move || {
                    barrier.wait();
                    for _ in 0..4 {
                        set.read(0).unwrap();
                    }
                });
            }
        });
        assert_eq!(set.refresh().unwrap(), 1, "frame {count}");
        assert_eq!(set.len(), count + 1);
        assert_eq!(set.refresh().unwrap(), 0);
        thread::scope(|scope| {
            for _ in 0..8 {
                let set = &set;
                scope.spawn(move || {
                    set.read(count).unwrap();
                });
            }
        });
    }
    writer.finish().unwrap();
    let written = indexed_frames(&path, &topology);
    for (index, frame) in written.iter().enumerate() {
        assert_frame(&set, index, &set.read(index).unwrap(), frame);
    }
}

/// A TRR file of [`chain`] in `directory`.
fn chain_trr(directory: &Directory, topology: &Arc<Topology>, count: u64) -> PathBuf {
    let path = directory.0.join("run.trr");
    write_structure_set(&path, &chain(topology, count), trr_options()).unwrap();
    path
}

fn trr_options() -> TrajectoryWriteOptions {
    TrajectoryWriteOptions::new(TrajectoryFormat::Trr).with_trr_options(
        TrrWriteOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero),
    )
}

/// Reads every geometry of `set` on `threads` threads at once, each in its
/// own order, so that some readers are opened for the reads.
fn read_at_once(
    set: &StructureSet<IndexedFile>,
    threads: usize,
) -> Vec<Vec<Result<Geometry, StructureSetReadError>>> {
    let barrier = std::sync::Barrier::new(threads);
    thread::scope(|scope| {
        let runs = (0..threads)
            .map(|offset| {
                let barrier = &barrier;
                scope.spawn(move || {
                    barrier.wait();
                    let mut reads = (0..set.len())
                        .map(|step| {
                            let index = (offset + step) % set.len();
                            (index, set.read(index))
                        })
                        .collect::<Vec<_>>();
                    reads.sort_by_key(|(index, _)| *index);
                    reads.into_iter().map(|(_, read)| read).collect()
                })
            })
            .collect::<Vec<_>>();
        runs.into_iter().map(|run| run.join().unwrap()).collect()
    })
}

#[test]
fn a_failed_read_leaves_the_reader_reading() {
    let topology = linear_carbon_topology(12);
    let directory = Directory::new();
    let path = chain_trr(&directory, &topology, 6);
    let frames = indexed_frames(&path, &topology);
    let set = open_indexed_structure_set(&path, Arc::clone(&topology), Default::default()).unwrap();
    // The last coordinates of frame 3 become NaN.
    let end = set.storage().frame_index().frame_offsets()[4];
    let mut file = OpenOptions::new().write(true).open(&path).unwrap();
    file.seek(SeekFrom::Start(end - 8)).unwrap();
    file.write_all(&[0xff; 8]).unwrap();
    drop(file);

    let failed = set.read(3).unwrap_err();
    assert!(
        matches!(
            failed,
            StructureSetReadError::Trajectory { geometry: 3, .. }
        ),
        "{failed:?}"
    );
    assert_eq!(failed.to_string(), "Cannot read geometry 3.");
    // The same reader, back in the pool, reads the other frames.
    for index in [2, 4, 0, 5, 1] {
        assert_frame(&set, index, &set.read(index).unwrap(), &frames[index]);
    }
    assert_eq!(set.read(3).unwrap_err(), failed);
}

/// The lost frames of a truncated file fail with one error, whether a kept
/// reader or a reader opened for the read reads them.
#[test]
fn a_truncated_file_fails_the_frames_it_lost_alike() {
    let topology = linear_carbon_topology(12);
    let directory = Directory::new();
    let path = chain_trr(&directory, &topology, 8);
    let mut set =
        open_indexed_structure_set(&path, Arc::clone(&topology), Default::default()).unwrap();
    let starts = set.storage().frame_index().frame_offsets().to_vec();
    let expected = (0..set.len())
        .map(|index| set.read(index).unwrap())
        .collect::<Vec<_>>();
    // Several readers are kept.
    read_at_once(&set, 4);

    for (cut, first_lost) in [(starts[5] + 10, 5), (0, 0)] {
        OpenOptions::new()
            .write(true)
            .open(&path)
            .unwrap()
            .set_len(cut)
            .unwrap();
        let lost = (first_lost..set.len())
            .map(|index| set.read(index).unwrap_err())
            .collect::<Vec<_>>();
        for (index, error) in (first_lost..).zip(&lost) {
            let StructureSetReadError::Trajectory {
                geometry,
                source: error,
            } = error
            else {
                panic!("frame {index}: {error:?}");
            };
            assert_eq!(*geometry, index);
            assert_eq!(error.kind(), TrajectoryIoErrorKind::IndexMismatch);
            assert_eq!(error.frame(), Some(index));
            let detail = error.detail().unwrap();
            assert!(detail.ends_with("Open the file again."), "{detail}");
            assert!(detail.contains(&format!("it is {cut} bytes")), "{detail}");
        }
        for reads in read_at_once(&set, 12) {
            for (index, read) in reads.into_iter().enumerate() {
                match read {
                    Ok(geometry) => {
                        assert!(index < first_lost, "frame {index}");
                        assert_eq!(geometry, expected[index]);
                    }
                    Err(error) => assert_eq!(error, lost[index - first_lost], "frame {index}"),
                }
            }
        }
        let refused = set.refresh().unwrap_err();
        assert_eq!(refused.kind(), TrajectoryIoErrorKind::IndexMismatch);
        assert!(refused.detail().unwrap().ends_with("Open the file again."));
        assert_eq!(set.len(), expected.len());
    }
}

/// A set keeps reading the file it opened when another file replaces it,
/// and its refresh says to open the new file.
#[test]
fn a_replaced_file_reads_as_before_and_refresh_fails() {
    let topology = linear_carbon_topology(12);
    let directory = Directory::new();
    let path = chain_trr(&directory, &topology, 6);
    let mut set =
        open_indexed_structure_set(&path, Arc::clone(&topology), Default::default()).unwrap();
    let expected = (0..set.len())
        .map(|index| set.read(index).unwrap())
        .collect::<Vec<_>>();
    let other = directory.0.join("other.trr");
    let mut reversed = frames_of(&chain(&topology, 9));
    reversed.reverse();
    let reversed = time_ordered(Arc::clone(&topology), reversed);
    write_structure_set(&other, &reversed, trr_options()).unwrap();
    fs::rename(&other, &path).unwrap();

    for reads in read_at_once(&set, 12) {
        for (index, read) in reads.into_iter().enumerate() {
            assert_eq!(read.unwrap(), expected[index], "frame {index}");
        }
    }
    let refused = set.refresh().unwrap_err();
    assert_eq!(refused.kind(), TrajectoryIoErrorKind::IndexMismatch);
    assert!(refused.detail().unwrap().ends_with("Open the file again."));
    assert_eq!(set.len(), expected.len());
    // Opened again, the set reads the new file.
    let reopened = open_indexed_structure_set(&path, topology, Default::default()).unwrap();
    assert_eq!(reopened.len(), 9);
}

#[test]
fn refresh_keeps_a_distribution_over_some_samples_or_none() {
    let topology = linear_carbon_topology(12);
    let all = chain(&topology, 6);
    let frames = frames_of(&all);
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc);
    let growing = TrajectoryOpenOptions::default().with_growing(true);
    for members in [vec![0, 2], vec![]] {
        let directory = Directory::new();
        let path = directory.0.join("run.xtc");
        let head = time_ordered(Arc::clone(&topology), frames[..3].iter().cloned());
        write_structure_set(&path, &head, options.clone()).unwrap();
        let mut set =
            open_indexed_structure_set(&path, Arc::clone(&topology), growing.clone()).unwrap();
        let partition = Partition::one_state(members.len());
        let distribution = Distribution::new(Members::Subset(members.clone()), partition).unwrap();
        set.set_distribution(Some(distribution.clone())).unwrap();
        let (mut writer, _) = append_trajectory(
            &path,
            Arc::clone(&topology),
            KeepFrames::First(3),
            options.clone(),
        )
        .unwrap();
        for frame in all.iter().skip(3) {
            writer
                .write_frame(FrameView::from_geometry(frame).unwrap())
                .unwrap();
        }
        writer.finish().unwrap();
        assert_eq!(set.refresh().unwrap(), 3);
        assert_eq!(set.len(), 6);
        assert_eq!(set.distribution(), Some(&distribution), "{members:?}");
        assert_eq!(set.states().count(), usize::from(!members.is_empty()));
    }
}

/// A host following a growing file takes the set out of its processed
/// reader, refreshes it, and processes it again.
#[test]
fn a_processed_growing_file_reads_the_appended_samples() {
    let topology = linear_carbon_topology(12);
    let all = chain(&topology, 7);
    let directory = Directory::new();
    let path = directory.0.join("run.xtc");
    let head = time_ordered(Arc::clone(&topology), frames_of(&all).into_iter().take(3));
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc);
    write_structure_set(&path, &head, options.clone()).unwrap();
    let growing = TrajectoryOpenOptions::default().with_growing(true);
    let set = open_indexed_structure_set(&path, Arc::clone(&topology), growing).unwrap();
    let ids = topology.atom_ids();
    let middle =
        AtomSelection::from_atoms(Arc::clone(&topology), ids[4..8].iter().copied()).unwrap();
    let steps = ProcessingSteps::default()
        .with_wrap(middle.clone())
        .with_center(middle);
    let mut processed = Processed::new(set, steps).unwrap();
    let before = (0..3)
        .map(|index| processed.read(index).unwrap())
        .collect::<Vec<_>>();

    let (mut writer, _) =
        append_trajectory(&path, Arc::clone(&topology), KeepFrames::First(3), options).unwrap();
    for frame in all.iter().skip(3) {
        writer
            .write_frame(FrameView::from_geometry(frame).unwrap())
            .unwrap();
    }
    writer.finish().unwrap();
    let (mut set, steps) = processed.into_parts();
    assert_eq!(set.refresh().unwrap(), 4);
    processed = Processed::new(set, steps.clone()).unwrap();
    assert_eq!(processed.len(), 7);

    let reopened = open_indexed_structure_set(&path, topology, Default::default()).unwrap();
    let expected = Processed::new(reopened, steps).unwrap();
    for index in 0..7 {
        let geometry = processed.read(index).unwrap();
        assert_eq!(geometry, expected.read(index).unwrap(), "geometry {index}");
        if let Some(read) = before.get(index) {
            assert_eq!(&geometry, read, "geometry {index}");
        }
    }
}

/// Wrapping and making whole leave geometries without a cell as they are.
#[test]
fn wrapping_a_file_without_cells_reads_the_samples_as_stored() {
    let (path, topology) = fixtures()[0].clone();
    let set = open_indexed_structure_set(&path, Arc::clone(&topology), Default::default()).unwrap();
    let steps = ProcessingSteps::default()
        .with_whole(true)
        .with_wrap(AtomSelection::all(Arc::clone(&topology)));
    let processed = Processed::new(&set, steps).unwrap();
    for index in 0..set.len() {
        let stored = set.read(index).unwrap();
        assert!(stored.cell().is_none(), "geometry {index}");
        let read = processed.read(index).unwrap();
        assert!(processed.skipped_periodic_steps(&read), "geometry {index}");
        assert_eq!(read, stored, "geometry {index}");
    }
    // Geometries with a cell are made whole.
    let chain_topology = linear_carbon_topology(12);
    let boxed = chain(&chain_topology, 2);
    let whole = Processed::new(&boxed, ProcessingSteps::default().with_whole(true)).unwrap();
    assert!(!whole.skipped_periodic_steps(&whole.read(1).unwrap()));
    // Without making whole or wrapping, nothing is skipped.
    let centered = Processed::new(
        &set,
        ProcessingSteps::default().with_center(AtomSelection::all(topology)),
    )
    .unwrap();
    assert!(!centered.skipped_periodic_steps(&centered.read(0).unwrap()));
}

#[test]
fn a_failed_refresh_changes_nothing() {
    let (fixture, topology) = fixtures()[3].clone();
    let directory = Directory::new();
    let path = directory.0.join("run.xtc");
    fs::copy(&fixture, &path).unwrap();
    let growing = TrajectoryOpenOptions::default().with_growing(true);
    let mut set = open_indexed_structure_set(&path, Arc::clone(&topology), growing).unwrap();
    let len = set.len();
    let order = set.order().cloned();
    // Another trajectory replaces the indexed bytes.
    let mut bytes = fs::read(&path).unwrap();
    let last = bytes.len() - 5;
    bytes[last] ^= 0x01;
    bytes.extend_from_slice(&fs::read(&fixture).unwrap());
    fs::write(&path, &bytes).unwrap();
    let refused = set.refresh().unwrap_err();
    assert_eq!(refused.kind(), TrajectoryIoErrorKind::IndexMismatch);
    assert!(refused.detail().unwrap().ends_with("Open the file again."));
    assert_eq!(set.len(), len);
    assert_eq!(set.order(), order.as_ref());
}

/// A set in memory and the same geometries in a file read alike as readers.
#[test]
fn a_set_in_memory_lends_what_a_file_reads() {
    let (path, topology) = fixtures()[3].clone();
    let file = open_indexed_structure_set(&path, topology, Default::default()).unwrap();
    let memory: StructureSet = file.load(0..file.len()).unwrap();
    let readers: [&dyn StructureSetRead; 2] = [&file, &memory];
    for index in 0..file.len() {
        let [file_sample, memory_sample] = readers.map(|reader| reader.read(index).unwrap());
        assert_eq!(file_sample, memory_sample);
        assert!(memory
            .peek(index)
            .unwrap()
            .payload()
            .shares_structure_with(&memory_sample));
    }
    assert_eq!(readers[0].order(), readers[1].order());
}

/// `set` with every atom moved by `dx` nanometres along x.
fn shifted(set: &StructureSet, dx: f64) -> StructureSet {
    let mut shifted = set.clone();
    let transform = RigidTransform::new(
        Matrix3::identity(),
        Quantity::new(Vector3::new(dx, 0.0, 0.0), NANOMETER),
    )
    .unwrap();
    for index in 0..shifted.len() {
        shifted
            .get_mut(index)
            .unwrap()
            .apply_rigid_transform(&transform)
            .unwrap();
    }
    shifted
}

/// The set reads whatever bytes the file holds at the indexed offsets: a
/// frame rewritten in place at the same size reads as it now is, until a
/// refresh finds the indexed bytes changed.
#[test]
fn a_frame_rewritten_in_place_reads_as_it_now_is() {
    let topology = linear_carbon_topology(12);
    let directory = Directory::new();
    let path = chain_trr(&directory, &topology, 6);
    let other = directory.0.join("other.trr");
    write_structure_set(&other, &shifted(&chain(&topology, 6), 0.05), trr_options()).unwrap();
    let rewritten = indexed_frames(&other, &topology);
    let mut set =
        open_indexed_structure_set(&path, Arc::clone(&topology), Default::default()).unwrap();
    let before = (0..set.len())
        .map(|index| set.read(index).unwrap())
        .collect::<Vec<_>>();
    read_at_once(&set, 4);

    // Frame 2 of the other file, at the same offsets.
    let starts = set.storage().frame_index().frame_offsets().to_vec();
    let bytes = fs::read(&other).unwrap();
    let frame = &bytes[starts[2] as usize..starts[3] as usize];
    let mut file = OpenOptions::new().write(true).open(&path).unwrap();
    file.seek(SeekFrom::Start(starts[2])).unwrap();
    file.write_all(frame).unwrap();
    drop(file);

    for reads in read_at_once(&set, 8) {
        for (index, read) in reads.into_iter().enumerate() {
            let geometry = read.unwrap();
            if index == 2 {
                assert_frame(&set, index, &geometry, &rewritten[index]);
                assert_ne!(geometry, before[index]);
            } else {
                assert_eq!(geometry, before[index], "frame {index}");
            }
        }
    }
    let refused = set.refresh().unwrap_err();
    assert_eq!(refused.kind(), TrajectoryIoErrorKind::IndexMismatch);
    assert!(refused.detail().unwrap().ends_with("Open the file again."));
}

/// A file replaced by one that begins with the indexed frames is taken by
/// the next refresh, and every reader reads the new file.
#[test]
fn refresh_takes_a_replacement_that_begins_with_the_indexed_frames() {
    let topology = linear_carbon_topology(12);
    let all = chain(&topology, 6);
    let frames = frames_of(&all);
    let directory = Directory::new();
    let path = directory.0.join("run.xtc");
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc);
    let head = time_ordered(Arc::clone(&topology), frames[..3].iter().cloned());
    write_structure_set(&path, &head, options.clone()).unwrap();
    let growing = TrajectoryOpenOptions::default().with_growing(true);
    let mut set = open_indexed_structure_set(&path, Arc::clone(&topology), growing).unwrap();
    read_at_once(&set, 4);

    let other = directory.0.join("other.xtc");
    write_structure_set(&other, &all, options).unwrap();
    assert!(fs::read(&other)
        .unwrap()
        .starts_with(&fs::read(&path).unwrap()));
    fs::rename(&other, &path).unwrap();

    assert_eq!(set.refresh().unwrap(), 3);
    let written = indexed_frames(&path, &topology);
    for reads in read_at_once(&set, 8) {
        for (index, read) in reads.into_iter().enumerate() {
            assert_frame(&set, index, &read.unwrap(), &written[index]);
        }
    }
}

/// Write options keeping what each format holds of [`chain`].
fn growing_options(format: TrajectoryFormat) -> TrajectoryWriteOptions {
    TrajectoryWriteOptions::new(format)
        .with_dcd_options(
            DcdWriteOptions::default()
                .with_cell_policy(DcdCellPolicy::Write)
                .with_step_sequence(0, 10)
                .with_frame_metadata(UnrepresentablePolicy::Omit),
        )
        .with_xyz_options(
            XyzWriteOptions::default()
                .with_cell(UnrepresentablePolicy::Omit)
                .with_frame_metadata(UnrepresentablePolicy::Omit),
        )
}

#[test]
fn xyz_and_dcd_sets_take_appended_frames_on_refresh() {
    let topology = linear_carbon_topology(12);
    let all = chain(&topology, 6);
    let frames = frames_of(&all);
    for (format, name) in [
        (TrajectoryFormat::Xyz, "run.xyz"),
        (TrajectoryFormat::Dcd, "run.dcd"),
    ] {
        let directory = Directory::new();
        let path = directory.0.join(name);
        let options = growing_options(format);
        let head = time_ordered(Arc::clone(&topology), frames[..3].iter().cloned());
        write_structure_set(&path, &head, options.clone()).unwrap();
        let growing = TrajectoryOpenOptions::default().with_growing(true);
        let mut set = open_indexed_structure_set(&path, Arc::clone(&topology), growing).unwrap();
        assert_eq!(set.len(), 3, "{format}");
        let before = read_at_once(&set, 4).remove(0);

        let (mut writer, _) =
            append_trajectory(&path, Arc::clone(&topology), KeepFrames::First(3), options).unwrap();
        for frame in all.iter().skip(3).take(2) {
            writer
                .write_frame(FrameView::from_geometry(frame).unwrap())
                .unwrap();
        }
        writer.flush().unwrap();
        assert_eq!(set.refresh().unwrap(), 2, "{format}");
        for frame in all.iter().skip(5) {
            writer
                .write_frame(FrameView::from_geometry(frame).unwrap())
                .unwrap();
        }
        writer.finish().unwrap();
        assert_eq!(set.refresh().unwrap(), 1, "{format}");
        assert_eq!(set.len(), 6, "{format}");

        let written = indexed_frames(&path, &topology);
        for reads in read_at_once(&set, 8) {
            for (index, read) in reads.into_iter().enumerate() {
                assert_frame(&set, index, &read.unwrap(), &written[index]);
            }
        }
        for (index, read) in before.into_iter().enumerate() {
            assert_eq!(set.read(index).unwrap(), read.unwrap(), "{format} {index}");
        }
    }
}

/// A DCD with fixed atoms fills each later frame's fixed atoms from frame 0,
/// read through a set, from many threads, and processed.
#[test]
fn a_dcd_with_fixed_atoms_reads_whole_frames_through_a_set() {
    let topology = topology(&["C", "H", "O"], &[(0, 1), (0, 2)]);
    let directory = Directory::new();
    let path = directory.0.join("fixed.dcd");
    fs::write(&path, fixed_atom_fixture(DcdEndian::Little)).unwrap();
    let frames = indexed_frames(&path, &topology);
    let set = open_indexed_structure_set(&path, Arc::clone(&topology), Default::default()).unwrap();
    assert_eq!(set.len(), 2);
    // Frame 1 first, so the fixed atoms do not come from a read of frame 0.
    for index in [1, 0, 1] {
        assert_frame(&set, index, &set.read(index).unwrap(), &frames[index]);
    }
    let points = set.read(1).unwrap().positions().values().value().to_vec();
    let expected = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 1.1, 1.2),
        Point3::new(0.2, 0.2, 0.2),
    ];
    for (point, expected) in points.iter().zip(expected) {
        assert!((*point - expected).norm() < 1.0e-6, "{point:?}");
    }
    for reads in read_at_once(&set, 6) {
        for (index, read) in reads.into_iter().enumerate() {
            assert_frame(&set, index, &read.unwrap(), &frames[index]);
        }
    }

    let memory = set.load(0..set.len()).unwrap();
    let steps = ProcessingSteps::default().with_center(AtomSelection::all(Arc::clone(&topology)));
    let from_file = Processed::new(&set, steps.clone()).unwrap();
    let from_memory = Processed::new(&memory, steps).unwrap();
    for index in [1, 0, 1] {
        assert_eq!(
            from_file.read(index).unwrap(),
            from_memory.read(index).unwrap(),
            "geometry {index}"
        );
    }
}

/// A set made from a reader already open reads as one opened from the path,
/// and a clone keeps its geometries and its file while the set follows the
/// file.
#[test]
fn a_set_from_an_open_reader_and_its_clone_read_without_opening_the_file() {
    let topology = linear_carbon_topology(12);
    let all = chain(&topology, 6);
    let directory = Directory::new();
    let path = directory.0.join("run.xtc");
    let head = time_ordered(Arc::clone(&topology), frames_of(&all).into_iter().take(3));
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc);
    write_structure_set(&path, &head, options.clone()).unwrap();

    let growing = TrajectoryOpenOptions::default().with_growing(true);
    let reader = open_indexed_trajectory(&path, Arc::clone(&topology), growing.clone()).unwrap();
    let mut set = reader.into_structure_set();
    let opened = open_indexed_structure_set(&path, Arc::clone(&topology), growing).unwrap();
    assert_eq!(set.len(), 3);
    assert_eq!(set.order(), opened.order());
    for index in 0..3 {
        assert_eq!(set.read(index).unwrap(), opened.read(index).unwrap());
    }
    let distribution = Distribution::new(Members::All, Partition::one_state(3)).unwrap();
    set.set_distribution(Some(distribution)).unwrap();
    let expected = (0..3)
        .map(|index| set.read(index).unwrap())
        .collect::<Vec<_>>();

    let mut copy = set.clone();
    assert_eq!(copy.order(), set.order());
    assert_eq!(copy.distribution(), set.distribution());
    let (mut writer, _) =
        append_trajectory(&path, Arc::clone(&topology), KeepFrames::First(3), options).unwrap();
    for frame in all.iter().skip(3) {
        writer
            .write_frame(FrameView::from_geometry(frame).unwrap())
            .unwrap();
    }
    writer.finish().unwrap();
    assert_eq!(set.refresh().unwrap(), 3);
    assert_eq!(set.len(), 6);
    assert_eq!(copy.len(), 3);
    assert_eq!(copy.storage().frame_index().frame_count(), 3);
    assert_eq!(copy.distribution().unwrap().members(), &Members::All);

    // Another file at the path: the clone reads the file it shares.
    let other = directory.0.join("other.xtc");
    write_structure_set(
        &other,
        &shifted(&head, 0.5),
        TrajectoryWriteOptions::new(TrajectoryFormat::Xtc),
    )
    .unwrap();
    fs::rename(&other, &path).unwrap();
    for (index, geometry) in expected.iter().enumerate() {
        assert_eq!(&copy.read(index).unwrap(), geometry, "geometry {index}");
    }
    // Refreshed, the clone follows the path on its own.
    let refused = copy.refresh().unwrap_err();
    assert_eq!(refused.kind(), TrajectoryIoErrorKind::IndexMismatch);
    assert_eq!(copy.len(), 3);
}

/// A clone and its set each take the frames appended to the file on
/// their own refresh, and the other keeps the geometries it had.
#[test]
fn a_clone_and_its_set_grow_apart() {
    let topology = linear_carbon_topology(12);
    let all = chain(&topology, 6);
    let directory = Directory::new();
    let path = directory.0.join("run.xtc");
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc);
    let first = |count: usize| {
        time_ordered(
            Arc::clone(&topology),
            frames_of(&all).into_iter().take(count),
        )
    };
    write_structure_set(&path, &first(2), options.clone()).unwrap();
    let append = |from: usize, to: usize| {
        let (mut writer, _) = append_trajectory(
            &path,
            Arc::clone(&topology),
            KeepFrames::First(from),
            options.clone(),
        )
        .unwrap();
        for frame in all.iter().take(to).skip(from) {
            writer
                .write_frame(FrameView::from_geometry(frame).unwrap())
                .unwrap();
        }
        writer.finish().unwrap();
    };
    let growing = TrajectoryOpenOptions::default().with_growing(true);
    let mut set = open_indexed_structure_set(&path, Arc::clone(&topology), growing).unwrap();
    let mut copy = set.clone();

    append(2, 4);
    assert_eq!(copy.refresh().unwrap(), 2);
    assert_eq!(copy.len(), 4);
    assert_eq!(set.len(), 2);
    assert_eq!(set.storage().frame_index().frame_count(), 2);

    append(4, 6);
    assert_eq!(set.refresh().unwrap(), 4);
    assert_eq!(set.len(), 6);
    assert_eq!(copy.len(), 4);
    assert_eq!(copy.storage().frame_index().frame_count(), 4);

    let frames = indexed_frames(&path, &topology);
    for (index, frame) in frames.iter().enumerate() {
        assert_frame(&set, index, &set.read(index).unwrap(), frame);
        if index < 4 {
            assert_frame(&copy, index, &copy.read(index).unwrap(), frame);
        }
    }
    assert!(copy.read(4).is_err());
}

/// Two frames of a chain in a cell: frame 0 on a regular line, frame 1
/// scattered far from it, so that frame 1 needs more scratch to read in
/// every format but DCD, whose frames all have one size.
fn uneven(topology: &Arc<Topology>) -> StructureSet {
    let atoms = topology.atom_count();
    let frames = (0..2_u32).map(|step| {
        let points = (0..atoms)
            .map(|atom| {
                let atom = atom as f64;
                if step == 0 {
                    Point3::new(0.1 * atom, 0.0, 0.0)
                } else {
                    Point3::new(
                        -123.456 + 37.7 * atom.sin(),
                        456.789 * atom.cos(),
                        -98.765 + 0.11 * atom * atom.sqrt(),
                    )
                }
            })
            .collect::<Vec<_>>();
        let cell = PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(3.0, 3.5, 4.0), NANOMETER),
            [true; 3],
        )
        .unwrap();
        let mut frame = Frame::new(Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap())
            .with_cell(cell);
        frame.set_step(Some(10 * u64::from(step)));
        frame
            .set_time(Some(Quantity::new(f64::from(step), PICOSECOND)))
            .unwrap();
        frame
    });
    time_ordered(Arc::clone(topology), frames)
}

/// Which reader reads geometry 1 in [`reads_under`].
#[derive(Debug, Clone, Copy)]
enum History {
    /// A fresh set.
    Fresh,
    /// The reader that read geometry 0 just before.
    ReadBefore,
    /// A pooled reader that read geometry 0 on another thread.
    Pooled,
}

/// A file as a set opens it, and as it reads it after the bytes were
/// rewritten in place.
struct Rewritten<'a> {
    path: &'a Path,
    opened: Vec<u8>,
    read: Vec<u8>,
}

/// Whether geometry 1 of `file` opens and reads under a scratch limit of
/// `bytes`, by the reader `history` names.
fn reads_under(
    file: &Rewritten<'_>,
    topology: &Arc<Topology>,
    bytes: usize,
    history: History,
) -> bool {
    fs::write(file.path, &file.opened).unwrap();
    let limits = TrajectoryIoLimits::default().with_max_scratch_bytes(bytes);
    let options = TrajectoryOpenOptions::default().with_limits(limits);
    let Ok(set) = open_indexed_structure_set(file.path, Arc::clone(topology), options) else {
        return false;
    };
    let mut rewrite = OpenOptions::new().write(true).open(file.path).unwrap();
    rewrite.write_all(&file.read).unwrap();
    drop(rewrite);
    match history {
        History::Fresh => {}
        History::ReadBefore => drop(set.read(0)),
        History::Pooled => thread::scope(|scope| {
            scope.spawn(|| drop(set.read(0)));
        }),
    }
    set.read(1).is_ok()
}

/// The least scratch limit under which [`reads_under`] reads geometry 1.
fn least_limit(file: &Rewritten<'_>, topology: &Arc<Topology>, history: History) -> usize {
    let (mut fails, mut reads) = (0, 1 << 20);
    assert!(reads_under(file, topology, reads, history), "{history:?}");
    while reads - fails > 1 {
        let middle = fails + (reads - fails) / 2;
        if reads_under(file, topology, middle, history) {
            reads = middle;
        } else {
            fails = middle;
        }
    }
    reads
}

/// `xyz` with frame 1's first atom line 200 bytes longer and the 200 lines
/// after it a byte shorter each, so that every frame keeps its bytes.
fn one_long_line(xyz: &[u8]) -> Vec<u8> {
    let lines = xyz
        .split_inclusive(|byte| *byte == b'\n')
        .collect::<Vec<_>>();
    let atoms = lines.len() / 2 - 2;
    let first = atoms + 4;
    let mut rewritten = lines[..first].concat();
    let long = lines[first];
    rewritten.extend_from_slice(&long[..long.len() - 1]);
    rewritten.extend(std::iter::repeat_n(b'0', 200));
    rewritten.push(b'\n');
    for line in &lines[first + 1..first + 201] {
        assert!(line[line.len() - 2].is_ascii_digit());
        rewritten.extend_from_slice(&line[..line.len() - 2]);
        rewritten.push(b'\n');
    }
    rewritten.extend(lines[first + 201..].concat());
    assert_eq!(rewritten.len(), xyz.len());
    rewritten
}

/// The scratch limit counts the arrays a read allocates, so a frame passes
/// or fails it alike on a fresh reader, on one that read before, and on a
/// pooled one, which gave the arrays of its last read away. Opening checks
/// every frame as indexed, so the XYZ file is also read after a line of
/// frame 1 grew in place.
#[test]
fn a_frame_meets_the_scratch_limit_alike_whatever_the_reader_read_before() {
    let topology = linear_carbon_topology(600);
    let trajectory = uneven(&topology);
    let directory = Directory::new();
    for (format, rewrite) in [
        (TrajectoryFormat::Xyz, false),
        (TrajectoryFormat::Xyz, true),
        (TrajectoryFormat::Dcd, false),
        (TrajectoryFormat::Trr, false),
        (TrajectoryFormat::Xtc, false),
    ] {
        let path = directory.0.join(format!(
            "uneven-{rewrite}.{}",
            format.to_string().to_ascii_lowercase()
        ));
        let options = growing_options(format).with_trr_options(
            TrrWriteOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero),
        );
        write_structure_set(&path, &trajectory, options).unwrap();
        let opened = fs::read(&path).unwrap();
        let read = if rewrite {
            one_long_line(&opened)
        } else {
            opened.clone()
        };
        let file = Rewritten {
            path: &path,
            opened,
            read,
        };
        let fresh = least_limit(&file, &topology, History::Fresh);
        assert!(fresh > 600 * size_of::<Point3>(), "{format}: {fresh}");
        for history in [History::ReadBefore, History::Pooled] {
            assert_eq!(
                least_limit(&file, &topology, history),
                fresh,
                "{format} {rewrite} {history:?}"
            );
        }
        // Just under it, the read fails for the limit.
        if reads_under(&file, &topology, fresh - 1, History::Fresh) {
            panic!("{format}: geometry 1 reads under {} bytes", fresh - 1);
        }
        fs::write(&path, &file.opened).unwrap();
        let limits = TrajectoryIoLimits::default().with_max_scratch_bytes(fresh - 1);
        if let Ok(set) = open_indexed_structure_set(
            &path,
            Arc::clone(&topology),
            TrajectoryOpenOptions::default().with_limits(limits),
        ) {
            let mut rewrite = OpenOptions::new().write(true).open(&path).unwrap();
            rewrite.write_all(&file.read).unwrap();
            drop(rewrite);
            let error = set.read(1).unwrap_err();
            let StructureSetReadError::Trajectory { geometry, source } = error else {
                panic!("{format}: {error:?}");
            };
            assert_eq!(geometry, 1);
            assert_eq!(
                source.kind(),
                TrajectoryIoErrorKind::LimitExceeded,
                "{format}"
            );
        }
    }
}

/// The distribution of three geometries a host saved, with weights, states
/// with energies and properties, and a temperature.
fn saved_distribution() -> Distribution {
    use kekule::properties::{PropertyKey, PropertyValue};
    use kekule::structure::State;
    use kekule::units::{KELVIN, KILOJOULE_PER_MOLE};

    let mut labelled = State::new()
        .with_energy(Quantity::new(-1.25, KILOJOULE_PER_MOLE))
        .unwrap();
    labelled
        .properties_mut()
        .insert(
            PropertyKey::new("label").unwrap(),
            PropertyValue::String("folded".into()),
        )
        .unwrap();
    let other = State::new()
        .with_energy(Quantity::new(0.1 + 0.2, KILOJOULE_PER_MOLE))
        .unwrap();
    Distribution::new(
        Members::All,
        Partition::new(vec![0, 1, 0], vec![2, 1]).unwrap(),
    )
    .unwrap()
    .with_weights(vec![0.1, 0.7, 1.0 / 3.0])
    .unwrap()
    .with_states(vec![labelled, other])
    .unwrap()
    .with_temperature(Quantity::new(298.15, KELVIN))
    .unwrap()
}

/// A saved distribution of three geometries binds to the file grown to five as
/// a distribution of its first three, and a refresh gives the same.
#[test]
fn a_saved_distribution_binds_to_a_grown_file_as_a_refresh_extends_it() {
    let topology = linear_carbon_topology(12);
    let all = chain(&topology, 5);
    let frames = frames_of(&all);
    let directory = Directory::new();
    let path = directory.0.join("run.xtc");
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc);
    let head = time_ordered(Arc::clone(&topology), frames[..3].iter().cloned());
    write_structure_set(&path, &head, options.clone()).unwrap();

    // The host saves the description of the three geometries and keeps
    // following the file with the set.
    let growing = TrajectoryOpenOptions::default().with_growing(true);
    let mut followed =
        open_indexed_structure_set(&path, Arc::clone(&topology), growing.clone()).unwrap();
    followed
        .set_distribution(Some(saved_distribution()))
        .unwrap();
    let saved = StructureSetDescription::of(&followed);
    let states = saved.distribution().unwrap().states();
    assert_eq!(states.len(), 2);
    assert_eq!(
        states[0].energy().unwrap().value().to_bits(),
        (-1.25_f64).to_bits()
    );
    assert_eq!(states[0].properties().len(), 1);

    let (mut writer, _) =
        append_trajectory(&path, Arc::clone(&topology), KeepFrames::First(3), options).unwrap();
    for frame in all.iter().skip(3) {
        writer
            .write_frame(FrameView::from_geometry(frame).unwrap())
            .unwrap();
    }
    writer.finish().unwrap();
    assert_eq!(followed.refresh().unwrap(), 2);

    // Opened again, the file has five geometries; the saved distribution of
    // three does not cover them, and extended it does.
    let mut reopened = open_indexed_structure_set(&path, Arc::clone(&topology), growing).unwrap();
    assert_eq!(reopened.len(), 5);
    let distribution = saved.distribution().unwrap().clone();
    assert!(reopened
        .set_distribution(Some(distribution.clone()))
        .is_err());
    let extended = distribution.clone().extended_to(5).unwrap();
    reopened.set_distribution(Some(extended.clone())).unwrap();
    assert_eq!(extended.members(), &Members::Subset(vec![0, 1, 2]));
    assert_eq!(extended.weights(), distribution.weights());
    assert_eq!(extended.partition(), distribution.partition());
    assert_eq!(extended.states(), distribution.states());
    assert_eq!(extended.temperature(), distribution.temperature());
    let bits = |weights: Option<&[f64]>| {
        weights
            .unwrap()
            .iter()
            .map(|weight| weight.to_bits())
            .collect::<Vec<_>>()
    };
    assert_eq!(bits(extended.weights()), bits(distribution.weights()));
    // The refresh and the public operation agree.
    assert_eq!(followed.distribution(), Some(&extended));
    assert_eq!(reopened.distribution(), followed.distribution());

    // The description extends alike, and takes the longer set's order.
    let grown = saved
        .clone()
        .extended_to(5, reopened.order().cloned())
        .unwrap();
    assert_eq!(grown.geometry_count(), 5);
    assert_eq!(grown.distribution(), Some(&extended));
    assert_eq!(grown.properties(), saved.properties());
    assert_eq!(grown.order(), reopened.order());
    let same = saved
        .clone()
        .extended_to(3, saved.order().cloned())
        .unwrap();
    assert_eq!(same.order(), saved.order());
    assert_eq!(same.distribution(), saved.distribution());
    // The order must cover the longer set, and the set must not be shorter.
    assert!(saved
        .clone()
        .extended_to(5, saved.order().cloned())
        .is_err());
    assert!(saved.extended_to(2, None).is_err());
    assert!(extended.extended_to(2).is_err());
}

/// The fingerprint of a file's first frames stays while the file grows or
/// its later frames are rewritten, and changes with a frame among them.
#[test]
fn the_fingerprint_of_the_first_frames_ignores_what_follows_them() {
    let topology = linear_carbon_topology(12);
    let all = chain(&topology, 5);
    let frames = frames_of(&all);
    let directory = Directory::new();
    let path = directory.0.join("run.xtc");
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc);
    let head = time_ordered(Arc::clone(&topology), frames[..3].iter().cloned());
    write_structure_set(&path, &head, options.clone()).unwrap();
    let three = open_indexed_trajectory(&path, Arc::clone(&topology), Default::default())
        .unwrap()
        .to_frame_index();
    let used = three.prefix_fingerprint(&path, 3).unwrap();
    assert!(used.same_content(three.fingerprint()));

    // The file grows by two frames.
    let (mut writer, _) = append_trajectory(
        &path,
        Arc::clone(&topology),
        KeepFrames::First(3),
        options.clone(),
    )
    .unwrap();
    for frame in all.iter().skip(3) {
        writer
            .write_frame(FrameView::from_geometry(frame).unwrap())
            .unwrap();
    }
    writer.finish().unwrap();
    let five = open_indexed_trajectory(&path, Arc::clone(&topology), Default::default())
        .unwrap()
        .to_frame_index();
    assert!(five
        .prefix_fingerprint(&path, 3)
        .unwrap()
        .same_content(&used));
    assert!(three
        .prefix_fingerprint(&path, 3)
        .unwrap()
        .same_content(&used));
    assert!(five
        .prefix_fingerprint(&path, 5)
        .unwrap()
        .same_content(five.fingerprint()));
    assert!(!five
        .prefix_fingerprint(&path, 4)
        .unwrap()
        .same_content(&used));
    assert!(used.is_prefix_of_path(&path).unwrap());
    // Frame 0 starts the file, so no frames cover no bytes.
    assert_eq!(five.prefix_fingerprint(&path, 0).unwrap().byte_count(), 0);

    // A restarted run rewrites the last frame: the first three are as they
    // were.
    let (mut writer, _) =
        append_trajectory(&path, Arc::clone(&topology), KeepFrames::First(4), options).unwrap();
    writer
        .write_frame(FrameView::from_geometry(shifted(&all, 0.25).get(4).unwrap()).unwrap())
        .unwrap();
    writer.finish().unwrap();
    let rewritten = open_indexed_trajectory(&path, Arc::clone(&topology), Default::default())
        .unwrap()
        .to_frame_index();
    assert!(!rewritten.fingerprint().same_content(five.fingerprint()));
    assert!(rewritten
        .prefix_fingerprint(&path, 3)
        .unwrap()
        .same_content(&used));

    // A byte inside frame 1 changes.
    let first_frame = five.prefix_fingerprint(&path, 1).unwrap();
    let (start, end) = (five.frame_offsets()[1], five.frame_offsets()[2]);
    let middle = start + (end - start) / 2;
    let mut file = OpenOptions::new()
        .read(true)
        .write(true)
        .open(&path)
        .unwrap();
    let mut byte = [0_u8];
    file.seek(SeekFrom::Start(middle)).unwrap();
    std::io::Read::read_exact(&mut file, &mut byte).unwrap();
    file.seek(SeekFrom::Start(middle)).unwrap();
    file.write_all(&[byte[0] ^ 0x5a]).unwrap();
    drop(file);
    assert!(!five
        .prefix_fingerprint(&path, 3)
        .unwrap()
        .same_content(&used));
    assert!(!used.is_prefix_of_path(&path).unwrap());
    // Frame 0 is unchanged.
    assert!(five
        .prefix_fingerprint(&path, 1)
        .unwrap()
        .same_content(&first_frame));

    let beyond = five.prefix_fingerprint(&path, 6).unwrap_err();
    assert_eq!(beyond.kind(), TrajectoryIoErrorKind::FrameIndexOutOfRange);
    OpenOptions::new()
        .write(true)
        .open(&path)
        .unwrap()
        .set_len(end - 1)
        .unwrap();
    let short = five.prefix_fingerprint(&path, 2).unwrap_err();
    assert_eq!(short.kind(), TrajectoryIoErrorKind::IndexMismatch);
    assert!(five.prefix_fingerprint(&path, 1).is_ok());
}
