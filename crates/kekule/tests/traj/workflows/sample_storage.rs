//! Sample sets read from indexed trajectory files: every read is the indexed
//! reader's frame, reads run from many threads at once, a refresh takes the
//! frames a writer appended, and loading gives the set in memory.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;
use std::thread;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::native::SampleSetDescription;
use kekule::structure::{
    Distribution, EnsemblePartition, Members, Order, Positions, Sample, SampleRead,
    SampleReadError, SampleSet, Trajectory, TrajectoryFrame,
};
use kekule::topology::Topology;
use kekule::traj::io::{
    append_trajectory, open_indexed_sample_set, open_indexed_sample_set_from_index,
    open_indexed_trajectory, read_sample_set, write_trajectory, KeepFrames, TrajectoryOpenOptions,
    TrajectoryWriteOptions,
};
use kekule::traj::{
    SeekableTrajectoryReader, TrajectoryFormat, TrajectoryReader, TrajectoryWriter,
};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};

use crate::support::{linear_carbon_topology, topology};

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
            "kekule-sample-storage-{}-{}",
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
fn chain(topology: &Arc<Topology>, count: u64) -> Trajectory {
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
        let mut frame =
            TrajectoryFrame::new(Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap())
                .with_cell(cell);
        frame.set_step(Some(10 * step));
        frame
            .set_time(Some(Quantity::new(0.5 * shift, PICOSECOND)))
            .unwrap();
        frame
    });
    Trajectory::from_items(Arc::clone(topology), frames).unwrap()
}

/// Every frame of the file at `path`, read through its indexed reader.
fn indexed_frames(path: &Path, topology: &Arc<Topology>) -> Vec<TrajectoryFrame> {
    let mut reader =
        open_indexed_trajectory(path, Arc::clone(topology), Default::default()).unwrap();
    let mut buffer = reader.frame_buffer();
    (0..reader.frame_count().unwrap())
        .map(|index| {
            reader.read_frame(index, &mut buffer).unwrap();
            buffer.frame_view().payload().clone()
        })
        .collect()
}

/// Checks that `sample`, sample `index` of `set`, is `frame` with its time
/// and step in the set's order.
fn assert_frame(set: &impl SampleRead, index: usize, sample: &Sample, frame: &TrajectoryFrame) {
    assert_eq!(sample.positions(), frame.positions(), "sample {index}");
    assert_eq!(sample.cell(), frame.cell(), "sample {index}");
    assert_eq!(sample.velocities(), frame.velocities(), "sample {index}");
    assert_eq!(sample.forces(), frame.forces(), "sample {index}");
    assert_eq!(sample.properties(), frame.properties(), "sample {index}");
    let Some(Order::Time(axis)) = set.order() else {
        panic!("a file's samples are in time order");
    };
    assert_eq!(
        axis.time(index)
            .map(|time| time.value_in(PICOSECOND).unwrap()),
        frame.time().map(|time| time.value_in(PICOSECOND).unwrap()),
        "sample {index}"
    );
    assert_eq!(axis.step(index), frame.step(), "sample {index}");
}

#[test]
fn every_read_is_the_indexed_frame_of_every_format() {
    for (path, topology) in fixtures() {
        let frames = indexed_frames(&path, &topology);
        let set =
            open_indexed_sample_set(&path, Arc::clone(&topology), Default::default()).unwrap();
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
            let sample = set.read(index).unwrap();
            assert_frame(&set, index, &sample, &frames[index]);
        }
        assert_eq!(
            set.read(set.len()),
            Err(SampleReadError::IndexOutOfRange {
                index: set.len(),
                len: set.len()
            })
        );

        // Opened again from the kept index, the set reads the same samples.
        let index = set.storage().frame_index().clone();
        let reopened = open_indexed_sample_set_from_index(
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
            open_indexed_sample_set(&path, Arc::clone(&topology), Default::default()).unwrap();
        let loaded = set.load(0..set.len()).unwrap();
        let read = read_sample_set(&path, Arc::clone(&topology), Default::default()).unwrap();
        assert_eq!(loaded.len(), read.len(), "{path:?}");
        assert_eq!(loaded.order(), read.order(), "{path:?}");
        for index in 0..read.len() {
            assert_eq!(
                loaded.get(index).unwrap().payload(),
                read.get(index).unwrap().payload(),
                "{path:?} sample {index}"
            );
        }
        let described = SampleSetDescription::of(&set);
        assert_eq!(described.sample_count(), set.len());
        assert_eq!(described.order(), set.order());
    }

    // Listed samples keep their order values and the distribution's
    // members among them, as selecting them in memory does.
    let topology = linear_carbon_topology(12);
    let directory = Directory::new();
    let path = directory.0.join("run.xtc");
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc);
    write_trajectory(&path, &chain(&topology, 8), options).unwrap();
    let mut set =
        open_indexed_sample_set(&path, Arc::clone(&topology), Default::default()).unwrap();
    let kept = (0..set.len())
        .filter(|index| index % 3 != 1)
        .collect::<Vec<_>>();
    let count = kept.len();
    let members = Members::Subset(kept);
    let weights = (0..count).map(|index| 1.0 + index as f64).collect();
    let distribution = Distribution::new(members, EnsemblePartition::one_state(count))
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
        Err(SampleReadError::IndexOutOfRange { .. })
    ));
    // The set in memory is the set's own: editing it leaves the file alone.
    memory.set_distribution(None).unwrap();
    assert!(set.distribution().is_some());
}

#[test]
fn reads_from_many_threads_at_once_equal_sequential_reads() {
    for (path, topology) in fixtures() {
        let set = open_indexed_sample_set(&path, topology, Default::default()).unwrap();
        let sequential = (0..set.len())
            .map(|index| set.read(index).unwrap())
            .collect::<Vec<_>>();
        let reader: &dyn SampleRead = &set;
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
    let frames = all
        .iter()
        .map(|frame| frame.payload().clone())
        .collect::<Vec<_>>();
    let directory = Directory::new();
    let path = directory.0.join("run.xtc");
    let head = Trajectory::from_items(Arc::clone(&topology), frames[..3].iter().cloned()).unwrap();
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc);
    write_trajectory(&path, &head, options.clone()).unwrap();

    let growing = TrajectoryOpenOptions::default().with_growing(true);
    let mut set = open_indexed_sample_set(&path, Arc::clone(&topology), growing).unwrap();
    assert_eq!(set.len(), 3);
    assert_eq!(set.refresh().unwrap(), 0);
    let distribution = Distribution::new(Members::All, EnsemblePartition::one_state(3)).unwrap();
    set.set_distribution(Some(distribution)).unwrap();
    let before = (0..3)
        .map(|index| set.read(index).unwrap())
        .collect::<Vec<_>>();

    // A restarted run continues the file after its three frames.
    let (mut writer, report) =
        append_trajectory(&path, Arc::clone(&topology), KeepFrames::First(3), options).unwrap();
    assert_eq!(report.kept_frames(), 3);
    for frame in all.iter().skip(3).take(2) {
        writer.write_frame(frame).unwrap();
    }
    writer.flush().unwrap();
    assert_eq!(set.refresh().unwrap(), 2);
    assert_eq!(set.len(), 5);
    for frame in all.iter().skip(5) {
        writer.write_frame(frame).unwrap();
    }
    writer.finish().unwrap();
    assert_eq!(set.refresh().unwrap(), frames.len() - 5);
    assert_eq!(set.len(), frames.len());

    let written = indexed_frames(&path, &topology);
    for (index, frame) in written.iter().enumerate() {
        assert_frame(&set, index, &set.read(index).unwrap(), frame);
    }
    for (index, sample) in before.iter().enumerate() {
        assert_eq!(&set.read(index).unwrap(), sample);
    }
    // The distribution over every sample covers the samples it covered.
    let distribution = set.distribution().unwrap();
    assert_eq!(distribution.members(), &Members::Subset(vec![0, 1, 2]));
    assert_eq!(set.storage().frame_index().frame_count(), frames.len());
    assert_eq!(
        open_indexed_sample_set(&path, topology, Default::default())
            .unwrap()
            .order(),
        set.order()
    );
}

/// Reads on many threads open readers of the grown file, which index more
/// frames than the set holds; a refresh still takes every appended frame.
#[test]
fn refresh_counts_from_the_set_after_concurrent_reads() {
    let topology = linear_carbon_topology(12);
    let all = chain(&topology, 24);
    let directory = Directory::new();
    let path = directory.0.join("run.xtc");
    let head = Trajectory::from_items(
        Arc::clone(&topology),
        all.iter().take(2).map(|frame| frame.payload().clone()),
    )
    .unwrap();
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc);
    write_trajectory(&path, &head, options.clone()).unwrap();
    let growing = TrajectoryOpenOptions::default().with_growing(true);
    let mut set = open_indexed_sample_set(&path, Arc::clone(&topology), growing).unwrap();
    let (mut writer, _) =
        append_trajectory(&path, Arc::clone(&topology), KeepFrames::First(2), options).unwrap();
    for (count, frame) in all.iter().enumerate().skip(2) {
        writer.write_frame(frame).unwrap();
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
    }
    writer.finish().unwrap();
}

#[test]
fn a_failed_refresh_changes_nothing() {
    let (fixture, topology) = fixtures()[3].clone();
    let directory = Directory::new();
    let path = directory.0.join("run.xtc");
    fs::copy(&fixture, &path).unwrap();
    let growing = TrajectoryOpenOptions::default().with_growing(true);
    let mut set = open_indexed_sample_set(&path, Arc::clone(&topology), growing).unwrap();
    let len = set.len();
    let order = set.order().cloned();
    // Another trajectory replaces the indexed bytes.
    let mut bytes = fs::read(&path).unwrap();
    let last = bytes.len() - 5;
    bytes[last] ^= 0x01;
    bytes.extend_from_slice(&fs::read(&fixture).unwrap());
    fs::write(&path, &bytes).unwrap();
    assert!(set.refresh().is_err());
    assert_eq!(set.len(), len);
    assert_eq!(set.order(), order.as_ref());
}

/// A set in memory and the same samples in a file read alike as readers.
#[test]
fn a_set_in_memory_lends_what_a_file_reads() {
    let (path, topology) = fixtures()[3].clone();
    let file = open_indexed_sample_set(&path, topology, Default::default()).unwrap();
    let memory: SampleSet = file.load(0..file.len()).unwrap();
    let readers: [&dyn SampleRead; 2] = [&file, &memory];
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
