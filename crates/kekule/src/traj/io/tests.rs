use super::limits::{index_hard_capacity, next_index_capacity};
use super::*;
use crate::traj::TrajectoryFormat;

#[test]
fn bounded_index_reservation_grows_logarithmically() {
    let entry_count = 100_000;
    let limits = TrajectoryIoLimits::default()
        .with_max_frames(entry_count)
        .with_max_index_entries(entry_count)
        .with_max_index_bytes(entry_count * std::mem::size_of::<u64>());
    let mut offsets = Vec::new();
    let mut growth_events = 0;
    for frame in 0..entry_count {
        assert_eq!(projected_index_limit(offsets.len(), &limits), None);
        let previous_capacity = offsets.capacity();
        reserve_index_for_push(
            &mut offsets,
            &limits,
            TrajectoryFormat::Xyz,
            "capacity-test.xyz",
            frame,
        )
        .unwrap();
        growth_events += usize::from(offsets.capacity() != previous_capacity);
        offsets.push(frame as u64);
    }
    assert_eq!(offsets.len(), entry_count);
    assert!(
        growth_events <= 16,
        "{growth_events} growth events are not logarithmic"
    );
}

#[test]
fn index_hard_capacity_uses_the_smallest_configured_bound() {
    let limits = TrajectoryIoLimits::default()
        .with_max_frames(200)
        .with_max_index_entries(150)
        .with_max_index_bytes(125 * std::mem::size_of::<u64>());
    assert_eq!(index_hard_capacity(&limits), 125);
    assert_eq!(next_index_capacity(64, 64, 125), Some(125));
    assert_eq!(next_index_capacity(125, 125, 125), None);
}

/// Cloning a set over a file reads nothing and shares the frame index;
/// reading a geometry of the clone reads the file.
#[test]
fn a_clone_of_a_set_over_a_file_reads_nothing() {
    use crate::math::Point3;
    use crate::structure::{Geometry, Positions, StructureSet, StructureSetRead};
    use crate::units::{Quantity, NANOMETER};

    let topology = crate::smiles::to_topology("CC").unwrap();
    let mut set = StructureSet::new(std::sync::Arc::clone(&topology));
    for shift in [0.0, 0.1, 0.2] {
        let points = vec![Point3::new(shift, 0.0, 0.0), Point3::new(0.15, shift, 0.0)];
        let positions = Positions::new(Quantity::new(points, NANOMETER)).unwrap();
        set.push(Geometry::new(positions), None).unwrap();
    }
    let path = std::env::temp_dir().join(format!(
        "kekule-clone-reads-nothing-{}.xyz",
        std::process::id()
    ));
    write_structure_set(
        &path,
        &set,
        TrajectoryWriteOptions::new(TrajectoryFormat::Xyz),
    )
    .unwrap();
    let opened =
        open_indexed_structure_set(&path, topology, TrajectoryOpenOptions::default()).unwrap();
    let reads = || super::shared_file::READS.with(std::cell::Cell::get);

    let before = reads();
    let copy = opened.clone();
    assert_eq!(reads(), before);
    assert_eq!(copy.len(), 3);
    assert_eq!(copy.order(), opened.order());
    assert_eq!(
        copy.storage().frame_index().frame_offsets().as_ptr(),
        opened.storage().frame_index().frame_offsets().as_ptr()
    );
    assert_eq!(copy.read(2).unwrap(), opened.read(2).unwrap());
    assert!(reads() > before);
    drop((copy, opened));
    std::fs::remove_file(&path).unwrap();
}

/// A read running when the file is released finishes, and the release
/// returns after it; reads after it fail.
#[test]
fn a_release_waits_for_the_reads_running() {
    use crate::math::Point3;
    use crate::structure::{
        Geometry, Positions, StructureSet, StructureSetRead, StructureSetReadError,
    };
    use crate::traj::TrajectoryIoErrorKind;
    use crate::units::{Quantity, NANOMETER};
    use std::sync::atomic::{AtomicBool, Ordering};
    use std::sync::{mpsc, Arc};
    use std::time::Duration;

    let topology = crate::smiles::to_topology("CC").unwrap();
    let mut set = StructureSet::new(Arc::clone(&topology));
    for shift in [0.0, 0.1, 0.2] {
        let points = vec![Point3::new(shift, 0.0, 0.0), Point3::new(0.15, shift, 0.0)];
        let positions = Positions::new(Quantity::new(points, NANOMETER)).unwrap();
        set.push(Geometry::new(positions), None).unwrap();
    }
    let path =
        std::env::temp_dir().join(format!("kekule-release-waits-{}.xyz", std::process::id()));
    write_structure_set(
        &path,
        &set,
        TrajectoryWriteOptions::new(TrajectoryFormat::Xyz),
    )
    .unwrap();
    let opened =
        open_indexed_structure_set(&path, topology, TrajectoryOpenOptions::default()).unwrap();
    let expected = opened.read(1).unwrap();

    let (started, paused) = mpsc::channel();
    let (resume, resumed_signal) = mpsc::channel::<()>();
    let resumed = Arc::new(AtomicBool::new(false));
    std::thread::scope(|scope| {
        let reading = {
            let (opened, resumed) = (&opened, Arc::clone(&resumed));
            scope.spawn(move || {
                let mut pause = Some((started, resumed_signal));
                super::shared_file::BEFORE_READ.with(|hook| {
                    *hook.borrow_mut() = Some(Box::new(move || {
                        if let Some((started, resumed_signal)) = pause.take() {
                            started.send(()).unwrap();
                            resumed_signal.recv().unwrap();
                            resumed.store(true, Ordering::SeqCst);
                        }
                    }));
                });
                opened.read(1)
            })
        };
        paused.recv_timeout(Duration::from_secs(60)).unwrap();
        let releasing = scope.spawn(|| {
            opened.storage().release();
            resumed.load(Ordering::SeqCst)
        });
        // The release has time to return, if it does not wait.
        std::thread::sleep(Duration::from_millis(200));
        resume.send(()).unwrap();
        assert_eq!(reading.join().unwrap().unwrap(), expected);
        assert!(
            releasing.join().unwrap(),
            "the release returned during the read"
        );
    });
    let Err(StructureSetReadError::Reader { source, .. }) = opened.read(0) else {
        panic!("a read after the release fails");
    };
    let error = source
        .downcast_ref::<crate::traj::TrajectoryIoError>()
        .unwrap();
    assert_eq!(error.kind(), TrajectoryIoErrorKind::Released);
    drop(opened);
    std::fs::remove_file(&path).unwrap();
}
