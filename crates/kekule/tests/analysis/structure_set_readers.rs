//! Reductions over any geometry reader, such as a set read from a file: equal
//! to the same reductions of the set in memory, for every geometry, a range
//! with a stride, and weights over part of the set.

use kekule::diagnostics::Countable;
use kekule::diagnostics::IndexOutOfRange;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::sync::Arc;

use kekule::alignment::{self, AlignedRmsdOptions, AlignmentError, AlignmentOptions, Weighting};
use kekule::clustering::{self, ClusteringError, Cut, Linkage};
use kekule::diagnostics::Cancellation;
use kekule::math::{Matrix3, Point3, RigidTransform, Vector3};
use kekule::parallel::{Run, ScopedThreads};
use kekule::statistics::{
    ContactOccupancyAccumulator, GroupContactOccupancyAccumulator, RmsfAccumulator,
    SampleObservations, StatisticsError, Weights,
};
use kekule::structure::{
    Geometry, Order, Positions, Structure, StructureSet, StructureSetRead, StructureSetReadError,
    TimeAxis,
};
use kekule::topology::{AtomSelection, InstanceAtomId, Topology};
use kekule::traj::io::trr::{TrrLambdaPolicy, TrrWriteOptions};
use kekule::traj::io::{
    open_indexed_structure_set, write_structure_set, IndexedFile, TrajectoryWriteOptions,
};
use kekule::traj::{FitReference, Processed, ProcessingSteps, TrajectoryFormat, TrajectoryIoError};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};

const ATOMS: usize = 6;

/// Structure `index`: one of two zig-zag shapes by parity, turned and moved
/// a little more with each index.
fn points(index: usize) -> Vec<Point3> {
    let angle = 0.4 * index as f64;
    let (sin, cos) = angle.sin_cos();
    let rotation = Matrix3::from_columns(
        Vector3::new(cos, sin, 0.0),
        Vector3::new(-sin, cos, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    );
    let transform = RigidTransform::new(
        rotation,
        Quantity::new(Vector3::new(angle, -0.5 * angle, 0.2), NANOMETER),
    )
    .unwrap();
    let fold = if index.is_multiple_of(2) { 0.25 } else { -0.25 };
    (0..ATOMS)
        .map(|atom| {
            let wobble = 0.01 * ((index * 7 + atom * 3) % 5) as f64;
            transform.transform_point(Point3::new(
                0.15 * atom as f64 + wobble,
                if atom.is_multiple_of(2) { 0.0 } else { 0.09 },
                if atom >= ATOMS / 2 { fold } else { wobble },
            ))
        })
        .collect()
}

/// Ten timed frames of a six-carbon chain.
fn trajectory() -> StructureSet {
    let topology = kekule::smiles::to_topology(&"C".repeat(ATOMS)).unwrap();
    let frames = (0..10).map(|index| {
        Geometry::new(Positions::new(Quantity::new(points(index), NANOMETER)).unwrap())
    });
    let mut set = StructureSet::from_geometries(topology, frames).unwrap();
    let times = (0..10).map(|index| 2.0 * index as f64).collect();
    let axis = TimeAxis::new(
        Some(Quantity::new(times, PICOSECOND)),
        Some((0..10).collect()),
    )
    .unwrap();
    set.set_order(Some(Order::Time(axis))).unwrap();
    set
}

/// A scratch file removed on drop.
struct Scratch(PathBuf);

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// [`trajectory`] written to a TRR file, read back as a set, and the same
/// set loaded into memory.
fn file_and_memory(scratch: &str) -> (Scratch, StructureSet<IndexedFile>, StructureSet) {
    let directory = std::env::temp_dir().join(format!(
        "kekule-geometry-readers-{scratch}-{}",
        std::process::id()
    ));
    std::fs::create_dir_all(&directory).unwrap();
    let path = directory.join("run.trr");
    let trajectory = trajectory();
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Trr).with_trr_options(
        TrrWriteOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero),
    );
    write_structure_set(&path, &trajectory, options).unwrap();
    let file = open_indexed_structure_set(&path, trajectory.shared_topology(), Default::default())
        .unwrap();
    let memory = StructureSet::from_reader(&file, 0..file.len()).unwrap();
    (Scratch(directory), file, memory)
}

/// Weights of geometries 1, 2, 4, 5, 7, and 8 of ten.
fn partial() -> Weights {
    let mut log_weights = vec![None; 10];
    for (geometry, weight) in [1, 2, 4, 5, 7, 8]
        .into_iter()
        .zip([1.0, 2.0, 3.0, 0.5, 1.5, 4.0])
    {
        log_weights[geometry] = Some(f64::ln(weight));
    }
    Weights::given(log_weights, None).unwrap()
}

fn selection(topology: &Arc<Topology>, atoms: &[usize]) -> AtomSelection {
    let ids = topology.atom_ids();
    AtomSelection::from_atoms(Arc::clone(topology), atoms.iter().map(|atom| ids[*atom])).unwrap()
}

fn pairs(topology: &Topology) -> Vec<(InstanceAtomId, InstanceAtomId)> {
    let ids = topology.atom_ids();
    vec![(ids[0], ids[5]), (ids[1], ids[4]), (ids[0], ids[3])]
}

fn groups(topology: &Arc<Topology>) -> (Vec<AtomSelection>, Vec<(usize, usize)>) {
    let groups = vec![
        selection(topology, &[0, 1]),
        selection(topology, &[4, 5]),
        selection(topology, &[2, 3]),
    ];
    (groups, vec![(0, 1), (2, 0), (1, 2)])
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

const CUTOFF: f64 = 0.45;

/// Checks that every statistic of `observations` of `reader` equals the
/// statistic of `expected`, a set in memory holding just those geometries,
/// with `weights`.
fn statistics_match(
    reader: &dyn StructureSetRead,
    observations: &SampleObservations,
    expected: &StructureSet,
    weights: Option<&Weights>,
) {
    let topology = Arc::clone(reader.topology());
    let cutoff = Quantity::new(CUTOFF, NANOMETER);
    let all = 0..observations.len();
    for atoms in [
        AtomSelection::all(Arc::clone(&topology)),
        selection(&topology, &[1, 4, 5]),
    ] {
        let mut accumulator = RmsfAccumulator::new(&atoms).unwrap();
        accumulator
            .observe_samples(reader, observations, all.clone())
            .unwrap();
        let read = accumulator.finish().unwrap();
        let held = expected.rmsf(&atoms, weights).unwrap();
        assert_eq!(bits(read.values().value()), bits(held.values().value()));
        assert_eq!(read.observation_count(), held.observation_count());
        assert_eq!(read.total_weight().to_bits(), held.total_weight().to_bits());
    }

    let mut contacts =
        ContactOccupancyAccumulator::new(&topology, pairs(&topology), cutoff).unwrap();
    contacts
        .observe_samples(reader, observations, all.clone())
        .unwrap();
    let read = contacts.finish().unwrap();
    let held = expected
        .contact_occupancy(pairs(&topology), cutoff, weights)
        .unwrap();
    assert_eq!(bits(read.contact_weights()), bits(held.contact_weights()));
    assert_eq!(read.observation_count(), held.observation_count());
    assert_eq!(read.total_weight().to_bits(), held.total_weight().to_bits());

    let (group_atoms, group_pairs) = groups(&topology);
    let mut contacts = GroupContactOccupancyAccumulator::new(
        &topology,
        group_atoms.clone(),
        group_pairs.clone(),
        cutoff,
    )
    .unwrap();
    contacts.observe_samples(reader, observations, all).unwrap();
    let read = contacts.finish().unwrap();
    let held = expected
        .group_contact_occupancy(group_atoms, group_pairs, cutoff, weights)
        .unwrap();
    assert_eq!(bits(read.contact_weights()), bits(held.contact_weights()));
    assert_eq!(read.observation_count(), held.observation_count());
    assert_eq!(read.total_weight().to_bits(), held.total_weight().to_bits());
}

#[test]
fn statistics_over_a_reader_equal_the_set_in_memory() {
    let (_scratch, file, memory) = file_and_memory("statistics");
    let atoms = AtomSelection::all(file.shared_topology());
    // Every geometry with weight one, as the set's own statistics.
    for reader in [&file as &dyn StructureSetRead, &memory] {
        let observations = SampleObservations::all(reader);
        assert_eq!(observations.samples(), (0..10).collect::<Vec<_>>());
        assert_eq!(observations.weights(), [1.0; 10]);
        statistics_match(reader, &observations, &memory, None);
        // A range with a stride.
        let listed = (1..9).step_by(3);
        let observations = SampleObservations::of(reader, listed.clone()).unwrap();
        assert_eq!(observations.samples(), [1, 4, 7]);
        statistics_match(reader, &observations, &memory.select(listed).unwrap(), None);
    }
    // The set's own statistics read the same.
    assert_eq!(
        bits(file.rmsf(&atoms, None).unwrap().values().value()),
        bits(memory.rmsf(&atoms, None).unwrap().values().value())
    );
    assert_eq!(
        SampleObservations::of(&file, [2, 10]),
        Err(StatisticsError::Read(
            StructureSetReadError::IndexOutOfRange(IndexOutOfRange::new(
                Countable::Geometry,
                10,
                10
            ))
        ))
    );

    // Every sample of weights over part of the set, with its probability,
    // as the set's own weighted statistics.
    let weights = partial();
    let observations = SampleObservations::weighted(&weights);
    assert_eq!(observations.samples(), [1, 2, 4, 5, 7, 8]);
    for reader in [&file as &dyn StructureSetRead, &memory] {
        statistics_match(reader, &observations, &memory, Some(&weights));
    }
    let listed = (1..9).step_by(3);
    let observations = SampleObservations::weighted_of(&weights, listed.clone()).unwrap();
    // Geometries 1, 4, and 7 are samples.
    assert_eq!(observations.samples(), [1, 4, 7]);
    for reader in [&file as &dyn StructureSetRead, &memory] {
        statistics_match(
            reader,
            &observations,
            &memory.select(listed.clone()).unwrap(),
            Some(&weights.select(listed.clone()).unwrap()),
        );
    }
    assert_eq!(
        bits(file.rmsf(&atoms, Some(&weights)).unwrap().values().value()),
        bits(
            memory
                .rmsf(&atoms, Some(&weights))
                .unwrap()
                .values()
                .value()
        )
    );
    let observations = SampleObservations::weighted_of(&weights, [8, 1, 4]).unwrap();
    assert_eq!(observations.samples(), [8, 1, 4]);
    let probabilities = [4.0, 1.0, 3.0].map(|weight: f64| weight / 8.0);
    for (weight, expected) in observations.weights().iter().zip(probabilities) {
        assert!((weight - expected).abs() < 1.0e-15);
    }
    // Non-samples are left out.
    assert_eq!(
        SampleObservations::weighted_of(&weights, [0, 3, 6])
            .unwrap()
            .len(),
        0
    );
    assert_eq!(
        SampleObservations::weighted_of(&weights, [2, 10]),
        Err(StatisticsError::Read(
            StructureSetReadError::IndexOutOfRange(IndexOutOfRange::new(
                Countable::Geometry,
                10,
                10
            ))
        ))
    );
}

#[test]
fn a_sample_listed_twice_is_observed_twice() {
    let (_scratch, file, memory) = file_and_memory("repeated");
    let listed = [4, 2, 4, 8, 2, 2];
    let selected = memory.select(listed).unwrap();
    for reader in [&file as &dyn StructureSetRead, &memory] {
        let observations = SampleObservations::of(reader, listed).unwrap();
        assert_eq!(observations.samples(), listed);
        assert_eq!(
            bits(observations.weights()),
            bits(SampleObservations::all(&selected).weights())
        );
        statistics_match(reader, &observations, &selected, None);
    }

    // Weighted, each listing of a sample is observed with its share.
    let weights = partial();
    let selected_weights = weights.select(listed).unwrap();
    let observations = SampleObservations::weighted_of(&weights, listed).unwrap();
    assert_eq!(observations.samples(), listed);
    assert_eq!(
        bits(observations.weights()),
        bits(SampleObservations::weighted(&selected_weights).weights())
    );
    for reader in [&file as &dyn StructureSetRead, &memory] {
        statistics_match(reader, &observations, &selected, Some(&selected_weights));
    }
}

#[test]
fn split_observations_merge_and_cancelled_ones_change_nothing() {
    let (_scratch, file, memory) = file_and_memory("merge");
    let atoms = AtomSelection::all(file.shared_topology());
    let observations = SampleObservations::all(&file);
    let mut first = RmsfAccumulator::new(&atoms).unwrap();
    first.observe_samples(&file, &observations, 0..4).unwrap();
    let mut second = RmsfAccumulator::new(&atoms).unwrap();
    second
        .observe_samples(&file, &observations, 4..observations.len())
        .unwrap();
    first.merge(&second).unwrap();
    let merged = first.finish().unwrap();
    let whole = memory.rmsf(&atoms, None).unwrap();
    assert_eq!(merged.observation_count(), whole.observation_count());
    for (merged, whole) in merged.values().value().iter().zip(whole.values().value()) {
        assert!((merged - whole).abs() < 1.0e-12);
    }

    let mut accumulator = RmsfAccumulator::new(&atoms).unwrap();
    accumulator
        .observe_samples(&file, &observations, 0..2)
        .unwrap();
    let cancelled = Cancellation::new(|| true);
    assert!(matches!(
        accumulator.observe_samples_with(
            &file,
            &observations,
            2..6,
            Run::new().with_cancellation(&cancelled)
        ),
        Err(StatisticsError::Cancelled(_))
    ));
    assert_eq!(accumulator.observation_count(), 2);
    assert_eq!(
        accumulator.observe_samples(&file, &observations, 8..11),
        Err(StatisticsError::ObservationsOutOfRange {
            positions: 8..11,
            len: 10
        })
    );
    assert_eq!(accumulator.observation_count(), 2);
    let running = Cancellation::new(|| false);
    accumulator
        .observe_samples_with(
            &file,
            &observations,
            2..10,
            Run::new().with_cancellation(&running),
        )
        .unwrap();
    assert_eq!(accumulator.observation_count(), 10);
}

#[test]
fn rmsd_of_listed_geometries_equals_the_set_in_memory() {
    let (_scratch, file, memory) = file_and_memory("rmsd");
    let topology = file.shared_topology();
    let fit = selection(&topology, &[0, 1, 2, 3]);
    let measured = selection(&topology, &[2, 3, 4, 5]);
    let pair_weights = [1.0, 2.0, 0.5, 3.0];
    let options = AlignmentOptions::default().with_weighting(Weighting::Explicit(&pair_weights));
    let aligned = AlignedRmsdOptions::default()
        .with_fit(options)
        .with_measurement_weighting(Weighting::Explicit(&pair_weights));
    let structure = Structure::new(
        Arc::clone(&topology),
        Positions::new(Quantity::new(points(11), NANOMETER)).unwrap(),
    )
    .unwrap();

    let held = memory.rmsd(3, &fit, options).unwrap();
    let held_aligned = memory.aligned_rmsd(3, &fit, &measured, aligned).unwrap();
    let held_structure = memory.rmsd(&structure, &fit, options).unwrap();
    // A set read from a file measures every geometry as the set in memory.
    assert_eq!(
        bits(file.rmsd(3, &fit, options).unwrap().value()),
        bits(held.value())
    );
    assert_eq!(
        bits(
            file.aligned_rmsd(3, &fit, &measured, aligned)
                .unwrap()
                .value()
        ),
        bits(held_aligned.value())
    );
    assert!(matches!(
        file.rmsd(10, &fit, options),
        Err(AlignmentError::ReferenceOutOfRange { .. })
    ));
    let listed = [9, 0, 3, 3, 6];
    for reader in [&file as &dyn StructureSetRead, &memory] {
        let read = alignment::geometry_rmsd(reader, 0..10, 3, &fit, options).unwrap();
        assert_eq!(bits(read.value()), bits(held.value()));
        let read =
            alignment::geometry_aligned_rmsd(reader, 0..10, 3, &fit, &measured, aligned).unwrap();
        assert_eq!(bits(read.value()), bits(held_aligned.value()));
        let read = alignment::geometry_rmsd(reader, listed, &structure, &fit, options).unwrap();
        assert_eq!(
            bits(read.value()),
            listed.map(|index| held_structure.value()[index].to_bits())
        );
        let cancelled = Cancellation::new(|| true);
        assert!(matches!(
            alignment::geometry_rmsd_with(
                reader,
                listed,
                3,
                &fit,
                options,
                Run::new().with_cancellation(&cancelled)
            ),
            Err(AlignmentError::Cancelled(_))
        ));
        // An endless list stops when cancelled, before or while it is listed.
        assert!(matches!(
            alignment::geometry_aligned_rmsd_with(
                reader,
                0..,
                3,
                &fit,
                &measured,
                aligned,
                Run::new().with_cancellation(&cancelled)
            ),
            Err(AlignmentError::Cancelled(_))
        ));
        let polls = std::sync::atomic::AtomicUsize::new(0);
        let later = Cancellation::new(move || {
            polls.fetch_add(1, std::sync::atomic::Ordering::Relaxed) >= 3
        });
        assert!(matches!(
            alignment::geometry_rmsd_with(
                reader,
                std::iter::repeat(0),
                3,
                &fit,
                options,
                Run::new().with_cancellation(&later)
            ),
            Err(AlignmentError::Cancelled(_))
        ));
        let none = alignment::geometry_rmsd(reader, [0_usize; 0], 3, &fit, options).unwrap();
        assert!(none.value().is_empty());
    }
    assert!(matches!(
        alignment::geometry_rmsd(&file, [0], 10, &fit, options),
        Err(AlignmentError::ReferenceOutOfRange { .. })
    ));
    assert!(matches!(
        alignment::geometry_rmsd(&file, [0, 10], 3, &fit, options),
        Err(AlignmentError::Read(error))
            if error == StructureSetReadError::IndexOutOfRange(IndexOutOfRange::new(Countable::Geometry, 10, 10))
    ));
}

#[test]
fn pairwise_fitted_rmsd_over_a_reader_equals_the_set_in_memory() {
    let (_scratch, file, memory) = file_and_memory("pairwise");
    let topology = file.shared_topology();
    let threads = ScopedThreads::new(NonZeroUsize::new(3).unwrap());
    let weights = partial();
    for atoms in [
        AtomSelection::all(Arc::clone(&topology)),
        selection(&topology, &[1, 2, 4, 5]),
    ] {
        let held = memory
            .pairwise_fitted_rmsd(&atoms, None, Default::default())
            .unwrap();
        assert_eq!(held.member_count(), 10);
        for run in [Run::default(), Run::default().with_executor(&threads)] {
            for reader in [&file as &dyn StructureSetRead, &memory] {
                let read = clustering::pairwise_fitted_rmsd_with(
                    reader,
                    &atoms,
                    None,
                    Default::default(),
                    run,
                )
                .unwrap();
                assert_eq!(read, held);
            }
        }
        // The file's geometries cluster as the set in memory does.
        for linkage in [Linkage::Single, Linkage::Average] {
            assert_eq!(
                file.cluster_hierarchical(&held, None, linkage, Cut::Count(2)),
                memory.cluster_hierarchical(&held, None, linkage, Cut::Count(2))
            );
        }
        assert_eq!(
            file.cluster_k_medoids(&held, None, 2),
            memory.cluster_k_medoids(&held, None, 2)
        );

        // With weights, the file fits and clusters its samples as the set in
        // memory does.
        let held = memory
            .pairwise_fitted_rmsd(&atoms, Some(&weights), Default::default())
            .unwrap();
        assert_eq!(held.member_count(), 6);
        for run in [Run::default(), Run::default().with_executor(&threads)] {
            assert_eq!(
                file.pairwise_fitted_rmsd_with(&atoms, Some(&weights), Default::default(), run)
                    .unwrap(),
                held
            );
            assert_eq!(
                clustering::pairwise_fitted_rmsd_with(
                    &file,
                    &atoms,
                    Some(&weights),
                    Default::default(),
                    run
                )
                .unwrap(),
                held
            );
        }
        for linkage in [Linkage::Single, Linkage::Average] {
            assert_eq!(
                file.cluster_hierarchical(&held, Some(&weights), linkage, Cut::Count(2)),
                memory.cluster_hierarchical(&held, Some(&weights), linkage, Cut::Count(2))
            );
        }
        assert_eq!(
            file.cluster_k_medoids(&held, Some(&weights), 2),
            memory.cluster_k_medoids(&held, Some(&weights), 2)
        );
    }
    let cancelled = Cancellation::new(|| true);
    let run = Run::default().with_cancellation(&cancelled);
    let atoms = AtomSelection::all(Arc::clone(&topology));
    assert!(matches!(
        clustering::pairwise_fitted_rmsd_with(&file, &atoms, None, Default::default(), run),
        Err(ClusteringError::Cancelled(_))
    ));
}

/// Several threads fit the processed geometries of one file at once, each on
/// threads of its own, and all get the matrix of the set in memory.
#[test]
fn pairwise_fitted_rmsd_over_processed_geometries_runs_from_several_threads() {
    let (_scratch, file, memory) = file_and_memory("processed-pairwise");
    let topology = file.shared_topology();
    let fit = selection(&topology, &[0, 1, 2, 3]);
    let steps = ProcessingSteps::default()
        .with_center(fit.clone())
        .with_fit(fit, FitReference::Geometry(2));
    let file = Processed::new(file, steps.clone()).unwrap();
    let memory = Processed::new(memory, steps).unwrap();
    let atoms = selection(&topology, &[1, 2, 4, 5]);
    let held = clustering::pairwise_fitted_rmsd(&memory, &atoms, None, Default::default()).unwrap();
    std::thread::scope(|scope| {
        let runs = (0..4)
            .map(|_| {
                scope.spawn(|| {
                    let threads = ScopedThreads::new(NonZeroUsize::new(3).unwrap());
                    let run = Run::default().with_executor(&threads);
                    clustering::pairwise_fitted_rmsd_with(
                        &file,
                        &atoms,
                        None,
                        Default::default(),
                        run,
                    )
                    .unwrap()
                })
            })
            .collect::<Vec<_>>();
        for run in runs {
            assert_eq!(run.join().unwrap(), held);
        }
    });
}

#[test]
fn a_reduction_over_a_cut_file_names_the_geometry_it_cannot_read() {
    let (scratch, file, _memory) = file_and_memory("cut");
    let starts = file.storage().frame_index().frame_offsets().to_vec();
    std::fs::OpenOptions::new()
        .write(true)
        .open(scratch.0.join("run.trr"))
        .unwrap()
        .set_len(starts[7] + 4)
        .unwrap();
    let atoms = AtomSelection::all(file.shared_topology());
    let observations = SampleObservations::all(&file);
    let mut accumulator = RmsfAccumulator::new(&atoms).unwrap();
    let error = accumulator
        .observe_samples_with(
            &file,
            &observations,
            0..10,
            Run::new().with_cancellation(&Cancellation::new(|| false)),
        )
        .unwrap_err();
    let StatisticsError::Read(StructureSetReadError::Reader { geometry, source }) = &error else {
        panic!("{error:?}");
    };
    assert_eq!(*geometry, 7);
    let source = source.downcast_ref::<TrajectoryIoError>().unwrap();
    assert_eq!(source.frame(), Some(7));
    assert_eq!(error.to_string(), "Cannot read geometry 7.");
}
