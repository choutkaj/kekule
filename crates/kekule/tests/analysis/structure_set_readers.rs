//! Reductions over any geometry reader, such as a set read from a file: equal
//! to the same reductions of the set in memory, for every geometry, a range
//! with a stride, and a distribution over part of the set.

use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::sync::Arc;

use kekule::alignment::{self, AlignedRmsdOptions, AlignmentError, AlignmentOptions, Weighting};
use kekule::clustering::{self, ClusteringError, ClusteringRun, Cut, Linkage};
use kekule::diagnostics::Cancellation;
use kekule::math::{Matrix3, Point3, RigidTransform, Vector3};
use kekule::parallel::ScopedThreads;
use kekule::statistics::{
    ContactOccupancyAccumulator, GroupContactOccupancyAccumulator, RmsfAccumulator,
    SampleObservations, StatisticsError,
};
use kekule::structure::{
    Distribution, Geometry, Members, Order, Partition, Positions, Structure, StructureSet,
    StructureSetRead, StructureSetReadError, TimeAxis,
};
use kekule::topology::{AtomSelection, InstanceAtomId, Topology};
use kekule::traj::io::trr::{TrrLambdaPolicy, TrrWriteOptions};
use kekule::traj::io::{
    open_indexed_structure_set, write_structure_set, IndexedFile, TrajectoryWriteOptions,
};
use kekule::traj::{FitReference, Processed, ProcessingSteps, TrajectoryFormat};
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
    let memory = file.load(0..file.len()).unwrap();
    (Scratch(directory), file, memory)
}

/// A distribution over geometries 1, 2, 4, 5, 7, and 8, weighted, in two
/// states.
fn partial() -> Distribution {
    Distribution::new(
        Members::Subset(vec![1, 2, 4, 5, 7, 8]),
        Partition::new(vec![0, 1, 0, 1, 1, 0], vec![0, 1]).unwrap(),
    )
    .unwrap()
    .with_weights(vec![1.0, 2.0, 3.0, 0.5, 1.5, 4.0])
    .unwrap()
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
/// statistic of `expected`, a set in memory holding just those geometries.
fn statistics_match(
    reader: &dyn StructureSetRead,
    observations: &SampleObservations,
    expected: &StructureSet,
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
        let held = expected.rmsf(&atoms).unwrap();
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
        .contact_occupancy(pairs(&topology), cutoff)
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
        .group_contact_occupancy(group_atoms, group_pairs, cutoff)
        .unwrap();
    assert_eq!(bits(read.contact_weights()), bits(held.contact_weights()));
    assert_eq!(read.observation_count(), held.observation_count());
    assert_eq!(read.total_weight().to_bits(), held.total_weight().to_bits());
}

#[test]
fn statistics_over_a_reader_equal_the_set_in_memory() {
    let (_scratch, mut file, mut memory) = file_and_memory("statistics");
    for distribution in [None, Some(partial())] {
        file.set_distribution(distribution.clone()).unwrap();
        memory.set_distribution(distribution).unwrap();
        for reader in [&file as &dyn StructureSetRead, &memory] {
            // Every geometry, or every member.
            statistics_match(reader, &SampleObservations::all(reader), &memory);
            // A range with a stride.
            let listed = (1..9).step_by(3);
            let observations = SampleObservations::of(reader, listed.clone()).unwrap();
            statistics_match(reader, &observations, &memory.select(listed).unwrap());
        }
        // The set's own statistics read the same.
        let atoms = AtomSelection::all(file.shared_topology());
        assert_eq!(
            bits(file.rmsf(&atoms).unwrap().values().value()),
            bits(memory.rmsf(&atoms).unwrap().values().value())
        );
    }
    let observations = SampleObservations::of(&file, [8, 1, 4]).unwrap();
    assert_eq!(observations.samples(), [8, 1, 4]);
    let probabilities = [4.0, 1.0, 3.0].map(|weight: f64| weight / 8.0);
    for (weight, expected) in observations.weights().iter().zip(probabilities) {
        assert!((weight - expected).abs() < 1.0e-15);
    }
    // Non-members are left out.
    assert_eq!(SampleObservations::of(&file, [0, 3, 6]).unwrap().len(), 0);
    assert_eq!(
        SampleObservations::of(&file, [2, 10]),
        Err(StatisticsError::Read(
            StructureSetReadError::IndexOutOfRange { index: 10, len: 10 }
        ))
    );
}

#[test]
fn a_sample_listed_twice_is_observed_twice() {
    let (_scratch, mut file, mut memory) = file_and_memory("repeated");
    let listed = [4, 2, 4, 8, 2, 2];
    for distribution in [None, Some(partial())] {
        file.set_distribution(distribution.clone()).unwrap();
        memory.set_distribution(distribution).unwrap();
        let selected = memory.select(listed).unwrap();
        for reader in [&file as &dyn StructureSetRead, &memory] {
            let observations = SampleObservations::of(reader, listed).unwrap();
            assert_eq!(observations.samples(), listed);
            assert_eq!(
                bits(observations.weights()),
                bits(SampleObservations::all(&selected).weights())
            );
            statistics_match(reader, &observations, &selected);
        }
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
    let whole = memory.rmsf(&atoms).unwrap();
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
        accumulator.observe_samples_cancellable(&file, &observations, 2..6, &cancelled),
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
        .observe_samples_cancellable(&file, &observations, 2..10, &running)
        .unwrap();
    assert_eq!(accumulator.observation_count(), 10);
}

#[test]
fn rmsd_of_listed_geometries_equals_the_set_in_memory() {
    let (_scratch, mut file, mut memory) = file_and_memory("rmsd");
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
    // A distribution does not restrict which geometries are measured.
    file.set_distribution(Some(partial())).unwrap();
    memory.set_distribution(Some(partial())).unwrap();

    let held = memory.rmsd(3, &fit, options).unwrap();
    let held_aligned = memory.aligned_rmsd(3, &fit, &measured, aligned).unwrap();
    let held_structure = memory.rmsd(&structure, &fit, options).unwrap();
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
            alignment::geometry_rmsd_cancellable(reader, listed, 3, &fit, options, &cancelled),
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
            if *error == StructureSetReadError::IndexOutOfRange { index: 10, len: 10 }
    ));
}

#[test]
fn pairwise_fitted_rmsd_over_a_reader_equals_the_set_in_memory() {
    let (_scratch, mut file, mut memory) = file_and_memory("pairwise");
    let topology = file.shared_topology();
    let threads = ScopedThreads::new(NonZeroUsize::new(3).unwrap());
    for distribution in [None, Some(partial())] {
        file.set_distribution(distribution.clone()).unwrap();
        memory.set_distribution(distribution).unwrap();
        for atoms in [
            AtomSelection::all(Arc::clone(&topology)),
            selection(&topology, &[1, 2, 4, 5]),
        ] {
            let held = memory
                .pairwise_fitted_rmsd(&atoms, Default::default())
                .unwrap();
            for run in [
                ClusteringRun::default(),
                ClusteringRun::default().with_executor(&threads),
            ] {
                for reader in [&file as &dyn StructureSetRead, &memory] {
                    let read =
                        clustering::pairwise_fitted_rmsd(reader, &atoms, Default::default(), run)
                            .unwrap();
                    assert_eq!(read, held);
                }
            }
            // The file's geometries cluster as the set in memory does.
            for linkage in [Linkage::Single, Linkage::Average] {
                assert_eq!(
                    file.cluster_hierarchical(&held, linkage, Cut::Count(2)),
                    memory.cluster_hierarchical(&held, linkage, Cut::Count(2))
                );
            }
            assert_eq!(
                file.cluster_k_medoids(&held, 2),
                memory.cluster_k_medoids(&held, 2)
            );
        }
    }
    let cancelled = Cancellation::new(|| true);
    let run = ClusteringRun::default().with_cancellation(&cancelled);
    let atoms = AtomSelection::all(Arc::clone(&topology));
    assert!(matches!(
        clustering::pairwise_fitted_rmsd(&file, &atoms, Default::default(), run),
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
    let held = clustering::pairwise_fitted_rmsd(
        &memory,
        &atoms,
        Default::default(),
        ClusteringRun::default(),
    )
    .unwrap();
    std::thread::scope(|scope| {
        let runs = (0..4)
            .map(|_| {
                scope.spawn(|| {
                    let threads = ScopedThreads::new(NonZeroUsize::new(3).unwrap());
                    let run = ClusteringRun::default().with_executor(&threads);
                    clustering::pairwise_fitted_rmsd(&file, &atoms, Default::default(), run)
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
        .observe_samples_cancellable(&file, &observations, 0..10, &Cancellation::new(|| false))
        .unwrap_err();
    let StatisticsError::Read(StructureSetReadError::Trajectory { geometry, source }) = &error
    else {
        panic!("{error:?}");
    };
    assert_eq!(*geometry, 7);
    assert_eq!(source.frame(), Some(7));
    assert_eq!(error.to_string(), "Cannot read geometry 7.");
}
