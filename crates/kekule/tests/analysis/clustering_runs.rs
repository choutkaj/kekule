//! Clustering runs on executors: every executor and thread count returns the
//! single-threaded matrix, partition, or error; progress counts up to the
//! total; and a cancelled run returns no result.

use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use kekule::alignment::{self, AlignmentError, AlignmentGeometry};
use kekule::clustering::{ClusteringError, ClusteringProgress, ClusteringRun};
use kekule::diagnostics::Cancellation;
use kekule::geometry::Point3;
use kekule::parallel::{Executor, ScopedThreads};
use kekule::structure::{Ensemble, EnsemblePartition, EnsembleSample, Positions};
use kekule::topology::{AtomSelection, Topology};
use kekule::units::{Quantity, ANGSTROM, CANONICAL_LENGTH_UNIT};

use crate::support::{Interleaved, Reversed};

/// Enough atoms that one fit makes a checkpoint and four a work item, so
/// small ensembles split.
const ATOMS: usize = 16_384;

/// Deterministic values in `[0, 1)` from a linear congruential generator.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64
    }
}

fn topology(atoms: usize) -> Arc<Topology> {
    kekule::smiles::to_topology(&"C".repeat(atoms)).unwrap()
}

fn sample(points: Vec<Point3>, weight: f64) -> EnsembleSample {
    EnsembleSample::new(
        Positions::new(Quantity::new(points, ANGSTROM)).unwrap(),
        weight,
    )
    .unwrap()
}

/// `distinct` random structures, then `copies` of them perturbed by a
/// hundredth of an angstrom, interleaved, with distinct weights. With
/// `collinear`, those members lie on a line instead, so no fit to them
/// succeeds.
fn ensemble(distinct: usize, copies: usize, collinear: &[usize]) -> Ensemble {
    let topology = topology(ATOMS);
    let mut random = Lcg(7);
    let bases = (0..distinct)
        .map(|_| {
            (0..ATOMS)
                .map(|_| {
                    Point3::new(
                        20.0 * random.next(),
                        20.0 * random.next(),
                        20.0 * random.next(),
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut samples = Vec::new();
    let mut next_base = 0;
    for index in 0..distinct + copies {
        let copy = index % 4 == 3 && index / 4 < copies;
        let points = if collinear.contains(&index) {
            (0..ATOMS)
                .map(|atom| Point3::new(atom as f64, 0.0, 0.0))
                .collect()
        } else if copy {
            bases[(index * 7) % distinct]
                .iter()
                .map(|point| {
                    Point3::new(
                        point.x + 0.01 * random.next(),
                        point.y + 0.01 * random.next(),
                        point.z,
                    )
                })
                .collect()
        } else {
            bases[next_base].clone()
        };
        next_base += usize::from(!copy);
        samples.push(sample(points, 1.0 + random.next()));
    }
    let partition = EnsemblePartition::state_per_sample(samples.len());
    Ensemble::from_samples(topology, samples, partition).unwrap()
}

/// `members` copies of `bases` random structures of `atoms` atoms, each
/// displaced by up to its own amplitude of zero to two angstroms, with
/// distinct weights.
fn jittered(atoms: usize, bases: usize, members: usize) -> Ensemble {
    let mut random = Lcg(3);
    let bases = (0..bases)
        .map(|_| {
            (0..atoms)
                .map(|_| {
                    Point3::new(
                        5.0 * random.next(),
                        5.0 * random.next(),
                        5.0 * random.next(),
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let samples = (0..members)
        .map(|index| {
            let amplitude = 2.0 * random.next();
            let points = bases[index % bases.len()]
                .iter()
                .map(|point| {
                    Point3::new(
                        point.x + amplitude * random.next(),
                        point.y + amplitude * random.next(),
                        point.z + amplitude * random.next(),
                    )
                })
                .collect();
            sample(points, 1.0 + random.next())
        })
        .collect::<Vec<_>>();
    let partition = EnsemblePartition::state_per_sample(samples.len());
    Ensemble::from_samples(topology(atoms), samples, partition).unwrap()
}

/// Leader clustering as documented, one fit at a time: members by
/// descending weight join the nearest earlier representative within
/// `radius`, the first on ties, or found a cluster.
fn naive_leaders(
    ensemble: &Ensemble,
    atoms: &AtomSelection,
    radius: Quantity<f64>,
) -> EnsemblePartition {
    let radius = radius.value_in(CANONICAL_LENGTH_UNIT).unwrap();
    let view = |member: usize| ensemble.get(member).unwrap().as_model_view();
    let mut assignments = vec![0; ensemble.len()];
    let mut representatives: Vec<usize> = Vec::new();
    for member in ensemble.ranked_by_weight() {
        let mut nearest: Option<(usize, f64)> = None;
        for (cluster, representative) in representatives.iter().enumerate() {
            let rmsd = alignment::fit(
                view(member),
                view(*representative),
                atoms,
                Default::default(),
            )
            .unwrap()
            .rmsd()
            .value_in(CANONICAL_LENGTH_UNIT)
            .unwrap();
            if rmsd <= radius && nearest.is_none_or(|(_, best)| rmsd < best) {
                nearest = Some((cluster, rmsd));
            }
        }
        assignments[member] = match nearest {
            Some((cluster, _)) => cluster,
            None => {
                representatives.push(member);
                representatives.len() - 1
            }
        };
    }
    EnsemblePartition::new(assignments, representatives).unwrap()
}

/// [`ScopedThreads`] that records the most work items one call ran.
struct Recording {
    threads: ScopedThreads,
    most: AtomicUsize,
}

impl Recording {
    fn new(threads: usize) -> Self {
        Self {
            threads: ScopedThreads::new(NonZeroUsize::new(threads).unwrap()),
            most: AtomicUsize::new(0),
        }
    }
}

impl Executor for Recording {
    fn parallelism(&self) -> NonZeroUsize {
        self.threads.parallelism()
    }

    fn run(&self, count: usize, body: &(dyn Fn(usize) + Sync)) {
        self.most.fetch_max(count, Ordering::Relaxed);
        self.threads.run(count, body);
    }
}

/// Every executor these tests compare with the calling thread alone.
fn executors() -> Vec<Box<dyn Executor>> {
    let mut executors: Vec<Box<dyn Executor>> = [2, 3, 8]
        .into_iter()
        .map(|threads| {
            Box::new(ScopedThreads::new(NonZeroUsize::new(threads).unwrap())) as Box<dyn Executor>
        })
        .collect();
    executors.push(Box::new(Reversed(5)));
    executors.push(Box::new(Interleaved));
    executors
}

fn radius() -> Quantity<f64> {
    Quantity::new(0.5, ANGSTROM)
}

/// Progress reports of one run, checked to be nondecreasing with a fixed
/// total.
#[derive(Default)]
struct Reports(Mutex<Vec<ClusteringProgress>>);

impl Reports {
    fn record(&self, progress: ClusteringProgress) {
        self.0.lock().unwrap().push(progress);
    }

    /// Takes the reports so far.
    fn checked(&self, total: usize) -> Vec<usize> {
        let reports = std::mem::take(&mut *self.0.lock().unwrap());
        assert!(reports.iter().all(|report| report.total() == total));
        let completed = reports
            .iter()
            .map(ClusteringProgress::completed)
            .collect::<Vec<_>>();
        assert!(
            completed.windows(2).all(|pair| pair[0] < pair[1]),
            "{completed:?}"
        );
        assert_eq!(completed.first(), Some(&0));
        assert!(completed.iter().all(|done| *done <= total));
        completed
    }
}

#[test]
fn pairwise_rmsd_is_the_serial_matrix_on_every_executor() {
    let ensemble = ensemble(5, 1, &[]);
    let atoms = AtomSelection::all(ensemble.shared_topology());
    let serial = ensemble
        .pairwise_fitted_rmsd(&atoms, Default::default())
        .unwrap();
    for executor in executors() {
        let reports = Reports::default();
        let progress = |progress| reports.record(progress);
        let run = ClusteringRun::default()
            .with_executor(executor.as_ref())
            .with_progress(&progress);
        let matrix = ensemble
            .pairwise_fitted_rmsd_with(&atoms, Default::default(), run)
            .unwrap();
        assert_eq!(matrix, serial);
        assert_eq!(reports.checked(15).last(), Some(&15));
    }

    let recording = Recording::new(4);
    let run = ClusteringRun::default().with_executor(&recording);
    assert_eq!(
        ensemble
            .pairwise_fitted_rmsd_with(&atoms, Default::default(), run)
            .unwrap(),
        serial
    );
    assert_eq!(
        recording.most.into_inner(),
        3,
        "15 pairs split into three work items of at least four"
    );
}

#[test]
fn radius_clustering_is_the_serial_partition_on_every_executor() {
    let ensemble = ensemble(9, 2, &[]);
    let atoms = AtomSelection::all(ensemble.shared_topology());
    let serial = ensemble
        .cluster_by_radius(radius(), &atoms, Default::default())
        .unwrap();
    assert_eq!(serial, naive_leaders(&ensemble, &atoms, radius()));
    assert_eq!(serial.representatives().len(), 9);
    for executor in executors() {
        let reports = Reports::default();
        let progress = |progress| reports.record(progress);
        let run = ClusteringRun::default()
            .with_executor(executor.as_ref())
            .with_progress(&progress);
        let partition = ensemble
            .cluster_by_radius_with(radius(), &atoms, Default::default(), run)
            .unwrap();
        assert_eq!(partition, serial);
        assert_eq!(reports.checked(11).last(), Some(&11));
    }

    let recording = Recording::new(2);
    let run = ClusteringRun::default().with_executor(&recording);
    assert_eq!(
        ensemble
            .cluster_by_radius_with(radius(), &atoms, Default::default(), run)
            .unwrap(),
        serial
    );
    assert_eq!(
        recording.most.into_inner(),
        2,
        "batches of members split in two"
    );
}

#[test]
fn radius_clustering_places_members_as_one_fit_at_a_time_would() {
    let ensemble = jittered(6, 7, 120);
    let atoms = AtomSelection::all(ensemble.shared_topology());
    let mut cluster_counts = Vec::new();
    for value in [0.0, 0.4, 0.8, 1.5, 100.0] {
        let radius = Quantity::new(value, ANGSTROM);
        let expected = naive_leaders(&ensemble, &atoms, radius);
        cluster_counts.push(expected.representatives().len());
        assert_eq!(
            ensemble
                .cluster_by_radius(radius, &atoms, Default::default())
                .unwrap(),
            expected,
            "{value} A"
        );
        for executor in executors() {
            let run = ClusteringRun::default().with_executor(executor.as_ref());
            assert_eq!(
                ensemble
                    .cluster_by_radius_with(radius, &atoms, Default::default(), run)
                    .unwrap(),
                expected,
                "{value} A"
            );
        }
    }
    assert_eq!(cluster_counts.first(), Some(&120));
    assert_eq!(cluster_counts.last(), Some(&1));
    assert!(
        cluster_counts[1..4]
            .iter()
            .all(|count| (2..120).contains(count)),
        "{cluster_counts:?}"
    );
}

#[test]
fn failed_fits_report_the_serial_error_on_every_executor() {
    let ensemble = ensemble(9, 0, &[4, 7]);
    let atoms = AtomSelection::all(ensemble.shared_topology());
    let matrix = ensemble.pairwise_fitted_rmsd(&atoms, Default::default());
    assert_eq!(
        matrix,
        Err(ClusteringError::PairFit {
            first: 0,
            second: 4,
            source: Box::new(AlignmentError::DegenerateGeometry {
                geometry: AlignmentGeometry::Reference,
            }),
        })
    );
    let partition = ensemble.cluster_by_radius(radius(), &atoms, Default::default());
    assert!(matches!(partition, Err(ClusteringError::Alignment(_))));
    for executor in executors() {
        let run = ClusteringRun::default().with_executor(executor.as_ref());
        assert_eq!(
            ensemble.pairwise_fitted_rmsd_with(&atoms, Default::default(), run),
            matrix
        );
        assert_eq!(
            ensemble.cluster_by_radius_with(radius(), &atoms, Default::default(), run),
            partition
        );
    }
}

fn cancelled<T>(result: Result<T, ClusteringError>) -> bool {
    matches!(result, Err(ClusteringError::Cancelled(_)))
}

#[test]
fn cancelled_runs_return_no_result() {
    let ensemble = ensemble(5, 1, &[]);
    let atoms = AtomSelection::all(ensemble.shared_topology());

    // Requested before the run: nothing is fitted.
    let requested = Cancellation::new(|| true);
    let reports = Reports::default();
    let progress = |progress| reports.record(progress);
    let run = ClusteringRun::default()
        .with_cancellation(&requested)
        .with_progress(&progress);
    assert!(cancelled(ensemble.pairwise_fitted_rmsd_with(
        &atoms,
        Default::default(),
        run
    )));
    assert_eq!(reports.checked(15), [0]);
    assert!(cancelled(ensemble.cluster_by_radius_with(
        radius(),
        &atoms,
        Default::default(),
        run
    )));
    assert!(cancelled(
        Ensemble::new(ensemble.shared_topology()).pairwise_fitted_rmsd_with(
            &atoms,
            Default::default(),
            run
        )
    ));

    // Requested partway, from the progress callback.
    for executor in executors() {
        let flag = Arc::new(AtomicBool::new(false));
        let cancellation = Cancellation::from_flag(Arc::clone(&flag));
        let reports = Reports::default();
        let progress = |progress: ClusteringProgress| {
            reports.record(progress);
            if progress.completed() > 0 {
                flag.store(true, Ordering::Relaxed);
            }
        };
        let run = ClusteringRun::default()
            .with_executor(executor.as_ref())
            .with_cancellation(&cancellation)
            .with_progress(&progress);
        assert!(cancelled(ensemble.pairwise_fitted_rmsd_with(
            &atoms,
            Default::default(),
            run
        )));
        let completed = reports.checked(15);
        assert!(*completed.last().unwrap() < 15, "{completed:?}");

        flag.store(false, Ordering::Relaxed);
        let reports = Reports::default();
        let progress = |progress: ClusteringProgress| {
            reports.record(progress);
            if progress.completed() > 0 {
                flag.store(true, Ordering::Relaxed);
            }
        };
        let run = run.with_progress(&progress);
        assert!(cancelled(ensemble.cluster_by_radius_with(
            radius(),
            &atoms,
            Default::default(),
            run
        )));
        let completed = reports.checked(6);
        assert!(*completed.last().unwrap() < 6, "{completed:?}");
    }
}

/// Runs every work item but the last.
struct Skipping;

impl Executor for Skipping {
    fn parallelism(&self) -> NonZeroUsize {
        NonZeroUsize::new(4).unwrap()
    }

    fn run(&self, count: usize, body: &(dyn Fn(usize) + Sync)) {
        (0..count.saturating_sub(1)).for_each(body);
    }
}

#[test]
#[should_panic(expected = "executor skipped work item 2")]
fn an_executor_that_skips_work_items_panics() {
    let ensemble = ensemble(5, 1, &[]);
    let atoms = AtomSelection::all(ensemble.shared_topology());
    let run = ClusteringRun::default().with_executor(&Skipping);
    let _ = ensemble.pairwise_fitted_rmsd_with(&atoms, Default::default(), run);
}

/// 1,000 samples of the first 352 alpha carbons of RCSB 1AKE, each atom
/// displaced by up to a sample's own amplitude of 0.2 to 1.5 angstroms.
fn adenylate_kinase_samples() -> (Ensemble, AtomSelection) {
    let model = kekule::mmcif::parse_str(
        include_str!("../fixtures/mmcif/1AKE.cif"),
        Default::default(),
    )
    .unwrap()
    .interpret(Default::default())
    .unwrap()
    .into_model();
    let alpha_carbons =
        AtomSelection::for_label_atom_names(model.shared_topology(), ["CA"]).unwrap();
    let first = alpha_carbons.atom_ids().take(352).collect::<Vec<_>>();
    let selected = AtomSelection::from_atoms(model.shared_topology(), first).unwrap();
    let model = model.subset(&selected).unwrap();
    let topology = model.shared_topology();
    let atoms = AtomSelection::all(Arc::clone(&topology));
    assert_eq!(atoms.len(), 352);
    let points = model.positions().values().value().to_vec();
    let unit = model.positions().values().unit();
    let angstrom = Quantity::new(1.0, ANGSTROM).value_in(unit).unwrap();
    let mut random = Lcg(11);
    let mut ensemble = Ensemble::new(topology);
    for _ in 0..1000 {
        let amplitude = angstrom * (0.2 + 1.3 * random.next());
        let mut jiggle = || amplitude * (2.0 * random.next() - 1.0);
        let moved = points
            .iter()
            .map(|point| Point3::new(point.x + jiggle(), point.y + jiggle(), point.z + jiggle()))
            .collect::<Vec<_>>();
        let positions = Positions::new(Quantity::new(moved, unit)).unwrap();
        ensemble
            .push_state(
                EnsembleSample::new(positions, 1.0).unwrap(),
                Default::default(),
            )
            .unwrap();
    }
    (ensemble, atoms)
}

#[test]
#[ignore = "timing of 1,000 samples of 352 atoms; run in release with --ignored --nocapture"]
fn clustering_timing_on_a_thousand_samples_of_352_atoms() {
    use kekule::clustering::{Cut, Linkage};
    use std::time::Instant;

    let (ensemble, atoms) = adenylate_kinase_samples();
    let available = std::thread::available_parallelism().map_or(1, NonZeroUsize::get);
    let serial = ensemble
        .pairwise_fitted_rmsd(&atoms, Default::default())
        .unwrap();
    let radius = Quantity::new(1.0, ANGSTROM);
    let leaders = ensemble
        .cluster_by_radius(radius, &atoms, Default::default())
        .unwrap();
    let mut report = format!(
        "1,000 samples of 352 atoms, {} leader clusters within 1 A:",
        leaders.representatives().len()
    );
    for threads in [1, 2, 4, 6, 8, available] {
        let executor = ScopedThreads::new(NonZeroUsize::new(threads).unwrap());
        let run = ClusteringRun::default().with_executor(&executor);
        let start = Instant::now();
        let matrix = ensemble
            .pairwise_fitted_rmsd_with(&atoms, Default::default(), run)
            .unwrap();
        let pairs = start.elapsed();
        assert_eq!(matrix, serial, "{threads} threads");
        let start = Instant::now();
        ensemble
            .cluster_hierarchical(&matrix, Linkage::Average, Cut::Count(10))
            .unwrap();
        let hierarchical = start.elapsed();
        let start = Instant::now();
        let partition = ensemble
            .cluster_by_radius_with(radius, &atoms, Default::default(), run)
            .unwrap();
        let leader = start.elapsed();
        assert_eq!(partition, leaders, "{threads} threads");
        report += &format!(
            "\n  {threads} threads: matrix {pairs:.2?}, hierarchical cut {hierarchical:.2?}, leader {leader:.2?}"
        );
    }
    eprintln!("{report}");
}
