//! Superposition, RMSD, weighted statistics, clustering, and Boltzmann
//! populations of sample sets: identical to those of the ensembles,
//! trajectories, and series they convert from, and restricted to the members
//! of a distribution over part of a set.

use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, OnceLock};

use kekule::alignment::{
    AlignedRmsdOptions, AlignmentError, AlignmentGeometry, AlignmentOptions, Weighting,
};
use kekule::clustering::{
    ClusterOrder, ClusteringError, ClusteringOptions, ClusteringProgress, ClusteringRun, Cut,
    Linkage,
};
use kekule::diagnostics::Cancellation;
use kekule::geometry::{Matrix3, Point3, RigidTransform, Vector3};
use kekule::parallel::ScopedThreads;
use kekule::statistics::{
    ContactOccupancyResult, GroupContactOccupancyResult, RmsfResult, StatisticsError,
};
use kekule::structure::{
    Distribution, Ensemble, EnsemblePartition, EnsembleSample, EnsembleState, Members, Model,
    Positions, Requirement, Sample, SampleSet, SampleSetError, Series, SeriesAxis, SeriesItem,
    Trajectory, TrajectoryFrame, Velocities,
};
use kekule::topology::{AtomSelection, InstanceAtomId, Topology};
use kekule::units::{
    Quantity, CANONICAL_VELOCITY_UNIT, DEGREE, KELVIN, KILOJOULE_PER_MOLE, NANOMETER, PICOSECOND,
};

const ATOMS: usize = 6;

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

/// One carbon chain, shared by every fixture so that their structures and
/// references share one layout.
fn topology() -> Arc<Topology> {
    static TOPOLOGY: OnceLock<Arc<Topology>> = OnceLock::new();
    Arc::clone(TOPOLOGY.get_or_init(another_topology))
}

/// A carbon chain of another layout than [`topology`]'s.
fn another_topology() -> Arc<Topology> {
    kekule::smiles::to_topology(&"C".repeat(ATOMS)).unwrap()
}

/// A rotation by `angle` radians about the body diagonal (Rodrigues'
/// formula), then a translation.
fn motion(angle: f64) -> RigidTransform {
    let (sin, cos) = angle.sin_cos();
    let axis = 1.0 / 3.0_f64.sqrt();
    let outer = (1.0 - cos) * axis * axis;
    let (diagonal, plus, minus) = (cos + outer, outer + sin * axis, outer - sin * axis);
    let rotation = Matrix3::from_columns(
        Vector3::new(diagonal, plus, minus),
        Vector3::new(minus, diagonal, plus),
        Vector3::new(plus, minus, diagonal),
    );
    RigidTransform::new(
        rotation,
        Quantity::new(Vector3::new(angle, -2.0 * angle, 0.5), NANOMETER),
    )
    .unwrap()
}

/// Structure `index`: one of two zig-zag shapes by parity, jittered by a
/// few hundredths of a nanometre and moved rigidly, so that structures of
/// one shape cluster together only after fitting.
fn points(index: usize) -> Vec<Point3> {
    let mut random = Lcg(17 + index as u64);
    let fold = if index.is_multiple_of(2) { 0.25 } else { -0.25 };
    let transform = motion(0.7 * index as f64);
    (0..ATOMS)
        .map(|atom| {
            let base = Point3::new(
                0.15 * atom as f64,
                if atom.is_multiple_of(2) { 0.0 } else { 0.09 },
                if atom >= ATOMS / 2 { fold } else { 0.0 },
            );
            let jitter = Vector3::new(random.next(), random.next(), random.next()) * 0.02;
            transform.transform_point(base + jitter)
        })
        .collect()
}

fn positions(points: Vec<Point3>) -> Positions {
    Positions::new(Quantity::new(points, NANOMETER)).unwrap()
}

/// Atoms on a line, to which no fit succeeds.
fn collinear() -> Vec<Point3> {
    (0..ATOMS)
        .map(|atom| Point3::new(0.1 * atom as f64, 0.0, 0.0))
        .collect()
}

const WEIGHTS: [f64; 8] = [1.0, 3.0, 0.5, 2.0, 1.5, 4.0, 0.25, 1.0];
/// The state of each sample of [`ensemble`], its first sample representing it.
const STATES: [usize; 8] = [0, 0, 1, 1, 2, 2, 2, 3];
const ENERGIES: [f64; 4] = [1.5, 0.0, 4.0, 0.7];

/// Eight weighted samples in four states, with state energies or without.
fn ensemble(energies: bool) -> Ensemble {
    let mut ensemble = Ensemble::new(topology());
    for (index, (weight, state)) in WEIGHTS.iter().zip(STATES).enumerate() {
        let sample = EnsembleSample::new(positions(points(index)), *weight).unwrap();
        if state == ensemble.state_count() {
            let mut payload = EnsembleState::new();
            if energies {
                payload = payload
                    .with_energy(Quantity::new(ENERGIES[state], KILOJOULE_PER_MOLE))
                    .unwrap();
            }
            ensemble.push_state(sample, payload).unwrap();
        } else {
            ensemble.push_sample(sample, state).unwrap();
        }
    }
    ensemble
}

/// Eight timed frames with velocities.
fn trajectory() -> Trajectory {
    let frames = (0..8).map(|index| {
        let mut frame = TrajectoryFrame::new(positions(points(index)));
        frame
            .set_time(Some(Quantity::new(2.0 * index as f64, PICOSECOND)))
            .unwrap();
        let velocities = (0..ATOMS)
            .map(|atom| Vector3::new(atom as f64, 1.0, -(index as f64)))
            .collect::<Vec<_>>();
        frame
            .set_velocities(Some(
                Velocities::new(Quantity::new(velocities, CANONICAL_VELOCITY_UNIT)).unwrap(),
            ))
            .unwrap();
        frame
    });
    Trajectory::from_items(topology(), frames).unwrap()
}

/// Six items of a scan of one angle.
fn series() -> Series {
    let items = (0..6).map(|index| {
        SeriesItem::at_coordinate(
            positions(points(index)),
            Quantity::new(30.0 * index as f64, DEGREE),
        )
        .unwrap()
    });
    Series::from_items(topology(), SeriesAxis::Coordinate(DEGREE), items).unwrap()
}

fn all(set_topology: Arc<Topology>) -> AtomSelection {
    AtomSelection::all(set_topology)
}

fn selection(topology: &Arc<Topology>, atoms: &[usize]) -> AtomSelection {
    let ids = topology.atom_ids();
    AtomSelection::from_atoms(Arc::clone(topology), atoms.iter().map(|atom| ids[*atom])).unwrap()
}

fn samples(set: &SampleSet) -> Vec<Sample> {
    set.iter().map(|sample| sample.payload().clone()).collect()
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn same_rmsf(new: &RmsfResult, old: &RmsfResult) {
    assert_eq!(bits(new.values().value()), bits(old.values().value()));
    assert_eq!(new.values().unit(), old.values().unit());
    assert_eq!(new.selection(), old.selection());
    assert_eq!(new.observation_count(), old.observation_count());
    assert_eq!(new.total_weight().to_bits(), old.total_weight().to_bits());
}

fn same_contacts(new: &ContactOccupancyResult, old: &ContactOccupancyResult) {
    assert_eq!(bits(new.contact_weights()), bits(old.contact_weights()));
    assert_eq!(new.pairs(), old.pairs());
    assert_eq!(new.cutoff(), old.cutoff());
    assert_eq!(new.observation_count(), old.observation_count());
    assert_eq!(new.total_weight().to_bits(), old.total_weight().to_bits());
}

fn same_group_contacts(new: &GroupContactOccupancyResult, old: &GroupContactOccupancyResult) {
    assert_eq!(bits(new.contact_weights()), bits(old.contact_weights()));
    assert_eq!(new.groups(), old.groups());
    assert_eq!(new.pairs(), old.pairs());
    assert_eq!(new.cutoff(), old.cutoff());
    assert_eq!(new.observation_count(), old.observation_count());
    assert_eq!(new.total_weight().to_bits(), old.total_weight().to_bits());
}

fn contact_pairs(topology: &Topology) -> Vec<(InstanceAtomId, InstanceAtomId)> {
    let ids = topology.atom_ids();
    vec![(ids[0], ids[5]), (ids[1], ids[4]), (ids[0], ids[1])]
}

/// The groups and group pairs of the group contact checks.
fn contact_groups(topology: &Arc<Topology>) -> (Vec<AtomSelection>, Vec<(usize, usize)>) {
    let groups = vec![
        selection(topology, &[0, 1]),
        selection(topology, &[4, 5]),
        selection(topology, &[2, 3]),
    ];
    (groups, vec![(0, 1), (2, 0), (1, 2)])
}

const CUTOFF: f64 = 0.62;

/// Checks every statistic of `set` against `ensemble`'s.
fn statistics_match_ensemble(set: &SampleSet, ensemble: &Ensemble) {
    let topology = ensemble.shared_topology();
    let cutoff = Quantity::new(CUTOFF, NANOMETER);
    for atoms in [all(Arc::clone(&topology)), selection(&topology, &[1, 4, 5])] {
        same_rmsf(&set.rmsf(&atoms).unwrap(), &ensemble.rmsf(&atoms).unwrap());
    }
    same_contacts(
        &set.contact_occupancy(contact_pairs(&topology), cutoff)
            .unwrap(),
        &ensemble
            .contact_occupancy(contact_pairs(&topology), cutoff)
            .unwrap(),
    );
    let (groups, pairs) = contact_groups(&topology);
    same_group_contacts(
        &set.group_contact_occupancy(groups.clone(), pairs.clone(), cutoff)
            .unwrap(),
        &ensemble
            .group_contact_occupancy(groups, pairs, cutoff)
            .unwrap(),
    );
}

/// Every alignment of a collection, as `(rmsd, aligned rmsd, superposition
/// report)`.
type AlignmentResults = (
    Result<Quantity<Vec<f64>>, AlignmentError>,
    Result<Quantity<Vec<f64>>, AlignmentError>,
    Result<kekule::alignment::SuperpositionReport, AlignmentError>,
);

/// Measures and then superposes any collection with one set of atoms,
/// weights, and options, against `reference` or else its item `index`.
macro_rules! alignments {
    ($collection:expr, $reference:expr, $index:expr) => {{
        let collection = &mut $collection;
        let topology = collection.shared_topology();
        let fit = selection(&topology, &[0, 1, 2, 3]);
        let measured = selection(&topology, &[2, 3, 4, 5]);
        let pair_weights = [1.0, 2.0, 0.5, 3.0];
        let options =
            AlignmentOptions::default().with_weighting(Weighting::Explicit(&pair_weights));
        let aligned = AlignedRmsdOptions::default()
            .with_fit(options)
            .with_measurement_weighting(Weighting::Explicit(&pair_weights));
        let reference: Option<&Model> = $reference;
        let results: AlignmentResults = match reference {
            Some(model) => (
                collection.rmsd(model, &fit, options),
                collection.aligned_rmsd(model, &fit, &measured, aligned),
                collection.superpose(model, &fit, options),
            ),
            None => (
                collection.rmsd($index, &fit, options),
                collection.aligned_rmsd($index, &fit, &measured, aligned),
                collection.superpose($index, &fit, options),
            ),
        };
        results
    }};
}

fn same_alignments(new: &AlignmentResults, old: &AlignmentResults) {
    let values = |result: &Result<Quantity<Vec<f64>>, AlignmentError>| {
        result
            .as_ref()
            .map(|values| bits(values.value()))
            .map_err(Clone::clone)
    };
    assert_eq!(values(&new.0), values(&old.0));
    assert_eq!(values(&new.1), values(&old.1));
    assert_eq!(new.2, old.2);
}

#[test]
fn sample_sets_align_as_ensembles_trajectories_and_series_do() {
    let model = Model::new(topology(), positions(points(11))).unwrap();
    // A sample of the set, an independent view, and a reference out of range.
    for (reference, index) in [(None, 3), (Some(&model), 0), (None, 9)] {
        let mut old = ensemble(true);
        let mut set = SampleSet::from_ensemble(old.clone());
        let new = alignments!(set, reference, index);
        same_alignments(&new, &alignments!(old, reference, index));
        assert_eq!(new.2.is_ok(), index < 8);
        assert_eq!(samples(&set), samples(&SampleSet::from_ensemble(old)));

        // Velocities rotate with the frames.
        let mut old = trajectory();
        let mut set = SampleSet::from_trajectory(old.clone()).unwrap();
        let new = alignments!(set, reference, index);
        same_alignments(&new, &alignments!(old, reference, index));
        let old = SampleSet::from_trajectory(old).unwrap();
        assert_eq!(samples(&set), samples(&old));
        assert_eq!(set.order(), old.order());

        let mut old = series();
        let mut set = SampleSet::from_series(old.clone());
        let new = alignments!(set, reference, index);
        same_alignments(&new, &alignments!(old, reference, index));
        let old = SampleSet::from_series(old);
        assert_eq!(samples(&set), samples(&old));
        assert_eq!(set.order(), old.order());
    }
}

#[test]
fn sample_set_alignment_covers_non_members_and_keeps_the_distribution() {
    // A failing sample outside the distribution still fails the fit, and a
    // fit keeps the order and the distribution.
    let mut set = SampleSet::from_trajectory(trajectory()).unwrap();
    let distribution = Distribution::new(
        Members::Subset(vec![1, 4]),
        EnsemblePartition::state_per_sample(2),
    )
    .unwrap()
    .with_weights(vec![2.0, 1.0])
    .unwrap();
    set.set_distribution(Some(distribution.clone())).unwrap();
    let atoms = all(set.shared_topology());
    let order = set.order().cloned();
    let report = set.superpose(0, &atoms, Default::default()).unwrap();
    assert_eq!(report.len(), 8);
    assert_eq!(set.distribution(), Some(&distribution));
    assert_eq!(set.order().cloned(), order);
    let rmsd = set.rmsd(0, &atoms, Default::default()).unwrap();
    assert_eq!(rmsd.value().len(), 8);

    let mut broken = set.clone();
    broken
        .get_mut(6)
        .unwrap()
        .set_positions(Quantity::new(collinear(), NANOMETER))
        .unwrap();
    let before = samples(&broken);
    assert_eq!(
        broken.superpose(6, &atoms, Default::default()),
        Err(AlignmentError::Item {
            index: 0,
            source: Box::new(AlignmentError::DegenerateGeometry {
                geometry: AlignmentGeometry::Reference,
            }),
        })
    );
    assert_eq!(samples(&broken), before);
}

#[test]
fn sample_set_statistics_equal_ensemble_and_trajectory_statistics() {
    for energies in [true, false] {
        let old = ensemble(energies);
        statistics_match_ensemble(&SampleSet::from_ensemble(old.clone()), &old);
    }

    // Without a distribution, every frame counts with weight one.
    let old = trajectory();
    let set = SampleSet::from_trajectory(old.clone()).unwrap();
    let topology = old.shared_topology();
    let cutoff = Quantity::new(CUTOFF, NANOMETER);
    let atoms = all(Arc::clone(&topology));
    let rmsf = set.rmsf(&atoms).unwrap();
    same_rmsf(&rmsf, &old.rmsf(&atoms).unwrap());
    assert_eq!(rmsf.total_weight(), 8.0);
    same_contacts(
        &set.contact_occupancy(contact_pairs(&topology), cutoff)
            .unwrap(),
        &old.contact_occupancy(contact_pairs(&topology), cutoff)
            .unwrap(),
    );
    let (groups, pairs) = contact_groups(&topology);
    same_group_contacts(
        &set.group_contact_occupancy(groups.clone(), pairs.clone(), cutoff)
            .unwrap(),
        &old.group_contact_occupancy(groups, pairs, cutoff).unwrap(),
    );

    // Errors are the old ones too.
    let empty = SampleSet::new(self::topology());
    assert_eq!(
        empty.rmsf(&all(empty.shared_topology())).unwrap_err(),
        StatisticsError::NoObservations
    );
    assert!(matches!(
        set.rmsf(&all(another_topology())),
        Err(StatisticsError::Selection(_))
    ));
}

/// The members of `set`'s distribution as an ensemble of their own, with
/// their weights and states, and without the velocities an ensemble cannot
/// hold.
fn members_ensemble(set: &SampleSet) -> Ensemble {
    let members = set
        .weighted_samples()
        .map(|(sample, _)| sample)
        .collect::<Vec<_>>();
    let mut members = set.select(members).unwrap();
    members.set_order(None).unwrap();
    for index in 0..members.len() {
        members.get_mut(index).unwrap().take_velocities();
    }
    members.into_ensemble().unwrap()
}

/// A time-ordered set of eight frames whose distribution covers frames 2,
/// 3, 5, 6, and 7, weighted, in two states with energies.
fn equilibrated() -> SampleSet {
    let mut set = SampleSet::from_trajectory(trajectory()).unwrap();
    let states = [0.0, 2.5].map(|energy| {
        EnsembleState::new()
            .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
            .unwrap()
    });
    let distribution = Distribution::new(
        Members::Subset(vec![2, 3, 5, 6, 7]),
        EnsemblePartition::new(vec![0, 1, 0, 1, 1], vec![0, 1]).unwrap(),
    )
    .unwrap()
    .with_weights(vec![1.0, 2.0, 3.0, 0.5, 1.5])
    .unwrap()
    .with_states(states.to_vec())
    .unwrap();
    set.set_distribution(Some(distribution)).unwrap();
    set
}

#[test]
fn statistics_over_part_of_a_time_ordered_set_observe_its_members_alone() {
    let set = equilibrated();
    let members = members_ensemble(&set);
    statistics_match_ensemble(&set, &members);
    let rmsf = set.rmsf(&all(set.shared_topology())).unwrap();
    assert_eq!(rmsf.observation_count(), 5);

    // Moving a non-member changes nothing.
    let mut moved = set.clone();
    moved
        .get_mut(0)
        .unwrap()
        .set_positions(Quantity::new(collinear(), NANOMETER))
        .unwrap();
    moved
        .get_mut(4)
        .unwrap()
        .apply_rigid_transform(&motion(2.0))
        .unwrap();
    statistics_match_ensemble(&moved, &members);

    // A distribution without weights observes each member with probability
    // one in five.
    let mut equal = set.clone();
    equal.distribution_mut().unwrap().set_weights(None).unwrap();
    let rmsf = equal.rmsf(&all(set.shared_topology())).unwrap();
    assert!((rmsf.total_weight() - 1.0).abs() < 1.0e-15);
    let cutoff = Quantity::new(CUTOFF, NANOMETER);
    let topology = set.shared_topology();
    let contacts = equal
        .contact_occupancy(contact_pairs(&topology), cutoff)
        .unwrap();
    let frames = equal.select([2, 3, 5, 6, 7]).unwrap();
    let mut frames = frames;
    frames.set_distribution(None).unwrap();
    let counts = frames
        .contact_occupancy(contact_pairs(&topology), cutoff)
        .unwrap();
    for (probability, count) in contacts
        .contact_weights()
        .iter()
        .zip(counts.contact_weights())
    {
        assert!((probability - count / 5.0).abs() < 1.0e-15);
    }
}

/// Every clustering of `set` equals `ensemble`'s, which holds `set`'s
/// members in member order.
fn clustering_matches_ensemble(set: &SampleSet, ensemble: &Ensemble, energies: bool) {
    let atoms = all(ensemble.shared_topology());
    // The two shapes are two clusters of one radius.
    let shapes = set
        .cluster_by_radius(Quantity::new(0.08, NANOMETER), &atoms, Default::default())
        .unwrap();
    assert_eq!(shapes.state_count(), 2);
    let matrix = set
        .pairwise_fitted_rmsd(&atoms, Default::default())
        .unwrap();
    assert_eq!(
        matrix,
        ensemble
            .pairwise_fitted_rmsd(&atoms, Default::default())
            .unwrap()
    );
    let mut orders = vec![ClusterOrder::Weight];
    if energies {
        orders.push(ClusterOrder::Energy);
    }
    for order in orders {
        for radius in [0.0, 0.05, 0.08, 0.2, 10.0] {
            let options = ClusteringOptions::default().with_order(order);
            let radius = Quantity::new(radius, NANOMETER);
            assert_eq!(
                set.cluster_by_radius(radius, &atoms, options).unwrap(),
                ensemble.cluster_by_radius(radius, &atoms, options).unwrap(),
            );
        }
    }
    for linkage in [Linkage::Single, Linkage::Complete, Linkage::Average] {
        for count in 1..=matrix.member_count() {
            let cut = Cut::Count(count);
            assert_eq!(
                set.cluster_hierarchical(&matrix, linkage, cut).unwrap(),
                ensemble
                    .cluster_hierarchical(&matrix, linkage, cut)
                    .unwrap(),
            );
        }
        let cut = Cut::Height(Quantity::new(0.06, NANOMETER));
        assert_eq!(
            set.cluster_hierarchical(&matrix, linkage, cut).unwrap(),
            ensemble
                .cluster_hierarchical(&matrix, linkage, cut)
                .unwrap(),
        );
    }
    for clusters in 1..=matrix.member_count() {
        assert_eq!(
            set.cluster_k_medoids(&matrix, clusters).unwrap(),
            ensemble.cluster_k_medoids(&matrix, clusters).unwrap(),
        );
    }
}

#[test]
fn sample_sets_cluster_as_ensembles_do() {
    for energies in [true, false] {
        let old = ensemble(energies);
        let set = SampleSet::from_ensemble(old.clone());
        clustering_matches_ensemble(&set, &old, energies);
    }
    // Without energies, ordering by energy is refused by both.
    let old = ensemble(false);
    let set = SampleSet::from_ensemble(old.clone());
    let atoms = all(old.shared_topology());
    let by_energy = ClusteringOptions::default().with_order(ClusterOrder::Energy);
    let radius = Quantity::new(0.1, NANOMETER);
    assert_eq!(
        set.cluster_by_radius(radius, &atoms, by_energy),
        Err(ClusteringError::SampleSet(Box::new(
            SampleSetError::Missing(Requirement::StateEnergies)
        )))
    );
    assert!(matches!(
        old.cluster_by_radius(radius, &atoms, by_energy),
        Err(ClusteringError::Ensemble(_))
    ));
    let matrix = set
        .pairwise_fitted_rmsd(&atoms, Default::default())
        .unwrap();
    assert_eq!(
        set.cluster_k_medoids(&matrix, 0),
        old.cluster_k_medoids(&matrix, 0)
    );
    let smaller = SampleSet::from_ensemble(old.select(0..3).unwrap());
    let mismatch = smaller.cluster_hierarchical(&matrix, Linkage::Single, Cut::Count(1));
    assert_eq!(
        mismatch,
        Err(ClusteringError::MatrixMismatch {
            expected: 3,
            actual: 8
        })
    );
    assert_eq!(
        mismatch.unwrap_err().to_string(),
        "The distance matrix covers 8 members, but there are 3 to cluster."
    );
    let refusal = set
        .cluster_by_radius(radius, &atoms, by_energy)
        .unwrap_err();
    assert_eq!(refusal.to_string(), "Cannot cluster the sample set.");
    assert_eq!(
        std::error::Error::source(&refusal).unwrap().to_string(),
        "The sample set needs an energy for every state."
    );

    // A trajectory's frames, as a set without a distribution or as an
    // ensemble of equally weighted frames, cluster alike.
    let frames = trajectory();
    let set = SampleSet::from_trajectory(frames.clone()).unwrap();
    clustering_matches_ensemble(&set, &frames.into_ensemble(), false);
    assert_eq!(
        set.cluster_by_radius(radius, &atoms, by_energy),
        Err(ClusteringError::SampleSet(Box::new(
            SampleSetError::Missing(Requirement::Distribution)
        )))
    );
}

#[test]
fn clustering_part_of_a_time_ordered_set_groups_its_members_by_position() {
    let set = equilibrated();
    clustering_matches_ensemble(&set, &members_ensemble(&set), true);
    let atoms = all(set.shared_topology());
    let matrix = set
        .pairwise_fitted_rmsd(&atoms, Default::default())
        .unwrap();
    assert_eq!(matrix.member_count(), 5);
    // Entry (0, 1) is between the first two members, samples 2 and 3.
    let fit = kekule::alignment::fit(
        set.get(2).unwrap().as_model_view(),
        set.get(3).unwrap().as_model_view(),
        &atoms,
        Default::default(),
    )
    .unwrap();
    let distance = matrix.get(0, 1).unwrap().into_value();
    assert!((distance - fit.rmsd().into_value()).abs() < 1.0e-12);

    // Even and odd samples have different shapes: two medoids group the
    // members by shape, each member by its position.
    let states = set.cluster_k_medoids(&matrix, 2).unwrap();
    let shapes = [2, 3, 5, 6, 7].map(|sample| sample % 2);
    for (position, state) in states.assignments().iter().enumerate() {
        let representative = states.representatives()[*state];
        assert_eq!(shapes[position], shapes[representative]);
    }

    // The partition installs as the distribution's grouping, keeping every
    // weight, member, and the order.
    let mut regrouped = set.clone();
    regrouped
        .distribution_mut()
        .unwrap()
        .regroup(states.clone())
        .unwrap();
    let distribution = regrouped.distribution().unwrap();
    assert_eq!(distribution.partition(), &states);
    assert_eq!(
        distribution.weights(),
        set.distribution().unwrap().weights()
    );
    assert_eq!(
        distribution.members(),
        &Members::Subset(vec![2, 3, 5, 6, 7])
    );
    assert_eq!(regrouped.order(), set.order());
    assert!(regrouped.get(4).unwrap().state().is_none());

    // Non-members are not clustered, even when no fit to them succeeds.
    let mut broken = set.clone();
    broken
        .get_mut(4)
        .unwrap()
        .set_positions(Quantity::new(collinear(), NANOMETER))
        .unwrap();
    assert_eq!(
        broken.pairwise_fitted_rmsd(&atoms, Default::default()),
        Ok(matrix.clone())
    );
    // A failing member is named by its sample index.
    broken
        .get_mut(5)
        .unwrap()
        .set_positions(Quantity::new(collinear(), NANOMETER))
        .unwrap();
    assert_eq!(
        broken.pairwise_fitted_rmsd(&atoms, Default::default()),
        Err(ClusteringError::PairFit {
            first: 2,
            second: 5,
            source: Box::new(AlignmentError::DegenerateGeometry {
                geometry: AlignmentGeometry::Reference,
            }),
        })
    );
    assert_eq!(
        broken
            .pairwise_fitted_rmsd(&atoms, Default::default())
            .unwrap_err()
            .to_string(),
        "Cannot fit samples 2 and 5."
    );
    // The heaviest member, sample 5, leads; sample 3, visited next, fails
    // to fit onto it.
    assert_eq!(
        broken.cluster_by_radius(Quantity::new(0.05, NANOMETER), &atoms, Default::default()),
        Err(ClusteringError::Alignment(Box::new(AlignmentError::Item {
            index: 3,
            source: Box::new(AlignmentError::DegenerateGeometry {
                geometry: AlignmentGeometry::Reference,
            }),
        })))
    );
}

#[test]
fn sample_set_clustering_runs_on_threads_report_progress_and_cancel() {
    let set = equilibrated();
    let atoms = all(set.shared_topology());
    let threads = ScopedThreads::new(NonZeroUsize::new(3).unwrap());
    let fitted = AtomicUsize::new(0);
    let progress = |progress: ClusteringProgress| {
        fitted.store(progress.completed(), Ordering::Relaxed);
    };
    let run = ClusteringRun::default()
        .with_executor(&threads)
        .with_progress(&progress);
    let matrix = set
        .pairwise_fitted_rmsd_with(&atoms, Default::default(), run)
        .unwrap();
    assert_eq!(
        matrix,
        set.pairwise_fitted_rmsd(&atoms, Default::default())
            .unwrap()
    );
    assert_eq!(fitted.load(Ordering::Relaxed), 10);
    let radius = Quantity::new(0.08, NANOMETER);
    assert_eq!(
        set.cluster_by_radius_with(radius, &atoms, Default::default(), run)
            .unwrap(),
        set.cluster_by_radius(radius, &atoms, Default::default())
            .unwrap()
    );
    assert_eq!(fitted.load(Ordering::Relaxed), 5);

    let cancellation = Cancellation::new(|| true);
    let cancelled = ClusteringRun::default().with_cancellation(&cancellation);
    assert!(matches!(
        set.pairwise_fitted_rmsd_with(&atoms, Default::default(), cancelled),
        Err(ClusteringError::Cancelled(_))
    ));
    assert!(matches!(
        set.cluster_by_radius_with(radius, &atoms, Default::default(), cancelled),
        Err(ClusteringError::Cancelled(_))
    ));
}

#[test]
fn a_set_without_a_distribution_installs_its_clusters_as_one() {
    let mut set = SampleSet::from_trajectory(trajectory()).unwrap();
    let atoms = all(set.shared_topology());
    let matrix = set
        .pairwise_fitted_rmsd(&atoms, Default::default())
        .unwrap();
    assert_eq!(matrix.member_count(), 8);
    let states = set
        .cluster_hierarchical(&matrix, Linkage::Average, Cut::Count(2))
        .unwrap();
    let distribution = Distribution::new(Members::All, states.clone()).unwrap();
    set.set_distribution(Some(distribution)).unwrap();
    assert_eq!(set.distribution().unwrap().partition(), &states);
    assert_eq!(set.states().len(), 2);
    // Frames of one shape share a state.
    for sample in set.iter() {
        let representative = sample.state().unwrap().representative();
        assert_eq!(sample.index() % 2, representative % 2);
    }
}

#[test]
fn boltzmann_populations_of_a_distribution_equal_an_ensemble_s() {
    for temperature in [1.0, 298.15, 5000.0] {
        let temperature = Quantity::new(temperature, KELVIN);
        let mut old = ensemble(true);
        let mut set = SampleSet::from_ensemble(old.clone());
        old.assign_boltzmann_populations(temperature).unwrap();
        set.distribution_mut()
            .unwrap()
            .assign_boltzmann_populations(temperature)
            .unwrap();
        let expected = SampleSet::from_ensemble(old);
        let (new, old) = (
            set.distribution().unwrap(),
            expected.distribution().unwrap(),
        );
        assert_eq!(bits(new.weights().unwrap()), bits(old.weights().unwrap()));
        assert_eq!(new, old);
        assert_eq!(new.temperature(), Some(temperature));
    }

    // Weights near the floating-point limits take the within-state path.
    let mut old = Ensemble::new(topology());
    for (index, (weight, energy)) in [(1.0e300, 0.0), (1.0e-300, 1.0), (2.0e-300, 1.0)]
        .into_iter()
        .enumerate()
    {
        let sample = EnsembleSample::new(positions(points(index)), weight).unwrap();
        if index < 2 {
            let state = EnsembleState::new()
                .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
                .unwrap();
            old.push_state(sample, state).unwrap();
        } else {
            old.push_sample(sample, 1).unwrap();
        }
    }
    let mut set = SampleSet::from_ensemble(old.clone());
    let temperature = Quantity::new(300.0, KELVIN);
    old.assign_boltzmann_populations(temperature).unwrap();
    set.distribution_mut()
        .unwrap()
        .assign_boltzmann_populations(temperature)
        .unwrap();
    assert_eq!(
        set.distribution(),
        SampleSet::from_ensemble(old).distribution()
    );

    // A distribution over part of a set without weights gains them.
    let mut part = equilibrated();
    part.distribution_mut().unwrap().set_weights(None).unwrap();
    let mut old = members_ensemble(&part);
    old.assign_boltzmann_populations(temperature).unwrap();
    part.distribution_mut()
        .unwrap()
        .assign_boltzmann_populations(temperature)
        .unwrap();
    assert_eq!(
        part.distribution().unwrap().weights(),
        Some(
            old.iter()
                .map(|sample| sample.weight())
                .collect::<Vec<_>>()
                .as_slice()
        )
    );
    assert_eq!(
        part.distribution().unwrap().temperature(),
        Some(temperature)
    );
}

#[test]
fn boltzmann_populations_refuse_without_energies_and_change_nothing() {
    let mut set = SampleSet::from_ensemble(ensemble(false));
    let before = set.distribution().cloned();
    let mut distribution = set.distribution_mut().unwrap();
    assert_eq!(
        distribution.assign_boltzmann_populations(Quantity::new(300.0, KELVIN)),
        Err(SampleSetError::Missing(Requirement::StateEnergies))
    );
    assert_eq!(
        distribution.assign_boltzmann_populations(Quantity::new(-1.0, KELVIN)),
        Err(SampleSetError::InvalidTemperature)
    );
    assert!(matches!(
        distribution.assign_boltzmann_populations(Quantity::new(300.0, NANOMETER)),
        Err(SampleSetError::Unit(_))
    ));
    assert_eq!(set.distribution().cloned(), before);
    assert_eq!(
        SampleSetError::Missing(Requirement::StateEnergies).to_string(),
        "The sample set needs an energy for every state."
    );

    // A state far above the others underflows like an ensemble's.
    let mut old = ensemble(false);
    for (state, energy) in [0.0, 1.0e6, 2.0, 3.0].into_iter().enumerate() {
        *old.state_mut(state).unwrap() = EnsembleState::new()
            .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
            .unwrap();
    }
    let mut set = SampleSet::from_ensemble(old.clone());
    let cold = Quantity::new(1.0, KELVIN);
    assert!(old.assign_boltzmann_populations(cold).is_err());
    let before = set.distribution().cloned();
    assert_eq!(
        set.distribution_mut()
            .unwrap()
            .assign_boltzmann_populations(cold),
        Err(SampleSetError::InvalidWeight { sample: 2 })
    );
    assert_eq!(set.distribution().cloned(), before);
}
