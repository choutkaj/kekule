//! Superposition, RMSD, weighted statistics, clustering, and Boltzmann
//! populations of structure sets: the values the retired ensemble, trajectory,
//! and series paths gave, and restricted to the members of a distribution
//! over part of a set.

use std::fmt::Debug;
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
use kekule::math::{Matrix3, Point3, RigidTransform, Vector3};
use kekule::parallel::ScopedThreads;
use kekule::statistics::{
    ContactOccupancyResult, GroupContactOccupancyResult, RmsfResult, StatisticsError,
};
use kekule::structure::{
    Distribution, Geometry, Members, Order, Partition, Positions, Requirement, State, Structure,
    StructureSet, StructureSetError, TimeAxis, Velocities,
};
use kekule::topology::{AtomSelection, InstanceAtomId, Topology};
use kekule::units::{
    Quantity, CANONICAL_VELOCITY_UNIT, DEGREE, KELVIN, KILOJOULE_PER_MOLE, NANOMETER, PICOSECOND,
};

use crate::fingerprint::Recorded;

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
/// The state of each geometry of [`ensemble`], its first geometry representing it.
const STATES: [usize; 8] = [0, 0, 1, 1, 2, 2, 2, 3];
const ENERGIES: [f64; 4] = [1.5, 0.0, 4.0, 0.7];

/// Unordered geometries of [`points`] `0..count`.
fn unordered(count: usize) -> StructureSet {
    let geometries = (0..count).map(|index| Geometry::new(positions(points(index))));
    StructureSet::from_geometries(topology(), geometries).unwrap()
}

/// Eight weighted geometries in four states, with state energies or without.
fn ensemble(energies: bool) -> StructureSet {
    let mut set = unordered(8);
    let states = ENERGIES.map(|energy| {
        if energies {
            State::new()
                .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
                .unwrap()
        } else {
            State::new()
        }
    });
    let partition = Partition::new(STATES.to_vec(), vec![0, 2, 4, 7]).unwrap();
    let distribution = Distribution::new(Members::All, partition)
        .unwrap()
        .with_weights(WEIGHTS.to_vec())
        .unwrap()
        .with_states(states.to_vec())
        .unwrap();
    set.set_distribution(Some(distribution)).unwrap();
    set
}

/// Eight frames two picoseconds apart, with velocities.
fn trajectory() -> StructureSet {
    let geometries = (0..8).map(|index| {
        let mut geometry = Geometry::new(positions(points(index)));
        let velocities = (0..ATOMS)
            .map(|atom| Vector3::new(atom as f64, 1.0, -(index as f64)))
            .collect::<Vec<_>>();
        geometry
            .set_velocities(Some(
                Velocities::new(Quantity::new(velocities, CANONICAL_VELOCITY_UNIT)).unwrap(),
            ))
            .unwrap();
        geometry
    });
    let mut set = StructureSet::from_geometries(topology(), geometries).unwrap();
    let times = (0..8).map(|index| 2.0 * index as f64).collect::<Vec<_>>();
    let axis = TimeAxis::new(Some(Quantity::new(times, PICOSECOND)), None).unwrap();
    set.set_order(Some(Order::Time(axis))).unwrap();
    set
}

/// Six geometries of a scan of one angle.
fn series() -> StructureSet {
    let mut set = unordered(6);
    let values = (0..6).map(|index| 30.0 * index as f64).collect();
    set.set_order(Some(Order::Coordinate {
        unit: DEGREE,
        values,
    }))
    .unwrap();
    set
}

/// What a set holds apart from its topology, whose layout identity differs
/// between runs.
fn contents(set: &StructureSet) -> impl Debug {
    (
        geometries(set),
        set.order().cloned(),
        set.distribution().cloned(),
        set.properties().clone(),
    )
}

#[test]
fn fixtures_are_the_sets_the_retired_collections_converted_to() {
    // The sets the retired ensemble, trajectory, and series fixtures
    // converted to, and the extreme ensemble of the Boltzmann checks.
    let mut recorded = Recorded::new(&[
        0x5fc83223b464ac24,
        0x774245ee1d22d7ad,
        0xd73ec02d570e765d,
        0x686d98249eb9038b,
        0xf74af6e81faa6f7a,
    ]);
    recorded.check(&contents(&ensemble(true)));
    recorded.check(&contents(&ensemble(false)));
    recorded.check(&contents(&trajectory()));
    recorded.check(&contents(&series()));
    recorded.check(&contents(&extreme()));
}

fn all(set_topology: Arc<Topology>) -> AtomSelection {
    AtomSelection::all(set_topology)
}

fn selection(topology: &Arc<Topology>, atoms: &[usize]) -> AtomSelection {
    let ids = topology.atom_ids();
    AtomSelection::from_atoms(Arc::clone(topology), atoms.iter().map(|atom| ids[*atom])).unwrap()
}

fn geometries(set: &StructureSet) -> Vec<Geometry> {
    set.iter()
        .map(|geometry| geometry.payload().clone())
        .collect()
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

/// The values of an RMSF result, apart from the selection it was asked for.
fn rmsf_values(result: &RmsfResult) -> impl Debug + PartialEq {
    (
        bits(result.values().value()),
        result.values().unit(),
        result.observation_count(),
        result.total_weight().to_bits(),
    )
}

fn contact_values(result: &ContactOccupancyResult) -> impl Debug + PartialEq {
    (
        bits(result.contact_weights()),
        result.pairs().to_vec(),
        result.cutoff(),
        result.observation_count(),
        result.total_weight().to_bits(),
    )
}

/// The values of a group contact result, apart from the groups it was
/// asked for.
fn group_contact_values(result: &GroupContactOccupancyResult) -> impl Debug + PartialEq {
    (
        bits(result.contact_weights()),
        result.pairs().to_vec(),
        result.cutoff(),
        result.observation_count(),
        result.total_weight().to_bits(),
    )
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

/// Every statistic of `set`: the RMSF over all atoms and over three, and
/// the contact and group contact occupancies.
fn statistics(set: &StructureSet) -> impl Debug + PartialEq {
    let topology = set.shared_topology();
    let cutoff = Quantity::new(CUTOFF, NANOMETER);
    let rmsf = [all(Arc::clone(&topology)), selection(&topology, &[1, 4, 5])]
        .map(|atoms| rmsf_values(&set.rmsf(&atoms).unwrap()));
    let contacts = contact_values(
        &set.contact_occupancy(contact_pairs(&topology), cutoff)
            .unwrap(),
    );
    let (groups, pairs) = contact_groups(&topology);
    let group_contacts =
        group_contact_values(&set.group_contact_occupancy(groups, pairs, cutoff).unwrap());
    (rmsf, contacts, group_contacts)
}

/// Every alignment of a set, as `(rmsd, aligned rmsd, superposition
/// report)`.
type AlignmentResults = (
    Result<Vec<u64>, AlignmentError>,
    Result<Vec<u64>, AlignmentError>,
    Result<kekule::alignment::SuperpositionReport, AlignmentError>,
);

/// Measures and then superposes a set with one set of atoms, weights, and
/// options, against `reference` or else its geometry `index`.
fn alignments(
    set: &mut StructureSet,
    reference: Option<&Structure>,
    index: usize,
) -> AlignmentResults {
    let topology = set.shared_topology();
    let fit = selection(&topology, &[0, 1, 2, 3]);
    let measured = selection(&topology, &[2, 3, 4, 5]);
    let pair_weights = [1.0, 2.0, 0.5, 3.0];
    let options = AlignmentOptions::default().with_weighting(Weighting::Explicit(&pair_weights));
    let aligned = AlignedRmsdOptions::default()
        .with_fit(options)
        .with_measurement_weighting(Weighting::Explicit(&pair_weights));
    let values = |result: Result<Quantity<Vec<f64>>, AlignmentError>| {
        result.map(|values| bits(values.value()))
    };
    match reference {
        Some(structure) => (
            values(set.rmsd(structure, &fit, options)),
            values(set.aligned_rmsd(structure, &fit, &measured, aligned)),
            set.superpose(structure, &fit, options),
        ),
        None => (
            values(set.rmsd(index, &fit, options)),
            values(set.aligned_rmsd(index, &fit, &measured, aligned)),
            set.superpose(index, &fit, options),
        ),
    }
}

#[test]
fn structure_set_alignments_keep_the_retired_values() {
    // Fingerprints of the values the retired paths gave.
    let mut recorded = Recorded::new(&[
        0xa708ddfed8a897b0,
        0x09d4654d2788ed2f,
        0xa708ddfed8a897b0,
        0xde6b72c6afe8f797,
        0x7d356d80e26864d5,
        0x1778f93bf2568fda,
        0xa46e16252c8a8818,
        0x8b60a38a15003ee9,
        0xa46e16252c8a8818,
        0x92e5d1e483871d4b,
        0xef9aa5844b998ecc,
        0x8ad435bcdd7707a3,
        0xd0250a17c77174b1,
        0x1a2abdde993e8a8a,
        0xd0250a17c77174b1,
        0x434c38a5bc6bdb2e,
        0x4f8e715505afaaa7,
        0x8d130dc2f8772296,
    ]);
    let structure = Structure::new(topology(), positions(points(11))).unwrap();
    // A geometry of the set, an independent view, and a reference out of range.
    for (reference, index) in [(None, 3), (Some(&structure), 0), (None, 9)] {
        // Velocities rotate with the frames, and the order is kept.
        for fixture in [ensemble(true), trajectory(), series()] {
            let mut set = fixture.clone();
            let results = alignments(&mut set, reference, index);
            assert_eq!(results.2.is_ok(), index < fixture.len());
            recorded.check(&results);
            recorded.check(&geometries(&set));
            assert_eq!(set.order(), fixture.order());
            assert_eq!(set.distribution(), fixture.distribution());
        }
    }
}

#[test]
fn structure_set_alignment_covers_non_members_and_keeps_the_distribution() {
    // A failing geometry outside the distribution still fails the fit, and a
    // fit keeps the order and the distribution.
    let mut set = trajectory();
    let distribution =
        Distribution::new(Members::Subset(vec![1, 4]), Partition::state_per_sample(2))
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
    let before = geometries(&broken);
    assert_eq!(
        broken.superpose(6, &atoms, Default::default()),
        Err(AlignmentError::Item {
            index: 0,
            source: Box::new(AlignmentError::DegenerateGeometry {
                geometry: AlignmentGeometry::Reference,
            }),
        })
    );
    assert_eq!(geometries(&broken), before);
}

#[test]
fn structure_set_statistics_keep_the_retired_values() {
    // Fingerprints of the values the retired paths gave.
    let mut recorded = Recorded::new(&[0x9c767ac7a8c90863, 0x9c767ac7a8c90863, 0x3a279110203f7b17]);
    for energies in [true, false] {
        recorded.check(&statistics(&ensemble(energies)));
    }

    // Without a distribution, every frame counts with weight one.
    let set = trajectory();
    recorded.check(&statistics(&set));
    let rmsf = set.rmsf(&all(set.shared_topology())).unwrap();
    assert_eq!(rmsf.total_weight(), 8.0);

    let empty = StructureSet::new(self::topology());
    assert_eq!(
        empty.rmsf(&all(empty.shared_topology())).unwrap_err(),
        StatisticsError::NoObservations
    );
    assert!(matches!(
        set.rmsf(&all(another_topology())),
        Err(StatisticsError::Selection(_))
    ));
}

/// The members of `set`'s distribution as a set of their own, unordered,
/// with their weights and states, and without velocities.
fn members(set: &StructureSet) -> StructureSet {
    let members = set
        .weighted_samples()
        .map(|(geometry, _)| geometry)
        .collect::<Vec<_>>();
    let mut members = set.select(members).unwrap();
    members.set_order(None).unwrap();
    for index in 0..members.len() {
        members.get_mut(index).unwrap().take_velocities();
    }
    members
}

/// A time-ordered set of eight frames whose distribution covers frames 2,
/// 3, 5, 6, and 7, weighted, in two states with energies.
fn equilibrated() -> StructureSet {
    let mut set = trajectory();
    let states = [0.0, 2.5].map(|energy| {
        State::new()
            .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
            .unwrap()
    });
    let distribution = Distribution::new(
        Members::Subset(vec![2, 3, 5, 6, 7]),
        Partition::new(vec![0, 1, 0, 1, 1], vec![0, 1]).unwrap(),
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
    // Fingerprints of the values the retired paths gave.
    let mut recorded = Recorded::new(&[0x5e9439675b22d875]);
    let set = equilibrated();
    let observed = statistics(&set);
    recorded.check(&observed);
    assert!(observed == statistics(&members(&set)));
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
    assert!(statistics(&moved) == observed);

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
    let mut frames = equal.select([2, 3, 5, 6, 7]).unwrap();
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

/// Every clustering of `set`: its fitted distance matrix, clusters of five
/// radii by weight and, with `energies`, by energy, hierarchical clusters of
/// every count and one height under three linkages, and k-medoids of every
/// count.
fn clusterings(set: &StructureSet, energies: bool) -> impl Debug + PartialEq {
    let atoms = all(set.shared_topology());
    // The two shapes are two clusters of one radius.
    let shapes = set
        .cluster_by_radius(Quantity::new(0.08, NANOMETER), &atoms, Default::default())
        .unwrap();
    assert_eq!(shapes.state_count(), 2);
    let matrix = set
        .pairwise_fitted_rmsd(&atoms, Default::default())
        .unwrap();
    let mut orders = vec![ClusterOrder::Weight];
    if energies {
        orders.push(ClusterOrder::Energy);
    }
    let mut radii = Vec::new();
    for order in orders {
        for radius in [0.0, 0.05, 0.08, 0.2, 10.0] {
            let options = ClusteringOptions::default().with_order(order);
            let radius = Quantity::new(radius, NANOMETER);
            radii.push(set.cluster_by_radius(radius, &atoms, options).unwrap());
        }
    }
    let mut hierarchical = Vec::new();
    for linkage in [Linkage::Single, Linkage::Complete, Linkage::Average] {
        for count in 1..=matrix.member_count() {
            let cut = Cut::Count(count);
            hierarchical.push(set.cluster_hierarchical(&matrix, linkage, cut).unwrap());
        }
        let cut = Cut::Height(Quantity::new(0.06, NANOMETER));
        hierarchical.push(set.cluster_hierarchical(&matrix, linkage, cut).unwrap());
    }
    let medoids = (1..=matrix.member_count())
        .map(|clusters| set.cluster_k_medoids(&matrix, clusters).unwrap())
        .collect::<Vec<_>>();
    (matrix, radii, hierarchical, medoids)
}

#[test]
fn structure_set_clusterings_keep_the_retired_values() {
    // Fingerprints of the values the retired paths gave.
    let mut recorded = Recorded::new(&[
        0x38518c62c0a59494,
        0xf515e1b2b996a0fe,
        0x7b21e643d2f3651e,
        0x2c116153dbaf17bd,
    ]);
    for energies in [true, false] {
        recorded.check(&clusterings(&ensemble(energies), energies));
    }
    // Without energies, ordering by energy is refused.
    let set = ensemble(false);
    let atoms = all(set.shared_topology());
    let by_energy = ClusteringOptions::default().with_order(ClusterOrder::Energy);
    let radius = Quantity::new(0.1, NANOMETER);
    assert_eq!(
        set.cluster_by_radius(radius, &atoms, by_energy),
        Err(ClusteringError::StructureSet(Box::new(
            StructureSetError::Missing(Requirement::StateEnergies)
        )))
    );
    let matrix = set
        .pairwise_fitted_rmsd(&atoms, Default::default())
        .unwrap();
    recorded.check(&set.cluster_k_medoids(&matrix, 0));
    let smaller = set.select(0..3).unwrap();
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
    assert_eq!(refusal.to_string(), "Cannot cluster the structure set.");
    assert_eq!(
        std::error::Error::source(&refusal).unwrap().to_string(),
        "The structure set needs an energy for every state."
    );

    // A trajectory's frames, as a set without a distribution, cluster as
    // equally weighted geometries did.
    let set = trajectory();
    recorded.check(&clusterings(&set, false));
    assert_eq!(
        set.cluster_by_radius(radius, &atoms, by_energy),
        Err(ClusteringError::StructureSet(Box::new(
            StructureSetError::Missing(Requirement::Distribution)
        )))
    );
}

#[test]
fn clustering_part_of_a_time_ordered_set_groups_its_members_by_position() {
    // Fingerprints of the values the retired paths gave.
    let mut recorded = Recorded::new(&[0x25b0a03d01643e37]);
    let set = equilibrated();
    let clustered = clusterings(&set, true);
    recorded.check(&clustered);
    assert!(clustered == clusterings(&members(&set), true));
    let atoms = all(set.shared_topology());
    let matrix = set
        .pairwise_fitted_rmsd(&atoms, Default::default())
        .unwrap();
    assert_eq!(matrix.member_count(), 5);
    // Entry (0, 1) is between the first two members, geometries 2 and 3.
    let fit = kekule::alignment::fit(
        set.get(2).unwrap().as_structure_view(),
        set.get(3).unwrap().as_structure_view(),
        &atoms,
        Default::default(),
    )
    .unwrap();
    let distance = matrix.get(0, 1).unwrap().into_value();
    assert!((distance - fit.rmsd().into_value()).abs() < 1.0e-12);

    // Even and odd geometries have different shapes: two medoids group the
    // members by shape, each member by its position.
    let states = set.cluster_k_medoids(&matrix, 2).unwrap();
    let shapes = [2, 3, 5, 6, 7].map(|geometry| geometry % 2);
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
    // A failing member is named by its geometry index.
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
    // The heaviest member, geometry 5, leads; geometry 3, visited next, fails
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
fn structure_set_clustering_runs_on_threads_report_progress_and_cancel() {
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
    let mut set = trajectory();
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
    for geometry in set.iter() {
        let representative = geometry.state().unwrap().representative();
        assert_eq!(geometry.index() % 2, representative % 2);
    }
}

/// Three geometries whose weights are near the floating-point limits, in two
/// states with energies.
fn extreme() -> StructureSet {
    let mut set = unordered(3);
    let states = [0.0, 1.0].map(|energy| {
        State::new()
            .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
            .unwrap()
    });
    let distribution = Distribution::new(
        Members::All,
        Partition::new(vec![0, 1, 1], vec![0, 1]).unwrap(),
    )
    .unwrap()
    .with_weights(vec![1.0e300, 1.0e-300, 2.0e-300])
    .unwrap()
    .with_states(states.to_vec())
    .unwrap();
    set.set_distribution(Some(distribution)).unwrap();
    set
}

#[test]
fn boltzmann_populations_keep_the_retired_values() {
    // Fingerprints of the values the retired paths gave.
    let mut recorded = Recorded::new(&[
        0x068a61a3e61a908f,
        0x85c659922dd3ad74,
        0xb077e016f001eac2,
        0x7f9b35233a28d61f,
        0xa88494baf97eb0ab,
    ]);
    for temperature in [1.0, 298.15, 5000.0] {
        let temperature = Quantity::new(temperature, KELVIN);
        let mut set = ensemble(true);
        set.distribution_mut()
            .unwrap()
            .assign_boltzmann_populations(temperature)
            .unwrap();
        let distribution = set.distribution().unwrap();
        recorded.check(&(bits(distribution.weights().unwrap()), distribution));
        assert_eq!(distribution.temperature(), Some(temperature));
    }

    // Weights near the floating-point limits take the within-state path.
    let mut set = extreme();
    let temperature = Quantity::new(300.0, KELVIN);
    set.distribution_mut()
        .unwrap()
        .assign_boltzmann_populations(temperature)
        .unwrap();
    let distribution = set.distribution().unwrap();
    recorded.check(&(bits(distribution.weights().unwrap()), distribution));

    // A distribution over part of a set without weights gains them.
    let mut part = equilibrated();
    part.distribution_mut().unwrap().set_weights(None).unwrap();
    part.distribution_mut()
        .unwrap()
        .assign_boltzmann_populations(temperature)
        .unwrap();
    recorded.check(&bits(part.distribution().unwrap().weights().unwrap()));
    assert_eq!(
        part.distribution().unwrap().temperature(),
        Some(temperature)
    );
}

#[test]
fn boltzmann_populations_refuse_without_energies_and_change_nothing() {
    let mut set = ensemble(false);
    let before = set.distribution().cloned();
    let mut distribution = set.distribution_mut().unwrap();
    assert_eq!(
        distribution.assign_boltzmann_populations(Quantity::new(300.0, KELVIN)),
        Err(StructureSetError::Missing(Requirement::StateEnergies))
    );
    assert_eq!(
        distribution.assign_boltzmann_populations(Quantity::new(-1.0, KELVIN)),
        Err(StructureSetError::InvalidTemperature)
    );
    assert!(matches!(
        distribution.assign_boltzmann_populations(Quantity::new(300.0, NANOMETER)),
        Err(StructureSetError::Unit(_))
    ));
    assert_eq!(set.distribution().cloned(), before);
    assert_eq!(
        StructureSetError::Missing(Requirement::StateEnergies).to_string(),
        "The structure set needs an energy for every state."
    );

    // A state far above the others underflows.
    let mut set = ensemble(false);
    let mut distribution = set.distribution_mut().unwrap();
    for (state, energy) in [0.0, 1.0e6, 2.0, 3.0].into_iter().enumerate() {
        *distribution.state_mut(state).unwrap() = State::new()
            .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
            .unwrap();
    }
    let cold = Quantity::new(1.0, KELVIN);
    let before = set.distribution().cloned();
    assert_eq!(
        set.distribution_mut()
            .unwrap()
            .assign_boltzmann_populations(cold),
        Err(StructureSetError::InvalidWeight { geometry: 2 })
    );
    assert_eq!(set.distribution().cloned(), before);
}
