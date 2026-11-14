//! Superposition, RMSD, statistics, clustering, and Boltzmann populations of
//! structure sets, unweighted over every geometry, and of ensembles, weighted
//! over the samples of their distributions, each checked against its
//! definition evaluated directly from the coordinates and weights.

use std::fmt::Debug;
use std::num::NonZeroUsize;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, OnceLock};

use kekule::alignment::{
    self, AlignedRmsdOptions, AlignmentError, AlignmentGeometry, AlignmentOptions, Weighting,
};
use kekule::clustering::{
    ClusterOrder, ClusteringError, ClusteringOptions, ClusteringProgress, ClusteringRun, Cut,
    DistanceMatrix, Linkage,
};
use kekule::diagnostics::Cancellation;
use kekule::math::{Matrix3, Point3, RigidTransform, Vector3};
use kekule::parallel::ScopedThreads;
use kekule::statistics::{
    ContactOccupancyResult, GroupContactOccupancyResult, RmsfResult, StatisticsError,
};
use kekule::structure::{
    Distribution, Ensemble, Geometry, Order, Partition, Positions, Requirement, Samples, StateData,
    Structure, StructureSet, StructureSetError, StructureView, TimeAxis, Velocities,
};
use kekule::topology::{AtomSelection, InstanceAtomId, Topology};
use kekule::units::{
    Quantity, CANONICAL_VELOCITY_UNIT, DEGREE, KELVIN, KILOJOULE_PER_MOLE, NANOMETER, PICOSECOND,
};

const ATOMS: usize = 6;

/// The molar gas constant in kJ/(mol K): the Avogadro constant times the
/// Boltzmann constant, both exact since the 2019 SI.
const GAS_CONSTANT: f64 = 6.022_140_76e23 * 1.380_649e-23 * 1.0e-3;

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

/// The velocity of `atom` in frame `index` of [`trajectory`].
fn velocity(index: usize, atom: usize) -> Vector3 {
    Vector3::new(atom as f64, 1.0, -(index as f64))
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

/// The natural logarithms of linear weights.
fn ln(weights: &[f64]) -> Vec<f64> {
    weights.iter().map(|weight| weight.ln()).collect()
}

/// Unordered geometries of [`points`] `0..count`.
fn unordered(count: usize) -> StructureSet {
    let geometries = (0..count).map(|index| Geometry::new(positions(points(index))));
    StructureSet::from_geometries(topology(), geometries).unwrap()
}

/// Eight weighted geometries in four states, with state energies or without.
fn ensemble(energies: bool) -> Ensemble {
    let states = ENERGIES.map(|energy| {
        if energies {
            StateData::new()
                .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
                .unwrap()
        } else {
            StateData::new()
        }
    });
    let partition = Partition::new(STATES.to_vec(), vec![0, 2, 4, 7]).unwrap();
    let distribution = Distribution::new(Samples::All, partition)
        .unwrap()
        .with_log_weights(ln(&WEIGHTS))
        .unwrap()
        .with_states(states.to_vec())
        .unwrap();
    Ensemble::new(unordered(8), distribution).unwrap()
}

/// Eight frames two picoseconds apart, with velocities.
fn trajectory() -> StructureSet {
    let geometries = (0..8).map(|index| {
        let mut geometry = Geometry::new(positions(points(index)));
        let velocities = (0..ATOMS)
            .map(|atom| velocity(index, atom))
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

/// The positions of geometry `index` of `set`, in nanometres.
fn nanometers(set: &StructureSet, index: usize) -> Vec<Point3> {
    let geometry = set.get(index).unwrap();
    let positions = geometry.positions().values();
    assert_eq!(positions.unit(), NANOMETER);
    positions.value().to_vec()
}

/// The velocities of geometry `index` of `set`, in the canonical unit.
fn velocities(set: &StructureSet, index: usize) -> Vec<Vector3> {
    let geometry = set.get(index).unwrap();
    let velocities = geometry.velocities().unwrap().values();
    assert_eq!(velocities.unit(), CANONICAL_VELOCITY_UNIT);
    velocities.value().to_vec()
}

fn assert_close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected} within {tolerance}"
    );
}

fn assert_points_close(actual: &[Point3], expected: &[Point3], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        assert!(
            (*actual - *expected).norm() <= tolerance,
            "{actual:?} != {expected:?}"
        );
    }
}

/// Asserts that `set` holds the geometries of [`points`] `0..count`.
fn assert_holds_points(set: &StructureSet, count: usize) {
    assert_eq!(set.len(), count);
    assert!(set.topology().shares_layout(&topology()));
    for index in 0..count {
        assert_eq!(nanometers(set, index), points(index));
    }
}

#[test]
fn fixtures_hold_their_geometries_order_and_distribution() {
    for energies in [true, false] {
        let ensemble = ensemble(energies);
        assert_holds_points(ensemble.set(), 8);
        assert_eq!(ensemble.set().order(), None);
        assert!(ensemble.set().iter().all(|g| g.velocities().is_none()));
        let distribution = ensemble.distribution();
        assert_eq!(distribution.samples(), &Samples::All);
        assert_eq!(distribution.partition().assignments(), STATES);
        assert_eq!(distribution.partition().representatives(), [0, 2, 4, 7]);
        assert_eq!(distribution.log_weights(), Some(&ln(&WEIGHTS)[..]));
        assert_eq!(distribution.temperature(), None);
        let expected =
            ENERGIES.map(|energy| energies.then(|| Quantity::new(energy, KILOJOULE_PER_MOLE)));
        for (state, energy) in ensemble.states().zip(expected) {
            assert_eq!(state.energy(), energy);
        }
    }

    let trajectory = trajectory();
    assert_holds_points(&trajectory, 8);
    let Some(Order::Time(axis)) = trajectory.order() else {
        panic!("{:?}", trajectory.order());
    };
    let times = axis.times().unwrap();
    assert_eq!(times.unit(), PICOSECOND);
    assert_eq!(*times.value(), [0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0]);
    for index in 0..8 {
        let expected = (0..ATOMS)
            .map(|atom| velocity(index, atom))
            .collect::<Vec<_>>();
        assert_eq!(velocities(&trajectory, index), expected);
    }

    let series = series();
    assert_holds_points(&series, 6);
    assert_eq!(
        series.order(),
        Some(&Order::Coordinate {
            unit: DEGREE,
            values: vec![0.0, 30.0, 60.0, 90.0, 120.0, 150.0],
        })
    );

    let extreme = extreme();
    assert_holds_points(extreme.set(), 3);
    let distribution = extreme.distribution();
    assert_eq!(distribution.partition().assignments(), [0, 1, 1]);
    assert_eq!(
        distribution.log_weights(),
        Some(&ln(&[1.0e300, 1.0e-300, 2.0e-300])[..])
    );
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

const CONTACT_PAIRS: [(usize, usize); 3] = [(0, 5), (1, 4), (0, 1)];

fn contact_pairs(topology: &Topology) -> Vec<(InstanceAtomId, InstanceAtomId)> {
    let ids = topology.atom_ids();
    CONTACT_PAIRS
        .iter()
        .map(|(first, second)| (ids[*first], ids[*second]))
        .collect()
}

const GROUPS: [[usize; 2]; 3] = [[0, 1], [4, 5], [2, 3]];
const GROUP_PAIRS: [(usize, usize); 3] = [(0, 1), (2, 0), (1, 2)];

/// The groups and group pairs of the group contact checks.
fn contact_groups(topology: &Arc<Topology>) -> (Vec<AtomSelection>, Vec<(usize, usize)>) {
    let groups = GROUPS
        .iter()
        .map(|atoms| selection(topology, atoms))
        .collect();
    (groups, GROUP_PAIRS.to_vec())
}

const CUTOFF: f64 = 0.62;
/// The atoms of the second RMSF selection.
const RMSF_ATOMS: [usize; 3] = [1, 4, 5];

/// Every statistic of a set or ensemble: the RMSF over all atoms and over
/// [`RMSF_ATOMS`], and the contact and group contact occupancies.
struct Statistics {
    rmsf: [RmsfResult; 2],
    contacts: ContactOccupancyResult,
    groups: GroupContactOccupancyResult,
}

impl Statistics {
    /// The statistics of `set`, every geometry with weight one.
    fn of_set(set: &StructureSet) -> Self {
        let topology = set.shared_topology();
        let cutoff = Quantity::new(CUTOFF, NANOMETER);
        let (groups, pairs) = contact_groups(&topology);
        Self {
            rmsf: [
                all(Arc::clone(&topology)),
                selection(&topology, &RMSF_ATOMS),
            ]
            .map(|atoms| set.rmsf(&atoms).unwrap()),
            contacts: set
                .contact_occupancy(contact_pairs(&topology), cutoff)
                .unwrap(),
            groups: set.group_contact_occupancy(groups, pairs, cutoff).unwrap(),
        }
    }

    /// The statistics of `ensemble`, its samples with their probabilities.
    fn of_ensemble(ensemble: &Ensemble) -> Self {
        let topology = ensemble.set().shared_topology();
        let cutoff = Quantity::new(CUTOFF, NANOMETER);
        let (groups, pairs) = contact_groups(&topology);
        Self {
            rmsf: [
                all(Arc::clone(&topology)),
                selection(&topology, &RMSF_ATOMS),
            ]
            .map(|atoms| ensemble.rmsf(&atoms).unwrap()),
            contacts: ensemble
                .contact_occupancy(contact_pairs(&topology), cutoff)
                .unwrap(),
            groups: ensemble
                .group_contact_occupancy(groups, pairs, cutoff)
                .unwrap(),
        }
    }

    /// Every value, bit for bit, for comparing two runs.
    fn values(&self) -> impl Debug + PartialEq {
        (
            self.rmsf.each_ref().map(rmsf_values),
            contact_values(&self.contacts),
            group_contact_values(&self.groups),
        )
    }

    /// Checks every statistic against its definition, evaluated over the
    /// geometries `observed` of `set` with their weights.
    fn check(&self, set: &StructureSet, observed: &[(usize, f64)]) {
        let total = observed.iter().map(|(_, weight)| weight).sum::<f64>();
        let atom = |geometry: usize, atom: usize| nanometers(set, geometry)[atom];
        let distance = |geometry: usize, first: usize, second: usize| {
            (atom(geometry, first) - atom(geometry, second)).norm()
        };
        for (result, atoms) in self
            .rmsf
            .iter()
            .zip([(0..ATOMS).collect::<Vec<_>>(), RMSF_ATOMS.to_vec()])
        {
            assert_eq!(result.observation_count(), observed.len());
            assert_close(result.total_weight(), total, 1.0e-12 * total);
            assert_eq!(result.values().unit(), NANOMETER);
            assert_eq!(result.values().value().len(), atoms.len());
            for (value, atom_index) in result.values().value().iter().zip(atoms) {
                // The weighted mean position, then the weighted mean squared
                // deviation from it.
                let mean =
                    observed
                        .iter()
                        .fold(Vector3::new(0.0, 0.0, 0.0), |sum, (geometry, weight)| {
                            let point = atom(*geometry, atom_index);
                            sum + Vector3::new(point.x, point.y, point.z) * (*weight / total)
                        });
                let variance = observed
                    .iter()
                    .map(|(geometry, weight)| {
                        let point = atom(*geometry, atom_index);
                        let deviation = Vector3::new(point.x, point.y, point.z) - mean;
                        weight * deviation.dot(deviation)
                    })
                    .sum::<f64>()
                    / total;
                assert_close(*value, variance.sqrt(), 1.0e-12);
            }
        }

        assert_eq!(self.contacts.observation_count(), observed.len());
        assert_close(self.contacts.total_weight(), total, 1.0e-12 * total);
        assert_eq!(self.contacts.cutoff(), Quantity::new(CUTOFF, NANOMETER));
        assert_eq!(self.contacts.pairs(), contact_pairs(set.topology()));
        for (weight, (first, second)) in self.contacts.contact_weights().iter().zip(CONTACT_PAIRS) {
            let expected = observed
                .iter()
                .filter(|(geometry, _)| distance(*geometry, first, second) <= CUTOFF)
                .map(|(_, weight)| weight)
                .sum::<f64>();
            assert_close(*weight, expected, 1.0e-12 * total);
        }

        assert_eq!(self.groups.observation_count(), observed.len());
        assert_close(self.groups.total_weight(), total, 1.0e-12 * total);
        assert_eq!(self.groups.pairs(), GROUP_PAIRS);
        for (weight, (first, second)) in self.groups.contact_weights().iter().zip(GROUP_PAIRS) {
            // Two groups touch when some atom of one is within the cutoff of
            // some atom of the other.
            let touching = |geometry: usize| {
                GROUPS[first].iter().any(|a| {
                    GROUPS[second]
                        .iter()
                        .any(|b| distance(geometry, *a, *b) <= CUTOFF)
                })
            };
            let expected = observed
                .iter()
                .filter(|(geometry, _)| touching(*geometry))
                .map(|(_, weight)| weight)
                .sum::<f64>();
            assert_close(*weight, expected, 1.0e-12 * total);
        }
    }
}

/// Each geometry with weight one.
fn unweighted(count: usize) -> Vec<(usize, f64)> {
    (0..count).map(|geometry| (geometry, 1.0)).collect()
}

/// The fit and measured atoms, pair weights, and options of the alignment
/// checks.
struct AlignmentCase {
    fit: AtomSelection,
    measured: AtomSelection,
}

const PAIR_WEIGHTS: [f64; 4] = [1.0, 2.0, 0.5, 3.0];

impl AlignmentCase {
    fn new(topology: &Arc<Topology>) -> Self {
        Self {
            fit: selection(topology, &[0, 1, 2, 3]),
            measured: selection(topology, &[2, 3, 4, 5]),
        }
    }

    fn options(&self) -> AlignmentOptions<'static> {
        AlignmentOptions::default().with_weighting(Weighting::Explicit(&PAIR_WEIGHTS))
    }

    fn aligned(&self) -> AlignedRmsdOptions<'static> {
        AlignedRmsdOptions::default()
            .with_fit(self.options())
            .with_measurement_weighting(Weighting::Explicit(&PAIR_WEIGHTS))
    }
}

/// The weighted RMSD over `atoms` between two lists of points, as stored:
/// `sqrt(sum(w_i |x_i - y_i|^2) / sum(w_i))`.
fn direct_rmsd(moving: &[Point3], reference: &[Point3], atoms: &[usize]) -> f64 {
    let total = PAIR_WEIGHTS.iter().sum::<f64>();
    let sum = atoms
        .iter()
        .zip(PAIR_WEIGHTS)
        .map(|(atom, weight)| {
            let difference = moving[*atom] - reference[*atom];
            weight * difference.dot(difference)
        })
        .sum::<f64>();
    (sum / total).sqrt()
}

/// Measures and superposes `fixture` against its geometry `index`, whose
/// points `reference` holds, or else against `reference` itself, and checks
/// every value against the definitions and against fits and measurements of
/// each geometry on its own.
fn check_alignments(fixture: &StructureSet, index: Option<usize>, reference: &Structure) {
    let topology = fixture.shared_topology();
    let case = AlignmentCase::new(&topology);
    let reference_points = reference.positions().values().value().to_vec();
    let mut set = fixture.clone();
    let (rmsd, aligned, report) = match index {
        Some(index) => (
            set.rmsd(index, &case.fit, case.options()).unwrap(),
            set.aligned_rmsd(index, &case.fit, &case.measured, case.aligned())
                .unwrap(),
            set.superpose(index, &case.fit, case.options()).unwrap(),
        ),
        None => (
            set.rmsd(reference, &case.fit, case.options()).unwrap(),
            set.aligned_rmsd(reference, &case.fit, &case.measured, case.aligned())
                .unwrap(),
            set.superpose(reference, &case.fit, case.options()).unwrap(),
        ),
    };
    let view = reference.as_structure_view();
    assert_eq!(rmsd.unit(), NANOMETER);
    assert_eq!(aligned.unit(), NANOMETER);
    assert_eq!(rmsd.value().len(), fixture.len());
    assert_eq!(aligned.value().len(), fixture.len());
    assert_eq!(report.len(), fixture.len());
    assert_eq!(report.reference_index(), index);
    for geometry in 0..fixture.len() {
        let before = nanometers(fixture, geometry);
        let moving = fixture.get(geometry).unwrap();
        let moving: StructureView<'_> = moving.as_structure_view();
        // The RMSD is measured as stored, without fitting.
        let direct = direct_rmsd(&before, &reference_points, &[0, 1, 2, 3]);
        assert_close(rmsd.value()[geometry], direct, 1.0e-12);
        assert_eq!(
            rmsd.value()[geometry],
            alignment::rmsd(moving, view, &case.fit, case.options())
                .unwrap()
                .into_value()
        );
        // The aligned RMSD is the free kernel's, geometry by geometry.
        assert_eq!(
            aligned.value()[geometry],
            alignment::aligned_rmsd(moving, view, &case.fit, &case.measured, case.aligned())
                .unwrap()
                .into_value()
        );
        // Superposing applies each geometry's own fit.
        let applied = report.alignment(geometry).unwrap();
        let fit = alignment::fit(moving, view, &case.fit, case.options()).unwrap();
        assert_eq!(applied, &fit);
        let transform = applied.transform();
        let after = nanometers(&set, geometry);
        let moved = before
            .iter()
            .map(|point| transform.transform_point(*point))
            .collect::<Vec<_>>();
        assert_points_close(&after, &moved, 1.0e-12);
        // A fit cannot be farther than the stored coordinates, and once
        // applied its RMSD is measured as stored; so is the aligned RMSD
        // over the measured atoms.
        let fitted = applied.rmsd().into_value();
        assert!(fitted <= direct + 1.0e-12, "{fitted} > {direct}");
        assert_close(
            direct_rmsd(&after, &reference_points, &[0, 1, 2, 3]),
            fitted,
            1.0e-12,
        );
        assert_close(
            direct_rmsd(&after, &reference_points, &[2, 3, 4, 5]),
            aligned.value()[geometry],
            1.0e-12,
        );
        // Velocities rotate with the frames.
        if fixture.get(geometry).unwrap().velocities().is_some() {
            let rotated = velocities(fixture, geometry)
                .into_iter()
                .map(|velocity| transform.transform_vector(velocity))
                .collect::<Vec<_>>();
            for (actual, expected) in velocities(&set, geometry).iter().zip(rotated) {
                assert!((*actual - expected).norm() <= 1.0e-12);
            }
        }
    }
    if let Some(index) = index {
        assert_eq!(rmsd.value()[index], 0.0);
        assert_close(aligned.value()[index], 0.0, 1.0e-12);
        assert_close(
            report.alignment(index).unwrap().rmsd().into_value(),
            0.0,
            1.0e-12,
        );
    }
    // The order is kept.
    assert_eq!(set.order(), fixture.order());
}

#[test]
fn structure_set_alignments_match_each_geometry_fitted_alone() {
    let structure = Structure::new(topology(), positions(points(11))).unwrap();
    // Velocities rotate with the frames, and the order is kept.
    for fixture in [ensemble(true).set().clone(), trajectory(), series()] {
        // A geometry of the set, and an independent view.
        let own = Structure::new(topology(), positions(points(3))).unwrap();
        check_alignments(&fixture, Some(3), &own);
        check_alignments(&fixture, None, &structure);

        // A reference out of range changes nothing.
        let mut set = fixture.clone();
        let case = AlignmentCase::new(&set.shared_topology());
        let out_of_range = AlignmentError::ReferenceOutOfRange {
            index: 9,
            len: fixture.len(),
        };
        assert_eq!(
            set.rmsd(9, &case.fit, case.options()),
            Err(out_of_range.clone())
        );
        assert_eq!(
            set.aligned_rmsd(9, &case.fit, &case.measured, case.aligned()),
            Err(out_of_range.clone())
        );
        assert_eq!(
            set.superpose(9, &case.fit, case.options()),
            Err(out_of_range)
        );
        assert_eq!(geometries(&set), geometries(&fixture));
    }
}

#[test]
fn structure_set_alignment_covers_every_geometry_and_keeps_the_order() {
    // A fit moves every geometry, keeps the order, and leaves a distribution
    // over part of the set fitting the set.
    let mut set = trajectory();
    let distribution =
        Distribution::new(Samples::Subset(vec![1, 4]), Partition::state_per_sample(2))
            .unwrap()
            .with_log_weights(ln(&[2.0, 1.0]))
            .unwrap();
    let atoms = all(set.shared_topology());
    let order = set.order().cloned();
    let report = set.superpose(0, &atoms, Default::default()).unwrap();
    assert_eq!(report.len(), 8);
    assert_eq!(set.order().cloned(), order);
    let rmsd = set.rmsd(0, &atoms, Default::default()).unwrap();
    assert_eq!(rmsd.value().len(), 8);
    let ensemble = Ensemble::new(set.clone(), distribution.clone()).unwrap();
    assert_eq!(ensemble.distribution(), &distribution);

    // A geometry no fit succeeds to fails the fit and changes nothing.
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

/// Each sample's geometry index with its probability: its weight divided by
/// the total.
fn probabilities(geometries: &[usize], weights: &[f64]) -> Vec<(usize, f64)> {
    let total = weights.iter().sum::<f64>();
    geometries
        .iter()
        .zip(weights)
        .map(|(geometry, weight)| (*geometry, weight / total))
        .collect()
}

#[test]
fn set_statistics_weigh_every_geometry_and_ensemble_statistics_their_samples() {
    // The ensembles observe each geometry with its probability; energies do
    // not enter.
    let samples = probabilities(&(0..8).collect::<Vec<_>>(), &WEIGHTS);
    let with_energies = Statistics::of_ensemble(&ensemble(true));
    with_energies.check(ensemble(true).set(), &samples);
    assert!(with_energies.values() == Statistics::of_ensemble(&ensemble(false)).values());

    // The set of an ensemble, on its own, weighs every geometry one.
    let set = ensemble(true).set().clone();
    let statistics = Statistics::of_set(&set);
    statistics.check(&set, &unweighted(8));
    assert!(statistics.values() != with_energies.values());

    // So does a trajectory, whatever its order and velocities.
    let set = trajectory();
    let statistics = Statistics::of_set(&set);
    statistics.check(&set, &unweighted(8));
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

/// The samples of `ensemble` as an ensemble of their own, unordered, with
/// their weights and states, and without velocities.
fn samples_alone(ensemble: &Ensemble) -> Ensemble {
    let samples = ensemble
        .distribution()
        .sample_probabilities()
        .unwrap()
        .into_iter()
        .map(|(geometry, _)| geometry)
        .collect::<Vec<_>>();
    let (mut set, distribution) = ensemble.select(samples).unwrap().into_parts();
    set.set_order(None).unwrap();
    for index in 0..set.len() {
        set.get_mut(index).unwrap().take_velocities();
    }
    Ensemble::new(set, distribution).unwrap()
}

/// The samples of [`equilibrated`], and their weights.
const EQUILIBRATED: [usize; 5] = [2, 3, 5, 6, 7];
const EQUILIBRATED_WEIGHTS: [f64; 5] = [1.0, 2.0, 3.0, 0.5, 1.5];

/// A time-ordered set of eight frames whose distribution covers frames 2,
/// 3, 5, 6, and 7, weighted, in two states with energies.
fn equilibrated() -> Ensemble {
    let states = [0.0, 2.5].map(|energy| {
        StateData::new()
            .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
            .unwrap()
    });
    let distribution = Distribution::new(
        Samples::Subset(EQUILIBRATED.to_vec()),
        Partition::new(vec![0, 1, 0, 1, 1], vec![0, 1]).unwrap(),
    )
    .unwrap()
    .with_log_weights(ln(&EQUILIBRATED_WEIGHTS))
    .unwrap()
    .with_states(states.to_vec())
    .unwrap();
    Ensemble::new(trajectory(), distribution).unwrap()
}

#[test]
fn statistics_over_part_of_a_time_ordered_set_observe_its_samples_alone() {
    let ensemble = equilibrated();
    let observed = Statistics::of_ensemble(&ensemble);
    observed.check(
        ensemble.set(),
        &probabilities(&EQUILIBRATED, &EQUILIBRATED_WEIGHTS),
    );
    assert!(observed.values() == Statistics::of_ensemble(&samples_alone(&ensemble)).values());
    let rmsf = ensemble
        .rmsf(&all(ensemble.set().shared_topology()))
        .unwrap();
    assert_eq!(rmsf.observation_count(), 5);

    // Moving a geometry that is not a sample changes nothing.
    let mut moved = ensemble.set().clone();
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
    let moved = Ensemble::new(moved, ensemble.distribution().clone()).unwrap();
    assert!(Statistics::of_ensemble(&moved).values() == observed.values());

    // A distribution without weights observes each sample with probability
    // one in five.
    let mut equal = ensemble.clone();
    equal.distribution_mut().set_log_weights(None).unwrap();
    let topology = ensemble.set().shared_topology();
    let rmsf = equal.rmsf(&all(Arc::clone(&topology))).unwrap();
    assert!((rmsf.total_weight() - 1.0).abs() < 1.0e-15);
    let cutoff = Quantity::new(CUTOFF, NANOMETER);
    let contacts = equal
        .contact_occupancy(contact_pairs(&topology), cutoff)
        .unwrap();
    let frames = ensemble.set().select(EQUILIBRATED).unwrap();
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
    Statistics::of_ensemble(&equal).check(ensemble.set(), &probabilities(&EQUILIBRATED, &[1.0; 5]));
}

/// What clustering a set or an ensemble needs: its members' geometries,
/// their linear weights, and the clustering methods.
trait Clustered {
    fn set(&self) -> &StructureSet;
    /// The geometry index of each member, in member order.
    fn members(&self) -> Vec<usize>;
    /// Each member's linear weight, up to a common factor.
    fn weights(&self) -> Vec<f64>;
    fn by_radius(
        &self,
        radius: f64,
        atoms: &AtomSelection,
        options: ClusteringOptions<'_>,
    ) -> Result<Partition, ClusteringError>;
    fn matrix(&self, atoms: &AtomSelection) -> DistanceMatrix;
    fn hierarchical(&self, matrix: &DistanceMatrix, linkage: Linkage, cut: Cut) -> Partition;
    fn k_medoids(
        &self,
        matrix: &DistanceMatrix,
        clusters: usize,
    ) -> Result<Partition, ClusteringError>;
}

impl Clustered for StructureSet {
    fn set(&self) -> &StructureSet {
        self
    }

    fn members(&self) -> Vec<usize> {
        (0..self.len()).collect()
    }

    fn weights(&self) -> Vec<f64> {
        vec![1.0; self.len()]
    }

    fn by_radius(
        &self,
        radius: f64,
        atoms: &AtomSelection,
        options: ClusteringOptions<'_>,
    ) -> Result<Partition, ClusteringError> {
        self.cluster_by_radius(Quantity::new(radius, NANOMETER), atoms, options)
    }

    fn matrix(&self, atoms: &AtomSelection) -> DistanceMatrix {
        self.pairwise_fitted_rmsd(atoms, Default::default())
            .unwrap()
    }

    fn hierarchical(&self, matrix: &DistanceMatrix, linkage: Linkage, cut: Cut) -> Partition {
        self.cluster_hierarchical(matrix, linkage, cut).unwrap()
    }

    fn k_medoids(
        &self,
        matrix: &DistanceMatrix,
        clusters: usize,
    ) -> Result<Partition, ClusteringError> {
        self.cluster_k_medoids(matrix, clusters)
    }
}

impl Clustered for Ensemble {
    fn set(&self) -> &StructureSet {
        Ensemble::set(self)
    }

    fn members(&self) -> Vec<usize> {
        self.distribution()
            .sample_probabilities()
            .unwrap()
            .into_iter()
            .map(|(geometry, _)| geometry)
            .collect()
    }

    fn weights(&self) -> Vec<f64> {
        self.distribution()
            .sample_probabilities()
            .unwrap()
            .into_iter()
            .map(|(_, probability)| probability)
            .collect()
    }

    fn by_radius(
        &self,
        radius: f64,
        atoms: &AtomSelection,
        options: ClusteringOptions<'_>,
    ) -> Result<Partition, ClusteringError> {
        self.cluster_by_radius(Quantity::new(radius, NANOMETER), atoms, options)
    }

    fn matrix(&self, atoms: &AtomSelection) -> DistanceMatrix {
        self.pairwise_fitted_rmsd(atoms, Default::default())
            .unwrap()
    }

    fn hierarchical(&self, matrix: &DistanceMatrix, linkage: Linkage, cut: Cut) -> Partition {
        self.cluster_hierarchical(matrix, linkage, cut).unwrap()
    }

    fn k_medoids(
        &self,
        matrix: &DistanceMatrix,
        clusters: usize,
    ) -> Result<Partition, ClusteringError> {
        self.cluster_k_medoids(matrix, clusters)
    }
}

/// The fitted RMSD in nanometres between two geometries of `set`, the first
/// moving.
fn fitted(set: &StructureSet, moving: usize, reference: usize, atoms: &AtomSelection) -> f64 {
    alignment::fit(
        set.get(moving).unwrap().as_structure_view(),
        set.get(reference).unwrap().as_structure_view(),
        atoms,
        Default::default(),
    )
    .unwrap()
    .rmsd()
    .into_value()
}

/// Leader clustering as documented, one fit at a time: members in
/// `visiting` order, by member position, join the nearest representative
/// within `radius`, the first on ties, or found a cluster.
fn naive_leaders(
    subject: &impl Clustered,
    visiting: &[usize],
    radius: f64,
    atoms: &AtomSelection,
) -> Partition {
    let members = subject.members();
    let mut assignments = vec![0; members.len()];
    let mut representatives: Vec<usize> = Vec::new();
    for &member in visiting {
        let mut nearest: Option<(usize, f64)> = None;
        for (cluster, representative) in representatives.iter().enumerate() {
            let rmsd = fitted(
                subject.set(),
                members[member],
                members[*representative],
                atoms,
            );
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
    Partition::new(assignments, representatives).unwrap()
}

/// Member positions by descending `keys`, ties in member order.
fn descending(keys: &[f64]) -> Vec<usize> {
    let mut order = (0..keys.len()).collect::<Vec<_>>();
    order.sort_by(|a, b| keys[*b].total_cmp(&keys[*a]));
    order
}

/// The weighted summed distance from `candidate` to the members of
/// `cluster`.
fn weighted_cost(
    matrix: &DistanceMatrix,
    weights: &[f64],
    cluster: &[usize],
    candidate: usize,
) -> f64 {
    cluster
        .iter()
        .map(|other| weights[*other] * matrix.get(candidate, *other).unwrap().into_value())
        .sum()
}

/// The member of `cluster` with the least weighted summed distance to the
/// others, the first on ties: the medoid that represents a cluster.
fn weighted_medoid(matrix: &DistanceMatrix, weights: &[f64], cluster: &[usize]) -> usize {
    let cost = |candidate: usize| weighted_cost(matrix, weights, cluster, candidate);
    let mut best = cluster[0];
    for &candidate in &cluster[1..] {
        if cost(candidate) < cost(best) {
            best = candidate;
        }
    }
    best
}

/// The samples of each state of `partition`.
fn clusters_of(partition: &Partition) -> Vec<Vec<usize>> {
    (0..partition.state_count())
        .map(|state| partition.samples(state).collect())
        .collect()
}

/// Members grouped by `shapes`, the shape of each member, as two states.
fn by_shape(partition: &Partition, shapes: &[usize]) {
    assert_eq!(partition.state_count(), 2);
    for (position, state) in partition.assignments().iter().enumerate() {
        let representative = partition.representatives()[*state];
        assert_eq!(shapes[position], shapes[representative], "{partition:?}");
    }
}

const RADII: [f64; 5] = [0.0, 0.05, 0.08, 0.2, 10.0];

/// Every clustering of `subject`: its fitted distance matrix, clusters of
/// [`RADII`] in each visiting order of `orders`, hierarchical clusters of
/// every count and one height under three linkages, and k-medoids of every
/// count; each checked against its definition. `shapes` is the shape of
/// each member, which two clusters separate.
fn clusterings(
    subject: &impl Clustered,
    shapes: &[usize],
    orders: &[(ClusterOrder, Vec<usize>)],
) -> impl Debug + PartialEq {
    let atoms = all(subject.set().shared_topology());
    let members = subject.members();
    let weights = subject.weights();
    let count = members.len();

    // Each entry is the fitted RMSD of its pair.
    let matrix = subject.matrix(&atoms);
    assert_eq!(matrix.member_count(), count);
    for first in 0..count {
        assert_eq!(matrix.get(first, first).unwrap().into_value(), 0.0);
        for second in first + 1..count {
            let fit = alignment::fit(
                subject
                    .set()
                    .get(members[first])
                    .unwrap()
                    .as_structure_view(),
                subject
                    .set()
                    .get(members[second])
                    .unwrap()
                    .as_structure_view(),
                &atoms,
                Default::default(),
            )
            .unwrap();
            assert_eq!(matrix.get(first, second), Some(fit.rmsd()));
            assert_eq!(matrix.get(second, first), Some(fit.rmsd()));
        }
    }
    let mut radii = Vec::new();
    for (order, visiting) in orders {
        for radius in RADII {
            let options = ClusteringOptions::default().with_order(*order);
            let partition = subject.by_radius(radius, &atoms, options).unwrap();
            assert_eq!(
                partition,
                naive_leaders(subject, visiting, radius, &atoms),
                "{order:?} {radius}"
            );
            radii.push(partition);
        }
        // Within no radius every member founds its own cluster, in visiting
        // order; within 0.08 nm the first of each shape leads its shape; within
        // 10 nm the first member leads one cluster.
        let alone = &radii[radii.len() - 5];
        assert_eq!(alone.representatives(), visiting);
        let shapes_at_radius = &radii[radii.len() - 3];
        by_shape(shapes_at_radius, shapes);
        assert_eq!(shapes_at_radius.representatives()[0], visiting[0]);
        let one = &radii[radii.len() - 1];
        assert_eq!(one.assignments(), vec![0; count]);
        assert_eq!(one.representatives(), [visiting[0]]);
    }

    let mut hierarchical = Vec::new();
    for linkage in [Linkage::Single, Linkage::Complete, Linkage::Average] {
        let dendrogram = matrix.dendrogram(linkage);
        let cuts = (1..=count)
            .map(Cut::Count)
            .chain([Cut::Height(Quantity::new(0.06, NANOMETER))]);
        for cut in cuts {
            let partition = subject.hierarchical(&matrix, linkage, cut);
            // The dendrogram's cut, each cluster represented by its weighted
            // medoid.
            assert_eq!(partition.assignments(), dendrogram.cut(cut).unwrap());
            if let Cut::Count(clusters) = cut {
                assert_eq!(partition.state_count(), clusters);
            }
            for (cluster, representative) in clusters_of(&partition)
                .iter()
                .zip(partition.representatives())
            {
                assert_eq!(*representative, weighted_medoid(&matrix, &weights, cluster));
            }
            if matches!(cut, Cut::Count(2)) {
                by_shape(&partition, shapes);
            }
            hierarchical.push(partition);
        }
    }

    let medoids = (1..=count)
        .map(|clusters| {
            let partition = subject.k_medoids(&matrix, clusters).unwrap();
            assert_eq!(partition.state_count(), clusters);
            let representatives = partition.representatives();
            // Each member joins its nearest medoid, and each medoid is the
            // weighted medoid of its cluster.
            for (member, state) in partition.assignments().iter().enumerate() {
                let distance = |state: usize| {
                    matrix
                        .get(member, representatives[state])
                        .unwrap()
                        .into_value()
                };
                assert!((0..clusters).all(|other| distance(*state) <= distance(other)));
            }
            for (cluster, representative) in clusters_of(&partition).iter().zip(representatives) {
                // Ties may keep either medoid: no swap lowers the total.
                let cost = |candidate: usize| weighted_cost(&matrix, &weights, cluster, candidate);
                assert!(cluster
                    .iter()
                    .all(|other| cost(*representative) <= cost(*other)));
            }
            partition
        })
        .collect::<Vec<_>>();
    // One medoid is the weighted medoid of every member, two split the
    // shapes, and one per member leaves each alone.
    let everyone = (0..count).collect::<Vec<_>>();
    assert_eq!(
        medoids[0].representatives(),
        [weighted_medoid(&matrix, &weights, &everyone)]
    );
    by_shape(&medoids[1], shapes);
    assert_eq!(medoids[count - 1].assignments(), everyone);
    assert_eq!(medoids[count - 1].representatives(), everyone);
    (matrix, radii, hierarchical, medoids)
}

/// Member positions by ascending energy of their states, ties in member
/// order.
fn by_energy(energies: &[f64]) -> Vec<usize> {
    let mut order = (0..energies.len()).collect::<Vec<_>>();
    order.sort_by(|a, b| energies[*a].total_cmp(&energies[*b]));
    order
}

#[test]
fn structure_set_and_ensemble_clusterings_follow_their_definitions() {
    // By weight 4, 3, 2, 1.5, then 1 twice in member order, 0.5, 0.25; by
    // state energy 0 (geometries 2, 3), 0.7 (7), 1.5 (0, 1), 4 (4, 5, 6).
    let weight_order = descending(&WEIGHTS);
    assert_eq!(weight_order, [5, 1, 3, 4, 0, 7, 2, 6]);
    let energies = STATES.map(|state| ENERGIES[state]);
    let energy_order = by_energy(&energies);
    assert_eq!(energy_order, [2, 3, 7, 0, 1, 4, 5, 6]);
    let parity = (0..8).map(|geometry| geometry % 2).collect::<Vec<_>>();
    clusterings(
        &ensemble(true),
        &parity,
        &[
            (ClusterOrder::Weight, weight_order.clone()),
            (ClusterOrder::Energy, energy_order),
        ],
    );
    clusterings(
        &ensemble(false),
        &parity,
        &[(ClusterOrder::Weight, weight_order)],
    );

    // Without energies, ordering by energy is refused.
    let set = ensemble(false);
    let atoms = all(set.set().shared_topology());
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
    assert_eq!(
        set.cluster_k_medoids(&matrix, 0),
        Err(ClusteringError::InvalidClusterCount {
            requested: 0,
            members: 8
        })
    );
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
    assert_eq!(refusal.to_string(), "Cannot cluster the geometries.");
    assert_eq!(
        std::error::Error::source(&refusal).unwrap().to_string(),
        "The structure set needs an energy for every state."
    );

    // The set of an ensemble, and a trajectory's frames, cluster every
    // geometry equally weighted, visited in geometry order.
    let geometry_order = (0..8).collect::<Vec<_>>();
    let own = clusterings(
        ensemble(true).set(),
        &parity,
        &[(ClusterOrder::Weight, geometry_order.clone())],
    );
    let set = trajectory();
    assert!(own == clusterings(&set, &parity, &[(ClusterOrder::Weight, geometry_order)]));
    // A set has no states, so no energies to order by.
    assert_eq!(
        set.cluster_by_radius(radius, &atoms, by_energy),
        Err(ClusteringError::StructureSet(Box::new(
            StructureSetError::Missing(Requirement::StateEnergies)
        )))
    );
}

#[test]
fn clustering_part_of_a_time_ordered_set_groups_its_samples_by_position() {
    let set = equilibrated();
    // Samples 2, 3, 5, 6, 7 weigh 1, 2, 3, 0.5, 1.5 and lie in states of
    // energies 0, 2.5, 0, 2.5, 2.5.
    let weight_order = descending(&EQUILIBRATED_WEIGHTS);
    assert_eq!(weight_order, [2, 1, 4, 0, 3]);
    let energy_order = by_energy(&[0.0, 2.5, 0.0, 2.5, 2.5]);
    assert_eq!(energy_order, [0, 2, 1, 3, 4]);
    let orders = [
        (ClusterOrder::Weight, weight_order),
        (ClusterOrder::Energy, energy_order),
    ];
    let parity = EQUILIBRATED.map(|geometry| geometry % 2);
    let clustered = clusterings(&set, &parity, &orders);
    assert!(clustered == clusterings(&samples_alone(&set), &parity, &orders));
    let atoms = all(set.set().shared_topology());
    let matrix = set
        .pairwise_fitted_rmsd(&atoms, Default::default())
        .unwrap();
    assert_eq!(matrix.member_count(), 5);
    // Entry (0, 1) is between the first two samples, geometries 2 and 3.
    let fit = alignment::fit(
        set.set().get(2).unwrap().as_structure_view(),
        set.set().get(3).unwrap().as_structure_view(),
        &atoms,
        Default::default(),
    )
    .unwrap();
    let distance = matrix.get(0, 1).unwrap().into_value();
    assert!((distance - fit.rmsd().into_value()).abs() < 1.0e-12);

    // Even and odd geometries have different shapes: two medoids group the
    // samples by shape, each sample by its position.
    let states = set.cluster_k_medoids(&matrix, 2).unwrap();
    let shapes = EQUILIBRATED.map(|geometry| geometry % 2);
    for (position, state) in states.assignments().iter().enumerate() {
        let representative = states.representatives()[*state];
        assert_eq!(shapes[position], shapes[representative]);
    }

    // The partition installs as the distribution's grouping, keeping every
    // weight, sample, and the order.
    let mut regrouped = set.clone();
    regrouped
        .distribution_mut()
        .regroup(states.clone())
        .unwrap();
    let distribution = regrouped.distribution();
    assert_eq!(distribution.partition(), &states);
    assert_eq!(distribution.log_weights(), set.distribution().log_weights());
    assert_eq!(
        distribution.samples(),
        &Samples::Subset(EQUILIBRATED.to_vec())
    );
    assert_eq!(regrouped.set().order(), set.set().order());
    assert!(distribution.state_of(4).is_none());

    // Geometries that are not samples are not clustered, even when no fit to
    // them succeeds.
    let broken = |geometries: &[usize]| {
        let mut broken = set.set().clone();
        for geometry in geometries {
            broken
                .get_mut(*geometry)
                .unwrap()
                .set_positions(Quantity::new(collinear(), NANOMETER))
                .unwrap();
        }
        Ensemble::new(broken, set.distribution().clone()).unwrap()
    };
    assert_eq!(
        broken(&[4]).pairwise_fitted_rmsd(&atoms, Default::default()),
        Ok(matrix.clone())
    );
    // A failing sample is named by its geometry index.
    let broken = broken(&[4, 5]);
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
        "Cannot fit geometries 2 and 5."
    );
    // The heaviest sample, geometry 5, leads; geometry 3, visited next, fails
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
fn ensemble_clustering_runs_on_threads_report_progress_and_cancel() {
    let set = equilibrated();
    let atoms = all(set.set().shared_topology());
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
    // Five samples make ten pairs.
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
fn a_set_installs_its_clusters_as_the_distribution_of_an_ensemble() {
    let set = trajectory();
    let atoms = all(set.shared_topology());
    let matrix = set
        .pairwise_fitted_rmsd(&atoms, Default::default())
        .unwrap();
    assert_eq!(matrix.member_count(), 8);
    let states = set
        .cluster_hierarchical(&matrix, Linkage::Average, Cut::Count(2))
        .unwrap();
    let distribution = Distribution::new(Samples::All, states.clone()).unwrap();
    let ensemble = Ensemble::new(set, distribution).unwrap();
    assert_eq!(ensemble.distribution().partition(), &states);
    assert_eq!(ensemble.states().len(), 2);
    // Frames of one shape share a state.
    for geometry in 0..8 {
        let state = ensemble.distribution().state_of(geometry).unwrap();
        let representative = ensemble.state(state).unwrap().representative();
        assert_eq!(geometry % 2, representative % 2);
    }
}

/// Three geometries whose weights lie near the floating-point limits, in two
/// states with energies.
fn extreme() -> Ensemble {
    let states = [0.0, 1.0].map(|energy| {
        StateData::new()
            .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
            .unwrap()
    });
    let distribution = Distribution::new(
        Samples::All,
        Partition::new(vec![0, 1, 1], vec![0, 1]).unwrap(),
    )
    .unwrap()
    .with_log_weights(ln(&[1.0e300, 1.0e-300, 2.0e-300]))
    .unwrap()
    .with_states(states.to_vec())
    .unwrap();
    Ensemble::new(unordered(3), distribution).unwrap()
}

/// Checks Boltzmann populations: each state's probability is
/// `exp(-E / RT)` over the sum of those factors, and its samples share it in
/// proportion to their `prior` linear weights; the total weight is kept.
fn check_boltzmann(ensemble: &Ensemble, kelvin: f64, energies: &[f64], prior_log_weights: &[f64]) {
    let distribution = ensemble.distribution();
    let partition = distribution.partition();
    let thermal = GAS_CONSTANT * kelvin;
    // In log space relative to the lowest energy, so nothing underflows.
    let minimum = energies.iter().copied().fold(f64::INFINITY, f64::min);
    let factors = energies
        .iter()
        .map(|energy| (-(energy - minimum) / thermal).exp())
        .collect::<Vec<_>>();
    let factor_total = factors.iter().sum::<f64>();
    let state_probabilities = factors
        .iter()
        .map(|factor| factor / factor_total)
        .collect::<Vec<_>>();
    for (actual, expected) in distribution
        .state_probabilities()
        .unwrap()
        .iter()
        .zip(&state_probabilities)
    {
        assert_close(*actual, *expected, 1.0e-12);
    }
    // Each sample's share of its state, relative to the heaviest sample of
    // the state, so that no share underflows.
    let shares = (0..prior_log_weights.len())
        .map(|position| {
            let state = partition.assignments()[position];
            let heaviest = partition
                .samples(state)
                .map(|other| prior_log_weights[other])
                .fold(f64::NEG_INFINITY, f64::max);
            let total = partition
                .samples(state)
                .map(|other| (prior_log_weights[other] - heaviest).exp())
                .sum::<f64>();
            (prior_log_weights[position] - heaviest).exp() / total
        })
        .collect::<Vec<_>>();
    for (position, probability) in distribution.probabilities().unwrap().iter().enumerate() {
        let state = partition.assignments()[position];
        let expected = state_probabilities[state] * shares[position];
        assert_close(*probability, expected, 1.0e-12);
    }
    // The total weight is kept.
    let log_sum = |values: &[f64]| {
        let maximum = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        maximum
            + values
                .iter()
                .map(|value| (value - maximum).exp())
                .sum::<f64>()
                .ln()
    };
    assert_close(
        log_sum(distribution.log_weights().unwrap()),
        log_sum(prior_log_weights),
        1.0e-9 * log_sum(prior_log_weights).abs().max(1.0),
    );
    assert_eq!(
        distribution.temperature(),
        Some(Quantity::new(kelvin, KELVIN))
    );
}

#[test]
fn boltzmann_populations_follow_the_energies_and_keep_shares_within_states() {
    for kelvin in [1.0, 298.15, 5000.0] {
        let mut ensemble = ensemble(true);
        ensemble
            .distribution_mut()
            .assign_boltzmann_populations(Quantity::new(kelvin, KELVIN))
            .unwrap();
        check_boltzmann(&ensemble, kelvin, &ENERGIES, &ln(&WEIGHTS));
    }
    // At room temperature the lowest state, 0 kJ/mol, is the most
    // populated, then 0.7, 1.5, and 4 kJ/mol.
    let mut ensemble = ensemble(true);
    ensemble
        .distribution_mut()
        .assign_boltzmann_populations(Quantity::new(298.15, KELVIN))
        .unwrap();
    assert_eq!(ensemble.distribution().ranked_by_population(), [1, 3, 0, 2]);

    // Weights near the floating-point limits keep their ratios within a
    // state: geometry 2 weighs twice geometry 1.
    let mut set = extreme();
    set.distribution_mut()
        .assign_boltzmann_populations(Quantity::new(300.0, KELVIN))
        .unwrap();
    check_boltzmann(
        &set,
        300.0,
        &[0.0, 1.0],
        &ln(&[1.0e300, 1.0e-300, 2.0e-300]),
    );
    let log_weights = set.distribution().log_weights().unwrap();
    assert_close(log_weights[2] - log_weights[1], 2.0_f64.ln(), 1.0e-9);

    // A distribution over part of a set without weights gains them, equal
    // within each state.
    let mut part = equilibrated();
    part.distribution_mut().set_log_weights(None).unwrap();
    part.distribution_mut()
        .assign_boltzmann_populations(Quantity::new(300.0, KELVIN))
        .unwrap();
    check_boltzmann(&part, 300.0, &[0.0, 2.5], &[0.0; 5]);
    let log_weights = part.distribution().log_weights().unwrap();
    assert_eq!(log_weights.len(), 5);
    assert_eq!(log_weights[0], log_weights[2]);
    assert_eq!(log_weights[1], log_weights[3]);
    assert_eq!(log_weights[1], log_weights[4]);
}

#[test]
fn boltzmann_populations_refuse_without_energies_and_change_nothing() {
    let mut ensemble = ensemble(false);
    let before = ensemble.distribution().clone();
    let mut distribution = ensemble.distribution_mut();
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
    assert_eq!(ensemble.distribution(), &before);
    assert_eq!(
        StructureSetError::Missing(Requirement::StateEnergies).to_string(),
        "The structure set needs an energy for every state."
    );

    // A state far above the others keeps a finite log population, lower by
    // its energy over RT, though its probability underflows.
    let mut ensemble = self::ensemble(false);
    let energies = [0.0, 1.0e6, 2.0, 3.0];
    let mut distribution = ensemble.distribution_mut();
    for (state, energy) in energies.into_iter().enumerate() {
        *distribution.state_mut(state).unwrap() = StateData::new()
            .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
            .unwrap();
    }
    distribution
        .assign_boltzmann_populations(Quantity::new(1.0, KELVIN))
        .unwrap();
    check_boltzmann(&ensemble, 1.0, &energies, &ln(&WEIGHTS));
    let distribution = ensemble.distribution();
    assert!(distribution
        .log_weights()
        .unwrap()
        .iter()
        .all(|w| w.is_finite()));
    let gap =
        ensemble.state(0).unwrap().log_population() - ensemble.state(1).unwrap().log_population();
    let expected = 1.0e6 / GAS_CONSTANT;
    assert_close(gap, expected, 1.0e-9 * expected);
    assert_eq!(distribution.state_probabilities().unwrap()[1], 0.0);
}
