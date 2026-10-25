//! Partial accumulators merged in order: the same result for one split
//! whatever runs the ranges, the serial result up to rounding, and
//! transactional failures.

use std::num::NonZeroUsize;
use std::ops::Range;
use std::sync::{Arc, Mutex};

use kekule::geometry::Point3;
use kekule::parallel::{Executor, ScopedThreads};
use kekule::statistics::{
    ContactOccupancyAccumulator, MergeMismatch, RmsfAccumulator, StatisticsError,
};
use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
use kekule::topology::{AtomSelection, InstanceAtomId, Topology};
use kekule::units::{Quantity, ANGSTROM, NANOMETER};

const ATOMS: usize = 8;
const OBSERVATIONS: usize = 240;

/// Common translations of the frames, and how far the RMSF of merged partials
/// may stray from the serial RMSF with each, relative to the serial RMSF of
/// about 0.2 nm. The largest errors seen are 7e-16 and 2e-13; rounding the
/// translated means costs about ε × 1000 nm / 0.2 nm.
const MERGED_RMSF_BOUNDS: [(f64, f64); 2] = [(0.0, 1e-14), (1000.0, 1e-11)];

fn chain_topology() -> Arc<Topology> {
    kekule::smiles::to_topology("CCCCCCCC").unwrap()
}

/// A xorshift generator of numbers in `[0, 1)`.
struct Numbers(u64);

impl Numbers {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Frames fluctuating by up to 0.4 nm about a lattice translated by
/// `offset`, with weights between 0.25 and 4.
fn weighted_frames(top: &Arc<Topology>, offset: f64) -> (Trajectory, Vec<f64>) {
    let mut numbers = Numbers(0x2545_f491_4f6c_dd1d);
    let mut frames = Vec::new();
    let mut weights = Vec::new();
    for _ in 0..OBSERVATIONS {
        let points = (0..ATOMS)
            .map(|atom| {
                let site = offset + 0.15 * atom as f64;
                Point3::new(
                    site + 0.4 * numbers.next(),
                    offset + 0.4 * numbers.next(),
                    offset + 0.4 * numbers.next(),
                )
            })
            .collect::<Vec<_>>();
        frames.push(TrajectoryFrame::new(
            Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
        ));
        weights.push(0.25 + 3.75 * numbers.next());
    }
    (
        Trajectory::from_items(Arc::clone(top), frames).unwrap(),
        weights,
    )
}

/// Splits of `0..OBSERVATIONS` into contiguous ranges, some of them empty or
/// of one observation.
fn splits() -> Vec<Vec<Range<usize>>> {
    let n = OBSERVATIONS;
    vec![
        vec![0..n],
        vec![0..n / 2, n / 2..n],
        vec![0..1, 1..17, 17..17, 17..n - 1, n - 1..n],
        (0..7).map(|k| k * n / 7..(k + 1) * n / 7).collect(),
        (0..n).map(|k| k..k + 1).collect(),
    ]
}

/// Observes each range into its own accumulator on `executor`, then merges
/// the partials in range order into a fresh accumulator.
fn merged<A: Send>(
    executor: &dyn Executor,
    ranges: &[Range<usize>],
    empty: impl Fn() -> A + Sync,
    observe: impl Fn(&mut A, usize) + Sync,
    merge: impl Fn(&mut A, &A),
) -> A {
    let partials = ranges.iter().map(|_| Mutex::new(None)).collect::<Vec<_>>();
    executor.run(ranges.len(), &|range| {
        let mut partial = empty();
        for index in ranges[range].clone() {
            observe(&mut partial, index);
        }
        *partials[range].lock().unwrap() = Some(partial);
    });
    let mut result = empty();
    for partial in partials {
        merge(&mut result, &partial.into_inner().unwrap().unwrap());
    }
    result
}

fn executors() -> Vec<ScopedThreads> {
    [1, 2, 3, 8]
        .map(|threads| ScopedThreads::new(NonZeroUsize::new(threads).unwrap()))
        .to_vec()
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

#[test]
fn merged_rmsf_depends_on_the_split_alone_and_matches_serial_observation() {
    let top = chain_topology();
    let all = AtomSelection::all(Arc::clone(&top));
    for (offset, bound) in MERGED_RMSF_BOUNDS {
        let (frames, weights) = weighted_frames(&top, offset);
        let mut serial = RmsfAccumulator::new(&all).unwrap();
        for (index, weight) in weights.iter().enumerate() {
            serial
                .observe(index, &frames.get(index).unwrap(), *weight)
                .unwrap();
        }
        let serial = serial.finish().unwrap();
        for ranges in splits() {
            let results = executors()
                .iter()
                .map(|executor| {
                    merged(
                        executor,
                        &ranges,
                        || RmsfAccumulator::new(&all).unwrap(),
                        |rmsf, index| {
                            rmsf.observe(index, &frames.get(index).unwrap(), weights[index])
                                .unwrap()
                        },
                        |rmsf, partial| rmsf.merge(partial).unwrap(),
                    )
                    .finish()
                    .unwrap()
                })
                .collect::<Vec<_>>();
            for result in &results {
                assert_eq!(
                    bits(result.values().value()),
                    bits(results[0].values().value())
                );
                assert_eq!(
                    result.total_weight().to_bits(),
                    results[0].total_weight().to_bits()
                );
            }
            let result = &results[0];
            assert_eq!(result.observation_count(), OBSERVATIONS);
            assert!(
                (result.total_weight() - serial.total_weight()).abs()
                    <= 1e-12 * serial.total_weight()
            );
            for (merged, serial) in result.values().value().iter().zip(serial.values().value()) {
                assert!(*serial > 0.1, "fluctuating atoms");
                assert!(
                    (merged - serial).abs() <= bound * serial,
                    "offset {offset}, {} ranges: {merged} != {serial}",
                    ranges.len()
                );
            }
        }
    }
}

#[test]
fn merged_contact_occupancy_depends_on_the_split_alone_and_matches_serial_observation() {
    let top = chain_topology();
    let ids = top.atom_ids();
    let pairs = (0..ATOMS)
        .flat_map(|a| (a + 1..ATOMS).map(move |b| (a, b)))
        .map(|(a, b)| (ids[a], ids[b]))
        .collect::<Vec<_>>();
    let cutoff = Quantity::new(0.45, NANOMETER);
    let (frames, weights) = weighted_frames(&top, 0.0);
    // Weights of one sum exactly; fractional weights sum up to rounding.
    for weights in [vec![1.0; OBSERVATIONS], weights] {
        let empty =
            || ContactOccupancyAccumulator::new(&top, pairs.iter().copied(), cutoff).unwrap();
        let mut serial = empty();
        for (index, weight) in weights.iter().enumerate() {
            serial
                .observe(index, &frames.get(index).unwrap(), *weight)
                .unwrap();
        }
        let serial = serial.finish().unwrap();
        assert!(serial
            .contacts()
            .any(|(_, occupancy)| occupancy > 0.0 && occupancy < 1.0));
        for ranges in splits() {
            let results = executors()
                .iter()
                .map(|executor| {
                    merged(
                        executor,
                        &ranges,
                        empty,
                        |contacts, index| {
                            contacts
                                .observe(index, &frames.get(index).unwrap(), weights[index])
                                .unwrap()
                        },
                        |contacts, partial| contacts.merge(partial).unwrap(),
                    )
                    .finish()
                    .unwrap()
                })
                .collect::<Vec<_>>();
            for result in &results {
                assert_eq!(
                    bits(result.contact_weights()),
                    bits(results[0].contact_weights())
                );
                assert_eq!(
                    result.total_weight().to_bits(),
                    results[0].total_weight().to_bits()
                );
            }
            let result = &results[0];
            assert_eq!(result.observation_count(), OBSERVATIONS);
            assert_eq!(result.pairs(), serial.pairs());
            let exact = weights.iter().all(|weight| *weight == 1.0);
            for (merged, serial) in result
                .contact_weights()
                .iter()
                .chain([&result.total_weight()])
                .zip(
                    serial
                        .contact_weights()
                        .iter()
                        .chain([&serial.total_weight()]),
                )
            {
                if exact {
                    assert_eq!(merged, serial);
                } else {
                    assert!((merged - serial).abs() <= 1e-12 * serial.max(1.0));
                }
            }
        }
    }
}

/// Squared distances decide contacts away from the cutoff; ties and near
/// ties must still follow the distance itself.
#[test]
fn contacts_near_the_cutoff_follow_the_distance() {
    let top = kekule::smiles::to_topology("CC").unwrap();
    let ids = top.atom_ids();
    let mut numbers = Numbers(0x9e37_79b9_7f4a_7c15);
    for cutoff in [0.0, 1e-170, 0.3, 0.45, 2.0, 1e160] {
        let accumulator = ContactOccupancyAccumulator::new(
            &top,
            [(ids[0], ids[1])],
            Quantity::new(cutoff, NANOMETER),
        )
        .unwrap();
        let mut directions = vec![(1.0, 0.0, 0.0), (0.0, 0.6, 0.8), (0.48, 0.6, 0.64)];
        for _ in 0..200 {
            let (x, y, z) = (
                numbers.next() - 0.5,
                numbers.next() - 0.5,
                numbers.next() - 0.5,
            );
            let length = x.hypot(y).hypot(z);
            directions.push((x / length, y / length, z / length));
        }
        for (x, y, z) in directions {
            for ulps in -8_i64..=8 {
                let scale = f64::from_bits((cutoff.to_bits() as i64 + ulps).max(0) as u64);
                let b = Point3::new(x * scale, y * scale, z * scale);
                let delta = b - Point3::new(0.0, 0.0, 0.0);
                let expected = delta.x.hypot(delta.y).hypot(delta.z) <= cutoff;
                let positions = Positions::new(Quantity::new(
                    vec![Point3::new(0.0, 0.0, 0.0), b],
                    NANOMETER,
                ))
                .unwrap();
                let frame =
                    Trajectory::from_items(Arc::clone(&top), [TrajectoryFrame::new(positions)])
                        .unwrap();
                let mut contacts = accumulator.clone();
                contacts.observe(0, &frame.get(0).unwrap(), 1.0).unwrap();
                let contacts = contacts.finish().unwrap();
                assert_eq!(
                    contacts.contact_weights()[0] == 1.0,
                    expected,
                    "cutoff {cutoff}, direction ({x}, {y}, {z}), {ulps} ulps"
                );
            }
        }
    }
}

fn frame(top: &Arc<Topology>, x: [f64; ATOMS]) -> Trajectory {
    let points = x.map(|x| Point3::new(x, 0.0, 0.0));
    Trajectory::from_items(
        Arc::clone(top),
        [TrajectoryFrame::new(
            Positions::new(Quantity::new(points, NANOMETER)).unwrap(),
        )],
    )
    .unwrap()
}

fn rmsf_of(
    top: &Arc<Topology>,
    selection: &AtomSelection,
    rows: &[([f64; ATOMS], f64)],
) -> RmsfAccumulator {
    let mut rmsf = RmsfAccumulator::new(selection).unwrap();
    for (index, (x, weight)) in rows.iter().enumerate() {
        rmsf.observe(index, &frame(top, *x).get(0).unwrap(), *weight)
            .unwrap();
    }
    rmsf
}

fn contacts_of(
    top: &Arc<Topology>,
    pairs: &[(InstanceAtomId, InstanceAtomId)],
    cutoff: f64,
    rows: &[([f64; ATOMS], f64)],
) -> ContactOccupancyAccumulator {
    let mut contacts = ContactOccupancyAccumulator::new(
        top,
        pairs.iter().copied(),
        Quantity::new(cutoff, NANOMETER),
    )
    .unwrap();
    for (index, (x, weight)) in rows.iter().enumerate() {
        contacts
            .observe(index, &frame(top, *x).get(0).unwrap(), *weight)
            .unwrap();
    }
    contacts
}

/// The state an accumulator would finish with, compared bit for bit.
fn rmsf_state(rmsf: &RmsfAccumulator) -> (usize, u64, Vec<u64>) {
    (
        rmsf.observation_count(),
        rmsf.total_weight().to_bits(),
        bits(rmsf.clone().finish().unwrap().values().value()),
    )
}

fn contact_state(contacts: &ContactOccupancyAccumulator) -> (usize, u64, Vec<u64>) {
    (
        contacts.observation_count(),
        contacts.total_weight().to_bits(),
        bits(contacts.clone().finish().unwrap().contact_weights()),
    )
}

#[test]
fn incompatible_or_overflowing_merges_fail_with_typed_errors_and_change_nothing() {
    let top = chain_topology();
    let ids = top.atom_ids();
    let all = AtomSelection::all(Arc::clone(&top));
    let first = AtomSelection::from_atoms(Arc::clone(&top), [ids[0]]).unwrap();
    let other_layout = chain_topology();
    let rows = [
        ([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], 1.0),
        ([0.1; ATOMS], 2.0),
    ];

    let mut rmsf = rmsf_of(&top, &all, &rows);
    let before = rmsf_state(&rmsf);
    for (other, mismatch) in [
        (rmsf_of(&top, &first, &rows), MergeMismatch::Selection),
        (
            rmsf_of(
                &other_layout,
                &AtomSelection::all(Arc::clone(&other_layout)),
                &rows,
            ),
            MergeMismatch::Topology,
        ),
    ] {
        assert_eq!(
            rmsf.merge(&other),
            Err(StatisticsError::MergeMismatch(mismatch))
        );
        assert_eq!(rmsf_state(&rmsf), before);
    }
    // An empty accumulator of another selection still does not merge.
    assert_eq!(
        rmsf.merge(&RmsfAccumulator::new(&first).unwrap()),
        Err(StatisticsError::MergeMismatch(MergeMismatch::Selection))
    );
    let heavy = rmsf_of(&top, &all, &[(rows[0].0, f64::MAX)]);
    let mut heavy_rmsf = heavy.clone();
    assert_eq!(
        heavy_rmsf.merge(&heavy),
        Err(StatisticsError::MergeWeightOverflow)
    );
    assert_eq!(rmsf_state(&heavy_rmsf), rmsf_state(&heavy));
    // Each partial holds one far position; their deviation squared overflows.
    let mut far = rmsf_of(
        &top,
        &all,
        &[([1e300, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 1.0)],
    );
    let far_before = rmsf_state(&far);
    let opposite = rmsf_of(
        &top,
        &all,
        &[([-1e300, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 1.0)],
    );
    assert_eq!(
        far.merge(&opposite),
        Err(StatisticsError::MergeNumericalFailure { atom: ids[0] })
    );
    assert_eq!(rmsf_state(&far), far_before);
    // A failed merge leaves the accumulator usable.
    rmsf.merge(&rmsf_of(&top, &all, &rows)).unwrap();
    assert_eq!(rmsf.observation_count(), 4);
    assert_eq!(rmsf.total_weight(), 6.0);

    let pairs = [(ids[0], ids[1]), (ids[2], ids[5])];
    let mut contacts = contacts_of(&top, &pairs, 0.15, &rows);
    let before = contact_state(&contacts);
    let other_pairs = [(ids[0], ids[1]), (ids[2], ids[6])];
    let reordered = [pairs[1], pairs[0]];
    let reversed = [(ids[1], ids[0]), pairs[1]];
    let other_ids = other_layout.atom_ids();
    for (other, mismatch) in [
        (
            contacts_of(&top, &other_pairs, 0.15, &rows),
            MergeMismatch::Pairs,
        ),
        (
            contacts_of(&top, &reordered, 0.15, &rows),
            MergeMismatch::Pairs,
        ),
        (
            contacts_of(&top, &reversed, 0.15, &rows),
            MergeMismatch::Pairs,
        ),
        (
            contacts_of(&top, &pairs[..1], 0.15, &rows),
            MergeMismatch::Pairs,
        ),
        (
            contacts_of(&top, &pairs, 0.25, &rows),
            MergeMismatch::Cutoff,
        ),
        (
            contacts_of(
                &other_layout,
                &[(other_ids[0], other_ids[1]), (other_ids[2], other_ids[5])],
                0.15,
                &rows,
            ),
            MergeMismatch::Topology,
        ),
    ] {
        assert_eq!(
            contacts.merge(&other),
            Err(StatisticsError::MergeMismatch(mismatch))
        );
        assert_eq!(contact_state(&contacts), before);
    }
    // The same cutoff in another unit is the same cutoff.
    let angstrom =
        ContactOccupancyAccumulator::new(&top, pairs, Quantity::new(1.5, ANGSTROM)).unwrap();
    contacts.merge(&angstrom).unwrap();
    assert_eq!(contact_state(&contacts), before);
    let heavy = contacts_of(&top, &pairs, 0.15, &[(rows[0].0, f64::MAX)]);
    let mut heavy_contacts = heavy.clone();
    assert_eq!(
        heavy_contacts.merge(&heavy),
        Err(StatisticsError::MergeWeightOverflow)
    );
    assert_eq!(contact_state(&heavy_contacts), contact_state(&heavy));
    contacts
        .merge(&contacts_of(&top, &pairs, 0.15, &rows))
        .unwrap();
    let contacts = contacts.finish().unwrap();
    assert_eq!(contacts.observation_count(), 4);
    assert_eq!(contacts.contact_weights(), &[6.0, 4.0]);
}

#[test]
fn merging_empty_accumulators_changes_nothing() {
    let top = chain_topology();
    let all = AtomSelection::all(Arc::clone(&top));
    let rows = [
        ([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7], 1.0),
        ([0.3; ATOMS], 3.0),
    ];
    let observed = rmsf_of(&top, &all, &rows);
    let mut rmsf = observed.clone();
    rmsf.merge(&RmsfAccumulator::new(&all).unwrap()).unwrap();
    assert_eq!(rmsf_state(&rmsf), rmsf_state(&observed));
    let mut empty = RmsfAccumulator::new(&all).unwrap();
    empty.merge(&RmsfAccumulator::new(&all).unwrap()).unwrap();
    assert_eq!(empty.observation_count(), 0);
    assert_eq!(
        empty.clone().finish().unwrap_err(),
        StatisticsError::NoObservations
    );
    empty.merge(&observed).unwrap();
    assert_eq!(rmsf_state(&empty), rmsf_state(&observed));

    let ids = top.atom_ids();
    let pairs = [(ids[0], ids[1])];
    let observed = contacts_of(&top, &pairs, 0.15, &rows);
    let mut empty = contacts_of(&top, &pairs, 0.15, &[]);
    empty.merge(&contacts_of(&top, &pairs, 0.15, &[])).unwrap();
    assert_eq!(empty.observation_count(), 0);
    empty.merge(&observed).unwrap();
    assert_eq!(contact_state(&empty), contact_state(&observed));
}

/// Observations far from overflowing update the moments in place; those near
/// it go through pending copies that can fail. Scaling every coordinate by a
/// power of two scales the RMSF exactly, so both routes must agree bit for
/// bit, for observation and for merging.
#[test]
fn moments_near_overflow_agree_with_moments_far_from_it() {
    let top = chain_topology();
    let all = AtomSelection::all(Arc::clone(&top));
    let (frames, weights) = weighted_frames(&top, 0.0);
    let scaled = |scale: f64| {
        Trajectory::from_items(
            Arc::clone(&top),
            frames.iter().map(|frame| {
                let points = frame
                    .positions()
                    .values()
                    .value()
                    .iter()
                    .map(|p| Point3::new(p.x * scale, p.y * scale, p.z * scale))
                    .collect::<Vec<_>>();
                TrajectoryFrame::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap())
            }),
        )
        .unwrap()
    };
    let reduce = |frames: &Trajectory, ranges: &[Range<usize>]| {
        let mut rmsf = RmsfAccumulator::new(&all).unwrap();
        for range in ranges {
            let mut partial = RmsfAccumulator::new(&all).unwrap();
            for index in range.clone() {
                partial
                    .observe(index, &frames.get(index).unwrap(), weights[index])
                    .unwrap();
            }
            rmsf.merge(&partial).unwrap();
        }
        rmsf.finish().unwrap()
    };
    // 2⁵⁰⁵ puts the squared coordinates times the total weight within a
    // factor of 48 of overflowing.
    let scale = 2.0_f64.powi(505);
    let near = scaled(scale);
    for ranges in splits() {
        let far = reduce(&frames, &ranges);
        let near = reduce(&near, &ranges);
        let unscaled = near
            .values()
            .value()
            .iter()
            .map(|value| value / scale)
            .collect::<Vec<_>>();
        assert_eq!(bits(&unscaled), bits(far.values().value()));
    }
}
