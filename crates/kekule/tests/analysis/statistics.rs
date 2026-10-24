use std::sync::Arc;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::statistics::{ContactOccupancyAccumulator, RmsfAccumulator, StatisticsError};
use kekule::structure::{
    Ensemble, EnsembleSample, EnsembleState, Positions, Trajectory, TrajectoryFrame,
};
use kekule::topology::{AtomSelection, InstanceAtomId, SelectionError, Topology};
use kekule::units::{Quantity, ANGSTROM, NANOMETER, PICOSECOND};

fn linear_carbon_topology() -> Arc<Topology> {
    kekule::smiles::to_topology("CCC").unwrap()
}

fn positions(points: [Point3; 3]) -> Positions {
    Positions::new(Quantity::new(points, NANOMETER)).unwrap()
}

/// One frame per row, with the three atoms on the x axis.
fn trajectory(top: &Arc<Topology>, rows: &[[f64; 3]]) -> Trajectory {
    Trajectory::from_items(
        top.clone(),
        rows.iter()
            .map(|row| TrajectoryFrame::new(positions(row.map(|x| Point3::new(x, 0.0, 0.0))))),
    )
    .unwrap()
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-8, "{a} != {b}");
}

/// Three structures with weights 1, 2, 1: atom 0 at x = 0, 1, 3; atom 1 at
/// x = 0.1 and y = 0, 0, 2; atom 2 fixed far away.
const WEIGHTS: [f64; 3] = [1.0, 2.0, 1.0];

fn weighted_structures() -> [Positions; 3] {
    [(0.0, 0.0), (1.0, 0.0), (3.0, 2.0)].map(|(x, y)| {
        positions([
            Point3::new(x, 0.0, 0.0),
            Point3::new(0.1, y, 0.0),
            Point3::new(5.0, 5.0, 5.0),
        ])
    })
}

/// Hand computation for [`weighted_structures`] with total weight 4:
/// atom 0 has mean 1.25 and variance (1.5625 + 2 × 0.0625 + 3.0625) / 4;
/// atom 1 has mean y 0.5 and variance (0.25 + 2 × 0.25 + 2.25) / 4.
fn expected_weighted_rmsf() -> [f64; 3] {
    [(4.75_f64 / 4.0).sqrt(), (3.0_f64 / 4.0).sqrt(), 0.0]
}

#[test]
fn weighted_accumulators_match_hand_computation() {
    let top = linear_carbon_topology();
    let ids = top.atom_ids();
    let structures =
        Trajectory::from_items(top.clone(), weighted_structures().map(TrajectoryFrame::new))
            .unwrap();
    let pairs = [(ids[0], ids[1]), (ids[1], ids[2])];
    let mut rmsf = RmsfAccumulator::new(&AtomSelection::all(top.clone())).unwrap();
    let mut contacts =
        ContactOccupancyAccumulator::new(&top, pairs, Quantity::new(1.0, NANOMETER)).unwrap();
    for (index, (frame, weight)) in structures.iter().zip(WEIGHTS).enumerate() {
        rmsf.observe(index, &frame, weight).unwrap();
        contacts.observe(index, &frame, weight).unwrap();
    }
    assert_eq!(rmsf.observation_count(), 3);
    assert_eq!(rmsf.total_weight(), 4.0);
    let rmsf = rmsf.finish().unwrap();
    for (value, expected) in rmsf.values().value().iter().zip(expected_weighted_rmsf()) {
        close(*value, expected);
    }
    assert_eq!(rmsf.total_weight(), 4.0);
    // Distances 0.1, 0.9, and √(2.9² + 2²) nm: in contact in the first two,
    // which carry weight 3 of 4.
    let contacts = contacts.finish().unwrap();
    assert_eq!(contacts.contact_weights(), &[3.0, 0.0]);
    assert_eq!(contacts.total_weight(), 4.0);
    assert_eq!(contacts.observation_count(), 3);
    assert_eq!(
        contacts.contacts().collect::<Vec<_>>(),
        [(pairs[0], 0.75), (pairs[1], 0.0)]
    );
}

#[test]
fn ensemble_statistics_use_sample_weights_not_states() {
    let top = linear_carbon_topology();
    let ids = top.atom_ids();
    let [first, second, third] = weighted_structures();
    // The two heavier samples share a state; populations do not enter twice.
    let mut ensemble = Ensemble::new(top.clone());
    ensemble
        .push_state(
            EnsembleSample::new(first, 10.0).unwrap(),
            EnsembleState::new(),
        )
        .unwrap();
    let state = ensemble
        .push_state(
            EnsembleSample::new(second, 20.0).unwrap(),
            EnsembleState::new(),
        )
        .unwrap();
    ensemble
        .push_sample(EnsembleSample::new(third, 10.0).unwrap(), state)
        .unwrap();

    let rmsf = ensemble.rmsf(&AtomSelection::all(top.clone())).unwrap();
    for (value, expected) in rmsf.values().value().iter().zip(expected_weighted_rmsf()) {
        close(*value, expected);
    }
    assert_eq!(rmsf.observation_count(), 3);
    // Samples are observed with their probabilities.
    close(rmsf.total_weight(), 1.0);

    let pairs = [(ids[0], ids[1]), (ids[1], ids[2])];
    let contacts = ensemble
        .contact_occupancy(pairs, Quantity::new(10.0, ANGSTROM))
        .unwrap();
    assert_eq!(contacts.contact_weights(), &[0.75, 0.0]);
    assert_eq!(
        contacts.contacts().collect::<Vec<_>>(),
        [(pairs[0], 0.75), (pairs[1], 0.0)]
    );

    // Weights near the floating-point limits are scaled before summing.
    let mut extreme = Ensemble::new(top.clone());
    for (sample, weight) in weighted_structures().into_iter().zip(WEIGHTS) {
        let sample = EnsembleSample::new(sample, weight * 1.0e307).unwrap();
        extreme.push_state(sample, EnsembleState::new()).unwrap();
    }
    let extreme = extreme.rmsf(&AtomSelection::all(top.clone())).unwrap();
    for (value, expected) in extreme
        .values()
        .value()
        .iter()
        .zip(expected_weighted_rmsf())
    {
        close(*value, expected);
    }

    assert!(matches!(
        Ensemble::new(top.clone()).rmsf(&AtomSelection::all(top.clone())),
        Err(StatisticsError::NoObservations)
    ));
    assert!(matches!(
        Ensemble::new(top.clone()).contact_occupancy(pairs, Quantity::new(1.0, NANOMETER)),
        Err(StatisticsError::NoObservations)
    ));
    assert!(matches!(
        ensemble.rmsf(&AtomSelection::all(linear_carbon_topology())),
        Err(StatisticsError::Selection(SelectionError::TopologyMismatch))
    ));
}

#[test]
fn equal_weights_reproduce_frame_counting() {
    let top = linear_carbon_topology();
    let ids = top.atom_ids();
    let pairs = [(ids[0], ids[1]), (ids[1], ids[2])];
    let all = AtomSelection::all(top.clone());
    let cutoff = Quantity::new(2.0, NANOMETER);
    let rows = [[0.0, 1.0, 20.0], [2.0, 5.0, 20.0], [4.0, 6.0, 20.0]];
    let t = trajectory(&top, &rows);
    let rmsf = t.rmsf(&all).unwrap();
    let contacts = t.contact_occupancy(pairs, cutoff).unwrap();
    // Frame counts and population variances (division by N).
    assert_eq!(rmsf.total_weight(), 3.0);
    assert_eq!(contacts.contact_weights(), &[2.0, 0.0]);
    close(rmsf.values().value()[0], (8.0_f64 / 3.0).sqrt());
    close(rmsf.values().value()[1], (14.0_f64 / 3.0).sqrt());

    // Any common weight gives the same statistics.
    let mut scaled_rmsf = RmsfAccumulator::new(&all).unwrap();
    let mut scaled_contacts = ContactOccupancyAccumulator::new(&top, pairs, cutoff).unwrap();
    for (index, frame) in t.iter().enumerate() {
        scaled_rmsf.observe(index, &frame, 0.3).unwrap();
        scaled_contacts.observe(index, &frame, 0.3).unwrap();
    }
    for (scaled, unit) in scaled_rmsf
        .finish()
        .unwrap()
        .values()
        .value()
        .iter()
        .zip(rmsf.values().value())
    {
        close(*scaled, *unit);
    }
    let scaled_contacts = scaled_contacts.finish().unwrap();
    for ((_, scaled), (_, unit)) in scaled_contacts.contacts().zip(contacts.contacts()) {
        close(scaled, unit);
    }

    // The projection to an equally weighted ensemble keeps every statistic.
    let ensemble = t.clone().into_ensemble();
    for (sample, frame) in ensemble
        .rmsf(&all)
        .unwrap()
        .values()
        .value()
        .iter()
        .zip(rmsf.values().value())
    {
        close(*sample, *frame);
    }
    for ((_, sample), (_, frame)) in ensemble
        .contact_occupancy(pairs, cutoff)
        .unwrap()
        .contacts()
        .zip(contacts.contacts())
    {
        close(sample, frame);
    }
}

#[test]
fn trajectory_statistics_agree_with_independent_computation() {
    let top = linear_carbon_topology();
    let ids = top.atom_ids();
    let pairs = [(ids[0], ids[1]), (ids[1], ids[2])];
    let all = AtomSelection::all(Arc::clone(&top));
    // Large common translation exposes cancellation in E[x²] - E[x]².
    let rows = [
        [1.0e9, 1.0e9 + 1.0, 1.0e9 + 20.0],
        [1.0e9 + 2.0, 1.0e9 + 5.0, 1.0e9 + 20.0],
        [1.0e9 + 4.0, 1.0e9 + 6.0, 1.0e9 + 20.0],
    ];
    let t = trajectory(&top, &rows);
    let before = format!("{t:?}");
    let rmsf = t.rmsf(&all).unwrap();
    let contacts = t
        .contact_occupancy(pairs, Quantity::new(2.0, NANOMETER))
        .unwrap();
    assert_eq!(rmsf.observation_count(), 3);
    assert_eq!(rmsf.selection(), &all);
    assert_eq!(contacts.pairs(), &pairs);
    assert!(std::ptr::eq(contacts.topology(), top.as_ref()));
    assert_eq!(contacts.contact_weights(), &[2.0, 0.0]);
    close(contacts.contacts().next().unwrap().1, 2.0 / 3.0);
    close(contacts.cutoff().value_in(ANGSTROM).unwrap(), 20.0);
    for (index, (atom, value)) in rmsf.atoms().enumerate() {
        assert_eq!(atom, ids[index]);
        let mut pairwise = 0.0;
        for a in &rows {
            for b in &rows {
                pairwise += (a[index] - b[index]).powi(2);
            }
        }
        close(value.into_value().powi(2), pairwise / (2.0 * 9.0));
    }
    close(rmsf.values().value()[0], (8.0_f64 / 3.0).sqrt());
    close(rmsf.values().value()[1], (14.0_f64 / 3.0).sqrt());
    close(rmsf.values().value()[2], 0.0);
    let subset = t
        .rmsf(&AtomSelection::from_atoms(Arc::clone(&top), [ids[2], ids[0]]).unwrap())
        .unwrap();
    assert_eq!(
        subset.atoms().map(|(atom, _)| atom).collect::<Vec<_>>(),
        [ids[0], ids[2]]
    );
    assert_eq!(subset.values().value(), &[rmsf.values().value()[0], 0.0]);
    close(
        rmsf.values().value_in(ANGSTROM).unwrap()[0],
        (8.0_f64 / 3.0).sqrt() * 10.0,
    );
    assert_eq!(
        t.contact_occupancy(pairs, Quantity::new(21.0, ANGSTROM))
            .unwrap()
            .contact_weights(),
        &[2.0, 0.0]
    );
    assert_eq!(format!("{t:?}"), before);
}

#[test]
fn failed_observations_are_atomic_and_retryable_for_both_accumulators() {
    let top = linear_carbon_topology();
    let ids = top.atom_ids();
    let pairs = [(ids[0], ids[1]), (ids[1], ids[2])];
    let good = trajectory(&top, &[[0.0, 1.0, 0.0]]);
    // RMSF can fail after updating the first pending atom; contact calculation
    // can fail after the first pending pair. Neither may publish that prefix.
    let bad = trajectory(&top, &[[0.5, 1.0e308, -1.0e308]]);
    let unrelated = trajectory(&linear_carbon_topology(), &[[0.0, 1.0, 0.0]]);
    let mut rmsf = RmsfAccumulator::new(&AtomSelection::all(Arc::clone(&top))).unwrap();
    let mut contacts =
        ContactOccupancyAccumulator::new(&top, pairs, Quantity::new(2.0, NANOMETER)).unwrap();
    rmsf.observe(10, &good.get(0).unwrap(), 1.0).unwrap();
    contacts.observe(10, &good.get(0).unwrap(), 1.0).unwrap();
    assert_eq!(
        rmsf.observe(11, &unrelated.get(0).unwrap(), 1.0),
        Err(StatisticsError::TopologyMismatch { index: 11 })
    );
    assert_eq!(
        contacts.observe(11, &unrelated.get(0).unwrap(), 1.0),
        Err(StatisticsError::TopologyMismatch { index: 11 })
    );
    assert!(matches!(
        rmsf.observe(11, &bad.get(0).unwrap(), 1.0),
        Err(StatisticsError::NumericalFailure { index: 11, .. })
    ));
    assert!(matches!(
        contacts.observe(11, &bad.get(0).unwrap(), 1.0),
        Err(StatisticsError::NumericalFailure { index: 11, .. })
    ));
    for weight in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            rmsf.observe(11, &good.get(0).unwrap(), weight),
            Err(StatisticsError::InvalidWeight { index: 11 })
        );
        assert_eq!(
            contacts.observe(11, &good.get(0).unwrap(), weight),
            Err(StatisticsError::InvalidWeight { index: 11 })
        );
    }
    let mut heavy_rmsf = RmsfAccumulator::new(&AtomSelection::all(Arc::clone(&top))).unwrap();
    let mut heavy_contacts =
        ContactOccupancyAccumulator::new(&top, pairs, Quantity::new(2.0, NANOMETER)).unwrap();
    heavy_rmsf
        .observe(0, &good.get(0).unwrap(), f64::MAX)
        .unwrap();
    heavy_contacts
        .observe(0, &good.get(0).unwrap(), f64::MAX)
        .unwrap();
    assert_eq!(
        heavy_rmsf.observe(1, &good.get(0).unwrap(), f64::MAX),
        Err(StatisticsError::WeightOverflow { index: 1 })
    );
    assert_eq!(
        heavy_contacts.observe(1, &good.get(0).unwrap(), f64::MAX),
        Err(StatisticsError::WeightOverflow { index: 1 })
    );
    assert_eq!(heavy_rmsf.observation_count(), 1);
    assert_eq!(heavy_contacts.total_weight(), f64::MAX);
    assert_eq!(rmsf.observation_count(), 1);
    assert_eq!(rmsf.total_weight(), 1.0);
    assert_eq!(contacts.observation_count(), 1);
    assert_eq!(contacts.total_weight(), 1.0);
    rmsf.observe(11, &good.get(0).unwrap(), 1.0).unwrap();
    contacts.observe(11, &good.get(0).unwrap(), 1.0).unwrap();
    let rmsf = rmsf.finish().unwrap();
    let contacts = contacts.finish().unwrap();
    assert_eq!(rmsf.observation_count(), 2);
    assert_eq!(rmsf.values().value(), &[0.0, 0.0, 0.0]);
    assert_eq!(contacts.observation_count(), 2);
    assert_eq!(contacts.contact_weights(), &[2.0, 2.0]);
}

#[test]
fn statistics_reject_empty_or_ambiguous_inputs_and_handle_single_observations() {
    let top = linear_carbon_topology();
    let ids = top.atom_ids();
    let all = AtomSelection::all(Arc::clone(&top));
    let empty = AtomSelection::from_atoms(Arc::clone(&top), []).unwrap();
    let cutoff = Quantity::new(1.0, NANOMETER);
    assert!(matches!(
        RmsfAccumulator::new(&empty),
        Err(StatisticsError::EmptySelection)
    ));
    assert!(matches!(
        RmsfAccumulator::new(&all).unwrap().finish(),
        Err(StatisticsError::NoObservations)
    ));
    assert!(matches!(
        ContactOccupancyAccumulator::new(&top, [], cutoff),
        Err(StatisticsError::EmptySelection)
    ));
    assert!(matches!(
        ContactOccupancyAccumulator::new(&top, [(ids[0], ids[1])], cutoff)
            .unwrap()
            .finish(),
        Err(StatisticsError::NoObservations)
    ));
    assert!(matches!(
        ContactOccupancyAccumulator::new(&top, [(ids[0], ids[0])], cutoff),
        Err(StatisticsError::SelfPair { pair: 0 })
    ));
    assert!(matches!(
        ContactOccupancyAccumulator::new(&top, [(ids[0], ids[1]), (ids[1], ids[0])], cutoff),
        Err(StatisticsError::DuplicatePair { pair: 1 })
    ));
    let invalid = InstanceAtomId::new(ids[0].instance(), kekule::molecule::AtomId::new(99));
    assert!(matches!(
        ContactOccupancyAccumulator::new(&top, [(ids[0], invalid)], cutoff),
        Err(StatisticsError::Selection(SelectionError::InvalidAtomId(_)))
    ));
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            ContactOccupancyAccumulator::new(
                &top,
                [(ids[0], ids[1])],
                Quantity::new(value, NANOMETER)
            ),
            Err(StatisticsError::InvalidCutoff)
        ));
    }
    assert!(matches!(
        ContactOccupancyAccumulator::new(&top, [(ids[0], ids[1])], Quantity::new(1.0, PICOSECOND)),
        Err(StatisticsError::Unit(_))
    ));
    let none = Trajectory::new(top.clone());
    assert!(matches!(
        none.rmsf(&all),
        Err(StatisticsError::NoObservations)
    ));
    assert!(matches!(
        none.contact_occupancy([(ids[0], ids[1])], cutoff),
        Err(StatisticsError::NoObservations)
    ));
    let one = trajectory(&top, &[[0.0, 0.0, 1.0]]);
    assert_eq!(one.rmsf(&all).unwrap().values().value(), &[0.0, 0.0, 0.0]);
    assert_eq!(
        one.contact_occupancy([(ids[0], ids[1])], Quantity::new(0.0, NANOMETER))
            .unwrap()
            .contact_weights(),
        &[1.0]
    );
    let foreign = AtomSelection::all(linear_carbon_topology());
    assert!(matches!(
        one.rmsf(&foreign),
        Err(StatisticsError::Selection(SelectionError::TopologyMismatch))
    ));
}

#[test]
fn stored_coordinate_statistics_do_not_implicitly_image_or_align_periodic_frames() {
    let top = linear_carbon_topology();
    let ids = top.atom_ids();
    let mut t = trajectory(&top, &[[0.0, 9.0, 2.0], [1.0, 10.0, 3.0]]);
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(10.0, 10.0, 10.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    for index in 0..t.len() {
        t.get_mut(index).unwrap().set_cell(Some(cell));
    }
    // A minimum-image contact would match; stored distance is nine nanometers.
    assert_eq!(
        t.contact_occupancy([(ids[0], ids[1])], Quantity::new(2.0, NANOMETER))
            .unwrap()
            .contact_weights(),
        &[0.0]
    );
    for value in t
        .rmsf(&AtomSelection::all(Arc::clone(&top)))
        .unwrap()
        .values()
        .value()
    {
        close(*value, 0.5);
    }
}
