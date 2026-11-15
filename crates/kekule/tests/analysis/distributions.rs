//! Distribution operations: probabilities, rankings, energies, and the free
//! energies of regrouped states.

use kekule::diagnostics::{Countable, LengthMismatch};
use kekule::math::Point3;
use kekule::molecule::{Atom, Element, MoleculeEditor};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{
    Distribution, Ensemble, Geometry, Partition, Positions, Samples, StateData, StructureSet,
    StructureSetError,
};
use kekule::topology::Topology;
use kekule::units::{Quantity, KELVIN, KILOJOULE_PER_MOLE, NANOMETER};

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

fn single_atom_topology() -> Topology {
    let mut editor = MoleculeEditor::new();
    editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    Topology::from_molecule(editor.finish().unwrap()).unwrap()
}

fn single_position(x: f64) -> Positions {
    Positions::new(Quantity::new([Point3::new(x, 0.0, 0.0)], NANOMETER)).unwrap()
}

/// Whether `actual` and `expected` agree to within a few units in the last
/// place of the larger, as values computed through `exp` and `ln` do.
fn close(actual: f64, expected: f64) -> bool {
    (actual - expected).abs() <= 4.0 * f64::EPSILON * actual.abs().max(expected.abs())
}

fn assert_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len(), "{actual:?} != {expected:?}");
    for (a, e) in actual.iter().zip(expected) {
        assert!(close(*a, *e), "{actual:?} != {expected:?}");
    }
}

/// Equally weighted conformers of a one-atom topology, one state each, with
/// state energies in kJ/mol.
fn energy_conformers(energies: &[f64]) -> Ensemble {
    let states = energies
        .iter()
        .map(|energy| {
            StateData::new()
                .with_energy(Quantity::new(*energy, KILOJOULE_PER_MOLE))
                .unwrap()
        })
        .collect();
    let set = StructureSet::from_geometries(
        single_atom_topology(),
        energies.iter().map(|_| Geometry::new(single_position(0.0))),
    )
    .unwrap();
    let distribution = Distribution::new(Samples::All, Partition::state_per_sample(energies.len()))
        .unwrap()
        .with_log_weights(vec![0.0; energies.len()])
        .unwrap()
        .with_states(states)
        .unwrap();
    Ensemble::new(set, distribution).unwrap()
}

/// The natural logarithm of each weight.
fn ln(weights: &[f64]) -> Vec<f64> {
    weights.iter().map(|weight| weight.ln()).collect()
}

/// The exponential of each log weight.
fn exponentials(log_weights: &[f64]) -> Vec<f64> {
    log_weights
        .iter()
        .map(|log_weight| log_weight.exp())
        .collect()
}

/// Geometries of a one-atom topology with these log weights, grouped by
/// `assignments` with `representatives`.
fn grouped(log_weights: &[f64], assignments: &[usize], representatives: &[usize]) -> Ensemble {
    let set = StructureSet::from_geometries(
        single_atom_topology(),
        log_weights
            .iter()
            .map(|_| Geometry::new(single_position(0.0))),
    )
    .unwrap();
    let partition = Partition::new(assignments.to_vec(), representatives.to_vec()).unwrap();
    let distribution = Distribution::new(Samples::All, partition)
        .unwrap()
        .with_log_weights(log_weights.to_vec())
        .unwrap();
    Ensemble::new(set, distribution).unwrap()
}

fn kj(value: f64) -> Option<Quantity<f64>> {
    Some(Quantity::new(value, KILOJOULE_PER_MOLE))
}

/// The energy of state `state` in kJ/mol.
fn kilojoules(ensemble: &Ensemble, state: usize) -> f64 {
    ensemble
        .state(state)
        .unwrap()
        .energy()
        .unwrap()
        .value_in(KILOJOULE_PER_MOLE)
        .unwrap()
}

/// Sets the energy of every state, in kJ/mol where given.
fn set_energies(set: &mut Ensemble, energies: &[Option<f64>]) {
    let mut distribution = set.distribution_mut();
    for (state, energy) in energies.iter().enumerate() {
        distribution
            .state_mut(state)
            .unwrap()
            .set_energy(energy.and_then(kj))
            .unwrap();
    }
}

/// The log weight of every geometry, each a sample.
fn log_weights(ensemble: &Ensemble) -> Vec<f64> {
    (0..ensemble.set().len())
        .map(|geometry| ensemble.distribution().log_weight(geometry).unwrap())
        .collect()
}

fn energies(ensemble: &Ensemble) -> Vec<Option<Quantity<f64>>> {
    ensemble.states().map(|state| state.energy()).collect()
}

fn temperature(ensemble: &Ensemble) -> Option<Quantity<f64>> {
    ensemble.distribution().temperature()
}

#[test]
fn regrouping_keeps_weights_and_gives_merged_states_their_free_energy() {
    // Four conformers, energies 0..3 kJ/mol, weights 1..4.
    let weights = ln(&[1.0, 2.0, 3.0, 4.0]);
    let mut set = energy_conformers(&[0.0, 1.0, 2.0, 3.0]);
    let mut distribution = set.distribution_mut();
    distribution.set_log_weights(Some(weights.clone())).unwrap();
    distribution
        .state_mut(1)
        .unwrap()
        .properties_mut()
        .insert(key("label"), PropertyValue::Int(1))
        .unwrap();
    distribution
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    let before = set.clone();

    // An unchanged grouping keeps everything.
    let unchanged = set.distribution().partition().clone();
    set.distribution_mut().regroup(unchanged).unwrap();
    assert_eq!(set.distribution(), before.distribution());
    assert_eq!(energies(&set), energies(&before));
    assert_eq!(
        set.state(1).unwrap().properties(),
        before.state(1).unwrap().properties()
    );

    // Merging conformers 0 and 1 under 1 gives their free energy at 300 K
    // and keeps conformer 1's properties; whole states keep their energy.
    let merged = Partition::new(vec![0, 0, 1, 2], vec![1, 2, 3]).unwrap();
    set.distribution_mut().regroup(merged).unwrap();
    assert_eq!(log_weights(&set), weights);
    let rt: f64 = 8.314_462_618_153_24e-3 * 300.0;
    let free = -rt * (1.0 + (-1.0 / rt).exp()).ln();
    assert!(close(kilojoules(&set, 0), free));
    assert_eq!(energies(&set)[1..], [kj(2.0), kj(3.0)]);
    assert_eq!(
        set.state(0).unwrap().properties().get(&key("label")),
        Some(&PropertyValue::Int(1))
    );
    // The merged state weighs 1 + 2.
    assert!(close(set.state(0).unwrap().log_population(), 3.0_f64.ln()));
    assert_eq!(temperature(&set), Some(Quantity::new(300.0, KELVIN)));
    // Geometries never move.
    assert_eq!(
        set.set()
            .iter()
            .map(|geometry| geometry.positions().clone())
            .collect::<Vec<_>>(),
        before
            .set()
            .iter()
            .map(|geometry| geometry.positions().clone())
            .collect::<Vec<_>>()
    );

    // Splitting the merged state loses its energy and properties; the state
    // that takes a whole old state with part of a split one has no energy.
    let split = Partition::new(vec![0, 1, 0, 2], vec![2, 1, 3]).unwrap();
    set.distribution_mut().regroup(split).unwrap();
    assert_eq!(log_weights(&set), weights);
    assert_eq!(energies(&set), [None, None, kj(3.0)]);
    assert!(set.state(1).unwrap().properties().is_empty());

    // A partition of another sample count changes nothing.
    let snapshot = set.clone();
    assert_eq!(
        set.distribution_mut().regroup(Partition::one_state(3)),
        Err(StructureSetError::PartitionMismatch(LengthMismatch::new(
            Countable::Sample,
            4,
            3
        )))
    );
    assert_eq!(set.distribution(), snapshot.distribution());

    // One state split into several loses its energy.
    let mut one = grouped(&[0.0; 4], &[0, 0, 0, 0], &[0]);
    set_energies(&mut one, &[Some(5.0)]);
    one.distribution_mut()
        .regroup(Partition::new(vec![0, 0, 1, 1], vec![0, 2]).unwrap())
        .unwrap();
    assert_eq!(energies(&one), [None, None]);
    assert_eq!(log_weights(&one), [0.0; 4]);
}

#[test]
fn merged_free_energies_match_boltzmann_populations() {
    // Five conformers; states {0, 1, 2} and {3, 4} after merging.
    let energies_kj = [4.0, 0.5, 9.0, 2.0, 2.0];
    let temperature = Quantity::new(310.0, KELVIN);
    let rt: f64 = 8.314_462_618_153_24e-3 * 310.0;
    let mut set = energy_conformers(&energies_kj);
    set.distribution_mut()
        .assign_boltzmann_populations(temperature)
        .unwrap();
    let source = set.clone();
    let partition = Partition::new(vec![0, 0, 0, 1, 1], vec![0, 3]).unwrap();
    set.distribution_mut().regroup(partition.clone()).unwrap();

    // The medoid of state 0 lies at 4 kJ/mol, but its free energy is set by
    // all members, and the states' free energy difference is the log ratio
    // of their populations.
    let free = |members: &[f64]| -rt * members.iter().map(|e| (-e / rt).exp()).sum::<f64>().ln();
    assert!((kilojoules(&set, 0) - free(&energies_kj[..3])).abs() < 1e-12);
    assert!((kilojoules(&set, 1) - (2.0 - rt * 2.0_f64.ln())).abs() < 1e-12);
    let log_ratio = set.state(0).unwrap().log_population() - set.state(1).unwrap().log_population();
    assert!((kilojoules(&set, 0) - kilojoules(&set, 1) + rt * log_ratio).abs() < 1e-12);

    // Reassigning populations from the free energies at the same
    // temperature changes no state population.
    let mut reweighted = set.clone();
    reweighted
        .distribution_mut()
        .assign_boltzmann_populations(temperature)
        .unwrap();
    for state in 0..2 {
        assert!(
            (reweighted.state(state).unwrap().log_population()
                - set.state(state).unwrap().log_population())
            .abs()
                < 1e-12
        );
    }

    // Merging in steps gives the free energy of merging at once.
    let mut stepwise = source.clone();
    stepwise
        .distribution_mut()
        .regroup(Partition::new(vec![0, 0, 1, 2, 2], vec![0, 2, 3]).unwrap())
        .unwrap();
    stepwise
        .distribution_mut()
        .regroup(partition.clone())
        .unwrap();
    for state in 0..2 {
        assert!((kilojoules(&stepwise, state) - kilojoules(&set, state)).abs() < 1e-12);
    }

    // A changed temperature keeps the free energy of the merge temperature.
    set.distribution_mut()
        .set_temperature(Some(Quantity::new(400.0, KELVIN)))
        .unwrap();
    assert!((kilojoules(&set, 0) - free(&energies_kj[..3])).abs() < 1e-12);

    // Without a temperature, or with a member lacking an energy, a merged
    // state has none; a state of one old state keeps its energy.
    let mut cold = source.clone();
    cold.distribution_mut().set_temperature(None).unwrap();
    cold.distribution_mut()
        .regroup(Partition::new(vec![0, 0, 0, 1, 2], vec![0, 3, 4]).unwrap())
        .unwrap();
    assert_eq!(energies(&cold), [None, kj(2.0), kj(2.0)]);
    let mut missing = source.clone();
    set_energies(
        &mut missing,
        &[Some(4.0), None, Some(9.0), Some(2.0), Some(2.0)],
    );
    missing.distribution_mut().regroup(partition).unwrap();
    assert_eq!(energies(&missing)[0], None);
    assert!((kilojoules(&missing, 1) - (2.0 - rt * 2.0_f64.ln())).abs() < 1e-12);

    // Energies far apart merge without overflow.
    let mut wide = energy_conformers(&[0.0, 1.0e6]);
    wide.distribution_mut()
        .set_temperature(Some(temperature))
        .unwrap();
    wide.distribution_mut()
        .regroup(Partition::one_state(2))
        .unwrap();
    assert_eq!(kilojoules(&wide, 0), 0.0);
}

#[test]
fn distribution_operations_keep_their_results() {
    // State 0: geometries 0 and 2 (weights 1 and 3), represented by 2;
    // state 1: geometries 1 and 3 (weights 2 and 4), represented by 1.
    let weights = [1.0_f64, 2.0, 3.0, 4.0].map(f64::ln);
    let mut set = grouped(&weights, &[0, 1, 0, 1], &[2, 1]);
    assert_eq!(set.set().len(), 4);
    let distribution = set.distribution();
    assert_eq!(distribution.samples(), &Samples::All);
    assert_eq!(distribution.log_weights(), Some(&weights[..]));
    assert_eq!(distribution.partition().assignments(), [0, 1, 0, 1]);
    assert_eq!(distribution.partition().representatives(), [2, 1]);
    assert_eq!(distribution.temperature(), None);
    // Weights of a total of 10.
    assert_close(
        &distribution.probabilities().unwrap(),
        &[0.1, 0.2, 0.3, 0.4],
    );
    assert_close(&distribution.state_probabilities().unwrap(), &[0.4, 0.6]);
    // Kish: 10² / (1 + 4 + 9 + 16) = 10 / 3.
    assert!(close(
        distribution.effective_sample_size().unwrap(),
        10.0 / 3.0
    ));
    assert_eq!(distribution.ranked_by_weight(), [3, 2, 1, 0]);
    assert_eq!(distribution.ranked_by_population(), [1, 0]);
    assert_eq!(distribution.covering_population(0.5).unwrap(), [1]);
    assert_eq!(distribution.covering_population(0.7).unwrap(), [1, 0]);
    assert_eq!(
        distribution.covering_population(0.0),
        Err(StructureSetError::InvalidFraction)
    );
    let collapsed = set.collapse();
    assert_eq!(collapsed.set().len(), 2);
    assert!(collapsed
        .set()
        .get(0)
        .unwrap()
        .shares_storage_with(&set.set().get(2).unwrap()));
    assert!(collapsed
        .set()
        .get(1)
        .unwrap()
        .shares_storage_with(&set.set().get(1).unwrap()));
    assert_eq!(
        collapsed.distribution().partition(),
        &Partition::state_per_sample(2)
    );
    assert_close(&log_weights(&collapsed), &[4.0_f64.ln(), 6.0_f64.ln()]);

    let mut energy_set = energy_conformers(&[3.0, 1.0, 2.0]);
    let state = |energy| StateData::new().with_energy(kj(energy).unwrap()).unwrap();
    let expected = Distribution::new(Samples::All, Partition::state_per_sample(3))
        .unwrap()
        .with_log_weights(vec![0.0; 3])
        .unwrap()
        .with_states(vec![state(3.0), state(1.0), state(2.0)])
        .unwrap();
    assert_eq!(energy_set.set().len(), 3);
    assert_eq!(energy_set.distribution(), &expected);
    let distribution = energy_set.distribution();
    assert_eq!(distribution.ranked_by_energy().unwrap(), [1, 2, 0]);
    // States 1 and 2 lie within 1 kJ/mol of the lowest; the bound is inclusive.
    let window = Quantity::new(1.0, KILOJOULE_PER_MOLE);
    assert_eq!(distribution.within_energy_window(window).unwrap(), [1, 2]);
    assert_eq!(
        distribution.within_energy_window(Quantity::new(-1.0, KILOJOULE_PER_MOLE)),
        Err(StructureSetError::InvalidEnergyWindow)
    );
    let temperature = Quantity::new(300.0, KELVIN);
    energy_set
        .distribution_mut()
        .set_temperature(Some(temperature))
        .unwrap();
    assert_eq!(energy_set.distribution().temperature(), Some(temperature));
    // A temperature changes nothing else.
    assert_eq!(
        energy_set.distribution(),
        &expected.clone().with_temperature(temperature).unwrap()
    );
    assert_eq!(
        energy_set
            .distribution_mut()
            .set_temperature(Some(Quantity::new(0.0, KELVIN))),
        Err(StructureSetError::InvalidTemperature)
    );

    // Regrouping merges whole states into their free energy at the
    // temperature, -RT ln(e^(-3/RT) + e^(-1/RT)); a state that is one old
    // state keeps its energy.
    let rt: f64 = 8.314_462_618_153_24e-3 * 300.0;
    let merged = 1.0 - rt * ((-2.0 / rt).exp() + 1.0).ln();
    let partition = Partition::new(vec![0, 0, 1], vec![1, 2]).unwrap();
    energy_set
        .distribution_mut()
        .regroup(partition.clone())
        .unwrap();
    assert_eq!(
        energy_set
            .states()
            .map(|state| state.energy())
            .collect::<Vec<_>>(),
        [merged, 2.0].map(|energy| Some(Quantity::new(energy, KILOJOULE_PER_MOLE)))
    );
    let regrouped = Distribution::new(Samples::All, partition)
        .unwrap()
        .with_log_weights(vec![0.0; 3])
        .unwrap()
        .with_states(vec![state(merged), state(2.0)])
        .unwrap()
        .with_temperature(temperature)
        .unwrap();
    assert_eq!(energy_set.distribution(), &regrouped);
    assert_eq!(
        energy_set
            .distribution_mut()
            .regroup(Partition::one_state(2)),
        Err(StructureSetError::PartitionMismatch(LengthMismatch::new(
            Countable::Sample,
            3,
            2
        )))
    );
    assert_eq!(energy_set.distribution(), &regrouped);

    // Normalizing changes the log weights alone, to log probabilities.
    let before = set.distribution().clone();
    set.distribution_mut().normalize_log_weights().unwrap();
    let distribution = set.distribution();
    assert_close(
        &exponentials(distribution.log_weights().unwrap()),
        &[0.1, 0.2, 0.3, 0.4],
    );
    assert_eq!(distribution.samples(), before.samples());
    assert_eq!(distribution.partition(), before.partition());
    assert_eq!(distribution.states(), before.states());
    assert_eq!(distribution.temperature(), None);
}
