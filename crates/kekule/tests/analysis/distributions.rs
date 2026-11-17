//! Distributions and their ensembles: weights, probabilities, rankings,
//! energies, and the free energies of regrouped states.

use std::sync::Arc;

use kekule::diagnostics::{Countable, LengthMismatch};
use kekule::math::{PeriodicCell, Point3, Vector3};
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{
    DenseArrayError, Distribution, Ensemble, Forces, Geometry, GeometryError, Order, Partition,
    Place, Positions, Requirement, Samples, StateData, Structure, StructureSet, StructureSetError,
    ThermodynamicState, Velocities, WeightModel,
};
use kekule::topology::{Topology, TopologyAtomIndex, TopologyBuilder};
use kekule::units::{Quantity, ANGSTROM, BAR, KELVIN, KILOJOULE_PER_MOLE, NANOMETER};

const ATOM0: TopologyAtomIndex = TopologyAtomIndex::new(0);

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

/// The canonical state at `kelvin`.
fn at(kelvin: f64) -> ThermodynamicState {
    ThermodynamicState::canonical(Quantity::new(kelvin, KELVIN)).unwrap()
}

#[test]
fn regrouping_keeps_weights_and_gives_merged_states_their_free_energy() {
    // Four conformers, energies 0..3 kJ/mol, weights 1..4.
    let weights = ln(&[1.0, 2.0, 3.0, 4.0]);
    let mut set = energy_conformers(&[0.0, 1.0, 2.0, 3.0]);
    let mut distribution = set.distribution_mut();
    distribution.set_log_weights(weights.clone()).unwrap();
    distribution
        .state_mut(1)
        .unwrap()
        .properties_mut()
        .insert(key("label"), PropertyValue::Int(1))
        .unwrap();
    distribution.set_target(Some(at(300.0))).unwrap();
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
    let temperature = at(310.0);
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
    set.distribution_mut().set_target(Some(at(400.0))).unwrap();
    assert!((kilojoules(&set, 0) - free(&energies_kj[..3])).abs() < 1e-12);

    // Without a temperature, or with a member lacking an energy, a merged
    // state has none; a state of one old state keeps its energy.
    let mut cold = source.clone();
    cold.distribution_mut().set_target(None).unwrap();
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
        .set_target(Some(temperature))
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
    assert_eq!(
        distribution.weight_model(),
        WeightModel::LogWeights(&weights)
    );
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
    energy_set
        .distribution_mut()
        .set_target(Some(at(300.0)))
        .unwrap();
    assert_eq!(
        energy_set.distribution().temperature(),
        Some(Quantity::new(300.0, KELVIN))
    );
    // A target changes nothing else.
    assert_eq!(
        energy_set.distribution(),
        &expected.clone().with_target(at(300.0)).unwrap()
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
        .with_target(at(300.0))
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
    assert!(matches!(
        distribution.weight_model(),
        WeightModel::LogWeights(_)
    ));
    assert_close(&exponentials(&log_weights(&set)), &[0.1, 0.2, 0.3, 0.4]);
    assert_eq!(distribution.samples(), before.samples());
    assert_eq!(distribution.partition(), before.partition());
    assert_eq!(distribution.states(), before.states());
    assert_eq!(distribution.target(), None);
}

#[test]
fn the_weight_model_says_how_the_weights_were_assigned() {
    // Boltzmann populations record their target.
    let room = at(298.15);
    let mut set = energy_conformers(&[0.0, 1.0, 2.0]);
    set.distribution_mut()
        .assign_boltzmann_populations(room)
        .unwrap();
    let boltzmann = log_weights(&set);
    assert_eq!(
        set.distribution().weight_model(),
        WeightModel::Boltzmann(&boltzmann)
    );
    assert_eq!(set.distribution().target(), Some(room));
    assert_eq!(temperature(&set), Some(Quantity::new(298.15, KELVIN)));

    // Edits to samples, states, and energies keep the tag: normalizing,
    // regrouping, selecting, removing, adding a sample, collapsing, and
    // setting the same target.
    let is_boltzmann =
        |set: &Ensemble| matches!(set.distribution().weight_model(), WeightModel::Boltzmann(_));
    set.distribution_mut().normalize_log_weights().unwrap();
    set.distribution_mut().set_target(Some(room)).unwrap();
    set.distribution_mut()
        .regroup(Partition::new(vec![0, 0, 1], vec![0, 2]).unwrap())
        .unwrap();
    assert!(is_boltzmann(&set));
    assert!(is_boltzmann(&set.select([2, 0]).unwrap()));
    let mut removed = set.clone();
    removed.remove(1).unwrap();
    assert!(is_boltzmann(&removed));
    let mut grown = set.clone();
    grown
        .push_sample(Geometry::new(single_position(0.0)), None, 0, Some(0.0))
        .unwrap();
    assert!(is_boltzmann(&grown));
    assert!(is_boltzmann(&set.collapse()));

    // Setting a weight or a target at another temperature leaves plain log
    // weights of the same values.
    let normalized = log_weights(&set);
    let mut warmer = set.clone();
    warmer
        .distribution_mut()
        .set_target(Some(at(310.0)))
        .unwrap();
    assert_eq!(
        warmer.distribution().weight_model(),
        WeightModel::LogWeights(&normalized)
    );
    let mut edited = set.clone();
    edited.distribution_mut().set_log_weight(0, 0.0).unwrap();
    assert!(matches!(
        edited.distribution().weight_model(),
        WeightModel::LogWeights(_)
    ));

    // Equal weights keep the target.
    set.distribution_mut().set_equal_weights();
    assert_eq!(set.distribution().weight_model(), WeightModel::Equal);
    assert_eq!(set.distribution().target(), Some(room));

    // Populations follow the target's temperature; a pressure is recorded.
    let isobaric = ThermodynamicState::isothermal_isobaric(
        Quantity::new(298.15, KELVIN),
        Quantity::new(1.0, BAR),
    )
    .unwrap();
    let mut npt = energy_conformers(&[0.0, 1.0, 2.0]);
    npt.distribution_mut()
        .assign_boltzmann_populations(isobaric)
        .unwrap();
    assert_eq!(log_weights(&npt), boltzmann);
    assert_eq!(npt.distribution().target(), Some(isobaric));
    // A pressure alone does not change Boltzmann populations.
    npt.distribution_mut().set_target(Some(room)).unwrap();
    assert!(is_boltzmann(&npt));
}

fn structure_fixture() -> (
    Structure,
    kekule::topology::InstanceAtomId,
    kekule::topology::InstanceBondId,
) {
    let mut editor = MoleculeEditor::new();
    let carbon = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let oxygen = editor
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    let bond = editor.add_bond(carbon, oxygen, BondOrder::Single).unwrap();
    let molecule = editor.finish().unwrap();
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(molecule).unwrap();
    let topology = Arc::new(builder.build().unwrap());
    let atom = kekule::topology::InstanceAtomId::new(instance, carbon);
    let positions = Positions::new(Quantity::new(
        [Point3::origin(), Point3::new(1.0, 0.0, 0.0)],
        NANOMETER,
    ))
    .unwrap();
    (
        Structure::new(topology, positions).unwrap(),
        atom,
        kekule::topology::InstanceBondId::new(instance, bond),
    )
}

/// Distinct, equally weighted geometries, each its own state.
fn conformers(
    topology: impl Into<Arc<kekule::topology::Topology>>,
    geometries: impl IntoIterator<Item = Geometry>,
) -> Result<Ensemble, StructureSetError> {
    let set = StructureSet::from_geometries(topology, geometries)?;
    let partition = Partition::state_per_sample(set.len());
    Ensemble::new(set, Distribution::new(Samples::All, partition)?)
}

/// Structures drawn from one distribution: one state of equal weights.
fn bundle(structures: impl IntoIterator<Item = Structure>) -> Result<Ensemble, StructureSetError> {
    Ok(Ensemble::uniform(StructureSet::from_structures(
        structures,
    )?))
}

/// Conformers of `topology` at the origin with these log weights, one state
/// each.
fn weighted(topology: Arc<kekule::topology::Topology>, log_weights: &[f64]) -> Ensemble {
    let geometries = log_weights
        .iter()
        .map(|_| Geometry::new(single_position(0.0)));
    let mut ensemble = conformers(topology, geometries).unwrap();
    ensemble
        .distribution_mut()
        .set_log_weights(log_weights.to_vec())
        .unwrap();
    ensemble
}

/// An empty set of `topology` with a distribution of no samples.
fn with_empty_distribution(topology: impl Into<Arc<kekule::topology::Topology>>) -> Ensemble {
    Ensemble::uniform(StructureSet::new(topology))
}

/// Whether every geometry of `a` shares its structure with the same geometry
/// of `b`.
fn shares_every_structure(a: &StructureSet, b: &StructureSet) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b.iter())
            .all(|(a, b)| a.payload().shares_storage_with(b.payload()))
}

#[test]
fn sample_log_weights_are_finite_and_relative() {
    let topology = Arc::new(single_atom_topology());
    let three = 3.0_f64.ln();
    let mut set = weighted(Arc::clone(&topology), &[three, 0.0]);
    // Any finite log weight, of either sign, is a positive weight.
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            set.distribution_mut().set_log_weight(1, invalid),
            Err(StructureSetError::NonFiniteLogWeight { geometry: 1 })
        );
    }
    assert_eq!(log_weights(&set), [three, 0.0]);

    // Weights are relative, so selection keeps them as stored.
    let resampled = set.select([0, 0, 1]).unwrap();
    assert_eq!(log_weights(&resampled), [three, three, 0.0]);

    // Normalized, each log weight is the logarithm of its probability.
    set.distribution_mut().normalize_log_weights().unwrap();
    assert_close(&exponentials(&log_weights(&set)), &[0.75, 0.25]);

    // Weights beyond the range of `f64`, e¹⁰⁰⁰ each, normalize in log
    // space.
    let mut huge = weighted(Arc::clone(&topology), &[1000.0, 1000.0]);
    assert_eq!(huge.distribution().probabilities().unwrap(), [0.5, 0.5]);
    huge.distribution_mut().normalize_log_weights().unwrap();
    let normalized = log_weights(&huge);
    assert_eq!(normalized[0], normalized[1]);
    // 1000 + ln 2 rounds to the spacing of doubles near 1000, 1.1e-13.
    assert!(
        (normalized[0] + 2.0_f64.ln()).abs() <= 1.2e-13,
        "{normalized:?}"
    );
    assert_eq!(huge.distribution().probabilities().unwrap(), [0.5, 0.5]);
    // A ratio of e¹⁰⁰⁰, beyond the range of linear weights, is kept: the
    // light sample's probability underflows to zero, its log weight does not.
    let mut extreme = weighted(Arc::clone(&topology), &[0.0, -1000.0]);
    extreme.distribution_mut().normalize_log_weights().unwrap();
    assert_eq!(log_weights(&extreme), [0.0, -1000.0]);
    assert_eq!(extreme.distribution().probabilities().unwrap(), [1.0, 0.0]);
    assert_eq!(
        with_empty_distribution(topology)
            .distribution_mut()
            .normalize_log_weights(),
        Err(StructureSetError::EmptyDistribution)
    );

    // Structures carry no statistical weight; they become equally weighted samples.
    let (structure, _, _) = structure_fixture();
    let set = bundle([structure.clone(), structure]).unwrap();
    assert_eq!(set.distribution().weight_model(), WeightModel::Equal);
    assert_eq!(log_weights(&set), [0.0, 0.0]);
}

#[test]
fn boltzmann_populations_follow_relative_state_energies_at_temperature() {
    // 1, 0, and 2 kcal/mol above the minimum at 298.15 K. Expected factors
    // exp(-dE/RT) evaluated independently to 40 digits (Python `decimal`,
    // R = N_A k_B exact in the 2019 SI).
    let mut set = energy_conformers(&[-100.0 + 4.184, -100.0, -100.0 + 8.368]);
    set.distribution_mut()
        .assign_boltzmann_populations(at(298.15))
        .unwrap();
    let expected = [0.184_924_695_965_256_4, 1.0, 0.034_197_143_177_842_51];
    let log_weights = log_weights(&set);
    for (log_weight, expected) in log_weights.iter().zip(expected) {
        let ratio = (log_weight - log_weights[1]).exp();
        assert!((ratio - expected).abs() <= 1.0e-15, "{ratio} != {expected}");
    }
    // The total weight, three samples of log weight 0, is kept.
    let total = exponentials(&log_weights).iter().sum::<f64>();
    assert!((total - 3.0).abs() <= 1.0e-15, "{total}");
    assert_eq!(temperature(&set), Some(Quantity::new(298.15, KELVIN)));

    // Missing energies change nothing.
    let before = self::log_weights(&set);
    let mut distribution = set.distribution_mut();
    distribution.state_mut(2).unwrap().set_energy(None).unwrap();
    assert_eq!(
        distribution.assign_boltzmann_populations(at(310.0)),
        Err(StructureSetError::Missing(Requirement::StateEnergies))
    );
    assert_eq!(self::log_weights(&set), before);
    assert_eq!(temperature(&set), Some(Quantity::new(298.15, KELVIN)));

    // A Boltzmann factor that underflows a linear weight is kept in log
    // space: 1000 MJ/mol above the lowest state, exp(-ΔE/RT) is zero in
    // `f64`, but its logarithm, -ΔE/RT, is not.
    let mut wide = energy_conformers(&[0.0, 1.0e6]);
    wide.distribution_mut()
        .assign_boltzmann_populations(at(298.15))
        .unwrap();
    let wide_log_weights = self::log_weights(&wide);
    // The total, two samples of log weight 0, stays with the lowest state.
    assert_eq!(wide_log_weights[0], 2.0_f64.ln());
    let thermal_energy = 8.314_462_618_153_24e-3 * 298.15;
    let gap = wide_log_weights[1] - wide_log_weights[0];
    // -ΔE/RT is about -4.0e5, where doubles are 5.8e-11 apart.
    assert!(
        (gap + 1.0e6 / thermal_energy).abs() <= 1.0e-9,
        "{gap} != {}",
        -1.0e6 / thermal_energy
    );
    assert_eq!(wide.distribution().probabilities().unwrap(), [1.0, 0.0]);
    assert_eq!(temperature(&wide), Some(Quantity::new(298.15, KELVIN)));
}

#[test]
fn boltzmann_populations_scale_the_samples_of_each_state_together() {
    // State 0 (0 kJ/mol) holds geometries of weights 1 and 3; state 1, RT ln 3
    // higher, one geometry of weight 2. Populations go from 4 : 2 to 3 : 1 of
    // the kept total 6, so 4.5 and 1.5: state 0's geometries scale by 4.5 / 4
    // and state 1's by 1.5 / 2.
    let kelvin = 300.0;
    let gap = 8.314_462_618_153_24e-3 * kelvin * 3.0_f64.ln();
    let mut set = grouped(&ln(&[1.0, 2.0, 3.0]), &[0, 1, 0], &[0, 1]);
    set_energies(&mut set, &[Some(0.0), Some(gap)]);
    set.distribution_mut()
        .assign_boltzmann_populations(at(kelvin))
        .unwrap();
    for (weight, expected) in exponentials(&log_weights(&set))
        .into_iter()
        .zip([1.125, 1.5, 3.375])
    {
        assert!(
            (weight - expected).abs() <= 1.0e-14,
            "{weight} != {expected}"
        );
    }
    assert_eq!(set.distribution().partition().assignments(), [0, 1, 0]);

    // States whose weights, e⁻¹⁰⁰⁰ and e¹⁰⁰⁰, differ by more than the
    // floating-point range get equal populations of half the total, which
    // is e¹⁰⁰⁰ to the last place.
    let mut spread = grouped(&[-1000.0, 1000.0], &[0, 1], &[0, 1]);
    set_energies(&mut spread, &[Some(0.0), Some(0.0)]);
    spread
        .distribution_mut()
        .assign_boltzmann_populations(at(kelvin))
        .unwrap();
    assert_eq!(log_weights(&spread), [1000.0 - 2.0_f64.ln(); 2]);
    assert_eq!(spread.distribution().probabilities().unwrap(), [0.5, 0.5]);

    // WeightModel one state changes no weight: its population is all there is.
    let three = ln(&[1.0, 2.0, 5.0]);
    let mut one = grouped(&three, &[0, 0, 0], &[1]);
    set_energies(&mut one, &[Some(7.0)]);
    one.distribution_mut()
        .assign_boltzmann_populations(at(298.15))
        .unwrap();
    assert_close(&log_weights(&one), &three);
    assert_eq!(temperature(&one), Some(Quantity::new(298.15, KELVIN)));

    // Geometries drawn from a simulation form a state without an energy, which
    // is never weighted a second time.
    let (structure, _, _) = structure_fixture();
    let mut drawn = bundle([structure.clone(), structure]).unwrap();
    let mut distribution = drawn.distribution_mut();
    assert_eq!(
        distribution.assign_boltzmann_populations(at(298.15)),
        Err(StructureSetError::Missing(Requirement::StateEnergies))
    );
    assert_eq!(distribution.target(), None);

    // A target is recorded or cleared explicitly.
    distribution.set_target(Some(at(298.15))).unwrap();
    assert_eq!(
        distribution.temperature(),
        Some(Quantity::new(298.15, KELVIN))
    );
    distribution.set_target(None).unwrap();
    assert_eq!(temperature(&drawn), None);
}

#[test]
fn regrouping_and_reweighting_a_shared_ensemble_copy_no_structure() {
    let shared = Arc::new(energy_conformers(&[0.0, 1.0, 2.0, 3.0]));
    let mut edited = Arc::clone(&shared);

    let mut distribution = Arc::make_mut(&mut edited).distribution_mut();
    distribution
        .regroup(Partition::new(vec![0, 0, 1, 1], vec![0, 2]).unwrap())
        .unwrap();
    distribution.set_log_weight(1, 3.0_f64.ln()).unwrap();
    distribution.normalize_log_weights().unwrap();
    distribution
        .state_mut(0)
        .unwrap()
        .set_energy(kj(0.0))
        .unwrap();
    distribution
        .state_mut(1)
        .unwrap()
        .set_energy(kj(1.0))
        .unwrap();
    distribution
        .assign_boltzmann_populations(at(300.0))
        .unwrap();

    assert!(shares_every_structure(edited.set(), shared.set()));
    // The held ensemble keeps its own grouping and weights.
    assert_eq!(shared.distribution().state_count(), 4);
    assert_eq!(log_weights(&shared), [0.0; 4]);
    assert_eq!(temperature(&shared), None);
    assert_eq!(edited.distribution().state_count(), 2);
    assert_ne!(log_weights(&edited), log_weights(&shared));

    // Selections and collapsing share the structures they keep.
    let selected = shared.select([3, 1]).unwrap();
    assert!(selected
        .set()
        .get(0)
        .unwrap()
        .payload()
        .shares_storage_with(shared.set().get(3).unwrap().payload()));
    assert!(shares_every_structure(
        shared.select_states(0..4).unwrap().set(),
        shared.set()
    ));
    assert!(shares_every_structure(
        shared.collapse().set(),
        shared.set()
    ));

    // A geometry taken from a set of the same topology stays shared.
    let mut other = with_empty_distribution(shared.set().shared_topology());
    let geometry = shared.set().get(2).unwrap().payload().clone();
    other
        .push_state(geometry, None, StateData::new(), None)
        .unwrap();
    assert!(other
        .set()
        .get(0)
        .unwrap()
        .payload()
        .shares_storage_with(shared.set().get(2).unwrap().payload()));
}

#[test]
fn collapse_keeps_one_representative_per_state_with_its_population() {
    let geometries = StructureSet::from_geometries(
        single_atom_topology(),
        [1.0, 2.0, 3.0, 4.0].map(|x| Geometry::new(single_position(x))),
    )
    .unwrap();
    let distribution = Distribution::new(
        Samples::All,
        Partition::new(vec![0, 1, 0, 1], vec![2, 1]).unwrap(),
    )
    .unwrap()
    .with_log_weights(ln(&[1.0, 2.0, 3.0, 4.0]))
    .unwrap();
    let mut set = Ensemble::new(geometries, distribution).unwrap();
    set_energies(&mut set, &[Some(0.0), None]);
    set.distribution_mut().set_target(Some(at(300.0))).unwrap();
    let source = set.clone();

    let collapsed = set.collapse();
    assert_eq!(collapsed.set().len(), 2);
    // Populations 1 + 3 and 2 + 4 become the representatives' weights.
    assert_close(&log_weights(&collapsed), &ln(&[4.0, 6.0]));
    assert_eq!(
        collapsed.set().get(0).unwrap().positions(),
        set.set().get(2).unwrap().positions()
    );
    assert_eq!(
        collapsed.set().get(1).unwrap().positions(),
        set.set().get(1).unwrap().positions()
    );
    assert_eq!(
        collapsed.distribution().partition(),
        &Partition::state_per_sample(2)
    );
    assert_eq!(energies(&collapsed), [kj(0.0), None]);
    assert_eq!(temperature(&collapsed), temperature(&set));
    for (collapsed, source) in collapsed
        .distribution()
        .state_probabilities()
        .unwrap()
        .into_iter()
        .zip(set.distribution().state_probabilities().unwrap())
    {
        assert!((collapsed - source).abs() <= 1.0e-15);
    }
    // The source is unchanged.
    assert_eq!(set.distribution(), source.distribution());
    let empty = with_empty_distribution(single_atom_topology()).collapse();
    assert_eq!(empty.set().len(), 0);
    assert_eq!(empty.distribution().state_count(), 0);
}

fn cell() -> PeriodicCell {
    PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(10.0, 10.0, 10.0), ANGSTROM),
        [true; 3],
    )
    .unwrap()
}

#[test]
fn member_views_project_borrowed_and_owned_structures_in_stable_order() {
    let topology = Arc::new(single_atom_topology());
    let score = key("score");
    let mut first = Geometry::new(single_position(1.0));
    first.set_cell(Some(cell()));
    let mut properties = first.properties_mut();
    properties
        .owner_mut()
        .insert(key("method"), PropertyValue::String("sampled".into()))
        .unwrap();
    properties
        .atoms_mut()
        .set_value(score.clone(), ATOM0, Some(PropertyValue::Int(7)))
        .unwrap();
    let second = Geometry::new(single_position(2.0));
    let mut ensemble = conformers(Arc::clone(&topology), [first, second]).unwrap();
    let log_weights = vec![0.25_f64.ln(), 0.75_f64.ln()];
    ensemble
        .distribution_mut()
        .set_log_weights(log_weights.clone())
        .unwrap();

    assert_eq!(
        ensemble
            .set()
            .iter()
            .map(|member| member.positions().values().value()[0].x)
            .collect::<Vec<_>>(),
        [1.0, 2.0]
    );
    assert_eq!(self::log_weights(&ensemble), log_weights);

    let member = ensemble.set().get(0).unwrap();
    assert!(Arc::ptr_eq(&member.shared_topology(), &topology));
    let borrowed = member.as_structure_view();
    assert!(Arc::ptr_eq(&borrowed.shared_topology(), &topology));
    assert_eq!(
        borrowed.positions().values().value().as_ptr(),
        member.positions().values().value().as_ptr()
    );
    let owned = member.to_structure();
    assert!(Arc::ptr_eq(&owned.shared_topology(), &topology));
    assert_eq!(owned.positions(), member.positions());
    assert_eq!(owned.cell(), member.cell());
    assert_eq!(owned.properties(), member.properties());

    ensemble
        .get_mut(0)
        .unwrap()
        .properties_mut()
        .atoms_mut()
        .set_value(score.clone(), ATOM0, Some(PropertyValue::Int(9)))
        .unwrap();
    assert_eq!(
        owned.properties().atoms().value(&score, ATOM0).unwrap(),
        Some(PropertyValue::Int(7))
    );
}

#[test]
fn ordered_geometries_join_a_distribution_without_copying_positions() {
    let (structure, _, _) = structure_fixture();
    let mut geometry = Geometry::new(structure.positions().clone());
    geometry.set_velocities(Some(Velocities::zeros(2))).unwrap();
    assert_eq!(
        geometry.set_forces(Some(Forces::zeros(3))),
        Err(GeometryError::Array(Box::new(
            DenseArrayError::CountMismatch(LengthMismatch::new(Countable::Row, 2, 3))
        )))
    );
    let mut set = StructureSet::from_geometries(structure.shared_topology(), [geometry]).unwrap();
    set.set_order(Some(Order::Step(vec![4]))).unwrap();
    let pointer = set.get(0).unwrap().positions().values().value().as_ptr();
    let ensemble = Ensemble::new(
        set,
        Distribution::new(Samples::All, Partition::one_state(1)).unwrap(),
    )
    .unwrap();
    let member = ensemble.set().get(0).unwrap();
    // Geometries join as equally weighted samples and keep their place.
    assert_eq!(ensemble.distribution().weight_model(), WeightModel::Equal);
    assert_eq!(ensemble.distribution().log_weight(0), Some(0.0));
    assert_eq!(member.place(), Some(Place::Step(4)));
    assert_eq!(member.positions().values().value().as_ptr(), pointer);
}
