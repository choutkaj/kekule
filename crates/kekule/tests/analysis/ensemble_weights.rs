//! The log-weight arithmetic of distributions and the bookkeeping of
//! ensembles, checked against values worked out by hand.

use std::sync::Arc;

use kekule::math::Point3;
use kekule::structure::{
    Distribution, Ensemble, Geometry, Partition, Positions, Samples, StateData, StructureSet,
    StructureSetError, ThermodynamicState, WeightModel,
};
use kekule::topology::Topology;
use kekule::units::{Quantity, KELVIN, KILOJOULE_PER_MOLE, NANOMETER};

/// `ln Σ exp(x)`, shifted by the largest value.
fn log_sum_exp(values: impl Iterator<Item = f64> + Clone) -> f64 {
    let maximum = values.clone().fold(f64::NEG_INFINITY, f64::max);
    maximum
        + values
            .map(|value| (value - maximum).exp())
            .sum::<f64>()
            .ln()
}

/// The molar gas constant in kJ/(mol K), exact since the 2019 SI.
const R: f64 = 8.314_462_618_153_24e-3;

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} differs from {expected}"
    );
}

fn assert_all_close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (actual, expected) in actual.iter().zip(expected) {
        assert_close(*actual, *expected);
    }
}

/// The canonical state at `kelvin`.
fn at(kelvin: f64) -> ThermodynamicState {
    ThermodynamicState::canonical(Quantity::new(kelvin, KELVIN)).unwrap()
}

fn water() -> Arc<Topology> {
    kekule::smiles::to_topology("O").unwrap()
}

/// A set of `count` water geometries, geometry `i` with its oxygen at `x = i` nm.
fn set_of(count: usize) -> StructureSet {
    let topology = water();
    let mut set = StructureSet::new(Arc::clone(&topology));
    for index in 0..count {
        let mut positions = vec![Point3::origin(); topology.atom_count()];
        positions[0].x = index as f64;
        let positions = Positions::new(Quantity::new(positions, NANOMETER)).unwrap();
        set.push(Geometry::new(positions), None).unwrap();
    }
    set
}

fn oxygen_x(ensemble: &Ensemble, index: usize) -> f64 {
    ensemble
        .set()
        .get(index)
        .unwrap()
        .payload()
        .positions()
        .values()
        .value()[0]
        .x
}

#[test]
fn probabilities_are_normalized_weights() -> Result<(), StructureSetError> {
    // Weights 1 and 3.
    let distribution = Distribution::new(Samples::All, Partition::state_per_sample(2))?
        .with_log_weights(vec![0.0, 3.0_f64.ln()])?;
    assert_all_close(&distribution.probabilities()?, &[0.25, 0.75]);
    // Only weight ratios matter: shifting every log weight changes nothing,
    // even far outside the range of a linear f64 weight.
    for shift in [-1000.0, 1000.0] {
        let shifted = Distribution::new(Samples::All, Partition::state_per_sample(2))?
            .with_log_weights(vec![shift, shift + 3.0_f64.ln()])?;
        assert_all_close(&shifted.probabilities()?, &[0.25, 0.75]);
    }
    // Without weights the samples weigh equally.
    let unweighted = Distribution::new(Samples::All, Partition::one_state(4))?;
    assert_eq!(unweighted.probabilities()?, vec![0.25; 4]);
    assert_eq!(unweighted.weight_model(), WeightModel::Equal);
    Ok(())
}

#[test]
fn normalizing_makes_log_weights_log_probabilities() -> Result<(), StructureSetError> {
    let mut ensemble = Ensemble::weighted(set_of(3), vec![10.0, 10.0 + 2.0_f64.ln(), 10.0])?;
    let before = ensemble.distribution().probabilities()?;
    ensemble.distribution_mut().normalize_log_weights()?;
    let WeightModel::LogWeights(log_weights) = ensemble.distribution().weight_model() else {
        panic!("normalized weights are log weights")
    };
    assert_close(log_sum_exp(log_weights.iter().copied()), 0.0);
    assert_all_close(log_weights, &[0.25_f64.ln(), 0.5_f64.ln(), 0.25_f64.ln()]);
    assert_all_close(&ensemble.distribution().probabilities()?, &before);
    Ok(())
}

#[test]
fn boltzmann_populations_follow_state_energies() -> Result<(), StructureSetError> {
    let temperature = 300.0;
    let rt = R * temperature;
    // State 1 lies RT ln 2 above state 0, so it holds half its population.
    let partition = Partition::new(vec![0, 0, 1], vec![0, 2]).unwrap();
    let states = vec![
        StateData::new()
            .with_energy(Quantity::new(0.0, KILOJOULE_PER_MOLE))
            .unwrap(),
        StateData::new()
            .with_energy(Quantity::new(rt * 2.0_f64.ln(), KILOJOULE_PER_MOLE))
            .unwrap(),
    ];
    // Within state 0, the second sample weighs three times the first.
    let mut ensemble = Ensemble::new(
        set_of(3),
        Distribution::new(Samples::All, partition)?
            .with_states(states)?
            .with_log_weights(vec![0.0, 3.0_f64.ln(), 0.0])?,
    )?;
    let WeightModel::LogWeights(given) = ensemble.distribution().weight_model() else {
        panic!("the weights are given")
    };
    let total_before = log_sum_exp(given.iter().copied());
    ensemble
        .distribution_mut()
        .assign_boltzmann_populations(at(temperature))?;
    let distribution = ensemble.distribution();
    assert_all_close(
        &distribution.state_probabilities()?,
        &[2.0 / 3.0, 1.0 / 3.0],
    );
    let WeightModel::Boltzmann(log_weights) = distribution.weight_model() else {
        panic!("the weights are Boltzmann populations")
    };
    assert_close(log_weights[1] - log_weights[0], 3.0_f64.ln());
    assert_close(log_sum_exp(log_weights.iter().copied()), total_before);
    assert_eq!(
        distribution
            .temperature()
            .map(|t| t.value_in(KELVIN).unwrap()),
        Some(temperature)
    );
    Ok(())
}

#[test]
fn a_remote_state_keeps_a_finite_log_population() -> Result<(), StructureSetError> {
    // 10⁵ kJ/mol above the ground state at 300 K: exp(-E/RT) ≈ e^-40000.
    let states = vec![
        StateData::new()
            .with_energy(Quantity::new(0.0, KILOJOULE_PER_MOLE))
            .unwrap(),
        StateData::new()
            .with_energy(Quantity::new(1.0e5, KILOJOULE_PER_MOLE))
            .unwrap(),
    ];
    let mut ensemble = Ensemble::new(
        set_of(2),
        Distribution::new(Samples::All, Partition::state_per_sample(2))?.with_states(states)?,
    )?;
    ensemble
        .distribution_mut()
        .assign_boltzmann_populations(at(300.0))?;
    let gap =
        ensemble.state(1).unwrap().log_population() - ensemble.state(0).unwrap().log_population();
    assert_close(gap, -1.0e5 / (R * 300.0));
    assert_eq!(
        ensemble.distribution().state_probabilities()?,
        vec![1.0, 0.0]
    );
    Ok(())
}

#[test]
fn an_unrepresentable_population_is_refused() -> Result<(), StructureSetError> {
    // 100 kJ/mol above the ground state at 10⁻³⁰⁸ K: E/RT ≈ 10³¹², beyond
    // any finite log weight.
    let states = vec![
        StateData::new()
            .with_energy(Quantity::new(0.0, KILOJOULE_PER_MOLE))
            .unwrap(),
        StateData::new()
            .with_energy(Quantity::new(100.0, KILOJOULE_PER_MOLE))
            .unwrap(),
    ];
    let mut ensemble = Ensemble::new(
        set_of(2),
        Distribution::new(Samples::All, Partition::state_per_sample(2))?
            .with_states(states)?
            .with_log_weights(vec![0.0, 1.0])?,
    )?;
    assert_eq!(
        ensemble
            .distribution_mut()
            .assign_boltzmann_populations(at(1.0e-308)),
        Err(StructureSetError::UnrepresentablePopulation { state: 1 })
    );
    // Failure changes nothing.
    let distribution = ensemble.distribution();
    assert_eq!(
        distribution.weight_model(),
        WeightModel::LogWeights(&[0.0, 1.0])
    );
    assert!(distribution.target().is_none());
    Ok(())
}

#[test]
fn collapse_gives_each_representative_its_state_population() -> Result<(), StructureSetError> {
    let partition = Partition::new(vec![0, 1, 0, 1], vec![0, 3]).unwrap();
    let ensemble = Ensemble::new(
        set_of(4),
        Distribution::new(Samples::All, partition)?.with_log_weights(vec![0.0, 1.0, 2.0, 3.0])?,
    )?;
    let collapsed = ensemble.collapse();
    assert_eq!(collapsed.set().len(), 2);
    assert_eq!(
        (oxygen_x(&collapsed, 0), oxygen_x(&collapsed, 1)),
        (0.0, 3.0)
    );
    let WeightModel::LogWeights(log_weights) = collapsed.distribution().weight_model() else {
        panic!("collapsed weights are log weights")
    };
    assert_close(log_weights[0], log_sum_exp([0.0, 2.0].into_iter()));
    assert_close(log_weights[1], log_sum_exp([1.0, 3.0].into_iter()));
    assert_all_close(
        &collapsed.distribution().state_probabilities()?,
        &ensemble.distribution().state_probabilities()?,
    );
    Ok(())
}

#[test]
fn edits_keep_the_distribution_on_its_set() -> Result<(), StructureSetError> {
    let mut ensemble = Ensemble::weighted(set_of(2), vec![0.0, 1.0])?;
    // A pushed geometry is not a sample; every geometry becomes the earlier ones.
    let pushed = ensemble.push(set_of(3).get(2).unwrap().payload().clone(), None)?;
    assert_eq!(pushed, 2);
    assert_eq!(
        ensemble.distribution().samples(),
        &Samples::Subset(vec![0, 1])
    );
    assert_eq!(ensemble.distribution().log_weight(2), None);
    // A weighted distribution needs the new sample's log weight.
    let geometry = set_of(1).get(0).unwrap().payload().clone();
    assert_eq!(
        ensemble.push_sample(geometry.clone(), None, 0, None),
        Err(StructureSetError::LogWeightRequired)
    );
    assert_eq!(ensemble.set().len(), 3);
    let index = ensemble.push_sample(geometry, None, 0, Some(-2.0))?;
    assert_eq!(ensemble.distribution().log_weight(index), Some(-2.0));
    // Removing a geometry renumbers the samples after it.
    ensemble.remove(2)?;
    assert_eq!(ensemble.distribution().samples(), &Samples::All);
    assert_eq!(
        ensemble.distribution().weight_model(),
        WeightModel::LogWeights(&[0.0, 1.0, -2.0])
    );
    // A new state's sample.
    let state = ensemble.push_state(
        set_of(1).get(0).unwrap().payload().clone(),
        None,
        StateData::new(),
        Some(4.0),
    )?;
    assert_eq!(state, 1);
    assert_eq!(
        ensemble.state(1).unwrap().geometries().collect::<Vec<_>>(),
        vec![3]
    );
    Ok(())
}
