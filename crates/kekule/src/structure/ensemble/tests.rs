//! The log-weight arithmetic of [`Distribution`] and the bookkeeping of
//! [`Ensemble`], checked against values worked out by hand.

use crate::diagnostics::Countable;
use crate::diagnostics::IndexOutOfRange;
use crate::diagnostics::LengthMismatch;
use std::sync::Arc;

use super::{log_sum_exp, Distribution, Ensemble, Samples, StateData};
use crate::math::Point3;
use crate::structure::{Geometry, Partition, Positions, StructureSet, StructureSetError};
use crate::topology::Topology;
use crate::units::{Quantity, KELVIN, KILOJOULE_PER_MOLE};

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

fn water() -> Arc<Topology> {
    crate::smiles::to_topology("O").unwrap()
}

/// A set of `count` water geometries, geometry `i` with its oxygen at `x = i` nm.
fn set_of(count: usize) -> StructureSet {
    let topology = water();
    let mut set = StructureSet::new(Arc::clone(&topology));
    for index in 0..count {
        let mut positions = vec![Point3::origin(); topology.atom_count()];
        positions[0].x = index as f64;
        let positions = Positions::new(Quantity::new(positions, crate::units::NANOMETER)).unwrap();
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
fn log_sum_exp_is_the_log_of_the_sum() {
    let values = [1.0_f64.ln(), 2.0_f64.ln(), 3.0_f64.ln()];
    assert_close(log_sum_exp(values.iter().copied()), 6.0_f64.ln());
    // exp(1000) overflows; the shifted sum does not.
    assert_close(
        log_sum_exp([1000.0, 1000.0].into_iter()),
        1000.0 + 2.0_f64.ln(),
    );
    assert_close(
        log_sum_exp([-1000.0, -1000.0].into_iter()),
        -1000.0 + 2.0_f64.ln(),
    );
    assert_eq!(log_sum_exp(std::iter::empty()), f64::NEG_INFINITY);
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
    assert_eq!(unweighted.log_weights(), None);
    Ok(())
}

#[test]
fn a_negligible_weight_has_probability_zero() -> Result<(), StructureSetError> {
    // e^-800 is below the smallest positive f64, yet the log weight is kept.
    let distribution = Distribution::new(Samples::All, Partition::state_per_sample(2))?
        .with_log_weights(vec![0.0, -800.0])?;
    assert_eq!(distribution.probabilities()?, vec![1.0, 0.0]);
    assert_eq!(distribution.log_weight(1), Some(-800.0));
    Ok(())
}

#[test]
fn sample_probabilities_name_geometries() -> Result<(), StructureSetError> {
    // Samples 1 and 3 of four geometries, weights 1 and 4.
    let mut distribution = Distribution::new(Samples::Subset(vec![1, 3]), Partition::one_state(2))?
        .with_log_weights(vec![0.0, 4.0_f64.ln()])?;
    distribution.bind(4)?;
    let pairs = distribution.sample_probabilities()?;
    assert_eq!(
        pairs.iter().map(|(g, _)| *g).collect::<Vec<_>>(),
        vec![1, 3]
    );
    assert_all_close(
        &pairs.iter().map(|(_, p)| *p).collect::<Vec<_>>(),
        &[0.2, 0.8],
    );
    assert_eq!(distribution.log_weight(0), None);
    assert_eq!(distribution.sample_position(3), Some(1));
    Ok(())
}

#[test]
fn effective_sample_size_is_kish() -> Result<(), StructureSetError> {
    let unweighted = Distribution::new(Samples::All, Partition::one_state(4))?;
    assert_eq!(unweighted.effective_sample_size()?, 4.0);
    // Weights 1, 1, 2: (Σw)² / Σw² = 16 / 6.
    let log_weights = [0.0, 0.0, 2.0_f64.ln()];
    for shift in [0.0, 700.0, -700.0] {
        let distribution = Distribution::new(Samples::All, Partition::one_state(3))?
            .with_log_weights(log_weights.iter().map(|w| w + shift).collect())?;
        assert_close(distribution.effective_sample_size()?, 16.0 / 6.0);
    }
    // Equal log weights: as many effective samples as samples.
    let equal = Distribution::new(Samples::All, Partition::one_state(5))?
        .with_log_weights(vec![-3.0; 5])?;
    assert_close(equal.effective_sample_size()?, 5.0);
    // A common shift too large for the unshifted sums to stay finite or
    // exact leaves the count unchanged.
    for offset in [1.0e16, 1.0e308, -1.0e308] {
        let shifted = Distribution::new(Samples::All, Partition::one_state(2))?
            .with_log_weights(vec![offset; 2])?;
        assert_eq!(shifted.effective_sample_size()?, 2.0);
    }
    Ok(())
}

#[test]
fn state_populations_sum_their_samples() -> Result<(), StructureSetError> {
    // States {0, 1} and {2}; weights 1, 1, 2, so each state holds 2 of 4.
    // Both log populations are ln 2 computed the same way, so they tie bit
    // for bit.
    let partition = Partition::new(vec![0, 0, 1], vec![0, 2]).unwrap();
    let ensemble = Ensemble::new(
        set_of(3),
        Distribution::new(Samples::All, partition)?.with_log_weights(vec![
            0.0,
            0.0,
            2.0_f64.ln(),
        ])?,
    )?;
    let distribution = ensemble.distribution();
    assert_eq!(distribution.state_probabilities()?, vec![0.5, 0.5]);
    assert_close(ensemble.state(0).unwrap().log_population(), 2.0_f64.ln());
    assert_close(ensemble.state(1).unwrap().log_population(), 2.0_f64.ln());
    assert_eq!(
        ensemble.state(0).unwrap().geometries().collect::<Vec<_>>(),
        vec![0, 1]
    );
    assert_eq!(ensemble.state(1).unwrap().representative(), 2);
    // Ties keep state and geometry order; the heaviest sample ranks first.
    assert_eq!(distribution.ranked_by_population(), vec![0, 1]);
    assert_eq!(distribution.ranked_by_weight(), vec![2, 0, 1]);
    assert_eq!(distribution.covering_population(0.5)?, vec![0]);
    assert_eq!(distribution.covering_population(0.51)?, vec![0, 1]);
    Ok(())
}

#[test]
fn normalizing_makes_log_weights_log_probabilities() -> Result<(), StructureSetError> {
    let mut ensemble = Ensemble::weighted(set_of(3), vec![10.0, 10.0 + 2.0_f64.ln(), 10.0])?;
    let before = ensemble.distribution().probabilities()?;
    ensemble.distribution_mut().normalize_log_weights()?;
    let log_weights = ensemble.distribution().log_weights().unwrap().to_vec();
    assert_close(log_sum_exp(log_weights.iter().copied()), 0.0);
    assert_all_close(&log_weights, &[0.25_f64.ln(), 0.5_f64.ln(), 0.25_f64.ln()]);
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
    let total_before = log_sum_exp(
        ensemble
            .distribution()
            .log_weights()
            .unwrap()
            .iter()
            .copied(),
    );
    ensemble
        .distribution_mut()
        .assign_boltzmann_populations(Quantity::new(temperature, KELVIN))?;
    let distribution = ensemble.distribution();
    assert_all_close(
        &distribution.state_probabilities()?,
        &[2.0 / 3.0, 1.0 / 3.0],
    );
    let log_weights = distribution.log_weights().unwrap();
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
        .assign_boltzmann_populations(Quantity::new(300.0, KELVIN))?;
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
            .assign_boltzmann_populations(Quantity::new(1.0e-308, KELVIN)),
        Err(StructureSetError::UnrepresentablePopulation { state: 1 })
    );
    // Failure changes nothing.
    let distribution = ensemble.distribution();
    assert_eq!(distribution.log_weights(), Some(&[0.0, 1.0][..]));
    assert!(distribution.temperature().is_none());
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
    let log_weights = collapsed.distribution().log_weights().unwrap();
    assert_close(log_weights[0], log_sum_exp([0.0, 2.0].into_iter()));
    assert_close(log_weights[1], log_sum_exp([1.0, 3.0].into_iter()));
    assert_all_close(
        &collapsed.distribution().state_probabilities()?,
        &ensemble.distribution().state_probabilities()?,
    );
    Ok(())
}

#[test]
fn select_remaps_samples_and_keeps_log_weights() -> Result<(), StructureSetError> {
    let mut distribution =
        Distribution::new(Samples::Subset(vec![0, 2, 3]), Partition::one_state(3))?
            .with_log_weights(vec![0.5, 1.5, 2.5])?;
    distribution.bind(4)?;
    let ensemble = Ensemble::new(set_of(4), distribution)?;
    // Keep geometries 3, 1, 2 in that order: 3 and 2 are samples, 1 is not.
    let selected = ensemble.select([3, 1, 2])?;
    assert_eq!(
        (0..3).map(|i| oxygen_x(&selected, i)).collect::<Vec<_>>(),
        vec![3.0, 1.0, 2.0]
    );
    let distribution = selected.distribution();
    assert_eq!(distribution.log_weight(0), Some(2.5));
    assert_eq!(distribution.log_weight(1), None);
    assert_eq!(distribution.log_weight(2), Some(1.5));
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
        ensemble.distribution().log_weights(),
        Some(&[0.0, 1.0, -2.0][..])
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

#[test]
fn invalid_log_weights_are_refused() -> Result<(), StructureSetError> {
    let distribution = Distribution::new(Samples::All, Partition::one_state(2))?;
    assert_eq!(
        distribution
            .clone()
            .with_log_weights(vec![0.0])
            .unwrap_err(),
        StructureSetError::LogWeightCountMismatch(LengthMismatch::new(Countable::LogWeight, 2, 1))
    );
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            distribution
                .clone()
                .with_log_weights(vec![0.0, bad])
                .unwrap_err(),
            StructureSetError::NonFiniteLogWeight { geometry: 1 }
        );
    }
    // Negative and large log weights are ordinary weights.
    assert!(distribution.with_log_weights(vec![-750.0, 750.0]).is_ok());
    // A subset naming a geometry past the set does not bind to it.
    let outside = Distribution::new(Samples::Subset(vec![0, 5]), Partition::one_state(2))?;
    assert!(matches!(
        Ensemble::new(set_of(3), outside),
        Err(StructureSetError::IndexOutOfRange(payload)) if payload == IndexOutOfRange::new(Countable::Geometry, 5, 3)
    ));
    Ok(())
}
