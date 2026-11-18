//! The log-weight arithmetic of [`Distribution`] and the bookkeeping of
//! [`Ensemble`], checked against values worked out by hand.

use crate::diagnostics::{Countable, IndexOutOfRange, LengthMismatch};
use std::sync::Arc;

use super::distribution::Weights;
use super::{Distribution, Ensemble, Samples};
use crate::math::Point3;
use crate::statistics::{CumulantOrder, Estimator};
use crate::structure::{Geometry, Partition, Positions, StructureSet, StructureSetError};
use crate::topology::Topology;
use crate::units::Quantity;

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

#[test]
fn a_large_common_offset_keeps_group_normalizations() -> Result<(), StructureSetError> {
    // One sample per state: each of the two has probability one half.
    let distribution = Distribution::new(Samples::All, Partition::state_per_sample(2))?
        .with_log_weights(vec![1.0e16; 2])?;
    let half = 0.5_f64.ln();
    for bin in distribution.group_log_populations(&[Some(0), Some(1)], 2)? {
        assert_close(bin.unwrap(), half);
    }
    assert_close(distribution.log_population(&[0]).unwrap(), half);
    assert_eq!(distribution.state_probabilities()?, vec![0.5, 0.5]);
    Ok(())
}

#[test]
fn large_reduced_biases_keep_exact_shares() -> Result<(), StructureSetError> {
    let cumulant = Weights::Reweighted {
        estimator: Estimator::Cumulant(CumulantOrder::Second),
        reduced: vec![0.0, 1.0e10, 1.0e10 + 1.0],
    };
    let distribution =
        Distribution::new(Samples::All, Partition::one_state(3))?.with_weights(cumulant)?;
    let probabilities = distribution.probabilities()?;
    let share = 1.0 / (1.0 + (-1.0_f64).exp());
    assert_eq!(probabilities[0], 0.0);
    assert_close(probabilities[1], 1.0 - share);
    assert_close(probabilities[2], share);
    // Log weights are log probabilities.
    assert_close(distribution.log_weight(2).unwrap(), share.ln());

    // Two states at a common offset of 1e16: relative values 0, 2 and 4, 2.
    let offset = 1.0e16;
    let reduced = [0.0, 2.0, 4.0, 2.0].map(|value| offset + value).to_vec();
    let partition = Partition::new(vec![0, 0, 1, 1], vec![0, 2]).unwrap();
    let distribution =
        Distribution::new(Samples::All, partition)?.with_weights(Weights::Reweighted {
            estimator: Estimator::Cumulant(CumulantOrder::Second),
            reduced,
        })?;
    // ln 2 + mean + variance / 2 relative to the offset: 1.5 and 3.5.
    let low = 1.0 / (1.0 + 2.0_f64.exp());
    let states = distribution.state_probabilities()?;
    assert_close(states[0], low);
    assert_close(states[1], 1.0 - low);
    assert_close(distribution.log_population(&[0, 1]).unwrap(), low.ln());
    let bins = distribution.group_log_populations(&[Some(0), Some(0), Some(1), Some(1)], 2)?;
    assert_close(bins[0].unwrap(), low.ln());
    assert_close(bins[1].unwrap(), (1.0 - low).ln());
    // Exact shares within each state.
    let probabilities = distribution.probabilities()?;
    let share = 1.0 / (1.0 + 2.0_f64.exp());
    assert_close(probabilities[0], low * share);
    assert_close(probabilities[1], low * (1.0 - share));
    assert_close(probabilities[2], (1.0 - low) * (1.0 - share));
    assert_close(probabilities[3], (1.0 - low) * share);
    Ok(())
}

#[test]
fn edits_clear_the_cached_populations() -> Result<(), StructureSetError> {
    let partition = Partition::new(vec![0, 0, 1], vec![0, 2]).unwrap();
    let mut ensemble = Ensemble::new(
        set_of(3),
        Distribution::new(Samples::All, partition)?.with_log_weights(vec![0.0; 3])?,
    )?;
    assert_all_close(
        &ensemble.distribution().state_probabilities()?,
        &[2.0 / 3.0, 1.0 / 3.0],
    );
    ensemble
        .distribution_mut()
        .set_log_weight(2, 2.0_f64.ln())?;
    assert_all_close(&ensemble.distribution().state_probabilities()?, &[0.5, 0.5]);
    ensemble
        .distribution_mut()
        .regroup(Partition::one_state(3))
        .unwrap();
    assert_all_close(&ensemble.distribution().state_probabilities()?, &[1.0]);
    ensemble.remove(2)?;
    assert_all_close(&ensemble.distribution().probabilities()?, &[0.5, 0.5]);
    ensemble
        .distribution_mut()
        .set_log_weights(vec![0.0, 3.0_f64.ln()])?;
    assert_all_close(&ensemble.distribution().probabilities()?, &[0.25, 0.75]);
    // A copy compares equal whether or not its populations are cached.
    let copy = ensemble.distribution().clone();
    assert_eq!(&copy, ensemble.distribution());
    Ok(())
}
