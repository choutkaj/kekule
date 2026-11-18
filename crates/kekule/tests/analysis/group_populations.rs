//! Group populations of distributions that are not reweighted: exact sums
//! of their log weights, normalized against the states.

use kekule::statistics::{CumulantOrder, Estimator};
use kekule::structure::{Distribution, Partition, Samples, StructureSetError};

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} differs from {expected}"
    );
}

/// Four samples in two states, weights 1, 3 | 2, 2.
fn weighted() -> Distribution {
    Distribution::new(
        Samples::All,
        Partition::new(vec![0, 0, 1, 1], vec![0, 2]).unwrap(),
    )
    .unwrap()
    .with_log_weights(vec![0.0, 3.0_f64.ln(), 2.0_f64.ln(), 2.0_f64.ln()])
    .unwrap()
}

#[test]
fn populations_of_unweighted_models_are_exact_sums() -> Result<(), StructureSetError> {
    let distribution = weighted();
    assert_eq!(distribution.estimator(), Estimator::Exponential);
    // States of total weight 4 and 4 of 8.
    assert_close(distribution.log_population(&[0, 1]).unwrap(), 0.5_f64.ln());
    assert_close(distribution.log_population(&[1]).unwrap(), (3.0 / 8.0_f64).ln());
    assert_close(distribution.log_population(&[0, 1, 2, 3]).unwrap(), 0.0);
    let states = distribution.state_probabilities()?;
    for (state, probability) in states.iter().enumerate() {
        let positions = distribution.partition().samples(state).collect::<Vec<_>>();
        assert_close(
            distribution.log_population(&positions).unwrap().exp(),
            *probability,
        );
    }
    // Shares within a state sum to the state's probability.
    let probabilities = distribution.probabilities()?;
    assert_close(probabilities[0] + probabilities[1], states[0]);
    Ok(())
}

#[test]
fn groups_are_normalized_among_themselves() -> Result<(), StructureSetError> {
    let distribution = weighted();
    // Groups that are the states give the state probabilities.
    let states = distribution.group_log_populations(&[Some(0), Some(0), Some(1), Some(1)], 2)?;
    let expected = distribution.state_probabilities()?;
    for (state, probability) in states.iter().zip(expected) {
        assert_close(state.unwrap().exp(), probability);
    }
    // Weights 1 | 3 | none, the third sample in no group, the fourth group empty.
    let bins = distribution.group_log_populations(&[Some(0), Some(1), None, Some(2)], 4)?;
    assert_close(bins[0].unwrap().exp(), 1.0 / 6.0);
    assert_close(bins[1].unwrap().exp(), 3.0 / 6.0);
    assert_close(bins[2].unwrap().exp(), 2.0 / 6.0);
    assert_eq!(bins[3], None);
    Ok(())
}

#[test]
fn groupings_must_fit_the_samples() {
    let distribution = weighted();
    assert!(matches!(
        distribution.group_log_populations(&[Some(0); 3], 1),
        Err(StructureSetError::GroupingMismatch(_))
    ));
    let error = distribution
        .group_log_populations(&[Some(0), Some(2), None, None], 2)
        .unwrap_err();
    assert!(matches!(error, StructureSetError::GroupOutOfRange(_)));
    assert_eq!(
        error.to_string(),
        "Group 2 does not exist. There are 2 groups."
    );
    assert_eq!(distribution.log_population(&[]), None);
    assert_eq!(distribution.log_population(&[4]), None);
    assert_eq!(distribution.group_effective_sample_size(&[]), None);
    assert_eq!(distribution.group_effective_sample_size(&[0, 9]), None);
}

#[test]
fn group_effective_sample_sizes_are_kish() -> Result<(), StructureSetError> {
    let distribution = weighted();
    // Weights 1 and 3: 16 / 10.
    assert_close(distribution.group_effective_sample_size(&[0, 1]).unwrap(), 1.6);
    assert_close(distribution.group_effective_sample_size(&[2, 3]).unwrap(), 2.0);
    // Weights 1, 3, 2, 2: 64 / 18.
    assert_close(distribution.effective_sample_size()?, 64.0 / 18.0);
    let equal = Distribution::new(Samples::All, Partition::one_state(5))?;
    assert_eq!(equal.group_effective_sample_size(&[0, 2, 4]), Some(3.0));
    Ok(())
}

#[test]
fn only_reweighted_samples_take_an_estimator() {
    let mut ensemble = kekule::structure::Ensemble::uniform(kekule::structure::StructureSet::new(
        kekule::smiles::to_topology("O").unwrap(),
    ));
    let error = ensemble
        .distribution_mut()
        .set_estimator(Estimator::Cumulant(CumulantOrder::Second))
        .unwrap_err();
    assert_eq!(error, StructureSetError::NotReweighted);
    assert_eq!(error.to_string(), "The distribution is not reweighted.");
}
