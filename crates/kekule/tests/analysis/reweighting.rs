//! Reweighting a biased run from its sampling record: reduced biases,
//! estimator-based populations with exact shares, refusals, and the
//! operations that keep or drop the reweighting.

use kekule::native::{self, NativeReadOptions, NativeWriteOptions};
use kekule::properties::{PropertyColumn, PropertyKey};
use kekule::statistics::{CumulantOrder, Estimator};
use kekule::structure::{
    Distribution, Ensemble, Geometry, Order, Partition, Positions, Requirement, Samples, Sampling,
    StateData, StructureSet, StructureSetError, ThermodynamicState, WeightModel,
};
use kekule::units::{Quantity, KELVIN, KILOCALORIE_PER_MOLE, KILOJOULE_PER_MOLE};

/// The molar gas constant in kJ/(mol K), exact since the 2019 SI.
const R: f64 = 8.314_462_618_153_24e-3;

const CUMULANT: Estimator = Estimator::Cumulant(CumulantOrder::Second);

fn assert_close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} differs from {expected}"
    );
}

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

fn kelvin(value: f64) -> Quantity<f64> {
    Quantity::new(value, KELVIN)
}

fn at(value: f64) -> ThermodynamicState {
    ThermodynamicState::canonical(kelvin(value)).unwrap()
}

fn energies(values: &[Option<f64>]) -> PropertyColumn {
    PropertyColumn::Real {
        unit: KILOJOULE_PER_MOLE,
        values: values.to_vec(),
    }
}

/// A run at 300 K from step 1000: geometries at `steps`, the potential
/// energy, and two bias columns, one of them in kcal/mol.
fn biased_run(
    steps: &[u64],
    potential: &[Option<f64>],
    boost: &[Option<f64>],
    dihedral: &[Option<f64>],
) -> StructureSet {
    let topology = kekule::smiles::to_topology("O").unwrap();
    let atoms = topology.atom_count();
    let geometries = steps.iter().map(|_| Geometry::new(Positions::zeros(atoms)));
    let mut set = StructureSet::from_geometries(topology, geometries).unwrap();
    set.set_order(Some(Order::Step(steps.to_vec()))).unwrap();
    let mut columns = set.columns_mut();
    columns
        .insert(key("potential"), energies(potential))
        .unwrap();
    columns.insert(key("boost"), energies(boost)).unwrap();
    columns
        .insert(
            key("dihedral"),
            PropertyColumn::Real {
                unit: KILOCALORIE_PER_MOLE,
                values: dihedral.to_vec(),
            },
        )
        .unwrap();
    let sampling = Sampling::new(at(300.0))
        .with_from_step(1_000)
        .with_potential_energy(key("potential"))
        .unwrap()
        .with_bias([key("boost"), key("dihedral")])
        .unwrap();
    set.set_sampling(Some(sampling)).unwrap();
    set
}

/// Six complete samples from step 1000 on, after one geometry before it.
fn six() -> StructureSet {
    let steps = [500, 1_000, 1_500, 2_000, 2_500, 3_000, 3_500];
    let potential = [-10.0, -12.0, -9.0, -11.5, -10.5, -8.0, -12.5].map(Some);
    let boost = [9.0, 0.0, 1.0, 2.0, 0.5, 3.0, 1.5].map(Some);
    let dihedral = [9.0, 0.1, 0.0, 0.3, 0.2, 0.0, 0.4].map(Some);
    biased_run(&steps, &potential, &boost, &dihedral)
}

fn reduced(distribution: &Distribution) -> Vec<f64> {
    match distribution.weight_model() {
        WeightModel::Reweighted { reduced, .. } => reduced.to_vec(),
        other => panic!("{other:?} is not reweighted"),
    }
}

/// `ln n + κ1 + κ2 / 2` of `values`, with population moments.
fn second_cumulant(values: &[f64]) -> f64 {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    let variance = values
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / n;
    n.ln() + mean + variance / 2.0
}

fn log_sum_exp(values: &[f64]) -> f64 {
    let maximum = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    maximum
        + values
            .iter()
            .map(|value| (value - maximum).exp())
            .sum::<f64>()
            .ln()
}

#[test]
fn reduced_biases_sum_the_bias_columns_at_the_sampled_temperature() {
    let ensemble = Ensemble::reweighted(six(), None, CUMULANT).unwrap();
    let distribution = ensemble.distribution();
    // The geometry before step 1000 is no sample.
    assert_eq!(
        distribution.samples(),
        &Samples::Subset(vec![1, 2, 3, 4, 5, 6])
    );
    assert_eq!(distribution.estimator(), CUMULANT);
    assert_eq!(distribution.target(), Some(at(300.0)));
    let beta = 1.0 / (R * 300.0);
    let boost = [0.0, 1.0, 2.0, 0.5, 3.0, 1.5];
    let dihedral = [0.1, 0.0, 0.3, 0.2, 0.0, 0.4];
    for ((reduced, boost), dihedral) in reduced(distribution).iter().zip(boost).zip(dihedral) {
        assert_close(*reduced, beta * (boost + 4.184 * dihedral));
    }
    assert_eq!(distribution.state_count(), 1);
    assert_eq!(ensemble.state(0).unwrap().representative(), 1);
}

#[test]
fn another_temperature_reads_the_potential_energy() {
    let ensemble =
        Ensemble::reweighted(six(), Some(kelvin(310.0)), Estimator::Exponential).unwrap();
    let distribution = ensemble.distribution();
    assert_eq!(distribution.target(), Some(at(310.0)));
    let (sampled, target) = (1.0 / (R * 300.0), 1.0 / (R * 310.0));
    let potential = [-12.0, -9.0, -11.5, -10.5, -8.0, -12.5];
    let bias = [0.0, 1.0, 2.0, 0.5, 3.0, 1.5]
        .iter()
        .zip([0.1, 0.0, 0.3, 0.2, 0.0, 0.4])
        .map(|(boost, dihedral)| boost + 4.184 * dihedral);
    for ((reduced, potential), bias) in reduced(distribution).iter().zip(potential).zip(bias) {
        assert_close(*reduced, sampled * (potential + bias) - target * potential);
    }
    // An isothermal-isobaric run keeps its pressure at the target.
    let mut set = six();
    let pressure = Quantity::new(1.0, kekule::units::BAR);
    let npt = ThermodynamicState::isothermal_isobaric(kelvin(300.0), pressure).unwrap();
    let sampling = set.sampling().unwrap().clone();
    let mut record = Sampling::new(npt).with_from_step(sampling.from_step());
    record = record
        .with_potential_energy(key("potential"))
        .unwrap()
        .with_bias(sampling.bias().iter().cloned())
        .unwrap();
    set.set_sampling(Some(record)).unwrap();
    let ensemble = Ensemble::reweighted(set, Some(kelvin(310.0)), CUMULANT).unwrap();
    let target = ensemble.distribution().target().unwrap();
    assert_eq!(target.pressure(), npt.pressure());
    assert_eq!(*target.temperature().value(), 310.0);
}

#[test]
fn every_estimator_shares_a_state_exactly() -> Result<(), StructureSetError> {
    for estimator in [
        Estimator::Exponential,
        Estimator::Cumulant(CumulantOrder::First),
        CUMULANT,
        Estimator::Cumulant(CumulantOrder::Third),
    ] {
        let mut ensemble = Ensemble::reweighted(six(), None, estimator)?;
        ensemble
            .distribution_mut()
            .regroup(Partition::new(vec![0, 0, 1, 1, 1, 0], vec![0, 2]).unwrap())?;
        let distribution = ensemble.distribution();
        let b = reduced(distribution);
        let first = [b[0], b[1], b[5]];
        let second = [b[2], b[3], b[4]];
        let populations = [
            estimator.log_sum(&first).unwrap(),
            estimator.log_sum(&second).unwrap(),
        ];
        let total = log_sum_exp(&populations);
        let states = distribution.state_probabilities()?;
        assert_close(states[0], (populations[0] - total).exp());
        assert_close(states[1], (populations[1] - total).exp());
        // Exact shares within each state: probabilities sum to the state's.
        let probabilities = distribution.probabilities()?;
        assert_close(
            probabilities[0] + probabilities[1] + probabilities[5],
            states[0],
        );
        assert_close(
            probabilities[0] / probabilities[1],
            (first[0] - first[1]).exp(),
        );
        // A group that is a state, alone or among bins, has its population.
        assert_close(
            distribution.log_population(&[0, 1, 5]).unwrap(),
            states[0].ln(),
        );
        let bins = distribution
            .group_log_populations(&[Some(0), Some(0), Some(1), Some(1), Some(1), Some(0)], 2)?;
        assert_close(bins[0].unwrap(), states[0].ln());
        assert_close(bins[1].unwrap(), states[1].ln());
        // The state view reads the same population.
        assert_close(ensemble.state(1).unwrap().log_population(), populations[1]);
        // Ranked by probability.
        let ranked = distribution.ranked_by_weight();
        let mut expected = (0..6).collect::<Vec<_>>();
        expected.sort_by(|&x, &y| probabilities[y].total_cmp(&probabilities[x]));
        assert_eq!(
            ranked,
            expected
                .iter()
                .map(|position| position + 1)
                .collect::<Vec<_>>()
        );
    }
    Ok(())
}

#[test]
fn the_second_cumulant_sums_a_state() -> Result<(), StructureSetError> {
    let ensemble = Ensemble::reweighted(six(), None, CUMULANT)?;
    let distribution = ensemble.distribution();
    let b = reduced(distribution);
    assert_close(
        ensemble.state(0).unwrap().log_population(),
        second_cumulant(&b),
    );
    // One state: every sample's probability is its exact share.
    let total = log_sum_exp(&b);
    for (probability, reduced) in distribution.probabilities()?.iter().zip(&b) {
        assert_close(*probability, (reduced - total).exp());
    }
    // A group of some samples, relative to the state population.
    assert_close(
        distribution.log_population(&[1, 3]).unwrap(),
        second_cumulant(&[b[1], b[3]]) - second_cumulant(&b),
    );
    Ok(())
}

#[test]
fn a_constant_reduced_bias_per_region_is_exact_under_every_estimator(
) -> Result<(), StructureSetError> {
    // Bias 1 kJ/mol on the first three samples and 4 kJ/mol on the other three.
    let steps = [1_000, 1_500, 2_000, 2_500, 3_000, 3_500];
    let potential = [Some(0.0); 6];
    let boost = [1.0, 1.0, 1.0, 4.0, 4.0, 4.0].map(Some);
    let set = biased_run(&steps, &potential, &boost, &[Some(0.0); 6]);
    let beta = 1.0 / (R * 300.0);
    let expected = 1.0 / (1.0 + (3.0 * beta).exp());
    for estimator in [
        Estimator::Exponential,
        Estimator::Cumulant(CumulantOrder::First),
        CUMULANT,
        Estimator::Cumulant(CumulantOrder::Third),
    ] {
        let distribution = Distribution::reweighted(
            &set,
            Samples::All,
            Partition::new(vec![0, 0, 0, 1, 1, 1], vec![0, 3]).unwrap(),
            None,
            estimator,
        )?;
        assert_close(distribution.state_probabilities()?[0], expected);
        let bins = distribution
            .group_log_populations(&[Some(0), Some(0), Some(0), Some(1), Some(1), Some(1)], 2)?;
        assert_close(bins[0].unwrap().exp(), expected);
        assert_close(
            distribution.log_population(&[3, 4, 5]).unwrap().exp(),
            1.0 - expected,
        );
    }
    Ok(())
}

#[test]
fn temperature_reweighting_recovers_a_two_level_system() -> Result<(), StructureSetError> {
    // At 300 K a level ε = RT ln 3 above the ground holds a quarter of an
    // unbiased run: three frames at 0 and one at ε.
    let epsilon = R * 300.0 * 3.0_f64.ln();
    let steps = [1_000, 1_500, 2_000, 2_500];
    let potential = [0.0, 0.0, 0.0, epsilon].map(Some);
    let set = biased_run(&steps, &potential, &[Some(0.0); 4], &[Some(0.0); 4]);
    for (temperature, estimator) in [(350.0, Estimator::Exponential), (250.0, CUMULANT)] {
        let mut ensemble = Ensemble::reweighted(set.clone(), Some(kelvin(temperature)), estimator)?;
        ensemble
            .distribution_mut()
            .regroup(Partition::new(vec![0, 0, 0, 1], vec![0, 3]).unwrap())?;
        let excited = (-epsilon / (R * temperature)).exp();
        assert_close(
            ensemble.distribution().state_probabilities()?[1],
            excited / (1.0 + excited),
        );
    }
    Ok(())
}

#[test]
fn geometries_that_cannot_be_samples_are_refused() {
    let steps = [500, 1_000, 1_500];
    let set = biased_run(
        &steps,
        &[Some(0.0), None, Some(1.0)],
        &[Some(0.0), Some(1.0), Some(2.0)],
        &[Some(0.0), Some(0.0), None],
    );
    let all = |set: &StructureSet, temperature| {
        Distribution::reweighted(
            set,
            Samples::All,
            Partition::one_state(3),
            temperature,
            CUMULANT,
        )
    };
    let error = all(&set, None).unwrap_err();
    assert_eq!(error, StructureSetError::BeforeSampling { geometry: 0 });
    assert_eq!(
        error.to_string(),
        "Geometry 0 lies before the first sampled step."
    );
    let later = |samples| {
        Distribution::reweighted(
            &set,
            Samples::Subset(samples),
            Partition::one_state(1),
            None,
            CUMULANT,
        )
    };
    let error = later(vec![2]).unwrap_err();
    assert_eq!(
        error,
        StructureSetError::MissingColumnValue {
            key: key("dihedral"),
            geometry: 2
        }
    );
    assert_eq!(
        error.to_string(),
        "Geometry 2 has no value in the column \"dihedral\"."
    );
    // The potential energy is read only at another temperature.
    assert!(later(vec![1]).is_ok());
    let warmer = Distribution::reweighted(
        &set,
        Samples::Subset(vec![1]),
        Partition::one_state(1),
        Some(kelvin(310.0)),
        CUMULANT,
    );
    assert_eq!(
        warmer.unwrap_err(),
        StructureSetError::MissingColumnValue {
            key: key("potential"),
            geometry: 1
        }
    );
    // An ensemble takes the geometries that can be samples.
    let ensemble = Ensemble::reweighted(set.clone(), None, CUMULANT).unwrap();
    assert_eq!(ensemble.distribution().samples(), &Samples::Subset(vec![1]));
    // None can be at 310 K: the reason the last one cannot.
    assert_eq!(
        Ensemble::reweighted(set.clone(), Some(kelvin(310.0)), CUMULANT).unwrap_err(),
        StructureSetError::MissingColumnValue {
            key: key("dihedral"),
            geometry: 2
        }
    );
    // A set without a record, or without the record's columns or steps.
    let mut unrecorded = set.clone();
    unrecorded.set_sampling(None).unwrap();
    assert_eq!(
        all(&unrecorded, None).unwrap_err(),
        StructureSetError::Missing(Requirement::Sampling)
    );
    let mut unordered = set.clone();
    unordered.set_order(None).unwrap();
    assert_eq!(
        all(&unordered, None).unwrap_err(),
        StructureSetError::Missing(Requirement::Steps)
    );
    let mut no_energy = set.clone();
    no_energy
        .set_sampling(Some(
            Sampling::new(at(300.0)).with_bias([key("boost")]).unwrap(),
        ))
        .unwrap();
    assert_eq!(
        all(&no_energy, Some(kelvin(310.0))).unwrap_err(),
        StructureSetError::Missing(Requirement::PotentialEnergy)
    );
    assert_eq!(
        all(&no_energy, Some(kelvin(-1.0))).unwrap_err(),
        StructureSetError::InvalidTemperature
    );
    // An unphysical bias.
    let huge = biased_run(&[1_000], &[Some(0.0)], &[Some(1e300)], &[Some(0.0)]);
    assert_eq!(
        Ensemble::reweighted(huge, None, CUMULANT).unwrap_err(),
        StructureSetError::InvalidReducedBias { geometry: 0 }
    );
}

#[test]
fn effective_sample_sizes_use_the_exact_weights() -> Result<(), StructureSetError> {
    let ensemble = Ensemble::reweighted(six(), None, CUMULANT)?;
    let distribution = ensemble.distribution();
    let b = reduced(distribution);
    let kish = |values: &[f64]| {
        let weights = values.iter().map(|value| value.exp()).collect::<Vec<_>>();
        weights.iter().sum::<f64>().powi(2) / weights.iter().map(|w| w * w).sum::<f64>()
    };
    assert_close(distribution.effective_sample_size()?, kish(&b));
    assert_close(
        distribution.group_effective_sample_size(&[0, 4]).unwrap(),
        kish(&[b[0], b[4]]),
    );
    Ok(())
}

#[test]
fn the_estimator_changes_and_the_target_takes_a_reweighting() -> Result<(), StructureSetError> {
    let mut ensemble = Ensemble::reweighted(six(), None, CUMULANT)?;
    let b = reduced(ensemble.distribution());
    ensemble
        .distribution_mut()
        .set_estimator(Estimator::Exponential)?;
    assert_eq!(ensemble.distribution().estimator(), Estimator::Exponential);
    assert_eq!(reduced(ensemble.distribution()), b);
    let error = ensemble
        .distribution_mut()
        .set_target(Some(at(310.0)))
        .unwrap_err();
    assert_eq!(error, StructureSetError::ReweightedTarget);
    assert_eq!(
        error.to_string(),
        "The distribution is reweighted to its target. Reweight it to change the target."
    );
    assert_eq!(ensemble.distribution().target(), Some(at(300.0)));
    // Reweighting again keeps the samples, the states, and their data.
    ensemble
        .distribution_mut()
        .regroup(Partition::new(vec![0, 1, 0, 1, 0, 1], vec![0, 1]).unwrap())?;
    ensemble
        .distribution_mut()
        .state_mut(1)
        .unwrap()
        .set_energy(Some(Quantity::new(1.0, KILOJOULE_PER_MOLE)))?;
    let partition = ensemble.distribution().partition().clone();
    ensemble.reweight(Some(kelvin(310.0)), CUMULANT)?;
    let distribution = ensemble.distribution();
    assert_eq!(distribution.target(), Some(at(310.0)));
    assert_eq!(distribution.partition(), &partition);
    assert_eq!(
        distribution.samples(),
        &Samples::Subset(vec![1, 2, 3, 4, 5, 6])
    );
    assert!(distribution.states()[1].energy().is_some());
    assert_ne!(reduced(distribution), b);
    // A failed reweighting changes nothing.
    let before = ensemble.distribution().clone();
    assert!(ensemble.reweight(Some(kelvin(0.0)), CUMULANT).is_err());
    assert_eq!(ensemble.distribution(), &before);
    // Equal weights reweight too.
    let (set, _) = ensemble.into_parts();
    let mut uniform = Ensemble::uniform(set);
    assert!(matches!(
        uniform.reweight(None, CUMULANT).unwrap_err(),
        StructureSetError::BeforeSampling { geometry: 0 }
    ));
    Ok(())
}

#[test]
fn edits_of_log_weights_keep_the_probabilities() -> Result<(), StructureSetError> {
    let mut ensemble = Ensemble::reweighted(six(), None, CUMULANT)?;
    ensemble
        .distribution_mut()
        .regroup(Partition::new(vec![0, 1, 0, 1, 0, 1], vec![0, 1]).unwrap())?;
    let probabilities = ensemble.distribution().probabilities()?;
    let mut normalized = ensemble.clone();
    normalized.distribution_mut().normalize_log_weights()?;
    let WeightModel::LogWeights(log_weights) = normalized.distribution().weight_model() else {
        panic!("normalized weights are log weights");
    };
    for (log_weight, probability) in log_weights.iter().zip(&probabilities) {
        assert_close(log_weight.exp(), *probability);
    }
    ensemble.distribution_mut().set_log_weight(1, 0.0)?;
    assert!(matches!(
        ensemble.distribution().weight_model(),
        WeightModel::LogWeights(_)
    ));
    assert_eq!(ensemble.distribution().estimator(), Estimator::Exponential);
    assert_close(
        ensemble.distribution().log_weight(2).unwrap()
            - ensemble.distribution().log_weight(4).unwrap(),
        (probabilities[1] / probabilities[3]).ln(),
    );
    Ok(())
}

#[test]
fn selecting_and_collapsing_keep_the_reweighting() -> Result<(), StructureSetError> {
    let mut ensemble = Ensemble::reweighted(six(), None, CUMULANT)?;
    ensemble
        .distribution_mut()
        .regroup(Partition::new(vec![0, 0, 1, 1, 1, 0], vec![0, 2]).unwrap())?;
    let b = reduced(ensemble.distribution());
    let selected = ensemble.select([6, 0, 3])?;
    let distribution = selected.distribution();
    assert_eq!(distribution.samples(), &Samples::Subset(vec![0, 2]));
    assert_eq!(reduced(distribution), [b[5], b[2]]);
    assert_eq!(distribution.estimator(), CUMULANT);
    assert_eq!(distribution.target(), Some(at(300.0)));
    let states = ensemble.select_states([1])?;
    assert_eq!(reduced(states.distribution()), [b[2], b[3], b[4]]);

    let collapsed = ensemble.collapse();
    let shares = ensemble.distribution().state_probabilities()?;
    let WeightModel::LogWeights(log_weights) = collapsed.distribution().weight_model() else {
        panic!("collapsed weights are log weights");
    };
    assert_close(
        log_weights[0] - log_weights[1],
        (shares[0] / shares[1]).ln(),
    );
    assert_eq!(collapsed.distribution().target(), Some(at(300.0)));
    assert!(collapsed.set().sampling().is_some());
    Ok(())
}

#[test]
fn a_new_sample_cannot_join_reweighted_samples() -> Result<(), StructureSetError> {
    let mut ensemble = Ensemble::reweighted(six(), None, CUMULANT)?;
    let geometry = ensemble.set().get(1).unwrap().payload().clone();
    let place = Some(kekule::structure::Place::Step(4_000));
    let error = ensemble
        .push_sample(geometry.clone(), place.clone(), 0, Some(0.0))
        .unwrap_err();
    assert_eq!(error, StructureSetError::ReweightedSample { geometry: 7 });
    assert_eq!(
        error.to_string(),
        "Geometry 7 cannot join the reweighted samples. Add it outside the distribution."
    );
    assert!(ensemble
        .push_state(geometry.clone(), place.clone(), StateData::new(), None)
        .is_err());
    assert_eq!(ensemble.set().len(), 7);
    assert_eq!(ensemble.push(geometry, place)?, 7);
    assert_eq!(
        ensemble.distribution().samples(),
        &Samples::Subset(vec![1, 2, 3, 4, 5, 6])
    );
    Ok(())
}

#[test]
fn reweighted_ensembles_round_trip_natively() -> Result<(), Box<dyn std::error::Error>> {
    for estimator in [
        Estimator::Exponential,
        CUMULANT,
        Estimator::Cumulant(CumulantOrder::Third),
    ] {
        let ensemble = Ensemble::reweighted(six(), Some(kelvin(310.0)), estimator)?;
        let bytes = native::write(&ensemble, NativeWriteOptions::default())?;
        let read = native::read(bytes.as_slice(), NativeReadOptions::default())?
            .into_ensemble()
            .unwrap();
        assert_eq!(read.distribution(), ensemble.distribution());
        assert_eq!(read.set().sampling(), ensemble.set().sampling());
    }
    Ok(())
}
