//! Weights of a structure set's geometries: equal, given, Boltzmann, and of a
//! biased run, with their probabilities, populations, rankings, effective
//! sample sizes, and the geometries a biased run leaves out.

use kekule::clustering::Clustering;
use kekule::diagnostics::{Countable, IndexOutOfRange};
use kekule::math::Point3;
use kekule::properties::{PropertyColumn, PropertyKey};
use kekule::statistics::{LeftOut, Unweighable, WeightSource, Weights, WeightsError};
use kekule::structure::{
    Geometry, Order, Positions, Requirement, Sampling, StructureSet, ThermodynamicState,
};
use kekule::units::{
    Quantity, CANONICAL_ENERGY_UNIT, KELVIN, KILOCALORIE_PER_MOLE, KILOJOULE_PER_MOLE, NANOMETER,
};

/// The molar gas constant in kJ/(mol K), exact since the 2019 SI.
const R: f64 = 8.314_462_618_153_24e-3;

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

/// Given weights of these linear weights.
fn linear(weights: &[f64]) -> Weights {
    Weights::given(
        weights.iter().map(|weight| Some(weight.ln())).collect(),
        None,
    )
    .unwrap()
}

/// `count` one-atom geometries along x.
fn set_of(count: usize) -> StructureSet {
    let topology = kekule::smiles::to_topology("O").unwrap();
    let geometries = (0..count).map(|index| {
        let point = Point3::new(0.1 * index as f64, 0.0, 0.0);
        Geometry::new(Positions::new(Quantity::new(vec![point], NANOMETER)).unwrap())
    });
    StructureSet::from_geometries(topology, geometries).unwrap()
}

/// A set with the energy column `energy`, in kJ/mol.
fn with_energies(values: &[Option<f64>]) -> StructureSet {
    let mut set = set_of(values.len());
    set.columns_mut()
        .insert(key("energy"), energies(values))
        .unwrap();
    set
}

#[test]
fn equal_weights_weigh_every_listed_geometry_the_same() {
    let weights = Weights::equal(4);
    assert_eq!(weights.len(), 4);
    assert_eq!(weights.sample_count(), 4);
    assert_eq!(weights.source(), &WeightSource::Equal);
    assert_eq!(weights.temperature(), None);
    assert_eq!(weights.probabilities().unwrap(), [0.25; 4]);
    // Kish's effective sample size is the sample count, exactly.
    assert_eq!(weights.effective_sample_size().unwrap(), 4.0);
    assert_eq!(weights.ranked(), [0, 1, 2, 3]);

    // A geometry listed twice is one sample; the others are none.
    let some = Weights::equal_over(5, [3, 1, 1]).unwrap();
    assert_eq!(some.samples().collect::<Vec<_>>(), [1, 3]);
    assert_eq!(some.log_weight(0), None);
    assert_eq!(some.probabilities().unwrap(), [0.0, 0.5, 0.0, 0.5, 0.0]);
    assert_eq!(some.effective_sample_size().unwrap(), 2.0);
    assert_eq!(
        Weights::equal_over(5, [5]),
        Err(WeightsError::IndexOutOfRange(IndexOutOfRange::new(
            Countable::Geometry,
            5,
            5
        )))
    );

    // Without samples there are no probabilities.
    for none in [Weights::equal(0), Weights::equal_over(3, []).unwrap()] {
        assert_eq!(none.probabilities(), Err(WeightsError::NoSamples));
        assert_eq!(none.effective_sample_size(), Err(WeightsError::NoSamples));
        assert_eq!(none.log_population(&[0]), None);
    }
    assert!(Weights::equal(0).is_empty());
    assert_eq!(
        WeightsError::NoSamples.to_string(),
        "The weights have no samples."
    );
}

#[test]
fn given_weights_are_finite_log_weights() {
    // Geometries 0, 2, and 3 weigh 1, 3, and 4; geometry 1 is no sample.
    let weights = Weights::given(
        vec![Some(0.0), None, Some(3.0_f64.ln()), Some(4.0_f64.ln())],
        Some(kelvin(300.0)),
    )
    .unwrap();
    assert_eq!(
        weights.source(),
        &WeightSource::Given {
            temperature: Some(kelvin(300.0))
        }
    );
    assert_eq!(weights.temperature(), Some(kelvin(300.0)));
    let probabilities = weights.probabilities().unwrap();
    for (probability, expected) in probabilities.iter().zip([0.125, 0.0, 0.375, 0.5]) {
        assert_close(*probability, expected);
    }
    // A group's population is the exact sum of its samples' weights; a
    // geometry listed twice counts twice, and one that is no sample adds
    // nothing.
    assert_close(weights.log_population(&[0, 2]).unwrap(), 0.5_f64.ln());
    assert_close(weights.log_population(&[2, 2, 1]).unwrap(), 0.75_f64.ln());
    assert_close(weights.log_population(&[0, 2, 3]).unwrap(), 0.0);
    assert_eq!(weights.log_population(&[1]), None);
    assert_eq!(weights.log_population(&[0, 4]), None);
    assert_eq!(weights.log_weight(9), None);

    // Weights of e^1000, beyond `f64::MAX`, keep their shares.
    let huge = Weights::given(vec![Some(1000.0), Some(1000.0 + 3.0_f64.ln())], None).unwrap();
    assert_close(huge.probabilities().unwrap()[1], 0.75);
    assert_close(huge.log_population(&[0]).unwrap(), 0.25_f64.ln());
    assert_close(huge.effective_sample_size().unwrap(), 16.0 / 10.0);

    for (log_weights, temperature, error, message) in [
        (
            vec![Some(0.0), Some(f64::NAN)],
            None,
            WeightsError::NonFiniteLogWeight { geometry: 1 },
            "The log weight of geometry 1 is not finite.",
        ),
        (
            vec![Some(f64::INFINITY)],
            None,
            WeightsError::NonFiniteLogWeight { geometry: 0 },
            "The log weight of geometry 0 is not finite.",
        ),
        (
            vec![Some(0.0)],
            Some(kelvin(-1.0)),
            WeightsError::InvalidTemperature,
            "The temperature is not finite and positive.",
        ),
    ] {
        let actual = Weights::given(log_weights, temperature).unwrap_err();
        assert_eq!(actual, error);
        assert_eq!(actual.to_string(), message);
    }
    assert!(matches!(
        Weights::given(
            vec![Some(0.0)],
            Some(Quantity::new(1.0, KILOJOULE_PER_MOLE))
        ),
        Err(WeightsError::Unit(_))
    ));
}

#[test]
fn effective_sample_sizes_are_kish() {
    // (1 + 2 + 3 + 4)^2 / (1 + 4 + 9 + 16).
    let weights = linear(&[1.0, 2.0, 3.0, 4.0]);
    assert_close(weights.effective_sample_size().unwrap(), 100.0 / 30.0);
    // One dominant sample leaves one effective sample.
    let dominated = Weights::given(vec![Some(0.0), Some(0.0), Some(800.0)], None).unwrap();
    assert_close(dominated.effective_sample_size().unwrap(), 1.0);
    // A large common offset changes nothing.
    let offset = Weights::given(
        [1.0_f64, 2.0, 3.0, 4.0]
            .iter()
            .map(|weight| Some(weight.ln() + 1e5))
            .collect(),
        None,
    )
    .unwrap();
    assert!((offset.effective_sample_size().unwrap() - 100.0 / 30.0).abs() < 1e-9);
}

#[test]
fn samples_rank_and_cover_by_weight() {
    // Weights 1, 3, 2, and 3 of 9: equal weights keep geometry order.
    let weights = linear(&[1.0, 3.0, 2.0, 3.0]);
    assert_eq!(weights.ranked(), [1, 3, 2, 0]);
    assert_eq!(weights.covering(0.3).unwrap(), [1]);
    assert_eq!(weights.covering(0.5).unwrap(), [1, 3]);
    assert_eq!(weights.covering(0.7).unwrap(), [1, 3, 2]);
    assert_eq!(weights.covering(1.0).unwrap(), [1, 3, 2, 0]);
    // Every sample, even one too light to change a floating-point total.
    let light = Weights::given(vec![Some(0.0), Some(-800.0)], None).unwrap();
    assert_eq!(light.covering(1.0).unwrap(), [0, 1]);
    for fraction in [0.0, -0.5, 1.5, f64::NAN] {
        assert_eq!(
            weights.covering(fraction),
            Err(WeightsError::InvalidFraction)
        );
    }
    assert_eq!(
        WeightsError::InvalidFraction.to_string(),
        "The population fraction is not greater than 0 and at most 1."
    );
    assert_eq!(
        Weights::equal_over(2, []).unwrap().covering(0.5),
        Err(WeightsError::NoSamples)
    );
    // Samples that are not geometries of the weights are not ranked.
    let some = weights.restricted([0, 2]).unwrap();
    assert_eq!(some.ranked(), [2, 0]);
}

#[test]
fn restricting_selecting_and_extending_keep_each_log_weight() {
    let weights = Weights::given(
        vec![Some(0.0), Some(1.0), None, Some(2.0)],
        Some(kelvin(310.0)),
    )
    .unwrap();
    // Restricted to some geometries, the others are no samples.
    let restricted = weights.restricted([3, 0, 2]).unwrap();
    assert_eq!(restricted.len(), 4);
    assert_eq!(restricted.samples().collect::<Vec<_>>(), [0, 3]);
    assert_eq!(restricted.log_weight(3), Some(2.0));
    assert_eq!(restricted.source(), weights.source());
    assert_eq!(
        weights.restricted([4]),
        Err(WeightsError::IndexOutOfRange(IndexOutOfRange::new(
            Countable::Geometry,
            4,
            4
        )))
    );
    // Selected, as a selected set's geometries: a geometry listed twice is
    // two samples.
    let selected = weights.select([3, 2, 3, 1]).unwrap();
    assert_eq!(selected.len(), 4);
    assert_eq!(
        (0..4).map(|g| selected.log_weight(g)).collect::<Vec<_>>(),
        [Some(2.0), None, Some(2.0), Some(1.0)]
    );
    assert_eq!(selected.temperature(), Some(kelvin(310.0)));
    let total = 2.0 * 2.0_f64.exp() + 1.0_f64.exp();
    assert_close(selected.probabilities().unwrap()[0], 2.0_f64.exp() / total);
    assert!(weights.select([7]).is_err());
    // Extended, the new geometries are no samples.
    let extended = weights.clone().extended_to(6).unwrap();
    assert_eq!(extended.len(), 6);
    assert_eq!(extended.sample_count(), 3);
    assert_eq!(extended.log_weight(5), None);
    let error = weights.extended_to(2).unwrap_err();
    assert_eq!(
        error,
        WeightsError::CoverageMismatch {
            covered: 4,
            geometries: 2
        }
    );
    assert_eq!(
        error.to_string(),
        "The weights cover 4 geometries, but the structure set has 2."
    );
}

/// Checks Boltzmann weights: each geometry's probability is
/// `exp(-E / RT)` over the sum of those factors.
fn check_boltzmann(weights: &Weights, kelvin_value: f64, energies: &[Option<f64>]) {
    let thermal = R * kelvin_value;
    let minimum = energies
        .iter()
        .flatten()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let factors = energies
        .iter()
        .map(|energy| energy.map_or(0.0, |energy| (-(energy - minimum) / thermal).exp()))
        .collect::<Vec<_>>();
    let total = factors.iter().sum::<f64>();
    for (probability, factor) in weights.probabilities().unwrap().iter().zip(&factors) {
        assert_close(*probability, factor / total);
    }
    // Relative to the lowest energy, which has log weight zero.
    for (geometry, energy) in energies.iter().enumerate() {
        assert_eq!(
            weights.log_weight(geometry).is_some(),
            energy.is_some(),
            "geometry {geometry}"
        );
        if let Some(energy) = energy {
            assert_close(
                weights.log_weight(geometry).unwrap(),
                -(energy - minimum) / thermal,
            );
        }
    }
    assert_eq!(weights.temperature(), Some(kelvin(kelvin_value)));
}

#[test]
fn boltzmann_weights_follow_the_energies() {
    let values = [1.5, 1.5, 0.0, 0.0, 4.0, 4.0, 4.0, 0.7].map(Some);
    let set = with_energies(&values);
    for temperature in [1.0, 298.15, 5000.0] {
        let weights = Weights::boltzmann(&set, &key("energy"), kelvin(temperature)).unwrap();
        check_boltzmann(&weights, temperature, &values);
        assert_eq!(
            weights.source(),
            &WeightSource::Boltzmann {
                temperature: kelvin(temperature),
                energy: key("energy")
            }
        );
    }
    // At room temperature the lowest energies weigh most, equal energies in
    // geometry order.
    let weights = Weights::boltzmann(&set, &key("energy"), kelvin(298.15)).unwrap();
    assert_eq!(weights.ranked(), [2, 3, 7, 0, 1, 4, 5, 6]);

    // A cluster's population is the exact sum of its geometries' factors,
    // so a cluster's free energy is -RT ln of that sum.
    let clusters = Clustering::new(
        [0, 0, 1, 1, 2, 2, 2, 3].map(Some).to_vec(),
        vec![0, 2, 4, 7],
    )
    .unwrap();
    let thermal = R * 298.15;
    let clustered: [&[f64]; 4] = [&[1.5, 1.5], &[0.0, 0.0], &[4.0, 4.0, 4.0], &[0.7]];
    let sums = clustered.map(|energies| {
        energies
            .iter()
            .map(|energy| (-energy / thermal).exp())
            .sum::<f64>()
    });
    let total = sums.iter().sum::<f64>();
    for (population, sum) in clusters.log_populations(&weights).unwrap().iter().zip(sums) {
        assert_close(population.unwrap().exp(), sum / total);
    }
    let free = clusters.free_energies(&weights).unwrap();
    assert_eq!(free.unit(), CANONICAL_ENERGY_UNIT);
    for (energy, sum) in free.value().iter().zip(sums) {
        assert!((energy.unwrap() - thermal * (sums[1] / sum).ln()).abs() < 1e-12);
    }

    // A geometry without an energy is no sample.
    let partial = with_energies(&[Some(0.0), None, Some(2.0)]);
    let weights = Weights::boltzmann(&partial, &key("energy"), kelvin(300.0)).unwrap();
    check_boltzmann(&weights, 300.0, &[Some(0.0), None, Some(2.0)]);
}

#[test]
fn a_remote_geometry_keeps_a_finite_log_weight() {
    // A geometry far above the others weighs less by its energy over RT,
    // though its probability underflows.
    let values = [0.0, 1.0e6, 2.0, 3.0].map(Some);
    let set = with_energies(&values);
    let weights = Weights::boltzmann(&set, &key("energy"), kelvin(1.0)).unwrap();
    let gap = weights.log_weight(0).unwrap() - weights.log_weight(1).unwrap();
    let expected = 1.0e6 / R;
    assert!((gap - expected).abs() <= 1.0e-9 * expected);
    assert_eq!(weights.probabilities().unwrap()[1], 0.0);
    assert!(weights.log_population(&[1]).unwrap().is_finite());
    // An energy so far above that its log weight overflows is refused.
    let beyond = with_energies(&[Some(0.0), Some(1.0e307)]);
    assert_eq!(
        Weights::boltzmann(&beyond, &key("energy"), kelvin(1.0)),
        Err(WeightsError::NonFiniteLogWeight { geometry: 1 })
    );
}

#[test]
fn boltzmann_weights_need_an_energy_column_and_a_temperature() {
    let mut set = with_energies(&[Some(0.0), Some(1.0)]);
    set.columns_mut()
        .insert(key("count"), PropertyColumn::Int(vec![Some(1), Some(2)]))
        .unwrap();
    set.columns_mut()
        .insert(
            key("length"),
            PropertyColumn::Real {
                unit: NANOMETER,
                values: vec![Some(1.0), Some(2.0)],
            },
        )
        .unwrap();
    for (column, error, message) in [
        (
            "missing",
            WeightsError::MissingColumn {
                key: key("missing"),
            },
            "The structure set has no column \"missing\".",
        ),
        (
            "count",
            WeightsError::NonEnergyColumn { key: key("count") },
            "The column \"count\" holds no energies.",
        ),
        (
            "length",
            WeightsError::NonEnergyColumn { key: key("length") },
            "The column \"length\" holds no energies.",
        ),
    ] {
        let actual = Weights::boltzmann(&set, &key(column), kelvin(300.0)).unwrap_err();
        assert_eq!(actual, error);
        assert_eq!(actual.to_string(), message);
    }
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            Weights::boltzmann(&set, &key("energy"), kelvin(temperature)),
            Err(WeightsError::InvalidTemperature)
        );
    }
    // Energies in kcal/mol are converted.
    set.columns_mut()
        .insert(
            key("kcal"),
            PropertyColumn::Real {
                unit: KILOCALORIE_PER_MOLE,
                values: vec![Some(0.0), Some(1.0)],
            },
        )
        .unwrap();
    let weights = Weights::boltzmann(&set, &key("kcal"), kelvin(300.0)).unwrap();
    assert_close(weights.log_weight(1).unwrap(), -4.184 / (R * 300.0));
}

type Bias<'a> = (&'a [Option<f64>], &'a [Option<f64>]);

/// A run at 300 K whose frames are samples from step `from_step`, with
/// geometries at `steps`, a potential energy, and the bias columns `a`, in
/// kJ/mol, and `b`, in kcal/mol, or no bias for `None`.
fn run(
    from_step: u64,
    steps: &[u64],
    potential: &[Option<f64>],
    bias: Option<Bias<'_>>,
) -> StructureSet {
    let mut set = set_of(steps.len());
    set.set_order(Some(Order::Step(steps.to_vec()))).unwrap();
    let mut columns = set.columns_mut();
    columns
        .insert(key("potential"), energies(potential))
        .unwrap();
    let mut sampling = Sampling::new(at(300.0))
        .with_from_step(from_step)
        .with_potential_energy(key("potential"))
        .unwrap();
    if let Some((a, b)) = bias {
        columns.insert(key("a"), energies(a)).unwrap();
        columns
            .insert(
                key("b"),
                PropertyColumn::Real {
                    unit: KILOCALORIE_PER_MOLE,
                    values: b.to_vec(),
                },
            )
            .unwrap();
        sampling = sampling.with_bias([key("a"), key("b")]).unwrap();
    }
    set.set_sampling(Some(sampling)).unwrap();
    set
}

fn reasons(left_out: &[LeftOut]) -> Vec<(usize, Unweighable)> {
    left_out
        .iter()
        .map(|left_out| (left_out.geometry(), left_out.reason().clone()))
        .collect()
}

#[test]
fn bias_weights_sum_the_bias_columns_at_the_sampled_temperature() {
    let steps = [500, 1_000, 1_500, 2_000, 2_500, 3_000, 3_500];
    let potential = [-10.0, -12.0, -9.0, -11.5, -10.5, -8.0, -12.5].map(Some);
    let a = [9.0, 0.0, 1.0, 2.0, 0.5, 3.0, 1.5].map(Some);
    let b = [9.0, 0.1, 0.0, 0.3, 0.2, 0.0, 0.4].map(Some);
    let set = run(1_000, &steps, &potential, Some((&a, &b)));
    let bias = Weights::from_bias(&set).unwrap();
    let weights = bias.weights();
    assert_eq!(weights.source(), &WeightSource::Bias { state: at(300.0) });
    assert_eq!(weights.temperature(), Some(kelvin(300.0)));
    // The geometry before step 1000 is no sample.
    assert_eq!(weights.samples().collect::<Vec<_>>(), [1, 2, 3, 4, 5, 6]);
    assert_eq!(reasons(bias.left_out()), [(0, Unweighable::BeforeSampling)]);
    // ln w = (ΔV_a + ΔV_b) / RT, with ΔV_b converted from kcal/mol.
    let beta = 1.0 / (R * 300.0);
    let mut reduced = Vec::new();
    for geometry in 1..7 {
        let expected = beta * (a[geometry].unwrap() + 4.184 * b[geometry].unwrap());
        assert_close(weights.log_weight(geometry).unwrap(), expected);
        reduced.push(expected);
    }
    // Kish's effective sample size of the exact weights.
    let linear = reduced.iter().map(|value| value.exp()).collect::<Vec<_>>();
    let kish = linear.iter().sum::<f64>().powi(2) / linear.iter().map(|w| w * w).sum::<f64>();
    assert_close(weights.effective_sample_size().unwrap(), kish);
    assert_eq!(bias.clone().into_weights(), *weights);

    // Without bias columns every sample weighs the same.
    let unbiased = run(1_000, &steps, &potential, None);
    let equal = Weights::from_bias(&unbiased).unwrap().into_weights();
    assert_eq!(equal.samples().collect::<Vec<_>>(), [1, 2, 3, 4, 5, 6]);
    assert!(equal.samples().all(|g| equal.log_weight(g) == Some(0.0)));
}

#[test]
fn geometries_a_biased_run_cannot_weigh_are_left_out_with_their_reason() {
    let set = run(
        1_000,
        &[500, 1_000, 1_500, 2_000, 2_500],
        &[Some(-1.0), Some(-2.0), None, Some(-3.0), Some(-4.0)],
        Some((
            &[Some(1.0), None, Some(2.0), Some(3.0), None],
            &[Some(0.5), None, Some(0.5), None, Some(0.5)],
        )),
    );
    let bias = Weights::from_bias(&set).unwrap();
    // The potential energy is not read: geometry 2 lacks one and is weighed.
    assert_eq!(bias.weights().samples().collect::<Vec<_>>(), [2]);
    // The first bias column without a value names the reason.
    assert_eq!(
        reasons(bias.left_out()),
        [
            (0, Unweighable::BeforeSampling),
            (1, Unweighable::MissingBias { key: key("a") }),
            (3, Unweighable::MissingBias { key: key("b") }),
            (4, Unweighable::MissingBias { key: key("a") }),
        ]
    );
}

#[test]
fn a_run_without_a_weighable_geometry_says_why() {
    // Every geometry lies before the first sampled step.
    let early = run(
        1_000,
        &[0, 500],
        &[Some(0.0); 2],
        Some((&[Some(1.0); 2], &[Some(0.0); 2])),
    );
    let error = Weights::from_bias(&early).unwrap_err();
    assert_eq!(
        error,
        WeightsError::NoWeighableSample(Unweighable::BeforeSampling)
    );
    assert_eq!(
        error.to_string(),
        "No geometry can be weighed. Every geometry lies before the first sampled step."
    );
    // The furthest reason wins: one before, one without a value in "b".
    let missing = run(
        1_000,
        &[0, 1_000],
        &[Some(0.0); 2],
        Some((&[Some(1.0); 2], &[Some(0.0), None])),
    );
    let error = Weights::from_bias(&missing).unwrap_err();
    assert_eq!(
        error,
        WeightsError::NoWeighableSample(Unweighable::MissingBias { key: key("b") })
    );
    assert_eq!(
        error.to_string(),
        "No geometry can be weighed. Geometries lack values in the bias column \"b\"."
    );
    // A set without geometries has no weights and nothing left out.
    let empty = run(0, &[], &[], None);
    let bias = Weights::from_bias(&empty).unwrap();
    assert!(bias.weights().is_empty() && bias.left_out().is_empty());
}

#[test]
fn a_biased_run_needs_its_record_steps_and_physical_biases() {
    let set = run(
        1_000,
        &[1_000, 1_500],
        &[Some(0.0); 2],
        Some((&[Some(1.0); 2], &[Some(0.0); 2])),
    );
    let mut unrecorded = set.clone();
    unrecorded.set_sampling(None).unwrap();
    let error = Weights::from_bias(&unrecorded).unwrap_err();
    assert_eq!(error, WeightsError::Missing(Requirement::Sampling));
    assert_eq!(
        error.to_string(),
        "The structure set needs a sampling record."
    );
    // A record from a later step needs steps to compare with.
    let mut unordered = set.clone();
    unordered.set_order(None).unwrap();
    let error = Weights::from_bias(&unordered).unwrap_err();
    assert_eq!(error, WeightsError::Missing(Requirement::Steps));
    assert_eq!(error.to_string(), "The structure set needs steps.");
    // From step 0 every geometry is sampled, with or without steps.
    let mut from_zero = unordered.clone();
    from_zero
        .set_sampling(Some(
            Sampling::new(at(300.0)).with_bias([key("a")]).unwrap(),
        ))
        .unwrap();
    assert_eq!(
        Weights::from_bias(&from_zero)
            .unwrap()
            .weights()
            .sample_count(),
        2
    );
    // An unphysical bias.
    let huge = run(
        0,
        &[0],
        &[Some(0.0)],
        Some((&[Some(1e308)], &[Some(1e308)])),
    );
    let error = Weights::from_bias(&huge).unwrap_err();
    assert_eq!(error, WeightsError::InvalidBias { geometry: 0 });
    assert_eq!(error.to_string(), "The bias of geometry 0 is out of range.");
}

/// A small deterministic generator for property-style cases.
struct Lcg(u64);

impl Lcg {
    fn next(&mut self, bound: u64) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) % bound
    }

    fn cell(&mut self) -> Option<f64> {
        (self.next(4) != 0).then(|| self.next(100) as f64 / 10.0)
    }
}

/// Over random runs, a geometry is a sample exactly when it lies at or after
/// the first sampled step and has every bias value, and every other
/// geometry is left out with the first reason that applies.
#[test]
fn bias_weights_agree_with_their_rule_on_random_runs() {
    let mut random = Lcg(7);
    for case in 0..200 {
        let len = 1 + random.next(8) as usize;
        let steps = (0..len as u64).map(|step| step * 10).collect::<Vec<_>>();
        let potential = (0..len).map(|_| random.cell()).collect::<Vec<_>>();
        let a = (0..len).map(|_| random.cell()).collect::<Vec<_>>();
        let b = (0..len).map(|_| random.cell()).collect::<Vec<_>>();
        let biased = random.next(4) != 0;
        let from_step = random.next(3) * 10;
        let set = run(
            from_step,
            &steps,
            &potential,
            biased.then_some((&a[..], &b[..])),
        );
        let expected = (0..len)
            .map(|geometry| {
                if steps[geometry] < from_step {
                    Some(Unweighable::BeforeSampling)
                } else if biased && a[geometry].is_none() {
                    Some(Unweighable::MissingBias { key: key("a") })
                } else if biased && b[geometry].is_none() {
                    Some(Unweighable::MissingBias { key: key("b") })
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        match Weights::from_bias(&set) {
            Ok(bias) => {
                let weights = bias.weights();
                for (geometry, reason) in expected.iter().enumerate() {
                    assert_eq!(
                        weights.log_weight(geometry).is_none(),
                        reason.is_some(),
                        "case {case}"
                    );
                }
                let left_out = expected
                    .iter()
                    .enumerate()
                    .filter_map(|(geometry, reason)| reason.clone().map(|r| (geometry, r)))
                    .collect::<Vec<_>>();
                assert_eq!(reasons(bias.left_out()), left_out, "case {case}");
                assert_eq!(weights.sample_count() + left_out.len(), len, "case {case}");
            }
            Err(WeightsError::NoWeighableSample(reason)) => {
                assert!(expected.iter().all(Option::is_some), "case {case}");
                // The reason of the first geometry that got furthest.
                let furthest = expected
                    .iter()
                    .flatten()
                    .find(|reason| matches!(reason, Unweighable::MissingBias { .. }))
                    .unwrap_or(&Unweighable::BeforeSampling);
                assert_eq!(&reason, furthest, "case {case}");
            }
            Err(error) => panic!("case {case}: {error:?}"),
        }
    }
}
