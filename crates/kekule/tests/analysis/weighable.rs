//! Which samples a reweighting can weigh: the query, the reasons, and the
//! reweighting that leaves the others out and renumbers the states.

use kekule::properties::{PropertyColumn, PropertyKey};
use kekule::statistics::{CumulantOrder, Estimator};
use kekule::structure::{
    Distribution, Ensemble, Geometry, Order, Partition, Positions, Requirement, Samples, Sampling,
    StateData, StructureSet, StructureSetError, ThermodynamicState, Unweighable, WeightModel,
};
use kekule::units::{Quantity, KELVIN, KILOJOULE_PER_MOLE};

/// Two bias columns' cells.
type Bias<'a> = (&'a [Option<f64>], &'a [Option<f64>]);

const CUMULANT: Estimator = Estimator::Cumulant(CumulantOrder::Second);

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

fn kelvin(value: f64) -> Quantity<f64> {
    Quantity::new(value, KELVIN)
}

fn energies(values: &[Option<f64>]) -> PropertyColumn {
    PropertyColumn::Real {
        unit: KILOJOULE_PER_MOLE,
        values: values.to_vec(),
    }
}

/// A run at 300 K whose frames are samples from step `from_step`, with
/// geometries at `steps`, a potential energy, and the bias columns `a` and
/// `b`, or no bias for `None`.
fn run(
    from_step: u64,
    steps: &[u64],
    potential: &[Option<f64>],
    bias: Option<Bias<'_>>,
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
    let mut sampling = Sampling::new(ThermodynamicState::canonical(kelvin(300.0)).unwrap())
        .with_from_step(from_step)
        .with_potential_energy(key("potential"))
        .unwrap();
    if let Some((a, b)) = bias {
        columns.insert(key("a"), energies(a)).unwrap();
        columns.insert(key("b"), energies(b)).unwrap();
        sampling = sampling.with_bias([key("a"), key("b")]).unwrap();
    }
    set.set_sampling(Some(sampling)).unwrap();
    set
}

fn reduced(distribution: &Distribution) -> Vec<f64> {
    match distribution.weight_model() {
        WeightModel::Reweighted { reduced, .. } => reduced.to_vec(),
        other => panic!("{other:?} is not reweighted"),
    }
}

fn reasons(weighable: &[kekule::structure::LeftOut]) -> Vec<(usize, Unweighable)> {
    weighable
        .iter()
        .map(|left_out| (left_out.geometry(), left_out.reason().clone()))
        .collect()
}

#[test]
fn every_sample_left_out_is_an_error_that_says_why() {
    // Every geometry lies before the first sampled step.
    let early = run(
        1_000,
        &[0, 500],
        &[Some(0.0), Some(0.0)],
        Some((&[Some(1.0), Some(1.0)], &[Some(0.0), Some(0.0)])),
    );
    let before = StructureSetError::NoWeighableSample(Unweighable::BeforeSampling);
    assert_eq!(
        Distribution::weighable(&early, &Samples::All, None).unwrap_err(),
        before
    );
    let error = Ensemble::reweighted(early.clone(), None, CUMULANT).unwrap_err();
    assert_eq!(error, before);
    assert_eq!(
        error.to_string(),
        "No sample can be reweighted. Every sample lies before the first sampled step."
    );
    // The furthest reason wins: one before, one without a value in "b".
    let missing = run(
        1_000,
        &[0, 1_000],
        &[Some(0.0), Some(0.0)],
        Some((&[Some(1.0), Some(1.0)], &[Some(0.0), None])),
    );
    let error = Distribution::weighable(&missing, &Samples::All, None).unwrap_err();
    assert_eq!(
        error,
        StructureSetError::NoWeighableSample(Unweighable::MissingBias { key: key("b") })
    );
    assert_eq!(
        error.to_string(),
        "No sample can be reweighted. Samples lack values in the bias column \"b\"."
    );
    // A failed reweighting changes nothing.
    let mut ensemble = Ensemble::uniform(missing);
    let before = ensemble.distribution().clone();
    assert!(ensemble.reweight(None, CUMULANT).is_err());
    assert_eq!(ensemble.distribution(), &before);
    // No proposed sample is no error.
    let empty = Distribution::weighable(ensemble.set(), &Samples::Subset(Vec::new()), None)
        .expect("no sample is proposed");
    assert!(empty.samples().is_empty() && empty.left_out().is_empty());
}

#[test]
fn mixed_samples_are_weighed_by_the_rule_reweighting_applies() -> Result<(), StructureSetError> {
    let set = run(
        1_000,
        &[500, 1_000, 1_500, 2_000, 2_500],
        &[Some(-1.0), Some(-2.0), None, Some(-3.0), Some(-4.0)],
        Some((
            &[Some(1.0), None, Some(2.0), Some(3.0), None],
            &[Some(0.5), None, Some(0.5), Some(0.5), Some(0.5)],
        )),
    );
    let weighable = Distribution::weighable(&set, &Samples::All, None)?;
    assert_eq!(weighable.samples(), [2, 3]);
    // The first bias column without a value names the reason.
    assert_eq!(
        reasons(weighable.left_out()),
        [
            (0, Unweighable::BeforeSampling),
            (1, Unweighable::MissingBias { key: key("a") }),
            (4, Unweighable::MissingBias { key: key("a") }),
        ]
    );
    // At another temperature the potential energy is read too.
    let warmer = Distribution::weighable(&set, &Samples::All, Some(kelvin(310.0)))?;
    assert_eq!(warmer.samples(), [3]);
    assert_eq!(
        reasons(warmer.left_out()),
        [
            (0, Unweighable::BeforeSampling),
            (1, Unweighable::MissingBias { key: key("a") }),
            (2, Unweighable::MissingPotentialEnergy),
            (4, Unweighable::MissingBias { key: key("a") }),
        ]
    );
    // The strict build refuses the first left-out sample, and takes the
    // weighable ones with the reduced biases the lenient one gives.
    assert_eq!(
        Distribution::reweighted(
            &set,
            Samples::Subset(vec![2, 4]),
            Partition::one_state(2),
            None,
            CUMULANT
        )
        .unwrap_err(),
        StructureSetError::MissingColumnValue {
            key: key("a"),
            geometry: 4
        }
    );
    assert_eq!(
        Distribution::reweighted(
            &set,
            Samples::Subset(vec![2]),
            Partition::one_state(1),
            Some(kelvin(310.0)),
            CUMULANT
        )
        .unwrap_err(),
        StructureSetError::MissingColumnValue {
            key: key("potential"),
            geometry: 2
        }
    );
    let strict = Distribution::reweighted(
        &set,
        Samples::Subset(weighable.samples().to_vec()),
        Partition::one_state(2),
        None,
        CUMULANT,
    )?;
    let mut ensemble = Ensemble::uniform(set);
    let reweighting = ensemble.reweight(None, CUMULANT)?;
    assert_eq!(reweighting.left_out(), weighable.left_out());
    assert_eq!(ensemble.distribution(), &strict);
    assert_eq!(reduced(ensemble.distribution()), reduced(&strict));
    Ok(())
}

#[test]
fn states_left_empty_go_and_the_others_are_renumbered() -> Result<(), StructureSetError> {
    // Bias values missing for geometries 1, 2, and 4.
    let set = run(
        0,
        &[0, 1, 2, 3, 4, 5],
        &[Some(0.0); 6],
        Some((
            &[Some(1.0), None, None, Some(2.0), None, Some(3.0)],
            &[Some(0.0); 6],
        )),
    );
    // State 0: 0, 3 (represented by 3); state 1: 1, 2; state 2: 4, 5
    // (represented by 4).
    let partition = Partition::new(vec![0, 1, 1, 0, 2, 2], vec![3, 1, 4]).unwrap();
    let mut ensemble = Ensemble::uniform(set);
    ensemble.distribution_mut().regroup(partition)?;
    let energy = Quantity::new(2.0, KILOJOULE_PER_MOLE);
    ensemble
        .distribution_mut()
        .state_mut(2)
        .unwrap()
        .set_energy(Some(energy))?;
    let reweighting = ensemble.reweight(Some(kelvin(300.0)), Estimator::Exponential)?;
    assert_eq!(reweighting.states(), [Some(0), None, Some(1)]);
    assert_eq!(
        reweighting
            .left_out()
            .iter()
            .map(|left_out| left_out.geometry())
            .collect::<Vec<_>>(),
        [1, 2, 4]
    );
    let distribution = ensemble.distribution();
    assert_eq!(distribution.samples(), &Samples::Subset(vec![0, 3, 5]));
    // State 0 keeps its representative; state 2's went, so its first kept
    // sample represents it, and it keeps its energy.
    assert_eq!(
        distribution.partition(),
        &Partition::new(vec![0, 0, 1], vec![1, 2]).unwrap()
    );
    assert_eq!(distribution.states()[0], StateData::new());
    assert_eq!(distribution.states()[1].energy(), Some(energy));
    // Reweighting again leaves nothing out and keeps every state.
    let again = ensemble.reweight(None, Estimator::Exponential)?;
    assert!(again.left_out().is_empty());
    assert_eq!(again.states(), [Some(0), Some(1)]);
    Ok(())
}

#[test]
fn an_unbiased_set_at_another_temperature_reads_the_potential_energy(
) -> Result<(), StructureSetError> {
    let set = run(0, &[0, 1, 2], &[Some(-1.0), None, Some(-3.0)], None);
    // At the sampled temperature every geometry weighs the same.
    let same = Distribution::weighable(&set, &Samples::All, None)?;
    assert_eq!(same.samples(), [0, 1, 2]);
    let ensemble = Ensemble::reweighted(set.clone(), None, CUMULANT)?;
    assert_eq!(reduced(ensemble.distribution()), [0.0; 3]);
    // At 350 K the geometry without a potential energy is left out.
    let warmer = Distribution::weighable(&set, &Samples::All, Some(kelvin(350.0)))?;
    assert_eq!(warmer.samples(), [0, 2]);
    assert_eq!(
        reasons(warmer.left_out()),
        [(1, Unweighable::MissingPotentialEnergy)]
    );
    let mut ensemble = Ensemble::uniform(set.clone());
    let reweighting = ensemble.reweight(Some(kelvin(350.0)), CUMULANT)?;
    assert_eq!(reweighting.left_out(), warmer.left_out());
    let r = 8.314_462_618_153_24e-3;
    let difference = 1.0 / (r * 300.0) - 1.0 / (r * 350.0);
    let b = reduced(ensemble.distribution());
    assert!((b[0] + difference).abs() < 1e-12);
    assert!((b[1] - difference * -3.0).abs() < 1e-12);
    // Without a potential-energy column, another temperature is refused.
    let mut unnamed = set;
    unnamed
        .set_sampling(Some(Sampling::new(
            ThermodynamicState::canonical(kelvin(300.0)).unwrap(),
        )))
        .unwrap();
    assert_eq!(
        Distribution::weighable(&unnamed, &Samples::All, Some(kelvin(350.0))).unwrap_err(),
        StructureSetError::Missing(Requirement::PotentialEnergy)
    );
    Ok(())
}

#[test]
fn proposed_samples_are_checked() {
    let set = run(0, &[0, 1], &[Some(0.0); 2], None);
    assert_eq!(
        Distribution::weighable(&set, &Samples::Subset(vec![1, 0]), None).unwrap_err(),
        StructureSetError::UnsortedSamples
    );
    assert!(matches!(
        Distribution::weighable(&set, &Samples::Subset(vec![0, 2]), None).unwrap_err(),
        StructureSetError::IndexOutOfRange(_)
    ));
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

/// Over random runs, partitions, and temperatures, the query and the
/// reweighting agree: the reweighted samples are the weighable ones, every
/// proposed sample is weighed or left out, and each kept state keeps its
/// kept samples in order.
#[test]
fn the_query_and_the_reweighting_agree_on_random_runs() {
    let mut random = Lcg(7);
    for case in 0..200 {
        let len = 1 + random.next(8) as usize;
        let steps = (0..len as u64).map(|step| step * 10).collect::<Vec<_>>();
        let potential = (0..len).map(|_| random.cell()).collect::<Vec<_>>();
        let a = (0..len).map(|_| random.cell()).collect::<Vec<_>>();
        let b = (0..len).map(|_| random.cell()).collect::<Vec<_>>();
        let biased = random.next(4) != 0;
        let set = run(
            random.next(3) * 10,
            &steps,
            &potential,
            biased.then_some((&a[..], &b[..])),
        );
        let temperature = (random.next(2) == 0).then(|| kelvin(320.0));
        let states = 1 + random.next(len as u64) as usize;
        let assignments = (0..len)
            .map(|sample| {
                if sample < states {
                    sample
                } else {
                    random.next(states as u64) as usize
                }
            })
            .collect::<Vec<_>>();
        let partition = Partition::new(assignments.clone(), (0..states).collect()).unwrap();
        let mut ensemble = Ensemble::uniform(set.clone());
        ensemble.distribution_mut().regroup(partition).unwrap();
        let query = Distribution::weighable(&set, &Samples::All, temperature);
        let outcome = ensemble.reweight(temperature, CUMULANT);
        match (query, outcome) {
            (Ok(weighable), Ok(reweighting)) => {
                assert_eq!(reweighting.left_out(), weighable.left_out(), "case {case}");
                let distribution = ensemble.distribution();
                let samples = (0..distribution.sample_count())
                    .map(|position| match distribution.samples() {
                        Samples::All => position,
                        Samples::Subset(indices) => indices[position],
                        _ => unreachable!(),
                    })
                    .collect::<Vec<_>>();
                assert_eq!(samples, weighable.samples(), "case {case}");
                assert_eq!(
                    samples.len() + weighable.left_out().len(),
                    len,
                    "case {case}"
                );
                // Each kept sample is in its old state's new index.
                for (position, &geometry) in samples.iter().enumerate() {
                    assert_eq!(
                        Some(distribution.partition().assignments()[position]),
                        reweighting.states()[assignments[geometry]],
                        "case {case}"
                    );
                }
                let kept = reweighting.states().iter().flatten().count();
                assert_eq!(kept, distribution.state_count(), "case {case}");
                assert!(reweighting.states().iter().flatten().copied().eq(0..kept));
            }
            (Err(query), Err(outcome)) => assert_eq!(query, outcome, "case {case}"),
            (query, outcome) => panic!("case {case}: {query:?} but {outcome:?}"),
        }
    }
}
