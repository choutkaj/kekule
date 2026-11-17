//! Ensembles over structure sets: selecting, removing, pushing, regrouping,
//! and collapsing keep the order and the distribution on the set.

use std::sync::Arc;

use kekule::diagnostics::{Countable, IndexOutOfRange, LengthMismatch};
use kekule::math::Point3;
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::*;
use kekule::topology::{TopologyAtomIndex, TopologyBuilder};
use kekule::units::{Quantity, KELVIN, KILOJOULE_PER_MOLE, NANOMETER};

const ATOM0: TopologyAtomIndex = TopologyAtomIndex::new(0);

/// Geometries of a one-atom topology at these x positions, unordered and
/// without a distribution.
fn set_at(xs: &[f64]) -> StructureSet {
    StructureSet::from_geometries(
        single_atom_topology(),
        xs.iter().map(|x| Geometry::new(single_position(*x))),
    )
    .unwrap()
}

fn xs(set: &StructureSet) -> Vec<f64> {
    set.iter()
        .map(|geometry| geometry.positions().values().value()[0].x)
        .collect()
}

fn distribution(
    samples: Samples,
    assignments: &[usize],
    representatives: &[usize],
) -> Distribution {
    let partition = Partition::new(assignments.to_vec(), representatives.to_vec()).unwrap();
    Distribution::new(samples, partition).unwrap()
}

/// The geometries at these x positions as `distribution` reads them.
fn ensemble_at(xs: &[f64], distribution: Distribution) -> Ensemble {
    Ensemble::new(set_at(xs), distribution).unwrap()
}

fn shares_every_geometry(a: &StructureSet, b: &StructureSet) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b.iter())
            .all(|(a, b)| a.payload().shares_storage_with(b.payload()))
}

/// The parts of two ensembles compared one by one, since a set is not
/// `PartialEq`; geometries must share their structures.
fn assert_same_ensemble(actual: &Ensemble, expected: &Ensemble) {
    assert!(shares_every_geometry(actual.set(), expected.set()));
    assert_eq!(actual.set().order(), expected.set().order());
    assert_eq!(actual.distribution(), expected.distribution());
    assert_eq!(actual.set().properties(), expected.set().properties());
}

#[test]
fn log_weights_are_finite_and_equal_when_absent() {
    let mut set = ensemble_at(
        &[0.0, 1.0, 2.0],
        distribution(Samples::Subset(vec![0, 2]), &[0, 0], &[0]),
    );
    let mut distribution = set.distribution_mut();
    assert_eq!(distribution.weight_model(), WeightModel::Equal);
    assert_eq!(distribution.log_weight(2), Some(0.0));
    assert_eq!(distribution.log_weight(1), None);
    assert_eq!(distribution.probabilities().unwrap(), [0.5, 0.5]);
    assert_eq!(distribution.effective_sample_size().unwrap(), 2.0);
    assert_eq!(
        distribution.set_log_weights(vec![0.0]),
        Err(StructureSetError::LogWeightCountMismatch(
            LengthMismatch::new(Countable::LogWeight, 2, 1)
        ))
    );
    assert_eq!(
        distribution.set_log_weights(vec![0.0, f64::NAN]),
        Err(StructureSetError::NonFiniteLogWeight { geometry: 2 })
    );
    assert_eq!(
        distribution.set_log_weight(1, 2.0),
        Err(StructureSetError::NotSample { geometry: 1 })
    );
    // Failed edits keep every sample equally weighted.
    assert_eq!(distribution.weight_model(), WeightModel::Equal);
    let three = 3.0_f64.ln();
    distribution.set_log_weight(2, three).unwrap();
    assert_eq!(
        distribution.weight_model(),
        WeightModel::LogWeights(&[0.0, three])
    );
    assert_eq!(distribution.ranked_by_weight(), [2, 0]);
    distribution.normalize_log_weights().unwrap();
    let WeightModel::LogWeights(normalized) = distribution.weight_model() else {
        panic!("normalized weights are log weights")
    };
    let normalized = normalized.to_vec();
    assert_close(&exponentials(&normalized), &[0.25, 0.75]);
    // Normalized weights sum to one, so the only state's log population is 0.
    let log_population = set.state(0).unwrap().log_population();
    assert!(
        log_population.abs() <= 4.0 * f64::EPSILON,
        "{log_population}"
    );

    // A new sample needs a log weight once the samples have log weights.
    assert_eq!(
        set.push_sample(Geometry::new(single_position(3.0)), None, 0, None),
        Err(StructureSetError::LogWeightRequired)
    );
    assert_eq!(
        set.push_sample(Geometry::new(single_position(3.0)), None, 1, Some(0.0)),
        Err(StructureSetError::StateOutOfRange(IndexOutOfRange::new(
            Countable::State,
            1,
            1
        )))
    );
    assert_eq!(
        set.push_sample(
            Geometry::new(single_position(3.0)),
            None,
            0,
            Some(f64::INFINITY)
        ),
        Err(StructureSetError::NonFiniteLogWeight { geometry: 3 })
    );
    assert_eq!(set.set().len(), 3);
    let half = 0.5_f64.ln();
    assert_eq!(
        set.push_sample(Geometry::new(single_position(3.0)), None, 0, Some(half)),
        Ok(3)
    );
    let state = StateData::new()
        .with_energy(Quantity::new(1.0, KILOJOULE_PER_MOLE))
        .unwrap();
    assert_eq!(
        set.push_state(Geometry::new(single_position(4.0)), None, state, Some(0.0)),
        Ok(1)
    );
    // A geometry that is not a sample keeps the distribution as it was.
    assert_eq!(set.push(Geometry::new(single_position(5.0)), None), Ok(5));
    let distribution = set.distribution();
    assert_eq!(distribution.samples(), &Samples::Subset(vec![0, 2, 3, 4]));
    assert_eq!(
        distribution.weight_model(),
        WeightModel::LogWeights(&[normalized[0], normalized[1], half, 0.0])
    );
    assert_eq!(distribution.partition().assignments(), [0, 0, 0, 1]);
    assert_eq!(distribution.partition().representatives(), [0, 3]);
    assert_eq!(
        distribution.ranked_by_energy(),
        Err(StructureSetError::Missing(Requirement::StateEnergies))
    );
    assert_eq!(set.state(1).unwrap().representative(), 4);
    assert_eq!(
        set.state(1).unwrap().energy(),
        Some(Quantity::new(1.0, KILOJOULE_PER_MOLE))
    );
    // Pushing a geometry that is not a sample into a distribution of every
    // geometry lists the others.
    let mut all = Ensemble::uniform(set_at(&[0.0]));
    all.push(Geometry::new(single_position(1.0)), None).unwrap();
    assert_eq!(all.distribution().samples(), &Samples::Subset(vec![0]));
}

#[test]
fn selecting_and_removing_geometries_remaps_the_samples() {
    // Five frames; frames 1-4 are an ensemble of two states.
    let mut frames = set_at(&[0.0, 1.0, 2.0, 3.0, 4.0]);
    frames
        .set_order(Some(Order::Step(vec![0, 1, 2, 3, 4])))
        .unwrap();
    let energies = [1.0, 2.0].map(|energy| {
        StateData::new()
            .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
            .unwrap()
    });
    let [one, two, three, four] = [1.0_f64, 2.0, 3.0, 4.0].map(f64::ln);
    let set = Ensemble::new(
        frames,
        distribution(Samples::Subset(vec![1, 2, 3, 4]), &[0, 1, 0, 1], &[0, 3])
            .with_log_weights(vec![one, two, three, four])
            .unwrap()
            .with_states(energies.to_vec())
            .unwrap(),
    )
    .unwrap();

    // Selecting 4, 0, 3 keeps part of each state and a geometry that is not
    // a sample.
    let selected = set.select([4, 0, 3]).unwrap();
    assert_eq!(xs(selected.set()), [4.0, 0.0, 3.0]);
    assert_eq!(selected.set().order(), Some(&Order::Step(vec![4, 0, 3])));
    let distribution = selected.distribution();
    assert_eq!(distribution.samples(), &Samples::Subset(vec![0, 2]));
    assert_eq!(
        distribution.weight_model(),
        WeightModel::LogWeights(&[four, three])
    );
    // States in state order: state 0 (geometry 3, a new representative),
    // then state 1 (geometry 4, its representative).
    assert_eq!(distribution.partition().assignments(), [1, 0]);
    assert_eq!(distribution.partition().representatives(), [1, 0]);
    assert_eq!(selected.state(0).unwrap().energy(), None);
    assert_eq!(selected.state(1).unwrap().energy(), None);

    // Selecting a whole state keeps its energy.
    let whole = set.select([1, 3]).unwrap();
    assert_eq!(whole.distribution().samples(), &Samples::All);
    assert_eq!(
        whole.state(0).unwrap().energy(),
        Some(Quantity::new(1.0, KILOJOULE_PER_MOLE))
    );

    // Selecting whole states keeps their samples and order values.
    let states = set.select_states([1]).unwrap();
    assert_eq!(xs(states.set()), [2.0, 4.0]);
    assert_eq!(states.set().order(), Some(&Order::Step(vec![2, 4])));
    assert_eq!(
        states.distribution().weight_model(),
        WeightModel::LogWeights(&[two, four])
    );
    assert_eq!(
        states.state(0).unwrap().energy(),
        Some(Quantity::new(2.0, KILOJOULE_PER_MOLE))
    );
    assert_eq!(
        set.select_states([2]).unwrap_err(),
        StructureSetError::StateOutOfRange(IndexOutOfRange::new(Countable::State, 2, 2))
    );

    // Removing a geometry that is not a sample shifts the samples down.
    let mut removed = set.clone();
    removed.remove(0).unwrap();
    assert_eq!(removed.distribution().samples(), &Samples::All);
    assert_eq!(removed.set().order(), Some(&Order::Step(vec![1, 2, 3, 4])));
    assert_eq!(removed.state(1).unwrap().representative(), 3);
    // The representative of a state with other samples stays.
    assert_eq!(
        removed.remove(0).unwrap_err(),
        StructureSetError::RepresentativeRemoval {
            geometry: 0,
            state: 0
        }
    );
    removed.remove(2).unwrap();
    removed.remove(1).unwrap();
    assert_eq!(xs(removed.set()), [1.0, 4.0]);
    assert_eq!(removed.distribution().partition().assignments(), [0, 1]);
    // Removing a sample alone in its state removes the state.
    removed.remove(0).unwrap();
    let distribution = removed.distribution();
    assert_eq!(distribution.state_count(), 1);
    assert_eq!(distribution.samples(), &Samples::All);
    assert_eq!(
        distribution.weight_model(),
        WeightModel::LogWeights(&[four])
    );
    assert_eq!(distribution.partition().representatives(), [0]);
    assert_eq!(
        removed.state(0).unwrap().energy(),
        Some(Quantity::new(2.0, KILOJOULE_PER_MOLE))
    );
}

#[test]
fn regrouping_reweighting_and_selecting_copy_no_structure() {
    let mut set = set_at(&[0.0, 1.0, 2.0]);
    set.set_order(Some(Order::Step(vec![0, 1, 2]))).unwrap();
    let shared = Arc::new(Ensemble::uniform(set));

    let mut regrouped = Arc::clone(&shared);
    let mut distribution = Arc::make_mut(&mut regrouped).distribution_mut();
    distribution
        .regroup(Partition::state_per_sample(3))
        .unwrap();
    let weights = [1.0_f64, 2.0, 3.0].map(f64::ln);
    distribution.set_log_weights(weights.to_vec()).unwrap();
    assert!(!Arc::ptr_eq(&regrouped, &shared));
    assert!(shares_every_geometry(regrouped.set(), shared.set()));
    assert_eq!(shared.distribution().state_count(), 1);

    let selected = shared.select([2, 0]).unwrap();
    assert!(selected
        .set()
        .get(0)
        .unwrap()
        .payload()
        .shares_storage_with(shared.set().get(2).unwrap().payload()));
    // Each state of one sample keeps its log weight as its log population.
    let collapsed = regrouped.collapse();
    assert_eq!(collapsed.set().order(), None);
    assert_eq!(
        collapsed.distribution().weight_model(),
        WeightModel::LogWeights(&weights)
    );
    assert!(shares_every_geometry(collapsed.set(), regrouped.set()));
    let states = regrouped.select_states([1]).unwrap();
    assert!(states
        .set()
        .get(0)
        .unwrap()
        .payload()
        .shares_storage_with(shared.set().get(1).unwrap().payload()));

    // Editing one geometry of a clone copies that geometry alone.
    let mut edited = (*shared).clone();
    edited
        .get_mut(1)
        .unwrap()
        .set_position_at(ATOM0, Quantity::new(Point3::new(9.0, 0.0, 0.0), NANOMETER))
        .unwrap();
    assert!(edited
        .set()
        .get(0)
        .unwrap()
        .payload()
        .shares_storage_with(shared.set().get(0).unwrap().payload()));
    assert!(!edited
        .set()
        .get(1)
        .unwrap()
        .shares_storage_with(&shared.set().get(1).unwrap()));
    assert_eq!(xs(shared.set()), [0.0, 1.0, 2.0]);
    assert_eq!(xs(edited.set()), [0.0, 9.0, 2.0]);

    // Replacing positions copies only what the shared geometry keeps.
    let mut moved = shared.set().clone();
    moved
        .replace_positions([3.0, 4.0, 5.0].map(single_position))
        .unwrap();
    assert_eq!(xs(&moved), [3.0, 4.0, 5.0]);
    assert_eq!(xs(shared.set()), [0.0, 1.0, 2.0]);
    assert_eq!(
        moved.replace_positions([single_position(0.0)]),
        Err(StructureSetError::GeometryCountMismatch(
            LengthMismatch::new(Countable::Geometry, 3, 1)
        ))
    );
}

#[test]
fn subsets_keep_the_order_and_the_distribution() {
    let (structure, atom, _) = structure_fixture();
    let mut set = StructureSet::from_structures([structure.clone(), structure]).unwrap();
    set.set_order(Some(Order::Step(vec![5, 6]))).unwrap();
    let state = StateData::new()
        .with_energy(Quantity::new(1.0, KILOJOULE_PER_MOLE))
        .unwrap();
    let mut state = state;
    state
        .properties_mut()
        .insert(key("label"), PropertyValue::String("a".into()))
        .unwrap();
    set.properties_mut()
        .insert(key("method"), PropertyValue::String("md".into()))
        .unwrap();
    set.get_mut(0)
        .unwrap()
        .set_velocities(Some(Velocities::zeros(2)))
        .unwrap();
    let ensemble = Ensemble::new(
        set.clone(),
        distribution(Samples::Subset(vec![1]), &[0], &[0])
            .with_log_weights(vec![2.0])
            .unwrap()
            .with_states(vec![state])
            .unwrap()
            .with_target(ThermodynamicState::canonical(Quantity::new(300.0, KELVIN)).unwrap())
            .unwrap(),
    )
    .unwrap();
    let selection =
        kekule::topology::AtomSelection::from_atoms(set.shared_topology(), [atom]).unwrap();
    let subset = ensemble.subset(&selection).unwrap();
    assert_eq!(subset.set().topology().atom_count(), 1);
    assert!(subset.set().properties().is_empty());
    assert_eq!(subset.set().order(), set.order());
    assert_eq!(subset.set().get(0).unwrap().velocities().unwrap().len(), 1);
    let state = subset.state(0).unwrap();
    assert_eq!(state.energy(), Some(Quantity::new(1.0, KILOJOULE_PER_MOLE)));
    assert!(state.properties().is_empty());
    let distribution = subset.distribution();
    assert_eq!(distribution.samples(), &Samples::Subset(vec![1]));
    assert_eq!(distribution.weight_model(), WeightModel::LogWeights(&[2.0]));
    assert_eq!(
        distribution.temperature(),
        Some(Quantity::new(300.0, KELVIN))
    );
    // The set's own subset is the same projection, without a distribution.
    let set_subset = set.subset(&selection).unwrap();
    assert_eq!(set_subset.topology().atom_count(), 1);
    assert_eq!(set_subset.order(), set.order());

    let mut perceived = set.clone();
    perceived.perceive().unwrap();
    assert!(perceived.topology().shares_layout(set.topology()));
    assert!(shares_every_geometry(&perceived, &set));
    assert_eq!(
        StructureSet::from_structures(std::iter::empty()).unwrap_err(),
        StructureSetError::EmptySource
    );
}

#[test]
fn duplicates_in_a_selection_become_separate_geometries() {
    let mut frames = set_at(&[0.0, 1.0]);
    frames.set_order(Some(Order::Step(vec![5, 6]))).unwrap();
    let two = 2.0_f64.ln();
    let set = Ensemble::new(
        frames,
        distribution(Samples::All, &[0, 1], &[0, 1])
            .with_log_weights(vec![0.0, two])
            .unwrap()
            .with_states(vec![
                StateData::new()
                    .with_energy(Quantity::new(1.0, KILOJOULE_PER_MOLE))
                    .unwrap(),
                StateData::new()
                    .with_energy(Quantity::new(2.0, KILOJOULE_PER_MOLE))
                    .unwrap(),
            ])
            .unwrap(),
    )
    .unwrap();
    let twice = set.select([1, 1]).unwrap();
    assert_eq!(xs(twice.set()), [1.0, 1.0]);
    assert_eq!(twice.set().order(), Some(&Order::Step(vec![6, 6])));
    let distribution = twice.distribution();
    assert_eq!(distribution.samples(), &Samples::All);
    assert_eq!(
        distribution.weight_model(),
        WeightModel::LogWeights(&[two, two])
    );
    assert_eq!(distribution.partition().assignments(), [0, 0]);
    assert_eq!(distribution.partition().representatives(), [0]);
    // Every sample of state 1 is selected, so it keeps its energy.
    assert_eq!(twice.state(0).unwrap().energy(), kj(2.0));
    assert!(twice
        .set()
        .get(0)
        .unwrap()
        .shares_storage_with(&twice.set().get(1).unwrap()));

    let states = set.select_states([0, 0]).unwrap();
    assert_eq!(xs(states.set()), [0.0, 0.0]);
    assert_eq!(states.set().order(), Some(&Order::Step(vec![5, 5])));
    let distribution = states.distribution();
    assert_eq!(
        distribution.weight_model(),
        WeightModel::LogWeights(&[0.0, 0.0])
    );
    assert_eq!(distribution.state_count(), 2);
    assert_eq!(distribution.partition().assignments(), [0, 1]);
    assert_eq!(distribution.partition().representatives(), [0, 1]);
    assert_eq!(states.state(1).unwrap().energy(), kj(1.0));
}

#[test]
fn collapse_keeps_population_shares_of_any_samples() {
    // Samples 0 and 2 of three geometries, one state, without weights.
    let mut frames = set_at(&[0.0, 1.0, 2.0]);
    frames.set_order(Some(Order::Step(vec![0, 1, 2]))).unwrap();
    let set = Ensemble::new(
        frames,
        distribution(Samples::Subset(vec![0, 2]), &[0, 0], &[1]),
    )
    .unwrap();
    let collapsed = set.collapse();
    assert_eq!(xs(collapsed.set()), [2.0]);
    assert_eq!(collapsed.set().order(), None);
    // Two samples of weight 1 make a population of 2.
    assert_eq!(
        collapsed.distribution().weight_model(),
        WeightModel::LogWeights(&[2.0_f64.ln()])
    );
    assert!(collapsed
        .set()
        .get(0)
        .unwrap()
        .shares_storage_with(&set.set().get(2).unwrap()));

    // Weights of e¹⁰⁰⁰, whose sum no linear weight holds, keep their ratios.
    let heavy = ensemble_at(
        &[0.0, 1.0, 2.0],
        distribution(Samples::All, &[0, 0, 1], &[0, 2])
            .with_log_weights(vec![1000.0; 3])
            .unwrap(),
    );
    assert_eq!(
        heavy.state(0).unwrap().log_population(),
        1000.0 + 2.0_f64.ln()
    );
    // 1000 + ln 2 rounds to the spacing of doubles near 1000, 1.1e-13.
    let shares = heavy.distribution().state_probabilities().unwrap();
    for (share, expected) in shares.iter().zip([2.0 / 3.0, 1.0 / 3.0]) {
        assert!((share - expected).abs() <= 1.0e-13, "{shares:?}");
    }
    let collapsed = heavy.collapse();
    assert_eq!(
        collapsed.distribution().weight_model(),
        WeightModel::LogWeights(&[1000.0 + 2.0_f64.ln(), 1000.0])
    );
}

#[test]
fn samples_pushed_into_a_distribution_of_every_geometry() {
    let mut frames = set_at(&[0.0]);
    frames.set_order(Some(Order::Step(vec![0]))).unwrap();
    let mut set = Ensemble::uniform(frames);
    assert_eq!(
        set.push_sample(
            Geometry::new(single_position(1.0)),
            Some(Place::Step(1)),
            0,
            None
        ),
        Ok(1)
    );
    assert_eq!(
        set.push_state(
            Geometry::new(single_position(2.0)),
            Some(Place::Step(2)),
            StateData::new(),
            None,
        ),
        Ok(1)
    );
    let distribution = set.distribution();
    assert_eq!(distribution.samples(), &Samples::All);
    assert_eq!(distribution.weight_model(), WeightModel::Equal);
    assert_eq!(distribution.partition().assignments(), [0, 0, 1]);
    assert_eq!(distribution.partition().representatives(), [0, 2]);
    assert_eq!(set.set().order(), Some(&Order::Step(vec![0, 1, 2])));

    // A geometry that does not bind, checked last, changes nothing.
    let before = set.clone();
    let wide = || Geometry::new(Positions::zeros(2));
    assert!(matches!(
        set.push_sample(wide(), Some(Place::Step(3)), 1, None),
        Err(StructureSetError::Geometry(_))
    ));
    assert!(matches!(
        set.push_state(wide(), Some(Place::Step(3)), StateData::new(), None),
        Err(StructureSetError::Geometry(_))
    ));
    assert!(matches!(
        set.push(wide(), Some(Place::Step(3))),
        Err(StructureSetError::Geometry(_))
    ));
    assert_same_ensemble(&set, &before);
}

fn single_atom_topology() -> kekule::topology::Topology {
    let mut editor = MoleculeEditor::new();
    editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    kekule::topology::Topology::from_molecule(editor.finish().unwrap()).unwrap()
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

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
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

/// The exponential of each log weight.
fn exponentials(log_weights: &[f64]) -> Vec<f64> {
    log_weights
        .iter()
        .map(|log_weight| log_weight.exp())
        .collect()
}

fn kj(value: f64) -> Option<Quantity<f64>> {
    Some(Quantity::new(value, KILOJOULE_PER_MOLE))
}
