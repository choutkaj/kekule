//! Sample sets: the order and the distribution, copy on write, and lossless
//! conversion from and to ensembles, trajectories, and series.

use super::*;
use crate::units::{DEGREE, PICOSECOND};

/// Samples of a one-atom topology at these x positions, unordered and
/// without a distribution.
fn set_at(xs: &[f64]) -> SampleSet {
    SampleSet::from_samples(
        single_atom_topology(),
        xs.iter().map(|x| Sample::new(single_position(*x))),
    )
    .unwrap()
}

fn xs(set: &SampleSet) -> Vec<f64> {
    set.iter()
        .map(|sample| sample.positions().values().value()[0].x)
        .collect()
}

fn ps(value: f64) -> Quantity<f64> {
    Quantity::new(value, PICOSECOND)
}

fn timed(time: f64, step: u64) -> Option<Place> {
    Some(Place::Time {
        time: Some(ps(time)),
        step: Some(step),
    })
}

fn time_axis(times: Vec<f64>, steps: Vec<u64>) -> Order {
    Order::Time(TimeAxis::new(Some(Quantity::new(times, PICOSECOND)), Some(steps)).unwrap())
}

fn distribution(
    members: Members,
    assignments: &[usize],
    representatives: &[usize],
) -> Distribution {
    let partition = EnsemblePartition::new(assignments.to_vec(), representatives.to_vec()).unwrap();
    Distribution::new(members, partition).unwrap()
}

fn velocities(x: f64) -> Velocities {
    Velocities::new(Quantity::new(
        [Vector3::new(x, 0.0, 0.0)],
        CANONICAL_VELOCITY,
    ))
    .unwrap()
}

const CANONICAL_VELOCITY: crate::units::Unit = crate::units::CANONICAL_VELOCITY_UNIT;

fn shares_every_sample(a: &SampleSet, b: &SampleSet) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b.iter())
            .all(|(a, b)| Arc::ptr_eq(a.payload().structure(), b.payload().structure()))
}

#[test]
fn requirements_read_as_they_are() {
    assert_eq!(Requirement::TimeOrder.to_string(), "needs a time order");
    assert_eq!(
        Requirement::StateEnergies.to_string(),
        "needs an energy for every state"
    );
    assert_eq!(
        SampleSetError::Missing(Requirement::Distribution).to_string(),
        "The sample set needs a distribution."
    );
    assert_eq!(
        SampleSetError::NotConvertible {
            target: ConversionTarget::Trajectory,
            content: Unconvertible::Occupancies { sample: 3 },
        }
        .to_string(),
        "A trajectory cannot hold the occupancies of sample 3."
    );
    let set = set_at(&[0.0]);
    assert_eq!(
        set.validate_monotonic_time(),
        Err(SampleSetError::Missing(Requirement::TimeOrder))
    );
    assert_eq!(
        set.select_states([0]).unwrap_err(),
        SampleSetError::Missing(Requirement::Distribution)
    );
    assert_eq!(
        set.collapse().unwrap_err(),
        SampleSetError::Missing(Requirement::Distribution)
    );
}

#[test]
fn an_order_has_one_value_per_sample_in_every_column() {
    let mut set = set_at(&[0.0, 1.0]);
    assert_eq!(
        set.set_order(Some(Order::Step(vec![0]))),
        Err(SampleSetError::OrderLengthMismatch {
            expected: 2,
            actual: 1
        })
    );
    assert_eq!(
        set.set_order(Some(Order::Coordinate {
            unit: DEGREE,
            values: vec![0.0, f64::NAN],
        })),
        Err(SampleSetError::NonFiniteOrderValue { sample: 1 })
    );
    assert_eq!(set.order(), None);
    assert_eq!(
        TimeAxis::new(Some(Quantity::new(vec![0.0], PICOSECOND)), Some(vec![0, 1])),
        Err(SampleSetError::TimeAxisMismatch { times: 1, steps: 2 })
    );
    // Times convert to picoseconds.
    let axis = TimeAxis::new(
        Some(Quantity::new(
            vec![1000.0, 2000.0],
            crate::units::FEMTOSECOND,
        )),
        None,
    )
    .unwrap();
    assert_eq!(axis.times().unwrap().value(), &[1.0, 2.0]);
    assert_eq!(axis.step(0), None);
    set.set_order(Some(Order::Time(axis))).unwrap();
    assert_eq!(
        set.get(1).unwrap().place(),
        Some(Place::Time {
            time: Some(ps(2.0)),
            step: None
        })
    );

    // A new sample states exactly the columns the axis has.
    let before = set.clone();
    let sample = || Sample::new(single_position(2.0));
    assert_eq!(
        set.push(sample(), None),
        Err(SampleSetError::PlaceMismatch { sample: 2 })
    );
    assert_eq!(
        set.push(
            sample(),
            Some(Place::Time {
                time: None,
                step: None
            })
        ),
        Err(SampleSetError::MissingTime { sample: 2 })
    );
    assert_eq!(
        set.push(sample(), timed(3.0, 7)),
        Err(SampleSetError::UnexpectedStep { sample: 2 })
    );
    assert_eq!(
        set.push(sample(), Some(Place::Step(2))),
        Err(SampleSetError::PlaceMismatch { sample: 2 })
    );
    assert_eq!(set.len(), before.len());
    assert_eq!(set.order(), before.order());
    let place = Some(Place::Time {
        time: Some(ps(3.0)),
        step: None,
    });
    assert_eq!(set.push(sample(), place), Ok(2));
    assert_eq!(set.get(2).unwrap().place(), place);

    // A coordinate converts to the order's unit.
    let mut scan = set_at(&[0.0]);
    scan.set_order(Some(Order::Coordinate {
        unit: NANOMETER,
        values: vec![0.1],
    }))
    .unwrap();
    scan.push(
        Sample::new(single_position(1.0)),
        Some(Place::Coordinate(Quantity::new(2.0, ANGSTROM))),
    )
    .unwrap();
    let Some(Order::Coordinate { unit, values }) = scan.order() else {
        panic!("the scan keeps its coordinate order");
    };
    assert_eq!(*unit, NANOMETER);
    assert!((values[1] - 0.2).abs() < 1e-12);
    assert!(matches!(
        scan.push(
            Sample::new(single_position(1.0)),
            Some(Place::Coordinate(Quantity::new(2.0, PICOSECOND))),
        ),
        Err(SampleSetError::Unit(_))
    ));
    // An unordered set takes no place.
    let mut unordered = set_at(&[0.0]);
    assert_eq!(
        unordered.push(Sample::new(single_position(1.0)), Some(Place::Step(1))),
        Err(SampleSetError::PlaceMismatch { sample: 1 })
    );
    assert_eq!(unordered.get(0).unwrap().place(), None);
}

#[test]
fn order_values_belong_to_the_set_not_to_its_samples() {
    let mut set = set_at(&[0.0, 1.0, 2.0]);
    set.set_order(Some(time_axis(vec![0.0, 2.0, 4.0], vec![0, 10, 20])))
        .unwrap();
    // Selecting keeps each sample's value, in the requested order.
    let reversed = set.select([2, 0]).unwrap();
    assert_eq!(xs(&reversed), [2.0, 0.0]);
    assert_eq!(reversed.get(0).unwrap().place(), timed(4.0, 20));
    assert_eq!(reversed.get(1).unwrap().place(), timed(0.0, 0));
    assert_eq!(
        reversed.validate_monotonic_time(),
        Err(SampleSetError::NonMonotonicTime { sample: 1 })
    );
    assert_eq!(set.validate_monotonic_time(), Ok(()));
    // A removed sample claims no place; the later samples keep theirs.
    let removed = set.remove(1).unwrap();
    assert_eq!(removed, Sample::new(single_position(1.0)));
    assert_eq!(set.get(1).unwrap().place(), timed(4.0, 20));
    set.set_place(1, timed(3.0, 30)).unwrap();
    assert_eq!(set.order(), Some(&time_axis(vec![0.0, 3.0], vec![0, 30])));
    assert_eq!(
        set.set_place(1, Some(Place::Step(1))),
        Err(SampleSetError::PlaceMismatch { sample: 1 })
    );
    // A time order without times cannot be checked for monotonic time.
    let mut untimed = set_at(&[0.0]);
    untimed
        .set_order(Some(Order::Time(TimeAxis::default())))
        .unwrap();
    assert_eq!(
        untimed.validate_monotonic_time(),
        Err(SampleSetError::Missing(Requirement::Times))
    );
    untimed
        .push(Sample::new(single_position(1.0)), None)
        .unwrap();
    assert_eq!(
        untimed.get(1).unwrap().place(),
        Some(Place::Time {
            time: None,
            step: None
        })
    );
}

#[test]
fn any_sample_may_carry_site_fields_and_dynamics() {
    let mut set = set_at(&[0.0]);
    set.set_order(Some(Order::Step(vec![0]))).unwrap();
    let mut sample = set.get_mut(0).unwrap();
    sample.set_occupancies(Some(vec![Some(0.5)])).unwrap();
    sample.set_velocities(Some(velocities(1.0))).unwrap();
    assert_eq!(
        sample.set_forces(Some(Forces::zeros(2))),
        Err(StructureError::Array(DenseArrayError::CountMismatch {
            expected: 1,
            actual: 2
        }))
    );
    let view = set.get(0).unwrap();
    assert_eq!(view.occupancy(ATOM0), Some(0.5));
    let StructureKind::Sample(context) = view.as_model_view().kind() else {
        panic!("a sample's view is a sample");
    };
    assert_eq!(context.place(), Some(Place::Step(0)));
    assert_eq!(context.membership(), None);
    assert_eq!(context.velocities(), Some(&velocities(1.0)));
    assert_eq!(context.forces(), None);
    assert_eq!(view.as_model_view().occupancy(ATOM0), Some(0.5));
    // A model keeps the structure and drops the sample's context.
    let model = view.to_model();
    assert_eq!(model.as_model_view().kind(), StructureKind::Model);
    assert_eq!(model.occupancy(ATOM0), Some(0.5));
    // A rigid motion rotates the velocities with the positions.
    use crate::geometry::{Matrix3, RigidTransform};
    let turn = Matrix3::from_columns(
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    );
    let transform = RigidTransform::from_canonical(turn, Vector3::zero()).unwrap();
    let mut sample = set.get_mut(0).unwrap();
    sample.apply_rigid_transform(&transform).unwrap();
    assert_eq!(
        sample.velocities().unwrap().values().value(),
        &[Vector3::new(0.0, 1.0, 0.0)]
    );
    assert!(set.get_mut(0).unwrap().take_velocities().is_some());
    assert_eq!(set.get(0).unwrap().velocities(), None);
}

#[test]
fn a_distribution_covers_all_samples_or_a_sorted_subset() {
    assert_eq!(
        Distribution::new(Members::Subset(vec![2, 1]), EnsemblePartition::one_state(2)),
        Err(SampleSetError::UnsortedMembers)
    );
    assert_eq!(
        Distribution::new(Members::Subset(vec![1, 2]), EnsemblePartition::one_state(3)),
        Err(SampleSetError::PartitionMismatch {
            expected: 2,
            actual: 3
        })
    );
    let mut set = set_at(&[0.0, 1.0, 2.0]);
    assert_eq!(
        set.set_distribution(Some(distribution(Members::All, &[0, 0], &[0]))),
        Err(SampleSetError::CoverageMismatch {
            covered: 2,
            samples: 3
        })
    );
    assert_eq!(
        set.set_distribution(Some(distribution(
            Members::Subset(vec![1, 3]),
            &[0, 0],
            &[0]
        ))),
        Err(SampleSetError::IndexOutOfRange { index: 3, len: 3 })
    );
    assert_eq!(set.distribution(), None);
    // A subset listing every sample is all of them.
    set.set_distribution(Some(distribution(
        Members::Subset(vec![0, 1, 2]),
        &[0, 0, 1],
        &[0, 2],
    )))
    .unwrap();
    assert_eq!(set.distribution().unwrap().members(), &Members::All);

    set.set_distribution(Some(distribution(
        Members::Subset(vec![1, 2]),
        &[0, 1],
        &[0, 1],
    )))
    .unwrap();
    let distribution = set.distribution().unwrap();
    assert_eq!(distribution.member_position(0), None);
    assert_eq!(distribution.member_position(2), Some(1));
    assert_eq!(distribution.state_of(2), Some(1));
    assert_eq!(set.get(0).unwrap().state().map(|state| state.index()), None);
    assert!(!set.get(0).unwrap().is_representative());
    assert!(set.get(1).unwrap().is_representative());
    assert_eq!(set.state(1).unwrap().representative(), 2);
    assert_eq!(set.state(1).unwrap().samples().collect::<Vec<_>>(), [2]);
    assert_eq!(set.states().len(), 2);
    assert!(set.state(2).is_none());
    // The view's context reports membership only for members.
    let StructureKind::Sample(context) = set.get(0).unwrap().as_model_view().kind() else {
        panic!("a sample's view is a sample");
    };
    assert_eq!(context.membership(), None);
    let StructureKind::Sample(context) = set.get(2).unwrap().as_model_view().kind() else {
        panic!("a sample's view is a sample");
    };
    assert_eq!(
        context.membership(),
        Some(Membership::new(1.0, 1, true, None, None))
    );
}

#[test]
fn weights_are_finite_positive_and_equal_when_absent() {
    let mut set = set_at(&[0.0, 1.0, 2.0]);
    set.set_distribution(Some(distribution(
        Members::Subset(vec![0, 2]),
        &[0, 0],
        &[0],
    )))
    .unwrap();
    let mut distribution = set.distribution_mut().unwrap();
    assert_eq!(distribution.weights(), None);
    assert_eq!(distribution.weight(2), Some(1.0));
    assert_eq!(distribution.weight(1), None);
    assert_eq!(distribution.probabilities().unwrap(), [0.5, 0.5]);
    assert_eq!(distribution.effective_sample_size().unwrap(), 2.0);
    assert_eq!(
        distribution.set_weights(Some(vec![1.0])),
        Err(SampleSetError::WeightCountMismatch {
            expected: 2,
            actual: 1
        })
    );
    assert_eq!(
        distribution.set_weights(Some(vec![1.0, -1.0])),
        Err(SampleSetError::InvalidWeight { sample: 2 })
    );
    assert_eq!(
        distribution.set_weight(1, 2.0),
        Err(SampleSetError::NotMember { sample: 1 })
    );
    distribution.set_weight(2, 3.0).unwrap();
    assert_eq!(distribution.weights(), Some(&[1.0, 3.0][..]));
    assert_eq!(distribution.ranked_by_weight(), [2, 0]);
    distribution.normalize_weights().unwrap();
    assert_eq!(distribution.weights(), Some(&[0.25, 0.75][..]));
    // Populations are relative to the largest weight.
    assert_eq!(set.state(0).unwrap().population(), 1.0 / 0.75);

    // A new member needs a weight once the members have weights.
    assert_eq!(
        set.push_member(Sample::new(single_position(3.0)), None, 0, None),
        Err(SampleSetError::WeightRequired)
    );
    assert_eq!(
        set.push_member(Sample::new(single_position(3.0)), None, 1, Some(1.0)),
        Err(SampleSetError::StateOutOfRange {
            state: 1,
            states: 1
        })
    );
    assert_eq!(set.len(), 3);
    assert_eq!(
        set.push_member(Sample::new(single_position(3.0)), None, 0, Some(0.5)),
        Ok(3)
    );
    let state = EnsembleState::new()
        .with_energy(Quantity::new(1.0, KILOJOULE_PER_MOLE))
        .unwrap();
    assert_eq!(
        set.push_state(Sample::new(single_position(4.0)), None, state, Some(1.0)),
        Ok(1)
    );
    // A non-member keeps the distribution as it was.
    assert_eq!(set.push(Sample::new(single_position(5.0)), None), Ok(5));
    let distribution = set.distribution().unwrap();
    assert_eq!(distribution.members(), &Members::Subset(vec![0, 2, 3, 4]));
    assert_eq!(distribution.weights(), Some(&[0.25, 0.75, 0.5, 1.0][..]));
    assert_eq!(distribution.partition().assignments(), [0, 0, 0, 1]);
    assert_eq!(distribution.partition().representatives(), [0, 3]);
    assert_eq!(
        distribution.ranked_by_energy(),
        Err(SampleSetError::Missing(Requirement::StateEnergies))
    );
    assert_eq!(set.state(1).unwrap().representative(), 4);
    assert_eq!(
        set.state(1).unwrap().energy(),
        Some(Quantity::new(1.0, KILOJOULE_PER_MOLE))
    );
    // Pushing a non-member into a distribution of every sample lists the
    // others.
    let mut all = set_at(&[0.0]);
    all.set_distribution(Some(distribution_of_one_state(1)))
        .unwrap();
    all.push(Sample::new(single_position(1.0)), None).unwrap();
    assert_eq!(
        all.distribution().unwrap().members(),
        &Members::Subset(vec![0])
    );
}

fn distribution_of_one_state(members: usize) -> Distribution {
    Distribution::new(Members::All, EnsemblePartition::one_state(members)).unwrap()
}

#[test]
fn selecting_and_removing_samples_remaps_the_members() {
    // Five frames; frames 1-4 are an ensemble of two states.
    let mut set = set_at(&[0.0, 1.0, 2.0, 3.0, 4.0]);
    set.set_order(Some(Order::Step(vec![0, 1, 2, 3, 4])))
        .unwrap();
    let energies = [1.0, 2.0].map(|energy| {
        EnsembleState::new()
            .with_energy(Quantity::new(energy, KILOJOULE_PER_MOLE))
            .unwrap()
    });
    set.set_distribution(Some(
        distribution(Members::Subset(vec![1, 2, 3, 4]), &[0, 1, 0, 1], &[0, 3])
            .with_weights(vec![1.0, 2.0, 3.0, 4.0])
            .unwrap()
            .with_states(energies.to_vec())
            .unwrap(),
    ))
    .unwrap();

    // Selecting 4, 0, 3 keeps part of each state and a non-member.
    let selected = set.select([4, 0, 3]).unwrap();
    assert_eq!(xs(&selected), [4.0, 0.0, 3.0]);
    assert_eq!(selected.order(), Some(&Order::Step(vec![4, 0, 3])));
    let distribution = selected.distribution().unwrap();
    assert_eq!(distribution.members(), &Members::Subset(vec![0, 2]));
    assert_eq!(distribution.weights(), Some(&[4.0, 3.0][..]));
    // States in state order: state 0 (sample 3, a new representative),
    // then state 1 (sample 4, its representative).
    assert_eq!(distribution.partition().assignments(), [1, 0]);
    assert_eq!(distribution.partition().representatives(), [1, 0]);
    assert_eq!(selected.state(0).unwrap().energy(), None);
    assert_eq!(selected.state(1).unwrap().energy(), None);

    // Selecting a whole state keeps its energy.
    let whole = set.select([1, 3]).unwrap();
    assert_eq!(whole.distribution().unwrap().members(), &Members::All);
    assert_eq!(
        whole.state(0).unwrap().energy(),
        Some(Quantity::new(1.0, KILOJOULE_PER_MOLE))
    );

    // Selecting whole states keeps their members and order values.
    let states = set.select_states([1]).unwrap();
    assert_eq!(xs(&states), [2.0, 4.0]);
    assert_eq!(states.order(), Some(&Order::Step(vec![2, 4])));
    assert_eq!(
        states.distribution().unwrap().weights(),
        Some(&[2.0, 4.0][..])
    );
    assert_eq!(
        states.state(0).unwrap().energy(),
        Some(Quantity::new(2.0, KILOJOULE_PER_MOLE))
    );
    assert_eq!(
        set.select_states([2]).unwrap_err(),
        SampleSetError::StateOutOfRange {
            state: 2,
            states: 2
        }
    );

    // Removing a non-member shifts the members down.
    let mut removed = set.clone();
    removed.remove(0).unwrap();
    assert_eq!(removed.distribution().unwrap().members(), &Members::All);
    assert_eq!(removed.order(), Some(&Order::Step(vec![1, 2, 3, 4])));
    assert_eq!(removed.state(1).unwrap().representative(), 3);
    // The representative of a state with other members stays.
    assert_eq!(
        removed.remove(0).unwrap_err(),
        SampleSetError::RepresentativeRemoval {
            sample: 0,
            state: 0
        }
    );
    removed.remove(2).unwrap();
    removed.remove(1).unwrap();
    assert_eq!(xs(&removed), [1.0, 4.0]);
    assert_eq!(
        removed.distribution().unwrap().partition().assignments(),
        [0, 1]
    );
    // Removing a member alone in its state removes the state.
    removed.remove(0).unwrap();
    let distribution = removed.distribution().unwrap();
    assert_eq!(distribution.state_count(), 1);
    assert_eq!(distribution.members(), &Members::All);
    assert_eq!(distribution.weights(), Some(&[4.0][..]));
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
    set.set_distribution(Some(distribution_of_one_state(3)))
        .unwrap();
    let shared = Arc::new(set);

    let mut regrouped = Arc::clone(&shared);
    let mut distribution = Arc::make_mut(&mut regrouped).distribution_mut().unwrap();
    distribution
        .regroup(EnsemblePartition::state_per_sample(3))
        .unwrap();
    distribution.set_weights(Some(vec![1.0, 2.0, 3.0])).unwrap();
    assert!(!Arc::ptr_eq(&regrouped, &shared));
    assert!(shares_every_sample(&regrouped, &shared));
    assert_eq!(shared.distribution().unwrap().state_count(), 1);

    let selected = shared.select([2, 0]).unwrap();
    assert!(Arc::ptr_eq(
        selected.get(0).unwrap().payload().structure(),
        shared.get(2).unwrap().payload().structure()
    ));
    let collapsed = regrouped.collapse().unwrap();
    assert_eq!(collapsed.order(), None);
    assert_eq!(
        collapsed.distribution().unwrap().weights(),
        Some(&[1.0 / 3.0, 2.0 / 3.0, 1.0][..])
    );
    assert!(shares_every_sample(&collapsed, &regrouped));
    let states = regrouped.select_states([1]).unwrap();
    assert!(Arc::ptr_eq(
        states.get(0).unwrap().payload().structure(),
        shared.get(1).unwrap().payload().structure()
    ));

    // Editing one sample of a clone copies that sample alone.
    let mut edited = (*shared).clone();
    edited
        .get_mut(1)
        .unwrap()
        .set_position_at(ATOM0, Quantity::new(Point3::new(9.0, 0.0, 0.0), NANOMETER))
        .unwrap();
    assert!(Arc::ptr_eq(
        edited.get(0).unwrap().payload().structure(),
        shared.get(0).unwrap().payload().structure()
    ));
    assert!(!edited
        .get(1)
        .unwrap()
        .shares_structure_with(&shared.get(1).unwrap()));
    assert_eq!(xs(&shared), [0.0, 1.0, 2.0]);
    assert_eq!(xs(&edited), [0.0, 9.0, 2.0]);

    // Replacing positions copies only what the shared sample keeps.
    let mut moved = (*shared).clone();
    moved
        .replace_positions([3.0, 4.0, 5.0].map(single_position))
        .unwrap();
    assert_eq!(xs(&moved), [3.0, 4.0, 5.0]);
    assert_eq!(xs(&shared), [0.0, 1.0, 2.0]);
    assert_eq!(
        moved.replace_positions([single_position(0.0)]),
        Err(SampleSetError::SampleCountMismatch {
            expected: 3,
            actual: 1
        })
    );
}

#[test]
fn distribution_operations_match_the_ensemble() {
    let ensemble = grouped(&[1.0, 2.0, 3.0, 4.0], &[0, 1, 0, 1], &[2, 1]);
    let mut set = SampleSet::from_ensemble(ensemble.clone());
    let distribution = set.distribution().unwrap();
    assert_eq!(
        distribution.probabilities(),
        Ok(ensemble.probabilities().unwrap())
    );
    assert_eq!(
        distribution.state_probabilities(),
        Ok(ensemble.state_probabilities().unwrap())
    );
    assert_eq!(
        distribution.effective_sample_size(),
        Ok(ensemble.effective_sample_size().unwrap())
    );
    assert_eq!(distribution.ranked_by_weight(), ensemble.ranked_by_weight());
    assert_eq!(
        distribution.ranked_by_population(),
        ensemble.ranked_by_population()
    );
    assert_eq!(
        distribution.covering_population(0.5),
        Ok(ensemble.covering_population(0.5).unwrap())
    );
    assert_eq!(
        distribution.covering_population(0.0),
        Err(SampleSetError::InvalidFraction)
    );
    let collapsed = set.collapse().unwrap().into_ensemble().unwrap();
    let expected = ensemble.collapse().unwrap();
    // The same populations, relative to the largest sample weight, 4.
    let expected_weights = weights(&expected)
        .into_iter()
        .map(|weight| weight / 4.0)
        .collect::<Vec<_>>();
    assert_eq!(weights(&collapsed), expected_weights);
    assert_eq!(representatives(&collapsed), representatives(&expected));

    let mut energetic = energy_ensemble(&[3.0, 1.0, 2.0]);
    let mut energy_set = SampleSet::from_ensemble(energetic.clone());
    let distribution = energy_set.distribution().unwrap();
    assert_eq!(
        distribution.ranked_by_energy(),
        Ok(energetic.ranked_by_energy().unwrap())
    );
    let window = Quantity::new(1.0, KILOJOULE_PER_MOLE);
    assert_eq!(
        distribution.within_energy_window(window),
        Ok(energetic.within_energy_window(window).unwrap())
    );
    assert_eq!(
        distribution.within_energy_window(Quantity::new(-1.0, KILOJOULE_PER_MOLE)),
        Err(SampleSetError::InvalidEnergyWindow)
    );
    let temperature = Quantity::new(300.0, KELVIN);
    energetic.set_temperature(Some(temperature)).unwrap();
    energy_set
        .distribution_mut()
        .unwrap()
        .set_temperature(Some(temperature))
        .unwrap();
    assert_eq!(
        energy_set.distribution().unwrap().temperature(),
        energetic.temperature()
    );
    assert_eq!(
        energy_set
            .distribution_mut()
            .unwrap()
            .set_temperature(Some(Quantity::new(0.0, KELVIN))),
        Err(SampleSetError::InvalidTemperature)
    );

    // Regrouping merges whole states with their representative's energy.
    let partition = EnsemblePartition::new(vec![0, 0, 1], vec![1, 2]).unwrap();
    energetic.regroup(partition.clone()).unwrap();
    energy_set
        .distribution_mut()
        .unwrap()
        .regroup(partition)
        .unwrap();
    assert_eq!(
        energies(&energy_set.to_ensemble().unwrap()),
        energies(&energetic)
    );
    assert_eq!(
        energy_set
            .distribution_mut()
            .unwrap()
            .regroup(EnsemblePartition::one_state(2)),
        Err(SampleSetError::PartitionMismatch {
            expected: 3,
            actual: 2
        })
    );
    set.distribution_mut().unwrap().normalize_weights().unwrap();
}

/// Fields of two ensembles compared one by one, since an ensemble is not
/// `PartialEq`.
fn assert_same_ensemble(actual: &Ensemble, expected: &Ensemble) {
    assert!(Arc::ptr_eq(
        &actual.shared_topology(),
        &expected.shared_topology()
    ));
    assert_eq!(actual.properties(), expected.properties());
    assert_eq!(actual.len(), expected.len());
    for (a, b) in actual.iter().zip(expected.iter()) {
        assert_eq!(a.payload(), b.payload());
        assert!(a.shares_structure_with(&b));
    }
    assert_eq!(actual.partition(), expected.partition());
    assert_eq!(actual.state_count(), expected.state_count());
    for (a, b) in actual.states().zip(expected.states()) {
        assert_eq!(a.payload(), b.payload());
    }
    assert_eq!(actual.temperature(), expected.temperature());
}

#[test]
fn ensembles_round_trip_sharing_every_structure() {
    let mut ensemble = grouped(&[1.0, 2.0, 3.0], &[0, 1, 0], &[2, 1]);
    ensemble.state_mut(0).unwrap().set_energy(kj(-1.5)).unwrap();
    ensemble
        .state_mut(1)
        .unwrap()
        .properties_mut()
        .insert(key("label"), PropertyValue::String("b".into()))
        .unwrap();
    ensemble
        .set_temperature(Some(Quantity::new(298.15, KELVIN)))
        .unwrap();
    ensemble
        .properties_mut()
        .insert(key("method"), PropertyValue::String("crest".into()))
        .unwrap();
    ensemble
        .get_mut(1)
        .unwrap()
        .set_occupancies(Some(vec![Some(0.25)]))
        .unwrap();

    let set = SampleSet::from_ensemble(ensemble.clone());
    assert_eq!(set.order(), None);
    assert_eq!(set.properties(), ensemble.properties());
    let distribution = set.distribution().unwrap();
    assert_eq!(distribution.members(), &Members::All);
    assert_eq!(distribution.weights(), Some(&[1.0, 2.0, 3.0][..]));
    assert_eq!(distribution.temperature(), ensemble.temperature());
    for (sample, original) in set.iter().zip(ensemble.iter()) {
        assert!(Arc::ptr_eq(
            sample.payload().structure(),
            &original.payload().structure
        ));
        // The view reports the same context as the ensemble sample's.
        assert_eq!(
            sample.as_model_view().kind(),
            original.as_model_view().kind()
        );
    }
    assert_same_ensemble(&set.to_ensemble().unwrap(), &ensemble);
    assert_same_ensemble(&set.into_ensemble().unwrap(), &ensemble);
    // An empty ensemble too.
    let empty = Ensemble::new(single_atom_topology());
    assert_same_ensemble(
        &SampleSet::from_ensemble(empty.clone())
            .to_ensemble()
            .unwrap(),
        &empty,
    );
}

#[test]
fn ensemble_conversion_refuses_what_an_ensemble_cannot_hold() {
    let ensemble = grouped(&[1.0, 2.0], &[0, 0], &[0]);
    let set = SampleSet::from_ensemble(ensemble);
    let not = |content| SampleSetError::NotConvertible {
        target: ConversionTarget::Ensemble,
        content,
    };
    assert_eq!(
        set_at(&[0.0]).to_ensemble().unwrap_err(),
        SampleSetError::Missing(Requirement::Distribution)
    );
    let mut ordered = set.clone();
    ordered.set_order(Some(Order::Step(vec![0, 1]))).unwrap();
    assert_eq!(
        ordered.to_ensemble().unwrap_err(),
        not(Unconvertible::Order)
    );
    let mut partial = set.clone();
    partial
        .push(Sample::new(single_position(1.0)), None)
        .unwrap();
    assert_eq!(
        partial.to_ensemble().unwrap_err(),
        not(Unconvertible::NonMembers)
    );
    let mut moving = set.clone();
    moving
        .get_mut(1)
        .unwrap()
        .set_velocities(Some(velocities(1.0)))
        .unwrap();
    assert_eq!(
        moving.into_ensemble().unwrap_err(),
        not(Unconvertible::Velocities { sample: 1 })
    );
    // Members without weights weigh one.
    let mut equal = set_at(&[0.0, 1.0]);
    equal
        .set_distribution(Some(distribution_of_one_state(2)))
        .unwrap();
    assert_eq!(weights(&equal.to_ensemble().unwrap()), [1.0, 1.0]);
}

fn assert_same_trajectory(actual: &Trajectory, expected: &Trajectory) {
    assert!(Arc::ptr_eq(
        &actual.shared_topology(),
        &expected.shared_topology()
    ));
    assert_eq!(actual.properties(), expected.properties());
    assert_eq!(actual.len(), expected.len());
    for (a, b) in actual.iter().zip(expected.iter()) {
        assert_eq!(a.payload(), b.payload());
    }
}

fn frame(x: f64, time: Option<f64>, step: Option<u64>) -> TrajectoryFrame {
    let mut frame = TrajectoryFrame::new(single_position(x));
    frame.set_time(time.map(ps)).unwrap();
    frame.set_step(step);
    frame
}

#[test]
fn trajectories_round_trip_with_their_times_steps_and_dynamics() {
    let topology = Arc::new(single_atom_topology());
    let mut frames = vec![
        frame(0.0, Some(0.0), Some(0)),
        frame(1.0, Some(2.0), Some(1000)),
    ];
    frames[0].set_velocities(Some(velocities(0.5))).unwrap();
    frames[1].set_forces(Some(Forces::zeros(1))).unwrap();
    frames[1] = frames[1].clone().with_cell(cell());
    frames[1]
        .properties_mut()
        .owner_mut()
        .insert(
            key("energy"),
            PropertyValue::real(-3.0, KILOJOULE_PER_MOLE).unwrap(),
        )
        .unwrap();
    let mut trajectory = Trajectory::from_items(Arc::clone(&topology), frames).unwrap();
    trajectory
        .properties_mut()
        .insert(key("engine"), PropertyValue::String("wiggle".into()))
        .unwrap();

    let set = SampleSet::from_trajectory(trajectory.clone()).unwrap();
    assert_eq!(set.order(), Some(&time_axis(vec![0.0, 2.0], vec![0, 1000])));
    assert_eq!(set.distribution(), None);
    assert_eq!(set.get(0).unwrap().velocities(), Some(&velocities(0.5)));
    for (sample, original) in set.iter().zip(trajectory.iter()) {
        assert_eq!(
            sample.as_model_view().kind(),
            original.as_model_view().kind()
        );
    }
    assert_same_trajectory(&set.to_trajectory().unwrap(), &trajectory);
    assert_same_trajectory(&set.into_trajectory().unwrap(), &trajectory);

    // Frames without times or steps give an axis without columns.
    let untimed = Trajectory::from_items(Arc::clone(&topology), [frame(0.0, None, None)]).unwrap();
    let set = SampleSet::from_trajectory(untimed.clone()).unwrap();
    assert_eq!(set.order(), Some(&Order::Time(TimeAxis::default())));
    assert_same_trajectory(&set.to_trajectory().unwrap(), &untimed);
    let empty = Trajectory::new(Arc::clone(&topology));
    assert_same_trajectory(
        &SampleSet::from_trajectory(empty.clone())
            .unwrap()
            .to_trajectory()
            .unwrap(),
        &empty,
    );

    // A column covers every frame or none.
    let partly = |frames: Vec<TrajectoryFrame>| {
        SampleSet::from_trajectory(Trajectory::from_items(Arc::clone(&topology), frames).unwrap())
            .unwrap_err()
    };
    assert_eq!(
        partly(vec![frame(0.0, Some(0.0), None), frame(1.0, None, None)]),
        SampleSetError::MissingTime { sample: 1 }
    );
    assert_eq!(
        partly(vec![frame(0.0, None, None), frame(1.0, Some(1.0), None)]),
        SampleSetError::UnexpectedTime { sample: 1 }
    );
    assert_eq!(
        partly(vec![frame(0.0, None, Some(0)), frame(1.0, None, None)]),
        SampleSetError::MissingStep { sample: 1 }
    );
    assert_eq!(
        partly(vec![frame(0.0, None, None), frame(1.0, None, Some(1))]),
        SampleSetError::UnexpectedStep { sample: 1 }
    );
}

#[test]
fn trajectory_conversion_refuses_what_a_trajectory_cannot_hold() {
    let not = |content| SampleSetError::NotConvertible {
        target: ConversionTarget::Trajectory,
        content,
    };
    let mut set = set_at(&[0.0, 1.0]);
    assert_eq!(
        set.to_trajectory().unwrap_err(),
        SampleSetError::Missing(Requirement::TimeOrder)
    );
    set.set_order(Some(Order::Step(vec![0, 1]))).unwrap();
    assert_eq!(
        set.to_trajectory().unwrap_err(),
        SampleSetError::Missing(Requirement::TimeOrder)
    );
    set.set_order(Some(time_axis(vec![0.0, 1.0], vec![0, 1])))
        .unwrap();
    let mut clustered = set.clone();
    clustered
        .set_distribution(Some(distribution_of_one_state(2)))
        .unwrap();
    assert_eq!(
        clustered.to_trajectory().unwrap_err(),
        not(Unconvertible::Distribution)
    );
    set.get_mut(1)
        .unwrap()
        .set_b_factors(Some(Quantity::new(vec![Some(1.0)], SQUARE_ANGSTROM)))
        .unwrap();
    assert_eq!(
        set.to_trajectory().unwrap_err(),
        not(Unconvertible::BFactors { sample: 1 })
    );
}

#[test]
fn series_round_trip_on_steps_and_coordinates() {
    let topology = Arc::new(single_atom_topology());
    let mut history = Series::new(Arc::clone(&topology), SeriesAxis::Step);
    for (step, x) in [(0, 0.0), (3, 1.0)] {
        history
            .push(SeriesItem::at_step(single_position(x), step))
            .unwrap();
    }
    history
        .properties_mut()
        .insert(key("method"), PropertyValue::String("lbfgs".into()))
        .unwrap();
    let set = SampleSet::from_series(history.clone());
    assert_eq!(set.order(), Some(&Order::Step(vec![0, 3])));
    let round = set.to_series().unwrap();
    assert_eq!(round.axis(), SeriesAxis::Step);
    assert_eq!(round.properties(), history.properties());
    for (a, b) in round.iter().zip(history.iter()) {
        assert_eq!(a.payload(), b.payload());
    }
    for (sample, original) in set.iter().zip(history.iter()) {
        assert_eq!(
            sample.as_model_view().kind(),
            original.as_model_view().kind()
        );
    }

    let mut scan = Series::new(Arc::clone(&topology), SeriesAxis::Coordinate(DEGREE));
    for (angle, x) in [(60.0, 0.0), (-60.0, 1.0)] {
        let item = SeriesItem::at_coordinate(single_position(x), Quantity::new(angle, DEGREE))
            .unwrap()
            .with_cell(cell());
        scan.push(item).unwrap();
    }
    let set = SampleSet::from_series(scan.clone());
    assert_eq!(
        set.order(),
        Some(&Order::Coordinate {
            unit: DEGREE,
            values: vec![60.0, -60.0],
        })
    );
    let round = set.into_series().unwrap();
    assert_eq!(round.axis(), SeriesAxis::Coordinate(DEGREE));
    for (a, b) in round.iter().zip(scan.iter()) {
        assert_eq!(a.payload(), b.payload());
    }

    let not = |content| SampleSetError::NotConvertible {
        target: ConversionTarget::Series,
        content,
    };
    let mut set = set_at(&[0.0]);
    assert_eq!(
        set.to_series().unwrap_err(),
        SampleSetError::Missing(Requirement::AnyOrder)
    );
    set.set_order(Some(Order::Time(TimeAxis::default())))
        .unwrap();
    assert_eq!(set.to_series().unwrap_err(), not(Unconvertible::Order));
    set.set_order(Some(Order::Step(vec![0]))).unwrap();
    set.get_mut(0)
        .unwrap()
        .set_forces(Some(Forces::zeros(1)))
        .unwrap();
    assert_eq!(
        set.to_series().unwrap_err(),
        not(Unconvertible::Forces { sample: 0 })
    );
}

#[test]
fn sample_set_views_write_as_the_old_views_do() {
    let mut ensemble = energy_ensemble(&[2.5, 1.0]);
    ensemble
        .push_sample(EnsembleSample::new(single_position(1.0), 3.0).unwrap(), 0)
        .unwrap();
    ensemble
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    let set = SampleSet::from_ensemble(ensemble.clone());
    let options = crate::sdf::SdfWriteOptions::default();
    let expected = crate::sdf::write(ensemble.iter(), options).unwrap();
    let views = set.iter().collect::<Vec<_>>();
    let written = crate::sdf::write(
        views.iter().map(crate::sdf::SdfRecordSource::model),
        options,
    )
    .unwrap();
    assert_eq!(written, expected);

    let mut scan = Series::new(single_atom_topology(), SeriesAxis::Coordinate(DEGREE));
    scan.push(
        SeriesItem::at_coordinate(single_position(0.0), Quantity::new(30.0, DEGREE)).unwrap(),
    )
    .unwrap();
    let set = SampleSet::from_series(scan.clone());
    let views = set.iter().collect::<Vec<_>>();
    assert_eq!(
        crate::sdf::write(
            views.iter().map(crate::sdf::SdfRecordSource::model),
            options
        )
        .unwrap(),
        crate::sdf::write(scan.iter(), options).unwrap()
    );
}

#[test]
fn subsets_keep_the_order_and_the_distribution() {
    let (model, atom, _) = model_fixture();
    let mut set = SampleSet::from_models([model.clone(), model]).unwrap();
    set.set_order(Some(Order::Step(vec![5, 6]))).unwrap();
    let state = EnsembleState::new()
        .with_energy(Quantity::new(1.0, KILOJOULE_PER_MOLE))
        .unwrap();
    let mut state = state;
    state
        .properties_mut()
        .insert(key("label"), PropertyValue::String("a".into()))
        .unwrap();
    set.set_distribution(Some(
        distribution(Members::Subset(vec![1]), &[0], &[0])
            .with_states(vec![state])
            .unwrap(),
    ))
    .unwrap();
    set.properties_mut()
        .insert(key("method"), PropertyValue::String("md".into()))
        .unwrap();
    set.get_mut(0)
        .unwrap()
        .set_velocities(Some(Velocities::zeros(2)))
        .unwrap();
    let selection =
        crate::topology::AtomSelection::from_atoms(set.shared_topology(), [atom]).unwrap();
    let subset = set.subset(&selection).unwrap();
    assert_eq!(subset.topology().atom_count(), 1);
    assert!(subset.properties().is_empty());
    assert_eq!(subset.order(), set.order());
    assert_eq!(subset.get(0).unwrap().velocities().unwrap().len(), 1);
    let state = subset.state(0).unwrap();
    assert_eq!(state.energy(), Some(Quantity::new(1.0, KILOJOULE_PER_MOLE)));
    assert!(state.properties().is_empty());
    assert_eq!(
        subset.distribution().unwrap().members(),
        &Members::Subset(vec![1])
    );

    let mut perceived = set.clone();
    perceived.perceive().unwrap();
    assert!(perceived.topology().shares_layout(set.topology()));
    assert!(shares_every_sample(&perceived, &set));
    assert_eq!(
        SampleSet::from_models(std::iter::empty()).unwrap_err(),
        SampleSetError::EmptySource
    );
}

/// The parts of two sets compared one by one, since a set is not
/// `PartialEq`; samples must share their structures.
fn assert_same_set(actual: &SampleSet, expected: &SampleSet) {
    assert!(shares_every_sample(actual, expected));
    assert_eq!(actual.order(), expected.order());
    assert_eq!(actual.distribution(), expected.distribution());
    assert_eq!(actual.properties(), expected.properties());
}

#[test]
fn replacing_a_sample_keeps_its_place_and_membership() {
    let mut set = set_at(&[0.0, 1.0]);
    set.set_order(Some(Order::Step(vec![3, 4]))).unwrap();
    set.set_distribution(Some(
        distribution(Members::Subset(vec![1]), &[0], &[0])
            .with_weights(vec![2.0])
            .unwrap(),
    ))
    .unwrap();
    let previous = set.replace(1, Sample::new(single_position(5.0))).unwrap();
    assert_eq!(previous, Sample::new(single_position(1.0)));
    assert_eq!(xs(&set), [0.0, 5.0]);
    assert_eq!(set.get(1).unwrap().place(), Some(Place::Step(4)));
    assert_eq!(set.get(1).unwrap().weight(), Some(2.0));
    assert!(set.get(1).unwrap().is_representative());

    let before = set.clone();
    // Two atoms do not bind to a one-atom topology.
    assert!(matches!(
        set.replace(0, Sample::new(Positions::zeros(2))),
        Err(SampleSetError::Structure(_))
    ));
    assert_eq!(
        set.replace(2, Sample::new(single_position(0.0))),
        Err(SampleSetError::IndexOutOfRange { index: 2, len: 2 })
    );
    assert_same_set(&set, &before);
}

#[test]
fn the_distribution_guard_keeps_the_members() {
    let mut set = set_at(&[0.0, 1.0]);
    set.set_distribution(Some(distribution_of_one_state(2)))
        .unwrap();
    // The replacement the guard cannot make is refused when installed.
    let foreign =
        Distribution::new(Members::Subset(vec![7]), EnsemblePartition::one_state(1)).unwrap();
    assert_eq!(
        set.set_distribution(Some(foreign)),
        Err(SampleSetError::IndexOutOfRange { index: 7, len: 2 })
    );
    let mut distribution = set.distribution_mut().unwrap();
    assert_eq!(distribution.member_count(), 2);
    distribution.set_weight(1, 3.0).unwrap();
    distribution
        .regroup(EnsemblePartition::state_per_sample(2))
        .unwrap();
    distribution
        .state_mut(1)
        .unwrap()
        .set_energy(kj(1.0))
        .unwrap();
    assert_eq!(set.distribution().unwrap().members(), &Members::All);
    assert_eq!(set.state(1).unwrap().representative_sample().index(), 1);
    assert_eq!(set.state(1).unwrap().energy(), kj(1.0));
}

#[test]
fn weighted_walks_join_members_and_samples() {
    let mut set = set_at(&[0.0, 1.0, 2.0]);
    // Without a distribution every sample counts once.
    assert_eq!(
        set.weighted_samples().collect::<Vec<_>>(),
        [(0, 1.0), (1, 1.0), (2, 1.0)]
    );
    set.set_distribution(Some(distribution(
        Members::Subset(vec![0, 2]),
        &[0, 0],
        &[0],
    )))
    .unwrap();
    let walk = set.weighted_samples();
    assert_eq!(walk.len(), 2);
    assert_eq!(walk.collect::<Vec<_>>(), [(0, 1.0), (2, 1.0)]);
    set.distribution_mut().unwrap().set_weight(2, 4.0).unwrap();
    let distribution = set.distribution().unwrap();
    assert_eq!(
        distribution.weighted_members().collect::<Vec<_>>(),
        [(0, 1.0), (2, 4.0)]
    );
    assert_eq!(
        set.weighted_samples().collect::<Vec<_>>(),
        distribution.weighted_members().collect::<Vec<_>>()
    );
    // Member positions index the weights and probabilities.
    let probabilities = distribution.probabilities().unwrap();
    for ((sample, weight), (position, probability)) in distribution
        .weighted_members()
        .zip(probabilities.iter().enumerate())
    {
        assert_eq!(distribution.member_position(sample), Some(position));
        assert_eq!(*probability, weight / 5.0);
    }
    assert_eq!(distribution.ranked_by_weight(), [2, 0]);
}

#[test]
fn duplicates_in_a_selection_become_separate_samples() {
    let mut set = set_at(&[0.0, 1.0]);
    set.set_order(Some(Order::Step(vec![5, 6]))).unwrap();
    set.set_distribution(Some(
        distribution(Members::All, &[0, 1], &[0, 1])
            .with_weights(vec![1.0, 2.0])
            .unwrap()
            .with_states(vec![
                EnsembleState::new()
                    .with_energy(Quantity::new(1.0, KILOJOULE_PER_MOLE))
                    .unwrap(),
                EnsembleState::new()
                    .with_energy(Quantity::new(2.0, KILOJOULE_PER_MOLE))
                    .unwrap(),
            ])
            .unwrap(),
    ))
    .unwrap();
    let twice = set.select([1, 1]).unwrap();
    assert_eq!(xs(&twice), [1.0, 1.0]);
    assert_eq!(twice.order(), Some(&Order::Step(vec![6, 6])));
    let distribution = twice.distribution().unwrap();
    assert_eq!(distribution.members(), &Members::All);
    assert_eq!(distribution.weights(), Some(&[2.0, 2.0][..]));
    assert_eq!(distribution.partition().assignments(), [0, 0]);
    assert_eq!(distribution.partition().representatives(), [0]);
    // Every member of state 1 is selected, so it keeps its energy.
    assert_eq!(twice.state(0).unwrap().energy(), kj(2.0));
    assert!(twice
        .get(0)
        .unwrap()
        .shares_structure_with(&twice.get(1).unwrap()));

    let states = set.select_states([0, 0]).unwrap();
    assert_eq!(xs(&states), [0.0, 0.0]);
    assert_eq!(states.order(), Some(&Order::Step(vec![5, 5])));
    let distribution = states.distribution().unwrap();
    assert_eq!(distribution.state_count(), 2);
    assert_eq!(distribution.partition().assignments(), [0, 1]);
    assert_eq!(distribution.partition().representatives(), [0, 1]);
    assert_eq!(states.state(1).unwrap().energy(), kj(1.0));
}

#[test]
fn collapse_keeps_population_shares_of_any_members() {
    // Members 0 and 2 of three samples, one state, without weights.
    let mut set = set_at(&[0.0, 1.0, 2.0]);
    set.set_order(Some(Order::Step(vec![0, 1, 2]))).unwrap();
    set.set_distribution(Some(distribution(
        Members::Subset(vec![0, 2]),
        &[0, 0],
        &[1],
    )))
    .unwrap();
    let collapsed = set.collapse().unwrap();
    assert_eq!(xs(&collapsed), [2.0]);
    assert_eq!(collapsed.order(), None);
    assert_eq!(
        collapsed.distribution().unwrap().weights(),
        Some(&[2.0][..])
    );
    assert!(collapsed
        .get(0)
        .unwrap()
        .shares_structure_with(&set.get(2).unwrap()));

    // Weights whose sum overflows keep their ratios.
    let mut heavy = set_at(&[0.0, 1.0, 2.0]);
    heavy
        .set_distribution(Some(
            distribution(Members::All, &[0, 0, 1], &[0, 2])
                .with_weights(vec![f64::MAX, f64::MAX, f64::MAX])
                .unwrap(),
        ))
        .unwrap();
    assert_eq!(heavy.state(0).unwrap().population(), 2.0);
    assert_eq!(
        heavy.distribution().unwrap().state_probabilities().unwrap(),
        [2.0 / 3.0, 1.0 / 3.0]
    );
    let collapsed = heavy.collapse().unwrap();
    assert_eq!(
        collapsed.distribution().unwrap().weights(),
        Some(&[2.0, 1.0][..])
    );
}

#[test]
fn members_pushed_into_a_distribution_of_every_sample() {
    let mut set = set_at(&[0.0]);
    set.set_order(Some(Order::Step(vec![0]))).unwrap();
    set.set_distribution(Some(distribution_of_one_state(1)))
        .unwrap();
    assert_eq!(
        set.push_member(
            Sample::new(single_position(1.0)),
            Some(Place::Step(1)),
            0,
            None
        ),
        Ok(1)
    );
    assert_eq!(
        set.push_state(
            Sample::new(single_position(2.0)),
            Some(Place::Step(2)),
            EnsembleState::new(),
            None,
        ),
        Ok(1)
    );
    let distribution = set.distribution().unwrap();
    assert_eq!(distribution.members(), &Members::All);
    assert_eq!(distribution.weights(), None);
    assert_eq!(distribution.partition().assignments(), [0, 0, 1]);
    assert_eq!(distribution.partition().representatives(), [0, 2]);
    assert_eq!(set.order(), Some(&Order::Step(vec![0, 1, 2])));

    // A sample that does not bind, checked last, changes nothing.
    let before = set.clone();
    let wide = || Sample::new(Positions::zeros(2));
    assert!(matches!(
        set.push_member(wide(), Some(Place::Step(3)), 1, None),
        Err(SampleSetError::Structure(_))
    ));
    assert!(matches!(
        set.push_state(wide(), Some(Place::Step(3)), EnsembleState::new(), None),
        Err(SampleSetError::Structure(_))
    ));
    assert!(matches!(
        set.push(wide(), Some(Place::Step(3))),
        Err(SampleSetError::Structure(_))
    ));
    assert_same_set(&set, &before);
    // Without a distribution there is nothing to join.
    let mut plain = set_at(&[0.0]);
    assert_eq!(
        plain.push_member(Sample::new(single_position(1.0)), None, 0, None),
        Err(SampleSetError::Missing(Requirement::Distribution))
    );
    assert_eq!(plain.len(), 1);
}

#[test]
fn places_are_set_only_on_an_order_and_inside_the_set() {
    let mut set = set_at(&[0.0]);
    assert_eq!(set.set_place(0, None), Ok(()));
    assert_eq!(
        set.set_place(0, Some(Place::Step(1))),
        Err(SampleSetError::PlaceMismatch { sample: 0 })
    );
    assert_eq!(
        set.set_place(1, None),
        Err(SampleSetError::IndexOutOfRange { index: 1, len: 1 })
    );
    set.set_order(Some(Order::Step(vec![0]))).unwrap();
    assert_eq!(
        set.set_place(3, Some(Place::Step(1))),
        Err(SampleSetError::IndexOutOfRange { index: 3, len: 1 })
    );
    assert_eq!(
        set.set_place(0, None),
        Err(SampleSetError::PlaceMismatch { sample: 0 })
    );
    assert_eq!(set.order(), Some(&Order::Step(vec![0])));
}

#[test]
fn empty_selections_and_distributions() {
    let mut set = set_at(&[0.0, 1.0]);
    set.set_order(Some(Order::Step(vec![0, 1]))).unwrap();
    set.set_distribution(Some(distribution(Members::Subset(vec![0]), &[0], &[0])))
        .unwrap();
    let none = set.select([]).unwrap();
    assert!(none.is_empty());
    assert_eq!(none.order(), Some(&Order::Step(Vec::new())));
    assert_eq!(none.distribution().unwrap().member_count(), 0);
    assert_eq!(none.states().len(), 0);
    let none = set.select_states([]).unwrap();
    assert!(none.is_empty());
    assert_eq!(none.distribution().unwrap().state_count(), 0);
    // A distribution without members collapses to an empty set.
    let collapsed = none.collapse().unwrap();
    assert!(collapsed.is_empty());
    assert_eq!(collapsed.distribution().unwrap().state_count(), 0);
    assert_eq!(
        none.distribution().unwrap().probabilities(),
        Err(SampleSetError::EmptyDistribution)
    );
    // Removing the last non-member leaves a distribution of every sample.
    set.remove(1).unwrap();
    assert_eq!(set.distribution().unwrap().members(), &Members::All);
    assert_eq!(set.order(), Some(&Order::Step(vec![0])));
    // Removing the last member leaves an empty distribution.
    set.remove(0).unwrap();
    assert!(set.is_empty());
    assert_eq!(set.distribution().unwrap().state_count(), 0);
    assert_eq!(set.distribution().unwrap().members(), &Members::All);
}

#[test]
fn regrouping_a_split_state_drops_its_energy() {
    let mut set = set_at(&[0.0, 1.0, 2.0]);
    set.set_distribution(Some(
        distribution(Members::All, &[0, 0, 1], &[0, 2])
            .with_states(vec![
                EnsembleState::new()
                    .with_energy(Quantity::new(1.0, KILOJOULE_PER_MOLE))
                    .unwrap(),
                EnsembleState::new()
                    .with_energy(Quantity::new(2.0, KILOJOULE_PER_MOLE))
                    .unwrap(),
            ])
            .unwrap(),
    ))
    .unwrap();
    let partition = EnsemblePartition::new(vec![0, 1, 2], vec![0, 1, 2]).unwrap();
    set.distribution_mut().unwrap().regroup(partition).unwrap();
    assert_eq!(set.states().len(), 3);
    assert_eq!(set.state(0).unwrap().energy(), None);
    assert_eq!(set.state(1).unwrap().energy(), None);
    // State 1 was whole and keeps its energy.
    assert_eq!(set.state(2).unwrap().energy(), kj(2.0));
}

#[test]
fn constructors_refuse_mismatched_samples_and_layouts() {
    assert!(matches!(
        SampleSet::from_samples(
            single_atom_topology(),
            [
                Sample::new(single_position(0.0)),
                Sample::new(Positions::zeros(2))
            ],
        ),
        Err(SampleSetError::Structure(StructureError::Array(
            DenseArrayError::CountMismatch {
                expected: 1,
                actual: 2
            }
        )))
    ));
    let (model, _, _) = model_fixture();
    let other = Model::new(Arc::new(single_atom_topology()), single_position(0.0)).unwrap();
    assert_eq!(
        SampleSet::from_models([model, other]).unwrap_err(),
        SampleSetError::TopologyMismatch
    );
}

#[test]
fn every_error_reads_as_a_sentence() {
    use SampleSetError as E;
    let requirements = [
        (Requirement::TimeOrder, "needs a time order"),
        (Requirement::AnyOrder, "needs an order"),
        (Requirement::Times, "needs times"),
        (Requirement::Distribution, "needs a distribution"),
        (Requirement::Weights, "needs weights"),
        (
            Requirement::StateEnergies,
            "needs an energy for every state",
        ),
        (Requirement::Temperature, "needs a temperature"),
        (Requirement::Cell, "needs a periodic cell on every sample"),
        (Requirement::InMemory, "needs its samples in memory"),
    ];
    for (requirement, text) in requirements {
        assert_eq!(requirement.to_string(), text);
        assert_eq!(
            E::Missing(requirement).to_string(),
            format!("The sample set {text}.")
        );
    }
    let unit = Quantity::new(1.0, ANGSTROM)
        .value_in(PICOSECOND)
        .unwrap_err();
    let not = |target, content| E::NotConvertible { target, content };
    let errors = [
        (E::EmptySource, "There are no models to collect."),
        (E::EmptyDistribution, "The distribution has no members."),
        (
            E::IndexOutOfRange { index: 3, len: 1 },
            "Sample 3 is outside the sample set of 1 sample.",
        ),
        (
            E::StateOutOfRange {
                state: 2,
                states: 2,
            },
            "State 2 is outside the distribution of 2 states.",
        ),
        (
            E::NotMember { sample: 4 },
            "Sample 4 is outside the distribution.",
        ),
        (
            E::TopologyMismatch,
            "The model has a different topology layout from the sample set.",
        ),
        (
            E::SampleCountMismatch {
                expected: 3,
                actual: 1,
            },
            "The sample set has 3 samples, but 1 replacement was given.",
        ),
        (
            E::OrderLengthMismatch {
                expected: 2,
                actual: 1,
            },
            "The order has 1 value, but the sample set has 2 samples.",
        ),
        (
            E::TimeAxisMismatch { times: 1, steps: 2 },
            "The time axis has 1 time and 2 steps.",
        ),
        (
            E::NonFiniteOrderValue { sample: 1 },
            "The order value of sample 1 is not finite.",
        ),
        (
            E::PlaceMismatch { sample: 1 },
            "The place of sample 1 is not on the order of the sample set.",
        ),
        (
            E::MissingTime { sample: 1 },
            "Sample 1 has no time, but the other samples have one.",
        ),
        (
            E::UnexpectedTime { sample: 1 },
            "Sample 1 has a time, but the other samples have none.",
        ),
        (
            E::MissingStep { sample: 1 },
            "Sample 1 has no step, but the other samples have one.",
        ),
        (
            E::UnexpectedStep { sample: 1 },
            "Sample 1 has a step, but the other samples have none.",
        ),
        (
            E::NonMonotonicTime { sample: 2 },
            "Sample 2 is earlier than the sample before it.",
        ),
        (
            E::UnsortedMembers,
            "The members of the distribution are not in increasing order.",
        ),
        (
            E::CoverageMismatch {
                covered: 2,
                samples: 3,
            },
            "The distribution covers 2 samples, but the sample set has 3.",
        ),
        (
            E::PartitionMismatch {
                expected: 2,
                actual: 3,
            },
            "The partition covers 3 samples, but the distribution has 2 members.",
        ),
        (
            E::WeightCountMismatch {
                expected: 2,
                actual: 1,
            },
            "The distribution has 2 members, but 1 weight was given.",
        ),
        (
            E::StateCountMismatch {
                expected: 1,
                actual: 2,
            },
            "The partition has 1 state, but 2 states were given.",
        ),
        (
            E::InvalidWeight { sample: 0 },
            "The weight of sample 0 is not finite and positive.",
        ),
        (
            E::WeightRequired,
            "The new member has no weight, but the other members have weights.",
        ),
        (
            E::RepresentativeRemoval {
                sample: 0,
                state: 1,
            },
            "Sample 0 represents state 1. Remove the state's other samples first.",
        ),
        (
            E::InvalidTemperature,
            "The temperature is not finite and positive.",
        ),
        (
            E::InvalidEnergyWindow,
            "The energy window is not finite and nonnegative.",
        ),
        (
            E::InvalidFraction,
            "The population fraction is not greater than 0 and at most 1.",
        ),
        (
            not(ConversionTarget::Ensemble, Unconvertible::Order),
            "An ensemble cannot hold the order of the sample set.",
        ),
        (
            not(ConversionTarget::Trajectory, Unconvertible::Distribution),
            "A trajectory cannot hold the distribution of the sample set.",
        ),
        (
            not(ConversionTarget::Ensemble, Unconvertible::NonMembers),
            "An ensemble cannot hold samples outside the distribution.",
        ),
        (
            not(
                ConversionTarget::Series,
                Unconvertible::BFactors { sample: 2 },
            ),
            "A series cannot hold the B-factors of sample 2.",
        ),
        (
            not(
                ConversionTarget::Series,
                Unconvertible::Velocities { sample: 2 },
            ),
            "A series cannot hold the velocities of sample 2.",
        ),
        (
            not(
                ConversionTarget::Series,
                Unconvertible::Forces { sample: 2 },
            ),
            "A series cannot hold the forces of sample 2.",
        ),
        (
            not(
                ConversionTarget::Trajectory,
                Unconvertible::Occupancies { sample: 2 },
            ),
            "A trajectory cannot hold the occupancies of sample 2.",
        ),
        (
            E::Structure(StructureError::NonFiniteTime),
            "A sample is invalid.",
        ),
        (
            E::Subset(crate::topology::TopologySubsetError::EmptySelection),
            "Cannot subset the sample set.",
        ),
    ];
    for (error, text) in errors {
        assert_eq!(error.to_string(), text, "{error:?}");
    }
    assert_eq!(E::from(unit).to_string(), unit.to_string());
    // Wrapped causes are the error's source.
    use std::error::Error;
    assert!(E::Structure(StructureError::NonFiniteTime)
        .source()
        .is_some());
    assert!(E::EmptySource.source().is_none());
}

#[test]
fn failed_and_empty_sample_edits_keep_the_structure_shared() {
    let shared = Sample::new(single_position(0.0)).with_cell(cell());
    let mut sample = shared.clone();
    assert!(sample
        .set_positions(Quantity::new(vec![Point3::origin(); 2], NANOMETER))
        .is_err());
    assert!(sample
        .set_position_at(
            TopologyAtomIndex::new(1),
            Quantity::new(Point3::origin(), NANOMETER)
        )
        .is_err());
    assert!(sample
        .set_occupancy(TopologyAtomIndex::new(1), Some(0.5))
        .is_err());
    assert!(sample.set_occupancies(Some(vec![Some(f64::NAN)])).is_err());
    assert!(sample
        .set_b_factors(Some(Quantity::new(vec![Some(1.0)], NANOMETER)))
        .is_err());
    assert!(sample.set_velocities(Some(Velocities::zeros(2))).is_err());
    sample.set_velocities(None).unwrap();
    sample.set_forces(None).unwrap();
    assert_eq!(sample.take_velocities(), None);
    sample.set_cell(Some(cell()));
    assert!(sample.shares_structure_with(&shared));
    // A successful edit copies this sample alone.
    sample.set_cell(None);
    assert!(!sample.shares_structure_with(&shared));
    assert_eq!(shared.cell(), Some(&cell()));

    // The same through a set's guard.
    let mut set = SampleSet::from_samples(single_atom_topology(), [shared.clone()]).unwrap();
    let before = set.clone();
    assert!(set
        .get_mut(0)
        .unwrap()
        .set_occupancies(Some(vec![Some(1.0), Some(1.0)]))
        .is_err());
    set.get_mut(0).unwrap().set_cell(Some(cell()));
    assert!(shares_every_sample(&set, &before));
}
