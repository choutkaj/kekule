//! Structure sets: the order, the distribution, views, and copy on write.

use super::*;
use crate::units::{DEGREE, PICOSECOND};

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
    let partition = Partition::new(assignments.to_vec(), representatives.to_vec()).unwrap();
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

fn shares_every_geometry(a: &StructureSet, b: &StructureSet) -> bool {
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
        StructureSetError::Missing(Requirement::Distribution).to_string(),
        "The structure set needs a distribution."
    );
    let set = set_at(&[0.0]);
    assert_eq!(
        set.validate_monotonic_time(),
        Err(StructureSetError::Missing(Requirement::TimeOrder))
    );
    assert_eq!(
        set.select_states([0]).unwrap_err(),
        StructureSetError::Missing(Requirement::Distribution)
    );
    assert_eq!(
        set.collapse().unwrap_err(),
        StructureSetError::Missing(Requirement::Distribution)
    );
}

#[test]
fn an_order_has_one_value_per_sample_in_every_column() {
    let mut set = set_at(&[0.0, 1.0]);
    assert_eq!(
        set.set_order(Some(Order::Step(vec![0]))),
        Err(StructureSetError::OrderLengthMismatch {
            expected: 2,
            actual: 1
        })
    );
    assert_eq!(
        set.set_order(Some(Order::Coordinate {
            unit: DEGREE,
            values: vec![0.0, f64::NAN],
        })),
        Err(StructureSetError::NonFiniteOrderValue { geometry: 1 })
    );
    assert_eq!(set.order(), None);
    assert_eq!(
        TimeAxis::new(Some(Quantity::new(vec![0.0], PICOSECOND)), Some(vec![0, 1])),
        Err(StructureSetError::TimeAxisMismatch { times: 1, steps: 2 })
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

    // A new geometry states exactly the columns the axis has.
    let before = set.clone();
    let geometry = || Geometry::new(single_position(2.0));
    assert_eq!(
        set.push(geometry(), None),
        Err(StructureSetError::PlaceMismatch { geometry: 2 })
    );
    assert_eq!(
        set.push(
            geometry(),
            Some(Place::Time {
                time: None,
                step: None
            })
        ),
        Err(StructureSetError::MissingTime { geometry: 2 })
    );
    assert_eq!(
        set.push(geometry(), timed(3.0, 7)),
        Err(StructureSetError::UnexpectedStep { geometry: 2 })
    );
    assert_eq!(
        set.push(geometry(), Some(Place::Step(2))),
        Err(StructureSetError::PlaceMismatch { geometry: 2 })
    );
    assert_eq!(set.len(), before.len());
    assert_eq!(set.order(), before.order());
    let place = Some(Place::Time {
        time: Some(ps(3.0)),
        step: None,
    });
    assert_eq!(set.push(geometry(), place), Ok(2));
    assert_eq!(set.get(2).unwrap().place(), place);
    // An axis of steps alone takes a step and no time.
    let mut stepped = set_at(&[0.0]);
    stepped
        .set_order(Some(Order::Time(
            TimeAxis::new(None, Some(vec![0])).unwrap(),
        )))
        .unwrap();
    assert_eq!(
        stepped.push(
            geometry(),
            Some(Place::Time {
                time: None,
                step: None
            })
        ),
        Err(StructureSetError::MissingStep { geometry: 1 })
    );
    assert_eq!(
        stepped.push(geometry(), timed(1.0, 1)),
        Err(StructureSetError::UnexpectedTime { geometry: 1 })
    );
    assert_eq!(stepped.len(), 1);

    // A coordinate converts to the order's unit.
    let mut scan = set_at(&[0.0]);
    scan.set_order(Some(Order::Coordinate {
        unit: NANOMETER,
        values: vec![0.1],
    }))
    .unwrap();
    scan.push(
        Geometry::new(single_position(1.0)),
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
            Geometry::new(single_position(1.0)),
            Some(Place::Coordinate(Quantity::new(2.0, PICOSECOND))),
        ),
        Err(StructureSetError::Unit(_))
    ));
    // An unordered set takes no place.
    let mut unordered = set_at(&[0.0]);
    assert_eq!(
        unordered.push(Geometry::new(single_position(1.0)), Some(Place::Step(1))),
        Err(StructureSetError::PlaceMismatch { geometry: 1 })
    );
    assert_eq!(unordered.get(0).unwrap().place(), None);
}

#[test]
fn order_values_belong_to_the_set_not_to_its_samples() {
    let mut set = set_at(&[0.0, 1.0, 2.0]);
    set.set_order(Some(time_axis(vec![0.0, 2.0, 4.0], vec![0, 10, 20])))
        .unwrap();
    // Selecting keeps each geometry's value, in the requested order.
    let reversed = set.select([2, 0]).unwrap();
    assert_eq!(xs(&reversed), [2.0, 0.0]);
    assert_eq!(reversed.get(0).unwrap().place(), timed(4.0, 20));
    assert_eq!(reversed.get(1).unwrap().place(), timed(0.0, 0));
    assert_eq!(
        reversed.validate_monotonic_time(),
        Err(StructureSetError::NonMonotonicTime { geometry: 1 })
    );
    assert_eq!(set.validate_monotonic_time(), Ok(()));
    // A removed geometry claims no place; the later geometries keep theirs.
    let removed = set.remove(1).unwrap();
    assert_eq!(removed, Geometry::new(single_position(1.0)));
    assert_eq!(set.get(1).unwrap().place(), timed(4.0, 20));
    set.set_place(1, timed(3.0, 30)).unwrap();
    assert_eq!(set.order(), Some(&time_axis(vec![0.0, 3.0], vec![0, 30])));
    assert_eq!(
        set.set_place(1, Some(Place::Step(1))),
        Err(StructureSetError::PlaceMismatch { geometry: 1 })
    );
    // A time order without times cannot be checked for monotonic time.
    let mut untimed = set_at(&[0.0]);
    untimed
        .set_order(Some(Order::Time(TimeAxis::default())))
        .unwrap();
    assert_eq!(
        untimed.validate_monotonic_time(),
        Err(StructureSetError::Missing(Requirement::Times))
    );
    untimed
        .push(Geometry::new(single_position(1.0)), None)
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
    let mut geometry = set.get_mut(0).unwrap();
    geometry.set_occupancies(Some(vec![Some(0.5)])).unwrap();
    geometry.set_velocities(Some(velocities(1.0))).unwrap();
    assert_eq!(
        geometry.set_forces(Some(Forces::zeros(2))),
        Err(GeometryError::Array(DenseArrayError::CountMismatch {
            expected: 1,
            actual: 2
        }))
    );
    let view = set.get(0).unwrap();
    assert_eq!(view.occupancy(ATOM0), Some(0.5));
    let StructureKind::InSet(context) = view.as_structure_view().kind() else {
        panic!("a geometry's view is a geometry");
    };
    assert_eq!(context.place(), Some(Place::Step(0)));
    assert_eq!(context.membership(), None);
    assert_eq!(context.velocities(), Some(&velocities(1.0)));
    assert_eq!(context.forces(), None);
    assert_eq!(view.as_structure_view().occupancy(ATOM0), Some(0.5));
    // A structure keeps the structure and drops the geometry's context.
    let structure = view.to_structure();
    assert_eq!(
        structure.as_structure_view().kind(),
        StructureKind::Standalone
    );
    assert_eq!(structure.occupancy(ATOM0), Some(0.5));
    // A rigid motion rotates the velocities with the positions.
    use crate::math::{Matrix3, RigidTransform};
    let turn = Matrix3::from_columns(
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    );
    let transform = RigidTransform::from_canonical(turn, Vector3::zero()).unwrap();
    let mut geometry = set.get_mut(0).unwrap();
    geometry.apply_rigid_transform(&transform).unwrap();
    assert_eq!(
        geometry.velocities().unwrap().values().value(),
        &[Vector3::new(0.0, 1.0, 0.0)]
    );
    assert!(set.get_mut(0).unwrap().take_velocities().is_some());
    assert_eq!(set.get(0).unwrap().velocities(), None);
}

#[test]
fn a_distribution_covers_all_samples_or_a_sorted_subset() {
    assert_eq!(
        Distribution::new(Members::Subset(vec![2, 1]), Partition::one_state(2)),
        Err(StructureSetError::UnsortedMembers)
    );
    assert_eq!(
        Distribution::new(Members::Subset(vec![1, 2]), Partition::one_state(3)),
        Err(StructureSetError::PartitionMismatch {
            expected: 2,
            actual: 3
        })
    );
    let mut set = set_at(&[0.0, 1.0, 2.0]);
    assert_eq!(
        set.set_distribution(Some(distribution(Members::All, &[0, 0], &[0]))),
        Err(StructureSetError::CoverageMismatch {
            covered: 2,
            geometries: 3
        })
    );
    assert_eq!(
        set.set_distribution(Some(distribution(
            Members::Subset(vec![1, 3]),
            &[0, 0],
            &[0]
        ))),
        Err(StructureSetError::IndexOutOfRange { index: 3, len: 3 })
    );
    assert_eq!(set.distribution(), None);
    // A subset listing every geometry is all of them.
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
    assert_eq!(set.state(1).unwrap().geometries().collect::<Vec<_>>(), [2]);
    assert_eq!(set.states().len(), 2);
    assert!(set.state(2).is_none());
    // The view's context reports membership only for members.
    let StructureKind::InSet(context) = set.get(0).unwrap().as_structure_view().kind() else {
        panic!("a geometry's view is a geometry");
    };
    assert_eq!(context.membership(), None);
    let StructureKind::InSet(context) = set.get(2).unwrap().as_structure_view().kind() else {
        panic!("a geometry's view is a geometry");
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
        Err(StructureSetError::WeightCountMismatch {
            expected: 2,
            actual: 1
        })
    );
    assert_eq!(
        distribution.set_weights(Some(vec![1.0, -1.0])),
        Err(StructureSetError::InvalidWeight { geometry: 2 })
    );
    assert_eq!(
        distribution.set_weight(1, 2.0),
        Err(StructureSetError::NotMember { geometry: 1 })
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
        set.push_member(Geometry::new(single_position(3.0)), None, 0, None),
        Err(StructureSetError::WeightRequired)
    );
    assert_eq!(
        set.push_member(Geometry::new(single_position(3.0)), None, 1, Some(1.0)),
        Err(StructureSetError::StateOutOfRange {
            state: 1,
            states: 1
        })
    );
    assert_eq!(set.len(), 3);
    assert_eq!(
        set.push_member(Geometry::new(single_position(3.0)), None, 0, Some(0.5)),
        Ok(3)
    );
    let state = State::new()
        .with_energy(Quantity::new(1.0, KILOJOULE_PER_MOLE))
        .unwrap();
    assert_eq!(
        set.push_state(Geometry::new(single_position(4.0)), None, state, Some(1.0)),
        Ok(1)
    );
    // A non-member keeps the distribution as it was.
    assert_eq!(set.push(Geometry::new(single_position(5.0)), None), Ok(5));
    let distribution = set.distribution().unwrap();
    assert_eq!(distribution.members(), &Members::Subset(vec![0, 2, 3, 4]));
    assert_eq!(distribution.weights(), Some(&[0.25, 0.75, 0.5, 1.0][..]));
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
    // Pushing a non-member into a distribution of every geometry lists the
    // others.
    let mut all = set_at(&[0.0]);
    all.set_distribution(Some(distribution_of_one_state(1)))
        .unwrap();
    all.push(Geometry::new(single_position(1.0)), None).unwrap();
    assert_eq!(
        all.distribution().unwrap().members(),
        &Members::Subset(vec![0])
    );
}

fn distribution_of_one_state(members: usize) -> Distribution {
    Distribution::new(Members::All, Partition::one_state(members)).unwrap()
}

#[test]
fn selecting_and_removing_samples_remaps_the_members() {
    // Five frames; frames 1-4 are an ensemble of two states.
    let mut set = set_at(&[0.0, 1.0, 2.0, 3.0, 4.0]);
    set.set_order(Some(Order::Step(vec![0, 1, 2, 3, 4])))
        .unwrap();
    let energies = [1.0, 2.0].map(|energy| {
        State::new()
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
    // States in state order: state 0 (geometry 3, a new representative),
    // then state 1 (geometry 4, its representative).
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
        StructureSetError::StateOutOfRange {
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
        StructureSetError::RepresentativeRemoval {
            geometry: 0,
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
        .regroup(Partition::state_per_sample(3))
        .unwrap();
    distribution.set_weights(Some(vec![1.0, 2.0, 3.0])).unwrap();
    assert!(!Arc::ptr_eq(&regrouped, &shared));
    assert!(shares_every_geometry(&regrouped, &shared));
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
    assert!(shares_every_geometry(&collapsed, &regrouped));
    let states = regrouped.select_states([1]).unwrap();
    assert!(Arc::ptr_eq(
        states.get(0).unwrap().payload().structure(),
        shared.get(1).unwrap().payload().structure()
    ));

    // Editing one geometry of a clone copies that geometry alone.
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

    // Replacing positions copies only what the shared geometry keeps.
    let mut moved = (*shared).clone();
    moved
        .replace_positions([3.0, 4.0, 5.0].map(single_position))
        .unwrap();
    assert_eq!(xs(&moved), [3.0, 4.0, 5.0]);
    assert_eq!(xs(&shared), [0.0, 1.0, 2.0]);
    assert_eq!(
        moved.replace_positions([single_position(0.0)]),
        Err(StructureSetError::GeometryCountMismatch {
            expected: 3,
            actual: 1
        })
    );
}

#[test]
fn distribution_operations_keep_their_results() {
    use crate::tests::debug_fingerprint as fingerprint;
    let mut set = grouped(&[1.0, 2.0, 3.0, 4.0], &[0, 1, 0, 1], &[2, 1]);
    assert_eq!(fingerprint(&set), 0x9661_9a29_dfa3_cb04);
    let distribution = set.distribution().unwrap();
    assert_eq!(
        fingerprint(&distribution.probabilities()),
        0x5264_938e_9878_2088
    );
    assert_eq!(
        fingerprint(&distribution.state_probabilities()),
        0xeeff_7097_2abe_e8d6
    );
    assert_eq!(
        fingerprint(&distribution.effective_sample_size()),
        0xb362_b541_7cb7_3645
    );
    assert_eq!(
        fingerprint(&distribution.ranked_by_weight()),
        0x2789_cca9_c741_8519
    );
    assert_eq!(
        fingerprint(&distribution.ranked_by_population()),
        0xc97a_7cef_68ff_64f6
    );
    assert_eq!(
        fingerprint(&distribution.covering_population(0.5)),
        0xd618_a48c_73fd_0377
    );
    assert_eq!(
        distribution.covering_population(0.0),
        Err(StructureSetError::InvalidFraction)
    );
    assert_eq!(fingerprint(&set.collapse().unwrap()), 0x111e_e5db_69d1_366a);

    let mut energy_set = energy_conformers(&[3.0, 1.0, 2.0]);
    assert_eq!(fingerprint(&energy_set), 0x5cc5_099f_3348_e72e);
    let distribution = energy_set.distribution().unwrap();
    assert_eq!(
        fingerprint(&distribution.ranked_by_energy()),
        0x1d17_fb39_6219_0729
    );
    let window = Quantity::new(1.0, KILOJOULE_PER_MOLE);
    assert_eq!(
        fingerprint(&distribution.within_energy_window(window)),
        0xdacb_19db_c0de_14a3
    );
    assert_eq!(
        distribution.within_energy_window(Quantity::new(-1.0, KILOJOULE_PER_MOLE)),
        Err(StructureSetError::InvalidEnergyWindow)
    );
    let temperature = Quantity::new(300.0, KELVIN);
    energy_set
        .distribution_mut()
        .unwrap()
        .set_temperature(Some(temperature))
        .unwrap();
    assert_eq!(
        energy_set.distribution().unwrap().temperature(),
        Some(temperature)
    );
    assert_eq!(
        fingerprint(&energy_set.distribution()),
        0xf47d_c1c2_0c11_71f4
    );
    assert_eq!(
        energy_set
            .distribution_mut()
            .unwrap()
            .set_temperature(Some(Quantity::new(0.0, KELVIN))),
        Err(StructureSetError::InvalidTemperature)
    );

    // Regrouping merges whole states with their representative's energy.
    let partition = Partition::new(vec![0, 0, 1], vec![1, 2]).unwrap();
    energy_set
        .distribution_mut()
        .unwrap()
        .regroup(partition)
        .unwrap();
    assert_eq!(
        energy_set
            .states()
            .map(|state| state.energy())
            .collect::<Vec<_>>(),
        [1.0, 2.0].map(|energy| Some(Quantity::new(energy, KILOJOULE_PER_MOLE)))
    );
    assert_eq!(
        fingerprint(&energy_set.distribution()),
        0x1270_e3e0_7f3f_2a5c
    );
    assert_eq!(
        energy_set
            .distribution_mut()
            .unwrap()
            .regroup(Partition::one_state(2)),
        Err(StructureSetError::PartitionMismatch {
            expected: 3,
            actual: 2
        })
    );
    set.distribution_mut().unwrap().normalize_weights().unwrap();
    assert_eq!(fingerprint(&set.distribution()), 0xb6fd_7ee1_eebd_78cd);
}

#[test]
fn sample_views_write_the_sdf_records_they_always_wrote() {
    use crate::tests::debug_fingerprint as fingerprint;
    let mut set = energy_conformers(&[2.5, 1.0]);
    set.push_member(Geometry::new(single_position(1.0)), None, 0, Some(3.0))
        .unwrap();
    set.distribution_mut()
        .unwrap()
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    assert_eq!(fingerprint(&set), 0x797f_f62d_d300_adc4);
    let options = crate::sdf::SdfWriteOptions::default();
    let views = set.iter().collect::<Vec<_>>();
    let written = crate::sdf::write(
        views.iter().map(crate::sdf::SdfRecordSource::structure),
        options,
    )
    .unwrap()
    .0;
    assert_eq!(fingerprint(&written), 0xfdd2_1ccf_3060_ad78);

    let mut scan = set_at(&[0.0]);
    scan.set_order(Some(Order::Coordinate {
        unit: DEGREE,
        values: vec![30.0],
    }))
    .unwrap();
    assert_eq!(fingerprint(&scan), 0xfcff_0753_c475_c131);
    let views = scan.iter().collect::<Vec<_>>();
    let written = crate::sdf::write(
        views.iter().map(crate::sdf::SdfRecordSource::structure),
        options,
    )
    .unwrap()
    .0;
    assert_eq!(fingerprint(&written), 0x3c58_9dc0_2772_0fe0);
}

#[test]
fn subsets_keep_the_order_and_the_distribution() {
    let (structure, atom, _) = structure_fixture();
    let mut set = StructureSet::from_structures([structure.clone(), structure]).unwrap();
    set.set_order(Some(Order::Step(vec![5, 6]))).unwrap();
    let state = State::new()
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
    assert!(shares_every_geometry(&perceived, &set));
    assert_eq!(
        StructureSet::from_structures(std::iter::empty()).unwrap_err(),
        StructureSetError::EmptySource
    );
}

/// The parts of two sets compared one by one, since a set is not
/// `PartialEq`; geometries must share their structures.
fn assert_same_set(actual: &StructureSet, expected: &StructureSet) {
    assert!(shares_every_geometry(actual, expected));
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
    let previous = set.replace(1, Geometry::new(single_position(5.0))).unwrap();
    assert_eq!(previous, Geometry::new(single_position(1.0)));
    assert_eq!(xs(&set), [0.0, 5.0]);
    assert_eq!(set.get(1).unwrap().place(), Some(Place::Step(4)));
    assert_eq!(set.get(1).unwrap().weight(), Some(2.0));
    assert!(set.get(1).unwrap().is_representative());

    let before = set.clone();
    // Two atoms do not bind to a one-atom topology.
    assert!(matches!(
        set.replace(0, Geometry::new(Positions::zeros(2))),
        Err(StructureSetError::Geometry(_))
    ));
    assert_eq!(
        set.replace(2, Geometry::new(single_position(0.0))),
        Err(StructureSetError::IndexOutOfRange { index: 2, len: 2 })
    );
    assert_same_set(&set, &before);
}

#[test]
fn the_distribution_guard_keeps_the_members() {
    let mut set = set_at(&[0.0, 1.0]);
    set.set_distribution(Some(distribution_of_one_state(2)))
        .unwrap();
    // The replacement the guard cannot make is refused when installed.
    let foreign = Distribution::new(Members::Subset(vec![7]), Partition::one_state(1)).unwrap();
    assert_eq!(
        set.set_distribution(Some(foreign)),
        Err(StructureSetError::IndexOutOfRange { index: 7, len: 2 })
    );
    let mut distribution = set.distribution_mut().unwrap();
    assert_eq!(distribution.member_count(), 2);
    distribution.set_weight(1, 3.0).unwrap();
    distribution
        .regroup(Partition::state_per_sample(2))
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
    // Without a distribution every geometry counts once.
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
    for ((geometry, weight), (position, probability)) in distribution
        .weighted_members()
        .zip(probabilities.iter().enumerate())
    {
        assert_eq!(distribution.member_position(geometry), Some(position));
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
                State::new()
                    .with_energy(Quantity::new(1.0, KILOJOULE_PER_MOLE))
                    .unwrap(),
                State::new()
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
    // Members 0 and 2 of three geometries, one state, without weights.
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
            State::new(),
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

    // A geometry that does not bind, checked last, changes nothing.
    let before = set.clone();
    let wide = || Geometry::new(Positions::zeros(2));
    assert!(matches!(
        set.push_member(wide(), Some(Place::Step(3)), 1, None),
        Err(StructureSetError::Geometry(_))
    ));
    assert!(matches!(
        set.push_state(wide(), Some(Place::Step(3)), State::new(), None),
        Err(StructureSetError::Geometry(_))
    ));
    assert!(matches!(
        set.push(wide(), Some(Place::Step(3))),
        Err(StructureSetError::Geometry(_))
    ));
    assert_same_set(&set, &before);
    // Without a distribution there is nothing to join.
    let mut plain = set_at(&[0.0]);
    assert_eq!(
        plain.push_member(Geometry::new(single_position(1.0)), None, 0, None),
        Err(StructureSetError::Missing(Requirement::Distribution))
    );
    assert_eq!(plain.len(), 1);
}

#[test]
fn places_are_set_only_on_an_order_and_inside_the_set() {
    let mut set = set_at(&[0.0]);
    assert_eq!(set.set_place(0, None), Ok(()));
    assert_eq!(
        set.set_place(0, Some(Place::Step(1))),
        Err(StructureSetError::PlaceMismatch { geometry: 0 })
    );
    assert_eq!(
        set.set_place(1, None),
        Err(StructureSetError::IndexOutOfRange { index: 1, len: 1 })
    );
    set.set_order(Some(Order::Step(vec![0]))).unwrap();
    assert_eq!(
        set.set_place(3, Some(Place::Step(1))),
        Err(StructureSetError::IndexOutOfRange { index: 3, len: 1 })
    );
    assert_eq!(
        set.set_place(0, None),
        Err(StructureSetError::PlaceMismatch { geometry: 0 })
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
        Err(StructureSetError::EmptyDistribution)
    );
    // Removing the last non-member leaves a distribution of every geometry.
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
                State::new()
                    .with_energy(Quantity::new(1.0, KILOJOULE_PER_MOLE))
                    .unwrap(),
                State::new()
                    .with_energy(Quantity::new(2.0, KILOJOULE_PER_MOLE))
                    .unwrap(),
            ])
            .unwrap(),
    ))
    .unwrap();
    let partition = Partition::new(vec![0, 1, 2], vec![0, 1, 2]).unwrap();
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
        StructureSet::from_geometries(
            single_atom_topology(),
            [
                Geometry::new(single_position(0.0)),
                Geometry::new(Positions::zeros(2))
            ],
        ),
        Err(StructureSetError::Geometry(GeometryError::Array(
            DenseArrayError::CountMismatch {
                expected: 1,
                actual: 2
            }
        )))
    ));
    let (structure, _, _) = structure_fixture();
    let other = Structure::new(Arc::new(single_atom_topology()), single_position(0.0)).unwrap();
    assert_eq!(
        StructureSet::from_structures([structure, other]).unwrap_err(),
        StructureSetError::TopologyMismatch
    );
}

#[test]
fn every_error_reads_as_a_sentence() {
    use StructureSetError as E;
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
        (Requirement::Cell, "needs a periodic cell on every geometry"),
        (Requirement::InMemory, "needs its geometries in memory"),
    ];
    for (requirement, text) in requirements {
        assert_eq!(requirement.to_string(), text);
        assert_eq!(
            E::Missing(requirement).to_string(),
            format!("The structure set {text}.")
        );
    }
    let unit = Quantity::new(1.0, ANGSTROM)
        .value_in(PICOSECOND)
        .unwrap_err();
    let errors = [
        (E::EmptySource, "There are no structures to collect."),
        (E::EmptyDistribution, "The distribution has no members."),
        (
            E::IndexOutOfRange { index: 3, len: 1 },
            "Geometry 3 is outside the structure set of 1 geometry.",
        ),
        (
            E::StateOutOfRange {
                state: 2,
                states: 2,
            },
            "State 2 is outside the distribution of 2 states.",
        ),
        (
            E::NotMember { geometry: 4 },
            "Geometry 4 is outside the distribution.",
        ),
        (
            E::TopologyMismatch,
            "The structure has a different topology layout from the structure set.",
        ),
        (
            E::GeometryCountMismatch {
                expected: 3,
                actual: 1,
            },
            "The structure set has 3 geometries, but 1 replacement was given.",
        ),
        (
            E::OrderLengthMismatch {
                expected: 2,
                actual: 1,
            },
            "The order has 1 value, but the structure set has 2 geometries.",
        ),
        (
            E::TimeAxisMismatch { times: 1, steps: 2 },
            "The time axis has 1 time and 2 steps.",
        ),
        (
            E::NonFiniteOrderValue { geometry: 1 },
            "The order value of geometry 1 is not finite.",
        ),
        (
            E::PlaceMismatch { geometry: 1 },
            "The place of geometry 1 is not on the order of the structure set.",
        ),
        (
            E::MissingTime { geometry: 1 },
            "Geometry 1 has no time, but the other geometries have one.",
        ),
        (
            E::UnexpectedTime { geometry: 1 },
            "Geometry 1 has a time, but the other geometries have none.",
        ),
        (
            E::MissingStep { geometry: 1 },
            "Geometry 1 has no step, but the other geometries have one.",
        ),
        (
            E::UnexpectedStep { geometry: 1 },
            "Geometry 1 has a step, but the other geometries have none.",
        ),
        (
            E::NonMonotonicTime { geometry: 2 },
            "Geometry 2 is earlier than the geometry before it.",
        ),
        (
            E::UnsortedMembers,
            "The members of the distribution are not in increasing order.",
        ),
        (
            E::CoverageMismatch {
                covered: 2,
                geometries: 3,
            },
            "The distribution covers 2 geometries, but the structure set has 3.",
        ),
        (
            E::PartitionMismatch {
                expected: 2,
                actual: 3,
            },
            "The partition covers 3 geometries, but the distribution has 2 members.",
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
            E::InvalidWeight { geometry: 0 },
            "The weight of geometry 0 is not finite and positive.",
        ),
        (
            E::WeightRequired,
            "The new member has no weight, but the other members have weights.",
        ),
        (
            E::RepresentativeRemoval {
                geometry: 0,
                state: 1,
            },
            "Geometry 0 represents state 1. Remove the state's other geometries first.",
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
            E::Geometry(GeometryError::NonFiniteTime),
            "A geometry is invalid.",
        ),
        (
            E::Subset(crate::topology::TopologySubsetError::EmptySelection),
            "Cannot subset the structure set.",
        ),
    ];
    for (error, text) in errors {
        assert_eq!(error.to_string(), text, "{error:?}");
    }
    assert_eq!(E::from(unit).to_string(), unit.to_string());
    // Wrapped causes are the error's source.
    use std::error::Error;
    assert!(E::Geometry(GeometryError::NonFiniteTime).source().is_some());
    assert!(E::EmptySource.source().is_none());
}

#[test]
fn failed_and_empty_sample_edits_keep_the_structure_shared() {
    let shared = Geometry::new(single_position(0.0)).with_cell(cell());
    let mut geometry = shared.clone();
    assert!(geometry
        .set_positions(Quantity::new(vec![Point3::origin(); 2], NANOMETER))
        .is_err());
    assert!(geometry
        .set_position_at(
            TopologyAtomIndex::new(1),
            Quantity::new(Point3::origin(), NANOMETER)
        )
        .is_err());
    assert!(geometry
        .set_occupancy(TopologyAtomIndex::new(1), Some(0.5))
        .is_err());
    assert!(geometry
        .set_occupancies(Some(vec![Some(f64::NAN)]))
        .is_err());
    assert!(geometry
        .set_b_factors(Some(Quantity::new(vec![Some(1.0)], NANOMETER)))
        .is_err());
    assert!(geometry.set_velocities(Some(Velocities::zeros(2))).is_err());
    geometry.set_velocities(None).unwrap();
    geometry.set_forces(None).unwrap();
    assert_eq!(geometry.take_velocities(), None);
    geometry.set_cell(Some(cell()));
    assert!(geometry.shares_structure_with(&shared));
    // A successful edit copies this geometry alone.
    geometry.set_cell(None);
    assert!(!geometry.shares_structure_with(&shared));
    assert_eq!(shared.cell(), Some(&cell()));

    // The same through a set's guard.
    let mut set = StructureSet::from_geometries(single_atom_topology(), [shared.clone()]).unwrap();
    let before = set.clone();
    assert!(set
        .get_mut(0)
        .unwrap()
        .set_occupancies(Some(vec![Some(1.0), Some(1.0)]))
        .is_err());
    set.get_mut(0).unwrap().set_cell(Some(cell()));
    assert!(shares_every_geometry(&set, &before));
}
