//! Structure sets and ensembles: the order, the distribution, views, and
//! copy on write.

use super::*;
use crate::diagnostics::{
    CapacityExceeded, Countable, IndexOutOfRange, LengthMismatch, TopologyMismatch,
};
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
    samples: Samples,
    assignments: &[usize],
    representatives: &[usize],
) -> Distribution {
    let partition = Partition::new(assignments.to_vec(), representatives.to_vec()).unwrap();
    Distribution::new(samples, partition).unwrap()
}

fn velocities(x: f64) -> Velocities {
    Velocities::new(Quantity::new(
        [Vector3::new(x, 0.0, 0.0)],
        CANONICAL_VELOCITY,
    ))
    .unwrap()
}

const CANONICAL_VELOCITY: crate::units::Unit = crate::units::CANONICAL_VELOCITY_UNIT;

/// The geometries at these x positions as `distribution` reads them.
fn ensemble_at(xs: &[f64], distribution: Distribution) -> Ensemble {
    Ensemble::new(set_at(xs), distribution).unwrap()
}

fn shares_every_geometry(a: &StructureSet, b: &StructureSet) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b.iter())
            .all(|(a, b)| Arc::ptr_eq(a.payload().data(), b.payload().data()))
}

#[test]
fn requirements_read_as_they_are() {
    assert_eq!(Requirement::TimeOrder.to_string(), "needs a time order");
    assert_eq!(
        Requirement::StateEnergies.to_string(),
        "needs an energy for every state"
    );
    assert_eq!(
        StructureSetError::Missing(Requirement::TimeOrder).to_string(),
        "The structure set needs a time order."
    );
    let set = set_at(&[0.0]);
    assert_eq!(
        set.validate_monotonic_time(),
        Err(StructureSetError::Missing(Requirement::TimeOrder))
    );
}

#[test]
fn an_order_has_one_value_per_geometry_in_every_column() {
    let mut set = set_at(&[0.0, 1.0]);
    assert_eq!(
        set.set_order(Some(Order::Step(vec![0]))),
        Err(StructureSetError::OrderLengthMismatch(LengthMismatch::new(
            Countable::OrderValue,
            2,
            1
        )))
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
fn order_values_belong_to_the_set_not_to_its_geometries() {
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
fn any_geometry_may_carry_site_fields_and_dynamics() {
    let mut set = set_at(&[0.0]);
    set.set_order(Some(Order::Step(vec![0]))).unwrap();
    let mut geometry = set.get_mut(0).unwrap();
    geometry.set_occupancies(Some(vec![Some(0.5)])).unwrap();
    geometry.set_velocities(Some(velocities(1.0))).unwrap();
    assert_eq!(
        geometry.set_forces(Some(Forces::zeros(2))),
        Err(GeometryError::Array(DenseArrayError::CountMismatch(
            LengthMismatch::new(Countable::Row, 1, 2)
        )))
    );
    let view = set.get(0).unwrap();
    assert_eq!(view.occupancy(ATOM0), Some(0.5));
    let StructureKind::InSet(context) = view.as_structure_view().kind() else {
        panic!("a geometry's view is a geometry");
    };
    assert_eq!(context.place(), Some(Place::Step(0)));
    assert_eq!(
        view.as_structure_view().velocities(),
        Some(&velocities(1.0))
    );
    assert_eq!(view.as_structure_view().forces(), None);
    assert_eq!(view.as_structure_view().occupancy(ATOM0), Some(0.5));
    // A structure keeps the geometry, velocities included, and drops its
    // context.
    let structure = view.to_structure();
    assert_eq!(
        structure.as_structure_view().kind(),
        StructureKind::Standalone
    );
    assert_eq!(structure.occupancy(ATOM0), Some(0.5));
    assert_eq!(structure.velocities(), Some(&velocities(1.0)));
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
fn a_distribution_covers_all_geometries_or_a_sorted_subset() {
    assert_eq!(
        Distribution::new(Samples::Subset(vec![2, 1]), Partition::one_state(2)),
        Err(StructureSetError::UnsortedSamples)
    );
    assert_eq!(
        Distribution::new(Samples::Subset(vec![1, 2]), Partition::one_state(3)),
        Err(StructureSetError::PartitionMismatch(LengthMismatch::new(
            Countable::Sample,
            2,
            3
        )))
    );
    assert_eq!(
        Ensemble::new(
            set_at(&[0.0, 1.0, 2.0]),
            distribution(Samples::All, &[0, 0], &[0])
        )
        .unwrap_err(),
        StructureSetError::CoverageMismatch {
            covered: 2,
            geometries: 3
        }
    );
    assert_eq!(
        Ensemble::new(
            set_at(&[0.0, 1.0, 2.0]),
            distribution(Samples::Subset(vec![1, 3]), &[0, 0], &[0])
        )
        .unwrap_err(),
        StructureSetError::IndexOutOfRange(IndexOutOfRange::new(Countable::Geometry, 3, 3))
    );
    // A subset listing every geometry is all of them.
    let all = ensemble_at(
        &[0.0, 1.0, 2.0],
        distribution(Samples::Subset(vec![0, 1, 2]), &[0, 0, 1], &[0, 2]),
    );
    assert_eq!(all.distribution().samples(), &Samples::All);

    let set = ensemble_at(
        &[0.0, 1.0, 2.0],
        distribution(Samples::Subset(vec![1, 2]), &[0, 1], &[0, 1]),
    );
    let distribution = set.distribution();
    assert_eq!(distribution.sample_position(0), None);
    assert_eq!(distribution.sample_position(2), Some(1));
    assert_eq!(distribution.state_of(2), Some(1));
    assert_eq!(distribution.state_of(0), None);
    assert_eq!(distribution.log_weight(0), None);
    assert_eq!(set.state(0).unwrap().representative(), 1);
    assert_eq!(set.state(1).unwrap().representative(), 2);
    assert_eq!(set.state(1).unwrap().geometries().collect::<Vec<_>>(), [2]);
    assert_eq!(set.states().len(), 2);
    assert!(set.state(2).is_none());
}

#[test]
fn log_weights_are_finite_and_equal_when_absent() {
    let mut set = ensemble_at(
        &[0.0, 1.0, 2.0],
        distribution(Samples::Subset(vec![0, 2]), &[0, 0], &[0]),
    );
    let mut distribution = set.distribution_mut();
    assert_eq!(distribution.log_weights(), None);
    assert_eq!(distribution.log_weight(2), Some(0.0));
    assert_eq!(distribution.log_weight(1), None);
    assert_eq!(distribution.probabilities().unwrap(), [0.5, 0.5]);
    assert_eq!(distribution.effective_sample_size().unwrap(), 2.0);
    assert_eq!(
        distribution.set_log_weights(Some(vec![0.0])),
        Err(StructureSetError::LogWeightCountMismatch(
            LengthMismatch::new(Countable::LogWeight, 2, 1)
        ))
    );
    assert_eq!(
        distribution.set_log_weights(Some(vec![0.0, f64::NAN])),
        Err(StructureSetError::NonFiniteLogWeight { geometry: 2 })
    );
    assert_eq!(
        distribution.set_log_weight(1, 2.0),
        Err(StructureSetError::NotSample { geometry: 1 })
    );
    // Failed edits keep every sample equally weighted.
    assert_eq!(distribution.log_weights(), None);
    let three = 3.0_f64.ln();
    distribution.set_log_weight(2, three).unwrap();
    assert_eq!(distribution.log_weights(), Some(&[0.0, three][..]));
    assert_eq!(distribution.ranked_by_weight(), [2, 0]);
    distribution.normalize_log_weights().unwrap();
    let normalized = distribution.log_weights().unwrap().to_vec();
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
        distribution.log_weights(),
        Some(&[normalized[0], normalized[1], half, 0.0][..])
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

fn distribution_of_one_state(samples: usize) -> Distribution {
    Distribution::new(Samples::All, Partition::one_state(samples)).unwrap()
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
    assert_eq!(distribution.log_weights(), Some(&[four, three][..]));
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
    assert_eq!(states.distribution().log_weights(), Some(&[two, four][..]));
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
    assert_eq!(distribution.log_weights(), Some(&[four][..]));
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
    distribution
        .set_log_weights(Some(weights.to_vec()))
        .unwrap();
    assert!(!Arc::ptr_eq(&regrouped, &shared));
    assert!(shares_every_geometry(regrouped.set(), shared.set()));
    assert_eq!(shared.distribution().state_count(), 1);

    let selected = shared.select([2, 0]).unwrap();
    assert!(Arc::ptr_eq(
        selected.set().get(0).unwrap().payload().data(),
        shared.set().get(2).unwrap().payload().data()
    ));
    // Each state of one sample keeps its log weight as its log population.
    let collapsed = regrouped.collapse();
    assert_eq!(collapsed.set().order(), None);
    assert_eq!(collapsed.distribution().log_weights(), Some(&weights[..]));
    assert!(shares_every_geometry(collapsed.set(), regrouped.set()));
    let states = regrouped.select_states([1]).unwrap();
    assert!(Arc::ptr_eq(
        states.set().get(0).unwrap().payload().data(),
        shared.set().get(1).unwrap().payload().data()
    ));

    // Editing one geometry of a clone copies that geometry alone.
    let mut edited = (*shared).clone();
    edited
        .get_mut(1)
        .unwrap()
        .set_position_at(ATOM0, Quantity::new(Point3::new(9.0, 0.0, 0.0), NANOMETER))
        .unwrap();
    assert!(Arc::ptr_eq(
        edited.set().get(0).unwrap().payload().data(),
        shared.set().get(0).unwrap().payload().data()
    ));
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
fn distribution_operations_keep_their_results() {
    // State 0: geometries 0 and 2 (weights 1 and 3), represented by 2;
    // state 1: geometries 1 and 3 (weights 2 and 4), represented by 1.
    let weights = [1.0_f64, 2.0, 3.0, 4.0].map(f64::ln);
    let mut set = grouped(&weights, &[0, 1, 0, 1], &[2, 1]);
    assert_eq!(set.set().len(), 4);
    let distribution = set.distribution();
    assert_eq!(distribution.samples(), &Samples::All);
    assert_eq!(distribution.log_weights(), Some(&weights[..]));
    assert_eq!(distribution.partition().assignments(), [0, 1, 0, 1]);
    assert_eq!(distribution.partition().representatives(), [2, 1]);
    assert_eq!(distribution.temperature(), None);
    // Weights of a total of 10.
    assert_close(
        &distribution.probabilities().unwrap(),
        &[0.1, 0.2, 0.3, 0.4],
    );
    assert_close(&distribution.state_probabilities().unwrap(), &[0.4, 0.6]);
    // Kish: 10² / (1 + 4 + 9 + 16) = 10 / 3.
    assert!(close(
        distribution.effective_sample_size().unwrap(),
        10.0 / 3.0
    ));
    assert_eq!(distribution.ranked_by_weight(), [3, 2, 1, 0]);
    assert_eq!(distribution.ranked_by_population(), [1, 0]);
    assert_eq!(distribution.covering_population(0.5).unwrap(), [1]);
    assert_eq!(distribution.covering_population(0.7).unwrap(), [1, 0]);
    assert_eq!(
        distribution.covering_population(0.0),
        Err(StructureSetError::InvalidFraction)
    );
    let collapsed = set.collapse();
    assert_eq!(collapsed.set().len(), 2);
    assert!(collapsed
        .set()
        .get(0)
        .unwrap()
        .shares_storage_with(&set.set().get(2).unwrap()));
    assert!(collapsed
        .set()
        .get(1)
        .unwrap()
        .shares_storage_with(&set.set().get(1).unwrap()));
    assert_eq!(
        collapsed.distribution().partition(),
        &Partition::state_per_sample(2)
    );
    assert_close(&log_weights(&collapsed), &[4.0_f64.ln(), 6.0_f64.ln()]);

    let mut energy_set = energy_conformers(&[3.0, 1.0, 2.0]);
    let state = |energy| StateData::new().with_energy(kj(energy).unwrap()).unwrap();
    let expected = Distribution::new(Samples::All, Partition::state_per_sample(3))
        .unwrap()
        .with_log_weights(vec![0.0; 3])
        .unwrap()
        .with_states(vec![state(3.0), state(1.0), state(2.0)])
        .unwrap();
    assert_eq!(energy_set.set().len(), 3);
    assert_eq!(energy_set.distribution(), &expected);
    let distribution = energy_set.distribution();
    assert_eq!(distribution.ranked_by_energy().unwrap(), [1, 2, 0]);
    // States 1 and 2 lie within 1 kJ/mol of the lowest; the bound is inclusive.
    let window = Quantity::new(1.0, KILOJOULE_PER_MOLE);
    assert_eq!(distribution.within_energy_window(window).unwrap(), [1, 2]);
    assert_eq!(
        distribution.within_energy_window(Quantity::new(-1.0, KILOJOULE_PER_MOLE)),
        Err(StructureSetError::InvalidEnergyWindow)
    );
    let temperature = Quantity::new(300.0, KELVIN);
    energy_set
        .distribution_mut()
        .set_temperature(Some(temperature))
        .unwrap();
    assert_eq!(energy_set.distribution().temperature(), Some(temperature));
    // A temperature changes nothing else.
    assert_eq!(
        energy_set.distribution(),
        &expected.clone().with_temperature(temperature).unwrap()
    );
    assert_eq!(
        energy_set
            .distribution_mut()
            .set_temperature(Some(Quantity::new(0.0, KELVIN))),
        Err(StructureSetError::InvalidTemperature)
    );

    // Regrouping merges whole states with their representative's energy.
    let partition = Partition::new(vec![0, 0, 1], vec![1, 2]).unwrap();
    energy_set
        .distribution_mut()
        .regroup(partition.clone())
        .unwrap();
    assert_eq!(
        energy_set
            .states()
            .map(|state| state.energy())
            .collect::<Vec<_>>(),
        [1.0, 2.0].map(|energy| Some(Quantity::new(energy, KILOJOULE_PER_MOLE)))
    );
    let regrouped = Distribution::new(Samples::All, partition)
        .unwrap()
        .with_log_weights(vec![0.0; 3])
        .unwrap()
        .with_states(vec![state(1.0), state(2.0)])
        .unwrap()
        .with_temperature(temperature)
        .unwrap();
    assert_eq!(energy_set.distribution(), &regrouped);
    assert_eq!(
        energy_set
            .distribution_mut()
            .regroup(Partition::one_state(2)),
        Err(StructureSetError::PartitionMismatch(LengthMismatch::new(
            Countable::Sample,
            3,
            2
        )))
    );
    assert_eq!(energy_set.distribution(), &regrouped);

    // Normalizing changes the log weights alone, to log probabilities.
    let before = set.distribution().clone();
    set.distribution_mut().normalize_log_weights().unwrap();
    let distribution = set.distribution();
    assert_close(
        &exponentials(distribution.log_weights().unwrap()),
        &[0.1, 0.2, 0.3, 0.4],
    );
    assert_eq!(distribution.samples(), before.samples());
    assert_eq!(distribution.partition(), before.partition());
    assert_eq!(distribution.states(), before.states());
    assert_eq!(distribution.temperature(), None);
}

/// A V2000 record of one carbon atom at `x` Å, followed by `data`.
fn carbon_record(x: &str, data: &str) -> String {
    format!(
        "\n  kekule            2D\n\n  1  0  0  0  0  0            999 V2000\n\
         {x:>10}    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n\
         {data}$$$$\n"
    )
}

#[test]
fn geometry_views_write_the_sdf_records_they_always_wrote() {
    let options = crate::sdf::SdfWriteOptions::default();
    // One carbon at x = 0 and 0.1 nm: flat, so 2D, with coordinates in Å
    // and no data fields.
    let set = set_at(&[0.0, 0.1]);
    let views = set.iter().collect::<Vec<_>>();
    let written = crate::sdf::write(
        views.iter().map(crate::sdf::SdfRecordSource::structure),
        options,
    )
    .unwrap()
    .0;
    assert_eq!(
        written,
        carbon_record("0.0000", "") + &carbon_record("1.0000", "")
    );

    // A geometry on a coordinate order writes its coordinate and its unit
    // as data fields.
    let mut scan = set_at(&[0.0]);
    scan.set_order(Some(Order::Coordinate {
        unit: DEGREE,
        values: vec![30.0],
    }))
    .unwrap();
    let views = scan.iter().collect::<Vec<_>>();
    let written = crate::sdf::write(
        views.iter().map(crate::sdf::SdfRecordSource::structure),
        options,
    )
    .unwrap()
    .0;
    assert_eq!(
        written,
        carbon_record(
            "0.0000",
            ">  <kekule.coordinate>\n30\n\n>  <kekule.coordinate.unit>\ndeg\n\n"
        )
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
            .with_temperature(Quantity::new(300.0, KELVIN))
            .unwrap(),
    )
    .unwrap();
    let selection =
        crate::topology::AtomSelection::from_atoms(set.shared_topology(), [atom]).unwrap();
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
    assert_eq!(distribution.log_weights(), Some(&[2.0][..]));
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

/// The parts of two ensembles compared one by one, since a set is not
/// `PartialEq`; geometries must share their structures.
fn assert_same_ensemble(actual: &Ensemble, expected: &Ensemble) {
    assert!(shares_every_geometry(actual.set(), expected.set()));
    assert_eq!(actual.set().order(), expected.set().order());
    assert_eq!(actual.distribution(), expected.distribution());
    assert_eq!(actual.set().properties(), expected.set().properties());
}

#[test]
fn replacing_a_geometry_keeps_its_place_and_sample() {
    let mut frames = set_at(&[0.0, 1.0]);
    frames.set_order(Some(Order::Step(vec![3, 4]))).unwrap();
    let two = 2.0_f64.ln();
    let mut set = Ensemble::new(
        frames,
        distribution(Samples::Subset(vec![1]), &[0], &[0])
            .with_log_weights(vec![two])
            .unwrap(),
    )
    .unwrap();
    let previous = set.replace(1, Geometry::new(single_position(5.0))).unwrap();
    assert_eq!(previous, Geometry::new(single_position(1.0)));
    assert_eq!(xs(set.set()), [0.0, 5.0]);
    assert_eq!(set.set().get(1).unwrap().place(), Some(Place::Step(4)));
    assert_eq!(set.distribution().log_weight(1), Some(two));
    assert_eq!(set.state(0).unwrap().representative(), 1);

    let before = set.clone();
    // Two atoms do not bind to a one-atom topology.
    assert!(matches!(
        set.replace(0, Geometry::new(Positions::zeros(2))),
        Err(StructureSetError::Geometry(_))
    ));
    assert_eq!(
        set.replace(2, Geometry::new(single_position(0.0))),
        Err(StructureSetError::IndexOutOfRange(IndexOutOfRange::new(
            Countable::Geometry,
            2,
            2
        )))
    );
    assert_same_ensemble(&set, &before);
}

#[test]
fn the_distribution_guard_keeps_the_samples() {
    // The replacement the guard cannot make is refused when paired.
    let foreign = Distribution::new(Samples::Subset(vec![7]), Partition::one_state(1)).unwrap();
    assert_eq!(
        Ensemble::new(set_at(&[0.0, 1.0]), foreign).unwrap_err(),
        StructureSetError::IndexOutOfRange(IndexOutOfRange::new(Countable::Geometry, 7, 2))
    );
    let mut set = ensemble_at(&[0.0, 1.0], distribution_of_one_state(2));
    let mut distribution = set.distribution_mut();
    assert_eq!(distribution.sample_count(), 2);
    distribution.set_log_weight(1, 3.0_f64.ln()).unwrap();
    distribution
        .regroup(Partition::state_per_sample(2))
        .unwrap();
    distribution
        .state_mut(1)
        .unwrap()
        .set_energy(kj(1.0))
        .unwrap();
    assert_eq!(set.distribution().samples(), &Samples::All);
    assert_eq!(set.state(1).unwrap().representative_sample().index(), 1);
    assert_eq!(set.state(1).unwrap().energy(), kj(1.0));
    assert_eq!(set.distribution().log_weight(1), Some(3.0_f64.ln()));
}

#[test]
fn weighted_observations_join_samples_and_probabilities() {
    use crate::statistics::SampleObservations;
    let set = set_at(&[0.0, 1.0, 2.0]);
    // A set's statistics count every geometry once.
    let all = SampleObservations::all(&set);
    assert_eq!(all.samples(), [0, 1, 2]);
    assert_eq!(all.weights(), [1.0, 1.0, 1.0]);
    let mut ensemble = Ensemble::new(
        set,
        distribution(Samples::Subset(vec![0, 2]), &[0, 0], &[0]),
    )
    .unwrap();
    // An ensemble's observe its samples with their probabilities.
    let observed = SampleObservations::weighted(&ensemble);
    assert_eq!(observed.samples(), [0, 2]);
    assert_close(observed.weights(), &[0.5, 0.5]);
    assert_eq!(
        ensemble.distribution().sample_probabilities().unwrap(),
        [(0, 0.5), (2, 0.5)]
    );
    ensemble
        .distribution_mut()
        .set_log_weight(2, 4.0_f64.ln())
        .unwrap();
    let distribution = ensemble.distribution();
    let joined = distribution.sample_probabilities().unwrap();
    assert_eq!(
        joined
            .iter()
            .map(|(geometry, _)| *geometry)
            .collect::<Vec<_>>(),
        [0, 2]
    );
    // Weights 1 and 4 of 5.
    let probabilities = joined
        .iter()
        .map(|(_, probability)| *probability)
        .collect::<Vec<_>>();
    assert_close(&probabilities, &[0.2, 0.8]);
    let observed = SampleObservations::weighted(&ensemble);
    assert_eq!(observed.samples(), [0, 2]);
    assert_close(observed.weights(), &probabilities);
    // Sample positions index the log weights and probabilities.
    let by_position = distribution.probabilities().unwrap();
    for (position, (geometry, probability)) in joined.iter().enumerate() {
        assert_eq!(distribution.sample_position(*geometry), Some(position));
        assert_eq!(*probability, by_position[position]);
        assert_eq!(
            distribution.log_weight(*geometry),
            Some(distribution.log_weights().unwrap()[position])
        );
    }
    assert_eq!(distribution.ranked_by_weight(), [2, 0]);
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
    assert_eq!(distribution.log_weights(), Some(&[two, two][..]));
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
    assert_eq!(distribution.log_weights(), Some(&[0.0, 0.0][..]));
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
        collapsed.distribution().log_weights(),
        Some(&[2.0_f64.ln()][..])
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
        collapsed.distribution().log_weights(),
        Some(&[1000.0 + 2.0_f64.ln(), 1000.0][..])
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
    assert_eq!(distribution.log_weights(), None);
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
        Err(StructureSetError::IndexOutOfRange(IndexOutOfRange::new(
            Countable::Geometry,
            1,
            1
        )))
    );
    set.set_order(Some(Order::Step(vec![0]))).unwrap();
    assert_eq!(
        set.set_place(3, Some(Place::Step(1))),
        Err(StructureSetError::IndexOutOfRange(IndexOutOfRange::new(
            Countable::Geometry,
            3,
            1
        )))
    );
    assert_eq!(
        set.set_place(0, None),
        Err(StructureSetError::PlaceMismatch { geometry: 0 })
    );
    assert_eq!(set.order(), Some(&Order::Step(vec![0])));
}

#[test]
fn empty_selections_and_distributions() {
    let mut frames = set_at(&[0.0, 1.0]);
    frames.set_order(Some(Order::Step(vec![0, 1]))).unwrap();
    let mut set =
        Ensemble::new(frames, distribution(Samples::Subset(vec![0]), &[0], &[0])).unwrap();
    let none = set.select([]).unwrap();
    assert!(none.set().is_empty());
    assert_eq!(none.set().order(), Some(&Order::Step(Vec::new())));
    assert_eq!(none.distribution().sample_count(), 0);
    assert_eq!(none.states().len(), 0);
    let none = set.select_states([]).unwrap();
    assert!(none.set().is_empty());
    assert_eq!(none.distribution().state_count(), 0);
    // A distribution without samples collapses to an empty ensemble.
    let collapsed = none.collapse();
    assert!(collapsed.set().is_empty());
    assert_eq!(collapsed.distribution().state_count(), 0);
    assert_eq!(
        none.distribution().probabilities(),
        Err(StructureSetError::EmptyDistribution)
    );
    // Removing the last geometry that is not a sample leaves a distribution
    // of every geometry.
    set.remove(1).unwrap();
    assert_eq!(set.distribution().samples(), &Samples::All);
    assert_eq!(set.set().order(), Some(&Order::Step(vec![0])));
    // Removing the last sample leaves an empty distribution.
    set.remove(0).unwrap();
    assert!(set.set().is_empty());
    assert_eq!(set.distribution().state_count(), 0);
    assert_eq!(set.distribution().samples(), &Samples::All);
}

#[test]
fn regrouping_a_split_state_drops_its_energy() {
    let mut set = ensemble_at(
        &[0.0, 1.0, 2.0],
        distribution(Samples::All, &[0, 0, 1], &[0, 2])
            .with_states(vec![
                StateData::new()
                    .with_energy(Quantity::new(1.0, KILOJOULE_PER_MOLE))
                    .unwrap(),
                StateData::new()
                    .with_energy(Quantity::new(2.0, KILOJOULE_PER_MOLE))
                    .unwrap(),
            ])
            .unwrap(),
    );
    let partition = Partition::new(vec![0, 1, 2], vec![0, 1, 2]).unwrap();
    set.distribution_mut().regroup(partition).unwrap();
    assert_eq!(set.states().len(), 3);
    assert_eq!(set.state(0).unwrap().energy(), None);
    assert_eq!(set.state(1).unwrap().energy(), None);
    // State 1 was whole and keeps its energy.
    assert_eq!(set.state(2).unwrap().energy(), kj(2.0));
}

#[test]
fn constructors_refuse_mismatched_geometries_and_layouts() {
    assert!(matches!(
        StructureSet::from_geometries(
            single_atom_topology(),
            [
                Geometry::new(single_position(0.0)),
                Geometry::new(Positions::zeros(2))
            ],
        ),
        Err(StructureSetError::Geometry(GeometryError::Array(
            DenseArrayError::CountMismatch(payload)
        ))) if payload == LengthMismatch::new(Countable::Row, 1, 2)
    ));
    let (structure, _, _) = structure_fixture();
    let other = Structure::new(Arc::new(single_atom_topology()), single_position(0.0)).unwrap();
    assert_eq!(
        StructureSet::from_structures([structure, other]).unwrap_err(),
        StructureSetError::TopologyMismatch(TopologyMismatch::Structure)
    );
}

#[test]
fn every_error_reads_as_a_sentence() {
    use StructureSetError as E;
    let requirements = [
        (Requirement::TimeOrder, "needs a time order"),
        (Requirement::AnyOrder, "needs an order"),
        (Requirement::Times, "needs times"),
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
        (E::EmptyDistribution, "The distribution has no samples."),
        (
            E::IndexOutOfRange(IndexOutOfRange::new(Countable::Geometry, 3, 1)),
            "Geometry 3 does not exist. There is 1 geometry.",
        ),
        (
            E::StateOutOfRange(IndexOutOfRange::new(Countable::State, 2, 2)),
            "State 2 does not exist. There are 2 states.",
        ),
        (
            E::NotSample { geometry: 4 },
            "Geometry 4 is not a sample of the distribution.",
        ),
        (
            E::TopologyMismatch(TopologyMismatch::Structure),
            "The structure has a different topology layout.",
        ),
        (
            E::GeometryCountMismatch(LengthMismatch::new(Countable::Geometry, 3, 1)),
            "There is 1 geometry instead of 3.",
        ),
        (
            E::OrderLengthMismatch(LengthMismatch::new(Countable::OrderValue, 2, 1)),
            "There is 1 order value instead of 2.",
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
            E::UnsortedSamples,
            "The samples of the distribution are not in increasing order.",
        ),
        (
            E::CoverageMismatch {
                covered: 2,
                geometries: 3,
            },
            "The distribution covers 2 geometries, but the structure set has 3.",
        ),
        (
            E::PartitionMismatch(LengthMismatch::new(Countable::Sample, 2, 3)),
            "There are 3 samples instead of 2.",
        ),
        (
            E::LogWeightCountMismatch(LengthMismatch::new(Countable::LogWeight, 2, 1)),
            "There is 1 log weight instead of 2.",
        ),
        (
            E::LogWeightCountMismatch(LengthMismatch::new(Countable::LogWeight, 1, 2)),
            "There are 2 log weights instead of 1.",
        ),
        (
            E::StateCountMismatch(LengthMismatch::new(Countable::State, 1, 2)),
            "There are 2 states instead of 1.",
        ),
        (
            E::NonFiniteLogWeight { geometry: 0 },
            "The log weight of geometry 0 is not finite.",
        ),
        (
            E::LogWeightRequired,
            "The new sample has no log weight, but the other samples have log weights.",
        ),
        (
            E::RepresentativeRemoval {
                geometry: 0,
                state: 1,
            },
            "Geometry 0 represents state 1. Remove the state's other samples first.",
        ),
        (
            E::InvalidTemperature,
            "The temperature is not finite and positive.",
        ),
        (
            E::UnrepresentablePopulation { state: 1 },
            "The Boltzmann population of state 1 is too small to represent at this temperature. Raise the temperature.",
        ),
        (
            E::InvalidEnergyWindow,
            "The energy window is not finite and nonnegative.",
        ),
        (
            E::InvalidFraction,
            "The population fraction is not greater than 0 and at most 1.",
        ),
        (E::NonFiniteTime, "The time is not finite."),
        (E::NonFiniteEnergy, "The state energy is not finite."),
        (
            E::Geometry(GeometryError::CapacityExceeded(CapacityExceeded::new(Countable::Atom))),
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
    assert!(
        E::Geometry(GeometryError::CapacityExceeded(CapacityExceeded::new(
            Countable::Atom
        )))
        .source()
        .is_some()
    );
    assert!(E::EmptySource.source().is_none());
}

#[test]
fn failed_and_empty_geometry_edits_keep_the_storage_shared() {
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
    assert!(geometry.shares_storage_with(&shared));
    // A successful edit copies this geometry alone.
    geometry.set_cell(None);
    assert!(!geometry.shares_storage_with(&shared));
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
