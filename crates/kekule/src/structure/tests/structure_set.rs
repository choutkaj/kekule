//! Structure sets: the order, views, and copy on write.

use super::*;
use crate::diagnostics::{Countable, IndexOutOfRange, LengthMismatch, TopologyMismatch};
use crate::units::{DEGREE, PICOSECOND};

/// Geometries of a one-atom topology at these x positions, unordered.
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

fn shares_every_geometry(a: &StructureSet, b: &StructureSet) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b.iter())
            .all(|(a, b)| Arc::ptr_eq(a.payload().data(), b.payload().data()))
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
fn constructors_refuse_mismatched_geometries_and_layouts() {
    assert!(matches!(
        StructureSet::from_geometries(
            single_atom_topology(),
            [
                Geometry::new(single_position(0.0)),
                Geometry::new(Positions::zeros(2))
            ],
        ),
        Err(StructureSetError::Geometry(GeometryError::Array(error)))
            if *error == DenseArrayError::CountMismatch(LengthMismatch::new(Countable::Row, 1, 2))
    ));
    let (structure, _, _) = structure_fixture();
    let other = Structure::new(Arc::new(single_atom_topology()), single_position(0.0)).unwrap();
    assert_eq!(
        StructureSet::from_structures([structure, other]).unwrap_err(),
        StructureSetError::TopologyMismatch(TopologyMismatch::Structure)
    );
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
