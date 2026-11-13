//! Structure sets on a step or a coordinate order: optimization histories,
//! relaxed scans, and paths.

use std::sync::Arc;

use kekule::alignment::AlignmentOptions;
use kekule::math::{PeriodicCell, Point3, Vector3};
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{
    AsStructureView, Geometry, GeometryError, Order, Place, Positions, Structure, StructureKind,
    StructureSet, StructureSetError,
};
use kekule::topology::{AtomSelection, Topology};
use kekule::units::{
    Quantity, ANGSTROM, DEGREE, KILOJOULE_PER_MOLE, NANOMETER, RADIAN, SQUARE_NANOMETER,
};

fn chain(atoms: usize) -> Arc<Topology> {
    let mut editor = MoleculeEditor::new();
    let mut previous = None;
    for _ in 0..atoms {
        let atom = editor
            .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
            .unwrap();
        if let Some(previous) = previous {
            editor.add_bond(previous, atom, BondOrder::Single).unwrap();
        }
        previous = Some(atom);
    }
    Arc::new(Topology::from_molecule(editor.finish().unwrap()).unwrap())
}

/// Positions of a bent three-carbon chain, its last atom moved by `shift`.
fn positions(shift: f64) -> Positions {
    Positions::new(Quantity::new(
        [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.15, 0.0, 0.0),
            Point3::new(0.2, 0.14 + shift, 0.0),
        ],
        NANOMETER,
    ))
    .unwrap()
}

fn energy() -> PropertyKey {
    PropertyKey::new("energy").unwrap()
}

fn with_energy(mut geometry: Geometry, value: f64) -> Geometry {
    geometry
        .properties_mut()
        .owner_mut()
        .insert(
            energy(),
            PropertyValue::real(value, KILOJOULE_PER_MOLE).unwrap(),
        )
        .unwrap();
    geometry
}

/// An empty set of three-carbon chains on `order`.
fn ordered(order: Order) -> StructureSet {
    let mut set = StructureSet::new(chain(3));
    set.set_order(Some(order)).unwrap();
    set
}

fn coordinates(unit: kekule::units::Unit) -> Order {
    Order::Coordinate {
        unit,
        values: Vec::new(),
    }
}

/// A history of four optimization steps with energies.
fn history() -> StructureSet {
    let mut set = ordered(Order::Step(Vec::new()));
    for step in 0..4 {
        let geometry = with_energy(
            Geometry::new(positions(0.01 * step as f64)),
            -100.0 - step as f64,
        );
        set.push(geometry, Some(Place::Step(step * 5))).unwrap();
    }
    set
}

fn places(set: &StructureSet) -> Vec<Option<Place>> {
    set.iter().map(|geometry| geometry.place()).collect()
}

fn steps(steps: &[u64]) -> Vec<Option<Place>> {
    steps.iter().map(|step| Some(Place::Step(*step))).collect()
}

/// The coordinate of a geometry on a coordinate order.
fn coordinate(place: Option<Place>) -> Quantity<f64> {
    match place {
        Some(Place::Coordinate(coordinate)) => coordinate,
        other => panic!("not a coordinate: {other:?}"),
    }
}

fn real(value: &PropertyValue) -> f64 {
    match value {
        PropertyValue::Real { value, .. } => *value,
        other => panic!("not a real value: {other:?}"),
    }
}

#[test]
fn a_history_keeps_its_geometries_steps_and_values_in_order() {
    let set = history();
    assert_eq!(set.order(), Some(&Order::Step(vec![0, 5, 10, 15])));
    assert_eq!(set.len(), 4);
    assert_eq!(places(&set), steps(&[0, 5, 10, 15]));
    let energies = set
        .iter()
        .map(|geometry| real(geometry.properties().owner().get(&energy()).unwrap()))
        .collect::<Vec<_>>();
    assert_eq!(energies, [-100.0, -101.0, -102.0, -103.0]);
    let geometry = set.get(2).unwrap();
    assert_eq!(geometry.index(), 2);
    assert_eq!(geometry.place(), Some(Place::Step(10)));
    assert!(set.get(4).is_none());
    assert_eq!(set.iter().len(), 4);
}

#[test]
fn a_scan_stores_coordinates_in_its_order_unit_in_any_order() {
    let mut scan = ordered(coordinates(DEGREE));
    // A scan swept forward and back, given in radians and degrees.
    for (index, target) in [
        Quantity::new(0.0, RADIAN),
        Quantity::new(std::f64::consts::FRAC_PI_2, RADIAN),
        Quantity::new(180.0, DEGREE),
        Quantity::new(90.0, DEGREE),
    ]
    .into_iter()
    .enumerate()
    {
        let geometry = Geometry::new(positions(0.01 * index as f64));
        scan.push(geometry, Some(Place::Coordinate(target)))
            .unwrap();
    }
    let degrees = scan
        .iter()
        .map(|geometry| {
            let coordinate = coordinate(geometry.place());
            assert_eq!(coordinate.unit(), DEGREE);
            *coordinate.value()
        })
        .collect::<Vec<_>>();
    assert_eq!(degrees[0], 0.0);
    assert!((degrees[1] - 90.0).abs() < 1e-12);
    assert_eq!(degrees[2], 180.0);
    assert_eq!(degrees[3], 90.0);
}

#[test]
fn geometries_must_lie_on_the_order_and_fit_its_topology() {
    let mut steps = ordered(Order::Step(Vec::new()));
    let mut scan = ordered(coordinates(ANGSTROM));
    let geometry = || Geometry::new(positions(0.0));
    let at = |value, unit| Some(Place::Coordinate(Quantity::new(value, unit)));
    assert_eq!(
        steps.push(geometry(), at(1.5, ANGSTROM)),
        Err(StructureSetError::PlaceMismatch { geometry: 0 })
    );
    assert_eq!(
        scan.push(geometry(), Some(Place::Step(0))),
        Err(StructureSetError::PlaceMismatch { geometry: 0 })
    );
    // An angle on a distance order.
    assert!(matches!(
        scan.push(geometry(), at(1.0, RADIAN)),
        Err(StructureSetError::Unit(_))
    ));
    assert_eq!(
        scan.push(geometry(), at(f64::NAN, ANGSTROM)),
        Err(StructureSetError::NonFiniteOrderValue { geometry: 0 })
    );
    // Too few atoms for the topology.
    assert!(matches!(
        steps.push(Geometry::new(Positions::zeros(2)), Some(Place::Step(0))),
        Err(StructureSetError::Geometry(GeometryError::Array(_)))
    ));
    assert!(steps.is_empty() && scan.is_empty());
    // An order has one value per geometry.
    let mut collected = StructureSet::from_geometries(chain(3), [geometry(), geometry()]).unwrap();
    assert_eq!(
        collected.set_order(Some(Order::Step(vec![0]))),
        Err(StructureSetError::OrderLengthMismatch {
            expected: 2,
            actual: 1
        })
    );
    assert_eq!(collected.order(), None);
}

#[test]
fn geometry_views_borrow_as_structures_with_their_place() {
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    let mut scan = ordered(coordinates(ANGSTROM));
    let geometry = with_energy(Geometry::new(positions(0.0)).with_cell(cell), -3.0);
    let place = Place::Coordinate(Quantity::new(1.25, ANGSTROM));
    scan.push(geometry, Some(place)).unwrap();
    let geometry = scan.get(0).unwrap();
    let view = geometry.as_structure_view();
    assert!(matches!(
        view.kind(),
        StructureKind::InSet(context)
            if context.place() == Some(place) && context.membership().is_none()
    ));
    assert_eq!(
        AsStructureView::as_structure_view(&geometry).atom_count(),
        3
    );
    assert_eq!(view.cell(), Some(&cell));
    assert!(view.occupancies().is_none() && view.b_factors().is_none());
    assert!(Arc::ptr_eq(
        &view.shared_topology(),
        &scan.shared_topology()
    ));
    // A geometry materialized as a structure keeps its geometry, not its place.
    let structure = geometry.to_structure();
    assert_eq!(
        structure.as_structure_view().kind(),
        StructureKind::Standalone
    );
    assert_eq!(structure.positions(), geometry.positions());
    assert_eq!(structure.cell(), Some(&cell));
    assert_eq!(structure.properties(), geometry.properties());
}

#[test]
fn sets_from_structures_take_a_step_order_and_require_one_layout() {
    let topology = chain(3);
    let structures = (0..3)
        .map(|index| {
            let mut structure =
                Structure::new(Arc::clone(&topology), positions(0.01 * index as f64)).unwrap();
            structure.set_occupancies(Some(vec![Some(1.0); 3])).unwrap();
            structure
                .set_b_factors(Some(Quantity::new(vec![Some(0.1); 3], SQUARE_NANOMETER)))
                .unwrap();
            structure
        })
        .collect::<Vec<_>>();
    let mut set = StructureSet::from_structures(structures.clone()).unwrap();
    assert_eq!(set.order(), None);
    set.set_order(Some(Order::Step(vec![0, 1, 2]))).unwrap();
    assert!(Arc::ptr_eq(&set.shared_topology(), &topology));
    for (index, geometry) in set.iter().enumerate() {
        assert_eq!(geometry.place(), Some(Place::Step(index as u64)));
        assert_eq!(geometry.positions(), structures[index].positions());
        // Geometries keep the structures' crystallographic site values.
        assert_eq!(geometry.occupancies(), structures[index].occupancies());
        assert_eq!(geometry.b_factors(), structures[index].b_factors());
    }
    // Structures read apart have equal layouts and form a set too.
    let apart = Structure::new(chain(3), positions(0.5)).unwrap();
    let set = StructureSet::from_structures([structures[0].clone(), apart]).unwrap();
    assert_eq!(set.len(), 2);
    let other = Structure::new(chain(4), Positions::zeros(4)).unwrap();
    assert_eq!(
        StructureSet::from_structures([structures[0].clone(), other]).unwrap_err(),
        StructureSetError::TopologyMismatch
    );
    assert_eq!(
        StructureSet::from_structures(Vec::new()).unwrap_err(),
        StructureSetError::EmptySource
    );
}

#[test]
fn geometries_move_on_the_order_and_are_selected_with_their_places() {
    let mut set = history();
    set.set_place(1, Some(Place::Step(6))).unwrap();
    assert_eq!(
        set.set_place(1, Some(Place::Coordinate(Quantity::new(1.0, ANGSTROM)))),
        Err(StructureSetError::PlaceMismatch { geometry: 1 })
    );
    let origin = Quantity::new([Point3::origin(); 3], NANOMETER);
    set.get_mut(1).unwrap().set_positions(origin).unwrap();
    assert_eq!(set.get(1).unwrap().place(), Some(Place::Step(6)));
    // A removed geometry claims no place.
    let removed = set.remove(0).unwrap();
    assert_eq!(removed.positions(), &positions(0.0));
    assert_eq!(set.order(), Some(&Order::Step(vec![6, 10, 15])));
    // A replacement takes the place of the geometry it replaces.
    let replaced = set.replace(0, Geometry::new(positions(0.2))).unwrap();
    assert_eq!(
        replaced.positions().values().value(),
        &[Point3::origin(); 3]
    );
    assert_eq!(set.get(0).unwrap().place(), Some(Place::Step(6)));
    set.set_place(0, Some(Place::Step(7))).unwrap();
    assert_eq!(
        set.remove(3).unwrap_err(),
        StructureSetError::IndexOutOfRange { index: 3, len: 3 }
    );

    // Selection keeps each geometry's step, in the order asked for.
    let selected = set.select([2, 0, 0]).unwrap();
    assert_eq!(places(&selected), steps(&[15, 7, 7]));
    assert_eq!(selected.order(), Some(&Order::Step(vec![15, 7, 7])));
    assert!(set.select([5]).is_err());

    // Replacing every position keeps the rest.
    let replacement = (0..3).map(|_| positions(0.3)).collect::<Vec<_>>();
    set.replace_positions(replacement).unwrap();
    assert!(set
        .iter()
        .all(|geometry| geometry.positions() == &positions(0.3)));
    assert_eq!(places(&set), steps(&[7, 10, 15]));
    assert_eq!(
        set.replace_positions([positions(0.0)]).unwrap_err(),
        StructureSetError::GeometryCountMismatch {
            expected: 3,
            actual: 1
        }
    );

    // A geometry moves along a coordinate order, in its unit.
    let mut scan = ordered(coordinates(DEGREE));
    let ten = Place::Coordinate(Quantity::new(10.0, DEGREE));
    scan.push(Geometry::new(positions(0.0)), Some(ten)).unwrap();
    let half_turn = Place::Coordinate(Quantity::new(std::f64::consts::PI, RADIAN));
    scan.set_place(0, Some(half_turn)).unwrap();
    assert!(matches!(
        scan.set_place(0, Some(Place::Coordinate(Quantity::new(1.0, ANGSTROM)))),
        Err(StructureSetError::Unit(_))
    ));
    assert_eq!(
        scan.set_place(0, Some(Place::Step(1))).unwrap_err(),
        StructureSetError::PlaceMismatch { geometry: 0 }
    );
    let coordinate = coordinate(scan.get(0).unwrap().place());
    assert_eq!(coordinate.unit(), DEGREE);
    assert!((coordinate.value() - 180.0).abs() < 1e-12);
}

#[test]
fn subsets_keep_places_and_drop_values_of_the_whole_structure() {
    let mut set = history();
    set.properties_mut()
        .insert(
            PropertyKey::new("method").unwrap(),
            PropertyValue::String("GFN2-xTB".into()),
        )
        .unwrap();
    let topology = set.shared_topology();
    let atoms = AtomSelection::from_atoms(
        Arc::clone(&topology),
        topology.atom_ids()[..2].iter().copied(),
    )
    .unwrap();
    let subset = set.subset(&atoms).unwrap();
    assert_eq!(subset.topology().atom_count(), 2);
    assert!(subset.properties().is_empty());
    assert_eq!(subset.order(), set.order());
    for (source, geometry) in set.iter().zip(&subset) {
        assert_eq!(geometry.place(), source.place());
        assert!(geometry.properties().owner().is_empty());
    }
    // Perception installs a snapshot of the same layout.
    let mut perceived = set.clone();
    perceived.perceive().unwrap();
    assert!(perceived.topology().shares_layout(set.topology()));
    assert_eq!(perceived.len(), set.len());
    assert_eq!(perceived.order(), set.order());
    let expected = set
        .iter()
        .map(|geometry| geometry.payload().clone())
        .collect::<Vec<_>>();
    assert_eq!(set.into_geometries(), expected);
}

#[test]
fn ordered_sets_superpose_and_measure_rmsd_geometry_by_geometry() {
    let mut set = history();
    let topology = set.shared_topology();
    let atoms = AtomSelection::all(Arc::clone(&topology));
    let rmsd = set.rmsd(0, &atoms, AlignmentOptions::default()).unwrap();
    assert_eq!(rmsd.value().len(), 4);
    assert_eq!(rmsd.value()[0], 0.0);
    assert!(rmsd.value()[3] > rmsd.value()[1]);
    let report = set
        .superpose(0, &atoms, AlignmentOptions::default())
        .unwrap();
    assert_eq!(report.len(), 4);
    assert_eq!(report.reference_index(), Some(0));
    // Superposition keeps every step and value.
    assert_eq!(places(&set), steps(&[0, 5, 10, 15]));
    assert!(set
        .iter()
        .all(|geometry| geometry.properties().owner().get(&energy()).is_some()));
    let aligned = set
        .aligned_rmsd(0, &atoms, &atoms, Default::default())
        .unwrap();
    for (direct, fitted) in set
        .rmsd(0, &atoms, AlignmentOptions::default())
        .unwrap()
        .value()
        .iter()
        .zip(aligned.value())
    {
        assert!((direct - fitted).abs() < 1e-9);
    }
}
