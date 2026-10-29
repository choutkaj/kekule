//! Ordered series of structures on a step or coordinate axis: optimization
//! histories, relaxed scans, and paths.

use std::sync::Arc;

use kekule::alignment::AlignmentOptions;
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{
    AsModelView, Model, Positions, Series, SeriesAxis, SeriesError, SeriesItem, StructureError,
    StructureKind,
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

fn with_energy(mut item: SeriesItem, value: f64) -> SeriesItem {
    item.properties_mut()
        .owner_mut()
        .insert(
            energy(),
            PropertyValue::real(value, KILOJOULE_PER_MOLE).unwrap(),
        )
        .unwrap();
    item
}

/// A history of four optimization steps with energies.
fn history() -> Series {
    let mut series = Series::new(chain(3), SeriesAxis::Step);
    for step in 0..4 {
        let item = SeriesItem::at_step(positions(0.01 * step as f64), step * 5);
        series
            .push(with_energy(item, -100.0 - step as f64))
            .unwrap();
    }
    series
}

fn real(value: &PropertyValue) -> f64 {
    match value {
        PropertyValue::Real { value, .. } => *value,
        other => panic!("not a real value: {other:?}"),
    }
}

#[test]
fn a_history_keeps_its_items_steps_and_values_in_order() {
    let series = history();
    assert_eq!(series.axis(), SeriesAxis::Step);
    assert_eq!(series.len(), 4);
    let steps = series.iter().map(|item| item.step()).collect::<Vec<_>>();
    assert_eq!(steps, [Some(0), Some(5), Some(10), Some(15)]);
    assert!(series.iter().all(|item| item.coordinate().is_none()));
    let energies = series
        .iter()
        .map(|item| real(item.properties().owner().get(&energy()).unwrap()))
        .collect::<Vec<_>>();
    assert_eq!(energies, [-100.0, -101.0, -102.0, -103.0]);
    let item = series.get(2).unwrap();
    assert_eq!(item.index(), 2);
    assert_eq!(item.payload().step(), Some(10));
    assert!(series.get(4).is_none());
    assert_eq!(series.iter().len(), 4);
}

#[test]
fn a_scan_stores_coordinates_in_its_axis_unit_in_any_order() {
    let mut scan = Series::new(chain(3), SeriesAxis::Coordinate(DEGREE));
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
        let item = SeriesItem::at_coordinate(positions(0.01 * index as f64), target).unwrap();
        assert_eq!(item.coordinate(), Some(target));
        scan.push(item).unwrap();
    }
    let degrees = scan
        .iter()
        .map(|item| {
            let coordinate = item.coordinate().unwrap();
            assert_eq!(coordinate.unit(), DEGREE);
            *coordinate.value()
        })
        .collect::<Vec<_>>();
    assert!((degrees[1] - 90.0).abs() < 1e-12);
    assert_eq!(degrees[2], 180.0);
    assert_eq!(degrees[3], 90.0);
    assert!(scan.iter().all(|item| item.step().is_none()));
}

#[test]
fn items_must_lie_on_the_series_axis_and_fit_its_topology() {
    let mut steps = Series::new(chain(3), SeriesAxis::Step);
    let mut scan = Series::new(chain(3), SeriesAxis::Coordinate(ANGSTROM));
    let coordinate = SeriesItem::at_coordinate(positions(0.0), Quantity::new(1.5, ANGSTROM));
    assert_eq!(
        steps.push(coordinate.unwrap()),
        Err(SeriesError::AxisMismatch {
            axis: SeriesAxis::Step
        })
    );
    assert_eq!(
        scan.push(SeriesItem::at_step(positions(0.0), 0)),
        Err(SeriesError::AxisMismatch {
            axis: SeriesAxis::Coordinate(ANGSTROM)
        })
    );
    // An angle on a distance axis.
    let angle = SeriesItem::at_coordinate(positions(0.0), Quantity::new(1.0, RADIAN)).unwrap();
    assert!(matches!(scan.push(angle), Err(SeriesError::Unit(_))));
    assert_eq!(
        SeriesItem::at_coordinate(positions(0.0), Quantity::new(f64::NAN, ANGSTROM)),
        Err(StructureError::NonFiniteCoordinate)
    );
    // Too few atoms for the topology.
    let short = Positions::zeros(2);
    assert!(matches!(
        steps.push(SeriesItem::at_step(short, 0)),
        Err(SeriesError::Structure(StructureError::Array(_)))
    ));
    assert!(steps.is_empty() && scan.is_empty());
    let from_items = Series::from_items(
        chain(3),
        SeriesAxis::Step,
        [
            SeriesItem::at_step(positions(0.0), 0),
            SeriesItem::at_coordinate(positions(0.0), Quantity::new(1.0, ANGSTROM)).unwrap(),
        ],
    );
    assert!(matches!(
        from_items,
        Err(SeriesError::AxisMismatch { .. })
    ));
}

#[test]
fn item_views_borrow_as_models_with_their_place() {
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(2.0, 2.0, 2.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    let mut scan = Series::new(chain(3), SeriesAxis::Coordinate(ANGSTROM));
    let item = SeriesItem::at_coordinate(positions(0.0), Quantity::new(1.25, ANGSTROM))
        .unwrap()
        .with_cell(cell);
    scan.push(with_energy(item, -3.0)).unwrap();
    let item = scan.get(0).unwrap();
    let view = item.as_model_view();
    assert!(matches!(
        view.kind(),
        StructureKind::Item {
            step: None,
            coordinate: Some(coordinate),
            ..
        } if coordinate == Quantity::new(1.25, ANGSTROM)
    ));
    assert_eq!(AsModelView::as_model_view(&item).atom_count(), 3);
    assert_eq!(view.cell(), Some(&cell));
    assert!(view.occupancies().is_none() && view.b_factors().is_none());
    assert!(Arc::ptr_eq(&view.shared_topology(), &scan.shared_topology()));
    // An item materialized as a model keeps its structure, not its place.
    let model = item.to_model();
    assert_eq!(model.as_model_view().kind(), StructureKind::Model);
    assert_eq!(model.positions(), item.positions());
    assert_eq!(model.cell(), Some(&cell));
    assert_eq!(model.properties(), item.properties());
}

#[test]
fn series_from_models_number_them_and_require_one_layout() {
    let topology = chain(3);
    let models = (0..3)
        .map(|index| {
            let mut model = Model::new(Arc::clone(&topology), positions(0.01 * index as f64))
                .unwrap();
            model.set_occupancies(Some(vec![Some(1.0); 3])).unwrap();
            model
                .set_b_factors(Some(Quantity::new(vec![Some(0.1); 3], SQUARE_NANOMETER)))
                .unwrap();
            model
        })
        .collect::<Vec<_>>();
    let series = Series::from_models(models.clone()).unwrap();
    assert_eq!(series.axis(), SeriesAxis::Step);
    assert!(Arc::ptr_eq(&series.shared_topology(), &topology));
    for (index, item) in series.iter().enumerate() {
        assert_eq!(item.step(), Some(index as u64));
        assert_eq!(item.positions(), models[index].positions());
        // Series items carry no crystallographic site values.
        assert!(item.as_model_view().occupancies().is_none());
    }
    // Models read apart have equal layouts and form a series too.
    let apart = Model::new(chain(3), positions(0.5)).unwrap();
    let series = Series::from_models([models[0].clone(), apart]).unwrap();
    assert_eq!(series.len(), 2);
    let other = Model::new(chain(4), Positions::zeros(4)).unwrap();
    assert_eq!(
        Series::from_models([models[0].clone(), other]).unwrap_err(),
        SeriesError::TopologyMismatch
    );
    assert_eq!(
        Series::from_models(Vec::new()).unwrap_err(),
        SeriesError::EmptySource
    );
}

#[test]
fn items_are_edited_on_the_axis_and_selected_with_their_places() {
    let mut series = history();
    {
        let mut item = series.get_mut(1).unwrap();
        item.set_step(6).unwrap();
        assert_eq!(
            item.set_coordinate(Quantity::new(1.0, ANGSTROM)),
            Err(SeriesError::AxisMismatch {
                axis: SeriesAxis::Step
            })
        );
        item.set_positions(Quantity::new([Point3::origin(); 3], NANOMETER))
            .unwrap();
    }
    assert_eq!(series.get(1).unwrap().step(), Some(6));
    let removed = series.remove(0).unwrap();
    assert_eq!(removed.step(), Some(0));
    let replaced = series
        .replace(0, SeriesItem::at_step(positions(0.2), 7))
        .unwrap();
    assert_eq!(replaced.step(), Some(6));
    assert_eq!(
        series
            .replace(0, SeriesItem::at_coordinate(positions(0.0), Quantity::new(1.0, ANGSTROM)).unwrap())
            .unwrap_err(),
        SeriesError::AxisMismatch {
            axis: SeriesAxis::Step
        }
    );
    assert_eq!(
        series.remove(3).unwrap_err(),
        SeriesError::IndexOutOfRange { index: 3, len: 3 }
    );

    // Selection keeps each item's step, in the order asked for.
    let selected = series.select([2, 0, 0]).unwrap();
    let steps = selected.iter().map(|item| item.step()).collect::<Vec<_>>();
    assert_eq!(steps, [Some(15), Some(7), Some(7)]);
    assert_eq!(selected.axis(), series.axis());
    assert!(series.select([5]).is_err());

    // Replacing every position keeps the rest.
    let replacement = (0..3).map(|_| positions(0.3)).collect::<Vec<_>>();
    series.replace_positions(replacement).unwrap();
    assert!(series.iter().all(|item| item.positions() == &positions(0.3)));
    assert_eq!(
        series.replace_positions([positions(0.0)]).unwrap_err(),
        SeriesError::ItemCountMismatch {
            expected: 3,
            actual: 1
        }
    );

    // A coordinate item moves along a coordinate axis, in its unit.
    let mut scan = Series::new(chain(3), SeriesAxis::Coordinate(DEGREE));
    scan.push(SeriesItem::at_coordinate(positions(0.0), Quantity::new(10.0, DEGREE)).unwrap())
        .unwrap();
    let mut item = scan.get_mut(0).unwrap();
    item.set_coordinate(Quantity::new(std::f64::consts::PI, RADIAN))
        .unwrap();
    assert!(matches!(
        item.set_coordinate(Quantity::new(1.0, ANGSTROM)),
        Err(SeriesError::Unit(_))
    ));
    assert_eq!(item.set_step(1).unwrap_err(), SeriesError::AxisMismatch {
        axis: SeriesAxis::Coordinate(DEGREE)
    });
    let coordinate = scan.get(0).unwrap().coordinate().unwrap();
    assert_eq!(coordinate.unit(), DEGREE);
    assert!((coordinate.value() - 180.0).abs() < 1e-12);
}

#[test]
fn subsets_keep_places_and_drop_values_of_the_whole_structure() {
    let mut series = history();
    series
        .properties_mut()
        .insert(PropertyKey::new("method").unwrap(), PropertyValue::String("GFN2-xTB".into()))
        .unwrap();
    let topology = series.shared_topology();
    let atoms = AtomSelection::from_atoms(
        Arc::clone(&topology),
        topology.atom_ids()[..2].iter().copied(),
    )
    .unwrap();
    let subset = series.subset(&atoms).unwrap();
    assert_eq!(subset.topology().atom_count(), 2);
    assert!(subset.properties().is_empty());
    for (source, item) in series.iter().zip(&subset) {
        assert_eq!(item.step(), source.step());
        assert!(item.properties().owner().is_empty());
    }
    // Perception installs a snapshot of the same layout.
    let mut perceived = series.clone();
    perceived.perceive().unwrap();
    assert!(perceived.topology().shares_layout(series.topology()));
    assert_eq!(perceived.len(), series.len());
    let (topology, properties, axis, items) = series.into_parts();
    assert_eq!(topology.atom_count(), 3);
    assert!(!properties.is_empty());
    assert_eq!(axis, SeriesAxis::Step);
    assert_eq!(items.len(), 4);
}

#[test]
fn series_superpose_and_measure_rmsd_item_by_item() {
    let mut series = history();
    let topology = series.shared_topology();
    let atoms = AtomSelection::all(Arc::clone(&topology));
    let rmsd = series.rmsd(0, &atoms, AlignmentOptions::default()).unwrap();
    assert_eq!(rmsd.value().len(), 4);
    assert_eq!(rmsd.value()[0], 0.0);
    assert!(rmsd.value()[3] > rmsd.value()[1]);
    let report = series
        .superpose(0, &atoms, AlignmentOptions::default())
        .unwrap();
    assert_eq!(report.len(), 4);
    assert_eq!(report.reference_index(), Some(0));
    // Superposition keeps every step and value.
    let steps = series.iter().map(|item| item.step()).collect::<Vec<_>>();
    assert_eq!(steps, [Some(0), Some(5), Some(10), Some(15)]);
    let aligned = series
        .aligned_rmsd(0, &atoms, &atoms, Default::default())
        .unwrap();
    for (direct, fitted) in series
        .rmsd(0, &atoms, AlignmentOptions::default())
        .unwrap()
        .value()
        .iter()
        .zip(aligned.value())
    {
        assert!((direct - fitted).abs() < 1e-9);
    }
}
