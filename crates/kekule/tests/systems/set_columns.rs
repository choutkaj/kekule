//! A set's columns: one row per geometry, kept with the geometries by every
//! operation, and read as numbers beside the order.

use std::sync::Arc;

use kekule::math::Point3;
use kekule::native::StructureSetDescription;
use kekule::properties::{
    PropertyColumn, PropertyKey, PropertyTable, PropertyValue, PropertyValueRef,
};
use kekule::structure::{
    ColumnRef, Ensemble, Geometry, GeometryIndex, GeometryView, Order, Place, Positions,
    Requirement, Structure, StructureSet, StructureSetError, StructureSetRead,
    StructureSetReadError, TimeAxis,
};
use kekule::topology::{AtomSelection, Topology, TopologyAtomIndex};
use kekule::traj::{Processed, ProcessingSteps};
use kekule::units::{
    Quantity, DEGREE, DIMENSIONLESS, KILOCALORIE_PER_MOLE, KILOJOULE_PER_MOLE, NANOMETER,
    PICOSECOND,
};

fn water() -> Arc<Topology> {
    kekule::smiles::to_topology("O").unwrap()
}

fn geometry(x: f64) -> Geometry {
    Geometry::new(Positions::new(Quantity::new(vec![Point3::new(x, 0.0, 0.0)], NANOMETER)).unwrap())
}

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

fn energy() -> PropertyKey {
    key("energy")
}

/// Four geometries whose energy column is the geometry's number, in kJ/mol.
fn numbered(count: usize) -> StructureSet {
    let mut set =
        StructureSet::from_geometries(water(), (0..count).map(|index| geometry(index as f64)))
            .unwrap();
    set.columns_mut()
        .insert(
            energy(),
            PropertyColumn::Real {
                unit: KILOJOULE_PER_MOLE,
                values: (0..count).map(|index| Some(index as f64)).collect(),
            },
        )
        .unwrap();
    set
}

/// The energy cells of `set`, missing ones as `None`.
fn energies(set: &StructureSet) -> Vec<Option<f64>> {
    set.column_values(&ColumnRef::Key(energy()))
        .unwrap()
        .into_value()
}

#[test]
fn every_constructor_gives_one_row_per_geometry() {
    let set = StructureSet::from_geometries(water(), [geometry(0.0), geometry(1.0)]).unwrap();
    assert_eq!(set.columns().len(), 2);
    assert!(!set.columns().has_data());

    let structures = [0.0, 1.0, 2.0]
        .map(|x| Structure::from_parts(water(), geometry(x)).unwrap())
        .to_vec();
    assert_eq!(
        StructureSet::from_structures(structures)
            .unwrap()
            .columns()
            .len(),
        3
    );
    assert_eq!(StructureSet::new(water()).columns().len(), 0);
}

#[test]
fn a_pushed_geometry_has_missing_cells_and_a_replaced_one_keeps_its_row() {
    let mut set = numbered(2);
    assert_eq!(set.push(geometry(5.0), None).unwrap(), 2);
    assert_eq!(energies(&set), [Some(0.0), Some(1.0), None]);
    set.replace(1, geometry(7.0)).unwrap();
    assert_eq!(energies(&set), [Some(0.0), Some(1.0), None]);
}

#[test]
fn a_removed_geometry_takes_its_row() {
    let mut set = numbered(4);
    set.remove(1).unwrap();
    assert_eq!(energies(&set), [Some(0.0), Some(2.0), Some(3.0)]);
    assert_eq!(set.columns().len(), 3);

    // A column whose last value goes is not stored.
    let mut set = StructureSet::from_geometries(water(), [geometry(0.0), geometry(1.0)]).unwrap();
    set.columns_mut()
        .set_value(
            key("flag"),
            GeometryIndex::new(0),
            Some(PropertyValue::Bool(true)),
        )
        .unwrap();
    set.remove(0).unwrap();
    assert!(set.columns().get(&key("flag")).is_none());
}

#[test]
fn a_selection_picks_rows_in_the_requested_order() {
    let set = numbered(4);
    let selected = set.select([3, 0, 3, 1]).unwrap();
    assert_eq!(
        energies(&selected),
        [Some(3.0), Some(0.0), Some(3.0), Some(1.0)]
    );

    let read = StructureSet::from_reader(&set, [2, 2, 0]).unwrap();
    assert_eq!(energies(&read), [Some(2.0), Some(2.0), Some(0.0)]);

    // A selection leaving out a column's only value drops the column.
    let mut sparse = numbered(3);
    sparse
        .columns_mut()
        .set_value(
            key("flag"),
            GeometryIndex::new(2),
            Some(PropertyValue::Bool(true)),
        )
        .unwrap();
    assert!(sparse
        .select([0, 1])
        .unwrap()
        .columns()
        .get(&key("flag"))
        .is_none());
}

#[test]
fn an_atom_subset_and_a_processed_reader_keep_the_columns() {
    let set = numbered(3);
    let atoms =
        AtomSelection::from_indices(set.shared_topology(), [TopologyAtomIndex::new(0)]).unwrap();
    let subset = set.subset(&atoms).unwrap();
    assert_eq!(subset.topology().atom_count(), 1);
    assert_eq!(subset.columns(), set.columns());

    let processed = Processed::new(&set, ProcessingSteps::default().with_subset(atoms)).unwrap();
    assert_eq!(processed.columns(), Some(set.columns()));
    let loaded = StructureSet::from_reader(&processed, [1, 2]).unwrap();
    assert_eq!(energies(&loaded), [Some(1.0), Some(2.0)]);
}

#[test]
fn ensemble_operations_keep_the_rows_of_their_geometries() {
    let ensemble = Ensemble::uniform(numbered(4));
    let selected = ensemble.select([2, 0]).unwrap();
    assert_eq!(energies(selected.set()), [Some(2.0), Some(0.0)]);
    let collapsed = ensemble.collapse();
    assert_eq!(energies(collapsed.set()), [Some(0.0)]);
    let mut grown = Ensemble::uniform(numbered(2));
    grown.push(geometry(9.0), None).unwrap();
    assert_eq!(energies(grown.set()), [Some(0.0), Some(1.0), None]);
    grown.remove(1).unwrap();
    assert_eq!(energies(grown.set()), [Some(0.0), None]);
    grown
        .columns_mut()
        .set_value(
            energy(),
            GeometryIndex::new(1),
            Some(PropertyValue::Real {
                value: 4.0,
                unit: KILOJOULE_PER_MOLE,
            }),
        )
        .unwrap();
    assert_eq!(energies(grown.set()), [Some(0.0), Some(4.0)]);
    assert!(grown
        .columns_mut()
        .insert(key("short"), PropertyColumn::Int(vec![Some(1)]))
        .is_err());
}

/// A reader of its own whose columns have one row too many.
#[derive(Debug)]
struct Overlong {
    set: StructureSet,
    columns: PropertyTable<GeometryIndex>,
}

impl StructureSetRead for Overlong {
    fn topology(&self) -> &Arc<Topology> {
        StructureSetRead::topology(&self.set)
    }

    fn len(&self) -> usize {
        self.set.len()
    }

    fn order(&self) -> Option<&Order> {
        self.set.order()
    }

    fn columns(&self) -> Option<&PropertyTable<GeometryIndex>> {
        Some(&self.columns)
    }

    fn peek(&self, index: usize) -> Option<GeometryView<'_, dyn StructureSetRead>> {
        self.set.peek(index)
    }

    fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        self.set.read(index)
    }
}

#[test]
fn a_reader_whose_columns_miss_rows_fails_to_load() {
    let reader = Overlong {
        set: numbered(2),
        columns: numbered(3).columns().clone(),
    };
    let error = StructureSet::from_reader(&reader, [0]).unwrap_err();
    assert_eq!(
        error,
        StructureSetReadError::ColumnRows {
            rows: 3,
            geometries: 2
        }
    );
    assert_eq!(
        error.to_string(),
        "The columns of the reader have 3 rows for 2 geometries."
    );
}

#[test]
fn column_values_read_the_order_and_numeric_columns() {
    let mut set = numbered(3);
    set.columns_mut()
        .insert(
            key("replica"),
            PropertyColumn::Int(vec![Some(4), None, Some(-1)]),
        )
        .unwrap();
    set.columns_mut()
        .insert(
            key("stage"),
            PropertyColumn::String(vec![Some("nvt".into()), None, None]),
        )
        .unwrap();

    let replica = set.column_values(&ColumnRef::Key(key("replica"))).unwrap();
    assert_eq!(replica.unit(), DIMENSIONLESS);
    assert_eq!(replica.value(), &[Some(4.0), None, Some(-1.0)]);
    let energy = set.column_values(&ColumnRef::Key(energy())).unwrap();
    assert_eq!(energy.unit(), KILOJOULE_PER_MOLE);
    assert_eq!(
        set.column_values(&ColumnRef::Key(key("stage"))),
        Err(StructureSetError::NonNumericColumn { key: key("stage") })
    );
    let absent = set.column_values(&ColumnRef::Key(key("bias"))).unwrap_err();
    assert_eq!(
        absent,
        StructureSetError::MissingColumn { key: key("bias") }
    );
    assert_eq!(
        absent.to_string(),
        "The structure set has no column \"bias\"."
    );

    // Without an order there are no order values.
    for (column, requirement) in [
        (ColumnRef::Time, Requirement::TimeOrder),
        (ColumnRef::Step, Requirement::Steps),
        (ColumnRef::Coordinate, Requirement::CoordinateOrder),
    ] {
        assert_eq!(set.column_values(&column), Err(requirement.into()));
    }
    assert_eq!(set.steps(), None);

    // A time order with times and steps.
    let times = Quantity::new(vec![0.0, 2.0, 4.0], PICOSECOND);
    set.set_order(Some(Order::Time(
        TimeAxis::new(Some(times), Some(vec![0, 1000, 2000])).unwrap(),
    )))
    .unwrap();
    let time = set.column_values(&ColumnRef::Time).unwrap();
    assert_eq!(time.unit(), PICOSECOND);
    assert_eq!(time.value(), &[Some(0.0), Some(2.0), Some(4.0)]);
    let steps = set.column_values(&ColumnRef::Step).unwrap();
    assert_eq!(steps.value(), &[Some(0.0), Some(1000.0), Some(2000.0)]);
    assert_eq!(set.steps(), Some(&[0, 1000, 2000][..]));
    assert_eq!(
        set.column_values(&ColumnRef::Coordinate),
        Err(Requirement::CoordinateOrder.into())
    );

    // A time order without times.
    set.set_order(Some(Order::Time(TimeAxis::default())))
        .unwrap();
    assert_eq!(
        set.column_values(&ColumnRef::Time),
        Err(Requirement::Times.into())
    );
    assert_eq!(set.steps(), None);

    // A step order and a coordinate order.
    set.set_order(Some(Order::Step(vec![5, 6, 7]))).unwrap();
    assert_eq!(set.steps(), Some(&[5, 6, 7][..]));
    set.set_order(Some(Order::Coordinate {
        unit: DEGREE,
        values: vec![-60.0, 0.0, 60.0],
    }))
    .unwrap();
    let coordinate = set.column_values(&ColumnRef::Coordinate).unwrap();
    assert_eq!(coordinate.unit(), DEGREE);
    assert_eq!(coordinate.value(), &[Some(-60.0), Some(0.0), Some(60.0)]);
    assert_eq!(
        set.column_values(&ColumnRef::Step),
        Err(Requirement::Steps.into())
    );
}

#[test]
fn columns_mut_preserves_the_row_count() {
    let mut set = numbered(2);
    let error = set
        .columns_mut()
        .insert(key("short"), PropertyColumn::Int(vec![Some(1)]))
        .unwrap_err();
    assert!(matches!(
        &error,
        StructureSetError::Column(source)
            if matches!(**source, kekule::properties::PropertyError::LengthMismatch(_))
    ));
    assert_eq!(
        error.to_string(),
        "Cannot change the columns of the structure set."
    );
    assert!(set
        .columns_mut()
        .set_value(key("x"), GeometryIndex::new(2), Some(PropertyValue::Int(1)))
        .is_err());
    assert_eq!(set.columns().len(), 2);
    set.columns_mut()
        .clear_value(energy(), GeometryIndex::new(0))
        .unwrap();
    assert_eq!(energies(&set), [None, Some(1.0)]);
    assert!(set.columns_mut().remove(&energy()).unwrap().is_some());
    assert!(!set.columns().has_data());
}

#[test]
fn a_description_keeps_the_columns_and_extends_them_with_missing_cells() {
    let mut set = numbered(2);
    let times = Quantity::new(vec![0.0, 1.0], PICOSECOND);
    set.set_order(Some(Order::Time(TimeAxis::new(Some(times), None).unwrap())))
        .unwrap();
    let description = StructureSetDescription::of(&set);
    assert_eq!(description.columns(), set.columns());

    let place = Place::Time {
        time: Some(Quantity::new(2.0, PICOSECOND)),
        step: None,
    };
    set.push(geometry(2.0), Some(place)).unwrap();
    let grown = description.extended_to(3, set.order().cloned()).unwrap();
    assert_eq!(grown.columns(), set.columns());

    let error = StructureSetDescription::new(water(), 3)
        .with_columns(numbered(2).columns().clone())
        .unwrap_err();
    assert!(matches!(error, StructureSetError::GeometryCountMismatch(_)));
}

/// `geometry(x)` with its own energy, in kcal/mol.
fn with_energy(x: f64, energy: f64) -> Geometry {
    let mut geometry = geometry(x);
    geometry
        .properties_mut()
        .owner_mut()
        .insert(
            self::energy(),
            PropertyValue::real(energy, KILOCALORIE_PER_MOLE).unwrap(),
        )
        .unwrap();
    geometry
}

fn close(left: &[Option<f64>], right: &[Option<f64>]) {
    assert_eq!(left.len(), right.len());
    for (left, right) in left.iter().zip(right) {
        match (left, right) {
            (Some(left), Some(right)) => assert!((left - right).abs() < 1.0e-9, "{left} {right}"),
            _ => assert_eq!(left, right),
        }
    }
}

#[test]
fn every_constructor_moves_owner_values_into_the_rows() {
    let owned = || [geometry(0.0), with_energy(1.0, 2.0), with_energy(2.0, 3.0)];

    // Pushed into a set whose column is in kJ/mol, the values convert to it.
    let mut set = numbered(1);
    for geometry in owned().into_iter().skip(1) {
        set.push(geometry, None).unwrap();
    }
    close(
        &energies(&set),
        &[Some(0.0), Some(2.0 * 4.184), Some(3.0 * 4.184)],
    );
    assert!(set
        .iter()
        .all(|geometry| geometry.properties().owner().is_empty()));

    // A new key makes a column in the first value's unit.
    let set = StructureSet::from_geometries(water(), owned()).unwrap();
    assert_eq!(
        set.columns().get(&energy()),
        Some(&PropertyColumn::Real {
            unit: KILOCALORIE_PER_MOLE,
            values: vec![None, Some(2.0), Some(3.0)],
        })
    );
    assert!(set
        .iter()
        .all(|geometry| geometry.properties().owner().is_empty()));

    let structures = owned().map(|geometry| Structure::from_parts(water(), geometry).unwrap());
    let from_structures = StructureSet::from_structures(structures).unwrap();
    assert_eq!(from_structures.columns(), set.columns());

    // A reader's geometries move their values in beside the reader's columns.
    let mut source =
        StructureSet::from_geometries(water(), [geometry(0.0), geometry(1.0)]).unwrap();
    source
        .columns_mut()
        .set_value(
            key("flag"),
            GeometryIndex::new(1),
            Some(PropertyValue::Bool(true)),
        )
        .unwrap();
    let loaded = StructureSet::from_reader(&Labelled(source), [1, 0]).unwrap();
    assert_eq!(
        loaded.columns().get(&key("flag")),
        Some(&PropertyColumn::Bool(vec![Some(true), None]))
    );
    assert_eq!(
        loaded.columns().get(&key("label")),
        Some(&PropertyColumn::Int(vec![Some(1), Some(0)]))
    );

    // Replacing merges the new geometry's values into the kept row and
    // returns the old geometry with the row it had.
    let mut set = StructureSet::from_geometries(water(), owned()).unwrap();
    let mut replacement = geometry(5.0);
    replacement
        .properties_mut()
        .owner_mut()
        .insert(key("flag"), PropertyValue::Bool(false))
        .unwrap();
    let replaced = set.replace(1, replacement).unwrap();
    assert_eq!(
        replaced.properties().owner().get(&energy()),
        Some(&PropertyValue::real(2.0, KILOCALORIE_PER_MOLE).unwrap())
    );
    assert!(replaced.properties().owner().get(&key("flag")).is_none());
    let row = set.row(1).unwrap();
    assert_eq!(row.get(&key("flag")), Some(PropertyValueRef::Bool(false)));
    assert!(row.get(&energy()).is_some());
}

#[test]
fn an_owner_value_of_another_type_than_its_column_fails_and_changes_nothing() {
    let mut set = numbered(2);
    let before = format!("{set:?}");
    let mut text = geometry(9.0);
    text.properties_mut()
        .owner_mut()
        .insert(energy(), PropertyValue::String("high".into()))
        .unwrap();
    assert_eq!(
        set.push(text.clone(), None),
        Err(StructureSetError::ColumnTypeMismatch {
            key: energy(),
            geometry: 2
        })
    );
    assert_eq!(
        set.replace(0, text.clone()),
        Err(StructureSetError::ColumnTypeMismatch {
            key: energy(),
            geometry: 0
        })
    );
    assert_eq!(format!("{set:?}"), before);
    let mut length = geometry(9.0);
    length
        .properties_mut()
        .owner_mut()
        .insert(energy(), PropertyValue::real(1.0, NANOMETER).unwrap())
        .unwrap();
    assert!(matches!(
        StructureSet::from_geometries(water(), [with_energy(0.0, 1.0), length]),
        Err(StructureSetError::ColumnTypeMismatch { geometry: 1, .. })
    ));
    assert_eq!(
        StructureSet::from_geometries(water(), [with_energy(0.0, 1.0), text])
            .unwrap_err()
            .to_string(),
        "The value \"energy\" of geometry 1 has another type or unit than the column \"energy\"."
    );

    // A reader's geometry fails at its index in the reader.
    let mut named = numbered(2);
    named
        .columns_mut()
        .insert(
            key("label"),
            PropertyColumn::String(vec![Some("first".into()), None]),
        )
        .unwrap();
    assert_eq!(
        StructureSet::from_reader(&Labelled(named), [1, 0]).unwrap_err(),
        StructureSetReadError::ColumnTypeMismatch {
            key: key("label"),
            geometry: 1
        }
    );
}

#[test]
fn geometries_leave_a_set_with_their_rows_and_come_back_to_them() {
    let set = numbered(3);
    let energy_of = |geometry: &Geometry| geometry.properties().owner().get(&energy()).cloned();
    let kilojoules = |value: f64| Some(PropertyValue::real(value, KILOJOULE_PER_MOLE).unwrap());
    let structure = set.get(2).unwrap().to_structure();
    assert_eq!(energy_of(&structure), kilojoules(2.0));
    let geometries = set.clone().into_geometries();
    for (index, geometry) in (0..).zip(&geometries) {
        assert_eq!(energy_of(geometry), kilojoules(f64::from(index)));
    }
    let again = StructureSet::from_geometries(water(), geometries).unwrap();
    assert_eq!(again.columns(), set.columns());
    let mut set = set;
    let removed = set.remove(0).unwrap();
    assert_eq!(energy_of(&removed), kilojoules(0.0));
}

#[test]
fn a_set_geometry_editor_edits_atom_tables_and_leaves_the_row() {
    let mut set = numbered(1);
    set.get_mut(0)
        .unwrap()
        .atom_properties_mut()
        .set_value(
            key("charge"),
            TopologyAtomIndex::new(0),
            Some(PropertyValue::Int(-1)),
        )
        .unwrap();
    let geometry = set.get(0).unwrap();
    assert_eq!(
        geometry
            .properties()
            .atoms()
            .value(&key("charge"), TopologyAtomIndex::new(0))
            .unwrap(),
        Some(PropertyValue::Int(-1))
    );
    assert!(geometry.properties().owner().is_empty());
    assert!(geometry.row().unwrap().get(&energy()).is_some());
}

/// A host's reader of a set's geometries, each read with its own value
/// `label`, its index, beside the set's columns.
#[derive(Debug)]
struct Labelled(StructureSet);

impl StructureSetRead for Labelled {
    fn topology(&self) -> &Arc<Topology> {
        StructureSetRead::topology(&self.0)
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn order(&self) -> Option<&Order> {
        self.0.order()
    }

    fn columns(&self) -> Option<&PropertyTable<GeometryIndex>> {
        Some(self.0.columns())
    }

    fn peek(&self, _: usize) -> Option<GeometryView<'_, dyn StructureSetRead>> {
        None
    }

    fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        let mut geometry = self.0.read(index)?;
        geometry
            .properties_mut()
            .owner_mut()
            .insert(key("label"), PropertyValue::Int(index as i64))
            .unwrap();
        Ok(geometry)
    }
}
