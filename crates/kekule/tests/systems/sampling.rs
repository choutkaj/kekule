//! The sampling record of a structure set: its columns stay real energies,
//! it cannot lose them, and every operation that keeps columns keeps it.

use std::sync::Arc;

use kekule::math::Point3;
use kekule::native::StructureSetDescription;
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{
    ColumnRef, Geometry, GeometryIndex, Order, Positions, Sampling, StructureSet,
    StructureSetError, StructureSetRead, ThermodynamicState,
};
use kekule::topology::{AtomSelection, Topology, TopologyAtomIndex};
use kekule::traj::{Processed, ProcessingSteps};
use kekule::units::{Quantity, KELVIN, KILOCALORIE_PER_MOLE, KILOJOULE_PER_MOLE, NANOMETER};

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

fn water() -> Arc<Topology> {
    kekule::smiles::to_topology("O").unwrap()
}

/// Four water geometries at steps 0, 10, 20, 30, with a potential energy
/// in every row and a bias from step 20 on.
fn biased_set() -> StructureSet {
    let topology = water();
    let geometries = (0..4).map(|index| {
        let points = vec![Point3::new(f64::from(index), 0.0, 0.0); topology.atom_count()];
        Geometry::new(Positions::new(Quantity::new(points, NANOMETER)).unwrap())
    });
    let mut set = StructureSet::from_geometries(Arc::clone(&topology), geometries).unwrap();
    set.set_order(Some(Order::Step(vec![0, 10, 20, 30])))
        .unwrap();
    let mut columns = set.columns_mut();
    columns
        .insert(
            key("potential_energy"),
            PropertyColumn::Real {
                unit: KILOJOULE_PER_MOLE,
                values: vec![Some(-10.0), Some(-11.0), Some(-12.0), Some(-13.0)],
            },
        )
        .unwrap();
    columns
        .insert(
            key("bias"),
            PropertyColumn::Real {
                unit: KILOCALORIE_PER_MOLE,
                values: vec![None, None, Some(1.0), Some(2.0)],
            },
        )
        .unwrap();
    columns
        .insert(
            key("stage"),
            PropertyColumn::String(vec![Some("md".into()); 4]),
        )
        .unwrap();
    set
}

fn sampling() -> Sampling {
    let state = ThermodynamicState::canonical(Quantity::new(300.0, KELVIN)).unwrap();
    Sampling::new(state)
        .with_from_step(20)
        .with_potential_energy(key("potential_energy"))
        .unwrap()
        .with_bias([key("bias")])
        .unwrap()
}

fn sampled_set() -> StructureSet {
    let mut set = biased_set();
    set.set_sampling(Some(sampling())).unwrap();
    set
}

#[test]
fn a_record_names_each_column_once() {
    let state = ThermodynamicState::canonical(Quantity::new(300.0, KELVIN)).unwrap();
    let error = Sampling::new(state)
        .with_bias([key("a"), key("a")])
        .unwrap_err();
    assert_eq!(
        error,
        StructureSetError::RepeatedSamplingColumn { key: key("a") }
    );
    let error = Sampling::new(state)
        .with_bias([key("a")])
        .unwrap()
        .with_potential_energy(key("a"))
        .unwrap_err();
    assert_eq!(
        error.to_string(),
        "The sampling record names the column \"a\" twice."
    );
    let sampling = sampling();
    assert_eq!(sampling.from_step(), 20);
    assert_eq!(sampling.method(), None);
    assert_eq!(
        sampling.columns().collect::<Vec<_>>(),
        [&key("potential_energy"), &key("bias")]
    );
}

#[test]
fn a_record_names_real_energy_columns_or_absent_ones() {
    let mut set = biased_set();
    let state = ThermodynamicState::canonical(Quantity::new(300.0, KELVIN)).unwrap();
    let error = set
        .set_sampling(Some(
            Sampling::new(state).with_bias([key("stage")]).unwrap(),
        ))
        .unwrap_err();
    assert_eq!(
        error,
        StructureSetError::SamplingColumn { key: key("stage") }
    );
    assert_eq!(error.to_string(), "The column \"stage\" holds no energies.");
    assert!(set.sampling().is_none());
    // A column without values is absent and reads as missing energies.
    let absent = Sampling::new(state).with_bias([key("later")]).unwrap();
    set.set_sampling(Some(absent)).unwrap();
    let values = set.column_values(&ColumnRef::Key(key("later"))).unwrap();
    assert_eq!(values.value(), &vec![None; 4]);
    assert_eq!(values.unit(), KILOJOULE_PER_MOLE);
    assert!(matches!(
        set.column_values(&ColumnRef::Key(key("nothing"))),
        Err(StructureSetError::MissingColumn { .. })
    ));
    set.set_sampling(None).unwrap();
    assert!(set.sampling().is_none());
}

#[test]
fn edits_keep_the_record_columns_real_energies() {
    let mut set = sampled_set();
    let before = set.columns().clone();
    let mut columns = set.columns_mut();
    assert_eq!(
        columns.remove(&key("bias")).unwrap_err(),
        StructureSetError::ColumnInUse { key: key("bias") }
    );
    assert!(matches!(
        columns.clear(),
        Err(StructureSetError::ColumnInUse { .. })
    ));
    let error = columns
        .insert(key("bias"), PropertyColumn::Int(vec![Some(1); 4]))
        .unwrap_err();
    assert_eq!(
        error,
        StructureSetError::SamplingColumn { key: key("bias") }
    );
    assert!(columns
        .set_value(
            key("potential_energy"),
            GeometryIndex::new(0),
            Some(PropertyValue::Real {
                value: 1.0,
                unit: NANOMETER
            }),
        )
        .is_err());
    assert!(columns
        .set_values(
            key("bias"),
            [(GeometryIndex::new(0), Some(PropertyValue::Int(3)))],
        )
        .is_err());
    // An absent named column takes energies only.
    let absent = sampling().with_bias([key("bias"), key("later")]).unwrap();
    set.set_sampling(Some(absent)).unwrap();
    let mut columns = set.columns_mut();
    assert!(columns
        .set_value(
            key("later"),
            GeometryIndex::new(1),
            Some(PropertyValue::Bool(true))
        )
        .is_err());
    columns
        .set_value(
            key("later"),
            GeometryIndex::new(1),
            Some(PropertyValue::Real {
                value: 4.0,
                unit: KILOJOULE_PER_MOLE,
            }),
        )
        .unwrap();
    // Clearing every cell is allowed; the column is then absent.
    columns
        .set_values(key("later"), [(GeometryIndex::new(1), None)])
        .unwrap();
    assert_eq!(set.columns(), &before);
    // Other columns are edited freely.
    assert!(set.columns_mut().remove(&key("stage")).unwrap().is_some());
    set.set_sampling(None).unwrap();
    set.columns_mut().clear().unwrap();
    assert!(!set.columns().has_data());
}

#[test]
fn operations_that_keep_columns_keep_the_record() {
    let set = sampled_set();
    let expected = Some(&sampling());
    assert_eq!(set.select([3, 1]).unwrap().sampling(), expected);
    assert_eq!(
        StructureSet::from_reader(&set, [2]).unwrap().sampling(),
        expected
    );
    assert_eq!(set.clone().sampling(), expected);
    let oxygen =
        AtomSelection::from_indices(set.shared_topology(), [TopologyAtomIndex::new(0)]).unwrap();
    assert_eq!(set.subset(&oxygen).unwrap().sampling(), expected);
    let processed = Processed::new(&set, ProcessingSteps::default()).unwrap();
    assert_eq!(processed.sampling(), expected);
    assert_eq!(
        StructureSet::from_reader(&processed, 0..4)
            .unwrap()
            .sampling(),
        expected
    );
    let description = StructureSetDescription::of(&set);
    assert_eq!(description.sampling(), expected);
    assert_eq!(
        description.extended_to(6, None).unwrap().sampling(),
        expected
    );

    // Selecting only frames before the sampling leaves the bias absent; it
    // reads as missing.
    let early = set.select([0, 1]).unwrap();
    assert!(early.columns().get(&key("bias")).is_none());
    assert_eq!(
        early
            .column_values(&ColumnRef::Key(key("bias")))
            .unwrap()
            .value(),
        &[None, None]
    );
    assert_eq!(early.sampling(), expected);
}

#[test]
fn a_description_checks_the_record_against_its_columns() {
    let set = biased_set();
    let description = StructureSetDescription::of(&set);
    let state = ThermodynamicState::canonical(Quantity::new(300.0, KELVIN)).unwrap();
    let stage = Sampling::new(state).with_bias([key("stage")]).unwrap();
    assert!(description
        .clone()
        .with_sampling(Some(stage.clone()))
        .is_err());
    let described = StructureSetDescription::new(water(), 4)
        .with_sampling(Some(stage))
        .unwrap();
    assert!(described.with_columns(set.columns().clone()).is_err());
    let described = description.with_sampling(Some(sampling())).unwrap();
    assert_eq!(described.sampling(), Some(&sampling()));
}
