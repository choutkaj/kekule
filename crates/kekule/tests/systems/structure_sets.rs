//! Structure sets: what an operation needs, error texts, places kept on
//! replacement, empty selections, reading another set, and copy on write.

use std::sync::Arc;

use kekule::diagnostics::{
    CapacityExceeded, Countable, IndexOutOfRange, LengthMismatch, TopologyMismatch,
};
use kekule::math::{Matrix3, PeriodicCell, Point3, RigidTransform, Vector3};
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::structure::{
    DenseArrayError, Forces, Geometry, GeometryError, Order, Place, Positions, Requirement,
    Structure, StructureKind, StructureSet, StructureSetError, StructureSetReadError, Velocities,
};
use kekule::topology::{
    InstanceAtomId, InstanceBondId, Topology, TopologyAtomIndex, TopologyBondIndex, TopologyBuilder,
};
use kekule::units::{
    Quantity, ANGSTROM, CANONICAL_VELOCITY_UNIT, DEGREE, NANOMETER, PICOSECOND, SQUARE_ANGSTROM,
};

const ATOM0: TopologyAtomIndex = TopologyAtomIndex::new(0);
const BOND0: TopologyBondIndex = TopologyBondIndex::new(0);

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

/// Carbon monoxide at the origin and 1 nm along x, its carbon, and its bond.
fn carbon_monoxide() -> (Structure, InstanceAtomId, InstanceBondId) {
    let mut editor = MoleculeEditor::new();
    let carbon = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let oxygen = editor
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    let bond = editor.add_bond(carbon, oxygen, BondOrder::Single).unwrap();
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(editor.finish().unwrap()).unwrap();
    let topology = Arc::new(builder.build().unwrap());
    let positions = Positions::new(Quantity::new(
        [Point3::origin(), Point3::new(1.0, 0.0, 0.0)],
        NANOMETER,
    ))
    .unwrap();
    (
        Structure::new(topology, positions).unwrap(),
        InstanceAtomId::new(instance, carbon),
        InstanceBondId::new(instance, bond),
    )
}

fn single_atom_topology() -> Topology {
    let mut editor = MoleculeEditor::new();
    editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    Topology::from_molecule(editor.finish().unwrap()).unwrap()
}

fn single_position(x: f64) -> Positions {
    Positions::new(Quantity::new([Point3::new(x, 0.0, 0.0)], NANOMETER)).unwrap()
}

/// Geometries of a one-atom topology at these x positions, unordered.
fn set_at(xs: &[f64]) -> StructureSet {
    StructureSet::from_geometries(
        single_atom_topology(),
        xs.iter().map(|x| Geometry::new(single_position(*x))),
    )
    .unwrap()
}

fn x(set: &StructureSet, geometry: usize) -> f64 {
    set.get(geometry).unwrap().positions().values().value()[0].x
}

fn xs(set: &StructureSet) -> Vec<f64> {
    (0..set.len()).map(|geometry| x(set, geometry)).collect()
}

fn velocities(x: f64) -> Velocities {
    Velocities::new(Quantity::new(
        [Vector3::new(x, 0.0, 0.0)],
        CANONICAL_VELOCITY_UNIT,
    ))
    .unwrap()
}

fn cell() -> PeriodicCell {
    PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(10.0, 10.0, 10.0), ANGSTROM),
        [true; 3],
    )
    .unwrap()
}

fn shares_every_geometry(a: &StructureSet, b: &StructureSet) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b.iter())
            .all(|(a, b)| a.payload().shares_storage_with(b.payload()))
}

#[test]
fn requirements_read_as_they_are() {
    assert_eq!(Requirement::TimeOrder.to_string(), "needs a time order");
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
fn any_geometry_may_carry_site_fields_and_dynamics() {
    let mut set = set_at(&[0.0]);
    set.set_order(Some(Order::Step(vec![0]))).unwrap();
    let mut geometry = set.get_mut(0).unwrap();
    geometry.set_occupancies(Some(vec![Some(0.5)])).unwrap();
    geometry.set_velocities(Some(velocities(1.0))).unwrap();
    assert_eq!(
        geometry.set_forces(Some(Forces::zeros(2))),
        Err(GeometryError::Array(Box::new(
            DenseArrayError::CountMismatch(LengthMismatch::new(Countable::Row, 1, 2))
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
    let turn = Matrix3::from_columns(
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    );
    let transform = RigidTransform::new(turn, Quantity::new(Vector3::zero(), NANOMETER)).unwrap();
    let mut geometry = set.get_mut(0).unwrap();
    geometry.apply_rigid_transform(&transform).unwrap();
    assert_eq!(
        geometry.velocities().unwrap().values().value(),
        &[Vector3::new(0.0, 1.0, 0.0)]
    );
    assert!(set.get_mut(0).unwrap().take_velocities().is_some());
    assert_eq!(set.get(0).unwrap().velocities(), None);
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
    let options = kekule::sdf::SdfWriteOptions::default();
    // One carbon at x = 0 and 0.1 nm: flat, so 2D, with coordinates in Å
    // and no data fields.
    let set = set_at(&[0.0, 0.1]);
    let views = set.iter().collect::<Vec<_>>();
    let written = kekule::sdf::write(
        views.iter().map(kekule::sdf::SdfRecordSource::structure),
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
    let written = kekule::sdf::write(
        views.iter().map(kekule::sdf::SdfRecordSource::structure),
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
fn replacing_a_geometry_keeps_its_place() {
    let mut set = set_at(&[0.0, 1.0]);
    set.set_order(Some(Order::Step(vec![3, 4]))).unwrap();
    let previous = set.replace(1, Geometry::new(single_position(5.0))).unwrap();
    assert_eq!(previous, Geometry::new(single_position(1.0)));
    assert_eq!(xs(&set), [0.0, 5.0]);
    assert_eq!(set.get(1).unwrap().place(), Some(Place::Step(4)));

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
    assert!(shares_every_geometry(&set, &before));
    assert_eq!(set.order(), before.order());
}

#[test]
fn empty_selections() {
    let mut frames = set_at(&[0.0, 1.0]);
    frames.set_order(Some(Order::Step(vec![0, 1]))).unwrap();
    let none = frames.select([]).unwrap();
    assert!(none.is_empty());
    assert_eq!(none.order(), Some(&Order::Step(Vec::new())));
    frames.remove(1).unwrap();
    frames.remove(0).unwrap();
    assert!(frames.is_empty());
    assert_eq!(frames.order(), Some(&Order::Step(Vec::new())));
}

#[test]
fn every_error_reads_as_a_sentence() {
    use StructureSetError as E;
    let requirements = [
        (Requirement::TimeOrder, "needs a time order"),
        (Requirement::AnyOrder, "needs an order"),
        (Requirement::Times, "needs times"),
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
        (
            E::IndexOutOfRange(IndexOutOfRange::new(Countable::Geometry, 3, 1)),
            "Geometry 3 does not exist. There is 1 geometry.",
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
            E::InvalidTemperature,
            "The temperature is not finite and positive.",
        ),
        (
            E::InvalidPressure,
            "The pressure is not finite and positive.",
        ),
        (E::NonFiniteTime, "The time is not finite."),
        (
            E::CoverageMismatch {
                covered: 2,
                geometries: 1,
            },
            "The description covers 2 geometries, but the structure set has 1.",
        ),
        (
            E::Geometry(GeometryError::CapacityExceeded(CapacityExceeded::new(
                Countable::Atom,
            ))),
            "A geometry is invalid.",
        ),
        (
            E::Subset(kekule::topology::TopologySubsetError::EmptySelection),
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
fn geometry_replacement_rejects_incompatible_geometries_without_changing_state() {
    let (structure, _, _) = carbon_monoxide();
    let topology = structure.shared_topology();
    let mut set = StructureSet::from_geometries(
        Arc::clone(&topology),
        [Geometry::new(structure.positions().clone())],
    )
    .unwrap();
    let original = set.get(0).unwrap().payload().clone();
    // Bound to another one-bond topology; matching bond rows rebind freely.
    let other = carbon_monoxide().0;
    let mut bound_elsewhere = StructureSet::from_geometries(
        Topology::from_molecules([other
            .topology()
            .instances()
            .next()
            .unwrap()
            .molecule()
            .clone()])
        .unwrap(),
        [Geometry::new(Positions::zeros(2))],
    )
    .unwrap()
    .into_geometries()
    .remove(0);
    bound_elsewhere
        .properties_mut()
        .bonds_mut()
        .set_value(key("bond"), BOND0, Some(PropertyValue::Int(1)))
        .unwrap();

    let replacement = Geometry::new(Positions::zeros(1));
    let expected = StructureSetError::Geometry(GeometryError::Array(Box::new(
        DenseArrayError::CountMismatch(LengthMismatch::new(Countable::Row, 2, 1)),
    )));
    assert_eq!(set.replace(0, replacement.clone()), Err(expected.clone()));
    assert_eq!(set.push(replacement, None), Err(expected));
    assert_eq!(set.len(), 1);
    let member = set.get(0).unwrap();
    assert_eq!(member.payload(), &original);
    assert_eq!(
        member.as_structure_view().positions().len(),
        topology.atom_count()
    );
    // Same bond count: a payload bound elsewhere with matching rows is accepted.
    set.push(bound_elsewhere, None).unwrap();
    assert_eq!(set.len(), 2);
    set.remove(1).unwrap();

    for index in [1, usize::MAX] {
        assert_eq!(
            set.replace(index, original.clone()),
            Err(StructureSetError::IndexOutOfRange(IndexOutOfRange::new(
                Countable::Geometry,
                index,
                1
            )))
        );
        assert!(set.get_mut(index).is_none());
    }
    let mut empty = StructureSet::new(topology);
    assert_eq!(
        empty.replace(0, original.clone()),
        Err(StructureSetError::IndexOutOfRange(IndexOutOfRange::new(
            Countable::Geometry,
            0,
            0
        )))
    );
    assert!(empty.is_empty());
    assert_eq!(set.replace(0, original.clone()).unwrap(), original);
}

#[test]
fn geometry_replacement_preserves_order_and_publishes_the_whole_geometry() {
    let (structure, _, _) = carbon_monoxide();
    let topology = structure.shared_topology();
    let original = Geometry::new(structure.positions().clone());
    let mut set =
        StructureSet::from_geometries(Arc::clone(&topology), [original.clone(), original]).unwrap();
    let score = key("score");
    set.properties_mut()
        .insert(score.clone(), PropertyValue::Int(9))
        .unwrap();
    let set_properties = set.properties().clone();
    let original = set.get(0).unwrap().payload().clone();
    let mut replacement = Geometry::new(Positions::zeros(2));
    replacement.set_cell(Some(cell()));
    replacement
        .properties_mut()
        .owner_mut()
        .insert(score.clone(), PropertyValue::Int(1))
        .unwrap();
    replacement
        .properties_mut()
        .atoms_mut()
        .set_value(score, ATOM0, Some(PropertyValue::Int(2)))
        .unwrap();

    assert_eq!(set.replace(1, replacement.clone()).unwrap(), original);
    assert_eq!(set.len(), 2);
    assert_eq!(set.properties(), &set_properties);
    assert_eq!(set.get(0).unwrap().payload(), &original);
    let member = set.get(1).unwrap();
    assert!(Arc::ptr_eq(&member.shared_topology(), &topology));
    assert_eq!(member.positions(), replacement.positions());
    assert_eq!(member.cell(), replacement.cell());
    assert!(member.properties().owner().is_empty());
    assert_eq!(
        &member.row().unwrap().to_properties(),
        replacement.properties().owner()
    );
    assert_eq!(
        member.properties().atoms(),
        replacement.properties().atoms()
    );
    // Binding allocated the bond rows the detached replacement lacked.
    assert_eq!(member.properties().bonds().len(), 1);
    assert_eq!(member.to_structure().positions(), replacement.positions());
}

#[test]
fn sets_read_chosen_geometries_of_another_set_without_consuming_it() {
    let (structure, carbon, _) = carbon_monoxide();
    let frames = (0..5).map(|index| {
        let mut geometry = Geometry::new(
            Positions::new(Quantity::new(
                [Point3::new(index as f64, 0.0, 0.0), Point3::origin()],
                NANOMETER,
            ))
            .unwrap(),
        );
        geometry.set_velocities(Some(Velocities::zeros(2))).unwrap();
        geometry
    });
    let mut trajectory =
        StructureSet::from_geometries(structure.shared_topology(), frames).unwrap();
    trajectory
        .properties_mut()
        .insert(key("run"), PropertyValue::Int(7))
        .unwrap();
    // The copy keeps the set properties.
    assert_eq!(
        StructureSet::from_reader(&trajectory, [2])
            .unwrap()
            .properties(),
        trajectory.properties()
    );

    // Leave out equilibration and take every second geometry.
    let sampled = StructureSet::from_reader(&trajectory, (1..5).step_by(2)).unwrap();
    assert_eq!(trajectory.len(), 5);
    let geometries = sampled
        .iter()
        .map(|geometry| geometry.positions().values().value()[0].x)
        .collect::<Vec<_>>();
    assert_eq!(geometries, [1.0, 3.0]);
    assert!(sampled
        .get(0)
        .unwrap()
        .shares_storage_with(&trajectory.get(1).unwrap()));

    let atoms =
        kekule::topology::AtomSelection::from_atoms(structure.shared_topology(), [carbon]).unwrap();
    let subset = StructureSet::from_reader(&trajectory, [4])
        .unwrap()
        .subset(&atoms)
        .unwrap();
    assert_eq!(subset.topology().atom_count(), 1);
    assert_eq!(subset.get(0).unwrap().positions().len(), 1);

    assert_eq!(
        StructureSet::from_reader(&trajectory, [5]).unwrap_err(),
        StructureSetReadError::IndexOutOfRange(IndexOutOfRange::new(Countable::Geometry, 5, 5))
    );
}

#[test]
fn editing_one_geometry_of_a_shared_set_copies_that_geometry_alone() {
    let original = StructureSet::from_geometries(
        single_atom_topology(),
        (0..3).map(|_| Geometry::new(single_position(0.0))),
    )
    .unwrap();
    let mut copy = original.clone();
    let moved = Quantity::new(Point3::new(5.0, 0.0, 0.0), NANOMETER);
    copy.get_mut(1)
        .unwrap()
        .set_position_at(ATOM0, moved)
        .unwrap();
    copy.get_mut(2)
        .unwrap()
        .set_b_factor(ATOM0, Some(Quantity::new(1.0, SQUARE_ANGSTROM)))
        .unwrap();

    assert_eq!(x(&copy, 1), 5.0);
    assert_eq!(x(&original, 1), 0.0);
    assert!(original.get(2).unwrap().b_factors().is_none());
    assert!(copy.get(2).unwrap().b_factors().is_some());
    for (geometry, shared) in [(0, true), (1, false), (2, false)] {
        assert_eq!(
            copy.get(geometry)
                .unwrap()
                .payload()
                .shares_storage_with(original.get(geometry).unwrap().payload()),
            shared,
            "geometry {geometry}"
        );
    }

    // Transforms and position replacement leave the original unchanged.
    let mut transformed = original.clone();
    let shift = RigidTransform::new(
        Matrix3::identity(),
        Quantity::new(Vector3::new(2.0, 0.0, 0.0), NANOMETER),
    )
    .unwrap();
    transformed
        .get_mut(0)
        .unwrap()
        .apply_rigid_transform(&shift)
        .unwrap();
    assert_eq!(x(&transformed, 0), 2.0);
    let mut replaced = original.clone();
    replaced
        .replace_positions((0..3).map(|index| single_position(f64::from(index) + 7.0)))
        .unwrap();
    assert_eq!(x(&replaced, 2), 9.0);
    assert_eq!(
        (0..3)
            .map(|geometry| x(&original, geometry))
            .collect::<Vec<_>>(),
        [0.0; 3]
    );

    // An edit that fails after copying changes neither set.
    let mut failed = original.clone();
    assert!(failed
        .get_mut(0)
        .unwrap()
        .set_positions(Quantity::new([Point3::origin(); 2], NANOMETER))
        .is_err());
    assert_eq!(
        failed.get(0).unwrap().payload(),
        original.get(0).unwrap().payload()
    );
}
