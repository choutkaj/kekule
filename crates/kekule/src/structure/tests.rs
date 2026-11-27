use super::*;
use crate::diagnostics::{Countable, IndexOutOfRange, LengthMismatch, TopologyMismatch};
use crate::math::{PeriodicCell, Point3, Vector3};
use crate::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use crate::properties::{PropertyKey, PropertyValue};
use crate::topology::{Topology, TopologyAtomIndex, TopologyBondIndex, TopologyBuilder};
use crate::units::{
    Quantity, ANGSTROM, KELVIN, NANOMETER, PICOSECOND, SQUARE_ANGSTROM, SQUARE_NANOMETER,
};
use std::sync::Arc;

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

const ATOM0: TopologyAtomIndex = TopologyAtomIndex::new(0);
const BOND0: TopologyBondIndex = TopologyBondIndex::new(0);

fn structure_fixture() -> (
    Structure,
    crate::topology::InstanceAtomId,
    crate::topology::InstanceBondId,
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
    let atom = crate::topology::InstanceAtomId::new(instance, carbon);
    let positions = Positions::new(Quantity::new(
        [Point3::origin(), Point3::new(1.0, 0.0, 0.0)],
        NANOMETER,
    ))
    .unwrap();
    (
        Structure::new(topology, positions).unwrap(),
        atom,
        crate::topology::InstanceBondId::new(instance, bond),
    )
}

fn single_atom_topology() -> crate::topology::Topology {
    let mut editor = MoleculeEditor::new();
    editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    crate::topology::Topology::from_molecule(editor.finish().unwrap()).unwrap()
}

fn single_position(x: f64) -> Positions {
    Positions::new(Quantity::new([Point3::new(x, 0.0, 0.0)], NANOMETER)).unwrap()
}

fn cell() -> PeriodicCell {
    PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(10.0, 10.0, 10.0), ANGSTROM),
        [true; 3],
    )
    .unwrap()
}

#[test]
fn canonical_structure_constructors_accept_owned_and_shared_topology() {
    let owned = Structure::new(single_atom_topology(), single_position(1.0)).unwrap();
    let shared = owned.shared_topology();
    let shared_structure = Structure::new(Arc::clone(&shared), single_position(2.0)).unwrap();
    assert!(Arc::ptr_eq(&shared_structure.shared_topology(), &shared));

    let complete = Structure::new(Arc::clone(&shared), single_position(3.0))
        .unwrap()
        .with_cell(cell());
    assert!(Arc::ptr_eq(&complete.shared_topology(), &shared));
    assert_eq!(complete.cell(), Some(&cell()));
    assert_eq!(complete.positions(), &single_position(3.0));
    assert_eq!(complete.occupancies(), None);
    assert!(complete.properties().is_empty());
}

#[test]
fn structure_builder_rejects_instance_position_mismatch_before_mutation() {
    let topology = single_atom_topology();
    let molecule = topology.instances().next().unwrap();
    let mut builder = StructureBuilder::new();
    let definition = builder
        .topology_mut()
        .add_molecule_definition(molecule.molecule().clone())
        .unwrap();
    let before = builder.clone();
    assert_eq!(
        builder.add_instance(definition, &Positions::zeros(0)),
        Err(StructureBuildError::InstancePositionCountMismatch(
            LengthMismatch::new(Countable::Position, 1, 0)
        ))
    );
    assert_eq!(builder, before);
    builder
        .add_instance(definition, &single_position(1.0))
        .unwrap();
    assert_eq!(builder.build().unwrap().atom_count(), 1);
}

#[test]
fn structure_builder_rejects_excess_replacement_positions_before_mutation() {
    let topology = single_atom_topology();
    let molecule = topology.instances().next().unwrap();
    let mut builder = StructureBuilder::new();
    builder
        .add_molecule(molecule.molecule().clone(), &single_position(1.0))
        .unwrap();
    let before = builder.clone();
    assert_eq!(
        builder
            .geometry_mut()
            .set_positions(Quantity::new(vec![Point3::origin(); 2], NANOMETER)),
        Err(GeometryError::Array(Box::new(
            DenseArrayError::CountMismatch(LengthMismatch::new(Countable::Row, 1, 2))
        )))
    );
    assert_eq!(builder, before);
    assert_eq!(builder.build().unwrap().positions(), &single_position(1.0));
}

#[test]
fn structures_reject_rows_of_other_dimensions() {
    let (structure, _, _) = structure_fixture();
    assert_eq!(
        Structure::new(structure.shared_topology(), single_position(1.0)).unwrap_err(),
        StructureError::Geometry(GeometryError::Array(Box::new(
            DenseArrayError::CountMismatch(LengthMismatch::new(Countable::Row, 2, 1))
        )))
    );
    // A geometry bound to a topology with one bond cannot bind to one with no
    // bonds: populated bond rows never resize silently.
    let two_atoms = Arc::new(
        crate::topology::Topology::from_molecules([
            single_atom_topology()
                .instances()
                .next()
                .unwrap()
                .molecule()
                .clone(),
            single_atom_topology()
                .instances()
                .next()
                .unwrap()
                .molecule()
                .clone(),
        ])
        .unwrap(),
    );
    let bound = StructureSet::from_structures([structure])
        .unwrap()
        .into_geometries()
        .pop()
        .unwrap();
    assert_eq!(bound.properties().bonds().len(), 1);
    assert_eq!(
        StructureSet::from_geometries(two_atoms, [bound]).unwrap_err(),
        StructureSetError::Geometry(GeometryError::BondCountMismatch(LengthMismatch::new(
            Countable::Row,
            0,
            1
        )))
    );
}

#[test]
fn geometry_position_replacement_checks_the_count_transactionally() {
    let (structure, _, _) = structure_fixture();
    let mut set = StructureSet::from_structures([structure.clone(), structure.clone()]).unwrap();
    assert_eq!(
        set.replace_positions([Positions::zeros(2)]),
        Err(StructureSetError::GeometryCountMismatch(
            LengthMismatch::new(Countable::Geometry, 2, 1)
        ))
    );
    assert_eq!(
        set.replace_positions([Positions::zeros(2), Positions::zeros(3)]),
        Err(StructureSetError::Geometry(GeometryError::Array(Box::new(
            DenseArrayError::CountMismatch(LengthMismatch::new(Countable::Row, 2, 3))
        ))))
    );
    assert!(set
        .iter()
        .all(|member| member.positions() == structure.positions()));
    set.replace_positions([Positions::zeros(2), Positions::zeros(2)])
        .unwrap();
    assert!(set
        .iter()
        .all(|member| member.positions() == &Positions::zeros(2)));
}

#[test]
fn canonical_structure_set_constructors_accept_owned_and_shared_topology() {
    let owned_new = StructureSet::new(single_atom_topology());
    let shared = owned_new.shared_topology();
    let shared_new = StructureSet::new(Arc::clone(&shared));
    assert!(Arc::ptr_eq(&shared_new.shared_topology(), &shared));

    let owned_members = StructureSet::from_geometries(
        single_atom_topology(),
        [Geometry::new(single_position(1.0))],
    )
    .unwrap();
    let shared_members = StructureSet::from_geometries(
        owned_members.shared_topology(),
        [Geometry::new(single_position(2.0))],
    )
    .unwrap();
    assert!(owned_members
        .topology()
        .same_layout(shared_members.topology()));
}

#[test]
fn sets_from_structures_accept_equal_layouts_and_refuse_others() {
    let structure = |smiles: &str| {
        let topology = crate::smiles::to_topology(smiles).unwrap();
        let atoms = topology.atom_count();
        Structure::new(topology, Positions::zeros(atoms)).unwrap()
    };
    // Published separately with equal contents: one layout.
    let set = StructureSet::from_structures([structure("CO"), structure("CO")]).unwrap();
    assert_eq!(set.len(), 2);
    // The same atoms and bonds in another dense order are another layout.
    for structures in [
        [structure("CO"), structure("OC")],
        [structure("CO"), structure("CC")],
    ] {
        assert_eq!(
            StructureSet::from_structures(structures).unwrap_err(),
            StructureSetError::TopologyMismatch(TopologyMismatch::Structure)
        );
    }
}

#[test]
fn structure_view_materialization_clones_geometry_and_shares_topology() {
    let mut source = Structure::new(single_atom_topology(), single_position(1.0)).unwrap();
    let topology = source.shared_topology();
    let score = key("score");
    let atom = topology.atom_ids()[0];
    source
        .geometry_mut()
        .properties_mut()
        .atoms_mut()
        .set_value(score.clone(), ATOM0, Some(PropertyValue::Int(4)))
        .unwrap();
    let owned = source.as_structure_view().to_structure();

    assert!(Arc::ptr_eq(&owned.shared_topology(), &topology));
    assert_eq!(owned, source);
    source
        .geometry_mut()
        .properties_mut()
        .atoms_mut()
        .set_value(score.clone(), ATOM0, Some(PropertyValue::Int(8)))
        .unwrap();
    source
        .set_position(atom, Quantity::new(Point3::new(9.0, 0.0, 0.0), NANOMETER))
        .unwrap();
    assert_eq!(
        owned.properties().atoms().value(&score, ATOM0).unwrap(),
        Some(PropertyValue::Int(4))
    );
    assert_eq!(owned.position(atom).unwrap().value().x, 1.0);
}

#[test]
fn occupancies_and_b_factors_are_typed_and_independent_of_generic_properties() {
    let (mut structure, _, _) = structure_fixture();
    structure
        .geometry_mut()
        .set_occupancy(ATOM0, Some(0.75))
        .unwrap();
    structure
        .geometry_mut()
        .set_b_factor(ATOM0, Some(Quantity::new(12.5, SQUARE_ANGSTROM)))
        .unwrap();
    // No reserved names: a generic property may use the same key.
    structure
        .geometry_mut()
        .properties_mut()
        .atoms_mut()
        .set_value(key("occupancy"), ATOM0, Some(PropertyValue::Int(1)))
        .unwrap();
    assert_eq!(structure.occupancy(ATOM0), Some(0.75));
    assert_eq!(structure.occupancies(), Some([Some(0.75), None].as_slice()));
    let b_factor = structure.b_factor(ATOM0).unwrap();
    assert_eq!(b_factor.unit(), SQUARE_NANOMETER);
    assert!((b_factor.into_value() - 0.125).abs() < 1.0e-12);
    assert_eq!(
        structure
            .properties()
            .atoms()
            .value(&key("occupancy"), ATOM0)
            .unwrap(),
        Some(PropertyValue::Int(1))
    );

    let before = structure.clone();
    assert_eq!(
        structure
            .geometry_mut()
            .set_occupancy(ATOM0, Some(f64::NAN)),
        Err(GeometryError::Array(Box::new(DenseArrayError::NonFinite {
            index: 0
        })))
    );
    assert!(matches!(
        structure
            .geometry_mut()
            .set_b_factor(ATOM0, Some(Quantity::new(1.0, KELVIN))),
        Err(GeometryError::Unit(_))
    ));
    assert_eq!(
        structure
            .geometry_mut()
            .set_occupancy(TopologyAtomIndex::new(2), Some(1.0)),
        Err(GeometryError::Array(Box::new(
            DenseArrayError::IndexOutOfRange(IndexOutOfRange::new(Countable::Row, 2, 2))
        )))
    );
    assert_eq!(
        structure
            .geometry_mut()
            .set_occupancies(Some(vec![Some(1.0)])),
        Err(GeometryError::Array(Box::new(
            DenseArrayError::CountMismatch(LengthMismatch::new(Countable::Row, 2, 1))
        )))
    );
    assert_eq!(structure, before);

    // Clearing every value compares equal to never recording any.
    let mut cleared = structure.clone();
    cleared.geometry_mut().set_occupancy(ATOM0, None).unwrap();
    cleared.geometry_mut().set_occupancies(None).unwrap();
    assert_eq!(cleared.occupancies(), None);
}

#[test]
fn structure_subset_projects_entity_state_and_drops_owner_properties() {
    let (mut structure, atom, bond) = structure_fixture();
    structure
        .geometry_mut()
        .set_occupancy(ATOM0, Some(0.5))
        .unwrap();
    let mut geometry = structure.geometry_mut();
    let mut properties = geometry.properties_mut();
    properties
        .owner_mut()
        .insert(key("energy"), PropertyValue::Int(3))
        .unwrap();
    properties
        .atoms_mut()
        .set_value(key("selected"), ATOM0, Some(PropertyValue::Bool(true)))
        .unwrap();
    properties
        .bonds_mut()
        .set_value(
            key("bond_selected"),
            BOND0,
            Some(PropertyValue::String("yes".into())),
        )
        .unwrap();
    let selection =
        crate::topology::AtomSelection::from_atoms(structure.shared_topology(), [atom]).unwrap();
    let subset = structure.subset(&selection).unwrap();
    assert!(subset.properties().owner().is_empty());
    assert_eq!(subset.occupancies(), Some([Some(0.5)].as_slice()));
    assert_eq!(
        subset
            .properties()
            .atoms()
            .value(&key("selected"), ATOM0)
            .unwrap(),
        Some(PropertyValue::Bool(true))
    );
    assert!(subset.properties().bonds().is_empty());
    let all_selection = crate::topology::AtomSelection::all(structure.shared_topology());
    let fully_retained = structure.subset(&all_selection).unwrap();
    let bond_index = fully_retained.topology().bond_ids()[0];
    assert_eq!(bond_index.bond(), bond.bond());
    assert_eq!(
        fully_retained
            .properties()
            .bonds()
            .value(&key("bond_selected"), BOND0)
            .unwrap(),
        Some(PropertyValue::String("yes".into()))
    );
}

/// A reader of three geometries on a step order that decodes a fresh
/// geometry on every read, even past its length, and counts its reads.
#[derive(Debug)]
struct DecodingReader {
    set: StructureSet,
    reads: std::sync::atomic::AtomicUsize,
}

impl DecodingReader {
    fn new() -> Self {
        let (structure, _, _) = structure_fixture();
        let mut set =
            StructureSet::from_geometries(structure.shared_topology(), (0..3).map(Self::decode))
                .unwrap();
        set.set_order(Some(Order::Step(vec![10, 11, 12]))).unwrap();
        Self {
            set,
            reads: std::sync::atomic::AtomicUsize::new(0),
        }
    }

    fn decode(index: usize) -> Geometry {
        let mut geometry = Geometry::new(
            Positions::new(Quantity::new(
                [Point3::new(index as f64, 0.0, 0.0), Point3::origin()],
                NANOMETER,
            ))
            .unwrap(),
        );
        geometry.set_velocities(Some(Velocities::zeros(2))).unwrap();
        geometry
            .set_occupancies(Some(vec![Some(0.5), None]))
            .unwrap();
        geometry
    }

    fn reads(&self) -> usize {
        self.reads.load(std::sync::atomic::Ordering::Relaxed)
    }
}

impl StructureSetRead for DecodingReader {
    fn topology(&self) -> &Arc<Topology> {
        StructureSetRead::topology(&self.set)
    }

    fn len(&self) -> usize {
        self.set.len()
    }

    fn order(&self) -> Option<&Order> {
        self.set.order()
    }

    fn peek(&self, _: usize) -> Option<GeometryView<'_, dyn StructureSetRead>> {
        None
    }

    fn read(&self, index: usize) -> Result<Geometry, StructureSetReadError> {
        self.reads
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let mut geometry = Self::decode(index);
        geometry.bind(self.set.topology()).unwrap();
        Ok(geometry)
    }
}

#[test]
fn sets_read_from_a_reader_are_faithful_copies_reading_each_geometry_once() {
    let reader = DecodingReader::new();

    // Each picked geometry whole, with its order value; a repeated index is
    // read once and shared.
    let set = StructureSet::from_reader(&reader, [2, 0, 2, 2]).unwrap();
    assert_eq!(reader.reads(), 2);
    assert_eq!(set.order(), Some(&Order::Step(vec![12, 10, 12, 12])));
    let first = set.get(0).unwrap();
    for index in [2, 3] {
        assert!(set.get(index).unwrap().shares_storage_with(&first));
    }
    assert!(!set.get(1).unwrap().shares_storage_with(&first));
    assert_eq!(first.positions().values().value()[0].x, 2.0);
    assert_eq!(set.get(1).unwrap().positions().values().value()[0].x, 0.0);
    assert!(first.velocities().is_some());
    assert_eq!(first.occupancies(), Some([Some(0.5), None].as_slice()));

    // Ascending indices are read in turn.
    let set = StructureSet::from_reader(&reader, [0, 1, 2]).unwrap();
    assert_eq!(reader.reads(), 5);
    assert_eq!(set.order(), reader.order());

    // Every index is checked before anything is read, so a reader that
    // reads past its length fails cleanly rather than selecting order
    // values it does not have.
    assert_eq!(
        StructureSet::from_reader(&reader, [1, 4]).unwrap_err(),
        StructureSetReadError::IndexOutOfRange(IndexOutOfRange::new(Countable::Geometry, 4, 3))
    );
    assert_eq!(reader.reads(), 5);

    // Geometries out of ascending order keep their times.
    let mut timed = StructureSet::new(StructureSetRead::topology(&reader.set).clone());
    timed
        .set_order(Some(Order::Time(
            TimeAxis::new(Some(Quantity::new(Vec::new(), PICOSECOND)), None).unwrap(),
        )))
        .unwrap();
    for index in 0..2 {
        let place = Place::Time {
            time: Some(Quantity::new(index as f64, PICOSECOND)),
            step: None,
        };
        timed
            .push(DecodingReader::decode(index), Some(place))
            .unwrap();
    }
    let reversed = StructureSet::from_reader(&timed, [1, 0]).unwrap();
    assert_eq!(
        reversed.validate_monotonic_time(),
        Err(StructureSetError::NonMonotonicTime { geometry: 1 })
    );
}

#[test]
fn clearing_velocities_and_forces_keeps_geometries_without_them_shared() {
    let reader = DecodingReader::new();
    let mut set = StructureSet::from_reader(&reader, [0, 1]).unwrap();
    set.get_mut(1).unwrap().take_velocities();
    let bare = set.get(1).unwrap().payload().clone();
    set.get_mut(0)
        .unwrap()
        .set_forces(Some(Forces::zeros(2)))
        .unwrap();
    set.clear_velocities_and_forces();
    for geometry in &set {
        assert!(geometry.velocities().is_none() && geometry.forces().is_none());
    }
    assert!(set.get(1).unwrap().shares_storage_with(&bare));
    assert_eq!(
        set.get(0).unwrap().occupancies(),
        Some([Some(0.5), None].as_slice())
    );
}

#[test]
fn occupancies_and_b_factors_follow_structures_into_geometries() {
    let (mut structure, _, _) = structure_fixture();
    structure
        .geometry_mut()
        .set_occupancy(ATOM0, Some(0.5))
        .unwrap();
    structure
        .geometry_mut()
        .set_b_factor(ATOM0, Some(Quantity::new(20.0, SQUARE_ANGSTROM)))
        .unwrap();
    structure.geometry_mut().set_cell(Some(cell()));

    // Structures become geometries with their occupancies, B-factors, cell, and rows.
    let set = StructureSet::from_structures([structure.clone()]).unwrap();
    let member = set.get(0).unwrap();
    assert_eq!(member.occupancies(), Some([Some(0.5), None].as_slice()));
    assert_eq!(
        member.b_factors().unwrap().value(),
        structure.b_factors().unwrap().value()
    );
    assert_eq!(member.cell(), Some(&cell()));
    assert_eq!(member.as_structure_view().occupancy(ATOM0), Some(0.5));
    assert_eq!(
        member
            .as_structure_view()
            .occupancy(TopologyAtomIndex::new(2)),
        None
    );
    assert_eq!(member.to_structure(), structure);
}

#[test]
fn structure_views_report_which_structure_they_borrow() {
    let topology = Arc::new(single_atom_topology());
    let structure = Structure::new(Arc::clone(&topology), single_position(0.0)).unwrap();
    assert_eq!(
        structure.as_structure_view().kind(),
        StructureKind::Standalone
    );

    let mut geometry = Geometry::new(single_position(0.0));
    geometry.set_velocities(Some(Velocities::zeros(1))).unwrap();
    let set = StructureSet::from_geometries(
        Arc::clone(&topology),
        [Geometry::new(single_position(0.0)), geometry],
    )
    .unwrap();
    let view = set.get(1).unwrap().as_structure_view();
    let StructureKind::InSet(context) = view.kind() else {
        panic!("a geometry of a set is in the set");
    };
    assert_eq!(context.place(), None);
    assert_eq!(view.velocities(), Some(&Velocities::zeros(1)));
    assert_eq!(view.forces(), None);
    let geometry = set.get(0).unwrap();
    assert_eq!(geometry.as_structure_view().velocities(), None);
    // Materializing a geometry as a structure drops its set context and
    // keeps its velocities.
    assert_eq!(
        geometry.to_structure().as_structure_view().kind(),
        StructureKind::Standalone
    );
    let structure = set.get(1).unwrap().to_structure();
    assert_eq!(structure.velocities(), Some(&Velocities::zeros(1)));
    assert_eq!(
        structure.as_structure_view().velocities(),
        Some(&Velocities::zeros(1))
    );

    let mut trajectory = StructureSet::new(topology);
    trajectory
        .set_order(Some(Order::Time(
            TimeAxis::new(None, Some(Vec::new())).unwrap(),
        )))
        .unwrap();
    let place = Place::Time {
        time: None,
        step: Some(4),
    };
    trajectory
        .push(Geometry::new(single_position(0.0)), Some(place))
        .unwrap();
    let frame = trajectory.get(0).unwrap();
    let StructureKind::InSet(context) = frame.as_structure_view().kind() else {
        panic!("a frame is a geometry");
    };
    assert_eq!(context.place(), Some(place));
    assert_eq!(frame.as_structure_view().forces(), None);
}

#[test]
fn rigid_transforms_move_structures_and_geometries_with_their_cells() {
    use crate::math::{Matrix3, RigidTransform};
    let topology = Arc::new(single_atom_topology());
    let box_cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(1.0, 2.0, 3.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    let mut structure = Structure::new(Arc::clone(&topology), single_position(1.0)).unwrap();
    structure.geometry_mut().set_cell(Some(box_cell));
    // A quarter turn about z followed by a 2 nm shift along x.
    let turn = Matrix3::from_columns(
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    );
    let transform = RigidTransform::from_canonical(turn, Vector3::new(2.0, 0.0, 0.0)).unwrap();
    let mut members = StructureSet::from_structures([structure.clone()]).unwrap();
    let mut trajectory = StructureSet::from_geometries(
        Arc::clone(&topology),
        [Geometry::new(single_position(1.0)).with_cell(box_cell)],
    )
    .unwrap();
    trajectory
        .set_order(Some(Order::Time(TimeAxis::default())))
        .unwrap();

    structure
        .geometry_mut()
        .apply_rigid_transform(&transform)
        .unwrap();
    for set in [&mut members, &mut trajectory] {
        set.get_mut(0)
            .unwrap()
            .apply_rigid_transform(&transform)
            .unwrap();
    }

    assert_eq!(
        structure.positions().values().value(),
        &[Point3::new(2.0, 1.0, 0.0)]
    );
    assert_eq!(
        structure.cell().unwrap().vectors().into_value(),
        [
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-2.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 3.0),
        ]
    );
    for view in [
        members.get(0).unwrap().as_structure_view(),
        trajectory.get(0).unwrap().as_structure_view(),
    ] {
        assert_eq!(view.positions(), structure.positions());
        assert_eq!(view.cell(), structure.cell());
    }
}

mod construction;
mod structure_set;
