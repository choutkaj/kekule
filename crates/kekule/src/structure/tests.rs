use super::*;
use crate::math::{PeriodicCell, Point3, Vector3};
use crate::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use crate::properties::{PropertyColumn, PropertyError, PropertyKey, PropertyValue};
use crate::topology::{Topology, TopologyAtomIndex, TopologyBondIndex, TopologyBuilder};
use crate::units::{
    Quantity, UnitError, ANGSTROM, KELVIN, KILOCALORIE_PER_MOLE, KILOJOULE_PER_MOLE, NANOMETER,
    PICOSECOND, SQUARE_ANGSTROM, SQUARE_NANOMETER,
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

/// Distinct, equally weighted geometries, each its own state.
fn conformers(
    topology: impl Into<Arc<crate::topology::Topology>>,
    geometries: impl IntoIterator<Item = Geometry>,
) -> Result<StructureSet, StructureSetError> {
    let mut set = StructureSet::from_geometries(topology, geometries)?;
    let partition = Partition::state_per_sample(set.len());
    set.set_distribution(Some(Distribution::new(Members::All, partition)?))?;
    Ok(set)
}

/// Structures drawn from one distribution: one state of equal weights.
fn bundle(
    structures: impl IntoIterator<Item = Structure>,
) -> Result<StructureSet, StructureSetError> {
    let mut set = StructureSet::from_structures(structures)?;
    let partition = Partition::one_state(set.len());
    set.set_distribution(Some(Distribution::new(Members::All, partition)?))?;
    Ok(set)
}

fn cell() -> PeriodicCell {
    PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(10.0, 10.0, 10.0), ANGSTROM),
        [true; 3],
    )
    .unwrap()
}

#[test]
fn canonical_model_constructors_accept_owned_and_shared_topology() {
    let owned = Structure::new(single_atom_topology(), single_position(1.0)).unwrap();
    let shared = owned.shared_topology();
    let shared_model = Structure::new(Arc::clone(&shared), single_position(2.0)).unwrap();
    assert!(Arc::ptr_eq(&shared_model.shared_topology(), &shared));

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
fn model_builder_rejects_instance_position_mismatch_before_mutation() {
    let topology = single_atom_topology();
    let molecule = topology.instances().next().unwrap();
    let mut builder = StructureBuilder::new();
    let definition = builder
        .add_molecule_definition(molecule.molecule().clone())
        .unwrap();
    let before = builder.clone();
    assert_eq!(
        builder.add_instance(definition, &Positions::zeros(0)),
        Err(StructureBuildError::InstancePositionCountMismatch {
            expected: 1,
            actual: 0
        })
    );
    assert_eq!(builder, before);
    builder
        .add_instance(definition, &single_position(1.0))
        .unwrap();
    assert_eq!(builder.build().unwrap().atom_count(), 1);
}

#[test]
fn model_builder_rejects_excess_replacement_positions_before_mutation() {
    let topology = single_atom_topology();
    let molecule = topology.instances().next().unwrap();
    let mut builder = StructureBuilder::new();
    builder
        .add_molecule(molecule.molecule().clone(), &single_position(1.0))
        .unwrap();
    let before = builder.clone();
    assert_eq!(
        builder.set_positions(Quantity::new(vec![Point3::origin(); 2], NANOMETER)),
        Err(GeometryError::Array(DenseArrayError::CountMismatch {
            expected: 1,
            actual: 2
        }))
    );
    assert_eq!(builder, before);
    assert_eq!(builder.build().unwrap().positions(), &single_position(1.0));
}

#[test]
fn structures_reject_rows_of_other_dimensions() {
    let (structure, _, _) = structure_fixture();
    assert_eq!(
        Structure::new(structure.shared_topology(), single_position(1.0)).unwrap_err(),
        StructureError::Geometry(GeometryError::Array(DenseArrayError::CountMismatch {
            expected: 2,
            actual: 1
        }))
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
    let bound = bundle([structure])
        .unwrap()
        .into_geometries()
        .pop()
        .unwrap();
    assert_eq!(bound.properties().bonds().len(), 1);
    assert_eq!(
        conformers(two_atoms, [bound]).unwrap_err(),
        StructureSetError::Geometry(GeometryError::BondCountMismatch {
            expected: 0,
            actual: 1
        })
    );
}

#[test]
fn sample_position_replacement_checks_the_count_transactionally() {
    let (structure, _, _) = structure_fixture();
    let mut set = bundle([structure.clone(), structure.clone()]).unwrap();
    assert_eq!(
        set.replace_positions([Positions::zeros(2)]),
        Err(StructureSetError::GeometryCountMismatch {
            expected: 2,
            actual: 1
        })
    );
    assert_eq!(
        set.replace_positions([Positions::zeros(2), Positions::zeros(3)]),
        Err(StructureSetError::Geometry(GeometryError::Array(
            DenseArrayError::CountMismatch {
                expected: 2,
                actual: 3
            }
        )))
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

    let owned_members = conformers(
        single_atom_topology(),
        [Geometry::new(single_position(1.0))],
    )
    .unwrap();
    let shared_members = conformers(
        owned_members.shared_topology(),
        [Geometry::new(single_position(2.0))],
    )
    .unwrap();
    assert!(owned_members
        .topology()
        .same_layout(shared_members.topology()));
}

#[test]
fn sets_from_models_accept_equal_layouts_and_refuse_others() {
    let structure = |smiles: &str| {
        let topology = crate::smiles::to_topology(smiles).unwrap();
        let atoms = topology.atom_count();
        Structure::new(topology, Positions::zeros(atoms)).unwrap()
    };
    // Published separately with equal contents: one layout.
    let set = bundle([structure("CO"), structure("CO")]).unwrap();
    assert_eq!(set.len(), 2);
    // The same atoms and bonds in another dense order are another layout.
    for structures in [
        [structure("CO"), structure("OC")],
        [structure("CO"), structure("CC")],
    ] {
        assert_eq!(
            bundle(structures).unwrap_err(),
            StructureSetError::TopologyMismatch
        );
    }
}

#[test]
fn member_views_project_borrowed_and_owned_models_in_stable_order() {
    let topology = Arc::new(single_atom_topology());
    let score = key("score");
    let mut first = Geometry::new(single_position(1.0));
    first.set_cell(Some(cell()));
    let mut properties = first.properties_mut();
    properties
        .owner_mut()
        .insert(key("method"), PropertyValue::String("sampled".into()))
        .unwrap();
    properties
        .atoms_mut()
        .set_value(score.clone(), ATOM0, Some(PropertyValue::Int(7)))
        .unwrap();
    let second = Geometry::new(single_position(2.0));
    let mut set = conformers(Arc::clone(&topology), [first, second]).unwrap();
    set.distribution_mut()
        .unwrap()
        .set_weights(Some(vec![0.25, 0.75]))
        .unwrap();

    assert_eq!(
        set.iter()
            .map(|member| member.positions().values().value()[0].x)
            .collect::<Vec<_>>(),
        [1.0, 2.0]
    );
    assert_eq!(weights(&set), [0.25, 0.75]);

    let member = set.get(0).unwrap();
    assert!(Arc::ptr_eq(&member.shared_topology(), &topology));
    let borrowed = member.as_structure_view();
    assert!(Arc::ptr_eq(&borrowed.shared_topology(), &topology));
    assert_eq!(
        borrowed.positions().values().value().as_ptr(),
        member.positions().values().value().as_ptr()
    );
    let owned = member.to_structure();
    assert!(Arc::ptr_eq(&owned.shared_topology(), &topology));
    assert_eq!(owned.positions(), member.positions());
    assert_eq!(owned.cell(), member.cell());
    assert_eq!(owned.properties(), member.properties());

    set.get_mut(0)
        .unwrap()
        .properties_mut()
        .atoms_mut()
        .set_value(score.clone(), ATOM0, Some(PropertyValue::Int(9)))
        .unwrap();
    assert_eq!(
        owned.properties().atoms().value(&score, ATOM0).unwrap(),
        Some(PropertyValue::Int(7))
    );
}

#[test]
fn detached_samples_get_bond_rows_only_when_bound() {
    let (structure, _, _) = structure_fixture();
    let label = key("label");
    let mut member = Geometry::new(structure.positions().clone());
    assert_eq!(member.properties().atoms().len(), 2);
    assert!(member.properties().bonds().is_empty());
    // Detached payloads cannot carry bond annotations yet.
    assert!(matches!(
        member.properties_mut().bonds_mut().set_value(
            label.clone(),
            BOND0,
            Some(PropertyValue::Int(1))
        ),
        Err(PropertyError::InvalidIndex { len: 0, index: 0 })
    ));
    member
        .properties_mut()
        .owner_mut()
        .insert(label.clone(), PropertyValue::Int(7))
        .unwrap();
    let mut set = conformers(structure.shared_topology(), [member]).unwrap();
    let stored = set.get(0).unwrap();
    assert_eq!(stored.properties().atoms().len(), 2);
    assert_eq!(stored.properties().bonds().len(), 1);
    assert_eq!(
        stored.properties().owner().get(&label),
        Some(&PropertyValue::Int(7))
    );

    let mut editor = set.get_mut(0).unwrap();
    let before = editor.properties().clone();
    assert!(editor
        .properties_mut()
        .bonds_mut()
        .insert(label.clone(), PropertyColumn::Int(vec![Some(8); 2]))
        .is_err());
    assert_eq!(editor.properties(), &before);
    editor
        .properties_mut()
        .bonds_mut()
        .insert(label.clone(), PropertyColumn::Int(vec![Some(8)]))
        .unwrap();
    assert_eq!(
        editor.properties().bonds().value(&label, BOND0).unwrap(),
        Some(PropertyValue::Int(8))
    );
}

#[test]
fn sample_replacement_rejects_incompatible_samples_without_changing_state() {
    let (structure, _, _) = structure_fixture();
    let topology = structure.shared_topology();
    let mut set = conformers(
        Arc::clone(&topology),
        [Geometry::new(structure.positions().clone())],
    )
    .unwrap();
    set.distribution_mut().unwrap().set_weight(0, 0.75).unwrap();
    let original = set.get(0).unwrap().payload().clone();
    // Bound to another one-bond topology; matching bond rows rebind freely.
    let other = structure_fixture().0;
    let mut bound_elsewhere = conformers(
        crate::topology::Topology::from_molecules([other
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
    let expected =
        StructureSetError::Geometry(GeometryError::Array(DenseArrayError::CountMismatch {
            expected: 2,
            actual: 1,
        }));
    assert_eq!(set.replace(0, replacement.clone()), Err(expected.clone()));
    assert_eq!(
        set.push_state(replacement, None, State::new(), Some(1.0)),
        Err(expected)
    );
    assert_eq!(set.len(), 1);
    let member = set.get(0).unwrap();
    assert_eq!(member.payload(), &original);
    assert_eq!(member.weight(), Some(0.75));
    assert_eq!(
        member.as_structure_view().positions().len(),
        topology.atom_count()
    );
    // Same bond count: a payload bound elsewhere with matching rows is accepted.
    set.push_member(bound_elsewhere, None, 0, Some(1.0))
        .unwrap();
    assert_eq!(set.len(), 2);
    set.remove(1).unwrap();

    for index in [1, usize::MAX] {
        assert_eq!(
            set.replace(index, original.clone()),
            Err(StructureSetError::IndexOutOfRange { index, len: 1 })
        );
        assert!(set.get_mut(index).is_none());
    }
    let mut empty = StructureSet::new(topology);
    assert_eq!(
        empty.replace(0, original.clone()),
        Err(StructureSetError::IndexOutOfRange { index: 0, len: 0 })
    );
    assert!(empty.is_empty());
    assert_eq!(set.replace(0, original.clone()).unwrap(), original);
}

#[test]
fn sample_replacement_preserves_order_and_membership_and_publishes_the_whole_sample() {
    let (structure, _, _) = structure_fixture();
    let topology = structure.shared_topology();
    let original = Geometry::new(structure.positions().clone());
    let mut set = conformers(Arc::clone(&topology), [original.clone(), original]).unwrap();
    set.distribution_mut()
        .unwrap()
        .set_weights(Some(vec![1.0, 0.5]))
        .unwrap();
    let original = set.get(0).unwrap().payload().clone();
    let score = key("score");
    set.properties_mut()
        .insert(score.clone(), PropertyValue::Int(9))
        .unwrap();
    let set_properties = set.properties().clone();
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
    // The replacement takes the replaced geometry's weight and state.
    assert_eq!(member.weight(), Some(0.5));
    assert_eq!(member.state().map(|state| state.index()), Some(1));
    assert_eq!(
        member.properties().owner(),
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
fn sample_editor_preserves_dimensions_and_validated_values() {
    let (structure, _, _) = structure_fixture();
    let mut set = bundle([structure]).unwrap();
    let mut distribution = set.distribution_mut().unwrap();
    distribution.set_weight(0, 0.5).unwrap();
    // Weights are finite and strictly positive.
    for weight in [f64::NAN, f64::INFINITY, -1.0, 0.0] {
        assert_eq!(
            distribution.set_weight(0, weight),
            Err(StructureSetError::InvalidWeight { geometry: 0 })
        );
        assert_eq!(distribution.weight(0), Some(0.5));
    }
    {
        let mut member = set.get_mut(0).unwrap();
        let properties = member.properties().clone();
        let detached = Geometry::new(Positions::zeros(2)).properties().clone();
        assert_eq!(
            member.set_properties(detached),
            Err(GeometryError::BondCountMismatch {
                expected: 1,
                actual: 0
            })
        );
        let wrong_atoms = Geometry::new(Positions::zeros(1)).properties().clone();
        assert_eq!(
            member.set_properties(wrong_atoms),
            Err(GeometryError::Array(DenseArrayError::CountMismatch {
                expected: 2,
                actual: 1
            }))
        );
        assert_eq!(member.properties(), &properties);
        assert_eq!(
            member.set_positions(Quantity::new(vec![Point3::origin(); 3], NANOMETER)),
            Err(GeometryError::Array(DenseArrayError::CountMismatch {
                expected: 2,
                actual: 3
            }))
        );
        assert_eq!(member.positions().len(), 2);
    }
    let member = set.get(0).unwrap();
    assert_eq!(member.weight(), Some(0.5));
    assert_eq!(member.as_structure_view().positions().len(), 2);
    assert_eq!(member.to_structure().properties().bonds().len(), 1);
}

#[test]
fn model_view_materialization_clones_structure_and_shares_topology() {
    let mut source = Structure::new(single_atom_topology(), single_position(1.0)).unwrap();
    let topology = source.shared_topology();
    let score = key("score");
    let atom = topology.atom_ids()[0];
    source
        .properties_mut()
        .atoms_mut()
        .set_value(score.clone(), ATOM0, Some(PropertyValue::Int(4)))
        .unwrap();
    let owned = source.as_structure_view().to_structure();

    assert!(Arc::ptr_eq(&owned.shared_topology(), &topology));
    assert_eq!(owned, source);
    source
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
    structure.set_occupancy(ATOM0, Some(0.75)).unwrap();
    structure
        .set_b_factor(ATOM0, Some(Quantity::new(12.5, SQUARE_ANGSTROM)))
        .unwrap();
    // No reserved names: a generic property may use the same key.
    structure
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
        structure.set_occupancy(ATOM0, Some(f64::NAN)),
        Err(GeometryError::Array(DenseArrayError::NonFinite {
            index: 0
        }))
    );
    assert!(matches!(
        structure.set_b_factor(ATOM0, Some(Quantity::new(1.0, KELVIN))),
        Err(GeometryError::Unit(_))
    ));
    assert_eq!(
        structure.set_occupancy(TopologyAtomIndex::new(2), Some(1.0)),
        Err(GeometryError::Array(DenseArrayError::IndexOutOfRange {
            index: 2,
            len: 2
        }))
    );
    assert_eq!(
        structure.set_occupancies(Some(vec![Some(1.0)])),
        Err(GeometryError::Array(DenseArrayError::CountMismatch {
            expected: 2,
            actual: 1
        }))
    );
    assert_eq!(structure, before);

    // Clearing every value compares equal to never recording any.
    let mut cleared = structure.clone();
    cleared.set_occupancy(ATOM0, None).unwrap();
    cleared.set_occupancies(None).unwrap();
    assert_eq!(cleared.occupancies(), None);
}

#[test]
fn model_subset_projects_entity_state_and_drops_owner_properties() {
    let (mut structure, atom, bond) = structure_fixture();
    structure.set_occupancy(ATOM0, Some(0.5)).unwrap();
    let mut properties = structure.properties_mut();
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

#[test]
fn set_and_sample_properties_are_separate() {
    let (mut structure, atom, _) = structure_fixture();
    let mut properties = structure.properties_mut();
    properties
        .atoms_mut()
        .set_value(key("member_atom"), ATOM0, Some(PropertyValue::Int(1)))
        .unwrap();
    properties
        .bonds_mut()
        .set_value(key("member_bond"), BOND0, Some(PropertyValue::Int(2)))
        .unwrap();
    let mut set = bundle([structure]).unwrap();
    set.properties_mut()
        .insert(key("collection"), PropertyValue::Bool(true))
        .unwrap();
    set.get_mut(0)
        .unwrap()
        .properties_mut()
        .owner_mut()
        .insert(key("member"), PropertyValue::Int(1))
        .unwrap();
    assert!(!set.properties().is_empty());
    let member = set.get(0).unwrap();
    assert!(!member.properties().owner().is_empty());
    assert_eq!(
        member
            .properties()
            .atoms()
            .value(&key("member_atom"), ATOM0)
            .unwrap(),
        Some(PropertyValue::Int(1))
    );
    assert_eq!(
        member
            .properties()
            .bonds()
            .value(&key("member_bond"), BOND0)
            .unwrap(),
        Some(PropertyValue::Int(2))
    );

    let cloned = set.clone();
    assert_eq!(cloned.properties(), set.properties());
    let selection =
        crate::topology::AtomSelection::from_atoms(set.shared_topology(), [atom]).unwrap();
    let subset = set.subset(&selection).unwrap();
    assert!(subset.properties().is_empty());
    assert!(subset.get(0).unwrap().properties().owner().is_empty());
}

#[test]
fn ordered_samples_join_a_distribution_without_copying_positions() {
    let (structure, _, _) = structure_fixture();
    let mut geometry = Geometry::new(structure.positions().clone());
    geometry.set_velocities(Some(Velocities::zeros(2))).unwrap();
    assert_eq!(
        geometry.set_forces(Some(Forces::zeros(3))),
        Err(GeometryError::Array(DenseArrayError::CountMismatch {
            expected: 2,
            actual: 3
        }))
    );
    let mut set = StructureSet::from_geometries(structure.shared_topology(), [geometry]).unwrap();
    set.set_order(Some(Order::Step(vec![4]))).unwrap();
    let pointer = set.get(0).unwrap().positions().values().value().as_ptr();
    set.set_distribution(Some(
        Distribution::new(Members::All, Partition::one_state(1)).unwrap(),
    ))
    .unwrap();
    let member = set.get(0).unwrap();
    // Geometries join as equally weighted members and keep their place.
    assert_eq!(member.weight(), Some(1.0));
    assert_eq!(member.place(), Some(Place::Step(4)));
    assert_eq!(member.positions().values().value().as_ptr(), pointer);
}

#[test]
fn sets_read_chosen_samples_of_another_set_without_consuming_it() {
    let (structure, carbon, _) = structure_fixture();
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
    let trajectory = StructureSet::from_geometries(structure.shared_topology(), frames).unwrap();

    // Leave out equilibration, take every second geometry, and weight them.
    let mut sampled = StructureSet::from_reader(&trajectory, (1..5).step_by(2)).unwrap();
    sampled
        .set_distribution(Some(
            Distribution::new(Members::All, Partition::one_state(2))
                .unwrap()
                .with_weights(vec![3.0, 1.0])
                .unwrap(),
        ))
        .unwrap();
    assert_eq!(trajectory.len(), 5);
    assert_eq!(sampled.distribution().unwrap().state_count(), 1);
    assert!(sampled.state(0).unwrap().energy().is_none());
    let geometries = sampled
        .iter()
        .map(|geometry| {
            (
                geometry.positions().values().value()[0].x,
                geometry.weight(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(geometries, [(1.0, Some(3.0)), (3.0, Some(1.0))]);
    assert!(sampled
        .get(0)
        .unwrap()
        .shares_structure_with(&trajectory.get(1).unwrap()));

    let atoms =
        crate::topology::AtomSelection::from_atoms(structure.shared_topology(), [carbon]).unwrap();
    let subset = StructureSet::from_reader(&trajectory, [4])
        .unwrap()
        .subset(&atoms)
        .unwrap();
    assert_eq!(subset.topology().atom_count(), 1);
    assert_eq!(subset.get(0).unwrap().positions().len(), 1);

    assert_eq!(
        StructureSet::from_reader(&trajectory, [5]).unwrap_err(),
        StructureSetReadError::IndexOutOfRange { index: 5, len: 5 }
    );
}

/// A reader of three geometries on a step order, weighted 1, 2, and 3, that
/// decodes a fresh geometry on every read, even past its length, and counts
/// its reads.
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
        set.set_distribution(Some(
            Distribution::new(Members::All, Partition::one_state(3))
                .unwrap()
                .with_weights(vec![1.0, 2.0, 3.0])
                .unwrap(),
        ))
        .unwrap();
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

    fn distribution(&self) -> Option<&Distribution> {
        self.set.distribution()
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
fn sets_read_from_a_reader_are_faithful_copies_reading_each_sample_once() {
    let reader = DecodingReader::new();

    // Each picked geometry whole, with its order value and weight; a repeated
    // index is read once and shared.
    let set = StructureSet::from_reader(&reader, [2, 0, 2, 2]).unwrap();
    assert_eq!(reader.reads(), 2);
    assert_eq!(set.order(), Some(&Order::Step(vec![12, 10, 12, 12])));
    assert_eq!(
        set.distribution().unwrap().weights(),
        Some([3.0, 1.0, 3.0, 3.0].as_slice())
    );
    let first = set.get(0).unwrap();
    for index in [2, 3] {
        assert!(set.get(index).unwrap().shares_structure_with(&first));
    }
    assert!(!set.get(1).unwrap().shares_structure_with(&first));
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
    // values and weights it does not have.
    assert_eq!(
        StructureSet::from_reader(&reader, [1, 4]).unwrap_err(),
        StructureSetReadError::IndexOutOfRange { index: 4, len: 3 }
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
fn clearing_velocities_and_forces_keeps_samples_without_them_shared() {
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
    assert!(set.get(1).unwrap().shares_structure_with(&bare));
    assert_eq!(
        set.get(0).unwrap().occupancies(),
        Some([Some(0.5), None].as_slice())
    );
}

#[test]
fn occupancies_and_b_factors_follow_models_into_samples() {
    let (mut structure, _, _) = structure_fixture();
    structure.set_occupancy(ATOM0, Some(0.5)).unwrap();
    structure
        .set_b_factor(ATOM0, Some(Quantity::new(20.0, SQUARE_ANGSTROM)))
        .unwrap();
    structure.set_cell(Some(cell()));

    // Structures become members with their occupancies, B-factors, cell, and rows.
    let set = bundle([structure.clone()]).unwrap();
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
fn model_atom_views_join_topology_identity_with_structure_state() {
    let (mut structure, carbon, _) = structure_fixture();
    structure.set_occupancy(ATOM0, Some(0.5)).unwrap();
    structure
        .set_b_factor(ATOM0, Some(Quantity::new(20.0, SQUARE_ANGSTROM)))
        .unwrap();
    structure
        .properties_mut()
        .atoms_mut()
        .set_value(key("charge"), ATOM0, Some(PropertyValue::Int(-1)))
        .unwrap();

    let atom = structure.atom(carbon).unwrap();
    assert_eq!(atom.id(), carbon);
    assert_eq!(atom.index(), ATOM0);
    assert_eq!(atom.element.symbol(), "C");
    assert_eq!(atom.position().into_value(), Point3::origin());
    assert_eq!(atom.occupancy(), Some(0.5));
    let b_factor = atom.b_factor().unwrap().value_in(SQUARE_ANGSTROM).unwrap();
    assert!((b_factor - 20.0).abs() < 1.0e-12);
    assert_eq!(
        atom.structure_property(&key("charge")),
        Some(PropertyValue::Int(-1))
    );
    // Static annotations stay on the topology view.
    assert_eq!(atom.property(&key("charge")), None);

    let oxygen = structure.atom_at(TopologyAtomIndex::new(1)).unwrap();
    assert_eq!(oxygen.element.symbol(), "O");
    assert_eq!(oxygen.occupancy(), None);
    assert_eq!(oxygen.b_factor(), None);
    assert_eq!(
        oxygen.neighbors().map(|atom| atom.id()).collect::<Vec<_>>(),
        [carbon]
    );
    assert_eq!(
        structure
            .atoms()
            .map(|atom| atom.index())
            .collect::<Vec<_>>(),
        [ATOM0, TopologyAtomIndex::new(1)]
    );
    assert!(structure.atom_at(TopologyAtomIndex::new(2)).is_none());
    let missing =
        crate::topology::InstanceAtomId::new(carbon.instance(), crate::molecule::AtomId::new(9));
    assert!(structure.atom(missing).is_none());
    let view = structure.as_structure_view();
    assert_eq!(view.atom(carbon).unwrap().position(), atom.position());
}

#[test]
fn member_weights_are_finite_positive_and_relative() {
    let topology = Arc::new(single_atom_topology());
    let mut set = weighted(Arc::clone(&topology), &[3.0, 1.0]);
    for weight in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            set.distribution_mut().unwrap().set_weight(1, weight),
            Err(StructureSetError::InvalidWeight { geometry: 1 })
        );
    }

    // Weights are relative, so selection keeps them as stored.
    let resampled = set.select([0, 0, 1]).unwrap();
    assert_eq!(weights(&resampled), [3.0, 3.0, 1.0]);

    set.distribution_mut().unwrap().normalize_weights().unwrap();
    assert_eq!(weights(&set), [0.75, 0.25]);
    // Scaling by the largest weight first keeps huge weights finite.
    let mut huge = weighted(Arc::clone(&topology), &[f64::MAX, f64::MAX]);
    huge.distribution_mut()
        .unwrap()
        .normalize_weights()
        .unwrap();
    assert_eq!(weights(&huge), [0.5, 0.5]);
    // A ratio that underflows is rejected without changing any weight.
    let mut extreme = weighted(Arc::clone(&topology), &[f64::MAX, f64::MIN_POSITIVE]);
    assert_eq!(
        extreme.distribution_mut().unwrap().normalize_weights(),
        Err(StructureSetError::InvalidWeight { geometry: 1 })
    );
    assert_eq!(extreme.get(1).unwrap().weight(), Some(f64::MIN_POSITIVE));
    assert_eq!(
        with_empty_distribution(topology)
            .distribution_mut()
            .unwrap()
            .normalize_weights(),
        Err(StructureSetError::EmptyDistribution)
    );

    // Structures carry no statistical weight; they become equally weighted members.
    let (structure, _, _) = structure_fixture();
    let set = bundle([structure.clone(), structure]).unwrap();
    assert_eq!(set.distribution().unwrap().weights(), None);
    assert_eq!(weights(&set), [1.0, 1.0]);
}

/// Unit-weight conformers of a one-atom topology, one state each, with
/// state energies in kJ/mol.
fn energy_conformers(energies: &[f64]) -> StructureSet {
    let states = energies
        .iter()
        .map(|energy| {
            State::new()
                .with_energy(Quantity::new(*energy, KILOJOULE_PER_MOLE))
                .unwrap()
        })
        .collect();
    let mut set = StructureSet::from_geometries(
        single_atom_topology(),
        energies.iter().map(|_| Geometry::new(single_position(0.0))),
    )
    .unwrap();
    let distribution = Distribution::new(Members::All, Partition::state_per_sample(energies.len()))
        .unwrap()
        .with_weights(vec![1.0; energies.len()])
        .unwrap()
        .with_states(states)
        .unwrap();
    set.set_distribution(Some(distribution)).unwrap();
    set
}

/// Geometries of a one-atom topology with these weights, grouped by
/// `assignments` with `representatives`.
fn grouped(weights: &[f64], assignments: &[usize], representatives: &[usize]) -> StructureSet {
    let mut set = StructureSet::from_geometries(
        single_atom_topology(),
        weights.iter().map(|_| Geometry::new(single_position(0.0))),
    )
    .unwrap();
    let partition = Partition::new(assignments.to_vec(), representatives.to_vec()).unwrap();
    let distribution = Distribution::new(Members::All, partition)
        .unwrap()
        .with_weights(weights.to_vec())
        .unwrap();
    set.set_distribution(Some(distribution)).unwrap();
    set
}

/// Conformers of `topology` at the origin with these weights, one state
/// each.
fn weighted(topology: Arc<crate::topology::Topology>, weights: &[f64]) -> StructureSet {
    let geometries = weights.iter().map(|_| Geometry::new(single_position(0.0)));
    let mut set = conformers(topology, geometries).unwrap();
    set.distribution_mut()
        .unwrap()
        .set_weights(Some(weights.to_vec()))
        .unwrap();
    set
}

/// An empty set of `topology` with a distribution of no members.
fn with_empty_distribution(topology: impl Into<Arc<crate::topology::Topology>>) -> StructureSet {
    let mut set = StructureSet::new(topology);
    set.set_distribution(Some(
        Distribution::new(Members::All, Partition::one_state(0)).unwrap(),
    ))
    .unwrap();
    set
}

fn kj(value: f64) -> Option<Quantity<f64>> {
    Some(Quantity::new(value, KILOJOULE_PER_MOLE))
}

/// Sets the energy of every state, in kJ/mol where given.
fn set_energies(set: &mut StructureSet, energies: &[Option<f64>]) {
    let mut distribution = set.distribution_mut().unwrap();
    for (state, energy) in energies.iter().enumerate() {
        distribution
            .state_mut(state)
            .unwrap()
            .set_energy(energy.and_then(kj))
            .unwrap();
    }
}

fn weights(set: &StructureSet) -> Vec<f64> {
    set.iter()
        .map(|geometry| geometry.weight().unwrap())
        .collect()
}

fn energies(set: &StructureSet) -> Vec<Option<Quantity<f64>>> {
    set.states().map(|state| state.energy()).collect()
}

fn representatives(set: &StructureSet) -> Vec<usize> {
    set.distribution()
        .unwrap()
        .partition()
        .representatives()
        .to_vec()
}

fn temperature(set: &StructureSet) -> Option<Quantity<f64>> {
    set.distribution().unwrap().temperature()
}

#[test]
fn state_energies_are_canonical_finite_and_never_on_samples() {
    let state = State::new();
    assert_eq!(state.energy(), None);
    let mut state = state
        .with_energy(Quantity::new(1.0, KILOCALORIE_PER_MOLE))
        .unwrap();
    assert_eq!(state.energy(), kj(4.184));

    // Invalid energies are rejected without changing the stored one.
    for value in [f64::NAN, f64::INFINITY] {
        assert_eq!(
            state.set_energy(Some(Quantity::new(value, KILOJOULE_PER_MOLE))),
            Err(GeometryError::NonFiniteEnergy)
        );
    }
    assert_eq!(
        state.set_energy(Some(Quantity::new(1.0, KELVIN))),
        Err(GeometryError::Unit(Box::new(
            UnitError::IncompatibleUnits {
                from: KELVIN,
                to: KILOJOULE_PER_MOLE,
            }
        )))
    );
    assert_eq!(state.energy(), kj(4.184));

    let mut set = with_empty_distribution(single_atom_topology());
    let index = set
        .push_state(Geometry::new(single_position(0.0)), None, state, Some(2.0))
        .unwrap();
    assert_eq!(index, 0);
    set.distribution_mut()
        .unwrap()
        .state_mut(0)
        .unwrap()
        .set_energy(kj(-3.0))
        .unwrap();
    // Energies describe the sampled structure, so subsets keep them.
    let all = crate::topology::AtomSelection::all(set.shared_topology());
    assert_eq!(energies(&set.subset(&all).unwrap()), [kj(-3.0)]);
    assert_eq!(set.get(0).unwrap().state().unwrap().energy(), kj(-3.0));
    let mut distribution = set.distribution_mut().unwrap();
    distribution.state_mut(0).unwrap().set_energy(None).unwrap();
    assert!(distribution.state_mut(1).is_none());
    assert_eq!(energies(&set), [None]);

    // Structures carry no energy.
    let (structure, _, _) = structure_fixture();
    assert_eq!(energies(&bundle([structure]).unwrap()), [None]);
}

#[test]
fn every_member_is_in_one_state_that_contains_its_representative() {
    let mut set = with_empty_distribution(single_atom_topology());
    assert_eq!(set.distribution().unwrap().state_count(), 0);
    let geometry = |x| Geometry::new(single_position(x));
    assert_eq!(
        set.push_member(geometry(0.0), None, 0, None),
        Err(StructureSetError::StateOutOfRange {
            state: 0,
            states: 0
        })
    );
    set.push_state(geometry(0.0), None, State::new(), None)
        .unwrap();
    set.push_state(geometry(1.0), None, State::new(), None)
        .unwrap();
    set.push_member(geometry(2.0), None, 0, None).unwrap();
    assert_eq!(
        set.distribution().unwrap().partition().assignments(),
        [0, 1, 0]
    );
    assert_eq!(representatives(&set), [0, 1]);
    let state = set.state(0).unwrap();
    assert_eq!(state.geometries().collect::<Vec<_>>(), [0, 2]);
    assert_eq!(state.representative_sample().index(), 0);
    assert!(set.get(0).unwrap().is_representative());
    assert!(!set.get(2).unwrap().is_representative());
    assert_eq!(set.get(2).unwrap().state().unwrap().index(), 0);

    // A representative of a state with other members stays.
    assert_eq!(
        set.remove(0),
        Err(StructureSetError::RepresentativeRemoval {
            geometry: 0,
            state: 0
        })
    );
    // Removing a state's only member removes the state; later indices shift.
    set.remove(1).unwrap();
    assert_eq!(
        set.distribution().unwrap().partition().assignments(),
        [0, 0]
    );
    assert_eq!(representatives(&set), [0]);
    set.remove(1).unwrap();
    set.remove(0).unwrap();
    assert!(set.is_empty());
    assert_eq!(set.distribution().unwrap().state_count(), 0);

    // A distribution of every geometry covers all of them.
    let mut two =
        StructureSet::from_geometries(single_atom_topology(), [geometry(0.0), geometry(1.0)])
            .unwrap();
    assert_eq!(
        two.set_distribution(Some(
            Distribution::new(Members::All, Partition::one_state(3)).unwrap()
        )),
        Err(StructureSetError::CoverageMismatch {
            covered: 3,
            geometries: 2
        })
    );
    assert_eq!(
        Partition::new(vec![0, 1], vec![1, 0]),
        Err(PartitionError::RepresentativeNotInState {
            state: 0,
            sample: 1
        })
    );
    assert_eq!(
        Partition::new(vec![0, 2], vec![0, 1]),
        Err(PartitionError::StateOutOfRange {
            sample: 1,
            state: 2,
            states: 2
        })
    );
    assert_eq!(
        Partition::new(vec![0], vec![0, 3]),
        Err(PartitionError::RepresentativeOutOfRange {
            state: 1,
            sample: 3,
            samples: 1
        })
    );
    assert_eq!(Partition::one_state(0).state_count(), 0);
    assert_eq!(Partition::one_state(2).representatives(), [0]);
}

#[test]
fn populations_are_derived_from_member_weights() {
    // Two states: geometries 0 and 2 (weights 1 and 3) and geometry 1 (weight 4).
    let mut set = grouped(&[1.0, 4.0, 3.0], &[0, 1, 0], &[2, 1]);
    // Populations are relative to the largest member weight, 4.
    assert_eq!(set.state(0).unwrap().population(), 1.0);
    assert_eq!(set.state(1).unwrap().population(), 1.0);
    let distribution = set.distribution().unwrap();
    assert_eq!(distribution.probabilities().unwrap(), [0.125, 0.5, 0.375]);
    assert_eq!(distribution.state_probabilities().unwrap(), [0.5, 0.5]);
    // Kish: (1 + 4 + 3)² / (1 + 16 + 9) = 64 / 26.
    assert_eq!(distribution.effective_sample_size().unwrap(), 64.0 / 26.0);

    set.distribution_mut().unwrap().set_weight(1, 12.0).unwrap();
    assert_eq!(set.state(1).unwrap().population(), 1.0);
    let distribution = set.distribution().unwrap();
    assert_eq!(distribution.state_probabilities().unwrap(), [0.25, 0.75]);
    assert_eq!(distribution.ranked_by_population(), [1, 0]);
    set.distribution_mut().unwrap().normalize_weights().unwrap();
    assert_eq!(weights(&set), [0.0625, 0.75, 0.1875]);
    assert_eq!(
        set.distribution().unwrap().state_probabilities().unwrap(),
        [0.25, 0.75]
    );

    let empty = with_empty_distribution(single_atom_topology());
    assert_eq!(
        empty.distribution().unwrap().state_probabilities(),
        Err(StructureSetError::EmptyDistribution)
    );
}

#[test]
fn boltzmann_populations_follow_relative_state_energies_at_temperature() {
    // 1, 0, and 2 kcal/mol above the minimum at 298.15 K. Expected factors
    // exp(-dE/RT) evaluated independently to 40 digits (Python `decimal`,
    // R = N_A k_B exact in the 2019 SI).
    let mut set = energy_conformers(&[-100.0 + 4.184, -100.0, -100.0 + 8.368]);
    set.distribution_mut()
        .unwrap()
        .assign_boltzmann_populations(Quantity::new(298.15, KELVIN))
        .unwrap();
    let expected = [0.184_924_695_965_256_4, 1.0, 0.034_197_143_177_842_51];
    let weights = weights(&set);
    for (weight, expected) in weights.iter().zip(expected) {
        let ratio = weight / weights[1];
        assert!((ratio - expected).abs() <= 1.0e-15, "{ratio} != {expected}");
    }
    // The total weight is kept.
    assert!((weights.iter().sum::<f64>() - 3.0).abs() <= 1.0e-15);
    assert_eq!(temperature(&set), Some(Quantity::new(298.15, KELVIN)));

    // Invalid temperatures and missing energies change nothing.
    let before = self::weights(&set);
    let mut distribution = set.distribution_mut().unwrap();
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            distribution.assign_boltzmann_populations(Quantity::new(temperature, KELVIN)),
            Err(StructureSetError::InvalidTemperature)
        );
    }
    assert_eq!(
        distribution.assign_boltzmann_populations(Quantity::new(1.0, KILOJOULE_PER_MOLE)),
        Err(StructureSetError::Unit(Box::new(
            UnitError::IncompatibleUnits {
                from: KILOJOULE_PER_MOLE,
                to: KELVIN,
            }
        )))
    );
    distribution.state_mut(2).unwrap().set_energy(None).unwrap();
    assert_eq!(
        distribution.assign_boltzmann_populations(Quantity::new(310.0, KELVIN)),
        Err(StructureSetError::Missing(Requirement::StateEnergies))
    );
    assert_eq!(self::weights(&set), before);
    assert_eq!(temperature(&set), Some(Quantity::new(298.15, KELVIN)));

    // A factor that underflows to zero is rejected transactionally.
    let mut wide = energy_conformers(&[0.0, 1.0e6]);
    assert_eq!(
        wide.distribution_mut()
            .unwrap()
            .assign_boltzmann_populations(Quantity::new(298.15, KELVIN)),
        Err(StructureSetError::InvalidWeight { geometry: 1 })
    );
    assert_eq!(self::weights(&wide), [1.0, 1.0]);
    assert_eq!(temperature(&wide), None);
}

#[test]
fn boltzmann_populations_scale_the_members_of_each_state_together() {
    // State 0 (0 kJ/mol) holds geometries of weights 1 and 3; state 1, RT ln 3
    // higher, one geometry of weight 2. Populations go from 4 : 2 to 3 : 1 of
    // the kept total 6, so 4.5 and 1.5: state 0's geometries scale by 4.5 / 4
    // and state 1's by 1.5 / 2.
    let kelvin = 300.0;
    let gap = 8.314_462_618_153_24e-3 * kelvin * 3.0_f64.ln();
    let mut set = grouped(&[1.0, 2.0, 3.0], &[0, 1, 0], &[0, 1]);
    set_energies(&mut set, &[Some(0.0), Some(gap)]);
    set.distribution_mut()
        .unwrap()
        .assign_boltzmann_populations(Quantity::new(kelvin, KELVIN))
        .unwrap();
    for (weight, expected) in weights(&set).into_iter().zip([1.125, 1.5, 3.375]) {
        assert!(
            (weight - expected).abs() <= 1.0e-14,
            "{weight} != {expected}"
        );
    }
    assert_eq!(
        set.distribution().unwrap().partition().assignments(),
        [0, 1, 0]
    );

    // States whose weights differ by more than the floating-point range
    // still get representable populations.
    let mut spread = grouped(&[1.0e-200, 1.0e200], &[0, 1], &[0, 1]);
    set_energies(&mut spread, &[Some(0.0), Some(0.0)]);
    spread
        .distribution_mut()
        .unwrap()
        .assign_boltzmann_populations(Quantity::new(kelvin, KELVIN))
        .unwrap();
    for weight in weights(&spread) {
        assert!((weight / 5.0e199 - 1.0).abs() <= 1.0e-12, "{weight}");
    }

    // Weighting one state changes nothing: its population is all there is.
    let mut one = grouped(&[1.0, 2.0, 5.0], &[0, 0, 0], &[1]);
    set_energies(&mut one, &[Some(7.0)]);
    one.distribution_mut()
        .unwrap()
        .assign_boltzmann_populations(Quantity::new(298.15, KELVIN))
        .unwrap();
    assert_eq!(weights(&one), [1.0, 2.0, 5.0]);
    assert_eq!(temperature(&one), Some(Quantity::new(298.15, KELVIN)));

    // Geometries drawn from a simulation form a state without an energy, which
    // is never weighted a second time.
    let (structure, _, _) = structure_fixture();
    let mut drawn = bundle([structure.clone(), structure]).unwrap();
    let mut distribution = drawn.distribution_mut().unwrap();
    assert_eq!(
        distribution.assign_boltzmann_populations(Quantity::new(298.15, KELVIN)),
        Err(StructureSetError::Missing(Requirement::StateEnergies))
    );
    assert_eq!(distribution.temperature(), None);

    // A temperature is recorded or cleared explicitly, and validated.
    distribution
        .set_temperature(Some(Quantity::new(298.15, KELVIN)))
        .unwrap();
    assert_eq!(
        distribution.temperature(),
        Some(Quantity::new(298.15, KELVIN))
    );
    assert_eq!(
        distribution.set_temperature(Some(Quantity::new(-1.0, KELVIN))),
        Err(StructureSetError::InvalidTemperature)
    );
    distribution.set_temperature(None).unwrap();
    assert_eq!(temperature(&drawn), None);
}

/// Whether every geometry of `a` shares its structure with the same geometry
/// of `b`.
fn shares_every_structure(a: &StructureSet, b: &StructureSet) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b.iter())
            .all(|(a, b)| a.payload().shares_structure_with(b.payload()))
}

fn x(set: &StructureSet, geometry: usize) -> f64 {
    set.get(geometry).unwrap().positions().values().value()[0].x
}

#[test]
fn regrouping_and_reweighting_a_shared_set_copy_no_structure() {
    let shared = Arc::new(energy_conformers(&[0.0, 1.0, 2.0, 3.0]));
    let mut edited = Arc::clone(&shared);

    let mut distribution = Arc::make_mut(&mut edited).distribution_mut().unwrap();
    distribution
        .regroup(Partition::new(vec![0, 0, 1, 1], vec![0, 2]).unwrap())
        .unwrap();
    distribution.set_weight(1, 3.0).unwrap();
    distribution.normalize_weights().unwrap();
    distribution
        .state_mut(0)
        .unwrap()
        .set_energy(kj(0.0))
        .unwrap();
    distribution
        .state_mut(1)
        .unwrap()
        .set_energy(kj(1.0))
        .unwrap();
    distribution
        .assign_boltzmann_populations(Quantity::new(300.0, KELVIN))
        .unwrap();

    assert!(shares_every_structure(&edited, &shared));
    // The held set keeps its own grouping and weights.
    assert_eq!(shared.distribution().unwrap().state_count(), 4);
    assert_eq!(weights(&shared), [1.0; 4]);
    assert_eq!(temperature(&shared), None);
    assert_eq!(edited.distribution().unwrap().state_count(), 2);
    assert_ne!(weights(&edited), weights(&shared));

    // Selections and collapsing share the structures they keep.
    let selected = shared.select([3, 1]).unwrap();
    assert!(selected
        .get(0)
        .unwrap()
        .payload()
        .shares_structure_with(shared.get(3).unwrap().payload()));
    assert!(shares_every_structure(
        &shared.select_states(0..4).unwrap(),
        &shared
    ));
    assert!(shares_every_structure(&shared.collapse().unwrap(), &shared));

    // A geometry taken from a set of the same topology stays shared.
    let mut other = with_empty_distribution(shared.shared_topology());
    let geometry = shared.get(2).unwrap().payload().clone();
    other
        .push_state(geometry, None, State::new(), None)
        .unwrap();
    assert!(other
        .get(0)
        .unwrap()
        .payload()
        .shares_structure_with(shared.get(2).unwrap().payload()));
}

#[test]
fn editing_one_sample_of_a_shared_set_copies_that_sample_alone() {
    let original = energy_conformers(&[0.0, 1.0, 2.0]);
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
                .shares_structure_with(original.get(geometry).unwrap().payload()),
            shared,
            "geometry {geometry}"
        );
    }

    // Transforms and position replacement leave the original unchanged.
    let mut transformed = original.clone();
    use crate::math::{Matrix3, RigidTransform};
    let shift =
        RigidTransform::from_canonical(Matrix3::identity(), Vector3::new(2.0, 0.0, 0.0)).unwrap();
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

#[test]
fn regrouping_keeps_weights_and_the_energies_of_whole_states() {
    // Four conformers, energies 0..3 kJ/mol, weights 1..4.
    let mut set = energy_conformers(&[0.0, 1.0, 2.0, 3.0]);
    let mut distribution = set.distribution_mut().unwrap();
    distribution
        .set_weights(Some(vec![1.0, 2.0, 3.0, 4.0]))
        .unwrap();
    distribution
        .state_mut(1)
        .unwrap()
        .properties_mut()
        .insert(key("label"), PropertyValue::Int(1))
        .unwrap();
    distribution
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    let before = set.clone();

    // An unchanged grouping keeps everything.
    let unchanged = set.distribution().unwrap().partition().clone();
    set.distribution_mut().unwrap().regroup(unchanged).unwrap();
    assert_eq!(energies(&set), energies(&before));
    assert_eq!(
        set.state(1).unwrap().properties(),
        before.state(1).unwrap().properties()
    );

    // Merging conformers 0 and 1 under 1 keeps conformer 1's state.
    let merged = Partition::new(vec![0, 0, 1, 2], vec![1, 2, 3]).unwrap();
    set.distribution_mut().unwrap().regroup(merged).unwrap();
    assert_eq!(weights(&set), [1.0, 2.0, 3.0, 4.0]);
    assert_eq!(energies(&set), [kj(1.0), kj(2.0), kj(3.0)]);
    assert_eq!(
        set.state(0).unwrap().properties().get(&key("label")),
        Some(&PropertyValue::Int(1))
    );
    // Relative to the largest member weight, 4.
    assert_eq!(set.state(0).unwrap().population(), 0.75);
    assert_eq!(temperature(&set), Some(Quantity::new(300.0, KELVIN)));
    // Geometries never move.
    assert_eq!(
        set.iter()
            .map(|geometry| geometry.positions().clone())
            .collect::<Vec<_>>(),
        before
            .iter()
            .map(|geometry| geometry.positions().clone())
            .collect::<Vec<_>>()
    );

    // Splitting the merged state loses its energy and properties; the state
    // that takes a whole old state with part of a split one keeps the whole
    // one's energy when its representative comes from it.
    let split = Partition::new(vec![0, 1, 0, 2], vec![2, 1, 3]).unwrap();
    set.distribution_mut().unwrap().regroup(split).unwrap();
    assert_eq!(weights(&set), [1.0, 2.0, 3.0, 4.0]);
    assert_eq!(energies(&set), [kj(2.0), None, kj(3.0)]);
    assert!(set.state(1).unwrap().properties().is_empty());

    // A partition of another member count changes nothing.
    let snapshot = set.clone();
    assert_eq!(
        set.distribution_mut()
            .unwrap()
            .regroup(Partition::one_state(3)),
        Err(StructureSetError::PartitionMismatch {
            expected: 4,
            actual: 3
        })
    );
    assert_eq!(set.distribution(), snapshot.distribution());

    // One state split into several loses its energy.
    let mut one = grouped(&[1.0, 1.0, 1.0, 1.0], &[0, 0, 0, 0], &[0]);
    set_energies(&mut one, &[Some(5.0)]);
    one.distribution_mut()
        .unwrap()
        .regroup(Partition::new(vec![0, 0, 1, 1], vec![0, 2]).unwrap())
        .unwrap();
    assert_eq!(energies(&one), [None, None]);
    assert_eq!(weights(&one), [1.0; 4]);
}

#[test]
fn collapse_keeps_one_representative_per_state_with_its_population() {
    let mut set = grouped(&[1.0, 2.0, 3.0, 4.0], &[0, 1, 0, 1], &[2, 1]);
    set.replace_positions([1.0, 2.0, 3.0, 4.0].map(single_position))
        .unwrap();
    set_energies(&mut set, &[Some(0.0), None]);
    set.distribution_mut()
        .unwrap()
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    let source = set.clone();

    let collapsed = set.collapse().unwrap();
    assert_eq!(collapsed.len(), 2);
    // Populations 4 and 6, relative to the largest member weight, 4.
    assert_eq!(weights(&collapsed), [1.0, 1.5]);
    assert_eq!(
        collapsed.get(0).unwrap().positions(),
        set.get(2).unwrap().positions()
    );
    assert_eq!(
        collapsed.get(1).unwrap().positions(),
        set.get(1).unwrap().positions()
    );
    assert_eq!(
        collapsed.distribution().unwrap().partition(),
        &Partition::state_per_sample(2)
    );
    assert_eq!(energies(&collapsed), [kj(0.0), None]);
    assert_eq!(temperature(&collapsed), temperature(&set));
    for (collapsed, source) in collapsed
        .distribution()
        .unwrap()
        .state_probabilities()
        .unwrap()
        .into_iter()
        .zip(set.distribution().unwrap().state_probabilities().unwrap())
    {
        assert!((collapsed - source).abs() <= 1.0e-15);
    }
    // The source is unchanged.
    assert_eq!(set.distribution(), source.distribution());
    assert_eq!(
        with_empty_distribution(single_atom_topology())
            .collapse()
            .unwrap()
            .len(),
        0
    );
}

#[test]
fn selection_copies_samples_and_whole_or_partial_states() {
    // State 0: geometries 0 and 1, represented by 1; state 1: geometry 2.
    let mut set = grouped(&[1.0, 2.0, 3.0], &[0, 0, 1], &[1, 2]);
    set_energies(&mut set, &[Some(1.0), Some(2.0)]);

    // Weights are relative, so selection keeps them as stored.
    let whole = set.select([2, 0, 1, 1]).unwrap();
    assert_eq!(weights(&whole), [3.0, 1.0, 2.0, 2.0]);
    assert_eq!(
        whole.distribution().unwrap().partition().assignments(),
        [1, 0, 0, 0]
    );
    assert_eq!(representatives(&whole), [2, 0]);
    assert_eq!(energies(&whole), [kj(1.0), kj(2.0)]);

    // A partly selected state has no energy.
    let partial = set.select([0, 2]).unwrap();
    assert_eq!(
        partial.distribution().unwrap().partition().assignments(),
        [0, 1]
    );
    assert_eq!(representatives(&partial), [0, 1]);
    assert_eq!(energies(&partial), [None, kj(2.0)]);

    // Whole states keep their members, representative, and energy.
    let states = set.select_states([1, 0]).unwrap();
    assert_eq!(weights(&states), [3.0, 1.0, 2.0]);
    assert_eq!(
        states.distribution().unwrap().partition().assignments(),
        [0, 1, 1]
    );
    assert_eq!(representatives(&states), [0, 2]);
    assert_eq!(energies(&states), [kj(2.0), kj(1.0)]);
    assert_eq!(
        set.select_states([2]).unwrap_err(),
        StructureSetError::StateOutOfRange {
            state: 2,
            states: 2
        }
    );
}

#[test]
fn population_ranking_and_coverage_select_the_most_populated_states() {
    let topology = Arc::new(single_atom_topology());
    let set = weighted(Arc::clone(&topology), &[1.0, 4.0, 2.0, 2.0, 1.0]);
    let distribution = set.distribution().unwrap();

    // Equal weights and populations keep geometry and state order.
    assert_eq!(distribution.ranked_by_weight(), [1, 2, 3, 0, 4]);
    assert_eq!(distribution.ranked_by_population(), [1, 2, 3, 0, 4]);
    assert_eq!(
        distribution.probabilities().unwrap(),
        [0.1, 0.4, 0.2, 0.2, 0.1]
    );
    assert_eq!(weights(&set), [1.0, 4.0, 2.0, 2.0, 1.0]);

    assert_eq!(distribution.covering_population(0.4).unwrap(), [1]);
    assert_eq!(distribution.covering_population(0.5).unwrap(), [1, 2]);
    assert_eq!(
        distribution.covering_population(0.85).unwrap(),
        [1, 2, 3, 0]
    );
    assert_eq!(
        distribution.covering_population(1.0).unwrap(),
        [1, 2, 3, 0, 4]
    );
    for fraction in [0.0, -0.5, 1.5, f64::NAN] {
        assert_eq!(
            distribution.covering_population(fraction),
            Err(StructureSetError::InvalidFraction)
        );
    }
    // States of several members rank by their summed weight.
    let pooled = grouped(&[3.0, 2.0, 2.0], &[0, 1, 1], &[0, 1]);
    let pooled = pooled.distribution().unwrap();
    assert_eq!(pooled.ranked_by_weight(), [0, 1, 2]);
    assert_eq!(pooled.ranked_by_population(), [1, 0]);
    assert_eq!(pooled.covering_population(0.5).unwrap(), [1]);

    // A full fraction keeps states too light to change the total.
    let light = weighted(Arc::clone(&topology), &[1.0, 1.0e-20]);
    let light = light.distribution().unwrap();
    assert_eq!(light.covering_population(1.0).unwrap(), [0, 1]);
    assert_eq!(light.covering_population(0.999).unwrap(), [0]);

    let skewed = weighted(Arc::clone(&topology), &[2.0, 1.0, 1.0]);
    assert_eq!(
        skewed
            .distribution()
            .unwrap()
            .effective_sample_size()
            .unwrap(),
        8.0 / 3.0
    );
    let uniform = weighted(Arc::clone(&topology), &[7.0; 4]);
    assert_eq!(
        uniform
            .distribution()
            .unwrap()
            .effective_sample_size()
            .unwrap(),
        4.0
    );

    let empty = with_empty_distribution(topology);
    let empty = empty.distribution().unwrap();
    assert_eq!(empty.ranked_by_weight(), Vec::<usize>::new());
    assert_eq!(empty.ranked_by_population(), Vec::<usize>::new());
    assert_eq!(
        empty.probabilities(),
        Err(StructureSetError::EmptyDistribution)
    );
    assert_eq!(
        empty.effective_sample_size(),
        Err(StructureSetError::EmptyDistribution)
    );
    assert_eq!(
        empty.covering_population(0.5),
        Err(StructureSetError::EmptyDistribution)
    );
}

#[test]
fn energy_ranking_and_windows_are_relative_to_the_lowest_state() {
    let mut set = energy_conformers(&[10.0, 0.0, 30.0, 4.184, 0.0]);
    let distribution = set.distribution().unwrap();
    assert_eq!(distribution.ranked_by_energy().unwrap(), [1, 4, 3, 0, 2]);
    // The window bound is inclusive and converted to the canonical unit.
    assert_eq!(
        distribution
            .within_energy_window(Quantity::new(1.0, KILOCALORIE_PER_MOLE))
            .unwrap(),
        [1, 3, 4]
    );
    assert_eq!(
        distribution
            .within_energy_window(Quantity::new(0.0, KILOJOULE_PER_MOLE))
            .unwrap(),
        [1, 4]
    );
    let window = distribution
        .within_energy_window(Quantity::new(10.0, KILOJOULE_PER_MOLE))
        .unwrap();
    assert_eq!(window, [0, 1, 3, 4]);
    assert_eq!(
        energies(&set.select_states(window).unwrap()),
        [10.0, 0.0, 4.184, 0.0].map(kj)
    );

    for value in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            distribution.within_energy_window(Quantity::new(value, KILOJOULE_PER_MOLE)),
            Err(StructureSetError::InvalidEnergyWindow)
        );
    }
    set_energies(&mut set, &[Some(10.0), Some(0.0), Some(30.0), None]);
    let distribution = set.distribution().unwrap();
    assert_eq!(
        distribution.ranked_by_energy(),
        Err(StructureSetError::Missing(Requirement::StateEnergies))
    );
    assert_eq!(
        distribution.within_energy_window(Quantity::new(1.0, KILOJOULE_PER_MOLE)),
        Err(StructureSetError::Missing(Requirement::StateEnergies))
    );
}

#[test]
fn model_views_report_which_structure_they_borrow() {
    let topology = Arc::new(single_atom_topology());
    let structure = Structure::new(Arc::clone(&topology), single_position(0.0)).unwrap();
    assert_eq!(
        structure.as_structure_view().kind(),
        StructureKind::Standalone
    );

    let mut set = energy_conformers(&[2.5]);
    set.push_member(Geometry::new(single_position(0.0)), None, 0, Some(3.0))
        .unwrap();
    set.distribution_mut()
        .unwrap()
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    let StructureKind::InSet(context) = set.get(1).unwrap().as_structure_view().kind() else {
        panic!("a member is a geometry");
    };
    assert_eq!(context.place(), None);
    assert_eq!(
        context.membership(),
        Some(Membership::new(
            3.0,
            0,
            false,
            kj(2.5),
            Some(Quantity::new(300.0, KELVIN)),
        ))
    );
    assert_eq!(context.velocities(), None);
    let geometry = set.get(0).unwrap();
    assert!(matches!(
        geometry.as_structure_view().kind(),
        StructureKind::InSet(context)
            if context.membership().is_some_and(Membership::is_representative)
    ));
    // Materializing a geometry as a structure drops its geometry state.
    assert_eq!(
        geometry.to_structure().as_structure_view().kind(),
        StructureKind::Standalone
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
    assert_eq!(context.membership(), None);
    assert_eq!(context.forces(), None);
}

#[test]
fn rigid_transforms_move_models_and_samples_with_their_cells() {
    use crate::math::{Matrix3, RigidTransform};
    let topology = Arc::new(single_atom_topology());
    let box_cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(1.0, 2.0, 3.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    let mut structure = Structure::new(Arc::clone(&topology), single_position(1.0)).unwrap();
    structure.set_cell(Some(box_cell));
    // A quarter turn about z followed by a 2 nm shift along x.
    let turn = Matrix3::from_columns(
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    );
    let transform = RigidTransform::from_canonical(turn, Vector3::new(2.0, 0.0, 0.0)).unwrap();
    let mut members = bundle([structure.clone()]).unwrap();
    let mut trajectory = StructureSet::from_geometries(
        Arc::clone(&topology),
        [Geometry::new(single_position(1.0)).with_cell(box_cell)],
    )
    .unwrap();
    trajectory
        .set_order(Some(Order::Time(TimeAxis::default())))
        .unwrap();

    structure.apply_rigid_transform(&transform).unwrap();
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

mod structure_set;
