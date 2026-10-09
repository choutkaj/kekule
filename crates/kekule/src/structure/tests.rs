use super::*;
use crate::core::{Atom, BondOrder, Element, MoleculeEditor};
use crate::geometry::{PeriodicCell, Point3, Vector3};
use crate::properties::{PropertyColumn, PropertyError, PropertyKey, PropertyValue};
use crate::topology::{TopologyAtomIndex, TopologyBondIndex, TopologyBuilder};
use crate::units::{Quantity, ANGSTROM, KELVIN, NANOMETER, SQUARE_ANGSTROM, SQUARE_NANOMETER};
use std::sync::Arc;

fn key(name: &str) -> PropertyKey {
    PropertyKey::new(name).unwrap()
}

const ATOM0: TopologyAtomIndex = TopologyAtomIndex::new(0);
const BOND0: TopologyBondIndex = TopologyBondIndex::new(0);

fn model_fixture() -> (
    Model,
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
        Model::new(topology, positions).unwrap(),
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
fn canonical_model_constructors_accept_owned_and_shared_topology_and_conformations() {
    let owned = Model::new(single_atom_topology(), single_position(1.0)).unwrap();
    let shared = owned.shared_topology();
    let shared_model = Model::new(Arc::clone(&shared), single_position(2.0)).unwrap();
    assert!(Arc::ptr_eq(&shared_model.shared_topology(), &shared));

    let conformation = Conformation::new(single_position(3.0)).with_cell(cell());
    let complete = Model::new(Arc::clone(&shared), conformation.clone()).unwrap();
    assert_eq!(complete.cell(), conformation.cell());
    assert_eq!(complete.positions(), conformation.positions());
    let (topology, recovered) = complete.into_parts();
    assert!(Arc::ptr_eq(&topology, &shared));
    assert_eq!(recovered.positions(), conformation.positions());
}

#[test]
fn model_builder_rejects_instance_position_mismatch_before_mutation() {
    let topology = single_atom_topology();
    let molecule = topology.molecules().next().unwrap();
    let mut builder = ModelBuilder::new();
    let definition = builder
        .add_molecule_definition(molecule.molecule().clone())
        .unwrap();
    let before = builder.clone();
    assert_eq!(
        builder.add_instance(definition, &Positions::zeros(0)),
        Err(ModelBuildError::InstancePositionCountMismatch {
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
    let molecule = topology.molecules().next().unwrap();
    let mut builder = ModelBuilder::new();
    builder
        .add_molecule(molecule.molecule().clone(), &single_position(1.0))
        .unwrap();
    let before = builder.clone();
    assert_eq!(
        builder
            .conformation_mut()
            .set_positions(Quantity::new(vec![Point3::origin(); 2], NANOMETER)),
        Err(ConformationError::Position(
            PositionError::PositionCountMismatch {
                expected: 1,
                actual: 2
            }
        ))
    );
    assert_eq!(builder, before);
    assert_eq!(builder.build().unwrap().positions(), &single_position(1.0));
}

#[test]
fn model_rejects_conformations_with_other_dimensions() {
    let (model, _, _) = model_fixture();
    assert_eq!(
        Model::new(model.shared_topology(), single_position(1.0)).unwrap_err(),
        ModelError::Conformation(ConformationError::AtomCountMismatch {
            expected: 2,
            actual: 1
        })
    );
    // A conformation bound to a topology with one bond cannot bind to one
    // with no bonds: populated bond rows never resize silently.
    let two_atoms = Arc::new(
        crate::topology::Topology::from_molecules([
            single_atom_topology()
                .molecules()
                .next()
                .unwrap()
                .molecule()
                .clone(),
            single_atom_topology()
                .molecules()
                .next()
                .unwrap()
                .molecule()
                .clone(),
        ])
        .unwrap(),
    );
    let bound = model.conformation().clone();
    assert_eq!(
        Model::new(two_atoms, bound).unwrap_err(),
        ModelError::Conformation(ConformationError::BondCountMismatch {
            expected: 0,
            actual: 1
        })
    );
}

#[test]
fn canonical_ensemble_constructors_accept_owned_and_shared_topology() {
    let owned_new = Ensemble::new(single_atom_topology());
    let shared = owned_new.shared_topology();
    let shared_new = Ensemble::new(Arc::clone(&shared));
    assert!(Arc::ptr_eq(&shared_new.shared_topology(), &shared));

    let owned_members = Ensemble::from_items(
        single_atom_topology(),
        [EnsembleMember::new(single_position(1.0), 1.0).unwrap()],
    )
    .unwrap();
    let shared_members = Ensemble::from_items(
        owned_members.shared_topology(),
        [EnsembleMember::new(single_position(2.0), 1.0).unwrap()],
    )
    .unwrap();
    assert!(owned_members
        .topology()
        .same_layout(shared_members.topology()));
}

#[test]
fn ensemble_member_views_project_borrowed_and_owned_models_in_stable_order() {
    let topology = Arc::new(single_atom_topology());
    let score = key("score");
    let mut first = EnsembleMember::new(single_position(1.0), 1.0).unwrap();
    let mut conformation = first.conformation_mut();
    conformation.set_cell(Some(cell()));
    let mut properties = conformation.properties_mut();
    properties
        .owner_mut()
        .insert(key("method"), PropertyValue::String("sampled".into()))
        .unwrap();
    properties
        .atoms_mut()
        .set_value(score.clone(), ATOM0, Some(PropertyValue::Int(7)))
        .unwrap();
    first.set_weight(0.25).unwrap();
    let mut second = EnsembleMember::new(single_position(2.0), 1.0).unwrap();
    second.set_weight(0.75).unwrap();
    let mut ensemble = Ensemble::from_items(Arc::clone(&topology), [first, second]).unwrap();

    assert_eq!(
        ensemble
            .iter()
            .map(|member| member.positions().values().value()[0].x)
            .collect::<Vec<_>>(),
        [1.0, 2.0]
    );
    assert_eq!(
        ensemble
            .iter()
            .map(|member| member.weight())
            .collect::<Vec<_>>(),
        [0.25, 0.75]
    );

    let member = ensemble.get(0).unwrap();
    assert!(Arc::ptr_eq(&member.shared_topology(), &topology));
    let borrowed = member.as_model_view();
    assert!(Arc::ptr_eq(&borrowed.shared_topology(), &topology));
    assert_eq!(
        borrowed.positions().values().value().as_ptr(),
        member.positions().values().value().as_ptr()
    );
    let owned = member.to_model();
    assert!(Arc::ptr_eq(&owned.shared_topology(), &topology));
    assert_eq!(owned.conformation(), member.conformation());

    ensemble
        .get_mut(0)
        .unwrap()
        .conformation_mut()
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
fn detached_members_get_bond_rows_only_when_bound() {
    let (model, _, _) = model_fixture();
    let label = key("label");
    let mut member = EnsembleMember::new(model.positions().clone(), 1.0).unwrap();
    assert_eq!(member.properties().atoms().len(), 2);
    assert!(member.properties().bonds().is_empty());
    // Detached payloads cannot carry bond annotations yet.
    assert!(matches!(
        member
            .conformation_mut()
            .properties_mut()
            .bonds_mut()
            .set_value(label.clone(), BOND0, Some(PropertyValue::Int(1))),
        Err(PropertyError::InvalidIndex { len: 0, index: 0 })
    ));
    member
        .conformation_mut()
        .properties_mut()
        .owner_mut()
        .insert(label.clone(), PropertyValue::Int(7))
        .unwrap();
    let mut ensemble = Ensemble::from_items(model.shared_topology(), [member]).unwrap();
    let stored = ensemble.get(0).unwrap();
    assert_eq!(stored.properties().atoms().len(), 2);
    assert_eq!(stored.properties().bonds().len(), 1);
    assert_eq!(
        stored.properties().owner().get(&label),
        Some(&PropertyValue::Int(7))
    );

    let mut editor = ensemble.get_mut(0).unwrap();
    let before = editor.properties().clone();
    let mut conformation = editor.conformation_mut();
    assert!(conformation
        .properties_mut()
        .bonds_mut()
        .insert(label.clone(), PropertyColumn::Int(vec![Some(8); 2]))
        .is_err());
    assert_eq!(conformation.properties(), &before);
    conformation
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
fn ensemble_replacement_rejects_incompatible_members_without_changing_state() {
    let (model, _, _) = model_fixture();
    let topology = model.shared_topology();
    let mut original = EnsembleMember::new(model.positions().clone(), 1.0).unwrap();
    original.set_weight(0.75).unwrap();
    let mut ensemble = Ensemble::from_items(Arc::clone(&topology), [original]).unwrap();
    let original = ensemble.get(0).unwrap().payload().clone();
    // Bound to another one-bond topology; matching bond rows rebind freely.
    let other = model_fixture().0;
    let mut bound_elsewhere = Ensemble::from_items(
        crate::topology::Topology::from_molecules([other
            .topology()
            .molecules()
            .next()
            .unwrap()
            .molecule()
            .clone()])
        .unwrap(),
        [EnsembleMember::new(Positions::zeros(2), 1.0).unwrap()],
    )
    .unwrap()
    .into_items()
    .remove(0);
    bound_elsewhere
        .conformation_mut()
        .properties_mut()
        .bonds_mut()
        .set_value(key("bond"), BOND0, Some(PropertyValue::Int(1)))
        .unwrap();

    let replacement = EnsembleMember::new(Positions::zeros(1), 1.0).unwrap();
    let expected = RealizationError::Conformation(ConformationError::AtomCountMismatch {
        expected: 2,
        actual: 1,
    });
    assert_eq!(
        ensemble.replace(0, replacement.clone()),
        Err(expected.clone())
    );
    assert_eq!(ensemble.push(replacement), Err(expected));
    assert_eq!(ensemble.len(), 1);
    let member = ensemble.get(0).unwrap();
    assert_eq!(member.payload(), &original);
    assert_eq!(
        member.as_model_view().positions().len(),
        topology.atom_count()
    );
    // Same bond count: a payload bound elsewhere with matching rows is accepted.
    ensemble.push(bound_elsewhere).unwrap();
    assert_eq!(ensemble.len(), 2);
    ensemble.remove(1).unwrap();

    for index in [1, usize::MAX] {
        assert_eq!(
            ensemble.replace(index, original.clone()),
            Err(RealizationError::IndexOutOfRange { index, len: 1 })
        );
        assert!(ensemble.get_mut(index).is_none());
    }
    let mut empty = Ensemble::new(topology);
    assert_eq!(
        empty.replace(0, original.clone()),
        Err(RealizationError::IndexOutOfRange { index: 0, len: 0 })
    );
    assert!(empty.is_empty());
    assert_eq!(ensemble.replace(0, original.clone()).unwrap(), original);
}

#[test]
fn ensemble_replacement_preserves_order_and_publishes_complete_member_state() {
    let (model, _, _) = model_fixture();
    let topology = model.shared_topology();
    let original = EnsembleMember::new(model.positions().clone(), 1.0).unwrap();
    let mut ensemble =
        Ensemble::from_items(Arc::clone(&topology), [original.clone(), original]).unwrap();
    let original = ensemble.get(0).unwrap().payload().clone();
    let score = key("score");
    ensemble
        .properties_mut()
        .insert(score.clone(), PropertyValue::Int(9))
        .unwrap();
    let collection_properties = ensemble.properties().clone();
    let mut replacement = EnsembleMember::new(Positions::zeros(2), 1.0).unwrap();
    replacement.conformation_mut().set_cell(Some(cell()));
    replacement.set_weight(0.5).unwrap();
    replacement
        .conformation_mut()
        .properties_mut()
        .owner_mut()
        .insert(score.clone(), PropertyValue::Int(1))
        .unwrap();
    replacement
        .conformation_mut()
        .properties_mut()
        .atoms_mut()
        .set_value(score, ATOM0, Some(PropertyValue::Int(2)))
        .unwrap();

    assert_eq!(ensemble.replace(1, replacement.clone()).unwrap(), original);
    assert_eq!(ensemble.len(), 2);
    assert_eq!(ensemble.properties(), &collection_properties);
    assert_eq!(ensemble.get(0).unwrap().payload(), &original);
    let member = ensemble.get(1).unwrap();
    assert!(Arc::ptr_eq(&member.shared_topology(), &topology));
    assert_eq!(member.positions(), replacement.positions());
    assert_eq!(member.cell(), replacement.cell());
    assert_eq!(member.weight(), replacement.weight());
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
    assert_eq!(member.to_model().positions(), replacement.positions());
}

#[test]
fn ensemble_member_editor_preserves_dimensions_and_validated_values() {
    let (model, _, _) = model_fixture();
    let mut ensemble = Ensemble::from_models([model]).unwrap();
    {
        let mut member = ensemble.get_mut(0).unwrap();
        member.set_weight(0.5).unwrap();
        // Weights are finite and strictly positive.
        for weight in [f64::NAN, f64::INFINITY, -1.0, 0.0] {
            assert_eq!(
                member.set_weight(weight),
                Err(ConformationError::InvalidWeight)
            );
            assert_eq!(member.weight(), 0.5);
        }
        let properties = member.properties().clone();
        let detached = Conformation::new(Positions::zeros(2)).properties().clone();
        assert_eq!(
            member.conformation_mut().set_properties(detached),
            Err(ConformationError::BondCountMismatch {
                expected: 1,
                actual: 0
            })
        );
        let wrong_atoms = Conformation::new(Positions::zeros(1)).properties().clone();
        assert_eq!(
            member.conformation_mut().set_properties(wrong_atoms),
            Err(ConformationError::AtomCountMismatch {
                expected: 2,
                actual: 1
            })
        );
        assert_eq!(member.properties(), &properties);
        assert_eq!(
            member
                .conformation_mut()
                .set_positions(Quantity::new(vec![Point3::origin(); 3], NANOMETER)),
            Err(ConformationError::Position(
                PositionError::PositionCountMismatch {
                    expected: 2,
                    actual: 3
                }
            ))
        );
        assert_eq!(member.positions().len(), 2);
    }
    let member = ensemble.get(0).unwrap();
    assert_eq!(member.weight(), 0.5);
    assert_eq!(member.as_model_view().positions().len(), 2);
    assert_eq!(member.to_model().properties().bonds().len(), 1);
}

#[test]
fn model_view_materialization_clones_realization_and_shares_topology() {
    let mut source = Model::new(single_atom_topology(), single_position(1.0)).unwrap();
    let topology = source.shared_topology();
    let score = key("score");
    let atom = topology.atom_ids()[0];
    source
        .conformation_mut()
        .properties_mut()
        .atoms_mut()
        .set_value(score.clone(), ATOM0, Some(PropertyValue::Int(4)))
        .unwrap();
    let owned = source.as_model_view().to_model();

    assert!(Arc::ptr_eq(&owned.shared_topology(), &topology));
    assert_eq!(owned.conformation(), source.conformation());
    source
        .conformation_mut()
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
    let (mut model, _, _) = model_fixture();
    let mut conformation = model.conformation_mut();
    conformation.set_occupancy(ATOM0, Some(0.75)).unwrap();
    conformation
        .set_b_factor(ATOM0, Some(Quantity::new(12.5, SQUARE_ANGSTROM)))
        .unwrap();
    // No reserved names: a generic property may use the same key.
    conformation
        .properties_mut()
        .atoms_mut()
        .set_value(key("occupancy"), ATOM0, Some(PropertyValue::Int(1)))
        .unwrap();
    assert_eq!(model.occupancy(ATOM0).unwrap(), Some(0.75));
    assert_eq!(model.occupancies(), Some([Some(0.75), None].as_slice()));
    let b_factor = model.b_factor(ATOM0).unwrap().unwrap();
    assert_eq!(b_factor.unit(), SQUARE_NANOMETER);
    assert!((b_factor.into_value() - 0.125).abs() < 1.0e-12);
    assert_eq!(
        model
            .properties()
            .atoms()
            .value(&key("occupancy"), ATOM0)
            .unwrap(),
        Some(PropertyValue::Int(1))
    );

    let before = model.conformation().clone();
    let mut conformation = model.conformation_mut();
    assert_eq!(
        conformation.set_occupancy(ATOM0, Some(f64::NAN)),
        Err(ConformationError::NonFiniteValue { index: 0 })
    );
    assert!(matches!(
        conformation.set_b_factor(ATOM0, Some(Quantity::new(1.0, KELVIN))),
        Err(ConformationError::Unit(_))
    ));
    assert_eq!(
        conformation.set_occupancy(TopologyAtomIndex::new(2), Some(1.0)),
        Err(ConformationError::InvalidAtomIndex { index: 2, len: 2 })
    );
    assert_eq!(
        conformation.set_occupancies(Some(vec![Some(1.0)])),
        Err(ConformationError::AtomCountMismatch {
            expected: 2,
            actual: 1
        })
    );
    assert_eq!(model.conformation(), &before);

    // Clearing every value compares equal to never recording any.
    let mut cleared = model.conformation().clone();
    cleared.set_occupancy(ATOM0, None).unwrap();
    cleared.set_occupancies(None).unwrap();
    assert_eq!(cleared.occupancies(), None);
}

#[test]
fn model_subset_projects_entity_state_and_drops_owner_properties() {
    let (mut model, atom, bond) = model_fixture();
    let mut conformation = model.conformation_mut();
    conformation.set_occupancy(ATOM0, Some(0.5)).unwrap();
    let mut properties = conformation.properties_mut();
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
        crate::topology::AtomSelection::from_atoms(&model.shared_topology(), [atom]).unwrap();
    let subset = model.subset(&selection).unwrap();
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
    let all_selection = crate::topology::AtomSelection::all(&model.shared_topology());
    let fully_retained = model.subset(&all_selection).unwrap();
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
fn ensemble_collection_and_member_properties_are_separate() {
    let (mut model, atom, _) = model_fixture();
    let mut conformation = model.conformation_mut();
    let mut properties = conformation.properties_mut();
    properties
        .atoms_mut()
        .set_value(key("member_atom"), ATOM0, Some(PropertyValue::Int(1)))
        .unwrap();
    properties
        .bonds_mut()
        .set_value(key("member_bond"), BOND0, Some(PropertyValue::Int(2)))
        .unwrap();
    let mut ensemble = Ensemble::from_models([model]).unwrap();
    ensemble
        .properties_mut()
        .insert(key("collection"), PropertyValue::Bool(true))
        .unwrap();
    ensemble
        .get_mut(0)
        .unwrap()
        .conformation_mut()
        .properties_mut()
        .owner_mut()
        .insert(key("member"), PropertyValue::Int(1))
        .unwrap();
    assert!(!ensemble.properties().is_empty());
    let member = ensemble.get(0).unwrap();
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

    let cloned = ensemble.clone();
    assert_eq!(cloned.properties(), ensemble.properties());
    let selection =
        crate::topology::AtomSelection::from_atoms(&ensemble.shared_topology(), [atom]).unwrap();
    let subset = ensemble.subset(&selection).unwrap();
    assert!(subset.properties().is_empty());
    assert!(subset.get(0).unwrap().properties().owner().is_empty());
}

#[test]
fn trajectories_reinterpret_as_ensembles_without_copying_conformations() {
    let (model, _, _) = model_fixture();
    let mut frame = TrajectoryFrame::new(model.conformation().clone());
    frame.set_step(Some(4));
    frame.set_velocities(Some(Velocities::zeros(2))).unwrap();
    assert_eq!(
        frame.set_forces(Some(Forces::zeros(3))),
        Err(ConformationError::AtomCountMismatch {
            expected: 2,
            actual: 3
        })
    );
    let trajectory = Trajectory::from_items(model.shared_topology(), [frame]).unwrap();
    let pointer = trajectory
        .get(0)
        .unwrap()
        .positions()
        .values()
        .value()
        .as_ptr();
    let ensemble = trajectory.into_ensemble();
    let member = ensemble.get(0).unwrap();
    // Frames become equally weighted members.
    assert_eq!(member.weight(), 1.0);
    assert_eq!(member.positions().values().value().as_ptr(), pointer);
}

#[test]
fn model_atom_views_join_topology_identity_with_realization_state() {
    let (mut model, carbon, _) = model_fixture();
    let mut conformation = model.conformation_mut();
    conformation.set_occupancy(ATOM0, Some(0.5)).unwrap();
    conformation
        .set_b_factor(ATOM0, Some(Quantity::new(20.0, SQUARE_ANGSTROM)))
        .unwrap();
    conformation
        .properties_mut()
        .atoms_mut()
        .set_value(key("charge"), ATOM0, Some(PropertyValue::Int(-1)))
        .unwrap();

    let atom = model.atom(carbon).unwrap();
    assert_eq!(atom.id(), carbon);
    assert_eq!(atom.index(), ATOM0);
    assert_eq!(atom.element.symbol(), "C");
    assert_eq!(atom.position().into_value(), Point3::origin());
    assert_eq!(atom.occupancy(), Some(0.5));
    let b_factor = atom.b_factor().unwrap().value_in(SQUARE_ANGSTROM).unwrap();
    assert!((b_factor - 20.0).abs() < 1.0e-12);
    assert_eq!(
        atom.realization_property(&key("charge")),
        Some(PropertyValue::Int(-1))
    );
    // Static annotations stay on the topology view.
    assert_eq!(atom.property(&key("charge")), None);

    let oxygen = model.atom_at(TopologyAtomIndex::new(1)).unwrap();
    assert_eq!(oxygen.element.symbol(), "O");
    assert_eq!(oxygen.occupancy(), None);
    assert_eq!(oxygen.b_factor(), None);
    assert_eq!(
        oxygen.neighbors().map(|atom| atom.id()).collect::<Vec<_>>(),
        [carbon]
    );
    assert_eq!(
        model.atoms().map(|atom| atom.index()).collect::<Vec<_>>(),
        [ATOM0, TopologyAtomIndex::new(1)]
    );
    assert!(model.atom_at(TopologyAtomIndex::new(2)).is_none());
    let missing =
        crate::topology::InstanceAtomId::new(carbon.molecule(), crate::core::AtomId::new(9));
    assert!(model.atom(missing).is_none());
    let view = model.as_model_view();
    assert_eq!(view.atom(carbon).unwrap().position(), atom.position());
}

#[test]
fn ensemble_members_always_carry_finite_positive_relative_weights() {
    for weight in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            EnsembleMember::new(single_position(0.0), weight),
            Err(ConformationError::InvalidWeight)
        );
    }
    let topology = Arc::new(single_atom_topology());
    let members =
        [3.0, 1.0].map(|weight| EnsembleMember::new(single_position(0.0), weight).unwrap());
    let mut ensemble = Ensemble::from_items(Arc::clone(&topology), members).unwrap();

    // Weights are relative, so selection keeps them as stored.
    let resampled = ensemble.select([0, 0, 1]).unwrap();
    assert_eq!(
        resampled
            .iter()
            .map(|member| member.weight())
            .collect::<Vec<_>>(),
        [3.0, 3.0, 1.0]
    );

    ensemble.normalize_weights().unwrap();
    assert_eq!(
        ensemble
            .iter()
            .map(|member| member.weight())
            .collect::<Vec<_>>(),
        [0.75, 0.25]
    );
    // Scaling by the largest weight first keeps huge weights finite.
    let huge = [f64::MAX, f64::MAX]
        .map(|weight| EnsembleMember::new(single_position(0.0), weight).unwrap());
    let mut huge = Ensemble::from_items(Arc::clone(&topology), huge).unwrap();
    huge.normalize_weights().unwrap();
    assert_eq!(
        huge.iter()
            .map(|member| member.weight())
            .collect::<Vec<_>>(),
        [0.5, 0.5]
    );
    // A ratio that underflows is rejected without changing any weight.
    let extreme = [f64::MAX, f64::MIN_POSITIVE]
        .map(|weight| EnsembleMember::new(single_position(0.0), weight).unwrap());
    let mut extreme = Ensemble::from_items(Arc::clone(&topology), extreme).unwrap();
    assert_eq!(
        extreme.normalize_weights(),
        Err(RealizationError::Conformation(
            ConformationError::InvalidWeight
        ))
    );
    assert_eq!(extreme.get(1).unwrap().weight(), f64::MIN_POSITIVE);
    assert_eq!(
        Ensemble::new(topology).normalize_weights(),
        Err(RealizationError::EmptySource)
    );

    // Models carry no statistical weight; they become equally weighted members.
    let (model, _, _) = model_fixture();
    let ensemble = Ensemble::from_models([model.clone(), model]).unwrap();
    assert!(ensemble.iter().all(|member| member.weight() == 1.0));
}
