use super::*;
use crate::geometry::{PeriodicCell, Point3, Vector3};
use crate::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use crate::properties::{PropertyColumn, PropertyError, PropertyKey, PropertyValue};
use crate::topology::{TopologyAtomIndex, TopologyBondIndex, TopologyBuilder};
use crate::units::{
    Quantity, UnitError, ANGSTROM, KELVIN, KILOCALORIE_PER_MOLE, KILOJOULE_PER_MOLE, NANOMETER,
    SQUARE_ANGSTROM, SQUARE_NANOMETER,
};
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

/// Distinct samples, each its own state.
fn conformers(
    topology: impl Into<Arc<crate::topology::Topology>>,
    samples: impl IntoIterator<Item = EnsembleSample>,
) -> Result<Ensemble, EnsembleError> {
    let samples = samples.into_iter().collect::<Vec<_>>();
    let partition = EnsemblePartition::state_per_sample(samples.len());
    Ensemble::from_samples(topology, samples, partition)
}

/// Models drawn from one distribution: one state.
fn bundle(models: impl IntoIterator<Item = Model>) -> Result<Ensemble, EnsembleError> {
    let models = models.into_iter().collect::<Vec<_>>();
    let partition = EnsemblePartition::one_state(models.len());
    Ensemble::from_models(models, partition)
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
    let owned = Model::new(single_atom_topology(), single_position(1.0)).unwrap();
    let shared = owned.shared_topology();
    let shared_model = Model::new(Arc::clone(&shared), single_position(2.0)).unwrap();
    assert!(Arc::ptr_eq(&shared_model.shared_topology(), &shared));

    let complete = Model::new(Arc::clone(&shared), single_position(3.0))
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
    let molecule = topology.instances().next().unwrap();
    let mut builder = ModelBuilder::new();
    builder
        .add_molecule(molecule.molecule().clone(), &single_position(1.0))
        .unwrap();
    let before = builder.clone();
    assert_eq!(
        builder.set_positions(Quantity::new(vec![Point3::origin(); 2], NANOMETER)),
        Err(StructureError::Array(DenseArrayError::CountMismatch {
            expected: 1,
            actual: 2
        }))
    );
    assert_eq!(builder, before);
    assert_eq!(builder.build().unwrap().positions(), &single_position(1.0));
}

#[test]
fn structures_reject_rows_of_other_dimensions() {
    let (model, _, _) = model_fixture();
    assert_eq!(
        Model::new(model.shared_topology(), single_position(1.0)).unwrap_err(),
        ModelError::Structure(StructureError::Array(DenseArrayError::CountMismatch {
            expected: 2,
            actual: 1
        }))
    );
    // A member bound to a topology with one bond cannot bind to one with no
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
    let bound = bundle([model]).unwrap().into_samples().pop().unwrap();
    assert_eq!(bound.properties().bonds().len(), 1);
    assert_eq!(
        conformers(two_atoms, [bound]).unwrap_err(),
        EnsembleError::Structure(StructureError::BondCountMismatch {
            expected: 0,
            actual: 1
        })
    );
}

#[test]
fn collection_position_replacement_checks_counts_per_collection_transactionally() {
    let (model, _, _) = model_fixture();
    let mut ensemble = bundle([model.clone(), model.clone()]).unwrap();
    assert_eq!(
        ensemble.replace_positions([Positions::zeros(2)]),
        Err(EnsembleError::SampleCountMismatch {
            expected: 2,
            actual: 1
        })
    );
    assert_eq!(
        ensemble.replace_positions([Positions::zeros(2), Positions::zeros(3)]),
        Err(EnsembleError::Structure(StructureError::Array(
            DenseArrayError::CountMismatch {
                expected: 2,
                actual: 3
            }
        )))
    );
    assert!(ensemble
        .iter()
        .all(|member| member.positions() == model.positions()));
    ensemble
        .replace_positions([Positions::zeros(2), Positions::zeros(2)])
        .unwrap();
    assert!(ensemble
        .iter()
        .all(|member| member.positions() == &Positions::zeros(2)));

    let frame = TrajectoryFrame::new(model.positions().clone());
    let mut trajectory = Trajectory::from_items(model.shared_topology(), [frame]).unwrap();
    assert_eq!(
        trajectory.replace_positions([]),
        Err(TrajectoryError::FrameCountMismatch {
            expected: 1,
            actual: 0
        })
    );
    assert_eq!(trajectory.get(0).unwrap().positions(), model.positions());
}

#[test]
fn canonical_ensemble_constructors_accept_owned_and_shared_topology() {
    let owned_new = Ensemble::new(single_atom_topology());
    let shared = owned_new.shared_topology();
    let shared_new = Ensemble::new(Arc::clone(&shared));
    assert!(Arc::ptr_eq(&shared_new.shared_topology(), &shared));

    let owned_members = conformers(
        single_atom_topology(),
        [EnsembleSample::new(single_position(1.0), 1.0).unwrap()],
    )
    .unwrap();
    let shared_members = conformers(
        owned_members.shared_topology(),
        [EnsembleSample::new(single_position(2.0), 1.0).unwrap()],
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
    let mut first = EnsembleSample::new(single_position(1.0), 1.0).unwrap();
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
    first.set_weight(0.25).unwrap();
    let mut second = EnsembleSample::new(single_position(2.0), 1.0).unwrap();
    second.set_weight(0.75).unwrap();
    let mut ensemble = conformers(Arc::clone(&topology), [first, second]).unwrap();

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
    assert_eq!(owned.positions(), member.positions());
    assert_eq!(owned.cell(), member.cell());
    assert_eq!(owned.properties(), member.properties());

    ensemble
        .get_mut(0)
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
fn detached_members_get_bond_rows_only_when_bound() {
    let (model, _, _) = model_fixture();
    let label = key("label");
    let mut member = EnsembleSample::new(model.positions().clone(), 1.0).unwrap();
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
    let mut ensemble = conformers(model.shared_topology(), [member]).unwrap();
    let stored = ensemble.get(0).unwrap();
    assert_eq!(stored.properties().atoms().len(), 2);
    assert_eq!(stored.properties().bonds().len(), 1);
    assert_eq!(
        stored.properties().owner().get(&label),
        Some(&PropertyValue::Int(7))
    );

    let mut editor = ensemble.get_mut(0).unwrap();
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
fn ensemble_replacement_rejects_incompatible_members_without_changing_state() {
    let (model, _, _) = model_fixture();
    let topology = model.shared_topology();
    let mut original = EnsembleSample::new(model.positions().clone(), 1.0).unwrap();
    original.set_weight(0.75).unwrap();
    let mut ensemble = conformers(Arc::clone(&topology), [original]).unwrap();
    let original = ensemble.get(0).unwrap().payload().clone();
    // Bound to another one-bond topology; matching bond rows rebind freely.
    let other = model_fixture().0;
    let mut bound_elsewhere = conformers(
        crate::topology::Topology::from_molecules([other
            .topology()
            .instances()
            .next()
            .unwrap()
            .molecule()
            .clone()])
        .unwrap(),
        [EnsembleSample::new(Positions::zeros(2), 1.0).unwrap()],
    )
    .unwrap()
    .into_samples()
    .remove(0);
    bound_elsewhere
        .properties_mut()
        .bonds_mut()
        .set_value(key("bond"), BOND0, Some(PropertyValue::Int(1)))
        .unwrap();

    let replacement = EnsembleSample::new(Positions::zeros(1), 1.0).unwrap();
    let expected =
        EnsembleError::Structure(StructureError::Array(DenseArrayError::CountMismatch {
            expected: 2,
            actual: 1,
        }));
    assert_eq!(
        ensemble.replace(0, replacement.clone()),
        Err(expected.clone())
    );
    assert_eq!(
        ensemble.push_state(replacement, EnsembleState::new()),
        Err(expected)
    );
    assert_eq!(ensemble.len(), 1);
    let member = ensemble.get(0).unwrap();
    assert_eq!(member.payload(), &original);
    assert_eq!(
        member.as_model_view().positions().len(),
        topology.atom_count()
    );
    // Same bond count: a payload bound elsewhere with matching rows is accepted.
    ensemble.push_sample(bound_elsewhere, 0).unwrap();
    assert_eq!(ensemble.len(), 2);
    ensemble.remove(1).unwrap();

    for index in [1, usize::MAX] {
        assert_eq!(
            ensemble.replace(index, original.clone()),
            Err(EnsembleError::IndexOutOfRange { index, len: 1 })
        );
        assert!(ensemble.get_mut(index).is_none());
    }
    let mut empty = Ensemble::new(topology);
    assert_eq!(
        empty.replace(0, original.clone()),
        Err(EnsembleError::IndexOutOfRange { index: 0, len: 0 })
    );
    assert!(empty.is_empty());
    assert_eq!(ensemble.replace(0, original.clone()).unwrap(), original);
}

#[test]
fn ensemble_replacement_preserves_order_and_publishes_complete_member_state() {
    let (model, _, _) = model_fixture();
    let topology = model.shared_topology();
    let original = EnsembleSample::new(model.positions().clone(), 1.0).unwrap();
    let mut ensemble = conformers(Arc::clone(&topology), [original.clone(), original]).unwrap();
    let original = ensemble.get(0).unwrap().payload().clone();
    let score = key("score");
    ensemble
        .properties_mut()
        .insert(score.clone(), PropertyValue::Int(9))
        .unwrap();
    let collection_properties = ensemble.properties().clone();
    let mut replacement = EnsembleSample::new(Positions::zeros(2), 1.0).unwrap();
    replacement.set_cell(Some(cell()));
    replacement.set_weight(0.5).unwrap();
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
    let mut ensemble = bundle([model]).unwrap();
    {
        let mut member = ensemble.get_mut(0).unwrap();
        member.set_weight(0.5).unwrap();
        // Weights are finite and strictly positive.
        for weight in [f64::NAN, f64::INFINITY, -1.0, 0.0] {
            assert_eq!(
                member.set_weight(weight),
                Err(StructureError::InvalidWeight)
            );
            assert_eq!(member.weight(), 0.5);
        }
        let properties = member.properties().clone();
        let detached = EnsembleSample::new(Positions::zeros(2), 1.0)
            .unwrap()
            .properties()
            .clone();
        assert_eq!(
            member.set_properties(detached),
            Err(StructureError::BondCountMismatch {
                expected: 1,
                actual: 0
            })
        );
        let wrong_atoms = EnsembleSample::new(Positions::zeros(1), 1.0)
            .unwrap()
            .properties()
            .clone();
        assert_eq!(
            member.set_properties(wrong_atoms),
            Err(StructureError::Array(DenseArrayError::CountMismatch {
                expected: 2,
                actual: 1
            }))
        );
        assert_eq!(member.properties(), &properties);
        assert_eq!(
            member.set_positions(Quantity::new(vec![Point3::origin(); 3], NANOMETER)),
            Err(StructureError::Array(DenseArrayError::CountMismatch {
                expected: 2,
                actual: 3
            }))
        );
        assert_eq!(member.positions().len(), 2);
    }
    let member = ensemble.get(0).unwrap();
    assert_eq!(member.weight(), 0.5);
    assert_eq!(member.as_model_view().positions().len(), 2);
    assert_eq!(member.to_model().properties().bonds().len(), 1);
}

#[test]
fn model_view_materialization_clones_structure_and_shares_topology() {
    let mut source = Model::new(single_atom_topology(), single_position(1.0)).unwrap();
    let topology = source.shared_topology();
    let score = key("score");
    let atom = topology.atom_ids()[0];
    source
        .properties_mut()
        .atoms_mut()
        .set_value(score.clone(), ATOM0, Some(PropertyValue::Int(4)))
        .unwrap();
    let owned = source.as_model_view().to_model();

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
    let (mut model, _, _) = model_fixture();
    model.set_occupancy(ATOM0, Some(0.75)).unwrap();
    model
        .set_b_factor(ATOM0, Some(Quantity::new(12.5, SQUARE_ANGSTROM)))
        .unwrap();
    // No reserved names: a generic property may use the same key.
    model
        .properties_mut()
        .atoms_mut()
        .set_value(key("occupancy"), ATOM0, Some(PropertyValue::Int(1)))
        .unwrap();
    assert_eq!(model.occupancy(ATOM0), Some(0.75));
    assert_eq!(model.occupancies(), Some([Some(0.75), None].as_slice()));
    let b_factor = model.b_factor(ATOM0).unwrap();
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

    let before = model.clone();
    assert_eq!(
        model.set_occupancy(ATOM0, Some(f64::NAN)),
        Err(StructureError::Array(DenseArrayError::NonFinite {
            index: 0
        }))
    );
    assert!(matches!(
        model.set_b_factor(ATOM0, Some(Quantity::new(1.0, KELVIN))),
        Err(StructureError::Unit(_))
    ));
    assert_eq!(
        model.set_occupancy(TopologyAtomIndex::new(2), Some(1.0)),
        Err(StructureError::Array(DenseArrayError::IndexOutOfRange {
            index: 2,
            len: 2
        }))
    );
    assert_eq!(
        model.set_occupancies(Some(vec![Some(1.0)])),
        Err(StructureError::Array(DenseArrayError::CountMismatch {
            expected: 2,
            actual: 1
        }))
    );
    assert_eq!(model, before);

    // Clearing every value compares equal to never recording any.
    let mut cleared = model.clone();
    cleared.set_occupancy(ATOM0, None).unwrap();
    cleared.set_occupancies(None).unwrap();
    assert_eq!(cleared.occupancies(), None);
}

#[test]
fn model_subset_projects_entity_state_and_drops_owner_properties() {
    let (mut model, atom, bond) = model_fixture();
    model.set_occupancy(ATOM0, Some(0.5)).unwrap();
    let mut properties = model.properties_mut();
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
        crate::topology::AtomSelection::from_atoms(model.shared_topology(), [atom]).unwrap();
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
    let all_selection = crate::topology::AtomSelection::all(model.shared_topology());
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
    let mut properties = model.properties_mut();
    properties
        .atoms_mut()
        .set_value(key("member_atom"), ATOM0, Some(PropertyValue::Int(1)))
        .unwrap();
    properties
        .bonds_mut()
        .set_value(key("member_bond"), BOND0, Some(PropertyValue::Int(2)))
        .unwrap();
    let mut ensemble = bundle([model]).unwrap();
    ensemble
        .properties_mut()
        .insert(key("collection"), PropertyValue::Bool(true))
        .unwrap();
    ensemble
        .get_mut(0)
        .unwrap()
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
        crate::topology::AtomSelection::from_atoms(ensemble.shared_topology(), [atom]).unwrap();
    let subset = ensemble.subset(&selection).unwrap();
    assert!(subset.properties().is_empty());
    assert!(subset.get(0).unwrap().properties().owner().is_empty());
}

#[test]
fn trajectories_reinterpret_as_ensembles_without_copying_positions() {
    let (model, _, _) = model_fixture();
    let mut frame = TrajectoryFrame::new(model.positions().clone());
    frame.set_step(Some(4));
    frame.set_velocities(Some(Velocities::zeros(2))).unwrap();
    assert_eq!(
        frame.set_forces(Some(Forces::zeros(3))),
        Err(StructureError::Array(DenseArrayError::CountMismatch {
            expected: 2,
            actual: 3
        }))
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
fn trajectories_sample_chosen_frames_into_ensembles_without_consuming() {
    let (model, carbon, _) = model_fixture();
    let frames = (0..5).map(|index| {
        let mut frame = TrajectoryFrame::new(
            Positions::new(Quantity::new(
                [Point3::new(index as f64, 0.0, 0.0), Point3::origin()],
                NANOMETER,
            ))
            .unwrap(),
        );
        frame.set_velocities(Some(Velocities::zeros(2))).unwrap();
        frame
    });
    let trajectory = Trajectory::from_items(model.shared_topology(), frames).unwrap();

    // Leave out equilibration, take every second frame, and weight them.
    let sampling = TrajectorySampling::new((1..5).step_by(2))
        .with_weights([3.0, 1.0])
        .unwrap();
    let ensemble = trajectory.to_ensemble(&sampling).unwrap();
    assert_eq!(trajectory.len(), 5);
    assert_eq!(ensemble.state_count(), 1);
    assert!(ensemble.state(0).unwrap().payload().energy().is_none());
    let samples = ensemble
        .iter()
        .map(|sample| (sample.positions().values().value()[0].x, sample.weight()))
        .collect::<Vec<_>>();
    assert_eq!(samples, [(1.0, 3.0), (3.0, 1.0)]);

    let atoms =
        crate::topology::AtomSelection::from_atoms(model.shared_topology(), [carbon]).unwrap();
    let subset = trajectory
        .to_ensemble(&TrajectorySampling::new([4]).with_atoms(atoms))
        .unwrap();
    assert_eq!(subset.topology().atom_count(), 1);
    assert_eq!(subset.get(0).unwrap().positions().len(), 1);

    assert_eq!(
        trajectory
            .to_ensemble(&TrajectorySampling::new([5]))
            .unwrap_err(),
        TrajectoryError::IndexOutOfRange { index: 5, len: 5 }
    );
}

#[test]
fn occupancies_and_b_factors_follow_models_and_members_but_never_frames() {
    let (mut model, _, _) = model_fixture();
    model.set_occupancy(ATOM0, Some(0.5)).unwrap();
    model
        .set_b_factor(ATOM0, Some(Quantity::new(20.0, SQUARE_ANGSTROM)))
        .unwrap();
    model.set_cell(Some(cell()));

    // Models become members with their occupancies, B-factors, cell, and rows.
    let ensemble = bundle([model.clone()]).unwrap();
    let member = ensemble.get(0).unwrap();
    assert_eq!(member.occupancies(), Some([Some(0.5), None].as_slice()));
    assert_eq!(
        member.b_factors().unwrap().value(),
        model.b_factors().unwrap().value()
    );
    assert_eq!(member.cell(), Some(&cell()));
    assert_eq!(member.as_model_view().occupancy(ATOM0), Some(0.5));
    assert_eq!(member.to_model(), model);

    // Frame views report no occupancies or B-factors but still check indices.
    let mut frame = TrajectoryFrame::new(model.positions().clone()).with_cell(cell());
    frame
        .properties_mut()
        .owner_mut()
        .insert(key("source"), PropertyValue::Int(1))
        .unwrap();
    let trajectory = Trajectory::from_items(model.shared_topology(), [frame]).unwrap();
    let view = trajectory.get(0).unwrap().as_model_view();
    assert_eq!(view.occupancies(), None);
    assert!(view.b_factors().is_none());
    assert_eq!(view.occupancy(ATOM0), None);
    assert_eq!(view.b_factor(ATOM0), None);
    assert_eq!(view.occupancy(TopologyAtomIndex::new(2)), None);
    assert_eq!(view.to_model().occupancies(), None);

    // Frames reinterpreted as members keep their cell and rows only.
    let ensemble = trajectory.into_ensemble();
    let member = ensemble.get(0).unwrap();
    assert_eq!(member.occupancies(), None);
    assert_eq!(member.cell(), Some(&cell()));
    assert_eq!(
        member.properties().owner().get(&key("source")),
        Some(&PropertyValue::Int(1))
    );
}

#[test]
fn model_atom_views_join_topology_identity_with_structure_state() {
    let (mut model, carbon, _) = model_fixture();
    model.set_occupancy(ATOM0, Some(0.5)).unwrap();
    model
        .set_b_factor(ATOM0, Some(Quantity::new(20.0, SQUARE_ANGSTROM)))
        .unwrap();
    model
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
        atom.structure_property(&key("charge")),
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
        crate::topology::InstanceAtomId::new(carbon.instance(), crate::molecule::AtomId::new(9));
    assert!(model.atom(missing).is_none());
    let view = model.as_model_view();
    assert_eq!(view.atom(carbon).unwrap().position(), atom.position());
}

#[test]
fn ensemble_members_always_carry_finite_positive_relative_weights() {
    for weight in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            EnsembleSample::new(single_position(0.0), weight),
            Err(StructureError::InvalidWeight)
        );
    }
    let topology = Arc::new(single_atom_topology());
    let members =
        [3.0, 1.0].map(|weight| EnsembleSample::new(single_position(0.0), weight).unwrap());
    let mut ensemble = conformers(Arc::clone(&topology), members).unwrap();

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
        .map(|weight| EnsembleSample::new(single_position(0.0), weight).unwrap());
    let mut huge = conformers(Arc::clone(&topology), huge).unwrap();
    huge.normalize_weights().unwrap();
    assert_eq!(
        huge.iter()
            .map(|member| member.weight())
            .collect::<Vec<_>>(),
        [0.5, 0.5]
    );
    // A ratio that underflows is rejected without changing any weight.
    let extreme = [f64::MAX, f64::MIN_POSITIVE]
        .map(|weight| EnsembleSample::new(single_position(0.0), weight).unwrap());
    let mut extreme = conformers(Arc::clone(&topology), extreme).unwrap();
    assert_eq!(
        extreme.normalize_weights(),
        Err(EnsembleError::Structure(StructureError::InvalidWeight))
    );
    assert_eq!(extreme.get(1).unwrap().weight(), f64::MIN_POSITIVE);
    assert_eq!(
        Ensemble::new(topology).normalize_weights(),
        Err(EnsembleError::EmptySource)
    );

    // Models carry no statistical weight; they become equally weighted members.
    let (model, _, _) = model_fixture();
    let ensemble = bundle([model.clone(), model]).unwrap();
    assert!(ensemble.iter().all(|member| member.weight() == 1.0));
}

/// Unit-weight conformers of a one-atom topology, one state each, with
/// state energies in kJ/mol.
fn energy_ensemble(energies: &[f64]) -> Ensemble {
    let mut ensemble = Ensemble::new(single_atom_topology());
    for energy in energies {
        let state = EnsembleState::new()
            .with_energy(Quantity::new(*energy, KILOJOULE_PER_MOLE))
            .unwrap();
        ensemble
            .push_state(
                EnsembleSample::new(single_position(0.0), 1.0).unwrap(),
                state,
            )
            .unwrap();
    }
    ensemble
}

/// Samples of a one-atom topology with these weights, grouped by
/// `assignments` with `representatives`.
fn grouped(weights: &[f64], assignments: &[usize], representatives: &[usize]) -> Ensemble {
    let samples = weights
        .iter()
        .map(|weight| EnsembleSample::new(single_position(0.0), *weight).unwrap());
    let partition = EnsemblePartition::new(assignments.to_vec(), representatives.to_vec()).unwrap();
    Ensemble::from_samples(single_atom_topology(), samples, partition).unwrap()
}

fn kj(value: f64) -> Option<Quantity<f64>> {
    Some(Quantity::new(value, KILOJOULE_PER_MOLE))
}

/// Sets the energy of every state, in kJ/mol where given.
fn set_energies(ensemble: &mut Ensemble, energies: &[Option<f64>]) {
    for (state, energy) in energies.iter().enumerate() {
        ensemble
            .state_mut(state)
            .unwrap()
            .set_energy(energy.and_then(kj))
            .unwrap();
    }
}

fn weights(ensemble: &Ensemble) -> Vec<f64> {
    ensemble.iter().map(|sample| sample.weight()).collect()
}

fn energies(ensemble: &Ensemble) -> Vec<Option<Quantity<f64>>> {
    ensemble.states().map(|state| state.energy()).collect()
}

fn representatives(ensemble: &Ensemble) -> Vec<usize> {
    ensemble.partition().representatives().to_vec()
}

#[test]
fn state_energies_are_canonical_finite_and_never_on_samples() {
    let state = EnsembleState::new();
    assert_eq!(state.energy(), None);
    let mut state = state
        .with_energy(Quantity::new(1.0, KILOCALORIE_PER_MOLE))
        .unwrap();
    assert_eq!(state.energy(), kj(4.184));

    // Invalid energies are rejected without changing the stored one.
    for value in [f64::NAN, f64::INFINITY] {
        assert_eq!(
            state.set_energy(Some(Quantity::new(value, KILOJOULE_PER_MOLE))),
            Err(StructureError::NonFiniteEnergy)
        );
    }
    assert_eq!(
        state.set_energy(Some(Quantity::new(1.0, KELVIN))),
        Err(StructureError::Unit(Box::new(
            UnitError::IncompatibleUnits {
                from: KELVIN,
                to: KILOJOULE_PER_MOLE,
            }
        )))
    );
    assert_eq!(state.energy(), kj(4.184));

    let mut ensemble = Ensemble::new(single_atom_topology());
    let index = ensemble
        .push_state(
            EnsembleSample::new(single_position(0.0), 2.0).unwrap(),
            state,
        )
        .unwrap();
    assert_eq!(index, 0);
    ensemble.state_mut(0).unwrap().set_energy(kj(-3.0)).unwrap();
    // Energies describe the sampled structure, so subsets keep them.
    let all = crate::topology::AtomSelection::all(ensemble.shared_topology());
    assert_eq!(energies(&ensemble.subset(&all).unwrap()), [kj(-3.0)]);
    assert_eq!(ensemble.get(0).unwrap().state().energy(), kj(-3.0));
    ensemble.state_mut(0).unwrap().set_energy(None).unwrap();
    assert_eq!(energies(&ensemble), [None]);
    assert!(ensemble.state_mut(1).is_none());

    // Models carry no energy.
    let (model, _, _) = model_fixture();
    assert_eq!(energies(&bundle([model]).unwrap()), [None]);
}

#[test]
fn every_sample_is_in_one_state_that_contains_its_representative() {
    let mut ensemble = Ensemble::new(single_atom_topology());
    assert_eq!(ensemble.state_count(), 0);
    let sample = |x| EnsembleSample::new(single_position(x), 1.0).unwrap();
    assert_eq!(
        ensemble.push_sample(sample(0.0), 0),
        Err(EnsembleError::StateOutOfRange {
            state: 0,
            states: 0
        })
    );
    ensemble
        .push_state(sample(0.0), EnsembleState::new())
        .unwrap();
    ensemble
        .push_state(sample(1.0), EnsembleState::new())
        .unwrap();
    ensemble.push_sample(sample(2.0), 0).unwrap();
    assert_eq!(ensemble.partition().assignments(), [0, 1, 0]);
    assert_eq!(representatives(&ensemble), [0, 1]);
    let state = ensemble.state(0).unwrap();
    assert_eq!(state.samples().collect::<Vec<_>>(), [0, 2]);
    assert_eq!(state.representative_sample().index(), 0);
    assert!(ensemble.get(0).unwrap().is_representative());
    assert!(!ensemble.get(2).unwrap().is_representative());
    assert_eq!(ensemble.get(2).unwrap().state().index(), 0);

    // A representative of a state with other samples stays.
    assert_eq!(
        ensemble.remove(0),
        Err(EnsembleError::RepresentativeRemoval {
            sample: 0,
            state: 0
        })
    );
    // Removing a state's only sample removes the state; later indices shift.
    ensemble.remove(1).unwrap();
    assert_eq!(ensemble.partition().assignments(), [0, 0]);
    assert_eq!(representatives(&ensemble), [0]);
    ensemble.remove(1).unwrap();
    ensemble.remove(0).unwrap();
    assert!(ensemble.is_empty());
    assert_eq!(ensemble.state_count(), 0);

    // Partitions cover the ensemble's samples.
    let samples = [sample(0.0), sample(1.0)];
    assert_eq!(
        Ensemble::from_samples(
            single_atom_topology(),
            samples,
            EnsemblePartition::one_state(3)
        )
        .unwrap_err(),
        EnsembleError::PartitionMismatch {
            expected: 2,
            actual: 3
        }
    );
    assert_eq!(
        EnsemblePartition::new(vec![0, 1], vec![1, 0]),
        Err(PartitionError::RepresentativeNotInState {
            state: 0,
            sample: 1
        })
    );
    assert_eq!(
        EnsemblePartition::new(vec![0, 2], vec![0, 1]),
        Err(PartitionError::StateOutOfRange {
            sample: 1,
            state: 2,
            states: 2
        })
    );
    assert_eq!(
        EnsemblePartition::new(vec![0], vec![0, 3]),
        Err(PartitionError::RepresentativeOutOfRange {
            state: 1,
            sample: 3,
            samples: 1
        })
    );
    assert_eq!(EnsemblePartition::one_state(0).state_count(), 0);
    assert_eq!(EnsemblePartition::one_state(2).representatives(), [0]);
}

#[test]
fn populations_are_derived_from_sample_weights() {
    // Two states: samples 0 and 2 (weights 1 and 3) and sample 1 (weight 4).
    let mut ensemble = grouped(&[1.0, 4.0, 3.0], &[0, 1, 0], &[2, 1]);
    assert_eq!(ensemble.state(0).unwrap().population(), 4.0);
    assert_eq!(ensemble.state(1).unwrap().population(), 4.0);
    assert_eq!(ensemble.probabilities().unwrap(), [0.125, 0.5, 0.375]);
    assert_eq!(ensemble.state_probabilities().unwrap(), [0.5, 0.5]);
    // Kish: (1 + 4 + 3)² / (1 + 16 + 9) = 64 / 26.
    assert_eq!(ensemble.effective_sample_size().unwrap(), 64.0 / 26.0);

    ensemble.get_mut(1).unwrap().set_weight(12.0).unwrap();
    assert_eq!(ensemble.state(1).unwrap().population(), 12.0);
    assert_eq!(ensemble.state_probabilities().unwrap(), [0.25, 0.75]);
    assert_eq!(ensemble.ranked_by_population(), [1, 0]);
    ensemble.normalize_weights().unwrap();
    assert_eq!(weights(&ensemble), [0.0625, 0.75, 0.1875]);
    assert_eq!(ensemble.state_probabilities().unwrap(), [0.25, 0.75]);

    let empty = Ensemble::new(single_atom_topology());
    assert_eq!(empty.state_probabilities(), Err(EnsembleError::EmptySource));
}

#[test]
fn boltzmann_populations_follow_relative_state_energies_at_temperature() {
    // 1, 0, and 2 kcal/mol above the minimum at 298.15 K. Expected factors
    // exp(-dE/RT) evaluated independently to 40 digits (Python `decimal`,
    // R = N_A k_B exact in the 2019 SI).
    let mut ensemble = energy_ensemble(&[-100.0 + 4.184, -100.0, -100.0 + 8.368]);
    ensemble
        .assign_boltzmann_populations(Quantity::new(298.15, KELVIN))
        .unwrap();
    let expected = [0.184_924_695_965_256_4, 1.0, 0.034_197_143_177_842_51];
    let weights = weights(&ensemble);
    for (weight, expected) in weights.iter().zip(expected) {
        let ratio = weight / weights[1];
        assert!((ratio - expected).abs() <= 1.0e-15, "{ratio} != {expected}");
    }
    // The total weight is kept.
    assert!((weights.iter().sum::<f64>() - 3.0).abs() <= 1.0e-15);
    assert_eq!(ensemble.temperature(), Some(Quantity::new(298.15, KELVIN)));

    // Invalid temperatures and missing energies change nothing.
    let before = self::weights(&ensemble);
    for temperature in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            ensemble.assign_boltzmann_populations(Quantity::new(temperature, KELVIN)),
            Err(EnsembleError::InvalidTemperature)
        );
    }
    assert_eq!(
        ensemble.assign_boltzmann_populations(Quantity::new(1.0, KILOJOULE_PER_MOLE)),
        Err(EnsembleError::Unit(Box::new(
            UnitError::IncompatibleUnits {
                from: KILOJOULE_PER_MOLE,
                to: KELVIN,
            }
        )))
    );
    ensemble.state_mut(2).unwrap().set_energy(None).unwrap();
    assert_eq!(
        ensemble.assign_boltzmann_populations(Quantity::new(310.0, KELVIN)),
        Err(EnsembleError::MissingEnergy { state: 2 })
    );
    assert_eq!(self::weights(&ensemble), before);
    assert_eq!(ensemble.temperature(), Some(Quantity::new(298.15, KELVIN)));

    // A factor that underflows to zero is rejected transactionally.
    let mut wide = energy_ensemble(&[0.0, 1.0e6]);
    assert_eq!(
        wide.assign_boltzmann_populations(Quantity::new(298.15, KELVIN)),
        Err(EnsembleError::Structure(StructureError::InvalidWeight))
    );
    assert_eq!(self::weights(&wide), [1.0, 1.0]);
    assert_eq!(wide.temperature(), None);
}

#[test]
fn boltzmann_populations_scale_the_samples_of_each_state_together() {
    // State 0 (0 kJ/mol) holds samples of weights 1 and 3; state 1, RT ln 3
    // higher, one sample of weight 2. Populations go from 4 : 2 to 3 : 1 of
    // the kept total 6, so 4.5 and 1.5: state 0's samples scale by 4.5 / 4
    // and state 1's by 1.5 / 2.
    let temperature = 300.0;
    let gap = 8.314_462_618_153_24e-3 * temperature * 3.0_f64.ln();
    let mut ensemble = grouped(&[1.0, 2.0, 3.0], &[0, 1, 0], &[0, 1]);
    set_energies(&mut ensemble, &[Some(0.0), Some(gap)]);
    ensemble
        .assign_boltzmann_populations(Quantity::new(temperature, KELVIN))
        .unwrap();
    for (weight, expected) in weights(&ensemble).into_iter().zip([1.125, 1.5, 3.375]) {
        assert!(
            (weight - expected).abs() <= 1.0e-14,
            "{weight} != {expected}"
        );
    }
    assert_eq!(ensemble.partition().assignments(), [0, 1, 0]);

    // States whose weights differ by more than the floating-point range
    // still get representable populations.
    let mut spread = grouped(&[1.0e-200, 1.0e200], &[0, 1], &[0, 1]);
    set_energies(&mut spread, &[Some(0.0), Some(0.0)]);
    spread
        .assign_boltzmann_populations(Quantity::new(temperature, KELVIN))
        .unwrap();
    for weight in weights(&spread) {
        assert!((weight / 5.0e199 - 1.0).abs() <= 1.0e-12, "{weight}");
    }

    // Weighting one state changes nothing: its population is all there is.
    let mut one = grouped(&[1.0, 2.0, 5.0], &[0, 0, 0], &[1]);
    set_energies(&mut one, &[Some(7.0)]);
    one.assign_boltzmann_populations(Quantity::new(298.15, KELVIN))
        .unwrap();
    assert_eq!(weights(&one), [1.0, 2.0, 5.0]);
    assert_eq!(one.temperature(), Some(Quantity::new(298.15, KELVIN)));

    // Samples drawn from a simulation form a state without an energy, which
    // is never weighted a second time.
    let (model, _, _) = model_fixture();
    let mut drawn = bundle([model.clone(), model]).unwrap();
    assert_eq!(
        drawn.assign_boltzmann_populations(Quantity::new(298.15, KELVIN)),
        Err(EnsembleError::MissingEnergy { state: 0 })
    );
    assert_eq!(drawn.temperature(), None);

    // A temperature is recorded or cleared explicitly, and validated.
    drawn
        .set_temperature(Some(Quantity::new(298.15, KELVIN)))
        .unwrap();
    assert_eq!(drawn.temperature(), Some(Quantity::new(298.15, KELVIN)));
    assert_eq!(
        drawn.set_temperature(Some(Quantity::new(-1.0, KELVIN))),
        Err(EnsembleError::InvalidTemperature)
    );
    drawn.set_temperature(None).unwrap();
    assert_eq!(drawn.temperature(), None);
}

/// Whether every sample of `a` shares its structure with the same sample
/// of `b`.
fn shares_every_structure(a: &Ensemble, b: &Ensemble) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b)
            .all(|(a, b)| a.payload().shares_structure_with(b.payload()))
}

fn x(ensemble: &Ensemble, sample: usize) -> f64 {
    ensemble.get(sample).unwrap().positions().values().value()[0].x
}

#[test]
fn regrouping_and_reweighting_a_shared_ensemble_copy_no_structure() {
    let shared = Arc::new(energy_ensemble(&[0.0, 1.0, 2.0, 3.0]));
    let mut edited = Arc::clone(&shared);

    let ensemble = Arc::make_mut(&mut edited);
    ensemble
        .regroup(EnsemblePartition::new(vec![0, 0, 1, 1], vec![0, 2]).unwrap())
        .unwrap();
    ensemble.get_mut(1).unwrap().set_weight(3.0).unwrap();
    ensemble.normalize_weights().unwrap();
    ensemble.state_mut(0).unwrap().set_energy(kj(0.0)).unwrap();
    ensemble.state_mut(1).unwrap().set_energy(kj(1.0)).unwrap();
    ensemble
        .assign_boltzmann_populations(Quantity::new(300.0, KELVIN))
        .unwrap();

    assert!(shares_every_structure(&edited, &shared));
    // The held ensemble keeps its own grouping and weights.
    assert_eq!(shared.state_count(), 4);
    assert_eq!(weights(&shared), [1.0; 4]);
    assert_eq!(shared.temperature(), None);
    assert_eq!(edited.state_count(), 2);
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

    // A sample taken from an ensemble of the same topology stays shared.
    let mut other = Ensemble::new(shared.shared_topology());
    let sample = shared.get(2).unwrap().payload().clone();
    other.push_state(sample, EnsembleState::new()).unwrap();
    assert!(other
        .get(0)
        .unwrap()
        .payload()
        .shares_structure_with(shared.get(2).unwrap().payload()));
}

#[test]
fn editing_one_sample_of_a_shared_ensemble_copies_that_sample_alone() {
    let original = energy_ensemble(&[0.0, 1.0, 2.0]);
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
    for (sample, shared) in [(0, true), (1, false), (2, false)] {
        assert_eq!(
            copy.get(sample)
                .unwrap()
                .payload()
                .shares_structure_with(original.get(sample).unwrap().payload()),
            shared,
            "sample {sample}"
        );
    }

    // Transforms and position replacement leave the original unchanged.
    let mut transformed = original.clone();
    use crate::geometry::{Matrix3, RigidTransform};
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
            .map(|sample| x(&original, sample))
            .collect::<Vec<_>>(),
        [0.0; 3]
    );

    // An edit that fails after copying changes neither ensemble.
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
    let mut ensemble = energy_ensemble(&[0.0, 1.0, 2.0, 3.0]);
    for (index, weight) in [1.0, 2.0, 3.0, 4.0].into_iter().enumerate() {
        ensemble.get_mut(index).unwrap().set_weight(weight).unwrap();
    }
    ensemble
        .state_mut(1)
        .unwrap()
        .properties_mut()
        .insert(key("label"), PropertyValue::Int(1))
        .unwrap();
    ensemble
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    let before = ensemble.clone();

    // An unchanged grouping keeps everything.
    ensemble.regroup(ensemble.partition().clone()).unwrap();
    assert_eq!(energies(&ensemble), energies(&before));
    assert_eq!(
        ensemble.state(1).unwrap().properties(),
        before.state(1).unwrap().properties()
    );

    // Merging conformers 0 and 1 under 1 keeps conformer 1's state.
    let merged = EnsemblePartition::new(vec![0, 0, 1, 2], vec![1, 2, 3]).unwrap();
    ensemble.regroup(merged).unwrap();
    assert_eq!(weights(&ensemble), [1.0, 2.0, 3.0, 4.0]);
    assert_eq!(energies(&ensemble), [kj(1.0), kj(2.0), kj(3.0)]);
    assert_eq!(
        ensemble.state(0).unwrap().properties().get(&key("label")),
        Some(&PropertyValue::Int(1))
    );
    assert_eq!(ensemble.state(0).unwrap().population(), 3.0);
    assert_eq!(ensemble.temperature(), Some(Quantity::new(300.0, KELVIN)));
    // Samples never move.
    assert_eq!(
        ensemble
            .iter()
            .map(|sample| sample.positions().clone())
            .collect::<Vec<_>>(),
        before
            .iter()
            .map(|sample| sample.positions().clone())
            .collect::<Vec<_>>()
    );

    // Splitting the merged state loses its energy and properties; the state
    // that takes a whole old state with part of a split one keeps the whole
    // one's energy when its representative comes from it.
    let split = EnsemblePartition::new(vec![0, 1, 0, 2], vec![2, 1, 3]).unwrap();
    ensemble.regroup(split).unwrap();
    assert_eq!(weights(&ensemble), [1.0, 2.0, 3.0, 4.0]);
    assert_eq!(energies(&ensemble), [kj(2.0), None, kj(3.0)]);
    assert!(ensemble.state(1).unwrap().properties().is_empty());

    // A partition of another sample count changes nothing.
    let snapshot = ensemble.clone();
    assert_eq!(
        ensemble.regroup(EnsemblePartition::one_state(3)),
        Err(EnsembleError::PartitionMismatch {
            expected: 4,
            actual: 3
        })
    );
    assert_eq!(ensemble.partition(), snapshot.partition());
    assert_eq!(energies(&ensemble), energies(&snapshot));

    // One state split into several loses its energy.
    let mut one = grouped(&[1.0, 1.0, 1.0, 1.0], &[0, 0, 0, 0], &[0]);
    set_energies(&mut one, &[Some(5.0)]);
    one.regroup(EnsemblePartition::new(vec![0, 0, 1, 1], vec![0, 2]).unwrap())
        .unwrap();
    assert_eq!(energies(&one), [None, None]);
    assert_eq!(weights(&one), [1.0; 4]);
}

#[test]
fn collapse_keeps_one_representative_per_state_with_its_population() {
    let mut ensemble = grouped(&[1.0, 2.0, 3.0, 4.0], &[0, 1, 0, 1], &[2, 1]);
    ensemble
        .replace_positions([1.0, 2.0, 3.0, 4.0].map(single_position))
        .unwrap();
    set_energies(&mut ensemble, &[Some(0.0), None]);
    ensemble
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    let source = ensemble.clone();

    let collapsed = ensemble.collapse().unwrap();
    assert_eq!(collapsed.len(), 2);
    assert_eq!(weights(&collapsed), [4.0, 6.0]);
    assert_eq!(
        collapsed.get(0).unwrap().positions(),
        ensemble.get(2).unwrap().positions()
    );
    assert_eq!(
        collapsed.get(1).unwrap().positions(),
        ensemble.get(1).unwrap().positions()
    );
    assert_eq!(
        collapsed.partition(),
        &EnsemblePartition::state_per_sample(2)
    );
    assert_eq!(energies(&collapsed), [kj(0.0), None]);
    assert_eq!(collapsed.temperature(), ensemble.temperature());
    for (collapsed, source) in collapsed
        .state_probabilities()
        .unwrap()
        .into_iter()
        .zip(ensemble.state_probabilities().unwrap())
    {
        assert!((collapsed - source).abs() <= 1.0e-15);
    }
    // The source is unchanged.
    assert_eq!(weights(&ensemble), weights(&source));
    assert_eq!(ensemble.partition(), source.partition());
    assert_eq!(
        Ensemble::new(single_atom_topology())
            .collapse()
            .unwrap()
            .len(),
        0
    );
}

#[test]
fn selection_copies_samples_and_whole_or_partial_states() {
    // State 0: samples 0 and 1, represented by 1; state 1: sample 2.
    let mut ensemble = grouped(&[1.0, 2.0, 3.0], &[0, 0, 1], &[1, 2]);
    set_energies(&mut ensemble, &[Some(1.0), Some(2.0)]);

    // Weights are relative, so selection keeps them as stored.
    let whole = ensemble.select([2, 0, 1, 1]).unwrap();
    assert_eq!(weights(&whole), [3.0, 1.0, 2.0, 2.0]);
    assert_eq!(whole.partition().assignments(), [1, 0, 0, 0]);
    assert_eq!(representatives(&whole), [2, 0]);
    assert_eq!(energies(&whole), [kj(1.0), kj(2.0)]);

    // A partly selected state has no energy.
    let partial = ensemble.select([0, 2]).unwrap();
    assert_eq!(partial.partition().assignments(), [0, 1]);
    assert_eq!(representatives(&partial), [0, 1]);
    assert_eq!(energies(&partial), [None, kj(2.0)]);

    // Whole states keep their samples, representative, and energy.
    let states = ensemble.select_states([1, 0]).unwrap();
    assert_eq!(weights(&states), [3.0, 1.0, 2.0]);
    assert_eq!(states.partition().assignments(), [0, 1, 1]);
    assert_eq!(representatives(&states), [0, 2]);
    assert_eq!(energies(&states), [kj(2.0), kj(1.0)]);
    assert_eq!(
        ensemble.select_states([2]).unwrap_err(),
        EnsembleError::StateOutOfRange {
            state: 2,
            states: 2
        }
    );
}

#[test]
fn population_ranking_and_coverage_select_the_most_populated_states() {
    let topology = Arc::new(single_atom_topology());
    let ensemble = conformers(
        Arc::clone(&topology),
        [1.0, 4.0, 2.0, 2.0, 1.0]
            .map(|weight| EnsembleSample::new(single_position(0.0), weight).unwrap()),
    )
    .unwrap();

    // Equal weights and populations keep sample and state order.
    assert_eq!(ensemble.ranked_by_weight(), [1, 2, 3, 0, 4]);
    assert_eq!(ensemble.ranked_by_population(), [1, 2, 3, 0, 4]);
    assert_eq!(ensemble.probabilities().unwrap(), [0.1, 0.4, 0.2, 0.2, 0.1]);
    assert_eq!(weights(&ensemble), [1.0, 4.0, 2.0, 2.0, 1.0]);

    assert_eq!(ensemble.covering_population(0.4).unwrap(), [1]);
    assert_eq!(ensemble.covering_population(0.5).unwrap(), [1, 2]);
    assert_eq!(ensemble.covering_population(0.85).unwrap(), [1, 2, 3, 0]);
    assert_eq!(ensemble.covering_population(1.0).unwrap(), [1, 2, 3, 0, 4]);
    for fraction in [0.0, -0.5, 1.5, f64::NAN] {
        assert_eq!(
            ensemble.covering_population(fraction),
            Err(EnsembleError::InvalidFraction)
        );
    }
    // States of several samples rank by their summed weight.
    let pooled = grouped(&[3.0, 2.0, 2.0], &[0, 1, 1], &[0, 1]);
    assert_eq!(pooled.ranked_by_weight(), [0, 1, 2]);
    assert_eq!(pooled.ranked_by_population(), [1, 0]);
    assert_eq!(pooled.covering_population(0.5).unwrap(), [1]);

    // A full fraction keeps states too light to change the total.
    let light = conformers(
        Arc::clone(&topology),
        [1.0, 1.0e-20].map(|weight| EnsembleSample::new(single_position(0.0), weight).unwrap()),
    )
    .unwrap();
    assert_eq!(light.covering_population(1.0).unwrap(), [0, 1]);
    assert_eq!(light.covering_population(0.999).unwrap(), [0]);

    let skewed = conformers(
        Arc::clone(&topology),
        [2.0, 1.0, 1.0].map(|weight| EnsembleSample::new(single_position(0.0), weight).unwrap()),
    )
    .unwrap();
    assert_eq!(skewed.effective_sample_size().unwrap(), 8.0 / 3.0);
    let uniform = conformers(
        Arc::clone(&topology),
        [7.0; 4].map(|weight| EnsembleSample::new(single_position(0.0), weight).unwrap()),
    )
    .unwrap();
    assert_eq!(uniform.effective_sample_size().unwrap(), 4.0);

    let empty = Ensemble::new(topology);
    assert_eq!(empty.ranked_by_weight(), Vec::<usize>::new());
    assert_eq!(empty.ranked_by_population(), Vec::<usize>::new());
    assert_eq!(empty.probabilities(), Err(EnsembleError::EmptySource));
    assert_eq!(
        empty.effective_sample_size(),
        Err(EnsembleError::EmptySource)
    );
    assert_eq!(
        empty.covering_population(0.5),
        Err(EnsembleError::EmptySource)
    );
}

#[test]
fn energy_ranking_and_windows_are_relative_to_the_lowest_state() {
    let mut ensemble = energy_ensemble(&[10.0, 0.0, 30.0, 4.184, 0.0]);
    assert_eq!(ensemble.ranked_by_energy().unwrap(), [1, 4, 3, 0, 2]);
    // The window bound is inclusive and converted to the canonical unit.
    assert_eq!(
        ensemble
            .within_energy_window(Quantity::new(1.0, KILOCALORIE_PER_MOLE))
            .unwrap(),
        [1, 3, 4]
    );
    assert_eq!(
        ensemble
            .within_energy_window(Quantity::new(0.0, KILOJOULE_PER_MOLE))
            .unwrap(),
        [1, 4]
    );
    let window = ensemble
        .within_energy_window(Quantity::new(10.0, KILOJOULE_PER_MOLE))
        .unwrap();
    assert_eq!(window, [0, 1, 3, 4]);
    assert_eq!(
        energies(&ensemble.select_states(window).unwrap()),
        [10.0, 0.0, 4.184, 0.0].map(kj)
    );

    for value in [-1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            ensemble.within_energy_window(Quantity::new(value, KILOJOULE_PER_MOLE)),
            Err(EnsembleError::InvalidEnergyWindow)
        );
    }
    ensemble.state_mut(3).unwrap().set_energy(None).unwrap();
    assert_eq!(
        ensemble.ranked_by_energy(),
        Err(EnsembleError::MissingEnergy { state: 3 })
    );
    assert_eq!(
        ensemble.within_energy_window(Quantity::new(1.0, KILOJOULE_PER_MOLE)),
        Err(EnsembleError::MissingEnergy { state: 3 })
    );
}

#[test]
fn model_views_report_which_structure_they_borrow() {
    let topology = Arc::new(single_atom_topology());
    let model = Model::new(Arc::clone(&topology), single_position(0.0)).unwrap();
    assert_eq!(model.as_model_view().kind(), StructureKind::Model);

    let mut ensemble = energy_ensemble(&[2.5]);
    ensemble
        .push_sample(EnsembleSample::new(single_position(0.0), 3.0).unwrap(), 0)
        .unwrap();
    ensemble
        .set_temperature(Some(Quantity::new(300.0, KELVIN)))
        .unwrap();
    assert_eq!(
        ensemble.get(1).unwrap().as_model_view().kind(),
        StructureKind::Sample {
            weight: 3.0,
            state: 0,
            representative: false,
            energy: kj(2.5).as_ref(),
            temperature: Some(&Quantity::new(300.0, KELVIN)),
        }
    );
    let sample = ensemble.get(0).unwrap();
    assert!(matches!(
        sample.as_model_view().kind(),
        StructureKind::Sample {
            representative: true,
            ..
        }
    ));
    // Materializing a sample as a model drops its sample state.
    assert_eq!(
        sample.to_model().as_model_view().kind(),
        StructureKind::Model
    );

    let mut trajectory = Trajectory::new(topology);
    let mut frame = TrajectoryFrame::new(single_position(0.0));
    frame.set_step(Some(4));
    trajectory.push(frame).unwrap();
    let frame = trajectory.get(0).unwrap();
    assert_eq!(
        frame.as_model_view().kind(),
        StructureKind::Frame {
            time: None,
            step: Some(4),
            velocities: None,
            forces: None,
        }
    );
}

#[test]
fn rigid_transforms_move_models_members_and_frames_with_their_cells() {
    use crate::geometry::{Matrix3, RigidTransform};
    let topology = Arc::new(single_atom_topology());
    let box_cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(1.0, 2.0, 3.0), NANOMETER),
        [true; 3],
    )
    .unwrap();
    let mut model = Model::new(Arc::clone(&topology), single_position(1.0)).unwrap();
    model.set_cell(Some(box_cell));
    // A quarter turn about z followed by a 2 nm shift along x.
    let turn = Matrix3::from_columns(
        Vector3::new(0.0, 1.0, 0.0),
        Vector3::new(-1.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 1.0),
    );
    let transform = RigidTransform::from_canonical(turn, Vector3::new(2.0, 0.0, 0.0)).unwrap();
    let mut ensemble = bundle([model.clone()]).unwrap();
    let mut trajectory = Trajectory::new(Arc::clone(&topology));
    trajectory
        .push(TrajectoryFrame::new(single_position(1.0)).with_cell(box_cell))
        .unwrap();

    model.apply_rigid_transform(&transform).unwrap();
    ensemble
        .get_mut(0)
        .unwrap()
        .apply_rigid_transform(&transform)
        .unwrap();
    trajectory
        .get_mut(0)
        .unwrap()
        .apply_rigid_transform(&transform)
        .unwrap();

    assert_eq!(
        model.positions().values().value(),
        &[Point3::new(2.0, 1.0, 0.0)]
    );
    assert_eq!(
        model.cell().unwrap().vectors().into_value(),
        [
            Vector3::new(0.0, 1.0, 0.0),
            Vector3::new(-2.0, 0.0, 0.0),
            Vector3::new(0.0, 0.0, 3.0),
        ]
    );
    for view in [
        ensemble.get(0).unwrap().as_model_view(),
        trajectory.get(0).unwrap().as_model_view(),
    ] {
        assert_eq!(view.positions(), model.positions());
        assert_eq!(view.cell(), model.cell());
    }
}
