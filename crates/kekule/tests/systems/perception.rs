use std::error::Error;

use crate::support::hierarchy_record;
use std::sync::{Arc, Weak};

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyColumn, PropertyKey, PropertyValue};
use kekule::structure::{Ensemble, EnsembleMember, Model, Positions};
use kekule::topology::{
    AtomSelection, AtomSiteMetadata, InstanceAtomId, MoleculeClass, MoleculeDefinitionId,
    ResidueClass, SelectionError, Topology, TopologyAtomIndex, TopologyBondIndex, TopologyBuilder,
};
use kekule::units::{Quantity, NANOMETER};
use kekule::{
    molecule::{Atom, BondOrder, Element, ImplicitHydrogens, MoleculeEditor},
    perception::Perception,
};
use kekule::{perception, smiles, stereo};

fn key() -> PropertyKey {
    PropertyKey::new("perception_regression").unwrap()
}

fn assert_default_perception(topology: &Topology) {
    for definition in topology.definitions() {
        let molecule = definition.molecule();
        let mut expected = molecule.clone();
        expected.clear_perception();
        expected.perceive().unwrap();
        assert_eq!(molecule.perception(), expected.perception());
        assert!(molecule.perception().valence().is_some());
        assert!(molecule.perception().rings().is_some());
        assert!(molecule.perception().aromaticity().is_some());
        assert!(molecule.perception().conjugation().is_some());
        assert!(molecule.perception().resonance().is_none());
        assert!(molecule.perception().cip().is_none());
    }
}

fn failing_topology() -> Arc<Topology> {
    let mut good = smiles::to_molecules("c1ccccc1[C@H](F)Cl")
        .unwrap()
        .remove(0);
    good.perceive().unwrap();
    stereo::cip::assign_cip_descriptors(&mut good, Default::default()).unwrap();
    perception::resonance::perceive_resonance(&mut good).unwrap();
    assert!(good
        .perception()
        .cip()
        .is_some_and(|cip| cip.cip_descriptors().len() != 0));
    assert!(good.perception().resonance().is_some());

    // Publication permits represented chemistry rejected by the default valence model.
    let mut editor = MoleculeEditor::new();
    let carbon =
        Atom::new(Element::from_symbol("C").unwrap()).with_hydrogens(ImplicitHydrogens::Fixed(5));
    editor.add_atom(carbon).unwrap();
    let mut bad = editor.finish().unwrap();
    perception::rings::perceive_rings(
        &mut bad,
        perception::rings::RingPerceptionOptions::membership_only(),
    )
    .expect("membership perception cannot fail");
    Arc::new(Topology::from_molecules([good, bad].clone()).unwrap())
}

fn installed(topology: &Topology) -> Vec<Perception> {
    topology
        .definitions()
        .map(|d| d.molecule().perception().clone())
        .collect()
}

#[test]
fn perceived_topology_preserves_edited_ids_reuse_hierarchy_and_annotations() {
    let mut editor = MoleculeEditor::new();
    let carbon = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let deleted = editor
        .add_atom(Atom::new(Element::from_symbol("H").unwrap()))
        .unwrap();
    let oxygen = editor
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    editor.delete_atom(deleted).unwrap();
    let deleted_bond = editor.add_bond(carbon, oxygen, BondOrder::Single).unwrap();
    editor.delete_bond(deleted_bond).unwrap();
    let bond = editor.add_bond(carbon, oxygen, BondOrder::Double).unwrap();
    let (mut molecule, ids) = editor.finish_with_correspondence().unwrap();
    // Publication renumbers the surviving atoms and bond densely.
    let (carbon, oxygen, bond) = (
        ids.atom(carbon).unwrap(),
        ids.atom(oxygen).unwrap(),
        ids.bond(bond).unwrap(),
    );
    molecule
        .properties_mut()
        .owner_mut()
        .insert(key(), PropertyValue::Int(1))
        .unwrap();
    molecule
        .properties_mut()
        .atoms_mut()
        .set_value(key(), carbon, Some(PropertyValue::Int(2)))
        .unwrap();
    molecule
        .properties_mut()
        .bonds_mut()
        .set_value(key(), bond, Some(PropertyValue::Int(3)))
        .unwrap();

    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(molecule.clone()).unwrap();
    let first = builder.add_instance(definition).unwrap();
    let second = builder.add_instance(definition).unwrap();
    builder
        .set_molecule_class(definition, MoleculeClass::Other)
        .unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "LIG", Some(1), None, None)
        .unwrap();
    builder
        .hierarchy_mut()
        .set_residue_class(residue, ResidueClass::Other)
        .unwrap();
    for instance in [first, second] {
        for atom in [carbon, oxygen] {
            builder
                .hierarchy_mut()
                .add_atom_site(
                    residue,
                    InstanceAtomId::new(instance, atom),
                    AtomSiteMetadata::default(),
                )
                .unwrap();
        }
    }
    builder
        .properties_mut()
        .owner_mut()
        .insert(key(), PropertyValue::Int(4))
        .unwrap();
    builder
        .properties_mut()
        .atoms_mut()
        .insert(key(), PropertyColumn::Int(vec![Some(5); 4]))
        .unwrap();
    let source = Arc::new(builder.build().unwrap());
    let target = Arc::new(source.perceived().unwrap());

    assert!(source.same_layout(&target));
    assert!(target.same_layout(&source));
    assert!(!Arc::ptr_eq(&source, &target));
    assert_eq!(target.atom_ids(), source.atom_ids());
    assert_eq!(target.bond_ids(), source.bond_ids());
    assert_eq!(hierarchy_record(&target), hierarchy_record(&source));
    assert_eq!(target.properties(), source.properties());
    assert_eq!(target.definition_count(), 1);
    assert_eq!(target.instance_count(), 2);
    let first = target.instance(first).unwrap();
    let second = target.instance(second).unwrap();
    assert!(std::ptr::eq(first.molecule(), second.molecule()));
    assert_eq!(first.class(), MoleculeClass::Other);
    assert_eq!(first.molecule(), &molecule);
    assert_eq!(first.molecule().properties(), molecule.properties());
    assert_eq!(
        first
            .molecule()
            .perception()
            .valence()
            .and_then(|valence| valence.inferred_hydrogen_count(carbon)),
        Some(2)
    );
    assert_eq!(
        source
            .definition(definition)
            .unwrap()
            .molecule()
            .perception(),
        &Perception::default()
    );
    assert_default_perception(&target);
}

#[test]
fn model_perception_changes_snapshot_preserving_structure_and_existing_bindings() {
    let topology = smiles::to_topology("c1ccccc1").unwrap();
    let selection =
        AtomSelection::from_atoms(Arc::clone(&topology), [topology.atom_ids()[0]]).unwrap();
    let positions = Positions::new(Quantity::new(
        (0..6)
            .map(|i| Point3::new(f64::from(i), 1.0, 2.0))
            .collect::<Vec<_>>(),
        NANOMETER,
    ))
    .unwrap();
    let mut model = Model::new(Arc::clone(&topology), positions).unwrap();
    let mut properties = model.properties_mut();
    properties
        .owner_mut()
        .insert(key(), PropertyValue::Int(7))
        .unwrap();
    properties
        .atoms_mut()
        .set_value(
            key(),
            TopologyAtomIndex::new(0),
            Some(PropertyValue::Int(8)),
        )
        .unwrap();
    properties
        .bonds_mut()
        .set_value(
            key(),
            TopologyBondIndex::new(0),
            Some(PropertyValue::Int(9)),
        )
        .unwrap();
    model.set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(10.0, 10.0, 10.0), NANOMETER),
            [true; 3],
        )
        .unwrap(),
    ));
    let original = model.clone();
    let positions_ptr = model.positions().values().value().as_ptr();
    model.perceive().unwrap();

    assert_default_perception(model.topology());
    assert!(topology.same_layout(model.topology()));
    assert!(!Arc::ptr_eq(&topology, &model.shared_topology()));
    assert!(Arc::ptr_eq(&topology, &original.shared_topology()));
    assert_eq!(model.positions(), original.positions());
    assert_eq!(positions_ptr, model.positions().values().value().as_ptr());
    assert_eq!(model.cell(), original.cell());
    assert_eq!(model.properties(), original.properties());
    // The new snapshot keeps the layout, so bindings made before perception
    // remain usable with the perceived model.
    assert!(model.topology().shares_layout(&topology));
    selection
        .ensure_compatible(&model.shared_topology())
        .unwrap();
    selection
        .ensure_compatible(&original.shared_topology())
        .unwrap();
    let sliced = model.subset(&selection).unwrap();
    assert_eq!(sliced.atom_count(), 1);
    let independent = smiles::to_topology("c1ccccc1").unwrap();
    assert_eq!(
        selection.ensure_compatible(&independent),
        Err(SelectionError::TopologyMismatch)
    );
    assert!(installed(&topology)
        .iter()
        .all(|p| p == &Perception::default()));
}

#[test]
fn reperception_replaces_even_a_uniquely_owned_snapshot_and_clears_cip() {
    let mut molecule = smiles::to_molecules("c1ccccc1[C@H](F)Cl")
        .unwrap()
        .remove(0);
    molecule.perceive().unwrap();
    stereo::cip::assign_cip_descriptors(&mut molecule, Default::default()).unwrap();
    assert!(molecule.perception().cip().is_some());
    let mut model =
        Model::from_molecule(molecule.clone(), &Positions::zeros(molecule.atom_count())).unwrap();
    let old = Arc::downgrade(&model.shared_topology());
    assert_eq!(old.strong_count(), 1);
    model.perceive().unwrap();
    assert!(!Weak::ptr_eq(
        &old,
        &Arc::downgrade(&model.shared_topology())
    ));
    assert_default_perception(model.topology());
    assert_eq!(
        model.topology().instances().next().unwrap().molecule(),
        &molecule
    );
}

#[test]
fn ensemble_perception_preserves_members_weights_and_collection_properties() {
    let topology = smiles::to_topology("c1ccccc1").unwrap();
    let mut member = EnsembleMember::new(Positions::zeros(6), 1.0).unwrap();
    member.set_weight(0.25).unwrap();
    member.set_cell(Some(
        PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(10.0, 10.0, 10.0), NANOMETER),
            [true; 3],
        )
        .unwrap(),
    ));
    let mut properties = member.properties_mut();
    properties
        .owner_mut()
        .insert(key(), PropertyValue::Int(10))
        .unwrap();
    properties
        .atoms_mut()
        .set_value(
            key(),
            TopologyAtomIndex::new(0),
            Some(PropertyValue::Int(11)),
        )
        .unwrap();
    let mut second = EnsembleMember::new(Positions::zeros(6), 1.0).unwrap();
    second.set_weight(0.75).unwrap();
    let mut ensemble = Ensemble::from_items(Arc::clone(&topology), [member, second]).unwrap();
    // Bond rows exist once the member is bound to the ensemble topology.
    ensemble
        .get_mut(0)
        .unwrap()
        .properties_mut()
        .bonds_mut()
        .insert(key(), PropertyColumn::Int(vec![Some(12); 6]))
        .unwrap();
    ensemble
        .properties_mut()
        .insert(key(), PropertyValue::Int(13))
        .unwrap();
    let before = ensemble.clone();
    let pointers = ensemble
        .iter()
        .map(|m| m.positions().values().value().as_ptr())
        .collect::<Vec<_>>();
    ensemble.perceive().unwrap();

    assert_default_perception(ensemble.topology());
    assert!(topology.same_layout(ensemble.topology()));
    assert!(!Arc::ptr_eq(&topology, &ensemble.shared_topology()));
    assert!(Arc::ptr_eq(&topology, &before.shared_topology()));
    assert_eq!(ensemble.properties(), before.properties());
    for ((member, original), pointer) in ensemble.iter().zip(before.iter()).zip(pointers) {
        assert!(Arc::ptr_eq(
            &member.shared_topology(),
            &ensemble.shared_topology()
        ));
        assert_eq!(member.positions(), original.positions());
        assert_eq!(member.positions().values().value().as_ptr(), pointer);
        assert_eq!(member.cell(), original.cell());
        assert_eq!(member.properties(), original.properties());
        assert_eq!(member.weight(), original.weight());
    }
    let mut empty = Ensemble::new(topology);
    empty.perceive().unwrap();
    assert_default_perception(empty.topology());
    assert_eq!(empty.iter().len(), 0);
}

#[test]
fn topology_model_and_ensemble_failure_preserve_complete_previous_state() {
    let topology = failing_topology();
    let previous = installed(&topology);
    let error = topology.perceived().unwrap_err();
    assert_eq!(error.definition, MoleculeDefinitionId::new(1));
    assert!(matches!(
        error.source,
        perception::PerceptionError::Valence(_)
    ));
    assert!(error.to_string().contains(&error.definition.to_string()));
    assert!(error.source().is_some());
    assert_eq!(installed(&topology), previous);

    let mut model = Model::new(
        Arc::clone(&topology),
        Positions::zeros(topology.atom_count()),
    )
    .unwrap();
    model
        .properties_mut()
        .owner_mut()
        .insert(key(), PropertyValue::Int(14))
        .unwrap();
    let before = model.clone();
    assert_eq!(model.perceive(), Err(error.clone()));
    assert_eq!(model, before);
    assert_eq!(installed(model.topology()), previous);

    let mut ensemble = Ensemble::from_items(
        Arc::clone(&topology),
        [EnsembleMember::new(Positions::zeros(topology.atom_count()), 1.0).unwrap()],
    )
    .unwrap();
    ensemble
        .properties_mut()
        .insert(key(), PropertyValue::Int(15))
        .unwrap();
    let before = ensemble.clone();
    assert_eq!(ensemble.perceive(), Err(error));
    assert!(Arc::ptr_eq(&topology, &ensemble.shared_topology()));
    assert_eq!(ensemble.properties(), before.properties());
    assert_eq!(
        ensemble.get(0).unwrap().to_model(),
        before.get(0).unwrap().to_model()
    );
    assert_eq!(installed(ensemble.topology()), previous);
}

#[test]
fn perceived_snapshots_share_one_layout_and_publications_create_new_ones() {
    use kekule::substructure::{PreparedTopologyTarget, SubstructureMatchOptions};
    use kekule::topology::edit::TopologyEditor;

    let source = smiles::to_topology("c1ccccc1.[Na+]").unwrap();
    let perceived = Arc::new(source.perceived().unwrap());
    let again = Arc::new(perceived.perceived().unwrap());
    assert!(perceived.shares_layout(&source) && again.shares_layout(&source));
    // Perception shares the layout storage instead of copying it.
    assert!(std::ptr::eq(
        source.atom_ids().as_ptr(),
        perceived.atom_ids().as_ptr()
    ));
    assert!(std::ptr::eq(source.properties(), perceived.properties()));

    // Publication always creates a new layout, even with equal contents.
    let rebuilt = smiles::to_topology("c1ccccc1.[Na+]").unwrap();
    assert!(rebuilt.same_layout(&source) && !rebuilt.shares_layout(&source));
    let republished = Topology::from_molecules(
        (source
            .instances()
            .map(|m| m.molecule().clone())
            .collect::<Vec<_>>())
        .clone(),
    )
    .unwrap();
    assert!(!republished.shares_layout(&source));
    let unchanged = TopologyEditor::from_topology(Arc::clone(&perceived))
        .finish()
        .unwrap();
    assert!(Arc::ptr_eq(&unchanged, &perceived));
    let mut editor = TopologyEditor::from_topology(Arc::clone(&perceived));
    let sodium = editor.atom_handle(perceived.atom_ids()[6]).unwrap();
    editor.delete_atom(sodium).unwrap();
    assert!(!editor.finish().unwrap().shares_layout(&source));

    // Selections from either snapshot combine, compare equal, and subset both.
    let ring =
        AtomSelection::from_atoms(Arc::clone(&source), source.atom_ids()[..6].iter().copied())
            .unwrap();
    let ion = AtomSelection::from_atoms(Arc::clone(&perceived), [perceived.atom_ids()[6]]).unwrap();
    assert_eq!(
        ring.union(&ion).unwrap(),
        AtomSelection::all(Arc::clone(&perceived))
    );
    assert_eq!(
        AtomSelection::all(Arc::clone(&source)),
        AtomSelection::all(Arc::clone(&perceived))
    );
    let subset = kekule::topology::Topology::subset(perceived.clone(), &ring).unwrap();
    assert_eq!(subset.topology().atom_count(), 6);
    // The correspondence keeps the receiver snapshot, not the selection's.
    assert!(std::ptr::eq(
        subset.correspondence().source_topology(),
        perceived.as_ref()
    ));
    assert_eq!(
        ring.union(&AtomSelection::all(Arc::clone(&rebuilt))),
        Err(SelectionError::TopologyMismatch)
    );

    // Perception-dependent matches keep their exact snapshot but select atoms
    // by layout.
    let query = kekule::smarts::parse_str("c", Default::default()).unwrap();
    let matches = PreparedTopologyTarget::new(&perceived)
        .all(&query, SubstructureMatchOptions::default())
        .unwrap();
    assert_eq!(matches.len(), 6);
    assert!(matches
        .iter()
        .all(|m| Arc::ptr_eq(&m.shared_topology(), &perceived)));
    assert_eq!(
        AtomSelection::from_topology_query_matches(Arc::clone(&source), &matches).unwrap(),
        ring
    );
    assert_eq!(
        AtomSelection::from_topology_query_matches(Arc::clone(&rebuilt), &matches),
        Err(SelectionError::TopologyMismatch)
    );
}

/// A cancellation requested from the `after`th poll on, counting every poll.
fn requested_after(
    after: usize,
) -> (
    kekule::diagnostics::Cancellation,
    Arc<std::sync::atomic::AtomicUsize>,
) {
    use std::sync::atomic::{AtomicUsize, Ordering};
    let polls = Arc::new(AtomicUsize::new(0));
    let counter = Arc::clone(&polls);
    let cancellation = kekule::diagnostics::Cancellation::new(move || {
        counter.fetch_add(1, Ordering::Relaxed) + 1 >= after
    });
    (cancellation, polls)
}

fn assert_unperceived(topology: &Topology) {
    assert!(topology.definitions().all(|definition| definition
        .molecule()
        .perception()
        .rings()
        .is_none()));
}

#[test]
fn cancellable_perception_stops_once_cancellation_is_requested() {
    use std::sync::atomic::Ordering;

    let document = kekule::mmcif::parse_str(
        include_str!("../fixtures/mmcif/1AKE.cif"),
        Default::default(),
    )
    .unwrap();
    let model = document.interpret(Default::default()).unwrap().into_model();
    let ensemble = document
        .interpret_ensemble(Default::default())
        .unwrap()
        .into_ensemble();
    assert_unperceived(model.topology());
    let mut expected = model.clone();
    expected.perceive().unwrap();

    // A check that never turns true changes nothing.
    let (never, polls) = requested_after(usize::MAX);
    let mut perceived = model.clone();
    perceived.perceive_cancellable(&never).unwrap();
    assert!(perceived.topology().shares_layout(model.topology()));
    for (left, right) in perceived
        .topology()
        .definitions()
        .zip(expected.topology().definitions())
    {
        assert_eq!(left.molecule().perception(), right.molecule().perception());
    }
    let polls = polls.load(Ordering::Relaxed);
    assert!(
        polls > model.topology().definitions().len(),
        "perception polled only {polls} times"
    );

    // A check that turns true partway stops perception and keeps the owner.
    for after in [1, polls / 2, polls] {
        let (cancellation, _) = requested_after(after);
        let mut cancelled = model.clone();
        let error = cancelled.perceive_cancellable(&cancellation).unwrap_err();
        assert!(
            matches!(error.source, perception::PerceptionError::Cancelled(_)),
            "poll {after}: {error}"
        );
        assert!(Arc::ptr_eq(
            &cancelled.shared_topology(),
            &model.shared_topology()
        ));
        assert_unperceived(cancelled.topology());

        let (cancellation, _) = requested_after(after);
        let mut cancelled = ensemble.clone();
        let error = cancelled.perceive_cancellable(&cancellation).unwrap_err();
        assert!(matches!(
            error.source,
            perception::PerceptionError::Cancelled(_)
        ));
        assert_unperceived(cancelled.topology());
    }
    let (never, _) = requested_after(usize::MAX);
    let mut perceived = ensemble.clone();
    perceived.perceive_cancellable(&never).unwrap();
    assert_default_perception(perceived.topology());
}
