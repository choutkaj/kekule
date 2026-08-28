use super::builder::{checked_future_len, checked_id};
use super::*;
use std::sync::Arc;

use crate::core::{BondOrder, Element};
use crate::properties::{PropertyColumn, PropertyKey, PropertyValue};
use crate::query;
use crate::substructure;
use crate::topology::AtomSiteMetadata;

fn perceived_molecule(smiles: &str) -> Molecule {
    let mut molecule = crate::tests::read_smiles(smiles).expect("SMILES should parse");
    molecule.perceive().expect("molecule should perceive");
    molecule
}

#[test]
fn hierarchy_errors_have_diagnostic_display_messages() {
    let error = TopologyHierarchyError::InconsistentResidueAtomSite {
        residue: ResidueId::new(2),
        site: AtomSiteId::new(7),
    };
    let message = error.to_string();
    assert!(message.contains("residue2"));
    assert!(message.contains("atom-site7"));
    assert!(message.contains("do not reference each other"));
}

fn tombstoned_molecule() -> (Molecule, AtomId, AtomId, BondId) {
    let mut graph = crate::core::MoleculeEditor::new();
    let carbon = graph
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .expect("atom identifier capacity");
    let tombstone = graph
        .add_atom(Atom::new(Element::from_symbol("H").unwrap()))
        .expect("atom identifier capacity");
    let oxygen = graph
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .expect("atom identifier capacity");
    graph.delete_atom(tombstone).unwrap();
    let deleted_bond = graph.add_bond(carbon, oxygen, BondOrder::Single).unwrap();
    graph.delete_bond(deleted_bond).unwrap();
    let bond = graph.add_bond(carbon, oxygen, BondOrder::Double).unwrap();
    (graph.finish().unwrap(), carbon, oxygen, bond)
}

fn topology_with_reused_definition() -> (Arc<Topology>, AtomId, AtomId, BondId) {
    let (molecule, carbon, oxygen, bond) = tombstoned_molecule();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    builder.add_instance(definition).unwrap();
    builder.add_instance(definition).unwrap();
    (Arc::new(builder.build().unwrap()), carbon, oxygen, bond)
}

fn topology_with_two_distinct_definitions(reverse: bool) -> Arc<Topology> {
    let carbon = perceived_molecule("C");
    let carbon_oxygen = perceived_molecule("CO");
    let molecules = if reverse {
        [&carbon_oxygen, &carbon]
    } else {
        [&carbon, &carbon_oxygen]
    };
    let mut builder = TopologyBuilder::new();
    for molecule in molecules {
        let definition = builder.add_molecule_definition(molecule).unwrap();
        builder.add_instance(definition).unwrap();
    }
    Arc::new(builder.build().unwrap())
}

#[test]
fn topology_reuses_definitions_and_preserves_qualified_dense_order() {
    let (topology, carbon, oxygen, bond) = topology_with_reused_definition();
    assert_eq!(topology.definition_count(), 1);
    assert_eq!(topology.instance_count(), 2);
    assert_eq!(topology.atom_count(), 4);
    assert_eq!(topology.bond_count(), 2);

    let first = MoleculeInstanceId::new(0);
    let second = MoleculeInstanceId::new(1);
    assert!(std::ptr::eq(
        topology.definition_for_instance(first).unwrap(),
        topology.definition_for_instance(second).unwrap()
    ));
    assert_eq!(
        topology.atom_ids(),
        &[
            InstanceAtomId::new(first, carbon),
            InstanceAtomId::new(first, oxygen),
            InstanceAtomId::new(second, carbon),
            InstanceAtomId::new(second, oxygen),
        ]
    );
    assert_eq!(
        topology.bond_ids(),
        &[
            InstanceBondId::new(first, bond),
            InstanceBondId::new(second, bond),
        ]
    );
    for (raw, atom) in topology.atom_ids().iter().copied().enumerate() {
        let index = topology.atom_index(atom).unwrap();
        assert_eq!(index.index(), raw);
        assert_eq!(topology.atom_id(index), Some(atom));
    }
    for (raw, bond) in topology.bond_ids().iter().copied().enumerate() {
        let index = topology.bond_index(bond).unwrap();
        assert_eq!(index.index(), raw);
        assert_eq!(topology.bond_id(index), Some(bond));
    }
}

#[test]
fn molecule_views_are_instance_qualified_and_share_definition_state() {
    let (topology, carbon, oxygen, bond) = topology_with_reused_definition();
    let molecules = topology.molecules().collect::<Vec<_>>();
    assert_eq!(molecules.len(), 2);
    assert_eq!(molecules[0].id(), MoleculeInstanceId::new(0));
    assert_eq!(molecules[1].id(), MoleculeInstanceId::new(1));
    assert_eq!(molecules[0].definition_id(), molecules[1].definition_id());
    assert!(std::ptr::eq(
        molecules[0].molecule(),
        molecules[1].molecule()
    ));
    assert_eq!(
        molecules[0].atoms().map(|(id, _)| id).collect::<Vec<_>>(),
        vec![
            InstanceAtomId::new(molecules[0].id(), carbon),
            InstanceAtomId::new(molecules[0].id(), oxygen),
        ]
    );
    assert_eq!(
        molecules[1].bonds().map(|(id, _)| id).collect::<Vec<_>>(),
        vec![InstanceBondId::new(molecules[1].id(), bond)]
    );
    assert_eq!(
        topology.molecule(molecules[1].id()).unwrap().id(),
        molecules[1].id()
    );
}

#[test]
fn builder_rejects_empty_topologies_and_unused_definitions() {
    assert!(matches!(
        TopologyBuilder::new().build(),
        Err(TopologyBuildError::NoMoleculeInstances)
    ));

    let molecule = perceived_molecule("O");
    let mut builder = TopologyBuilder::new();
    let used = builder.add_molecule_definition(&molecule).unwrap();
    builder.add_instance(used).unwrap();
    let unused = builder.add_molecule_definition(&molecule).unwrap();
    assert!(matches!(
        builder.build(),
        Err(TopologyBuildError::UnusedMoleculeDefinition(id)) if id == unused
    ));
}

#[test]
fn builder_add_molecule_is_the_concise_single_instance_path() {
    let molecule = perceived_molecule("CO");
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(&molecule).unwrap();
    let topology = builder.build().unwrap();
    assert_eq!(instance, MoleculeInstanceId::new(0));
    assert_eq!(topology.definition_count(), 1);
    assert_eq!(topology.instance_count(), 1);
    assert_eq!(topology.molecule(instance).unwrap().molecule(), &molecule);
}

#[test]
fn topology_directly_owns_its_layout_collections() {
    let (topology, ..) = topology_with_reused_definition();
    assert_eq!(topology.definitions.len(), 1);
    assert_eq!(topology.instances.len(), 2);
    assert_eq!(topology.instance_atoms.len(), 4);
    assert_eq!(topology.instance_bonds.len(), 2);
    assert_eq!(topology.atom_indices.len(), 4);
    assert_eq!(topology.bond_indices.len(), 2);
    for &atom in &topology.instance_atoms {
        assert_eq!(
            topology.instance_atoms[topology.atom_indices[&atom].index()],
            atom
        );
    }
    for &bond in &topology.instance_bonds {
        assert_eq!(
            topology.instance_bonds[topology.bond_indices[&bond].index()],
            bond
        );
    }

    let debug = format!("{topology:?}");
    assert!(debug.contains("instance_atoms"));
    assert!(debug.contains("instance_bonds"));
}

#[test]
fn shared_allocation_is_exact_while_layout_can_match() {
    let (topology, ..) = topology_with_reused_definition();
    let clone = Arc::clone(&topology);
    let (independent, ..) = topology_with_reused_definition();

    assert!(Arc::ptr_eq(&topology, &clone));
    assert!(topology.same_layout(&clone));
    assert!(!Arc::ptr_eq(&topology, &independent));
    assert!(topology.same_layout(&independent));
    assert_eq!(clone.atom_ids(), topology.atom_ids());
}

#[test]
fn layout_equality_does_not_reorder_definitions_instances_or_dense_state() {
    let forward = topology_with_two_distinct_definitions(false);
    let reverse = topology_with_two_distinct_definitions(true);

    assert!(!Arc::ptr_eq(&forward, &reverse));
    assert!(!forward.same_layout(&reverse));
}

#[test]
fn topology_properties_cover_every_domain_and_do_not_change_layout_identity() {
    let (molecule, carbon, oxygen, bond) = tombstoned_molecule();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    let instance = builder.add_instance(definition).unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "LIG", Some(1), None, None)
        .unwrap();
    builder
        .hierarchy_mut()
        .add_atom_site(
            residue,
            InstanceAtomId::new(instance, carbon),
            AtomSiteMetadata {
                label_atom_id: Some("C1".into()),
                ..AtomSiteMetadata::default()
            },
        )
        .unwrap();

    let owner_key = PropertyKey::new("source").unwrap();
    let value_key = PropertyKey::new("tag").unwrap();
    builder
        .insert_property(owner_key.clone(), PropertyValue::String("test".into()))
        .unwrap();
    fn insert_tag(table: &mut crate::properties::PropertyTable, key: &PropertyKey) {
        table
            .insert(key.clone(), PropertyColumn::Int(vec![Some(1); table.len()]))
            .unwrap();
    }
    insert_tag(builder.molecule_instance_properties_mut(), &value_key);
    insert_tag(builder.atom_properties_mut(), &value_key);
    insert_tag(builder.bond_properties_mut(), &value_key);
    insert_tag(builder.chain_properties_mut(), &value_key);
    insert_tag(builder.residue_properties_mut(), &value_key);
    insert_tag(builder.atom_site_properties_mut(), &value_key);
    let enriched = builder.build().unwrap();

    assert_eq!(
        enriched.properties().get(&owner_key),
        Some(&PropertyValue::String("test".into()))
    );
    assert_eq!(enriched.molecule_instance_properties().len(), 1);
    assert_eq!(enriched.atom_properties().len(), 2);
    assert_eq!(enriched.bond_properties().len(), 1);
    assert_eq!(enriched.chain_properties().len(), 1);
    assert_eq!(enriched.residue_properties().len(), 1);
    assert_eq!(enriched.atom_site_properties().len(), 1);
    assert_eq!(
        enriched
            .molecule_instance_property(instance, &value_key)
            .unwrap(),
        Some(PropertyValue::Int(1))
    );
    assert_eq!(
        enriched
            .molecule(instance)
            .unwrap()
            .property(&value_key)
            .unwrap(),
        Some(PropertyValue::Int(1))
    );
    assert_eq!(
        enriched
            .atom_property(InstanceAtomId::new(instance, oxygen), &value_key)
            .unwrap(),
        Some(PropertyValue::Int(1))
    );
    assert_eq!(
        enriched
            .bond_property(InstanceBondId::new(instance, bond), &value_key)
            .unwrap(),
        Some(PropertyValue::Int(1))
    );

    let mut plain_builder = TopologyBuilder::new();
    let definition = plain_builder.add_molecule_definition(&molecule).unwrap();
    let instance = plain_builder.add_instance(definition).unwrap();
    let chain = plain_builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = plain_builder
        .hierarchy_mut()
        .add_residue(chain, "LIG", Some(1), None, None)
        .unwrap();
    plain_builder
        .hierarchy_mut()
        .add_atom_site(
            residue,
            InstanceAtomId::new(instance, carbon),
            AtomSiteMetadata {
                label_atom_id: Some("C1".into()),
                ..AtomSiteMetadata::default()
            },
        )
        .unwrap();
    let plain = plain_builder.build().unwrap();
    assert!(enriched.same_layout(&plain));
}

#[test]
fn builder_is_transactional_and_does_not_intern_equal_definitions() {
    let (molecule, ..) = tombstoned_molecule();
    let mut builder = TopologyBuilder::new();
    let first = builder.add_molecule_definition(&molecule).unwrap();
    assert_eq!(
        builder.add_instance(MoleculeDefinitionId::new(99)),
        Err(TopologyBuildError::InvalidMoleculeDefinitionId(
            MoleculeDefinitionId::new(99)
        ))
    );
    assert!(builder.instances.is_empty());
    builder.add_instance(first).unwrap();
    let second = builder.add_molecule_definition(&molecule).unwrap();
    builder.add_instance(second).unwrap();
    let topology = Arc::new(builder.build().unwrap());
    assert_eq!(topology.definition_count(), 2);
    assert_eq!(topology.instance_count(), 2);
    assert_eq!(
        checked_future_len(usize::MAX, 1, TopologyIdKind::Atom),
        Err(TopologyBuildError::IdentifierCapacityExceeded(
            TopologyIdKind::Atom
        ))
    );
}

#[cfg(target_pointer_width = "64")]
#[test]
fn every_topology_identifier_space_checks_its_boundary() {
    let max_slot = usize::try_from(u64::from(u32::MAX)).expect("64-bit usize");
    let first_unsupported_slot = usize::try_from(u64::from(u32::MAX) + 1).expect("64-bit usize");

    assert_eq!(
        checked_id::<MoleculeDefinitionId>(max_slot, TopologyIdKind::MoleculeDefinition,),
        Ok(MoleculeDefinitionId::new(u32::MAX))
    );
    assert_eq!(
        checked_id::<MoleculeDefinitionId>(
            first_unsupported_slot,
            TopologyIdKind::MoleculeDefinition,
        ),
        Err(TopologyBuildError::IdentifierCapacityExceeded(
            TopologyIdKind::MoleculeDefinition
        ))
    );
    assert_eq!(
        checked_id::<MoleculeInstanceId>(first_unsupported_slot, TopologyIdKind::MoleculeInstance,),
        Err(TopologyBuildError::IdentifierCapacityExceeded(
            TopologyIdKind::MoleculeInstance
        ))
    );
    assert_eq!(
        checked_id::<TopologyAtomIndex>(first_unsupported_slot, TopologyIdKind::Atom,),
        Err(TopologyBuildError::IdentifierCapacityExceeded(
            TopologyIdKind::Atom
        ))
    );
    assert_eq!(
        checked_id::<TopologyBondIndex>(first_unsupported_slot, TopologyIdKind::Bond,),
        Err(TopologyBuildError::IdentifierCapacityExceeded(
            TopologyIdKind::Bond
        ))
    );
    assert_eq!(
        checked_future_len(first_unsupported_slot, 1, TopologyIdKind::Atom),
        Err(TopologyBuildError::IdentifierCapacityExceeded(
            TopologyIdKind::Atom
        ))
    );
}

#[test]
fn selections_distinguish_instances_elements_and_queries() {
    let ethane = perceived_molecule("CC");
    let water = perceived_molecule("O");
    let mut builder = TopologyBuilder::new();
    let ethane_definition = builder.add_molecule_definition(&ethane).unwrap();
    let water_definition = builder.add_molecule_definition(&water).unwrap();
    let ethane_instance = builder.add_instance(ethane_definition).unwrap();
    let water_instance = builder.add_instance(water_definition).unwrap();
    let topology = Arc::new(builder.build().unwrap());

    let selected =
        AtomSelection::for_instances(&topology, [ethane_instance, water_instance]).unwrap();
    assert_eq!(selected.indices().len(), 3);
    let ethane_selection = AtomSelection::for_instances(&topology, [ethane_instance]).unwrap();
    let water_selection = AtomSelection::for_instances(&topology, [water_instance]).unwrap();
    assert_eq!(ethane_selection.indices().len(), 2);
    assert_eq!(water_selection.indices().len(), 1);
    let oxygen =
        AtomSelection::for_elements(&topology, [Element::from_symbol("O").unwrap()]).unwrap();
    assert_eq!(oxygen.indices().len(), 1);

    let query = query::parse_smarts("O").unwrap();
    let matches = substructure::find_substructure_matches(&water, &query).unwrap();
    let from_query =
        AtomSelection::from_query_matches(&topology, water_instance, &matches).unwrap();
    assert_eq!(
        from_query.semantic_ids(&topology).unwrap(),
        oxygen.semantic_ids(&topology).unwrap()
    );

    let mut independent_builder = TopologyBuilder::new();
    let definition = independent_builder
        .add_molecule_definition(&ethane)
        .unwrap();
    independent_builder.add_instance(definition).unwrap();
    let independent = Arc::new(independent_builder.build().unwrap());
    assert_eq!(
        selected.ensure_compatible(&independent),
        Err(SelectionError::TopologyMismatch)
    );
}

#[test]
fn topology_global_hierarchy_crosses_molecule_boundaries() {
    let mut macro_builder = crate::core::MoleculeEditor::new();
    let atom = macro_builder
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .expect("atom identifier capacity");
    let molecule = macro_builder.finish().unwrap();

    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    let first = builder.add_instance(definition).unwrap();
    let small = perceived_molecule("O");
    let small_definition = builder.add_molecule_definition(&small).unwrap();
    let small_instance = builder.add_instance(small_definition).unwrap();
    let second = builder.add_instance(definition).unwrap();
    let chain = builder
        .hierarchy_mut()
        .add_chain("A", Some("AUTH".into()))
        .unwrap();
    let first_residue = builder
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), Some("10".into()), None)
        .unwrap();
    let second_residue = builder
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(2), Some("11".into()), None)
        .unwrap();
    let first_atom = InstanceAtomId::new(first, atom);
    let second_atom = InstanceAtomId::new(second, atom);
    let first_site = builder
        .hierarchy_mut()
        .add_atom_site(first_residue, first_atom, AtomSiteMetadata::default())
        .unwrap();
    let second_site = builder
        .hierarchy_mut()
        .add_atom_site(second_residue, second_atom, AtomSiteMetadata::default())
        .unwrap();
    let topology = Arc::new(builder.build().unwrap());
    assert_eq!(topology.definition_count(), 2);
    assert_eq!(topology.hierarchy().chains().count(), 1);
    let first_molecule = topology.molecule(first).unwrap();
    assert_eq!(first_molecule.molecule(), &molecule);
    assert_eq!(
        first_molecule
            .chains()
            .map(ChainView::id)
            .collect::<Vec<_>>(),
        vec![chain]
    );
    assert!(topology
        .molecule(small_instance)
        .unwrap()
        .chains()
        .next()
        .is_none());

    assert_eq!(
        topology.chains().map(ChainView::id).collect::<Vec<_>>(),
        vec![chain]
    );
    assert_eq!(
        topology.residues().map(ResidueView::id).collect::<Vec<_>>(),
        vec![first_residue, second_residue]
    );
    assert_eq!(
        topology
            .atom_sites()
            .map(AtomSiteView::id)
            .collect::<Vec<_>>(),
        vec![first_site, second_site]
    );

    let chain_view = topology.chain(chain).unwrap();
    assert_eq!(chain_view.residues().count(), 2);

    assert_eq!(topology.atom_for_site(first_site).unwrap(), first_atom);
    assert_eq!(
        topology
            .atom_site_for_atom(first_atom)
            .unwrap()
            .unwrap()
            .id(),
        first_site
    );
    assert_eq!(
        topology.residue_for_atom(first_atom).unwrap().unwrap().id(),
        first_residue
    );
    assert_eq!(
        topology.chain_for_atom(first_atom).unwrap().unwrap().id(),
        chain
    );
    assert_eq!(
        topology.residue_for_site(first_site).unwrap().id(),
        first_residue
    );
    assert_eq!(
        topology.chain_for_residue(first_residue).unwrap().id(),
        chain
    );

    let small_atom = topology
        .atom_ids()
        .iter()
        .copied()
        .find(|id| id.molecule() == small_instance)
        .unwrap();
    assert!(topology.atom_site_for_atom(small_atom).unwrap().is_none());
    assert!(topology.residue_for_atom(small_atom).unwrap().is_none());
    assert!(topology.chain_for_atom(small_atom).unwrap().is_none());

    assert_eq!(
        AtomSelection::for_chains(&topology, [chain])
            .unwrap()
            .semantic_ids(&topology)
            .unwrap(),
        vec![first_atom, second_atom]
    );
    assert_eq!(
        AtomSelection::for_residues(&topology, [second_residue])
            .unwrap()
            .semantic_ids(&topology)
            .unwrap(),
        vec![second_atom]
    );
    assert_eq!(
        AtomSelection::for_atom_sites(&topology, [first_site])
            .unwrap()
            .semantic_ids(&topology)
            .unwrap(),
        vec![first_atom]
    );
    assert_eq!(
        AtomSelection::for_chain_label(&topology, "A")
            .unwrap()
            .semantic_ids(&topology)
            .unwrap(),
        vec![first_atom, second_atom]
    );

    let chain_selection = AtomSelection::for_chains(&topology, [chain]).unwrap();
    let subset = topology.subset(&chain_selection).unwrap();
    assert_eq!(subset.topology().instance_count(), 2);
    assert_eq!(subset.topology().chains().count(), 1);
    assert_eq!(subset.topology().residues().count(), 2);
    assert_eq!(subset.topology().atom_sites().count(), 2);
    assert_eq!(
        subset.correspondence().source_atom_indices(),
        [TopologyAtomIndex::new(0), TopologyAtomIndex::new(2)]
    );
    assert!(subset.correspondence().target_atom(small_atom).is_none());
}

#[test]
fn induced_subset_splits_molecules_and_filters_hierarchy_deterministically() {
    let mut editor = crate::core::MoleculeEditor::new();
    let first = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let middle = editor
        .add_atom(Atom::new(Element::from_symbol("N").unwrap()))
        .unwrap();
    let tombstone = editor
        .add_atom(Atom::new(Element::from_symbol("H").unwrap()))
        .unwrap();
    editor.delete_atom(tombstone).unwrap();
    let last = editor
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    editor.add_bond(first, middle, BondOrder::Single).unwrap();
    editor.add_bond(middle, last, BondOrder::Double).unwrap();
    let molecule = editor.finish().unwrap();
    assert_eq!(last.raw(), 3);

    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(&molecule).unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    for (sequence, atom) in [(1, first), (2, middle), (3, last)] {
        let residue = builder
            .hierarchy_mut()
            .add_residue(chain, "RES", Some(sequence), None, None)
            .unwrap();
        builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, atom),
                AtomSiteMetadata::default(),
            )
            .unwrap();
    }
    let owner_key = PropertyKey::new("owner_note").unwrap();
    let value_key = PropertyKey::new("source_index").unwrap();
    builder
        .insert_property(owner_key.clone(), PropertyValue::Bool(true))
        .unwrap();
    builder
        .molecule_instance_properties_mut()
        .insert(value_key.clone(), PropertyColumn::Int(vec![Some(10)]))
        .unwrap();
    builder
        .atom_properties_mut()
        .insert(
            value_key.clone(),
            PropertyColumn::Int(vec![Some(1), Some(2), Some(3)]),
        )
        .unwrap();
    builder
        .bond_properties_mut()
        .insert(
            value_key.clone(),
            PropertyColumn::Int(vec![Some(4), Some(5)]),
        )
        .unwrap();
    builder
        .chain_properties_mut()
        .insert(value_key.clone(), PropertyColumn::Int(vec![Some(6)]))
        .unwrap();
    builder
        .residue_properties_mut()
        .insert(
            value_key.clone(),
            PropertyColumn::Int(vec![Some(7), Some(8), Some(9)]),
        )
        .unwrap();
    builder
        .atom_site_properties_mut()
        .insert(
            value_key.clone(),
            PropertyColumn::Int(vec![Some(10), Some(11), Some(12)]),
        )
        .unwrap();
    let source = Arc::new(builder.build().unwrap());
    let whole_selection = AtomSelection::from_atoms(
        &source,
        [first, middle, last].map(|atom| InstanceAtomId::new(instance, atom)),
    )
    .unwrap();
    let whole = source.subset(&whole_selection).unwrap();
    assert_eq!(
        whole
            .topology()
            .molecule_instance_properties()
            .get(&value_key),
        Some(&PropertyColumn::Int(vec![Some(10)]))
    );

    let partial_selection = AtomSelection::from_atoms(
        &source,
        [first, middle].map(|atom| InstanceAtomId::new(instance, atom)),
    )
    .unwrap();
    let partial = source.subset(&partial_selection).unwrap();
    assert!(!partial.topology().molecule_instance_properties().has_data());

    let selection = AtomSelection::from_atoms(
        &source,
        [
            InstanceAtomId::new(instance, first),
            InstanceAtomId::new(instance, last),
        ],
    )
    .unwrap();

    let subset = source.subset(&selection).unwrap();
    assert_eq!(subset.topology().instance_count(), 2);
    assert_eq!(subset.topology().atom_count(), 2);
    assert_eq!(subset.topology().bond_count(), 0);
    assert_eq!(subset.topology().chains().count(), 1);
    assert_eq!(subset.topology().residues().count(), 2);
    assert_eq!(subset.topology().atom_sites().count(), 2);
    assert_eq!(
        subset
            .topology()
            .residues()
            .map(|residue| residue.label_seq_id())
            .collect::<Vec<_>>(),
        [Some(1), Some(3)]
    );
    assert_eq!(
        subset.correspondence().source_atom_indices(),
        [TopologyAtomIndex::new(0), TopologyAtomIndex::new(2)]
    );
    assert!(subset
        .correspondence()
        .target_atom(InstanceAtomId::new(instance, middle))
        .is_none());
    let target_last = subset
        .correspondence()
        .target_atom(InstanceAtomId::new(instance, last))
        .expect("selected tombstone-separated atom is mapped");
    assert_eq!(subset.topology().atom_index(target_last).unwrap().raw(), 1);
    let projected = subset.topology();
    assert!(projected.properties().get(&owner_key).is_none());
    assert!(!projected.molecule_instance_properties().has_data());
    assert_eq!(
        projected.atom_properties().get(&value_key).unwrap(),
        &PropertyColumn::Int(vec![Some(1), Some(3)])
    );
    assert!(!projected.bond_properties().has_data());
    assert_eq!(
        projected.chain_properties().get(&value_key).unwrap(),
        &PropertyColumn::Int(vec![Some(6)])
    );
    assert_eq!(
        projected.residue_properties().get(&value_key).unwrap(),
        &PropertyColumn::Int(vec![Some(7), Some(9)])
    );
    assert_eq!(
        projected.atom_site_properties().get(&value_key).unwrap(),
        &PropertyColumn::Int(vec![Some(10), Some(12)])
    );
}

#[test]
fn retaining_one_whole_instance_projects_exactly_one_instance_property_row() {
    let mut editor = crate::core::MoleculeEditor::new();
    editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let molecule = editor.finish().unwrap();
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule).unwrap();
    builder.add_instance(definition).unwrap();
    let second = builder.add_instance(definition).unwrap();
    let key = PropertyKey::new("instance_score").unwrap();
    builder
        .molecule_instance_properties_mut()
        .insert(key.clone(), PropertyColumn::Int(vec![Some(10), Some(20)]))
        .unwrap();
    let source = Arc::new(builder.build().unwrap());

    let retained = transform::retain_instances(&source, [second]).unwrap();
    assert_eq!(retained.instance_count(), 1);
    assert_eq!(
        retained.molecule_instance_properties().get(&key),
        Some(&PropertyColumn::Int(vec![Some(20)]))
    );
}
