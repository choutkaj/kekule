use std::sync::Arc;

use kekule::{
    core::Element,
    mmcif, smiles,
    structure::measure,
    topology::{
        AtomSelection, AtomSiteMetadata, HierarchyLookupError, InstanceAtomId, ResidueClass,
        SelectionError, Topology, TopologyBuilder,
    },
    units::{Quantity, ANGSTROM},
};

fn topology() -> Arc<Topology> {
    let molecule = smiles::to_molecules("CC").unwrap().pop().unwrap();
    let local = molecule.atoms().map(|(id, _)| id).collect::<Vec<_>>();
    let mut builder = TopologyBuilder::new();
    let mut atoms = Vec::new();
    for _ in 0..3 {
        let instance = builder.add_molecule(&molecule).unwrap();
        atoms.extend(local.iter().map(|id| InstanceAtomId::new(instance, *id)));
    }
    let chain = builder
        .hierarchy_mut()
        .add_chain("L", Some("A".into()))
        .unwrap();
    let first = builder
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), Some("7".into()), None)
        .unwrap();
    let second = builder
        .hierarchy_mut()
        .add_residue(chain, "ALA", Some(2), Some("7".into()), Some("B".into()))
        .unwrap();
    builder
        .set_residue_class(first, ResidueClass::AminoAcid)
        .unwrap();
    builder
        .set_residue_class(second, ResidueClass::AminoAcid)
        .unwrap();
    // Each residue crosses molecule instances. Two atoms have no hierarchy.
    for (index, residue, label, author) in [
        (0, first, "CA", "alpha"),
        (2, first, "N", "nitrogen"),
        (1, second, "CA", "alpha"),
        (3, second, "CA", "other"),
    ] {
        builder
            .hierarchy_mut()
            .add_atom_site(
                residue,
                atoms[index],
                AtomSiteMetadata {
                    label_atom_id: Some(label.into()),
                    auth_atom_id: Some(author.into()),
                    ..Default::default()
                },
            )
            .unwrap();
    }
    // One author chain identifier can describe multiple label chains.
    builder
        .hierarchy_mut()
        .add_chain("M", Some("A".into()))
        .unwrap();
    Arc::new(builder.build().unwrap())
}

#[test]
fn set_operations_preserve_dense_order_and_exact_snapshot_identity() {
    let top = topology();
    let ids = top.atom_ids();
    let a = AtomSelection::from_atoms(&top, [ids[4], ids[2], ids[0], ids[2]]).unwrap();
    let b = AtomSelection::from_atoms(&top, [ids[5], ids[2], ids[1]]).unwrap();
    assert_eq!(
        a.union(&b).unwrap().atom_ids().collect::<Vec<_>>(),
        [ids[0], ids[1], ids[2], ids[4], ids[5]]
    );
    assert_eq!(
        a.intersection(&b).unwrap().atom_ids().collect::<Vec<_>>(),
        [ids[2]]
    );
    assert_eq!(
        a.difference(&b).unwrap().atom_ids().collect::<Vec<_>>(),
        [ids[0], ids[4]]
    );
    let empty = AtomSelection::from_atoms(&top, []).unwrap();
    assert_eq!(a.union(&empty).unwrap(), a);
    assert_eq!(empty.union(&a).unwrap(), a);
    assert_eq!(a.intersection(&empty).unwrap(), empty);
    assert_eq!(a.difference(&a).unwrap(), empty);
    assert_eq!(empty.difference(&a).unwrap(), empty);
    assert_eq!(a.atom_ids().len(), 3);
    let other = AtomSelection::all(&topology());
    assert_eq!(a.union(&other), Err(SelectionError::TopologyMismatch));
    assert_eq!(
        empty.intersection(&other),
        Err(SelectionError::TopologyMismatch)
    );
    assert_eq!(a.difference(&other), Err(SelectionError::TopologyMismatch));
}

#[test]
fn residue_expansion_crosses_instances_and_preserves_unassigned_atoms() {
    let top = topology();
    let ids = top.atom_ids();
    let selected = AtomSelection::from_atoms(&top, [ids[0], ids[4]]).unwrap();
    let expanded = selected.expand_to_residues();
    assert_eq!(
        expanded.atom_ids().collect::<Vec<_>>(),
        [ids[0], ids[2], ids[4]]
    );
    assert_eq!(expanded.expand_to_residues(), expanded);
    assert!(std::ptr::eq(expanded.topology(), top.as_ref()));
    let empty = AtomSelection::from_atoms(&top, []).unwrap();
    assert_eq!(empty.expand_to_residues(), empty);
    let bare = AtomSelection::from_atoms(&top, [ids[5]]).unwrap();
    assert_eq!(bare.expand_to_residues(), bare);
}

#[test]
fn lookup_respects_identifier_namespace_scope_insertion_and_ambiguity() {
    let top = topology();
    let chain = top.chain_by_label("L").unwrap();
    assert!(matches!(
        top.chain_by_author("A"),
        Err(HierarchyLookupError::Ambiguous { matches: 2, .. })
    ));
    assert!(matches!(
        top.chain_by_author("L"),
        Err(HierarchyLookupError::NotFound { .. })
    ));
    let first = chain.residue_by_author("7", None).unwrap();
    let second = chain.residue_by_author("7", Some("B")).unwrap();
    assert_eq!(first.id(), chain.residue_by_label(1).unwrap().id());
    assert_eq!(second.id(), chain.residue_by_label(2).unwrap().id());
    assert_eq!(
        top.residue_by_author("A", "7", None).unwrap().id(),
        first.id()
    );
    assert!(matches!(
        chain.residue_by_author("7", Some("C")),
        Err(HierarchyLookupError::NotFound { .. })
    ));
    assert_eq!(
        first.atom_site_by_label("CA").unwrap().atom(),
        top.atom_ids()[0]
    );
    assert_eq!(
        first.atom_site_by_author("alpha").unwrap().atom(),
        top.atom_ids()[0]
    );
    assert!(matches!(
        first.atom_site_by_author("CA"),
        Err(HierarchyLookupError::NotFound { .. })
    ));
    assert!(matches!(
        second.atom_site_by_label("CA"),
        Err(HierarchyLookupError::Ambiguous { matches: 2, .. })
    ));
    let ca = AtomSelection::for_label_atom_names(&top, ["CA", "missing"]).unwrap();
    assert_eq!(
        ca.atom_ids().collect::<Vec<_>>(),
        [top.atom_ids()[0], top.atom_ids()[1], top.atom_ids()[3]]
    );
    let authors = AtomSelection::for_author_atom_names(&top, ["alpha"]).unwrap();
    assert_eq!(
        authors.atom_ids().collect::<Vec<_>>(),
        [top.atom_ids()[0], top.atom_ids()[1]]
    );
    assert!(AtomSelection::for_author_atom_names(&top, ["CA"])
        .unwrap()
        .indices()
        .is_empty());
    let amino = AtomSelection::for_residue_classes(&top, [ResidueClass::AminoAcid]).unwrap();
    assert_eq!(ca.intersection(&amino).unwrap(), ca);
}

#[test]
fn duplicate_residue_addresses_and_label_chains_are_not_silently_resolved() {
    let top = topology();
    let chain = top.chain_by_label("L").unwrap().id();
    let mut builder = Arc::try_unwrap(top).unwrap().into_builder();
    builder
        .hierarchy_mut()
        .add_residue(chain, "ALA", Some(1), Some("7".into()), None)
        .unwrap();
    builder.hierarchy_mut().add_chain("L", None).unwrap();
    let top = builder.build().unwrap();
    assert!(matches!(
        top.chain_by_label("L"),
        Err(HierarchyLookupError::Ambiguous { .. })
    ));
    assert!(matches!(
        top.chain(chain).unwrap().residue_by_label(1),
        Err(HierarchyLookupError::Ambiguous { .. })
    ));
    assert!(matches!(
        top.residue_by_author("A", "7", None),
        Err(HierarchyLookupError::Ambiguous { .. })
    ));
}

#[test]
fn original_1ake_pocket_uses_checked_composition_and_retains_source_associations() {
    let model = mmcif::parse_str(include_str!("fixtures/mmcif/1AKE.cif"))
        .unwrap()
        .interpret()
        .unwrap()
        .into_parts()
        .0;
    let top = model.shared_topology();
    let residue = top.residue_by_author("A", "215", None).unwrap();
    assert_eq!(residue.name(), "AP5");
    let ligand = AtomSelection::for_residues(&top, [residue.id()]).unwrap();
    assert_eq!(ligand.atom_ids().len(), 57);
    let hydrogen = AtomSelection::for_elements(&top, [Element::from_symbol("H").unwrap()]).unwrap();
    let candidates = AtomSelection::all(&top).difference(&hydrogen).unwrap();
    let nearby = measure::within(
        model.view(),
        &candidates,
        &ligand,
        Quantity::new(5.0, ANGSTROM),
    )
    .unwrap();
    let selected = nearby.expand_to_residues();
    let pocket = model.slice(&selected).unwrap();
    assert_eq!(pocket.atom_count(), 428);
    assert_eq!(pocket.residues().count(), 69);
    for (source, target) in selected.atom_ids().zip(pocket.atom_ids()) {
        assert_eq!(
            model.position(source).unwrap(),
            pocket.position(*target).unwrap()
        );
        assert_eq!(
            model
                .residue_for_atom(source)
                .unwrap()
                .unwrap()
                .author_seq_id(),
            pocket
                .residue_for_atom(*target)
                .unwrap()
                .unwrap()
                .author_seq_id()
        );
    }
}
