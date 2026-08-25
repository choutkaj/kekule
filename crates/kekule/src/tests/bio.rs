use super::*;
use crate::bio::{Hierarchy, HierarchyError, SmcraAtomSiteMetadata, SmcraChainId};
use crate::topology::{InstanceAtomId, MoleculeInstanceId, TopologyBuildError, TopologyBuilder};

#[test]
fn hierarchy_adds_chains_residues_and_atom_sites_in_order() {
    let mut hierarchy = Hierarchy::new();
    let chain = hierarchy.add_chain("A", Some("authA".to_owned())).unwrap();
    let residue = hierarchy
        .add_residue(
            chain,
            "GLY",
            Some(10),
            Some("42".to_owned()),
            Some("A".to_owned()),
        )
        .unwrap();
    let metadata = SmcraAtomSiteMetadata {
        type_symbol: Some("C".to_owned()),
        label_asym_id: Some("A".to_owned()),
        auth_asym_id: Some("authA".to_owned()),
        label_atom_id: Some("CA".to_owned()),
        auth_atom_id: Some("CAY".to_owned()),
    };
    let atom = InstanceAtomId::new(MoleculeInstanceId::new(3), AtomId::new(7));
    let site = hierarchy
        .add_atom_site(residue, atom, metadata.clone())
        .unwrap();
    let second_chain = hierarchy.add_chain("B", None).unwrap();

    assert_eq!(
        hierarchy.chains().map(|(id, _)| id).collect::<Vec<_>>(),
        vec![chain, second_chain]
    );
    assert_eq!(hierarchy.chain(chain).unwrap().residues(), &[residue]);
    assert_eq!(hierarchy.residue(residue).unwrap().atom_sites(), &[site]);
    assert_eq!(
        hierarchy.atom_site_for_atom(atom).unwrap().metadata(),
        &metadata
    );
}

#[test]
fn hierarchy_rejects_missing_parents_and_duplicate_atom_placement() {
    let mut hierarchy = Hierarchy::new();
    let chain = hierarchy.add_chain("A", None).unwrap();
    assert_eq!(
        hierarchy
            .add_residue(SmcraChainId::new(99), "GLY", None, None, None)
            .unwrap_err(),
        HierarchyError::InvalidChainId(SmcraChainId::new(99))
    );
    let residue = hierarchy
        .add_residue(chain, "GLY", None, None, None)
        .unwrap();
    let atom = InstanceAtomId::new(MoleculeInstanceId::new(0), AtomId::new(2));
    hierarchy
        .add_atom_site(residue, atom, SmcraAtomSiteMetadata::default())
        .unwrap();
    assert_eq!(
        hierarchy
            .add_atom_site(residue, atom, SmcraAtomSiteMetadata::default())
            .unwrap_err(),
        HierarchyError::DuplicateAtomPlacement(atom)
    );
}

#[test]
fn topology_publication_validates_hierarchy_references() {
    let mut editor = MoleculeEditor::new();
    let atom = editor.add_atom(carbon()).unwrap();
    let molecule = editor.finish().unwrap();
    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(&molecule).unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), None, None)
        .unwrap();
    builder
        .hierarchy_mut()
        .add_atom_site(
            residue,
            InstanceAtomId::new(instance, AtomId::new(99)),
            SmcraAtomSiteMetadata::default(),
        )
        .unwrap();
    assert!(matches!(
        builder.build(),
        Err(TopologyBuildError::InvalidHierarchy(_))
    ));

    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(&molecule).unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), None, None)
        .unwrap();
    builder
        .hierarchy_mut()
        .add_atom_site(
            residue,
            InstanceAtomId::new(instance, atom),
            SmcraAtomSiteMetadata::default(),
        )
        .unwrap();
    assert_eq!(builder.build().unwrap().hierarchy().atom_sites().count(), 1);
}
