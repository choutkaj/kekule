use crate::molecule::{AtomId, MoleculeEditor};
use crate::test_fixtures::carbon;
use crate::topology::{AtomSiteId, AtomSiteMetadata, ChainId, HierarchyBuilder, HierarchyError};
use crate::topology::{InstanceAtomId, MoleculeInstanceId, TopologyBuildError, TopologyBuilder};
use crate::topology::{
    RawAtomSites, RawChain, RawHierarchyError, RawResidues, ResidueClass, ResidueId, ResidueLabels,
};

#[test]
fn hierarchy_adds_chains_residues_and_atom_sites_in_order() {
    let mut hierarchy = HierarchyBuilder::new();
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
    let metadata = AtomSiteMetadata::default()
        .with_type_symbol(Some("C".to_owned()))
        .with_label_asym_id(Some("A".to_owned()))
        .with_auth_asym_id(Some("authA".to_owned()))
        .with_label_atom_id(Some("CA".to_owned()))
        .with_auth_atom_id(Some("CAY".to_owned()));
    let atom = InstanceAtomId::new(MoleculeInstanceId::new(3), AtomId::new(7));
    let site = hierarchy
        .add_atom_site(residue, atom, metadata.clone())
        .unwrap();
    let second_chain = hierarchy.add_chain("B", None).unwrap();

    assert_eq!(
        hierarchy.storage().chain_ids().collect::<Vec<_>>(),
        vec![chain, second_chain]
    );
    assert_eq!(hierarchy.storage().chain_residues(chain), &[residue]);
    assert_eq!(hierarchy.storage().residue_sites(residue), &[site]);
    assert_eq!(hierarchy.atom_site_for_atom(atom), Some(site));
    assert_eq!(hierarchy.storage().site_metadata(site), &metadata);
}

#[test]
fn hierarchy_rejects_missing_parents_and_duplicate_atom_placement() {
    let mut hierarchy = HierarchyBuilder::new();
    let chain = hierarchy.add_chain("A", None).unwrap();
    assert_eq!(
        hierarchy
            .add_residue(ChainId::new(99), "GLY", None, None, None)
            .unwrap_err(),
        HierarchyError::InvalidChainId(ChainId::new(99))
    );
    let residue = hierarchy
        .add_residue(chain, "GLY", None, None, None)
        .unwrap();
    let atom = InstanceAtomId::new(MoleculeInstanceId::new(0), AtomId::new(2));
    hierarchy
        .add_atom_site(residue, atom, AtomSiteMetadata::default())
        .unwrap();
    assert_eq!(
        hierarchy
            .add_atom_site(residue, atom, AtomSiteMetadata::default())
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
    let instance = builder.add_molecule(molecule.clone()).unwrap();
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
            AtomSiteMetadata::default(),
        )
        .unwrap();
    assert!(matches!(
        builder.build(),
        Err(TopologyBuildError::MissingHierarchyAtom { .. })
    ));

    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(molecule.clone()).unwrap();
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
            AtomSiteMetadata::default(),
        )
        .unwrap();
    assert_eq!(builder.build().unwrap().atom_sites().count(), 1);
}

/// One chain of one residue of one site, with `residues` and `sites`
/// changed by the caller.
fn raw_parts(
    residues: impl FnOnce(&mut RawResidues),
    sites: impl FnOnce(&mut RawAtomSites),
) -> Result<HierarchyBuilder, RawHierarchyError> {
    let chains = vec![RawChain {
        label_id: "A".to_owned(),
        author_id: None,
        residues: vec![ResidueId::new(0)],
    }];
    let mut raw_residues = RawResidues {
        label_keys: vec![0],
        labels: vec![ResidueLabels::new("GLY".to_owned(), None, None, None)],
        classes: vec![ResidueClass::AminoAcid],
        explicit: vec![false],
        site_starts: vec![0, 1],
        sites: vec![AtomSiteId::new(0)],
    };
    let mut raw_sites = RawAtomSites {
        atoms: vec![InstanceAtomId::new(
            MoleculeInstanceId::new(0),
            AtomId::new(0),
        )],
        metadata_keys: vec![0],
        metadata: vec![AtomSiteMetadata::default()],
    };
    residues(&mut raw_residues);
    sites(&mut raw_sites);
    HierarchyBuilder::from_raw_parts(chains, raw_residues, raw_sites)
}

#[test]
fn raw_parts_reject_columns_that_do_not_match() {
    assert!(raw_parts(|_| {}, |_| {}).is_ok());
    let residue_cases: [fn(&mut RawResidues); 5] = [
        |residues| residues.label_keys[0] = 1,
        |residues| residues.explicit.clear(),
        |residues| residues.site_starts = vec![0, 2],
        |residues| residues.site_starts = vec![1, 1],
        |residues| residues.site_starts = vec![0],
    ];
    for (index, case) in residue_cases.into_iter().enumerate() {
        assert!(
            matches!(
                raw_parts(case, |_| {}),
                Err(RawHierarchyError::ResidueColumns)
            ),
            "residue case {index}"
        );
    }
    let site_cases: [fn(&mut RawAtomSites); 2] = [
        |sites| sites.metadata_keys[0] = 1,
        |sites| sites.metadata_keys.clear(),
    ];
    for (index, case) in site_cases.into_iter().enumerate() {
        assert!(
            matches!(
                raw_parts(|_| {}, case),
                Err(RawHierarchyError::AtomSiteColumns)
            ),
            "site case {index}"
        );
    }
}
