use super::*;

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
    let site = hierarchy
        .add_atom_site(residue, AtomId::new(7), metadata.clone())
        .unwrap();
    let second_chain = hierarchy.add_chain("B", None).unwrap();

    assert_eq!(
        hierarchy.chains().map(|(id, _)| id).collect::<Vec<_>>(),
        vec![chain, second_chain]
    );
    assert_eq!(hierarchy.chain(chain).unwrap().residues(), &[residue]);
    assert_eq!(hierarchy.residue(residue).unwrap().atom_sites(), &[site]);
    assert_eq!(
        hierarchy
            .atom_site_for_atom(AtomId::new(7))
            .unwrap()
            .metadata(),
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
    let atom = AtomId::new(2);
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
fn molecule_publication_validates_hierarchy_references() {
    let mut editor = MoleculeEditor::new();
    let atom = editor.add_atom(carbon()).unwrap();
    let chain = editor.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = editor
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), None, None)
        .unwrap();
    editor
        .add_atom_site(residue, atom, SmcraAtomSiteMetadata::default())
        .unwrap();
    let molecule = editor.finish().unwrap();
    assert_eq!(
        molecule
            .hierarchy()
            .atom_site_for_atom(atom)
            .unwrap()
            .atom(),
        atom
    );

    let mut edited = molecule.edit();
    edited.delete_atom(atom).unwrap();
    assert_eq!(edited.finish(), Err(MoleculePublicationError::EmptyGraph));
}

#[test]
fn deleting_annotated_atom_removes_site_without_renumbering_hierarchy() {
    let mut editor = MoleculeEditor::new();
    let annotated = editor.add_atom(carbon()).unwrap();
    let retained = editor.add_atom(carbon()).unwrap();
    editor
        .add_bond(annotated, retained, BondOrder::Single)
        .unwrap();
    let chain = editor.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = editor
        .hierarchy_mut()
        .add_residue(chain, "LIG", Some(1), None, None)
        .unwrap();
    let removed_site = editor
        .add_atom_site(residue, annotated, SmcraAtomSiteMetadata::default())
        .unwrap();
    let molecule = editor.finish().unwrap();

    let mut editor = molecule.edit();
    editor.delete_atom(annotated).unwrap();
    let molecule = editor.finish().expect("remaining atom is connected");

    assert_eq!(molecule.atom_count(), 1);
    assert!(molecule.atom(annotated).is_err());
    assert!(molecule.hierarchy().atom_site(removed_site).is_err());
    assert!(molecule.hierarchy().atom_site_for_atom(annotated).is_none());
    assert!(molecule
        .hierarchy()
        .residue(residue)
        .unwrap()
        .atom_sites()
        .is_empty());
    assert_eq!(
        molecule.hierarchy().chain(chain).unwrap().residues(),
        &[residue]
    );

    let mut editor = molecule.edit();
    let replacement_site = editor
        .add_atom_site(residue, retained, SmcraAtomSiteMetadata::default())
        .unwrap();
    assert!(replacement_site.raw() > removed_site.raw());
    let molecule = editor.finish().expect("replacement site publishes");
    assert_eq!(
        molecule
            .hierarchy()
            .atom_site_for_atom(retained)
            .unwrap()
            .id(),
        replacement_site
    );
}

#[test]
fn deleting_annotated_bridge_atom_still_fails_transactionally() {
    let mut editor = MoleculeEditor::new();
    let left = editor.add_atom(carbon()).unwrap();
    let bridge = editor.add_atom(carbon()).unwrap();
    let right = editor.add_atom(carbon()).unwrap();
    editor.add_bond(left, bridge, BondOrder::Single).unwrap();
    editor.add_bond(bridge, right, BondOrder::Single).unwrap();
    let chain = editor.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = editor
        .hierarchy_mut()
        .add_residue(chain, "LIG", Some(1), None, None)
        .unwrap();
    let site = editor
        .add_atom_site(residue, bridge, SmcraAtomSiteMetadata::default())
        .unwrap();
    let molecule = editor.finish().unwrap();
    let original = molecule.clone();

    let mut editor = molecule.edit();
    editor.delete_atom(bridge).unwrap();
    assert!(matches!(
        editor.finish(),
        Err(MoleculePublicationError::DisconnectedGraph(_))
    ));

    assert_eq!(molecule, original);
    assert_eq!(
        molecule
            .hierarchy()
            .atom_site_for_atom(bridge)
            .unwrap()
            .id(),
        site
    );
}
