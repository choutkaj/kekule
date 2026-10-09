use kekule::core::{Atom, Element, MoleculeEditor};
use kekule::topology::{
    AtomSite, AtomSiteId, AtomSiteMetadata, AtomSiteView, Chain, ChainId, ChainView, Hierarchy,
    HierarchyError, HierarchyIdKind, InstanceAtomId, MoleculeClass, Residue, ResidueId,
    ResidueView, TopologyBuilder, TopologyEditor,
};

fn atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).unwrap())
}

#[test]
fn canonical_topology_hierarchy_types_and_views_are_public() {
    let mut editor = MoleculeEditor::new();
    let atom = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let molecule = editor.finish().unwrap();

    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(molecule.clone()).unwrap();
    let chain: ChainId = builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue: ResidueId = builder
        .hierarchy_mut()
        .add_residue(chain, "GLY", Some(1), None, None)
        .unwrap();
    let site: AtomSiteId = builder
        .hierarchy_mut()
        .add_atom_site(
            residue,
            InstanceAtomId::new(instance, atom),
            AtomSiteMetadata::default(),
        )
        .unwrap();
    let topology = builder.build().unwrap();

    let hierarchy: &Hierarchy = topology.hierarchy();
    let stored_chain: &Chain = hierarchy.chain(chain).unwrap();
    let stored_residue: &Residue = hierarchy.residue(residue).unwrap();
    let stored_site: &AtomSite = hierarchy.atom_site(site).unwrap();
    assert_eq!(stored_chain.residues(), &[residue]);
    assert_eq!(stored_residue.atom_sites(), &[site]);
    assert_eq!(stored_site.atom(), InstanceAtomId::new(instance, atom));

    let chain_view: ChainView<'_> = topology.chain(chain).unwrap();
    let residue_view: ResidueView<'_> = chain_view.residues().next().unwrap();
    let site_view: AtomSiteView<'_> = residue_view.atom_sites().next().unwrap();
    assert_eq!(site_view.atom().id(), InstanceAtomId::new(instance, atom));

    assert_eq!(
        hierarchy.chain(ChainId::new(99)).unwrap_err(),
        HierarchyError::InvalidChainId(ChainId::new(99))
    );
    assert_eq!(HierarchyIdKind::AtomSite.to_string(), "atom-site");
}

#[test]
fn editor_hierarchy_membership_stays_aligned_during_bulk_assembly_and_site_edits() {
    let mut editor = TopologyEditor::new();
    let chain = editor.add_chain("A", None).unwrap();
    let mut residues = Vec::new();
    let mut sites = Vec::new();
    for _ in 0..512 {
        let atom = editor.add_atom(atom("O")).unwrap();
        let residue = editor.add_residue(chain, "HOH", None, None, None).unwrap();
        sites.push(
            editor
                .add_atom_site(residue, atom, AtomSiteMetadata::default())
                .unwrap(),
        );
        residues.push(residue);
    }
    editor.set_atom_site_residue(sites[0], residues[1]).unwrap();
    assert!(editor.residue(residues[0]).is_err());
    editor.delete_atom_site(sites[1]).unwrap();
    assert!(editor.residue(residues[1]).is_ok());
    for &site in sites.iter().skip(2).step_by(2) {
        editor.delete_atom_site(site).unwrap();
    }
    editor.validate().unwrap();
    let mut cleared = editor.clone();
    let topology = editor.finish().unwrap();
    assert_eq!(topology.atom_count(), 512);
    assert_eq!(topology.hierarchy().atom_sites().count(), 256);
    assert_eq!(topology.residues().count(), 256);
    assert_eq!(
        topology
            .molecules()
            .filter(|m| m.class() == MoleculeClass::Water)
            .count(),
        256
    );
    let mut imported = topology.edit();
    let site = imported.atom_sites().next().unwrap().0;
    imported.delete_atom_site(site).unwrap();
    assert_eq!(
        imported.finish().unwrap().hierarchy().atom_sites().count(),
        255
    );
    cleared.clear();
    cleared.add_atom(atom("C")).unwrap();
    assert!(cleared.finish().unwrap().hierarchy().is_empty());
}
