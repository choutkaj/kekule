use kekule::core::{Atom, Element, MoleculeEditor};
use kekule::topology::{
    AtomSite, AtomSiteId, AtomSiteMetadata, AtomSiteView, Chain, ChainId, ChainView, Hierarchy,
    HierarchyError, HierarchyIdKind, InstanceAtomId, Residue, ResidueId, ResidueView,
    TopologyBuilder,
};

#[test]
fn canonical_topology_hierarchy_types_and_views_are_public() {
    let mut editor = MoleculeEditor::new();
    let atom = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let molecule = editor.finish().unwrap();

    let mut builder = TopologyBuilder::new();
    let instance = builder.add_molecule(&molecule).unwrap();
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
    assert_eq!(site_view.atom(), InstanceAtomId::new(instance, atom));

    assert_eq!(
        hierarchy.chain(ChainId::new(99)).unwrap_err(),
        HierarchyError::InvalidChainId(ChainId::new(99))
    );
    assert_eq!(HierarchyIdKind::AtomSite.to_string(), "atom-site");
}
