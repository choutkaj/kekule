use kekule::molecule::{Atom, Element, MoleculeEditor};
use kekule::topology::{
    edit::TopologyEditor, AtomSiteId, AtomSiteMetadata, AtomSiteView, ChainId, ChainView,
    HierarchyIdKind, InstanceAtomId, MoleculeClass, ResidueId, ResidueView, TopologyBuilder,
};

fn atom(symbol: &str) -> Atom {
    Atom::new(Element::from_symbol(symbol).unwrap())
}

#[test]
fn hierarchy_views_read_what_the_builder_staged() {
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
    let staged = builder.hierarchy();
    assert_eq!(
        (
            staged.chain_count(),
            staged.residue_count(),
            staged.atom_site_count()
        ),
        (1, 1, 1)
    );
    let topology = builder.build().unwrap();

    let chain_view: ChainView<'_> = topology.chain(chain).unwrap();
    let residue_view: ResidueView<'_> = topology.residue(residue).unwrap();
    let site_view: AtomSiteView<'_> = topology.atom_site(site).unwrap();
    assert_eq!(
        chain_view.residues().map(|r| r.id()).collect::<Vec<_>>(),
        [residue]
    );
    assert_eq!(
        residue_view
            .atom_sites()
            .map(|s| s.id())
            .collect::<Vec<_>>(),
        [site]
    );
    assert_eq!(site_view.atom().id(), InstanceAtomId::new(instance, atom));
    assert_eq!(site_view.residue().chain().id(), chain);
    // Every hierarchy iterator knows its length.
    assert_eq!(topology.chains().len(), 1);
    assert_eq!(topology.residues().len(), 1);
    assert_eq!(topology.atom_sites().len(), 1);
    let occurrence = topology.instance(instance).unwrap();
    assert_eq!(
        (
            occurrence.chains().len(),
            occurrence.residues().len(),
            occurrence.atom_sites().len()
        ),
        (1, 1, 1)
    );

    assert!(topology.chain(ChainId::new(99)).is_none());
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
    assert!(editor.residue(residues[0]).is_none());
    editor.delete_atom_site(sites[1]).unwrap();
    assert!(editor.residue(residues[1]).is_some());
    for &site in sites.iter().skip(2).step_by(2) {
        editor.delete_atom_site(site).unwrap();
    }
    editor.validate().unwrap();
    let mut cleared = editor.clone();
    let topology = editor.finish().unwrap();
    assert_eq!(topology.atom_count(), 512);
    assert_eq!(topology.atom_sites().len(), 256);
    assert_eq!(topology.residues().count(), 256);
    assert_eq!(
        topology
            .instances()
            .filter(|m| m.class() == MoleculeClass::Water)
            .count(),
        256
    );
    let mut imported = topology.to_editor();
    let site = imported.atom_sites().next().unwrap().0;
    imported.delete_atom_site(site).unwrap();
    assert_eq!(imported.finish().unwrap().atom_sites().len(), 255);
    cleared.clear();
    cleared.add_atom(atom("C")).unwrap();
    let cleared = cleared.finish().unwrap();
    assert_eq!(
        (
            cleared.chains().len(),
            cleared.residues().len(),
            cleared.atom_sites().len()
        ),
        (0, 0, 0)
    );
}
