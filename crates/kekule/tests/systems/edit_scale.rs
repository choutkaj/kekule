//! Edits of large systems copy and perceive only what they change, and report
//! dense source-to-published maps that agree with the handle lookups.

use std::sync::Arc;

use kekule::geometry::Point3;
use kekule::molecule::{Atom, BondOrder, Element};
use kekule::structure::{Model, ModelEditor, Positions};
use kekule::topology::edit::{EditAtomId, EditBondId, EditCorrespondence, TopologyEditor};
use kekule::topology::{
    AtomSiteMetadata, ChainId, ResidueId, Topology, TopologyAtomIndex, TopologyBondIndex,
    TopologyBuilder,
};
use kekule::units::{Quantity, ANGSTROM};

/// Ethanol, benzene, two waters and a sodium ion, each a residue; the first
/// three in chain A, the rest in chain B.
fn system() -> Arc<Topology> {
    let mut builder = TopologyBuilder::new();
    let names = ["EOH", "BNZ", "HOH", "HOH", "NA"];
    let molecules = kekule::smiles::to_molecules("CCO.c1ccccc1.O.O.[Na+]").unwrap();
    let mut instances = Vec::new();
    for molecule in molecules {
        instances.push(builder.add_molecule(molecule).unwrap());
    }
    let chains = [
        builder.hierarchy_mut().add_chain("A", None).unwrap(),
        builder.hierarchy_mut().add_chain("B", None).unwrap(),
    ];
    let atoms = builder.atom_ids().to_vec();
    for (index, (&instance, name)) in instances.iter().zip(names).enumerate() {
        let chain = chains[usize::from(index >= 3)];
        let residue = builder
            .hierarchy_mut()
            .add_residue(chain, name, Some(index as i32 + 1), None, None)
            .unwrap();
        for &atom in atoms.iter().filter(|atom| atom.instance() == instance) {
            let metadata = AtomSiteMetadata::default()
                .with_label_atom_id(Some(format!("X{}", atom.atom().index())));
            builder
                .hierarchy_mut()
                .add_atom_site(residue, atom, metadata)
                .unwrap();
        }
    }
    Arc::new(builder.build().unwrap())
}

/// The draft handle of every source atom and bond, in dense order.
type Handles = (Vec<Option<EditAtomId>>, Vec<Option<EditBondId>>);

/// Checks the dense maps against the handle lookups: `handles` holds the draft
/// handle of every source atom and bond, taken before publication.
fn assert_dense_maps_agree(
    source: &Topology,
    published: &Topology,
    handles: &Handles,
    correspondence: &EditCorrespondence,
) {
    let atoms = correspondence.atoms();
    assert_eq!(atoms.to_published().len(), source.atom_count());
    assert_eq!(atoms.to_source().len(), published.atom_count());
    for (index, handle) in handles.0.iter().enumerate() {
        let source_index = TopologyAtomIndex::new(index as u32);
        let expected = handle
            .and_then(|handle| correspondence.atom(handle))
            .map(|(_, index)| index);
        assert_eq!(atoms.to_published()[index], expected);
        assert_eq!(atoms.published(source_index), expected);
        if let Some(target) = expected {
            assert_eq!(atoms.source(target), Some(source_index));
            let (id, _) = correspondence.atom(handle.unwrap()).unwrap();
            assert_eq!(published.atom_id(target), Some(id));
            assert_eq!(
                published.atom_at(target).unwrap().atom(),
                source.atom_at(source_index).unwrap().atom()
            );
        }
    }
    for (index, source_index) in atoms.to_source().iter().enumerate() {
        if let Some(source_index) = source_index {
            assert_eq!(
                atoms.to_published()[source_index.index()],
                Some(TopologyAtomIndex::new(index as u32))
            );
        }
    }
    let bonds = correspondence.bonds();
    assert_eq!(bonds.to_published().len(), source.bond_count());
    assert_eq!(bonds.to_source().len(), published.bond_count());
    for (index, handle) in handles.1.iter().enumerate() {
        let expected = handle
            .and_then(|handle| correspondence.bond(handle))
            .map(|(_, index)| index);
        assert_eq!(bonds.to_published()[index], expected);
        if let Some(target) = expected {
            assert_eq!(
                bonds.source(target),
                Some(TopologyBondIndex::new(index as u32))
            );
        }
    }
    for (index, source_index) in bonds.to_source().iter().enumerate() {
        if let Some(source_index) = source_index {
            assert_eq!(
                bonds.published(*source_index),
                Some(TopologyBondIndex::new(index as u32))
            );
        }
    }
    let residues = correspondence.residues();
    assert_eq!(residues.to_published().len(), source.residues().len());
    assert_eq!(residues.to_source().len(), published.residues().len());
    for (index, source_residue) in residues.to_source().iter().enumerate() {
        let residue = published.residue(ResidueId::new(index as u32)).unwrap();
        if let Some(source_residue) = source_residue {
            let before = source.residue(*source_residue).unwrap();
            assert_eq!(residue.name(), before.name());
            assert_eq!(residue.label_seq_id(), before.label_seq_id());
            assert_eq!(
                residues.published(*source_residue),
                Some(ResidueId::new(index as u32))
            );
        }
    }
    let chains = correspondence.chains();
    assert_eq!(chains.to_published().len(), source.chains().len());
    assert_eq!(chains.to_source().len(), published.chains().len());
    for (index, source_chain) in chains.to_source().iter().enumerate() {
        let chain = published.chain(ChainId::new(index as u32)).unwrap();
        if let Some(source_chain) = source_chain {
            assert_eq!(
                chain.label_id(),
                source.chain(*source_chain).unwrap().label_id()
            );
        }
    }
}

fn source_handles(source: &Topology, editor: &TopologyEditor) -> Handles {
    (
        source
            .atom_ids()
            .iter()
            .map(|&atom| editor.atom_handle(atom))
            .collect(),
        source
            .bond_ids()
            .iter()
            .map(|&bond| editor.bond_handle(bond))
            .collect(),
    )
}

#[test]
fn dense_maps_agree_with_handle_lookups_through_structural_edits() {
    let source = system();
    let atom = |index: usize| source.atom_ids()[index];
    type Edit = fn(&mut TopologyEditor, &dyn Fn(usize) -> kekule::topology::InstanceAtomId);
    let edits: [(&str, Edit); 6] = [
        ("delete a terminal atom", |editor, atom| {
            let handle = editor.atom_handle(atom(2)).unwrap();
            editor.delete_atom(handle).unwrap();
        }),
        (
            "delete a bridging atom, splitting ethanol",
            |editor, atom| {
                let handle = editor.atom_handle(atom(1)).unwrap();
                editor.delete_atom(handle).unwrap();
            },
        ),
        ("add an atom bonded to benzene", |editor, atom| {
            let ring = editor.atom_handle(atom(3)).unwrap();
            let added = editor
                .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
                .unwrap();
            editor.add_bond(ring, added, BondOrder::Single).unwrap();
        }),
        ("merge two waters and add a molecule", |editor, atom| {
            let first = editor.atom_handle(atom(9)).unwrap();
            let second = editor.atom_handle(atom(10)).unwrap();
            editor.add_bond(first, second, BondOrder::Single).unwrap();
            editor
                .add_molecule(kekule::smiles::to_molecules("CC").unwrap().remove(0))
                .unwrap();
        }),
        ("delete the ion, emptying its residue", |editor, atom| {
            let ion = editor.atom_handle(atom(11)).unwrap();
            editor.delete_atom(ion).unwrap();
        }),
        ("delete a whole chain's molecules", |editor, atom| {
            for index in 9..12 {
                let handle = editor.atom_handle(atom(index)).unwrap();
                editor.delete_atom(handle).unwrap();
            }
        }),
    ];
    for (name, edit) in edits {
        let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
        edit(&mut editor, &atom);
        let handles = source_handles(&source, &editor);
        let (published, correspondence) = editor.finish_with_correspondence().unwrap();
        assert!(!published.shares_layout(&source), "{name}");
        assert_dense_maps_agree(&source, &published, &handles, &correspondence);
    }
}

#[test]
fn dense_maps_name_deleted_residues_and_chains() {
    let source = system();
    let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
    for index in 10..12 {
        let handle = editor.atom_handle(source.atom_ids()[index]).unwrap();
        editor.delete_atom(handle).unwrap();
    }
    let handles = source_handles(&source, &editor);
    let (published, correspondence) = editor.finish_with_correspondence().unwrap();
    assert_dense_maps_agree(&source, &published, &handles, &correspondence);
    // Chain B held only the second water and the ion.
    let residues = correspondence.residues();
    assert_eq!(
        residues.to_published(),
        [
            Some(ResidueId::new(0)),
            Some(ResidueId::new(1)),
            Some(ResidueId::new(2)),
            None,
            None
        ]
    );
    assert_eq!(residues.to_source().len(), 3);
    assert_eq!(
        correspondence.chains().to_published(),
        [Some(ChainId::new(0)), None]
    );
}

#[test]
fn unedited_and_sourceless_drafts_report_identity_and_empty_maps() {
    let source = system();
    let editor = TopologyEditor::from_topology(Arc::clone(&source));
    let handles = source_handles(&source, &editor);
    let (published, correspondence) = editor.finish_with_correspondence().unwrap();
    assert!(Arc::ptr_eq(&published, &source));
    assert_dense_maps_agree(&source, &published, &handles, &correspondence);
    let identity = (0..source.atom_count() as u32)
        .map(|index| Some(TopologyAtomIndex::new(index)))
        .collect::<Vec<_>>();
    assert_eq!(correspondence.atoms().to_published(), identity);
    assert_eq!(correspondence.atoms().to_source(), identity);

    let mut editor = TopologyEditor::new();
    editor
        .add_molecule(kekule::smiles::to_molecules("CO").unwrap().remove(0))
        .unwrap();
    let (published, correspondence) = editor.finish_with_correspondence().unwrap();
    assert!(correspondence.atoms().to_published().is_empty());
    assert_eq!(correspondence.atoms().to_source(), [None, None]);
    assert_eq!(
        correspondence.bonds().to_source().len(),
        published.bond_count()
    );

    // A cleared draft keeps its source but none of its entities.
    let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
    editor.clear();
    editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    let (published, correspondence) = editor.finish_with_correspondence().unwrap();
    assert_eq!(published.atom_count(), 1);
    assert_eq!(
        correspondence.atoms().to_published(),
        vec![None; source.atom_count()]
    );
    assert_eq!(correspondence.atoms().to_source(), [None]);
}

#[test]
fn model_edits_remap_columns_through_the_dense_maps() {
    let source = system();
    let positions = Positions::new(Quantity::new(
        (0..source.atom_count())
            .map(|index| Point3::new(index as f64, 0.0, 0.0))
            .collect::<Vec<_>>(),
        ANGSTROM,
    ))
    .unwrap();
    let model = Model::new(Arc::clone(&source), positions).unwrap();
    let mut editor: ModelEditor = model.to_editor();
    let ring = editor.topology().atom_handle(source.atom_ids()[3]).unwrap();
    editor.delete_atom(ring).unwrap();
    let ethanol = editor.topology().atom_handle(source.atom_ids()[0]).unwrap();
    let added = editor
        .add_atom(
            Atom::new(Element::from_symbol("N").unwrap()),
            Quantity::new(Point3::new(-1.0, 0.0, 0.0), ANGSTROM),
        )
        .unwrap();
    editor.add_bond(ethanol, added, BondOrder::Single).unwrap();
    let (published, correspondence) = editor.finish_with_correspondence().unwrap();
    // Each published atom's position is its source atom's, read by index.
    for (index, source_index) in correspondence.atoms().to_source().iter().enumerate() {
        let position = published
            .position(published.topology().atom_ids()[index])
            .unwrap();
        let expected = match source_index {
            Some(source_index) => model.position(source.atom_ids()[source_index.index()]),
            None => Some(Quantity::new(Point3::new(-1.0, 0.0, 0.0), ANGSTROM)),
        }
        .unwrap();
        let (got, want) = (
            position.value_in(ANGSTROM).unwrap(),
            expected.value_in(ANGSTROM).unwrap(),
        );
        assert!((got - want).norm() < 1e-12);
    }
    assert_eq!(correspondence.atoms().to_published()[3], None);
}

#[test]
fn publication_shares_untouched_definitions_with_a_shared_source() {
    let source = system();
    let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
    let handle = editor.atom_handle(source.atom_ids()[2]).unwrap();
    editor.delete_atom(handle).unwrap();
    let published = editor.finish().unwrap();
    // The source is still shared, yet only the edited ethanol is new.
    let shared = |instance: usize| {
        std::ptr::eq(
            source.instances().nth(instance).unwrap().molecule(),
            published.instances().nth(instance).unwrap().molecule(),
        )
    };
    assert!(!shared(0));
    assert!((1..5).all(shared));
}

#[test]
fn perception_after_an_edit_perceives_only_the_edited_definitions() {
    let perceived = Arc::new(system().perceived().unwrap());
    // Perceiving a perceived snapshot again perceives nothing.
    let again = perceived.perceived().unwrap();
    assert!(perceived
        .definitions()
        .zip(again.definitions())
        .all(|(a, b)| std::ptr::eq(a.molecule(), b.molecule())));

    let mut editor = TopologyEditor::from_topology(Arc::clone(&perceived));
    let ring = editor.atom_handle(perceived.atom_ids()[3]).unwrap();
    let added = editor
        .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
        .unwrap();
    editor.add_bond(ring, added, BondOrder::Single).unwrap();
    let edited = editor.finish().unwrap();
    let benzene = |topology: &Topology| topology.instances().nth(1).unwrap().molecule().clone();
    assert!(benzene(&edited).perception().aromaticity().is_none());
    let reperceived = edited.perceived().unwrap();
    for (index, (before, after)) in edited.instances().zip(reperceived.instances()).enumerate() {
        assert_eq!(
            std::ptr::eq(before.molecule(), after.molecule()),
            index != 1,
            "instance {index}"
        );
    }
    // The result equals perceiving every definition from scratch.
    for definition in reperceived.definitions() {
        let mut expected = definition.molecule().clone();
        expected.clear_perception();
        expected.perceive().unwrap();
        assert_eq!(definition.molecule().perception(), expected.perception());
    }
    let mut toluene = benzene(&reperceived);
    toluene.clear_perception();
    toluene.perceive().unwrap();
    assert_eq!(benzene(&reperceived).perception(), toluene.perception());
}

#[test]
fn untouched_residues_keep_their_source_atoms_in_a_draft() {
    // Ethanol named C1, C2, O3 in one residue: placement names the next
    // carbon after its source atoms before any edit touches the residue.
    let mut builder = TopologyBuilder::new();
    let instance = builder
        .add_molecule(kekule::smiles::to_molecules("CCO").unwrap().remove(0))
        .unwrap();
    let chain = builder.hierarchy_mut().add_chain("A", None).unwrap();
    let residue = builder
        .hierarchy_mut()
        .add_residue(chain, "EOH", Some(1), None, None)
        .unwrap();
    let atoms = builder.atom_ids().to_vec();
    for (atom, name) in atoms
        .iter()
        .filter(|atom| atom.instance() == instance)
        .zip(["C1", "C2", "O3"])
    {
        let metadata = AtomSiteMetadata::default().with_label_atom_id(Some(name.to_owned()));
        builder
            .hierarchy_mut()
            .add_atom_site(residue, *atom, metadata)
            .unwrap();
    }
    let topology = Arc::new(builder.build().unwrap());
    let model = Model::new(Arc::clone(&topology), Positions::zeros(3)).unwrap();
    let editor = model.into_editor();
    let handle = editor.topology().residue_handle(residue).unwrap();
    let carbon = Element::from_symbol("C").unwrap();
    assert_eq!(editor.next_atom_name(handle, carbon).unwrap(), "C3");
    let oxygen = Element::from_symbol("O").unwrap();
    assert_eq!(editor.next_atom_name(handle, oxygen).unwrap(), "O4");
}
