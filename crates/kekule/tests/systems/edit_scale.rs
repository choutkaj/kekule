//! Edits of large systems copy and perceive only what they change, and report
//! dense source-to-published maps that agree with the handle lookups.

use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;

use kekule::math::{Point3, Vector3};
use kekule::molecule::{Atom, BondOrder, Element};
use kekule::structure::placement::AttachmentSite;
use kekule::structure::{Positions, Structure, StructureEditor};
use kekule::topology::edit::{
    EditAtomId, EditBondId, EditCorrespondence, EditResidueId, TopologyEditor,
};
use kekule::topology::{
    AtomSiteMetadata, AtomView, ChainId, InstanceAtomId, ResidueId, Topology, TopologyAtomIndex,
    TopologyBondIndex, TopologyBuilder,
};
use kekule::units::{Quantity, ANGSTROM, SQUARE_ANGSTROM};

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
    assert!(published.shares_layout(&source));
    assert_dense_maps_agree(&source, &published, &handles, &correspondence);
    let identity = (0..source.atom_count() as u32)
        .map(|index| Some(TopologyAtomIndex::new(index)))
        .collect::<Vec<_>>();
    assert_eq!(correspondence.atoms().to_published(), identity);
    assert_eq!(correspondence.atoms().to_source(), identity);

    let mut editor = StructureEditor::new();
    editor
        .add_molecule(
            kekule::smiles::to_molecules("CO").unwrap().remove(0),
            &Positions::zeros(2),
        )
        .unwrap();
    let (published, correspondence) = editor.finish_with_correspondence().unwrap();
    assert!(correspondence.atoms().to_published().is_empty());
    assert_eq!(correspondence.atoms().to_source(), [None, None]);
    assert_eq!(
        correspondence.bonds().to_source().len(),
        published.topology().bond_count()
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
fn structure_edits_remap_columns_through_the_dense_maps() {
    let source = system();
    let positions = Positions::new(Quantity::new(
        (0..source.atom_count())
            .map(|index| Point3::new(index as f64, 0.0, 0.0))
            .collect::<Vec<_>>(),
        ANGSTROM,
    ))
    .unwrap();
    let structure = Structure::new(Arc::clone(&source), positions).unwrap();
    let mut editor: StructureEditor = structure.to_editor();
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
            Some(source_index) => structure.position(source.atom_ids()[source_index.index()]),
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
    let structure = Structure::new(Arc::clone(&topology), Positions::zeros(3)).unwrap();
    let editor = structure.into_editor();
    let handle = editor.topology().residue_handle(residue).unwrap();
    let carbon = Element::from_symbol("C").unwrap();
    assert_eq!(editor.next_atom_name(handle, carbon).unwrap(), "C3");
    let oxygen = Element::from_symbol("O").unwrap();
    assert_eq!(editor.next_atom_name(handle, oxygen).unwrap(), "O4");
}

/// A deep description of a structure read through its public views, so that a
/// change anywhere in its shared state shows.
fn describe(structure: &Structure) -> String {
    use std::fmt::Write;
    let topology = structure.topology();
    let view = structure.as_structure_view();
    let mut text = String::new();
    for atom in topology.atoms() {
        writeln!(
            text,
            "atom {} {:?} {:?} occupancy {:?} b {:?}",
            atom.id(),
            atom.atom(),
            structure
                .position(atom.id())
                .map(|p| p.value_in(ANGSTROM).unwrap()),
            view.occupancy(atom.index()),
            view.b_factor(atom.index())
                .map(|b| b.value_in(SQUARE_ANGSTROM).unwrap()),
        )
        .unwrap();
    }
    for bond in topology.bonds() {
        writeln!(text, "bond {:?} {:?}", bond.id(), bond.bond()).unwrap();
    }
    for instance in topology.instances() {
        writeln!(text, "instance {} {:?}", instance.id(), instance.class()).unwrap();
    }
    for chain in topology.chains() {
        writeln!(text, "chain {:?} {:?}", chain.label_id(), chain.author_id()).unwrap();
    }
    for residue in topology.residues() {
        writeln!(
            text,
            "residue {} {:?} {:?} {:?} {:?}",
            residue.name(),
            residue.class_assignment(),
            residue.label_comp_id(),
            residue.label_seq_id(),
            residue
                .atom_sites()
                .map(|site| site.id())
                .collect::<Vec<_>>(),
        )
        .unwrap();
    }
    for site in topology.atom_sites() {
        writeln!(text, "site {} {:?}", site.atom().id(), site.metadata()).unwrap();
    }
    text
}

#[test]
fn edits_through_shared_state_leave_the_source_structure_unchanged() {
    let source = system();
    let positions = Positions::new(Quantity::new(
        (0..source.atom_count())
            .map(|index| Point3::new(index as f64, 1.0, 2.0))
            .collect::<Vec<_>>(),
        ANGSTROM,
    ))
    .unwrap();
    let mut seeded = Structure::new(Arc::clone(&source), positions)
        .unwrap()
        .to_editor();
    for atom in source.atom_ids() {
        let handle = seeded.topology().atom_handle(*atom).unwrap();
        seeded.set_occupancy(handle, Some(0.5)).unwrap();
        seeded
            .set_b_factor(handle, Some(Quantity::new(20.0, SQUARE_ANGSTROM)))
            .unwrap();
    }
    let structure = seeded.finish().unwrap();
    let before = describe(&structure);
    let atom = |editor: &StructureEditor, index: usize| {
        editor
            .topology()
            .atom_handle(structure.topology().atom_ids()[index])
            .unwrap()
    };
    type Edit = fn(&mut StructureEditor, &dyn Fn(&StructureEditor, usize) -> EditAtomId);
    let edits: [(&str, Edit); 12] = [
        ("delete an atom", |editor, atom| {
            editor.delete_atom(atom(editor, 2)).unwrap();
        }),
        ("attach an atom", |editor, atom| {
            let anchor = atom(editor, 4);
            let added = editor
                .add_atom(
                    Atom::new(Element::from_symbol("C").unwrap()),
                    Quantity::new(Point3::new(9.0, 9.0, 9.0), ANGSTROM),
                )
                .unwrap();
            editor.add_bond(anchor, added, BondOrder::Single).unwrap();
        }),
        ("change an element", |editor, atom| {
            let handle = atom(editor, 0);
            let changed = editor
                .topology()
                .atom(handle)
                .unwrap()
                .clone()
                .with_element(Element::from_symbol("N").unwrap());
            editor.replace_atom(handle, changed).unwrap();
        }),
        ("change a charge", |editor, atom| {
            let handle = atom(editor, 9);
            let changed = editor
                .topology()
                .atom(handle)
                .unwrap()
                .clone()
                .with_formal_charge(-1);
            editor.replace_atom(handle, changed).unwrap();
        }),
        ("move an atom and change its annotations", |editor, atom| {
            let handle = atom(editor, 5);
            editor
                .set_position(handle, Quantity::new(Point3::new(-5.0, 0.0, 0.0), ANGSTROM))
                .unwrap();
            editor.set_occupancy(handle, None).unwrap();
            editor
                .set_b_factor(handle, Some(Quantity::new(1.0, SQUARE_ANGSTROM)))
                .unwrap();
        }),
        ("merge two molecules", |editor, atom| {
            let (ethanol, water) = (atom(editor, 2), atom(editor, 9));
            editor.add_bond(ethanol, water, BondOrder::Single).unwrap();
        }),
        ("split a molecule", |editor, _| {
            let bond = editor
                .topology()
                .bond_handle(editor.topology().source_topology().unwrap().bond_ids()[0])
                .unwrap();
            editor.delete_bond(bond).unwrap();
        }),
        ("rename a residue's component", |editor, _| {
            let residue = editor.topology().residue_handle(ResidueId::new(1)).unwrap();
            editor
                .set_residue_component_ids(residue, Some("TOL".into()), None)
                .unwrap();
        }),
        ("rename an atom site", |editor, atom| {
            let site = editor
                .topology()
                .atom_site_for_atom(atom(editor, 1))
                .unwrap()
                .unwrap();
            let metadata = AtomSiteMetadata::default().with_label_atom_id(Some("Q".into()));
            editor.set_atom_site_metadata(site, metadata).unwrap();
        }),
        ("delete a residue", |editor, _| {
            let residue = editor.topology().residue_handle(ResidueId::new(2)).unwrap();
            editor.delete_residue(residue).unwrap();
        }),
        ("delete a chain", |editor, _| {
            let chain = editor.topology().chain_handle(ChainId::new(1)).unwrap();
            editor.delete_chain(chain).unwrap();
        }),
        ("add hierarchy", |editor, atom| {
            let sodium = atom(editor, 11);
            let site = editor
                .topology()
                .atom_site_for_atom(sodium)
                .unwrap()
                .unwrap();
            editor.delete_atom_site(site).unwrap();
            let chain = editor.add_chain("C", None).unwrap();
            let residue = editor
                .add_residue(chain, "ION", Some(1), None, None)
                .unwrap();
            editor
                .add_atom_site(residue, sodium, AtomSiteMetadata::default())
                .unwrap();
        }),
    ];
    for (name, edit) in edits {
        let mut editor = structure.to_editor();
        edit(&mut editor, &atom);
        assert_eq!(
            describe(&structure),
            before,
            "{name}: the draft changed its source"
        );
        let published = editor.finish().unwrap();
        assert_ne!(
            describe(&published),
            before,
            "{name}: the edit changed nothing"
        );
        assert_eq!(
            describe(&structure),
            before,
            "{name}: publication changed its source"
        );
    }
}

#[test]
fn chained_edits_keep_every_snapshot_and_their_labels() {
    use kekule::topology::AtomSiteId;
    // Each edit starts from the previous publication, as an undo history
    // does; every snapshot keeps its own site names and residue labels.
    let sites = system().atom_count();
    let mut history = vec![Structure::new(system(), Positions::zeros(sites)).unwrap()];
    for round in 0..40 {
        let mut editor = history.last().unwrap().to_editor();
        let site = editor
            .topology()
            .atom_site_handle(AtomSiteId::new((round % sites) as u32))
            .unwrap();
        let residue = editor
            .topology()
            .residue_handle(ResidueId::new(round as u32 % 5))
            .unwrap();
        let metadata = AtomSiteMetadata::default().with_label_atom_id(Some(format!("R{round}")));
        editor.set_atom_site_metadata(site, metadata).unwrap();
        editor
            .set_residue_component_ids(residue, Some(format!("L{round}")), None)
            .unwrap();
        history.push(editor.finish().unwrap());
    }
    let site_name = |structure: &Structure, site: usize| {
        structure
            .topology()
            .atom_site(AtomSiteId::new(site as u32))
            .unwrap()
            .metadata()
            .label_atom_id()
            .map(str::to_owned)
    };
    for round in 0..40 {
        let (before, after) = (&history[round], &history[round + 1]);
        let site = round % sites;
        assert_eq!(site_name(after, site), Some(format!("R{round}")));
        let earlier = (round >= sites).then(|| format!("R{}", round - sites));
        assert_eq!(
            site_name(before, site),
            earlier.or_else(|| site_name(&history[0], site))
        );
        let residue = after
            .topology()
            .residue(ResidueId::new(round as u32 % 5))
            .unwrap();
        assert_eq!(residue.label_comp_id(), Some(format!("L{round}").as_str()));
    }
    assert_eq!(site_name(&history[0], 0).as_deref(), Some("X0"));
}

#[test]
fn appended_perceived_copies_keep_their_perception_profile() {
    // A host builds a large structure by appending perceived copies and never
    // perceives the result: perceiving it, or an edit of it, perceives only
    // what changed.
    let copy = Structure::new(
        Arc::new(system().perceived().unwrap()),
        Positions::zeros(12),
    )
    .unwrap();
    let mut merged = StructureEditor::new();
    for _ in 0..3 {
        merged.append_structure(&copy).unwrap();
    }
    // An added molecule has no known perception and is perceived.
    merged
        .add_molecule(
            kekule::smiles::to_molecules("C").unwrap().remove(0),
            &Positions::zeros(1),
        )
        .unwrap();
    let merged = merged.finish().unwrap();
    let shared = |before: &Topology, after: &Topology| {
        before
            .definitions()
            .zip(after.definitions())
            .map(|(a, b)| std::ptr::eq(a.molecule(), b.molecule()))
            .collect::<Vec<_>>()
    };
    let perceived = merged.topology().perceived().unwrap();
    let mut expected = vec![true; merged.topology().definition_count()];
    *expected.last_mut().unwrap() = false;
    assert_eq!(shared(merged.topology(), &perceived), expected);

    let mut editor = merged.to_editor();
    let ring = editor
        .topology()
        .atom_handle(merged.topology().atom_ids()[12 + 3])
        .unwrap();
    editor.delete_atom(ring).unwrap();
    let edited = editor.finish().unwrap();
    let reperceived = edited.topology().perceived().unwrap();
    let fresh = shared(edited.topology(), &reperceived)
        .iter()
        .filter(|shared| !**shared)
        .count();
    // The edited benzene's new definition and the added methane.
    assert_eq!(fresh, 2);
}

/// Deleting every glycine cuts the chains of adenylate kinase into fragments.
/// Each fragment publishes as the molecule the same deletion gives when the
/// rest of its chain goes too, which leaves the chain one molecule.
#[test]
fn split_fragments_equal_the_molecules_they_publish_as_alone() {
    use kekule::topology::MoleculeInstanceId;

    let structure = kekule::mmcif::parse_str(
        include_str!("../fixtures/mmcif/1AKE.cif"),
        Default::default(),
    )
    .unwrap()
    .interpret(Default::default())
    .unwrap()
    .into_structure();
    let topology = structure.topology();
    let glycines = topology
        .atoms()
        .filter(|atom| {
            atom.residue()
                .is_some_and(|residue| residue.name() == "GLY")
        })
        .map(|atom| atom.id())
        .collect::<BTreeSet<_>>();
    let deleted = |atoms: &BTreeSet<InstanceAtomId>| {
        let mut editor = structure.to_editor();
        let handles = topology
            .atom_ids()
            .iter()
            .map(|&atom| (atom, editor.topology().atom_handle(atom).unwrap()))
            .collect::<Vec<_>>();
        editor
            .delete_atoms(
                handles
                    .iter()
                    .filter(|(atom, _)| atoms.contains(atom))
                    .map(|&(_, handle)| handle),
            )
            .unwrap();
        let (edited, map) = editor.finish_with_correspondence().unwrap();
        // The source atoms of each published instance.
        let mut sources = BTreeMap::<MoleculeInstanceId, Vec<InstanceAtomId>>::new();
        for (atom, handle) in handles {
            if let Some((published, _)) = map.atom(handle) {
                sources.entry(published.instance()).or_default().push(atom);
            }
        }
        (edited, sources)
    };
    let (split, fragments) = deleted(&glycines);
    let mut compared = BTreeMap::<MoleculeInstanceId, usize>::new();
    for (&published, atoms) in &fragments {
        let source = atoms[0].instance();
        if !glycines.iter().any(|atom| atom.instance() == source) {
            continue;
        }
        let rest = topology
            .atom_ids()
            .iter()
            .copied()
            .filter(|atom| atom.instance() == source && !atoms.contains(atom));
        let (alone, instances) = deleted(&rest.chain(glycines.iter().copied()).collect());
        let (&instance, _) = instances
            .iter()
            .find(|(_, sources)| *sources == atoms)
            .expect("the fragment is one molecule alone");
        assert_eq!(
            split.topology().instance(published).unwrap().molecule(),
            alone.topology().instance(instance).unwrap().molecule(),
        );
        *compared.entry(source).or_default() += 1;
    }
    // Both chains split into many fragments.
    assert_eq!(compared.len(), 2, "{compared:?}");
    assert!(
        compared.values().all(|&fragments| fragments > 10),
        "{compared:?}"
    );
}

/// Deleting a sulfur and a hydrogen cuts a molecule in two. Each fragment
/// equals a molecule built by hand: the atoms with the hydrogens their cut
/// bonds give them, the bonds, both stereocenters, the one whose explicit
/// hydrogen went with that hydrogen as its implicit one, the stereo group
/// split between them, and the atom properties.
#[test]
fn split_fragments_equal_hand_built_molecules() {
    use kekule::molecule::{ImplicitHydrogens, Molecule, MoleculeEditor};
    use kekule::properties::{PropertyKey, PropertyValue};
    use kekule::stereo::{StereoCarrier, StereoElement, StereoElementKind, StereoGroup};
    use kekule::stereo::{StereoGroupKind, TetrahedralStereo};
    use std::collections::HashMap;

    // Atoms by SMILES position: H0 C1 F2 Cl3 C4 S5 C6 C7 H8 Br9 O10.
    let parsed = kekule::smiles::to_molecules("[H][C@](F)(Cl)[CH2]SC[C@@]([H])(Br)O")
        .unwrap()
        .remove(0);
    let atoms = parsed.atom_ids().collect::<Vec<_>>();
    assert_eq!(atoms.len(), 11);
    let hydrogens = |index: usize| parsed.atom(atoms[index]).unwrap().hydrogens();
    assert_eq!(hydrogens(4), ImplicitHydrogens::Fixed(2));
    assert_eq!(hydrogens(7), ImplicitHydrogens::Fixed(0));
    let tag = PropertyKey::new("tag").unwrap();
    let mut editor = parsed.clone().into_editor();
    for (index, &atom) in atoms.iter().enumerate() {
        editor
            .properties_mut()
            .atoms_mut()
            .set_value(tag.clone(), atom, Some(PropertyValue::Int(index as i64)))
            .unwrap();
    }
    let centers = editor.graph().stereo_element_ids().collect::<Vec<_>>();
    assert_eq!(centers.len(), 2);
    editor
        .add_stereo_group(StereoGroup {
            kind: StereoGroupKind::Absolute,
            members: centers,
        })
        .unwrap();
    let source: Molecule = editor.finish().unwrap();
    let index_of = |atom| atoms.iter().position(|&id| id == atom).unwrap();

    // The fragments by hand, atoms in source order.
    let expected = [&[0, 1, 2, 3, 4][..], &[6, 7, 9, 10]].map(|members| {
        let mut editor = MoleculeEditor::new();
        let mut local = HashMap::new();
        for &index in members {
            let atom = source.atom(atoms[index]).unwrap().clone();
            // C4 loses its bond to S5 and C7 its bond to H8: each fixed
            // count takes the valence lost.
            let atom = match index {
                4 => atom.with_hydrogens(ImplicitHydrogens::Fixed(3)),
                7 => atom.with_hydrogens(ImplicitHydrogens::Fixed(1)),
                _ => atom,
            };
            let added = editor.add_atom(atom).unwrap();
            editor
                .properties_mut()
                .atoms_mut()
                .set_value(tag.clone(), added, Some(PropertyValue::Int(index as i64)))
                .unwrap();
            local.insert(index, added);
        }
        for (_, bond) in source.bonds() {
            let ends = (
                local.get(&index_of(bond.a())),
                local.get(&index_of(bond.b())),
            );
            if let (Some(&a), Some(&b)) = ends {
                editor.add_bond(a, b, bond.order()).unwrap();
            }
        }
        let mut members = Vec::new();
        for (_, element) in source.stereo_elements() {
            let StereoElementKind::Tetrahedral(stereo) = &element.kind else {
                panic!("only tetrahedral centers");
            };
            let Some(&center) = local.get(&index_of(stereo.center)) else {
                continue;
            };
            let carriers = stereo.carriers.map(|carrier| match carrier {
                StereoCarrier::Atom(atom) if index_of(atom) == 8 => StereoCarrier::ImplicitHydrogen,
                StereoCarrier::Atom(atom) => StereoCarrier::Atom(local[&index_of(atom)]),
                other => other,
            });
            let kind = StereoElementKind::Tetrahedral(TetrahedralStereo {
                center,
                carriers,
                orientation: stereo.orientation,
            });
            members.push(editor.add_stereo_element(StereoElement::new(kind)).unwrap());
        }
        assert_eq!(members.len(), 1);
        editor
            .add_stereo_group(StereoGroup {
                kind: StereoGroupKind::Absolute,
                members,
            })
            .unwrap();
        editor.finish().unwrap()
    });

    let topology = Topology::from_molecules([source]).unwrap();
    let structure = Structure::new(topology, Positions::zeros(11)).unwrap();
    let mut editor = structure.to_editor();
    let deleted = [(5, "S"), (8, "H")].map(|(index, symbol)| {
        let handle = editor
            .topology()
            .atom_handle(structure.topology().atom_ids()[index])
            .unwrap();
        let atom = editor.topology().atom(handle).unwrap();
        assert_eq!(atom.element().symbol(), symbol);
        handle
    });
    editor.delete_atoms(deleted).unwrap();
    let split = editor.finish().unwrap();
    let published = split
        .topology()
        .instances()
        .map(|instance| instance.molecule().clone())
        .collect::<Vec<_>>();
    assert_eq!(published.len(), 2);
    for expected in &expected {
        let molecule = published
            .iter()
            .find(|molecule| molecule.atom_count() == expected.atom_count())
            .unwrap();
        assert_eq!(molecule, expected);
        assert_eq!(molecule.stereo_groups().count(), 1);
        assert_eq!(molecule.properties(), expected.properties());
    }
}

/// KRas G12C with GDP, Mg2+, its inhibitor and waters (PDB 9UOH),
/// perceived as a host perceives an import.
fn perceived_kras() -> Structure {
    let mut structure = kekule::mmcif::parse_str(
        include_str!("../fixtures/mmcif/9UOH.cif"),
        Default::default(),
    )
    .unwrap()
    .interpret(Default::default())
    .unwrap()
    .into_structure();
    structure.perceive().unwrap();
    structure
}

/// The publication order of a draft's atoms by the rule the editor
/// documents, from a scan of every live atom: handle order, except that an
/// added atom placed in a residue that holds earlier atoms follows the last
/// of them, in creation order.
fn documented_order(editor: &TopologyEditor) -> Vec<EditAtomId> {
    let source = editor
        .source_topology()
        .map_or_else(BTreeSet::new, |source| {
            source
                .atom_ids()
                .iter()
                .filter_map(|&atom| editor.atom_handle(atom))
                .collect()
        });
    let atoms = editor.atom_ids().collect::<Vec<_>>();
    let mut following = vec![Vec::new(); atoms.len()];
    let mut stays = vec![true; atoms.len()];
    let mut last = BTreeMap::<EditResidueId, usize>::new();
    for (position, &atom) in atoms.iter().enumerate() {
        let Some(site) = editor.atom_site_for_atom(atom).unwrap() else {
            continue;
        };
        let residue = editor.atom_site(site).unwrap().residue();
        match last.get(&residue) {
            Some(&anchor) if !source.contains(&atom) => {
                following[anchor].push(position);
                stays[position] = false;
            }
            _ => {
                last.insert(residue, position);
            }
        }
    }
    let mut order = Vec::with_capacity(atoms.len());
    for (position, &atom) in atoms.iter().enumerate() {
        if stays[position] {
            order.push(atom);
            order.extend(following[position].iter().map(|&behind| atoms[behind]));
        }
    }
    order
}

/// What a published atom should be, read from the draft by handle: the atom,
/// its position, and its chain, residue and site metadata.
fn drafted_atom(editor: &StructureEditor, atom: EditAtomId) -> String {
    let topology = editor.topology();
    let site = topology.atom_site_for_atom(atom).unwrap().map(|site| {
        let site = topology.atom_site(site).unwrap();
        let residue = topology.residue(site.residue()).unwrap();
        let chain = topology.chain(residue.chain()).unwrap();
        format!(
            "{} {} {:?} {:?}",
            chain.label_id(),
            residue.name(),
            residue.label_seq_id(),
            site.metadata()
        )
    });
    let position = editor.position(atom).unwrap().value_in(ANGSTROM).unwrap();
    format!("{:?} {position:?} {site:?}", topology.atom(atom).unwrap())
}

/// The same description of a published atom.
fn published_atom(structure: &Structure, index: usize) -> String {
    let topology = structure.topology();
    let atom = topology
        .atom_at(TopologyAtomIndex::new(index as u32))
        .unwrap();
    let site = atom.atom_site().map(|site| {
        let residue = site.residue();
        format!(
            "{} {} {:?} {:?}",
            residue.chain().label_id(),
            residue.name(),
            residue.label_seq_id(),
            site.metadata()
        )
    });
    let position = structure
        .position(atom.id())
        .unwrap()
        .value_in(ANGSTROM)
        .unwrap();
    format!("{:?} {position:?} {site:?}", atom.atom())
}

/// Publishes `editor`, checking the result against the draft: atoms in the
/// documented order, each with its draft atom, position and site, bonds
/// between the atoms their draft bonds join, dense maps that agree with the
/// handle lookups, and perception equal to perceiving from scratch.
fn publish_checked(editor: StructureEditor, edit: &str) -> Structure {
    let order = documented_order(editor.topology());
    let expected = order
        .iter()
        .map(|&atom| drafted_atom(&editor, atom))
        .collect::<Vec<_>>();
    let bonds = editor.topology().bonds().collect::<Vec<_>>();
    let handles = editor
        .topology()
        .source_topology()
        .map(|source| source_handles(source, editor.topology()));
    let (published, correspondence) = editor.finish_with_correspondence().unwrap();
    let index = |atom| correspondence.atom(atom).map(|(_, index)| index.index());
    let published_order = order.iter().map(|&atom| index(atom)).collect::<Vec<_>>();
    let dense = (0..order.len()).map(Some).collect::<Vec<_>>();
    assert_eq!(published_order, dense, "{edit}: the publication order");
    for (dense, expected) in expected.iter().enumerate() {
        assert_eq!(
            &published_atom(&published, dense),
            expected,
            "{edit}: atom {dense}"
        );
    }
    let topology = published.topology();
    assert_eq!(topology.bond_count(), bonds.len(), "{edit}: bonds");
    for (handle, bond) in bonds {
        let (_, published_bond) = correspondence.bond(handle).unwrap();
        let view = topology.bond_at(published_bond).unwrap();
        let mut ends = view.atoms().map(|atom| atom.index().index());
        let mut drafted = [index(bond.a()).unwrap(), index(bond.b()).unwrap()];
        ends.sort_unstable();
        drafted.sort_unstable();
        assert_eq!(ends, drafted, "{edit}: bond ends");
        assert_eq!(view.bond().order(), bond.order, "{edit}: bond order");
    }
    // The dense maps agree with the handle lookups and with each other.
    if let Some((atom_handles, bond_handles)) = handles {
        let atoms = correspondence.atoms();
        let by_handle = atom_handles
            .iter()
            .map(|handle| handle.and_then(|handle| correspondence.atom(handle)))
            .map(|found| found.map(|(_, index)| index))
            .collect::<Vec<_>>();
        assert_eq!(atoms.to_published(), by_handle, "{edit}: atom map");
        let bonds = correspondence.bonds();
        let by_handle = bond_handles
            .iter()
            .map(|handle| handle.and_then(|handle| correspondence.bond(handle)))
            .map(|found| found.map(|(_, index)| index))
            .collect::<Vec<_>>();
        assert_eq!(bonds.to_published(), by_handle, "{edit}: bond map");
        for (published, source) in atoms.to_source().iter().enumerate() {
            if let Some(source) = source {
                let back = atoms.to_published()[source.index()].map(|index| index.index());
                assert_eq!(back, Some(published), "{edit}: atom map back");
            }
        }
        let sources = atoms.to_source().iter().flatten().count();
        let survivors = atoms.to_published().iter().flatten().count();
        assert_eq!(sources, survivors, "{edit}: atom maps");
        assert_eq!(atoms.to_source().len(), topology.atom_count());
        assert_eq!(bonds.to_source().len(), topology.bond_count());
    }
    let perceived = topology.perceived().unwrap();
    for definition in perceived.definitions() {
        let mut fresh = definition.molecule().clone();
        fresh.clear_perception();
        fresh.perceive().unwrap();
        assert_eq!(
            definition.molecule().perception(),
            fresh.perception(),
            "{edit}: perception"
        );
    }
    published
}

/// Bonds a carbon to `parent` at its free position nearest a direction, as a
/// drag in a sketch does; the carbon joins the residue of `parent`.
fn attach_carbon(editor: &mut StructureEditor, parent: EditAtomId) -> EditAtomId {
    let direction = AttachmentSite::Toward(Vector3::new(1.0, 0.0, 0.0));
    let carbon = Atom::new(Element::from_symbol("C").unwrap());
    editor
        .attach_atom(parent, carbon, BondOrder::Single, direction)
        .unwrap()
        .atom()
}

/// Adds a carbon at `x` Å as a molecule and residue of its own.
fn place_carbon(editor: &mut StructureEditor, x: f64) -> EditAtomId {
    let carbon = Atom::new(Element::from_symbol("C").unwrap());
    let point = Quantity::new(Point3::new(x, 50.0, 50.0), ANGSTROM);
    editor.place_atom(carbon, point).unwrap()
}

/// The residue of an atom that has a site.
fn residue_of(editor: &StructureEditor, atom: EditAtomId) -> EditResidueId {
    let topology = editor.topology();
    let site = topology.atom_site_for_atom(atom).unwrap().unwrap();
    topology.atom_site(site).unwrap().residue()
}

const N_END: usize = 0;
const C_END: usize = 1;
const LIGAND_RING: usize = 2;
const PHE_RING: usize = 3;
const LIGAND_END: usize = 4;
const CYS_SG: usize = 5;
const WATER: usize = 6;
const CA_40: usize = 7;
const C_40: usize = 8;
const O_41: usize = 9;
const CB_41: usize = 10;
const N_END_CA: usize = 11;

/// The atoms the edits below start from, by the indices above: both ends of
/// chain A, ring atoms of the inhibitor and of a phenylalanine, the
/// inhibitor's end, a cysteine's sulfur, a water, atoms inside the chain,
/// and the first residue's alpha carbon.
fn kras_atoms(topology: &Topology) -> [InstanceAtomId; 12] {
    let find = |test: &dyn Fn(&AtomView<'_>) -> bool| {
        topology.atoms().find(|atom| test(atom)).unwrap().id()
    };
    let chain_a = topology
        .residues()
        .filter(|residue| residue.chain().label_id() == "A")
        .filter_map(|residue| residue.label_seq_id())
        .collect::<Vec<_>>();
    let site = |seq: i32, name: &str| {
        find(&|atom| {
            atom.atom_site().is_some_and(|site| {
                let residue = site.residue();
                residue.chain().label_id() == "A"
                    && residue.label_seq_id() == Some(seq)
                    && site.metadata().label_atom_id() == Some(name)
            })
        })
    };
    let named = |residue: &str, name: &str| {
        find(&|atom| {
            atom.atom_site().is_some_and(|site| {
                site.residue().name() == residue && site.metadata().label_atom_id() == Some(name)
            })
        })
    };
    let inhibitor = |test: &dyn Fn(&AtomView<'_>) -> bool| {
        find(&|atom| atom.residue().is_some_and(|r| r.name() == "A1L9E") && test(atom))
    };
    [
        site(chain_a[0], "N"),
        site(*chain_a.last().unwrap(), "CA"),
        inhibitor(&|atom| {
            atom.is_aromatic() == Some(true) && atom.implicit_hydrogens().unwrap_or(0) > 0
        }),
        named("PHE", "CZ"),
        inhibitor(&|atom| atom.neighbors().count() == 1),
        named("CYS", "SG"),
        named("HOH", "O"),
        site(40, "CA"),
        site(40, "C"),
        site(41, "O"),
        site(41, "CB"),
        site(chain_a[0], "CA"),
    ]
}

type KrasEdit = fn(&mut StructureEditor, &[EditAtomId; 12]);

/// Edits of every kind: moves, deletions, additions, element and charge
/// changes, bond changes, and attachments at chain ends, on rings and
/// across residues.
const KRAS_EDITS: [(&str, KrasEdit); 17] = [
    ("move an atom", |editor, atoms| {
        let point = Quantity::new(Point3::new(1.0, 2.0, 3.0), ANGSTROM);
        editor.set_position(atoms[PHE_RING], point).unwrap();
    }),
    ("delete the inhibitor's end", |editor, atoms| {
        editor.delete_atom(atoms[LIGAND_END]).unwrap();
    }),
    ("delete an atom inside a chain", |editor, atoms| {
        editor.delete_atom(atoms[CA_40]).unwrap();
    }),
    ("add an atom without a site", |editor, atoms| {
        let point = Quantity::new(Point3::new(9.0, 9.0, 9.0), ANGSTROM);
        let carbon = Atom::new(Element::from_symbol("C").unwrap());
        let added = editor.add_atom(carbon, point).unwrap();
        editor
            .add_bond(atoms[LIGAND_RING], added, BondOrder::Single)
            .unwrap();
    }),
    ("change an element", |editor, atoms| {
        let atom = editor.topology().atom(atoms[PHE_RING]).unwrap().clone();
        let nitrogen = atom.with_element(Element::from_symbol("N").unwrap());
        editor.replace_atom(atoms[PHE_RING], nitrogen).unwrap();
    }),
    ("change a charge", |editor, atoms| {
        let atom = editor.topology().atom(atoms[WATER]).unwrap().clone();
        editor
            .replace_atom(atoms[WATER], atom.with_formal_charge(-1))
            .unwrap();
    }),
    ("bond the inhibitor to a cysteine", |editor, atoms| {
        editor
            .add_bond(atoms[LIGAND_END], atoms[CYS_SG], BondOrder::Single)
            .unwrap();
    }),
    ("cut a chain and change a bond order", |editor, atoms| {
        let topology = editor.topology();
        let cut = topology.bond_between(atoms[CA_40], atoms[C_40]).unwrap();
        let carbonyl = topology.incident_bonds(atoms[O_41]).unwrap().next();
        let carbonyl = carbonyl.unwrap();
        editor.delete_bond(cut.unwrap()).unwrap();
        editor.set_bond_order(carbonyl, BondOrder::Single).unwrap();
    }),
    ("attach at both ends of a chain", |editor, atoms| {
        attach_carbon(editor, atoms[N_END]);
        attach_carbon(editor, atoms[C_END]);
    }),
    ("attach to rings", |editor, atoms| {
        attach_carbon(editor, atoms[LIGAND_RING]);
        attach_carbon(editor, atoms[PHE_RING]);
    }),
    ("attach across residues in turn", |editor, atoms| {
        // Residues gain atoms out of creation order, and a new residue's
        // atom anchors the atoms attached to it.
        let added = attach_carbon(editor, atoms[N_END]);
        attach_carbon(editor, atoms[WATER]);
        attach_carbon(editor, added);
        let placed = place_carbon(editor, 30.0);
        attach_carbon(editor, atoms[CB_41]);
        attach_carbon(editor, placed);
        attach_carbon(editor, atoms[N_END]);
        attach_carbon(editor, placed);
    }),
    (
        "attach after moving an atom to an earlier residue",
        |editor, atoms| {
            // The moved atom is now the earlier residue's last source atom, and
            // the residue loses its own last atom: added atoms follow the moved
            // atom.
            let residue = residue_of(editor, atoms[CA_40]);
            let site = editor.topology().atom_site_for_atom(atoms[CB_41]).unwrap();
            editor
                .set_atom_site_residue(site.unwrap(), residue)
                .unwrap();
            let own = editor
                .topology()
                .atom_ids()
                .filter(|&atom| atom != atoms[CB_41])
                .filter(|&atom| {
                    let site = editor.topology().atom_site_for_atom(atom).unwrap();
                    site.is_some() && residue_of(editor, atom) == residue
                })
                .last()
                .unwrap();
            editor.delete_atom(own).unwrap();
            attach_carbon(editor, atoms[CA_40]);
            attach_carbon(editor, atoms[CB_41]);
        },
    ),
    (
        "attach to a placed atom, then delete it",
        |editor, atoms| {
            // The first atom attached to the placed one then begins the new
            // residue.
            let placed = place_carbon(editor, 30.0);
            attach_carbon(editor, placed);
            attach_carbon(editor, atoms[LIGAND_RING]);
            attach_carbon(editor, placed);
            editor.delete_atom(placed).unwrap();
        },
    ),
    (
        "attach, then move the atom to another residue",
        |editor, atoms| {
            let added = attach_carbon(editor, atoms[N_END]);
            attach_carbon(editor, atoms[WATER]);
            let residue = residue_of(editor, atoms[WATER]);
            let site = editor.topology().atom_site_for_atom(added).unwrap();
            editor
                .set_atom_site_residue(site.unwrap(), residue)
                .unwrap();
        },
    ),
    ("attach after deleting the first atom", |editor, atoms| {
        editor.delete_atom(atoms[N_END]).unwrap();
        attach_carbon(editor, atoms[N_END_CA]);
    }),
    ("attach and settle the hydrogens", |editor, atoms| {
        let added = attach_carbon(editor, atoms[LIGAND_RING]);
        editor.settle_atoms([added, atoms[LIGAND_RING]]).unwrap();
    }),
    ("attach after clearing the draft", |editor, _| {
        editor.clear();
        let (first, second) = (place_carbon(editor, 0.0), place_carbon(editor, 9.0));
        attach_carbon(editor, first);
        attach_carbon(editor, second);
        attach_carbon(editor, first);
    }),
];

#[test]
fn publication_keeps_the_documented_atom_order_through_every_kind_of_edit() {
    let kras = perceived_kras();
    let atoms = kras_atoms(kras.topology());
    for (name, edit) in KRAS_EDITS {
        let mut editor = kras.to_editor();
        let handles = atoms.map(|atom| editor.topology().atom_handle(atom).unwrap());
        edit(&mut editor, &handles);
        // Each attachment moves atoms out of handle order.
        let handle_order = editor.topology().atom_ids().collect::<Vec<_>>();
        let reordered = documented_order(editor.topology()) != handle_order;
        assert_eq!(reordered, name.starts_with("attach"), "{name}: reordered");
        let published = publish_checked(editor, name);
        // The next edit starts from the published order.
        let mut again = published.to_editor();
        let first = again.topology().atom_ids().next().unwrap();
        again.delete_atom(first).unwrap();
        let next = again.topology().atom_ids().next().unwrap();
        attach_carbon(&mut again, next);
        publish_checked(again, name);
    }
}

#[test]
fn a_draft_without_a_source_publishes_in_the_documented_order() {
    let mut editor = StructureEditor::new();
    editor.append_structure(&perceived_kras()).unwrap();
    let placed = [
        place_carbon(&mut editor, 0.0),
        place_carbon(&mut editor, 9.0),
    ];
    let first = editor.topology().atom_ids().next().unwrap();
    for parent in [placed[0], first, placed[1], placed[0]] {
        attach_carbon(&mut editor, parent);
    }
    publish_checked(editor, "a draft without a source");
}
