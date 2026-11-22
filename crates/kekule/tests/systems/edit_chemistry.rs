//! Edits keep represented chemistry: hydrogen counts follow lost and gained
//! bond valence, and stereo carriers follow deleted hydrogens.
use kekule::math::{Point3, Vector3};
use kekule::molecule::{
    Atom, AtomId, AtomRadical, BondOrder, Element, ImplicitHydrogens, Molecule,
};
use kekule::perception::valence::open_valences;
use kekule::structure::placement::AttachmentSite;
use kekule::structure::StructureEditor;
use kekule::topology::edit::{EditAtomId, TopologyEditor};
use kekule::topology::{AtomSelection, Topology, TopologyAtomIndex, TopologyBuilder};
use kekule::units::{Quantity, NANOMETER};
use std::sync::Arc;

fn molecules(text: &str) -> Vec<Molecule> {
    kekule::smiles::to_molecules(text).unwrap()
}

fn topology(text: &str) -> Arc<Topology> {
    let mut builder = TopologyBuilder::new();
    for molecule in molecules(text) {
        builder.add_molecule(molecule).unwrap();
    }
    Arc::new(builder.build().unwrap())
}

fn handle(editor: &TopologyEditor, index: usize) -> EditAtomId {
    editor.atom_ids().nth(index).unwrap()
}

fn hydrogens(topology: &Topology, index: u32) -> ImplicitHydrogens {
    topology
        .atom_at(TopologyAtomIndex::new(index))
        .unwrap()
        .atom()
        .hydrogens()
}

#[test]
fn deleting_a_substituent_gives_its_fixed_neighbor_a_hydrogen() {
    // Toluene with every count fixed, as a structure file states it: the
    // ring carbon that loses the methyl becomes a CH, so the result is benzene.
    let source = topology("[CH3][C]1=[CH][CH]=[CH][CH]=[CH]1");
    let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
    editor.delete_atom(handle(&editor, 0)).unwrap();
    let edited = editor.finish().unwrap().perceived().unwrap();
    assert_eq!(hydrogens(&edited, 0), ImplicitHydrogens::Fixed(1));
    for index in 0..6 {
        let atom = edited.atom_at(TopologyAtomIndex::new(index)).unwrap();
        assert_eq!(atom.total_hydrogens(), Some(1));
    }
}

#[test]
fn deleting_an_explicit_hydrogen_keeps_the_count_of_a_fixed_parent() {
    // Water from an XYZ file fixes the oxygen at zero implicit hydrogens.
    let xyz = "3\nwater\nO 0.000 0.000 0.000\nH 0.957 0.000 0.000\nH -0.240 0.927 0.000\n";
    let structure = kekule::xyz::parse_str(xyz, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
        .into_structure();
    assert_eq!(
        hydrogens(structure.topology(), 0),
        ImplicitHydrogens::Fixed(0)
    );
    let mut editor = structure.to_editor();
    let hydrogen = editor.topology().atom_ids().nth(1).unwrap();
    editor.delete_atom(hydrogen).unwrap();
    let edited = editor.finish().unwrap();
    let oxygen = edited
        .topology()
        .atom_at(TopologyAtomIndex::new(0))
        .unwrap();
    assert_eq!(oxygen.atom().hydrogens(), ImplicitHydrogens::Fixed(1));
    assert_eq!(oxygen.total_hydrogens(), Some(2));
}

#[test]
fn lowering_an_order_and_deleting_a_bond_raise_fixed_counts() {
    let source = topology("[CH2]=[CH2]");
    let mut editor = TopologyEditor::from_topology(source);
    let bond = editor.bond_ids().next().unwrap();
    editor.set_bond_order(bond, BondOrder::Single).unwrap();
    assert_eq!(
        editor.atom(handle(&editor, 0)).unwrap().hydrogens(),
        ImplicitHydrogens::Fixed(3)
    );
    editor.delete_bond(bond).unwrap();
    let edited = editor.finish().unwrap();
    assert_eq!(edited.instance_count(), 2);
    assert_eq!(hydrogens(&edited, 0), ImplicitHydrogens::Fixed(4));
    assert_eq!(hydrogens(&edited, 1), ImplicitHydrogens::Fixed(4));
}

#[test]
fn adding_a_bond_takes_hydrogens_from_fixed_counts() {
    let source = topology("[CH4].[CH4]");
    let mut editor = TopologyEditor::from_topology(source);
    let (a, b) = (handle(&editor, 0), handle(&editor, 1));
    let bond = editor.add_bond(a, b, BondOrder::Single).unwrap();
    assert_eq!(
        editor.atom(a).unwrap().hydrogens(),
        ImplicitHydrogens::Fixed(3)
    );
    editor.set_bond_order(bond, BondOrder::Double).unwrap();
    assert_eq!(
        editor.atom(b).unwrap().hydrogens(),
        ImplicitHydrogens::Fixed(2)
    );
    // Deleting the bond undoes both edits.
    editor.delete_bond(bond).unwrap();
    assert_eq!(
        editor.atom(a).unwrap().hydrogens(),
        ImplicitHydrogens::Fixed(4)
    );
    let edited = editor.finish().unwrap();
    assert_eq!(hydrogens(&edited, 1), ImplicitHydrogens::Fixed(4));
}

#[test]
fn rewiring_a_bond_moves_the_hydrogen_it_took() {
    let source = topology("[CH3][CH3].[CH4]");
    let mut editor = TopologyEditor::from_topology(source);
    let bond = editor.bond_ids().next().unwrap();
    let (a, c) = (handle(&editor, 0), handle(&editor, 2));
    editor.set_bond_endpoints(bond, a, c).unwrap();
    let edited = editor.finish().unwrap();
    assert_eq!(hydrogens(&edited, 0), ImplicitHydrogens::Fixed(3));
    assert_eq!(hydrogens(&edited, 1), ImplicitHydrogens::Fixed(4));
    assert_eq!(hydrogens(&edited, 2), ImplicitHydrogens::Fixed(3));
}

#[test]
fn a_subset_caps_fixed_atoms_at_its_cuts() {
    let source = topology("[CH3][CH2][OH]");
    let selection =
        AtomSelection::from_indices(Arc::clone(&source), [0, 1].map(TopologyAtomIndex::new))
            .unwrap();
    let subset = Topology::subset(source, &selection).unwrap();
    assert_eq!(hydrogens(subset.topology(), 1), ImplicitHydrogens::Fixed(3));
}

#[test]
fn open_valences_name_fixed_atoms_short_of_their_valence() {
    // A molecule draft edits the graph as stated, so the carbon that loses
    // a neighbor keeps two hydrogens and an open valence.
    let mut editor = molecules("[CH3][CH2][CH3]").pop().unwrap().to_editor();
    editor.delete_atom(AtomId::new(2)).unwrap();
    let radical = editor.finish().unwrap();
    let open = open_valences(&radical);
    assert_eq!(open.len(), 1);
    assert_eq!(open[0].atom(), AtomId::new(1));
    assert_eq!(open[0].missing(), 1);
    // A declared radical accounts for the missing valence.
    let declared = molecules("[CH3][CH2]").pop().unwrap();
    assert_eq!(open_valences(&declared), Vec::new());
    // Inferred counts never leave an open valence.
    assert!(open_valences(&molecules("CC").pop().unwrap()).is_empty());
}

fn canonical(molecule: &Molecule) -> String {
    kekule::smiles::write(molecule, kekule::smiles::SmilesWriteOptions::canonical()).unwrap()
}

/// The molecule of `text` with explicit hydrogen `hydrogen` deleted by a
/// topology edit, and the same molecule with its hydrogens collapsed.
fn deleted_and_collapsed(text: &str, hydrogen: usize) -> (Molecule, Molecule) {
    let source = topology(text);
    let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
    editor.delete_atom(handle(&editor, hydrogen)).unwrap();
    let edited = editor.finish().unwrap();
    let deleted = edited.instances().next().unwrap().molecule().clone();
    let mut collapsed = molecules(text).pop().unwrap();
    collapsed.perceive().unwrap();
    collapsed.remove_hydrogens().unwrap();
    (deleted, collapsed)
}

#[test]
fn deleting_a_hydrogen_carrier_keeps_tetrahedral_stereo() {
    let (deleted, collapsed) = deleted_and_collapsed("[H][C@](F)(Cl)Br", 0);
    assert_eq!(deleted.stereo_elements().count(), 1);
    assert_eq!(canonical(&deleted), canonical(&collapsed));
}

#[test]
fn deleting_a_hydrogen_carrier_keeps_double_bond_stereo() {
    let (deleted, collapsed) = deleted_and_collapsed("F/C=C(/[H])Cl", 3);
    assert_eq!(deleted.stereo_elements().count(), 1);
    assert_eq!(canonical(&deleted), canonical(&collapsed));
}

#[test]
fn closing_a_ring_at_a_center_puts_the_new_neighbor_where_its_hydrogen_was() {
    let source = topology("F[C@@H](Cl)CCC");
    let mut editor = TopologyEditor::from_topology(source);
    let (center, end) = (handle(&editor, 1), handle(&editor, 5));
    editor.add_bond(center, end, BondOrder::Single).unwrap();
    let edited = editor.finish().unwrap();
    let mut ring = edited.instances().next().unwrap().molecule().clone();
    assert_eq!(ring.stereo_elements().count(), 1);
    ring.perceive().unwrap();
    let mut expected = molecules("F[C@@]1(Cl)CCC1").pop().unwrap();
    expected.perceive().unwrap();
    assert_eq!(canonical(&ring), canonical(&expected));
}

/// Cisplatin with every count fixed: Pt 0, Cl 0, N 3.
const CISPLATIN: &str = "[Cl][Pt]([Cl])([NH3])[NH3]";

fn fixed_counts(topology: &Topology) -> Vec<ImplicitHydrogens> {
    (0..topology.atom_count() as u32)
        .map(|index| hydrogens(topology, index))
        .collect()
}

#[test]
fn only_donor_bonds_to_a_metal_keep_hydrogen_counts() {
    // The ammines donate to platinum, so their bonds change no count. The
    // chlorides' bonds are covalent: a chloride cut off takes a hydrogen.
    let source = topology(CISPLATIN);
    let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
    let bonds = editor.bonds().collect::<Vec<_>>();
    let ending = |index| {
        bonds
            .iter()
            .find(|(_, bond)| {
                bond.a() == handle(&editor, index) || bond.b() == handle(&editor, index)
            })
            .unwrap()
            .0
    };
    let (chloride, amine) = (ending(0), ending(3));
    editor.delete_bond(chloride).unwrap();
    editor.set_bond_order(amine, BondOrder::Dative).unwrap();
    let edited = editor.finish().unwrap();
    assert_eq!(hydrogens(&edited, 0), ImplicitHydrogens::Fixed(1));
    assert_eq!(hydrogens(&edited, 1), ImplicitHydrogens::Fixed(0));
    assert_eq!(hydrogens(&edited, 3), ImplicitHydrogens::Fixed(3));
    let ligands = [
        ImplicitHydrogens::Fixed(1),
        ImplicitHydrogens::Fixed(1),
        ImplicitHydrogens::Fixed(3),
        ImplicitHydrogens::Fixed(3),
    ];
    let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
    editor.delete_atom(handle(&editor, 1)).unwrap();
    assert_eq!(fixed_counts(&editor.finish().unwrap()), ligands);
    let without_metal = AtomSelection::from_predicate(Arc::clone(&source), |atom| {
        atom.atom().element().symbol() != "Pt"
    });
    let subset = Topology::subset(Arc::clone(&source), &without_metal).unwrap();
    assert_eq!(fixed_counts(subset.topology()), ligands);
}

#[test]
fn a_carbon_cut_from_lithium_takes_a_hydrogen() {
    let source = topology("[CH3][Li]");
    let mut editor = TopologyEditor::from_topology(source);
    editor.delete_atom(handle(&editor, 1)).unwrap();
    let edited = editor.finish().unwrap();
    assert_eq!(fixed_counts(&edited), [ImplicitHydrogens::Fixed(4)]);
}

#[test]
fn a_kekule_swap_keeps_every_count_in_any_order() {
    // Doubling the single bonds first leaves the substituted carbon owing a
    // hydrogen until a double bond becomes single.
    let source = topology("[CH3][C]1=[CH][CH]=[CH][CH]=[CH]1");
    let before = fixed_counts(&source);
    let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
    let ring = editor
        .bonds()
        .filter(|(_, bond)| bond.a() != handle(&editor, 0) && bond.b() != handle(&editor, 0))
        .map(|(id, bond)| (id, bond.order))
        .collect::<Vec<_>>();
    assert_eq!(ring.len(), 6);
    for &(id, order) in &ring {
        if order == BondOrder::Single {
            editor.set_bond_order(id, BondOrder::Double).unwrap();
        }
    }
    for &(id, order) in &ring {
        if order == BondOrder::Double {
            editor.set_bond_order(id, BondOrder::Single).unwrap();
        }
    }
    let edited = editor.finish().unwrap();
    assert_eq!(fixed_counts(&edited), before);
}

#[test]
fn a_bond_past_the_last_hydrogen_fails_publication() {
    let source = topology("[O]=[C]=[O]");
    let mut editor = TopologyEditor::from_topology(source);
    let (left, right) = (handle(&editor, 0), handle(&editor, 2));
    editor.add_bond(left, right, BondOrder::Single).unwrap();
    let error = editor.finish().unwrap_err();
    assert!(error.to_string().contains("bond valence"), "{error}");
}

#[test]
fn a_radical_neighbor_takes_the_hydrogen_its_valence_has_room_for() {
    // Propyl radical without its middle carbon: methane and a methyl
    // radical.
    let source = topology("[CH3][CH2][CH2]");
    let radical = source.atom_at(TopologyAtomIndex::new(2)).unwrap();
    assert!(radical.atom().radical().is_some());
    let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
    editor.delete_atom(handle(&editor, 1)).unwrap();
    let edited = editor.finish().unwrap();
    assert_eq!(hydrogens(&edited, 0), ImplicitHydrogens::Fixed(4));
    assert_eq!(hydrogens(&edited, 1), ImplicitHydrogens::Fixed(3));
    let methyl = edited.atom_at(TopologyAtomIndex::new(1)).unwrap();
    assert_eq!(methyl.atom().radical(), AtomRadical::new(1, None));
}

#[test]
fn a_bond_between_radicals_takes_their_electrons() {
    let source = topology("[CH3].[CH3]");
    let mut editor = TopologyEditor::from_topology(source);
    let (left, right) = (handle(&editor, 0), handle(&editor, 1));
    editor.add_bond(left, right, BondOrder::Single).unwrap();
    let edited = editor.finish().unwrap();
    assert_eq!(fixed_counts(&edited), [ImplicitHydrogens::Fixed(3); 2]);
    assert!(edited.atoms().all(|atom| atom.atom().radical().is_none()));
}

#[test]
fn a_declared_radical_keeps_its_hydrogen_count_when_a_hydrogen_goes() {
    // Methane with one explicit hydrogen, made a methyl radical by declaring
    // the radical and then deleting that hydrogen.
    let source = topology("[H][CH3]");
    let mut editor = TopologyEditor::from_topology(source);
    let carbon = handle(&editor, 1);
    let radical = editor
        .atom(carbon)
        .unwrap()
        .clone()
        .with_radical(AtomRadical::new(1, None));
    editor.replace_atom(carbon, radical).unwrap();
    editor.delete_atom(handle(&editor, 0)).unwrap();
    let edited = editor.finish().unwrap();
    assert_eq!(fixed_counts(&edited), [ImplicitHydrogens::Fixed(3)]);
    let methyl = edited.atom_at(TopologyAtomIndex::new(0)).unwrap();
    assert!(methyl.atom().radical().is_some());
}

#[test]
fn a_replaced_atom_owes_what_its_new_state_still_lacks() {
    // Methylating pyridine's nitrogen gives it a fourth bond, which only
    // a +1 charge allows; a new isotope does not.
    let source = topology("[CH]1=[CH][CH]=[N][CH]=[CH]1.[CH4]");
    let methylated = |charge: i8, isotope: Option<u16>| {
        let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
        let (nitrogen, carbon) = (handle(&editor, 3), handle(&editor, 6));
        editor
            .add_bond(nitrogen, carbon, BondOrder::Single)
            .unwrap();
        let mut atom = editor
            .atom(nitrogen)
            .unwrap()
            .clone()
            .with_formal_charge(charge);
        if let Some(mass) = isotope {
            atom = atom.with_isotope(Some(std::num::NonZeroU16::new(mass).unwrap()));
        }
        editor.replace_atom(nitrogen, atom).unwrap();
        editor.finish()
    };
    let error = methylated(0, Some(15)).unwrap_err();
    assert!(error.to_string().contains("charge"), "{error}");
    let pyridinium = methylated(1, None).unwrap();
    assert_eq!(hydrogens(&pyridinium, 3), ImplicitHydrogens::Fixed(0));
    assert_eq!(hydrogens(&pyridinium, 6), ImplicitHydrogens::Fixed(3));
}

#[test]
fn attaching_to_a_radical_takes_its_electron() {
    let mut editor = StructureEditor::new();
    let element = Element::from_symbol("C").unwrap();
    let methyl = Atom::new(element)
        .with_hydrogens(ImplicitHydrogens::Fixed(3))
        .with_radical(AtomRadical::new(1, None));
    let parent = editor
        .place_atom(methyl, Quantity::new(Point3::origin(), NANOMETER))
        .unwrap();
    editor
        .attach_atom(
            parent,
            Atom::new(element).with_hydrogens(ImplicitHydrogens::Fixed(3)),
            BondOrder::Single,
            AttachmentSite::Toward(Vector3::new(1.0, 0.0, 0.0)),
        )
        .unwrap();
    let ethane = editor.finish().unwrap();
    assert_eq!(
        fixed_counts(ethane.topology()),
        [ImplicitHydrogens::Fixed(3); 2]
    );
    assert!(ethane
        .topology()
        .atoms()
        .all(|atom| atom.atom().radical().is_none()));
}

#[test]
fn deleting_a_deuterium_carrier_drops_the_stereo_it_told_apart() {
    // An implicit protium in the deuterium's place would give the centre two
    // protium atoms, so the stereo goes.
    let source = topology("[2H][C@]([H])(F)Cl");
    let mut editor = TopologyEditor::from_topology(source);
    editor.delete_atom(handle(&editor, 0)).unwrap();
    let edited = editor.finish().unwrap();
    let molecule = edited.instances().next().unwrap().molecule().clone();
    assert_eq!(molecule.stereo_elements().count(), 0);
}
