//! Edits keep represented chemistry: hydrogen counts follow lost and gained
//! bond valence, and stereo carriers follow deleted hydrogens.
use kekule::math::{Point3, Vector3};
use kekule::molecule::{
    Atom, AtomId, AtomRadical, BondOrder, Element, ImplicitHydrogens, Molecule,
};
use kekule::perception::valence::open_valences;
use kekule::structure::placement::AttachmentSite;
use kekule::structure::StructureEditor;
use kekule::topology::edit::{EditAtomId, EditBond, EditBondId, TopologyEditor};
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

fn fixed_counts(topology: &Topology) -> Vec<ImplicitHydrogens> {
    (0..topology.atom_count() as u32)
        .map(|index| hydrogens(topology, index))
        .collect()
}

/// The published draft, which must perceive.
fn published(editor: TopologyEditor) -> Topology {
    let topology = editor.finish().unwrap();
    topology.perceived().unwrap()
}

/// The bond between atoms `a` and `b` of `editor`.
fn bond_of(editor: &TopologyEditor, a: usize, b: usize) -> EditBondId {
    let (a, b) = (handle(editor, a), handle(editor, b));
    editor
        .bonds()
        .find(|(_, bond)| (bond.a() == a && bond.b() == b) || (bond.a() == b && bond.b() == a))
        .unwrap()
        .0
}

/// Cisplatin with every count fixed: Pt 0, Cl 0, N 3, the ammines dative.
fn cisplatin() -> Arc<Topology> {
    let mut editor = TopologyEditor::from_topology(topology("[Cl][Pt]([Cl])([NH3])[NH3]"));
    for nitrogen in [3, 4] {
        let bond = bond_of(&editor, 1, nitrogen);
        editor.set_bond_order(bond, BondOrder::Dative).unwrap();
    }
    Arc::new(published(editor))
}

#[test]
fn a_dative_bond_to_a_metal_keeps_hydrogen_counts() {
    // The ammines coordinate platinum through dative bonds, which trade no
    // hydrogen. The chlorides' bonds are single: a chloride cut off takes a
    // hydrogen.
    let source = cisplatin();
    let ammine = [ImplicitHydrogens::Fixed(3); 2];
    assert_eq!(fixed_counts(&source)[3..], ammine);
    let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
    editor.delete_bond(bond_of(&editor, 0, 1)).unwrap();
    let edited = published(editor);
    assert_eq!(hydrogens(&edited, 0), ImplicitHydrogens::Fixed(1));
    assert_eq!(hydrogens(&edited, 1), ImplicitHydrogens::Fixed(0));
    let ligands = [
        ImplicitHydrogens::Fixed(1),
        ImplicitHydrogens::Fixed(1),
        ImplicitHydrogens::Fixed(3),
        ImplicitHydrogens::Fixed(3),
    ];
    let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
    editor.delete_atom(handle(&editor, 1)).unwrap();
    assert_eq!(fixed_counts(&published(editor)), ligands);
    let without_metal = AtomSelection::from_predicate(Arc::clone(&source), |atom| {
        atom.atom().element().symbol() != "Pt"
    });
    let subset = Topology::subset(Arc::clone(&source), &without_metal).unwrap();
    assert_eq!(fixed_counts(subset.topology()), ligands);
    subset.topology().perceived().unwrap();
}

#[test]
fn a_single_bond_to_a_metal_trades_a_hydrogen() {
    // A Grignard reagent, sodium methoxide, and a platinum amide: each
    // single bond to the metal takes a hydrogen from the non-metal, and the
    // metal's count stays.
    for (text, [non_metal, metal], counts) in [
        (
            "[CH4].[Mg][Br]",
            [0, 1],
            [3, 0, 0].map(ImplicitHydrogens::Fixed).to_vec(),
        ),
        (
            "[CH3][OH].[Na]",
            [1, 2],
            [3, 0, 0].map(ImplicitHydrogens::Fixed).to_vec(),
        ),
        (
            "[NH3].[Pt]",
            [0, 1],
            [2, 0].map(ImplicitHydrogens::Fixed).to_vec(),
        ),
    ] {
        let mut editor = TopologyEditor::from_topology(topology(text));
        let (a, b) = (handle(&editor, non_metal), handle(&editor, metal));
        let bond = editor.add_bond(a, b, BondOrder::Single).unwrap();
        let bonded = editor.clone();
        assert_eq!(fixed_counts(&published(bonded)), counts, "{text}");
        // Deleting the bond gives the hydrogen back.
        editor.delete_bond(bond).unwrap();
        assert_eq!(
            fixed_counts(&published(editor)),
            fixed_counts(&topology(text)),
            "{text}"
        );
    }
    // A dative bond, as an ammine's, takes none.
    let mut editor = TopologyEditor::from_topology(topology("[NH3].[Pt]"));
    let (a, b) = (handle(&editor, 0), handle(&editor, 1));
    editor.add_bond(a, b, BondOrder::Dative).unwrap();
    assert_eq!(
        hydrogens(&published(editor), 0),
        ImplicitHydrogens::Fixed(3)
    );
}

#[test]
fn attaching_lithium_to_ethane_gives_ethyllithium() {
    let mut editor = StructureEditor::new();
    let carbon = Element::from_symbol("C").unwrap();
    let methyl = |count| Atom::new(carbon).with_hydrogens(ImplicitHydrogens::Fixed(count));
    let first = editor
        .place_atom(methyl(4), Quantity::new(Point3::origin(), NANOMETER))
        .unwrap();
    let toward = |x| AttachmentSite::Toward(Vector3::new(x, 0.0, 0.0));
    editor
        .attach_atom(first, methyl(3), BondOrder::Single, toward(1.0))
        .unwrap();
    let lithium =
        Atom::new(Element::from_symbol("Li").unwrap()).with_hydrogens(ImplicitHydrogens::Fixed(0));
    editor
        .attach_atom(first, lithium, BondOrder::Single, toward(-1.0))
        .unwrap();
    let ethyllithium = editor.finish().unwrap();
    let topology = ethyllithium.topology().perceived().unwrap();
    assert_eq!(
        fixed_counts(&topology),
        [2, 3, 0].map(ImplicitHydrogens::Fixed)
    );
    let mut editor = TopologyEditor::from_topology(Arc::new(topology));
    editor.delete_atom(handle(&editor, 2)).unwrap();
    assert_eq!(
        fixed_counts(&published(editor)),
        [ImplicitHydrogens::Fixed(3); 2]
    );
}

#[test]
fn a_carbon_cut_from_lithium_takes_a_hydrogen() {
    let source = topology("[CH3][Li]");
    let mut editor = TopologyEditor::from_topology(source);
    editor.delete_atom(handle(&editor, 1)).unwrap();
    assert_eq!(
        fixed_counts(&published(editor)),
        [ImplicitHydrogens::Fixed(4)]
    );
}

#[test]
fn a_lost_bond_fills_only_to_the_smallest_valence_left() {
    // Each case deletes atom `deleted` and expects the count of atom `atom`
    // of what remains.
    for (text, deleted, atom, count) in [
        // Dimethyl sulfoxide less its oxygen is dimethyl sulfide.
        ("[CH3][S](=[O])[CH3]", 2, 1, 0),
        // Trimethylphosphine oxide less its oxygen is trimethylphosphine.
        ("[O]=[P]([CH3])([CH3])[CH3]", 0, 0, 0),
        // Methyl phosphate less its double-bonded oxygen is a phosphite.
        ("[CH3][O][P](=[O])([OH])[OH]", 3, 2, 0),
        // An H-phosphonate keeps its P-H as an implicit hydrogen.
        ("[H][P](=[O])([O][CH3])[O][CH3]", 0, 0, 1),
        // A thiol less its carbon is hydrogen sulfide.
        ("[CH3][SH]", 0, 0, 2),
        // Methylammonium less its carbon is ammonium.
        ("[CH3][NH3+]", 0, 0, 4),
        // Ethane less a carbon is methane.
        ("[CH3][CH3]", 0, 0, 4),
    ] {
        let mut editor = TopologyEditor::from_topology(topology(text));
        editor.delete_atom(handle(&editor, deleted)).unwrap();
        let edited = published(editor);
        assert_eq!(
            hydrogens(&edited, atom),
            ImplicitHydrogens::Fixed(count),
            "{text}"
        );
    }
}

#[test]
fn deleting_atoms_and_taking_the_rest_as_a_subset_agree() {
    for (text, deleted) in [
        // A chloride-bridged platinum dimer less one platinum.
        ("[Cl][Pt]1([Cl])[Cl][Pt]([Cl])([Cl])[Cl]1", vec![1]),
        // Sulfur hexafluoride less three fluorines is SF3H.
        ("[S]([F])([F])([F])([F])([F])[F]", vec![1, 2, 3]),
        ("[CH3][S](=[O])[CH3]", vec![0, 2]),
        ("[CH3][C]1=[CH][CH]=[CH][CH]=[CH]1", vec![0]),
        ("[CH3][CH2][OH]", vec![2]),
        ("[CH3][O][P](=[O])([OH])[OH]", vec![0, 3]),
    ] {
        let source = topology(text);
        let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
        let handles = deleted.iter().map(|&index| handle(&editor, index));
        editor.delete_atoms(handles.collect::<Vec<_>>()).unwrap();
        let edited = published(editor);
        let kept = (0..source.atom_count() as u32)
            .filter(|index| !deleted.contains(&(*index as usize)))
            .map(TopologyAtomIndex::new);
        let selection = AtomSelection::from_indices(Arc::clone(&source), kept).unwrap();
        let subset = Topology::subset(Arc::clone(&source), &selection).unwrap();
        let subset = subset.topology().perceived().unwrap();
        let counts = |topology: &Topology| {
            let mut counts = topology
                .atoms()
                .map(|atom| {
                    (
                        atom.atom().element().atomic_number(),
                        atom.atom().hydrogens(),
                    )
                })
                .collect::<Vec<_>>();
            counts.sort_by_key(|&(number, hydrogens)| (number, format!("{hydrogens:?}")));
            counts
        };
        assert_eq!(counts(&edited), counts(&subset), "{text}");
    }
    let mut editor = TopologyEditor::from_topology(topology("[S]([F])([F])([F])([F])([F])[F]"));
    let fluorines = (1..4)
        .map(|index| handle(&editor, index))
        .collect::<Vec<_>>();
    editor.delete_atoms(fluorines).unwrap();
    assert_eq!(
        hydrogens(&published(editor), 0),
        ImplicitHydrogens::Fixed(1)
    );
}

#[test]
fn replacing_a_bond_settles_its_endpoints_once() {
    // Ethene's double bond becomes a single bond from one of its carbons to
    // methanol's oxygen: methane and dimethyl ether.
    let source = topology("[CH2]=[CH2].[CH3][OH]");
    let mut editor = TopologyEditor::from_topology(source);
    let double = bond_of(&editor, 0, 1);
    let (carbon, oxygen) = (handle(&editor, 0), handle(&editor, 3));
    editor
        .replace_bond(double, EditBond::new(carbon, oxygen, BondOrder::Single))
        .unwrap();
    let edited = published(editor);
    assert_eq!(
        fixed_counts(&edited),
        [3, 4, 3, 0].map(ImplicitHydrogens::Fixed)
    );
}

/// Seeded xorshift.
struct Random(u64);

impl Random {
    fn below(&mut self, bound: usize) -> usize {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 % bound as u64) as usize
    }
}

#[test]
fn random_edits_of_molecules_with_metals_publish_only_what_perceives() {
    let sources = [
        "[CH3][Li]",
        "[CH3][Mg][Br].[CH4]",
        "[CH3][OH].[Na].[Cl]",
        "[NH3].[Pt]([Cl])[Cl].[CH3][CH2][SH]",
        "[CH3][S](=[O])[CH3].[Zn]",
        "[O]=[P]([CH3])([CH3])[CH3].[K]",
        "[CH2]=[CH][CH]=[CH2].[Fe]",
    ];
    let orders = [
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Dative,
        BondOrder::Zero,
    ];
    let mut random = Random(0x9e37_79b9_7f4a_7c15);
    let mut published_count = 0;
    for text in sources {
        let source = topology(text);
        for _ in 0..60 {
            let mut editor = TopologyEditor::from_topology(Arc::clone(&source));
            for _ in 0..1 + random.below(4) {
                let atoms = editor.atom_ids().collect::<Vec<_>>();
                let bonds = editor.bonds().map(|(id, _)| id).collect::<Vec<_>>();
                if atoms.len() < 2 {
                    break;
                }
                let atom = |random: &mut Random| atoms[random.below(atoms.len())];
                let order = orders[random.below(orders.len())];
                // Edits the editor refuses, such as a duplicate bond, leave
                // the draft as it was.
                let _ = match random.below(5) {
                    0 => editor.delete_atom(atom(&mut random)).map(|_| ()),
                    1 if !bonds.is_empty() => editor
                        .delete_bond(bonds[random.below(bonds.len())])
                        .map(|_| ()),
                    2 => {
                        let (a, b) = (atom(&mut random), atom(&mut random));
                        editor.add_bond(a, b, order).map(|_| ())
                    }
                    3 if !bonds.is_empty() => {
                        editor.set_bond_order(bonds[random.below(bonds.len())], order)
                    }
                    4 if !bonds.is_empty() => {
                        let bond = bonds[random.below(bonds.len())];
                        let (a, b) = (atom(&mut random), atom(&mut random));
                        editor
                            .replace_bond(bond, EditBond::new(a, b, order))
                            .map(|_| ())
                    }
                    _ => Ok(()),
                };
            }
            if let Ok(topology) = editor.finish() {
                published_count += 1;
                if let Err(error) = topology.perceived() {
                    panic!("{text}: {error}");
                }
            }
        }
    }
    assert!(published_count > 100, "{published_count}");
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

#[test]
fn a_deleted_hydrogen_carrier_leaves_its_parent_count_as_settled() {
    // An inferred parent stays inferred and perceives the hydrogen that
    // takes the stereo place.
    let mut editor = TopologyEditor::from_topology(topology("[H][C@](F)(Cl)Br"));
    let carbon = handle(&editor, 1);
    let inferred = editor
        .atom(carbon)
        .unwrap()
        .clone()
        .with_hydrogens(ImplicitHydrogens::Inferred);
    editor.replace_atom(carbon, inferred).unwrap();
    let source = Arc::new(published(editor));
    let mut editor = TopologyEditor::from_topology(source);
    editor.delete_atom(handle(&editor, 0)).unwrap();
    let edited = published(editor);
    assert_eq!(hydrogens(&edited, 0), ImplicitHydrogens::Inferred);
    let center = edited.atom_at(TopologyAtomIndex::new(0)).unwrap();
    assert_eq!(center.total_hydrogens(), Some(1));
    let molecule = edited.instances().next().unwrap().molecule();
    assert_eq!(molecule.stereo_elements().count(), 1);

    // A fixed radical centre gains no hydrogen, so the stereo the deleted
    // hydrogen told apart goes.
    let mut editor = TopologyEditor::from_topology(topology("[H][C@](F)(Cl)Br"));
    let carbon = handle(&editor, 1);
    let radical = editor
        .atom(carbon)
        .unwrap()
        .clone()
        .with_radical(AtomRadical::new(1, None));
    editor.replace_atom(carbon, radical).unwrap();
    editor.delete_atom(handle(&editor, 0)).unwrap();
    let edited = published(editor);
    assert_eq!(hydrogens(&edited, 0), ImplicitHydrogens::Fixed(0));
    let molecule = edited.instances().next().unwrap().molecule();
    assert_eq!(molecule.stereo_elements().count(), 0);
}
