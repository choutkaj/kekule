//! Structures edited after reading keep their hierarchy identifiers, atom
//! order, and connection records consistent for the writers.
use kekule::diagnostics::UnrepresentablePolicy;
use kekule::hydrogens::AddHydrogensOptions;
use kekule::macromolecular::StructureInterpretation;
use kekule::math::Point3;
use kekule::mmcif;
use kekule::pdb::{self, PdbSource, PdbWriteOptions};
use kekule::properties::{PropertyColumn, PropertyKey};
use kekule::structure::StructureEditor;
use kekule::topology::edit::EditAtomId;
use kekule::topology::InstanceAtomId;
use kekule::units::{Quantity, ANGSTROM};

/// A `HETATM` record of atom `name` of residue `residue` `sequence` of chain
/// `chain` at `x` Å.
fn hetatm(serial: usize, name: &str, residue: &str, chain: char, sequence: i32, x: f64) -> String {
    let element = name.trim_end_matches(char::is_numeric);
    let name = if element.len() == 2 {
        format!("{name:<4}")
    } else {
        format!(" {name:<3}")
    };
    format!(
        "HETATM{serial:>5} {name} {residue:>3} {chain}{sequence:>4}    {x:>8.3}{:>8.3}{:>8.3}  1.00 20.00          {element:>2}",
        0.0, 0.0
    )
}

fn entry(records: &[String]) -> String {
    let mut lines = records.to_vec();
    lines.push("END".to_owned());
    lines.join("\n") + "\n"
}

fn read(text: &str) -> StructureInterpretation {
    pdb::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

fn records<'a>(text: &'a str, name: &str) -> Vec<&'a str> {
    text.lines().filter(|line| line.starts_with(name)).collect()
}

fn atom(editor: &StructureEditor, index: usize) -> EditAtomId {
    editor.topology().atom_ids().nth(index).unwrap()
}

fn omit_hydrogens() -> PdbWriteOptions {
    PdbWriteOptions::default().with_implicit_hydrogens(UnrepresentablePolicy::Omit)
}

#[test]
fn added_hydrogens_join_their_residue_in_atom_order() {
    let source = read(&entry(&[
        hetatm(1, "O", "HOH", 'A', 1, 0.0),
        hetatm(2, "O", "HOH", 'A', 2, 3.0),
    ]))
    .into_structure();
    let mut editor = source.to_editor();
    let oxygen = atom(&editor, 0);
    editor
        .add_hydrogens([oxygen], AddHydrogensOptions::default())
        .unwrap();
    let edited = editor.finish().unwrap();
    let (written, _) = pdb::write(&edited, omit_hydrogens()).unwrap();
    let residues = records(&written, "HETATM")
        .iter()
        .map(|line| line[22..26].trim().to_owned())
        .collect::<Vec<_>>();
    assert_eq!(residues, ["1", "1", "1", "2"]);
}

#[test]
fn a_renamed_chain_is_written_under_its_new_identifiers() {
    let source = read(&entry(&[hetatm(1, "O", "HOH", 'A', 1, 0.0)])).into_structure();
    let mut editor = source.to_editor();
    let (chain, _) = editor.topology().chains().next().unwrap();
    editor
        .set_chain_identifiers(chain, "B", Some("B".to_owned()))
        .unwrap();
    let edited = editor.finish().unwrap();
    let site = edited.topology().atom_sites().next().unwrap();
    assert_eq!(site.metadata().auth_asym_id(), Some("B"));
    let (written, _) = pdb::write(&edited, omit_hydrogens()).unwrap();
    assert_eq!(&records(&written, "HETATM")[0][21..22], "B");
}

#[test]
fn a_site_moved_to_another_chain_takes_its_identifiers() {
    let source = read(&entry(&[
        hetatm(1, "NA", "NA", 'A', 1, 0.0),
        hetatm(2, "CL", "CL", 'B', 1, 5.0),
    ]))
    .into_structure();
    let mut editor = source.to_editor();
    let sodium = editor
        .topology()
        .atom_site_for_atom(atom(&editor, 0))
        .unwrap()
        .unwrap();
    let chloride = editor
        .topology()
        .atom_site_for_atom(atom(&editor, 1))
        .unwrap()
        .unwrap();
    let target = editor.topology().atom_site(chloride).unwrap().residue();
    editor.set_atom_site_residue(sodium, target).unwrap();
    let edited = editor.finish().unwrap();
    for site in edited.topology().atom_sites() {
        assert_eq!(site.metadata().auth_asym_id(), Some("B"));
        assert_eq!(site.residue().chain().label_id(), "B");
    }
}

/// Zinc ions 1 and 3 coordinate the waters 2 and 5; water 4 is free.
fn zinc_entry() -> String {
    let link = |zinc: i32, water: i32| {
        let gap = " ".repeat(15);
        format!("LINK        ZN    ZN A{zinc:>4} {gap} O   HOH A{water:>4}     1555   1555  2.10")
    };
    entry(&[
        link(1, 2),
        link(3, 5),
        hetatm(1, "ZN", "ZN", 'A', 1, 0.0),
        hetatm(2, "O", "HOH", 'A', 2, 2.1),
        hetatm(3, "ZN", "ZN", 'A', 3, 10.0),
        hetatm(4, "O", "HOH", 'A', 4, 20.0),
        hetatm(5, "O", "HOH", 'A', 5, 12.1),
    ])
}

#[test]
fn coordination_follows_its_atoms_through_edits() {
    let source = read(&zinc_entry());
    assert_eq!(source.report().coordination().len(), 2);
    let mut editor = source.structure().to_editor();
    let first = [atom(&editor, 0), atom(&editor, 1)];
    editor.delete_atoms(first).unwrap();
    let edited = editor.finish().unwrap();
    let (written, _) = pdb::write(
        PdbSource::structure(&edited).with_reports(std::slice::from_ref(source.report())),
        omit_hydrogens(),
    )
    .unwrap();
    let links = records(&written, "LINK");
    assert_eq!(links.len(), 1, "{written}");
    assert_eq!(&links[0][22..26], "   3");
    assert_eq!(&links[0][52..56], "   5");
}

#[test]
fn appending_a_chain_under_a_taken_identifier_is_reported() {
    let water = read(&entry(&[hetatm(1, "O", "HOH", 'A', 1, 0.0)])).into_structure();
    let ion = read(&entry(&[
        hetatm(1, "NA", "NA", 'A', 1, 5.0),
        hetatm(2, "CL", "CL", 'C', 2, 9.0),
    ]))
    .into_structure();
    let mut editor = water.to_editor();
    let appended = editor.append_structure(&ion).unwrap();
    let colliding = appended.report().colliding_chains();
    assert_eq!(colliding.len(), 1);
    let chain = ion.topology().chain(colliding[0]).unwrap();
    assert_eq!(chain.label_id(), "A");
}

fn two_waters() -> StructureEditor {
    read(&entry(&[
        hetatm(1, "O", "HOH", 'A', 1, 0.0),
        hetatm(2, "O", "HOH", 'A', 2, 3.0),
    ]))
    .into_structure()
    .to_editor()
}

#[test]
fn a_column_in_editor_order_follows_its_atoms_through_publication() {
    let mut editor = two_waters();
    let oxygen = atom(&editor, 0);
    editor
        .add_hydrogens([oxygen], AddHydrogensOptions::default())
        .unwrap();
    let numbers = editor
        .topology()
        .atoms()
        .map(|(_, atom)| Some(i64::from(atom.element().atomic_number())))
        .collect::<Vec<_>>();
    let key = PropertyKey::new("element").unwrap();
    editor
        .properties_mut()
        .atoms_mut()
        .insert(key.clone(), PropertyColumn::Int(numbers))
        .unwrap();
    let edited = editor.finish().unwrap();
    let column = edited.properties().atoms().get(&key).unwrap().clone();
    let expected = edited
        .topology()
        .atoms()
        .map(|atom| Some(i64::from(atom.atom().element().atomic_number())))
        .collect::<Vec<_>>();
    assert_eq!(column, PropertyColumn::Int(expected));
}

#[test]
fn positions_keep_handle_order_when_an_added_atom_joins_a_residue() {
    let mut editor = two_waters();
    let positions = editor.positions();
    let unit = positions.values().unit();
    let before = positions.values().value().to_vec();
    let oxygen = atom(&editor, 0);
    editor
        .add_hydrogens([oxygen], AddHydrogensOptions::default())
        .unwrap();
    let mut points = before.clone();
    points.extend([Point3::new(0.09, 0.0, 0.0), Point3::new(-0.02, 0.09, 0.0)]);
    editor.set_positions(Quantity::new(points, unit)).unwrap();
    let second = atom(&editor, 1);
    let at = editor.position(second).unwrap().value_in(unit).unwrap();
    assert_eq!(at, before[1]);
}

#[test]
fn a_followed_report_keeps_coordination_through_a_chain_rename() {
    let source = read(&zinc_entry());
    let mut editor = source.structure().to_editor();
    let (chain, _) = editor.topology().chains().next().unwrap();
    editor
        .set_chain_identifiers(chain, "B", Some("B".to_owned()))
        .unwrap();
    let (edited, correspondence) = editor.finish_with_correspondence().unwrap();
    let followed =
        source
            .report()
            .follow_edit(source.topology(), edited.topology(), &correspondence);
    let (written, _) = pdb::write(
        PdbSource::structure(&edited).with_reports(std::slice::from_ref(&followed)),
        omit_hydrogens(),
    )
    .unwrap();
    let links = records(&written, "LINK");
    assert_eq!(links.len(), 2, "{written}");
    assert!(links.iter().all(|link| &link[21..22] == "B"), "{written}");
}

/// A `HETATM` record of residue `sequence` with a blank chain and segment
/// `segment`.
fn segment_hetatm(
    serial: usize,
    name: &str,
    residue: &str,
    segment: &str,
    sequence: i32,
    x: f64,
) -> String {
    let record = hetatm(serial, name, residue, ' ', sequence, x);
    format!("{}{segment:<4}{}", &record[..72], &record[76..])
}

#[test]
fn coordination_in_segments_that_share_residue_numbers_is_written() {
    // Each segment holds zinc 1, coordinated by the oxygen of water 2.
    let source = read(&entry(&[
        segment_hetatm(1, "ZN", "ZN", "PROA", 1, 0.0),
        segment_hetatm(2, "O", "HOH", "PROA", 2, 2.1),
        segment_hetatm(3, "ZN", "ZN", "PROB", 1, 10.0),
        segment_hetatm(4, "O", "HOH", "PROB", 2, 12.1),
        "CONECT    1    2".to_owned(),
        "CONECT    3    4".to_owned(),
    ]));
    assert_eq!(source.report().coordination().len(), 2);
    let (written, report) = mmcif::write([&source], Default::default()).unwrap();
    assert!(report.omissions().is_empty(), "{:?}", report.omissions());
    assert_eq!(written.matches(" metalc ").count(), 2, "{written}");
}

#[test]
fn a_report_follows_its_structure_through_an_append() {
    let water = read(&entry(&[hetatm(1, "O", "HOH", 'W', 1, 30.0)])).into_structure();
    let zinc = read(&zinc_entry());
    let mut editor = water.to_editor();
    let append = editor.append_structure(zinc.structure()).unwrap();
    let (edited, correspondence) = editor.finish_with_correspondence().unwrap();
    let followed = zinc
        .report()
        .follow_append(&append, edited.topology(), &correspondence);
    let (written, _) = pdb::write(
        PdbSource::structure(&edited).with_reports(std::slice::from_ref(&followed)),
        omit_hydrogens(),
    )
    .unwrap();
    assert_eq!(records(&written, "LINK").len(), 2, "{written}");
}

/// An atom by author chain, component, author sequence number, and name.
type AtomKey = (String, String, String, String);

/// Each atom's element, formal charge, hydrogen count, and bonded atoms.
type AtomChemistry = (String, i8, Option<usize>, Vec<AtomKey>);

fn atom_key(structure: &kekule::structure::Structure, atom: InstanceAtomId) -> AtomKey {
    let site = structure
        .topology()
        .atom(atom)
        .and_then(|atom| atom.atom_site())
        .unwrap();
    let residue = site.residue();
    let chain = residue.chain();
    (
        chain.author_id().unwrap_or(chain.label_id()).to_owned(),
        residue.name().to_owned(),
        residue.author_seq_id().unwrap_or_default().to_owned(),
        site.metadata()
            .label_atom_id()
            .or(site.metadata().auth_atom_id())
            .unwrap_or_default()
            .to_owned(),
    )
}

/// The chemistry of every atom by its key.
fn atom_chemistry(
    structure: &kekule::structure::Structure,
) -> std::collections::BTreeMap<AtomKey, AtomChemistry> {
    structure
        .topology()
        .atoms()
        .map(|atom| {
            let mut neighbors = atom
                .neighbors()
                .map(|neighbor| atom_key(structure, neighbor.id()))
                .collect::<Vec<_>>();
            neighbors.sort();
            (
                atom_key(structure, atom.id()),
                (
                    atom.atom().element().symbol().to_owned(),
                    atom.atom().formal_charge(),
                    atom.total_hydrogens(),
                    neighbors,
                ),
            )
        })
        .collect()
}

/// The keys of the atoms whose chemistry differs between `left` and `right`,
/// with both sides.
fn chemistry_changes(
    left: &kekule::structure::Structure,
    right: &kekule::structure::Structure,
) -> Vec<(AtomKey, Option<AtomChemistry>, Option<AtomChemistry>)> {
    let (left, right) = (atom_chemistry(left), atom_chemistry(right));
    left.keys()
        .chain(right.keys())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .filter(|key| left.get(*key) != right.get(*key))
        .map(|key| (key.clone(), left.get(key).cloned(), right.get(key).cloned()))
        .collect()
}

fn pdb_fixture(text: &str) -> StructureInterpretation {
    pdb::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

/// A source with hydrogens added to every atom, and its report followed
/// through the edit.
fn with_hydrogens(
    source: &StructureInterpretation,
) -> (
    kekule::structure::Structure,
    kekule::macromolecular::StructureInterpretationReport,
) {
    let mut editor = source.structure().to_editor();
    let atoms = editor.topology().atom_ids().collect::<Vec<_>>();
    let added = editor
        .add_hydrogens(atoms, AddHydrogensOptions::default())
        .unwrap();
    assert!(added.added_count() > 0);
    let (edited, correspondence) = editor.finish_with_correspondence().unwrap();
    let followed =
        source
            .report()
            .follow_edit(source.topology(), edited.topology(), &correspondence);
    (edited, followed)
}

fn reread(written: &str) -> StructureInterpretation {
    mmcif::parse_str(written, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

#[test]
fn a_followed_report_writes_mmcif_with_added_hydrogens() {
    // Rubredoxin's iron is coordinated by four cysteines; its waters have
    // no hydrogens.
    let (edited, followed) = with_hydrogens(&pdb_fixture(include_str!("../fixtures/pdb/1IRO.pdb")));
    let (written, report) =
        mmcif::write(
            [mmcif::MmcifBlockSource::structure(&edited)
                .with_reports(std::slice::from_ref(&followed))],
            Default::default(),
        )
        .unwrap();
    assert!(report.omissions().is_empty(), "{:?}", report.omissions());
    assert_eq!(written.matches(" metalc ").count(), 4, "{written}");
    let reread = reread(&written);
    assert_eq!(reread.report().coordination().len(), 4);
    let changes = chemistry_changes(&edited, reread.structure());
    assert!(
        changes.is_empty(),
        "{:#?}",
        &changes[..changes.len().min(4)]
    );
}

#[test]
fn added_hydrogens_keep_their_bonds_through_mmcif() {
    // Added hydrogens take their templates' names, never one a template
    // places on another atom, as `H2` on a peptide's `N`.
    for text in [
        include_str!("../fixtures/pdb/1IRO.pdb"),
        include_str!("../fixtures/pdb/1AKE.pdb"),
    ] {
        let (edited, _) = with_hydrogens(&pdb_fixture(text));
        let (written, _) = mmcif::write([&edited], Default::default()).unwrap();
        let reread = reread(&written);
        let changes = chemistry_changes(&edited, reread.structure());
        assert!(
            changes.is_empty(),
            "{:#?}",
            &changes[..changes.len().min(4)]
        );
    }
}

#[test]
fn a_followed_report_writes_mmcif_with_a_molecule_added_outside_the_hierarchy() {
    let source = read(&zinc_entry());
    let mut editor = source.structure().to_editor();
    let mut sodium = kekule::molecule::MoleculeEditor::new();
    sodium
        .add_atom(
            kekule::molecule::Atom::new(kekule::molecule::Element::from_symbol("Na").unwrap())
                .with_formal_charge(1),
        )
        .unwrap();
    editor
        .add_molecule(
            sodium.finish().unwrap(),
            &kekule::structure::Positions::zeros(1),
        )
        .unwrap();
    let (edited, correspondence) = editor.finish_with_correspondence().unwrap();
    let followed =
        source
            .report()
            .follow_edit(source.topology(), edited.topology(), &correspondence);
    let (written, _) =
        mmcif::write(
            [mmcif::MmcifBlockSource::structure(&edited)
                .with_reports(std::slice::from_ref(&followed))],
            Default::default(),
        )
        .unwrap();
    // The zincs keep their coordination; the sodium is an asym of its own.
    assert_eq!(written.matches(" metalc ").count(), 2, "{written}");
    let reread = mmcif::parse_str(&written, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    assert_eq!(reread.report().coordination().len(), 2);
    assert_eq!(reread.structure().atom_count(), edited.atom_count());
    let chains = reread
        .structure()
        .topology()
        .chains()
        .map(|chain| chain.label_id().to_owned())
        .collect::<Vec<_>>();
    assert!(chains.last().unwrap().starts_with('M'), "{written}");
}

/// Reads back the block written from 1IRO with hydrogens whose second `HA`
/// is renamed `H2`, placed `between` its `CA` and the `N` the template bonds
/// an `H2` to, or left in place; returns its issues on `H2` and whether the
/// read changed any atom's chemistry.
fn reread_ha_as_h2(between: bool) -> (Vec<String>, bool) {
    let (edited, _) = with_hydrogens(&pdb_fixture(include_str!("../fixtures/pdb/1IRO.pdb")));
    let mut editor = edited.to_editor();
    let (site, hydrogen, residue, metadata) = editor
        .topology()
        .atom_sites()
        .filter(|(_, site)| site.metadata().auth_atom_id() == Some("HA"))
        .map(|(id, site)| (id, site.atom(), site.residue(), site.metadata().clone()))
        .nth(1)
        .unwrap();
    editor
        .set_atom_site_metadata(
            site,
            metadata
                .with_label_atom_id(Some("H2".to_owned()))
                .with_auth_atom_id(Some("H2".to_owned())),
        )
        .unwrap();
    if between {
        let point = |editor: &StructureEditor, name: &str| {
            let atom = editor
                .topology()
                .atom_sites()
                .find(|(_, site)| {
                    site.residue() == residue && site.metadata().auth_atom_id() == Some(name)
                })
                .map(|(_, site)| site.atom())
                .unwrap();
            editor.position(atom).unwrap().value_in(ANGSTROM).unwrap()
        };
        let (nitrogen, carbon, at) = (
            point(&editor, "N"),
            point(&editor, "CA"),
            point(&editor, "H2"),
        );
        // 1.0 Å from both: 0.68 Å off the middle of the N-CA bond, towards
        // where the hydrogen was.
        let middle = [
            (nitrogen.x + carbon.x) / 2.0,
            (nitrogen.y + carbon.y) / 2.0,
            (nitrogen.z + carbon.z) / 2.0,
        ];
        let axis = [
            carbon.x - nitrogen.x,
            carbon.y - nitrogen.y,
            carbon.z - nitrogen.z,
        ];
        let toward = [at.x - middle[0], at.y - middle[1], at.z - middle[2]];
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let along = dot(toward, axis) / dot(axis, axis);
        let off = [0, 1, 2].map(|index| toward[index] - along * axis[index]);
        let length = dot(off, off).sqrt();
        let half = dot(axis, axis).sqrt() / 2.0;
        let reach = (1.0 - half * half).sqrt();
        let moved = [0, 1, 2].map(|index| middle[index] + off[index] / length * reach);
        editor
            .set_position(
                hydrogen,
                Quantity::new(Point3::new(moved[0], moved[1], moved[2]), ANGSTROM),
            )
            .unwrap();
    }
    let renamed = editor.finish().unwrap();
    let (written, _) = mmcif::write([&renamed], Default::default()).unwrap();
    let reread = reread(&written);
    let issues = reread
        .report()
        .issues()
        .iter()
        .filter_map(|issue| match issue {
            kekule::macromolecular::StructureInterpretIssue::HydrogenWithSeveralBonds {
                atom_name,
                bonds,
                ..
            } => Some(format!("{atom_name} has {bonds} bonds")),
            kekule::macromolecular::StructureInterpretIssue::HydrogenBeyondTemplateBond {
                atom_name,
                anchor_name,
                ..
            } => Some(format!("{atom_name} is far from {anchor_name}")),
            _ => None,
        })
        .collect();
    let changed = !chemistry_changes(&renamed, reread.structure()).is_empty();
    (issues, changed)
}

#[test]
fn a_hydrogen_named_after_another_atom_reads_back_reported() {
    // An `HA` renamed `H2` keeps its bond to `CA`; the template bonds an
    // `H2` to `N`, 2.1 Å away, which is left out.
    assert_eq!(
        reread_ha_as_h2(false),
        (vec!["H2 is far from N".to_owned()], false)
    );
    // Within bonding distance of `N` as well, it takes both bonds.
    assert_eq!(
        reread_ha_as_h2(true),
        (vec!["H2 has 2 bonds".to_owned()], true)
    );
}

#[test]
fn added_hydrogens_on_a_residue_after_a_gap_write_without_a_report() {
    // KRas starts at residue 14, whose N has one fixed hydrogen for the gap
    // before it; an added hydrogen takes its place.
    let source = mmcif::parse_str(
        include_str!("../fixtures/mmcif/6OIM.cif"),
        Default::default(),
    )
    .unwrap()
    .interpret(Default::default())
    .unwrap();
    let (edited, _) = with_hydrogens(&source);
    let (written, report) = mmcif::write([&edited], Default::default()).unwrap();
    assert!(report.omissions().is_empty(), "{:?}", report.omissions());
    let changes = chemistry_changes(&edited, reread(&written).structure());
    assert!(
        changes.is_empty(),
        "{:#?}",
        &changes[..changes.len().min(4)]
    );
}

#[test]
fn a_hydrogen_added_beside_one_under_an_older_name_takes_the_other_name() {
    // A leucine's CB carries `1HB`, the older name of `HB2`; the hydrogen
    // added beside it is `HB3`, so the two read back as two template atoms.
    let (edited, _) = with_hydrogens(&pdb_fixture(include_str!("../fixtures/pdb/1AKE.pdb")));
    let mut editor = edited.to_editor();
    let leucine = |editor: &StructureEditor, name: &str| {
        let residue = editor
            .topology()
            .atom_sites()
            .find(|(_, site)| editor.topology().residue(site.residue()).unwrap().name() == "LEU")
            .map(|(_, site)| site.residue())
            .unwrap();
        editor
            .topology()
            .atom_sites()
            .find(|(_, site)| {
                site.residue() == residue && site.metadata().auth_atom_id() == Some(name)
            })
            .map(|(id, site)| (id, site.atom(), site.metadata().clone()))
            .unwrap()
    };
    let (site, _, metadata) = leucine(&editor, "HB2");
    editor
        .set_atom_site_metadata(
            site,
            metadata
                .with_label_atom_id(Some("1HB".to_owned()))
                .with_auth_atom_id(Some("1HB".to_owned())),
        )
        .unwrap();
    let (_, hydrogen, _) = leucine(&editor, "HB3");
    editor.delete_atoms([hydrogen]).unwrap();
    let (_, carbon, _) = leucine(&editor, "CB");
    editor
        .add_hydrogens([carbon], AddHydrogensOptions::default())
        .unwrap();
    // The new hydrogen is `HB3`.
    leucine(&editor, "HB3");
    let renamed = editor.finish().unwrap();
    let (written, _) = mmcif::write([&renamed], Default::default()).unwrap();
    let changes = chemistry_changes(&renamed, reread(&written).structure());
    assert!(
        changes.is_empty(),
        "{:#?}",
        &changes[..changes.len().min(4)]
    );
}

#[test]
fn a_hydrogen_far_from_the_atom_its_name_places_it_on_keeps_its_place() {
    // An isoleucine's `HG12` renamed `H2`, which its template bonds to `N`
    // 3.9 Å away, stays on `CG1`, and `N` keeps its charge.
    let (edited, _) = with_hydrogens(&pdb_fixture(include_str!("../fixtures/pdb/1AKE.pdb")));
    let mut editor = edited.to_editor();
    let (site, metadata) = editor
        .topology()
        .atom_sites()
        .find(|(_, site)| {
            editor.topology().residue(site.residue()).unwrap().name() == "ILE"
                && site.metadata().auth_atom_id() == Some("HG12")
        })
        .map(|(id, site)| (id, site.metadata().clone()))
        .unwrap();
    editor
        .set_atom_site_metadata(
            site,
            metadata
                .with_label_atom_id(Some("H2".to_owned()))
                .with_auth_atom_id(Some("H2".to_owned())),
        )
        .unwrap();
    let renamed = editor.finish().unwrap();
    let (written, _) = mmcif::write([&renamed], Default::default()).unwrap();
    let reread = reread(&written);
    let changes = chemistry_changes(&renamed, reread.structure());
    assert!(
        changes.is_empty(),
        "{:#?}",
        &changes[..changes.len().min(4)]
    );
    let far = reread
        .report()
        .issues()
        .iter()
        .filter_map(|issue| match issue {
            kekule::macromolecular::StructureInterpretIssue::HydrogenBeyondTemplateBond {
                atom_name,
                anchor_name,
                ..
            } => Some((atom_name.as_str(), anchor_name.as_str())),
            _ => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(far, [("H2", "N")]);
}
