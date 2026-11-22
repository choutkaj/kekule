//! Structures edited after reading keep their hierarchy identifiers, atom
//! order, and connection records consistent for the writers.
use kekule::diagnostics::UnrepresentablePolicy;
use kekule::hydrogens::AddHydrogensOptions;
use kekule::macromolecular::StructureInterpretation;
use kekule::pdb::{self, PdbSource, PdbWriteOptions};
use kekule::structure::StructureEditor;
use kekule::topology::edit::EditAtomId;

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
