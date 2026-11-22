//! mmCIF round trips keep each chain's author name and residues, and group
//! the instances of one chemical component into one entity.

use kekule::macromolecular::StructureInterpretation;
use kekule::mmcif::{self, MmcifWriteOptions};
use kekule::molecule::{Atom, Element, MoleculeEditor};
use kekule::pdb;
use kekule::structure::{Positions, Structure};

/// A coordinate record of atom `name` of residue `residue` `sequence`, with
/// chain column `chain` and segment `segment`, at `position` Å.
struct Site<'a> {
    record: &'a str,
    name: &'a str,
    element: &'a str,
    residue: &'a str,
    chain: char,
    sequence: i32,
    segment: &'a str,
    position: [f64; 3],
}

impl Site<'_> {
    fn line(&self, serial: usize) -> String {
        let Self {
            record,
            name,
            element,
            residue,
            chain,
            sequence,
            segment,
            position: [x, y, z],
        } = self;
        let name = if element.len() == 2 {
            format!("{name:<4}")
        } else {
            format!(" {name:<3}")
        };
        format!(
            "{record:<6}{serial:>5} {name} {residue:>3} {chain}{sequence:>4}    {x:>8.3}{y:>8.3}{z:>8.3}  1.00 20.00      {segment:<4}{element:>2}"
        )
    }
}

/// A one-atom residue with a blank chain column in segment `segment`.
fn ion<'a>(residue: &'a str, segment: &'a str, sequence: i32, position: [f64; 3]) -> Site<'a> {
    let element = if residue == "HOH" { "O" } else { residue };
    Site {
        record: "HETATM",
        name: element,
        element,
        residue,
        chain: ' ',
        sequence,
        segment,
        position,
    }
}

fn read(sites: &[Site<'_>]) -> StructureInterpretation {
    let mut lines = (1..)
        .zip(sites)
        .map(|(serial, site)| site.line(serial))
        .collect::<Vec<_>>();
    lines.push("END".to_owned());
    pdb::parse_str(&(lines.join("\n") + "\n"), Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

fn read_back(text: &str) -> Structure {
    mmcif::parse_str(text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
        .into_structure()
}

/// Each chain name with its residues, by component and author sequence
/// number, in order. A chain's name is its author ID, or its label where it
/// has none; chains of one name are one author chain.
fn author_chains(structure: &Structure) -> Vec<(String, Vec<(String, String)>)> {
    let mut chains = Vec::<(String, Vec<(String, String)>)>::new();
    for chain in structure.topology().chains() {
        let name = chain.author_id().unwrap_or(chain.label_id()).to_owned();
        let residues = chain.residues().map(|residue| {
            (
                residue.name().to_owned(),
                residue.author_seq_id().unwrap_or_default().to_owned(),
            )
        });
        match chains.iter_mut().find(|(known, _)| *known == name) {
            Some((_, known)) => known.extend(residues),
            None => chains.push((name, residues.collect())),
        }
    }
    chains
}

/// The `_struct_asym` rows of a written block: each asym with its entity.
fn asym_entities(text: &str) -> Vec<(String, String)> {
    let document = mmcif::parse_str(text, Default::default()).unwrap();
    let table = document.blocks()[0]
        .loop_with_tag("_struct_asym.id")
        .unwrap();
    (0..table.row_count())
        .map(|row| {
            ["_struct_asym.id", "_struct_asym.entity_id"]
                .map(|tag| table.value(row, tag).unwrap().text().to_owned())
                .into()
        })
        .collect()
}

/// A solvated peptide as a simulation's PDB writes it: blank chain columns,
/// with segments naming the chains. Segment `PROA` holds a glycine dipeptide
/// and a sodium ion, `ION1` two sodium ions, `ION2` three chloride ions, and
/// `WAT1` and `WAT2` waters.
fn solvated_peptide() -> StructureInterpretation {
    let glycine = |name, element, sequence, position| Site {
        record: "ATOM",
        name,
        element,
        residue: "GLY",
        chain: ' ',
        sequence,
        segment: "PROA",
        position,
    };
    read(&[
        glycine("N", "N", 1, [0.000, 0.000, 0.0]),
        glycine("CA", "C", 1, [1.458, 0.000, 0.0]),
        glycine("C", "C", 1, [2.009, 1.420, 0.0]),
        glycine("O", "O", 1, [1.251, 2.390, 0.0]),
        glycine("N", "N", 2, [3.332, 1.536, 0.0]),
        glycine("CA", "C", 2, [3.988, 2.836, 0.0]),
        glycine("C", "C", 2, [5.488, 2.736, 0.0]),
        glycine("O", "O", 2, [6.100, 1.680, 0.0]),
        ion("NA", "PROA", 3, [15.0, 0.0, 0.0]),
        ion("NA", "ION1", 1, [0.0, 10.0, 0.0]),
        ion("NA", "ION1", 2, [5.0, 10.0, 0.0]),
        ion("CL", "ION2", 1, [0.0, 20.0, 0.0]),
        ion("CL", "ION2", 2, [5.0, 20.0, 0.0]),
        ion("CL", "ION2", 3, [10.0, 20.0, 0.0]),
        ion("HOH", "WAT1", 1, [0.0, 30.0, 0.0]),
        ion("HOH", "WAT1", 2, [5.0, 30.0, 0.0]),
        ion("HOH", "WAT2", 1, [0.0, 40.0, 0.0]),
    ])
}

#[test]
fn chains_without_author_ids_read_back_under_their_own_names() {
    let source = solvated_peptide();
    let chains = author_chains(source.structure());
    let residues = |name: &str, sequences: &[&str]| {
        sequences
            .iter()
            .map(|sequence| (name.to_owned(), (*sequence).to_owned()))
            .collect::<Vec<_>>()
    };
    let mut peptide = residues("GLY", &["1", "2"]);
    peptide.extend(residues("NA", &["3"]));
    assert_eq!(
        chains,
        [
            ("PROA".to_owned(), peptide),
            ("ION1".to_owned(), residues("NA", &["1", "2"])),
            ("ION2".to_owned(), residues("CL", &["1", "2", "3"])),
            ("WAT1".to_owned(), residues("HOH", &["1", "2"])),
            ("WAT2".to_owned(), residues("HOH", &["1"])),
        ]
    );
    // With the source report and without it, each ion is an asym of its
    // own, written under the name of the chain that holds it.
    for written in [
        mmcif::write([&source], MmcifWriteOptions::default()),
        mmcif::write([source.structure()], MmcifWriteOptions::default()),
    ] {
        let (text, report) = written.unwrap();
        assert!(report.omissions().is_empty(), "{:?}", report.omissions());
        assert_eq!(
            asym_entities(&text)
                .into_iter()
                .map(|(asym, _)| asym)
                .collect::<Vec<_>>(),
            ["PROA", "PROA_2", "ION1", "ION1_2", "ION2", "ION2_2", "ION2_3", "WAT1", "WAT2"]
        );
        assert_eq!(author_chains(&read_back(&text)), chains, "{text}");
    }
}

#[test]
fn instances_of_one_component_are_one_entity() {
    let source = solvated_peptide();
    // The peptide is entity 1; every sodium is entity 2, every chloride 3,
    // and every water 4.
    for (written, entities) in [
        (
            mmcif::write([&source], MmcifWriteOptions::default()),
            ["K1", "K2", "K3", "K4"],
        ),
        (
            mmcif::write([source.structure()], MmcifWriteOptions::default()),
            ["1", "2", "3", "4"],
        ),
    ] {
        let (text, _) = written.unwrap();
        let [peptide, sodium, chloride, water] = entities;
        assert_eq!(
            asym_entities(&text)
                .into_iter()
                .map(|(_, entity)| entity)
                .collect::<Vec<_>>(),
            [peptide, sodium, sodium, sodium, chloride, chloride, chloride, water, water]
        );
        let document = mmcif::parse_str(&text, Default::default()).unwrap();
        let table = document.blocks()[0].loop_with_tag("_entity.id").unwrap();
        assert_eq!(
            (0..table.row_count())
                .map(|row| ["_entity.id", "_entity.type"].map(|tag| table
                    .value(row, tag)
                    .unwrap()
                    .text()
                    .to_owned()))
                .collect::<Vec<_>>(),
            [
                [peptide, "polymer"],
                [sodium, "non-polymer"],
                [chloride, "non-polymer"],
                [water, "water"],
            ]
            .map(|row| row.map(str::to_owned))
        );
        // The reader keeps one molecule per ion and water.
        let reread = read_back(&text);
        assert_eq!(
            reread.topology().instances().count(),
            source.structure().topology().instances().count()
        );
    }
}

#[test]
fn a_chain_without_an_author_id_takes_a_name_no_other_chain_states() {
    // The ions of segment ION1 have no author chain, and chain B states
    // ION1 as its author ID.
    let source = read(&[
        ion("NA", "ION1", 1, [0.0, 0.0, 0.0]),
        ion("NA", "ION1", 2, [5.0, 0.0, 0.0]),
        Site {
            chain: 'B',
            ..ion("CL", "", 1, [10.0, 0.0, 0.0])
        },
    ])
    .into_structure();
    let mut editor = source.to_editor();
    let (chain, _) = editor.topology().chains().nth(1).unwrap();
    editor
        .set_chain_identifiers(chain, "B", Some("ION1".to_owned()))
        .unwrap();
    let edited = editor.finish().unwrap();
    let (text, _) = mmcif::write([&edited], MmcifWriteOptions::default()).unwrap();
    let sodium =
        [("NA", "1"), ("NA", "2")].map(|(name, sequence)| (name.to_owned(), sequence.to_owned()));
    assert_eq!(
        author_chains(&read_back(&text)),
        [
            ("ION1-2".to_owned(), sodium.to_vec()),
            ("ION1".to_owned(), vec![("CL".to_owned(), "1".to_owned())]),
        ],
        "{text}"
    );
}

#[test]
fn a_molecule_outside_the_hierarchy_takes_a_name_no_chain_states() {
    // Chain B states M2 as its author ID, the name the writer would give
    // the sodium ion added outside the hierarchy, instance 1.
    let source = read(&[Site {
        chain: 'B',
        ..ion("CL", "", 1, [10.0, 0.0, 0.0])
    }])
    .into_structure();
    let mut editor = source.to_editor();
    let (chain, _) = editor.topology().chains().next().unwrap();
    editor
        .set_chain_identifiers(chain, "B", Some("M2".to_owned()))
        .unwrap();
    let mut sodium = MoleculeEditor::new();
    sodium
        .add_atom(Atom::new(Element::from_symbol("Na").unwrap()).with_formal_charge(1))
        .unwrap();
    editor
        .add_molecule(sodium.finish().unwrap(), &Positions::zeros(1))
        .unwrap();
    let edited = editor.finish().unwrap();
    let (text, _) = mmcif::write([&edited], MmcifWriteOptions::default()).unwrap();
    let names = read_back(&text)
        .topology()
        .chains()
        .map(|chain| {
            (
                chain.label_id().to_owned(),
                chain.author_id().map(str::to_owned),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(
        names,
        [
            ("B".to_owned(), Some("M2".to_owned())),
            ("M2".to_owned(), Some("M2-2".to_owned())),
        ],
        "{text}"
    );
}
