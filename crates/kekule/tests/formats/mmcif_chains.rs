//! mmCIF round trips keep each chain's author name and residues, and group
//! the instances of one chemical component into one entity.

use kekule::macromolecular::StructureInterpretation;
use kekule::math::Point3;
use kekule::mmcif::{self, MmcifWriteOptions};
use kekule::molecule::{Atom, AtomId, BondOrder, Element, MoleculeEditor};
use kekule::pdb;
use kekule::structure::{Positions, Structure, StructureBuilder, StructureEditor};
use kekule::topology::{AtomSiteMetadata, InstanceAtomId};
use kekule::units::{Quantity, ANGSTROM};

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

#[test]
fn a_site_without_an_author_id_takes_the_one_its_chain_states() {
    // Two waters of chain A, one stating author chain A, one stating none.
    let mut builder = StructureBuilder::new();
    let chain = builder
        .topology_mut()
        .hierarchy_mut()
        .add_chain("A", None)
        .unwrap();
    for (sequence, author) in [(1, Some("A".to_owned())), (2, None)] {
        let mut water = MoleculeEditor::new();
        water
            .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
            .unwrap();
        let position = Point3::new(5.0 * f64::from(sequence), 0.0, 0.0);
        let positions = Positions::new(Quantity::new(vec![position], ANGSTROM)).unwrap();
        let instance = builder
            .add_molecule(water.finish().unwrap(), &positions)
            .unwrap();
        let hierarchy = builder.topology_mut().hierarchy_mut();
        let residue = hierarchy
            .add_residue(chain, "HOH", None, Some(sequence.to_string()), None)
            .unwrap();
        hierarchy
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, AtomId::new(0)),
                AtomSiteMetadata::default()
                    .with_label_atom_id(Some("O".to_owned()))
                    .with_auth_asym_id(author),
            )
            .unwrap();
    }
    let structure = builder.build().unwrap();
    let (text, _) = mmcif::write([&structure], MmcifWriteOptions::default()).unwrap();
    let water = [("HOH", "1"), ("HOH", "2")]
        .map(|(name, sequence)| (name.to_owned(), sequence.to_owned()))
        .to_vec();
    assert_eq!(
        author_chains(&read_back(&text)),
        [("A".to_owned(), water)],
        "{text}"
    );
}

fn mmcif_fixture(name: &str) -> StructureInterpretation {
    let path = format!("{}/tests/fixtures/mmcif/{name}", env!("CARGO_MANIFEST_DIR"));
    mmcif::parse_str(&std::fs::read_to_string(path).unwrap(), Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

/// Adds residue `component` of one atom `name`, of `element` with `charge`,
/// to the chain labelled `chain`, a new one of author chain `A` where none
/// is, at `position` Å, numbered `sequence`.
fn add_ion(
    editor: &mut StructureEditor,
    chain: &str,
    (component, name, element, charge): (&str, &str, &str, i8),
    sequence: &str,
    position: [f64; 3],
) {
    let mut ion = MoleculeEditor::new();
    ion.add_atom(Atom::new(Element::from_symbol(element).unwrap()).with_formal_charge(charge))
        .unwrap();
    let [x, y, z] = position;
    let positions = Positions::new(Quantity::new(vec![Point3::new(x, y, z)], ANGSTROM)).unwrap();
    let molecule = editor
        .add_molecule(ion.finish().unwrap(), &positions)
        .unwrap();
    let atom = molecule.atoms()[&AtomId::new(0)];
    let found = editor
        .topology()
        .chains()
        .find(|(_, candidate)| candidate.label_id() == chain)
        .map(|(id, _)| id);
    let chain = found.unwrap_or_else(|| editor.add_chain(chain, Some("A".to_owned())).unwrap());
    let residue = editor
        .add_residue(chain, component, None, Some(sequence.to_owned()), None)
        .unwrap();
    editor
        .add_atom_site(
            residue,
            atom,
            AtomSiteMetadata::default().with_label_atom_id(Some(name.to_owned())),
        )
        .unwrap();
}

/// The block written from `structure`, with `report` followed from `source`
/// when given.
fn write_edit(
    source: &StructureInterpretation,
    editor: StructureEditor,
    with_report: bool,
) -> String {
    let (edited, correspondence) = editor.finish_with_correspondence().unwrap();
    let followed =
        source
            .report()
            .follow_edit(source.topology(), edited.topology(), &correspondence);
    let block = if with_report {
        mmcif::MmcifBlockSource::structure(&edited).with_reports(std::slice::from_ref(&followed))
    } else {
        mmcif::MmcifBlockSource::structure(&edited)
    };
    let (text, _) = mmcif::write([block], MmcifWriteOptions::default()).unwrap();
    text
}

/// The asym and entity of each atom-site row of component `component`.
fn component_rows(text: &str, component: &str) -> Vec<(String, String)> {
    let document = mmcif::parse_str(text, Default::default()).unwrap();
    let table = document.blocks()[0].loop_with_tag("_atom_site.id").unwrap();
    (0..table.row_count())
        .filter(|&row| table.value(row, "_atom_site.label_comp_id").unwrap().text() == component)
        .map(|row| {
            ["_atom_site.label_asym_id", "_atom_site.label_entity_id"]
                .map(|tag| table.value(row, tag).unwrap().text().to_owned())
                .into()
        })
        .collect()
}

#[test]
fn one_component_name_of_other_chemistry_is_another_entity() {
    // Ligand LIG is a carbon in chains A and B and a nitrogen in chain C.
    let ligand = |chain, element| Site {
        chain,
        name: if element == "C" { "C1" } else { "N1" },
        element,
        residue: "LIG",
        ..ion("NA", "", 1, [0.0, 0.0, 0.0])
    };
    let source = read(&[
        Site {
            position: [0.0, 0.0, 0.0],
            ..ligand('A', "C")
        },
        Site {
            position: [10.0, 0.0, 0.0],
            ..ligand('B', "C")
        },
        Site {
            position: [20.0, 0.0, 0.0],
            ..ligand('C', "N")
        },
    ]);
    for written in [
        mmcif::write([&source], MmcifWriteOptions::default()),
        mmcif::write([source.structure()], MmcifWriteOptions::default()),
    ] {
        let (text, _) = written.unwrap();
        let entities = asym_entities(&text);
        assert_eq!(entities[0].1, entities[1].1, "{text}");
        assert_ne!(entities[0].1, entities[2].1, "{text}");
    }
}

#[test]
fn an_ion_added_to_a_ligand_chain_is_an_instance_of_its_own() {
    // Rubredoxin's iron is asym B of entity 2.
    let source = mmcif_fixture("1IRO.cif");
    for with_report in [true, false] {
        let mut editor = source.structure().to_editor();
        add_ion(
            &mut editor,
            "B",
            ("NA", "NA", "Na", 1),
            "200",
            [30.0, 0.0, 0.0],
        );
        let text = write_edit(&source, editor, with_report);
        let [(iron_asym, iron_entity)] = component_rows(&text, "FE").try_into().unwrap();
        let [(sodium_asym, sodium_entity)] = component_rows(&text, "NA").try_into().unwrap();
        assert_eq!(iron_asym, "B");
        assert_ne!(sodium_asym, iron_asym, "{text}");
        assert_ne!(sodium_entity, iron_entity, "{text}");
    }
}

#[test]
fn a_ligand_added_to_a_polymer_chain_is_an_instance_of_its_own() {
    let source = mmcif_fixture("1IRO.cif");
    for with_report in [true, false] {
        let mut editor = source.structure().to_editor();
        add_ion(
            &mut editor,
            "A",
            ("NA", "NA", "Na", 1),
            "200",
            [30.0, 0.0, 0.0],
        );
        let text = write_edit(&source, editor, with_report);
        let [(sodium_asym, _)] = component_rows(&text, "NA").try_into().unwrap();
        assert_ne!(sodium_asym, "A", "{text}");
        let reread = read_back(&text);
        let chain = reread
            .topology()
            .chains()
            .find(|chain| chain.label_id() == sodium_asym)
            .unwrap();
        assert_eq!(chain.author_id(), Some("A"));
    }
}

#[test]
fn added_instances_of_source_entities_join_them() {
    // KRas: magnesium is entity 2 and water entity 5. A new chain Z holds
    // another magnesium and water, and a water lies outside the hierarchy.
    let source = mmcif_fixture("6OIM.cif");
    let mut editor = source.structure().to_editor();
    add_ion(
        &mut editor,
        "Z",
        ("MG", "MG", "Mg", 2),
        "300",
        [60.0, 0.0, 0.0],
    );
    add_ion(
        &mut editor,
        "Z",
        ("HOH", "O", "O", 0),
        "900",
        [70.0, 0.0, 0.0],
    );
    let mut water = MoleculeEditor::new();
    let oxygen = water
        .add_atom(Atom::new(Element::from_symbol("O").unwrap()))
        .unwrap();
    for _ in 0..2 {
        let hydrogen = water
            .add_atom(Atom::new(Element::from_symbol("H").unwrap()))
            .unwrap();
        water.add_bond(oxygen, hydrogen, BondOrder::Single).unwrap();
    }
    let points = [[80.0, 0.0, 0.0], [80.96, 0.0, 0.0], [79.76, 0.93, 0.0]]
        .map(|[x, y, z]| Point3::new(x, y, z));
    let positions = Positions::new(Quantity::new(points.to_vec(), ANGSTROM)).unwrap();
    editor
        .add_molecule(water.finish().unwrap(), &positions)
        .unwrap();
    let text = write_edit(&source, editor, true);
    let entities = |component| {
        component_rows(&text, component)
            .into_iter()
            .map(|(_, entity)| entity)
            .collect::<std::collections::BTreeSet<_>>()
    };
    assert_eq!(entities("MG"), ["2".to_owned()].into(), "{text}");
    assert_eq!(entities("HOH"), ["5".to_owned()].into(), "{text}");
    assert_eq!(component_rows(&text, "HOH").len(), 207 + 1 + 3);
}

/// A curated benchmark entry.
fn curated(name: &str) -> StructureInterpretation {
    use std::io::Read;
    let path = format!(
        "{}/../../benchmarks/datasets/bio/curated/pdb/{name}",
        env!("CARGO_MANIFEST_DIR")
    );
    let mut text = String::new();
    flate2::read::GzDecoder::new(std::fs::File::open(path).unwrap())
        .read_to_string(&mut text)
        .unwrap();
    mmcif::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

/// The `_entity.type` of each entity of a written block.
fn entity_types(text: &str) -> Vec<String> {
    let document = mmcif::parse_str(text, Default::default()).unwrap();
    let table = document.blocks()[0].loop_with_tag("_entity.id").unwrap();
    (0..table.row_count())
        .map(|row| table.value(row, "_entity.type").unwrap().text().to_owned())
        .collect()
}

#[test]
fn entries_write_without_their_reports() {
    // 7LPM ends its chain with a leucine a break leaves alone; 7X9L holds
    // a sucrose, a branched entity of two saccharide residues.
    for (name, types) in [
        (
            "pdb_00007lpm_xyz_v1-3.cif.gz",
            &["polymer", "non-polymer", "non-polymer", "water"][..],
        ),
        (
            "pdb_00007x9l_xyz_v1-1.cif.gz",
            &[
                "polymer",
                "branched",
                "non-polymer",
                "non-polymer",
                "non-polymer",
                "non-polymer",
                "water",
            ],
        ),
    ] {
        let source = curated(name);
        let (text, report) =
            mmcif::write([source.structure()], MmcifWriteOptions::default()).unwrap();
        assert!(report.omissions().is_empty(), "{:?}", report.omissions());
        assert_eq!(entity_types(&text), types, "{name}");
        let reread = read_back(&text);
        assert_eq!(reread.atom_count(), source.structure().atom_count());
        assert_eq!(
            reread.topology().bonds().count(),
            source.structure().topology().bonds().count()
        );
    }
}

#[test]
fn waters_whose_numbers_wrap_are_written_as_two_asyms() {
    // Chain W numbers its third water 1 again.
    let water = |sequence, x| Site {
        chain: 'W',
        ..ion("HOH", "", sequence, [x, 0.0, 0.0])
    };
    let source = read(&[water(1, 0.0), water(2, 5.0), water(1, 10.0)]);
    for written in [
        mmcif::write([&source], MmcifWriteOptions::default()),
        mmcif::write([source.structure()], MmcifWriteOptions::default()),
    ] {
        let (text, _) = written.unwrap();
        let asyms = asym_entities(&text);
        assert_eq!(asyms.len(), 2, "{text}");
        assert_eq!(asyms[0].1, asyms[1].1);
        let residues = |sequences: &[&str]| {
            sequences
                .iter()
                .map(|sequence| ("HOH".to_owned(), (*sequence).to_owned()))
                .collect::<Vec<_>>()
        };
        assert_eq!(
            author_chains(&read_back(&text)),
            [("W".to_owned(), residues(&["1", "2", "1"]))]
        );
    }
}
