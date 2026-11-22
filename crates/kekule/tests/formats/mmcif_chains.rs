//! mmCIF round trips keep each chain's author name and residues, and group
//! the instances of one chemical component into one entity.

use kekule::diagnostics::FormatErrorDetail;
use kekule::macromolecular::{
    EntityKind, ModelSelection, StructureInterpretOptions, StructureInterpretation,
};
use kekule::math::Point3;
use kekule::mmcif::{
    self, MmcifBlockSource, MmcifEntityClassifications, MmcifWriteErrorKind, MmcifWriteOptions,
};
use kekule::molecule::{Atom, AtomId, BondOrder, Element, MoleculeEditor};
use kekule::pdb;
use kekule::structure::{Positions, Structure, StructureBuilder, StructureEditor};
use kekule::topology::{AtomSiteMetadata, InstanceAtomId, MoleculeClass};
use kekule::units::{Quantity, ANGSTROM};
use std::collections::{BTreeMap, BTreeSet};

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
            .collect::<BTreeSet<_>>()
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
        .interpret(StructureInterpretOptions::default().with_model_selection(ModelSelection::First))
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

/// Adds a glycine dipeptide, residues 1001 and 1002, to the chain labelled
/// `chain`, a new one of author chain `A` where none is.
fn add_dipeptide(editor: &mut StructureEditor, chain: &str) {
    let element = |symbol| Atom::new(Element::from_symbol(symbol).unwrap());
    let mut peptide = MoleculeEditor::new();
    // The N terminus is charged, as readers give it.
    let atoms = ["N", "C", "C", "O", "N", "C", "C", "O"]
        .map(element)
        .map(|atom| {
            let charge = i8::from(peptide.atom_count() == 0);
            peptide.add_atom(atom.with_formal_charge(charge)).unwrap()
        });
    for (left, right, order) in [
        (0, 1, BondOrder::Single),
        (1, 2, BondOrder::Single),
        (2, 3, BondOrder::Double),
        (2, 4, BondOrder::Single),
        (4, 5, BondOrder::Single),
        (5, 6, BondOrder::Single),
        (6, 7, BondOrder::Double),
    ] {
        peptide.add_bond(atoms[left], atoms[right], order).unwrap();
    }
    let points = [
        [0.000, 0.000],
        [1.458, 0.000],
        [2.009, 1.420],
        [1.251, 2.390],
        [3.332, 1.536],
        [3.988, 2.836],
        [5.488, 2.736],
        [6.100, 1.680],
    ]
    .map(|[x, y]| Point3::new(x + 60.0, y, 0.0));
    let positions = Positions::new(Quantity::new(points.to_vec(), ANGSTROM)).unwrap();
    let molecule = editor
        .add_molecule(peptide.finish().unwrap(), &positions)
        .unwrap();
    let found = editor
        .topology()
        .chains()
        .find(|(_, candidate)| candidate.label_id() == chain)
        .map(|(id, _)| id);
    let chain = found.unwrap_or_else(|| editor.add_chain(chain, Some("A".to_owned())).unwrap());
    let names = ["N", "CA", "C", "O"];
    for (residue, sequence) in [(0u32, "1001"), (1, "1002")] {
        let residue_id = editor
            .add_residue(chain, "GLY", None, Some(sequence.to_owned()), None)
            .unwrap();
        for (offset, name) in (0u32..).zip(names) {
            let atom = molecule.atoms()[&AtomId::new(4 * residue + offset)];
            editor
                .add_atom_site(
                    residue_id,
                    atom,
                    AtomSiteMetadata::default().with_label_atom_id(Some(name.to_owned())),
                )
                .unwrap();
        }
    }
}

/// Each polymer entity's `_entity_poly_seq` monomers, in entity order.
fn polymer_sequences(text: &str) -> Vec<Vec<String>> {
    let document = mmcif::parse_str(text, Default::default()).unwrap();
    let table = document.blocks()[0]
        .loop_with_tag("_entity_poly_seq.entity_id")
        .unwrap();
    let mut sequences = Vec::<(String, Vec<String>)>::new();
    for row in 0..table.row_count() {
        let [entity, monomer] = ["_entity_poly_seq.entity_id", "_entity_poly_seq.mon_id"]
            .map(|tag| table.value(row, tag).unwrap().text().to_owned());
        match sequences.iter_mut().find(|(known, _)| *known == entity) {
            Some((_, monomers)) => monomers.push(monomer),
            None => sequences.push((entity, vec![monomer])),
        }
    }
    sequences
        .into_iter()
        .map(|(_, monomers)| monomers)
        .collect()
}

/// The formal charge and hydrogens of each atom of residues 1001 and 1002,
/// by residue and atom name.
fn dipeptide_chemistry(structure: &Structure) -> Vec<(String, String, i8, Option<usize>)> {
    let mut atoms = structure
        .topology()
        .atom_sites()
        .filter(|site| matches!(site.residue().author_seq_id(), Some("1001" | "1002")))
        .map(|site| {
            (
                site.residue().author_seq_id().unwrap().to_owned(),
                site.metadata().label_atom_id().unwrap().to_owned(),
                site.atom().atom().formal_charge(),
                site.atom().total_hydrogens(),
            )
        })
        .collect::<Vec<_>>();
    atoms.sort();
    atoms
}

#[test]
fn a_peptide_added_to_a_polymer_chain_takes_its_own_sequence() {
    // Rubredoxin's 54 monomers are the source sequence of chain A. A
    // dipeptide added to chain A reads back as the same one added to a
    // chain of its own.
    let source = mmcif_fixture("1IRO.cif");
    let written = |chain| {
        let mut editor = source.structure().to_editor();
        add_dipeptide(&mut editor, chain);
        write_edit(&source, editor, true)
    };
    let (in_chain, alone) = (written("A"), written("P"));
    let sequences = polymer_sequences(&in_chain);
    assert_eq!(sequences.len(), 2, "{in_chain}");
    assert_eq!(sequences[0].len(), 54);
    assert_eq!(sequences[1], ["GLY", "GLY"]);
    assert_eq!(
        dipeptide_chemistry(&read_back(&in_chain)),
        dipeptide_chemistry(&read_back(&alone))
    );
}

/// One carbon atom per component, as a residue of chain A, bonded as
/// `bonds` says.
fn carbohydrate_structure(components: &[&str], bonds: &[(usize, usize)]) -> Structure {
    let mut graph = MoleculeEditor::new();
    let atoms = components
        .iter()
        .map(|_| {
            graph
                .add_atom(Atom::new(Element::from_symbol("C").unwrap()))
                .unwrap()
        })
        .collect::<Vec<_>>();
    for &(left, right) in bonds {
        graph
            .add_bond(atoms[left], atoms[right], BondOrder::Single)
            .unwrap();
    }
    let mut builder = StructureBuilder::new();
    let instance = builder
        .add_molecule(graph.finish().unwrap(), &Positions::zeros(atoms.len()))
        .unwrap();
    let hierarchy = builder.topology_mut().hierarchy_mut();
    let chain = hierarchy.add_chain("A", None).unwrap();
    for (index, (&component, atom)) in components.iter().zip(atoms).enumerate() {
        let sequence = i32::try_from(index).unwrap() + 1;
        let residue = hierarchy
            .add_residue(chain, component, Some(sequence), None, None)
            .unwrap();
        hierarchy
            .add_atom_site(
                residue,
                InstanceAtomId::new(instance, atom),
                AtomSiteMetadata::default().with_label_atom_id(Some("C1".to_owned())),
            )
            .unwrap();
    }
    builder.build().unwrap()
}

#[test]
fn saccharides_of_several_residues_are_branched_entities() {
    // PDBx writes every oligosaccharide as a branched entity, linear or not.
    let linear = carbohydrate_structure(&["GLC", "GAL", "GLC"], &[(0, 1), (1, 2)]);
    let graph_branched =
        carbohydrate_structure(&["GLC", "GAL", "GLC", "MAN"], &[(0, 1), (0, 2), (0, 3)]);
    for structure in [&linear, &graph_branched] {
        let (written, _) = mmcif::write([structure], MmcifWriteOptions::default()).unwrap();
        assert_eq!(entity_types(&written), ["branched"]);
    }
}

/// A trisaccharide whose block states it as a branched entity.
const BRANCHED_CARBOHYDRATE: &str = r#"
data_branched_carbohydrate
loop_
_entity.id
_entity.type
1 branched
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
HETATM 1 C C1 GLC A 1 1 0.0 0.0 0.0
HETATM 2 C C1 GAL A 1 2 3.0 0.0 0.0
HETATM 3 C C1 GLC A 1 3 6.0 0.0 0.0
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_label_seq_id
link-1 covale A C1 1 A C1 2
link-2 covale A C1 2 A C1 3
"#;

#[test]
fn a_source_branched_saccharide_stays_branched() {
    let interpreted = mmcif::parse_str(BRANCHED_CARBOHYDRATE, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    let structure = interpreted.structure();
    let molecule = structure.topology().instances().next().unwrap();
    assert_eq!(molecule.class(), MoleculeClass::Carbohydrate);
    assert_eq!(molecule.residues().count(), 3);
    for written in [
        mmcif::write([structure], MmcifWriteOptions::default()),
        mmcif::write([&interpreted], MmcifWriteOptions::default()),
    ] {
        let (text, _) = written.unwrap();
        assert_eq!(entity_types(&text), ["branched"]);
    }
}

/// The block written without a report from the structure read back from
/// the PDB entry written from `source`.
fn via_pdb(source: &StructureInterpretation) -> String {
    let read = pdb_round_trip(source);
    let (text, _) = mmcif::write([read.structure()], MmcifWriteOptions::default()).unwrap();
    text
}

/// The block written with its report from the interpretation of the PDB
/// entry written from `source`.
fn via_pdb_with_report(source: &StructureInterpretation) -> String {
    let (text, _) = mmcif::write([&pdb_round_trip(source)], MmcifWriteOptions::default()).unwrap();
    text
}

/// The interpretation of the PDB entry written from `source`.
fn pdb_round_trip(source: &StructureInterpretation) -> StructureInterpretation {
    let (text, _) = pdb::write(source, Default::default()).unwrap();
    pdb::parse_str(&text, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap()
}

/// The values of `tags` in each atom-site row, with the row's entity type
/// last.
fn site_values<const N: usize>(text: &str, tags: [&str; N]) -> Vec<Vec<String>> {
    let document = mmcif::parse_str(text, Default::default()).unwrap();
    let block = &document.blocks()[0];
    let entities = block.loop_with_tag("_entity.id").unwrap();
    let types = (0..entities.row_count())
        .map(|row| {
            ["_entity.id", "_entity.type"].map(|tag| entities.value(row, tag).unwrap().text())
        })
        .map(|[id, kind]| (id.to_owned(), kind.to_owned()))
        .collect::<BTreeMap<_, _>>();
    let table = block.loop_with_tag("_atom_site.id").unwrap();
    (0..table.row_count())
        .map(|row| {
            let mut values = tags
                .iter()
                .map(|tag| table.value(row, tag).unwrap().text().to_owned())
                .collect::<Vec<_>>();
            let entity = table
                .value(row, "_atom_site.label_entity_id")
                .unwrap()
                .text();
            values.push(types[entity].clone());
            values
        })
        .collect()
}

#[test]
fn covalently_bound_ligands_are_instances_of_their_own() {
    // Sotorasib (MOV) binds cysteine 12 of KRas, and retinal (RET) a lysine
    // of its opsin; neither is a monomer of the chain.
    for (source, ligand) in [
        (mmcif_fixture("6OIM.cif"), "MOV"),
        (curated("pdb_00006g7l_xyz_v1-4.cif.gz"), "RET"),
    ] {
        let text = via_pdb(&source);
        let rows = site_values(
            &text,
            ["_atom_site.label_comp_id", "_atom_site.label_asym_id"],
        );
        let polymer_asyms = rows
            .iter()
            .filter(|row| row[2] == "polymer")
            .map(|row| row[1].clone())
            .collect::<BTreeSet<_>>();
        let ligand_asyms = rows
            .iter()
            .filter(|row| row[0] == ligand)
            .map(|row| (row[1].clone(), row[2].clone()))
            .collect::<BTreeSet<_>>();
        let [(asym, kind)] = ligand_asyms
            .into_iter()
            .collect::<Vec<_>>()
            .try_into()
            .unwrap();
        assert_eq!(kind, "non-polymer", "{ligand}");
        assert!(!polymer_asyms.contains(&asym), "{ligand}");
    }
}

/// The asyms of a block that hold saccharides `names`: each one's entity
/// type and components.
fn glycans(text: &str, names: &[&str]) -> Vec<(String, Vec<String>)> {
    let mut asyms = BTreeMap::<String, (String, BTreeSet<(String, String)>)>::new();
    for row in site_values(
        text,
        [
            "_atom_site.label_comp_id",
            "_atom_site.label_asym_id",
            "_atom_site.auth_seq_id",
        ],
    ) {
        if names.contains(&row[0].as_str()) {
            asyms
                .entry(row[1].clone())
                .or_insert_with(|| (row[3].clone(), BTreeSet::new()))
                .1
                .insert((row[2].clone(), row[0].clone()));
        }
    }
    let mut glycans = asyms
        .into_values()
        .map(|(kind, residues)| {
            let mut components = residues
                .into_iter()
                .map(|(_, name)| name)
                .collect::<Vec<_>>();
            components.sort();
            (kind, components)
        })
        .collect::<Vec<_>>();
    glycans.sort();
    glycans
}

#[test]
fn glycans_bound_to_a_protein_are_branched_entities() {
    // 7UKO's N-linked glycans of two and more saccharides are branched
    // entities, and its single N-acetylglucosamines non-polymer ones; 7X9L's
    // sucrose, bound to no protein, is branched too.
    for (name, saccharides) in [
        ("pdb_00007uko_xyz_v1-4.cif.gz", &["NAG", "BMA", "MAN"][..]),
        ("pdb_00007x9l_xyz_v1-1.cif.gz", &["GLC", "FRU"][..]),
    ] {
        let source = curated(name);
        let (direct, _) = mmcif::write([&source], MmcifWriteOptions::default()).unwrap();
        let expected = glycans(&direct, saccharides);
        assert!(expected.iter().any(|(kind, _)| kind == "branched"));
        assert_eq!(glycans(&via_pdb(&source), saccharides), expected, "{name}");
        assert_eq!(
            glycans(&via_pdb_with_report(&source), saccharides),
            expected,
            "{name} with the report"
        );
    }
}

#[test]
fn instances_of_a_named_component_are_one_entity_however_much_they_model() {
    // 6G7L models seven monoolein lipids (LI1), each in part.
    let text = via_pdb(&curated("pdb_00006g7l_xyz_v1-4.cif.gz"));
    let entities = component_rows(&text, "LI1")
        .into_iter()
        .map(|(_, entity)| entity)
        .collect::<BTreeSet<_>>();
    assert_eq!(entities.len(), 1, "{entities:?}");
}

#[test]
fn written_polymer_residues_state_their_places() {
    let directory = format!(
        "{}/../../benchmarks/datasets/bio/curated/pdb",
        env!("CARGO_MANIFEST_DIR")
    );
    for entry in std::fs::read_dir(directory).unwrap() {
        let name = entry.unwrap().file_name().into_string().unwrap();
        // 6GCY's residues pair two groups of alternate locations, which the
        // default selection cannot combine.
        if name.starts_with("pdb_00006gcy") {
            continue;
        }
        let text = via_pdb(&curated(&name));
        let unplaced = site_values(
            &text,
            [
                "_atom_site.label_asym_id",
                "_atom_site.label_comp_id",
                "_atom_site.label_seq_id",
            ],
        )
        .into_iter()
        .filter(|row| row[3] == "polymer" && row[2] == ".")
        .map(|row| (row[0].clone(), row[1].clone()))
        .collect::<BTreeSet<_>>();
        assert!(unplaced.is_empty(), "{name}: {unplaced:?}");
    }
}

/// An atom's site: its author chain, component, author number, insertion
/// code, and name; blank for an atom outside the hierarchy.
type SiteKey = [String; 5];

fn site_key(atom: kekule::topology::AtomView<'_>) -> SiteKey {
    let Some(site) = atom.atom_site() else {
        return SiteKey::default();
    };
    let metadata = site.metadata();
    let residue = site.residue();
    let chain = residue.chain();
    [
        metadata
            .auth_asym_id()
            .or(chain.author_id())
            .unwrap_or(chain.label_id()),
        residue.name(),
        residue.author_seq_id().unwrap_or_default(),
        residue.insertion_code().unwrap_or_default(),
        metadata
            .label_atom_id()
            .or(metadata.auth_atom_id())
            .unwrap_or_default(),
    ]
    .map(str::to_owned)
}

/// A structure's bonds by their atoms' sites and orders, each with the
/// number of bonds that share them.
fn site_bonds(structure: &Structure) -> BTreeMap<(SiteKey, SiteKey, String), usize> {
    let mut bonds = BTreeMap::new();
    for bond in structure.topology().bonds() {
        let [mut left, mut right] = bond.atoms().map(site_key);
        if left > right {
            std::mem::swap(&mut left, &mut right);
        }
        *bonds
            .entry((left, right, format!("{:?}", bond.bond().order())))
            .or_insert(0) += 1;
    }
    bonds
}

#[test]
fn every_bond_survives_a_write_without_the_report() {
    // Strand B of 7B3K is a peptide of D-amino acids, and strand B of 1K1H
    // a DNA of modified nucleotides: no residue class names their monomers.
    let directory = format!(
        "{}/../../benchmarks/datasets/bio/curated/pdb",
        env!("CARGO_MANIFEST_DIR")
    );
    let mut entries = std::fs::read_dir(directory)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().into_string().unwrap())
        // 6GCY's residues pair two groups of alternate locations, which the
        // default selection cannot combine.
        .filter(|name| !name.starts_with("pdb_00006gcy"))
        .map(|name| (curated(&name), name))
        .collect::<Vec<_>>();
    for name in ["7B3K-B.cif", "1K1H.cif"] {
        entries.push((mmcif_fixture(name), name.to_owned()));
    }
    for (source, name) in entries {
        let (text, _) = mmcif::write([source.structure()], MmcifWriteOptions::default()).unwrap();
        let (before, after) = (
            site_bonds(source.structure()),
            site_bonds(&read_back(&text)),
        );
        let lost = before
            .iter()
            .filter(|(bond, count)| after.get(*bond) != Some(count))
            .take(5)
            .collect::<Vec<_>>();
        let gained = after
            .iter()
            .filter(|(bond, count)| before.get(*bond) != Some(count))
            .take(5)
            .collect::<Vec<_>>();
        assert!(
            lost.is_empty() && gained.is_empty(),
            "{name}: lost {lost:?}, gained {gained:?}"
        );
        // Only polymer residues state places.
        let placed = site_values(
            &text,
            ["_atom_site.label_comp_id", "_atom_site.label_seq_id"],
        )
        .into_iter()
        .filter(|row| row[2] != "polymer" && row[1] != ".")
        .map(|row| row[0].clone())
        .collect::<BTreeSet<_>>();
        assert!(placed.is_empty(), "{name}: {placed:?}");
    }
}

#[test]
fn chains_of_monomers_no_residue_class_names_are_polymers() {
    for (name, monomers) in [
        (
            "7B3K-B.cif",
            &["DAR", "DPR", "2TL", "DLE", "DHI", "DSG"][..],
        ),
        ("1K1H.cif", &["DC", "CMR", "RMP", "SMP"][..]),
    ] {
        let source = mmcif_fixture(name);
        let (text, _) = mmcif::write([source.structure()], MmcifWriteOptions::default()).unwrap();
        let kinds = site_values(&text, ["_atom_site.label_comp_id"])
            .into_iter()
            .filter(|row| monomers.contains(&row[0].as_str()))
            .map(|row| row[1].clone())
            .collect::<BTreeSet<_>>();
        assert_eq!(kinds, ["polymer".to_owned()].into(), "{name}");
    }
}

#[test]
fn a_chain_the_caller_classifies_as_non_polymer_keeps_its_links() {
    let source = mmcif_fixture("7B3K-B.cif");
    let structure = source.structure();
    let mut classifications = MmcifEntityClassifications::new();
    for instance in structure.topology().instances() {
        classifications
            .insert(instance.id(), EntityKind::NonPolymer)
            .unwrap();
    }
    let (text, _) = mmcif::write(
        [MmcifBlockSource::structure(structure).with_classifications(&classifications)],
        MmcifWriteOptions::default(),
    )
    .unwrap();
    assert_eq!(site_bonds(&read_back(&text)), site_bonds(structure));
    let places = site_values(&text, ["_atom_site.label_seq_id"])
        .into_iter()
        .map(|row| row[0].clone())
        .collect::<BTreeSet<_>>();
    assert_eq!(places, [".".to_owned()].into());
}

#[test]
fn a_residue_in_a_polymer_asym_without_a_place_fails_the_write() {
    // A free glycine states asym A, rubredoxin's chain, but no place in its
    // sequence.
    let source = mmcif_fixture("1IRO.cif");
    let mut editor = source.structure().to_editor();
    let element = |symbol| Atom::new(Element::from_symbol(symbol).unwrap());
    let mut glycine = MoleculeEditor::new();
    let atoms = [("N", 1), ("C", 0), ("C", 0), ("O", 0), ("O", -1)].map(|(symbol, charge)| {
        glycine
            .add_atom(element(symbol).with_formal_charge(charge))
            .unwrap()
    });
    for (left, right, order) in [
        (0, 1, BondOrder::Single),
        (1, 2, BondOrder::Single),
        (2, 3, BondOrder::Double),
        (2, 4, BondOrder::Single),
    ] {
        glycine.add_bond(atoms[left], atoms[right], order).unwrap();
    }
    let points = [
        [0.0, 0.0],
        [1.458, 0.0],
        [2.009, 1.42],
        [1.251, 2.39],
        [3.25, 1.5],
    ]
    .map(|[x, y]| Point3::new(x + 60.0, y, 0.0));
    let positions = Positions::new(Quantity::new(points.to_vec(), ANGSTROM)).unwrap();
    let molecule = editor
        .add_molecule(glycine.finish().unwrap(), &positions)
        .unwrap();
    let chain = editor
        .topology()
        .chains()
        .find(|(_, chain)| chain.label_id() == "A")
        .map(|(id, _)| id)
        .unwrap();
    let residue = editor
        .add_residue(chain, "GLY", None, Some("1001".to_owned()), None)
        .unwrap();
    for (index, name) in (0u32..).zip(["N", "CA", "C", "O", "OXT"]) {
        editor
            .add_atom_site(
                residue,
                molecule.atoms()[&AtomId::new(index)],
                AtomSiteMetadata::default()
                    .with_label_atom_id(Some(name.to_owned()))
                    .with_label_asym_id(Some("A".to_owned())),
            )
            .unwrap();
    }
    let edited = editor.finish().unwrap();
    let error = mmcif::write([&edited], MmcifWriteOptions::default()).unwrap_err();
    assert!(
        matches!(
            error.detail(),
            Some(FormatErrorDetail::MmcifWrite(
                MmcifWriteErrorKind::ConflictingAsymEntityClassifications { .. }
            ))
        ),
        "{error:?}"
    );
}

/// A block whose asym A holds a glycine and a ligand, with a ligand asym
/// and a water.
const GLYCINE_AND_LIGAND: &str = "data_mixed
loop_
_entity.id
_entity.type
1 polymer
2 non-polymer
3 water
loop_
_struct_asym.id
_struct_asym.entity_id
A 1
L 2
W 3
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
ATOM 1 N N GLY A 1 1 0.0 0.0 0.0 1
HETATM 2 C C1 LIG A 1 3 9.0 0.0 0.0 1
HETATM 3 C C1 LIG L 2 . 2.0 0.0 0.0 1
HETATM 4 O O HOH W 3 . 3.0 0.0 0.0 1
";

#[test]
fn a_partial_classification_still_conflicts_within_one_asym() {
    // The caller classifies the glycine a polymer; the ligand in its asym
    // stays non-polymer, so the asym conflicts.
    let interpreted = mmcif::parse_str(GLYCINE_AND_LIGAND, Default::default())
        .unwrap()
        .interpret(Default::default())
        .unwrap();
    let structure = interpreted.structure();
    let mut instances = structure
        .topology()
        .instances()
        .map(|instance| instance.id());
    let classified = instances.next().unwrap();
    instances.next().unwrap();
    let mut classifications = MmcifEntityClassifications::new();
    classifications
        .insert(classified, EntityKind::Polymer)
        .unwrap();
    let error = mmcif::write(
        [MmcifBlockSource::structure(structure).with_classifications(&classifications)],
        MmcifWriteOptions::default(),
    )
    .unwrap_err();
    assert!(
        matches!(
            error.detail(),
            Some(FormatErrorDetail::MmcifWrite(
                MmcifWriteErrorKind::ConflictingAsymEntityClassifications { .. }
            ))
        ),
        "{error:?}"
    );
}
