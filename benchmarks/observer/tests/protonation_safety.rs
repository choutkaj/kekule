//! Protonation safety net: every curated structure that Kekule's tests and
//! benchmarks carry, and a heme with every hydrogen, reads as it is and
//! with its hydrogens stripped, and both structures perceive. Without
//! hydrogens, the standard residues take the net charge their pH 7 states
//! give, counted here from the residues' atoms alone. As it is, a structure
//! keeps its net charge through a PDB round trip.

use std::collections::{BTreeMap, BTreeSet};
use std::io::Read;
use std::path::{Path, PathBuf};

use kekule::macromolecular::{
    ModelSelection, SourceResiduePosition, StructureInterpretIssue, StructureInterpretOptions,
    StructureInterpretation,
};
use kekule::topology::{InstanceAtomId, ResidueView, Topology};
use kekule::{mmcif, pdb};

/// A structure file and whether it is PDB rather than mmCIF.
struct Entry {
    path: PathBuf,
    pdb: bool,
}

/// Heme with every hydrogen of its template, at the template's ideal
/// coordinates, as an mmCIF block: its carboxylic acids are protonated.
fn protonated_heme() -> String {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let template =
        std::fs::read_to_string(root.join("../../crates/kekule/data/components/HEM.cif"))
            .expect("bundled heme template");
    let header = template
        .lines()
        .filter(|line| line.starts_with("_chem_comp_atom."))
        .collect::<Vec<_>>();
    let column = |name: &str| {
        header
            .iter()
            .position(|line| line.trim_end() == name)
            .expect("template column")
    };
    let (atom, symbol, x) = (
        column("_chem_comp_atom.atom_id"),
        column("_chem_comp_atom.type_symbol"),
        column("_chem_comp_atom.pdbx_model_Cartn_x_ideal"),
    );
    let mut text = "data_HEM\nloop_\n_atom_site.group_PDB\n_atom_site.id\n\
                    _atom_site.type_symbol\n_atom_site.label_atom_id\n\
                    _atom_site.label_comp_id\n_atom_site.label_asym_id\n\
                    _atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n\
                    _atom_site.Cartn_z\n_atom_site.auth_seq_id\n"
        .to_owned();
    let rows = template
        .lines()
        .map(tokens)
        .filter(|fields| fields.len() > x + 2 && fields[0] == "HEM");
    for (serial, fields) in (1..).zip(rows) {
        text.push_str(&format!(
            "HETATM {serial} {} \"{}\" HEM A . {} {} {} 1\n",
            fields[symbol],
            fields[atom],
            fields[x],
            fields[x + 1],
            fields[x + 2]
        ));
    }
    text
}

fn entries() -> Vec<Entry> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut entries = Vec::new();
    for (directory, pdb) in [
        ("../datasets/bio/curated/pdb", false),
        ("../../crates/kekule/tests/fixtures/mmcif", false),
        ("../../crates/kekule/tests/fixtures/pdb", true),
    ] {
        let mut paths = std::fs::read_dir(root.join(directory))
            .expect("readable fixture directory")
            .map(|entry| entry.expect("directory entry").path())
            .filter(|path| {
                let name = path.to_string_lossy();
                name.ends_with(".cif.gz") || name.ends_with(".cif") || name.ends_with(".pdb")
            })
            .collect::<Vec<_>>();
        paths.sort();
        entries.extend(paths.into_iter().map(|path| Entry { path, pdb }));
    }
    entries
}

fn text(path: &Path) -> String {
    if path.to_string_lossy().ends_with(".gz") {
        let mut text = String::new();
        flate2::read::GzDecoder::new(std::fs::File::open(path).expect("structure file"))
            .read_to_string(&mut text)
            .expect("gzip text");
        text
    } else {
        std::fs::read_to_string(path).expect("structure file")
    }
}

/// The whitespace-separated tokens of an mmCIF data line, quotes removed.
fn tokens(line: &str) -> Vec<&str> {
    let mut tokens = Vec::new();
    let bytes = line.as_bytes();
    let mut index = 0;
    while index < bytes.len() {
        if bytes[index].is_ascii_whitespace() {
            index += 1;
            continue;
        }
        let quote = bytes[index];
        if quote == b'\'' || quote == b'"' {
            let start = index + 1;
            let mut end = start;
            while end < bytes.len()
                && !(bytes[end] == quote
                    && bytes
                        .get(end + 1)
                        .is_none_or(|next| next.is_ascii_whitespace()))
            {
                end += 1;
            }
            tokens.push(&line[start..end.min(bytes.len())]);
            index = end + 1;
        } else {
            let start = index;
            while index < bytes.len() && !bytes[index].is_ascii_whitespace() {
                index += 1;
            }
            tokens.push(&line[start..index]);
        }
    }
    tokens
}

fn is_hydrogen_symbol(symbol: &str) -> bool {
    matches!(symbol.trim(), "H" | "D" | "T")
}

/// The mmCIF text without its hydrogen atom-site rows.
fn strip_mmcif(text: &str) -> String {
    let lines = text.lines().collect::<Vec<_>>();
    let mut kept = Vec::with_capacity(lines.len());
    let mut index = 0;
    while index < lines.len() {
        kept.push(lines[index]);
        if lines[index].trim() != "loop_" {
            index += 1;
            continue;
        }
        index += 1;
        let mut header = Vec::new();
        while index < lines.len() && lines[index].starts_with('_') {
            header.push(lines[index].trim());
            kept.push(lines[index]);
            index += 1;
        }
        let symbol = header
            .iter()
            .position(|item| *item == "_atom_site.type_symbol");
        while index < lines.len()
            && !lines[index].starts_with('_')
            && !lines[index].starts_with('#')
            && lines[index].trim() != "loop_"
            && !lines[index].starts_with("data_")
        {
            let hydrogen = symbol.is_some_and(|column| {
                tokens(lines[index])
                    .get(column)
                    .is_some_and(|value| is_hydrogen_symbol(value))
            });
            if !hydrogen {
                kept.push(lines[index]);
            }
            index += 1;
        }
    }
    kept.join("\n") + "\n"
}

/// The PDB text without its hydrogen coordinate records.
fn strip_pdb(text: &str) -> String {
    text.lines()
        .filter(|line| {
            let coordinate = line.starts_with("ATOM  ") || line.starts_with("HETATM");
            !(coordinate && line.get(76..78).is_some_and(is_hydrogen_symbol))
        })
        .collect::<Vec<_>>()
        .join("\n")
        + "\n"
}

/// The first model of the entry.
fn read(entry: &Entry, text: &str) -> Result<StructureInterpretation, String> {
    let first = || StructureInterpretOptions::default().with_model_selection(ModelSelection::First);
    let result = if entry.pdb {
        pdb::parse_str(text, Default::default()).and_then(|document| document.interpret(first()))
    } else {
        mmcif::parse_str(text, Default::default()).and_then(|document| document.interpret(first()))
    };
    result.map_err(|error| error.to_string())
}

/// The residues, as chain and label sequence number, for which the source
/// states a charge that differs from the derived one. The source's charge
/// holds there, so they have no pH 7 expectation.
fn stated_charges(result: &StructureInterpretation) -> BTreeSet<(String, i32)> {
    result
        .report()
        .issues()
        .iter()
        .filter_map(|issue| match issue {
            StructureInterpretIssue::StatedChargeKept { residue, .. }
            | StructureInterpretIssue::StatedChargeRejected { residue, .. }
            | StructureInterpretIssue::StatedChargeResonance { residue, .. } => {
                match residue.position() {
                    SourceResiduePosition::Label(sequence) => {
                        Some((residue.asym_id().to_owned(), *sequence))
                    }
                    _ => None,
                }
            }
            _ => None,
        })
        .collect()
}

const AMINO_ACIDS: [&str; 20] = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET",
    "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
];
const NUCLEOTIDES: [&str; 8] = ["A", "C", "G", "U", "DA", "DC", "DG", "DT"];

/// One atom of a residue: its charge, whether it coordinates a metal, and
/// whether it bonds to another residue.
struct ResidueAtom {
    charge: i8,
    coordinated: bool,
    linked: bool,
}

fn residue_atoms(
    residue: ResidueView<'_>,
    coordinated: &BTreeSet<InstanceAtomId>,
) -> BTreeMap<String, ResidueAtom> {
    residue
        .atom_sites()
        .map(|site| {
            let atom = site.atom();
            let name = site
                .metadata()
                .label_atom_id()
                .or(site.metadata().auth_atom_id())
                .unwrap_or_default()
                .to_owned();
            let linked = atom.neighbors().any(|other| {
                other
                    .atom_site()
                    .is_some_and(|other| other.residue().id() != residue.id())
            });
            (
                name,
                ResidueAtom {
                    charge: atom.atom().formal_charge(),
                    coordinated: coordinated.contains(&atom.id()),
                    linked,
                },
            )
        })
        .collect()
}

/// The net charge of a standard residue without hydrogens at pH 7, from
/// its atoms: acids anionic, lysine and arginine cationic, histidine
/// neutral, a true N terminus an ammonium and a C terminus a carboxylate.
/// Metal coordination deprotonates a cysteine or tyrosine and takes the
/// place of an amine's added proton.
fn expected_charge(residue: ResidueView<'_>, atoms: &BTreeMap<String, ResidueAtom>) -> i32 {
    let has = |name: &str| atoms.contains_key(name);
    let coordinated = |name: &str| atoms.get(name).is_some_and(|atom| atom.coordinated);
    let name = residue.name();
    let mut charge = 0;
    if AMINO_ACIDS.contains(&name) {
        charge += match name {
            "ASP" if has("OD2") => -1,
            "GLU" if has("OE2") => -1,
            "LYS" if has("NZ") && !coordinated("NZ") => 1,
            "ARG" if has("NH2") => 1,
            "CYS" if coordinated("SG") => -1,
            "TYR" if coordinated("OH") => -1,
            _ => 0,
        };
        let free_n = atoms
            .get("N")
            .is_some_and(|atom| !atom.linked && !atom.coordinated);
        if residue.label_seq_id() == Some(1) && free_n {
            charge += 1;
        }
        if has("OXT") {
            charge -= 1;
        }
    } else {
        for oxygen in ["OP2", "OP3"] {
            if atoms.get(oxygen).is_some_and(|atom| !atom.linked) {
                charge -= 1;
            }
        }
    }
    charge
}

/// The standard residues whose charge differs from their pH 7 one, as
/// `chain residue name: charge (expected)`.
fn charge_mismatches(result: &StructureInterpretation) -> Vec<String> {
    let topology: &Topology = result.structure().topology();
    let coordinated = result
        .report()
        .coordination()
        .iter()
        .flat_map(|record| record.atoms())
        .collect::<BTreeSet<_>>();
    let stated = stated_charges(result);
    let mut mismatches = Vec::new();
    for chain in topology.chains() {
        for residue in chain.residues() {
            if !AMINO_ACIDS.contains(&residue.name()) && !NUCLEOTIDES.contains(&residue.name()) {
                continue;
            }
            let key = (
                chain.label_id().to_owned(),
                residue.label_seq_id().unwrap_or(0),
            );
            if stated.contains(&key) {
                continue;
            }
            let atoms = residue_atoms(residue, &coordinated);
            let charge = atoms
                .values()
                .map(|atom| i32::from(atom.charge))
                .sum::<i32>();
            let expected = expected_charge(residue, &atoms);
            if charge != expected {
                mismatches.push(format!(
                    "{} {} {}: {charge} ({expected})",
                    chain.label_id(),
                    residue.author_seq_id().unwrap_or_default(),
                    residue.name()
                ));
            }
        }
    }
    mismatches
}

/// The net formal charge of a structure.
fn net_charge(result: &StructureInterpretation) -> i32 {
    result
        .structure()
        .topology()
        .atoms()
        .map(|atom| i32::from(atom.atom().formal_charge()))
        .sum()
}

/// How a PDB round trip changes the net charge, if it does.
fn round_trip_charge(result: &StructureInterpretation) -> Option<String> {
    let (text, _) = match pdb::write(result, Default::default()) {
        Ok(written) => written,
        // PDB has no room for five-character component names.
        Err(error) if error.to_string().contains("residue-name columns") => return None,
        Err(error) => return Some(format!("does not write as PDB: {error}")),
    };
    let entry = Entry {
        path: PathBuf::new(),
        pdb: true,
    };
    let reread = match read(&entry, &text) {
        Ok(reread) => reread,
        Err(error) => return Some(format!("does not read back from PDB: {error}")),
    };
    let (before, after) = (net_charge(result), net_charge(&reread));
    (before != after).then(|| format!("net charge {before} reads back from PDB as {after}"))
}

#[test]
fn curated_structures_perceive_with_and_without_hydrogens() {
    let mut failures = Vec::new();
    let mut sources = entries()
        .into_iter()
        .map(|entry| {
            let label = entry
                .path
                .file_name()
                .expect("file name")
                .to_string_lossy()
                .into_owned();
            let source = text(&entry.path);
            (label, entry, source)
        })
        .collect::<Vec<_>>();
    sources.push((
        "protonated heme".to_owned(),
        Entry {
            path: PathBuf::new(),
            pdb: false,
        },
        protonated_heme(),
    ));
    for (label, entry, source) in sources {
        let stripped = if entry.pdb {
            strip_pdb(&source)
        } else {
            strip_mmcif(&source)
        };
        for (version, text) in [("as is", &source), ("without hydrogens", &stripped)] {
            let result = match read(&entry, text) {
                Ok(result) => result,
                // Alternate-location selection rejects some entries before
                // protonation; that is no protonation failure.
                Err(error) if error.contains("alternate-location combination") => continue,
                Err(error) => {
                    failures.push(format!("{label} {version}: does not read: {error}"));
                    continue;
                }
            };
            let mut structure = result.structure().clone();
            if let Err(error) = structure.perceive() {
                failures.push(format!("{label} {version}: does not perceive: {error}"));
            }
            if version == "as is" {
                if let Some(failure) = round_trip_charge(&result) {
                    failures.push(format!("{label}: {failure}"));
                }
            }
            if version == "without hydrogens" {
                let mismatches = charge_mismatches(&result);
                if !mismatches.is_empty() {
                    failures.push(format!(
                        "{label}: standard residues off their pH 7 charge: {}",
                        mismatches.join(", ")
                    ));
                }
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
