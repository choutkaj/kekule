//! Fixture access shared by the integration test binaries.
#![allow(dead_code)]

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use kekule::core::{
    Atom, AtomId, BondOrder, Element, Molecule, MoleculeEditor, StereoCarrier, StereoElement,
    StereoElementKind, StereoGroup,
};
use kekule::{sdf, smiles, stereo};

/// Absolute path of a file under `tests/fixtures`.
pub fn fixture(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(relative)
}

/// Absolute path of a file under the repository's `benchmarks/corpora`.
pub fn corpus(relative: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../benchmarks/corpora")
        .join(relative)
}

/// Interprets one connected molecule from SMILES or from the first record of
/// an SDF file, written as `sdf:<path>` relative to `directory`. Perception is
/// not run.
pub fn read_molecule(input: &str, directory: &Path) -> Molecule {
    let mut molecules = match input.strip_prefix("sdf:") {
        Some(file) => {
            let text = std::fs::read_to_string(directory.join(file))
                .unwrap_or_else(|error| panic!("{file}: {error}"));
            let document = sdf::parse_str(&text).unwrap_or_else(|error| panic!("{file}: {error}"));
            document.records()[0]
                .interpret()
                .unwrap_or_else(|error| panic!("{file}: {error}"))
                .into_molecules()
        }
        None => smiles::to_molecules(input).unwrap_or_else(|error| panic!("{input}: {error}")),
    };
    assert_eq!(molecules.len(), 1, "{input} must be one connected molecule");
    molecules.pop().unwrap()
}

/// One row of `fixtures/perception/aromaticity.tsv`.
pub struct AromaticityCase {
    pub label: String,
    pub input: String,
    pub aromatic_atoms: Vec<usize>,
    pub nonaromatic_bonds: Vec<(usize, usize)>,
}

impl AromaticityCase {
    pub fn molecule(&self) -> Molecule {
        read_molecule(&self.input, &fixture("perception"))
    }
}

pub fn aromaticity_cases() -> Vec<AromaticityCase> {
    let path = fixture("perception/aromaticity.tsv");
    std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()))
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            assert_eq!(fields.len(), 4, "malformed aromaticity row: {line}");
            AromaticityCase {
                label: fields[0].to_owned(),
                input: fields[1].to_owned(),
                aromatic_atoms: index_ranges(fields[2]),
                nonaromatic_bonds: index_pairs(fields[3]),
            }
        })
        .collect()
}

/// Parses `0-5,7` into `[0, 1, 2, 3, 4, 5, 7]`; `-` is empty.
fn index_ranges(field: &str) -> Vec<usize> {
    if field == "-" {
        return Vec::new();
    }
    field
        .split(',')
        .flat_map(|part| {
            let (start, end) = part.split_once('-').unwrap_or((part, part));
            start.parse::<usize>().unwrap()..=end.parse::<usize>().unwrap()
        })
        .collect()
}

/// Parses `3:12,4:6` into `[(3, 12), (4, 6)]`; `-` is empty.
fn index_pairs(field: &str) -> Vec<(usize, usize)> {
    if field == "-" {
        return Vec::new();
    }
    field
        .split(',')
        .map(|pair| {
            let (left, right) = pair.split_once(':').unwrap();
            (left.parse().unwrap(), right.parse().unwrap())
        })
        .collect()
}

/// Rebuilds `molecule` with a seeded atom permutation, reversed bond order and
/// swapped bond endpoints, remapping every stereo element and group. Odd seeds
/// also leave deleted draft slots behind so sparse storage cannot leak into
/// numbering-dependent results. Returns the copy and, for every source atom
/// index, the corresponding atom of the copy. Perception is not copied.
pub fn renumbered(molecule: &Molecule, seed: u64) -> (Molecule, Vec<AtomId>) {
    let mut order = molecule.atom_ids().collect::<Vec<_>>();
    let mut state = seed;
    for index in (1..order.len()).rev() {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        order.swap(index, (state >> 33) as usize % (index + 1));
    }
    let mut editor = MoleculeEditor::new();
    if seed % 2 == 1 {
        let carbon = Atom::new(Element::from_symbol("C").unwrap());
        let left = editor.add_atom(carbon.clone()).unwrap();
        let right = editor.add_atom(carbon).unwrap();
        editor.add_bond(left, right, BondOrder::Single).unwrap();
        editor.delete_atom(left).unwrap();
        editor.delete_atom(right).unwrap();
    }
    let mut atoms = BTreeMap::new();
    for old in order {
        let new = editor
            .add_atom(molecule.atom(old).unwrap().clone())
            .unwrap();
        atoms.insert(old, new);
    }
    let mut bonds = BTreeMap::new();
    for (old, bond) in molecule.bonds().collect::<Vec<_>>().into_iter().rev() {
        let new = editor
            .add_bond(atoms[&bond.b()], atoms[&bond.a()], bond.order)
            .unwrap();
        bonds.insert(old, new);
    }
    let carrier = |carrier: StereoCarrier| match carrier {
        StereoCarrier::Atom(id) => StereoCarrier::Atom(atoms[&id]),
        other => other,
    };
    let mut elements = BTreeMap::new();
    for (old, element) in molecule.stereo_elements() {
        let mut kind = element.kind.clone();
        match &mut kind {
            StereoElementKind::Tetrahedral(value) => {
                value.center = atoms[&value.center];
                value.carriers = value.carriers.iter().copied().map(carrier).collect();
            }
            StereoElementKind::DoubleBond(value) => {
                value.bond = bonds[&value.bond];
                value.left = atoms[&value.left];
                value.right = atoms[&value.right];
                value.left_carrier = carrier(value.left_carrier);
                value.right_carrier = carrier(value.right_carrier);
            }
            StereoElementKind::Axis(value) => {
                value.axis = bonds[&value.axis];
                value.carriers = value.carriers.iter().copied().map(carrier).collect();
            }
        }
        let new = editor.add_stereo_element(StereoElement::new(kind)).unwrap();
        elements.insert(old, new);
    }
    for (_, group) in molecule.stereo_groups() {
        editor
            .add_stereo_group(StereoGroup {
                kind: group.kind,
                members: group.members.iter().map(|id| elements[id]).collect(),
            })
            .unwrap();
    }
    let (copy, published) = editor.finish_with_correspondence().unwrap();
    let mapping = molecule
        .atom_ids()
        .map(|old| published.atom(atoms[&old]).unwrap())
        .collect();
    (copy, mapping)
}

/// Where a stereo descriptor applies: a center atom, or the endpoints of a
/// stereogenic double bond or axis as sorted atom indices.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum StereoFocus {
    Atom(usize),
    Bond(usize, usize),
}

impl StereoFocus {
    pub fn mapped(self, mapping: &[AtomId]) -> Self {
        self.try_map(|atom| Some(mapping[atom].index())).unwrap()
    }

    /// Maps every focus atom index, or returns `None` if one has no image.
    pub fn try_map(self, atom: impl Fn(usize) -> Option<usize>) -> Option<Self> {
        Some(match self {
            Self::Atom(center) => Self::Atom(atom(center)?),
            Self::Bond(left, right) => Self::bond(atom(left)?, atom(right)?),
        })
    }

    fn bond(left: usize, right: usize) -> Self {
        Self::Bond(left.min(right), left.max(right))
    }
}

/// Assigns CIP descriptors to a perceived molecule and keys them by focus.
pub fn cip_labels(molecule: &mut Molecule) -> BTreeMap<StereoFocus, String> {
    try_cip_labels(molecule).unwrap_or_else(|error| panic!("CIP assignment failed: {error}"))
}

/// Like [`cip_labels`], but reports a failed assignment instead of panicking.
pub fn try_cip_labels(molecule: &mut Molecule) -> Result<BTreeMap<StereoFocus, String>, String> {
    let report = stereo::assign_cip_descriptors(molecule).map_err(|error| format!("{error:?}"))?;
    Ok(report
        .assigned
        .iter()
        .map(|assignment| {
            let element = molecule.stereo_element(assignment.element).unwrap();
            let focus = match &element.kind {
                StereoElementKind::Tetrahedral(value) => StereoFocus::Atom(value.center.index()),
                StereoElementKind::DoubleBond(value) => {
                    StereoFocus::bond(value.left.index(), value.right.index())
                }
                StereoElementKind::Axis(value) => {
                    let axis = molecule.bond(value.axis).unwrap();
                    StereoFocus::bond(axis.a().index(), axis.b().index())
                }
            };
            (focus, format!("{:?}", assignment.descriptor))
        })
        .collect())
}
