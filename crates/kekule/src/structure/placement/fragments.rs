//! Named fragments with ideal 3D coordinates and an attachment atom.
use super::PlacementError;
use crate::geometry::Point3;
use crate::molecule::{AtomId, BondOrder, Molecule};
use crate::units::{Quantity, CANONICAL_LENGTH_UNIT};
use std::collections::BTreeSet;
use std::sync::OnceLock;

/// The library: PubChem 3D conformers of the fragments' capped parent
/// compounds, written by `tools/fragments/build_fragments.py`; each record
/// names its fragment and gives its kind, attachment, and leaving atom.
const LIBRARY: &str = include_str!("fragments.sdf");

/// What a [`Fragment`] is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum FragmentKind {
    /// A ring or ring system.
    Ring,
    /// A substituent group.
    Group,
}

/// A named fragment: a capped molecule with ideal 3D coordinates, the atom
/// that bonds when it is attached, and the leaving atom that caps that bond.
///
/// The molecule is the fragment's parent compound with explicit hydrogens,
/// as PubChem's computed 3D conformer gives it (phenyl's is benzene, nitro's
/// nitromethane). Attaching it removes the leaving atom and, when that is a
/// heavy atom, the branch beyond it; placing it alone keeps the whole capped
/// molecule. See [`fragments`] for the library.
///
/// Atoms are in a meaningful order, which placement names follow. A ring's
/// atoms come first in its conventional numbering (IUPAC locants: pyridine
/// N1, C2, …, C6; naphthalene C1, …, C4, C4a, C5, …, C8, C8a), so a placed
/// benzene's carbons are C1 to C6 around the ring and C1 and C4 are para. A
/// group's attachment atom comes first, then its other heavy atoms. Each
/// hydrogen follows, in the order of the atoms that carry them.
#[derive(Debug, Clone, PartialEq)]
pub struct Fragment {
    name: &'static str,
    kind: FragmentKind,
    pubchem_cid: u64,
    molecule: Molecule,
    /// Canonical positions in molecule atom order.
    positions: Vec<Point3>,
    attachment: AtomId,
    leaving: AtomId,
}

impl Fragment {
    /// The fragment's name, as a palette shows it (`"pyridine"`,
    /// `"tert-butyl"`).
    pub const fn name(&self) -> &'static str {
        self.name
    }
    /// Whether it is a ring or a group.
    pub const fn kind(&self) -> FragmentKind {
        self.kind
    }
    /// The PubChem compound the coordinates come from.
    pub const fn pubchem_cid(&self) -> u64 {
        self.pubchem_cid
    }
    /// The capped molecule, with explicit hydrogens.
    pub const fn molecule(&self) -> &Molecule {
        &self.molecule
    }
    /// The ideal position of one atom of [`Self::molecule`], or `None` when
    /// it names none.
    pub fn position(&self, atom: AtomId) -> Option<Quantity<Point3>> {
        let index = self.molecule.atom_ids().position(|id| id == atom)?;
        Some(Quantity::new(self.positions[index], CANONICAL_LENGTH_UNIT))
    }
    /// The atom that bonds to the atom the fragment is attached to.
    pub const fn attachment_atom(&self) -> AtomId {
        self.attachment
    }
    /// The atom bonded to the attachment atom that attaching removes, with
    /// any branch beyond it; usually a hydrogen.
    pub const fn leaving_atom(&self) -> AtomId {
        self.leaving
    }
    /// The bond order between the attachment and leaving atoms, which the
    /// new bond takes.
    pub fn attachment_order(&self) -> BondOrder {
        let bond = self
            .molecule
            .bond_between(self.attachment, self.leaving)
            .expect("the leaving atom is bonded to the attachment atom");
        self.molecule.bond(bond).expect("live bond").order
    }
    /// The same fragment attached through `atom` instead, which must carry an
    /// explicit hydrogen; that hydrogen leaves (for example pyridine through
    /// a carbon next to the nitrogen).
    ///
    /// # Errors
    ///
    /// [`PlacementError::FragmentAtomWithoutHydrogen`] when `atom` is not an atom
    /// of the fragment with a hydrogen.
    pub fn with_attachment(&self, atom: AtomId) -> Result<Self, PlacementError> {
        let leaving = self
            .molecule
            .neighbors(atom)
            .into_iter()
            .flatten()
            .filter(|neighbor| {
                self.molecule
                    .atom(*neighbor)
                    .is_some_and(|value| value.element.atomic_number() == 1)
            })
            .min()
            .ok_or(PlacementError::FragmentAtomWithoutHydrogen {
                fragment: self.name,
                atom,
            })?;
        Ok(Self {
            attachment: atom,
            leaving,
            ..self.clone()
        })
    }
    /// The atoms an attachment keeps: all but the leaving atom and the branch
    /// beyond it, in molecule order, with their positions.
    pub(crate) fn kept_atoms(&self) -> Vec<(AtomId, Point3)> {
        let mut removed = BTreeSet::from([self.leaving]);
        let mut pending = vec![self.leaving];
        while let Some(atom) = pending.pop() {
            for neighbor in self.molecule.neighbors(atom).into_iter().flatten() {
                if neighbor != self.attachment && removed.insert(neighbor) {
                    pending.push(neighbor);
                }
            }
        }
        self.molecule
            .atom_ids()
            .zip(self.positions.iter().copied())
            .filter(|(atom, _)| !removed.contains(atom))
            .collect()
    }
    /// Canonical positions in molecule atom order.
    pub(crate) fn points(&self) -> &[Point3] {
        &self.positions
    }
}

/// The fragment library, rings first and then groups:
///
/// - rings: benzene, pyridine (attached at C4), pyrimidine (C2), pyrazine,
///   pyrrole (N1), furan (C2), thiophene (C2), imidazole (N1), pyrazole
///   (N1), oxazole (C2), thiazole (C2), cyclopropane, cyclobutane,
///   cyclopentane, cyclohexane (chair), piperidine (N1), piperazine (N1),
///   morpholine (N4), tetrahydrofuran (C2), naphthalene (C2), indole (C3),
///   quinoline (C4), benzimidazole (C2);
/// - groups: methyl, ethyl, isopropyl, tert-butyl, vinyl, ethynyl, phenyl,
///   benzyl, hydroxyl, methoxy, amino, methylamino, dimethylamino, formyl,
///   acetyl, carboxyl, methoxycarbonyl, carbamoyl, acetamido (through N),
///   nitro, cyano, trifluoromethyl, methylsulfonyl, sulfonamide
///   (`-SO2NH2`), thiol.
///
/// Coordinates are PubChem's computed 3D conformers of the capped parent
/// compounds (MMFF94s-optimized, from OpenEye OMEGA).
pub fn fragments() -> &'static [Fragment] {
    static LIBRARY_FRAGMENTS: OnceLock<Vec<Fragment>> = OnceLock::new();
    LIBRARY_FRAGMENTS.get_or_init(load)
}

/// The library fragment named `name`, compared without ASCII case.
///
/// # Errors
///
/// [`PlacementError::UnknownFragment`] when no fragment has that name.
pub fn fragment(name: &str) -> Result<&'static Fragment, PlacementError> {
    fragments()
        .iter()
        .find(|fragment| fragment.name.eq_ignore_ascii_case(name))
        .ok_or_else(|| PlacementError::UnknownFragment {
            name: name.to_owned(),
        })
}

fn load() -> Vec<Fragment> {
    let document = crate::io::parse_sdf_document(LIBRARY, Default::default())
        .expect("the bundled fragment library parses");
    document
        .records()
        .iter()
        .map(|record| {
            let field = |name: &str| {
                record
                    .data_fields()
                    .iter()
                    .find(|field| field.name() == name)
                    .map(|field| field.value().trim().to_owned())
                    .expect("bundled fragments carry their fields")
            };
            let model = record
                .interpret(Default::default())
                .expect("bundled fragments interpret")
                .into_model();
            let topology = model.topology();
            let atom_at = |field_value: String| {
                let index: usize = field_value.parse().expect("an atom number");
                topology
                    .atom_at(crate::topology::TopologyAtomIndex::new(
                        u32::try_from(index - 1).expect("small"),
                    ))
                    .expect("a bundled atom")
                    .id()
                    .atom()
            };
            let attachment = atom_at(field("attachment"));
            let leaving = atom_at(field("leaving"));
            let molecule = topology
                .definitions()
                .next()
                .expect("one molecule")
                .molecule()
                .clone();
            // Stored positions are in the canonical length unit.
            let points = model.positions().values().value().to_vec();
            let mut positions = vec![Point3::origin(); points.len()];
            for view in topology.atoms() {
                let order = molecule
                    .atom_ids()
                    .position(|id| id == view.id().atom())
                    .expect("one molecule");
                positions[order] = points[view.index().index()];
            }
            Fragment {
                name: name_of(record.title()),
                kind: match field("kind").as_str() {
                    "ring" => FragmentKind::Ring,
                    _ => FragmentKind::Group,
                },
                pubchem_cid: field("PUBCHEM_COMPOUND_CID").parse().expect("a CID"),
                molecule,
                positions,
                attachment,
                leaving,
            }
        })
        .collect()
}

/// The library's static name for a record title.
fn name_of(title: &str) -> &'static str {
    const NAMES: [&str; 48] = [
        "benzene",
        "pyridine",
        "pyrimidine",
        "pyrazine",
        "pyrrole",
        "furan",
        "thiophene",
        "imidazole",
        "pyrazole",
        "oxazole",
        "thiazole",
        "cyclopropane",
        "cyclobutane",
        "cyclopentane",
        "cyclohexane",
        "piperidine",
        "piperazine",
        "morpholine",
        "tetrahydrofuran",
        "naphthalene",
        "indole",
        "quinoline",
        "benzimidazole",
        "methyl",
        "ethyl",
        "isopropyl",
        "tert-butyl",
        "vinyl",
        "ethynyl",
        "phenyl",
        "benzyl",
        "hydroxyl",
        "methoxy",
        "amino",
        "methylamino",
        "dimethylamino",
        "formyl",
        "acetyl",
        "carboxyl",
        "methoxycarbonyl",
        "carbamoyl",
        "acetamido",
        "nitro",
        "cyano",
        "trifluoromethyl",
        "methylsulfonyl",
        "sulfonamide",
        "thiol",
    ];
    NAMES
        .iter()
        .find(|name| **name == title.trim())
        .expect("bundled fragment names are listed")
}
