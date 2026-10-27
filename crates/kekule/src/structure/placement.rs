//! Builder geometry: coordinates for the atoms an editor adds.
//!
//! Kekule never invents geometry implicitly; placement is an operation a host
//! asks for, on a [`ModelEditor`] draft, as part of the edit it publishes:
//!
//! - [`ModelEditor::add_hydrogens`] makes the implicit hydrogens of chosen
//!   atoms explicit at ideal positions, and [`ModelEditor::remove_hydrogens`]
//!   collapses their explicit hydrogens back into implicit counts, each
//!   reporting per atom what changed ([`PlacedHydrogensReport`],
//!   [`RemovedHydrogensReport`]).
//! - [`ModelEditor::atom_geometry`] perceives an atom's [`AtomGeometry`] and
//!   [`ModelEditor::free_directions`] gives where its next substituents go.
//! - [`ModelEditor::attach_atom`] bonds a new atom to an existing one, in place
//!   of an explicit hydrogen or at the free position nearest a direction
//!   ([`AttachmentSite`]); [`ModelEditor::place_atom`] starts a new molecule
//!   at a point.
//! - [`ModelEditor::attach_fragment`] attaches a named [`Fragment`] of the
//!   library ([`fragments`], [`fragment`]) the same ways, turned to the
//!   clearest of a few staggered torsions about the new bond, and
//!   [`ModelEditor::place_fragment`] places one as a new capped molecule at a
//!   point.
//! - [`ideal_bond_length`] gives the length placement uses for an element pair
//!   and bond order.
//!
//! # Ideal geometry
//!
//! An atom's shape is its number of electron domains: its neighbors, its
//! implicit hydrogens, and the lone pairs that its element's valence
//! electrons, formal charge, radical electrons, and bond orders leave
//! ([`AtomShape`]). A neutral three-coordinate nitrogen with single bonds to a
//! boron, carbon, nitrogen, or oxygen that has a multiple bond (an amide,
//! aniline, or pyrrole nitrogen) is conjugated and planar. Elements outside
//! the s and p blocks have no lone pairs here, so their shape follows their
//! neighbors. Up to six domains are supported.
//!
//! New sites keep the ideal angles of the shape toward the atom's existing
//! neighbors. With one neighbor, the first new site is anti to that
//! neighbor's main substituent, a heavy atom preferring the one with the most
//! heavy neighbors, so chains grow anti and substituents stay staggered;
//! trigonal atoms stay in the plane of that substituent, as across a double
//! bond. Lone pairs of a trigonal bipyramid are equatorial; those of an
//! octahedron are trans to each other or to a neighbor. Placement never moves
//! an existing atom; relaxing the result is a separate step.
//!
//! # Names
//!
//! An atom a placement adds joins the residue of the atom it bonds to, named
//! by its element's uppercase symbol and the next free number in that
//! residue, as chemical components name theirs (`C21`, `H31`): one more than
//! the largest number following that symbol in a name of the residue. An
//! atom without a residue gives its new neighbors none. A new molecule gets a
//! chain and a residue `UNL` of its own. [`ModelEditor::next_atom_name`] and
//! [`ModelEditor::rename_atom`] apply the same rule to existing atoms. Names
//! are label and author atom names alike, and the host may rename.
//!
//! ```
//! use kekule::{smiles, structure::{Model, Positions}};
//! use kekule::hydrogens::AddHydrogensOptions;
//! use kekule::units::{Quantity, NANOMETER};
//! use kekule::geometry::Point3;
//!
//! // Ethanol's heavy atoms, then its hydrogens placed in the same draft.
//! let topology = smiles::to_topology("CCO")?;
//! let positions = Positions::new(Quantity::new(
//!     [Point3::origin(), Point3::new(0.153, 0.0, 0.0), Point3::new(0.2, 0.135, 0.0)],
//!     NANOMETER,
//! ))?;
//! let mut editor = Model::new(topology, positions)?.into_editor();
//! let atoms = editor.topology().atom_ids().collect::<Vec<_>>();
//! let report = editor.add_hydrogens(atoms, AddHydrogensOptions::default())?;
//! assert_eq!(report.added_count(), 6);
//! assert_eq!(editor.finish()?.atom_count(), 9);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
mod fragments;
mod lengths;
mod names;
mod sites;

pub use fragments::{fragment, fragments, Fragment, FragmentKind};

pub(crate) use names::{name_prefix, next_free_number, NameAllocator};
pub(crate) use sites::{free_sites, unit, Rotation, SiteQuery};

use super::ModelEditError;
use crate::algorithms::{AddedHydrogenOrigin, RetainedHydrogenReason};
use crate::diagnostics::LimitExceeded;
use crate::molecule::{BondOrder, Element, ImplicitHydrogens};
use crate::topology::edit::{EditAtomId, EditBondId};
use crate::units::{Quantity, ANGSTROM};
use std::fmt;

#[cfg(doc)]
use super::ModelEditor;

/// The ideal length of a bond of `order` between atoms of elements `a` and
/// `b`: a typical length for common pairs (C–H 1.09 Å, C=O 1.22 Å, …), and
/// otherwise the sum of the elements' Cordero covalent radii, shortened for
/// multiple bonds by Pauling's relation `0.71 Å · log10(order)`. Zero-order
/// and dative bonds take the single-bond length.
///
/// Returns `None` when an element has no covalent radius.
///
/// ```
/// use kekule::{molecule::{BondOrder, Element}, structure::placement::ideal_bond_length};
/// use kekule::units::ANGSTROM;
/// let carbon = Element::from_symbol("C").unwrap();
/// let length = ideal_bond_length(carbon, carbon, BondOrder::Double).unwrap();
/// assert_eq!(length.into_unit(ANGSTROM)?.into_value(), 1.34);
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
pub fn ideal_bond_length(a: Element, b: Element, order: BondOrder) -> Option<Quantity<f64>> {
    lengths::ideal_length_angstrom(a, b, order).map(|length| Quantity::new(length, ANGSTROM))
}

/// The arrangement of an atom's electron domains.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum AtomShape {
    /// One domain, as of a hydrogen.
    Terminal,
    /// Two domains at 180°.
    Linear,
    /// Three coplanar domains at 120°.
    TrigonalPlanar,
    /// Four domains at 109.47°.
    Tetrahedral,
    /// Five domains: two axial, three equatorial.
    TrigonalBipyramidal,
    /// Six domains at 90°.
    Octahedral,
}

impl AtomShape {
    /// The number of domains the shape arranges.
    pub const fn domains(self) -> usize {
        match self {
            Self::Terminal => 1,
            Self::Linear => 2,
            Self::TrigonalPlanar => 3,
            Self::Tetrahedral => 4,
            Self::TrigonalBipyramidal => 5,
            Self::Octahedral => 6,
        }
    }

    pub(crate) const fn from_domains(domains: usize) -> Option<Self> {
        Some(match domains {
            1 => Self::Terminal,
            2 => Self::Linear,
            3 => Self::TrigonalPlanar,
            4 => Self::Tetrahedral,
            5 => Self::TrigonalBipyramidal,
            6 => Self::Octahedral,
            _ => return None,
        })
    }
}

/// An atom's perceived ideal geometry; see the [module](self) for the rules.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomGeometry {
    pub(crate) shape: Option<AtomShape>,
    pub(crate) neighbors: usize,
    pub(crate) implicit_hydrogens: usize,
    pub(crate) lone_pairs: usize,
}

impl AtomGeometry {
    /// The arrangement of the atom's domains, or `None` for an atom with none
    /// (an isolated atom without hydrogens or lone pairs).
    pub const fn shape(&self) -> Option<AtomShape> {
        self.shape
    }
    /// Graph neighbors, explicit hydrogens included.
    pub const fn neighbor_count(&self) -> usize {
        self.neighbors
    }
    /// Implicit hydrogens, which take free sites when made explicit.
    pub const fn implicit_hydrogens(&self) -> usize {
        self.implicit_hydrogens
    }
    /// Lone pairs that take sites of the shape; a conjugated nitrogen's lone
    /// pair takes none.
    pub const fn lone_pairs(&self) -> usize {
        self.lone_pairs
    }
    /// Sites left for new substituents: the shape's domains less neighbors
    /// and lone pairs. Implicit hydrogens count as free.
    pub fn free_site_count(&self) -> usize {
        self.shape.map_or(0, |shape| {
            shape
                .domains()
                .saturating_sub(self.neighbors + self.lone_pairs)
        })
    }
}

/// Where [`ModelEditor::attach_atom`] puts the new atom.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum AttachmentSite {
    /// In place of this explicit hydrogen of the parent, which is deleted;
    /// the new atom lies along the parent–hydrogen direction.
    Hydrogen(EditAtomId),
    /// At the parent's free position nearest this direction, such as a drag
    /// in a viewport; the direction need not be normalized.
    Toward(crate::geometry::Vector3),
}

/// The atom [`ModelEditor::attach_atom`] added.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AttachedAtom {
    pub(crate) atom: EditAtomId,
    pub(crate) bond: EditBondId,
    pub(crate) replaced: Option<EditAtomId>,
}

impl AttachedAtom {
    /// The new atom.
    pub const fn atom(&self) -> EditAtomId {
        self.atom
    }
    /// The new bond to the parent.
    pub const fn bond(&self) -> EditBondId {
        self.bond
    }
    /// The explicit hydrogen deleted for it, under
    /// [`AttachmentSite::Hydrogen`].
    pub const fn replaced_hydrogen(&self) -> Option<EditAtomId> {
        self.replaced
    }
}

/// The atoms [`ModelEditor::attach_fragment`] added.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AttachedFragment {
    pub(crate) atoms: Vec<(crate::molecule::AtomId, EditAtomId)>,
    pub(crate) attachment: EditAtomId,
    pub(crate) bond: EditBondId,
    pub(crate) replaced: Option<EditAtomId>,
}

impl AttachedFragment {
    /// The new atoms with the fragment molecule's atoms they copy, in
    /// molecule order; the leaving atom and its branch are absent.
    pub fn atoms(&self) -> &[(crate::molecule::AtomId, EditAtomId)] {
        &self.atoms
    }
    /// The new atom bonded to the parent.
    pub const fn attachment(&self) -> EditAtomId {
        self.attachment
    }
    /// The new bond to the parent.
    pub const fn bond(&self) -> EditBondId {
        self.bond
    }
    /// The explicit hydrogen deleted for it, under
    /// [`AttachmentSite::Hydrogen`].
    pub const fn replaced_hydrogen(&self) -> Option<EditAtomId> {
        self.replaced
    }
}

/// The atoms [`ModelEditor::place_fragment`] added.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PlacedFragment {
    pub(crate) atoms: Vec<(crate::molecule::AtomId, EditAtomId)>,
}

impl PlacedFragment {
    /// The new atoms with the fragment molecule's atoms they copy, in
    /// molecule order.
    pub fn atoms(&self) -> &[(crate::molecule::AtomId, EditAtomId)] {
        &self.atoms
    }
}

/// The hydrogens [`ModelEditor::add_hydrogens`] placed on one atom.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PlacedHydrogens {
    pub(crate) atom: EditAtomId,
    pub(crate) hydrogens: Vec<EditAtomId>,
    pub(crate) origin: AddedHydrogenOrigin,
}

impl PlacedHydrogens {
    /// The atom whose hydrogens were placed.
    pub const fn atom(&self) -> EditAtomId {
        self.atom
    }
    /// The new hydrogens, in the order of the atom's free sites.
    pub fn hydrogens(&self) -> &[EditAtomId] {
        &self.hydrogens
    }
    /// Whether they came from a fixed or an inferred count.
    pub const fn origin(&self) -> AddedHydrogenOrigin {
        self.origin
    }
}

/// What [`ModelEditor::add_hydrogens`] placed: one entry per distinct atom
/// asked for, in request order, including atoms that needed none.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PlacedHydrogensReport {
    pub(crate) atoms: Vec<PlacedHydrogens>,
}

impl PlacedHydrogensReport {
    /// Per-atom placements.
    pub fn atoms(&self) -> &[PlacedHydrogens] {
        &self.atoms
    }
    /// The number of hydrogens added.
    pub fn added_count(&self) -> usize {
        self.atoms.iter().map(|atom| atom.hydrogens.len()).sum()
    }
}

/// The explicit hydrogens [`ModelEditor::remove_hydrogens`] collapsed on one
/// atom.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemovedHydrogens {
    pub(crate) atom: EditAtomId,
    pub(crate) removed: Vec<EditAtomId>,
    pub(crate) retained: Vec<(EditAtomId, RetainedHydrogenReason)>,
    pub(crate) hydrogens: ImplicitHydrogens,
    pub(crate) implicit_hydrogens: usize,
}

impl RemovedHydrogens {
    /// The atom whose hydrogens were collapsed.
    pub const fn atom(&self) -> EditAtomId {
        self.atom
    }
    /// The deleted hydrogens, which no longer name draft atoms.
    pub fn removed(&self) -> &[EditAtomId] {
        &self.removed
    }
    /// Explicit hydrogens kept, with why, in handle order.
    pub fn retained(&self) -> &[(EditAtomId, RetainedHydrogenReason)] {
        &self.retained
    }
    /// The atom's resulting declaration: inference is kept only when it
    /// reproduces the complete implicit count, which is otherwise fixed.
    pub const fn declaration(&self) -> ImplicitHydrogens {
        self.hydrogens
    }
    /// The atom's complete implicit hydrogen count afterwards.
    pub const fn implicit_hydrogens(&self) -> usize {
        self.implicit_hydrogens
    }
}

/// What [`ModelEditor::remove_hydrogens`] collapsed: one entry per distinct
/// atom asked for, in request order, except hydrogens asked for that were
/// collapsed into another atom asked for, which its entry lists as removed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemovedHydrogensReport {
    pub(crate) atoms: Vec<RemovedHydrogens>,
}

impl RemovedHydrogensReport {
    /// Per-atom outcomes.
    pub fn atoms(&self) -> &[RemovedHydrogens] {
        &self.atoms
    }
    /// The number of hydrogens deleted.
    pub fn removed_count(&self) -> usize {
        self.atoms.iter().map(|atom| atom.removed.len()).sum()
    }
}

/// A placement was rejected; the draft is unchanged.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PlacementError {
    /// Editing the draft failed, for example on a deleted handle; see the
    /// source.
    Edit(Box<ModelEditError>),
    /// The atom's represented valence exceeds what its element allows, so its
    /// hydrogens cannot be completed.
    InvalidValence {
        /// The atom.
        atom: EditAtomId,
    },
    /// The atom has no free site for a new substituent.
    NoFreePosition {
        /// The atom.
        atom: EditAtomId,
    },
    /// The atom would arrange more than six domains.
    UnsupportedGeometry {
        /// The atom.
        atom: EditAtomId,
        /// Its domain count.
        domains: usize,
    },
    /// The atom's neighbors coincide with it or lie so that no orientation
    /// follows from them.
    DegenerateGeometry {
        /// The atom.
        atom: EditAtomId,
    },
    /// An element of the bond has no covalent radius.
    UnknownBondLength {
        /// One element of the bond.
        first: Element,
        /// The other element.
        second: Element,
    },
    /// The atom named as the hydrogen to replace is not an explicit hydrogen
    /// singly bonded to the parent.
    NotReplaceableHydrogen {
        /// The parent.
        atom: EditAtomId,
        /// The atom named.
        hydrogen: EditAtomId,
    },
    /// The direction is zero or not finite.
    InvalidDirection,
    /// The atom's represented stereo names an implicit hydrogen that the
    /// requested count cannot materialize as exactly one atom, or an axis
    /// names one.
    StereoHydrogen {
        /// The atom.
        atom: EditAtomId,
    },
    /// Collapsing would give the atom more implicit hydrogens than a count
    /// can hold.
    HydrogenCountOverflow {
        /// The atom.
        atom: EditAtomId,
        /// The count.
        count: usize,
    },
    /// More hydrogens would be added than the options allow.
    LimitExceeded(LimitExceeded),
    /// No library fragment has this name.
    UnknownFragment {
        /// The name asked for.
        name: String,
    },
    /// The fragment atom asked to attach through carries no hydrogen.
    FragmentAtomWithoutHydrogen {
        /// The fragment.
        fragment: &'static str,
        /// The atom of its molecule.
        atom: crate::molecule::AtomId,
    },
}

impl fmt::Display for PlacementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Edit(_) => f.write_str("cannot edit the draft"),
            Self::InvalidValence { atom } => write!(f, "{atom:?} exceeds its valence"),
            Self::NoFreePosition { atom } => write!(f, "{atom:?} has no free position"),
            Self::UnsupportedGeometry { atom, domains } => {
                write!(
                    f,
                    "{atom:?} would arrange {domains} domains; at most 6 are supported"
                )
            }
            Self::DegenerateGeometry { atom } => {
                write!(f, "the neighbors of {atom:?} give it no orientation")
            }
            Self::UnknownBondLength { first, second } => write!(
                f,
                "no bond length between {} and {}",
                first.symbol(),
                second.symbol()
            ),
            Self::NotReplaceableHydrogen { atom, hydrogen } => write!(
                f,
                "{hydrogen:?} is not an explicit hydrogen singly bonded to {atom:?}"
            ),
            Self::InvalidDirection => f.write_str("direction is zero or not finite"),
            Self::StereoHydrogen { atom } => write!(
                f,
                "the represented stereo of {atom:?} cannot follow its hydrogens"
            ),
            Self::HydrogenCountOverflow { atom, count } => {
                write!(f, "{atom:?} would have {count} implicit hydrogens")
            }
            Self::LimitExceeded(_) => f.write_str("too many hydrogens to add"),
            Self::UnknownFragment { name } => write!(f, "no fragment named `{name}`"),
            Self::FragmentAtomWithoutHydrogen { fragment, atom } => {
                write!(f, "{atom:?} of fragment {fragment} carries no hydrogen")
            }
        }
    }
}

impl std::error::Error for PlacementError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Edit(error) => Some(error.as_ref()),
            Self::LimitExceeded(error) => Some(error),
            _ => None,
        }
    }
}

impl From<ModelEditError> for PlacementError {
    fn from(error: ModelEditError) -> Self {
        Self::Edit(Box::new(error))
    }
}

impl From<crate::topology::edit::TopologyEditError> for PlacementError {
    fn from(error: crate::topology::edit::TopologyEditError) -> Self {
        Self::Edit(Box::new(error.into()))
    }
}

/// The electron-domain count of an atom from its chemistry.
#[derive(Debug, Clone, Copy)]
pub(crate) struct DomainInput {
    pub(crate) element: Element,
    pub(crate) formal_charge: i8,
    pub(crate) radical_electrons: usize,
    /// Neighbors (any bond order) plus implicit hydrogens.
    pub(crate) substituents: usize,
    /// Bond orders plus implicit hydrogens.
    pub(crate) bond_valence: usize,
    /// Zero-order and dative bonds, which may take a lone pair.
    pub(crate) coordinate_bonds: usize,
    /// A three-coordinate nitrogen bonded to an atom with a multiple bond.
    pub(crate) conjugated_nitrogen: bool,
}

/// The domain count, which may be zero or exceed six, and the lone pairs
/// that take sites among them.
pub(crate) fn perceive_domains(input: DomainInput) -> (usize, usize) {
    let lone_pairs = valence_electrons(input.element.atomic_number()).map_or(0, |valence| {
        let free = i64::from(valence)
            - i64::from(input.formal_charge)
            - i64::try_from(input.bond_valence).unwrap_or(i64::MAX / 4)
            - i64::try_from(input.radical_electrons).unwrap_or(i64::MAX / 4);
        let pairs = usize::try_from(free.max(0) / 2).expect("non-negative");
        // A coordinate bond donates one of the atom's own pairs.
        pairs.saturating_sub(input.coordinate_bonds)
    });
    if input.conjugated_nitrogen && lone_pairs == 1 && input.substituents == 3 {
        return (3, 0);
    }
    (input.substituents + lone_pairs, lone_pairs)
}

/// Valence electrons of s- and p-block elements; `None` elsewhere.
const fn valence_electrons(atomic_number: u8) -> Option<u8> {
    Some(match atomic_number {
        1 | 2 => atomic_number,
        3..=10 => atomic_number - 2,
        11..=18 => atomic_number - 10,
        19 | 37 | 55 | 87 => 1,
        20 | 38 | 56 | 88 => 2,
        31..=36 => atomic_number - 28,
        49..=54 => atomic_number - 46,
        81..=86 => atomic_number - 78,
        113..=118 => atomic_number - 110,
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn domains(symbol: &str, charge: i8, substituents: usize, valence: usize) -> (usize, usize) {
        perceive_domains(DomainInput {
            element: Element::from_symbol(symbol).unwrap(),
            formal_charge: charge,
            radical_electrons: 0,
            substituents,
            bond_valence: valence,
            coordinate_bonds: 0,
            conjugated_nitrogen: false,
        })
    }

    #[test]
    fn domain_counts_follow_valence_shell_electron_pairs() {
        assert_eq!(domains("C", 0, 4, 4), (4, 0)); // methane
        assert_eq!(domains("C", 0, 3, 4), (3, 0)); // alkene carbon
        assert_eq!(domains("C", 0, 2, 4), (2, 0)); // alkyne carbon
        assert_eq!(domains("N", 0, 1, 3), (2, 1)); // nitrile nitrogen
        assert_eq!(domains("O", 0, 2, 2), (4, 2)); // water
        assert_eq!(domains("O", 0, 1, 2), (3, 2)); // carbonyl oxygen
        assert_eq!(domains("N", 1, 4, 4), (4, 0)); // ammonium
        assert_eq!(domains("N", 1, 3, 4), (3, 0)); // nitro nitrogen
        assert_eq!(domains("C", 1, 3, 3), (3, 0)); // carbocation
        assert_eq!(domains("C", -1, 3, 3), (4, 1)); // carbanion
        assert_eq!(domains("S", 0, 4, 6), (4, 0)); // sulfonyl
        assert_eq!(domains("S", 0, 3, 4), (4, 1)); // sulfoxide
        assert_eq!(domains("P", 0, 5, 5), (5, 0)); // phosphorus pentachloride
        assert_eq!(domains("S", 0, 6, 6), (6, 0)); // sulfur hexafluoride
        assert_eq!(domains("S", 0, 4, 4), (5, 1)); // sulfur tetrafluoride
        assert_eq!(domains("Xe", 0, 4, 4), (6, 2)); // xenon tetrafluoride
        assert_eq!(domains("Xe", 0, 2, 2), (5, 3)); // xenon difluoride
        assert_eq!(domains("Fe", 0, 6, 6), (6, 0)); // no lone pairs off the main group
    }
}
