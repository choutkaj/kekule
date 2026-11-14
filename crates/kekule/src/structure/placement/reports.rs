//! What placement edits report: attachment sites, placed and removed atoms
//! and hydrogens, and the placement error type.

use crate::algorithms::{AddedHydrogenOrigin, RetainedHydrogenReason};
use crate::diagnostics::LimitExceeded;
use crate::molecule::{Element, ImplicitHydrogens};
use crate::structure::StructureEditError;
use crate::topology::edit::{EditAtomId, EditBondId};
use std::fmt;

#[cfg(doc)]
use crate::structure::StructureEditor;

/// Where [`StructureEditor::attach_atom`] puts the new atom.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum AttachmentSite {
    /// In place of this explicit hydrogen of the parent, which is deleted;
    /// the new atom lies along the parent–hydrogen direction.
    Hydrogen(EditAtomId),
    /// At the parent's free position nearest this direction, such as a drag
    /// in a viewport; the direction need not be normalized.
    Toward(crate::math::Vector3),
}

/// The atom [`StructureEditor::attach_atom`] added.
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

/// The atoms [`StructureEditor::attach_fragment`] added.
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

/// The atoms [`StructureEditor::place_fragment`] added.
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

/// The hydrogens [`StructureEditor::add_hydrogens`] placed on one atom.
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

/// What [`StructureEditor::add_hydrogens`] placed: one entry per distinct atom
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

/// The explicit hydrogens [`StructureEditor::remove_hydrogens`] collapsed on one
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

/// What [`StructureEditor::remove_hydrogens`] collapsed: one entry per distinct
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
    Edit(Box<StructureEditError>),
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
        use crate::diagnostics::Capitalized;
        match self {
            Self::Edit(error) => error.fmt(f),
            Self::InvalidValence { atom } => {
                write!(f, "{} exceeds its valence.", Capitalized(atom))
            }
            Self::NoFreePosition { atom } => {
                write!(f, "{} has no free position.", Capitalized(atom))
            }
            Self::UnsupportedGeometry { atom, domains } => write!(
                f,
                "{} would arrange {domains} electron domains. At most 6 are supported.",
                Capitalized(atom)
            ),
            Self::DegenerateGeometry { atom } => {
                write!(f, "The neighbors of {atom} give it no orientation.")
            }
            Self::UnknownBondLength { first, second } => write!(
                f,
                "No bond length is known between {} and {}.",
                first.symbol(),
                second.symbol()
            ),
            Self::NotReplaceableHydrogen { atom, hydrogen } => write!(
                f,
                "{} is not an explicit hydrogen singly bonded to {atom}.",
                Capitalized(hydrogen)
            ),
            Self::InvalidDirection => f.write_str("The direction is zero or not finite."),
            Self::StereoHydrogen { atom } => write!(
                f,
                "Cannot place the hydrogens of {atom}. Its stereo names an implicit hydrogen."
            ),
            Self::HydrogenCountOverflow { atom, count } => write!(
                f,
                "{} would have {count} implicit hydrogens, more than a count can hold.",
                Capitalized(atom)
            ),
            Self::LimitExceeded(limit) => {
                write!(f, "Hydrogen placement exceeded its {}.", limit.describe())
            }
            Self::UnknownFragment { name } => write!(f, "No fragment is named {name:?}."),
            Self::FragmentAtomWithoutHydrogen { fragment, atom } => write!(
                f,
                "{} of fragment {fragment:?} carries no hydrogen.",
                Capitalized(atom)
            ),
        }
    }
}

impl std::error::Error for PlacementError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Edit(error) => error.source(),
            _ => None,
        }
    }
}

impl From<StructureEditError> for PlacementError {
    fn from(error: StructureEditError) -> Self {
        Self::Edit(Box::new(error))
    }
}

impl From<crate::topology::edit::TopologyEditError> for PlacementError {
    fn from(error: crate::topology::edit::TopologyEditError) -> Self {
        Self::Edit(Box::new(error.into()))
    }
}
