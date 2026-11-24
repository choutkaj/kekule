//! What placement edits report: attachment sites, placed and removed atoms
//! and hydrogens, and the placement error type.

use super::BuilderValenceError;
use crate::algorithms::AddedHydrogenOrigin;
use crate::diagnostics::LimitExceeded;
use crate::molecule::Element;
use crate::structure::StructureEditError;
use crate::topology::edit::{EditAtomId, EditBondId};
use crate::topology::editor::RemovedHydrogens;
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
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AttachedAtom {
    pub(crate) atom: EditAtomId,
    pub(crate) bond: EditBondId,
    pub(crate) replaced: Vec<EditAtomId>,
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
    /// The explicit hydrogens of the parent deleted for it, which no longer
    /// name draft atoms: under [`AttachmentSite::Hydrogen`] that hydrogen
    /// first, then those a multiple bond takes.
    pub fn replaced_hydrogens(&self) -> &[EditAtomId] {
        &self.replaced
    }
}

/// The atoms [`StructureEditor::attach_fragment`] added.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AttachedFragment {
    pub(crate) atoms: Vec<(crate::molecule::AtomId, EditAtomId)>,
    pub(crate) attachment: EditAtomId,
    pub(crate) bond: EditBondId,
    pub(crate) replaced: Vec<EditAtomId>,
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
    /// The explicit hydrogens of the parent deleted for it, which no longer
    /// name draft atoms: under [`AttachmentSite::Hydrogen`] that hydrogen
    /// first, then those a multiple bond takes.
    pub fn replaced_hydrogens(&self) -> &[EditAtomId] {
        &self.replaced
    }
}

/// The atoms [`StructureEditor::fuse_fragment`] added and removed.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FusedFragment {
    pub(crate) atoms: Vec<(crate::molecule::AtomId, EditAtomId)>,
    pub(crate) fused: [(crate::molecule::AtomId, EditAtomId); 2],
    pub(crate) bond: EditBondId,
    pub(crate) removed: Vec<EditAtomId>,
}

impl FusedFragment {
    /// The new atoms with the fragment molecule's atoms they copy, in
    /// molecule order; the fragment's two fused atoms and their hydrogens
    /// are absent.
    pub fn atoms(&self) -> &[(crate::molecule::AtomId, EditAtomId)] {
        &self.atoms
    }
    /// The fragment's fusion-bond atoms with the target atoms they merged
    /// into.
    pub const fn fused_atoms(&self) -> [(crate::molecule::AtomId, EditAtomId); 2] {
        self.fused
    }
    /// The target bond, now shared by the two rings.
    pub const fn bond(&self) -> EditBondId {
        self.bond
    }
    /// The explicit hydrogens of the target atoms that were deleted, which
    /// no longer name draft atoms.
    pub fn removed_hydrogens(&self) -> &[EditAtomId] {
        &self.removed
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
        self.atoms.iter().map(|atom| atom.removed().len()).sum()
    }
}

/// What [`StructureEditor::settle_atoms`] did to one atom.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SettledAtom {
    pub(crate) atom: EditAtomId,
    pub(crate) added: Vec<EditAtomId>,
    pub(crate) removed: Vec<EditAtomId>,
    pub(crate) moved: Vec<EditAtomId>,
    pub(crate) turned: Vec<EditAtomId>,
    pub(crate) angle_in_ring: bool,
}

impl SettledAtom {
    /// The settled atom.
    pub const fn atom(&self) -> EditAtomId {
        self.atom
    }
    /// Its new hydrogens.
    pub fn added_hydrogens(&self) -> &[EditAtomId] {
        &self.added
    }
    /// Its deleted hydrogens, which no longer name draft atoms.
    pub fn removed_hydrogens(&self) -> &[EditAtomId] {
        &self.removed
    }
    /// Its kept hydrogens that moved to ideal positions.
    pub fn moved_hydrogens(&self) -> &[EditAtomId] {
        &self.moved
    }
    /// The atoms of the side that turned about the atom to give its two
    /// neighbors the angle of its shape; empty when none turned.
    pub fn turned_atoms(&self) -> &[EditAtomId] {
        &self.turned
    }
    /// Whether the angle between the atom's two neighbors is off its shape's
    /// but both lie in a ring with the atom, so no side turned; cleaning up
    /// the geometry fixes it.
    pub const fn angle_in_ring(&self) -> bool {
        self.angle_in_ring
    }
}

/// What [`StructureEditor::settle_atoms`] did: one entry per distinct atom
/// asked for, in request order, except hydrogens that belong to another
/// atom.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SettledAtomsReport {
    pub(crate) atoms: Vec<SettledAtom>,
}

impl SettledAtomsReport {
    /// Per-atom outcomes.
    pub fn atoms(&self) -> &[SettledAtom] {
        &self.atoms
    }
    /// The number of hydrogens added.
    pub fn added_count(&self) -> usize {
        self.atoms.iter().map(|atom| atom.added.len()).sum()
    }
    /// The number of hydrogens deleted.
    pub fn removed_count(&self) -> usize {
        self.atoms.iter().map(|atom| atom.removed.len()).sum()
    }
}

/// What [`StructureEditor::settle_bonds`] did to one bond.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SettledBond {
    pub(crate) bond: EditBondId,
    pub(crate) moved: Vec<EditAtomId>,
    pub(crate) in_ring: bool,
}

impl SettledBond {
    /// The settled bond.
    pub const fn bond(&self) -> EditBondId {
        self.bond
    }
    /// The atoms of the side that moved along the bond to give it its ideal
    /// length; empty when none moved.
    pub fn moved_atoms(&self) -> &[EditAtomId] {
        &self.moved
    }
    /// Whether the bond's length is off its ideal one but the bond lies in a
    /// ring, so no side moved; cleaning up the geometry fixes it.
    pub const fn in_ring(&self) -> bool {
        self.in_ring
    }
}

/// What [`StructureEditor::settle_bonds`] did: one entry per distinct bond
/// asked for, in request order.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SettledBondsReport {
    pub(crate) bonds: Vec<SettledBond>,
}

impl SettledBondsReport {
    /// Per-bond outcomes.
    pub fn bonds(&self) -> &[SettledBond] {
        &self.bonds
    }
}

/// A placement was rejected; the draft is unchanged.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PlacementError {
    /// Editing the draft failed, for example on a deleted handle; see the
    /// source.
    Edit(StructureEditError),
    /// The atom's element, charge, or valence is outside the builder's
    /// states, so [`StructureEditor::settle_atoms`] cannot give it hydrogens;
    /// see [`builder_hydrogen_count`](super::builder_hydrogen_count).
    BuilderValence {
        /// The atom.
        atom: EditAtomId,
        /// What is outside the builder's states.
        error: BuilderValenceError,
    },
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
    /// The fragment asked to fuse is a group; only rings fuse onto a bond.
    NotRingFragment {
        /// The fragment.
        fragment: &'static str,
    },
    /// An atom of the target bond is a hydrogen.
    FusionAtHydrogen {
        /// The bond.
        bond: EditBondId,
    },
    /// The target bond is neither single nor double.
    FusionBondOrder {
        /// The bond.
        bond: EditBondId,
        /// Its order.
        order: crate::molecule::BondOrder,
    },
    /// The target bond already belongs to two rings of up to eight atoms,
    /// as the shared bond of naphthalene does.
    SharedFusionBond {
        /// The bond.
        bond: EditBondId,
    },
    /// An atom of the target bond carries no hydrogen, explicit or
    /// implicit, for the fused ring's bond to replace.
    FusionAtomWithoutHydrogen {
        /// The bond.
        bond: EditBondId,
        /// The atom.
        atom: EditAtomId,
    },
    /// The target bond is a double bond whose represented stereo puts the
    /// fused ring's two bonds trans, which a small ring cannot be.
    TransFusionBond {
        /// The bond.
        bond: EditBondId,
    },
}

impl fmt::Display for PlacementError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use crate::diagnostics::Capitalized;
        match self {
            Self::Edit(error) => error.fmt(f),
            Self::BuilderValence { atom, error } => error.write_about(f, Capitalized(atom)),
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
            Self::NotRingFragment { fragment } => {
                write!(f, "Fragment {fragment:?} is a group, not a ring.")
            }
            Self::FusionAtHydrogen { bond } => {
                write!(f, "Cannot fuse a ring onto {bond}. It ends at a hydrogen.")
            }
            Self::FusionBondOrder { bond, order } => {
                use crate::molecule::BondOrder;
                let order = match order {
                    BondOrder::Zero => "a zero-order",
                    BondOrder::Single => "a single",
                    BondOrder::Double => "a double",
                    BondOrder::Triple => "a triple",
                    BondOrder::Quadruple => "a quadruple",
                    BondOrder::Dative => "a dative",
                };
                write!(f, "Cannot fuse a ring onto {bond}. It is {order} bond.")
            }
            Self::SharedFusionBond { bond } => write!(
                f,
                "Cannot fuse a ring onto {bond}. It is already shared by two rings."
            ),
            Self::TransFusionBond { bond } => write!(
                f,
                "Cannot fuse a ring onto {bond}. Its stereo puts the ring bonds trans."
            ),
            Self::FusionAtomWithoutHydrogen { bond, atom } => write!(
                f,
                "Cannot fuse a ring onto {bond}. {} carries no hydrogen.",
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
        Self::Edit(error)
    }
}

impl From<crate::topology::edit::TopologyEditError> for PlacementError {
    fn from(error: crate::topology::edit::TopologyEditError) -> Self {
        Self::Edit(error.into())
    }
}
