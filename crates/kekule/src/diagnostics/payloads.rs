//! Error payloads that several areas share, so that each fact has one type
//! and one message.

use std::fmt;

use super::{counted, Capitalized};
use crate::molecule::AtomId;

/// A kind of numbered item that a shared error payload counts or indexes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Countable {
    /// Atoms.
    Atom,
    /// Bonds.
    Bond,
    /// Stereo elements of a molecule.
    StereoElement,
    /// Stereo groups of a molecule.
    StereoGroup,
    /// Molecule definitions of a topology.
    MoleculeDefinition,
    /// Molecule instances of a topology.
    MoleculeInstance,
    /// Chains of a hierarchy.
    Chain,
    /// Residues of a hierarchy.
    Residue,
    /// Atom sites of a hierarchy.
    AtomSite,
    /// Rows of a dense per-atom array or property table.
    Row,
    /// Positions, one per atom.
    Position,
    /// Geometries of a structure set.
    Geometry,
    /// Structures of a collection.
    Structure,
    /// States of a distribution.
    State,
    /// Samples of a distribution.
    Sample,
    /// Weights, one per atom or sample.
    Weight,
    /// Log weights, one per sample of a distribution.
    LogWeight,
    /// Values of an order, one per geometry of a structure set.
    OrderValue,
    /// Resonance contributors.
    ResonanceContributor,
    /// Solute radii, one per atom.
    SoluteRadius,
    /// Ions to place.
    Ion,
    /// Water molecules to place.
    Water,
    /// Interpretation reports, one per structure written.
    Report,
    /// Values of a property column.
    Value,
}

impl Countable {
    const fn singular(self) -> &'static str {
        match self {
            Self::Atom => "atom",
            Self::Bond => "bond",
            Self::StereoElement => "stereo element",
            Self::StereoGroup => "stereo group",
            Self::MoleculeDefinition => "molecule definition",
            Self::MoleculeInstance => "molecule instance",
            Self::Chain => "chain",
            Self::Residue => "residue",
            Self::AtomSite => "atom site",
            Self::Row => "row",
            Self::Position => "position",
            Self::Geometry => "geometry",
            Self::Structure => "structure",
            Self::State => "state",
            Self::Sample => "sample",
            Self::Weight => "weight",
            Self::LogWeight => "log weight",
            Self::OrderValue => "order value",
            Self::ResonanceContributor => "resonance contributor",
            Self::SoluteRadius => "solute radius",
            Self::Ion => "ion",
            Self::Water => "water molecule",
            Self::Report => "interpretation report",
            Self::Value => "value",
        }
    }

    const fn plural(self) -> &'static str {
        match self {
            Self::Atom => "atoms",
            Self::Bond => "bonds",
            Self::StereoElement => "stereo elements",
            Self::StereoGroup => "stereo groups",
            Self::MoleculeDefinition => "molecule definitions",
            Self::MoleculeInstance => "molecule instances",
            Self::Chain => "chains",
            Self::Residue => "residues",
            Self::AtomSite => "atom sites",
            Self::Row => "rows",
            Self::Position => "positions",
            Self::Geometry => "geometries",
            Self::Structure => "structures",
            Self::State => "states",
            Self::Sample => "samples",
            Self::Weight => "weights",
            Self::LogWeight => "log weights",
            Self::OrderValue => "order values",
            Self::ResonanceContributor => "resonance contributors",
            Self::SoluteRadius => "solute radii",
            Self::Ion => "ions",
            Self::Water => "water molecules",
            Self::Report => "interpretation reports",
            Self::Value => "values",
        }
    }

    /// `count` of this item with its verb, such as `is 1 weight` or `are 4
    /// weights`.
    const fn there(self, count: usize) -> There {
        There { item: self, count }
    }
}

struct There {
    item: Countable,
    count: usize,
}

impl fmt::Display for There {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let verb = if self.count == 1 { "is" } else { "are" };
        let items = counted(self.count, self.item.singular(), self.item.plural());
        write!(formatter, "{verb} {items}")
    }
}

/// Writes the singular noun, such as `geometry`.
impl fmt::Display for Countable {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.singular())
    }
}

/// An index names no item of a collection.
///
/// Writes `Geometry 7 does not exist. There are 5 geometries.`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct IndexOutOfRange {
    item: Countable,
    index: usize,
    len: usize,
}

impl IndexOutOfRange {
    /// Reports that `index` names no `item` of a collection of `len` items.
    pub const fn new(item: Countable, index: usize, len: usize) -> Self {
        Self { item, index, len }
    }

    /// What the collection holds.
    pub const fn item(&self) -> Countable {
        self.item
    }

    /// The requested index.
    pub const fn index(&self) -> usize {
        self.index
    }

    /// The number of items in the collection.
    pub const fn item_count(&self) -> usize {
        self.len
    }
}

impl fmt::Display for IndexOutOfRange {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "{} {} does not exist. There {}.",
            Capitalized(self.item),
            self.index,
            self.item.there(self.len)
        )
    }
}

impl std::error::Error for IndexOutOfRange {}

/// A number of items differs from the number required.
///
/// Writes `There are 4 weights instead of 5.`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct LengthMismatch {
    item: Countable,
    expected: usize,
    actual: usize,
}

impl LengthMismatch {
    /// Reports `actual` items of `item` where `expected` are required.
    pub const fn new(item: Countable, expected: usize, actual: usize) -> Self {
        Self {
            item,
            expected,
            actual,
        }
    }

    /// What was counted.
    pub const fn item(&self) -> Countable {
        self.item
    }

    /// The required number.
    pub const fn expected(&self) -> usize {
        self.expected
    }

    /// The supplied number.
    pub const fn actual(&self) -> usize {
        self.actual
    }
}

impl fmt::Display for LengthMismatch {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "There {} instead of {}.",
            self.item.there(self.actual),
            self.expected
        )
    }
}

impl std::error::Error for LengthMismatch {}

/// More items than their fixed-width identifiers or addressable storage can
/// hold.
///
/// Writes `There are too many atoms.` Unlike a [`LimitExceeded`], no option
/// raises this bound.
///
/// [`LimitExceeded`]: super::LimitExceeded
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct CapacityExceeded {
    item: Countable,
}

impl CapacityExceeded {
    /// Reports more of `item` than can be held.
    pub const fn new(item: Countable) -> Self {
        Self { item }
    }

    /// What there are too many of.
    pub const fn item(&self) -> Countable {
        self.item
    }
}

impl fmt::Display for CapacityExceeded {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "There are too many {}.", self.item.plural())
    }
}

impl std::error::Error for CapacityExceeded {}

/// An atom's implicit hydrogen count is inferred but not perceived, so an
/// operation that needs every hydrogen cannot run.
///
/// Missing valence perception is not zero hydrogens. Perceive valence, for
/// example with `Molecule::perceive`, before the operation.
///
/// Writes `The hydrogen count of atom 3 is unknown. Perceive valence first.`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct UnknownHydrogenCount {
    atom: AtomId,
}

impl UnknownHydrogenCount {
    /// Reports that the hydrogen count of `atom` is unknown.
    pub const fn new(atom: AtomId) -> Self {
        Self { atom }
    }

    /// The atom whose implicit hydrogen count is unknown.
    pub const fn atom(&self) -> AtomId {
        self.atom
    }
}

impl fmt::Display for UnknownHydrogenCount {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "The hydrogen count of {} is unknown. Perceive valence first.",
            self.atom
        )
    }
}

impl std::error::Error for UnknownHydrogenCount {}

/// What has a topology layout other than the one an operation requires.
///
/// Layouts are compared with `Topology::shares_layout`: a topology built
/// independently with equal contents does not share the layout. Each error
/// variant carrying this payload documents what the subject is compared
/// against.
///
/// Writes `The atom selection has a different topology layout.`
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum TopologyMismatch {
    /// An atom selection.
    Selection,
    /// The structure, view, or frame an operation reads.
    Structure,
    /// The moving structure of an alignment.
    Moving,
    /// The reference structure of an alignment or of imaging.
    Reference,
    /// The frame with this index.
    Frame(usize),
    /// The observation with this index.
    Observation(usize),
    /// Solute radii.
    SoluteRadii,
}

impl fmt::Display for TopologyMismatch {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Selection => formatter.write_str("The atom selection has"),
            Self::Structure => formatter.write_str("The structure has"),
            Self::Moving => formatter.write_str("The moving structure has"),
            Self::Reference => formatter.write_str("The reference structure has"),
            Self::Frame(frame) => write!(formatter, "Frame {frame} has"),
            Self::Observation(index) => write!(formatter, "Observation {index} has"),
            Self::SoluteRadii => formatter.write_str("The solute radii have"),
        }?;
        formatter.write_str(" a different topology layout.")
    }
}

impl std::error::Error for TopologyMismatch {}
