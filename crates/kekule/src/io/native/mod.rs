//! Kekule's native binary serialization: exact, versioned persistence of one
//! molecule, topology, structure, structure set, or ensemble, and of a
//! structure or structure set with its structure interpretation reports.
//!
//! The byte layout is specified in `docs/native-format.md`; the encoder and
//! decoder here implement that specification and nothing else.

mod crc;
mod decode;
mod description;
mod encode;
pub(crate) mod entry_points;

use std::fmt;

use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail};
use crate::io::{
    StructureInterpretation, StructureInterpretationReport, StructureSetInterpretation,
};
use crate::molecule::Molecule;
use crate::structure::{Ensemble, Structure, StructureSet};
use crate::topology::Topology;

pub(crate) use decode::read_native;
pub use description::{EnsembleDescription, StructureSetDescription};
pub(crate) use encode::{write_native, write_native_to};

/// The eight bytes every native input starts with.
///
/// The non-ASCII first byte and the CR LF, EOF, and LF bytes detect
/// text-mode transfer damage, as in PNG.
pub const NATIVE_MAGIC: [u8; 8] = [0x89, b'K', b'E', b'K', 0x0D, 0x0A, 0x1A, 0x0A];

/// The native format version this crate writes and reads.
///
/// Before Kekule 1.0 the format may change incompatibly, but every
/// incompatible change increments this version, and readers reject versions
/// they do not implement.
pub const NATIVE_FORMAT_VERSION: u32 = 6;

/// Default [`NativeReadOptions::max_input_bytes`]: 64 GiB.
pub const DEFAULT_NATIVE_MAX_INPUT_BYTES: usize = 64 << 30;
/// Default [`NativeReadOptions::max_atoms`].
pub const DEFAULT_NATIVE_MAX_ATOMS: usize = 100_000_000;
/// Default [`NativeReadOptions::max_bonds`].
pub const DEFAULT_NATIVE_MAX_BONDS: usize = 200_000_000;
/// Default [`NativeReadOptions::max_items`].
pub const DEFAULT_NATIVE_MAX_ITEMS: usize = 100_000_000;
/// Default [`NativeReadOptions::max_string_bytes`]: 1 MiB.
pub const DEFAULT_NATIVE_MAX_STRING_BYTES: usize = 1 << 20;

/// The kind of object a native input holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum NativeKind {
    /// A [`Molecule`].
    Molecule,
    /// A [`Topology`].
    Topology,
    /// A [`Structure`].
    Structure,
    /// A [`Structure`] with the report of its interpretation, a
    /// [`StructureInterpretation`].
    StructureInterpretation,
    /// A [`StructureSet`], with its geometries or with only their number
    /// ([`StructureSetDescription`]).
    StructureSet,
    /// A [`StructureSet`] with one interpretation report per geometry, a
    /// [`StructureSetInterpretation`].
    StructureSetInterpretation,
    /// An [`Ensemble`], with its set's geometries or with only their number
    /// ([`EnsembleDescription`]).
    Ensemble,
}

impl NativeKind {
    pub(super) const fn code(self) -> u8 {
        match self {
            Self::Molecule => 1,
            Self::Topology => 2,
            Self::Structure => 3,
            Self::StructureInterpretation => 7,
            Self::StructureSet => 9,
            Self::StructureSetInterpretation => 10,
            Self::Ensemble => 11,
        }
    }

    pub(super) const fn from_code(code: u8) -> Option<Self> {
        Some(match code {
            1 => Self::Molecule,
            2 => Self::Topology,
            3 => Self::Structure,
            7 => Self::StructureInterpretation,
            9 => Self::StructureSet,
            10 => Self::StructureSetInterpretation,
            11 => Self::Ensemble,
            _ => return None,
        })
    }

    /// The kind with its indefinite article, such as "an interpreted structure".
    fn with_article(self) -> String {
        let article = match self {
            Self::StructureInterpretation | Self::StructureSetInterpretation | Self::Ensemble => {
                "an"
            }
            _ => "a",
        };
        format!("{article} {self}")
    }
}

impl fmt::Display for NativeKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Molecule => "molecule",
            Self::Topology => "topology",
            Self::Structure => "structure",
            Self::StructureInterpretation => "interpreted structure",
            Self::StructureSet => "structure set",
            Self::StructureSetInterpretation => "interpreted structure set",
            Self::Ensemble => "ensemble",
        })
    }
}

/// The object a native writer encodes, borrowed from its owner.
///
/// Converts from a borrowed [`Molecule`], [`Topology`], [`Structure`],
/// [`StructureSet`], [`StructureSetDescription`], [`Ensemble`],
/// [`EnsembleDescription`], [`StructureInterpretation`], or
/// [`StructureSetInterpretation`]. [`Self::interpreted_structure`] and
/// [`Self::interpreted_structure_set`] pair a structure or structure set with
/// reports held apart from it.
#[derive(Debug, Clone, Copy)]
pub struct NativeSource<'a> {
    pub(super) content: NativeContent<'a>,
}

#[derive(Debug, Clone, Copy)]
pub(super) enum NativeContent<'a> {
    Molecule(&'a Molecule),
    Topology(&'a Topology),
    Structure(&'a Structure),
    InterpretedStructure(&'a Structure, &'a StructureInterpretationReport),
    StructureSet(&'a StructureSet),
    InterpretedStructureSet(&'a StructureSet, &'a [StructureInterpretationReport]),
    StructureSetDescription(&'a StructureSetDescription),
    Ensemble(&'a Ensemble),
    EnsembleDescription(&'a EnsembleDescription),
}

impl<'a> NativeSource<'a> {
    /// A structure with the report of its interpretation.
    ///
    /// The report is stored as it is. Its IDs name atoms, chains, and
    /// residues of the structure it was interpreted with.
    pub const fn interpreted_structure(
        structure: &'a Structure,
        report: &'a StructureInterpretationReport,
    ) -> Self {
        Self {
            content: NativeContent::InterpretedStructure(structure, report),
        }
    }

    /// A structure set with one interpretation report per geometry, in geometry
    /// order, as [`StructureSetInterpretation::reports`] gives them.
    ///
    /// Writing fails with [`DiagnosticKind::Inconsistent`] unless there is
    /// one report per geometry.
    pub const fn interpreted_structure_set(
        set: &'a StructureSet,
        reports: &'a [StructureInterpretationReport],
    ) -> Self {
        Self {
            content: NativeContent::InterpretedStructureSet(set, reports),
        }
    }

    /// The kind of object this source encodes.
    pub const fn kind(&self) -> NativeKind {
        match self.content {
            NativeContent::Molecule(_) => NativeKind::Molecule,
            NativeContent::Topology(_) => NativeKind::Topology,
            NativeContent::Structure(_) => NativeKind::Structure,
            NativeContent::InterpretedStructure(..) => NativeKind::StructureInterpretation,
            NativeContent::StructureSet(_) | NativeContent::StructureSetDescription(_) => {
                NativeKind::StructureSet
            }
            NativeContent::InterpretedStructureSet(..) => NativeKind::StructureSetInterpretation,
            NativeContent::Ensemble(_) | NativeContent::EnsembleDescription(_) => {
                NativeKind::Ensemble
            }
        }
    }
}

macro_rules! native_source_from {
    ($($variant:ident),*) => {$(
        impl<'a> From<&'a $variant> for NativeSource<'a> {
            fn from(value: &'a $variant) -> Self {
                Self {
                    content: NativeContent::$variant(value),
                }
            }
        }
    )*};
}

native_source_from!(
    Molecule,
    Topology,
    Structure,
    StructureSet,
    StructureSetDescription,
    Ensemble,
    EnsembleDescription
);

impl<'a> From<&'a StructureInterpretation> for NativeSource<'a> {
    fn from(value: &'a StructureInterpretation) -> Self {
        Self::interpreted_structure(value.structure(), value.report())
    }
}

impl<'a> From<&'a StructureSetInterpretation> for NativeSource<'a> {
    fn from(value: &'a StructureSetInterpretation) -> Self {
        Self::interpreted_structure_set(value.structure_set(), value.reports())
    }
}

/// One object decoded from native input.
///
/// A decoded topology, including the shared topology of a structure, is a
/// new publication with a fresh layout identity and the stored semantic IDs.
#[derive(Debug, Clone)]
#[non_exhaustive]
pub enum NativeObject {
    /// A molecule.
    Molecule(Box<Molecule>),
    /// A topology.
    Topology(Topology),
    /// A structure with its topology.
    Structure(Structure),
    /// A structure with its topology and the report of its interpretation.
    StructureInterpretation(Box<StructureInterpretation>),
    /// A structure set with its topology and geometries.
    StructureSet(StructureSet),
    /// A structure set with its topology and one interpretation report per
    /// geometry.
    StructureSetInterpretation(Box<StructureSetInterpretation>),
    /// A structure set record that stores the number of its geometries, not the
    /// geometries: its topology, properties, and order.
    StructureSetDescription(Box<StructureSetDescription>),
    /// An ensemble with its set's topology and geometries.
    Ensemble(Box<Ensemble>),
    /// An ensemble record that stores the number of its geometries, not the
    /// geometries: its set's topology, properties, and order, and its
    /// distribution.
    EnsembleDescription(Box<EnsembleDescription>),
}

impl NativeObject {
    /// The kind of this object.
    pub const fn kind(&self) -> NativeKind {
        match self {
            Self::Molecule(_) => NativeKind::Molecule,
            Self::Topology(_) => NativeKind::Topology,
            Self::Structure(_) => NativeKind::Structure,
            Self::StructureInterpretation(_) => NativeKind::StructureInterpretation,
            Self::StructureSet(_) | Self::StructureSetDescription(_) => NativeKind::StructureSet,
            Self::StructureSetInterpretation(_) => NativeKind::StructureSetInterpretation,
            Self::Ensemble(_) | Self::EnsembleDescription(_) => NativeKind::Ensemble,
        }
    }

    fn kind_error(&self, expected: NativeKind) -> FormatError {
        let found = self.kind();
        FormatError::interpret(
            Format::Native,
            DiagnosticKind::Inconsistent,
            format!(
                "The input holds {}, not {}.",
                found.with_article(),
                expected.with_article()
            ),
        )
        .with_detail(FormatErrorDetail::NativeKind { expected, found })
    }

    /// The decoded molecule.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeKind`] when the object is of another kind.
    pub fn into_molecule(self) -> Result<Molecule, FormatError> {
        match self {
            Self::Molecule(value) => Ok(*value),
            other => Err(other.kind_error(NativeKind::Molecule)),
        }
    }

    /// The decoded topology.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeKind`] when the object is of another kind.
    pub fn into_topology(self) -> Result<Topology, FormatError> {
        match self {
            Self::Topology(value) => Ok(value),
            other => Err(other.kind_error(NativeKind::Topology)),
        }
    }

    /// The decoded structure, of a structure or of an interpreted structure,
    /// whose report it drops.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeKind`] when the object is of another kind.
    pub fn into_structure(self) -> Result<Structure, FormatError> {
        match self {
            Self::Structure(value) => Ok(value),
            Self::StructureInterpretation(value) => Ok(value.into_structure()),
            other => Err(other.kind_error(NativeKind::Structure)),
        }
    }

    /// The decoded structure with the report of its interpretation.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeKind`] when the object is of another kind,
    /// such as a structure stored without a report.
    pub fn into_interpretation(self) -> Result<StructureInterpretation, FormatError> {
        match self {
            Self::StructureInterpretation(value) => Ok(*value),
            other => Err(other.kind_error(NativeKind::StructureInterpretation)),
        }
    }

    /// The decoded structure set, of a structure set, of an interpreted
    /// structure set, whose reports it drops, or of an ensemble, whose
    /// distribution it drops.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeKind`] when the object is of another kind,
    /// and with [`FormatErrorDetail::NativeExternalGeometries`] for a record
    /// that stores its geometries elsewhere.
    pub fn into_structure_set(self) -> Result<StructureSet, FormatError> {
        match self {
            Self::StructureSet(value) => Ok(value),
            Self::StructureSetInterpretation(value) => Ok(value.into_structure_set()),
            Self::Ensemble(ensemble) => Ok(ensemble.into_parts().0),
            Self::StructureSetDescription(description) => Err(external_geometries(&description)),
            Self::EnsembleDescription(description) => Err(external_geometries(description.set())),
            other => Err(other.kind_error(NativeKind::StructureSet)),
        }
    }

    /// The decoded ensemble.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeKind`] when the object is of another kind,
    /// such as a structure set stored without a distribution, and with
    /// [`FormatErrorDetail::NativeExternalGeometries`] for a record that
    /// stores its geometries elsewhere.
    pub fn into_ensemble(self) -> Result<Ensemble, FormatError> {
        match self {
            Self::Ensemble(ensemble) => Ok(*ensemble),
            Self::EnsembleDescription(description) => Err(external_geometries(description.set())),
            other => Err(other.kind_error(NativeKind::Ensemble)),
        }
    }

    /// What the decoded ensemble states besides its geometries: of a record
    /// that stores only their number, or of any ensemble, whose geometries
    /// it drops.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeKind`] when the object is no ensemble.
    pub fn into_ensemble_description(self) -> Result<EnsembleDescription, FormatError> {
        match self {
            Self::EnsembleDescription(description) => Ok(*description),
            Self::Ensemble(ensemble) => Ok(EnsembleDescription::of(&ensemble)),
            other => Err(other.kind_error(NativeKind::Ensemble)),
        }
    }

    /// The decoded structure set with one interpretation report per geometry.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeKind`] when the object is of another kind,
    /// such as a structure set stored without reports.
    pub fn into_structure_set_interpretation(
        self,
    ) -> Result<StructureSetInterpretation, FormatError> {
        match self {
            Self::StructureSetInterpretation(value) => Ok(*value),
            other => Err(other.kind_error(NativeKind::StructureSetInterpretation)),
        }
    }

    /// What the decoded structure set states besides its geometries: of a record
    /// that stores only their number, or of any structure set, whose geometries
    /// and reports it drops, or of an ensemble, whose distribution it also
    /// drops.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeKind`] when the object is no structure set.
    pub fn into_structure_set_description(self) -> Result<StructureSetDescription, FormatError> {
        match self {
            Self::StructureSetDescription(value) => Ok(*value),
            Self::StructureSet(set) => Ok(StructureSetDescription::of(&set)),
            Self::StructureSetInterpretation(value) => {
                Ok(StructureSetDescription::of(value.structure_set()))
            }
            Self::Ensemble(ensemble) => Ok(StructureSetDescription::of(ensemble.set())),
            Self::EnsembleDescription(description) => Ok(description.into_parts().0),
            other => Err(other.kind_error(NativeKind::StructureSet)),
        }
    }
}

fn external_geometries(description: &StructureSetDescription) -> FormatError {
    let geometries = description.geometry_count();
    FormatError::interpret(
        Format::Native,
        DiagnosticKind::Inconsistent,
        format!(
            "The input's structure set stores its {} elsewhere.",
            crate::diagnostics::counted(geometries, "geometry", "geometries")
        ),
    )
    .with_detail(FormatErrorDetail::NativeExternalGeometries { geometries })
}

/// Options of the native writer.
///
/// The default stores installed perception, so a reader can reinstall it
/// instead of perceiving again.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NativeWriteOptions {
    perception: bool,
}

impl Default for NativeWriteOptions {
    fn default() -> Self {
        Self { perception: true }
    }
}

impl NativeWriteOptions {
    /// Whether installed perception is stored with its provenance.
    pub const fn perception(&self) -> bool {
        self.perception
    }

    /// Sets [`Self::perception`]. Omitting perception stores represented
    /// chemistry only; readers then publish molecules without perception.
    #[must_use]
    pub const fn with_perception(mut self, perception: bool) -> Self {
        self.perception = perception;
        self
    }
}

/// Options and resource limits of the native reader.
///
/// Declared counts are checked against these limits before anything is
/// allocated; exceeding one is a [`FormatError`] of kind
/// [`DiagnosticKind::Limit`] carrying
/// [`LimitExceeded`](crate::diagnostics::LimitExceeded).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NativeReadOptions {
    perception: bool,
    max_input_bytes: usize,
    max_atoms: usize,
    max_bonds: usize,
    max_items: usize,
    max_string_bytes: usize,
}

impl Default for NativeReadOptions {
    fn default() -> Self {
        Self {
            perception: true,
            max_input_bytes: DEFAULT_NATIVE_MAX_INPUT_BYTES,
            max_atoms: DEFAULT_NATIVE_MAX_ATOMS,
            max_bonds: DEFAULT_NATIVE_MAX_BONDS,
            max_items: DEFAULT_NATIVE_MAX_ITEMS,
            max_string_bytes: DEFAULT_NATIVE_MAX_STRING_BYTES,
        }
    }
}

fn nonzero_option(name: &str, value: usize) -> Result<(), FormatError> {
    if value == 0 {
        return Err(FormatError::parse(
            Format::Native,
            DiagnosticKind::Limit,
            format!(
                "The {} is 0. Set a positive limit.",
                match name {
                    "max_input_bytes" => "input size limit",
                    "max_atoms" => "atom limit",
                    "max_bonds" => "bond limit",
                    "max_items" => "item limit",
                    "max_string_bytes" => "string size limit",
                    _ => "limit",
                }
            ),
        ));
    }
    Ok(())
}

macro_rules! native_limit {
    ($(#[$meta:meta])* $field:ident, $setter:ident) => {
        $(#[$meta])*
        pub const fn $field(&self) -> usize {
            self.$field
        }

        #[doc = concat!("Sets [`Self::", stringify!($field), "`].")]
        ///
        /// # Errors
        ///
        /// [`DiagnosticKind::Limit`] when the bound is zero.
        pub fn $setter(mut self, $field: usize) -> Result<Self, FormatError> {
            nonzero_option(stringify!($field), $field)?;
            self.$field = $field;
            Ok(self)
        }
    };
}

impl NativeReadOptions {
    /// Whether stored perception is installed, after validation against the
    /// decoded chemistry through
    /// [`Molecule::install_perception`](crate::molecule::Molecule::install_perception),
    /// rather than discarded.
    pub const fn perception(&self) -> bool {
        self.perception
    }

    /// Sets [`Self::perception`].
    #[must_use]
    pub const fn with_perception(mut self, perception: bool) -> Self {
        self.perception = perception;
        self
    }

    native_limit!(
        /// Maximum bytes read, trailer included.
        max_input_bytes,
        with_max_input_bytes
    );
    native_limit!(
        /// Maximum atoms of one molecule or topology, counting every instance
        /// of a reused definition.
        max_atoms,
        with_max_atoms
    );
    native_limit!(
        /// Maximum bonds of one molecule or topology, counting every instance
        /// of a reused definition.
        max_bonds,
        with_max_bonds
    );
    native_limit!(
        /// Maximum geometries of a structure set, counting those a record stores
        /// only the number of, and its distribution's states, and maximum
        /// entries of any other stored list that atoms or bonds do not
        /// bound: property columns and owner values, hierarchy chains and
        /// residues, and the residues listed by all chains together.
        max_items,
        with_max_items
    );
    native_limit!(
        /// Maximum bytes of one stored string.
        max_string_bytes,
        with_max_string_bytes
    );
}
