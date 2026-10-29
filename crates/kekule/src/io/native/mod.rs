//! Kekule's native binary serialization: exact, versioned persistence of one
//! molecule, topology, model, ensemble, trajectory, or series.
//!
//! The byte layout is specified in `docs/native-format.md`; the encoder and
//! decoder here implement that specification and nothing else.

mod crc;
mod decode;
mod encode;

use std::fmt;

use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail};
use crate::molecule::Molecule;
use crate::structure::{Ensemble, Model, Series, Trajectory};
use crate::topology::Topology;

pub(crate) use decode::read_native;
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
pub const NATIVE_FORMAT_VERSION: u32 = 2;

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
    /// A [`Model`].
    Model,
    /// An [`Ensemble`].
    Ensemble,
    /// A [`Trajectory`].
    Trajectory,
    /// A [`Series`].
    Series,
}

impl NativeKind {
    pub(super) const fn code(self) -> u8 {
        match self {
            Self::Molecule => 1,
            Self::Topology => 2,
            Self::Model => 3,
            Self::Ensemble => 4,
            Self::Trajectory => 5,
            Self::Series => 6,
        }
    }

    pub(super) const fn from_code(code: u8) -> Option<Self> {
        Some(match code {
            1 => Self::Molecule,
            2 => Self::Topology,
            3 => Self::Model,
            4 => Self::Ensemble,
            5 => Self::Trajectory,
            6 => Self::Series,
            _ => return None,
        })
    }
}

impl fmt::Display for NativeKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Molecule => "molecule",
            Self::Topology => "topology",
            Self::Model => "model",
            Self::Ensemble => "ensemble",
            Self::Trajectory => "trajectory",
            Self::Series => "series",
        })
    }
}

/// The object a native writer encodes, borrowed from its owner.
///
/// Converts from a borrowed [`Molecule`], [`Topology`], [`Model`],
/// [`Ensemble`], [`Trajectory`], or [`Series`].
#[derive(Debug, Clone, Copy)]
pub struct NativeSource<'a> {
    pub(super) content: NativeContent<'a>,
}

#[derive(Debug, Clone, Copy)]
pub(super) enum NativeContent<'a> {
    Molecule(&'a Molecule),
    Topology(&'a Topology),
    Model(&'a Model),
    Ensemble(&'a Ensemble),
    Trajectory(&'a Trajectory),
    Series(&'a Series),
}

impl NativeSource<'_> {
    /// The kind of object this source encodes.
    pub const fn kind(&self) -> NativeKind {
        match self.content {
            NativeContent::Molecule(_) => NativeKind::Molecule,
            NativeContent::Topology(_) => NativeKind::Topology,
            NativeContent::Model(_) => NativeKind::Model,
            NativeContent::Ensemble(_) => NativeKind::Ensemble,
            NativeContent::Trajectory(_) => NativeKind::Trajectory,
            NativeContent::Series(_) => NativeKind::Series,
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

native_source_from!(Molecule, Topology, Model, Ensemble, Trajectory, Series);

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
    /// A model with its topology.
    Model(Model),
    /// An ensemble with its topology.
    Ensemble(Ensemble),
    /// A trajectory with its topology.
    Trajectory(Trajectory),
    /// A series with its topology.
    Series(Series),
}

macro_rules! native_projection {
    ($(#[$meta:meta])* $name:ident, $variant:ident) => {
        $(#[$meta])*
        ///
        /// # Errors
        ///
        /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
        /// [`FormatErrorDetail::NativeKind`] when the object is of another
        /// kind.
        pub fn $name(self) -> Result<$variant, FormatError> {
            match self {
                Self::$variant(value) => Ok(value),
                other => Err(other.kind_error(NativeKind::$variant)),
            }
        }
    };
}

impl NativeObject {
    /// The kind of this object.
    pub const fn kind(&self) -> NativeKind {
        match self {
            Self::Molecule(_) => NativeKind::Molecule,
            Self::Topology(_) => NativeKind::Topology,
            Self::Model(_) => NativeKind::Model,
            Self::Ensemble(_) => NativeKind::Ensemble,
            Self::Trajectory(_) => NativeKind::Trajectory,
            Self::Series(_) => NativeKind::Series,
        }
    }

    fn kind_error(&self, expected: NativeKind) -> FormatError {
        let found = self.kind();
        FormatError::interpret(
            Format::Native,
            DiagnosticKind::Inconsistent,
            format!("native input holds a {found}, not a {expected}"),
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

    native_projection!(
        /// The decoded topology.
        into_topology,
        Topology
    );
    native_projection!(
        /// The decoded model.
        into_model,
        Model
    );
    native_projection!(
        /// The decoded ensemble.
        into_ensemble,
        Ensemble
    );
    native_projection!(
        /// The decoded trajectory.
        into_trajectory,
        Trajectory
    );
    native_projection!(
        /// The decoded series.
        into_series,
        Series
    );
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
            format!("native option {name} must be greater than zero"),
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
        /// Maximum ensemble samples or states, trajectory frames, or series
        /// items, and maximum entries of any other stored list that atoms or
        /// bonds do not bound: property columns and owner values, hierarchy
        /// chains and residues, and the residues listed by all chains
        /// together.
        max_items,
        with_max_items
    );
    native_limit!(
        /// Maximum bytes of one stored string.
        max_string_bytes,
        with_max_string_bytes
    );
}
