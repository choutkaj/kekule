//! Kekule's native binary serialization: exact, versioned persistence of one
//! molecule, topology, model, or sample set, and of a model or sample set
//! with its structure interpretation reports.
//!
//! The byte layout is specified in `docs/native-format.md`; the encoder and
//! decoder here implement that specification and nothing else.

mod crc;
mod decode;
mod description;
mod encode;

use std::fmt;

use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail};
use crate::io::{
    StructureEnsembleInterpretation, StructureInterpretation, StructureInterpretationReport,
    StructureSampleSetInterpretation,
};
use crate::molecule::Molecule;
use crate::structure::{
    ConversionTarget, Ensemble, Model, SampleSet, SampleSetError, Series, Trajectory,
};
use crate::topology::Topology;

pub(crate) use decode::read_native;
pub use description::SampleSetDescription;
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
pub const NATIVE_FORMAT_VERSION: u32 = 5;

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
    /// A [`Model`] with the report of its interpretation, a
    /// [`StructureInterpretation`].
    StructureInterpretation,
    /// A [`SampleSet`], with its samples or with only their number
    /// ([`SampleSetDescription`]).
    SampleSet,
    /// A [`SampleSet`] with one interpretation report per sample, a
    /// [`StructureSampleSetInterpretation`].
    SampleSetInterpretation,
}

impl NativeKind {
    pub(super) const fn code(self) -> u8 {
        match self {
            Self::Molecule => 1,
            Self::Topology => 2,
            Self::Model => 3,
            Self::StructureInterpretation => 7,
            Self::SampleSet => 9,
            Self::SampleSetInterpretation => 10,
        }
    }

    pub(super) const fn from_code(code: u8) -> Option<Self> {
        Some(match code {
            1 => Self::Molecule,
            2 => Self::Topology,
            3 => Self::Model,
            7 => Self::StructureInterpretation,
            9 => Self::SampleSet,
            10 => Self::SampleSetInterpretation,
            _ => return None,
        })
    }

    /// The kind with its indefinite article, such as "an interpreted model".
    fn with_article(self) -> String {
        let article = match self {
            Self::StructureInterpretation | Self::SampleSetInterpretation => "an",
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
            Self::Model => "model",
            Self::StructureInterpretation => "interpreted model",
            Self::SampleSet => "sample set",
            Self::SampleSetInterpretation => "interpreted sample set",
        })
    }
}

/// The object a native writer encodes, borrowed from its owner.
///
/// Converts from a borrowed [`Molecule`], [`Topology`], [`Model`],
/// [`SampleSet`], [`SampleSetDescription`], [`StructureInterpretation`], or
/// [`StructureSampleSetInterpretation`]. An [`Ensemble`], [`Trajectory`], or
/// [`Series`] is written as the sample set it converts to
/// ([`SampleSet::from_ensemble`], [`SampleSet::from_trajectory`],
/// [`SampleSet::from_series`]) without copying it, and a
/// [`StructureEnsembleInterpretation`] as an interpreted sample set.
/// [`Self::interpreted_model`], [`Self::interpreted_sample_set`], and
/// [`Self::interpreted_ensemble`] pair a model, sample set, or ensemble with
/// reports held apart from it.
#[derive(Debug, Clone, Copy)]
pub struct NativeSource<'a> {
    pub(super) content: NativeContent<'a>,
}

#[derive(Debug, Clone, Copy)]
pub(super) enum NativeContent<'a> {
    Molecule(&'a Molecule),
    Topology(&'a Topology),
    Model(&'a Model),
    InterpretedModel(&'a Model, &'a StructureInterpretationReport),
    SampleSet(&'a SampleSet),
    InterpretedSampleSet(&'a SampleSet, &'a [StructureInterpretationReport]),
    SampleSetDescription(&'a SampleSetDescription),
    Ensemble(&'a Ensemble),
    InterpretedEnsemble(&'a Ensemble, &'a [StructureInterpretationReport]),
    Trajectory(&'a Trajectory),
    Series(&'a Series),
}

impl<'a> NativeSource<'a> {
    /// A model with the report of its interpretation.
    ///
    /// The report is stored as it is. Its IDs name atoms, chains, and
    /// residues of the model it was interpreted with.
    pub const fn interpreted_model(
        model: &'a Model,
        report: &'a StructureInterpretationReport,
    ) -> Self {
        Self {
            content: NativeContent::InterpretedModel(model, report),
        }
    }

    /// A sample set with one interpretation report per sample, in sample
    /// order, as [`StructureSampleSetInterpretation::reports`] gives them.
    ///
    /// Writing fails with [`DiagnosticKind::Inconsistent`] unless there is
    /// one report per sample.
    pub const fn interpreted_sample_set(
        set: &'a SampleSet,
        reports: &'a [StructureInterpretationReport],
    ) -> Self {
        Self {
            content: NativeContent::InterpretedSampleSet(set, reports),
        }
    }

    /// An ensemble with one interpretation report per sample, in sample
    /// order, as [`StructureEnsembleInterpretation::reports`] gives them,
    /// written as the interpreted sample set it converts to.
    ///
    /// Writing fails with [`DiagnosticKind::Inconsistent`] unless there is
    /// one report per sample.
    pub const fn interpreted_ensemble(
        ensemble: &'a Ensemble,
        reports: &'a [StructureInterpretationReport],
    ) -> Self {
        Self {
            content: NativeContent::InterpretedEnsemble(ensemble, reports),
        }
    }

    /// The kind of object this source encodes.
    pub const fn kind(&self) -> NativeKind {
        match self.content {
            NativeContent::Molecule(_) => NativeKind::Molecule,
            NativeContent::Topology(_) => NativeKind::Topology,
            NativeContent::Model(_) => NativeKind::Model,
            NativeContent::InterpretedModel(..) => NativeKind::StructureInterpretation,
            NativeContent::SampleSet(_)
            | NativeContent::SampleSetDescription(_)
            | NativeContent::Ensemble(_)
            | NativeContent::Trajectory(_)
            | NativeContent::Series(_) => NativeKind::SampleSet,
            NativeContent::InterpretedSampleSet(..) | NativeContent::InterpretedEnsemble(..) => {
                NativeKind::SampleSetInterpretation
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
    Model,
    SampleSet,
    SampleSetDescription,
    Ensemble,
    Trajectory,
    Series
);

impl<'a> From<&'a StructureInterpretation> for NativeSource<'a> {
    fn from(value: &'a StructureInterpretation) -> Self {
        Self::interpreted_model(value.model(), value.report())
    }
}

impl<'a> From<&'a StructureSampleSetInterpretation> for NativeSource<'a> {
    fn from(value: &'a StructureSampleSetInterpretation) -> Self {
        Self::interpreted_sample_set(value.sample_set(), value.reports())
    }
}

impl<'a> From<&'a StructureEnsembleInterpretation> for NativeSource<'a> {
    fn from(value: &'a StructureEnsembleInterpretation) -> Self {
        Self::interpreted_ensemble(value.ensemble(), value.reports())
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
    /// A model with its topology.
    Model(Model),
    /// A model with its topology and the report of its interpretation.
    StructureInterpretation(Box<StructureInterpretation>),
    /// A sample set with its topology and samples.
    SampleSet(SampleSet),
    /// A sample set with its topology and one interpretation report per
    /// sample.
    SampleSetInterpretation(Box<StructureSampleSetInterpretation>),
    /// A sample set record that stores the number of its samples, not the
    /// samples: its topology, properties, order, and distribution.
    SampleSetDescription(Box<SampleSetDescription>),
}

impl NativeObject {
    /// The kind of this object.
    pub const fn kind(&self) -> NativeKind {
        match self {
            Self::Molecule(_) => NativeKind::Molecule,
            Self::Topology(_) => NativeKind::Topology,
            Self::Model(_) => NativeKind::Model,
            Self::StructureInterpretation(_) => NativeKind::StructureInterpretation,
            Self::SampleSet(_) | Self::SampleSetDescription(_) => NativeKind::SampleSet,
            Self::SampleSetInterpretation(_) => NativeKind::SampleSetInterpretation,
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

    /// The decoded model, of a model or of an interpreted model, whose
    /// report it drops.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeKind`] when the object is of another kind.
    pub fn into_model(self) -> Result<Model, FormatError> {
        match self {
            Self::Model(value) => Ok(value),
            Self::StructureInterpretation(value) => Ok(value.into_model()),
            other => Err(other.kind_error(NativeKind::Model)),
        }
    }

    /// The decoded model with the report of its interpretation.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeKind`] when the object is of another kind,
    /// such as a model stored without a report.
    pub fn into_interpretation(self) -> Result<StructureInterpretation, FormatError> {
        match self {
            Self::StructureInterpretation(value) => Ok(*value),
            other => Err(other.kind_error(NativeKind::StructureInterpretation)),
        }
    }

    /// The decoded sample set, of a sample set or of an interpreted sample
    /// set, whose reports it drops.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeKind`] when the object is of another kind,
    /// and with [`FormatErrorDetail::NativeExternalSamples`] for a record
    /// that stores its samples elsewhere.
    pub fn into_sample_set(self) -> Result<SampleSet, FormatError> {
        match self {
            Self::SampleSet(value) => Ok(value),
            Self::SampleSetInterpretation(value) => Ok(value.into_sample_set()),
            Self::SampleSetDescription(description) => Err(external_samples(&description)),
            other => Err(other.kind_error(NativeKind::SampleSet)),
        }
    }

    /// The decoded sample set with one interpretation report per sample.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeKind`] when the object is of another kind,
    /// such as a sample set stored without reports.
    pub fn into_sample_set_interpretation(
        self,
    ) -> Result<StructureSampleSetInterpretation, FormatError> {
        match self {
            Self::SampleSetInterpretation(value) => Ok(*value),
            other => Err(other.kind_error(NativeKind::SampleSetInterpretation)),
        }
    }

    /// What the decoded sample set states besides its samples: of a record
    /// that stores only their number, or of any sample set, whose samples
    /// and reports it drops.
    ///
    /// # Errors
    ///
    /// A [`FormatError`] of kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeKind`] when the object is no sample set.
    pub fn into_sample_set_description(self) -> Result<SampleSetDescription, FormatError> {
        match self {
            Self::SampleSetDescription(value) => Ok(*value),
            Self::SampleSet(set) => Ok(SampleSetDescription::of(&set)),
            Self::SampleSetInterpretation(value) => {
                Ok(SampleSetDescription::of(value.sample_set()))
            }
            other => Err(other.kind_error(NativeKind::SampleSet)),
        }
    }

    /// The decoded sample set as an ensemble, of a sample set or of an
    /// interpreted sample set, whose reports it drops; see
    /// [`SampleSet::into_ensemble`].
    ///
    /// # Errors
    ///
    /// The errors of [`Self::into_sample_set`], and a [`FormatError`] of
    /// kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeSampleSetShape`] when the set is not an
    /// ensemble, such as an ordered one; its source says why.
    pub fn into_ensemble(self) -> Result<Ensemble, FormatError> {
        self.into_sample_set()?
            .into_ensemble()
            .map_err(|error| shape_error(ConversionTarget::Ensemble, error))
    }

    /// The decoded sample set as an ensemble with one interpretation report
    /// per sample; see [`SampleSet::into_ensemble`].
    ///
    /// # Errors
    ///
    /// The errors of [`Self::into_sample_set_interpretation`], and those of
    /// [`Self::into_ensemble`] for a set that is not an ensemble.
    pub fn into_ensemble_interpretation(
        self,
    ) -> Result<StructureEnsembleInterpretation, FormatError> {
        let (set, reports) = self.into_sample_set_interpretation()?.into_parts();
        let ensemble = set
            .into_ensemble()
            .map_err(|error| shape_error(ConversionTarget::Ensemble, error))?;
        Ok(StructureEnsembleInterpretation { ensemble, reports })
    }

    /// The decoded sample set as a trajectory; see
    /// [`SampleSet::into_trajectory`].
    ///
    /// # Errors
    ///
    /// The errors of [`Self::into_sample_set`], and a [`FormatError`] of
    /// kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeSampleSetShape`] when the set is not a
    /// trajectory, such as one without a time order.
    pub fn into_trajectory(self) -> Result<Trajectory, FormatError> {
        self.into_sample_set()?
            .into_trajectory()
            .map_err(|error| shape_error(ConversionTarget::Trajectory, error))
    }

    /// The decoded sample set as a series; see [`SampleSet::into_series`].
    ///
    /// # Errors
    ///
    /// The errors of [`Self::into_sample_set`], and a [`FormatError`] of
    /// kind [`DiagnosticKind::Inconsistent`] with
    /// [`FormatErrorDetail::NativeSampleSetShape`] when the set is not a
    /// series, such as one without a step or coordinate order.
    pub fn into_series(self) -> Result<Series, FormatError> {
        self.into_sample_set()?
            .into_series()
            .map_err(|error| shape_error(ConversionTarget::Series, error))
    }
}

fn external_samples(description: &SampleSetDescription) -> FormatError {
    let samples = description.sample_count();
    FormatError::interpret(
        Format::Native,
        DiagnosticKind::Inconsistent,
        format!(
            "The input's sample set stores its {} elsewhere.",
            crate::diagnostics::counted(samples, "sample", "samples")
        ),
    )
    .with_detail(FormatErrorDetail::NativeExternalSamples { samples })
}

fn shape_error(target: ConversionTarget, error: SampleSetError) -> FormatError {
    FormatError::interpret(
        Format::Native,
        DiagnosticKind::Inconsistent,
        format!("The input's sample set is not {target}."),
    )
    .with_detail(FormatErrorDetail::NativeSampleSetShape { target })
    .with_source(error)
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
        /// Maximum samples of a sample set, counting those a record stores
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
