//! Public writer API and orchestration of planning, preparation, and emission.
use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail, StructureField};
use crate::geometry::PeriodicCell;
use crate::io::mmcif_cell::{cell_parameters, CellParameterError, CellParameters};
use crate::io::structure_interpret::{
    EntityKind, StructureEnsembleInterpretation, StructureInterpretation,
    StructureInterpretationReport,
};
use crate::io::write_error::{check_multi_model_ensemble, check_structure_state, Omissions};
use crate::io::{UnrepresentablePolicy, WriteOmission, WriteReport};
use crate::molecule::BondOrder;
use crate::structure::{AsModelView, Ensemble, Model, ModelView};
use crate::topology::{InstanceAtomId, InstanceBondId, MoleculeClass, MoleculeInstanceId};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::io::Write;

mod emit;
mod entities;
mod prepare;
use emit::{render_model_to, write_atom_rows, write_block_end, write_block_start};
use entities::{
    ensemble_entity_plan, entity_plan_from_report, generic_entity_plan,
    normalize_entity_classifications,
};
use prepare::{prepare_model, validate_ensemble_sample};

const MAX_COORDINATE_PRECISION: usize = 15;

/// Options of the mmCIF writer.
///
/// The default writes a block named `model` with coordinates to 3 decimal
/// places and fails on ensemble metadata, a periodic cell, or fixed
/// implicit-hydrogen counts mmCIF cannot carry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifWriteOptions {
    block_name: String,
    coordinate_precision: usize,
    ensemble_metadata: UnrepresentablePolicy,
    frame_metadata: UnrepresentablePolicy,
    cell: UnrepresentablePolicy,
    implicit_hydrogens: UnrepresentablePolicy,
}

impl MmcifWriteOptions {
    /// The name of the written data block, after `data_`; `model` by default.
    /// Several sources are written as blocks `{name}_1`, `{name}_2`, and so
    /// on, numbered by position; a source with its own
    /// [`MmcifBlockSource::block_name`] takes that name instead.
    pub fn block_name(&self) -> &str {
        &self.block_name
    }

    /// Sets [`Self::block_name`].
    ///
    /// # Errors
    ///
    /// A write [`FormatError`] with [`MmcifWriteIssue::InvalidBlockName`]
    /// unless the name is a nonempty run of ASCII letters, digits, `_`, `-`
    /// and `.`.
    pub fn with_block_name(mut self, block_name: impl Into<String>) -> Result<Self, FormatError> {
        let block_name = block_name.into();
        validate_block_name(&block_name)?;
        self.block_name = block_name;
        Ok(self)
    }

    /// Decimal places of written Cartesian coordinates, in Å; 3 by default and
    /// at most 15.
    pub const fn coordinate_precision(&self) -> usize {
        self.coordinate_precision
    }

    /// Sets [`Self::coordinate_precision`].
    ///
    /// # Errors
    ///
    /// A write [`FormatError`] with
    /// [`MmcifWriteIssue::CoordinatePrecisionTooLarge`] above 15 digits.
    pub fn with_coordinate_precision(
        mut self,
        coordinate_precision: usize,
    ) -> Result<Self, FormatError> {
        if coordinate_precision > MAX_COORDINATE_PRECISION {
            return Err(fail(MmcifWriteIssue::CoordinatePrecisionTooLarge(
                coordinate_precision,
            )));
        }
        self.coordinate_precision = coordinate_precision;
        Ok(self)
    }

    /// Ensemble sample weights, states, state energies, and the temperature,
    /// which mmCIF cannot carry: by default an ensemble with unequal weights,
    /// more than one state, any energy, or a temperature fails to write, and
    /// so does a model block of one ensemble sample, whose weight is relative
    /// to samples the block does not hold.
    pub const fn ensemble_metadata(&self) -> UnrepresentablePolicy {
        self.ensemble_metadata
    }

    /// Sets [`Self::ensemble_metadata`].
    #[must_use]
    pub fn with_ensemble_metadata(mut self, ensemble_metadata: UnrepresentablePolicy) -> Self {
        self.ensemble_metadata = ensemble_metadata;
        self
    }

    /// A trajectory frame's time, step, velocities, and forces, which mmCIF
    /// cannot carry: by default a model block of a frame with any of them
    /// fails to write.
    pub const fn frame_metadata(&self) -> UnrepresentablePolicy {
        self.frame_metadata
    }

    /// Sets [`Self::frame_metadata`].
    #[must_use]
    pub fn with_frame_metadata(mut self, frame_metadata: UnrepresentablePolicy) -> Self {
        self.frame_metadata = frame_metadata;
        self
    }

    /// A periodic cell mmCIF cannot carry: one not periodic along all three
    /// axes, not in the standard crystallographic orientation, or differing
    /// between ensemble samples. By default such a cell fails to write;
    /// representable cells are always written as `_cell` in space group P 1.
    pub const fn cell(&self) -> UnrepresentablePolicy {
        self.cell
    }

    /// Sets [`Self::cell`].
    #[must_use]
    pub fn with_cell(mut self, cell: UnrepresentablePolicy) -> Self {
        self.cell = cell;
        self
    }

    /// Fixed implicit-hydrogen counts, which mmCIF cannot carry: readers
    /// infer implicit hydrogens from valence. By default an atom with one
    /// fails to write; [`UnrepresentablePolicy::Omit`] writes it without the
    /// count and reports [`WriteOmission::ImplicitHydrogens`]. Structures
    /// whose files list every atom, such as XYZ interpretations, fix their
    /// counts at zero.
    pub const fn implicit_hydrogens(&self) -> UnrepresentablePolicy {
        self.implicit_hydrogens
    }

    /// Sets [`Self::implicit_hydrogens`].
    #[must_use]
    pub fn with_implicit_hydrogens(mut self, implicit_hydrogens: UnrepresentablePolicy) -> Self {
        self.implicit_hydrogens = implicit_hydrogens;
        self
    }
}

impl Default for MmcifWriteOptions {
    fn default() -> Self {
        Self {
            block_name: "model".to_owned(),
            coordinate_precision: 3,
            ensemble_metadata: UnrepresentablePolicy::Fail,
            frame_metadata: UnrepresentablePolicy::Fail,
            cell: UnrepresentablePolicy::Fail,
            implicit_hydrogens: UnrepresentablePolicy::Fail,
        }
    }
}

/// A failure specific to the mmCIF writer: the typed
/// [`detail`](crate::diagnostics::FormatError::detail) of a write-stage
/// [`FormatError`], as
/// [`FormatErrorDetail::MmcifWrite`](crate::diagnostics::FormatErrorDetail::MmcifWrite).
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MmcifWriteIssue {
    /// [`MmcifWriteOptions::block_name`] or
    /// [`MmcifBlockSource::block_name`] is empty or not a valid data block
    /// name.
    InvalidBlockName(String),
    /// Two written data blocks would have this name, compared
    /// case-insensitively as mmCIF block names are.
    DuplicateBlockName(String),
    /// [`MmcifWriteOptions::coordinate_precision`] exceeds the supported
    /// maximum.
    CoordinatePrecisionTooLarge(usize),
    /// Explicit entity semantics leave this molecule instance unclassified.
    MissingEntityClassification(MoleculeInstanceId),
    /// Entity semantics classify this molecule instance more than once.
    DuplicateEntityClassification(MoleculeInstanceId),
    /// The instances one mmCIF entity would cover have different
    /// classifications.
    ConflictingEntityClassifications {
        /// The molecule instance.
        molecule: MoleculeInstanceId,
        /// The conflicting classifications.
        classifications: Vec<EntityKind>,
    },
    /// The classification names an entity type the writer does not support.
    UnsupportedEntityClassification {
        /// The molecule instance.
        molecule: MoleculeInstanceId,
        /// The source entity type.
        classification: String,
    },
    /// Canonical topology classification does not determine one mmCIF entity
    /// kind for this molecule instance.
    UnresolvedCanonicalEntityClassification {
        /// The molecule instance.
        molecule: MoleculeInstanceId,
        /// Its canonical class.
        classification: MoleculeClass,
    },
    /// A structural instance has several source entity IDs.
    ConflictingAsymEntityIds {
        /// The structural instance (`label_asym_id`).
        asym_id: String,
        /// The conflicting entity IDs.
        entity_ids: Vec<String>,
    },
    /// A structural instance has several entity classifications.
    ConflictingAsymEntityClassifications {
        /// The structural instance (`label_asym_id`).
        asym_id: String,
        /// The conflicting classifications.
        classifications: Vec<EntityKind>,
    },
    /// A source entity has several classifications.
    ConflictingSourceEntityClassifications {
        /// The source entity ID.
        entity_id: String,
        /// The conflicting classifications.
        classifications: Vec<EntityKind>,
    },
    /// Entity semantics name a molecule instance the model does not have.
    UnknownClassifiedMolecule(MoleculeInstanceId),
    /// Two structural instances share one ID.
    DuplicateAsymId(String),
    /// A hierarchy atom has no biomolecular atom site.
    MissingAtomSite(InstanceAtomId),
    /// An atom appears in more than one biomolecular atom site.
    DuplicateAtomSite(InstanceAtomId),
    /// An atom's atom-site field disagrees with the rest of its residue or
    /// chain.
    InconsistentAtomSite {
        /// The atom.
        atom: InstanceAtomId,
        /// The inconsistent field.
        field: MmcifAtomField,
    },
    /// An atom duplicates another atom's identity within one residue.
    DuplicateAtomIdentity(InstanceAtomId),
    /// Interpretation reports carry no atom-level provenance for this atom.
    MissingAtomProvenance(InstanceAtomId),
    /// Interpretation reports carry atom-level provenance for this atom twice.
    DuplicateAtomProvenance(InstanceAtomId),
    /// Interpretation reports carry provenance for an atom the model does not
    /// have.
    UnknownAtomProvenance(InstanceAtomId),
    /// An atom has state the writer does not represent.
    UnsupportedAtomField {
        /// The atom.
        atom: InstanceAtomId,
        /// The unsupported field.
        field: MmcifAtomField,
    },
    /// An atom's formal charge is outside the PDBx/mmCIF range -8..=8.
    FormalChargeOutOfRange {
        /// The atom.
        atom: InstanceAtomId,
        /// Its formal charge.
        charge: i8,
    },
    /// A molecule instance carries stereochemistry the writer does not
    /// represent.
    UnsupportedStereo(MoleculeInstanceId),
    /// A bond has an order mmCIF connectivity cannot represent.
    UnsupportedBondOrder {
        /// The bond.
        bond: InstanceBondId,
        /// Its order.
        order: BondOrder,
    },
    /// A `_struct_conn` partner cannot select this atom unambiguously.
    AmbiguousConnectionSelector(InstanceAtomId),
    /// A text value cannot be written as a single mmCIF token, for example
    /// because it is empty or spans lines.
    UnsupportedTextValue {
        /// The mmCIF item the value belongs to, such as
        /// `_atom_site.label_atom_id`.
        item: String,
    },
    /// The ensemble has no samples.
    EmptyEnsemble,
    /// The number of interpretation reports does not match the number of
    /// models.
    ReportCountMismatch {
        /// The number of models.
        expected: usize,
        /// The number of reports.
        actual: usize,
    },
    /// The 1-based ensemble sample `sample` differs from the first sample in
    /// state one mmCIF block shares between models.
    IncompatibleEnsembleSample {
        /// The 1-based sample number.
        sample: usize,
        /// The differing state.
        field: MmcifSampleField,
    },
    /// The 1-based ensemble sample `sample` has state mmCIF cannot
    /// represent, with [`UnrepresentablePolicy::Fail`]. A temperature is
    /// reported for the first sample.
    UnrepresentableEnsembleSample {
        /// The 1-based sample number.
        sample: usize,
        /// [`MmcifSampleField::Weight`] for unequal weights,
        /// [`MmcifSampleField::State`] for more than one state,
        /// [`MmcifSampleField::Energy`] for a state energy, or
        /// [`MmcifSampleField::Temperature`] for the ensemble temperature.
        field: MmcifSampleField,
    },
    /// A periodic cell mmCIF cannot represent, with
    /// [`UnrepresentablePolicy::Fail`].
    UnrepresentableCell {
        /// Why the cell cannot be written.
        reason: MmcifCellReason,
    },
}

/// Atom state named by an [`MmcifWriteIssue`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MmcifAtomField {
    /// The atom's isotope.
    Isotope,
    /// The atom's radical electrons.
    Radical,
    /// The atom's implicit hydrogens.
    Hydrogens,
    /// The atom's atom-map number.
    AtomMap,
    /// `_atom_site.label_asym_id`.
    LabelAsymId,
    /// The source `asym_id` recorded in provenance.
    SourceAsymId,
    /// `_atom_site.type_symbol`.
    TypeSymbol,
}

impl fmt::Display for MmcifAtomField {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Isotope => "isotope",
            Self::Radical => "radical",
            Self::Hydrogens => "hydrogens",
            Self::AtomMap => "atom_map",
            Self::LabelAsymId => "label_asym_id",
            Self::SourceAsymId => "source asym_id",
            Self::TypeSymbol => "type_symbol",
        })
    }
}

/// Ensemble sample state named by an [`MmcifWriteIssue`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MmcifSampleField {
    /// The `_entity` representation.
    Entities,
    /// The `_struct_asym` representation.
    StructuralAsymmetry,
    /// The `_struct_conn` and `_chem_comp_bond` representation.
    Connectivity,
    /// The `_atom_site` identity columns.
    AtomSites,
    /// The sample's relative weight.
    Weight,
    /// The energy of the sample's state.
    Energy,
    /// The sample's state.
    State,
    /// The temperature the ensemble's populations describe.
    Temperature,
}

impl fmt::Display for MmcifSampleField {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Entities => "entity representation",
            Self::StructuralAsymmetry => "structural asymmetry representation",
            Self::Connectivity => "connectivity representation",
            Self::AtomSites => "atom-site representation",
            Self::Weight => "weight",
            Self::Energy => "state energy",
            Self::State => "state",
            Self::Temperature => "temperature",
        })
    }
}

/// Why mmCIF cannot represent a periodic cell.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MmcifCellReason {
    /// Ensemble samples have different cells; one block has one `_cell`.
    SampleCellsDiffer,
    /// The cell is not periodic along all three axes.
    PartiallyPeriodic,
    /// The cell does not have `a` along x and `b` in the xy plane, as mmCIF
    /// coordinates assume.
    NonstandardOrientation,
}

impl fmt::Display for MmcifCellReason {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::SampleCellsDiffer => "ensemble samples have different cells",
            Self::PartiallyPeriodic => "an mmCIF cell is periodic along all three axes",
            Self::NonstandardOrientation => {
                "an mmCIF cell has a along x and b in the xy plane, as its coordinates assume"
            }
        })
    }
}

impl MmcifWriteIssue {
    /// The diagnostic kind of a write error with this detail.
    pub const fn kind(&self) -> DiagnosticKind {
        match self {
            Self::InvalidBlockName(_)
            | Self::DuplicateBlockName(_)
            | Self::CoordinatePrecisionTooLarge(_) => DiagnosticKind::InvalidValue,
            Self::UnknownClassifiedMolecule(_) | Self::UnknownAtomProvenance(_) => {
                DiagnosticKind::InvalidReference
            }
            Self::UnsupportedEntityClassification { .. }
            | Self::UnresolvedCanonicalEntityClassification { .. }
            | Self::UnsupportedAtomField { .. }
            | Self::FormalChargeOutOfRange { .. }
            | Self::UnsupportedStereo(_)
            | Self::UnsupportedBondOrder { .. }
            | Self::AmbiguousConnectionSelector(_)
            | Self::UnsupportedTextValue { .. }
            | Self::UnrepresentableEnsembleSample { .. }
            | Self::UnrepresentableCell { .. } => DiagnosticKind::Unsupported,
            Self::MissingEntityClassification(_)
            | Self::DuplicateEntityClassification(_)
            | Self::ConflictingEntityClassifications { .. }
            | Self::ConflictingAsymEntityIds { .. }
            | Self::ConflictingAsymEntityClassifications { .. }
            | Self::ConflictingSourceEntityClassifications { .. }
            | Self::DuplicateAsymId(_)
            | Self::MissingAtomSite(_)
            | Self::DuplicateAtomSite(_)
            | Self::InconsistentAtomSite { .. }
            | Self::DuplicateAtomIdentity(_)
            | Self::MissingAtomProvenance(_)
            | Self::DuplicateAtomProvenance(_)
            | Self::EmptyEnsemble
            | Self::ReportCountMismatch { .. }
            | Self::IncompatibleEnsembleSample { .. } => DiagnosticKind::Inconsistent,
        }
    }
}

impl fmt::Display for MmcifWriteIssue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBlockName(name) => {
                write!(f, "invalid mmCIF data block name `{name}`")
            }
            Self::DuplicateBlockName(name) => {
                write!(f, "duplicate mmCIF data block name `{name}`")
            }
            Self::CoordinatePrecisionTooLarge(precision) => write!(
                f,
                "mmCIF coordinate precision {precision} exceeds the supported maximum of {MAX_COORDINATE_PRECISION}"
            ),
            Self::MissingEntityClassification(molecule) => write!(
                f,
                "{molecule} has no explicit mmCIF entity classification"
            ),
            Self::DuplicateEntityClassification(molecule) => write!(
                f,
                "mmCIF entity semantics classify {molecule} more than once"
            ),
            Self::ConflictingEntityClassifications {
                molecule,
                classifications,
            } => write!(
                f,
                "{molecule} has conflicting mmCIF entity classifications {classifications:?}"
            ),
            Self::UnsupportedEntityClassification {
                molecule,
                classification,
            } => write!(
                f,
                "{molecule} has unsupported mmCIF entity classification `{classification}`"
            ),
            Self::UnresolvedCanonicalEntityClassification {
                molecule,
                classification,
            } => write!(
                f,
                "cannot derive an unambiguous mmCIF entity kind for {molecule} from canonical class {classification:?}"
            ),
            Self::ConflictingAsymEntityIds {
                asym_id,
                entity_ids,
            } => write!(
                f,
                "mmCIF structural instance `{asym_id}` has conflicting source entity IDs {entity_ids:?}"
            ),
            Self::ConflictingAsymEntityClassifications {
                asym_id,
                classifications,
            } => write!(
                f,
                "mmCIF structural instance `{asym_id}` has conflicting entity classifications {classifications:?}"
            ),
            Self::ConflictingSourceEntityClassifications {
                entity_id,
                classifications,
            } => write!(
                f,
                "source mmCIF entity `{entity_id}` has conflicting classifications {classifications:?}"
            ),
            Self::UnknownClassifiedMolecule(molecule) => write!(
                f,
                "mmCIF entity semantics reference unknown {molecule}"
            ),
            Self::DuplicateAsymId(id) => {
                write!(f, "duplicate mmCIF structural-instance ID `{id}`")
            }
            Self::MissingAtomSite(atom) => write!(f, "{atom} has no biomolecular atom site"),
            Self::DuplicateAtomSite(atom) => {
                write!(f, "{atom} appears in more than one biomolecular atom site")
            }
            Self::InconsistentAtomSite { atom, field } => {
                write!(f, "{atom} has inconsistent atom-site {field}")
            }
            Self::DuplicateAtomIdentity(atom) => write!(
                f,
                "{atom} duplicates an mmCIF atom identity within one residue"
            ),
            Self::MissingAtomProvenance(atom) => {
                write!(f, "{atom} has no atom-level mmCIF source provenance")
            }
            Self::DuplicateAtomProvenance(atom) => {
                write!(f, "{atom} has duplicate atom-level mmCIF source provenance")
            }
            Self::UnknownAtomProvenance(atom) => {
                write!(f, "mmCIF source provenance references unknown {atom}")
            }
            Self::UnsupportedAtomField { atom, field } => {
                write!(f, "{atom} has unsupported atom field `{field}`")
            }
            Self::FormalChargeOutOfRange { atom, charge } => write!(
                f,
                "{atom} formal charge {charge} is outside the PDBx/mmCIF range -8..=8"
            ),
            Self::UnsupportedStereo(molecule) => write!(
                f,
                "{molecule} contains stereochemistry not represented by the foundational mmCIF writer"
            ),
            Self::UnsupportedBondOrder { bond, order } => {
                write!(f, "{bond} has unsupported mmCIF bond order {order:?}")
            }
            Self::AmbiguousConnectionSelector(atom) => write!(
                f,
                "{atom} cannot be selected unambiguously by an mmCIF struct_conn partner"
            ),
            Self::UnsupportedTextValue { item } => write!(
                f,
                "{item} contains a text value that cannot be emitted as a single mmCIF token"
            ),
            Self::EmptyEnsemble => f.write_str("cannot write an empty ensemble as mmCIF"),
            Self::ReportCountMismatch { expected, actual } => write!(
                f,
                "mmCIF writing requires {expected} interpretation reports, but received {actual}"
            ),
            Self::IncompatibleEnsembleSample { sample, field } => write!(
                f,
                "ensemble sample {sample} has incompatible mmCIF {field}"
            ),
            Self::UnrepresentableEnsembleSample { sample, field } => write!(
                f,
                "mmCIF cannot represent the {field} of ensemble sample {sample}"
            ),
            Self::UnrepresentableCell { reason } => {
                write!(f, "mmCIF cannot represent the periodic cell: {reason}")
            }
        }
    }
}

/// Fails with [`MmcifWriteIssue::InvalidBlockName`] unless `name` is a
/// nonempty run of ASCII letters, digits, `_`, `-` and `.`.
fn validate_block_name(name: &str) -> Result<(), FormatError> {
    if name.is_empty()
        || !name
            .chars()
            .all(|character| character.is_ascii_alphanumeric() || "_-.".contains(character))
    {
        return Err(fail(MmcifWriteIssue::InvalidBlockName(name.to_owned())));
    }
    Ok(())
}

/// The write error carrying `issue` as its detail.
pub(super) fn fail(issue: MmcifWriteIssue) -> FormatError {
    FormatError::write(Format::Mmcif, issue.kind(), issue.to_string())
        .with_detail(FormatErrorDetail::MmcifWrite(issue))
}

/// The source model is internally inconsistent or incomplete.
pub(super) fn invalid_model(message: impl Into<String>) -> FormatError {
    FormatError::write(Format::Mmcif, DiagnosticKind::Inconsistent, message)
}

/// The output sink failed.
pub(super) fn sink(error: std::io::Error) -> FormatError {
    FormatError::sink(Format::Mmcif, error)
}

/// Expert mmCIF entity-kind overrides for canonical molecule instances.
///
/// Ordinary writing derives kinds from canonical topology classification.
/// Entries in this collection take precedence where exact source distinctions
/// or an expert format-specific choice is required.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MmcifEntityClassifications {
    kinds: BTreeMap<MoleculeInstanceId, EntityKind>,
}

impl MmcifEntityClassifications {
    /// An empty collection; every instance keeps its derived kind.
    pub const fn new() -> Self {
        Self {
            kinds: BTreeMap::new(),
        }
    }

    /// Assigns one explicit mmCIF entity kind to a molecule instance.
    ///
    /// # Errors
    ///
    /// A write [`FormatError`] with
    /// [`MmcifWriteIssue::DuplicateEntityClassification`] when the instance
    /// already has a kind; the collection is unchanged.
    pub fn insert(
        &mut self,
        molecule: MoleculeInstanceId,
        kind: EntityKind,
    ) -> Result<(), FormatError> {
        if self.kinds.contains_key(&molecule) {
            return Err(fail(MmcifWriteIssue::DuplicateEntityClassification(
                molecule,
            )));
        }
        self.kinds.insert(molecule, kind);
        Ok(())
    }

    /// Returns the assigned kind for one molecule instance.
    pub fn get(&self, molecule: MoleculeInstanceId) -> Option<&EntityKind> {
        self.kinds.get(&molecule)
    }

    /// Iterates over explicitly classified molecule instances in ID order.
    pub fn iter(&self) -> impl ExactSizeIterator<Item = (MoleculeInstanceId, &EntityKind)> {
        self.kinds.iter().map(|(&molecule, kind)| (molecule, kind))
    }

    /// Returns the number of explicitly classified molecule instances.
    pub fn len(&self) -> usize {
        self.kinds.len()
    }

    /// Returns whether no molecule instance has been classified.
    pub fn is_empty(&self) -> bool {
        self.kinds.is_empty()
    }
}

/// One mmCIF data block to write: a model or an ensemble, plus the entity
/// semantics used for its `_entity` and asymmetry assignments.
///
/// Convert from `&Model`, `&Ensemble`, `&StructureInterpretation`, or
/// `&StructureEnsembleInterpretation`; interpretations keep their source reports.
/// [`Self::model`] accepts any borrowed model view; an ensemble sample's
/// weight and state energy and a trajectory frame's time, step, velocities, and
/// forces fail to write unless the write options omit them.
///
/// A block is named by [`MmcifWriteOptions::block_name`] unless
/// [`Self::with_block_name`] names it.
#[derive(Debug, Clone, Copy)]
pub struct MmcifBlockSource<'a> {
    content: BlockContent<'a>,
    entities: EntitySemantics<'a>,
    block_name: Option<&'a str>,
}

#[derive(Debug, Clone, Copy)]
enum BlockContent<'a> {
    Model(ModelView<'a>),
    Ensemble(&'a Ensemble),
}

#[derive(Debug, Clone, Copy)]
enum EntitySemantics<'a> {
    Canonical,
    Classifications(&'a MmcifEntityClassifications),
    Reports(&'a [StructureInterpretationReport]),
}

impl<'a> MmcifBlockSource<'a> {
    /// One model block with entity kinds derived from canonical topology
    /// classification.
    pub fn model(model: &'a (impl AsModelView + ?Sized)) -> Self {
        Self {
            content: BlockContent::Model(model.as_model_view()),
            entities: EntitySemantics::Canonical,
            block_name: None,
        }
    }

    /// One multi-model block with one model per ensemble sample.
    pub fn ensemble(ensemble: &'a Ensemble) -> Self {
        Self {
            content: BlockContent::Ensemble(ensemble),
            entities: EntitySemantics::Canonical,
            block_name: None,
        }
    }

    /// Overrides derived entity kinds for the listed molecule instances;
    /// omitted instances keep canonical classification.
    ///
    /// One mmCIF entity is assigned to each populated hierarchy chain (and one
    /// to each hierarchy-free instance), so instances touched by the same chain
    /// must have the same classification.
    #[must_use]
    pub fn with_classifications(mut self, classifications: &'a MmcifEntityClassifications) -> Self {
        self.entities = EntitySemantics::Classifications(classifications);
        self
    }

    /// Preserves interpreted mmCIF entity and asymmetry semantics: one report
    /// for a model block, one per sample for an ensemble block.
    ///
    /// Atom-level provenance keeps one source entity and structural asymmetry
    /// consistent even when it spans several molecule instances; conflicting
    /// source identity is rejected. An auth-only source is normalized by
    /// copying author identifiers into the required label fields.
    #[must_use]
    pub fn with_reports(mut self, reports: &'a [StructureInterpretationReport]) -> Self {
        self.entities = EntitySemantics::Reports(reports);
        self
    }

    /// The name of this block, after `data_`, or `None` when the write
    /// options name it.
    pub const fn block_name(&self) -> Option<&'a str> {
        self.block_name
    }

    /// Names this block `block_name`, overriding
    /// [`MmcifWriteOptions::block_name`] and its numbering. Reading the file
    /// back gives the name as [`MmcifBlock::name`](crate::mmcif::MmcifBlock::name).
    ///
    /// Blocks of one write must have different names, compared
    /// case-insensitively, including the names the options give unnamed
    /// blocks; a clash fails the write with
    /// [`MmcifWriteIssue::DuplicateBlockName`] before anything is written.
    ///
    /// # Errors
    ///
    /// A write [`FormatError`] with [`MmcifWriteIssue::InvalidBlockName`]
    /// unless the name is a nonempty run of ASCII letters, digits, `_`, `-`
    /// and `.`.
    pub fn with_block_name(mut self, block_name: &'a str) -> Result<Self, FormatError> {
        validate_block_name(block_name)?;
        self.block_name = Some(block_name);
        Ok(self)
    }
}

impl<'a> From<&'a Model> for MmcifBlockSource<'a> {
    fn from(model: &'a Model) -> Self {
        Self::model(model)
    }
}

impl<'a> From<&'a Ensemble> for MmcifBlockSource<'a> {
    fn from(ensemble: &'a Ensemble) -> Self {
        Self::ensemble(ensemble)
    }
}

impl<'a> From<&'a StructureInterpretation> for MmcifBlockSource<'a> {
    fn from(interpretation: &'a StructureInterpretation) -> Self {
        Self::model(interpretation.model())
            .with_reports(std::slice::from_ref(interpretation.report()))
    }
}

impl<'a> From<&'a StructureEnsembleInterpretation> for MmcifBlockSource<'a> {
    fn from(interpretation: &'a StructureEnsembleInterpretation) -> Self {
        Self::ensemble(interpretation.ensemble()).with_reports(interpretation.reports())
    }
}

/// Writes data blocks to a string. See [`write_mmcif_to`].
pub(crate) fn write_mmcif<'a, B: Into<MmcifBlockSource<'a>>>(
    blocks: impl IntoIterator<Item = B>,
    options: MmcifWriteOptions,
) -> Result<(String, WriteReport), FormatError> {
    let mut output = Vec::new();
    let report = write_mmcif_to(&mut output, blocks, options)?;
    Ok((
        String::from_utf8(output).expect("mmCIF writer emits UTF-8"),
        report,
    ))
}

/// Writes one deterministic data block per source in input order.
///
/// A block with its own [`MmcifBlockSource::block_name`] takes it. Otherwise a
/// single block is named [`MmcifWriteOptions::block_name`], and among several
/// blocks the one at 1-based position `n` is named `{block_name}_n`. Rows
/// follow dense atom order.
pub(crate) fn write_mmcif_to<'a, B: Into<MmcifBlockSource<'a>>>(
    writer: &mut impl Write,
    blocks: impl IntoIterator<Item = B>,
    options: MmcifWriteOptions,
) -> Result<WriteReport, FormatError> {
    let mut omissions = Omissions::default();
    let blocks = blocks.into_iter().map(Into::into).collect::<Vec<_>>();
    let numbered = blocks.len() > 1;
    let names = blocks
        .iter()
        .enumerate()
        .map(|(index, block)| match block.block_name {
            Some(name) => name.to_owned(),
            None if numbered => format!("{}_{}", options.block_name, index + 1),
            None => options.block_name.clone(),
        })
        .collect::<Vec<_>>();
    let mut seen = BTreeSet::new();
    if let Some(duplicate) = names
        .iter()
        .find(|name| !seen.insert(name.to_ascii_lowercase()))
    {
        return Err(fail(MmcifWriteIssue::DuplicateBlockName(duplicate.clone())));
    }
    for (block, block_name) in blocks.into_iter().zip(names) {
        let block_options = MmcifWriteOptions {
            block_name,
            ..options.clone()
        };
        match block.content {
            BlockContent::Model(view) => {
                write_model_block_to(writer, view, block.entities, &block_options, &mut omissions)?
            }
            BlockContent::Ensemble(ensemble) => write_ensemble_block_to(
                writer,
                ensemble,
                block.entities,
                &block_options,
                &mut omissions,
            )?,
        }
    }
    Ok(omissions.into_report())
}

fn write_model_block_to(
    writer: &mut impl Write,
    view: ModelView<'_>,
    entities: EntitySemantics<'_>,
    options: &MmcifWriteOptions,
    omissions: &mut Omissions,
) -> Result<(), FormatError> {
    check_structure_state(
        Format::Mmcif,
        view.kind(),
        Some(options.ensemble_metadata),
        options.frame_metadata,
    )?;
    let plan = match entities {
        EntitySemantics::Reports(reports) => {
            let [report] = reports else {
                return Err(fail(MmcifWriteIssue::ReportCountMismatch {
                    expected: 1,
                    actual: reports.len(),
                }));
            };
            entity_plan_from_report(view, report)?
        }
        EntitySemantics::Classifications(classifications) => {
            let normalized = normalize_entity_classifications(
                view,
                classifications
                    .iter()
                    .map(|(molecule, kind)| (molecule, vec![kind.clone()])),
            )?;
            generic_entity_plan(view, &normalized)?
        }
        EntitySemantics::Canonical => {
            let normalized = normalize_entity_classifications(view, std::iter::empty())?;
            generic_entity_plan(view, &normalized)?
        }
    };
    let cell = block_cell([view.cell()], options.cell)?;
    let prepared = prepare_model(view, plan, options.implicit_hydrogens)?;
    render_model_to(writer, &prepared, cell, options)?;
    // Nothing is reported for a block that did not reach the sink.
    omissions.record_structure_state(
        view.kind(),
        Some(options.ensemble_metadata),
        options.frame_metadata,
    );
    if cell.is_none() && view.cell().is_some() {
        omissions.record(WriteOmission::Cell { structures: 1 });
    }
    omissions.record(WriteOmission::ImplicitHydrogens {
        atoms: prepared.omitted_hydrogens,
    });
    Ok(())
}

/// The `_cell` parameters of a block whose models have these cells, or none
/// when they have none or the caller omits a cell mmCIF cannot represent.
fn block_cell<'a>(
    cells: impl IntoIterator<Item = Option<&'a PeriodicCell>>,
    policy: UnrepresentablePolicy,
) -> Result<Option<CellParameters>, FormatError> {
    let mut cells = cells.into_iter();
    let first = cells.next().flatten();
    let parameters = if cells.any(|cell| cell != first) {
        Err(MmcifCellReason::SampleCellsDiffer)
    } else {
        match first.map(cell_parameters).transpose() {
            Ok(parameters) => return Ok(parameters),
            Err(CellParameterError::PartiallyPeriodic) => Err(MmcifCellReason::PartiallyPeriodic),
            Err(CellParameterError::NonstandardOrientation) => {
                Err(MmcifCellReason::NonstandardOrientation)
            }
        }
    };
    match (parameters, policy) {
        (Err(reason), UnrepresentablePolicy::Fail) => {
            Err(fail(MmcifWriteIssue::UnrepresentableCell { reason }))
        }
        (Err(_), UnrepresentablePolicy::Omit) => Ok(None),
        (Ok(parameters), _) => Ok(parameters),
    }
}

fn write_ensemble_block_to(
    writer: &mut impl Write,
    ensemble: &Ensemble,
    entities: EntitySemantics<'_>,
    options: &MmcifWriteOptions,
    omissions: &mut Omissions,
) -> Result<(), FormatError> {
    let (classifications, reports) = match entities {
        EntitySemantics::Canonical => (None, None),
        EntitySemantics::Classifications(classifications) => (Some(classifications), None),
        EntitySemantics::Reports(reports) => {
            if reports.len() != ensemble.len() {
                return Err(fail(MmcifWriteIssue::ReportCountMismatch {
                    expected: ensemble.len(),
                    actual: reports.len(),
                }));
            }
            (None, Some(reports))
        }
    };
    // Omissions are recorded only once the block has reached the sink.
    let mut ensemble_omissions = Omissions::default();
    check_multi_model_ensemble(ensemble, options.ensemble_metadata, &mut ensemble_omissions)
        .map_err(|(sample, field)| {
            let field = match field {
                StructureField::Weight => MmcifSampleField::Weight,
                StructureField::Energy => MmcifSampleField::Energy,
                StructureField::State => MmcifSampleField::State,
                _ => MmcifSampleField::Temperature,
            };
            fail(MmcifWriteIssue::UnrepresentableEnsembleSample { sample, field })
        })?;
    let cell = block_cell(
        ensemble.iter().map(|sample| sample.as_model_view().cell()),
        options.cell,
    )?;
    let mut samples = ensemble.iter().enumerate();
    let (_, first_sample) = samples
        .next()
        .ok_or_else(|| fail(MmcifWriteIssue::EmptyEnsemble))?;
    let first_view = first_sample.as_model_view();
    let first_plan = ensemble_entity_plan(first_view, classifications, reports.map(|all| &all[0]))?;
    let first = prepare_model(first_view, first_plan, options.implicit_hydrogens)?;

    write_block_start(writer, &first, cell, options)?;
    let mut atom_serial = 1u64;
    write_atom_rows(writer, &first.atoms, 1, &mut atom_serial, options)?;

    for (index, sample) in samples {
        let view = sample.as_model_view();
        let plan = ensemble_entity_plan(view, classifications, reports.map(|all| &all[index]))?;
        let candidate = prepare_model(view, plan, options.implicit_hydrogens)?;
        validate_ensemble_sample(&first, &candidate, index + 1)?;
        write_atom_rows(
            writer,
            &candidate.atoms,
            index + 1,
            &mut atom_serial,
            options,
        )?;
    }
    write_block_end(writer, &first)?;
    for omission in ensemble_omissions.into_report().omissions() {
        omissions.record(*omission);
    }
    if cell.is_none() {
        omissions.record(WriteOmission::Cell {
            structures: ensemble
                .iter()
                .filter(|sample| sample.as_model_view().cell().is_some())
                .count(),
        });
    }
    // Samples share one topology, so the block omits each atom's count once.
    omissions.record(WriteOmission::ImplicitHydrogens {
        atoms: first.omitted_hydrogens,
    });
    Ok(())
}

fn one_based_serial(raw: u32) -> u64 {
    u64::from(raw) + 1
}
