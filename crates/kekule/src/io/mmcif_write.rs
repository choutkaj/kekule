//! Public writer API and orchestration of planning, preparation, and emission.
use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail};
use crate::geometry::PeriodicCell;
use crate::io::mmcif_cell::{cell_parameters, CellParameterError, CellParameters};
use crate::io::structure_interpret::{
    EntityKind, StructureEnsembleInterpretation, StructureInterpretation,
    StructureInterpretationReport,
};
use crate::io::write_error::check_structure_state;
use crate::io::UnrepresentablePolicy;
use crate::molecule::BondOrder;
use crate::structure::{AsModelView, Ensemble, Model, ModelView};
use crate::topology::{InstanceAtomId, InstanceBondId, MoleculeClass, MoleculeInstanceId};
use std::collections::BTreeMap;
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
use prepare::{prepare_model, validate_ensemble_member};

const MAX_COORDINATE_PRECISION: usize = 15;

/// Options of the mmCIF writer.
///
/// The default writes a block named `model` with coordinates to 3 decimal
/// places and fails on ensemble metadata or a periodic cell mmCIF cannot
/// carry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifWriteOptions {
    block_name: String,
    coordinate_precision: usize,
    ensemble_metadata: UnrepresentablePolicy,
    frame_metadata: UnrepresentablePolicy,
    cell: UnrepresentablePolicy,
}

impl MmcifWriteOptions {
    /// The name of the written data block, after `data_`; `model` by default.
    /// Several sources are written as blocks `{name}_1`, `{name}_2`, and so on.
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
        if block_name.is_empty()
            || !block_name
                .chars()
                .all(|character| character.is_ascii_alphanumeric() || "_-.".contains(character))
        {
            return Err(fail(MmcifWriteIssue::InvalidBlockName(block_name)));
        }
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

    /// Ensemble member weights and energies, which mmCIF cannot carry: by
    /// default an ensemble with unequal weights or any energy fails to write,
    /// and so does a model block of one ensemble member, whose weight is
    /// relative to members the block does not hold.
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
    /// between ensemble members. By default such a cell fails to write;
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
}

impl Default for MmcifWriteOptions {
    fn default() -> Self {
        Self {
            block_name: "model".to_owned(),
            coordinate_precision: 3,
            ensemble_metadata: UnrepresentablePolicy::Fail,
            frame_metadata: UnrepresentablePolicy::Fail,
            cell: UnrepresentablePolicy::Fail,
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
    /// [`MmcifWriteOptions::block_name`] is empty or not a valid data block
    /// name.
    InvalidBlockName(String),
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
    /// The ensemble has no members.
    EmptyEnsemble,
    /// The number of interpretation reports does not match the number of
    /// models.
    ReportCountMismatch {
        /// The number of models.
        expected: usize,
        /// The number of reports.
        actual: usize,
    },
    /// The 1-based ensemble member `member` differs from the first member in
    /// state one mmCIF block shares between models.
    IncompatibleEnsembleMember {
        /// The 1-based member number.
        member: usize,
        /// The differing state.
        field: MmcifMemberField,
    },
    /// The 1-based ensemble member `member` has state mmCIF cannot represent,
    /// with [`UnrepresentablePolicy::Fail`].
    UnrepresentableEnsembleMember {
        /// The 1-based member number.
        member: usize,
        /// [`MmcifMemberField::Weight`] for unequal weights or
        /// [`MmcifMemberField::Energy`] for a member energy.
        field: MmcifMemberField,
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

/// Ensemble member state named by an [`MmcifWriteIssue`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MmcifMemberField {
    /// The `_entity` representation.
    Entities,
    /// The `_struct_asym` representation.
    StructuralAsymmetry,
    /// The `_struct_conn` and `_chem_comp_bond` representation.
    Connectivity,
    /// The `_atom_site` identity columns.
    AtomSites,
    /// The member's relative weight.
    Weight,
    /// The member's energy.
    Energy,
}

impl fmt::Display for MmcifMemberField {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Entities => "entity representation",
            Self::StructuralAsymmetry => "structural asymmetry representation",
            Self::Connectivity => "connectivity representation",
            Self::AtomSites => "atom-site representation",
            Self::Weight => "weight",
            Self::Energy => "energy",
        })
    }
}

/// Why mmCIF cannot represent a periodic cell.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MmcifCellReason {
    /// Ensemble members have different cells; one block has one `_cell`.
    MemberCellsDiffer,
    /// The cell is not periodic along all three axes.
    PartiallyPeriodic,
    /// The cell does not have `a` along x and `b` in the xy plane, as mmCIF
    /// coordinates assume.
    NonstandardOrientation,
}

impl fmt::Display for MmcifCellReason {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::MemberCellsDiffer => "ensemble members have different cells",
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
            Self::InvalidBlockName(_) | Self::CoordinatePrecisionTooLarge(_) => {
                DiagnosticKind::InvalidValue
            }
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
            | Self::UnrepresentableEnsembleMember { .. }
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
            | Self::IncompatibleEnsembleMember { .. } => DiagnosticKind::Inconsistent,
        }
    }
}

impl fmt::Display for MmcifWriteIssue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBlockName(name) => {
                write!(f, "invalid mmCIF data block name `{name}`")
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
            Self::IncompatibleEnsembleMember { member, field } => write!(
                f,
                "ensemble member {member} has incompatible mmCIF {field}"
            ),
            Self::UnrepresentableEnsembleMember { member, field } => write!(
                f,
                "mmCIF cannot represent the {field} of ensemble member {member}"
            ),
            Self::UnrepresentableCell { reason } => {
                write!(f, "mmCIF cannot represent the periodic cell: {reason}")
            }
        }
    }
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
/// [`Self::model`] accepts any borrowed model view; an ensemble member's
/// weight and energy and a trajectory frame's time, step, velocities, and
/// forces fail to write unless the write options omit them.
#[derive(Debug, Clone, Copy)]
pub struct MmcifBlockSource<'a> {
    content: BlockContent<'a>,
    entities: EntitySemantics<'a>,
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
        }
    }

    /// One multi-model block with one model per ensemble member.
    pub fn ensemble(ensemble: &'a Ensemble) -> Self {
        Self {
            content: BlockContent::Ensemble(ensemble),
            entities: EntitySemantics::Canonical,
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
    /// for a model block, one per member for an ensemble block.
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
) -> Result<String, FormatError> {
    let mut output = Vec::new();
    write_mmcif_to(&mut output, blocks, options)?;
    Ok(String::from_utf8(output).expect("mmCIF writer emits UTF-8"))
}

/// Writes one deterministic data block per source in input order.
///
/// A single block is named [`MmcifWriteOptions::block_name`]; several blocks
/// are named `{block_name}_1`, `{block_name}_2`, and so on. Rows follow
/// dense atom order.
pub(crate) fn write_mmcif_to<'a, B: Into<MmcifBlockSource<'a>>>(
    writer: &mut impl Write,
    blocks: impl IntoIterator<Item = B>,
    options: MmcifWriteOptions,
) -> Result<(), FormatError> {
    let blocks = blocks.into_iter().map(Into::into).collect::<Vec<_>>();
    let numbered = blocks.len() > 1;
    for (index, block) in blocks.into_iter().enumerate() {
        let block_options = MmcifWriteOptions {
            block_name: if numbered {
                format!("{}_{}", options.block_name, index + 1)
            } else {
                options.block_name.clone()
            },
            ..options.clone()
        };
        match block.content {
            BlockContent::Model(view) => {
                write_model_block_to(writer, view, block.entities, &block_options)?
            }
            BlockContent::Ensemble(ensemble) => {
                write_ensemble_block_to(writer, ensemble, block.entities, &block_options)?
            }
        }
    }
    Ok(())
}

fn write_model_block_to(
    writer: &mut impl Write,
    view: ModelView<'_>,
    entities: EntitySemantics<'_>,
    options: &MmcifWriteOptions,
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
    let prepared = prepare_model(view, plan)?;
    render_model_to(writer, &prepared, cell, options)
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
        Err(MmcifCellReason::MemberCellsDiffer)
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
    if options.ensemble_metadata == UnrepresentablePolicy::Fail {
        // Interpretation gives every member weight 1.0 and no energy, so
        // only equal weights and absent energies survive a round trip.
        let weight = ensemble.get(0).map(|member| member.weight());
        for (index, member) in ensemble.iter().enumerate() {
            let field = if member.energy().is_some() {
                MmcifMemberField::Energy
            } else if Some(member.weight()) != weight {
                MmcifMemberField::Weight
            } else {
                continue;
            };
            return Err(fail(MmcifWriteIssue::UnrepresentableEnsembleMember {
                member: index + 1,
                field,
            }));
        }
    }
    let cell = block_cell(
        ensemble.iter().map(|member| member.as_model_view().cell()),
        options.cell,
    )?;
    let mut members = ensemble.iter().enumerate();
    let (_, first_member) = members
        .next()
        .ok_or_else(|| fail(MmcifWriteIssue::EmptyEnsemble))?;
    let first_view = first_member.as_model_view();
    let first_plan = ensemble_entity_plan(first_view, classifications, reports.map(|all| &all[0]))?;
    let first = prepare_model(first_view, first_plan)?;

    write_block_start(writer, &first, cell, options)?;
    let mut atom_serial = 1u64;
    write_atom_rows(writer, &first.atoms, 1, &mut atom_serial, options)?;

    for (index, member) in members {
        let view = member.as_model_view();
        let plan = ensemble_entity_plan(view, classifications, reports.map(|all| &all[index]))?;
        let candidate = prepare_model(view, plan)?;
        validate_ensemble_member(&first, &candidate, index + 1)?;
        write_atom_rows(
            writer,
            &candidate.atoms,
            index + 1,
            &mut atom_serial,
            options,
        )?;
    }
    write_block_end(writer, &first)
}

fn one_based_serial(raw: u32) -> u64 {
    u64::from(raw) + 1
}
