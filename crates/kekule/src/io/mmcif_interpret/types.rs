use std::fmt;

use crate::structure::{Ensemble, EnsembleError, Model};
use crate::topology::{InstanceAtomId, MoleculeInstanceId};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MmcifAltLocPolicy {
    HighestOccupancy,
    SelectLabel(String),
    ErrorOnAlternateLocations,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MmcifModelSelection {
    RequireSingle,
    Select(String),
    First,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifInterpretOptions {
    pub strict_entity_metadata: bool,
    pub altloc_policy: MmcifAltLocPolicy,
    pub model_selection: MmcifModelSelection,
}

impl Default for MmcifInterpretOptions {
    fn default() -> Self {
        Self {
            strict_entity_metadata: false,
            altloc_policy: MmcifAltLocPolicy::HighestOccupancy,
            model_selection: MmcifModelSelection::RequireSingle,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum MmcifEntityKind {
    Polymer,
    Branched,
    NonPolymer,
    Water,
    Other(String),
}

impl MmcifEntityKind {
    pub(super) fn from_mmcif(value: &str) -> Self {
        match value.to_ascii_lowercase().as_str() {
            "polymer" => Self::Polymer,
            "branched" => Self::Branched,
            "non-polymer" => Self::NonPolymer,
            "water" => Self::Water,
            _ => Self::Other(value.to_owned()),
        }
    }

    pub(super) fn is_macro(&self) -> bool {
        matches!(self, Self::Polymer | Self::Branched)
    }
}

/// Non-fatal scientific interpretation issues retained in the mmCIF report.
///
/// Connection issues preserve source identity and distinguish unresolved
/// selectors from ambiguous selectors. Neither case creates a topology bond.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MmcifInterpretIssue {
    EntityTypeInferred {
        asym_id: String,
        kind: MmcifEntityKind,
    },
    ConnectionIgnored {
        connection_type: String,
    },
    /// One declared connection partner matched no selected-model atom.
    ConnectionUnresolved {
        connection_id: Option<String>,
        connection_type: String,
        partner: u8,
        source_line: Option<usize>,
        reason: MmcifConnectionResolutionReason,
    },
    /// One declared connection partner matched more than one selected-model atom.
    ConnectionAmbiguous {
        connection_id: Option<String>,
        connection_type: String,
        partner: u8,
        source_line: Option<usize>,
        candidates: usize,
        reason: MmcifConnectionResolutionReason,
    },
    CoordinateModelIgnored {
        model_id: String,
        atom_site_rows: usize,
    },
    AlternateLocationOmitted {
        atom_name: String,
        alt_id: Option<String>,
    },
    ConnectivityCandidatesInferred {
        atom_count: usize,
        candidate_count: usize,
    },
}

/// Why a `_struct_conn` partner could not be resolved uniquely.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MmcifConnectionResolutionReason {
    /// The record requests a symmetry mate outside the asymmetric unit.
    UnsupportedSymmetry { symmetry: String },
    /// No label or author atom-identity selector was supplied.
    MissingSelector,
    /// Alias fields for the same selector supplied conflicting values.
    ConflictingSelectorValues { selector: &'static str },
    /// Label and author selector families matched disjoint atom sets.
    ConflictingLabelAndAuthorSelectors,
    /// An explicitly named alternate location existed before selection but was omitted.
    AlternateLocationOmitted { alternate_location: String },
    /// No selected-model atom satisfied every supplied selector field.
    NoMatchingAtom,
    /// More than one selected-model atom satisfied every supplied selector field.
    MultipleMatchingAtoms,
}

impl fmt::Display for MmcifConnectionResolutionReason {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedSymmetry { symmetry } => {
                write!(formatter, "unsupported symmetry mate `{symmetry}`")
            }
            Self::MissingSelector => formatter.write_str("missing partner atom selector"),
            Self::ConflictingSelectorValues { selector } => {
                write!(formatter, "conflicting values for {selector}")
            }
            Self::ConflictingLabelAndAuthorSelectors => {
                formatter.write_str("label and author selectors identify different atoms")
            }
            Self::AlternateLocationOmitted { alternate_location } => write!(
                formatter,
                "alternate location `{alternate_location}` was omitted by selection policy"
            ),
            Self::NoMatchingAtom => formatter.write_str("no atom matches the partner selector"),
            Self::MultipleMatchingAtoms => {
                formatter.write_str("multiple atoms match the partner selector")
            }
        }
    }
}

/// Stable source correspondence for one interpreted mmCIF atom.
///
/// Residue, occurrence, asymmetry, and selected alternate-location fields are
/// retained independently of derived molecule-instance and dense ordering.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifAtomProvenance {
    pub(crate) atom: InstanceAtomId,
    pub(crate) type_symbol: String,
    pub(crate) source_line: usize,
    pub(crate) atom_site_id: Option<String>,
    pub(crate) atom_name: String,
    pub(crate) auth_atom_name: Option<String>,
    pub(crate) component_id: String,
    pub(crate) auth_component_id: Option<String>,
    pub(crate) asym_id: String,
    pub(crate) auth_asym_id: Option<String>,
    pub(crate) entity_id: Option<String>,
    pub(crate) label_sequence_id: Option<i32>,
    pub(crate) author_sequence_id: Option<String>,
    pub(crate) insertion_code: Option<String>,
    pub(crate) occurrence: Option<usize>,
    pub(crate) selected_alternate_location: Option<String>,
}

impl MmcifAtomProvenance {
    pub const fn atom(&self) -> InstanceAtomId {
        self.atom
    }

    pub const fn source_line(&self) -> usize {
        self.source_line
    }

    pub fn atom_site_id(&self) -> Option<&str> {
        self.atom_site_id.as_deref()
    }

    pub fn atom_name(&self) -> &str {
        &self.atom_name
    }

    pub fn auth_atom_name(&self) -> Option<&str> {
        self.auth_atom_name.as_deref()
    }

    pub fn component_id(&self) -> &str {
        &self.component_id
    }

    pub fn auth_component_id(&self) -> Option<&str> {
        self.auth_component_id.as_deref()
    }

    /// Returns `_atom_site.label_asym_id`, falling back to
    /// `_atom_site.auth_asym_id` when the label identifier is absent.
    pub fn asym_id(&self) -> &str {
        &self.asym_id
    }

    /// Returns the independently preserved author asymmetry identifier.
    pub fn auth_asym_id(&self) -> Option<&str> {
        self.auth_asym_id.as_deref()
    }

    pub fn entity_id(&self) -> Option<&str> {
        self.entity_id.as_deref()
    }

    pub const fn label_sequence_id(&self) -> Option<i32> {
        self.label_sequence_id
    }

    pub fn author_sequence_id(&self) -> Option<&str> {
        self.author_sequence_id.as_deref()
    }

    pub fn insertion_code(&self) -> Option<&str> {
        self.insertion_code.as_deref()
    }

    /// Returns the zero-based source occurrence discriminator used when both
    /// label and author sequence identifiers are absent.
    ///
    /// The discriminator is scoped to one coordinate model, asymmetry
    /// identifier, and component identifier.
    pub const fn occurrence(&self) -> Option<usize> {
        self.occurrence
    }

    /// Returns the alternate-location label selected for this source atom.
    ///
    /// This is coordinate-row provenance and is not topology atom identity.
    pub fn selected_alternate_location(&self) -> Option<&str> {
        self.selected_alternate_location.as_deref()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifInstanceProvenance {
    pub(crate) molecule: MoleculeInstanceId,
    pub(crate) coordinate_model_id: String,
    pub(crate) asym_ids: Vec<String>,
    pub(crate) entity_ids: Vec<String>,
    pub(crate) entity_kinds: Vec<MmcifEntityKind>,
    pub(crate) atoms: Vec<MmcifAtomProvenance>,
}

impl MmcifInstanceProvenance {
    pub const fn molecule(&self) -> MoleculeInstanceId {
        self.molecule
    }

    pub fn coordinate_model_id(&self) -> &str {
        &self.coordinate_model_id
    }

    pub fn asym_ids(&self) -> &[String] {
        &self.asym_ids
    }

    pub fn entity_ids(&self) -> &[String] {
        &self.entity_ids
    }

    pub fn entity_kinds(&self) -> &[MmcifEntityKind] {
        &self.entity_kinds
    }

    pub fn atoms(&self) -> &[MmcifAtomProvenance] {
        &self.atoms
    }
}

/// Structured record of one mmCIF interpretation.
///
/// Applied connection counts include only uniquely resolved local covalent
/// records. [`Self::issues`] retains source-aware unresolved and ambiguous
/// partner diagnostics.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MmcifInterpretationReport {
    pub(crate) data_block: String,
    pub(crate) entity_definitions: usize,
    pub(crate) coordinate_models: usize,
    pub(crate) selected_model: Option<String>,
    pub(crate) ignored_coordinate_models: Vec<String>,
    pub(crate) macromolecules: usize,
    pub(crate) small_molecules: usize,
    pub(crate) solvent_molecules: usize,
    pub(crate) applied_connections: usize,
    pub(crate) connectivity_candidates: usize,
    pub(crate) instances: Vec<MmcifInstanceProvenance>,
    pub(crate) issues: Vec<MmcifInterpretIssue>,
}

impl MmcifInterpretationReport {
    pub fn data_block(&self) -> &str {
        &self.data_block
    }

    pub const fn entity_definitions(&self) -> usize {
        self.entity_definitions
    }

    pub const fn coordinate_models(&self) -> usize {
        self.coordinate_models
    }

    pub fn selected_model(&self) -> Option<&str> {
        self.selected_model.as_deref()
    }

    pub fn ignored_coordinate_models(&self) -> &[String] {
        &self.ignored_coordinate_models
    }

    pub const fn macromolecules(&self) -> usize {
        self.macromolecules
    }

    pub const fn small_molecules(&self) -> usize {
        self.small_molecules
    }

    pub const fn solvent_molecules(&self) -> usize {
        self.solvent_molecules
    }

    /// Returns the number of uniquely resolved covalent connection records.
    pub const fn applied_connections(&self) -> usize {
        self.applied_connections
    }

    pub const fn connectivity_candidates(&self) -> usize {
        self.connectivity_candidates
    }

    pub fn instances(&self) -> &[MmcifInstanceProvenance] {
        &self.instances
    }

    /// Returns source-aware non-fatal interpretation issues.
    ///
    /// Unresolved and ambiguous `_struct_conn` partners are reported here and
    /// never create molecule merges or topology bonds.
    pub fn issues(&self) -> &[MmcifInterpretIssue] {
        &self.issues
    }
}

#[derive(Debug, Clone, PartialEq)]
/// One final canonical model and its format-specific mmCIF interpretation report.
pub struct MmcifInterpretation {
    pub(super) model: Model,
    pub(super) report: MmcifInterpretationReport,
}

impl MmcifInterpretation {
    pub fn model(&self) -> &Model {
        &self.model
    }

    pub fn report(&self) -> &MmcifInterpretationReport {
        &self.report
    }

    pub fn topology(&self) -> &crate::topology::Topology {
        self.model.topology()
    }

    pub fn to_model(self) -> Model {
        self.model
    }

    pub fn to_parts(self) -> (Model, MmcifInterpretationReport) {
        (self.model, self.report)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifInterpretError {
    pub(crate) line: Option<usize>,
    pub(crate) message: String,
}

impl MmcifInterpretError {
    pub(super) fn new(line: Option<usize>, message: impl Into<String>) -> Self {
        Self {
            line,
            message: message.into(),
        }
    }

    pub const fn line(&self) -> Option<usize> {
        self.line
    }

    pub fn message(&self) -> &str {
        &self.message
    }
}

impl fmt::Display for MmcifInterpretError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.line {
            Some(line) => write!(
                f,
                "mmCIF interpretation error at line {line}: {}",
                self.message
            ),
            None => write!(f, "mmCIF interpretation error: {}", self.message),
        }
    }
}

impl std::error::Error for MmcifInterpretError {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifEnsembleInterpretOptions {
    pub strict_entity_metadata: bool,
    pub altloc_policy: MmcifAltLocPolicy,
    /// `None` selects all coordinate models in source order. `Some` must
    /// contain at least one model ID.
    pub model_ids: Option<Vec<String>>,
}

impl Default for MmcifEnsembleInterpretOptions {
    fn default() -> Self {
        Self {
            strict_entity_metadata: false,
            altloc_policy: MmcifAltLocPolicy::HighestOccupancy,
            model_ids: None,
        }
    }
}

#[derive(Debug, Clone)]
/// One shared-topology ensemble and one mmCIF report per coordinate model.
pub struct MmcifEnsembleInterpretation {
    pub(super) ensemble: Ensemble,
    pub(super) reports: Vec<MmcifInterpretationReport>,
}

impl MmcifEnsembleInterpretation {
    pub fn ensemble(&self) -> &Ensemble {
        &self.ensemble
    }

    pub fn reports(&self) -> &[MmcifInterpretationReport] {
        &self.reports
    }

    pub fn to_parts(self) -> (Ensemble, Vec<MmcifInterpretationReport>) {
        (self.ensemble, self.reports)
    }
}

#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum MmcifEnsembleInterpretError {
    NoCoordinateModels,
    EmptyModelSelection,
    MultipleAtomSiteDataBlocks,
    DuplicateRequestedModel(String),
    UnknownRequestedModel(String),
    Model {
        model_id: String,
        error: MmcifInterpretError,
    },
    InconsistentTopology {
        model_id: String,
    },
    InconsistentAtomSet {
        model_id: String,
    },
    InconsistentDenseAtomOrder {
        model_id: String,
    },
    Ensemble(Box<EnsembleError>),
    Position(crate::structure::PositionError),
}

impl fmt::Display for MmcifEnsembleInterpretError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoCoordinateModels => {
                formatter.write_str("mmCIF document contains no coordinate models")
            }
            Self::EmptyModelSelection => {
                formatter.write_str("explicit mmCIF coordinate model selection is empty")
            }
            Self::MultipleAtomSiteDataBlocks => formatter
                .write_str("mmCIF document has atom-site data in more than one data block"),
            Self::DuplicateRequestedModel(model) => {
                write!(formatter, "coordinate model `{model}` was requested more than once")
            }
            Self::UnknownRequestedModel(model) => {
                write!(formatter, "coordinate model `{model}` is not present")
            }
            Self::Model { model_id, error } => {
                write!(formatter, "cannot interpret coordinate model `{model_id}`: {error}")
            }
            Self::InconsistentTopology { model_id } => write!(
                formatter,
                "coordinate model `{model_id}` has inconsistent molecule partition, chemistry, connectivity, or hierarchy"
            ),
            Self::InconsistentAtomSet { model_id } => write!(
                formatter,
                "coordinate model `{model_id}` has an inconsistent atom identity set"
            ),
            Self::InconsistentDenseAtomOrder { model_id } => write!(
                formatter,
                "coordinate model `{model_id}` has an inconsistent dense atom identity order"
            ),
            Self::Ensemble(error) => write!(formatter, "cannot assemble mmCIF ensemble: {error}"),
            Self::Position(error) => {
                write!(formatter, "cannot transfer mmCIF member positions: {error}")
            }
        }
    }
}

impl std::error::Error for MmcifEnsembleInterpretError {}
