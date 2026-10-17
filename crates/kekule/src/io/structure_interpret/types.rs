use std::fmt;
use std::sync::Arc;

use super::altloc::{AltLocDecision, AltLocPolicy, SourceResidueId};
use crate::components::ComponentLibrary;
use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail, SourceLocation};
use crate::molecule::Molecule;
use crate::structure::{Ensemble, Model};
use crate::topology::{InstanceAtomId, MoleculeInstanceId, Topology};

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ModelSelection {
    RequireSingle,
    Select(String),
    First,
}

/// How interpretation resolves bonds the file does not state.
///
/// Each residue takes its bonds from the block's `_chem_comp_bond` rows for
/// its component when present, otherwise from a component template: the
/// caller's [`components`](Self::components) first, then
/// [`ComponentLibrary::bundled`]. A residue that neither covers has its bonds
/// inferred from covalent radii when
/// [`infer_bonds_from_distance`](Self::infer_bonds_from_distance) is set, and
/// fails otherwise.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConnectivityOptions {
    components: Option<Arc<ComponentLibrary>>,
    infer_bonds_from_distance: bool,
}

impl ConnectivityOptions {
    /// Templates consulted before the bundled ones, such as
    /// [`ComponentLibrary::amber`] for structures with AMBER residue names.
    pub fn components(&self) -> Option<&Arc<ComponentLibrary>> {
        self.components.as_ref()
    }

    /// Sets [`Self::components`].
    #[must_use]
    pub fn with_components(mut self, components: Option<Arc<ComponentLibrary>>) -> Self {
        self.components = components;
        self
    }

    /// Infer the bonds of residues without a template or bond table from
    /// covalent radii. Inferred bonds are single and reported. On by default.
    pub const fn infer_bonds_from_distance(&self) -> bool {
        self.infer_bonds_from_distance
    }

    /// Sets [`Self::infer_bonds_from_distance`].
    #[must_use]
    pub fn with_infer_bonds_from_distance(mut self, infer_bonds_from_distance: bool) -> Self {
        self.infer_bonds_from_distance = infer_bonds_from_distance;
        self
    }
}

impl Default for ConnectivityOptions {
    fn default() -> Self {
        Self {
            components: None,
            infer_bonds_from_distance: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructureInterpretOptions {
    strict_entity_metadata: bool,
    altloc_policy: AltLocPolicy,
    model_selection: ModelSelection,
    connectivity: ConnectivityOptions,
}

impl StructureInterpretOptions {
    pub const fn strict_entity_metadata(&self) -> bool {
        self.strict_entity_metadata
    }

    /// Sets [`Self::strict_entity_metadata`].
    #[must_use]
    pub fn with_strict_entity_metadata(mut self, strict_entity_metadata: bool) -> Self {
        self.strict_entity_metadata = strict_entity_metadata;
        self
    }

    pub fn altloc_policy(&self) -> &AltLocPolicy {
        &self.altloc_policy
    }

    /// Sets [`Self::altloc_policy`].
    #[must_use]
    pub fn with_altloc_policy(mut self, altloc_policy: AltLocPolicy) -> Self {
        self.altloc_policy = altloc_policy;
        self
    }

    pub fn model_selection(&self) -> &ModelSelection {
        &self.model_selection
    }

    /// Sets [`Self::model_selection`].
    #[must_use]
    pub fn with_model_selection(mut self, model_selection: ModelSelection) -> Self {
        self.model_selection = model_selection;
        self
    }

    pub fn connectivity(&self) -> &ConnectivityOptions {
        &self.connectivity
    }

    /// Sets [`Self::connectivity`].
    #[must_use]
    pub fn with_connectivity(mut self, connectivity: ConnectivityOptions) -> Self {
        self.connectivity = connectivity;
        self
    }
}

impl Default for StructureInterpretOptions {
    fn default() -> Self {
        Self {
            strict_entity_metadata: false,
            altloc_policy: AltLocPolicy::HighestOccupancy,
            model_selection: ModelSelection::RequireSingle,
            connectivity: ConnectivityOptions::default(),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum EntityKind {
    Polymer,
    Branched,
    NonPolymer,
    Water,
    Other(String),
}

impl EntityKind {
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

/// Non-fatal scientific interpretation issues retained in the report.
///
/// Connection issues preserve source identity and distinguish unresolved
/// selectors from ambiguous selectors. Neither case creates a topology bond.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StructureInterpretIssue {
    EntityTypeInferred {
        asym_id: String,
        kind: EntityKind,
    },
    /// Atom records without an element symbol had their elements inferred
    /// from the alignment of their atom names, as PDB convention places them.
    ElementsInferredFromAtomNames {
        source_lines: Vec<usize>,
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
        reason: ConnectionResolutionReason,
    },
    /// One declared connection partner matched more than one selected-model atom.
    ConnectionAmbiguous {
        connection_id: Option<String>,
        connection_type: String,
        partner: u8,
        source_line: Option<usize>,
        candidates: usize,
        reason: ConnectionResolutionReason,
    },
    CoordinateModelIgnored {
        model_id: String,
        atom_site_rows: usize,
    },
    AlternateLocationOmitted {
        atom_name: String,
        alt_id: Option<String>,
        residue: SourceResidueId,
        source_line: usize,
    },
    /// An alternative lacks sites present in another alternative of the same
    /// chemical component. It is excluded rather than completed from other labels.
    AlternateLocationIncomplete {
        residue: SourceResidueId,
        alternative: String,
        missing_atoms: Vec<String>,
    },
    AlternateLocationOccupancyMissing {
        residue: SourceResidueId,
        source_lines: Vec<usize>,
    },
    /// A source-declared combination cannot represent this residue (for example,
    /// it omits every label or includes overlapping sites). Other declared
    /// combinations can still be considered; the rejected one is never repaired.
    AlternateLocationCombinationRejected {
        residue: SourceResidueId,
        alternative: String,
        reason: String,
    },
    /// Equal heuristic scores, resolved by lexical alternative identifier.
    AlternateLocationOccupancyTie {
        residues: Vec<SourceResidueId>,
        alternatives: Vec<String>,
    },
    /// Atom pairs within covalent bonding distance that no evidence bonds,
    /// such as close contacts between residues. They are counted, not bonded.
    ConnectivityCandidatesInferred {
        atom_count: usize,
        candidate_count: usize,
    },
    /// A residue without template or bond table had its bonds inferred from
    /// covalent radii; inferred bonds are single.
    ResidueBondsInferred {
        residue: SourceResidueId,
        component_id: String,
        bonds: usize,
    },
    /// A hydrogen that no bond source names (the residue template, the
    /// block's own component bond table, or the source's listed bonds) was
    /// bonded to the nearest heavy atom of its residue.
    HydrogenAttachedByDistance {
        residue: SourceResidueId,
        atom_name: String,
        heavy_atom_name: String,
    },
    /// Residues adjacent in sequence, or joined by a branched entity's link,
    /// whose linking atoms lie beyond bonding distance are not bonded.
    ChainBreak {
        residue: SourceResidueId,
        next: SourceResidueId,
    },
    /// The model's hydrogens differ from the reference template: a protonation
    /// variant was selected, or heteroatom charges were adjusted by one per
    /// added or removed hydrogen.
    ProtonationFromHydrogens {
        residue: SourceResidueId,
        component_id: String,
        variant: Option<String>,
        charged_atoms: Vec<(String, i8)>,
    },
    /// The resolved atoms of one residue lie in more than one molecule, for
    /// example because a bridging atom is missing from the model.
    ResidueFragmented {
        residue: SourceResidueId,
        fragments: usize,
    },
    /// The source states unit-cell parameters (mmCIF `_cell`, PDB `CRYST1`)
    /// that are not the model's periodic cell, so the model has none.
    CellOmitted {
        reason: CellOmissionReason,
        location: SourceLocation,
    },
}

/// Why stated unit-cell parameters are not a model's periodic cell.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CellOmissionReason {
    /// Some of the six lengths and angles are missing or unknown.
    Incomplete,
    /// The coordinates are the asymmetric unit of this space group, whose
    /// translates alone are not the crystal.
    SpaceGroup { name: String },
    /// PDB's 1 Å cube, recorded for structures without a crystal lattice
    /// such as NMR and EM models.
    Placeholder,
}

/// The evidence a bond was resolved from.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BondEvidence {
    /// A `_chem_comp_bond` row of the interpreted block.
    ComponentTable,
    /// A component template; `variant` names a protonation variant.
    Template {
        component_id: String,
        variant: Option<String>,
    },
    /// Standard linkage between polymer residues adjacent in sequence.
    PolymerLink,
    /// A `_pdbx_entity_branch_link` row.
    BranchLink,
    /// A connection record: an mmCIF `_struct_conn` row, whose type
    /// `connection_type` names, or a PDB `LINK`, `SSBOND`, or `CONECT` record,
    /// named by its record type. Metal coordination yields a dative bond.
    Connection {
        connection_id: Option<String>,
        connection_type: String,
    },
    /// Covalent radii, for a residue without template or bond table.
    Distance,
    /// A hydrogen that no template, component bond table, or listed bond
    /// names, bonded to the nearest heavy atom of its residue.
    HydrogenDistance,
}

/// The evidence for one interpreted bond.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BondProvenance {
    pub(crate) atoms: [InstanceAtomId; 2],
    pub(crate) evidence: BondEvidence,
}

impl BondProvenance {
    pub const fn atoms(&self) -> [InstanceAtomId; 2] {
        self.atoms
    }

    pub fn evidence(&self) -> &BondEvidence {
        &self.evidence
    }
}

/// A connection-partner selector whose alias fields disagree.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ConflictingSelector {
    /// The label alternate-location aliases.
    LabelAlternateLocationAliases,
    /// The label and author alternate locations.
    LabelAndAuthorAlternateLocations,
}

impl fmt::Display for ConflictingSelector {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::LabelAlternateLocationAliases => "label alternate-location aliases",
            Self::LabelAndAuthorAlternateLocations => "label and author alternate locations",
        })
    }
}

/// Why a connection record's partner could not be resolved uniquely.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ConnectionResolutionReason {
    /// The record requests a symmetry mate outside the asymmetric unit.
    UnsupportedSymmetry { symmetry: String },
    /// No label or author atom-identity selector was supplied.
    MissingSelector,
    /// Alias fields for the same selector supplied conflicting values.
    ConflictingSelectorValues {
        /// The selector whose fields conflict.
        selector: ConflictingSelector,
    },
    /// Label and author selector families matched disjoint atom sets.
    ConflictingLabelAndAuthorSelectors,
    /// An explicitly named alternate location existed before selection but was omitted.
    AlternateLocationOmitted { alternate_location: String },
    /// No selected-model atom satisfied every supplied selector field.
    NoMatchingAtom,
    /// More than one selected-model atom satisfied every supplied selector field.
    MultipleMatchingAtoms,
}

impl fmt::Display for ConnectionResolutionReason {
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

/// Stable source correspondence for one interpreted atom.
///
/// Residue, occurrence, asymmetry, and selected alternate-location fields are
/// retained independently of derived molecule-instance and dense ordering.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomProvenance {
    pub(crate) atom: InstanceAtomId,
    pub(crate) type_symbol: String,
    pub(crate) source_line: usize,
    pub(crate) source_order: usize,
    pub(crate) atom_site_id: Option<String>,
    pub(crate) label_atom_name: Option<String>,
    pub(crate) atom_name: String,
    pub(crate) auth_atom_name: Option<String>,
    pub(crate) label_component_id: Option<String>,
    pub(crate) component_id: String,
    pub(crate) auth_component_id: Option<String>,
    pub(crate) label_asym_id: Option<String>,
    pub(crate) asym_id: String,
    pub(crate) auth_asym_id: Option<String>,
    pub(crate) entity_id: Option<String>,
    pub(crate) entity_kind: EntityKind,
    pub(crate) residue_key: String,
    pub(crate) label_sequence_id: Option<i32>,
    pub(crate) author_sequence_id: Option<String>,
    pub(crate) insertion_code: Option<String>,
    pub(crate) occurrence: Option<usize>,
    pub(crate) selected_alternate_location: Option<String>,
}

impl AtomProvenance {
    pub const fn atom(&self) -> InstanceAtomId {
        self.atom
    }

    pub const fn source_line(&self) -> usize {
        self.source_line
    }

    pub fn atom_site_id(&self) -> Option<&str> {
        self.atom_site_id.as_deref()
    }

    /// Returns `_atom_site.label_atom_id`, falling back to
    /// `_atom_site.auth_atom_id` when the label identifier is absent.
    pub fn atom_name(&self) -> &str {
        &self.atom_name
    }

    /// Returns the source `_atom_site.label_atom_id`, without applying an
    /// author-identifier fallback.
    pub fn label_atom_name(&self) -> Option<&str> {
        self.label_atom_name.as_deref()
    }

    pub fn auth_atom_name(&self) -> Option<&str> {
        self.auth_atom_name.as_deref()
    }

    /// Returns `_atom_site.label_comp_id`, falling back to
    /// `_atom_site.auth_comp_id` when the label identifier is absent.
    pub fn component_id(&self) -> &str {
        &self.component_id
    }

    /// Returns the source `_atom_site.label_comp_id`, without applying an
    /// author-identifier fallback.
    pub fn label_component_id(&self) -> Option<&str> {
        self.label_component_id.as_deref()
    }

    pub fn auth_component_id(&self) -> Option<&str> {
        self.auth_component_id.as_deref()
    }

    /// Returns `_atom_site.label_asym_id`, falling back to
    /// `_atom_site.auth_asym_id` when the label identifier is absent.
    pub fn asym_id(&self) -> &str {
        &self.asym_id
    }

    /// Returns the source `_atom_site.label_asym_id`, without applying an
    /// author-identifier fallback.
    pub fn label_asym_id(&self) -> Option<&str> {
        self.label_asym_id.as_deref()
    }

    /// Returns the independently preserved author asymmetry identifier.
    pub fn auth_asym_id(&self) -> Option<&str> {
        self.auth_asym_id.as_deref()
    }

    pub fn entity_id(&self) -> Option<&str> {
        self.entity_id.as_deref()
    }

    /// Returns the source or explicitly inferred entity kind for this atom.
    pub const fn entity_kind(&self) -> &EntityKind {
        &self.entity_kind
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
pub struct InstanceProvenance {
    pub(crate) molecule: MoleculeInstanceId,
    pub(crate) coordinate_model_id: String,
    pub(crate) asym_ids: Vec<String>,
    pub(crate) entity_ids: Vec<String>,
    pub(crate) entity_kinds: Vec<EntityKind>,
    pub(crate) atoms: Vec<AtomProvenance>,
    pub(crate) bonds: Vec<BondProvenance>,
}

impl InstanceProvenance {
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

    pub fn entity_kinds(&self) -> &[EntityKind] {
        &self.entity_kinds
    }

    pub fn atoms(&self) -> &[AtomProvenance] {
        &self.atoms
    }

    /// The evidence for every bond of the molecule, in bond order.
    pub fn bonds(&self) -> &[BondProvenance] {
        &self.bonds
    }
}

/// Structured record of one structure interpretation.
///
/// Applied connection counts include only uniquely resolved local covalent and
/// metal-coordination records. [`Self::issues`] retains source-aware unresolved and ambiguous
/// partner diagnostics.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct StructureInterpretationReport {
    pub(crate) block_name: String,
    pub(crate) entity_definitions: usize,
    pub(crate) coordinate_models: usize,
    pub(crate) selected_model: Option<String>,
    pub(crate) ignored_coordinate_models: Vec<String>,
    pub(crate) macromolecules: usize,
    pub(crate) small_molecules: usize,
    pub(crate) solvent_molecules: usize,
    pub(crate) applied_connections: usize,
    pub(crate) connectivity_candidates: usize,
    pub(crate) instances: Vec<InstanceProvenance>,
    pub(crate) alternate_locations: Vec<AltLocDecision>,
    pub(crate) issues: Vec<StructureInterpretIssue>,
}

impl StructureInterpretationReport {
    /// The interpreted mmCIF data block's name, or the PDB entry's ID code
    /// from its `HEADER` record (empty without one).
    pub fn block_name(&self) -> &str {
        &self.block_name
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

    pub fn instances(&self) -> &[InstanceProvenance] {
        &self.instances
    }

    /// Decisions for the selected source coordinate model only. The source
    /// document retains every alternative, including omitted rows and metadata.
    pub fn alternate_locations(&self) -> &[AltLocDecision] {
        &self.alternate_locations
    }

    /// Returns source-aware non-fatal interpretation issues.
    ///
    /// Unresolved and ambiguous connection partners are reported here and
    /// never create molecule merges or topology bonds.
    pub fn issues(&self) -> &[StructureInterpretIssue] {
        &self.issues
    }
}

#[derive(Debug, Clone, PartialEq)]
/// One final canonical model and its interpretation report.
pub struct StructureInterpretation {
    pub(super) model: Model,
    pub(super) report: StructureInterpretationReport,
}

impl StructureInterpretation {
    pub fn model(&self) -> &Model {
        &self.model
    }

    pub fn report(&self) -> &StructureInterpretationReport {
        &self.report
    }

    pub fn topology(&self) -> &crate::topology::Topology {
        self.model.topology()
    }

    /// Borrows molecule definitions once per occurrence in authoritative
    /// topology instance order.
    pub fn molecules(&self) -> impl ExactSizeIterator<Item = &Molecule> + DoubleEndedIterator {
        self.model
            .topology()
            .molecules()
            .map(|occurrence| occurrence.molecule())
    }

    pub fn into_molecules(self) -> Vec<Molecule> {
        self.model
            .topology()
            .molecules()
            .map(|occurrence| occurrence.molecule().clone())
            .collect()
    }

    /// Retains shared ownership of the model's exact topology allocation.
    pub fn into_topology(self) -> Arc<Topology> {
        self.model.shared_topology()
    }

    pub fn into_model(self) -> Model {
        self.model
    }

    pub fn into_parts(self) -> (Model, StructureInterpretationReport) {
        (self.model, self.report)
    }
}

/// A macromolecular interpretation failure, located at the source line it
/// concerns when it concerns one. The pipeline is shared, so errors start as
/// mmCIF errors and PDB entry points relabel them.
pub(in crate::io) fn structure_error(
    kind: DiagnosticKind,
    line: Option<usize>,
    message: impl Into<String>,
) -> FormatError {
    FormatError::interpret(Format::Mmcif, kind, message)
        .with_location(line.map(SourceLocation::at_line))
}

/// An ensemble-level failure with typed `detail`.
pub(in crate::io) fn ensemble_error(
    kind: DiagnosticKind,
    message: impl Into<String>,
    detail: FormatErrorDetail,
) -> FormatError {
    FormatError::interpret(Format::Mmcif, kind, message).with_detail(detail)
}

/// Converts deposited coordinate models to members, applying `altloc_policy`
/// within each model. Alternate locations never generate additional members or
/// statistical weights. A single coordinate model yields a single member.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructureEnsembleInterpretOptions {
    strict_entity_metadata: bool,
    altloc_policy: AltLocPolicy,
    connectivity: ConnectivityOptions,
    model_ids: Option<Vec<String>>,
}

impl StructureEnsembleInterpretOptions {
    pub const fn strict_entity_metadata(&self) -> bool {
        self.strict_entity_metadata
    }

    /// Sets [`Self::strict_entity_metadata`].
    #[must_use]
    pub fn with_strict_entity_metadata(mut self, strict_entity_metadata: bool) -> Self {
        self.strict_entity_metadata = strict_entity_metadata;
        self
    }

    pub fn altloc_policy(&self) -> &AltLocPolicy {
        &self.altloc_policy
    }

    /// Sets [`Self::altloc_policy`].
    #[must_use]
    pub fn with_altloc_policy(mut self, altloc_policy: AltLocPolicy) -> Self {
        self.altloc_policy = altloc_policy;
        self
    }

    pub fn connectivity(&self) -> &ConnectivityOptions {
        &self.connectivity
    }

    /// Sets [`Self::connectivity`].
    #[must_use]
    pub fn with_connectivity(mut self, connectivity: ConnectivityOptions) -> Self {
        self.connectivity = connectivity;
        self
    }

    /// `None` selects all coordinate models in source order. `Some` must
    /// contain at least one model ID.
    pub fn model_ids(&self) -> Option<&[String]> {
        self.model_ids.as_deref()
    }

    /// Sets [`Self::model_ids`].
    ///
    /// # Errors
    ///
    /// An interpretation [`FormatError`] with [`FormatErrorDetail::EmptyModelSelection`](crate::diagnostics::FormatErrorDetail::EmptyModelSelection) for `Some` of an empty list.
    pub fn with_model_ids(mut self, model_ids: Option<Vec<String>>) -> Result<Self, FormatError> {
        if model_ids.as_ref().is_some_and(Vec::is_empty) {
            return Err(super::ensemble::empty_model_selection());
        }
        self.model_ids = model_ids;
        Ok(self)
    }
}

impl Default for StructureEnsembleInterpretOptions {
    fn default() -> Self {
        Self {
            strict_entity_metadata: false,
            altloc_policy: AltLocPolicy::HighestOccupancy,
            connectivity: ConnectivityOptions::default(),
            model_ids: None,
        }
    }
}

#[derive(Debug, Clone)]
/// One shared-topology ensemble and one report per selected conformation.
pub struct StructureEnsembleInterpretation {
    pub(super) ensemble: Ensemble,
    pub(super) reports: Vec<StructureInterpretationReport>,
}

impl StructureEnsembleInterpretation {
    pub fn ensemble(&self) -> &Ensemble {
        &self.ensemble
    }

    pub fn reports(&self) -> &[StructureInterpretationReport] {
        &self.reports
    }

    pub fn topology(&self) -> &Topology {
        self.ensemble.topology()
    }

    pub fn into_ensemble(self) -> Ensemble {
        self.ensemble
    }

    pub fn into_parts(self) -> (Ensemble, Vec<StructureInterpretationReport>) {
        (self.ensemble, self.reports)
    }
}
