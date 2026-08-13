use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::sync::Arc;

use crate::bio::{MacroMolecule, SmcraAtomSiteMetadata, SmcraHierarchy};
use crate::core::{Atom, AtomId, BondOrder, Conformer, ConformerId, Element, Molecule};
use crate::geometry::Point3;
use crate::small::model::SmallMolecule;
use crate::structure::{
    AtomData, Ensemble, EnsembleError, EnsembleMember, Model, ModelBuilder, Positions,
};
use crate::topology::{
    InstanceAtomId, MoleculeInstanceId, MoleculeInstanceMetadata, MoleculeRole, TopologyMapping,
};
use crate::units::{Quantity, ANGSTROM};

use super::{MmcifDataBlock, MmcifDocument, MmcifLoopTable, MmcifValue};

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
    fn from_mmcif(value: &str) -> Self {
        match value.to_ascii_lowercase().as_str() {
            "polymer" => Self::Polymer,
            "branched" => Self::Branched,
            "non-polymer" => Self::NonPolymer,
            "water" => Self::Water,
            _ => Self::Other(value.to_owned()),
        }
    }

    fn is_macro(&self) -> bool {
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
    pub(crate) source_line: usize,
    pub(crate) atom_site_id: Option<String>,
    pub(crate) atom_name: String,
    pub(crate) component_id: String,
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

    pub fn component_id(&self) -> &str {
        &self.component_id
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

    /// Returns the alternate-location identity selected for this source atom.
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
    pub(crate) template_bonds_pending: usize,
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

    pub const fn template_bonds_pending(&self) -> usize {
        self.template_bonds_pending
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
pub struct MmcifInterpretation {
    model: Model,
    report: MmcifInterpretationReport,
}

impl MmcifInterpretation {
    pub fn into_parts(self) -> (Model, MmcifInterpretationReport) {
        (self.model, self.report)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifInterpretError {
    pub(crate) line: Option<usize>,
    pub(crate) message: String,
}

impl MmcifInterpretError {
    fn new(line: Option<usize>, message: impl Into<String>) -> Self {
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

pub fn interpret_mmcif(
    document: &MmcifDocument,
    options: MmcifInterpretOptions,
) -> Result<MmcifInterpretation, MmcifInterpretError> {
    let blocks = document
        .blocks()
        .iter()
        .filter(|block| block.loop_with_tag("_atom_site.type_symbol").is_some())
        .collect::<Vec<_>>();
    if blocks.is_empty() {
        return Err(MmcifInterpretError::new(
            None,
            "document has no atom-site loop",
        ));
    }
    if blocks.len() > 1 {
        return Err(MmcifInterpretError::new(
            None,
            "document has atom-site data in more than one data block",
        ));
    }
    interpret_block(blocks[0], options)
}

fn interpret_block(
    block: &MmcifDataBlock,
    options: MmcifInterpretOptions,
) -> Result<MmcifInterpretation, MmcifInterpretError> {
    let entities = read_entity_types(block)?;
    let asym_entities = read_asym_entities(block)?;
    let atom_table = block
        .loop_with_tag("_atom_site.type_symbol")
        .expect("selected block has atom-site data");
    if atom_table.row_count() == 0 {
        return Err(MmcifInterpretError::new(
            None,
            "atom-site loop contains no rows",
        ));
    }
    let mut report = MmcifInterpretationReport {
        data_block: block.name().to_owned(),
        entity_definitions: entities.len(),
        ..MmcifInterpretationReport::default()
    };
    let rows = read_atom_rows(atom_table, &entities, &asym_entities, &options, &mut report)?;
    let selected = select_alt_locations(&rows, &options.altloc_policy, &mut report)?;
    let selected = select_coordinate_model(selected, &options.model_selection, &mut report)?;
    let selected_model = report
        .selected_model
        .clone()
        .ok_or_else(|| MmcifInterpretError::new(None, "coordinate model selection was lost"))?;
    let mut union = InstanceUnion::new(selected.iter().map(|row| row.instance_key.clone()));
    let connections = read_connections(
        block,
        &selected,
        &rows,
        &selected_model,
        &mut union,
        &mut report,
    )?;
    let polymer_asym_order = polymer_asym_order(block);
    let groups = group_rows(selected, &mut union, &polymer_asym_order);
    let mut builder = ModelBuilder::new();
    let mut qualified_atom_data = Vec::new();
    for group in groups {
        let built = build_molecule(group, &connections, &mut report)?;
        match built {
            BuiltMolecule::Macro {
                molecule,
                conformer,
                metadata,
                provenance,
            } => {
                let id = builder
                    .add_macro_molecule_with_metadata_unchecked_connectedness(
                        &molecule, conformer, metadata,
                    )
                    .map_err(graph_error)?;
                let (provenance, atom_data) = provenance.qualify(id);
                report.instances.push(provenance);
                qualified_atom_data.extend(atom_data);
            }
            BuiltMolecule::Small {
                molecule,
                conformer,
                metadata,
                provenance,
            } => {
                let id = builder
                    .add_small_molecule_with_metadata_unchecked_connectedness(
                        &molecule, conformer, metadata,
                    )
                    .map_err(graph_error)?;
                let (provenance, atom_data) = provenance.qualify(id);
                report.instances.push(provenance);
                qualified_atom_data.extend(atom_data);
            }
        }
    }
    let mut model = builder.build().map_err(graph_error)?;
    let mut occupancies = vec![None; model.topology().atom_count()];
    let mut b_factors = vec![None; model.topology().atom_count()];
    for (atom, occupancy, b_factor) in qualified_atom_data {
        let index = model
            .topology()
            .atom_index(atom)
            .expect("interpreted atom has a dense topology index");
        occupancies[index.index()] = occupancy;
        b_factors[index.index()] = b_factor;
    }
    let topology = model.shared_topology();
    let atom_data = AtomData::from_columns(&topology, Some(occupancies), Some(b_factors))
        .map_err(graph_error)?;
    model.set_atom_data(atom_data).map_err(graph_error)?;
    report.macromolecules = model
        .topology()
        .instances()
        .filter(|(id, _)| {
            model
                .topology()
                .definition_for_instance(*id)
                .is_ok_and(|definition| definition.macro_molecule().is_some())
        })
        .count();
    report.small_molecules = model
        .topology()
        .instances()
        .filter(|(id, _)| {
            model
                .topology()
                .definition_for_instance(*id)
                .is_ok_and(|definition| definition.small_molecule().is_some())
        })
        .count();
    report.solvent_molecules = model
        .topology()
        .instances()
        .filter(|(_, molecule)| molecule.has_role(MoleculeRole::Solvent))
        .count();
    Ok(MmcifInterpretation { model, report })
}

fn coordinate_model_ids(block: &MmcifDataBlock) -> Result<Vec<String>, MmcifInterpretError> {
    let table = block
        .loop_with_tag("_atom_site.type_symbol")
        .ok_or_else(|| MmcifInterpretError::new(None, "data block has no atom-site loop"))?;
    let mut seen = BTreeSet::new();
    let mut models = Vec::new();
    for row in 0..table.row_count() {
        let model = optional(table, row, "_atom_site.pdbx_PDB_model_num")
            .unwrap_or("1")
            .to_owned();
        if seen.insert(model.clone()) {
            models.push(model);
        }
    }
    Ok(models)
}

fn read_entity_types(
    block: &MmcifDataBlock,
) -> Result<BTreeMap<String, MmcifEntityKind>, MmcifInterpretError> {
    let mut entities = BTreeMap::new();
    if let Some(table) = block.loop_with_tag("_entity.id") {
        for row in 0..table.row_count() {
            let id = required(table, row, "_entity.id")?;
            let kind = required(table, row, "_entity.type")?;
            if entities
                .insert(id.to_owned(), MmcifEntityKind::from_mmcif(kind))
                .is_some()
            {
                return Err(row_error(table, row, format!("duplicate entity `{id}`")));
            }
        }
    } else if let (Some(id), Some(kind)) = (
        block.item("_entity.id").and_then(MmcifValue::optional_text),
        block
            .item("_entity.type")
            .and_then(MmcifValue::optional_text),
    ) {
        entities.insert(id.to_owned(), MmcifEntityKind::from_mmcif(kind));
    }
    Ok(entities)
}

fn read_asym_entities(
    block: &MmcifDataBlock,
) -> Result<BTreeMap<String, String>, MmcifInterpretError> {
    let mut instances = BTreeMap::new();
    if let Some(table) = block.loop_with_tag("_struct_asym.id") {
        for row in 0..table.row_count() {
            let id = required(table, row, "_struct_asym.id")?;
            let entity = required(table, row, "_struct_asym.entity_id")?;
            if instances.insert(id.to_owned(), entity.to_owned()).is_some() {
                return Err(row_error(
                    table,
                    row,
                    format!("duplicate structural instance `{id}`"),
                ));
            }
        }
    } else if let (Some(id), Some(entity)) = (
        block
            .item("_struct_asym.id")
            .and_then(MmcifValue::optional_text),
        block
            .item("_struct_asym.entity_id")
            .and_then(MmcifValue::optional_text),
    ) {
        instances.insert(id.to_owned(), entity.to_owned());
    }
    Ok(instances)
}

#[derive(Debug, Clone)]
struct AtomRow {
    line: usize,
    row_index: usize,
    model_id: String,
    entity_id: Option<String>,
    kind: MmcifEntityKind,
    instance_key: String,
    label_asym_id: Option<String>,
    asym_id: String,
    auth_asym_id: Option<String>,
    residue_key: String,
    label_seq_id: Option<i32>,
    auth_seq_id: Option<String>,
    insertion_code: Option<String>,
    occurrence: Option<usize>,
    label_comp_id: Option<String>,
    comp_id: String,
    auth_comp_id: Option<String>,
    label_atom_name: Option<String>,
    atom_name: String,
    auth_atom_name: Option<String>,
    atom_site_id: Option<String>,
    alt_id: Option<String>,
    occupancy: Option<f64>,
    b_factor: Option<f64>,
    point: Option<Point3>,
    element: Element,
    formal_charge: i8,
}

impl AtomRow {
    fn atom_key(&self) -> String {
        format!("{}|{}|{}", self.asym_id, self.residue_key, self.atom_name)
    }
}

#[derive(Debug, Default)]
struct OccurrenceState {
    occurrence: usize,
    seen: BTreeMap<String, BTreeSet<Option<String>>>,
}

fn read_atom_rows(
    table: &MmcifLoopTable,
    entities: &BTreeMap<String, MmcifEntityKind>,
    asym_entities: &BTreeMap<String, String>,
    options: &MmcifInterpretOptions,
    report: &mut MmcifInterpretationReport,
) -> Result<Vec<AtomRow>, MmcifInterpretError> {
    let mut rows = Vec::with_capacity(table.row_count());
    let mut occurrences = BTreeMap::<(String, String, String), OccurrenceState>::new();
    let mut inferred = BTreeSet::new();
    for row in 0..table.row_count() {
        let type_symbol = required(table, row, "_atom_site.type_symbol")?;
        let type_value = table
            .value(row, "_atom_site.type_symbol")
            .expect("required");
        let element = Element::from_symbol(&canonical_mmcif_element_symbol(type_symbol))
            .ok_or_else(|| {
                MmcifInterpretError::new(
                    Some(type_value.line()),
                    format!("unknown atom-site element `{type_symbol}`"),
                )
            })?;
        let label_asym_id = optional(table, row, "_atom_site.label_asym_id").map(str::to_owned);
        let asym_id = label_asym_id
            .as_deref()
            .or_else(|| optional(table, row, "_atom_site.auth_asym_id"))
            .ok_or_else(|| row_error(table, row, "missing atom-site chain identifier"))?
            .to_owned();
        let auth_asym_id = optional(table, row, "_atom_site.auth_asym_id").map(str::to_owned);
        let label_comp_id = optional(table, row, "_atom_site.label_comp_id").map(str::to_owned);
        let comp_id = label_comp_id
            .as_deref()
            .or_else(|| optional(table, row, "_atom_site.auth_comp_id"))
            .ok_or_else(|| row_error(table, row, "missing atom-site component identifier"))?
            .to_owned();
        let label_atom_name = optional(table, row, "_atom_site.label_atom_id").map(str::to_owned);
        let atom_name = label_atom_name
            .as_deref()
            .or_else(|| optional(table, row, "_atom_site.auth_atom_id"))
            .ok_or_else(|| row_error(table, row, "missing atom-site atom identifier"))?
            .to_owned();
        let model_id = optional(table, row, "_atom_site.pdbx_PDB_model_num")
            .unwrap_or("1")
            .to_owned();
        let entity_id = optional(table, row, "_atom_site.label_entity_id")
            .map(str::to_owned)
            .or_else(|| asym_entities.get(&asym_id).cloned());
        let group_pdb = optional(table, row, "_atom_site.group_PDB").map(str::to_owned);
        let kind = entity_id
            .as_ref()
            .and_then(|entity| entities.get(entity))
            .cloned()
            .unwrap_or_else(|| infer_entity_kind(group_pdb.as_deref(), &comp_id));
        if entity_id
            .as_ref()
            .and_then(|entity| entities.get(entity))
            .is_none()
        {
            if options.strict_entity_metadata {
                return Err(row_error(
                    table,
                    row,
                    format!("missing entity type for structural instance `{asym_id}`"),
                ));
            }
            if inferred.insert(asym_id.clone()) {
                report.issues.push(MmcifInterpretIssue::EntityTypeInferred {
                    asym_id: asym_id.clone(),
                    kind: kind.clone(),
                });
            }
        }
        let label_seq_id = optional_i32(table, row, "_atom_site.label_seq_id")?;
        let auth_seq_id = optional(table, row, "_atom_site.auth_seq_id").map(str::to_owned);
        let insertion_code =
            optional(table, row, "_atom_site.pdbx_PDB_ins_code").map(str::to_owned);
        let alt_id = optional(table, row, "_atom_site.label_alt_id").map(str::to_owned);
        let (residue_key, occurrence) = if let Some(sequence) = label_seq_id {
            (
                format!(
                    "label:{sequence}:{}",
                    insertion_code.as_deref().unwrap_or("")
                ),
                None,
            )
        } else if let Some(sequence) = &auth_seq_id {
            (
                format!(
                    "auth:{sequence}:{}",
                    insertion_code.as_deref().unwrap_or("")
                ),
                None,
            )
        } else {
            let state = occurrences
                .entry((model_id.clone(), asym_id.clone(), comp_id.clone()))
                .or_default();
            let prior = state.seen.get(&atom_name);
            let repeats = prior.is_some_and(|labels| {
                alt_id.is_none() || labels.contains(&None) || labels.contains(&alt_id)
            });
            if repeats {
                state.occurrence += 1;
                state.seen.clear();
            }
            state
                .seen
                .entry(atom_name.clone())
                .or_default()
                .insert(alt_id.clone());
            (
                format!("occurrence:{}", state.occurrence),
                Some(state.occurrence),
            )
        };
        let instance_key = if kind.is_macro() {
            format!("macro:{asym_id}")
        } else {
            format!("small:{asym_id}:{residue_key}")
        };
        let formal_charge =
            optional_i8(table, row, "_atom_site.pdbx_formal_charge")?.unwrap_or_default();
        let point = optional_point(table, row)?;
        rows.push(AtomRow {
            line: type_value.line(),
            row_index: row,
            model_id,
            entity_id,
            kind,
            instance_key,
            label_asym_id,
            asym_id,
            auth_asym_id,
            residue_key,
            label_seq_id,
            auth_seq_id,
            insertion_code,
            occurrence,
            label_comp_id,
            comp_id,
            auth_comp_id: optional(table, row, "_atom_site.auth_comp_id").map(str::to_owned),
            label_atom_name,
            atom_name,
            auth_atom_name: optional(table, row, "_atom_site.auth_atom_id").map(str::to_owned),
            atom_site_id: optional(table, row, "_atom_site.id").map(str::to_owned),
            alt_id,
            occupancy: optional_f64(table, row, "_atom_site.occupancy")?,
            b_factor: optional_f64(table, row, "_atom_site.B_iso_or_equiv")?,
            point,
            element,
            formal_charge,
        });
    }
    Ok(rows)
}

fn infer_entity_kind(group_pdb: Option<&str>, comp_id: &str) -> MmcifEntityKind {
    if ["HOH", "WAT", "DOD"]
        .iter()
        .any(|water| comp_id.eq_ignore_ascii_case(water))
    {
        MmcifEntityKind::Water
    } else if group_pdb.is_some_and(|group| group.eq_ignore_ascii_case("ATOM")) {
        MmcifEntityKind::Polymer
    } else {
        MmcifEntityKind::NonPolymer
    }
}

fn canonical_mmcif_element_symbol(symbol: &str) -> String {
    let mut chars = symbol.chars();
    let Some(first) = chars.next() else {
        return String::new();
    };
    let mut canonical = first.to_ascii_uppercase().to_string();
    canonical.extend(chars.flat_map(char::to_lowercase));
    canonical
}

fn select_alt_locations(
    rows: &[AtomRow],
    policy: &MmcifAltLocPolicy,
    report: &mut MmcifInterpretationReport,
) -> Result<Vec<AtomRow>, MmcifInterpretError> {
    let mut grouped = BTreeMap::<(String, String, String), Vec<&AtomRow>>::new();
    for row in rows {
        grouped
            .entry((
                row.instance_key.clone(),
                row.atom_key(),
                row.model_id.clone(),
            ))
            .or_default()
            .push(row);
    }
    let mut selected = Vec::new();
    for (_, mut candidates) in grouped {
        candidates.sort_by_key(|row| row.row_index);
        let mut identities = BTreeSet::new();
        if let Some(duplicate) = candidates
            .iter()
            .find(|row| !identities.insert(row.alt_id.clone()))
        {
            return Err(MmcifInterpretError::new(
                Some(duplicate.line),
                format!(
                    "atom `{}` has duplicate records for one alternate location",
                    duplicate.atom_name
                ),
            ));
        }
        let labels = candidates
            .iter()
            .filter_map(|row| row.alt_id.clone())
            .collect::<BTreeSet<_>>();
        if candidates.len() > 1
            && !labels.is_empty()
            && matches!(policy, MmcifAltLocPolicy::ErrorOnAlternateLocations)
        {
            return Err(MmcifInterpretError::new(
                Some(candidates[0].line),
                format!("atom `{}` has alternate locations", candidates[0].atom_name),
            ));
        }
        let chosen = match policy {
            MmcifAltLocPolicy::HighestOccupancy => candidates
                .iter()
                .max_by(|left, right| {
                    left.occupancy
                        .unwrap_or(0.0)
                        .partial_cmp(&right.occupancy.unwrap_or(0.0))
                        .unwrap_or(Ordering::Equal)
                        .then_with(|| right.alt_id.cmp(&left.alt_id))
                })
                .map(|row| (**row).clone()),
            MmcifAltLocPolicy::SelectLabel(label) => candidates
                .iter()
                .find(|row| row.alt_id.as_deref() == Some(label.as_str()))
                .map(|row| (**row).clone())
                .or_else(|| {
                    candidates
                        .iter()
                        .find(|row| row.alt_id.is_none())
                        .map(|row| (**row).clone())
                }),
            MmcifAltLocPolicy::ErrorOnAlternateLocations => {
                candidates.first().map(|row| (**row).clone())
            }
        };
        let Some(chosen) = chosen else {
            return Err(MmcifInterpretError::new(
                None,
                "requested alternate-location label is unavailable",
            ));
        };
        for omitted in candidates
            .iter()
            .filter(|candidate| candidate.row_index != chosen.row_index)
        {
            report
                .issues
                .push(MmcifInterpretIssue::AlternateLocationOmitted {
                    atom_name: omitted.atom_name.clone(),
                    alt_id: omitted.alt_id.clone(),
                });
        }
        selected.push(chosen);
    }
    selected.sort_by_key(|row| row.row_index);
    Ok(selected)
}

fn select_coordinate_model(
    rows: Vec<AtomRow>,
    selection: &MmcifModelSelection,
    report: &mut MmcifInterpretationReport,
) -> Result<Vec<AtomRow>, MmcifInterpretError> {
    let mut model_ids = Vec::new();
    let mut seen = BTreeSet::new();
    for row in &rows {
        if seen.insert(row.model_id.clone()) {
            model_ids.push(row.model_id.clone());
        }
    }
    report.coordinate_models = model_ids.len();
    let selected = match selection {
        MmcifModelSelection::RequireSingle if model_ids.len() == 1 => model_ids[0].clone(),
        MmcifModelSelection::RequireSingle => {
            return Err(MmcifInterpretError::new(
                None,
                format!(
                    "coordinate data contains {} models; select one explicitly",
                    model_ids.len()
                ),
            ));
        }
        MmcifModelSelection::Select(id) if seen.contains(id) => id.clone(),
        MmcifModelSelection::Select(id) => {
            return Err(MmcifInterpretError::new(
                None,
                format!("coordinate model `{id}` is unavailable"),
            ));
        }
        MmcifModelSelection::First => model_ids
            .first()
            .cloned()
            .ok_or_else(|| MmcifInterpretError::new(None, "coordinate data contains no models"))?,
    };
    report.selected_model = Some(selected.clone());
    for ignored in model_ids.iter().filter(|id| **id != selected) {
        let atom_site_rows = rows.iter().filter(|row| row.model_id == *ignored).count();
        report.ignored_coordinate_models.push(ignored.clone());
        report
            .issues
            .push(MmcifInterpretIssue::CoordinateModelIgnored {
                model_id: ignored.clone(),
                atom_site_rows,
            });
    }
    let selected_rows = rows
        .into_iter()
        .filter(|row| row.model_id == selected)
        .collect::<Vec<_>>();
    if let Some(row) = selected_rows.iter().find(|row| row.point.is_none()) {
        return Err(MmcifInterpretError::new(
            Some(row.line),
            format!(
                "selected coordinate model `{selected}` has no complete position for atom `{}`",
                row.atom_name
            ),
        ));
    }
    Ok(selected_rows)
}

#[derive(Debug, Clone)]
struct DeclaredConnection {
    left_atom: String,
    right_atom: String,
    order: BondOrder,
}

fn read_connections(
    block: &MmcifDataBlock,
    rows: &[AtomRow],
    all_rows: &[AtomRow],
    selected_model: &str,
    union: &mut InstanceUnion,
    report: &mut MmcifInterpretationReport,
) -> Result<Vec<DeclaredConnection>, MmcifInterpretError> {
    let Some(table) = block.loop_with_tag("_struct_conn.conn_type_id") else {
        return Ok(Vec::new());
    };
    let mut connections = Vec::new();
    for row in 0..table.row_count() {
        let kind = required(table, row, "_struct_conn.conn_type_id")?.to_owned();
        let connection_id = optional(table, row, "_struct_conn.id").map(str::to_owned);
        let source_line = table
            .row(row)
            .and_then(|values| values.first())
            .map(MmcifValue::line);
        if !is_covalent_connection(&kind) {
            report.issues.push(MmcifInterpretIssue::ConnectionIgnored {
                connection_type: kind,
            });
            continue;
        }
        let order = connection_bond_order(table, row)?;
        let left = connection_partner(table, row, 1, rows, all_rows, selected_model)?;
        let right = connection_partner(table, row, 2, rows, all_rows, selected_model)?;
        let left = report_connection_partner_resolution(
            left,
            connection_id.as_deref(),
            &kind,
            1,
            source_line,
            report,
        );
        let right = report_connection_partner_resolution(
            right,
            connection_id.as_deref(),
            &kind,
            2,
            source_line,
            report,
        );
        let (Some(left), Some(right)) = (left, right) else {
            continue;
        };
        union.union(&left.instance_key, &right.instance_key);
        connections.push(DeclaredConnection {
            left_atom: left.atom_key(),
            right_atom: right.atom_key(),
            order,
        });
        report.applied_connections += 1;
    }
    Ok(connections)
}

fn report_connection_partner_resolution<'a>(
    resolution: ConnectionPartnerResolution<'a>,
    connection_id: Option<&str>,
    connection_type: &str,
    partner: u8,
    source_line: Option<usize>,
    report: &mut MmcifInterpretationReport,
) -> Option<&'a AtomRow> {
    match resolution {
        ConnectionPartnerResolution::Resolved(atom) => Some(atom),
        ConnectionPartnerResolution::Unresolved(reason) => {
            report
                .issues
                .push(MmcifInterpretIssue::ConnectionUnresolved {
                    connection_id: connection_id.map(str::to_owned),
                    connection_type: connection_type.to_owned(),
                    partner,
                    source_line,
                    reason,
                });
            None
        }
        ConnectionPartnerResolution::Ambiguous { candidates, reason } => {
            report
                .issues
                .push(MmcifInterpretIssue::ConnectionAmbiguous {
                    connection_id: connection_id.map(str::to_owned),
                    connection_type: connection_type.to_owned(),
                    partner,
                    source_line,
                    candidates,
                    reason,
                });
            None
        }
    }
}

fn connection_bond_order(
    table: &MmcifLoopTable,
    row: usize,
) -> Result<BondOrder, MmcifInterpretError> {
    let Some(order) = optional(table, row, "_struct_conn.pdbx_value_order") else {
        return Ok(BondOrder::Single);
    };
    match order.to_ascii_lowercase().as_str() {
        "sing" => Ok(BondOrder::Single),
        "doub" => Ok(BondOrder::Double),
        "trip" => Ok(BondOrder::Triple),
        "quad" => Ok(BondOrder::Quadruple),
        _ => Err(row_error(
            table,
            row,
            format!("unsupported struct_conn bond order `{order}`"),
        )),
    }
}

fn is_covalent_connection(kind: &str) -> bool {
    let kind = kind.to_ascii_lowercase();
    kind.starts_with("covale") || kind == "disulf" || kind == "modres"
}

#[derive(Debug, Default)]
struct LabelConnectionSelector {
    asym_id: Option<String>,
    component_id: Option<String>,
    sequence_id: Option<i32>,
    atom_id: Option<String>,
    alternate_location: Option<String>,
}

impl LabelConnectionSelector {
    fn is_empty(&self) -> bool {
        self.asym_id.is_none()
            && self.component_id.is_none()
            && self.sequence_id.is_none()
            && self.atom_id.is_none()
            && self.alternate_location.is_none()
    }

    fn matches(&self, candidate: &AtomRow) -> bool {
        self.asym_id
            .as_deref()
            .is_none_or(|expected| candidate.label_asym_id.as_deref() == Some(expected))
            && self
                .component_id
                .as_deref()
                .is_none_or(|expected| candidate.label_comp_id.as_deref() == Some(expected))
            && self
                .sequence_id
                .is_none_or(|expected| candidate.label_seq_id == Some(expected))
            && self
                .atom_id
                .as_deref()
                .is_none_or(|expected| candidate.label_atom_name.as_deref() == Some(expected))
            && self
                .alternate_location
                .as_deref()
                .is_none_or(|expected| candidate.alt_id.as_deref() == Some(expected))
    }
}

#[derive(Debug, Default)]
struct AuthorConnectionSelector {
    asym_id: Option<String>,
    component_id: Option<String>,
    sequence_id: Option<String>,
    atom_id: Option<String>,
    insertion_code: Option<String>,
    alternate_location: Option<String>,
}

impl AuthorConnectionSelector {
    fn is_empty(&self) -> bool {
        self.asym_id.is_none()
            && self.component_id.is_none()
            && self.sequence_id.is_none()
            && self.atom_id.is_none()
            && self.insertion_code.is_none()
            && self.alternate_location.is_none()
    }

    fn matches(&self, candidate: &AtomRow) -> bool {
        self.asym_id
            .as_deref()
            .is_none_or(|expected| candidate.auth_asym_id.as_deref() == Some(expected))
            && self
                .component_id
                .as_deref()
                .is_none_or(|expected| candidate.auth_comp_id.as_deref() == Some(expected))
            && self
                .sequence_id
                .as_deref()
                .is_none_or(|expected| candidate.auth_seq_id.as_deref() == Some(expected))
            && self
                .atom_id
                .as_deref()
                .is_none_or(|expected| candidate.auth_atom_name.as_deref() == Some(expected))
            && self
                .insertion_code
                .as_deref()
                .is_none_or(|expected| candidate.insertion_code.as_deref() == Some(expected))
            && self
                .alternate_location
                .as_deref()
                .is_none_or(|expected| candidate.alt_id.as_deref() == Some(expected))
    }
}

#[derive(Debug, Default)]
struct ConnectionPartnerSelector {
    label: LabelConnectionSelector,
    author: AuthorConnectionSelector,
    conflict: Option<MmcifConnectionResolutionReason>,
}

impl ConnectionPartnerSelector {
    fn is_empty(&self) -> bool {
        self.label.is_empty() && self.author.is_empty()
    }

    fn matches(&self, candidate: &AtomRow) -> bool {
        self.label.matches(candidate) && self.author.matches(candidate)
    }

    fn explicit_alternate_location(&self) -> Option<&str> {
        self.label
            .alternate_location
            .as_deref()
            .or(self.author.alternate_location.as_deref())
    }
}

#[derive(Debug)]
enum ConnectionPartnerResolution<'a> {
    Resolved(&'a AtomRow),
    Unresolved(MmcifConnectionResolutionReason),
    Ambiguous {
        candidates: usize,
        reason: MmcifConnectionResolutionReason,
    },
}

fn connection_partner<'a>(
    table: &MmcifLoopTable,
    row: usize,
    partner: u8,
    selected_rows: &'a [AtomRow],
    all_rows: &[AtomRow],
    selected_model: &str,
) -> Result<ConnectionPartnerResolution<'a>, MmcifInterpretError> {
    let symmetry_tag = format!("_struct_conn.ptnr{partner}_symmetry");
    if let Some(symmetry) = optional(table, row, &symmetry_tag) {
        if symmetry != "1_555" {
            return Ok(ConnectionPartnerResolution::Unresolved(
                MmcifConnectionResolutionReason::UnsupportedSymmetry {
                    symmetry: symmetry.to_owned(),
                },
            ));
        }
    }

    let label_alt_tag = format!("_struct_conn.ptnr{partner}_label_alt_id");
    let pdbx_label_alt_tag = format!("_struct_conn.pdbx_ptnr{partner}_label_alt_id");
    let label_alt = optional(table, row, &label_alt_tag);
    let pdbx_label_alt = optional(table, row, &pdbx_label_alt_tag);
    let (alternate_location, conflict) = match (label_alt, pdbx_label_alt) {
        (Some(left), Some(right)) if left != right => (
            Some(left.to_owned()),
            Some(MmcifConnectionResolutionReason::ConflictingSelectorValues {
                selector: "label alternate-location aliases",
            }),
        ),
        (Some(value), _) | (_, Some(value)) => (Some(value.to_owned()), None),
        (None, None) => (None, None),
    };

    let mut selector = ConnectionPartnerSelector {
        label: LabelConnectionSelector {
            asym_id: optional(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_label_asym_id"),
            )
            .map(str::to_owned),
            component_id: optional(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_label_comp_id"),
            )
            .map(str::to_owned),
            sequence_id: optional_i32(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_label_seq_id"),
            )?,
            atom_id: optional(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_label_atom_id"),
            )
            .map(str::to_owned),
            alternate_location,
        },
        author: AuthorConnectionSelector {
            asym_id: optional(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_auth_asym_id"),
            )
            .map(str::to_owned),
            component_id: optional(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_auth_comp_id"),
            )
            .map(str::to_owned),
            sequence_id: optional(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_auth_seq_id"),
            )
            .map(str::to_owned),
            atom_id: optional(
                table,
                row,
                &format!("_struct_conn.ptnr{partner}_auth_atom_id"),
            )
            .map(str::to_owned),
            insertion_code: optional(
                table,
                row,
                &format!("_struct_conn.pdbx_ptnr{partner}_PDB_ins_code"),
            )
            .map(str::to_owned),
            alternate_location: optional(
                table,
                row,
                &format!("_struct_conn.pdbx_ptnr{partner}_auth_alt_id"),
            )
            .map(str::to_owned),
        },
        conflict,
    };
    if let (Some(label), Some(author)) = (
        selector.label.alternate_location.as_deref(),
        selector.author.alternate_location.as_deref(),
    ) {
        if label != author {
            selector.conflict = Some(MmcifConnectionResolutionReason::ConflictingSelectorValues {
                selector: "label and author alternate locations",
            });
        }
    }
    if let Some(reason) = selector.conflict.take() {
        return Ok(ConnectionPartnerResolution::Unresolved(reason));
    }
    if selector.is_empty() {
        return Ok(ConnectionPartnerResolution::Unresolved(
            MmcifConnectionResolutionReason::MissingSelector,
        ));
    }

    let candidates = selected_rows
        .iter()
        .filter(|candidate| selector.matches(candidate))
        .collect::<Vec<_>>();
    match candidates.as_slice() {
        [candidate] => Ok(ConnectionPartnerResolution::Resolved(candidate)),
        [] => {
            if let Some(alternate_location) = selector.explicit_alternate_location() {
                if all_rows.iter().any(|candidate| {
                    candidate.model_id == selected_model && selector.matches(candidate)
                }) {
                    return Ok(ConnectionPartnerResolution::Unresolved(
                        MmcifConnectionResolutionReason::AlternateLocationOmitted {
                            alternate_location: alternate_location.to_owned(),
                        },
                    ));
                }
            }
            if !selector.label.is_empty()
                && !selector.author.is_empty()
                && selected_rows
                    .iter()
                    .any(|candidate| selector.label.matches(candidate))
                && selected_rows
                    .iter()
                    .any(|candidate| selector.author.matches(candidate))
            {
                return Ok(ConnectionPartnerResolution::Unresolved(
                    MmcifConnectionResolutionReason::ConflictingLabelAndAuthorSelectors,
                ));
            }
            Ok(ConnectionPartnerResolution::Unresolved(
                MmcifConnectionResolutionReason::NoMatchingAtom,
            ))
        }
        candidates => Ok(ConnectionPartnerResolution::Ambiguous {
            candidates: candidates.len(),
            reason: MmcifConnectionResolutionReason::MultipleMatchingAtoms,
        }),
    }
}

#[derive(Debug)]
struct MoleculeGroup {
    rows: Vec<AtomRow>,
    kinds: BTreeSet<MmcifEntityKind>,
    instance_keys: BTreeSet<String>,
}

fn polymer_asym_order(block: &MmcifDataBlock) -> BTreeMap<String, usize> {
    let mut order = BTreeMap::new();
    let Some(table) = block.loop_with_tag("_pdbx_poly_seq_scheme.asym_id") else {
        return order;
    };
    for row in 0..table.row_count() {
        if let Some(asym_id) = optional(table, row, "_pdbx_poly_seq_scheme.asym_id") {
            let next = order.len();
            order.entry(asym_id.to_owned()).or_insert(next);
        }
    }
    order
}

fn group_rows(
    rows: Vec<AtomRow>,
    union: &mut InstanceUnion,
    polymer_asym_order: &BTreeMap<String, usize>,
) -> Vec<MoleculeGroup> {
    let mut group_indices = BTreeMap::new();
    let mut groups = Vec::new();
    for row in rows {
        let root = union.find(&row.instance_key);
        let index = *group_indices.entry(root).or_insert_with(|| {
            groups.push(MoleculeGroup {
                rows: Vec::new(),
                kinds: BTreeSet::new(),
                instance_keys: BTreeSet::new(),
            });
            groups.len() - 1
        });
        let group = &mut groups[index];
        group.kinds.insert(row.kind.clone());
        group.instance_keys.insert(row.instance_key.clone());
        group.rows.push(row);
    }
    for group in &mut groups {
        group.rows.sort_by_key(|row| {
            polymer_asym_order
                .get(&row.asym_id)
                .copied()
                .unwrap_or(usize::MAX)
        });
    }
    groups.sort_by_key(|group| {
        group
            .rows
            .iter()
            .filter_map(|row| polymer_asym_order.get(&row.asym_id).copied())
            .min()
            .unwrap_or(usize::MAX)
    });
    groups
}

enum BuiltMolecule {
    Small {
        molecule: SmallMolecule,
        conformer: ConformerId,
        metadata: MoleculeInstanceMetadata,
        provenance: BuiltMoleculeProvenance,
    },
    Macro {
        molecule: MacroMolecule,
        conformer: ConformerId,
        metadata: MoleculeInstanceMetadata,
        provenance: BuiltMoleculeProvenance,
    },
}

struct BuiltMoleculeProvenance {
    coordinate_model_id: String,
    asym_ids: Vec<String>,
    entity_ids: Vec<String>,
    entity_kinds: Vec<MmcifEntityKind>,
    atoms: Vec<BuiltAtomProvenance>,
}

struct BuiltAtomProvenance {
    atom: AtomId,
    source_line: usize,
    atom_site_id: Option<String>,
    atom_name: String,
    component_id: String,
    asym_id: String,
    auth_asym_id: Option<String>,
    entity_id: Option<String>,
    label_sequence_id: Option<i32>,
    author_sequence_id: Option<String>,
    insertion_code: Option<String>,
    occurrence: Option<usize>,
    selected_alternate_location: Option<String>,
    occupancy: Option<f64>,
    b_factor: Option<f64>,
}

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
pub struct MmcifEnsembleInterpretation {
    ensemble: Ensemble,
    reports: Vec<MmcifInterpretationReport>,
}

impl MmcifEnsembleInterpretation {
    pub fn into_parts(self) -> (Ensemble, Vec<MmcifInterpretationReport>) {
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

/// Interprets explicitly selected or all coordinate models as one
/// shared-topology non-temporal ensemble.
pub fn interpret_mmcif_ensemble(
    document: &MmcifDocument,
    options: MmcifEnsembleInterpretOptions,
) -> Result<MmcifEnsembleInterpretation, MmcifEnsembleInterpretError> {
    if options.model_ids.as_ref().is_some_and(Vec::is_empty) {
        return Err(MmcifEnsembleInterpretError::EmptyModelSelection);
    }
    let blocks = document
        .blocks()
        .iter()
        .filter(|block| block.loop_with_tag("_atom_site.type_symbol").is_some())
        .collect::<Vec<_>>();
    if blocks.is_empty() {
        return Err(MmcifEnsembleInterpretError::NoCoordinateModels);
    }
    if blocks.len() > 1 {
        return Err(MmcifEnsembleInterpretError::MultipleAtomSiteDataBlocks);
    }
    let block = blocks[0];
    let available =
        coordinate_model_ids(block).map_err(|error| MmcifEnsembleInterpretError::Model {
            model_id: "<model inventory>".to_owned(),
            error,
        })?;
    if available.is_empty() {
        return Err(MmcifEnsembleInterpretError::NoCoordinateModels);
    }
    let selected = options.model_ids.unwrap_or_else(|| available.clone());
    let mut seen = BTreeSet::new();
    for model in &selected {
        if !seen.insert(model.clone()) {
            return Err(MmcifEnsembleInterpretError::DuplicateRequestedModel(
                model.clone(),
            ));
        }
        if !available.contains(model) {
            return Err(MmcifEnsembleInterpretError::UnknownRequestedModel(
                model.clone(),
            ));
        }
    }

    let mut interpreted = Vec::with_capacity(selected.len());
    for model_id in &selected {
        let interpretation = interpret_mmcif(
            document,
            MmcifInterpretOptions {
                strict_entity_metadata: options.strict_entity_metadata,
                altloc_policy: options.altloc_policy.clone(),
                model_selection: MmcifModelSelection::Select(model_id.clone()),
            },
        )
        .map_err(|error| MmcifEnsembleInterpretError::Model {
            model_id: model_id.clone(),
            error,
        })?;
        interpreted.push(interpretation);
    }

    let first = interpreted
        .first()
        .ok_or(MmcifEnsembleInterpretError::EmptyModelSelection)?;
    let shared_topology = first.model.shared_topology();
    let shared_atom_identity = provenance_identity(&first.report);
    let mut ensemble = Ensemble::new(Arc::clone(&shared_topology));
    let mut reports = Vec::with_capacity(interpreted.len());
    for interpretation in interpreted {
        let (model, report) = interpretation.into_parts();
        let model_id = report.selected_model().unwrap_or("<unknown>").to_owned();
        let atom_identity = provenance_identity(&report);
        if atom_identity != shared_atom_identity {
            let error = if atom_identity.sorted_atoms() != shared_atom_identity.sorted_atoms() {
                MmcifEnsembleInterpretError::InconsistentAtomSet { model_id }
            } else {
                MmcifEnsembleInterpretError::InconsistentDenseAtomOrder { model_id }
            };
            return Err(error);
        }
        if !shared_topology.same_layout(model.topology()) {
            return Err(MmcifEnsembleInterpretError::InconsistentTopology { model_id });
        }
        if shared_topology.atom_ids() != model.topology().atom_ids()
            || shared_topology.bond_ids() != model.topology().bond_ids()
        {
            return Err(MmcifEnsembleInterpretError::InconsistentDenseAtomOrder { model_id });
        }
        let model_topology = model.shared_topology();
        TopologyMapping::between_identical_layouts(&model_topology, &shared_topology).map_err(
            |_| MmcifEnsembleInterpretError::InconsistentTopology {
                model_id: model_id.clone(),
            },
        )?;
        let positions = Positions::new(&shared_topology, model.positions())
            .map_err(MmcifEnsembleInterpretError::Position)?;
        let mut member = EnsembleMember::new(positions);
        member.set_cell(model.cell().copied());
        let atom_data = AtomData::from_columns(
            &shared_topology,
            model.atom_data().occupancies().map(<[Option<f64>]>::to_vec),
            model.atom_data().b_factors().map(<[Option<f64>]>::to_vec),
        )
        .map_err(|error| MmcifEnsembleInterpretError::Model {
            model_id: model_id.clone(),
            error: graph_error(error),
        })?;
        member
            .set_atom_data(atom_data)
            .map_err(|error| MmcifEnsembleInterpretError::Ensemble(Box::new(error)))?;
        ensemble
            .push(member)
            .map_err(|error| MmcifEnsembleInterpretError::Ensemble(Box::new(error)))?;
        reports.push(report);
    }
    Ok(MmcifEnsembleInterpretation { ensemble, reports })
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct ProvenanceAtomIdentity {
    atom_name: String,
    component_id: String,
    asym_id: String,
    auth_asym_id: Option<String>,
    entity_id: Option<String>,
    label_sequence_id: Option<i32>,
    author_sequence_id: Option<String>,
    insertion_code: Option<String>,
    occurrence: Option<usize>,
    selected_alternate_location: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ProvenanceIdentity {
    atoms: Vec<ProvenanceAtomIdentity>,
}

impl ProvenanceIdentity {
    fn sorted_atoms(&self) -> Vec<ProvenanceAtomIdentity> {
        let mut atoms = self.atoms.clone();
        atoms.sort_unstable();
        atoms
    }
}

type QualifiedAtomData = Vec<(InstanceAtomId, Option<f64>, Option<f64>)>;

fn provenance_identity(report: &MmcifInterpretationReport) -> ProvenanceIdentity {
    ProvenanceIdentity {
        atoms: report
            .instances
            .iter()
            .flat_map(|instance| {
                instance.atoms.iter().map(|atom| ProvenanceAtomIdentity {
                    atom_name: atom.atom_name.clone(),
                    component_id: atom.component_id.clone(),
                    asym_id: atom.asym_id.clone(),
                    auth_asym_id: atom.auth_asym_id.clone(),
                    entity_id: atom.entity_id.clone(),
                    label_sequence_id: atom.label_sequence_id,
                    author_sequence_id: atom.author_sequence_id.clone(),
                    insertion_code: atom.insertion_code.clone(),
                    occurrence: atom.occurrence,
                    selected_alternate_location: atom.selected_alternate_location.clone(),
                })
            })
            .collect(),
    }
}

impl BuiltMoleculeProvenance {
    fn qualify(self, molecule: MoleculeInstanceId) -> (MmcifInstanceProvenance, QualifiedAtomData) {
        let mut atom_data = Vec::with_capacity(self.atoms.len());
        let provenance = MmcifInstanceProvenance {
            molecule,
            coordinate_model_id: self.coordinate_model_id,
            asym_ids: self.asym_ids,
            entity_ids: self.entity_ids,
            entity_kinds: self.entity_kinds,
            atoms: self
                .atoms
                .into_iter()
                .map(|atom| {
                    let qualified = InstanceAtomId::new(molecule, atom.atom);
                    atom_data.push((qualified, atom.occupancy, atom.b_factor));
                    MmcifAtomProvenance {
                        atom: qualified,
                        source_line: atom.source_line,
                        atom_site_id: atom.atom_site_id,
                        atom_name: atom.atom_name,
                        component_id: atom.component_id,
                        asym_id: atom.asym_id,
                        auth_asym_id: atom.auth_asym_id,
                        entity_id: atom.entity_id,
                        label_sequence_id: atom.label_sequence_id,
                        author_sequence_id: atom.author_sequence_id,
                        insertion_code: atom.insertion_code,
                        occurrence: atom.occurrence,
                        selected_alternate_location: atom.selected_alternate_location,
                    }
                })
                .collect(),
        };
        (provenance, atom_data)
    }
}

fn build_molecule(
    group: MoleculeGroup,
    connections: &[DeclaredConnection],
    report: &mut MmcifInterpretationReport,
) -> Result<BuiltMolecule, MmcifInterpretError> {
    let is_macro = group.kinds.iter().any(MmcifEntityKind::is_macro);
    let mut graph = Molecule::new();
    let mut atoms = BTreeMap::new();
    let mut representative = Vec::<(String, AtomRow)>::new();
    let mut seen_atoms = BTreeMap::<String, usize>::new();
    for row in &group.rows {
        let key = row.atom_key();
        if let Some(&index) = seen_atoms.get(&key) {
            let prior = &representative[index].1;
            if prior.element != row.element
                || prior.formal_charge != row.formal_charge
                || prior.comp_id != row.comp_id
                || prior.entity_id != row.entity_id
            {
                return Err(MmcifInterpretError::new(
                    Some(row.line),
                    format!(
                        "atom `{}` has inconsistent topology payload across coordinate models",
                        row.atom_name
                    ),
                ));
            }
        } else {
            seen_atoms.insert(key.clone(), representative.len());
            representative.push((key, row.clone()));
        }
    }
    for (key, row) in &representative {
        let mut atom = Atom::new(row.element);
        atom.formal_charge = row.formal_charge;
        atoms.insert(key.clone(), graph.add_atom(atom).map_err(graph_error)?);
    }
    let model_id = group
        .rows
        .first()
        .map(|row| row.model_id.clone())
        .ok_or_else(|| MmcifInterpretError::new(None, "empty molecule group"))?;
    let mut conformer = Conformer::new(ANGSTROM).expect("angstrom is a length unit");
    for row in &group.rows {
        let point = row.point.ok_or_else(|| {
            MmcifInterpretError::new(
                Some(row.line),
                format!(
                    "missing selected-model position for atom `{}`",
                    row.atom_name
                ),
            )
        })?;
        conformer
            .set_position(atoms[&row.atom_key()], Quantity::new(point, ANGSTROM))
            .expect("matching coordinate units");
    }
    for connection in connections {
        let Some(&left) = atoms.get(&connection.left_atom) else {
            continue;
        };
        let Some(&right) = atoms.get(&connection.right_atom) else {
            continue;
        };
        if graph
            .bond_between(left, right)
            .map_err(graph_error)?
            .is_none()
        {
            graph
                .add_bond(left, right, connection.order)
                .map_err(graph_error)?;
        } else {
            let existing = graph
                .bond_between(left, right)
                .map_err(graph_error)?
                .expect("existing bond was found");
            if graph.bond(existing).map_err(graph_error)?.order != connection.order {
                return Err(MmcifInterpretError::new(
                    None,
                    "duplicate struct_conn records assign conflicting bond orders",
                ));
            }
        }
    }
    let connectivity_candidates = infer_covalent_bonds(&graph, &representative, &atoms)?;
    if connectivity_candidates > 0 {
        report.connectivity_candidates += connectivity_candidates;
        report
            .issues
            .push(MmcifInterpretIssue::ConnectivityCandidatesInferred {
                atom_count: graph.atom_count(),
                candidate_count: connectivity_candidates,
            });
    }
    let asym_ids = representative
        .iter()
        .map(|(_, row)| row)
        .map(|row| row.asym_id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let entity_ids = representative
        .iter()
        .map(|(_, row)| row)
        .filter_map(|row| row.entity_id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    let entity_kinds = group.kinds.iter().cloned().collect::<Vec<_>>();
    let atom_provenance = representative
        .iter()
        .map(|(key, row)| BuiltAtomProvenance {
            atom: atoms[key],
            source_line: row.line,
            atom_site_id: row.atom_site_id.clone(),
            atom_name: row.atom_name.clone(),
            component_id: row.comp_id.clone(),
            asym_id: row.asym_id.clone(),
            auth_asym_id: row.auth_asym_id.clone(),
            entity_id: row.entity_id.clone(),
            label_sequence_id: row.label_seq_id,
            author_sequence_id: row.auth_seq_id.clone(),
            insertion_code: row.insertion_code.clone(),
            occurrence: row.occurrence,
            selected_alternate_location: row.alt_id.clone(),
            occupancy: row.occupancy,
            b_factor: row.b_factor,
        })
        .collect();
    let provenance = BuiltMoleculeProvenance {
        coordinate_model_id: model_id,
        asym_ids,
        entity_ids,
        entity_kinds,
        atoms: atom_provenance,
    };
    if graph.atom_count() > 1 {
        report.template_bonds_pending += 1;
    }

    let mut metadata = MoleculeInstanceMetadata::default();
    for kind in &group.kinds {
        match kind {
            MmcifEntityKind::Polymer => {
                metadata.insert_role(MoleculeRole::Polymer);
            }
            MmcifEntityKind::Branched => {
                metadata.insert_role(MoleculeRole::Branched);
            }
            MmcifEntityKind::NonPolymer => {
                metadata.insert_role(MoleculeRole::NonPolymer);
            }
            MmcifEntityKind::Water => {
                metadata.insert_role(MoleculeRole::Solvent);
            }
            MmcifEntityKind::Other(_) => {}
        }
    }
    if graph.atom_count() == 1
        && graph
            .atoms()
            .next()
            .is_some_and(|(_, atom)| atom.formal_charge != 0)
    {
        metadata.insert_role(MoleculeRole::Ion);
    }
    let conformer = graph
        .add_conformer(conformer)
        .expect("interpreted coordinates reference live atoms");
    if is_macro {
        let hierarchy = build_hierarchy(&graph, &representative, &atoms)?;
        Ok(BuiltMolecule::Macro {
            molecule: MacroMolecule::try_from_parts_unchecked_connectedness(graph, hierarchy)
                .map_err(graph_error)?,
            conformer,
            metadata,
            provenance,
        })
    } else {
        let molecule = SmallMolecule::from_graph_unchecked_connectedness(graph);
        Ok(BuiltMolecule::Small {
            molecule,
            conformer,
            metadata,
            provenance,
        })
    }
}

const COVALENT_BOND_CELL_ANGSTROM: f64 = 2.1;
const COVALENT_BOND_TOLERANCE_ANGSTROM: f64 = 0.45;
const MIN_COVALENT_BOND_DISTANCE_SQUARED: f64 = 0.16;
const FALLBACK_COVALENT_RADIUS_ANGSTROM: f64 = 0.77;

fn infer_covalent_bonds(
    graph: &Molecule,
    representative: &[(String, AtomRow)],
    atoms: &BTreeMap<String, AtomId>,
) -> Result<usize, MmcifInterpretError> {
    let mut cells = BTreeMap::<[i64; 3], Vec<usize>>::new();
    let mut candidates = 0usize;

    for (right_index, (right_key, right_row)) in representative.iter().enumerate() {
        let right_point = right_row
            .point
            .expect("selected mmCIF atom rows have complete positions");
        let right_cell = covalent_bond_cell(right_point, right_row.line)?;
        for offset_x in -1_i64..=1 {
            for offset_y in -1_i64..=1 {
                for offset_z in -1_i64..=1 {
                    let neighbor = covalent_bond_neighbor(
                        right_cell,
                        [offset_x, offset_y, offset_z],
                        right_row.line,
                    )?;
                    let Some(left_indexes) = cells.get(&neighbor) else {
                        continue;
                    };
                    for &left_index in left_indexes {
                        let (left_key, left_row) = &representative[left_index];
                        let left_point = left_row
                            .point
                            .expect("selected mmCIF atom rows have complete positions");
                        let distance_squared = point_distance_squared(left_point, right_point);
                        if distance_squared <= MIN_COVALENT_BOND_DISTANCE_SQUARED {
                            continue;
                        }
                        let left_radius = left_row
                            .element
                            .covalent_radius_angstrom()
                            .unwrap_or(FALLBACK_COVALENT_RADIUS_ANGSTROM);
                        let right_radius = right_row
                            .element
                            .covalent_radius_angstrom()
                            .unwrap_or(FALLBACK_COVALENT_RADIUS_ANGSTROM);
                        let cutoff =
                            (left_radius + right_radius + COVALENT_BOND_TOLERANCE_ANGSTROM)
                                .min(COVALENT_BOND_CELL_ANGSTROM);
                        if distance_squared > cutoff * cutoff {
                            continue;
                        }

                        let left = atoms[left_key];
                        let right = atoms[right_key];
                        if graph
                            .bond_between(left, right)
                            .map_err(graph_error)?
                            .is_none()
                        {
                            candidates += 1;
                        }
                    }
                }
            }
        }
        cells.entry(right_cell).or_default().push(right_index);
    }

    Ok(candidates)
}

fn covalent_bond_cell(point: Point3, line: usize) -> Result<[i64; 3], MmcifInterpretError> {
    Ok([
        covalent_bond_cell_axis(point.x, "x", line)?,
        covalent_bond_cell_axis(point.y, "y", line)?,
        covalent_bond_cell_axis(point.z, "z", line)?,
    ])
}

fn covalent_bond_cell_axis(
    coordinate: f64,
    axis: &str,
    line: usize,
) -> Result<i64, MmcifInterpretError> {
    let cell = (coordinate / COVALENT_BOND_CELL_ANGSTROM).floor();
    if !cell.is_finite() || cell <= i64::MIN as f64 || cell >= i64::MAX as f64 {
        return Err(MmcifInterpretError::new(
            Some(line),
            format!(
                "_atom_site.Cartn_{axis} coordinate is outside the supported covalent-connectivity diagnostic cell range"
            ),
        ));
    }
    Ok(cell as i64)
}

fn covalent_bond_neighbor(
    cell: [i64; 3],
    offset: [i64; 3],
    line: usize,
) -> Result<[i64; 3], MmcifInterpretError> {
    let checked = |axis: usize| {
        cell[axis].checked_add(offset[axis]).ok_or_else(|| {
            MmcifInterpretError::new(
                Some(line),
                "atom-site coordinate exceeds the covalent-connectivity diagnostic neighbor range",
            )
        })
    };
    Ok([checked(0)?, checked(1)?, checked(2)?])
}

fn point_distance_squared(left: Point3, right: Point3) -> f64 {
    let dx = left.x - right.x;
    let dy = left.y - right.y;
    let dz = left.z - right.z;
    dx * dx + dy * dy + dz * dz
}

fn build_hierarchy(
    graph: &Molecule,
    representative: &[(String, AtomRow)],
    atoms: &BTreeMap<String, AtomId>,
) -> Result<SmcraHierarchy, MmcifInterpretError> {
    let mut hierarchy = SmcraHierarchy::new();
    let mut chains = BTreeMap::new();
    let mut residues = BTreeMap::new();
    for (key, row) in representative {
        let chain = if let Some(chain) = chains.get(&row.asym_id) {
            *chain
        } else {
            let chain = hierarchy
                .add_chain(row.asym_id.clone(), row.auth_asym_id.clone())
                .map_err(hierarchy_error)?;
            chains.insert(row.asym_id.clone(), chain);
            chain
        };
        let residue_key = (row.asym_id.clone(), row.residue_key.clone());
        let residue = if let Some(residue) = residues.get(&residue_key) {
            *residue
        } else {
            let residue = hierarchy
                .add_residue(
                    chain,
                    row.comp_id.clone(),
                    row.label_seq_id,
                    row.auth_seq_id.clone(),
                    row.insertion_code.clone(),
                )
                .map_err(hierarchy_error)?;
            let record = &mut hierarchy.residues[residue.index()];
            record.label_comp_id = Some(row.comp_id.clone());
            record.author_comp_id = row.auth_comp_id.clone();
            residues.insert(residue_key, residue);
            residue
        };
        let atom = atoms[key];
        graph.atom(atom).map_err(graph_error)?;
        hierarchy
            .add_atom_site(
                residue,
                atom,
                SmcraAtomSiteMetadata {
                    type_symbol: Some(row.element.symbol().to_owned()),
                    label_asym_id: Some(row.asym_id.clone()),
                    auth_asym_id: row.auth_asym_id.clone(),
                    label_atom_id: Some(row.atom_name.clone()),
                    auth_atom_id: row.auth_atom_name.clone(),
                },
            )
            .map_err(hierarchy_error)?;
    }
    Ok(hierarchy)
}

fn graph_error(error: impl fmt::Display) -> MmcifInterpretError {
    MmcifInterpretError::new(None, error.to_string())
}

fn hierarchy_error(error: impl fmt::Display) -> MmcifInterpretError {
    MmcifInterpretError::new(None, error.to_string())
}

#[derive(Debug)]
struct InstanceUnion {
    parent: BTreeMap<String, String>,
}

impl InstanceUnion {
    fn new(keys: impl IntoIterator<Item = String>) -> Self {
        let parent = keys.into_iter().map(|key| (key.clone(), key)).collect();
        Self { parent }
    }

    fn find(&mut self, key: &str) -> String {
        let mut current = key.to_owned();
        let mut path = Vec::new();
        loop {
            let parent = self
                .parent
                .get(&current)
                .cloned()
                .unwrap_or_else(|| current.clone());
            if parent == current {
                break;
            }
            path.push(current);
            current = parent;
        }
        self.parent
            .entry(current.clone())
            .or_insert_with(|| current.clone());
        for node in path {
            self.parent.insert(node, current.clone());
        }
        current
    }

    fn union(&mut self, left: &str, right: &str) {
        let left_root = self.find(left);
        let right_root = self.find(right);
        if left_root != right_root {
            let (root, child) = if left_root < right_root {
                (left_root, right_root)
            } else {
                (right_root, left_root)
            };
            self.parent.insert(child, root);
        }
    }
}

fn required<'a>(
    table: &'a MmcifLoopTable,
    row: usize,
    tag: &str,
) -> Result<&'a str, MmcifInterpretError> {
    let value = table
        .value(row, tag)
        .ok_or_else(|| row_error(table, row, format!("missing required {tag}")))?;
    value.optional_text().ok_or_else(|| {
        MmcifInterpretError::new(Some(value.line()), format!("missing required {tag}"))
    })
}

fn optional<'a>(table: &'a MmcifLoopTable, row: usize, tag: &str) -> Option<&'a str> {
    table.value(row, tag).and_then(MmcifValue::optional_text)
}

fn optional_f64(
    table: &MmcifLoopTable,
    row: usize,
    tag: &str,
) -> Result<Option<f64>, MmcifInterpretError> {
    optional(table, row, tag)
        .map(|value| {
            let parsed = value
                .parse::<f64>()
                .map_err(|_| row_error(table, row, format!("invalid float {tag}")))?;
            if !parsed.is_finite() {
                return Err(row_error(table, row, format!("non-finite float {tag}")));
            }
            Ok(parsed)
        })
        .transpose()
}

fn optional_i32(
    table: &MmcifLoopTable,
    row: usize,
    tag: &str,
) -> Result<Option<i32>, MmcifInterpretError> {
    optional(table, row, tag)
        .map(|value| {
            value
                .parse::<i32>()
                .map_err(|_| row_error(table, row, format!("invalid integer {tag}")))
        })
        .transpose()
}

fn optional_i8(
    table: &MmcifLoopTable,
    row: usize,
    tag: &str,
) -> Result<Option<i8>, MmcifInterpretError> {
    optional(table, row, tag)
        .map(|value| {
            value
                .parse::<i8>()
                .map_err(|_| row_error(table, row, format!("invalid integer {tag}")))
        })
        .transpose()
}

fn optional_point(
    table: &MmcifLoopTable,
    row: usize,
) -> Result<Option<Point3>, MmcifInterpretError> {
    let x = optional_f64(table, row, "_atom_site.Cartn_x")?;
    let y = optional_f64(table, row, "_atom_site.Cartn_y")?;
    let z = optional_f64(table, row, "_atom_site.Cartn_z")?;
    match (x, y, z) {
        (Some(x), Some(y), Some(z)) => Ok(Some(Point3::new(x, y, z))),
        (None, None, None) => Ok(None),
        _ => Err(row_error(
            table,
            row,
            "partial atom-site coordinate triplet",
        )),
    }
}

fn row_error(
    table: &MmcifLoopTable,
    row: usize,
    message: impl Into<String>,
) -> MmcifInterpretError {
    MmcifInterpretError::new(
        table
            .row(row)
            .and_then(|row| row.first())
            .map(MmcifValue::line),
        message,
    )
}
