//! Public writer API and orchestration of planning, preparation, and emission.
use crate::diagnostics::{
    counted, Capitalized, DiagnosticKind, Format, FormatError, FormatErrorDetail, GeometryField,
};
use crate::io::mmcif_cell::{cell_parameters, CellParameterError, CellParameters};
use crate::io::structure_interpret::{
    EntityKind, StructureInterpretation, StructureInterpretationReport, StructureSetInterpretation,
};
use crate::io::write_error::{
    check_geometry_state, record_multi_model_structure_set, record_structure_set_properties,
    unrepresentable_state, Omissions, StatePolicies,
};
use crate::io::{UnrepresentablePolicy, WriteOmission, WriteReport};
use crate::math::PeriodicCell;
use crate::molecule::BondOrder;
use crate::structure::{AsStructureView, GeometryView, Structure, StructureSet, StructureView};
use crate::topology::{InstanceAtomId, InstanceBondId, MoleculeClass, MoleculeInstanceId};
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::io::Write;

mod emit;
mod entities;
mod implied;
mod prepare;
use emit::{render_model_to, write_atom_rows, write_block_end, write_block_start};
use entities::{
    entity_plan_from_report, generic_entity_plan, multi_model_entity_plan,
    normalize_entity_classifications,
};
use prepare::{prepare_model, validate_geometry};

const MAX_COORDINATE_PRECISION: usize = 15;

/// Options of the mmCIF writer.
///
/// The default writes a block named `model` with coordinates to 3 decimal
/// places and fails on geometry state, properties, a periodic cell, or fixed
/// implicit-hydrogen counts mmCIF cannot carry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifWriteOptions {
    block_name: String,
    coordinate_precision: usize,
    geometry_metadata: UnrepresentablePolicy,
    properties: UnrepresentablePolicy,
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

    /// The state of geometries that mmCIF cannot carry, whether of a structure set
    /// or a trajectory frame. An mmCIF reader gives one state of equally
    /// weighted geometries without an order, represented by the first geometry,
    /// so by default a set fails to write when a geometry is outside its
    /// distribution, weights are unequal, it has more than one state, a
    /// state energy, state properties, another representative, a
    /// temperature, or an order, a time order without times or steps
    /// included, or a geometry has velocities or forces; so does a structure
    /// block of one geometry with a membership, a place (a frame's time and
    /// MD step among them), velocities, or forces. [`UnrepresentablePolicy::Omit`] writes
    /// the models without them and reports a
    /// [`WriteOmission::GeometryState`] for each field.
    pub const fn geometry_metadata(&self) -> UnrepresentablePolicy {
        self.geometry_metadata
    }

    /// Sets [`Self::geometry_metadata`].
    #[must_use]
    pub fn with_geometry_metadata(mut self, geometry_metadata: UnrepresentablePolicy) -> Self {
        self.geometry_metadata = geometry_metadata;
        self
    }

    /// A structure's per-atom, per-bond, and own properties and the
    /// properties of a structure set or of a geometry's set, which mmCIF cannot
    /// carry: by default a structure or set with any fails to write, and
    /// under [`UnrepresentablePolicy::Omit`] the report names each
    /// ([`WriteOmission::GeometryState`]).
    pub const fn properties(&self) -> UnrepresentablePolicy {
        self.properties
    }

    /// Sets [`Self::properties`].
    #[must_use]
    pub fn with_properties(mut self, properties: UnrepresentablePolicy) -> Self {
        self.properties = properties;
        self
    }

    /// The policies for the state of a structure block.
    const fn state_policies(&self) -> StatePolicies {
        StatePolicies::with_sites(self.geometry_metadata, self.properties)
    }

    /// A periodic cell mmCIF cannot carry: one not periodic along all three
    /// axes, not in the standard crystallographic orientation, or differing
    /// between geometries. By default such a cell fails to write;
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
            geometry_metadata: UnrepresentablePolicy::Fail,
            properties: UnrepresentablePolicy::Fail,
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
    /// Entity semantics name a molecule instance the structure does not have.
    UnknownClassifiedMolecule(MoleculeInstanceId),
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
    /// Interpretation reports carry provenance for an atom the structure does
    /// not have.
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
    /// The number of interpretation reports does not match the number of
    /// structures written.
    ReportCountMismatch {
        /// The number of structures.
        expected: usize,
        /// The number of reports.
        actual: usize,
    },
    /// The 1-based geometry `geometry` differs from the first geometry in state
    /// one mmCIF block shares between models.
    IncompatibleGeometry {
        /// The 1-based geometry number.
        geometry: usize,
        /// The differing state.
        field: MmcifGeometryField,
    },
    /// A periodic cell mmCIF cannot represent, with
    /// [`UnrepresentablePolicy::Fail`].
    UnrepresentableCell {
        /// Why the cell cannot be written.
        reason: MmcifCellReason,
    },
    /// The structure set has no geometries.
    EmptyStructureSet,
    /// The 1-based geometry `geometry` of a structure set has state mmCIF cannot
    /// represent, with [`MmcifWriteOptions::geometry_metadata`] set to fail. A
    /// temperature is reported for the first geometry.
    UnrepresentableGeometry {
        /// The 1-based geometry number.
        geometry: usize,
        /// The state; see [`MmcifWriteOptions::geometry_metadata`].
        field: GeometryField,
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

impl MmcifAtomField {
    /// The verb for the field as the subject: "have" for a plural field,
    /// "has" otherwise.
    const fn verb(self) -> &'static str {
        match self {
            Self::Radical | Self::Hydrogens => "have",
            _ => "has",
        }
    }
}

impl fmt::Display for MmcifAtomField {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Isotope => "isotope",
            Self::Radical => "radical electrons",
            Self::Hydrogens => "implicit hydrogens",
            Self::AtomMap => "atom-map number",
            Self::LabelAsymId => "label_asym_id",
            Self::SourceAsymId => "source asym_id",
            Self::TypeSymbol => "type_symbol",
        })
    }
}

/// Geometry state named by [`MmcifWriteIssue::IncompatibleGeometry`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MmcifGeometryField {
    /// The `_entity` representation.
    Entities,
    /// The `_struct_asym` representation.
    StructuralAsymmetry,
    /// The `_struct_conn` and `_chem_comp_bond` representation.
    Connectivity,
    /// The `_atom_site` identity columns.
    AtomSites,
}

impl fmt::Display for MmcifGeometryField {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Entities => "entity representation",
            Self::StructuralAsymmetry => "structural asymmetry representation",
            Self::Connectivity => "connectivity representation",
            Self::AtomSites => "atom-site representation",
        })
    }
}

/// Why mmCIF cannot represent a periodic cell.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MmcifCellReason {
    /// The geometries have different cells; one block has one `_cell`.
    GeometryCellsDiffer,
    /// The cell is not periodic along all three axes.
    PartiallyPeriodic,
    /// The cell does not have `a` along x and `b` in the xy plane, as mmCIF
    /// coordinates assume.
    NonstandardOrientation,
}

/// Writes the reason as a sentence, such as `The cell is not periodic along
/// all three axes.`
impl fmt::Display for MmcifCellReason {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::GeometryCellsDiffer => "The geometries have different cells.",
            Self::PartiallyPeriodic => "The cell is not periodic along all three axes.",
            Self::NonstandardOrientation => {
                "Cell vector a is not along x, or b is not in the xy plane."
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
            | Self::UnrepresentableGeometry { .. }
            | Self::UnrepresentableCell { .. } => DiagnosticKind::Unsupported,
            Self::MissingEntityClassification(_)
            | Self::DuplicateEntityClassification(_)
            | Self::ConflictingEntityClassifications { .. }
            | Self::ConflictingAsymEntityIds { .. }
            | Self::ConflictingAsymEntityClassifications { .. }
            | Self::ConflictingSourceEntityClassifications { .. }
            | Self::MissingAtomSite(_)
            | Self::DuplicateAtomSite(_)
            | Self::InconsistentAtomSite { .. }
            | Self::DuplicateAtomIdentity(_)
            | Self::MissingAtomProvenance(_)
            | Self::DuplicateAtomProvenance(_)
            | Self::EmptyStructureSet
            | Self::ReportCountMismatch { .. }
            | Self::IncompatibleGeometry { .. } => DiagnosticKind::Inconsistent,
        }
    }
}

impl fmt::Display for MmcifWriteIssue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidBlockName(name) => write!(
                f,
                "Data block name {name:?} is invalid. Use letters, digits, \"_\", \"-\" and \".\"."
            ),
            Self::DuplicateBlockName(name) => {
                write!(f, "Data block name {name:?} is used twice.")
            }
            Self::CoordinatePrecisionTooLarge(precision) => write!(
                f,
                "Coordinate precision {precision} exceeds the maximum of {MAX_COORDINATE_PRECISION} decimal places."
            ),
            Self::MissingEntityClassification(molecule) => write!(
                f,
                "{} has no entity classification.",
                Capitalized(molecule)
            ),
            Self::DuplicateEntityClassification(molecule) => write!(
                f,
                "{} is classified more than once.",
                Capitalized(molecule)
            ),
            Self::ConflictingEntityClassifications {
                molecule,
                classifications,
            } => write!(
                f,
                "The entity of {molecule} covers instances classified as {}.",
                entity_kinds(classifications)
            ),
            Self::UnsupportedEntityClassification {
                molecule,
                classification,
            } => write!(
                f,
                "{} has unsupported entity type {classification:?}.",
                Capitalized(molecule)
            ),
            Self::UnresolvedCanonicalEntityClassification {
                molecule,
                classification,
            } => write!(
                f,
                "The {} class of {molecule} does not determine one entity type.",
                class_name(*classification)
            ),
            Self::ConflictingAsymEntityIds {
                asym_id,
                entity_ids,
            } => write!(
                f,
                "Structural instance {asym_id:?} has several entity IDs: {}.",
                quoted_list(entity_ids)
            ),
            Self::ConflictingAsymEntityClassifications {
                asym_id,
                classifications,
            } => write!(
                f,
                "Structural instance {asym_id:?} has several entity classifications: {}.",
                entity_kinds(classifications)
            ),
            Self::ConflictingSourceEntityClassifications {
                entity_id,
                classifications,
            } => write!(
                f,
                "Source entity {entity_id:?} has several classifications: {}.",
                entity_kinds(classifications)
            ),
            Self::UnknownClassifiedMolecule(molecule) => write!(
                f,
                "The entity classifications name {molecule}. The structure has no such instance."
            ),
            Self::MissingAtomSite(atom) => {
                write!(f, "{} has no atom site.", Capitalized(atom))
            }
            Self::DuplicateAtomSite(atom) => {
                write!(f, "{} is in more than one atom site.", Capitalized(atom))
            }
            Self::InconsistentAtomSite { atom, field } => write!(
                f,
                "The {field} of {atom} differs from the rest of its residue or chain."
            ),
            Self::DuplicateAtomIdentity(atom) => write!(
                f,
                "{} repeats the identity of another atom of its residue.",
                Capitalized(atom)
            ),
            Self::MissingAtomProvenance(atom) => write!(
                f,
                "The interpretation reports have no source record of {atom}."
            ),
            Self::DuplicateAtomProvenance(atom) => write!(
                f,
                "The interpretation reports have two source records of {atom}."
            ),
            Self::UnknownAtomProvenance(atom) => write!(
                f,
                "The interpretation reports have a source record of {atom}. The structure has no such atom."
            ),
            Self::UnsupportedAtomField { atom, field } => write!(
                f,
                "The {field} of {atom} {} no mmCIF representation.",
                field.verb()
            ),
            Self::FormalChargeOutOfRange { atom, charge } => write!(
                f,
                "{} has formal charge {charge}, outside mmCIF's range of -8 to 8.",
                Capitalized(atom)
            ),
            Self::UnsupportedStereo(molecule) => write!(
                f,
                "The stereochemistry of {molecule} has no mmCIF representation."
            ),
            Self::UnsupportedBondOrder { bond, order } => write!(
                f,
                "The {} order of {bond} has no mmCIF representation.",
                bond_order_name(*order)
            ),
            Self::AmbiguousConnectionSelector(atom) => write!(
                f,
                "{} has no unique _struct_conn selector.",
                Capitalized(atom)
            ),
            Self::UnsupportedTextValue { item } => {
                write!(f, "A value of {item} does not fit one mmCIF token.")
            }
            Self::ReportCountMismatch { expected, actual } => write!(
                f,
                "The block has {} and {}.",
                counted(*expected, "structure", "structures"),
                counted(*actual, "interpretation report", "interpretation reports")
            ),
            Self::IncompatibleGeometry { geometry, field } => write!(
                f,
                "Geometry {geometry} has a different {field} from the first geometry."
            ),
            Self::UnrepresentableCell { reason } => reason.fmt(f),
            Self::EmptyStructureSet => f.write_str("The structure set has no geometries."),
            Self::UnrepresentableGeometry { geometry, field } => {
                f.write_str(&unrepresentable_state("mmCIF", *field, Some(*geometry)))
            }
        }
    }
}

/// Names entity kinds in a sentence, such as `polymer and water`.
fn entity_kinds(kinds: &[EntityKind]) -> String {
    let names = kinds
        .iter()
        .map(|kind| match kind {
            EntityKind::Polymer => "polymer".to_owned(),
            EntityKind::Branched => "branched".to_owned(),
            EntityKind::NonPolymer => "non-polymer".to_owned(),
            EntityKind::Water => "water".to_owned(),
            EntityKind::Other(name) => format!("{name:?}"),
        })
        .collect::<Vec<_>>();
    and_list(&names)
}

/// Quotes values in a sentence, such as `"1" and "2"`.
fn quoted_list(values: &[String]) -> String {
    and_list(
        &values
            .iter()
            .map(|value| format!("{value:?}"))
            .collect::<Vec<_>>(),
    )
}

/// Joins items as `a`, `a and b`, or `a, b and c`.
fn and_list(items: &[String]) -> String {
    match items {
        [] => String::new(),
        [only] => only.clone(),
        [rest @ .., last] => format!("{} and {last}", rest.join(", ")),
    }
}

/// The name of a molecule class in a sentence.
fn class_name(class: MoleculeClass) -> &'static str {
    match class {
        MoleculeClass::Protein => "protein",
        MoleculeClass::Dna => "DNA",
        MoleculeClass::Rna => "RNA",
        MoleculeClass::Carbohydrate => "carbohydrate",
        MoleculeClass::Water => "water",
        MoleculeClass::Ion => "ion",
        MoleculeClass::SmallMolecule => "small-molecule",
        MoleculeClass::Other => "other",
    }
}

/// The name of a bond order in a sentence, such as `double` in `the double
/// order of bond 3`.
pub(super) fn bond_order_name(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Zero => "zero",
        BondOrder::Single => "single",
        BondOrder::Double => "double",
        BondOrder::Triple => "triple",
        BondOrder::Quadruple => "quadruple",
        BondOrder::Dative => "dative",
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

/// The source structure is internally inconsistent or incomplete.
pub(super) fn invalid_structure(message: impl Into<String>) -> FormatError {
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

/// One mmCIF data block to write: a structure or a structure set, plus the
/// entity semantics used for its `_entity` and asymmetry assignments.
///
/// Convert from `&Structure`, `&StructureSet`, [`GeometryView`],
/// `&StructureInterpretation`, or `&StructureSetInterpretation`;
/// interpretations keep their source reports. [`Self::structure`] accepts any
/// borrowed structure view, such as a geometry or a trajectory frame; a
/// geometry's place, velocities, and forces fail to write unless
/// the write options omit them.
///
/// A block is named by [`MmcifWriteOptions::block_name`] unless
/// [`Self::with_block_name`] names it.
///
/// Entity kinds follow canonical classification; hierarchy alone does not
/// establish polymer status. The writer states only the bonds its reader
/// cannot derive, and states a residue's bonds whole or not at all, since one
/// stated internal bond makes the reader take all of that residue's bonds
/// from the source.
///
/// Every hierarchy chain is a structural instance of its own, even when
/// chains share a label, as the copies
/// [`StructureEditor::append_structure`](crate::structure::StructureEditor::append_structure)
/// appends do. The first chain with a label keeps it as `label_asym_id`; each
/// later one is written as the label with a copy serial, `A-2`, `A-3`, and so
/// on, skipping IDs another chain or asym uses, while `auth_asym_id` keeps
/// the author chain name, which falls back to the label. Reading the block
/// back gives the same chains, residues, atoms, and bonds, in the same order
/// and under the same author names, with the written asym IDs as their
/// labels. The hierarchy itself keeps labels as appended, since a label
/// names a source's instance, not one of the combined structure's.
#[derive(Debug, Clone, Copy)]
pub struct MmcifBlockSource<'a> {
    content: BlockContent<'a>,
    entities: EntitySemantics<'a>,
    block_name: Option<&'a str>,
}

#[derive(Debug, Clone, Copy)]
enum BlockContent<'a> {
    Structure(StructureView<'a>),
    StructureSet(&'a StructureSet),
}

#[derive(Debug, Clone, Copy)]
enum EntitySemantics<'a> {
    Canonical,
    Classifications(&'a MmcifEntityClassifications),
    Reports(&'a [StructureInterpretationReport]),
}

impl<'a> MmcifBlockSource<'a> {
    /// One structure block with entity kinds derived from canonical topology
    /// classification.
    pub fn structure(structure: &'a (impl AsStructureView + ?Sized)) -> Self {
        Self {
            content: BlockContent::Structure(structure.as_structure_view()),
            entities: EntitySemantics::Canonical,
            block_name: None,
        }
    }

    /// One multi-model block with one model per geometry of a structure set, in
    /// geometry order. The state mmCIF cannot carry follows
    /// [`MmcifWriteOptions::geometry_metadata`].
    pub fn structure_set(set: &'a StructureSet) -> Self {
        Self {
            content: BlockContent::StructureSet(set),
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
    /// for a structure block, one per geometry for a structure set block.
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

impl<'a> From<&'a Structure> for MmcifBlockSource<'a> {
    fn from(structure: &'a Structure) -> Self {
        Self::structure(structure)
    }
}

impl<'a> From<&'a StructureSet> for MmcifBlockSource<'a> {
    fn from(set: &'a StructureSet) -> Self {
        Self::structure_set(set)
    }
}

impl<'a, R: ?Sized> From<GeometryView<'a, R>> for MmcifBlockSource<'a> {
    fn from(geometry: GeometryView<'a, R>) -> Self {
        Self {
            content: BlockContent::Structure(geometry.as_structure_view()),
            entities: EntitySemantics::Canonical,
            block_name: None,
        }
    }
}

impl<'a> From<&'a StructureSetInterpretation> for MmcifBlockSource<'a> {
    fn from(interpretation: &'a StructureSetInterpretation) -> Self {
        Self::structure_set(interpretation.structure_set()).with_reports(interpretation.reports())
    }
}

impl<'a> From<&'a StructureInterpretation> for MmcifBlockSource<'a> {
    fn from(interpretation: &'a StructureInterpretation) -> Self {
        Self::structure(interpretation.structure())
            .with_reports(std::slice::from_ref(interpretation.report()))
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
/// follow dense atom order, and `_struct_conn` states the bonds the reader
/// does not derive (see the `implied` module).
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
            BlockContent::Structure(view) => {
                write_model_block_to(writer, view, block.entities, &block_options, &mut omissions)?
            }
            BlockContent::StructureSet(set) => {
                let views = set
                    .iter()
                    .map(|geometry| geometry.as_structure_view())
                    .collect::<Vec<_>>();
                let record_state = |state: &mut Omissions| {
                    let unrepresentable = |(geometry, field)| {
                        fail(MmcifWriteIssue::UnrepresentableGeometry { geometry, field })
                    };
                    record_multi_model_structure_set(set, options.geometry_metadata, state)
                        .map_err(unrepresentable)?;
                    record_structure_set_properties(set, options.properties, state)
                        .map_err(unrepresentable)
                };
                write_multi_model_block_to(
                    writer,
                    &views,
                    record_state,
                    MmcifWriteIssue::EmptyStructureSet,
                    block.entities,
                    &block_options,
                    &mut omissions,
                )?;
            }
        }
    }
    Ok(omissions.into_report())
}

fn write_model_block_to(
    writer: &mut impl Write,
    view: StructureView<'_>,
    entities: EntitySemantics<'_>,
    options: &MmcifWriteOptions,
    omissions: &mut Omissions,
) -> Result<(), FormatError> {
    check_geometry_state(Format::Mmcif, view, options.state_policies())?;
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
    let prepared = prepare_model(
        view,
        plan,
        options.implicit_hydrogens,
        options.coordinate_precision,
    )?;
    render_model_to(writer, &prepared, cell, options)?;
    // Nothing is reported for a block that did not reach the sink.
    omissions.record_geometry_state(view, options.state_policies());
    if cell.is_none() && view.cell().is_some() {
        omissions.record(WriteOmission::Cell { structures: 1 });
    }
    omissions.record(WriteOmission::ImplicitHydrogens {
        atoms: prepared.omitted_hydrogens,
    });
    omissions.record(WriteOmission::Coordination {
        records: prepared.omitted_coordination,
    });
    omissions.record(WriteOmission::CoordinationBonds {
        bonds: prepared.coordination_bonds,
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
        Err(MmcifCellReason::GeometryCellsDiffer)
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

/// Writes one block of `views` as models `1..`. It checks the report count,
/// then the multi-model state the block omits, which `record_state` records
/// or rejects, then the cell, and fails with `empty` without views. The
/// state's omissions are recorded once the block has reached the sink.
fn write_multi_model_block_to(
    writer: &mut impl Write,
    views: &[StructureView<'_>],
    record_state: impl FnOnce(&mut Omissions) -> Result<(), FormatError>,
    empty: MmcifWriteIssue,
    entities: EntitySemantics<'_>,
    options: &MmcifWriteOptions,
    omissions: &mut Omissions,
) -> Result<(), FormatError> {
    let (classifications, reports) = match entities {
        EntitySemantics::Canonical => (None, None),
        EntitySemantics::Classifications(classifications) => (Some(classifications), None),
        EntitySemantics::Reports(reports) => {
            if reports.len() != views.len() {
                return Err(fail(MmcifWriteIssue::ReportCountMismatch {
                    expected: views.len(),
                    actual: reports.len(),
                }));
            }
            (None, Some(reports))
        }
    };
    let mut state = Omissions::default();
    record_state(&mut state)?;
    let cell = block_cell(views.iter().map(|view| view.cell()), options.cell)?;
    let mut geometries = views.iter().copied().enumerate();
    let (_, first_view) = geometries.next().ok_or_else(|| fail(empty))?;
    let first_plan =
        multi_model_entity_plan(first_view, classifications, reports.map(|all| &all[0]))?;
    let mut first = prepare_model(
        first_view,
        first_plan,
        options.implicit_hydrogens,
        options.coordinate_precision,
    )?;

    write_block_start(writer, &first, cell, options)?;
    let mut atom_serial = 1u64;
    write_atom_rows(writer, &first.atoms, 1, &mut atom_serial, options)?;

    for (index, view) in geometries {
        let plan = multi_model_entity_plan(view, classifications, reports.map(|all| &all[index]))?;
        let candidate = prepare_model(
            view,
            plan,
            options.implicit_hydrogens,
            options.coordinate_precision,
        )?;
        validate_geometry(&first, &candidate, index + 1)?;
        // A bond the reader derives in one model but not another is stated.
        for (connection, geometry) in first.connections.iter_mut().zip(&candidate.connections) {
            if connection.implied != geometry.implied {
                connection.implied = None;
            }
        }
        write_atom_rows(
            writer,
            &candidate.atoms,
            index + 1,
            &mut atom_serial,
            options,
        )?;
    }
    write_block_end(writer, &first)?;
    for omission in state.into_report().omissions() {
        omissions.record(*omission);
    }
    if cell.is_none() {
        omissions.record(WriteOmission::Cell {
            structures: views.iter().filter(|view| view.cell().is_some()).count(),
        });
    }
    // Geometries share one topology, so the block omits each atom's count once.
    omissions.record(WriteOmission::ImplicitHydrogens {
        atoms: first.omitted_hydrogens,
    });
    omissions.record(WriteOmission::Coordination {
        records: first.omitted_coordination,
    });
    omissions.record(WriteOmission::CoordinationBonds {
        bonds: first.coordination_bonds,
    });
    Ok(())
}

fn one_based_serial(raw: u32) -> u64 {
    u64::from(raw) + 1
}
