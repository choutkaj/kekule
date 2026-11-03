//! Public options, reports, issues, and source provenance of macromolecular
//! structure interpretation.
use std::fmt;
use std::sync::Arc;

use super::altloc::{AltLocDecision, AltLocPolicy, SourceResidueId};
use super::sequences::PolymerEntity;
use crate::components::ComponentLibrary;
use crate::diagnostics::{
    counted, listed, Cancellation, Cancelled, Capitalized, Count, DiagnosticKind, Format,
    FormatError, FormatErrorDetail, SourceLocation,
};
use crate::molecule::Molecule;
use crate::structure::{Ensemble, Model};
use crate::topology::{InstanceAtomId, MoleculeInstanceId, Topology};

/// Which coordinate model of a source one interpretation reads.
///
/// Model IDs are mmCIF `_atom_site.pdbx_PDB_model_num` values or PDB `MODEL`
/// serial numbers as text; a source without them has the one model `"1"`.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum ModelSelection {
    /// The only coordinate model; a source with several fails with
    /// [`DiagnosticKind::Inconsistent`].
    RequireSingle,
    /// The coordinate model with this ID.
    Select(String),
    /// The first coordinate model in source order.
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

/// How interpretation treats structural instances without a declared entity type.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
#[non_exhaustive]
pub enum EntityMetadataPolicy {
    /// Infer the entity kind from the residues and report the inference.
    #[default]
    InferMissing,
    /// Fail interpretation; PDB entries, which declare no entity types,
    /// always fail.
    RequireDeclared,
}

impl EntityMetadataPolicy {
    const fn requires_declared(self) -> bool {
        match self {
            Self::InferMissing => false,
            Self::RequireDeclared => true,
        }
    }
}

/// Policy for interpreting one coordinate model of an mmCIF data block or a
/// PDB entry as a [`Model`].
///
/// The default requires a single coordinate model, keeps the
/// highest-occupancy alternate locations, infers missing entity types, uses
/// the default [`ConnectivityOptions`], and cannot be cancelled.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructureInterpretOptions {
    entity_metadata: EntityMetadataPolicy,
    altloc_policy: AltLocPolicy,
    model_selection: ModelSelection,
    connectivity: ConnectivityOptions,
    cancellation: Option<Cancellation>,
}

impl StructureInterpretOptions {
    /// Whether missing entity types are inferred or rejected.
    pub const fn entity_metadata(&self) -> EntityMetadataPolicy {
        self.entity_metadata
    }

    /// Sets [`Self::entity_metadata`].
    #[must_use]
    pub fn with_entity_metadata(mut self, entity_metadata: EntityMetadataPolicy) -> Self {
        self.entity_metadata = entity_metadata;
        self
    }

    pub(crate) const fn requires_entity_metadata(&self) -> bool {
        self.entity_metadata.requires_declared()
    }

    /// Which alternate locations the selected model keeps;
    /// [`AltLocPolicy::HighestOccupancy`] by default.
    pub fn altloc_policy(&self) -> &AltLocPolicy {
        &self.altloc_policy
    }

    /// Sets [`Self::altloc_policy`].
    #[must_use]
    pub fn with_altloc_policy(mut self, altloc_policy: AltLocPolicy) -> Self {
        self.altloc_policy = altloc_policy;
        self
    }

    /// Which coordinate model is interpreted; [`ModelSelection::RequireSingle`]
    /// by default.
    pub fn model_selection(&self) -> &ModelSelection {
        &self.model_selection
    }

    /// Sets [`Self::model_selection`].
    #[must_use]
    pub fn with_model_selection(mut self, model_selection: ModelSelection) -> Self {
        self.model_selection = model_selection;
        self
    }

    /// How bonds the source does not state are resolved.
    pub fn connectivity(&self) -> &ConnectivityOptions {
        &self.connectivity
    }

    /// Sets [`Self::connectivity`].
    #[must_use]
    pub fn with_connectivity(mut self, connectivity: ConnectivityOptions) -> Self {
        self.connectivity = connectivity;
        self
    }

    /// The cancellation interpretation polls, if any; `None` by default.
    ///
    /// Interpretation polls it between its stages and periodically inside
    /// them, and once it is requested fails with a [`FormatError`] of kind
    /// [`DiagnosticKind::Cancelled`] instead of returning a partial result.
    /// It does not change what a completed interpretation returns.
    pub fn cancellation(&self) -> Option<&Cancellation> {
        self.cancellation.as_ref()
    }

    /// Sets [`Self::cancellation`].
    #[must_use]
    pub fn with_cancellation(mut self, cancellation: Option<Cancellation>) -> Self {
        self.cancellation = cancellation;
        self
    }
}

impl Default for StructureInterpretOptions {
    fn default() -> Self {
        Self {
            entity_metadata: EntityMetadataPolicy::InferMissing,
            altloc_policy: AltLocPolicy::HighestOccupancy,
            model_selection: ModelSelection::RequireSingle,
            connectivity: ConnectivityOptions::default(),
            cancellation: None,
        }
    }
}

/// The kind of an mmCIF entity, the `_entity.type` of a structural instance.
///
/// Matching of source values is case-insensitive.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
#[non_exhaustive]
pub enum EntityKind {
    /// `polymer`: a polypeptide or polynucleotide chain.
    Polymer,
    /// `branched`: a branched oligosaccharide.
    Branched,
    /// `non-polymer`: a ligand, ion, or other small molecule.
    NonPolymer,
    /// `water`.
    Water,
    /// Any other `_entity.type` value, as written in the source.
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

/// Writes the `_entity.type` value, such as `non-polymer`.
impl fmt::Display for EntityKind {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Polymer => "polymer",
            Self::Branched => "branched",
            Self::NonPolymer => "non-polymer",
            Self::Water => "water",
            Self::Other(value) => value,
        })
    }
}

/// Non-fatal scientific interpretation issues retained in the report.
///
/// Connection issues preserve source identity and distinguish unresolved
/// selectors from ambiguous selectors. Neither case creates a topology bond.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StructureInterpretIssue {
    /// A structural instance has no declared entity type, so its kind was
    /// inferred from its residues under [`EntityMetadataPolicy::InferMissing`].
    EntityTypeInferred {
        /// The structural instance (`label_asym_id`, or the PDB chain).
        asym_id: String,
        /// The inferred kind.
        kind: EntityKind,
    },
    /// Atom records without an element symbol had their elements inferred
    /// from the alignment of their atom names, as PDB convention places them.
    ElementsInferredFromAtomNames {
        /// The 1-based source lines of the affected records.
        source_lines: Vec<usize>,
    },
    /// Connection records of a non-covalent type, such as the hydrogen bonds
    /// a file lists, were not interpreted as bonds. Each such type is
    /// reported once.
    ConnectionsIgnored {
        /// The source connection type, such as `hydrog`.
        connection_type: String,
        /// The number of records of this type.
        count: usize,
    },
    /// One declared connection partner matched no selected-model atom.
    ConnectionUnresolved {
        /// The source connection ID, if the record has one.
        connection_id: Option<String>,
        /// The source connection type or record name.
        connection_type: String,
        /// Which partner of the record, `1` or `2`.
        partner: u8,
        /// The 1-based source line of the record, if known.
        source_line: Option<usize>,
        /// Why the partner did not resolve.
        reason: ConnectionResolutionReason,
    },
    /// One declared connection partner matched more than one selected-model atom.
    ConnectionAmbiguous {
        /// The source connection ID, if the record has one.
        connection_id: Option<String>,
        /// The source connection type or record name.
        connection_type: String,
        /// Which partner of the record, `1` or `2`.
        partner: u8,
        /// The 1-based source line of the record, if known.
        source_line: Option<usize>,
        /// The number of matching atoms.
        candidates: usize,
        /// Why the partner is ambiguous.
        reason: ConnectionResolutionReason,
    },
    /// A coordinate model other than the selected one was not interpreted.
    CoordinateModelIgnored {
        /// The ignored model's ID.
        model_id: String,
        /// The number of its atom-site rows.
        atom_site_rows: usize,
    },
    /// One atom-site row of an alternative that alternate-location selection
    /// did not keep.
    AlternateLocationOmitted {
        /// The atom name of the row.
        atom_name: String,
        /// The row's alternate-location label, if any.
        alt_id: Option<String>,
        /// The residue of the row.
        residue: SourceResidueId,
        /// The 1-based source line of the row.
        source_line: usize,
    },
    /// An alternative lacks sites present in another alternative of the same
    /// chemical component. It is excluded rather than completed from other labels.
    AlternateLocationIncomplete {
        /// The residue.
        residue: SourceResidueId,
        /// The excluded alternative.
        alternative: String,
        /// Names of the atoms the alternative lacks.
        missing_atoms: Vec<String>,
    },
    /// Alternate-location rows of a residue state no occupancy.
    AlternateLocationOccupancyMissing {
        /// The residue.
        residue: SourceResidueId,
        /// The 1-based source lines of the rows without occupancy.
        source_lines: Vec<usize>,
    },
    /// A source-declared combination cannot represent this residue (for example,
    /// it omits every label or includes overlapping sites). Other declared
    /// combinations can still be considered; the rejected one is never repaired.
    AlternateLocationCombinationRejected {
        /// The residue.
        residue: SourceResidueId,
        /// The rejected source-declared combination.
        alternative: String,
        /// Why it was rejected.
        reason: String,
    },
    /// Equal heuristic scores, resolved by lexical alternative identifier.
    AlternateLocationOccupancyTie {
        /// The residues selected together.
        residues: Vec<SourceResidueId>,
        /// The tied alternatives.
        alternatives: Vec<String>,
    },
    /// Atom pairs within covalent bonding distance that no evidence bonds,
    /// such as close contacts between residues. They are counted, not bonded.
    ConnectivityCandidatesInferred {
        /// The number of atoms in the molecule searched.
        atom_count: usize,
        /// The number of unbonded close pairs found.
        candidate_count: usize,
    },
    /// A residue without template or bond table had its bonds inferred from
    /// covalent radii; inferred bonds are single.
    ResidueBondsInferred {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The number of inferred bonds.
        bonds: usize,
    },
    /// A hydrogen that no bond source names (the residue template, the
    /// block's own component bond table, or the source's listed bonds) was
    /// bonded to the nearest heavy atom of its residue.
    HydrogenAttachedByDistance {
        /// The residue.
        residue: SourceResidueId,
        /// The hydrogen's atom name.
        atom_name: String,
        /// The atom name of the heavy atom it was bonded to.
        heavy_atom_name: String,
    },
    /// Residues adjacent in sequence, or joined by a branched entity's link,
    /// whose linking atoms lie beyond bonding distance are not bonded.
    ChainBreak {
        /// The first residue.
        residue: SourceResidueId,
        /// The residue that would be linked to it.
        next: SourceResidueId,
    },
    /// A residue that no template gives charges had atoms with more bonds
    /// than a neutral atom allows. Each took the formal charge its bonds
    /// imply, as the methylated N7 of 7-methylguanosine takes +1.
    ChargesFromBonds {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The charged atoms, by atom name, with their formal charge in e.
        charged_atoms: Vec<(String, i8)>,
    },
    /// The model's hydrogens differ from the reference template: a protonation
    /// variant was selected, or heteroatom charges were adjusted by one per
    /// added or removed hydrogen.
    ProtonationFromHydrogens {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The selected protonation variant, if any.
        variant: Option<String>,
        /// The atoms whose formal charge was adjusted, by atom name, with
        /// their resulting formal charge in e.
        charged_atoms: Vec<(String, i8)>,
    },
    /// The resolved atoms of one residue lie in more than one molecule, for
    /// example because a bridging atom is missing from the model.
    ResidueFragmented {
        /// The residue.
        residue: SourceResidueId,
        /// The number of molecules its atoms lie in.
        fragments: usize,
    },
    /// The source states unit-cell parameters (mmCIF `_cell`, PDB `CRYST1`)
    /// that are not the model's periodic cell, so the model has none.
    CellOmitted {
        /// Why the cell was omitted.
        reason: CellOmissionReason,
        /// The source location of the cell parameters.
        location: SourceLocation,
    },
}

/// Writes the issue as one or two sentences, such as `Residue 12 of chain
/// "A" in model "1" (component "7MG") has no template charges. Atom "N7"
/// (+1) takes charges from its bonds.`
impl fmt::Display for StructureInterpretIssue {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EntityTypeInferred { asym_id, kind } => write!(
                formatter,
                "Chain {asym_id:?} has no entity type. Its residues give it type \"{kind}\"."
            ),
            Self::ElementsInferredFromAtomNames { source_lines } => match source_lines.as_slice()
            {
                [line] => write!(
                    formatter,
                    "The atom record at line {line} has no element symbol. Its element comes \
                     from its atom name."
                ),
                lines => write!(
                    formatter,
                    "{}{} have no element symbol. Their elements come from their atom names.",
                    counted(lines.len(), "atom record", "atom records"),
                    FirstLine(lines),
                ),
            },
            Self::ConnectionsIgnored {
                connection_type,
                count,
            } => {
                let (records, verb) = if *count == 1 {
                    ("connection record", "is not read as a bond")
                } else {
                    ("connection records", "are not read as bonds")
                };
                write!(
                    formatter,
                    "{} {records} of type {connection_type:?} {verb}.",
                    Count(*count)
                )
            }
            Self::ConnectionUnresolved {
                connection_id,
                connection_type,
                partner,
                source_line,
                reason,
            } => write!(
                formatter,
                "Partner {partner} of {} is not found. {reason}",
                Connection {
                    id: connection_id.as_deref(),
                    kind: connection_type,
                    line: *source_line,
                }
            ),
            Self::ConnectionAmbiguous {
                connection_id,
                connection_type,
                partner,
                source_line,
                candidates,
                reason,
            } => {
                write!(
                    formatter,
                    "Partner {partner} of {} matches {} atoms, not one.",
                    Connection {
                        id: connection_id.as_deref(),
                        kind: connection_type,
                        line: *source_line,
                    },
                    Count(*candidates)
                )?;
                if *reason != ConnectionResolutionReason::MultipleMatchingAtoms {
                    write!(formatter, " {reason}")?;
                }
                Ok(())
            }
            Self::CoordinateModelIgnored {
                model_id,
                atom_site_rows,
            } => write!(
                formatter,
                "Coordinate model {model_id:?} with {} is not read.",
                counted(*atom_site_rows, "atom site", "atom sites")
            ),
            Self::AlternateLocationOmitted {
                atom_name,
                alt_id,
                residue,
                source_line,
            } => {
                match alt_id {
                    Some(alt_id) => write!(
                        formatter,
                        "Alternate location {alt_id:?} of atom {atom_name:?}"
                    )?,
                    None => write!(
                        formatter,
                        "Atom {atom_name:?} without an alternate location"
                    )?,
                }
                write!(
                    formatter,
                    " at line {source_line} is not read. The atom is in {residue}."
                )
            }
            Self::AlternateLocationIncomplete {
                residue,
                alternative,
                missing_atoms,
            } => write!(
                formatter,
                "Alternate location {alternative:?} of {residue} lacks {} {}. It is not read.",
                if missing_atoms.len() == 1 {
                    "atom"
                } else {
                    "atoms"
                },
                listed(&Quoted::all(missing_atoms))
            ),
            Self::AlternateLocationOccupancyMissing {
                residue,
                source_lines,
            } => write!(
                formatter,
                "{} has {} without occupancy{}.",
                Capitalized(residue),
                counted(
                    source_lines.len(),
                    "alternate-location atom site",
                    "alternate-location atom sites"
                ),
                FirstLine(source_lines),
            ),
            Self::AlternateLocationCombinationRejected {
                residue,
                alternative,
                reason,
            } => write!(
                formatter,
                "The alternate-location combination {alternative:?} of {residue} is not read. {reason}"
            ),
            Self::AlternateLocationOccupancyTie {
                residues,
                alternatives,
            } => write!(
                formatter,
                "Alternate locations {} of {} have equal occupancy. The first by label is kept.",
                listed(&Quoted::all(alternatives)),
                listed(residues)
            ),
            Self::ConnectivityCandidatesInferred {
                atom_count,
                candidate_count,
            } => write!(
                formatter,
                "A molecule of {} has {} within bonding distance.",
                counted(*atom_count, "atom", "atoms"),
                counted(
                    *candidate_count,
                    "unbonded atom pair",
                    "unbonded atom pairs"
                )
            ),
            Self::ResidueBondsInferred {
                residue,
                component_id,
                bonds,
            } => write!(
                formatter,
                "{} (component {component_id:?}) has no template or bond table. Its bonds come \
                 from atom distances: {}.",
                Capitalized(residue),
                counted(*bonds, "single bond", "single bonds")
            ),
            Self::HydrogenAttachedByDistance {
                residue,
                atom_name,
                heavy_atom_name,
            } => write!(
                formatter,
                "Hydrogen {atom_name:?} of {residue} has no stated bond. It is bonded to \
                 {heavy_atom_name:?}, the nearest heavy atom."
            ),
            Self::ChainBreak { residue, next } => write!(
                formatter,
                "{} and {next} are not bonded. Their linking atoms are too far apart.",
                Capitalized(residue)
            ),
            Self::ChargesFromBonds {
                residue,
                component_id,
                charged_atoms,
            } => {
                let (atoms, verb) = if charged_atoms.len() == 1 {
                    ("Atom", "takes charges from its bonds")
                } else {
                    ("Atoms", "take charges from their bonds")
                };
                write!(
                    formatter,
                    "{} (component {component_id:?}) has no template charges. {atoms} {} {verb}.",
                    Capitalized(residue),
                    listed(&Charged::all(charged_atoms))
                )
            }
            Self::ProtonationFromHydrogens {
                residue,
                component_id,
                variant,
                charged_atoms,
            } => {
                write!(
                    formatter,
                    "The hydrogens of {residue} (component {component_id:?}) differ from its \
                     template."
                )?;
                if let Some(variant) = variant {
                    write!(formatter, " It is read as variant {variant:?}.")?;
                }
                if !charged_atoms.is_empty() {
                    let (atoms, verb) = if charged_atoms.len() == 1 {
                        ("Atom", "takes its charge")
                    } else {
                        ("Atoms", "take their charges")
                    };
                    write!(
                        formatter,
                        " {atoms} {} {verb} from the hydrogens.",
                        listed(&Charged::all(charged_atoms))
                    )?;
                }
                Ok(())
            }
            Self::ResidueFragmented { residue, fragments } => write!(
                formatter,
                "The atoms of {residue} lie in {}.",
                counted(*fragments, "molecule", "molecules")
            ),
            Self::CellOmitted { reason, location } => {
                write!(formatter, "The unit cell at {location} is not read. {reason}")
            }
        }
    }
}

/// A connection record as a noun phrase, such as `connection "covale1" of
/// type "covale" at line 40`.
struct Connection<'a> {
    id: Option<&'a str>,
    kind: &'a str,
    line: Option<usize>,
}

impl fmt::Display for Connection<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.id {
            Some(id) => write!(formatter, "connection {id:?} of type {:?}", self.kind)?,
            None => write!(formatter, "a connection of type {:?}", self.kind)?,
        }
        if let Some(line) = self.line {
            write!(formatter, " at line {line}")?;
        }
        Ok(())
    }
}

/// Where several source lines start, such as ` from line 12`, or the one
/// line, ` at line 12`.
struct FirstLine<'a>(&'a [usize]);

impl fmt::Display for FirstLine<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            [] => Ok(()),
            [line] => write!(formatter, " at line {line}"),
            [first, ..] => write!(formatter, " from line {first}"),
        }
    }
}

/// Text from the input in double quotes.
struct Quoted<'a>(&'a str);

impl<'a> Quoted<'a> {
    fn all(texts: &'a [String]) -> Vec<Self> {
        texts.iter().map(|text| Self(text)).collect()
    }
}

impl fmt::Display for Quoted<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{:?}", self.0)
    }
}

/// An atom name with its formal charge, such as `"N7" (+1)`.
struct Charged<'a>(&'a str, i8);

impl<'a> Charged<'a> {
    fn all(atoms: &'a [(String, i8)]) -> Vec<Self> {
        atoms
            .iter()
            .map(|(name, charge)| Self(name, *charge))
            .collect()
    }
}

impl fmt::Display for Charged<'_> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.1 {
            0 => write!(formatter, "{:?} (0)", self.0),
            charge => write!(formatter, "{:?} ({charge:+})", self.0),
        }
    }
}

/// Why stated unit-cell parameters are not a model's periodic cell.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CellOmissionReason {
    /// Some of the six lengths and angles are missing or unknown.
    Incomplete,
    /// The coordinates are the asymmetric unit of this space group, whose
    /// translates alone are not the crystal.
    SpaceGroup {
        /// The space-group symbol as stated.
        name: String,
    },
    /// PDB's 1 Å cube, recorded for structures without a crystal lattice
    /// such as NMR and EM models.
    Placeholder,
}

/// Writes the reason as a sentence, such as `Some cell lengths and angles
/// are missing.`
impl fmt::Display for CellOmissionReason {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Incomplete => formatter.write_str("Some cell lengths and angles are missing."),
            Self::SpaceGroup { name } => write!(
                formatter,
                "The coordinates are the asymmetric unit of space group {name:?}."
            ),
            Self::Placeholder => formatter.write_str("The cell is the 1 Å placeholder cube."),
        }
    }
}

/// The evidence a bond was resolved from.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BondEvidence {
    /// A `_chem_comp_bond` row of the interpreted block.
    ComponentTable,
    /// A component template; `variant` names a protonation variant.
    Template {
        /// The template's component ID.
        component_id: String,
        /// The protonation variant, if one was selected.
        variant: Option<String>,
    },
    /// Standard linkage between polymer residues adjacent in sequence.
    PolymerLink,
    /// A `_pdbx_entity_branch_link` row.
    BranchLink,
    /// A connection record: an mmCIF `_struct_conn` row, whose type
    /// `connection_type` names, or a PDB `LINK`, `SSBOND`, or `CONECT` record,
    /// named by its record type. Metal coordination within one residue
    /// yields a dative bond; between residues it is no bond but a
    /// [`CoordinationProvenance`].
    Connection {
        /// The source connection ID, if the record has one.
        connection_id: Option<String>,
        /// The connection type or record name.
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
    /// The bonded atoms.
    pub const fn atoms(&self) -> [InstanceAtomId; 2] {
        self.atoms
    }

    /// The evidence the bond was resolved from.
    pub fn evidence(&self) -> &BondEvidence {
        &self.evidence
    }
}

/// Metal coordination a connection record states between atoms of different
/// residues, such as an mmCIF `metalc` row or a PDB `LINK` from a metal to a
/// non-metal.
///
/// Coordination between residues is source connection metadata, not a bond:
/// covalent connectedness alone decides molecule boundaries, so a metal ion,
/// its ligands, and the protein that coordinates it stay separate molecules.
/// The mmCIF writer states these records again as `metalc` rows when it is
/// given the report. Coordination within one residue, such as heme's Fe-N
/// bonds, stays a dative bond.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoordinationProvenance {
    pub(crate) atoms: [InstanceAtomId; 2],
    pub(crate) connection_id: Option<String>,
    pub(crate) connection_type: String,
    pub(crate) source_line: Option<usize>,
}

impl CoordinationProvenance {
    /// The two partners, in the order the record names them.
    pub const fn atoms(&self) -> [InstanceAtomId; 2] {
        self.atoms
    }

    /// The source connection ID, if the record has one.
    pub fn connection_id(&self) -> Option<&str> {
        self.connection_id.as_deref()
    }

    /// The source connection type, such as `metalc`, or the record name, such
    /// as `LINK` or `CONECT`.
    pub fn connection_type(&self) -> &str {
        &self.connection_type
    }

    /// The 1-based source line of the record, if known.
    pub const fn source_line(&self) -> Option<usize> {
        self.source_line
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
    UnsupportedSymmetry {
        /// The symmetry operator as stated.
        symmetry: String,
    },
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
    AlternateLocationOmitted {
        /// The named alternate-location label.
        alternate_location: String,
    },
    /// No selected-model atom satisfied every supplied selector field.
    NoMatchingAtom,
    /// More than one selected-model atom satisfied every supplied selector field.
    MultipleMatchingAtoms,
}

impl fmt::Display for ConnectionResolutionReason {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedSymmetry { symmetry } => write!(
                formatter,
                "The partner is symmetry mate {symmetry:?}. Symmetry mates are not supported."
            ),
            Self::MissingSelector => formatter.write_str("The record names no partner atom."),
            Self::ConflictingSelectorValues { selector } => {
                write!(formatter, "The {selector} disagree.")
            }
            Self::ConflictingLabelAndAuthorSelectors => {
                formatter.write_str("The label and author selectors name different atoms.")
            }
            Self::AlternateLocationOmitted { alternate_location } => write!(
                formatter,
                "Alternate location {alternate_location:?} is not in the selected structure."
            ),
            Self::NoMatchingAtom => formatter.write_str("No atom matches the partner selector."),
            Self::MultipleMatchingAtoms => {
                formatter.write_str("Several atoms match the partner selector.")
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
    /// The interpreted atom.
    pub const fn atom(&self) -> InstanceAtomId {
        self.atom
    }

    /// The 1-based source line of the atom's coordinate row or record.
    pub const fn source_line(&self) -> usize {
        self.source_line
    }

    /// `_atom_site.id`, or the PDB atom serial number; `None` when absent.
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

    /// `_atom_site.auth_atom_id`, or the PDB atom name; `None` when absent.
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

    /// `_atom_site.auth_comp_id`, or the PDB residue name; `None` when absent.
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

    /// `_atom_site.label_entity_id`; always `None` for PDB entries, which declare
    /// no entities.
    pub fn entity_id(&self) -> Option<&str> {
        self.entity_id.as_deref()
    }

    /// Returns the source or explicitly inferred entity kind for this atom.
    pub const fn entity_kind(&self) -> &EntityKind {
        &self.entity_kind
    }

    /// `_atom_site.label_seq_id`; always `None` for PDB entries.
    pub const fn label_sequence_id(&self) -> Option<i32> {
        self.label_sequence_id
    }

    /// `_atom_site.auth_seq_id`, or the PDB residue sequence number, as text;
    /// `None` when absent.
    pub fn author_sequence_id(&self) -> Option<&str> {
        self.author_sequence_id.as_deref()
    }

    /// `_atom_site.pdbx_PDB_ins_code`, or the PDB insertion code; `None` when
    /// absent.
    pub fn insertion_code(&self) -> Option<&str> {
        self.insertion_code.as_deref()
    }

    /// Returns the zero-based source occurrence discriminator used when both
    /// label and author sequence identifiers are absent.
    ///
    /// The discriminator is scoped to one coordinate model, asymmetry
    /// identifier, and component identifier. For a PDB entry it instead
    /// numbers the repeats of a residue identity that recurs in a chain, from
    /// `1`, and is `None` for its first occurrence.
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

/// Source correspondence for one interpreted molecule instance: its
/// structural instances, entities, atoms, and bond evidence.
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
    /// The molecule instance.
    pub const fn instance(&self) -> MoleculeInstanceId {
        self.molecule
    }

    /// The ID of the coordinate model the instance was read from.
    pub fn coordinate_model_id(&self) -> &str {
        &self.coordinate_model_id
    }

    /// The structural instances (`label_asym_id`, or PDB chains) the molecule
    /// spans, sorted and distinct.
    pub fn asym_ids(&self) -> &[String] {
        &self.asym_ids
    }

    /// The source entity IDs the molecule spans, sorted and distinct; empty for
    /// PDB entries.
    pub fn entity_ids(&self) -> &[String] {
        &self.entity_ids
    }

    /// The entity kinds of the molecule's residues, sorted and distinct.
    pub fn entity_kinds(&self) -> &[EntityKind] {
        &self.entity_kinds
    }

    /// Source correspondence of every atom of the molecule.
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
    pub(crate) coordination: Vec<CoordinationProvenance>,
    pub(crate) alternate_locations: Vec<AltLocDecision>,
    pub(crate) issues: Vec<StructureInterpretIssue>,
    pub(crate) polymer_entities: Vec<PolymerEntity>,
}

impl StructureInterpretationReport {
    /// The interpreted mmCIF data block's name, or the PDB entry's ID code
    /// from its `HEADER` record (empty without one).
    pub fn block_name(&self) -> &str {
        &self.block_name
    }

    /// The number of `_entity` definitions in the block; `0` for PDB entries.
    pub const fn entity_definitions(&self) -> usize {
        self.entity_definitions
    }

    /// The number of coordinate models in the source.
    pub const fn coordinate_models(&self) -> usize {
        self.coordinate_models
    }

    /// The ID of the interpreted coordinate model.
    pub fn selected_model(&self) -> Option<&str> {
        self.selected_model.as_deref()
    }

    /// The IDs of the source's other coordinate models, in source order.
    pub fn ignored_coordinate_models(&self) -> &[String] {
        &self.ignored_coordinate_models
    }

    /// The number of molecule instances with a polymer or branched entity.
    pub const fn macromolecules(&self) -> usize {
        self.macromolecules
    }

    /// The number of molecule instances that are not macromolecules, water
    /// included.
    pub const fn small_molecules(&self) -> usize {
        self.small_molecules
    }

    /// The number of molecule instances with a water entity.
    pub const fn solvent_molecules(&self) -> usize {
        self.solvent_molecules
    }

    /// Returns the number of uniquely resolved covalent and
    /// metal-coordination connection records.
    pub const fn applied_connections(&self) -> usize {
        self.applied_connections
    }

    /// The number of close atom pairs no evidence bonds, summed over molecules;
    /// see [`StructureInterpretIssue::ConnectivityCandidatesInferred`].
    pub const fn connectivity_candidates(&self) -> usize {
        self.connectivity_candidates
    }

    /// Source correspondence of every molecule instance, in topology instance
    /// order.
    pub fn instances(&self) -> &[InstanceProvenance] {
        &self.instances
    }

    /// Metal coordination the source states between residues, in record
    /// order; see [`CoordinationProvenance`].
    pub fn coordination(&self) -> &[CoordinationProvenance] {
        &self.coordination
    }

    /// Decisions for the selected source coordinate model only. The source
    /// document retains every alternative, including omitted rows and metadata.
    pub fn alternate_locations(&self) -> &[AltLocDecision] {
        &self.alternate_locations
    }

    /// The polymer entities the block states, in the block's order, with
    /// their canonical sequences and the model's residues at each position;
    /// empty for PDB entries, which state no entities. Entities whose
    /// monomers neither `_entity_poly_seq` nor the atom sites give are left
    /// out. [`crate::sequence`] gives the sequences of the chains of any
    /// topology.
    pub fn polymer_entities(&self) -> &[PolymerEntity] {
        &self.polymer_entities
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
    /// The interpreted model.
    pub fn model(&self) -> &Model {
        &self.model
    }

    /// How the model was read from its source.
    pub fn report(&self) -> &StructureInterpretationReport {
        &self.report
    }

    /// The topology of [`Self::model`].
    pub fn topology(&self) -> &crate::topology::Topology {
        self.model.topology()
    }

    /// Borrows molecule definitions once per occurrence in authoritative
    /// topology instance order.
    pub fn molecules(&self) -> impl ExactSizeIterator<Item = &Molecule> + DoubleEndedIterator {
        self.model
            .topology()
            .instances()
            .map(|occurrence| occurrence.molecule())
    }

    /// Clones each molecule definition once per occurrence, in topology instance
    /// order.
    pub fn into_molecules(self) -> Vec<Molecule> {
        self.model
            .topology()
            .instances()
            .map(|occurrence| occurrence.molecule().clone())
            .collect()
    }

    /// Retains shared ownership of the model's exact topology allocation.
    pub fn into_topology(self) -> Arc<Topology> {
        self.model.shared_topology()
    }

    /// The interpreted model, discarding the report.
    pub fn into_model(self) -> Model {
        self.model
    }

    /// The model and its report.
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

/// Interpretation stopped because its cancellation was requested.
pub(in crate::io) fn cancelled(cancelled: Cancelled) -> FormatError {
    FormatError::interpret(Format::Mmcif, DiagnosticKind::Cancelled, "").with_source(cancelled)
}

/// Polls interpretation's cancellation, if any.
#[inline]
pub(in crate::io) fn poll(cancellation: Option<&Cancellation>) -> Result<(), FormatError> {
    Cancellation::poll(cancellation).map_err(cancelled)
}

/// Polls interpretation's cancellation on every few steps of a cheap loop.
#[inline]
pub(in crate::io) fn poll_every(
    cancellation: Option<&Cancellation>,
    step: usize,
) -> Result<(), FormatError> {
    Cancellation::poll_every(cancellation, step).map_err(cancelled)
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
    entity_metadata: EntityMetadataPolicy,
    altloc_policy: AltLocPolicy,
    connectivity: ConnectivityOptions,
    model_ids: Option<Vec<String>>,
    cancellation: Option<Cancellation>,
}

impl StructureEnsembleInterpretOptions {
    /// Whether missing entity types are inferred or rejected.
    pub const fn entity_metadata(&self) -> EntityMetadataPolicy {
        self.entity_metadata
    }

    /// Sets [`Self::entity_metadata`].
    #[must_use]
    pub fn with_entity_metadata(mut self, entity_metadata: EntityMetadataPolicy) -> Self {
        self.entity_metadata = entity_metadata;
        self
    }

    pub(crate) const fn requires_entity_metadata(&self) -> bool {
        self.entity_metadata.requires_declared()
    }

    /// Which alternate locations each member keeps;
    /// [`AltLocPolicy::HighestOccupancy`] by default.
    pub fn altloc_policy(&self) -> &AltLocPolicy {
        &self.altloc_policy
    }

    /// Sets [`Self::altloc_policy`].
    #[must_use]
    pub fn with_altloc_policy(mut self, altloc_policy: AltLocPolicy) -> Self {
        self.altloc_policy = altloc_policy;
        self
    }

    /// How bonds the source does not state are resolved.
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

    /// The cancellation interpretation polls, if any; `None` by default.
    ///
    /// Interpretation polls it between its stages and periodically inside
    /// them, for every member, and once it is requested fails with a
    /// [`FormatError`] of kind [`DiagnosticKind::Cancelled`] instead of
    /// returning a partial ensemble. It does not change what a completed
    /// interpretation returns.
    pub fn cancellation(&self) -> Option<&Cancellation> {
        self.cancellation.as_ref()
    }

    /// Sets [`Self::cancellation`].
    #[must_use]
    pub fn with_cancellation(mut self, cancellation: Option<Cancellation>) -> Self {
        self.cancellation = cancellation;
        self
    }
}

impl Default for StructureEnsembleInterpretOptions {
    fn default() -> Self {
        Self {
            entity_metadata: EntityMetadataPolicy::InferMissing,
            altloc_policy: AltLocPolicy::HighestOccupancy,
            connectivity: ConnectivityOptions::default(),
            model_ids: None,
            cancellation: None,
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
    /// The interpreted ensemble, one member per selected model or conformation.
    pub fn ensemble(&self) -> &Ensemble {
        &self.ensemble
    }

    /// One report per member, in member order.
    pub fn reports(&self) -> &[StructureInterpretationReport] {
        &self.reports
    }

    /// The topology all members share.
    pub fn topology(&self) -> &Topology {
        self.ensemble.topology()
    }

    /// The interpreted ensemble, discarding the reports.
    pub fn into_ensemble(self) -> Ensemble {
        self.ensemble
    }

    /// The ensemble and its per-member reports.
    pub fn into_parts(self) -> (Ensemble, Vec<StructureInterpretationReport>) {
        (self.ensemble, self.reports)
    }
}

#[cfg(test)]
mod tests {
    use super::super::altloc::SourceResiduePosition;
    use super::*;

    fn residue(sequence: i32) -> SourceResidueId {
        SourceResidueId::from_parts("1", "A", SourceResiduePosition::Label(sequence), None)
    }

    fn names(names: &[&str]) -> Vec<String> {
        names.iter().map(|name| (*name).to_owned()).collect()
    }

    #[test]
    fn issues_read_as_sentences() {
        use StructureInterpretIssue as Issue;
        let cases = [
            (
                Issue::EntityTypeInferred {
                    asym_id: "B".into(),
                    kind: EntityKind::NonPolymer,
                },
                r#"Chain "B" has no entity type. Its residues give it type "non-polymer"."#,
            ),
            (
                Issue::ElementsInferredFromAtomNames {
                    source_lines: vec![12],
                },
                "The atom record at line 12 has no element symbol. Its element comes from its \
                 atom name.",
            ),
            (
                Issue::ElementsInferredFromAtomNames {
                    source_lines: vec![12, 13, 1500],
                },
                "3 atom records from line 12 have no element symbol. Their elements come from \
                 their atom names.",
            ),
            (
                Issue::ConnectionsIgnored {
                    connection_type: "hydrog".into(),
                    count: 1_204,
                },
                r#"1,204 connection records of type "hydrog" are not read as bonds."#,
            ),
            (
                Issue::ConnectionsIgnored {
                    connection_type: "hydrog".into(),
                    count: 1,
                },
                r#"1 connection record of type "hydrog" is not read as a bond."#,
            ),
            (
                Issue::ConnectionUnresolved {
                    connection_id: Some("covale1".into()),
                    connection_type: "covale".into(),
                    partner: 2,
                    source_line: Some(40),
                    reason: ConnectionResolutionReason::NoMatchingAtom,
                },
                r#"Partner 2 of connection "covale1" of type "covale" at line 40 is not found. No atom matches the partner selector."#,
            ),
            (
                Issue::ConnectionAmbiguous {
                    connection_id: None,
                    connection_type: "LINK".into(),
                    partner: 1,
                    source_line: None,
                    candidates: 2,
                    reason: ConnectionResolutionReason::MultipleMatchingAtoms,
                },
                r#"Partner 1 of a connection of type "LINK" matches 2 atoms, not one."#,
            ),
            (
                Issue::CoordinateModelIgnored {
                    model_id: "2".into(),
                    atom_site_rows: 1_632,
                },
                r#"Coordinate model "2" with 1,632 atom sites is not read."#,
            ),
            (
                Issue::AlternateLocationOmitted {
                    atom_name: "CG".into(),
                    alt_id: Some("B".into()),
                    residue: residue(12),
                    source_line: 300,
                },
                r#"Alternate location "B" of atom "CG" at line 300 is not read. The atom is in residue 12 of chain "A" in model "1"."#,
            ),
            (
                Issue::AlternateLocationIncomplete {
                    residue: residue(12),
                    alternative: "B".into(),
                    missing_atoms: names(&["CD", "NE"]),
                },
                r#"Alternate location "B" of residue 12 of chain "A" in model "1" lacks atoms "CD" and "NE". It is not read."#,
            ),
            (
                Issue::AlternateLocationOccupancyMissing {
                    residue: residue(12),
                    source_lines: vec![300],
                },
                r#"Residue 12 of chain "A" in model "1" has 1 alternate-location atom site without occupancy at line 300."#,
            ),
            (
                Issue::AlternateLocationCombinationRejected {
                    residue: residue(12),
                    alternative: "1".into(),
                    reason: "The combination has no alternate location of this residue.".into(),
                },
                r#"The alternate-location combination "1" of residue 12 of chain "A" in model "1" is not read. The combination has no alternate location of this residue."#,
            ),
            (
                Issue::AlternateLocationOccupancyTie {
                    residues: vec![residue(12), residue(13)],
                    alternatives: names(&["A", "B"]),
                },
                r#"Alternate locations "A" and "B" of residue 12 of chain "A" in model "1" and residue 13 of chain "A" in model "1" have equal occupancy. The first by label is kept."#,
            ),
            (
                Issue::ConnectivityCandidatesInferred {
                    atom_count: 2_410,
                    candidate_count: 3,
                },
                "A molecule of 2,410 atoms has 3 unbonded atom pairs within bonding distance.",
            ),
            (
                Issue::ResidueBondsInferred {
                    residue: residue(401),
                    component_id: "LIG".into(),
                    bonds: 24,
                },
                r#"Residue 401 of chain "A" in model "1" (component "LIG") has no template or bond table. Its bonds come from atom distances: 24 single bonds."#,
            ),
            (
                Issue::HydrogenAttachedByDistance {
                    residue: residue(401),
                    atom_name: "H12".into(),
                    heavy_atom_name: "C5".into(),
                },
                r#"Hydrogen "H12" of residue 401 of chain "A" in model "1" has no stated bond. It is bonded to "C5", the nearest heavy atom."#,
            ),
            (
                Issue::ChainBreak {
                    residue: residue(12),
                    next: residue(13),
                },
                r#"Residue 12 of chain "A" in model "1" and residue 13 of chain "A" in model "1" are not bonded. Their linking atoms are too far apart."#,
            ),
            (
                Issue::ChargesFromBonds {
                    residue: residue(5),
                    component_id: "7MG".into(),
                    charged_atoms: vec![("N7".into(), 1)],
                },
                r#"Residue 5 of chain "A" in model "1" (component "7MG") has no template charges. Atom "N7" (+1) takes charges from its bonds."#,
            ),
            (
                Issue::ChargesFromBonds {
                    residue: residue(5),
                    component_id: "XYZ".into(),
                    charged_atoms: vec![("N1".into(), 1), ("B2".into(), -1)],
                },
                r#"Residue 5 of chain "A" in model "1" (component "XYZ") has no template charges. Atoms "N1" (+1) and "B2" (-1) take charges from their bonds."#,
            ),
            (
                Issue::ProtonationFromHydrogens {
                    residue: residue(57),
                    component_id: "HIS".into(),
                    variant: Some("HIS_LSN".into()),
                    charged_atoms: Vec::new(),
                },
                r#"The hydrogens of residue 57 of chain "A" in model "1" (component "HIS") differ from its template. It is read as variant "HIS_LSN"."#,
            ),
            (
                Issue::ProtonationFromHydrogens {
                    residue: residue(15),
                    component_id: "LYS".into(),
                    variant: None,
                    charged_atoms: vec![("NZ".into(), 0)],
                },
                r#"The hydrogens of residue 15 of chain "A" in model "1" (component "LYS") differ from its template. Atom "NZ" (0) takes its charge from the hydrogens."#,
            ),
            (
                Issue::ResidueFragmented {
                    residue: residue(7),
                    fragments: 2,
                },
                r#"The atoms of residue 7 of chain "A" in model "1" lie in 2 molecules."#,
            ),
            (
                Issue::CellOmitted {
                    reason: CellOmissionReason::SpaceGroup {
                        name: "P 21 21 21".into(),
                    },
                    location: SourceLocation::at_line(8),
                },
                r#"The unit cell at line 8 is not read. The coordinates are the asymmetric unit of space group "P 21 21 21"."#,
            ),
        ];
        for (issue, text) in cases {
            assert_eq!(issue.to_string(), text, "{issue:?}");
        }
    }
}
