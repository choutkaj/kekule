//! Public options, reports, issues, and source provenance of macromolecular
//! structure interpretation.
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::sync::Arc;

use super::altloc::{AltLocDecision, AltLocPolicy, SourceResidueId};
use super::assemblies::Assembly;
use super::sequences::PolymerEntity;
use crate::components::ComponentLibrary;
use crate::diagnostics::{
    counted, listed, Cancellation, Cancelled, Capitalized, Count, DiagnosticKind, Format,
    FormatError, FormatErrorDetail, SourceLocation,
};
use crate::molecule::{BondOrder, Molecule};
use crate::structure::{Structure, StructureSet};
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
/// [`ComponentLibrary::bundled`] when the residue's atoms fit it, then the
/// bundled template under CHARMM atom names, such as an isoleucine with `CD`,
/// then the [`ComponentLibrary::amber`] or [`ComponentLibrary::charmm`]
/// template of its name, such as AMBER `CYX`: the first that names all its
/// heavy atoms. A residue that none covers takes the bonds connection
/// records list within it, and otherwise has its bonds inferred from
/// covalent radii when
/// [`infer_bonds_from_distance`](Self::infer_bonds_from_distance) is set, and
/// fails otherwise. Listed bonds keep the orders their records state; PDB
/// `CONECT` records state orders only in a bonded group of atoms that lists
/// some bond more than once, and the report names residues whose bond orders
/// no record states.
///
/// A template always applies to a residue whose atoms fit it. Connection
/// records within the residue are checked against it: a bond the template
/// lacks is kept, and an order a record states is kept, and the report
/// names both. Stated orders outrank hydrogens: a residue takes only a
/// protonation variant that gives every stated order. When none gives them
/// with the residue's hydrogens, its charges follow from both, and the
/// report names it. A residue with a heavy atom that its
/// template does not list, or lists with another element, takes the bonds
/// its records state when they bond each of its heavy atoms, and the report
/// names it; otherwise it is an error, never a silently separate molecule.
///
/// Formal charges come from the source when stated, otherwise from the
/// template, also when the source's own bond table supplies the bonds.
/// Hydrogens present in the model are atoms, and the remaining hydrogens
/// follow from valence, so a template's protonation is never asserted as a
/// fixed hydrogen count. When a residue's model hydrogens on a heteroatom
/// differ from its template, a variant template that matches them exactly is
/// selected; otherwise each added or removed hydrogen on a heteroatom whose
/// template bonds are all single, or on a ring nitrogen, changes its charge by
/// one, and any other difference fails interpretation. Hydrogens absent from
/// carbon are left to valence.
///
/// Hydrogen bonds, salt bridges, and links to symmetry mates are never bonds.
/// Formats without residues (XYZ, Molfile, SDF, SMILES) have no
/// inter-residue coordination, so their dative bonds are bonds.
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
    /// covalent radii, and disulfides between cysteine `SG` atoms within
    /// bonding distance that carry no hydrogen and no stated link. Inferred
    /// bonds are single and reported. On by default.
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
/// PDB entry as a [`Structure`].
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
    /// A residue whose atoms fit no dictionary template, or whose name none
    /// defines, was read with a force-field template of that name, which
    /// names all its heavy atoms: the dictionary template under CHARMM atom
    /// names, with the dictionary's hydrogens and charges, as for
    /// isoleucine's `CD`, else the AMBER or CHARMM template, as AMBER `CYX`
    /// names a cysteine in a disulfide.
    ForceFieldTemplate {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The force field whose template it was read with: `AMBER` or
        /// `CHARMM`.
        force_field: String,
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
        /// Onium centres, by atom name, that took no charge because their
        /// bond orders are unknown, with the hydrogens each took to its next
        /// allowed valence. Each may instead carry a charge or a double bond.
        neutral_centres: Vec<(String, u8)>,
    },
    /// A residue that no template or bond table covers took its bonds from
    /// connection records that state no bond orders, such as PDB `CONECT`
    /// records in a group that lists no bond twice, or `LINK` records. Its
    /// bonds are single whatever its chemistry, as for
    /// [`Self::ResidueBondsInferred`]. A PDB ligand with only single bonds
    /// is reported too: `CONECT` records cannot state a single bond.
    ResidueBondOrdersUnstated {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The number of its bonds between heavy atoms without a stated
        /// order.
        bonds: usize,
        /// Onium centres, by atom name, that took no charge because their
        /// bond orders are unknown, with the hydrogens each took to its next
        /// allowed valence. Each may instead carry a charge or a double bond.
        neutral_centres: Vec<(String, u8)>,
    },
    /// The atoms of a residue, by name and element, fit no template of its
    /// name, and connection records bond each of its heavy atoms. It takes
    /// the bonds and orders they state, as a residue without a template
    /// does, and its template gives it neither bonds nor charges. A residue
    /// whose records leave a heavy atom out fails to read.
    ResidueFitsNoTemplate {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
    },
    /// Connection records state bond orders within a residue that no
    /// protonation variant gives together with the residue's hydrogens. The
    /// stated orders are kept, and its charges follow from them and its
    /// hydrogens. A histidine ring stated as NE2-protonated with a hydrogen
    /// on ND1 reads as zwitterionic.
    StatedOrdersDisagreeWithHydrogens {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
    },
    /// A connection record bonds two atoms of a residue that its template
    /// does not bond. The bond is kept, with the record's order.
    BondNotInTemplate {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The atom names of the two atoms, in record order.
        atoms: [String; 2],
    },
    /// A connection record states an order for a bond within a residue that
    /// its template gives another order. The stated order is kept.
    BondOrderDiffersFromTemplate {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The atom names of the two atoms, in record order.
        atoms: [String; 2],
        /// The stated order.
        order: BondOrder,
        /// The template's order, or `None` for an aromatic or delocalized
        /// template bond.
        template_order: Option<BondOrder>,
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
    /// Two cysteines whose `SG` atoms lie within bonding distance, without
    /// hydrogens or a record linking them, were bonded as a disulfide.
    DisulfideInferred {
        /// The two cysteines.
        residues: [SourceResidueId; 2],
    },
    /// Residues adjacent in sequence, or joined by a branched entity's link,
    /// whose linking atoms lie beyond bonding distance are not bonded. Their
    /// polymer linking atoms are reported as [`Self::PolymerLinkMissing`].
    ChainBreak {
        /// The first residue.
        residue: SourceResidueId,
        /// The residue that would be linked to it.
        next: SourceResidueId,
    },
    /// Residues adjacent in sequence whose standard polymer link is not
    /// made, because a connection record bonds one of its linking atoms to
    /// another atom, as a beta-aspartyl isopeptide bonds the next residue's
    /// `N` to `CG`.
    PolymerLinkSkipped {
        /// The first residue.
        residue: SourceResidueId,
        /// The residue that would be linked to it.
        next: SourceResidueId,
    },
    /// Residues adjacent in sequence, one of which lacks an atom their
    /// polymer link joins, such as a residue modelled without its `C`, are
    /// not bonded.
    LinkingAtomMissing {
        /// The first residue.
        residue: SourceResidueId,
        /// The residue that would be linked to it.
        next: SourceResidueId,
    },
    /// A residue that no template gives charges had atoms with more bonds
    /// than a neutral atom allows, or onium centres that a neutral atom
    /// completes only with a hydrogen above its default valence. Each took
    /// the formal charge its bonds imply, as the methylated N7 of
    /// 7-methylguanosine and the sulfur of trimethylsulfonium take +1. A
    /// singly bonded terminal chalcogen or boron neighbor of an onium centre
    /// took -1, as the oxygen of DMSO drawn with a single bond.
    ChargesFromBonds {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The charged atoms, by atom name, with their formal charge in e.
        charged_atoms: Vec<(String, i8)>,
    },
    /// An alkali, alkaline-earth, or p-block metal atom that no template
    /// gives a charge has a valence its element and charge do not allow, as
    /// a tin atom with three bonds or an unbonded neutral sodium atom. It
    /// takes no hydrides; the source likely omits a charge, a bond, or a
    /// hydrogen. Supplying the component as context gives it its charge.
    MetalValenceIncomplete {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The atom's name.
        atom_name: String,
        /// Its valence: bond orders and radical electrons.
        valence: u8,
    },
    /// Non-metal atoms of a residue that no template gives charges bond
    /// metals at a bond valence their element and charge do not allow, but
    /// one they allow when some of those bonds are not counted, as the
    /// phosphorus of a phosphine on palladium or a sulfide bridging three iron
    /// atoms. Those bonds are donor bonds stated as covalent ones. The atoms
    /// keep their charge and take no hydrogens.
    MetalDonors {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The donor atoms, by atom name.
        atom_names: Vec<String>,
    },
    /// The structure's hydrogens differ from the reference template: a
    /// protonation variant was selected, or heteroatom charges were adjusted by
    /// one per added or removed hydrogen. A cationic nitrogen with its
    /// neutral count of hydrogens, as lysine's `NZ` with two, is the neutral
    /// base.
    ProtonationFromHydrogens {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The selected protonation variant, if any.
        variant: Option<String>,
        /// The atoms whose formal charge the hydrogens set to a nonzero
        /// value, by atom name, with that charge in e.
        charged_atoms: Vec<(String, i8)>,
    },
    /// Residues without hydrogens on their heteroatoms took an assumed
    /// protonation. Strong acids take their state at pH 7 in every residue
    /// a dictionary template covers: carboxylic, phosphoric, phosphonic,
    /// and sulfonic acids and sulfates are deprotonated. Polymer residues
    /// also take the ammonium of a true N terminus and their least charged
    /// variant, as neutral histidine with its hydrogen on `NE2`. Ligands
    /// keep the template's state elsewhere, amines included. Reported once
    /// per structure.
    ProtonationAssumed {
        /// Component IDs in order, each with its number of residues.
        components: Vec<(String, usize)>,
    },
    /// Heteroatoms of a residue with model hydrogens lack hydrogens that no
    /// charge explains, as an amine, a phenol, a thiol, or a hydroxyl group
    /// without its hydrogens. They keep their template charge and take
    /// hydrogens by valence. Strong acids that lack their hydrogen are
    /// deprotonated instead, and a cationic nitrogen that keeps its neutral
    /// count is the neutral base.
    HydrogensMissing {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The heavy atoms lacking hydrogens, by atom name.
        atoms: Vec<String>,
    },
    /// Hydrogens of a residue whose source occupancy is zero, while its
    /// heavy atoms have occupancy, as placed but unobserved hydrogens. They
    /// are not read, so they decide no protonation.
    ZeroOccupancyHydrogens {
        /// The residue.
        residue: SourceResidueId,
        /// The hydrogens, by atom name.
        atoms: Vec<String>,
    },
    /// Template heteroatoms that coordinate a metal in another residue, as
    /// zinc-bound cysteine does, have one hydrogen fewer than the template:
    /// a thiolate, a carboxylate, or an imidazole nitrogen without its
    /// hydrogen. Alcohols and water keep their hydrogens.
    HydrogenReplacedByCoordination {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The coordinating atoms, by atom name.
        atoms: Vec<String>,
    },
    /// An atom's stated formal charge differs from the one its template and
    /// hydrogens give, and some hydrogen count holds it, so it is kept with
    /// that count. In a residue without hydrogens the charge sets the count,
    /// as a neutral aspartate `OD2` takes one.
    StatedChargeKept {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The atom name.
        atom_name: String,
        /// The stated charge in e.
        stated: i8,
        /// The charge the template and hydrogens give, in e.
        derived: i8,
    },
    /// An atom's stated formal charge needs fewer hydrogens than none, or
    /// fewer than the model gives it, as -2 on a hydroxyl oxygen does, so
    /// the charge its template and hydrogens give replaces it.
    StatedChargeRejected {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The atom name.
        atom_name: String,
        /// The stated charge in e.
        stated: i8,
        /// The charge kept, in e.
        kept: i8,
    },
    /// Stated charges give another resonance form of their groups than the
    /// template, as `OD1` -1 of an aspartate or `NH1` +1 of an arginine
    /// do. The double bonds move to that form, and the charges are kept.
    StatedChargeResonance {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The atoms whose charge the form moves, by atom name, with their
        /// formal charge in e.
        charged_atoms: Vec<(String, i8)>,
    },
    /// A polymer residue's linking atom, a peptide's `N` or `C` or a
    /// nucleotide's `P` or `O3'`, has no bond to its neighbour in the
    /// chain: the model lacks that residue, as before the first modelled
    /// residue of an entity that its sequence continues, or places it
    /// beyond bonding distance. The atom keeps the hydrogens it has when
    /// bonded, so an `N` after a gap is an amide with one hydrogen, not an
    /// ammonium.
    PolymerLinkMissing {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The linking atom's name.
        atom_name: String,
    },
    /// A residue's template states no charge for some atoms, as the
    /// dictionary does for the iron of heme and of iron-sulfur clusters.
    /// They are read as neutral.
    TemplateChargeUnknown {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The atoms without a template charge, by atom name.
        atoms: Vec<String>,
    },
    /// A residue of one atom, or one holding a metal, has no template, and
    /// neither the source nor its component tables state its charges. Its
    /// atoms are read as neutral and take hydrogens by valence.
    TemplateMissing {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
    },
    /// Heavy atoms of a residue's template that the residue lacks, other
    /// than atoms that leave when it links. Their bonded atoms take
    /// hydrogens by valence in their place.
    AtomsMissing {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The missing atoms, by atom name.
        atoms: Vec<String>,
    },
    /// A terminal atom that leaves when the residue links, such as a
    /// carboxyl `OXT` or a phosphate `OP3`, is missing although nothing
    /// links in its place and the chain ends there. Its bonded atom takes no
    /// hydrogen for it, so a carboxyl terminus without `OXT` does not become
    /// an aldehyde.
    TerminalAtomMissing {
        /// The residue.
        residue: SourceResidueId,
        /// Its component ID.
        component_id: String,
        /// The missing atom's name.
        atom_name: String,
        /// The name of the atom it bonds to in the template.
        anchor: String,
    },
    /// The resolved atoms of one residue lie in more than one molecule, for
    /// example because a bridging atom is missing from the structure.
    ResidueFragmented {
        /// The residue.
        residue: SourceResidueId,
        /// The number of molecules its atoms lie in.
        fragments: usize,
    },
    /// The source states unit-cell parameters (mmCIF `_cell`, PDB `CRYST1`)
    /// that are not the structure's periodic cell, so the structure has none.
    CellOmitted {
        /// Why the cell was omitted.
        reason: CellOmissionReason,
        /// The source location of the cell parameters.
        location: SourceLocation,
    },
    /// The source states a biological assembly that cannot be read, so the
    /// report does not list it.
    AssemblyOmitted {
        /// The assembly's identifier, or `None` when the source gives none.
        assembly: Option<String>,
        /// Why the assembly was omitted.
        reason: AssemblyOmissionReason,
        /// The 1-based line of the assembly's first record, when known.
        line: Option<usize>,
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
            Self::ForceFieldTemplate {
                residue,
                component_id,
                force_field,
            } => write!(
                formatter,
                "{} (component {component_id:?}) fits no dictionary template. It is read with \
                 the {force_field} template of that name.",
                Capitalized(residue)
            ),
            Self::ResidueBondsInferred {
                residue,
                component_id,
                bonds,
                neutral_centres,
            } => {
                write!(
                    formatter,
                    "{} (component {component_id:?}) has no template or bond table. Its bonds \
                     come from atom distances: {}.",
                    Capitalized(residue),
                    counted(*bonds, "single bond", "single bonds")
                )?;
                write_neutral_centres(formatter, neutral_centres)
            }
            Self::ResidueBondOrdersUnstated {
                residue,
                component_id,
                bonds,
                neutral_centres,
            } => {
                write!(
                    formatter,
                    "{} (component {component_id:?}) is bonded by records that state no orders, \
                     and no template gives them. Its {} between heavy atoms {} read as single.",
                    Capitalized(residue),
                    counted(*bonds, "bond", "bonds"),
                    if *bonds == 1 { "is" } else { "are" }
                )?;
                write_neutral_centres(formatter, neutral_centres)
            }
            Self::ResidueFitsNoTemplate {
                residue,
                component_id,
            } => write!(
                formatter,
                "{} (component {component_id:?}) fits no template of its name. Its bonds are \
                 read from its connection records.",
                Capitalized(residue)
            ),
            Self::StatedOrdersDisagreeWithHydrogens {
                residue,
                component_id,
            } => write!(
                formatter,
                "{} (component {component_id:?}) has stated bond orders that disagree with its \
                 hydrogens. Its charges are derived from the stated orders and the hydrogens.",
                Capitalized(residue)
            ),
            Self::BondNotInTemplate {
                residue,
                component_id,
                atoms: [left, right],
            } => write!(
                formatter,
                "A connection record bonds atoms {left:?} and {right:?} of {residue} (component \
                 {component_id:?}). Its template does not bond them."
            ),
            Self::BondOrderDiffersFromTemplate {
                residue,
                component_id,
                atoms: [left, right],
                order,
                template_order,
            } => {
                write!(
                    formatter,
                    "A connection record states a {order} bond between atoms {left:?} and \
                     {right:?} of {residue} (component {component_id:?}). Its template has "
                )?;
                match template_order {
                    Some(order) => write!(formatter, "a {order} bond."),
                    None => formatter.write_str("an aromatic bond."),
                }
            }
            Self::HydrogenAttachedByDistance {
                residue,
                atom_name,
                heavy_atom_name,
            } => write!(
                formatter,
                "Hydrogen {atom_name:?} of {residue} has no stated bond. It is bonded to \
                 {heavy_atom_name:?}, the nearest heavy atom."
            ),
            Self::DisulfideInferred {
                residues: [left, right],
            } => write!(
                formatter,
                "{} and {right} have no stated disulfide. Their SG atoms are within bonding \
                 distance, so they are bonded.",
                Capitalized(left)
            ),
            Self::ChainBreak { residue, next } => write!(
                formatter,
                "{} and {next} are not bonded. Their linking atoms are too far apart.",
                Capitalized(residue)
            ),
            Self::PolymerLinkSkipped { residue, next } => write!(
                formatter,
                "{} and {next} have no standard polymer link. A connection record bonds one of \
                 their linking atoms to another atom.",
                Capitalized(residue)
            ),
            Self::LinkingAtomMissing { residue, next } => write!(
                formatter,
                "{} and {next} are not bonded. One of their linking atoms is missing.",
                Capitalized(residue)
            ),
            Self::PolymerLinkMissing {
                residue,
                component_id,
                atom_name,
            } => write!(
                formatter,
                "Atom {atom_name:?} of {residue} (component {component_id:?}) has no bond to \
                 its neighbour in the chain, which the model lacks or places too far away. It \
                 keeps the hydrogens it has when bonded."
            ),
            Self::StatedChargeResonance {
                residue,
                component_id,
                charged_atoms,
            } => write!(
                formatter,
                "The stated {} of {residue} (component {component_id:?}) give another \
                 resonance form than the template. The double bonds move to it.",
                if charged_atoms.len() == 1 {
                    format!("charge of {}", listed(&Charged::all(charged_atoms)))
                } else {
                    format!("charges of {}", listed(&Charged::all(charged_atoms)))
                }
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
            Self::MetalValenceIncomplete {
                residue,
                component_id,
                atom_name,
                valence,
            } => write!(
                formatter,
                "Atom {atom_name:?} of {residue} (component {component_id:?}) is a metal atom \
                 with valence {valence}, which its element does not allow. It takes no hydrogens."
            ),
            Self::MetalDonors {
                residue,
                component_id,
                atom_names,
            } => {
                let (atoms, verb) = if atom_names.len() == 1 {
                    ("Atom", "bonds metals as a donor and takes")
                } else {
                    ("Atoms", "bond metals as donors and take")
                };
                write!(
                    formatter,
                    "{} (component {component_id:?}) has no template charges. {atoms} {} {verb} \
                     no hydrogens.",
                    Capitalized(residue),
                    listed(&Quoted::all(atom_names))
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
            Self::ProtonationAssumed { components } => {
                let residues = components.iter().map(|(_, count)| count).sum::<usize>();
                let (have, take) = if residues == 1 {
                    ("has", "It takes its")
                } else {
                    ("have", "They take their")
                };
                let counts = components
                    .iter()
                    .map(|(component, count)| format!("{component:?} ({count})"))
                    .collect::<Vec<_>>();
                write!(
                    formatter,
                    "{} {have} no hydrogens on heteroatoms. {take} assumed protonation: \
                     {}. Acid groups take their state at pH 7, and polymer residues their \
                     usual one. Ligands take nothing else.",
                    counted(residues, "residue", "residues"),
                    listed(&counts)
                )
            }
            Self::HydrogensMissing {
                residue,
                component_id,
                atoms,
            } => {
                let (noun, pronoun) = if atoms.len() == 1 {
                    ("atom", "It takes")
                } else {
                    ("atoms", "They take")
                };
                write!(
                    formatter,
                    "{} (component {component_id:?}) lacks hydrogens on {noun} {}. {pronoun} \
                     hydrogens by valence.",
                    Capitalized(residue),
                    listed(&Quoted::all(atoms))
                )
            }
            Self::ZeroOccupancyHydrogens { residue, atoms } => {
                let (noun, verb) = if atoms.len() == 1 {
                    ("Hydrogen", "has zero occupancy. It is")
                } else {
                    ("Hydrogens", "have zero occupancy. They are")
                };
                write!(
                    formatter,
                    "{noun} {} of {residue} {verb} not read.",
                    listed(&Quoted::all(atoms))
                )
            }
            Self::HydrogenReplacedByCoordination {
                residue,
                component_id,
                atoms,
            } => {
                let (noun, verb) = if atoms.len() == 1 {
                    ("Atom", "coordinates a metal in place of a hydrogen")
                } else {
                    ("Atoms", "coordinate metals in place of hydrogens")
                };
                write!(
                    formatter,
                    "{noun} {} of {residue} (component {component_id:?}) {verb}.",
                    listed(&Quoted::all(atoms))
                )
            }
            Self::StatedChargeKept {
                residue,
                component_id,
                atom_name,
                stated,
                derived,
            } => write!(
                formatter,
                "Atom {atom_name:?} of {residue} (component {component_id:?}) states charge \
                 {}. Its template and hydrogens give {}. The stated charge is kept.",
                Signed(*stated),
                Signed(*derived)
            ),
            Self::StatedChargeRejected {
                residue,
                component_id,
                atom_name,
                stated,
                kept,
            } => write!(
                formatter,
                "Atom {atom_name:?} of {residue} (component {component_id:?}) states charge \
                 {}, which no hydrogen count it can have holds. It takes charge {} from its \
                 template and hydrogens.",
                Signed(*stated),
                Signed(*kept)
            ),
            Self::TemplateChargeUnknown {
                residue,
                component_id,
                atoms,
            } => {
                let (noun, pronoun) = if atoms.len() == 1 {
                    ("atom", "It is")
                } else {
                    ("atoms", "They are")
                };
                write!(
                    formatter,
                    "The template of {residue} (component {component_id:?}) states no charge \
                     for {noun} {}. {pronoun} read as neutral.",
                    listed(&Quoted::all(atoms))
                )
            }
            Self::TemplateMissing {
                residue,
                component_id,
            } => write!(
                formatter,
                "{} (component {component_id:?}) has no template, and the source states no \
                 charges for it. Its atoms are read as neutral.",
                Capitalized(residue)
            ),
            Self::AtomsMissing {
                residue,
                component_id,
                atoms,
            } => write!(
                formatter,
                "{} (component {component_id:?}) lacks template {} {}. Hydrogens take their \
                 place.",
                Capitalized(residue),
                if atoms.len() == 1 { "atom" } else { "atoms" },
                listed(&Quoted::all(atoms))
            ),
            Self::TerminalAtomMissing {
                residue,
                component_id,
                atom_name,
                anchor,
            } => write!(
                formatter,
                "{} (component {component_id:?}) lacks terminal atom {atom_name:?}. Atom \
                 {anchor:?} takes no hydrogen in its place.",
                Capitalized(residue)
            ),
            Self::ResidueFragmented { residue, fragments } => write!(
                formatter,
                "The atoms of {residue} lie in {}.",
                counted(*fragments, "molecule", "molecules")
            ),
            Self::CellOmitted { reason, location } => {
                write!(formatter, "The unit cell at {location} is not read. {reason}")
            }
            Self::AssemblyOmitted {
                assembly,
                reason,
                line,
            } => {
                match assembly {
                    Some(id) => write!(formatter, "Assembly {id:?}")?,
                    None => formatter.write_str("An assembly")?,
                }
                if let Some(line) = line {
                    write!(formatter, " at line {line}")?;
                }
                write!(formatter, " is not read. {reason}")
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

/// Appends the onium centres that took no charge, and their hydrogens, to an
/// issue about unknown bond orders.
fn write_neutral_centres(
    formatter: &mut fmt::Formatter<'_>,
    centres: &[(String, u8)],
) -> fmt::Result {
    if centres.is_empty() {
        return Ok(());
    }
    let atoms = centres
        .iter()
        .map(|(name, count)| {
            format!(
                "{name:?} took {}",
                counted(usize::from(*count), "hydrogen", "hydrogens")
            )
        })
        .collect::<Vec<_>>();
    let noun = if centres.len() == 1 { "Atom" } else { "Atoms" };
    write!(
        formatter,
        " Its onium centres took no charge. {noun} {}.",
        listed(&atoms)
    )
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

/// A formal charge in e with its sign, such as `+1`; zero is `0`.
struct Signed(i8);

impl fmt::Display for Signed {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.0 {
            0 => formatter.write_str("0"),
            charge => write!(formatter, "{charge:+}"),
        }
    }
}

/// Why stated unit-cell parameters are not a structure's periodic cell.
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

/// Why a stated biological assembly is not listed.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AssemblyOmissionReason {
    /// The assembly has no identifier.
    Unidentified,
    /// It states no chains, or chains without operators.
    NoGenerators,
    /// An operator expression cannot be read, or names more than 100,000
    /// operations.
    Expression {
        /// The expression as stated.
        expression: String,
    },
    /// An expression names an operator the source does not state.
    UnknownOperator {
        /// The operator's identifier.
        operator: String,
    },
    /// An operator's matrix or vector is incomplete or not numbers.
    UnreadableOperator {
        /// The operator's identifier.
        operator: String,
    },
    /// An operator's matrix is not a rotation, as a mirror or a scaling is
    /// not.
    NotRotation {
        /// The operator's identifier.
        operator: String,
    },
    /// With it, the source's assemblies would use their operators more
    /// than 100,000 times, counting each operator of a product.
    TooManyOperations,
}

/// Writes the reason as a sentence, such as `Operator "61" is not stated.`
impl fmt::Display for AssemblyOmissionReason {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unidentified => formatter.write_str("It has no identifier."),
            Self::NoGenerators => formatter.write_str("Its chains or operators are not stated."),
            Self::Expression { expression } => {
                write!(
                    formatter,
                    "Operator expression {expression:?} cannot be read."
                )
            }
            Self::UnknownOperator { operator } => {
                write!(formatter, "Operator {operator:?} is not stated.")
            }
            Self::UnreadableOperator { operator } => write!(
                formatter,
                "Operator {operator:?} has missing or unreadable numbers."
            ),
            Self::NotRotation { operator } => {
                write!(formatter, "Operator {operator:?} is not a rotation.")
            }
            Self::TooManyOperations => formatter.write_str(
                "With it, the file's assemblies would apply more than 100,000 operators.",
            ),
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
    /// A disulfide between the `SG` atoms of two cysteines within bonding
    /// distance, which no record states.
    Disulfide,
    /// Covalent radii, for a residue without template or bond table.
    Distance,
    /// A hydrogen that no template, component bond table, or listed bond
    /// names, bonded to the nearest heavy atom of its residue.
    HydrogenDistance,
    /// A terminal hydrogen bonded to the atom its conventional name places
    /// it on: `H1`-`H3`, `HT1`-`HT3` or `HN1`-`HN3` on a peptide's `N`,
    /// `HO5'` or `H5T` on a nucleotide's `O5'`, and `H3T` or `HO3*` on its
    /// `O3'`.
    HydrogenName,
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
/// The mmCIF and PDB writers state these records again, as `metalc` rows
/// and `LINK` records, when they are given the report. Coordination within
/// one residue, such as heme's Fe-N bonds, stays a dative bond.
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

/// The coordination records of a report that match a topology; see
/// [`StructureInterpretationReport::coordination_in`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CoordinationMatch {
    /// The partners of each matching record, in record order.
    pub(crate) records: Vec<[InstanceAtomId; 2]>,
    /// The records that do not match.
    pub(crate) unmatched: usize,
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
    pub(crate) assemblies: Vec<Assembly>,
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

    /// The coordination records that still match `topology`, as writers
    /// state them again. A record matches when each partner is a topology
    /// atom whose atom site keeps the atom names, type symbol, and component
    /// the report recorded for it; a structure edited since it was read can
    /// hold other atoms under the recorded IDs, or none.
    pub(crate) fn coordination_in(&self, topology: &Topology) -> CoordinationMatch {
        let partners = self
            .coordination
            .iter()
            .flat_map(|record| record.atoms)
            .collect::<HashSet<_>>();
        let sources = self
            .instances
            .iter()
            .flat_map(|instance| &instance.atoms)
            .filter(|atom| partners.contains(&atom.atom))
            .map(|atom| (atom.atom, atom))
            .collect::<HashMap<_, _>>();
        let matches = |atom: InstanceAtomId| {
            let (Some(source), Some(site)) = (sources.get(&atom), topology.atom_site_of(atom))
            else {
                return false;
            };
            let metadata = site.metadata();
            metadata.label_atom_id() == source.label_atom_name.as_deref()
                && metadata.auth_atom_id() == source.auth_atom_name.as_deref()
                && metadata.type_symbol() == Some(source.type_symbol.as_str())
                && site.residue().name() == source.component_id
        };
        let (records, unmatched): (Vec<_>, Vec<_>) = self
            .coordination
            .iter()
            .map(|record| record.atoms)
            .partition(|atoms| atoms.iter().all(|&atom| matches(atom)));
        CoordinationMatch {
            records,
            unmatched: unmatched.len(),
        }
    }

    /// Decisions for the selected source coordinate model only. The source
    /// document retains every alternative, including omitted rows and metadata.
    pub fn alternate_locations(&self) -> &[AltLocDecision] {
        &self.alternate_locations
    }

    /// The polymer entities the block states, in the block's order, with
    /// their canonical sequences and the structure's residues at each position.
    /// A PDB entry has one for each chain whose `SEQRES` sequence places its
    /// residues. Entities whose
    /// monomers neither `_entity_poly_seq` nor the atom sites give are left
    /// out. [`crate::sequence`] gives the sequences of the chains of any
    /// topology.
    pub fn polymer_entities(&self) -> &[PolymerEntity] {
        &self.polymer_entities
    }

    /// The biological assemblies the source states, in source order: mmCIF
    /// `_pdbx_struct_assembly`, PDB `REMARK 350` biomolecules. Assemblies
    /// that cannot be read are left out with
    /// [`StructureInterpretIssue::AssemblyOmitted`].
    /// [`Assembly::to_structure`] builds one.
    pub fn assemblies(&self) -> &[Assembly] {
        &self.assemblies
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
/// One final canonical structure and its interpretation report.
pub struct StructureInterpretation {
    pub(crate) structure: Structure,
    pub(crate) report: StructureInterpretationReport,
}

impl StructureInterpretation {
    /// The interpreted structure.
    pub fn structure(&self) -> &Structure {
        &self.structure
    }

    /// How the structure was read from its source.
    pub fn report(&self) -> &StructureInterpretationReport {
        &self.report
    }

    /// The topology of [`Self::structure`].
    pub fn topology(&self) -> &crate::topology::Topology {
        self.structure.topology()
    }

    /// Borrows molecule definitions once per instance in authoritative
    /// topology instance order.
    pub fn molecules(&self) -> impl ExactSizeIterator<Item = &Molecule> + DoubleEndedIterator {
        self.structure
            .topology()
            .instances()
            .map(|instance| instance.molecule())
    }

    /// Clones each molecule definition once per instance, in topology instance
    /// order.
    pub fn into_molecules(self) -> Vec<Molecule> {
        self.structure
            .topology()
            .instances()
            .map(|instance| instance.molecule().clone())
            .collect()
    }

    /// Retains shared ownership of the structure's exact topology allocation.
    pub fn into_topology(self) -> Arc<Topology> {
        self.structure.shared_topology()
    }

    /// The interpreted structure, discarding the report.
    pub fn into_structure(self) -> Structure {
        self.structure
    }

    /// The structure and its report.
    pub fn into_parts(self) -> (Structure, StructureInterpretationReport) {
        (self.structure, self.report)
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

/// A structure-set-level failure with typed `detail`.
pub(in crate::io) fn structure_set_error(
    kind: DiagnosticKind,
    message: impl Into<String>,
    detail: FormatErrorDetail,
) -> FormatError {
    FormatError::interpret(Format::Mmcif, kind, message).with_detail(detail)
}

/// Converts deposited coordinate models to geometries, applying `altloc_policy`
/// within each model. Alternate locations never generate additional geometries or
/// statistical weights. A single coordinate model yields a single geometry.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StructureSetInterpretOptions {
    entity_metadata: EntityMetadataPolicy,
    altloc_policy: AltLocPolicy,
    connectivity: ConnectivityOptions,
    model_ids: Option<Vec<String>>,
    cancellation: Option<Cancellation>,
}

impl StructureSetInterpretOptions {
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
            return Err(super::structure_set::empty_model_selection());
        }
        self.model_ids = model_ids;
        Ok(self)
    }

    /// The cancellation interpretation polls, if any; `None` by default.
    ///
    /// Interpretation polls it between its stages and periodically inside
    /// them, for every geometry, and once it is requested fails with a
    /// [`FormatError`] of kind [`DiagnosticKind::Cancelled`] instead of
    /// returning a partial structure set. It does not change what a completed
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

impl Default for StructureSetInterpretOptions {
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

/// One shared-topology structure set and one report per geometry, in geometry
/// order.
#[derive(Debug, Clone)]
pub struct StructureSetInterpretation {
    pub(crate) set: StructureSet,
    pub(crate) reports: Vec<StructureInterpretationReport>,
}

impl StructureSetInterpretation {
    /// The interpreted structure set, one geometry per selected model.
    pub fn structure_set(&self) -> &StructureSet {
        &self.set
    }

    /// One report per geometry, in geometry order.
    pub fn reports(&self) -> &[StructureInterpretationReport] {
        &self.reports
    }

    /// The topology all geometries share.
    pub fn topology(&self) -> &Topology {
        self.set.topology()
    }

    /// The interpreted structure set, discarding the reports.
    pub fn into_structure_set(self) -> StructureSet {
        self.set
    }

    /// The structure set and its per-geometry reports.
    pub fn into_parts(self) -> (StructureSet, Vec<StructureInterpretationReport>) {
        (self.set, self.reports)
    }
}

#[cfg(test)]
mod tests;
