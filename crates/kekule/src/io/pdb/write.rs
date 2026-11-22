//! PDB writing of structures and structure sets.
//!
//! A written entry holds these records, in this order:
//! - `SEQRES` for chains whose first or last polymer residue bonds to a
//!   residue the structure lacks, with an `UNK` residue past that end;
//! - `LINK` for the metal coordination of interpretation reports, when the
//!   source has them ([`PdbSource::with_reports`]);
//! - `CRYST1` for a representable cell;
//! - `MODEL` and `ENDMDL` around each geometry of a structure set, none
//!   for a single structure;
//! - `ATOM` and `HETATM` in dense atom order, with a `TER` after each
//!   chain's last `ATOM` record;
//! - `CONECT` for the bonds a reader cannot derive from the Chemical
//!   Component Dictionary;
//! - `END`.
//!
//! A chain identifier longer than the one-character chain column, such as
//! the `SOL1` and `ION1` chains [`crate::solvation::solvate`] adds, is
//! written as a segment identifier under [`PdbChainIdPolicy::SegmentId`].
use crate::diagnostics::{Countable, LengthMismatch};
use crate::fixed_u32_id::FixedId;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt;
use std::io::Write;
use std::sync::Arc;

use crate::diagnostics::{
    AtomField, Capitalized, CellReason, Count, DiagnosticKind, Format, FormatError,
    FormatErrorDetail, GeometryField,
};
use crate::io::element_symbol::structure_element_symbol;
use crate::io::macromolecular::connectivity::{
    hydrogen_sites, is_coordination, link_sides, link_valence, reader_protonation, reader_template,
    standard_link, within_bonding_distance, AtomEvidence, PolymerLink, Protonation,
    ResidueEvidence, ResidueTemplate, StandardLink,
};
use crate::io::macromolecular::{
    is_metal, StructureInterpretation, StructureInterpretationReport, StructureSetInterpretation,
};
use crate::io::mmcif::write::{bond_order_name, cell_values};
use crate::io::write_policy::{
    check_geometry_state, record_multi_model_structure_set, record_structure_set_properties,
    unrepresentable_state, Carries, Omissions, StatePolicies,
};
use crate::io::{Component, ComponentLibrary, UnrepresentablePolicy, WriteOmission, WriteReport};
use crate::math::{PeriodicCell, PeriodicCellError, Point3};
use crate::molecule::{BondOrder, ImplicitHydrogens};
use crate::structure::{AsStructureView, GeometryView, Structure, StructureSet, StructureView};
use crate::topology::{
    InstanceAtomId, InstanceBondId, MoleculeClass, MoleculeInstanceId, ResidueClass, ResidueId,
    Topology, TopologyAtomIndex,
};
use crate::units::{ANGSTROM, SQUARE_ANGSTROM};

/// PDB carries occupancies and B factors, and a cell when [`entry_cell`]
/// can write it.
const PDB_CARRIES: Carries = Carries {
    cell: true,
    sites: true,
    ..Carries::POSITIONS
};

/// The largest atom serial number five hybrid-36 columns hold.
const MAX_ATOM_SERIAL: i64 = 87_440_031;
/// The smallest residue sequence number four columns hold.
const MIN_RESIDUE_NUMBER: i64 = -999;
/// The largest residue sequence number four hybrid-36 columns hold.
const MAX_RESIDUE_NUMBER: i64 = 2_436_111;

/// How the PDB writer writes a chain identifier longer than PDB's
/// one-character chain column.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Hash)]
pub enum PdbChainIdPolicy {
    /// Fail with [`PdbWriteErrorKind::ChainIdTooLong`].
    #[default]
    Fail,
    /// Write an identifier of up to four printable ASCII characters, such as
    /// a chain [`crate::solvation::solvate`] adds, in the segment identifier
    /// columns 73-76 and leave the chain column blank, as CHARMM and NAMD
    /// name segments. Kekule's PDB reader, like CHARMM and VMD, reads the
    /// segment of an atom with a blank chain column as its chain, so the
    /// chain keeps its identifier; readers that take chains from the chain
    /// column alone see a blank chain. [`WriteOmission::ChainIds`] counts the
    /// chains written so; a longer identifier still fails.
    SegmentId,
}

/// Options of the PDB writer.
///
/// The default fails on every state PDB cannot carry. Identifiers that do
/// not fit their columns fail, as omitting or shortening them would change
/// which atoms a reader identifies, except chain identifiers that
/// [`PdbChainIdPolicy::SegmentId`] writes as segment identifiers.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct PdbWriteOptions {
    geometry_metadata: UnrepresentablePolicy,
    properties: UnrepresentablePolicy,
    cell: UnrepresentablePolicy,
    implicit_hydrogens: UnrepresentablePolicy,
    chemistry: UnrepresentablePolicy,
    chain_ids: PdbChainIdPolicy,
    components: Option<Arc<ComponentLibrary>>,
}

impl PdbWriteOptions {
    /// The state of geometries that PDB cannot carry, whether of a structure,
    /// a structure set, or a trajectory frame. A PDB reader gives geometries
    /// without an order, so by default a set fails to write when a geometry
    /// has a place in its order (a frame's time and MD step among them),
    /// velocities, or forces, or when its time order has neither times nor
    /// steps; so does a single geometry with a place, or a geometry or
    /// structure with velocities or forces.
    /// [`UnrepresentablePolicy::Omit`] writes the models without them and
    /// reports a [`WriteOmission::GeometryState`] for each field.
    pub const fn geometry_metadata(&self) -> UnrepresentablePolicy {
        self.geometry_metadata
    }

    /// Sets [`Self::geometry_metadata`].
    #[must_use]
    pub const fn with_geometry_metadata(
        mut self,
        geometry_metadata: UnrepresentablePolicy,
    ) -> Self {
        self.geometry_metadata = geometry_metadata;
        self
    }

    /// A structure's per-atom, per-bond, and own properties and the
    /// properties and columns of a structure set or of a geometry's set,
    /// which PDB cannot carry: by default a structure or set with any fails to write, and
    /// under [`UnrepresentablePolicy::Omit`] the report names each
    /// ([`WriteOmission::GeometryState`]).
    pub const fn properties(&self) -> UnrepresentablePolicy {
        self.properties
    }

    /// Sets [`Self::properties`].
    #[must_use]
    pub const fn with_properties(mut self, properties: UnrepresentablePolicy) -> Self {
        self.properties = properties;
        self
    }

    /// The policies for the state of one structure.
    const fn state_policies(&self) -> StatePolicies {
        StatePolicies {
            carries: PDB_CARRIES,
            geometry_metadata: self.geometry_metadata,
            properties: self.properties,
            cell: self.cell,
        }
    }

    /// A periodic cell `CRYST1` cannot carry as the structure's cell: one not
    /// periodic along all three axes, not in the standard crystallographic
    /// orientation, differing between geometries, or written as PDB's
    /// 1 Å placeholder cube. By default such a cell fails to write;
    /// representable cells are written as `CRYST1` in space group P 1.
    pub const fn cell(&self) -> UnrepresentablePolicy {
        self.cell
    }

    /// Sets [`Self::cell`].
    #[must_use]
    pub const fn with_cell(mut self, cell: UnrepresentablePolicy) -> Self {
        self.cell = cell;
        self
    }

    /// Fixed implicit-hydrogen counts, which PDB cannot carry: readers infer
    /// implicit hydrogens from valence. A count readers fix again is written
    /// as it is: that of an atom whose residue lacks a leaving atom bonded to
    /// it, such as the carbonyl carbon of a residue without `OXT`, or of a
    /// polymer linking atom whose neighbour the chain lacks, such as the `N`
    /// after a gap. `SEQRES` states the residues past a chain end that is
    /// such a gap. By default
    /// another atom with one fails to write; [`UnrepresentablePolicy::Omit`]
    /// writes it without the count and reports
    /// [`WriteOmission::ImplicitHydrogens`]. Structures whose files list
    /// every atom, such as XYZ interpretations, fix their counts at zero.
    pub const fn implicit_hydrogens(&self) -> UnrepresentablePolicy {
        self.implicit_hydrogens
    }

    /// Sets [`Self::implicit_hydrogens`].
    #[must_use]
    pub const fn with_implicit_hydrogens(
        mut self,
        implicit_hydrogens: UnrepresentablePolicy,
    ) -> Self {
        self.implicit_hydrogens = implicit_hydrogens;
        self
    }

    /// Other represented chemistry PDB cannot carry: isotopes other than
    /// deuterium and tritium, radicals, atom maps, represented stereo, bonds
    /// `CONECT` cannot state (see [`PdbWriteErrorKind::UnsupportedBondOrder`]),
    /// and residue bonds that differ from the dictionary readers apply (see
    /// [`PdbWriteErrorKind::DictionaryBondConflict`]). By default it fails to
    /// write; [`UnrepresentablePolicy::Omit`] writes the records without it
    /// and reports each kind left out. Fixed implicit-hydrogen counts have
    /// their own policy, [`Self::implicit_hydrogens`].
    pub const fn chemistry(&self) -> UnrepresentablePolicy {
        self.chemistry
    }

    /// Sets [`Self::chemistry`].
    #[must_use]
    pub const fn with_chemistry(mut self, chemistry: UnrepresentablePolicy) -> Self {
        self.chemistry = chemistry;
        self
    }

    /// Chain identifiers longer than PDB's one-character chain column: by
    /// default they fail; [`PdbChainIdPolicy::SegmentId`] writes those of up
    /// to four characters as segment identifiers.
    pub const fn chain_ids(&self) -> PdbChainIdPolicy {
        self.chain_ids
    }

    /// Sets [`Self::chain_ids`].
    #[must_use]
    pub const fn with_chain_ids(mut self, chain_ids: PdbChainIdPolicy) -> Self {
        self.chain_ids = chain_ids;
        self
    }

    /// The templates readers of the entry consult before the bundled ones,
    /// as [`ConnectivityOptions::components`](crate::macromolecular::ConnectivityOptions::components)
    /// gives them. The writer leaves out the bonds and fixed hydrogen
    /// counts that reading with them gives back, and checks residue bonds
    /// against them. `None` by default: readers use the bundled templates.
    pub fn components(&self) -> Option<&Arc<ComponentLibrary>> {
        self.components.as_ref()
    }

    /// Sets [`Self::components`].
    #[must_use]
    pub fn with_components(mut self, components: Option<Arc<ComponentLibrary>>) -> Self {
        self.components = components;
        self
    }
}

/// A failure specific to the PDB writer: the typed
/// [`detail`](crate::diagnostics::FormatError::detail) of a write-stage
/// [`FormatError`], as
/// [`FormatErrorDetail::PdbWrite`](crate::diagnostics::FormatErrorDetail::PdbWrite).
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PdbWriteErrorKind {
    /// The chain identifier is longer than PDB's one column, and under
    /// [`PdbChainIdPolicy::SegmentId`] than the four segment identifier
    /// columns or not printable ASCII. [`PdbChainIdPolicy::SegmentId`]
    /// writes identifiers of up to four printable characters.
    ChainIdTooLong {
        /// An atom of the chain.
        atom: InstanceAtomId,
        /// The identifier.
        chain_id: String,
    },
    /// The residue name is longer than PDB's three columns.
    ResidueNameTooLong {
        /// An atom of the residue.
        atom: InstanceAtomId,
        /// The name.
        name: String,
    },
    /// The atom name is empty, longer than PDB's four columns, or not
    /// printable ASCII without whitespace.
    InvalidAtomName {
        /// The atom.
        atom: InstanceAtomId,
        /// The name.
        name: String,
    },
    /// Two atoms of one residue have the same name, so a reader could not
    /// tell them apart.
    DuplicateAtomName {
        /// The second atom with the name.
        atom: InstanceAtomId,
        /// The name.
        name: String,
    },
    /// The insertion code is longer than PDB's one column.
    InsertionCodeTooLong {
        /// An atom of the residue.
        atom: InstanceAtomId,
        /// The code.
        code: String,
    },
    /// The residue has no sequence number, or one that is not an integer
    /// from -999 to 2,436,111, the range of four hybrid-36 columns.
    InvalidResidueNumber {
        /// An atom of the residue.
        atom: InstanceAtomId,
        /// The author sequence number, if any.
        number: Option<String>,
    },
    /// The residue's atoms are not consecutive in dense atom order, which
    /// writers keep, so a reader would split the residue.
    NonContiguousResidue {
        /// The first atom that resumes the residue.
        atom: InstanceAtomId,
    },
    /// More atom and `TER` records than five hybrid-36 columns number.
    TooManyAtoms {
        /// The serial numbers the entry would need.
        serials: usize,
    },
    /// A coordinate does not fit the `%8.3f` columns, -999.999 to 9999.999 Å.
    CoordinateOutOfRange {
        /// The atom.
        atom: InstanceAtomId,
    },
    /// An occupancy does not fit the `%6.2f` columns, -99.99 to 999.99.
    OccupancyOutOfRange {
        /// The atom.
        atom: InstanceAtomId,
    },
    /// A B factor does not fit the `%6.2f` columns, -99.99 to 999.99 Å².
    BFactorOutOfRange {
        /// The atom.
        atom: InstanceAtomId,
    },
    /// A formal charge outside PDB's one-digit range -9..=9.
    FormalChargeOutOfRange {
        /// The atom.
        atom: InstanceAtomId,
        /// Its charge.
        charge: i8,
    },
    /// An atom has state PDB cannot carry, with
    /// [`PdbWriteOptions::implicit_hydrogens`] (for
    /// [`AtomField::Hydrogens`]) or [`PdbWriteOptions::chemistry`] set to
    /// fail.
    UnsupportedAtomField {
        /// The atom.
        atom: InstanceAtomId,
        /// The field.
        field: AtomField,
    },
    /// A molecule instance carries represented stereo, with
    /// [`PdbWriteOptions::chemistry`] set to fail.
    UnsupportedStereo(MoleculeInstanceId),
    /// A bond `CONECT` cannot state, with [`PdbWriteOptions::chemistry`] set
    /// to fail. Readers take a pair listed once as single, twice as double,
    /// and three times as triple, and a pair of a metal and a non-metal as
    /// dative whatever its listing; zero-order and quadruple bonds, dative
    /// bonds between other pairs, and other orders between a metal and a
    /// non-metal have no listing.
    UnsupportedBondOrder {
        /// The bond.
        bond: InstanceBondId,
        /// Its order.
        order: BondOrder,
    },
    /// A periodic cell `CRYST1` cannot carry, with
    /// [`PdbWriteOptions::cell`] set to fail.
    UnrepresentableCell {
        /// Why the cell cannot be written.
        reason: CellReason,
    },
    /// The bonds between two named heavy atoms of a residue the bundled
    /// [`ComponentLibrary`] covers differ from its template, which a reader
    /// applies to the residue whatever `CONECT` lists, with
    /// [`PdbWriteOptions::chemistry`] set to fail. A source with its own
    /// component bond table, such as an mmCIF `_chem_comp_bond` category that
    /// states a heme's iron-nitrogen bonds as single where the dictionary has
    /// coordination, can hold such bonds.
    DictionaryBondConflict {
        /// The two atoms.
        atoms: [InstanceAtomId; 2],
        /// The structure's bond order between them, if they are bonded.
        order: Option<BondOrder>,
    },
    /// The 1-based geometry `geometry` of a structure set has state PDB cannot
    /// carry, with [`PdbWriteOptions::geometry_metadata`] set to fail. A
    /// bare time order is reported for the first geometry.
    UnrepresentableGeometry {
        /// The 1-based geometry number.
        geometry: usize,
        /// The state; see [`PdbWriteOptions::geometry_metadata`].
        field: GeometryField,
    },
    /// The structure set has no geometries.
    EmptyStructureSet,
    /// The number of interpretation reports does not match the number of
    /// structures: one for a structure, one per geometry or item.
    ReportCountMismatch(LengthMismatch),
}

impl PdbWriteErrorKind {
    /// The diagnostic kind of a write error with this detail:
    /// [`DiagnosticKind::Inconsistent`] for an empty structure set and for
    /// reports that do not match the structures, and
    /// [`DiagnosticKind::Unsupported`] otherwise.
    pub const fn kind(&self) -> DiagnosticKind {
        match self {
            Self::EmptyStructureSet | Self::ReportCountMismatch(_) => DiagnosticKind::Inconsistent,
            _ => DiagnosticKind::Unsupported,
        }
    }
}

impl fmt::Display for PdbWriteErrorKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ChainIdTooLong { atom, chain_id } => write!(
                f,
                "Chain ID {chain_id:?} of {atom} does not fit PDB's chain columns."
            ),
            Self::ResidueNameTooLong { atom, name } => write!(
                f,
                "Residue name {name:?} of {atom} is longer than PDB's three characters."
            ),
            Self::InvalidAtomName { atom, name } => write!(
                f,
                "Atom name {name:?} of {atom} does not fit PDB's four atom-name columns."
            ),
            Self::DuplicateAtomName { atom, name } => write!(
                f,
                "{} repeats atom name {name:?} in its residue.",
                Capitalized(atom)
            ),
            Self::InsertionCodeTooLong { atom, code } => write!(
                f,
                "Insertion code {code:?} of {atom} is longer than PDB's one character."
            ),
            Self::InvalidResidueNumber { atom, number } => match number {
                Some(number) => write!(
                    f,
                    "Residue number {number:?} of {atom} is not a PDB residue number."
                ),
                None => write!(f, "The residue of {atom} has no sequence number."),
            },
            Self::NonContiguousResidue { atom } => write!(
                f,
                "{} resumes its residue after atoms of another residue.",
                Capitalized(atom)
            ),
            Self::TooManyAtoms { serials } => write!(
                f,
                "The entry needs {} serial numbers. PDB has {}.",
                Count(*serials),
                Count(MAX_ATOM_SERIAL as usize)
            ),
            Self::CoordinateOutOfRange { atom } => write!(
                f,
                "A coordinate of {atom} is outside PDB's range of -999.999 to 9999.999 Å."
            ),
            Self::OccupancyOutOfRange { atom } => write!(
                f,
                "The occupancy of {atom} is outside PDB's range of -99.99 to 999.99."
            ),
            Self::BFactorOutOfRange { atom } => write!(
                f,
                "The B factor of {atom} is outside PDB's range of -99.99 to 999.99 Å²."
            ),
            Self::FormalChargeOutOfRange { atom, charge } => write!(
                f,
                "{} has formal charge {charge}, outside PDB's range of -9 to 9.",
                Capitalized(atom)
            ),
            Self::UnsupportedAtomField { atom, field } => write!(
                f,
                "The {field} of {atom} {} no PDB representation.",
                field.verb()
            ),
            Self::UnsupportedStereo(molecule) => write!(
                f,
                "The stereochemistry of {molecule} has no PDB representation."
            ),
            Self::UnsupportedBondOrder { bond, order } => write!(
                f,
                "The {} order of {bond} has no CONECT representation.",
                bond_order_name(*order)
            ),
            Self::UnrepresentableCell { reason } => reason.fmt(f),
            Self::DictionaryBondConflict { atoms, order } => match order {
                Some(order) => write!(
                    f,
                    "The {} bond between {} and {} differs from the component dictionary.",
                    bond_order_name(*order),
                    atoms[0],
                    atoms[1]
                ),
                None => write!(
                    f,
                    "{} and {} are not bonded, but the component dictionary bonds them.",
                    Capitalized(atoms[0]),
                    atoms[1]
                ),
            },
            Self::UnrepresentableGeometry { geometry, field } => {
                f.write_str(&unrepresentable_state("PDB", *field, Some(*geometry)))
            }
            Self::EmptyStructureSet => f.write_str("The structure set has no geometries."),
            Self::ReportCountMismatch(error) => error.fmt(f),
        }
    }
}

fn fail(issue: PdbWriteErrorKind) -> FormatError {
    FormatError::write(Format::Pdb, issue.kind(), issue.to_string())
        .with_detail(FormatErrorDetail::PdbWrite(issue))
}

/// One PDB entry to write: a structure or a structure set, and the
/// interpretation reports whose metal coordination it states.
///
/// Convert from `&Structure`, `&StructureSet`, [`GeometryView`],
/// `&StructureInterpretation`, or `&StructureSetInterpretation`;
/// interpretations keep their source reports. [`Self::structure`] accepts
/// any borrowed structure view, such as a geometry of a structure set or a
/// trajectory frame, whose own state fails to write unless the write
/// options omit it.
#[derive(Debug, Clone, Copy)]
pub struct PdbSource<'a> {
    content: Content<'a>,
    reports: Option<&'a [StructureInterpretationReport]>,
}

#[derive(Debug, Clone, Copy)]
enum Content<'a> {
    Structure(StructureView<'a>),
    StructureSet(&'a StructureSet),
}

impl<'a> PdbSource<'a> {
    /// One structure, written without `MODEL` records.
    pub fn structure(structure: &'a (impl AsStructureView + ?Sized)) -> Self {
        Self {
            content: Content::Structure(structure.as_structure_view()),
            reports: None,
        }
    }

    /// A structure set, one `MODEL` per geometry in geometry order. The state PDB
    /// cannot carry follows [`PdbWriteOptions::geometry_metadata`], and the
    /// set's and geometries' properties [`PdbWriteOptions::properties`].
    pub fn structure_set(set: &'a StructureSet) -> Self {
        Self {
            content: Content::StructureSet(set),
            reports: None,
        }
    }

    /// States the metal coordination of interpretation reports as `LINK`
    /// records: one report for a structure, one per geometry for a structure
    /// set. Coordination between residues is no bond of the
    /// structure, so without reports the entry has no `LINK` records.
    ///
    /// Each record names both partners as their atom records do, with blank
    /// alternate locations, the identity symmetry operator `1555`, and the
    /// partners' distance in the first structure. `LINK` records belong to
    /// the entry rather than to a model, so a structure set states the
    /// coordination every report states. Each pair of atoms is stated once.
    /// A record is left out and counted as [`WriteOmission::Coordination`]
    /// when:
    /// - the structure no longer matches the report: a partner is missing,
    ///   or its atom site has other names or another component than the
    ///   report recorded;
    /// - another report of the structure set does not state it;
    /// - a partner's chain is written as a segment identifier, which the
    ///   `LINK` chain column cannot hold;
    /// - the partners are not one metal and one non-metal, such as two metal
    ///   ions an mmCIF `metalc` row joins: readers take such a `LINK` for a
    ///   bond;
    /// - a partner's identifiers also name another atom: readers take a blank
    ///   chain or insertion code in a `LINK` for any, so residue `10` also
    ///   names residue `10A`.
    #[must_use]
    pub fn with_reports(mut self, reports: &'a [StructureInterpretationReport]) -> Self {
        self.reports = Some(reports);
        self
    }
}

impl<'a> From<&'a Structure> for PdbSource<'a> {
    fn from(structure: &'a Structure) -> Self {
        Self::structure(structure)
    }
}

impl<'a> From<&'a StructureSet> for PdbSource<'a> {
    fn from(set: &'a StructureSet) -> Self {
        Self::structure_set(set)
    }
}

impl<'a, R: ?Sized> From<GeometryView<'a, R>> for PdbSource<'a> {
    fn from(geometry: GeometryView<'a, R>) -> Self {
        Self {
            content: Content::Structure(geometry.as_structure_view()),
            reports: None,
        }
    }
}

impl<'a> From<&'a StructureSetInterpretation> for PdbSource<'a> {
    fn from(interpretation: &'a StructureSetInterpretation) -> Self {
        Self::structure_set(interpretation.structure_set()).with_reports(interpretation.reports())
    }
}

impl<'a> From<&'a StructureInterpretation> for PdbSource<'a> {
    fn from(interpretation: &'a StructureInterpretation) -> Self {
        Self::structure(interpretation.structure())
            .with_reports(std::slice::from_ref(interpretation.report()))
    }
}

/// Writes one entry to a string. See [`write_pdb_to`].
pub(crate) fn write_pdb<'a>(
    source: impl Into<PdbSource<'a>>,
    options: PdbWriteOptions,
) -> Result<(String, WriteReport), FormatError> {
    let mut output = Vec::new();
    let report = write_pdb_to(&mut output, source, options)?;
    Ok((
        String::from_utf8(output).expect("the PDB writer emits UTF-8"),
        report,
    ))
}

/// Validates the whole entry, then streams its records into `writer`.
pub(crate) fn write_pdb_to<'a>(
    writer: &mut impl Write,
    source: impl Into<PdbSource<'a>>,
    options: PdbWriteOptions,
) -> Result<WriteReport, FormatError> {
    let mut omissions = Omissions::default();
    let source = source.into();
    let (views, models) = match source.content {
        Content::Structure(view) => {
            omissions.extend(check_geometry_state(
                Format::Pdb,
                view,
                options.state_policies(),
            )?);
            (vec![view], false)
        }
        Content::StructureSet(set) => {
            if set.is_empty() {
                return Err(fail(PdbWriteErrorKind::EmptyStructureSet));
            }
            let unrepresentable = |(geometry, field)| {
                fail(PdbWriteErrorKind::UnrepresentableGeometry { geometry, field })
            };
            record_multi_model_structure_set(set, options.geometry_metadata, &mut omissions)
                .map_err(unrepresentable)?;
            record_structure_set_properties(set, options.properties, &mut omissions)
                .map_err(unrepresentable)?;
            let views = set
                .iter()
                .map(|geometry| geometry.as_structure_view())
                .collect();
            (views, true)
        }
    };
    let cell = entry_cell(views.iter().map(|view| view.cell()), options.cell)?;
    if cell.is_none() {
        omissions.record(WriteOmission::Cell {
            structures: views.iter().filter(|view| view.cell().is_some()).count(),
        });
    }
    let positions = views
        .iter()
        .map(|view| written_positions(*view))
        .collect::<Vec<_>>();
    let mut plan = Plan::new(views[0].topology(), &positions, &options, &mut omissions)?;

    let links = match source.reports {
        Some(reports) if reports.len() != views.len() => {
            return Err(fail(PdbWriteErrorKind::ReportCountMismatch(
                LengthMismatch::new(Countable::Report, views.len(), reports.len()),
            )));
        }
        Some(reports) => plan.links(&views, reports, &mut omissions),
        None => Vec::new(),
    };
    plan.check_reader(
        views[0].topology(),
        views[0],
        &links,
        &options,
        &mut omissions,
    )?;
    // Check every structure's values before the first byte is written.
    let rows = views
        .iter()
        .map(|view| structure_values(*view, &plan))
        .collect::<Result<Vec<_>, _>>()?;
    let io = |error: std::io::Error| FormatError::sink(Format::Pdb, error);
    plan.write_seqres(writer).map_err(io)?;
    for link in &links {
        plan.write_link(writer, *link, &rows[0]).map_err(io)?;
    }
    if let Some([a, b, c, alpha, beta, gamma]) = cell {
        writeln!(
            writer,
            "CRYST1{a:>9.3}{b:>9.3}{c:>9.3}{alpha:>7.2}{beta:>7.2}{gamma:>7.2} {:<11}{:>4}",
            "P 1", 1
        )
        .map_err(io)?;
    }
    for (index, rows) in rows.iter().enumerate() {
        if models {
            writeln!(writer, "MODEL     {:>4}", index + 1).map_err(io)?;
        }
        plan.write_atoms(writer, rows).map_err(io)?;
        if models {
            writeln!(writer, "ENDMDL").map_err(io)?;
        }
    }
    plan.write_conect(writer).map_err(io)?;
    writeln!(writer, "END").map_err(io)?;
    Ok(omissions.into_report())
}

/// The `CRYST1` parameters of an entry whose structures have these cells, or
/// none when they have none or the caller omits one PDB cannot represent.
fn entry_cell<'a>(
    cells: impl IntoIterator<Item = Option<&'a PeriodicCell>>,
    policy: UnrepresentablePolicy,
) -> Result<Option<[f64; 6]>, FormatError> {
    let mut cells = cells.into_iter();
    let first = cells.next().flatten();
    let parameters = if cells.any(|cell| cell != first) {
        Err(CellReason::GeometryCellsDiffer)
    } else {
        match first.map(|cell| cell_values(*cell)).transpose() {
            Ok(Some(parameters)) if is_written_placeholder(parameters) => {
                Err(CellReason::Placeholder)
            }
            Ok(Some(parameters))
                if parameters[..3]
                    .iter()
                    .any(|&length| rounded(length, 3) > 99_999.999)
                    || parameters[3..]
                        .iter()
                        .any(|&angle| rounded(angle, 2) > 9_999.99) =>
            {
                Err(CellReason::ParametersOutOfRange)
            }
            Ok(parameters) => return Ok(parameters),
            Err(PeriodicCellError::PartiallyPeriodic) => Err(CellReason::PartiallyPeriodic),
            Err(_) => Err(CellReason::NonstandardOrientation),
        }
    };
    match (parameters, policy) {
        (Err(reason), UnrepresentablePolicy::Fail) => {
            Err(fail(PdbWriteErrorKind::UnrepresentableCell { reason }))
        }
        (Err(_), _) => Ok(None),
        (Ok(parameters), _) => Ok(parameters),
    }
}

fn rounded(value: f64, places: i32) -> f64 {
    let scale = 10f64.powi(places);
    (value * scale).round() / scale
}

fn is_written_placeholder(parameters: [f64; 6]) -> bool {
    parameters[..3]
        .iter()
        .all(|&length| rounded(length, 3) == 1.0)
        && parameters[3..]
            .iter()
            .all(|&angle| rounded(angle, 2) == 90.0)
}

/// The identity columns of one atom record, shared by every model.
#[derive(Debug, Clone)]
struct AtomRecord {
    atom: InstanceAtomId,
    hetatm: bool,
    serial: i64,
    /// Columns 13-16, aligned.
    name: String,
    residue_name: String,
    chain: char,
    /// Columns 73-76: a chain identifier the chain column cannot hold, or
    /// empty.
    segment: String,
    sequence: String,
    insertion_code: char,
    element: String,
    charge: String,
    /// Whether a `TER` record, numbered `serial + 1`, follows this atom.
    ter: bool,
}

/// The per-structure columns of one atom record.
#[derive(Debug, Clone, Copy)]
struct SiteValues {
    position: Point3,
    occupancy: Option<f64>,
    b_factor: Option<f64>,
}

/// One residue's identity, for grouping and polymer adjacency.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum ResidueKey {
    Hierarchy(ResidueId),
    /// The synthetic residue of an instance without atom sites.
    Instance(MoleculeInstanceId),
}

/// The residue facts bond listing needs for one atom.
#[derive(Debug, Clone)]
struct AtomResidue {
    key: ResidueKey,
    polymer: bool,
    class: ResidueClass,
}

/// One atom's identity columns before validation.
struct Identity {
    residue: AtomResidue,
    name: String,
    residue_name: String,
    chain_id: String,
    sequence: Option<String>,
    insertion_code: Option<String>,
}

/// Atom records by dense index, grouped by atom name, residue name, and
/// residue number.
type RecordIdentities<'a> = HashMap<(&'a str, &'a str, &'a str), Vec<usize>>;

/// Everything about the entry that every structure shares: atom records and
/// `CONECT` listings.
#[derive(Debug)]
struct Plan<'o> {
    records: Vec<AtomRecord>,
    /// Each record's residue.
    residues: Vec<AtomResidue>,
    /// For each atom with listed bonds, its partners by dense index, each
    /// listed once per step of bond order.
    conect: BTreeMap<usize, Vec<usize>>,
    /// The residues other than water that readers apply a template to.
    templated: Vec<Templated<'o>>,
    /// The atoms with fixed implicit hydrogen counts.
    fixed_atoms: Vec<(InstanceAtomId, ImplicitHydrogens)>,
    /// The residue names of the chains `SEQRES` states, by chain.
    seqres: Vec<(char, Vec<String>)>,
}

/// Encodes `value` in `width` columns as a decimal or, beyond it, as
/// hybrid-36, which continues with base-36 numbers beginning with an
/// upper-case and then a lower-case letter. `None` when it does not fit.
fn hybrid36(value: i64, width: u32) -> Option<String> {
    let decimal = 10i64.pow(width);
    let span = 36i64.pow(width - 1);
    let digits = |mut value: i64, alphabet: &[u8; 36]| {
        let mut text = vec![b'0'; width as usize];
        for slot in text.iter_mut().rev() {
            *slot = alphabet[usize::try_from(value % 36).expect("a base-36 digit")];
            value /= 36;
        }
        String::from_utf8(text).expect("ASCII digits")
    };
    if value < decimal {
        let text = value.to_string();
        return (text.len() <= width as usize).then_some(text);
    }
    let upper = value - decimal + 10 * span;
    if upper < 36 * span {
        return Some(digits(upper, b"0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"));
    }
    let lower = value - decimal - 16 * span;
    (lower < 36 * span).then(|| digits(lower, b"0123456789abcdefghijklmnopqrstuvwxyz"))
}

/// Columns 13-16: four-character names fill them, names of two-letter
/// elements start in column 13, and others in column 14.
fn aligned_name(name: &str, element: &str) -> String {
    if name.len() >= 4 || element.len() == 2 {
        format!("{name:<4}")
    } else {
        format!(" {name:<3}")
    }
}

/// Columns 79-80, such as `2+` or `1-`; blank for a neutral atom.
fn charge_field(charge: i8) -> String {
    match charge {
        0 => String::new(),
        charge if charge > 0 => format!("{charge}+"),
        charge => format!("{}-", -i16::from(charge)),
    }
}

fn single_character(text: &str) -> Option<char> {
    let mut characters = text.chars();
    match (characters.next(), characters.next()) {
        (None, _) => Some(' '),
        (Some(character), None) => Some(character),
        _ => None,
    }
}

impl<'o> Plan<'o> {
    /// `positions` are each structure's coordinates as written, by dense
    /// index.
    fn new(
        topology: &Topology,
        positions: &[Vec<Option<Point3>>],
        options: &'o PdbWriteOptions,

        omissions: &mut Omissions,
    ) -> Result<Self, FormatError> {
        let strict = options.chemistry == UnrepresentablePolicy::Fail;
        let mut fixed_atoms = Vec::new();
        let (mut isotopes, mut radicals, mut atom_maps) = (0, 0, 0);
        let mut stereo = BTreeSet::new();
        let mut synthetic = BTreeMap::new();
        let mut records = Vec::with_capacity(topology.atom_count());
        let mut residues = Vec::with_capacity(topology.atom_count());
        let mut names = BTreeSet::new();
        let mut segments = BTreeSet::new();
        for &atom_id in topology.atom_ids() {
            let view = topology
                .atom(atom_id)
                .expect("dense atom IDs name topology atoms");
            let atom = view.atom();
            let instance = view.instance();
            let molecule = instance.definition().molecule();
            if molecule.stereo_elements().next().is_some()
                || molecule.stereo_groups().next().is_some()
            {
                if strict {
                    return Err(fail(PdbWriteErrorKind::UnsupportedStereo(instance.id())));
                }
                stereo.insert(instance.id());
            }
            let isotope = atom.isotope.map(std::num::NonZero::get);
            let fields = [
                (
                    structure_element_symbol(atom.element, isotope).is_none(),
                    AtomField::Isotope,
                ),
                (atom.radical.is_some(), AtomField::Radical),
                (atom.atom_map.is_some(), AtomField::AtomMap),
            ];
            for (present, field) in fields {
                if !present {
                    continue;
                }
                if strict {
                    return Err(fail(PdbWriteErrorKind::UnsupportedAtomField {
                        atom: atom_id,
                        field,
                    }));
                }
                match field {
                    AtomField::Isotope => isotopes += 1,
                    AtomField::Radical => radicals += 1,
                    _ => atom_maps += 1,
                }
            }
            if atom.hydrogens != ImplicitHydrogens::default() {
                fixed_atoms.push((atom_id, atom.hydrogens));
            }
            if !(-9..=9).contains(&atom.formal_charge) {
                return Err(fail(PdbWriteErrorKind::FormalChargeOutOfRange {
                    atom: atom_id,
                    charge: atom.formal_charge,
                }));
            }
            let element =
                structure_element_symbol(atom.element, atom.isotope.map(std::num::NonZero::get))
                    .unwrap_or_else(|| atom.element.symbol());
            let identity = Self::identity(topology, atom_id, &mut synthetic);
            let (chain, segment) = match single_character(&identity.chain_id) {
                Some(chain) => (chain, String::new()),
                None if options.chain_ids == PdbChainIdPolicy::SegmentId
                    && identity.chain_id.len() <= 4
                    && identity
                        .chain_id
                        .bytes()
                        .all(|byte| byte.is_ascii_graphic()) =>
                {
                    segments.insert(identity.chain_id.clone());
                    (' ', identity.chain_id.clone())
                }
                None => {
                    return Err(fail(PdbWriteErrorKind::ChainIdTooLong {
                        atom: atom_id,
                        chain_id: identity.chain_id,
                    }))
                }
            };
            if identity.residue_name.chars().count() > 3 {
                return Err(fail(PdbWriteErrorKind::ResidueNameTooLong {
                    atom: atom_id,
                    name: identity.residue_name,
                }));
            }
            let name = identity.name;
            if name.is_empty()
                || name.len() > 4
                || !name.bytes().all(|byte| byte.is_ascii_graphic())
            {
                return Err(fail(PdbWriteErrorKind::InvalidAtomName {
                    atom: atom_id,
                    name,
                }));
            }
            let insertion_code = single_character(identity.insertion_code.as_deref().unwrap_or(""))
                .ok_or_else(|| {
                    fail(PdbWriteErrorKind::InsertionCodeTooLong {
                        atom: atom_id,
                        code: identity.insertion_code.clone().unwrap_or_default(),
                    })
                })?;
            let number = identity
                .sequence
                .as_deref()
                .and_then(|text| text.trim().parse::<i64>().ok())
                .filter(|value| (MIN_RESIDUE_NUMBER..=MAX_RESIDUE_NUMBER).contains(value));
            let Some(number) = number else {
                return Err(fail(PdbWriteErrorKind::InvalidResidueNumber {
                    atom: atom_id,
                    number: identity.sequence,
                }));
            };
            if !names.insert((identity.residue.key.clone(), name.clone())) {
                return Err(fail(PdbWriteErrorKind::DuplicateAtomName {
                    atom: atom_id,
                    name,
                }));
            }
            let standard = identity.residue.polymer
                && matches!(
                    identity.residue.class,
                    ResidueClass::AminoAcid
                        | ResidueClass::DnaNucleotide
                        | ResidueClass::RnaNucleotide
                );
            records.push(AtomRecord {
                atom: atom_id,
                hetatm: !standard,
                serial: 0,
                name: aligned_name(&name, element),
                residue_name: identity.residue_name,
                chain,
                segment,
                sequence: hybrid36(number, 4).expect("validated residue number range"),
                insertion_code,
                element: element.to_ascii_uppercase(),
                charge: charge_field(atom.formal_charge),
                ter: false,
            });
            residues.push(identity.residue);
        }

        // A reader groups consecutive records into residues, so each
        // residue's atoms must be consecutive.
        let mut left = BTreeSet::new();
        for index in 1..residues.len() {
            let (previous, current) = (&residues[index - 1].key, &residues[index].key);
            if previous != current {
                left.insert(previous.clone());
                if left.contains(current) {
                    return Err(fail(PdbWriteErrorKind::NonContiguousResidue {
                        atom: records[index].atom,
                    }));
                }
            }
        }

        // `TER` follows the last `ATOM` record of each chain and takes the
        // next serial number, so ligands bonded to the polymer, written
        // after it, read as no part of it.
        let mut last_polymer_atom = BTreeMap::new();
        for (index, record) in records.iter().enumerate() {
            if !record.hetatm {
                last_polymer_atom.insert((record.chain, record.segment.as_str()), index);
            }
        }
        let last_polymer_atom = last_polymer_atom.into_values().collect::<Vec<_>>();
        for &index in &last_polymer_atom {
            records[index].ter = true;
        }
        let serials = records.len() + last_polymer_atom.len();
        if i64::try_from(serials).map_or(true, |count| count > MAX_ATOM_SERIAL) {
            return Err(fail(PdbWriteErrorKind::TooManyAtoms { serials }));
        }
        let mut serial = 0;
        for record in &mut records {
            serial += 1;
            record.serial = serial;
            if record.ter {
                serial += 1;
            }
        }

        let library = options.components().map(|library| &**library);
        let templates = residue_templates(library, topology, &records, &residues);
        let templated = templated_residues(&records, &residues, &templates);
        let ConectPlan {
            listings: conect,
            unstated,
            coordination_bonds,
        } = Self::conect(topology, positions, &records, &residues, &templates, strict)?;
        omissions.record(WriteOmission::ChainIds {
            chains: segments.len(),
        });
        omissions.record(WriteOmission::Isotopes { atoms: isotopes });
        omissions.record(WriteOmission::Radicals { atoms: radicals });
        omissions.record(WriteOmission::AtomMaps { atoms: atom_maps });
        omissions.record(WriteOmission::Stereo {
            molecules: stereo.len(),
        });
        omissions.record(WriteOmission::Bonds { bonds: unstated });
        omissions.record(WriteOmission::CoordinationBonds {
            bonds: coordination_bonds,
        });
        drop(templates);
        Ok(Self {
            records,
            residues,
            conect,
            templated,
            fixed_atoms,
            seqres: Vec::new(),
        })
    }

    /// The source identity of one atom: its hierarchy atom site, or a
    /// synthetic `MOL` or `HOH` residue, numbered in instance order, for an
    /// instance without atom sites. Missing atom names are the element
    /// symbol in upper case and the 1-based molecule-local atom number.
    fn identity(
        topology: &Topology,
        atom_id: InstanceAtomId,
        synthetic: &mut BTreeMap<MoleculeInstanceId, usize>,
    ) -> Identity {
        let view = topology
            .atom(atom_id)
            .expect("dense atom IDs name topology atoms");
        let element = view.atom().element;
        let instance = view.instance();
        let polymer = matches!(
            instance.class(),
            MoleculeClass::Protein | MoleculeClass::Dna | MoleculeClass::Rna
        );
        let generated_name = || {
            format!(
                "{}{}",
                element.symbol().to_ascii_uppercase(),
                atom_id.atom().raw() + 1
            )
        };
        let Some(site) = topology.atom_site_of(atom_id) else {
            let next = synthetic.len() + 1;
            let number = *synthetic.entry(instance.id()).or_insert(next);
            let water = instance.class() == MoleculeClass::Water;
            return Identity {
                residue: AtomResidue {
                    key: ResidueKey::Instance(instance.id()),
                    polymer: false,
                    class: if water {
                        ResidueClass::Water
                    } else {
                        ResidueClass::Other
                    },
                },
                name: generated_name(),
                residue_name: if water { "HOH" } else { "MOL" }.to_owned(),
                chain_id: String::new(),
                sequence: Some(number.to_string()),
                insertion_code: None,
            };
        };
        let residue = site.residue();
        let chain = residue.chain();
        let metadata = site.metadata();
        // A sequenced residue is a polymer's even when a chain break leaves
        // it a molecule of its own.
        let polymer = polymer || residue.label_seq_id().is_some();
        Identity {
            residue: AtomResidue {
                key: ResidueKey::Hierarchy(residue.id()),
                polymer,
                class: residue.class(),
            },
            name: metadata
                .auth_atom_id()
                .or(metadata.label_atom_id())
                .map_or_else(generated_name, str::to_owned),
            residue_name: residue
                .author_comp_id()
                .unwrap_or(residue.name())
                .to_owned(),
            chain_id: metadata
                .auth_asym_id()
                .or(chain.author_id())
                .unwrap_or(chain.label_id())
                .to_owned(),
            sequence: residue
                .author_seq_id()
                .map(str::to_owned)
                .or_else(|| residue.label_seq_id().map(|value| value.to_string())),
            insertion_code: residue.insertion_code().map(str::to_owned),
        }
    }

    /// `CONECT` lists every bond a reader cannot derive: all but bonds
    /// within one standard residue that a default template covers, the
    /// standard links between standard residues adjacent in one polymer
    /// chain that a reader rebuilds in every structure, and bonds within
    /// water.
    fn conect(
        topology: &Topology,
        positions: &[Vec<Option<Point3>>],
        records: &[AtomRecord],
        residues: &[AtomResidue],
        templates: &BTreeMap<&ResidueKey, ResidueTemplate<'_>>,

        strict: bool,
    ) -> Result<ConectPlan, FormatError> {
        let chain_positions = residue_chain_positions(topology);
        let mut named = BTreeMap::<&ResidueKey, BTreeMap<&str, usize>>::new();
        for (index, (record, residue)) in records.iter().zip(residues).enumerate() {
            named
                .entry(&residue.key)
                .or_default()
                .insert(record.name.trim(), index);
        }
        let links = LinkContext {
            topology,
            positions,
            records,
            residues,
            chain_positions: &chain_positions,
            named: &named,
        };
        let mut conect = BTreeMap::<usize, Vec<usize>>::new();
        let mut unstated = 0;
        let mut coordination_bonds = 0;
        for bond in topology.bonds() {
            let [left, right] = bond.atoms().map(|atom| atom.index().index());
            let standard = !records[left].hetatm && !records[right].hetatm;
            let same_residue = residues[left].key == residues[right].key;
            // A residue a default template covers takes its template's bonds
            // and orders; its listings only supplement them.
            let templated = same_residue && templates.contains_key(&residues[left].key);
            let derived = if same_residue {
                residues[left].class == ResidueClass::Water || (standard && templated)
            } else {
                standard
                    && adjacent(&chain_positions, &residues[left].key, &residues[right].key)
                    && links.rebuilds(left, right)
            };
            if derived {
                continue;
            }
            let elements = bond.atoms().map(|atom| atom.atom().element);
            let coordination = is_metal(elements[0]) != is_metal(elements[1]);
            if bond.order == BondOrder::Dative && coordination && !same_residue {
                // Readers record this listing as coordination, not a bond.
                coordination_bonds += 1;
            }
            let listings = match (bond.order, coordination) {
                (BondOrder::Dative, true) | (BondOrder::Single, false) => 1,
                (BondOrder::Double, false) => 2,
                (BondOrder::Triple, false) => 3,
                // A template supplies bonds a listing would misstate, such
                // as a heme's single iron-nitrogen bonds, which a listing
                // would make dative.
                _ if templated => continue,
                (order, _) => {
                    if strict {
                        return Err(fail(PdbWriteErrorKind::UnsupportedBondOrder {
                            bond: bond.id(),
                            order,
                        }));
                    }
                    unstated += 1;
                    continue;
                }
            };
            for _ in 0..listings {
                conect.entry(left).or_default().push(right);
                conect.entry(right).or_default().push(left);
            }
        }
        for partners in conect.values_mut() {
            partners.sort_unstable();
        }
        Ok(ConectPlan {
            listings: conect,
            unstated,
            coordination_bonds,
        })
    }

    fn write_atoms(&self, writer: &mut impl Write, values: &[SiteValues]) -> std::io::Result<()> {
        let optional = |value: Option<f64>| {
            value.map_or_else(|| "      ".to_owned(), |value| format!("{value:>6.2}"))
        };
        for (record, values) in self.records.iter().zip(values) {
            let line = format!(
                "{:<6}{:>5} {} {:>3} {}{:>4}{}   {:>8.3}{:>8.3}{:>8.3}{}{}      {:<4}{:>2}{:<2}",
                if record.hetatm { "HETATM" } else { "ATOM" },
                hybrid36(record.serial, 5).expect("validated serial range"),
                record.name,
                record.residue_name,
                record.chain,
                record.sequence,
                record.insertion_code,
                values.position.x,
                values.position.y,
                values.position.z,
                optional(values.occupancy),
                optional(values.b_factor),
                record.segment,
                record.element,
                record.charge,
            );
            // Blank trailing columns are left out, as readers allow.
            writeln!(writer, "{}", line.trim_end())?;
            if record.ter {
                let line = format!(
                    "TER   {:>5}      {:>3} {}{:>4}{}",
                    hybrid36(record.serial + 1, 5).expect("validated serial range"),
                    record.residue_name,
                    record.chain,
                    record.sequence,
                    record.insertion_code,
                );
                writeln!(writer, "{}", line.trim_end())?;
            }
        }
        Ok(())
    }

    /// The coordination of the reports that the entry states, as pairs of
    /// dense indices in the first report's record order; see
    /// [`PdbSource::with_reports`] for what is left out.
    fn links(
        &self,
        views: &[StructureView<'_>],
        reports: &[StructureInterpretationReport],
        omissions: &mut Omissions,
    ) -> Vec<[usize; 2]> {
        let pair = |[left, right]: [InstanceAtomId; 2]| {
            if left <= right {
                [left, right]
            } else {
                [right, left]
            }
        };
        let stated = reports
            .iter()
            .flat_map(|report| report.coordination())
            .map(|record| pair(record.atoms()))
            .collect::<BTreeSet<_>>();
        let mut matched = views
            .iter()
            .zip(reports)
            .map(|(view, report)| report.coordination_in(view.topology()).records);
        let first = matched.next().unwrap_or_default();
        let others = matched
            .map(|records| records.into_iter().map(pair).collect::<BTreeSet<_>>())
            .collect::<Vec<_>>();
        let topology = views[0].topology();
        let mut identities = None;
        let mut written = BTreeSet::new();
        let mut links = Vec::new();
        for atoms in first {
            let key = pair(atoms);
            if written.contains(&key) || !others.iter().all(|records| records.contains(&key)) {
                continue;
            }
            let indices = atoms.map(|atom| {
                topology
                    .atom_index(atom)
                    .expect("matched partners are topology atoms")
                    .index()
            });
            // The `LINK` chain column cannot hold a segment identifier.
            if indices
                .iter()
                .any(|&index| !self.records[index].segment.is_empty())
            {
                continue;
            }
            // Readers take a `LINK` for coordination only between a metal and
            // a non-metal; between other pairs it would read back as a bond.
            let metals = atoms.map(|atom| {
                topology
                    .atom(atom)
                    .is_some_and(|view| is_metal(view.atom().element))
            });
            if metals[0] == metals[1] {
                continue;
            }
            if !indices
                .iter()
                .all(|&index| self.selects_one(index, &mut identities))
            {
                continue;
            }
            written.insert(key);
            links.push(indices);
        }
        omissions.record(WriteOmission::Coordination {
            records: stated.len() - written.len(),
        });
        links
    }

    /// Whether a `LINK` partner naming the atom record at `index` selects it
    /// alone. Readers take a blank chain or insertion code in a `LINK` for
    /// any, so a residue `10` also matches residue `10A`. `identities` groups
    /// the records by atom name, residue name, and number, built on first
    /// use.
    fn selects_one<'a>(
        &'a self,
        index: usize,
        identities: &mut Option<RecordIdentities<'a>>,
    ) -> bool {
        let identity = |record: &'a AtomRecord| {
            (
                record.name.trim(),
                record.residue_name.as_str(),
                record.sequence.as_str(),
            )
        };
        let identities = identities.get_or_insert_with(|| {
            let mut groups = HashMap::<_, Vec<usize>>::new();
            for (index, record) in self.records.iter().enumerate() {
                groups.entry(identity(record)).or_default().push(index);
            }
            groups
        });
        let partner = &self.records[index];
        identities[&identity(partner)]
            .iter()
            .filter(|&&other| {
                let candidate = &self.records[other];
                (partner.chain == ' ' || candidate.chain == partner.chain)
                    && (partner.insertion_code == ' '
                        || candidate.insertion_code == partner.insertion_code)
            })
            .count()
            == 1
    }

    /// One `LINK` record in the columns of PDB format version 3.3, with the
    /// partners' distance in `values`.
    fn write_link(
        &self,
        writer: &mut impl Write,
        [left, right]: [usize; 2],
        values: &[SiteValues],
    ) -> std::io::Result<()> {
        let partner = |index: usize| {
            let record = &self.records[index];
            format!(
                "{} {:>3} {}{:>4}{}",
                record.name,
                record.residue_name,
                record.chain,
                record.sequence,
                record.insertion_code
            )
        };
        let distance = (values[left].position - values[right].position).norm();
        // Columns 74-78 hold distances up to 99.99 Å.
        let distance = if rounded(distance, 2) <= 99.99 {
            format!("{distance:>5.2}")
        } else {
            String::new()
        };
        let line = format!(
            "LINK        {}               {}  {:>6} {:>6} {distance}",
            partner(left),
            partner(right),
            "1555",
            "1555",
        );
        writeln!(writer, "{}", line.trim_end())
    }

    fn write_conect(&self, writer: &mut impl Write) -> std::io::Result<()> {
        let serial =
            |index: usize| hybrid36(self.records[index].serial, 5).expect("validated serial range");
        for (&atom, partners) in &self.conect {
            for chunk in partners.chunks(4) {
                write!(writer, "CONECT{:>5}", serial(atom))?;
                for &partner in chunk {
                    write!(writer, "{:>5}", serial(partner))?;
                }
                writeln!(writer)?;
            }
        }
        Ok(())
    }
}

/// Counts, or with `strict` fails on, residues the bundled dictionary covers
/// with bonds between named
/// heavy atoms that neither its template nor any of its protonation variants,
/// one of which a reader selects, would give it: a missing bond or a different
/// order (an extra bond is listed instead). Coordination is dative, aromatic
/// and delocalized template bonds accept either localized order, and
/// hydrogens, whose names variants change, are not compared.
/// The `CONECT` listings of an entry and the bonds they do not carry as
/// bonds.
struct ConectPlan {
    /// Each atom's listed partners by dense index, a partner once per
    /// listing.
    listings: BTreeMap<usize, Vec<usize>>,
    /// Bonds `CONECT` cannot state, left out.
    unstated: usize,
    /// Dative bonds between residues, listed as coordination that readers
    /// record rather than bond.
    coordination_bonds: usize,
}

/// The template a reader applies to each residue (see
/// `reader_template`), with the caller's `library` first, when its atoms
/// fit one. A reader takes the bonds of a residue that fits none from its
/// `CONECT` records, so all of them are listed.
fn residue_templates<'r, 'c>(
    library: Option<&'c ComponentLibrary>,
    topology: &Topology,
    records: &'r [AtomRecord],
    residues: &'r [AtomResidue],
) -> BTreeMap<&'r ResidueKey, ResidueTemplate<'c>> {
    let mut members = BTreeMap::<&ResidueKey, Vec<usize>>::new();
    for (index, residue) in residues.iter().enumerate() {
        members.entry(&residue.key).or_default().push(index);
    }
    let element = |index: usize| {
        topology
            .atom_at(TopologyAtomIndex::from_index(index))
            .expect("records follow dense atoms")
            .atom()
            .element
    };
    members
        .into_iter()
        .filter_map(|(key, atoms)| {
            let id = &records[atoms[0]].residue_name;
            let atoms = atoms
                .iter()
                .map(|&index| (records[index].name.trim(), element(index)));
            let template = reader_template(library, id, atoms)?;
            template.fits.then_some((key, template))
        })
        .collect()
}

/// One residue other than water that readers apply a template to.
#[derive(Debug)]
struct Templated<'c> {
    key: ResidueKey,
    component: &'c Component,
    variants: Vec<&'c Component>,
    /// The residue's atoms by name, as dense indices.
    atoms: BTreeMap<String, usize>,
}

/// The residues other than water that readers apply a template to, with
/// their atoms by name.
fn templated_residues<'c>(
    records: &[AtomRecord],
    residues: &[AtomResidue],
    templates: &BTreeMap<&ResidueKey, ResidueTemplate<'c>>,
) -> Vec<Templated<'c>> {
    let mut named = BTreeMap::<&ResidueKey, BTreeMap<String, usize>>::new();
    for (index, (record, residue)) in records.iter().zip(residues).enumerate() {
        if residue.class != ResidueClass::Water && templates.contains_key(&residue.key) {
            named
                .entry(&residue.key)
                .or_default()
                .insert(record.name.trim().to_owned(), index);
        }
    }
    named
        .into_iter()
        .map(|(key, atoms)| {
            let template = &templates[key];
            Templated {
                key: key.clone(),
                component: template.component,
                variants: template.variants.clone(),
                atoms,
            }
        })
        .collect()
}

/// The charge readers take as stated for an atom record: its charge
/// columns, blank for a neutral atom.
fn stated_charge(charge: i8) -> Option<i8> {
    (charge != 0).then_some(charge)
}

impl<'o> Plan<'o> {
    /// Mirrors the readers' protonation states for the templated residues,
    /// given the coordination the entry states as `LINK` records (`links`,
    /// as dense indices) and the first structure's occupancies, and checks
    /// that readers give the structure back:
    /// - fixed implicit hydrogen counts readers do not fix again fail, or
    ///   are counted as [`WriteOmission::ImplicitHydrogens`] under
    ///   [`UnrepresentablePolicy::Omit`];
    /// - residue bonds that differ from the template or variant readers
    ///   choose fail as [`PdbWriteErrorKind::DictionaryBondConflict`], or
    ///   are counted as [`WriteOmission::DictionaryBonds`].
    ///
    /// A chain whose first or last polymer residue bonds to a residue the
    /// structure lacks, so that its linking atom has the fixed hydrogens of
    /// a gap, gets `SEQRES` records with an `UNK` residue on that side.
    fn check_reader(
        &mut self,
        topology: &Topology,
        view: StructureView<'_>,
        links: &[[usize; 2]],
        options: &PdbWriteOptions,
        omissions: &mut Omissions,
    ) -> Result<(), FormatError> {
        let mut coordinated = BTreeMap::<usize, usize>::new();
        for &index in links.iter().flatten() {
            *coordinated.entry(index).or_default() += 1;
        }
        let atom = |index: usize| {
            topology
                .atom_at(TopologyAtomIndex::from_index(index))
                .expect("records follow dense atoms")
        };
        let hydrogen = |index: usize| atom(index).atom().element.atomic_number() == 1;
        // Readers number the residues of each chain up to its `TER` in
        // record order, water aside; a chain of one residue is no polymer.
        let ends = self
            .records
            .iter()
            .enumerate()
            .filter(|(_, record)| record.ter)
            .map(|(index, record)| ((record.chain, record.segment.as_str()), index))
            .collect::<BTreeMap<_, _>>();
        let mut chains = BTreeMap::<(char, &str), Vec<&ResidueKey>>::new();
        for (index, (record, residue)) in self.records.iter().zip(&self.residues).enumerate() {
            let end = ends.get(&(record.chain, record.segment.as_str()));
            if residue.class == ResidueClass::Water || end.is_none_or(|&end| index > end) {
                continue;
            }
            let chain = chains
                .entry((record.chain, record.segment.as_str()))
                .or_default();
            if chain.last() != Some(&&residue.key) {
                chain.push(&residue.key);
            }
        }
        let mut place = BTreeMap::new();
        for (&chain, keys) in &chains {
            if keys.len() > 1 {
                for (position, &key) in keys.iter().enumerate() {
                    place.insert(key, (chain, position, keys.len()));
                }
            }
        }
        let mut pads = BTreeMap::<(char, &str), [bool; 2]>::new();
        let mut states = Vec::with_capacity(self.templated.len());
        for templated in &self.templated {
            let members = templated.atoms.values().copied().collect::<BTreeSet<_>>();
            let mut heavy = BTreeMap::new();
            for (name, &index) in &templated.atoms {
                let Some((_, template_atom)) = templated.component.atom(name) else {
                    continue;
                };
                if hydrogen(index) {
                    continue;
                }
                let mut evidence = AtomEvidence {
                    coordination: coordinated.get(&index).copied().unwrap_or(0),
                    stated: stated_charge(atom(index).atom().formal_charge),
                    ..AtomEvidence::default()
                };
                let mut sites = Vec::new();
                for bond in atom(index).bonds() {
                    let Some(other) = bond
                        .atoms()
                        .into_iter()
                        .map(|atom| atom.index().index())
                        .find(|&other| other != index)
                    else {
                        continue;
                    };
                    if hydrogen(other) {
                        let isotope = atom(other).atom().isotope;
                        sites.push((
                            isotope.is_some_and(|isotope| isotope.get() == 2),
                            view.occupancy(TopologyAtomIndex::from_index(other)),
                            None,
                        ));
                    } else if !members.contains(&other) {
                        evidence.links += link_valence(bond.bond().order);
                    }
                }
                evidence.hydrogens = hydrogen_sites(sites);
                heavy.insert(template_atom.name(), evidence);
            }
            let place = place.get(&templated.key).copied();
            let evidence = ResidueEvidence {
                heavy,
                has_hydrogens: members.iter().any(|&index| hydrogen(index)),
                dictionary: templated.component.force_field().is_none(),
                polymer: place.is_some(),
            };
            let ends = place.map_or([false, false], |(_, position, length)| {
                [position > 0, position + 1 < length]
            });
            let read = |ends| {
                reader_protonation(
                    templated.component,
                    templated.variants.iter().copied(),
                    evidence.clone(),
                    ends,
                )
            };
            let mut state = read(ends);
            // A chain end whose linking atom has a gap's fixed hydrogens is
            // a gap: `SEQRES` states residues past it.
            if let Some((chain, position, length)) = place.filter(|&((chain, _), ..)| chain != ' ')
            {
                let sides = link_sides(templated.component.linking());
                for side in 0..2 {
                    let end = if side == 0 {
                        position == 0
                    } else {
                        position + 1 == length
                    };
                    let open = |state: &Option<Protonation<'_>>| {
                        sides[side].iter().any(|name| {
                            let fixed = templated.atoms.get(*name).is_some_and(|&index| {
                                matches!(atom(index).atom().hydrogens, ImplicitHydrogens::Fixed(_))
                            });
                            fixed
                                && state.as_ref().is_some_and(|state| {
                                    state
                                        .atoms
                                        .iter()
                                        .any(|atom| atom.name == *name && atom.fixed)
                                })
                        })
                    };
                    if !end || open(&state) {
                        continue;
                    }
                    let mut padded = ends;
                    padded[side] = true;
                    let candidate = read(padded);
                    if open(&candidate) {
                        pads.entry(chain).or_default()[side] = true;
                        state = candidate;
                    }
                }
            }
            states.push(state);
        }

        // Fixed hydrogen counts have no PDB field; those a reader fixes
        // again are written as they are.
        let mut refixed = BTreeMap::new();
        for (templated, state) in self.templated.iter().zip(&states) {
            let Some(state) = state else {
                continue;
            };
            for atom_state in state.atoms.iter().filter(|state| state.fixed) {
                let Some(&index) = templated.atoms.get(&atom_state.name) else {
                    continue;
                };
                let explicit = atom(index)
                    .bonds()
                    .filter(|bond| {
                        bond.atoms().iter().any(|atom| {
                            atom.index().index() != index && hydrogen(atom.index().index())
                        })
                    })
                    .count();
                refixed.insert(
                    index,
                    u8::try_from(atom_state.hydrogens.saturating_sub(explicit)).unwrap_or(u8::MAX),
                );
            }
        }
        let mut omitted_hydrogens = 0;
        for &(atom_id, hydrogens) in &self.fixed_atoms {
            let index = topology
                .atom_index(atom_id)
                .expect("dense atoms name topology atoms")
                .index();
            if refixed
                .get(&index)
                .map(|&count| ImplicitHydrogens::Fixed(count))
                == Some(hydrogens)
            {
                continue;
            }
            if options.implicit_hydrogens == UnrepresentablePolicy::Fail {
                return Err(fail(PdbWriteErrorKind::UnsupportedAtomField {
                    atom: atom_id,
                    field: AtomField::Hydrogens,
                }));
            }
            omitted_hydrogens += 1;
        }
        omissions.record(WriteOmission::ImplicitHydrogens {
            atoms: omitted_hydrogens,
        });

        let strict = options.chemistry == UnrepresentablePolicy::Fail;
        let mut conflicting = 0;
        for (templated, state) in self.templated.iter().zip(&states) {
            // Readers fail loudly on hydrogens that fit no candidate.
            let Some(state) = state else {
                continue;
            };
            let names = templated
                .atoms
                .iter()
                .map(|(name, &index)| (name.as_str(), index))
                .collect::<BTreeMap<_, _>>();
            let Some((left, right, order)) = template_conflict(topology, state, &names) else {
                continue;
            };
            if strict {
                return Err(fail(PdbWriteErrorKind::DictionaryBondConflict {
                    atoms: [self.records[left].atom, self.records[right].atom],
                    order,
                }));
            }
            conflicting += 1;
        }
        omissions.record(WriteOmission::DictionaryBonds {
            residues: conflicting,
        });

        let mut seqres = Vec::new();
        for (chain, [before, after]) in pads {
            let keys = &chains[&chain];
            let mut names = Vec::with_capacity(keys.len() + 2);
            if before {
                names.push("UNK".to_owned());
            }
            for key in keys {
                let index = self
                    .residues
                    .iter()
                    .position(|residue| &residue.key == *key)
                    .expect("chains list record residues");
                names.push(self.records[index].residue_name.clone());
            }
            if after {
                names.push("UNK".to_owned());
            }
            seqres.push((chain.0, names));
        }
        self.seqres = seqres;
        Ok(())
    }

    /// `SEQRES` records for the chains [`Self::check_reader`] pads.
    fn write_seqres(&self, writer: &mut impl Write) -> std::io::Result<()> {
        for (chain, names) in &self.seqres {
            for (serial, line) in names.chunks(13).enumerate() {
                let residues = line
                    .iter()
                    .map(|name| format!("{name:>3}"))
                    .collect::<Vec<_>>()
                    .join(" ");
                writeln!(
                    writer,
                    "SEQRES {:>3} {chain} {:>4}  {residues}",
                    serial + 1,
                    names.len()
                )?;
            }
        }
        Ok(())
    }
}

/// The first heavy-atom bond of the chosen candidate that the residue,
/// whose atoms are named by dense index, does not have with the order a
/// reader would give it.
fn template_conflict(
    topology: &Topology,
    state: &Protonation<'_>,
    atoms: &BTreeMap<&str, usize>,
) -> Option<(usize, usize, Option<BondOrder>)> {
    let template = state.candidate;
    for bond in template.bonds() {
        let [a, b] = bond.atoms().map(|atom| &template.atoms()[atom]);
        if a.element().atomic_number() == 1 || b.element().atomic_number() == 1 {
            continue;
        }
        let (Some(&left), Some(&right)) = (atoms.get(a.name()), atoms.get(b.name())) else {
            continue;
        };
        let order = topology
            .atom_at(TopologyAtomIndex::from_index(left))
            .expect("records follow dense atoms")
            .bonds()
            .find(|candidate| {
                candidate
                    .atoms()
                    .iter()
                    .any(|atom| atom.index().index() == right)
            })
            .map(|candidate| candidate.bond().order);
        let agrees = match (order, bond.metal_coordination()) {
            (None, _) => false,
            (Some(order), true) => order == BondOrder::Dative,
            (Some(order), false) => match state.bond_order(bond) {
                Some(expected) => order == expected,
                None => matches!(order, BondOrder::Single | BondOrder::Double),
            },
        };
        if !agrees {
            return Some((left, right, order));
        }
    }
    None
}

/// Each hierarchy residue's chain and position within it.
fn residue_chain_positions(topology: &Topology) -> BTreeMap<ResidueId, (usize, usize)> {
    let mut positions = BTreeMap::new();
    for (chain_index, chain) in topology.chains().enumerate() {
        for (position, residue) in chain.residues().enumerate() {
            positions.insert(residue.id(), (chain_index, position));
        }
    }
    positions
}

/// What a reader needs to decide whether it rebuilds a standard polymer
/// link between adjacent residues.
struct LinkContext<'a> {
    topology: &'a Topology,
    positions: &'a [Vec<Option<Point3>>],
    records: &'a [AtomRecord],
    residues: &'a [AtomResidue],
    chain_positions: &'a BTreeMap<ResidueId, (usize, usize)>,
    named: &'a BTreeMap<&'a ResidueKey, BTreeMap<&'a str, usize>>,
}

impl LinkContext<'_> {
    /// Whether a reader rebuilds the bond between dense atoms `a` and `b`
    /// of adjacent residues as their standard polymer link, by the reader's
    /// own rule, in every structure written.
    fn rebuilds(&self, a: usize, b: usize) -> bool {
        let position = |atom: usize| match &self.residues[atom].key {
            ResidueKey::Hierarchy(id) => self.chain_positions.get(id).map(|&(_, at)| at),
            ResidueKey::Instance(_) => None,
        };
        let (Some(at_a), Some(at_b)) = (position(a), position(b)) else {
            return false;
        };
        let (earlier, later) = if at_a < at_b { (a, b) } else { (b, a) };
        let library = ComponentLibrary::bundled();
        let linking = |atom: usize| {
            library.get(&self.records[atom].residue_name).map_or(
                crate::components::ComponentLinking::None,
                Component::linking,
            )
        };
        let Some(link) = PolymerLink::undeclared(linking(earlier), linking(later)) else {
            return false;
        };
        let named = |atom: usize| {
            let names = &self.named[&self.residues[atom].key];
            move |name: &str| names.get(name).copied()
        };
        let element = |atom: usize| {
            self.topology
                .atom_at(TopologyAtomIndex::from_index(atom))
                .expect("records follow dense atoms")
                .atom()
                .element
        };
        let outcome = standard_link(
            link,
            named(earlier),
            named(later),
            |atom: usize, partner: usize| {
                let bonded = self
                    .topology
                    .atom_at(TopologyAtomIndex::from_index(atom))
                    .expect("records follow dense atoms")
                    .bonds()
                    .any(|bond| {
                        let [x, y] = bond.atoms().map(|end| end.index().index());
                        let other = if x == atom { y } else { x };
                        other != partner
                            && self.residues[other].key != self.residues[atom].key
                            && !is_coordination(bond.order, [element(atom), element(other)])
                    });
                Ok::<_, std::convert::Infallible>(bonded)
            },
            |left, right| {
                self.positions.iter().all(|positions| {
                    matches!(
                        (positions[left], positions[right]),
                        (Some(l), Some(r))
                            if within_bonding_distance((element(left), l), (element(right), r))
                    )
                })
            },
        );
        matches!(outcome, Ok(StandardLink::Linked(pair)) if pair == [earlier, later])
    }
}

/// Each atom's position in ångström as written, to three decimal places,
/// or `None` when it has none in that unit.
fn written_positions(view: StructureView<'_>) -> Vec<Option<Point3>> {
    let positions = view.positions();
    (0..positions.len())
        .map(|index| {
            let position = positions.get(index)?;
            position.value_in(ANGSTROM).ok().map(|position| {
                Point3::new(
                    rounded(position.x, 3),
                    rounded(position.y, 3),
                    rounded(position.z, 3),
                )
            })
        })
        .collect()
}

fn adjacent(
    positions: &BTreeMap<ResidueId, (usize, usize)>,
    left: &ResidueKey,
    right: &ResidueKey,
) -> bool {
    let (ResidueKey::Hierarchy(left), ResidueKey::Hierarchy(right)) = (left, right) else {
        return false;
    };
    match (positions.get(left), positions.get(right)) {
        (Some(&(left_chain, left)), Some(&(right_chain, right))) => {
            left_chain == right_chain && left.abs_diff(right) == 1
        }
        _ => false,
    }
}

/// One structure's coordinates, occupancies, and B factors in dense order,
/// each checked against its columns after rounding.
fn structure_values(view: StructureView<'_>, plan: &Plan) -> Result<Vec<SiteValues>, FormatError> {
    let invalid = |what: &'static str, atom: InstanceAtomId| {
        move |error: crate::units::UnitError| {
            FormatError::write(
                Format::Pdb,
                DiagnosticKind::Inconsistent,
                format!("The {what} of {atom} is invalid."),
            )
            .with_source(error)
        }
    };
    let fits = |value: f64, low: f64, high: f64, places: i32| {
        (low..=high).contains(&rounded(value, places))
    };
    let mut rows = Vec::with_capacity(plan.records.len());
    for (index, record) in plan.records.iter().enumerate() {
        let dense = TopologyAtomIndex::from_index(index);
        let position = view
            .positions()
            .get(index)
            .expect("a structure has one position per topology atom")
            .value_in(ANGSTROM)
            .map_err(invalid("position", record.atom))?;
        if ![position.x, position.y, position.z]
            .into_iter()
            .all(|value| fits(value, -999.999, 9999.999, 3))
        {
            return Err(fail(PdbWriteErrorKind::CoordinateOutOfRange {
                atom: record.atom,
            }));
        }
        let occupancy = view.occupancy(dense);
        if occupancy.is_some_and(|value| !fits(value, -99.99, 999.99, 2)) {
            return Err(fail(PdbWriteErrorKind::OccupancyOutOfRange {
                atom: record.atom,
            }));
        }
        let b_factor = view
            .b_factor(dense)
            .map(|value| value.value_in(SQUARE_ANGSTROM))
            .transpose()
            .map_err(invalid("B factor", record.atom))?;
        if b_factor.is_some_and(|value| !fits(value, -99.99, 999.99, 2)) {
            return Err(fail(PdbWriteErrorKind::BFactorOutOfRange {
                atom: record.atom,
            }));
        }
        rows.push(SiteValues {
            position,
            occupancy,
            b_factor,
        });
    }
    Ok(rows)
}

#[cfg(test)]
mod tests;
