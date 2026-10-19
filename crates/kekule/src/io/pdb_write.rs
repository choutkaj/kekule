//! PDB writing of models, ensembles, and collection items.
//!
//! A written entry holds an optional `CRYST1` record, one `MODEL` per
//! ensemble member (none for a single structure), `ATOM` and `HETATM`
//! records in dense atom order with a `TER` record after each chain's
//! polymer, `CONECT` records for the bonds a reader cannot derive from the
//! Chemical Component Dictionary, and `END`.
use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::io::Write;

use crate::diagnostics::{DiagnosticKind, Format, FormatError, FormatErrorDetail, StructureField};
use crate::geometry::{PeriodicCell, Point3};
use crate::io::element_symbol::structure_element_symbol;
use crate::io::mmcif_cell::{cell_parameters, CellParameterError, CellParameters};
use crate::io::structure_interpret::is_metal;
use crate::io::write_error::{check_structure_state, Omissions};
use crate::io::{
    ComponentBondOrder, ComponentLibrary, UnrepresentablePolicy, WriteOmission, WriteReport,
};
use crate::molecule::{BondOrder, ImplicitHydrogens};
use crate::structure::{AsModelView, Ensemble, Model, ModelView};
use crate::topology::{
    InstanceAtomId, InstanceBondId, MoleculeClass, MoleculeInstanceId, ResidueClass, ResidueId,
    Topology, TopologyAtomIndex,
};
use crate::units::{ANGSTROM, SQUARE_ANGSTROM};

/// The largest atom serial number five hybrid-36 columns hold.
const MAX_ATOM_SERIAL: i64 = 87_440_031;
/// The smallest residue sequence number four columns hold.
const MIN_RESIDUE_NUMBER: i64 = -999;
/// The largest residue sequence number four hybrid-36 columns hold.
const MAX_RESIDUE_NUMBER: i64 = 2_436_111;

/// Options of the PDB writer.
///
/// The default fails on every state PDB cannot carry. Identifiers that do
/// not fit their columns always fail: omitting or shortening them would
/// change which atoms a reader identifies.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PdbWriteOptions {
    ensemble_metadata: UnrepresentablePolicy,
    frame_metadata: UnrepresentablePolicy,
    cell: UnrepresentablePolicy,
    implicit_hydrogens: UnrepresentablePolicy,
    chemistry: UnrepresentablePolicy,
}

impl PdbWriteOptions {
    /// Ensemble member weights and energies, which PDB cannot carry: by
    /// default an ensemble with unequal weights or any energy fails to write,
    /// and so does a single ensemble member, whose weight is relative to
    /// members the entry does not hold.
    pub const fn ensemble_metadata(&self) -> UnrepresentablePolicy {
        self.ensemble_metadata
    }

    /// Sets [`Self::ensemble_metadata`].
    #[must_use]
    pub const fn with_ensemble_metadata(
        mut self,
        ensemble_metadata: UnrepresentablePolicy,
    ) -> Self {
        self.ensemble_metadata = ensemble_metadata;
        self
    }

    /// A trajectory frame's time, step, velocities, and forces, which PDB
    /// cannot carry: by default a frame with any of them fails to write.
    pub const fn frame_metadata(&self) -> UnrepresentablePolicy {
        self.frame_metadata
    }

    /// Sets [`Self::frame_metadata`].
    #[must_use]
    pub const fn with_frame_metadata(mut self, frame_metadata: UnrepresentablePolicy) -> Self {
        self.frame_metadata = frame_metadata;
        self
    }

    /// A periodic cell `CRYST1` cannot carry as the model's cell: one not
    /// periodic along all three axes, not in the standard crystallographic
    /// orientation, differing between ensemble members, or written as PDB's
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
    /// implicit hydrogens from valence. By default an atom with one fails to
    /// write; [`UnrepresentablePolicy::Omit`] writes it without the count
    /// and reports [`WriteOmission::ImplicitHydrogens`]. Structures whose
    /// files list every atom, such as XYZ interpretations, fix their counts
    /// at zero.
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
    /// `CONECT` cannot state (see [`PdbWriteIssue::UnsupportedBondOrder`]),
    /// and residue bonds that differ from the dictionary readers apply (see
    /// [`PdbWriteIssue::DictionaryBondConflict`]). By default it fails to
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
}

/// Atom state named by [`PdbWriteIssue::UnsupportedAtomField`].
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PdbAtomField {
    /// An isotope other than deuterium (`D`) or tritium (`T`).
    Isotope,
    /// Radical electrons.
    Radical,
    /// A fixed implicit hydrogen count; PDB readers infer implicit
    /// hydrogens from valence.
    Hydrogens,
    /// An atom-map number.
    AtomMap,
}

impl fmt::Display for PdbAtomField {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Isotope => "isotope",
            Self::Radical => "radical",
            Self::Hydrogens => "fixed implicit hydrogen count",
            Self::AtomMap => "atom map",
        })
    }
}

/// Why `CRYST1` cannot carry a periodic cell as the model's cell.
#[non_exhaustive]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PdbCellReason {
    /// Ensemble members have different cells; an entry has one `CRYST1`.
    MemberCellsDiffer,
    /// The cell is not periodic along all three axes.
    PartiallyPeriodic,
    /// The cell does not have `a` along x and `b` in the xy plane, as PDB
    /// coordinates assume.
    NonstandardOrientation,
    /// The written parameters are PDB's 1 Å placeholder cube, which readers
    /// take to mean no crystal lattice.
    Placeholder,
    /// A length or angle does not fit its `CRYST1` columns.
    ParametersOutOfRange,
}

impl fmt::Display for PdbCellReason {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::MemberCellsDiffer => "ensemble members have different cells",
            Self::PartiallyPeriodic => "a CRYST1 cell is periodic along all three axes",
            Self::NonstandardOrientation => {
                "a CRYST1 cell has a along x and b in the xy plane, as its coordinates assume"
            }
            Self::Placeholder => "the cell would be written as the 1 Å placeholder cube",
            Self::ParametersOutOfRange => "a cell parameter does not fit its CRYST1 columns",
        })
    }
}

/// A failure specific to the PDB writer: the typed
/// [`detail`](crate::diagnostics::FormatError::detail) of a write-stage
/// [`FormatError`], as
/// [`FormatErrorDetail::PdbWrite`](crate::diagnostics::FormatErrorDetail::PdbWrite).
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PdbWriteIssue {
    /// The chain identifier is longer than PDB's one column.
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
    /// [`PdbAtomField::Hydrogens`]) or [`PdbWriteOptions::chemistry`] set to
    /// fail.
    UnsupportedAtomField {
        /// The atom.
        atom: InstanceAtomId,
        /// The field.
        field: PdbAtomField,
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
    /// The 1-based ensemble member `member` has a weight or energy PDB cannot
    /// carry, with [`PdbWriteOptions::ensemble_metadata`] set to fail.
    UnrepresentableEnsembleMember {
        /// The 1-based member number.
        member: usize,
        /// [`StructureField::Weight`] or [`StructureField::Energy`].
        field: StructureField,
    },
    /// A periodic cell `CRYST1` cannot carry, with
    /// [`PdbWriteOptions::cell`] set to fail.
    UnrepresentableCell {
        /// Why the cell cannot be written.
        reason: PdbCellReason,
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
    /// The ensemble has no members.
    EmptyEnsemble,
}

impl PdbWriteIssue {
    /// The diagnostic kind of a write error with this detail:
    /// [`DiagnosticKind::Inconsistent`] for an empty ensemble and
    /// [`DiagnosticKind::Unsupported`] otherwise.
    pub const fn kind(&self) -> DiagnosticKind {
        match self {
            Self::EmptyEnsemble => DiagnosticKind::Inconsistent,
            _ => DiagnosticKind::Unsupported,
        }
    }
}

impl fmt::Display for PdbWriteIssue {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ChainIdTooLong { atom, chain_id } => write!(
                f,
                "chain `{chain_id}` of {atom} is longer than PDB's one-character chain identifier"
            ),
            Self::ResidueNameTooLong { atom, name } => write!(
                f,
                "residue name `{name}` of {atom} is longer than PDB's three characters"
            ),
            Self::InvalidAtomName { atom, name } => write!(
                f,
                "atom name `{name}` of {atom} does not fit PDB's four atom-name columns"
            ),
            Self::DuplicateAtomName { atom, name } => {
                write!(f, "{atom} repeats atom name `{name}` within its residue")
            }
            Self::InsertionCodeTooLong { atom, code } => write!(
                f,
                "insertion code `{code}` of {atom} is longer than PDB's one character"
            ),
            Self::InvalidResidueNumber { atom, number } => match number {
                Some(number) => write!(
                    f,
                    "residue number `{number}` of {atom} is not a PDB residue sequence number"
                ),
                None => write!(f, "the residue of {atom} has no sequence number"),
            },
            Self::NonContiguousResidue { atom } => write!(
                f,
                "{atom} resumes a residue whose atoms are not consecutive in dense order"
            ),
            Self::TooManyAtoms { serials } => write!(
                f,
                "the entry needs {serials} serial numbers, more than PDB's {MAX_ATOM_SERIAL}"
            ),
            Self::CoordinateOutOfRange { atom } => {
                write!(f, "a coordinate of {atom} does not fit PDB's columns")
            }
            Self::OccupancyOutOfRange { atom } => {
                write!(f, "the occupancy of {atom} does not fit PDB's columns")
            }
            Self::BFactorOutOfRange { atom } => {
                write!(f, "the B factor of {atom} does not fit PDB's columns")
            }
            Self::FormalChargeOutOfRange { atom, charge } => write!(
                f,
                "{atom} formal charge {charge} is outside PDB's range -9..=9"
            ),
            Self::UnsupportedAtomField { atom, field } => {
                write!(f, "PDB cannot represent the {field} of {atom}")
            }
            Self::UnsupportedStereo(molecule) => {
                write!(f, "PDB cannot represent the stereo assertions of {molecule}")
            }
            Self::UnsupportedBondOrder { bond, order } => {
                write!(f, "CONECT records cannot state {bond} as {order:?}")
            }
            Self::UnrepresentableEnsembleMember { member, field } => write!(
                f,
                "PDB cannot represent the {field} of ensemble member {member}"
            ),
            Self::UnrepresentableCell { reason } => {
                write!(f, "CRYST1 cannot represent the periodic cell: {reason}")
            }
            Self::DictionaryBondConflict { atoms, order } => write!(
                f,
                "the bond between {} and {} ({order:?}) differs from the component dictionary a PDB reader applies",
                atoms[0], atoms[1]
            ),
            Self::EmptyEnsemble => f.write_str("cannot write an empty ensemble as PDB"),
        }
    }
}

fn fail(issue: PdbWriteIssue) -> FormatError {
    FormatError::write(Format::Pdb, issue.kind(), issue.to_string())
        .with_detail(FormatErrorDetail::PdbWrite(issue))
}

/// One PDB entry to write: a structure or an ensemble.
///
/// Convert from `&Model` or `&Ensemble`; [`Self::structure`] accepts any
/// borrowed model view, such as an ensemble member or a trajectory frame,
/// whose own state fails to write unless the write options omit it.
#[derive(Debug, Clone, Copy)]
pub struct PdbSource<'a> {
    content: Content<'a>,
}

#[derive(Debug, Clone, Copy)]
enum Content<'a> {
    Structure(ModelView<'a>),
    Ensemble(&'a Ensemble),
}

impl<'a> PdbSource<'a> {
    /// One structure, written without `MODEL` records.
    pub fn structure(structure: &'a (impl AsModelView + ?Sized)) -> Self {
        Self {
            content: Content::Structure(structure.as_model_view()),
        }
    }

    /// An ensemble, one `MODEL` per member in member order.
    pub fn ensemble(ensemble: &'a Ensemble) -> Self {
        Self {
            content: Content::Ensemble(ensemble),
        }
    }
}

impl<'a> From<&'a Model> for PdbSource<'a> {
    fn from(model: &'a Model) -> Self {
        Self::structure(model)
    }
}

impl<'a> From<&'a Ensemble> for PdbSource<'a> {
    fn from(ensemble: &'a Ensemble) -> Self {
        Self::ensemble(ensemble)
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
    let (views, models) = match source.into().content {
        Content::Structure(view) => {
            check_structure_state(
                Format::Pdb,
                view.kind(),
                Some(options.ensemble_metadata),
                options.frame_metadata,
            )?;
            omissions.record_structure_state(
                view.kind(),
                Some(options.ensemble_metadata),
                options.frame_metadata,
            );
            (vec![view], false)
        }
        Content::Ensemble(ensemble) => {
            if ensemble.is_empty() {
                return Err(fail(PdbWriteIssue::EmptyEnsemble));
            }
            check_member_metadata(ensemble, options.ensemble_metadata, &mut omissions)?;
            let views = ensemble
                .iter()
                .map(|member| member.as_model_view())
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
    let plan = Plan::new(views[0].topology(), &options, &mut omissions)?;
    // Check every structure's values before the first byte is written.
    let rows = views
        .iter()
        .map(|view| structure_values(*view, &plan))
        .collect::<Result<Vec<_>, _>>()?;
    let io = |error: std::io::Error| FormatError::sink(Format::Pdb, error);
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

fn check_member_metadata(
    ensemble: &Ensemble,
    policy: UnrepresentablePolicy,
    omissions: &mut Omissions,
) -> Result<(), FormatError> {
    // Interpretation gives every member weight 1.0 and no energy, so only
    // equal weights and absent energies survive a round trip.
    let weight = ensemble.get(0).map(|member| member.weight());
    if policy != UnrepresentablePolicy::Fail {
        let unequal = ensemble
            .iter()
            .any(|member| Some(member.weight()) != weight);
        if unequal {
            omissions.record(WriteOmission::StructureState {
                field: StructureField::Weight,
                structures: ensemble.len(),
            });
        }
        omissions.record(WriteOmission::StructureState {
            field: StructureField::Energy,
            structures: ensemble
                .iter()
                .filter(|member| member.energy().is_some())
                .count(),
        });
        return Ok(());
    }
    for (index, member) in ensemble.iter().enumerate() {
        let field = if member.energy().is_some() {
            StructureField::Energy
        } else if Some(member.weight()) != weight {
            StructureField::Weight
        } else {
            continue;
        };
        return Err(fail(PdbWriteIssue::UnrepresentableEnsembleMember {
            member: index + 1,
            field,
        }));
    }
    Ok(())
}

/// The `CRYST1` parameters of an entry whose structures have these cells, or
/// none when they have none or the caller omits one PDB cannot represent.
fn entry_cell<'a>(
    cells: impl IntoIterator<Item = Option<&'a PeriodicCell>>,
    policy: UnrepresentablePolicy,
) -> Result<Option<CellParameters>, FormatError> {
    let mut cells = cells.into_iter();
    let first = cells.next().flatten();
    let parameters = if cells.any(|cell| cell != first) {
        Err(PdbCellReason::MemberCellsDiffer)
    } else {
        match first.map(cell_parameters).transpose() {
            Ok(Some(parameters)) if is_written_placeholder(parameters) => {
                Err(PdbCellReason::Placeholder)
            }
            Ok(Some(parameters))
                if parameters[..3]
                    .iter()
                    .any(|&length| rounded(length, 3) > 99_999.999)
                    || parameters[3..]
                        .iter()
                        .any(|&angle| rounded(angle, 2) > 9_999.99) =>
            {
                Err(PdbCellReason::ParametersOutOfRange)
            }
            Ok(parameters) => return Ok(parameters),
            Err(CellParameterError::PartiallyPeriodic) => Err(PdbCellReason::PartiallyPeriodic),
            Err(CellParameterError::NonstandardOrientation) => {
                Err(PdbCellReason::NonstandardOrientation)
            }
        }
    };
    match (parameters, policy) {
        (Err(reason), UnrepresentablePolicy::Fail) => {
            Err(fail(PdbWriteIssue::UnrepresentableCell { reason }))
        }
        (Err(_), _) => Ok(None),
        (Ok(parameters), _) => Ok(parameters),
    }
}

fn rounded(value: f64, places: i32) -> f64 {
    let scale = 10f64.powi(places);
    (value * scale).round() / scale
}

fn is_written_placeholder(parameters: CellParameters) -> bool {
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

/// Everything about the entry that every structure shares: atom records and
/// `CONECT` listings.
#[derive(Debug)]
struct Plan {
    records: Vec<AtomRecord>,
    /// For each atom with listed bonds, its partners by dense index, each
    /// listed once per step of bond order.
    conect: BTreeMap<usize, Vec<usize>>,
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

impl Plan {
    fn new(
        topology: &Topology,
        options: &PdbWriteOptions,
        omissions: &mut Omissions,
    ) -> Result<Self, FormatError> {
        let strict = options.chemistry == UnrepresentablePolicy::Fail;
        let mut omitted_hydrogens = 0;
        let (mut isotopes, mut radicals, mut atom_maps) = (0, 0, 0);
        let mut stereo = BTreeSet::new();
        let mut synthetic = BTreeMap::new();
        let mut records = Vec::with_capacity(topology.atom_count());
        let mut residues = Vec::with_capacity(topology.atom_count());
        let mut names = BTreeSet::new();
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
                    return Err(fail(PdbWriteIssue::UnsupportedStereo(instance.id())));
                }
                stereo.insert(instance.id());
            }
            let isotope = atom.isotope.map(std::num::NonZero::get);
            let fields = [
                (
                    structure_element_symbol(atom.element, isotope).is_none(),
                    PdbAtomField::Isotope,
                ),
                (atom.radical.is_some(), PdbAtomField::Radical),
                (atom.atom_map.is_some(), PdbAtomField::AtomMap),
            ];
            for (present, field) in fields {
                if !present {
                    continue;
                }
                if strict {
                    return Err(fail(PdbWriteIssue::UnsupportedAtomField {
                        atom: atom_id,
                        field,
                    }));
                }
                match field {
                    PdbAtomField::Isotope => isotopes += 1,
                    PdbAtomField::Radical => radicals += 1,
                    _ => atom_maps += 1,
                }
            }
            if atom.hydrogens != ImplicitHydrogens::default() {
                if options.implicit_hydrogens == UnrepresentablePolicy::Fail {
                    return Err(fail(PdbWriteIssue::UnsupportedAtomField {
                        atom: atom_id,
                        field: PdbAtomField::Hydrogens,
                    }));
                }
                omitted_hydrogens += 1;
            }
            if !(-9..=9).contains(&atom.formal_charge) {
                return Err(fail(PdbWriteIssue::FormalChargeOutOfRange {
                    atom: atom_id,
                    charge: atom.formal_charge,
                }));
            }
            let element =
                structure_element_symbol(atom.element, atom.isotope.map(std::num::NonZero::get))
                    .unwrap_or_else(|| atom.element.symbol());
            let identity = Self::identity(topology, atom_id, &mut synthetic);
            let chain = single_character(&identity.chain_id).ok_or_else(|| {
                fail(PdbWriteIssue::ChainIdTooLong {
                    atom: atom_id,
                    chain_id: identity.chain_id.clone(),
                })
            })?;
            if identity.residue_name.chars().count() > 3 {
                return Err(fail(PdbWriteIssue::ResidueNameTooLong {
                    atom: atom_id,
                    name: identity.residue_name,
                }));
            }
            let name = identity.name;
            if name.is_empty()
                || name.len() > 4
                || !name.bytes().all(|byte| byte.is_ascii_graphic())
            {
                return Err(fail(PdbWriteIssue::InvalidAtomName {
                    atom: atom_id,
                    name,
                }));
            }
            let insertion_code = single_character(identity.insertion_code.as_deref().unwrap_or(""))
                .ok_or_else(|| {
                    fail(PdbWriteIssue::InsertionCodeTooLong {
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
                return Err(fail(PdbWriteIssue::InvalidResidueNumber {
                    atom: atom_id,
                    number: identity.sequence,
                }));
            };
            if !names.insert((identity.residue.key.clone(), name.clone())) {
                return Err(fail(PdbWriteIssue::DuplicateAtomName {
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
                    return Err(fail(PdbWriteIssue::NonContiguousResidue {
                        atom: records[index].atom,
                    }));
                }
            }
        }

        // `TER` follows the last polymer atom of each chain and takes the
        // next serial number.
        let mut last_polymer_atom = BTreeMap::new();
        for (index, (residue, record)) in residues.iter().zip(&records).enumerate() {
            if residue.polymer {
                last_polymer_atom.insert(record.chain, index);
            }
        }
        for &index in last_polymer_atom.values() {
            records[index].ter = true;
        }
        let serials = records.len() + last_polymer_atom.len();
        if i64::try_from(serials).map_or(true, |count| count > MAX_ATOM_SERIAL) {
            return Err(fail(PdbWriteIssue::TooManyAtoms { serials }));
        }
        let mut serial = 0;
        for record in &mut records {
            serial += 1;
            record.serial = serial;
            if record.ter {
                serial += 1;
            }
        }

        let conflicting = check_dictionary_bonds(topology, &records, &residues, strict)?;
        let (conect, unstated) = Self::conect(topology, &records, &residues, strict)?;
        omissions.record(WriteOmission::ImplicitHydrogens {
            atoms: omitted_hydrogens,
        });
        omissions.record(WriteOmission::Isotopes { atoms: isotopes });
        omissions.record(WriteOmission::Radicals { atoms: radicals });
        omissions.record(WriteOmission::AtomMaps { atoms: atom_maps });
        omissions.record(WriteOmission::Stereo {
            molecules: stereo.len(),
        });
        omissions.record(WriteOmission::Bonds { bonds: unstated });
        omissions.record(WriteOmission::DictionaryBonds {
            residues: conflicting,
        });
        Ok(Self { records, conect })
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
        let hierarchy = topology.hierarchy();
        let Some(site) = hierarchy.atom_site_for_atom(atom_id) else {
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
        let residue = hierarchy
            .require_residue(site.residue())
            .expect("published atom sites reference live residues");
        let chain = hierarchy
            .require_chain(residue.chain())
            .expect("published residues reference live chains");
        let metadata = site.metadata();
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
    /// within one standard residue that the bundled component dictionary
    /// covers, linkage between standard residues adjacent in one polymer
    /// chain, and bonds within water.
    fn conect(
        topology: &Topology,
        records: &[AtomRecord],
        residues: &[AtomResidue],
        strict: bool,
    ) -> Result<(BTreeMap<usize, Vec<usize>>, usize), FormatError> {
        let chain_positions = residue_chain_positions(topology);
        let mut conect = BTreeMap::<usize, Vec<usize>>::new();
        let mut unstated = 0;
        for bond in topology.bonds() {
            let [left, right] = bond.atoms().map(|atom| atom.index().index());
            let standard = !records[left].hetatm && !records[right].hetatm;
            let same_residue = residues[left].key == residues[right].key;
            // A residue the bundled dictionary covers takes its template's
            // bonds and orders; its listings only supplement them.
            let templated = same_residue
                && ComponentLibrary::bundled()
                    .get(&records[left].residue_name)
                    .is_some();
            let derived = if same_residue {
                residues[left].class == ResidueClass::Water || (standard && templated)
            } else {
                standard && adjacent(&chain_positions, &residues[left].key, &residues[right].key)
            };
            if derived {
                continue;
            }
            let elements = bond.atoms().map(|atom| atom.atom().element);
            let coordination = is_metal(elements[0]) != is_metal(elements[1]);
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
                        return Err(fail(PdbWriteIssue::UnsupportedBondOrder {
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
        Ok((conect, unstated))
    }

    fn write_atoms(&self, writer: &mut impl Write, values: &[SiteValues]) -> std::io::Result<()> {
        let optional = |value: Option<f64>| {
            value.map_or_else(|| "      ".to_owned(), |value| format!("{value:>6.2}"))
        };
        for (record, values) in self.records.iter().zip(values) {
            let line = format!(
                "{:<6}{:>5} {} {:>3} {}{:>4}{}   {:>8.3}{:>8.3}{:>8.3}{}{}          {:>2}{:<2}",
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
fn check_dictionary_bonds(
    topology: &Topology,
    records: &[AtomRecord],
    residues: &[AtomResidue],
    strict: bool,
) -> Result<usize, FormatError> {
    let mut conflicting = 0;
    let mut named = BTreeMap::<&ResidueKey, BTreeMap<&str, usize>>::new();
    for (index, (record, residue)) in records.iter().zip(residues).enumerate() {
        named
            .entry(&residue.key)
            .or_default()
            .insert(record.name.trim(), index);
    }
    let library = ComponentLibrary::bundled();
    for atoms in named.values() {
        let Some(&first) = atoms.values().next() else {
            continue;
        };
        let residue_name = &records[first].residue_name;
        if residues[first].class == ResidueClass::Water {
            continue;
        }
        let Some(template) = library.get(residue_name) else {
            continue;
        };
        let Some((left, right, order)) = template_conflict(topology, template, atoms) else {
            continue;
        };
        if library
            .variants(residue_name)
            .any(|variant| template_conflict(topology, variant, atoms).is_none())
        {
            continue;
        }
        if strict {
            return Err(fail(PdbWriteIssue::DictionaryBondConflict {
                atoms: [records[left].atom, records[right].atom],
                order,
            }));
        }
        conflicting += 1;
    }
    Ok(conflicting)
}

/// The first heavy-atom bond of `template` that the residue, whose atoms are
/// named by dense index, does not have with the order a reader would give it.
fn template_conflict(
    topology: &Topology,
    template: &crate::io::Component,
    atoms: &BTreeMap<&str, usize>,
) -> Option<(usize, usize, Option<BondOrder>)> {
    for bond in template.bonds() {
        let [a, b] = bond.atoms().map(|atom| &template.atoms()[atom]);
        if a.element().atomic_number() == 1 || b.element().atomic_number() == 1 {
            continue;
        }
        let (Some(&left), Some(&right)) = (atoms.get(a.name()), atoms.get(b.name())) else {
            continue;
        };
        let order = topology
            .atom_at(TopologyAtomIndex::new(
                u32::try_from(left).expect("dense indices fit u32"),
            ))
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
            (Some(order), false) => match bond.order() {
                ComponentBondOrder::Single => order == BondOrder::Single,
                ComponentBondOrder::Double | ComponentBondOrder::Pi => order == BondOrder::Double,
                ComponentBondOrder::Triple => order == BondOrder::Triple,
                ComponentBondOrder::Quadruple => order == BondOrder::Quadruple,
                ComponentBondOrder::Aromatic | ComponentBondOrder::Delocalized => {
                    matches!(order, BondOrder::Single | BondOrder::Double)
                }
                _ => true,
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
    for (chain_index, (_, chain)) in topology.hierarchy().chains().enumerate() {
        for (position, &residue) in chain.residues().iter().enumerate() {
            positions.insert(residue, (chain_index, position));
        }
    }
    positions
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
fn structure_values(view: ModelView<'_>, plan: &Plan) -> Result<Vec<SiteValues>, FormatError> {
    let inconsistent = |error: crate::units::UnitError| {
        FormatError::write(
            Format::Pdb,
            DiagnosticKind::Inconsistent,
            "inconsistent source model",
        )
        .with_source(error)
    };
    let fits = |value: f64, low: f64, high: f64, places: i32| {
        (low..=high).contains(&rounded(value, places))
    };
    let mut rows = Vec::with_capacity(plan.records.len());
    for (index, record) in plan.records.iter().enumerate() {
        let dense = TopologyAtomIndex::new(u32::try_from(index).expect("dense indices fit u32"));
        let position = view
            .positions()
            .get(index)
            .expect("a structure has one position per topology atom")
            .value_in(ANGSTROM)
            .map_err(inconsistent)?;
        if ![position.x, position.y, position.z]
            .into_iter()
            .all(|value| fits(value, -999.999, 9999.999, 3))
        {
            return Err(fail(PdbWriteIssue::CoordinateOutOfRange {
                atom: record.atom,
            }));
        }
        let occupancy = view.occupancy(dense);
        if occupancy.is_some_and(|value| !fits(value, -99.99, 999.99, 2)) {
            return Err(fail(PdbWriteIssue::OccupancyOutOfRange {
                atom: record.atom,
            }));
        }
        let b_factor = view
            .b_factor(dense)
            .map(|value| value.value_in(SQUARE_ANGSTROM))
            .transpose()
            .map_err(inconsistent)?;
        if b_factor.is_some_and(|value| !fits(value, -99.99, 999.99, 2)) {
            return Err(fail(PdbWriteIssue::BFactorOutOfRange { atom: record.atom }));
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
mod tests {
    use super::*;

    #[test]
    fn hybrid36_encodes_what_the_reader_decodes() {
        for (value, width, text) in [
            (99_999, 5, "99999"),
            (100_000, 5, "A0000"),
            (100_001, 5, "A0001"),
            (87_440_031, 5, "zzzzz"),
            (9_999, 4, "9999"),
            (10_000, 4, "A000"),
            (-1, 4, "-1"),
            (-999, 4, "-999"),
            (2_436_111, 4, "zzzz"),
        ] {
            assert_eq!(hybrid36(value, width).as_deref(), Some(text), "{value}");
        }
        assert_eq!(hybrid36(87_440_032, 5), None);
        assert_eq!(hybrid36(-1000, 4), None);
    }

    #[test]
    fn names_and_charges_follow_the_column_conventions() {
        assert_eq!(aligned_name("CA", "C"), " CA ");
        assert_eq!(aligned_name("FE", "FE"), "FE  ");
        assert_eq!(aligned_name("HG21", "H"), "HG21");
        assert_eq!(charge_field(-2), "2-");
        assert_eq!(charge_field(1), "1+");
        assert_eq!(charge_field(0), "");
    }
}
