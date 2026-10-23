//! DSSP 4-compatible protein secondary-structure assignment.
//!
//! [`assign`] is a read-only analysis of one
//! [`Model`](crate::structure::Model) coordinate snapshot that returns a
//! [`DsspResult`]: per-residue [`DsspResidue`] values (summary
//! [`DsspSecondaryStructure`], backbone torsions, hydrogen-bond and beta
//! partners, helix flags), aggregate [`DsspStatistics`], and a [`DsspReport`]
//! of skipped residues and work counts. [`DsspOptions`] and [`DsspLimits`]
//! configure the run; failures are [`DsspError`].
//!
//! Assignment does not parse structures, choose coordinate models or alternate
//! locations, repair hierarchy data, or install derived state in the model.
//! Geometry uses single-precision coordinates and arithmetic, matching DSSP 4;
//! model coordinates are unchanged and result values are exposed as `f64`.
//!
//! # Example
//!
//! ```
//! use kekule::dssp::{self, DsspChainBreak, DsspOptions, DsspSecondaryStructure};
//! use kekule::pdb;
//!
//! // A glycine dipeptide: too short for any secondary structure.
//! let entry = concat!(
//!     "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00  0.00           N\n",
//!     "ATOM      2  CA  GLY A   1       1.458   0.000   0.000  1.00  0.00           C\n",
//!     "ATOM      3  C   GLY A   1       2.009   1.420   0.000  1.00  0.00           C\n",
//!     "ATOM      4  O   GLY A   1       1.251   2.390   0.000  1.00  0.00           O\n",
//!     "ATOM      5  N   GLY A   2       3.332   1.536   0.000  1.00  0.00           N\n",
//!     "ATOM      6  CA  GLY A   2       3.988   2.836   0.000  1.00  0.00           C\n",
//!     "ATOM      7  C   GLY A   2       5.488   2.736   0.000  1.00  0.00           C\n",
//!     "ATOM      8  O   GLY A   2       6.100   1.680   0.000  1.00  0.00           O\n",
//!     "END\n",
//! );
//! let model = pdb::parse_str(entry, Default::default())?
//!     .interpret(Default::default())?
//!     .into_model();
//!
//! let result = dssp::assign(&model, DsspOptions::default())?;
//! let residues = result.residues().collect::<Vec<_>>();
//! assert_eq!(residues.len(), 2);
//! assert_eq!(residues[0].chain_break(), DsspChainBreak::NewChain);
//! assert_eq!(residues[1].chain_break(), DsspChainBreak::None);
//! assert_eq!(residues[0].secondary_structure(), DsspSecondaryStructure::Loop);
//! assert_eq!(residues[0].source().author_sequence_id(), Some("1"));
//! assert!(residues[0].phi().is_none());
//! assert_eq!(result.statistics().chain_segments(), 1);
//! assert_eq!(
//!     result.statistics().secondary_structure_count(DsspSecondaryStructure::Loop),
//!     2
//! );
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

use crate::diagnostics::{LimitExceeded, Resource};
use std::collections::BTreeMap;
use std::fmt;

use crate::structure::AsModelView;
use crate::topology::{MoleculeInstanceId, ResidueId};
use crate::units::{Quantity, DEGREE, KILOCALORIE_PER_MOLE, KILOJOULE_PER_MOLE};

mod kernel;

/// Assign DSSP secondary structure to the protein residues in `model`.
///
/// The returned value is a snapshot. Updating the model coordinates does not
/// update an existing result; callers must run assignment again explicitly.
/// Amino-acid residues and otherwise unclassified residues with a complete
/// backbone are eligible without source sequence labels. Each chain uses label
/// numbering when complete, otherwise numeric author numbering when complete,
/// otherwise its topology hierarchy order. Insertion codes break numbered ties.
///
/// `model` is any model view: a [`Model`](crate::structure::Model), an
/// ensemble sample, a trajectory frame, or a borrowed
/// [`ModelView`](crate::structure::ModelView).
///
/// # Errors
///
/// - [`DsspError::NoAnalyzableProteinResidues`] when no residue has a complete
///   `N`/`CA`/`C`/`O` backbone.
/// - [`DsspError::AmbiguousBackboneAtom`] when a residue repeats a backbone
///   atom name.
/// - [`DsspError::CoordinateOutOfRange`], [`DsspError::DegenerateBackboneGeometry`]
///   or [`DsspError::NonFiniteGeometry`] for backbone geometry DSSP cannot
///   evaluate.
/// - [`DsspError::LimitExceeded`] when the work exceeds a bound of
///   [`DsspOptions::limits`].
pub fn assign(
    model: &(impl AsModelView + ?Sized),
    options: DsspOptions,
) -> Result<DsspResult, DsspError> {
    kernel::assign(model.as_model_view(), options)
}

/// DSSP 4 summary assignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum DsspSecondaryStructure {
    /// No assigned structure, code `' '`.
    Loop,
    /// Alpha helix, `H`.
    AlphaHelix,
    /// Isolated beta bridge, `B`.
    BetaBridge,
    /// Extended strand in a beta ladder, `E`.
    ExtendedStrand,
    /// 3-10 helix, `G`.
    Helix3_10,
    /// Pi helix, `I`.
    PiHelix,
    /// Polyproline-II helix, `P`.
    PolyprolineIIHelix,
    /// Hydrogen-bonded turn, `T`.
    Turn,
    /// Bend, `S`.
    Bend,
}

impl DsspSecondaryStructure {
    /// Every assignment, in declaration order.
    pub const ALL: [Self; 9] = [
        Self::Loop,
        Self::AlphaHelix,
        Self::BetaBridge,
        Self::ExtendedStrand,
        Self::Helix3_10,
        Self::PiHelix,
        Self::PolyprolineIIHelix,
        Self::Turn,
        Self::Bend,
    ];

    /// Return the canonical one-character DSSP 4 code.
    pub const fn code(self) -> char {
        match self {
            Self::Loop => ' ',
            Self::AlphaHelix => 'H',
            Self::BetaBridge => 'B',
            Self::ExtendedStrand => 'E',
            Self::Helix3_10 => 'G',
            Self::PiHelix => 'I',
            Self::PolyprolineIIHelix => 'P',
            Self::Turn => 'T',
            Self::Bend => 'S',
        }
    }
}

/// Parses a one-character DSSP 4 code; `'-'` is accepted as an alias of
/// [`DsspSecondaryStructure::Loop`]. Any other character fails with
/// [`DsspCodeError`].
impl TryFrom<char> for DsspSecondaryStructure {
    type Error = DsspCodeError;

    fn try_from(value: char) -> Result<Self, Self::Error> {
        match value {
            ' ' | '-' => Ok(Self::Loop),
            'H' => Ok(Self::AlphaHelix),
            'B' => Ok(Self::BetaBridge),
            'E' => Ok(Self::ExtendedStrand),
            'G' => Ok(Self::Helix3_10),
            'I' => Ok(Self::PiHelix),
            'P' => Ok(Self::PolyprolineIIHelix),
            'T' => Ok(Self::Turn),
            'S' => Ok(Self::Bend),
            code => Err(DsspCodeError { code }),
        }
    }
}

/// A character is not a DSSP secondary-structure code.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DsspCodeError {
    code: char,
}

impl DsspCodeError {
    /// The rejected character.
    pub const fn code(self) -> char {
        self.code
    }
}

impl fmt::Display for DsspCodeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "unknown DSSP secondary-structure code {:?}", self.code)
    }
}

impl std::error::Error for DsspCodeError {}

/// Source hierarchy labels retained as metadata, not result identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DsspResidueSource {
    pub(crate) residue_name: String,
    pub(crate) chain_label_id: String,
    pub(crate) chain_author_id: Option<String>,
    pub(crate) label_sequence_id: Option<i32>,
    pub(crate) author_sequence_id: Option<String>,
    pub(crate) insertion_code: Option<String>,
}

impl DsspResidueSource {
    /// Residue (component) name, such as `ALA`.
    pub fn residue_name(&self) -> &str {
        &self.residue_name
    }

    /// Label chain identifier of the residue's hierarchy chain.
    pub fn chain_label_id(&self) -> &str {
        &self.chain_label_id
    }

    /// Author chain identifier; `None` when the hierarchy carries none.
    pub fn chain_author_id(&self) -> Option<&str> {
        self.chain_author_id.as_deref()
    }

    /// Label sequence number; `None` when the hierarchy carries none, as for
    /// PDB input.
    pub const fn label_sequence_id(&self) -> Option<i32> {
        self.label_sequence_id
    }

    /// Author sequence number as written; `None` when the hierarchy carries
    /// none.
    pub fn author_sequence_id(&self) -> Option<&str> {
        self.author_sequence_id.as_deref()
    }

    /// Insertion code; `None` when the residue has none.
    pub fn insertion_code(&self) -> Option<&str> {
        self.insertion_code.as_deref()
    }
}

/// Whether a residue starts a chain segment in the analyzed snapshot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum DsspChainBreak {
    /// The residue continues the segment of the preceding analyzed residue.
    None,
    /// The residue is the first analyzed residue of its chain.
    NewChain,
    /// The residue starts a new segment within its chain: the preceding
    /// residue's `C` is more than 2.5 Å from this residue's `N`, or a skipped
    /// or non-peptide residue lies between them.
    Gap,
}

/// Position of a residue inside one DSSP helix pattern.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum DsspHelixPosition {
    /// The residue is not part of this helix pattern (DSSP blank).
    None,
    /// The residue starts the pattern (DSSP `>`).
    Start,
    /// The residue ends the pattern (DSSP `<`).
    End,
    /// The residue ends one occurrence and starts another (DSSP `X`).
    StartAndEnd,
    /// The residue lies inside the pattern (DSSP digit).
    Middle,
}

/// One of the two strongest DSSP electrostatic partners.
///
/// DSSP retains the two strongest negative interactions for reporting even
/// when they do not satisfy the strict `energy < -0.5 kcal/mol` assignment
/// cutoff.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DsspHydrogenBond {
    pub(crate) partner: ResidueId,
    /// DSSP electrostatic energy in kcal/mol, the unit DSSP computes in.
    pub(crate) energy_kcal_per_mol: f64,
}

impl DsspHydrogenBond {
    /// The partner residue; it is always present in the same [`DsspResult`].
    pub const fn partner(&self) -> ResidueId {
        self.partner
    }

    /// DSSP electrostatic energy, in the canonical kJ/mol.
    pub fn energy(&self) -> Quantity<f64> {
        Quantity::new(self.energy_kcal_per_mol, KILOCALORIE_PER_MOLE)
            .into_unit(KILOJOULE_PER_MOLE)
            .expect("kcal/mol converts to kJ/mol")
    }
}

/// Whether a beta bridge pairs strands running the same way or opposite ways.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum DsspBridgeDirection {
    /// Both strands run in the same direction.
    Parallel,
    /// The strands run in opposite directions.
    Antiparallel,
}

/// One beta partner produced by DSSP bridge/ladder construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DsspBetaPartner {
    pub(crate) partner: ResidueId,
    pub(crate) ladder: usize,
    pub(crate) direction: DsspBridgeDirection,
}

impl DsspBetaPartner {
    /// The bridge partner residue; it is always present in the same
    /// [`DsspResult`].
    pub const fn partner(&self) -> ResidueId {
        self.partner
    }

    /// Zero-based deterministic ladder identifier.
    pub const fn ladder(&self) -> usize {
        self.ladder
    }

    /// Whether the bridge is parallel or antiparallel.
    pub const fn direction(&self) -> DsspBridgeDirection {
        self.direction
    }
}

/// DSSP values for one analyzed residue.
#[derive(Debug, Clone, PartialEq)]
pub struct DsspResidue {
    key: ResidueId,
    source: DsspResidueSource,
    secondary_structure: DsspSecondaryStructure,
    chain_break: DsspChainBreak,
    phi: Option<f64>,
    psi: Option<f64>,
    omega: Option<f64>,
    alpha: Option<f64>,
    kappa: Option<f64>,
    tco: Option<f64>,
    acceptors: [Option<DsspHydrogenBond>; 2],
    donors: [Option<DsspHydrogenBond>; 2],
    beta_partners: [Option<DsspBetaPartner>; 2],
    sheet: Option<usize>,
    strand: Option<usize>,
    helix_positions: [DsspHelixPosition; 4],
}

impl DsspResidue {
    /// The analyzed hierarchy residue of the input model.
    pub const fn key(&self) -> ResidueId {
        self.key
    }

    /// Source hierarchy labels of the residue.
    pub fn source(&self) -> &DsspResidueSource {
        &self.source
    }

    /// DSSP 4 summary assignment.
    pub const fn secondary_structure(&self) -> DsspSecondaryStructure {
        self.secondary_structure
    }

    /// Whether the residue starts a chain or a new segment after a gap.
    pub const fn chain_break(&self) -> DsspChainBreak {
        self.chain_break
    }

    /// Backbone torsion `C(i-1)-N-CA-C`, in degrees within (-180°, 180°].
    ///
    /// `None` for the first residue of a segment or when the torsion is
    /// geometrically undefined.
    pub fn phi(&self) -> Option<Quantity<f64>> {
        self.phi.map(|angle| Quantity::new(angle, DEGREE))
    }

    /// Backbone torsion `N-CA-C-N(i+1)`, in degrees.
    ///
    /// `None` for the last residue of a segment or when the torsion is
    /// geometrically undefined.
    pub fn psi(&self) -> Option<Quantity<f64>> {
        self.psi.map(|angle| Quantity::new(angle, DEGREE))
    }

    /// Peptide torsion `CA-C-N(i+1)-CA(i+1)`, in degrees.
    ///
    /// `None` for the last residue of a segment or when the torsion is
    /// geometrically undefined.
    pub fn omega(&self) -> Option<Quantity<f64>> {
        self.omega.map(|angle| Quantity::new(angle, DEGREE))
    }

    /// Virtual torsion `CA(i-1)-CA(i)-CA(i+1)-CA(i+2)`, in degrees.
    ///
    /// `None` unless the residue has one preceding and two following
    /// residues in its segment, or when the torsion is geometrically
    /// undefined.
    pub fn alpha(&self) -> Option<Quantity<f64>> {
        self.alpha.map(|angle| Quantity::new(angle, DEGREE))
    }

    /// Virtual bend angle between `CA(i-2)->CA(i)` and `CA(i)->CA(i+2)`, in
    /// degrees within [0°, 180°].
    ///
    /// `None` unless the residue has two residues on each side in its
    /// segment whose label sequence numbers, when both are known, span
    /// exactly four.
    pub fn kappa(&self) -> Option<Quantity<f64>> {
        self.kappa.map(|angle| Quantity::new(angle, DEGREE))
    }

    /// Cosine of the angle between this residue's and the preceding residue's
    /// `C->O` carbonyl vectors, in [-1, 1].
    ///
    /// `None` for the first residue of a segment or a zero-length carbonyl.
    pub const fn tco(&self) -> Option<f64> {
        self.tco
    }

    /// Up to two strongest acceptors for this residue's N-H donor.
    pub const fn acceptors(&self) -> &[Option<DsspHydrogenBond>; 2] {
        &self.acceptors
    }

    /// Up to two strongest N-H donors to this residue's carbonyl oxygen.
    pub const fn donors(&self) -> &[Option<DsspHydrogenBond>; 2] {
        &self.donors
    }

    /// Up to two beta-bridge partners, in DSSP's ladder order.
    pub const fn beta_partners(&self) -> &[Option<DsspBetaPartner>; 2] {
        &self.beta_partners
    }

    /// One-based sheet identifier, matching DSSP's sheet construction; `None`
    /// when the residue is in no ladder.
    pub const fn sheet(&self) -> Option<usize> {
        self.sheet
    }

    /// One-based deterministic strand identifier; `None` when the residue is
    /// in no ladder.
    pub const fn strand(&self) -> Option<usize> {
        self.strand
    }

    /// Positions in the 3-10 (`i+3`), alpha (`i+4`), pi (`i+5`), and
    /// polyproline-II patterns, in that order.
    pub const fn helix_positions(&self) -> &[DsspHelixPosition; 4] {
        &self.helix_positions
    }
}

/// Aggregate counts for one DSSP assignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DsspStatistics {
    analyzed_residues: usize,
    chain_segments: usize,
    hydrogen_bonds: usize,
    beta_bridges: usize,
    intra_chain_beta_bridges: usize,
    inter_chain_beta_bridges: usize,
    secondary_structure_counts: BTreeMap<DsspSecondaryStructure, usize>,
}

impl DsspStatistics {
    /// Number of residues in [`DsspResult::residues`].
    pub const fn analyzed_residues(&self) -> usize {
        self.analyzed_residues
    }

    /// Number of continuous backbone segments, counting every
    /// [`DsspChainBreak::NewChain`] and [`DsspChainBreak::Gap`] start.
    pub const fn chain_segments(&self) -> usize {
        self.chain_segments
    }

    /// Number of donor hydrogen bonds below DSSP's -0.5 kcal/mol cutoff.
    pub const fn hydrogen_bonds(&self) -> usize {
        self.hydrogen_bonds
    }

    /// Number of beta bridges before ladder merging.
    pub const fn beta_bridges(&self) -> usize {
        self.beta_bridges
    }

    /// Number of beta bridges between residues of one chain.
    pub const fn intra_chain_beta_bridges(&self) -> usize {
        self.intra_chain_beta_bridges
    }

    /// Number of beta bridges between residues of different chains.
    pub const fn inter_chain_beta_bridges(&self) -> usize {
        self.inter_chain_beta_bridges
    }

    /// Number of analyzed residues with this summary assignment; zero when
    /// none has it.
    pub fn secondary_structure_count(&self, structure: DsspSecondaryStructure) -> usize {
        self.secondary_structure_counts
            .get(&structure)
            .copied()
            .unwrap_or_default()
    }

    /// Residue counts for every [`DsspSecondaryStructure::ALL`] entry,
    /// including zero counts.
    pub fn secondary_structure_counts(&self) -> &BTreeMap<DsspSecondaryStructure, usize> {
        &self.secondary_structure_counts
    }
}

/// One of the four backbone atoms DSSP reads from a residue.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[non_exhaustive]
pub enum DsspBackboneAtom {
    /// The amide nitrogen, `N`.
    N,
    /// The alpha carbon, `CA`.
    Ca,
    /// The carbonyl carbon, `C`.
    C,
    /// The carbonyl oxygen, `O`.
    O,
}

impl DsspBackboneAtom {
    /// The backbone atoms in DSSP's order.
    pub const ALL: [Self; 4] = [Self::N, Self::Ca, Self::C, Self::O];

    /// The PDB atom name, such as `CA`.
    pub const fn name(self) -> &'static str {
        match self {
            Self::N => "N",
            Self::Ca => "CA",
            Self::C => "C",
            Self::O => "O",
        }
    }
}

impl fmt::Display for DsspBackboneAtom {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.name())
    }
}

/// Reason an otherwise hierarchy-visible residue was not analyzed.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum DsspSkipReason {
    /// The residue lacks these backbone atoms.
    MissingBackboneAtoms {
        /// The missing atoms, in DSSP's order.
        missing: Vec<DsspBackboneAtom>,
    },
}

/// A peptide residue that was not analyzed; it also breaks its chain segment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DsspSkippedResidue {
    pub(crate) key: ResidueId,
    pub(crate) source: DsspResidueSource,
    pub(crate) reason: DsspSkipReason,
}

impl DsspSkippedResidue {
    /// The skipped hierarchy residue of the input model.
    pub const fn key(&self) -> ResidueId {
        self.key
    }

    /// Source hierarchy labels of the residue.
    pub fn source(&self) -> &DsspResidueSource {
        &self.source
    }

    /// Why the residue was skipped.
    pub fn reason(&self) -> &DsspSkipReason {
        &self.reason
    }
}

/// Non-fatal extraction and work accounting for one assignment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DsspReport {
    ignored_instances: Vec<MoleculeInstanceId>,
    non_peptide_residues: usize,
    skipped_residues: Vec<DsspSkippedResidue>,
    detected_gaps: usize,
    reconstructed_amide_hydrogens: usize,
    candidate_hydrogen_bond_pairs: usize,
    generated_ladders: usize,
}

impl DsspReport {
    /// Molecule instances with no atom in the topology hierarchy, in instance
    /// order; their atoms are not analyzed.
    pub fn ignored_instances(&self) -> &[MoleculeInstanceId] {
        &self.ignored_instances
    }

    /// Hierarchy residues ignored as non-peptide: residues classified as
    /// neither amino acid nor other, and unclassified residues without any
    /// backbone atom.
    pub const fn non_peptide_residues(&self) -> usize {
        self.non_peptide_residues
    }

    /// Peptide residues skipped for an incomplete backbone, in analysis order.
    pub fn skipped_residues(&self) -> &[DsspSkippedResidue] {
        &self.skipped_residues
    }

    /// Number of residues marked [`DsspChainBreak::Gap`].
    pub const fn detected_gaps(&self) -> usize {
        self.detected_gaps
    }

    /// Number of amide hydrogens placed from the preceding residue's carbonyl;
    /// prolines and the first analyzed residue get none.
    pub const fn reconstructed_amide_hydrogens(&self) -> usize {
        self.reconstructed_amide_hydrogens
    }

    /// Number of residue pairs with `CA` atoms within 9 Å, the pairs whose
    /// hydrogen-bond energies were evaluated.
    pub const fn candidate_hydrogen_bond_pairs(&self) -> usize {
        self.candidate_hydrogen_bond_pairs
    }

    /// Number of provisional ladders before DSSP bulge merging.
    pub const fn generated_ladders(&self) -> usize {
        self.generated_ladders
    }
}

/// Bounded work limits for DSSP analysis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DsspLimits {
    max_residues: usize,
    max_candidate_pairs: usize,
    max_ladders: usize,
}

impl DsspLimits {
    /// Maximum residues with a complete backbone analyzed in one run.
    pub const fn max_residues(&self) -> usize {
        self.max_residues
    }

    /// Sets [`Self::max_residues`].
    #[must_use]
    pub const fn with_max_residues(mut self, max_residues: usize) -> Self {
        self.max_residues = max_residues;
        self
    }

    /// Maximum residue pairs with `CA` atoms within 9 Å evaluated for
    /// hydrogen bonds; see [`DsspReport::candidate_hydrogen_bond_pairs`].
    pub const fn max_candidate_pairs(&self) -> usize {
        self.max_candidate_pairs
    }

    /// Sets [`Self::max_candidate_pairs`].
    #[must_use]
    pub const fn with_max_candidate_pairs(mut self, max_candidate_pairs: usize) -> Self {
        self.max_candidate_pairs = max_candidate_pairs;
        self
    }

    /// Maximum provisional ladders retained before DSSP bulge merging.
    pub const fn max_ladders(&self) -> usize {
        self.max_ladders
    }

    /// Sets [`Self::max_ladders`].
    #[must_use]
    pub const fn with_max_ladders(mut self, max_ladders: usize) -> Self {
        self.max_ladders = max_ladders;
        self
    }
}

impl DsspLimits {
    /// The [`Default`] bound of [`Self::max_residues`].
    pub const DEFAULT_MAX_RESIDUES: usize = 1_000_000;
    /// The [`Default`] bound of [`Self::max_candidate_pairs`].
    pub const DEFAULT_MAX_CANDIDATE_PAIRS: usize = 50_000_000;
    /// The [`Default`] bound of [`Self::max_ladders`].
    pub const DEFAULT_MAX_LADDERS: usize = 1_000_000;
}

impl Default for DsspLimits {
    fn default() -> Self {
        Self {
            max_residues: Self::DEFAULT_MAX_RESIDUES,
            max_candidate_pairs: Self::DEFAULT_MAX_CANDIDATE_PAIRS,
            max_ladders: Self::DEFAULT_MAX_LADDERS,
        }
    }
}

/// DSSP 4 assignment options.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DsspOptions {
    min_polyproline_stretch: usize,
    limits: DsspLimits,
}

impl DsspOptions {
    /// DSSP4 accepts two or three qualifying residues; the reference default is three.
    pub const fn min_polyproline_stretch(&self) -> usize {
        self.min_polyproline_stretch
    }

    /// Sets [`Self::min_polyproline_stretch`].
    ///
    /// # Errors
    ///
    /// [`DsspError::InvalidPolyprolineStretch`] unless the stretch is 2 or 3.
    pub fn with_min_polyproline_stretch(
        mut self,
        min_polyproline_stretch: usize,
    ) -> Result<Self, DsspError> {
        if !matches!(min_polyproline_stretch, 2 | 3) {
            return Err(DsspError::InvalidPolyprolineStretch {
                value: min_polyproline_stretch,
            });
        }
        self.min_polyproline_stretch = min_polyproline_stretch;
        Ok(self)
    }

    /// Work bounds; exceeding one fails with [`DsspError::LimitExceeded`].
    pub const fn limits(&self) -> DsspLimits {
        self.limits
    }

    /// Sets [`Self::limits`].
    #[must_use]
    pub const fn with_limits(mut self, limits: DsspLimits) -> Self {
        self.limits = limits;
        self
    }
}

impl Default for DsspOptions {
    fn default() -> Self {
        Self {
            min_polyproline_stretch: 3,
            limits: DsspLimits::default(),
        }
    }
}

/// Structured DSSP assignment failure.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum DsspError {
    /// [`DsspOptions::min_polyproline_stretch`] is neither 2 nor 3.
    InvalidPolyprolineStretch {
        /// The configured value.
        value: usize,
    },
    /// A residue has more than one atom with this backbone name.
    AmbiguousBackboneAtom {
        /// The residue.
        residue: ResidueId,
        /// The repeated backbone atom.
        atom: DsspBackboneAtom,
    },
    /// The residue's preceding carbonyl bond is degenerate, so its amide
    /// hydrogen cannot be placed.
    DegenerateBackboneGeometry {
        /// The residue.
        residue: ResidueId,
    },
    /// A hydrogen-bond energy involving this donor residue is not finite.
    NonFiniteGeometry {
        /// The donor residue.
        residue: ResidueId,
    },
    /// A backbone coordinate of the residue is outside DSSP's numeric range.
    CoordinateOutOfRange {
        /// The residue.
        residue: ResidueId,
    },
    /// No residue has a complete N/CA/C/O backbone.
    NoAnalyzableProteinResidues,
    /// The model exceeds a bound of [`DsspLimits`].
    LimitExceeded(LimitExceeded),
}

impl fmt::Display for DsspError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPolyprolineStretch { value } => {
                write!(f, "DSSP polyproline-II stretch must be 2 or 3, got {value}")
            }
            Self::AmbiguousBackboneAtom { residue, atom } => {
                write!(f, "{residue} has more than one `{atom}` backbone atom")
            }
            Self::DegenerateBackboneGeometry { residue } => {
                write!(
                    f,
                    "{residue} has degenerate geometry for preceding carbonyl bond"
                )
            }
            Self::NonFiniteGeometry { residue } => {
                write!(f, "{residue} produced non-finite hydrogen-bond energy")
            }
            Self::CoordinateOutOfRange { residue } => {
                write!(
                    f,
                    "{residue} has backbone coordinate outside DSSP's numeric range"
                )
            }
            Self::NoAnalyzableProteinResidues => {
                f.write_str("model contains no residue with a complete N/CA/C/O backbone")
            }
            Self::LimitExceeded(limit) => write!(f, "DSSP {limit}"),
        }
    }
}

impl std::error::Error for DsspError {}

/// Complete read-only DSSP result for one coordinate snapshot.
#[derive(Debug, Clone, PartialEq)]
pub struct DsspResult {
    residues: Vec<DsspResidue>,
    lookup: BTreeMap<ResidueId, usize>,
    statistics: DsspStatistics,
    report: DsspReport,
}

impl DsspResult {
    /// Analyzed residues in DSSP order: chains by label identifier (shorter
    /// first), residues in the chain order described on [`assign`].
    pub fn residues(&self) -> impl ExactSizeIterator<Item = &DsspResidue> {
        self.residues.iter()
    }

    /// The values of one analyzed residue; `None` when `key` was not analyzed
    /// (non-peptide, skipped, or not in the model).
    pub fn get(&self, key: ResidueId) -> Option<&DsspResidue> {
        self.lookup
            .get(&key)
            .and_then(|index| self.residues.get(*index))
    }

    /// Aggregate counts of the assignment.
    pub fn statistics(&self) -> &DsspStatistics {
        &self.statistics
    }

    /// Non-fatal extraction issues and work counts.
    pub fn report(&self) -> &DsspReport {
        &self.report
    }
}
