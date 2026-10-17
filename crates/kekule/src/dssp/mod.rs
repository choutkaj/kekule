//! DSSP 4-compatible protein secondary-structure assignment.
//!
//! Assignment is a read-only analysis of one
//! [`Model`](crate::structure::Model) coordinate snapshot.
//! It does not parse structures, choose coordinate models or alternate
//! locations, repair hierarchy data, or install derived state in the model.
//! Geometry uses single-precision coordinates and arithmetic, matching DSSP 4;
//! model coordinates are unchanged and result values are exposed as `f64`.

use crate::diagnostics::{LimitExceeded, Resource};
use std::collections::BTreeMap;
use std::fmt;

use crate::structure::AsModelView;
use crate::topology::{MoleculeInstanceId, ResidueId};

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
/// ensemble member, a trajectory frame, or a borrowed
/// [`ModelView`](crate::structure::ModelView).
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
    Loop,
    AlphaHelix,
    BetaBridge,
    ExtendedStrand,
    Helix3_10,
    PiHelix,
    PolyprolineIIHelix,
    Turn,
    Bend,
}

impl DsspSecondaryStructure {
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
    pub fn residue_name(&self) -> &str {
        &self.residue_name
    }

    pub fn chain_label_id(&self) -> &str {
        &self.chain_label_id
    }

    pub fn chain_author_id(&self) -> Option<&str> {
        self.chain_author_id.as_deref()
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
}

/// Whether a residue starts a chain segment in the analyzed snapshot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum DsspChainBreak {
    None,
    NewChain,
    Gap,
}

/// Position of a residue inside one DSSP helix pattern.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum DsspHelixPosition {
    None,
    Start,
    End,
    StartAndEnd,
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
    pub(crate) energy_kcal_per_mol: f64,
}

impl DsspHydrogenBond {
    pub const fn partner(&self) -> ResidueId {
        self.partner
    }

    /// DSSP electrostatic energy in kcal/mol.
    pub const fn energy_kcal_per_mol(&self) -> f64 {
        self.energy_kcal_per_mol
    }
}

/// Whether a beta bridge pairs strands running the same way or opposite ways.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum DsspBridgeDirection {
    Parallel,
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
    pub const fn partner(&self) -> ResidueId {
        self.partner
    }

    /// Zero-based deterministic ladder identifier.
    pub const fn ladder(&self) -> usize {
        self.ladder
    }

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
    phi_degrees: Option<f64>,
    psi_degrees: Option<f64>,
    omega_degrees: Option<f64>,
    alpha_degrees: Option<f64>,
    kappa_degrees: Option<f64>,
    tco: Option<f64>,
    acceptors: [Option<DsspHydrogenBond>; 2],
    donors: [Option<DsspHydrogenBond>; 2],
    beta_partners: [Option<DsspBetaPartner>; 2],
    sheet: Option<usize>,
    strand: Option<usize>,
    helix_positions: [DsspHelixPosition; 4],
}

impl DsspResidue {
    pub const fn key(&self) -> ResidueId {
        self.key
    }

    pub fn source(&self) -> &DsspResidueSource {
        &self.source
    }

    pub const fn secondary_structure(&self) -> DsspSecondaryStructure {
        self.secondary_structure
    }

    pub const fn chain_break(&self) -> DsspChainBreak {
        self.chain_break
    }

    pub const fn phi_degrees(&self) -> Option<f64> {
        self.phi_degrees
    }

    pub const fn psi_degrees(&self) -> Option<f64> {
        self.psi_degrees
    }

    pub const fn omega_degrees(&self) -> Option<f64> {
        self.omega_degrees
    }

    pub const fn alpha_degrees(&self) -> Option<f64> {
        self.alpha_degrees
    }

    pub const fn kappa_degrees(&self) -> Option<f64> {
        self.kappa_degrees
    }

    /// Cosine of the angle between consecutive peptide carbonyl vectors.
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

    pub const fn beta_partners(&self) -> &[Option<DsspBetaPartner>; 2] {
        &self.beta_partners
    }

    /// One-based sheet identifier, matching DSSP's sheet construction.
    pub const fn sheet(&self) -> Option<usize> {
        self.sheet
    }

    /// One-based deterministic strand identifier.
    pub const fn strand(&self) -> Option<usize> {
        self.strand
    }

    /// Positions for 3-10, alpha, pi, and polyproline-II helices.
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
    pub const fn analyzed_residues(&self) -> usize {
        self.analyzed_residues
    }

    pub const fn chain_segments(&self) -> usize {
        self.chain_segments
    }

    pub const fn hydrogen_bonds(&self) -> usize {
        self.hydrogen_bonds
    }

    pub const fn beta_bridges(&self) -> usize {
        self.beta_bridges
    }

    pub const fn intra_chain_beta_bridges(&self) -> usize {
        self.intra_chain_beta_bridges
    }

    pub const fn inter_chain_beta_bridges(&self) -> usize {
        self.inter_chain_beta_bridges
    }

    pub fn secondary_structure_count(&self, structure: DsspSecondaryStructure) -> usize {
        self.secondary_structure_counts
            .get(&structure)
            .copied()
            .unwrap_or_default()
    }

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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DsspSkippedResidue {
    pub(crate) key: ResidueId,
    pub(crate) source: DsspResidueSource,
    pub(crate) reason: DsspSkipReason,
}

impl DsspSkippedResidue {
    pub const fn key(&self) -> ResidueId {
        self.key
    }

    pub fn source(&self) -> &DsspResidueSource {
        &self.source
    }

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
    pub fn ignored_instances(&self) -> &[MoleculeInstanceId] {
        &self.ignored_instances
    }

    pub const fn non_peptide_residues(&self) -> usize {
        self.non_peptide_residues
    }

    pub fn skipped_residues(&self) -> &[DsspSkippedResidue] {
        &self.skipped_residues
    }

    pub const fn detected_gaps(&self) -> usize {
        self.detected_gaps
    }

    pub const fn reconstructed_amide_hydrogens(&self) -> usize {
        self.reconstructed_amide_hydrogens
    }

    pub const fn candidate_hydrogen_bond_pairs(&self) -> usize {
        self.candidate_hydrogen_bond_pairs
    }

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
    pub const fn max_residues(&self) -> usize {
        self.max_residues
    }

    /// Sets [`Self::max_residues`].
    #[must_use]
    pub const fn with_max_residues(mut self, max_residues: usize) -> Self {
        self.max_residues = max_residues;
        self
    }

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
    pub fn residues(&self) -> impl ExactSizeIterator<Item = &DsspResidue> {
        self.residues.iter()
    }

    pub fn get(&self, key: ResidueId) -> Option<&DsspResidue> {
        self.lookup
            .get(&key)
            .and_then(|index| self.residues.get(*index))
    }

    pub fn statistics(&self) -> &DsspStatistics {
        &self.statistics
    }

    pub fn report(&self) -> &DsspReport {
        &self.report
    }
}
