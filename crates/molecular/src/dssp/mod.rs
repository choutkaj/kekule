//! DSSP 4-compatible protein secondary-structure assignment.
//!
//! Assignment is a read-only analysis of one [`Model`] coordinate snapshot.
//! It does not parse structures, choose coordinate models or alternate
//! locations, repair hierarchy data, or install derived state in the model.

use std::collections::BTreeMap;
use std::fmt;

use crate::bio::SmcraResidueId;
use crate::modeling::{Model, MoleculeInstanceId};

mod kernel;

/// Assign DSSP secondary structure to the protein residues in `model`.
///
/// The returned value is a snapshot. Updating the model coordinates does not
/// update an existing result; callers must run assignment again explicitly.
pub fn assign(model: &Model, options: DsspOptions) -> Result<DsspResult, DsspError> {
    kernel::assign(model, options)
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

/// A canonical residue identity qualified by its model molecule instance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DsspResidueKey {
    molecule: MoleculeInstanceId,
    residue: SmcraResidueId,
}

impl DsspResidueKey {
    pub const fn new(molecule: MoleculeInstanceId, residue: SmcraResidueId) -> Self {
        Self { molecule, residue }
    }

    pub const fn molecule(self) -> MoleculeInstanceId {
        self.molecule
    }

    pub const fn residue(self) -> SmcraResidueId {
        self.residue
    }
}

impl fmt::Display for DsspResidueKey {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:residue{}", self.molecule, self.residue.raw())
    }
}

/// Source hierarchy labels retained as metadata, not result identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DsspResidueSource {
    pub residue_name: String,
    pub chain_label_id: String,
    pub chain_author_id: Option<String>,
    pub label_sequence_id: Option<i32>,
    pub author_sequence_id: Option<String>,
    pub insertion_code: Option<String>,
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
    pub partner: DsspResidueKey,
    /// DSSP electrostatic energy in kcal/mol.
    pub energy_kcal_per_mol: f64,
}

/// One beta partner produced by DSSP bridge/ladder construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DsspBetaPartner {
    pub partner: DsspResidueKey,
    /// Zero-based deterministic ladder identifier.
    pub ladder: usize,
    pub parallel: bool,
}

/// DSSP values for one analyzed residue.
#[derive(Debug, Clone, PartialEq)]
pub struct DsspResidue {
    key: DsspResidueKey,
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
    pub const fn key(&self) -> DsspResidueKey {
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

/// Reason an otherwise hierarchy-visible residue was not analyzed.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum DsspSkipReason {
    MissingBackboneAtoms { missing: Vec<&'static str> },
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DsspSkippedResidue {
    pub key: DsspResidueKey,
    pub source: DsspResidueSource,
    pub reason: DsspSkipReason,
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
    pub max_residues: usize,
    pub max_candidate_pairs: usize,
    /// Maximum provisional ladders retained before DSSP bulge merging.
    pub max_ladders: usize,
}

impl Default for DsspLimits {
    fn default() -> Self {
        Self {
            max_residues: 1_000_000,
            max_candidate_pairs: 50_000_000,
            max_ladders: 1_000_000,
        }
    }
}

/// DSSP 4 assignment options.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DsspOptions {
    /// DSSP4 accepts two or three qualifying residues; the reference default is three.
    pub min_polyproline_stretch: usize,
    pub limits: DsspLimits,
}

impl Default for DsspOptions {
    fn default() -> Self {
        Self {
            min_polyproline_stretch: 3,
            limits: DsspLimits::default(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[non_exhaustive]
pub enum DsspResource {
    Residues,
    CandidatePairs,
    Ladders,
}

/// Structured DSSP assignment failure.
#[derive(Debug, Clone, PartialEq, Eq)]
#[non_exhaustive]
pub enum DsspError {
    InvalidPolyprolineStretch {
        value: usize,
    },
    AmbiguousBackboneAtom {
        residue: DsspResidueKey,
        atom_name: &'static str,
    },
    DegenerateBackboneGeometry {
        residue: DsspResidueKey,
        quantity: &'static str,
    },
    NonFiniteGeometry {
        residue: DsspResidueKey,
        quantity: &'static str,
    },
    CoordinateOutOfRange {
        residue: DsspResidueKey,
        quantity: &'static str,
    },
    InvalidHierarchy {
        molecule: MoleculeInstanceId,
        message: String,
    },
    NoAnalyzableProteinResidues,
    ResourceLimitExceeded {
        resource: DsspResource,
        limit: usize,
    },
}

impl fmt::Display for DsspError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPolyprolineStretch { value } => {
                write!(f, "DSSP polyproline-II stretch must be 2 or 3, got {value}")
            }
            Self::AmbiguousBackboneAtom { residue, atom_name } => {
                write!(f, "{residue} has more than one `{atom_name}` backbone atom")
            }
            Self::DegenerateBackboneGeometry { residue, quantity } => {
                write!(f, "{residue} has degenerate geometry for {quantity}")
            }
            Self::NonFiniteGeometry { residue, quantity } => {
                write!(f, "{residue} produced non-finite {quantity}")
            }
            Self::CoordinateOutOfRange { residue, quantity } => {
                write!(f, "{residue} has {quantity} outside DSSP's numeric range")
            }
            Self::InvalidHierarchy { molecule, message } => {
                write!(f, "invalid hierarchy in {molecule}: {message}")
            }
            Self::NoAnalyzableProteinResidues => {
                f.write_str("model contains no residue with a complete N/CA/C/O backbone")
            }
            Self::ResourceLimitExceeded { resource, limit } => {
                write!(f, "DSSP {resource:?} resource limit {limit} exceeded")
            }
        }
    }
}

impl std::error::Error for DsspError {}

/// Complete read-only DSSP result for one coordinate snapshot.
#[derive(Debug, Clone, PartialEq)]
pub struct DsspResult {
    residues: Vec<DsspResidue>,
    lookup: BTreeMap<DsspResidueKey, usize>,
    statistics: DsspStatistics,
    report: DsspReport,
}

impl DsspResult {
    pub fn residues(&self) -> impl ExactSizeIterator<Item = &DsspResidue> {
        self.residues.iter()
    }

    pub fn get(&self, key: DsspResidueKey) -> Option<&DsspResidue> {
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
