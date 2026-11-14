//! Secondary-structure assignment: the entry points, the reusable plan, its
//! options and limits, and the per-residue result, report, and error types.

use crate::diagnostics::{Cancelled, LimitExceeded};
use std::collections::BTreeMap;
use std::fmt;

use crate::parallel::Run;
use std::sync::Arc;

use crate::structure::AsStructureView;
use crate::topology::{MoleculeInstanceId, ResidueId, Topology};
use crate::units::{Quantity, DEGREE, KILOCALORIE_PER_MOLE, KILOJOULE_PER_MOLE};

/// Assign DSSP secondary structure to the protein residues in `structure`.
///
/// The returned value is a snapshot. Updating the structure coordinates does not
/// update an existing result; callers must run assignment again explicitly.
/// Amino-acid residues and otherwise unclassified residues with a complete
/// backbone are eligible without source sequence labels. Each chain uses label
/// numbering when complete, otherwise numeric author numbering when complete,
/// otherwise its topology hierarchy order. Insertion codes break numbered ties.
///
/// `structure` is any structure view: a [`Structure`](crate::structure::Structure), a
/// geometry of a structure set, a trajectory frame, or a borrowed
/// [`StructureView`](crate::structure::StructureView).
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
    structure: &(impl AsStructureView + ?Sized),
    options: DsspOptions,
) -> Result<DsspResult, DsspError> {
    assign_with(structure, options, Run::default())
}

/// [`assign`] under a [`Run`]. Each stage splits its atom sites or residues
/// into contiguous ranges on the run's executor; only the short steps
/// between stages, placing residues in chain segments and assembling beta
/// ladders and sheets, run on the calling thread. Progress counts the
/// stages done, of five.
///
/// ```
/// use kekule::dssp::{self, DsspOptions};
/// use kekule::parallel::{Executor, Run};
/// use std::num::NonZeroUsize;
///
/// /// Runs every work item on the calling thread, last first.
/// struct Reversed;
///
/// impl Executor for Reversed {
///     fn parallelism(&self) -> NonZeroUsize {
///         NonZeroUsize::new(4).unwrap()
///     }
///
///     fn run(&self, count: usize, body: &(dyn Fn(usize) + Sync)) {
///         (0..count).rev().for_each(body);
///     }
/// }
///
/// let text = std::fs::read_to_string(concat!(
///     env!("CARGO_MANIFEST_DIR"),
///     "/tests/fixtures/mmcif/1AKE.cif"
/// ))?;
/// let structure = kekule::mmcif::parse_str(&text, Default::default())?
///     .interpret(Default::default())?
///     .into_structure();
/// let options = DsspOptions::default();
/// assert_eq!(
///     dssp::assign_with(&structure, options, Run::new().with_executor(&Reversed)),
///     dssp::assign(&structure, options)
/// );
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
///
/// # Errors
///
/// As [`assign`], and [`DsspError::Cancelled`] once `run`'s cancellation is
/// requested.
pub fn assign_with(
    structure: &(impl AsStructureView + ?Sized),
    options: DsspOptions,
    run: Run<'_>,
) -> Result<DsspResult, DsspError> {
    let structure = structure.as_structure_view();
    run.check().map_err(DsspError::Cancelled)?;
    let prepared = super::kernel::Prepared::new(structure.topology(), run.executor());
    super::kernel::assign(&prepared, structure, options, run)
}

/// DSSP assignment prepared once for every structure of one topology
/// layout, such as each frame of a trajectory.
///
/// Preparing reads everything assignment needs from the topology: which
/// molecule instances the hierarchy leaves out, the chains and residues in
/// DSSP order, each residue's backbone atoms, and its source labels. Each
/// [`Self::assign`] then reads only backbone coordinates, and returns what
/// [`assign`] returns for the same structure and options, errors included.
///
/// The plan binds its topology's layout: it assigns any structure, ensemble
/// geometry, or trajectory frame whose topology
/// [shares that layout](Topology::shares_layout).
///
/// ```
/// use kekule::dssp::{self, DsspOptions, DsspPlan};
///
/// let text = std::fs::read_to_string(concat!(
///     env!("CARGO_MANIFEST_DIR"),
///     "/tests/fixtures/mmcif/1AKE.cif"
/// ))?;
/// let structure = kekule::mmcif::parse_str(&text, Default::default())?
///     .interpret(Default::default())?
///     .into_structure();
/// let plan = DsspPlan::new(structure.shared_topology(), DsspOptions::default());
/// assert_eq!(plan.assign(&structure), dssp::assign(&structure, DsspOptions::default()));
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Clone)]
pub struct DsspPlan {
    topology: Arc<Topology>,
    options: DsspOptions,
    prepared: super::kernel::Prepared,
}

impl fmt::Debug for DsspPlan {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("DsspPlan")
            .field("atom_count", &self.topology.atom_count())
            .field("peptide_residues", &self.prepared.peptide_residues())
            .field("options", &self.options)
            .finish_non_exhaustive()
    }
}

impl DsspPlan {
    /// Prepares assignment of `topology`'s structures with `options`, on
    /// the calling thread.
    pub fn new(topology: impl Into<Arc<Topology>>, options: DsspOptions) -> Self {
        let topology = topology.into();
        let prepared = super::kernel::Prepared::new(&topology, Run::new().executor());
        Self {
            topology,
            options,
            prepared,
        }
    }

    /// The topology whose layout every assigned structure must share.
    pub fn topology(&self) -> &Topology {
        &self.topology
    }

    /// The options every assignment uses.
    pub const fn options(&self) -> DsspOptions {
        self.options
    }

    /// Assigns `structure`, as [`assign`] does.
    ///
    /// # Errors
    ///
    /// [`DsspError::TopologyMismatch`] when `structure`'s topology does not
    /// share the plan's layout, otherwise the errors of [`assign`].
    pub fn assign(
        &self,
        structure: &(impl AsStructureView + ?Sized),
    ) -> Result<DsspResult, DsspError> {
        self.assign_with(structure, Run::default())
    }

    /// [`Self::assign`] under a [`Run`], as [`assign_with`] runs.
    ///
    /// # Errors
    ///
    /// As [`Self::assign`], and [`DsspError::Cancelled`] once `run`'s
    /// cancellation is requested.
    pub fn assign_with(
        &self,
        structure: &(impl AsStructureView + ?Sized),
        run: Run<'_>,
    ) -> Result<DsspResult, DsspError> {
        let structure = structure.as_structure_view();
        if !structure.topology().shares_layout(&self.topology) {
            return Err(DsspError::TopologyMismatch);
        }
        run.check().map_err(DsspError::Cancelled)?;
        super::kernel::assign(&self.prepared, structure, self.options, run)
    }
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
/// [`DsspError::InvalidCode`].
impl TryFrom<char> for DsspSecondaryStructure {
    type Error = DsspError;

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
            code => Err(DsspError::InvalidCode { code }),
        }
    }
}

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
    pub(super) key: ResidueId,
    pub(super) source: DsspResidueSource,
    pub(super) secondary_structure: DsspSecondaryStructure,
    pub(super) chain_break: DsspChainBreak,
    pub(super) phi: Option<f64>,
    pub(super) psi: Option<f64>,
    pub(super) omega: Option<f64>,
    pub(super) alpha: Option<f64>,
    pub(super) kappa: Option<f64>,
    pub(super) tco: Option<f64>,
    pub(super) acceptors: [Option<DsspHydrogenBond>; 2],
    pub(super) donors: [Option<DsspHydrogenBond>; 2],
    pub(super) beta_partners: [Option<DsspBetaPartner>; 2],
    pub(super) sheet: Option<usize>,
    pub(super) strand: Option<usize>,
    pub(super) helix_positions: [DsspHelixPosition; 4],
}

impl DsspResidue {
    /// The analyzed hierarchy residue of the input structure.
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
    pub(super) analyzed_residues: usize,
    pub(super) chain_segments: usize,
    pub(super) hydrogen_bonds: usize,
    pub(super) beta_bridges: usize,
    pub(super) intra_chain_beta_bridges: usize,
    pub(super) inter_chain_beta_bridges: usize,
    pub(super) secondary_structure_counts: BTreeMap<DsspSecondaryStructure, usize>,
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
    /// The skipped hierarchy residue of the input structure.
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
    pub(super) ignored_instances: Vec<MoleculeInstanceId>,
    pub(super) non_peptide_residues: usize,
    pub(super) skipped_residues: Vec<DsspSkippedResidue>,
    pub(super) detected_gaps: usize,
    pub(super) reconstructed_amide_hydrogens: usize,
    pub(super) candidate_hydrogen_bond_pairs: usize,
    pub(super) generated_ladders: usize,
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
    pub(super) max_residues: usize,
    pub(super) max_candidate_pairs: usize,
    pub(super) max_ladders: usize,
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
    pub(super) min_polyproline_stretch: usize,
    pub(super) limits: DsspLimits,
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
    /// A character is not a DSSP secondary-structure code.
    InvalidCode {
        /// The rejected character.
        code: char,
    },
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
    /// The structure exceeds a bound of [`DsspLimits`].
    LimitExceeded(LimitExceeded),
    /// The structure does not share the layout of the topology a
    /// [`DsspPlan`] was prepared for.
    TopologyMismatch,
    /// The run's [`Cancellation`](crate::diagnostics::Cancellation) stopped
    /// [`assign_with`] or [`DsspPlan::assign_with`].
    Cancelled(Cancelled),
}

impl fmt::Display for DsspError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidCode { code } => write!(
                f,
                "The character \"{}\" is not a DSSP secondary-structure code.",
                code.escape_debug()
            ),
            Self::InvalidPolyprolineStretch { value } => write!(
                f,
                "The polyproline-II stretch is {value}. Use 2 or 3."
            ),
            Self::AmbiguousBackboneAtom { residue, atom } => write!(
                f,
                "{} has more than one {atom} atom.",
                crate::diagnostics::Capitalized(residue)
            ),
            Self::DegenerateBackboneGeometry { residue } => write!(
                f,
                "The preceding carbonyl bond of {residue} is degenerate. Its amide hydrogen cannot be placed."
            ),
            Self::NonFiniteGeometry { residue } => write!(
                f,
                "A hydrogen-bond energy of donor {residue} is not finite."
            ),
            Self::CoordinateOutOfRange { residue } => write!(
                f,
                "A backbone coordinate of {residue} is outside the DSSP range."
            ),
            Self::NoAnalyzableProteinResidues => {
                f.write_str("No residue has a complete N, CA, C and O backbone.")
            }
            Self::LimitExceeded(limit) => {
                write!(f, "DSSP exceeded its {}.", limit.describe())
            }
            Self::TopologyMismatch => {
                f.write_str("The structure's topology layout differs from the DSSP plan's.")
            }
            Self::Cancelled(_) => f.write_str("DSSP assignment stopped."),
        }
    }
}

impl std::error::Error for DsspError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Cancelled(cancelled) => Some(cancelled),
            _ => None,
        }
    }
}

/// Complete read-only DSSP result for one coordinate snapshot.
#[derive(Debug, Clone, PartialEq)]
pub struct DsspResult {
    pub(super) residues: Vec<DsspResidue>,
    /// Each residue's key and position in `residues`, by key.
    pub(super) lookup: Vec<(ResidueId, usize)>,
    pub(super) statistics: DsspStatistics,
    pub(super) report: DsspReport,
}

impl DsspResult {
    /// Analyzed residues in DSSP order: chains by label identifier (shorter
    /// first), residues in the chain order described on [`assign`].
    pub fn residues(&self) -> impl ExactSizeIterator<Item = &DsspResidue> {
        self.residues.iter()
    }

    /// The values of one analyzed residue; `None` when `key` was not analyzed
    /// (non-peptide, skipped, or not in the structure).
    pub fn get(&self, key: ResidueId) -> Option<&DsspResidue> {
        self.lookup
            .binary_search_by_key(&key, |&(key, _)| key)
            .ok()
            .map(|found| &self.residues[self.lookup[found].1])
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
