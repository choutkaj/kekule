use std::collections::{BTreeMap, VecDeque};
use std::ops::Range;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::geometry::Point3;
use crate::structure::ModelView;
use crate::topology::{Chain, Hierarchy, Residue, ResidueClass, ResidueId};
use crate::units::ANGSTROM;

use super::*;
use crate::parallel::map_ranges;

/// Fewest residues or residue pairs a stage gives one thread.
const MIN_CHUNK: usize = 256;

const MAX_PEPTIDE_BOND_ANGSTROMS: f64 = 2.5;
const MAX_CA_PAIR_DISTANCE_ANGSTROMS: f64 = 9.0;
const MIN_ELECTROSTATIC_DISTANCE_ANGSTROMS: f64 = 0.5;
const MIN_HBOND_ENERGY_KCAL_PER_MOL: f64 = -9.9;
const MAX_HBOND_ENERGY_KCAL_PER_MOL: f64 = -0.5;
const DSSP_COUPLING_KCAL_ANGSTROM_PER_MOL: f64 = -27.888_f32 as f64;

#[derive(Debug, Clone, Copy, PartialEq)]
// DSSP 4.6.1 evaluates points and angular geometry in float. Keep that
// compatibility boundary local; model coordinates and public results are f64.
struct Vec3 {
    x: f32,
    y: f32,
    z: f32,
}

impl From<Point3> for Vec3 {
    fn from(value: Point3) -> Self {
        Self {
            x: value.x as f32,
            y: value.y as f32,
            z: value.z as f32,
        }
    }
}

impl Vec3 {
    fn sub(self, other: Self) -> Self {
        Self {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    fn dot(self, other: Self) -> f32 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn cross(self, other: Self) -> Self {
        Self {
            x: self.y * other.z - other.y * self.z,
            y: self.z * other.x - other.z * self.x,
            z: self.x * other.y - other.x * self.y,
        }
    }

    fn norm_squared(self) -> f32 {
        self.dot(self)
    }

    fn is_dssp_representable(self) -> bool {
        [self.x, self.y, self.z].into_iter().all(|value| {
            let cell = (value / MAX_CA_PAIR_DISTANCE_ANGSTROMS as f32).floor();
            value.is_finite() && cell > i64::MIN as f32 && cell < i64::MAX as f32
        })
    }
}

fn dssp_distance_squared(first: Vec3, second: Vec3) -> f64 {
    let dx = first.x - second.x;
    let dy = first.y - second.y;
    let dz = first.z - second.z;
    ((dx * dx + dy * dy) + dz * dz) as f64
}

fn dssp_distance(first: Vec3, second: Vec3) -> f64 {
    (dssp_distance_squared(first, second) as f32).sqrt() as f64
}

fn dssp_reconstructed_hydrogen(n: Vec3, preceding_c: Vec3, preceding_o: Vec3) -> Option<Vec3> {
    let dx = preceding_c.x - preceding_o.x;
    let dy = preceding_c.y - preceding_o.y;
    let dz = preceding_c.z - preceding_o.z;
    let length = ((dx * dx + dy * dy) + dz * dz).sqrt();
    (length.is_finite() && length > f32::EPSILON).then(|| Vec3 {
        x: n.x + dx / length,
        y: n.y + dy / length,
        z: n.z + dz / length,
    })
}

#[derive(Debug, Clone, Copy, PartialEq)]
struct BondSlot {
    partner: Option<usize>,
    energy: f64,
}

impl Default for BondSlot {
    fn default() -> Self {
        Self {
            partner: None,
            energy: 0.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct BetaSlot {
    partner: usize,
    ladder: usize,
    parallel: bool,
}

#[derive(Debug, Clone)]
struct BackboneResidue {
    key: ResidueId,
    source: DsspResidueSource,
    chain: usize,
    n: Vec3,
    ca: Vec3,
    c: Vec3,
    o: Vec3,
    h: Option<Vec3>,
    is_proline: bool,
    segment: usize,
    prev: Option<usize>,
    next: Option<usize>,
    chain_break: DsspChainBreak,
    phi: Option<f64>,
    psi: Option<f64>,
    omega: Option<f64>,
    alpha: Option<f64>,
    kappa: Option<f64>,
    tco: Option<f64>,
    acceptors: [BondSlot; 2],
    donors: [BondSlot; 2],
    beta_partners: [Option<BetaSlot>; 2],
    sheet: Option<usize>,
    strand: Option<usize>,
    helix_positions: [DsspHelixPosition; 4],
    secondary_structure: DsspSecondaryStructure,
}

#[derive(Debug, Clone, Copy)]
struct ChainFragment<'a> {
    hierarchy: &'a Hierarchy,
    chain: &'a Chain,
    discovery_order: usize,
}

#[derive(Debug, Clone, Copy)]
struct ChainResidue<'a> {
    hierarchy: &'a Hierarchy,
    chain: &'a Chain,
    residue: ResidueId,
    discovery_order: usize,
}

#[derive(Debug)]
struct LogicalChain<'a> {
    label_id: String,
    residues: Vec<ChainResidue<'a>>,
    discovery_order: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BridgeType {
    Parallel,
    Antiparallel,
}

#[derive(Debug, Clone)]
struct Bridge {
    bridge_type: BridgeType,
    i: Vec<usize>,
    j: VecDeque<usize>,
    sheet: usize,
    ladder: usize,
}

#[derive(Debug, Clone, Copy, Default)]
struct BetaCounts {
    bridges: usize,
    intra_chain: usize,
    inter_chain: usize,
    ladders: usize,
}

pub(super) fn assign(
    model: ModelView<'_>,
    options: DsspOptions,
    executor: &dyn Executor,
) -> Result<DsspResult, DsspError> {
    let (mut residues, mut report, chain_segments) =
        extract_backbones(model, options.limits, executor)?;
    if residues.is_empty() {
        return Err(DsspError::NoAnalyzableProteinResidues);
    }

    let candidates = analyze_backbone(&mut residues, options.limits.max_candidate_pairs, executor)?;
    report.candidate_hydrogen_bond_pairs = candidates.len();
    let beta_counts = calculate_beta_sheets(
        &mut residues,
        &candidates,
        options.limits.max_ladders,
        executor,
    )?;
    report.generated_ladders = beta_counts.ladders;
    calculate_helices_turns_and_bends(&mut residues);
    calculate_polyproline_helices(&mut residues, options.min_polyproline_stretch);

    let hydrogen_bonds = residues
        .iter()
        .flat_map(|residue| residue.donors)
        .filter(|bond| bond.partner.is_some() && bond.energy < MAX_HBOND_ENERGY_KCAL_PER_MOL)
        .count();
    let mut secondary_structure_counts = DsspSecondaryStructure::ALL
        .into_iter()
        .map(|structure| (structure, 0))
        .collect::<BTreeMap<_, _>>();
    for residue in &residues {
        *secondary_structure_counts
            .entry(residue.secondary_structure)
            .or_default() += 1;
    }

    let keys = residues
        .iter()
        .map(|residue| residue.key)
        .collect::<Vec<_>>();
    let public_residues = residues
        .into_iter()
        .map(|residue| public_residue(residue, &keys))
        .collect::<Vec<_>>();
    let lookup = public_residues
        .iter()
        .enumerate()
        .map(|(index, residue)| (residue.key, index))
        .collect();
    let statistics = DsspStatistics {
        analyzed_residues: public_residues.len(),
        chain_segments,
        hydrogen_bonds,
        beta_bridges: beta_counts.bridges,
        intra_chain_beta_bridges: beta_counts.intra_chain,
        inter_chain_beta_bridges: beta_counts.inter_chain,
        secondary_structure_counts,
    };
    Ok(DsspResult {
        residues: public_residues,
        lookup,
        statistics,
        report,
    })
}

fn extract_backbones(
    model: ModelView<'_>,
    limits: DsspLimits,
    executor: &dyn Executor,
) -> Result<(Vec<BackboneResidue>, DsspReport, usize), DsspError> {
    let mut report = DsspReport {
        ignored_instances: Vec::new(),
        non_peptide_residues: 0,
        skipped_residues: Vec::new(),
        detected_gaps: 0,
        reconstructed_amide_hydrogens: 0,
        candidate_hydrogen_bond_pairs: 0,
        generated_ladders: 0,
    };
    let mut chain_segments = 0;
    let mut fragments = Vec::new();

    let hierarchy = model.topology().hierarchy();
    let mut represented_instances = hierarchy
        .atom_sites()
        .map(|(_, site)| site.atom().instance())
        .collect::<Vec<_>>();
    represented_instances.sort_unstable();
    represented_instances.dedup();
    for molecule_id in model.topology().instances().map(|molecule| molecule.id()) {
        if represented_instances.binary_search(&molecule_id).is_err() {
            report.ignored_instances.push(molecule_id);
        }
    }
    for (_, chain) in hierarchy.chains() {
        fragments.push(ChainFragment {
            hierarchy,
            chain,
            discovery_order: fragments.len(),
        });
    }

    // Hierarchy chains are already topology-global and may span molecule
    // instances, so no fragment reassembly is required here.
    let mut logical_chains = Vec::<LogicalChain<'_>>::new();
    let mut residue_discovery_order = 0usize;
    for fragment in fragments {
        logical_chains.push(LogicalChain {
            label_id: fragment.chain.label_id().to_owned(),
            residues: Vec::new(),
            discovery_order: fragment.discovery_order,
        });
        let logical = logical_chains
            .last_mut()
            .expect("just inserted logical chain");
        for &residue in fragment.chain.residues() {
            logical.residues.push(ChainResidue {
                hierarchy: fragment.hierarchy,
                chain: fragment.chain,
                residue,
                discovery_order: residue_discovery_order,
            });
            residue_discovery_order += 1;
        }
    }

    // mkdssp builds one global residue table in label-asym order. Stable
    // shortlex ordering matches generated mmCIF label-asym series (A..Z,
    // AA..) while preserving discovery order for repeated logical chains.
    logical_chains.sort_by(|left, right| {
        left.label_id
            .len()
            .cmp(&right.label_id.len())
            .then_with(|| left.label_id.cmp(&right.label_id))
            .then_with(|| left.discovery_order.cmp(&right.discovery_order))
    });
    // Ordering a chain and reading a residue's backbone depend on nothing
    // else; segments, limits and errors are then decided in DSSP order.
    let chains = map_ranges(executor, logical_chains.len(), 1, |range| {
        logical_chains[range]
            .iter()
            .map(|logical| order_chain(&logical.residues))
            .collect::<Vec<_>>()
    })
    .into_iter()
    .flatten()
    .collect::<Vec<_>>();
    let entries = chains.iter().flatten().collect::<Vec<_>>();
    let mut backbones = map_ranges(executor, entries.len(), MIN_CHUNK, |range| {
        entries[range]
            .iter()
            .map(|entry| read_backbone(model, entry))
            .collect::<Vec<_>>()
    })
    .into_iter()
    .flatten();
    let mut residues = Vec::with_capacity(entries.len().min(limits.max_residues));
    let mut degenerate = None;
    for (chain_index, chain) in chains.iter().enumerate() {
        extract_chain(
            chain_index,
            backbones.by_ref().take(chain.len()),
            &mut residues,
            &mut report,
            &mut chain_segments,
            &mut degenerate,
            limits,
        )?;
    }
    // Hydrogen placement follows extraction in DSSP, so any extraction
    // failure takes precedence.
    if let Some(residue) = degenerate {
        return Err(DsspError::DegenerateBackboneGeometry { residue });
    }
    Ok((residues, report, chain_segments))
}

/// The residues of one logical chain in DSSP order.
fn order_chain<'a>(residues: &[ChainResidue<'a>]) -> Vec<ChainResidue<'a>> {
    // Choose one ordering convention for the complete chain. Mixing label
    // and author numbering in a pairwise comparator can be non-transitive.
    // Native hierarchies need not carry either source numbering convention.
    let ordering_residues = residues
        .iter()
        .filter(|entry| {
            entry
                .hierarchy
                .residue(entry.residue)
                .is_some_and(|residue| {
                    residue.class() == ResidueClass::AminoAcid
                        || residue.class() == ResidueClass::Other
                            && residue.atom_sites().iter().any(|site| {
                                entry.hierarchy.atom_site(*site).is_some_and(|site| {
                                    site.metadata()
                                        .label_atom_id
                                        .as_deref()
                                        .or(site.metadata().auth_atom_id.as_deref())
                                        .and_then(|name| backbone_index(name.trim()))
                                        .is_some()
                                })
                            })
                })
        })
        .collect::<Vec<_>>();
    let label_numbered = ordering_residues.iter().all(|entry| {
        entry
            .hierarchy
            .residue(entry.residue)
            .is_some_and(|residue| residue.label_seq_id().is_some())
    });
    let author_numbered = ordering_residues.iter().all(|entry| {
        entry
            .hierarchy
            .residue(entry.residue)
            .is_some_and(|residue| {
                residue
                    .author_seq_id()
                    .and_then(|value| value.parse::<i64>().ok())
                    .is_some()
            })
    });
    let mut ordered = residues.to_vec();
    ordered.sort_by(|left, right| {
        let left_residue = left
            .hierarchy
            .residue(left.residue)
            .expect("validated logical chain residue");
        let right_residue = right
            .hierarchy
            .residue(right.residue)
            .expect("validated logical chain residue");
        let sequence_order = if label_numbered {
            left_residue
                .label_seq_id()
                .cmp(&right_residue.label_seq_id())
        } else if author_numbered {
            let number = |residue: &Residue| {
                residue
                    .author_seq_id()
                    .and_then(|value| value.parse::<i64>().ok())
            };
            number(left_residue).cmp(&number(right_residue))
        } else {
            return left.discovery_order.cmp(&right.discovery_order);
        };
        sequence_order
            .then_with(|| {
                left_residue
                    .insertion_code()
                    .cmp(&right_residue.insertion_code())
            })
            .then_with(|| left.discovery_order.cmp(&right.discovery_order))
    });
    ordered
}

/// The backbone atoms of one peptide residue, as its atom sites give them.
#[derive(Debug)]
struct ResidueBackbone {
    key: ResidueId,
    source: DsspResidueSource,
    /// `N`, `CA`, `C`, `O`, or the first repeated one.
    atoms: Result<[Option<Vec3>; 4], DsspBackboneAtom>,
    any_backbone_atom: bool,
    amino_acid: bool,
    is_proline: bool,
}

/// Reads a residue's backbone atoms; `None` for a residue classified as
/// neither amino acid nor other.
fn read_backbone(model: ModelView<'_>, entry: &ChainResidue<'_>) -> Option<ResidueBackbone> {
    let hierarchy = entry.hierarchy;
    let residue = hierarchy
        .residue(entry.residue)
        .expect("chain residues come from the model hierarchy");
    if !matches!(
        residue.class(),
        ResidueClass::AminoAcid | ResidueClass::Other
    ) {
        return None;
    }
    let mut backbone = [None; 4];
    let mut any_backbone_atom = false;
    let mut repeated = None;
    for &site_id in &residue.atom_sites {
        let site = hierarchy
            .atom_site(site_id)
            .expect("residue atom sites come from the model hierarchy");
        let atom_name = site
            .metadata
            .label_atom_id
            .as_deref()
            .or(site.metadata.auth_atom_id.as_deref())
            .map(str::trim);
        let Some(backbone_index) = atom_name.and_then(backbone_index) else {
            continue;
        };
        any_backbone_atom = true;
        if backbone[backbone_index].is_some() {
            repeated = Some(DsspBackboneAtom::ALL[backbone_index]);
            break;
        }
        backbone[backbone_index] = Some(Vec3::from(
            model
                .position(site.atom)
                .expect("hierarchy atoms have model positions")
                .value_in(ANGSTROM)
                .expect("positions are lengths"),
        ));
    }
    Some(ResidueBackbone {
        key: entry.residue,
        source: residue_source(entry.chain, residue),
        atoms: repeated.map_or(Ok(backbone), Err),
        any_backbone_atom,
        amino_acid: residue.class() == ResidueClass::AminoAcid,
        is_proline: residue.name.eq_ignore_ascii_case("PRO"),
    })
}

fn extract_chain(
    chain_index: usize,
    backbones: impl Iterator<Item = Option<ResidueBackbone>>,
    residues: &mut Vec<BackboneResidue>,
    report: &mut DsspReport,
    chain_segments: &mut usize,
    degenerate: &mut Option<ResidueId>,
    limits: DsspLimits,
) -> Result<(), DsspError> {
    let mut previous = None;
    let mut forced_break = false;

    for backbone in backbones {
        let Some(backbone) = backbone else {
            report.non_peptide_residues += 1;
            forced_break = true;
            continue;
        };
        let key = backbone.key;
        let atoms = backbone
            .atoms
            .map_err(|atom| DsspError::AmbiguousBackboneAtom { residue: key, atom })?;
        let missing = DsspBackboneAtom::ALL
            .into_iter()
            .enumerate()
            .filter_map(|(index, atom)| atoms[index].is_none().then_some(atom))
            .collect::<Vec<_>>();
        if !missing.is_empty() {
            if !backbone.any_backbone_atom && !backbone.amino_acid {
                report.non_peptide_residues += 1;
            } else {
                report.skipped_residues.push(DsspSkippedResidue {
                    key,
                    source: backbone.source,
                    reason: DsspSkipReason::MissingBackboneAtoms { missing },
                });
                forced_break = true;
            }
            continue;
        }
        if residues.len() >= limits.max_residues {
            return Err(DsspError::LimitExceeded(LimitExceeded::first_beyond(
                Resource::Residues,
                limits.max_residues,
            )));
        }

        let [n, ca, c, o] = atoms.map(|atom| atom.expect("checked complete backbone"));
        if [n, ca, c, o]
            .into_iter()
            .any(|point| !point.is_dssp_representable())
        {
            return Err(DsspError::CoordinateOutOfRange { residue: key });
        }
        let index = residues.len();
        let continuous = previous.is_some_and(|prior: usize| {
            !forced_break && dssp_distance(residues[prior].c, n) <= MAX_PEPTIDE_BOND_ANGSTROMS
        });
        let (segment, chain_break, prev) = if continuous {
            let prior = previous.expect("continuous chain has previous residue");
            (residues[prior].segment, DsspChainBreak::None, Some(prior))
        } else {
            *chain_segments += 1;
            let chain_break = if previous.is_some() {
                report.detected_gaps += 1;
                DsspChainBreak::Gap
            } else {
                DsspChainBreak::NewChain
            };
            (*chain_segments, chain_break, None)
        };
        if let Some(prior) = prev {
            residues[prior].next = Some(index);
        }
        // DSSP assigns each non-initial amide hydrogen from the immediately
        // preceding complete residue in its internal table, even across a
        // numbered gap or chain boundary. Chain-continuity checks still gate
        // torsions and local secondary-structure patterns.
        let h = match residues.last() {
            Some(preceding) if !backbone.is_proline => {
                let h = dssp_reconstructed_hydrogen(n, preceding.c, preceding.o);
                match h {
                    Some(_) => report.reconstructed_amide_hydrogens += 1,
                    None => {
                        degenerate.get_or_insert(key);
                    }
                }
                h
            }
            _ => None,
        };
        residues.push(BackboneResidue {
            key,
            source: backbone.source,
            chain: chain_index,
            n,
            ca,
            c,
            o,
            h,
            is_proline: backbone.is_proline,
            segment,
            prev,
            next: None,
            chain_break,
            phi: None,
            psi: None,
            omega: None,
            alpha: None,
            kappa: None,
            tco: None,
            acceptors: [BondSlot::default(); 2],
            donors: [BondSlot::default(); 2],
            beta_partners: [None; 2],
            sheet: None,
            strand: None,
            helix_positions: [DsspHelixPosition::None; 4],
            secondary_structure: DsspSecondaryStructure::Loop,
        });
        previous = Some(index);
        forced_break = false;
    }
    Ok(())
}

const BACKBONE_NAMES: [&str; 4] = ["N", "CA", "C", "O"];

fn backbone_index(name: &str) -> Option<usize> {
    BACKBONE_NAMES
        .iter()
        .position(|candidate| *candidate == name)
}

fn residue_source(chain: &Chain, residue: &Residue) -> DsspResidueSource {
    DsspResidueSource {
        residue_name: residue.name.clone(),
        chain_label_id: chain.label_id.clone(),
        chain_author_id: chain.author_id.clone(),
        label_sequence_id: residue.label_seq_id,
        author_sequence_id: residue.author_seq_id.clone(),
        insertion_code: residue.insertion_code.clone(),
    }
}

/// The backbone geometry of one residue.
#[derive(Debug, Clone, Copy)]
struct ResidueGeometry {
    phi: Option<f64>,
    psi: Option<f64>,
    omega: Option<f64>,
    alpha: Option<f64>,
    kappa: Option<f64>,
    tco: Option<f64>,
}

/// One range of residues' part of [`analyze_backbone`].
#[derive(Debug)]
struct RangeAnalysis {
    geometry: Vec<ResidueGeometry>,
    /// The candidate pairs whose first residue is in the range.
    pairs: Vec<(usize, usize)>,
    /// For each pair, the energies of the bonds from the first residue's
    /// N-H to the second's C=O and back.
    energies: Vec<[Option<f64>; 2]>,
}

/// Computes every residue's backbone geometry, the candidate pairs and
/// their hydrogen-bond energies, installs the geometry and the bonds, and
/// returns the candidate pairs in order.
fn analyze_backbone(
    residues: &mut [BackboneResidue],
    max_candidate_pairs: usize,
    executor: &dyn Executor,
) -> Result<Vec<(usize, usize)>, DsspError> {
    let grid = CellGrid::new(residues);
    let found = AtomicUsize::new(0);
    let shared = &*residues;
    let ranges = map_ranges(executor, shared.len(), MIN_CHUNK, |range| {
        analyze_range(shared, &grid, range, &found, max_candidate_pairs)
    });
    // DSSP finds every candidate before it evaluates any energy, so the
    // candidate limit takes precedence over an earlier range's energy error.
    if let Some(Err(limit)) = ranges
        .iter()
        .find(|range| matches!(range, Err(DsspError::LimitExceeded(_))))
    {
        return Err(limit.clone());
    }
    // Bonds are inserted in candidate order, which decides DSSP's ties.
    let mut candidates = Vec::with_capacity(found.into_inner());
    let mut index = 0;
    for range in ranges {
        let range = range?;
        for geometry in range.geometry {
            let residue = &mut residues[index];
            residue.phi = geometry.phi;
            residue.psi = geometry.psi;
            residue.omega = geometry.omega;
            residue.alpha = geometry.alpha;
            residue.kappa = geometry.kappa;
            residue.tco = geometry.tco;
            index += 1;
        }
        for (&(first, second), [forward, backward]) in range.pairs.iter().zip(range.energies) {
            if let Some(energy) = forward {
                insert_hydrogen_bond(residues, first, second, energy);
            }
            if let Some(energy) = backward {
                insert_hydrogen_bond(residues, second, first, energy);
            }
        }
        candidates.extend(range.pairs);
    }
    Ok(candidates)
}

fn analyze_range(
    residues: &[BackboneResidue],
    grid: &CellGrid,
    range: Range<usize>,
    found: &AtomicUsize,
    max_candidate_pairs: usize,
) -> Result<RangeAnalysis, DsspError> {
    let geometry = range
        .clone()
        .map(|index| residue_geometry(residues, index))
        .collect();
    let pairs = candidate_pairs(residues, grid, range, found, max_candidate_pairs)?;
    let energies = pairs
        .iter()
        .map(|&(first, second)| {
            let forward = hydrogen_bond_energy(residues, first, second)?;
            let backward = if second != first + 1 {
                hydrogen_bond_energy(residues, second, first)?
            } else {
                None
            };
            Ok([forward, backward])
        })
        .collect::<Result<_, DsspError>>()?;
    Ok(RangeAnalysis {
        geometry,
        pairs,
        energies,
    })
}

fn residue_geometry(residues: &[BackboneResidue], index: usize) -> ResidueGeometry {
    let residue = &residues[index];
    let mut geometry = ResidueGeometry {
        phi: None,
        psi: None,
        omega: None,
        alpha: None,
        kappa: None,
        tco: None,
    };
    if let Some(previous) = residue.prev.map(|previous| &residues[previous]) {
        geometry.phi = dihedral(previous.c, residue.n, residue.ca, residue.c);
        geometry.tco =
            cosine_between(residue.c.sub(residue.o), previous.c.sub(previous.o)).map(f64::from);
    }
    if let Some(next) = residue.next.map(|next| &residues[next]) {
        geometry.psi = dihedral(residue.n, residue.ca, residue.c, next.n);
        geometry.omega = dihedral(residue.ca, residue.c, next.n, next.ca);
    }
    let (Some(previous), Some(next)) = (residue.prev, residue.next) else {
        return geometry;
    };
    let Some(next_next) = residues[next].next else {
        return geometry;
    };
    geometry.alpha = dihedral(
        residues[previous].ca,
        residue.ca,
        residues[next].ca,
        residues[next_next].ca,
    );
    let Some(previous_previous) = residues[previous].prev else {
        return geometry;
    };
    let sequence_is_contiguous = match (
        residues[previous_previous].source.label_sequence_id,
        residues[next_next].source.label_sequence_id,
    ) {
        (Some(first), Some(last)) => first.checked_add(4) == Some(last),
        _ => true,
    };
    if sequence_is_contiguous {
        geometry.kappa = angle_degrees(
            residue.ca.sub(residues[previous_previous].ca),
            residues[next_next].ca.sub(residue.ca),
        );
    }
    geometry
}

fn dihedral(p1: Vec3, p2: Vec3, p3: Vec3, p4: Vec3) -> Option<f64> {
    let v12 = p1.sub(p2);
    let v43 = p4.sub(p3);
    let z = p2.sub(p3);
    let p = z.cross(v12);
    let x = z.cross(v43);
    let y = z.cross(x);
    let x_norm_squared = x.norm_squared();
    let y_norm_squared = y.norm_squared();
    if !x_norm_squared.is_finite()
        || !y_norm_squared.is_finite()
        || x_norm_squared <= 0.0
        || y_norm_squared <= 0.0
    {
        return None;
    }
    let u = p.dot(x) / x_norm_squared.sqrt();
    let v = p.dot(y) / y_norm_squared.sqrt();
    (u.is_finite() && v.is_finite() && (u != 0.0 || v != 0.0))
        .then(|| f64::from(v.atan2(u) * 180.0 / std::f32::consts::PI))
}

fn cosine_between(first: Vec3, second: Vec3) -> Option<f32> {
    let denominator = (first.norm_squared() * second.norm_squared()).sqrt();
    (denominator.is_finite() && denominator > 0.0)
        .then(|| (first.dot(second) / denominator).clamp(-1.0, 1.0))
}

fn angle_degrees(first: Vec3, second: Vec3) -> Option<f64> {
    cosine_between(first, second).map(|cosine| {
        let sine = (1.0 - cosine * cosine).sqrt();
        f64::from(sine.atan2(cosine) * (180.0 / std::f64::consts::PI) as f32)
    })
}

/// The residues binned by the [`spatial_cell`] of their `CA`.
#[derive(Debug)]
enum CellGrid {
    /// Every cell of the residues' box, padded by one cell on each side, in
    /// x, y, z order: residue `members[starts[slot]..starts[slot + 1]]`
    /// lie in cell `slot`.
    Dense {
        origin: (i64, i64, i64),
        dims: (usize, usize, usize),
        starts: Vec<usize>,
        members: Vec<usize>,
    },
    /// Occupied cells only, for a box with many more cells than residues:
    /// `cells[k]` holds residue `members[k]`, sorted by cell.
    Sparse {
        cells: Vec<(i64, i64, i64)>,
        members: Vec<usize>,
    },
}

impl CellGrid {
    /// A box is stored densely when it has at most this many cells per
    /// residue, plus [`Self::DENSE_SPARE_CELLS`].
    const DENSE_CELLS_PER_RESIDUE: usize = 8;
    const DENSE_SPARE_CELLS: usize = 4096;

    fn new(residues: &[BackboneResidue]) -> Self {
        let cells = residues
            .iter()
            .map(|residue| spatial_cell(residue.ca))
            .collect::<Vec<_>>();
        let low = cells
            .iter()
            .fold((i64::MAX, i64::MAX, i64::MAX), |low, cell| {
                (low.0.min(cell.0), low.1.min(cell.1), low.2.min(cell.2))
            });
        let high = cells
            .iter()
            .fold((i64::MIN, i64::MIN, i64::MIN), |high, cell| {
                (high.0.max(cell.0), high.1.max(cell.1), high.2.max(cell.2))
            });
        // Representable cells lie well inside i64, so the padded spans fit
        // in i128; a span that does not fit in usize is never dense.
        let span = |low: i64, high: i64| usize::try_from(i128::from(high) - i128::from(low) + 3);
        let budget = residues
            .len()
            .saturating_mul(Self::DENSE_CELLS_PER_RESIDUE)
            .saturating_add(Self::DENSE_SPARE_CELLS);
        let dense = match (
            span(low.0, high.0),
            span(low.1, high.1),
            span(low.2, high.2),
        ) {
            (Ok(x), Ok(y), Ok(z)) => x
                .checked_mul(y)
                .and_then(|area| area.checked_mul(z))
                .filter(|&count| count <= budget)
                .map(|_| (x, y, z)),
            _ => None,
        };
        let Some(dims) = dense else {
            let mut sorted = cells.into_iter().zip(0..).collect::<Vec<_>>();
            sorted.sort_unstable();
            let (cells, members) = sorted.into_iter().unzip();
            return Self::Sparse { cells, members };
        };
        let origin = (low.0 - 1, low.1 - 1, low.2 - 1);
        let slot = |(x, y, z): (i64, i64, i64)| {
            (((x - origin.0) as usize * dims.1) + (y - origin.1) as usize) * dims.2
                + (z - origin.2) as usize
        };
        let mut starts = vec![0; dims.0 * dims.1 * dims.2 + 1];
        for &cell in &cells {
            starts[slot(cell) + 1] += 1;
        }
        for index in 1..starts.len() {
            starts[index] += starts[index - 1];
        }
        let mut next = starts.clone();
        let mut members = vec![0; cells.len()];
        for (index, &cell) in cells.iter().enumerate() {
            members[next[slot(cell)]] = index;
            next[slot(cell)] += 1;
        }
        Self::Dense {
            origin,
            dims,
            starts,
            members,
        }
    }

    /// The residues in cells `(x, y, z - 1)` through `(x, y, z + 1)`, where
    /// `(x, y, z)` is within one cell of a residue's cell.
    fn column(&self, x: i64, y: i64, z: i64) -> &[usize] {
        match self {
            Self::Dense {
                origin,
                dims,
                starts,
                members,
            } => {
                let base = (((x - origin.0) as usize * dims.1) + (y - origin.1) as usize) * dims.2;
                let first = base + (z - 1 - origin.2) as usize;
                &members[starts[first]..starts[first + 3]]
            }
            Self::Sparse { cells, members } => {
                let start = cells.partition_point(|&cell| cell < (x, y, z - 1));
                let end = cells.partition_point(|&cell| cell <= (x, y, z + 1));
                &members[start..end]
            }
        }
    }
}

/// The candidate pairs whose first residue is in `range`: residues with
/// `CA` atoms within 9 Å, in order. `found` counts pairs across ranges.
fn candidate_pairs(
    residues: &[BackboneResidue],
    grid: &CellGrid,
    range: Range<usize>,
    found: &AtomicUsize,
    max_candidate_pairs: usize,
) -> Result<Vec<(usize, usize)>, DsspError> {
    let maximum_squared = MAX_CA_PAIR_DISTANCE_ANGSTROMS.powi(2);
    let mut pairs = Vec::new();
    let mut partners = Vec::new();
    for first in range {
        let (x, y, z) = spatial_cell(residues[first].ca);
        partners.clear();
        for dx in -1..=1 {
            for dy in -1..=1 {
                partners.extend(grid.column(x + dx, y + dy, z).iter().filter_map(|&second| {
                    (second > first
                        && dssp_distance_squared(residues[first].ca, residues[second].ca)
                            <= maximum_squared)
                        .then_some(second)
                }));
            }
        }
        if found.fetch_add(partners.len(), Ordering::Relaxed) + partners.len() > max_candidate_pairs
        {
            return Err(DsspError::LimitExceeded(LimitExceeded::first_beyond(
                Resource::CandidatePairs,
                max_candidate_pairs,
            )));
        }
        partners.sort_unstable();
        pairs.extend(partners.iter().map(|&second| (first, second)));
    }
    Ok(pairs)
}

fn spatial_cell(point: Vec3) -> (i64, i64, i64) {
    (
        (point.x / MAX_CA_PAIR_DISTANCE_ANGSTROMS as f32).floor() as i64,
        (point.y / MAX_CA_PAIR_DISTANCE_ANGSTROMS as f32).floor() as i64,
        (point.z / MAX_CA_PAIR_DISTANCE_ANGSTROMS as f32).floor() as i64,
    )
}

/// The DSSP energy of the bond from `donor`'s N-H to `acceptor`'s C=O;
/// `None` when the donor is a proline or has no amide hydrogen.
fn hydrogen_bond_energy(
    residues: &[BackboneResidue],
    donor: usize,
    acceptor: usize,
) -> Result<Option<f64>, DsspError> {
    if residues[donor].is_proline {
        return Ok(None);
    }
    let Some(hydrogen) = residues[donor].h else {
        return Ok(None);
    };
    let distances = [
        dssp_distance(hydrogen, residues[acceptor].o),
        dssp_distance(hydrogen, residues[acceptor].c),
        dssp_distance(residues[donor].n, residues[acceptor].c),
        dssp_distance(residues[donor].n, residues[acceptor].o),
    ];
    if distances.iter().any(|distance| !distance.is_finite()) {
        return Err(DsspError::NonFiniteGeometry {
            residue: residues[donor].key,
        });
    }
    let mut energy = if distances
        .iter()
        .any(|distance| *distance < MIN_ELECTROSTATIC_DISTANCE_ANGSTROMS)
    {
        MIN_HBOND_ENERGY_KCAL_PER_MOL
    } else {
        DSSP_COUPLING_KCAL_ANGSTROM_PER_MOL / distances[0]
            - DSSP_COUPLING_KCAL_ANGSTROM_PER_MOL / distances[1]
            + DSSP_COUPLING_KCAL_ANGSTROM_PER_MOL / distances[2]
            - DSSP_COUPLING_KCAL_ANGSTROM_PER_MOL / distances[3]
    };
    energy = (energy * 1000.0).round() / 1000.0;
    Ok(Some(energy.max(MIN_HBOND_ENERGY_KCAL_PER_MOL)))
}

fn insert_hydrogen_bond(
    residues: &mut [BackboneResidue],
    donor: usize,
    acceptor: usize,
    energy: f64,
) {
    insert_bond(&mut residues[donor].acceptors, acceptor, energy);
    insert_bond(&mut residues[acceptor].donors, donor, energy);
}

fn insert_bond(slots: &mut [BondSlot; 2], partner: usize, energy: f64) {
    if energy < slots[0].energy {
        slots[1] = slots[0];
        slots[0] = BondSlot {
            partner: Some(partner),
            energy,
        };
    } else if energy < slots[1].energy {
        slots[1] = BondSlot {
            partner: Some(partner),
            energy,
        };
    }
}

fn test_bond(residues: &[BackboneResidue], donor: usize, acceptor: usize) -> bool {
    residues[donor]
        .acceptors
        .iter()
        .any(|bond| bond.partner == Some(acceptor) && bond.energy < MAX_HBOND_ENERGY_KCAL_PER_MOL)
}

fn test_bridge(residues: &[BackboneResidue], first: usize, second: usize) -> Option<BridgeType> {
    let a = residues[first].prev?;
    let c = residues[first].next?;
    let d = residues[second].prev?;
    let f = residues[second].next?;
    if residues[a].segment != residues[c].segment || residues[d].segment != residues[f].segment {
        return None;
    }
    if (test_bond(residues, c, second) && test_bond(residues, second, a))
        || (test_bond(residues, f, first) && test_bond(residues, first, d))
    {
        Some(BridgeType::Parallel)
    } else if (test_bond(residues, c, d) && test_bond(residues, f, a))
        || (test_bond(residues, second, first) && test_bond(residues, first, second))
    {
        Some(BridgeType::Antiparallel)
    } else {
        None
    }
}

fn calculate_beta_sheets(
    residues: &mut [BackboneResidue],
    candidates: &[(usize, usize)],
    max_ladders: usize,
    executor: &dyn Executor,
) -> Result<BetaCounts, DsspError> {
    let shared = &*residues;
    let found = map_ranges(executor, candidates.len(), MIN_CHUNK, |range| {
        candidates[range]
            .iter()
            .filter_map(|&(first, second)| {
                test_bridge(shared, first, second).map(|bridge_type| (first, second, bridge_type))
            })
            .collect::<Vec<_>>()
    });
    let mut bridges = Vec::<Bridge>::new();
    // A bridge extends the ladder whose last pair is its predecessor. Every
    // candidate pair joins one ladder, so at most one ladder of each type
    // ends in a given pair; this index finds it.
    let mut ends = BTreeMap::<(bool, usize, usize), usize>::new();
    let mut counts = BetaCounts::default();
    for (first, second, bridge_type) in found.into_iter().flatten() {
        counts.bridges += 1;
        if residues[first].chain == residues[second].chain {
            counts.intra_chain += 1;
        } else {
            counts.inter_chain += 1;
        }
        let parallel = bridge_type == BridgeType::Parallel;
        let predecessor = first.checked_sub(1).zip(if parallel {
            second.checked_sub(1)
        } else {
            second.checked_add(1)
        });
        let extended = predecessor
            .and_then(|(i, j)| ends.remove(&(parallel, i, j)))
            .inspect(|&ladder| {
                let bridge = &mut bridges[ladder];
                bridge.i.push(first);
                if parallel {
                    bridge.j.push_back(second);
                } else {
                    bridge.j.push_front(second);
                }
            });
        let ladder = match extended {
            Some(ladder) => ladder,
            None => {
                if bridges.len() >= max_ladders {
                    return Err(DsspError::LimitExceeded(LimitExceeded::first_beyond(
                        Resource::Ladders,
                        max_ladders,
                    )));
                }
                bridges.push(Bridge {
                    bridge_type,
                    i: vec![first],
                    j: VecDeque::from([second]),
                    sheet: 0,
                    ladder: 0,
                });
                bridges.len() - 1
            }
        };
        ends.insert((parallel, first, second), ladder);
    }

    bridges.sort_by(|left, right| {
        let left_residue = &residues[left.i[0]];
        let right_residue = &residues[right.i[0]];
        (
            &left_residue.source.chain_label_id,
            left_residue.key.raw(),
            left.i[0],
        )
            .cmp(&(
                &right_residue.source.chain_label_id,
                right_residue.key.raw(),
                right.i[0],
            ))
    });
    merge_bulges(residues, &mut bridges);
    if bridges.len() > max_ladders {
        return Err(DsspError::LimitExceeded(LimitExceeded::new(
            Resource::Ladders,
            max_ladders,
            bridges.len(),
        )));
    }
    counts.ladders = bridges.len();
    assign_sheets(&mut bridges, residues.len());
    assign_beta_results(residues, &bridges);
    Ok(counts)
}

/// Merges ladders across beta bulges, as DSSP does: each ladder in order
/// absorbs every later ladder, in order, that continues it across a bulge,
/// growing as it goes.
fn merge_bulges(residues: &[BackboneResidue], bridges: &mut Vec<Bridge>) {
    // A ladder continues another only if its first strand starts at most
    // five residues after the other's first strand ends, so later ladders
    // are found by where their first strand starts, which merging never
    // changes.
    let mut starts = BTreeMap::<usize, Vec<usize>>::new();
    for (index, bridge) in bridges.iter().enumerate() {
        starts.entry(bridge.i[0]).or_default().push(index);
    }
    let mut merged = vec![false; bridges.len()];
    for first in 0..bridges.len() {
        if merged[first] {
            continue;
        }
        let mut examined = first;
        loop {
            let end = *bridges[first].i.last().expect("bridge has first strand");
            let next = starts
                .range(end..=end.saturating_add(5))
                .filter_map(|(_, ladders)| {
                    let later = ladders.partition_point(|&ladder| ladder <= examined);
                    ladders[later..]
                        .iter()
                        .copied()
                        .find(|&ladder| !merged[ladder])
                })
                .min();
            let Some(second) = next else {
                break;
            };
            examined = second;
            if !bulge_joins(residues, &bridges[first], &bridges[second]) {
                continue;
            }
            merged[second] = true;
            let other_i = std::mem::take(&mut bridges[second].i);
            let other_j = std::mem::take(&mut bridges[second].j);
            bridges[first].i.extend(other_i);
            match bridges[first].bridge_type {
                BridgeType::Parallel => bridges[first].j.extend(other_j),
                BridgeType::Antiparallel => {
                    for value in other_j.into_iter().rev() {
                        bridges[first].j.push_front(value);
                    }
                }
            }
        }
    }
    let mut merged = merged.into_iter();
    bridges.retain(|_| !merged.next().expect("one flag per ladder"));
}

/// Whether ladder `second` continues ladder `first` across a beta bulge.
fn bulge_joins(residues: &[BackboneResidue], first: &Bridge, second: &Bridge) -> bool {
    let ibi = first.i[0];
    let iei = *first.i.last().expect("bridge has first strand");
    let jbi = first.j[0];
    let jei = *first.j.back().expect("bridge has second strand");
    let ibj = second.i[0];
    let iej = *second.i.last().expect("bridge has first strand");
    let jbj = second.j[0];
    let jej = *second.j.back().expect("bridge has second strand");

    let incompatible = first.bridge_type != second.bridge_type
        || residues[ibi].segment != residues[iej].segment
        || residues[jbi.min(jbj)].segment != residues[jei.max(jej)].segment
        || ibj.checked_sub(iei).is_none_or(|gap| gap >= 6)
        || (iei >= ibj && ibi <= iej);
    if incompatible {
        return false;
    }
    let first_gap = ibj - iei;
    match first.bridge_type {
        BridgeType::Parallel => jbj
            .checked_sub(jei)
            .is_some_and(|second_gap| (second_gap < 6 && first_gap < 3) || second_gap < 3),
        BridgeType::Antiparallel => jbi
            .checked_sub(jej)
            .is_some_and(|second_gap| (second_gap < 6 && first_gap < 3) || second_gap < 3),
    }
}

/// Groups ladders that share a residue into sheets. Sheets are numbered in
/// the order of their first ladder, and ladders in sheet order, each sheet's
/// ladders in order.
fn assign_sheets(bridges: &mut [Bridge], residue_count: usize) {
    fn root(parents: &mut [usize], mut ladder: usize) -> usize {
        while parents[ladder] != ladder {
            parents[ladder] = parents[parents[ladder]];
            ladder = parents[ladder];
        }
        ladder
    }
    let mut parents = (0..bridges.len()).collect::<Vec<_>>();
    let mut owner = vec![None; residue_count];
    for (ladder, bridge) in bridges.iter().enumerate() {
        for &residue in bridge.i.iter().chain(&bridge.j) {
            match owner[residue] {
                None => owner[residue] = Some(ladder),
                Some(other) => {
                    let (left, right) = (root(&mut parents, ladder), root(&mut parents, other));
                    parents[left.max(right)] = left.min(right);
                }
            }
        }
    }
    let mut sheets = Vec::<Vec<usize>>::new();
    let mut sheet_of = vec![usize::MAX; bridges.len()];
    for ladder in 0..bridges.len() {
        let first = root(&mut parents, ladder);
        if sheet_of[first] == usize::MAX {
            sheet_of[first] = sheets.len();
            sheets.push(Vec::new());
        }
        sheets[sheet_of[first]].push(ladder);
    }
    let mut ladder_number = 0;
    for (sheet, ladders) in sheets.into_iter().enumerate() {
        for ladder in ladders {
            bridges[ladder].sheet = sheet + 1;
            bridges[ladder].ladder = ladder_number;
            ladder_number += 1;
        }
    }
}

fn assign_beta_results(residues: &mut [BackboneResidue], bridges: &[Bridge]) {
    for bridge in bridges {
        let first_slot = bridge
            .i
            .iter()
            .any(|index| residues[*index].beta_partners[0].is_some())
            as usize;
        let second_slot = bridge
            .j
            .iter()
            .any(|index| residues[*index].beta_partners[0].is_some())
            as usize;
        let parallel = bridge.bridge_type == BridgeType::Parallel;
        let first_partners = if parallel {
            bridge.j.iter().copied().collect::<Vec<_>>()
        } else {
            bridge.j.iter().rev().copied().collect::<Vec<_>>()
        };
        let second_partners = if parallel {
            bridge.i.clone()
        } else {
            bridge.i.iter().rev().copied().collect()
        };
        for (&residue, partner) in bridge.i.iter().zip(first_partners) {
            residues[residue].beta_partners[first_slot] = Some(BetaSlot {
                partner,
                ladder: bridge.ladder,
                parallel,
            });
        }
        for (&residue, partner) in bridge.j.iter().zip(second_partners) {
            residues[residue].beta_partners[second_slot] = Some(BetaSlot {
                partner,
                ladder: bridge.ladder,
                parallel,
            });
        }

        let structure = if bridge.i.len() > 1 {
            DsspSecondaryStructure::ExtendedStrand
        } else {
            DsspSecondaryStructure::BetaBridge
        };
        let first_start = bridge.i[0];
        let first_end = *bridge.i.last().expect("bridge strand");
        for residue in &mut residues[first_start..=first_end] {
            if residue.secondary_structure != DsspSecondaryStructure::ExtendedStrand {
                residue.secondary_structure = structure;
            }
            residue.sheet = Some(bridge.sheet);
        }
        let second_start = bridge.j[0];
        let second_end = *bridge.j.back().expect("bridge strand");
        for residue in &mut residues[second_start..=second_end] {
            if residue.secondary_structure != DsspSecondaryStructure::ExtendedStrand {
                residue.secondary_structure = structure;
            }
            residue.sheet = Some(bridge.sheet);
        }
    }

    // Strands are numbered sheet by sheet, each sheet's in residue order.
    let mut sheets = vec![Vec::new(); bridges.len() + 1];
    for (index, residue) in residues.iter().enumerate() {
        if let Some(sheet) = residue.sheet {
            sheets[sheet].push(index);
        }
    }
    let mut strand = 0;
    for members in sheets {
        let mut previous: Option<usize> = None;
        for index in members {
            if previous.is_none_or(|prior| {
                prior + 1 != index || residues[prior].segment != residues[index].segment
            }) {
                strand += 1;
            }
            residues[index].strand = Some(strand);
            previous = Some(index);
        }
    }
}

fn calculate_helices_turns_and_bends(residues: &mut [BackboneResidue]) {
    for helix_index in 0..3 {
        let stride = helix_index + 3;
        for start in 0..residues.len().saturating_sub(stride) {
            let end = start + stride;
            if residues[start].segment != residues[end].segment || !test_bond(residues, end, start)
            {
                continue;
            }
            residues[end].helix_positions[helix_index] = DsspHelixPosition::End;
            for residue in residues.iter_mut().take(end).skip(start + 1) {
                if residue.helix_positions[helix_index] == DsspHelixPosition::None {
                    residue.helix_positions[helix_index] = DsspHelixPosition::Middle;
                }
            }
            residues[start].helix_positions[helix_index] =
                if residues[start].helix_positions[helix_index] == DsspHelixPosition::End {
                    DsspHelixPosition::StartAndEnd
                } else {
                    DsspHelixPosition::Start
                };
        }
    }
    for start in 1..residues.len().saturating_sub(4) {
        if is_helix_start(&residues[start], 1) && is_helix_start(&residues[start - 1], 1) {
            for residue in residues.iter_mut().skip(start).take(4) {
                residue.secondary_structure = DsspSecondaryStructure::AlphaHelix;
            }
        }
    }
    for start in 1..residues.len().saturating_sub(3) {
        if is_helix_start(&residues[start], 0) && is_helix_start(&residues[start - 1], 0) {
            let empty = residues[start..=start + 2].iter().all(|residue| {
                matches!(
                    residue.secondary_structure,
                    DsspSecondaryStructure::Loop | DsspSecondaryStructure::Helix3_10
                )
            });
            if empty {
                for residue in residues.iter_mut().skip(start).take(3) {
                    residue.secondary_structure = DsspSecondaryStructure::Helix3_10;
                }
            }
        }
    }
    for start in 1..residues.len().saturating_sub(5) {
        if is_helix_start(&residues[start], 2) && is_helix_start(&residues[start - 1], 2) {
            let empty = residues[start..=start + 4].iter().all(|residue| {
                matches!(
                    residue.secondary_structure,
                    DsspSecondaryStructure::Loop
                        | DsspSecondaryStructure::PiHelix
                        | DsspSecondaryStructure::AlphaHelix
                )
            });
            if empty {
                for residue in residues.iter_mut().skip(start).take(5) {
                    residue.secondary_structure = DsspSecondaryStructure::PiHelix;
                }
            }
        }
    }
    for index in 1..residues.len().saturating_sub(1) {
        if residues[index].secondary_structure != DsspSecondaryStructure::Loop {
            continue;
        }
        let mut turn = false;
        for helix_index in 0..3 {
            let stride = helix_index + 3;
            for offset in 1..stride {
                if index >= offset && is_helix_start(&residues[index - offset], helix_index) {
                    turn = true;
                    break;
                }
            }
            if turn {
                break;
            }
        }
        if turn {
            residues[index].secondary_structure = DsspSecondaryStructure::Turn;
        } else if residues[index].kappa.is_some_and(|kappa| kappa > 70.0) {
            residues[index].secondary_structure = DsspSecondaryStructure::Bend;
        }
    }
}

fn is_helix_start(residue: &BackboneResidue, helix_index: usize) -> bool {
    matches!(
        residue.helix_positions[helix_index],
        DsspHelixPosition::Start | DsspHelixPosition::StartAndEnd
    )
}

fn calculate_polyproline_helices(residues: &mut [BackboneResidue], stretch: usize) {
    const EPSILON: f64 = 29.0;
    const PHI_MIN: f64 = -75.0 - EPSILON;
    const PHI_MAX: f64 = -75.0 + EPSILON;
    const PSI_MIN: f64 = 145.0 - EPSILON;
    const PSI_MAX: f64 = 145.0 + EPSILON;

    for start in 1..residues.len().saturating_sub(stretch) {
        let qualifies = residues[start..start + stretch].iter().all(|residue| {
            residue
                .phi
                .is_some_and(|phi| (PHI_MIN..=PHI_MAX).contains(&phi))
                && residue
                    .psi
                    .is_some_and(|psi| (PSI_MIN..=PSI_MAX).contains(&psi))
        });
        if !qualifies {
            continue;
        }
        let prior = residues[start].helix_positions[3];
        residues[start].helix_positions[3] = match (stretch, prior) {
            (2, DsspHelixPosition::End) => DsspHelixPosition::Middle,
            (3, DsspHelixPosition::End) => DsspHelixPosition::StartAndEnd,
            (_, DsspHelixPosition::None) => DsspHelixPosition::Start,
            _ => prior,
        };
        for residue in residues.iter_mut().skip(start + 1).take(stretch - 2) {
            residue.helix_positions[3] = DsspHelixPosition::Middle;
        }
        residues[start + stretch - 1].helix_positions[3] = DsspHelixPosition::End;
        for residue in residues.iter_mut().skip(start).take(stretch) {
            if residue.secondary_structure == DsspSecondaryStructure::Loop {
                residue.secondary_structure = DsspSecondaryStructure::PolyprolineIIHelix;
            }
        }
    }
}

fn public_residue(residue: BackboneResidue, keys: &[ResidueId]) -> DsspResidue {
    let convert_bond = |bond: BondSlot| {
        bond.partner.map(|partner| DsspHydrogenBond {
            partner: keys[partner],
            energy_kcal_per_mol: bond.energy,
        })
    };
    let convert_beta = |partner: Option<BetaSlot>| {
        partner.map(|partner| DsspBetaPartner {
            partner: keys[partner.partner],
            ladder: partner.ladder,
            direction: if partner.parallel {
                DsspBridgeDirection::Parallel
            } else {
                DsspBridgeDirection::Antiparallel
            },
        })
    };
    DsspResidue {
        key: residue.key,
        source: residue.source,
        secondary_structure: residue.secondary_structure,
        chain_break: residue.chain_break,
        phi: residue.phi,
        psi: residue.psi,
        omega: residue.omega,
        alpha: residue.alpha,
        kappa: residue.kappa,
        tco: residue.tco,
        acceptors: residue.acceptors.map(convert_bond),
        donors: residue.donors.map(convert_bond),
        beta_partners: residue.beta_partners.map(convert_beta),
        sheet: residue.sheet,
        strand: residue.strand,
        helix_positions: residue.helix_positions,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn serial() -> ScopedThreads {
        ScopedThreads::new(NonZeroUsize::MIN)
    }

    fn backbone(count: usize) -> Vec<BackboneResidue> {
        let mut residues = (0..count)
            .zip(0..=u32::MAX)
            .map(|(index, raw)| BackboneResidue {
                key: crate::topology::ResidueId::new(raw),
                source: DsspResidueSource {
                    residue_name: "ALA".to_owned(),
                    chain_label_id: "A".to_owned(),
                    chain_author_id: Some("A".to_owned()),
                    label_sequence_id: Some(index as i32 + 1),
                    author_sequence_id: Some((index + 1).to_string()),
                    insertion_code: None,
                },
                chain: 0,
                n: Vec3 {
                    x: index as f32,
                    y: 0.0,
                    z: 0.0,
                },
                ca: Vec3 {
                    x: index as f32,
                    y: 1.0,
                    z: 0.0,
                },
                c: Vec3 {
                    x: index as f32,
                    y: 1.0,
                    z: 1.0,
                },
                o: Vec3 {
                    x: index as f32,
                    y: 1.0,
                    z: 2.0,
                },
                h: None,
                is_proline: false,
                segment: 1,
                prev: index.checked_sub(1),
                next: (index + 1 < count).then_some(index + 1),
                chain_break: if index == 0 {
                    DsspChainBreak::NewChain
                } else {
                    DsspChainBreak::None
                },
                phi: None,
                psi: None,
                omega: None,
                alpha: None,
                kappa: None,
                tco: None,
                acceptors: [BondSlot::default(); 2],
                donors: [BondSlot::default(); 2],
                beta_partners: [None; 2],
                sheet: None,
                strand: None,
                helix_positions: [DsspHelixPosition::None; 4],
                secondary_structure: DsspSecondaryStructure::Loop,
            })
            .collect::<Vec<_>>();
        if let Some(last) = residues.last_mut() {
            last.next = None;
        }
        residues
    }

    fn assert_helix(stride: usize, structure: DsspSecondaryStructure) {
        let mut residues = backbone(stride + 3);
        insert_bond(&mut residues[stride].acceptors, 0, -1.0);
        insert_bond(&mut residues[stride + 1].acceptors, 1, -1.0);
        calculate_helices_turns_and_bends(&mut residues);
        assert!(residues[1..=stride]
            .iter()
            .all(|residue| residue.secondary_structure == structure));
    }

    #[test]
    fn dssp_angles_with_overflowing_intermediates_are_undefined() {
        let points = [
            Point3::new(1.0e10, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 1.0e10, 0.0),
            Point3::new(0.0, 1.0e10, 1.0e10),
        ]
        .map(Vec3::from);
        assert_eq!(dihedral(points[0], points[1], points[2], points[3]), None);
        assert_eq!(cosine_between(points[0], points[3]), None);
        assert_eq!(angle_degrees(points[0], points[3]), None);
    }

    #[test]
    fn dssp_torsion_respects_single_precision_at_output_rounding_boundary() {
        // Supplied PDB 1A2L, chain B, phi of residue 67. DSSP 4.6.1
        // reports -52.0 degrees; evaluating the decimal coordinates in f64
        // instead crosses the one-decimal output boundary to -51.9.
        let points = [
            Point3::new(69.913, 96.863, 62.518),
            Point3::new(70.576, 95.880, 63.129),
            Point3::new(70.017, 94.541, 63.318),
            Point3::new(68.636, 94.556, 63.965),
        ]
        .map(Vec3::from);
        let angle = dihedral(points[0], points[1], points[2], points[3]).unwrap();
        assert_eq!(format!("{angle:.1}"), "-52.0");
    }

    #[test]
    fn dssp_dihedral_matches_reference_sign_convention() {
        let angle = dihedral(
            Vec3 {
                x: 1.0,
                y: 0.0,
                z: 0.0,
            },
            Vec3 {
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            Vec3 {
                x: 0.0,
                y: 1.0,
                z: 0.0,
            },
            Vec3 {
                x: 0.0,
                y: 1.0,
                z: 1.0,
            },
        )
        .expect("defined dihedral");
        assert!((angle - -90.0).abs() < 1.0e-12);
    }

    #[test]
    fn bond_slots_keep_two_strictly_strongest_candidates() {
        let mut slots = [BondSlot::default(); 2];
        insert_bond(&mut slots, 1, -0.6);
        insert_bond(&mut slots, 2, -1.2);
        insert_bond(&mut slots, 3, -0.8);
        insert_bond(&mut slots, 4, -0.8);
        assert_eq!(slots[0].partner, Some(2));
        assert_eq!(slots[1].partner, Some(3));
    }

    #[test]
    fn consecutive_turns_assign_three_ten_alpha_and_pi_helices() {
        assert_helix(3, DsspSecondaryStructure::Helix3_10);
        assert_helix(4, DsspSecondaryStructure::AlphaHelix);
        assert_helix(5, DsspSecondaryStructure::PiHelix);
    }

    #[test]
    fn dssp4_polyproline_window_and_default_stretch_assign_p() {
        let mut residues = backbone(7);
        for residue in &mut residues[1..=3] {
            residue.phi = Some(-75.0);
            residue.psi = Some(145.0);
        }
        calculate_polyproline_helices(&mut residues, 3);
        assert!(residues[1..=3].iter().all(|residue| {
            residue.secondary_structure == DsspSecondaryStructure::PolyprolineIIHelix
        }));
        assert_eq!(residues[1].helix_positions[3], DsspHelixPosition::Start);
        assert_eq!(residues[2].helix_positions[3], DsspHelixPosition::Middle);
        assert_eq!(residues[3].helix_positions[3], DsspHelixPosition::End);
    }

    #[test]
    fn two_residue_polyproline_mode_includes_the_last_complete_window() {
        let mut residues = backbone(7);
        for residue in &mut residues[4..=5] {
            residue.phi = Some(-75.0);
            residue.psi = Some(145.0);
        }
        calculate_polyproline_helices(&mut residues, 2);
        assert!(residues[4..=5].iter().all(|residue| {
            residue.secondary_structure == DsspSecondaryStructure::PolyprolineIIHelix
        }));
        assert_eq!(residues[4].helix_positions[3], DsspHelixPosition::Start);
        assert_eq!(residues[5].helix_positions[3], DsspHelixPosition::End);
    }

    #[test]
    fn ladder_limit_bounds_provisional_topology_before_bulge_merging() {
        let mut residues = backbone(10);
        insert_bond(&mut residues[3].acceptors, 5, -1.0);
        insert_bond(&mut residues[5].acceptors, 1, -1.0);
        insert_bond(&mut residues[5].acceptors, 7, -1.0);
        insert_bond(&mut residues[7].acceptors, 3, -1.0);

        assert!(matches!(
            calculate_beta_sheets(&mut residues, &[(2, 5), (4, 7)], 1, &serial()),
            Err(DsspError::LimitExceeded(limit))
                if limit == LimitExceeded::first_beyond(Resource::Ladders, 1)
        ));
    }

    #[test]
    fn parallel_and_antiparallel_bridge_formulas_are_distinct() {
        let mut parallel = backbone(8);
        insert_bond(&mut parallel[3].acceptors, 5, -1.0);
        insert_bond(&mut parallel[5].acceptors, 1, -1.0);
        let counts = calculate_beta_sheets(&mut parallel, &[(2, 5)], 1, &serial()).expect("bridge");
        assert_eq!(counts.bridges, 1);
        assert_eq!(
            parallel[2].secondary_structure,
            DsspSecondaryStructure::BetaBridge
        );
        assert!(parallel[2].beta_partners[0]
            .is_some_and(|partner| partner.parallel && partner.partner == 5));

        let mut antiparallel = backbone(8);
        insert_bond(&mut antiparallel[2].acceptors, 5, -1.0);
        insert_bond(&mut antiparallel[5].acceptors, 2, -1.0);
        calculate_beta_sheets(&mut antiparallel, &[(2, 5)], 1, &serial()).expect("bridge");
        assert!(antiparallel[2].beta_partners[0]
            .is_some_and(|partner| !partner.parallel && partner.partner == 5));
        assert_eq!(
            calculate_beta_sheets(&mut backbone(8), &[], 0, &serial())
                .expect("no bridge")
                .ladders,
            0
        );
    }
}
