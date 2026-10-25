use std::collections::{BTreeMap, VecDeque};
use std::ops::Range;
use std::sync::atomic::{AtomicUsize, Ordering};

use crate::geometry::Point3;
use crate::structure::ModelView;
use crate::topology::{Chain, ChainId, Hierarchy, Residue, ResidueClass, ResidueId, Topology};
use crate::units::ANGSTROM;

use super::*;
use crate::parallel::{map_ranges, map_weighted_ranges};

/// Fewest residues a stage gives one thread.
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

/// One residue of DSSP's residue table: its backbone, its place in the
/// table, and the geometry and hydrogen bonds that analysis gives it.
#[derive(Debug, Clone)]
struct BackboneResidue {
    key: ResidueId,
    /// The hierarchy chain, for the published source labels.
    source_chain: ChainId,
    label_sequence_id: Option<i32>,
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
    geometry: ResidueGeometry,
    acceptors: [BondSlot; 2],
    donors: [BondSlot; 2],
}

/// A residue's beta-sheet assignment.
#[derive(Debug, Clone, Copy, Default)]
struct BetaResidue {
    partners: [Option<BetaSlot>; 2],
    sheet: Option<usize>,
    strand: Option<usize>,
}

/// A residue's secondary structure and helix positions.
type Classified = (DsspSecondaryStructure, [DsspHelixPosition; 4]);

#[derive(Debug, Clone, Copy)]
struct ChainResidue<'a> {
    hierarchy: &'a Hierarchy,
    chain: &'a Chain,
    residue: ResidueId,
    discovery_order: usize,
    /// Whether the residue is classified as amino acid or other, the
    /// classes DSSP reads.
    peptide: bool,
}

#[derive(Debug)]
struct LogicalChain<'a> {
    chain: &'a Chain,
    discovery_order: usize,
    /// The discovery order of the chain's first residue.
    first_residue: usize,
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

/// Assigns `model`, whose topology `prepared` was prepared from.
pub(super) fn assign(
    prepared: &Prepared,
    model: ModelView<'_>,
    options: DsspOptions,
    executor: &dyn Executor,
) -> Result<DsspResult, DsspError> {
    let (mut residues, mut report, chain_segments) =
        extract_backbones(prepared, model, options.limits, executor)?;
    if residues.is_empty() {
        return Err(DsspError::NoAnalyzableProteinResidues);
    }

    let candidates = analyze_backbone(&mut residues, options.limits.max_candidate_pairs, executor)?;
    report.candidate_hydrogen_bond_pairs = candidates.len();
    let hierarchy = model.topology().hierarchy();
    let sheets = calculate_beta_sheets(
        hierarchy,
        &residues,
        &candidates,
        options.limits.max_ladders,
        executor,
    )?;
    drop(candidates);
    report.generated_ladders = sheets.counts.ladders;

    // Each range classifies its residues and publishes them.
    let published = map_ranges(executor, residues.len(), MIN_CHUNK, |range| {
        let classified = classify_range(
            &residues,
            &sheets.structure,
            options.min_polyproline_stretch,
            range.clone(),
        );
        let mut counts = [0; DsspSecondaryStructure::ALL.len()];
        let mut hydrogen_bonds = 0;
        let published = range
            .zip(classified)
            .map(|(index, classified)| {
                hydrogen_bonds += residues[index]
                    .donors
                    .iter()
                    .filter(|bond| {
                        bond.partner.is_some() && bond.energy < MAX_HBOND_ENERGY_KCAL_PER_MOL
                    })
                    .count();
                // `ALL` lists the assignments in declaration order.
                counts[classified.0 as usize] += 1;
                public_residue(
                    hierarchy,
                    &residues,
                    index,
                    sheets.residues[index],
                    classified,
                )
            })
            .collect::<Vec<_>>();
        (published, counts, hydrogen_bonds)
    });
    let mut public_residues = Vec::with_capacity(residues.len());
    let mut counts = [0; DsspSecondaryStructure::ALL.len()];
    let mut hydrogen_bonds = 0;
    for (published, range_counts, range_bonds) in published {
        public_residues.extend(published);
        for (count, range_count) in counts.iter_mut().zip(range_counts) {
            *count += range_count;
        }
        hydrogen_bonds += range_bonds;
    }
    debug_assert_eq!(prepared.lookup.len(), public_residues.len());
    let lookup = prepared.lookup.clone();
    let statistics = DsspStatistics {
        analyzed_residues: public_residues.len(),
        chain_segments,
        hydrogen_bonds,
        beta_bridges: sheets.counts.bridges,
        intra_chain_beta_bridges: sheets.counts.intra_chain,
        inter_chain_beta_bridges: sheets.counts.inter_chain,
        secondary_structure_counts: DsspSecondaryStructure::ALL
            .into_iter()
            .zip(counts)
            .collect(),
    };
    Ok(DsspResult {
        residues: public_residues,
        lookup,
        statistics,
        report,
    })
}

/// What assignment reads from a topology, prepared once for every structure
/// of its layout.
#[derive(Debug, Clone)]
pub(super) struct Prepared {
    ignored_instances: Vec<MoleculeInstanceId>,
    /// Each logical chain's residues in DSSP order, as positions in
    /// `residues`, or `None` for one classified as neither amino acid nor
    /// other.
    chains: Vec<Vec<Option<usize>>>,
    /// The peptide residues, in DSSP order.
    residues: Vec<PreparedResidue>,
    /// The key and table position of each residue with a complete backbone,
    /// by key: an assignment that succeeds places exactly these.
    lookup: Vec<(ResidueId, usize)>,
}

/// A peptide residue's backbone atoms and labels.
#[derive(Debug, Clone)]
struct PreparedResidue {
    key: ResidueId,
    source_chain: ChainId,
    label_sequence_id: Option<i32>,
    /// The dense indices of `N`, `CA`, `C`, `O`, or the first repeated
    /// atom.
    atoms: Result<[Option<usize>; 4], DsspBackboneAtom>,
    any_backbone_atom: bool,
    amino_acid: bool,
    is_proline: bool,
}

impl Prepared {
    /// Reads `topology`'s molecule instances, orders its chains and finds
    /// each peptide residue's backbone atoms.
    pub(super) fn new(topology: &Topology, executor: &dyn Executor) -> Self {
        let hierarchy = topology.hierarchy();
        let mut represented = vec![false; topology.instance_count()];
        for (_, site) in hierarchy.atom_sites() {
            if let Some(flag) = represented.get_mut(site.atom().instance().index()) {
                *flag = true;
            }
        }
        let ignored_instances = topology
            .instances()
            .map(|molecule| molecule.id())
            .filter(|id| represented.get(id.index()) != Some(&true))
            .collect();

        // Hierarchy chains are already topology-global and may span molecule
        // instances, so no fragment reassembly is required here.
        let mut first_residue = 0;
        let mut logical_chains = hierarchy
            .chains()
            .enumerate()
            .map(|(discovery_order, (_, chain))| {
                let logical = LogicalChain {
                    chain,
                    discovery_order,
                    first_residue,
                };
                first_residue += chain.residues().len();
                logical
            })
            .collect::<Vec<_>>();
        // mkdssp builds one global residue table in label-asym order. Stable
        // shortlex ordering matches generated mmCIF label-asym series (A..Z,
        // AA..) while preserving discovery order for repeated logical chains.
        logical_chains.sort_by(|left, right| {
            let (left_label, right_label) = (left.chain.label_id(), right.chain.label_id());
            left_label
                .len()
                .cmp(&right_label.len())
                .then_with(|| left_label.cmp(right_label))
                .then_with(|| left.discovery_order.cmp(&right.discovery_order))
        });
        // Ordering a chain and reading a residue's atoms depend on nothing
        // else. Chains differ widely in size, and solvent residues are cheap,
        // so each work item takes a share of the atom sites, not of the
        // chains.
        let sizes = logical_chains
            .iter()
            .map(|logical| {
                logical
                    .chain
                    .residues()
                    .iter()
                    .filter_map(|&residue| hierarchy.residue(residue))
                    .map(|residue| residue.atom_sites().len())
                    .sum::<usize>()
                    .max(1)
            })
            .collect::<Vec<_>>();
        let prepared = map_weighted_ranges(executor, &sizes, |range| {
            logical_chains[range]
                .iter()
                .map(|logical| {
                    let mut residues = Vec::new();
                    let entries = order_chain(hierarchy, logical)
                        .iter()
                        .map(|entry| {
                            entry.peptide.then(|| {
                                residues.push(prepare_residue(topology, entry));
                                residues.len() - 1
                            })
                        })
                        .collect::<Vec<_>>();
                    (entries, residues)
                })
                .collect::<Vec<_>>()
        });
        let mut chains = Vec::with_capacity(logical_chains.len());
        let mut residues = Vec::new();
        for (entries, chain_residues) in prepared.into_iter().flatten() {
            let offset = residues.len();
            chains.push(
                entries
                    .into_iter()
                    .map(|entry| entry.map(|index| offset + index))
                    .collect(),
            );
            residues.extend(chain_residues);
        }
        let mut lookup = residues
            .iter()
            .filter(|residue| {
                residue
                    .atoms
                    .is_ok_and(|atoms| atoms.iter().all(Option::is_some))
            })
            .zip(0..)
            .map(|(residue, index)| (residue.key, index))
            .collect::<Vec<_>>();
        lookup.sort_unstable();
        Self {
            ignored_instances,
            chains,
            residues,
            lookup,
        }
    }

    /// How many residues assignment reads.
    pub(super) fn peptide_residues(&self) -> usize {
        self.residues.len()
    }
}

/// The residues of one logical chain in DSSP order.
fn order_chain<'a>(hierarchy: &'a Hierarchy, logical: &LogicalChain<'a>) -> Vec<ChainResidue<'a>> {
    let mut ordered = logical
        .chain
        .residues()
        .iter()
        .zip(logical.first_residue..)
        .map(|(&residue, discovery_order)| ChainResidue {
            hierarchy,
            chain: logical.chain,
            residue,
            discovery_order,
            peptide: hierarchy.residue(residue).is_some_and(|residue| {
                matches!(
                    residue.class(),
                    ResidueClass::AminoAcid | ResidueClass::Other
                )
            }),
        })
        .collect::<Vec<_>>();
    // Choose one ordering convention for the complete chain. Mixing label
    // and author numbering in a pairwise comparator can be non-transitive.
    // Native hierarchies need not carry either source numbering convention.
    let ordering_residues = ordered
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

/// Finds a peptide residue's backbone atoms by name.
fn prepare_residue(topology: &Topology, entry: &ChainResidue<'_>) -> PreparedResidue {
    let hierarchy = entry.hierarchy;
    let residue = hierarchy
        .residue(entry.residue)
        .expect("chain residues come from the model hierarchy");
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
        backbone[backbone_index] = Some(
            topology
                .atom_index(site.atom)
                .expect("hierarchy atoms are topology atoms")
                .index(),
        );
    }
    PreparedResidue {
        key: entry.residue,
        source_chain: entry.chain.id,
        label_sequence_id: residue.label_seq_id,
        atoms: repeated.map_or(Ok(backbone), Err),
        any_backbone_atom,
        amino_acid: residue.class() == ResidueClass::AminoAcid,
        is_proline: residue.name.eq_ignore_ascii_case("PRO"),
    }
}

/// A peptide residue's backbone in one structure.
#[derive(Debug)]
struct ResidueBackbone<'p> {
    prepared: &'p PreparedResidue,
    /// `N`, `CA`, `C`, `O`, or the first repeated one.
    atoms: Result<[Option<Vec3>; 4], DsspBackboneAtom>,
    /// Whether DSSP can index every backbone atom found.
    representable: bool,
}

/// Reads a prepared residue's backbone coordinates.
fn read_backbone<'p>(model: ModelView<'_>, prepared: &'p PreparedResidue) -> ResidueBackbone<'p> {
    let positions = model.positions();
    let atoms = prepared.atoms.map(|atoms| {
        atoms.map(|atom| {
            atom.map(|atom| {
                Vec3::from(
                    positions
                        .get(atom)
                        .expect("hierarchy atoms have model positions")
                        .value_in(ANGSTROM)
                        .expect("positions are lengths"),
                )
            })
        })
    });
    ResidueBackbone {
        prepared,
        atoms,
        representable: atoms.is_ok_and(|atoms| {
            atoms
                .iter()
                .flatten()
                .all(|point| point.is_dssp_representable())
        }),
    }
}

/// Builds DSSP's residue table of `model` from its prepared topology.
fn extract_backbones(
    prepared: &Prepared,
    model: ModelView<'_>,
    limits: DsspLimits,
    executor: &dyn Executor,
) -> Result<(Vec<BackboneResidue>, DsspReport, usize), DsspError> {
    let hierarchy = model.topology().hierarchy();
    let mut report = DsspReport {
        ignored_instances: prepared.ignored_instances.clone(),
        non_peptide_residues: 0,
        skipped_residues: Vec::new(),
        detected_gaps: 0,
        reconstructed_amide_hydrogens: 0,
        candidate_hydrogen_bond_pairs: 0,
        generated_ladders: 0,
    };
    let mut chain_segments = 0;

    // Reading a residue's backbone depends on nothing else; segments,
    // limits and errors are then decided in DSSP order.
    let backbones = map_ranges(executor, prepared.residues.len(), MIN_CHUNK, |range| {
        prepared.residues[range]
            .iter()
            .map(|residue| read_backbone(model, residue))
            .collect::<Vec<_>>()
    })
    .into_iter()
    .flatten()
    .collect::<Vec<_>>();
    let mut residues = Vec::with_capacity(prepared.lookup.len().min(limits.max_residues));
    let mut degenerate = None;
    for (chain_index, chain) in prepared.chains.iter().enumerate() {
        extract_chain(
            chain_index,
            chain
                .iter()
                .map(|residue| residue.map(|index| &backbones[index])),
            hierarchy,
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

/// Appends the complete residues of one chain to the residue table, in DSSP
/// order, and reports the others.
#[allow(clippy::too_many_arguments)]
fn extract_chain<'p>(
    chain_index: usize,
    backbones: impl Iterator<Item = Option<&'p ResidueBackbone<'p>>>,
    hierarchy: &Hierarchy,
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
        let prepared = backbone.prepared;
        let key = prepared.key;
        let atoms = backbone
            .atoms
            .map_err(|atom| DsspError::AmbiguousBackboneAtom { residue: key, atom })?;
        let missing = DsspBackboneAtom::ALL
            .into_iter()
            .enumerate()
            .filter_map(|(index, atom)| atoms[index].is_none().then_some(atom))
            .collect::<Vec<_>>();
        if !missing.is_empty() {
            if !prepared.any_backbone_atom && !prepared.amino_acid {
                report.non_peptide_residues += 1;
            } else {
                report.skipped_residues.push(DsspSkippedResidue {
                    key,
                    source: source(hierarchy, prepared.source_chain, key),
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
        if !backbone.representable {
            return Err(DsspError::CoordinateOutOfRange { residue: key });
        }

        let [n, ca, c, o] = atoms.map(|atom| atom.expect("checked complete backbone"));
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
            Some(preceding) if !prepared.is_proline => {
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
            source_chain: prepared.source_chain,
            label_sequence_id: prepared.label_sequence_id,
            chain: chain_index,
            n,
            ca,
            c,
            o,
            h,
            is_proline: prepared.is_proline,
            segment,
            prev,
            next: None,
            chain_break,
            geometry: ResidueGeometry::default(),
            acceptors: [BondSlot::default(); 2],
            donors: [BondSlot::default(); 2],
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

/// The source labels of residue `key` of chain `chain`.
fn source(hierarchy: &Hierarchy, chain: ChainId, key: ResidueId) -> DsspResidueSource {
    residue_source(
        hierarchy
            .chain(chain)
            .expect("analyzed residues come from the model hierarchy"),
        hierarchy
            .residue(key)
            .expect("analyzed residues come from the model hierarchy"),
    )
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
#[derive(Debug, Clone, Copy, Default)]
struct ResidueGeometry {
    phi: Option<f64>,
    psi: Option<f64>,
    omega: Option<f64>,
    alpha: Option<f64>,
    kappa: Option<f64>,
    tco: Option<f64>,
}

/// One range of residues' part of [`analyze_backbone`].
#[derive(Debug, Default)]
struct RangeAnalysis {
    /// The first residue of the range.
    start: usize,
    geometry: Vec<ResidueGeometry>,
    /// Residue `start + k`'s candidate partners after it are
    /// `later[later_starts[k]..later_starts[k + 1]]`, in increasing order.
    later_starts: Vec<usize>,
    later: Vec<usize>,
    /// For each pair with a later partner, the energies of the bonds from
    /// the first residue's N-H to the second's C=O and back.
    energies: Vec<[Option<f64>; 2]>,
    /// The same pairs by second residue: those ending in residue
    /// `start + k` are `earlier[earlier_starts[k]..earlier_starts[k + 1]]`,
    /// as the first residue and the energies, in increasing order of the
    /// first.
    earlier_starts: Vec<usize>,
    earlier: Vec<(usize, [Option<f64>; 2])>,
}

/// The candidate pairs: residues whose `CA` atoms lie within 9 Å, with the
/// energies of the hydrogen bonds between them, by range of residues.
#[derive(Debug)]
struct CandidatePairs {
    ranges: Vec<RangeAnalysis>,
    pairs: usize,
}

impl CandidatePairs {
    /// How many pairs there are.
    fn len(&self) -> usize {
        self.pairs
    }

    /// The range holding `residue`, and the residue's place in it.
    fn locate(&self, residue: usize) -> (&RangeAnalysis, usize) {
        let range = &self.ranges[self.ranges.partition_point(|range| range.start <= residue) - 1];
        (range, residue - range.start)
    }

    /// The partners after `residue`, in increasing order, with the energies
    /// of the bonds from its N-H to each partner's C=O and back.
    fn later(&self, residue: usize) -> (&[usize], &[[Option<f64>; 2]]) {
        let (range, local) = self.locate(residue);
        let span = range.later_starts[local]..range.later_starts[local + 1];
        (&range.later[span.clone()], &range.energies[span])
    }

    /// The pairs ending in `residue`, by range of their first residue in
    /// order: each slice holds a range's first residues, in increasing
    /// order, with the energies of the bonds from the first's N-H to
    /// `residue`'s C=O and back.
    fn earlier(&self, residue: usize) -> impl Iterator<Item = &[(usize, [Option<f64>; 2])]> {
        self.ranges
            .iter()
            .take_while(move |range| range.start < residue)
            .map(move |range| {
                let local = residue - range.start;
                &range.earlier[range.earlier_starts[local]..range.earlier_starts[local + 1]]
            })
    }
}

/// Computes every residue's backbone geometry, the candidate pairs and
/// their hydrogen-bond energies, installs the geometry and the bonds, and
/// returns the candidate pairs.
fn analyze_backbone(
    residues: &mut [BackboneResidue],
    max_candidate_pairs: usize,
    executor: &dyn Executor,
) -> Result<CandidatePairs, DsspError> {
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
    let candidates = CandidatePairs {
        ranges: ranges.into_iter().collect::<Result<_, _>>()?,
        pairs: found.into_inner(),
    };
    // DSSP inserts the bonds pair by pair in candidate order, which decides
    // its ties. A residue's own bonds then come from its pairs with earlier
    // residues and then from those with later ones, each in increasing
    // order of the partner, so every residue fills its slots on its own.
    for (index, residue) in residues.iter_mut().enumerate() {
        let (part, local) = candidates.locate(index);
        residue.geometry = part.geometry[local];
        for &(earlier, [forward, backward]) in candidates.earlier(index).flatten() {
            if let Some(energy) = forward {
                insert_bond(&mut residue.donors, earlier, energy);
            }
            if let Some(energy) = backward {
                insert_bond(&mut residue.acceptors, earlier, energy);
            }
        }
        let (later, energies) = candidates.later(index);
        for (&later, &[forward, backward]) in later.iter().zip(energies) {
            if let Some(energy) = forward {
                insert_bond(&mut residue.acceptors, later, energy);
            }
            if let Some(energy) = backward {
                insert_bond(&mut residue.donors, later, energy);
            }
        }
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
    let mut analysis = candidate_pairs(residues, grid, range.clone(), found, max_candidate_pairs)?;
    analysis.geometry = range
        .clone()
        .map(|index| residue_geometry(residues, index))
        .collect();
    analysis.energies = Vec::with_capacity(analysis.later.len());
    for (local, first) in range.enumerate() {
        let span = analysis.later_starts[local]..analysis.later_starts[local + 1];
        for &second in &analysis.later[span] {
            let forward = hydrogen_bond_energy(residues, first, second)?;
            let backward = if second != first + 1 {
                hydrogen_bond_energy(residues, second, first)?
            } else {
                None
            };
            analysis.energies.push([forward, backward]);
        }
    }
    // A counting sort by second residue keeps each one's pairs in order of
    // the first.
    let start = analysis.start;
    let mut starts = vec![0; residues.len() - start + 1];
    for &second in &analysis.later {
        starts[second - start + 1] += 1;
    }
    for index in 1..starts.len() {
        starts[index] += starts[index - 1];
    }
    let mut next = starts.clone();
    let mut earlier = vec![(0, [None; 2]); analysis.later.len()];
    for (local, ends) in analysis.later_starts.windows(2).enumerate() {
        for (&second, &energies) in analysis.later[ends[0]..ends[1]]
            .iter()
            .zip(&analysis.energies[ends[0]..ends[1]])
        {
            let slot = &mut next[second - start];
            earlier[*slot] = (start + local, energies);
            *slot += 1;
        }
    }
    analysis.earlier_starts = starts;
    analysis.earlier = earlier;
    Ok(analysis)
}

fn residue_geometry(residues: &[BackboneResidue], index: usize) -> ResidueGeometry {
    let residue = &residues[index];
    let mut geometry = ResidueGeometry::default();
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
        residues[previous_previous].label_sequence_id,
        residues[next_next].label_sequence_id,
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

/// The residues binned by the [`spatial_cell`] of their `CA`, with the
/// `CA` coordinates stored in cell order so that a column scan reads them
/// contiguously.
#[derive(Debug)]
struct CellGrid {
    cells: CellIndex,
    /// The residues in cell order.
    members: Vec<usize>,
    /// The `CA` of each of `members`.
    points: Vec<Vec3>,
}

/// Where each cell's residues lie in [`CellGrid::members`].
#[derive(Debug)]
enum CellIndex {
    /// Every cell of the residues' box, padded by one cell on each side, in
    /// x, y, z order: members `starts[slot]..starts[slot + 1]` lie in cell
    /// `slot`.
    Dense {
        origin: (i64, i64, i64),
        dims: (usize, usize, usize),
        starts: Vec<usize>,
    },
    /// Occupied cells only, for a box with many more cells than residues:
    /// member `k` lies in `cells[k]`, sorted by cell.
    Sparse { cells: Vec<(i64, i64, i64)> },
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
        let (cells, members) = match dense {
            None => {
                let mut sorted = cells.into_iter().zip(0..).collect::<Vec<_>>();
                sorted.sort_unstable();
                let (cells, members) = sorted.into_iter().unzip();
                (CellIndex::Sparse { cells }, members)
            }
            Some(dims) => {
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
                (
                    CellIndex::Dense {
                        origin,
                        dims,
                        starts,
                    },
                    members,
                )
            }
        };
        let points = members.iter().map(|&member| residues[member].ca).collect();
        Self {
            cells,
            members,
            points,
        }
    }

    /// The residues in cells `(x, y, z - 1)` through `(x, y, z + 1)`, where
    /// `(x, y, z)` is within one cell of a residue's cell, and their `CA`.
    fn column(&self, x: i64, y: i64, z: i64) -> (&[usize], &[Vec3]) {
        let span = match &self.cells {
            CellIndex::Dense {
                origin,
                dims,
                starts,
            } => {
                let base = (((x - origin.0) as usize * dims.1) + (y - origin.1) as usize) * dims.2;
                let first = base + (z - 1 - origin.2) as usize;
                starts[first]..starts[first + 3]
            }
            CellIndex::Sparse { cells } => {
                cells.partition_point(|&cell| cell < (x, y, z - 1))
                    ..cells.partition_point(|&cell| cell <= (x, y, z + 1))
            }
        };
        (&self.members[span.clone()], &self.points[span])
    }
}

/// The candidate partners after each residue in `range`: residues with
/// `CA` atoms within 9 Å, in increasing order. `found` counts the pairs
/// across ranges.
fn candidate_pairs(
    residues: &[BackboneResidue],
    grid: &CellGrid,
    range: Range<usize>,
    found: &AtomicUsize,
    max_candidate_pairs: usize,
) -> Result<RangeAnalysis, DsspError> {
    let maximum_squared = MAX_CA_PAIR_DISTANCE_ANGSTROMS.powi(2);
    let mut analysis = RangeAnalysis {
        start: range.start,
        later_starts: vec![0],
        ..RangeAnalysis::default()
    };
    let mut partners = Vec::new();
    for first in range {
        let ca = residues[first].ca;
        let (x, y, z) = spatial_cell(ca);
        partners.clear();
        for dx in -1..=1 {
            for dy in -1..=1 {
                let (members, points) = grid.column(x + dx, y + dy, z);
                partners.extend(members.iter().zip(points).filter_map(|(&second, &point)| {
                    (second > first && dssp_distance_squared(ca, point) <= maximum_squared)
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
        analysis.later.extend_from_slice(&partners);
        analysis.later_starts.push(analysis.later.len());
    }
    Ok(analysis)
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

/// The beta-sheet assignment of every residue.
#[derive(Debug)]
struct BetaSheets {
    residues: Vec<BetaResidue>,
    /// `E`, `B`, or loop.
    structure: Vec<DsspSecondaryStructure>,
    counts: BetaCounts,
}

fn calculate_beta_sheets(
    hierarchy: &Hierarchy,
    residues: &[BackboneResidue],
    candidates: &CandidatePairs,
    max_ladders: usize,
    executor: &dyn Executor,
) -> Result<BetaSheets, DsspError> {
    let found = map_ranges(executor, residues.len(), MIN_CHUNK, |range| {
        range
            .flat_map(|first| {
                candidates.later(first).0.iter().filter_map(move |&second| {
                    test_bridge(residues, first, second)
                        .map(|bridge_type| (first, second, bridge_type))
                })
            })
            .collect::<Vec<_>>()
    });
    let mut bridges = Vec::<Bridge>::new();
    // A bridge extends the ladder whose last pair is its predecessor, one
    // residue back on its first strand. Bridges come in candidate order, so
    // the predecessor is a bridge of the previous first residue, and no
    // other bridge has the same one. `previous` and `current` hold the
    // second residue, type and ladder of the bridges of residue `row - 1`
    // and `row`.
    let (mut previous, mut current) = (Vec::<(usize, BridgeType, usize)>::new(), Vec::new());
    let mut row = None;
    let mut counts = BetaCounts::default();
    for (first, second, bridge_type) in found.into_iter().flatten() {
        if row != Some(first) {
            if row.is_some_and(|row| row + 1 == first) {
                std::mem::swap(&mut previous, &mut current);
            } else {
                previous.clear();
            }
            current.clear();
            row = Some(first);
        }
        counts.bridges += 1;
        if residues[first].chain == residues[second].chain {
            counts.intra_chain += 1;
        } else {
            counts.inter_chain += 1;
        }
        let parallel = bridge_type == BridgeType::Parallel;
        let predecessor = if parallel {
            second.checked_sub(1)
        } else {
            second.checked_add(1)
        };
        let extended = predecessor
            .and_then(|j| {
                previous
                    .iter()
                    .find(|&&(second, other_type, _)| second == j && other_type == bridge_type)
            })
            .map(|&(_, _, ladder)| {
                let bridge = &mut bridges[ladder];
                bridge.i.push(first);
                if parallel {
                    bridge.j.push_back(second);
                } else {
                    bridge.j.push_front(second);
                }
                ladder
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
        current.push((second, bridge_type, ladder));
    }

    let chain_label = |residue: &BackboneResidue| {
        hierarchy
            .chain(residue.source_chain)
            .expect("analyzed residues come from the model hierarchy")
            .label_id()
    };
    bridges.sort_by(|left, right| {
        let left_residue = &residues[left.i[0]];
        let right_residue = &residues[right.i[0]];
        (chain_label(left_residue), left_residue.key.raw(), left.i[0]).cmp(&(
            chain_label(right_residue),
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
    let mut sheets = BetaSheets {
        residues: vec![BetaResidue::default(); residues.len()],
        structure: vec![DsspSecondaryStructure::Loop; residues.len()],
        counts,
    };
    assign_beta_results(residues, &bridges, &mut sheets);
    Ok(sheets)
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

fn assign_beta_results(residues: &[BackboneResidue], bridges: &[Bridge], sheets: &mut BetaSheets) {
    let beta = &mut sheets.residues;
    for bridge in bridges {
        let first_slot = bridge
            .i
            .iter()
            .any(|index| beta[*index].partners[0].is_some()) as usize;
        let second_slot = bridge
            .j
            .iter()
            .any(|index| beta[*index].partners[0].is_some()) as usize;
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
            beta[residue].partners[first_slot] = Some(BetaSlot {
                partner,
                ladder: bridge.ladder,
                parallel,
            });
        }
        for (&residue, partner) in bridge.j.iter().zip(second_partners) {
            beta[residue].partners[second_slot] = Some(BetaSlot {
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
        let first = bridge.i[0]..=*bridge.i.last().expect("bridge strand");
        let second = bridge.j[0]..=*bridge.j.back().expect("bridge strand");
        for strand in [first, second] {
            for index in strand {
                if sheets.structure[index] != DsspSecondaryStructure::ExtendedStrand {
                    sheets.structure[index] = structure;
                }
                beta[index].sheet = Some(bridge.sheet);
            }
        }
    }

    // Strands are numbered sheet by sheet, each sheet's in residue order.
    let mut members = vec![Vec::new(); bridges.len() + 1];
    for (index, residue) in beta.iter().enumerate() {
        if let Some(sheet) = residue.sheet {
            members[sheet].push(index);
        }
    }
    let mut strand = 0;
    for members in members {
        let mut previous: Option<usize> = None;
        for index in members {
            if previous.is_none_or(|prior| {
                prior + 1 != index || residues[prior].segment != residues[index].segment
            }) {
                strand += 1;
            }
            beta[index].strand = Some(strand);
            previous = Some(index);
        }
    }
}

/// How far, in residues, classification looks around a range: an
/// assignment depends on residues at most 10 before and 10 after it.
///
/// Whether a residue starts a helix depends only on its own turn, which
/// reaches 5 residues ahead; a helix position also on the turns starting up
/// to 5 residues before. The alpha-helix pass reads the starts of a residue
/// and the one before, from up to 3 residues back, so it reaches 4 each way.
/// The 3-10 pass reads the alpha assignments up to 2 residues each way, and
/// the pi pass the 3-10 assignments up to 4 each way. Turns, bends and
/// polyproline helices reach no further.
const CLASSIFY_HALO: usize = 10;

/// Helices, turns, bends and polyproline helices of the residues in
/// `range`, over the beta assignment `structure`.
///
/// The range is classified on a window reaching [`CLASSIFY_HALO`] residues
/// beyond it, which gives the range what DSSP gives it over the whole
/// table.
fn classify_range(
    residues: &[BackboneResidue],
    structure: &[DsspSecondaryStructure],
    min_polyproline_stretch: usize,
    range: Range<usize>,
) -> Vec<Classified> {
    let window =
        range.start.saturating_sub(CLASSIFY_HALO)..(range.end + CLASSIFY_HALO).min(residues.len());
    let mut classified = structure[window.clone()]
        .iter()
        .map(|&structure| (structure, [DsspHelixPosition::None; 4]))
        .collect::<Vec<_>>();
    calculate_helices_turns_and_bends(residues, window.start, &mut classified);
    calculate_polyproline_helices(
        residues,
        window.start,
        &mut classified,
        min_polyproline_stretch,
    );
    classified.drain(..range.start - window.start);
    classified.truncate(range.len());
    classified
}

/// The positions `global` of the residue table at which a pass over a
/// window of the table starting at residue `start` and holding `len`
/// residues examines a pattern reaching `before` residues back and `reach`
/// ahead.
fn window_starts(
    global: Range<usize>,
    start: usize,
    len: usize,
    before: usize,
    reach: usize,
) -> Range<usize> {
    global.start.max(start + before)..global.end.min((start + len).saturating_sub(reach))
}

/// DSSP's helix, turn and bend passes over the residues of the table from
/// `start` on that `window` holds, as they assign them over the whole table
/// to residues far enough inside the window; see [`CLASSIFY_HALO`].
fn calculate_helices_turns_and_bends(
    residues: &[BackboneResidue],
    start: usize,
    window: &mut [Classified],
) {
    let len = residues.len();
    let window_len = window.len();
    let starts = |global: Range<usize>, before, reach| {
        window_starts(global, start, window_len, before, reach)
    };
    let starts_helix = |window: &[Classified], index: usize, helix_index: usize| {
        is_helix_start(window[index - start].1[helix_index])
    };
    for helix_index in 0..3 {
        let stride = helix_index + 3;
        for first in starts(0..len.saturating_sub(stride), 0, stride) {
            let last = first + stride;
            if residues[first].segment != residues[last].segment
                || !test_bond(residues, last, first)
            {
                continue;
            }
            window[last - start].1[helix_index] = DsspHelixPosition::End;
            for (_, positions) in &mut window[first + 1 - start..last - start] {
                if positions[helix_index] == DsspHelixPosition::None {
                    positions[helix_index] = DsspHelixPosition::Middle;
                }
            }
            let positions = &mut window[first - start].1;
            positions[helix_index] = if positions[helix_index] == DsspHelixPosition::End {
                DsspHelixPosition::StartAndEnd
            } else {
                DsspHelixPosition::Start
            };
        }
    }
    for first in starts(1..len.saturating_sub(4), 1, 3) {
        if starts_helix(window, first, 1) && starts_helix(window, first - 1, 1) {
            for (structure, _) in &mut window[first - start..first + 4 - start] {
                *structure = DsspSecondaryStructure::AlphaHelix;
            }
        }
    }
    for first in starts(1..len.saturating_sub(3), 1, 2) {
        if starts_helix(window, first, 0) && starts_helix(window, first - 1, 0) {
            let helix = &mut window[first - start..first + 3 - start];
            let empty = helix.iter().all(|(structure, _)| {
                matches!(
                    structure,
                    DsspSecondaryStructure::Loop | DsspSecondaryStructure::Helix3_10
                )
            });
            if empty {
                for (structure, _) in helix {
                    *structure = DsspSecondaryStructure::Helix3_10;
                }
            }
        }
    }
    for first in starts(1..len.saturating_sub(5), 1, 4) {
        if starts_helix(window, first, 2) && starts_helix(window, first - 1, 2) {
            let helix = &mut window[first - start..first + 5 - start];
            let empty = helix.iter().all(|(structure, _)| {
                matches!(
                    structure,
                    DsspSecondaryStructure::Loop
                        | DsspSecondaryStructure::PiHelix
                        | DsspSecondaryStructure::AlphaHelix
                )
            });
            if empty {
                for (structure, _) in helix {
                    *structure = DsspSecondaryStructure::PiHelix;
                }
            }
        }
    }
    for index in starts(1..len.saturating_sub(1), 0, 0) {
        if window[index - start].0 != DsspSecondaryStructure::Loop {
            continue;
        }
        let turn = (0..3).any(|helix_index| {
            (1..helix_index + 3).any(|offset| {
                index >= start + offset && starts_helix(window, index - offset, helix_index)
            })
        });
        if turn {
            window[index - start].0 = DsspSecondaryStructure::Turn;
        } else if residues[index]
            .geometry
            .kappa
            .is_some_and(|kappa| kappa > 70.0)
        {
            window[index - start].0 = DsspSecondaryStructure::Bend;
        }
    }
}

fn is_helix_start(position: DsspHelixPosition) -> bool {
    matches!(
        position,
        DsspHelixPosition::Start | DsspHelixPosition::StartAndEnd
    )
}

/// DSSP's polyproline pass over a window, as
/// [`calculate_helices_turns_and_bends`].
fn calculate_polyproline_helices(
    residues: &[BackboneResidue],
    start: usize,
    window: &mut [Classified],
    stretch: usize,
) {
    const EPSILON: f64 = 29.0;
    const PHI_MIN: f64 = -75.0 - EPSILON;
    const PHI_MAX: f64 = -75.0 + EPSILON;
    const PSI_MIN: f64 = 145.0 - EPSILON;
    const PSI_MAX: f64 = 145.0 + EPSILON;

    let global = 1..residues.len().saturating_sub(stretch);
    for first in window_starts(global, start, window.len(), 0, stretch - 1) {
        let qualifies = residues[first..first + stretch].iter().all(|residue| {
            residue
                .geometry
                .phi
                .is_some_and(|phi| (PHI_MIN..=PHI_MAX).contains(&phi))
                && residue
                    .geometry
                    .psi
                    .is_some_and(|psi| (PSI_MIN..=PSI_MAX).contains(&psi))
        });
        if !qualifies {
            continue;
        }
        let helix = &mut window[first - start..first + stretch - start];
        let prior = helix[0].1[3];
        helix[0].1[3] = match (stretch, prior) {
            (2, DsspHelixPosition::End) => DsspHelixPosition::Middle,
            (3, DsspHelixPosition::End) => DsspHelixPosition::StartAndEnd,
            (_, DsspHelixPosition::None) => DsspHelixPosition::Start,
            _ => prior,
        };
        for (_, positions) in &mut helix[1..stretch - 1] {
            positions[3] = DsspHelixPosition::Middle;
        }
        helix[stretch - 1].1[3] = DsspHelixPosition::End;
        for (structure, _) in helix {
            if *structure == DsspSecondaryStructure::Loop {
                *structure = DsspSecondaryStructure::PolyprolineIIHelix;
            }
        }
    }
}

fn public_residue(
    hierarchy: &Hierarchy,
    residues: &[BackboneResidue],
    index: usize,
    beta: BetaResidue,
    (secondary_structure, helix_positions): Classified,
) -> DsspResidue {
    let convert_bond = |bond: BondSlot| {
        bond.partner.map(|partner| DsspHydrogenBond {
            partner: residues[partner].key,
            energy_kcal_per_mol: bond.energy,
        })
    };
    let convert_beta = |partner: Option<BetaSlot>| {
        partner.map(|partner| DsspBetaPartner {
            partner: residues[partner.partner].key,
            ladder: partner.ladder,
            direction: if partner.parallel {
                DsspBridgeDirection::Parallel
            } else {
                DsspBridgeDirection::Antiparallel
            },
        })
    };
    let residue = &residues[index];
    DsspResidue {
        key: residue.key,
        source: source(hierarchy, residue.source_chain, residue.key),
        secondary_structure,
        chain_break: residue.chain_break,
        phi: residue.geometry.phi,
        psi: residue.geometry.psi,
        omega: residue.geometry.omega,
        alpha: residue.geometry.alpha,
        kappa: residue.geometry.kappa,
        tco: residue.geometry.tco,
        acceptors: residue.acceptors.map(convert_bond),
        donors: residue.donors.map(convert_bond),
        beta_partners: beta.partners.map(convert_beta),
        sheet: beta.sheet,
        strand: beta.strand,
        helix_positions,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn serial() -> ScopedThreads {
        ScopedThreads::new(NonZeroUsize::MIN)
    }

    fn backbone(count: usize) -> Vec<BackboneResidue> {
        (0..count)
            .zip(0..=u32::MAX)
            .map(|(index, raw)| BackboneResidue {
                key: ResidueId::new(raw),
                source_chain: ChainId::new(0),
                label_sequence_id: None,
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
                geometry: ResidueGeometry::default(),
                acceptors: [BondSlot::default(); 2],
                donors: [BondSlot::default(); 2],
            })
            .collect()
    }

    /// Candidate pairs given by their first and second residues, without
    /// energies.
    fn candidates(count: usize, pairs: &[(usize, usize)]) -> CandidatePairs {
        let mut analysis = RangeAnalysis {
            later_starts: vec![0],
            ..RangeAnalysis::default()
        };
        for residue in 0..count {
            let later = pairs
                .iter()
                .filter(|&&(first, _)| first == residue)
                .map(|&(_, second)| second);
            analysis.later.extend(later);
            analysis.later_starts.push(analysis.later.len());
        }
        analysis.energies = vec![[None; 2]; analysis.later.len()];
        analysis.earlier_starts = vec![0; count + 1];
        CandidatePairs {
            ranges: vec![analysis],
            pairs: pairs.len(),
        }
    }

    fn classify_alone(residues: &[BackboneResidue], stretch: usize) -> Vec<Classified> {
        let structure = vec![DsspSecondaryStructure::Loop; residues.len()];
        classify_range(residues, &structure, stretch, 0..residues.len())
    }

    fn assert_helix(stride: usize, structure: DsspSecondaryStructure) {
        let mut residues = backbone(stride + 3);
        insert_bond(&mut residues[stride].acceptors, 0, -1.0);
        insert_bond(&mut residues[stride + 1].acceptors, 1, -1.0);
        let classified = classify_alone(&residues, 3);
        assert!(classified[1..=stride]
            .iter()
            .all(|&(assigned, _)| assigned == structure));
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
            residue.geometry.phi = Some(-75.0);
            residue.geometry.psi = Some(145.0);
        }
        let classified = classify_alone(&residues, 3);
        assert!(classified[1..=3]
            .iter()
            .all(|&(structure, _)| { structure == DsspSecondaryStructure::PolyprolineIIHelix }));
        assert_eq!(classified[1].1[3], DsspHelixPosition::Start);
        assert_eq!(classified[2].1[3], DsspHelixPosition::Middle);
        assert_eq!(classified[3].1[3], DsspHelixPosition::End);
    }

    #[test]
    fn two_residue_polyproline_mode_includes_the_last_complete_window() {
        let mut residues = backbone(7);
        for residue in &mut residues[4..=5] {
            residue.geometry.phi = Some(-75.0);
            residue.geometry.psi = Some(145.0);
        }
        let classified = classify_alone(&residues, 2);
        assert!(classified[4..=5]
            .iter()
            .all(|&(structure, _)| { structure == DsspSecondaryStructure::PolyprolineIIHelix }));
        assert_eq!(classified[4].1[3], DsspHelixPosition::Start);
        assert_eq!(classified[5].1[3], DsspHelixPosition::End);
    }

    /// A deterministic pseudo-random sequence.
    struct Lcg(u64);

    impl Lcg {
        fn below(&mut self, bound: u64) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            (self.0 >> 33) % bound
        }
    }

    /// Every range classified alone matches the whole table; a halo of 9
    /// residues fails here.
    #[test]
    fn classifying_ranges_matches_classifying_the_whole_table() {
        let count = 240;
        for seed in 0..64 {
            let mut random = Lcg(seed);
            let mut residues = backbone(count);
            let mut segment = 1;
            for residue in &mut residues {
                if random.below(60) == 0 {
                    segment += 1;
                }
                residue.segment = segment;
                residue.geometry.kappa = Some(random.below(180) as f64);
                residue.geometry.phi = Some(random.below(100) as f64 - 125.0);
                residue.geometry.psi = Some(random.below(100) as f64 + 95.0);
            }
            // Dense turns of every stride, so helices overlap and compete.
            let density = 2 + seed % 3;
            for (donor, residue) in residues.iter_mut().enumerate().skip(3) {
                for stride in 3..=5 {
                    if donor >= stride && random.below(density) == 0 {
                        insert_bond(&mut residue.acceptors, donor - stride, -1.0);
                    }
                }
            }
            let beta = (0..count)
                .map(|_| match random.below(16) {
                    0 => DsspSecondaryStructure::ExtendedStrand,
                    1 => DsspSecondaryStructure::BetaBridge,
                    _ => DsspSecondaryStructure::Loop,
                })
                .collect::<Vec<_>>();
            for stretch in [2, 3] {
                let whole = classify_range(&residues, &beta, stretch, 0..count);
                for size in [1, 2, 5, 16, 33] {
                    let pieces = (0..count)
                        .step_by(size)
                        .flat_map(|start| {
                            let range = start..(start + size).min(count);
                            classify_range(&residues, &beta, stretch, range)
                        })
                        .collect::<Vec<_>>();
                    assert_eq!(pieces, whole, "seed {seed}, stretch {stretch}, size {size}");
                }
            }
        }
    }

    #[test]
    fn ladder_limit_bounds_provisional_topology_before_bulge_merging() {
        let mut residues = backbone(10);
        insert_bond(&mut residues[3].acceptors, 5, -1.0);
        insert_bond(&mut residues[5].acceptors, 1, -1.0);
        insert_bond(&mut residues[5].acceptors, 7, -1.0);
        insert_bond(&mut residues[7].acceptors, 3, -1.0);

        assert!(matches!(
            calculate_beta_sheets(&Hierarchy::default(), &residues, &candidates(10, &[(2, 5), (4, 7)]), 1, &serial()),
            Err(DsspError::LimitExceeded(limit))
                if limit == LimitExceeded::first_beyond(Resource::Ladders, 1)
        ));
    }

    #[test]
    fn parallel_and_antiparallel_bridge_formulas_are_distinct() {
        let mut parallel = backbone(8);
        insert_bond(&mut parallel[3].acceptors, 5, -1.0);
        insert_bond(&mut parallel[5].acceptors, 1, -1.0);
        let sheets = calculate_beta_sheets(
            &Hierarchy::default(),
            &parallel,
            &candidates(8, &[(2, 5)]),
            1,
            &serial(),
        )
        .expect("bridge");
        assert_eq!(sheets.counts.bridges, 1);
        assert_eq!(sheets.structure[2], DsspSecondaryStructure::BetaBridge);
        assert!(sheets.residues[2].partners[0]
            .is_some_and(|partner| partner.parallel && partner.partner == 5));

        let mut antiparallel = backbone(8);
        insert_bond(&mut antiparallel[2].acceptors, 5, -1.0);
        insert_bond(&mut antiparallel[5].acceptors, 2, -1.0);
        let sheets = calculate_beta_sheets(
            &Hierarchy::default(),
            &antiparallel,
            &candidates(8, &[(2, 5)]),
            1,
            &serial(),
        )
        .expect("bridge");
        assert!(sheets.residues[2].partners[0]
            .is_some_and(|partner| !partner.parallel && partner.partner == 5));
        assert_eq!(
            calculate_beta_sheets(
                &Hierarchy::default(),
                &backbone(8),
                &candidates(8, &[]),
                0,
                &serial()
            )
            .expect("no bridge")
            .counts
            .ladders,
            0
        );
    }
}
