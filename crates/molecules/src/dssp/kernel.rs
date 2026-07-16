use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::bio::{SmcraChain, SmcraChainId, SmcraHierarchy, SmcraResidue};
use crate::core::Point3;
use crate::modeling::{InstanceAtomId, Model, MoleculeInstanceId};

use super::*;

const MAX_PEPTIDE_BOND_ANGSTROMS: f64 = 2.5;
const MAX_CA_PAIR_DISTANCE_ANGSTROMS: f64 = 9.0;
const MIN_ELECTROSTATIC_DISTANCE_ANGSTROMS: f64 = 0.5;
const MIN_HBOND_ENERGY_KCAL_PER_MOL: f64 = -9.9;
const MAX_HBOND_ENERGY_KCAL_PER_MOL: f64 = -0.5;
const DSSP_COUPLING_KCAL_ANGSTROM_PER_MOL: f64 = -27.888_f32 as f64;

#[derive(Debug, Clone, Copy, PartialEq)]
struct Vec3 {
    x: f64,
    y: f64,
    z: f64,
}

impl From<Point3> for Vec3 {
    fn from(value: Point3) -> Self {
        Self {
            x: value.x,
            y: value.y,
            z: value.z,
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

    fn dot(self, other: Self) -> f64 {
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn cross(self, other: Self) -> Self {
        Self {
            x: self.y * other.z - other.y * self.z,
            y: self.z * other.x - other.z * self.x,
            z: self.x * other.y - other.x * self.y,
        }
    }

    fn norm_squared(self) -> f64 {
        self.dot(self)
    }

    fn is_dssp_representable(self) -> bool {
        [self.x, self.y, self.z].into_iter().all(|value| {
            let value = value as f32;
            let cell = (value / MAX_CA_PAIR_DISTANCE_ANGSTROMS as f32).floor();
            value.is_finite() && cell > i64::MIN as f32 && cell < i64::MAX as f32
        })
    }
}

// DSSP 4.6.1 stores points and evaluates point distances in `float`. These
// helpers preserve that narrow compatibility boundary while the public values
// and the rest of the molecules coordinate model remain `f64`.
fn dssp_distance_squared(first: Vec3, second: Vec3) -> f64 {
    let dx = (first.x as f32) - (second.x as f32);
    let dy = (first.y as f32) - (second.y as f32);
    let dz = (first.z as f32) - (second.z as f32);
    ((dx * dx + dy * dy) + dz * dz) as f64
}

fn dssp_distance(first: Vec3, second: Vec3) -> f64 {
    (dssp_distance_squared(first, second) as f32).sqrt() as f64
}

fn dssp_reconstructed_hydrogen(n: Vec3, preceding_c: Vec3, preceding_o: Vec3) -> Option<Vec3> {
    let dx = (preceding_c.x as f32) - (preceding_o.x as f32);
    let dy = (preceding_c.y as f32) - (preceding_o.y as f32);
    let dz = (preceding_c.z as f32) - (preceding_o.z as f32);
    let length = ((dx * dx + dy * dy) + dz * dz).sqrt();
    (length.is_finite() && length > f32::EPSILON).then(|| Vec3 {
        x: ((n.x as f32) + dx / length) as f64,
        y: ((n.y as f32) + dy / length) as f64,
        z: ((n.z as f32) + dz / length) as f64,
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
    key: DsspResidueKey,
    source: DsspResidueSource,
    chain: SmcraChainId,
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

pub(super) fn assign(model: &Model, options: DsspOptions) -> Result<DsspResult, DsspError> {
    if !matches!(options.min_polyproline_stretch, 2 | 3) {
        return Err(DsspError::InvalidPolyprolineStretch {
            value: options.min_polyproline_stretch,
        });
    }

    let (mut residues, mut report, chain_segments) = extract_backbones(model, options.limits)?;
    if residues.is_empty() {
        return Err(DsspError::NoAnalyzableProteinResidues);
    }

    calculate_geometry(&mut residues, &mut report)?;
    let candidates = candidate_pairs(&residues, options.limits.max_candidate_pairs)?;
    report.candidate_hydrogen_bond_pairs = candidates.len();
    calculate_hydrogen_bonds(&mut residues, &candidates)?;
    let beta_counts =
        calculate_beta_sheets(&mut residues, &candidates, options.limits.max_ladders)?;
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

    let public_residues = residues
        .iter()
        .map(|residue| public_residue(residue, &residues))
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
    model: &Model,
    limits: DsspLimits,
) -> Result<(Vec<BackboneResidue>, DsspReport, usize), DsspError> {
    let mut residues = Vec::new();
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

    for (molecule_id, instance) in model.topology().molecules() {
        let Some(macro_molecule) = instance.macro_molecule() else {
            report.ignored_instances.push(molecule_id);
            continue;
        };
        let hierarchy = macro_molecule.hierarchy();
        for (_, hierarchy_model) in hierarchy.models() {
            for &chain_id in &hierarchy_model.chains {
                let chain =
                    hierarchy
                        .chain(chain_id)
                        .map_err(|error| DsspError::InvalidHierarchy {
                            molecule: molecule_id,
                            message: error.to_string(),
                        })?;
                extract_chain(
                    model,
                    molecule_id,
                    hierarchy,
                    chain,
                    &mut residues,
                    &mut report,
                    &mut chain_segments,
                    limits,
                )?;
            }
        }
    }
    Ok((residues, report, chain_segments))
}

#[allow(clippy::too_many_arguments)]
fn extract_chain(
    model: &Model,
    molecule_id: MoleculeInstanceId,
    hierarchy: &SmcraHierarchy,
    chain: &SmcraChain,
    residues: &mut Vec<BackboneResidue>,
    report: &mut DsspReport,
    chain_segments: &mut usize,
    limits: DsspLimits,
) -> Result<(), DsspError> {
    let mut previous = None;
    let mut forced_break = false;

    for &residue_id in &chain.residues {
        let residue =
            hierarchy
                .residue(residue_id)
                .map_err(|error| DsspError::InvalidHierarchy {
                    molecule: molecule_id,
                    message: error.to_string(),
                })?;
        let key = DsspResidueKey::new(molecule_id, residue_id);
        let source = residue_source(chain, residue);
        let mut backbone = [None; 4];
        let mut seen_backbone = 0;

        for &site_id in &residue.atom_sites {
            let site =
                hierarchy
                    .atom_site(site_id)
                    .map_err(|error| DsspError::InvalidHierarchy {
                        molecule: molecule_id,
                        message: error.to_string(),
                    })?;
            let atom_name = site
                .metadata
                .label_atom_id
                .as_deref()
                .or(site.metadata.auth_atom_id.as_deref())
                .map(str::trim);
            let Some(backbone_index) = atom_name.and_then(backbone_index) else {
                continue;
            };
            seen_backbone += 1;
            if backbone[backbone_index].is_some() {
                return Err(DsspError::AmbiguousBackboneAtom {
                    residue: key,
                    atom_name: BACKBONE_NAMES[backbone_index],
                });
            }
            let atom = InstanceAtomId::new(molecule_id, site.atom);
            backbone[backbone_index] = Some(Vec3::from(
                model
                    .position(atom)
                    .map_err(|error| DsspError::InvalidHierarchy {
                        molecule: molecule_id,
                        message: error.to_string(),
                    })?
                    .into_value(),
            ));
        }

        let missing = BACKBONE_NAMES
            .iter()
            .enumerate()
            .filter_map(|(index, name)| backbone[index].is_none().then_some(*name))
            .collect::<Vec<_>>();
        if !missing.is_empty() {
            if seen_backbone == 0 && !is_standard_amino_acid(&residue.name) {
                report.non_peptide_residues += 1;
            } else {
                report.skipped_residues.push(DsspSkippedResidue {
                    key,
                    source,
                    reason: DsspSkipReason::MissingBackboneAtoms { missing },
                });
                forced_break = true;
            }
            continue;
        }
        if residues.len() >= limits.max_residues {
            return Err(DsspError::ResourceLimitExceeded {
                resource: DsspResource::Residues,
                limit: limits.max_residues,
            });
        }

        let n = backbone[0].expect("checked complete backbone");
        let ca = backbone[1].expect("checked complete backbone");
        let c = backbone[2].expect("checked complete backbone");
        let o = backbone[3].expect("checked complete backbone");
        if [n, ca, c, o]
            .into_iter()
            .any(|point| !point.is_dssp_representable())
        {
            return Err(DsspError::CoordinateOutOfRange {
                residue: key,
                quantity: "backbone coordinate",
            });
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
        residues.push(BackboneResidue {
            key,
            source,
            chain: chain.id,
            n,
            ca,
            c,
            o,
            h: None,
            is_proline: residue.name.eq_ignore_ascii_case("PRO"),
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

fn residue_source(chain: &SmcraChain, residue: &SmcraResidue) -> DsspResidueSource {
    DsspResidueSource {
        residue_name: residue.name.clone(),
        chain_label_id: chain.label_id.clone(),
        chain_author_id: chain.author_id.clone(),
        label_sequence_id: residue.label_seq_id,
        author_sequence_id: residue.author_seq_id.clone(),
        insertion_code: residue.insertion_code.clone(),
    }
}

fn is_standard_amino_acid(name: &str) -> bool {
    matches!(
        name.to_ascii_uppercase().as_str(),
        "ALA"
            | "ARG"
            | "ASN"
            | "ASP"
            | "CYS"
            | "GLN"
            | "GLU"
            | "GLY"
            | "HIS"
            | "ILE"
            | "LEU"
            | "LYS"
            | "MET"
            | "PHE"
            | "PRO"
            | "SER"
            | "THR"
            | "TRP"
            | "TYR"
            | "VAL"
    )
}

fn calculate_geometry(
    residues: &mut [BackboneResidue],
    report: &mut DsspReport,
) -> Result<(), DsspError> {
    for index in 0..residues.len() {
        // DSSP assigns each non-initial amide hydrogen from the immediately
        // preceding complete residue in its internal table, even across a
        // numbered gap or chain boundary. Chain-continuity checks still gate
        // torsions and local secondary-structure patterns below.
        if index > 0 && !residues[index].is_proline {
            let previous = index - 1;
            let hydrogen = dssp_reconstructed_hydrogen(
                residues[index].n,
                residues[previous].c,
                residues[previous].o,
            );
            if hydrogen.is_none() {
                return Err(DsspError::DegenerateBackboneGeometry {
                    residue: residues[index].key,
                    quantity: "preceding carbonyl bond",
                });
            }
            residues[index].h = hydrogen;
            report.reconstructed_amide_hydrogens += 1;
        }
        if let Some(previous) = residues[index].prev {
            residues[index].phi = dihedral(
                residues[previous].c,
                residues[index].n,
                residues[index].ca,
                residues[index].c,
            );
            residues[index].tco = cosine_between(
                residues[index].c.sub(residues[index].o),
                residues[previous].c.sub(residues[previous].o),
            );
        }
        if let Some(next) = residues[index].next {
            residues[index].psi = dihedral(
                residues[index].n,
                residues[index].ca,
                residues[index].c,
                residues[next].n,
            );
            residues[index].omega = dihedral(
                residues[index].ca,
                residues[index].c,
                residues[next].n,
                residues[next].ca,
            );
        }
        if let Some(previous) = residues[index].prev {
            if let Some(next) = residues[index].next {
                if let Some(next_next) = residues[next].next {
                    residues[index].alpha = dihedral(
                        residues[previous].ca,
                        residues[index].ca,
                        residues[next].ca,
                        residues[next_next].ca,
                    );
                }
            }
        }
    }

    for index in 0..residues.len() {
        let Some(previous) = residues[index].prev else {
            continue;
        };
        let Some(previous_previous) = residues[previous].prev else {
            continue;
        };
        let Some(next) = residues[index].next else {
            continue;
        };
        let Some(next_next) = residues[next].next else {
            continue;
        };
        let sequence_is_contiguous = match (
            residues[previous_previous].source.label_sequence_id,
            residues[next_next].source.label_sequence_id,
        ) {
            (Some(first), Some(last)) => first.checked_add(4) == Some(last),
            _ => true,
        };
        if sequence_is_contiguous {
            residues[index].kappa = angle_degrees(
                residues[index].ca.sub(residues[previous_previous].ca),
                residues[next_next].ca.sub(residues[index].ca),
            );
        }
    }
    Ok(())
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
    if x_norm_squared <= f64::EPSILON || y_norm_squared <= f64::EPSILON {
        return None;
    }
    let u = p.dot(x) / x_norm_squared.sqrt();
    let v = p.dot(y) / y_norm_squared.sqrt();
    (u != 0.0 || v != 0.0).then(|| v.atan2(u).to_degrees())
}

fn cosine_between(first: Vec3, second: Vec3) -> Option<f64> {
    let denominator = (first.norm_squared() * second.norm_squared()).sqrt();
    (denominator > f64::EPSILON).then(|| (first.dot(second) / denominator).clamp(-1.0, 1.0))
}

fn angle_degrees(first: Vec3, second: Vec3) -> Option<f64> {
    cosine_between(first, second).map(|cosine| cosine.acos().to_degrees())
}

fn candidate_pairs(
    residues: &[BackboneResidue],
    max_candidate_pairs: usize,
) -> Result<Vec<(usize, usize)>, DsspError> {
    let mut cells = BTreeMap::<(i64, i64, i64), Vec<usize>>::new();
    let mut pairs = Vec::new();
    let maximum_squared = MAX_CA_PAIR_DISTANCE_ANGSTROMS.powi(2);
    for (current, residue) in residues.iter().enumerate() {
        let cell = spatial_cell(residue.ca);
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let neighbor = (cell.0 + dx, cell.1 + dy, cell.2 + dz);
                    let Some(indices) = cells.get(&neighbor) else {
                        continue;
                    };
                    for &prior in indices {
                        if dssp_distance_squared(residues[prior].ca, residue.ca) <= maximum_squared
                        {
                            if pairs.len() >= max_candidate_pairs {
                                return Err(DsspError::ResourceLimitExceeded {
                                    resource: DsspResource::CandidatePairs,
                                    limit: max_candidate_pairs,
                                });
                            }
                            pairs.push((prior, current));
                        }
                    }
                }
            }
        }
        cells.entry(cell).or_default().push(current);
    }
    pairs.sort_unstable();
    Ok(pairs)
}

fn spatial_cell(point: Vec3) -> (i64, i64, i64) {
    (
        ((point.x as f32) / MAX_CA_PAIR_DISTANCE_ANGSTROMS as f32).floor() as i64,
        ((point.y as f32) / MAX_CA_PAIR_DISTANCE_ANGSTROMS as f32).floor() as i64,
        ((point.z as f32) / MAX_CA_PAIR_DISTANCE_ANGSTROMS as f32).floor() as i64,
    )
}

fn calculate_hydrogen_bonds(
    residues: &mut [BackboneResidue],
    candidates: &[(usize, usize)],
) -> Result<(), DsspError> {
    for &(first, second) in candidates {
        calculate_hydrogen_bond(residues, first, second)?;
        if second != first + 1 {
            calculate_hydrogen_bond(residues, second, first)?;
        }
    }
    Ok(())
}

fn calculate_hydrogen_bond(
    residues: &mut [BackboneResidue],
    donor: usize,
    acceptor: usize,
) -> Result<(), DsspError> {
    if residues[donor].is_proline {
        return Ok(());
    }
    let Some(hydrogen) = residues[donor].h else {
        return Ok(());
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
            quantity: "hydrogen-bond energy",
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
    energy = energy.max(MIN_HBOND_ENERGY_KCAL_PER_MOL);

    insert_bond(&mut residues[donor].acceptors, acceptor, energy);
    insert_bond(&mut residues[acceptor].donors, donor, energy);
    Ok(())
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
) -> Result<BetaCounts, DsspError> {
    let mut bridges = Vec::<Bridge>::new();
    let mut counts = BetaCounts::default();
    for &(first, second) in candidates {
        let Some(bridge_type) = test_bridge(residues, first, second) else {
            continue;
        };
        counts.bridges += 1;
        if residues[first].key.molecule() == residues[second].key.molecule()
            && residues[first].chain == residues[second].chain
        {
            counts.intra_chain += 1;
        } else {
            counts.inter_chain += 1;
        }
        let mut extended = false;
        for bridge in &mut bridges {
            if bridge.bridge_type != bridge_type
                || bridge
                    .i
                    .last()
                    .copied()
                    .and_then(|index| index.checked_add(1))
                    != Some(first)
            {
                continue;
            }
            let matches = match bridge_type {
                BridgeType::Parallel => {
                    bridge
                        .j
                        .back()
                        .copied()
                        .and_then(|index| index.checked_add(1))
                        == Some(second)
                }
                BridgeType::Antiparallel => {
                    bridge
                        .j
                        .front()
                        .copied()
                        .and_then(|index| index.checked_sub(1))
                        == Some(second)
                }
            };
            if matches {
                bridge.i.push(first);
                match bridge_type {
                    BridgeType::Parallel => bridge.j.push_back(second),
                    BridgeType::Antiparallel => bridge.j.push_front(second),
                }
                extended = true;
                break;
            }
        }
        if !extended {
            if bridges.len() >= max_ladders {
                return Err(DsspError::ResourceLimitExceeded {
                    resource: DsspResource::Ladders,
                    limit: max_ladders,
                });
            }
            bridges.push(Bridge {
                bridge_type,
                i: vec![first],
                j: VecDeque::from([second]),
                sheet: 0,
                ladder: 0,
            });
        }
    }

    bridges.sort_by_key(|bridge| {
        let residue = &residues[bridge.i[0]];
        (
            residue.key.molecule().raw(),
            residue.source.chain_label_id.clone(),
            bridge.i[0],
        )
    });
    merge_bulges(residues, &mut bridges);
    if bridges.len() > max_ladders {
        return Err(DsspError::ResourceLimitExceeded {
            resource: DsspResource::Ladders,
            limit: max_ladders,
        });
    }
    counts.ladders = bridges.len();
    assign_sheets(&mut bridges);
    assign_beta_results(residues, &bridges);
    Ok(counts)
}

fn merge_bulges(residues: &[BackboneResidue], bridges: &mut Vec<Bridge>) {
    let mut first = 0;
    while first < bridges.len() {
        let mut second = first + 1;
        while second < bridges.len() {
            let ibi = bridges[first].i[0];
            let iei = *bridges[first].i.last().expect("bridge has first strand");
            let jbi = bridges[first].j[0];
            let jei = *bridges[first].j.back().expect("bridge has second strand");
            let ibj = bridges[second].i[0];
            let iej = *bridges[second].i.last().expect("bridge has first strand");
            let jbj = bridges[second].j[0];
            let jej = *bridges[second].j.back().expect("bridge has second strand");

            let incompatible = bridges[first].bridge_type != bridges[second].bridge_type
                || residues[ibi].segment != residues[iej].segment
                || residues[jbi.min(jbj)].segment != residues[jei.max(jej)].segment
                || ibj.checked_sub(iei).is_none_or(|gap| gap >= 6)
                || (iei >= ibj && ibi <= iej);
            if incompatible {
                second += 1;
                continue;
            }
            let first_gap = ibj - iei;
            let bulge = match bridges[first].bridge_type {
                BridgeType::Parallel => jbj
                    .checked_sub(jei)
                    .is_some_and(|second_gap| (second_gap < 6 && first_gap < 3) || second_gap < 3),
                BridgeType::Antiparallel => jbi
                    .checked_sub(jej)
                    .is_some_and(|second_gap| (second_gap < 6 && first_gap < 3) || second_gap < 3),
            };
            if !bulge {
                second += 1;
                continue;
            }
            let other = bridges.remove(second);
            bridges[first].i.extend(other.i);
            match bridges[first].bridge_type {
                BridgeType::Parallel => bridges[first].j.extend(other.j),
                BridgeType::Antiparallel => {
                    for value in other.j.into_iter().rev() {
                        bridges[first].j.push_front(value);
                    }
                }
            }
        }
        first += 1;
    }
}

fn assign_sheets(bridges: &mut [Bridge]) {
    let mut remaining = (0..bridges.len()).collect::<BTreeSet<_>>();
    let mut sheet = 1;
    let mut ladder = 0;
    while let Some(first) = remaining.pop_first() {
        let mut component = BTreeSet::from([first]);
        loop {
            let linked = remaining.iter().copied().find(|candidate| {
                component
                    .iter()
                    .any(|existing| bridges_linked(&bridges[*existing], &bridges[*candidate]))
            });
            let Some(linked) = linked else {
                break;
            };
            remaining.remove(&linked);
            component.insert(linked);
        }
        for index in component {
            bridges[index].sheet = sheet;
            bridges[index].ladder = ladder;
            ladder += 1;
        }
        sheet += 1;
    }
}

fn bridges_linked(first: &Bridge, second: &Bridge) -> bool {
    first
        .i
        .iter()
        .any(|value| second.i.contains(value) || second.j.contains(value))
        || first
            .j
            .iter()
            .any(|value| second.i.contains(value) || second.j.contains(value))
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

    let maximum_sheet = residues.iter().filter_map(|residue| residue.sheet).max();
    let mut strand = 0;
    for sheet in 1..=maximum_sheet.unwrap_or_default() {
        let mut previous: Option<usize> = None;
        for index in 0..residues.len() {
            if residues[index].sheet != Some(sheet) {
                continue;
            }
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

fn public_residue(residue: &BackboneResidue, all: &[BackboneResidue]) -> DsspResidue {
    let convert_bond = |bond: BondSlot| {
        bond.partner.map(|partner| DsspHydrogenBond {
            partner: all[partner].key,
            energy_kcal_per_mol: bond.energy,
        })
    };
    let convert_beta = |partner: Option<BetaSlot>| {
        partner.map(|partner| DsspBetaPartner {
            partner: all[partner.partner].key,
            ladder: partner.ladder,
            parallel: partner.parallel,
        })
    };
    DsspResidue {
        key: residue.key,
        source: residue.source.clone(),
        secondary_structure: residue.secondary_structure,
        chain_break: residue.chain_break,
        phi_degrees: residue.phi,
        psi_degrees: residue.psi,
        omega_degrees: residue.omega,
        alpha_degrees: residue.alpha,
        kappa_degrees: residue.kappa,
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

    fn backbone(count: usize) -> Vec<BackboneResidue> {
        let mut residues = (0..count)
            .map(|index| BackboneResidue {
                key: DsspResidueKey::new(
                    MoleculeInstanceId::new(0),
                    crate::bio::SmcraResidueId::new(index as u32),
                ),
                source: DsspResidueSource {
                    residue_name: "ALA".to_owned(),
                    chain_label_id: "A".to_owned(),
                    chain_author_id: Some("A".to_owned()),
                    label_sequence_id: Some(index as i32 + 1),
                    author_sequence_id: Some((index + 1).to_string()),
                    insertion_code: None,
                },
                chain: SmcraChainId::new(0),
                n: Vec3 {
                    x: index as f64,
                    y: 0.0,
                    z: 0.0,
                },
                ca: Vec3 {
                    x: index as f64,
                    y: 1.0,
                    z: 0.0,
                },
                c: Vec3 {
                    x: index as f64,
                    y: 1.0,
                    z: 1.0,
                },
                o: Vec3 {
                    x: index as f64,
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
            calculate_beta_sheets(&mut residues, &[(2, 5), (4, 7)], 1),
            Err(DsspError::ResourceLimitExceeded {
                resource: DsspResource::Ladders,
                limit: 1,
            })
        ));
    }

    #[test]
    fn parallel_and_antiparallel_bridge_formulas_are_distinct() {
        let mut parallel = backbone(8);
        insert_bond(&mut parallel[3].acceptors, 5, -1.0);
        insert_bond(&mut parallel[5].acceptors, 1, -1.0);
        let counts = calculate_beta_sheets(&mut parallel, &[(2, 5)], 1).expect("bridge");
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
        calculate_beta_sheets(&mut antiparallel, &[(2, 5)], 1).expect("bridge");
        assert!(antiparallel[2].beta_partners[0]
            .is_some_and(|partner| !partner.parallel && partner.partner == 5));
        assert_eq!(
            calculate_beta_sheets(&mut backbone(8), &[], 0)
                .expect("no bridge")
                .ladders,
            0
        );
    }
}
