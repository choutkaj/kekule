use std::num::NonZeroUsize;

use super::*;
use crate::parallel::ScopedThreads;

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

/// Candidate pairs found through the grid match every pair of `CA`
/// atoms within 9 Å, by DSSP's float distance and in increasing order,
/// with pairs on 9 Å cell boundaries and exactly 9 Å apart.
#[test]
fn candidate_pairs_match_every_close_pair() {
    for seed in 0..8 {
        let mut random = Lcg(seed);
        let count = 400;
        let mut residues = backbone(count);
        let spread = [30, 60, 200, 4000][seed as usize % 4];
        for index in 0..count {
            let ca = if index % 5 == 1 {
                // Exactly 9 Å from the previous residue's `CA`.
                let previous = residues[index - 1].ca;
                Vec3 {
                    x: previous.x + 9.0,
                    ..previous
                }
            } else if index % 5 == 2 {
                // On a multiple of the 9 Å cell edge.
                Vec3 {
                    x: 9.0 * random.below(8) as f32,
                    y: 9.0 * random.below(8) as f32 - 36.0,
                    z: random.below(10) as f32 * 0.5,
                }
            } else {
                let mut axis = || random.below(spread * 100) as f32 * 0.01 - spread as f32 * 0.5;
                Vec3 {
                    x: axis(),
                    y: axis(),
                    z: axis(),
                }
            };
            residues[index].ca = ca;
        }
        let grid = ca_grid(&residues);
        let found = AtomicUsize::new(0);
        let analysis =
            candidate_pairs(&residues, &grid, 0..count, &found, usize::MAX, None).unwrap();
        let maximum_squared = MAX_CA_PAIR_DISTANCE_ANGSTROMS.powi(2);
        let mut later = Vec::new();
        let mut later_starts = vec![0];
        for first in 0..count {
            later.extend((first + 1..count).filter(|&second| {
                dssp_distance_squared(residues[first].ca, residues[second].ca) <= maximum_squared
            }));
            later_starts.push(later.len());
        }
        assert!(!later.is_empty(), "seed {seed}");
        assert_eq!(analysis.later, later, "seed {seed}");
        assert_eq!(analysis.later_starts, later_starts, "seed {seed}");
        assert_eq!(found.into_inner(), later.len(), "seed {seed}");
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
        calculate_beta_sheets(&Hierarchy::default(), &residues, &candidates(10, &[(2, 5), (4, 7)]), 1, &serial(), None),
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
        None,
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
        None,
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
            &serial(),
            None,
        )
        .expect("no bridge")
        .counts
        .ladders,
        0
    );
}
