//! Focused regressions: feature classification and every interaction type at
//! the edges of its criteria, on toy molecules placed by hand.

use std::sync::Arc;

use super::*;
use crate::geometry::{Point3, Vector3};
use crate::structure::{Model, Positions};
use crate::topology::{AtomSelection, Topology, TopologyAtomIndex};
use crate::units::{Quantity, ANGSTROM, DEGREE, NANOMETER};

/// Just inside and just outside a threshold.
const INSIDE: f64 = 1.0 - 1e-9;
const OUTSIDE: f64 = 1.0 + 1e-9;
/// Angles just inside and outside a threshold, in degrees.
const NUDGE: f64 = 1e-6;

fn model(smiles: &str, points: &[Point3]) -> (Arc<Topology>, Model) {
    let topology = crate::smiles::to_topology(smiles).expect("valid SMILES");
    assert_eq!(
        topology.atom_count(),
        points.len(),
        "one point per atom of {smiles}"
    );
    let positions = Positions::new(Quantity::new(points.to_vec(), ANGSTROM)).expect("finite");
    let model = Model::new(Arc::clone(&topology), positions).expect("bound");
    (topology, model)
}

fn select(topology: &Arc<Topology>, atoms: impl IntoIterator<Item = usize>) -> AtomSelection {
    AtomSelection::from_indices(
        Arc::clone(topology),
        atoms
            .into_iter()
            .map(|atom| TopologyAtomIndex::new(atom as u32)),
    )
    .expect("valid atoms")
}

fn find_with(
    smiles: &str,
    points: &[Point3],
    first: impl IntoIterator<Item = usize>,
    second: impl IntoIterator<Item = usize>,
    criteria: &InteractionCriteria,
) -> Interactions {
    let (topology, model) = model(smiles, points);
    let features = InteractionFeatures::new(Arc::clone(&topology));
    let first = select(&topology, first);
    let second = select(&topology, second);
    find_interactions_within(&features, &model, &first, &second, criteria).expect("search")
}

fn count(smiles: &str, points: &[Point3], split: usize, kind: InteractionKind) -> usize {
    let all = points.len();
    find_with(
        smiles,
        points,
        0..split,
        split..all,
        &InteractionCriteria::default(),
    )
    .of_kind(kind)
    .count()
}

fn p(x: f64, y: f64, z: f64) -> Point3 {
    Point3::new(x, y, z)
}

/// `length` from `from` in the xy plane at `degrees` from +x.
fn polar(from: Point3, length: f64, degrees: f64) -> Point3 {
    let radians = degrees.to_radians();
    from + Vector3::new(length * radians.cos(), length * radians.sin(), 0.0)
}

/// A regular hexagon of radius 1.39 Å around `center`, in the plane spanned
/// by the unit vectors `u` and `v`.
fn hexagon(center: Point3, u: Vector3, v: Vector3) -> Vec<Point3> {
    (0..6)
        .map(|k| {
            let angle = f64::from(k) * std::f64::consts::PI / 3.0;
            center + u * (1.39 * angle.cos()) + v * (1.39 * angle.sin())
        })
        .collect()
}

fn degrees(value: Quantity<f64>) -> f64 {
    value.value_in(DEGREE).unwrap()
}

fn angstrom(value: Quantity<f64>) -> f64 {
    value.value_in(ANGSTROM).unwrap()
}

#[test]
fn classifies_features() {
    let features =
        |smiles: &str| InteractionFeatures::new(crate::smiles::to_topology(smiles).unwrap());
    let groups = |smiles: &str| {
        features(smiles)
            .charged_groups()
            .map(|group| (group.sign(), group.kind(), group.charged_atoms().len()))
            .collect::<Vec<_>>()
    };
    use ChargeSign::{Negative, Positive};
    use ChargedGroupKind as K;
    assert_eq!(groups("CC(=O)O"), [(Negative, K::Carboxylate, 2)]);
    assert_eq!(groups("CC(=O)[O-]"), [(Negative, K::Carboxylate, 2)]);
    assert_eq!(groups("COP(=O)(O)O"), [(Negative, K::Phosphate, 3)]);
    assert_eq!(groups("CS(=O)(=O)O"), [(Negative, K::Sulfonate, 3)]);
    assert_eq!(groups("CCN"), [(Positive, K::Amine, 1)]);
    assert_eq!(groups("CC[NH3+]"), [(Positive, K::Amine, 1)]);
    assert_eq!(groups("C[N+](C)(C)C"), [(Positive, K::Amine, 1)]);
    assert_eq!(groups("CC(=N)N"), [(Positive, K::Amidine, 2)]);
    assert_eq!(groups("CNC(=N)N"), [(Positive, K::Guanidine, 3)]);
    assert_eq!(groups("C[S+](C)C"), [(Positive, K::FormalCharge, 1)]);
    // Amides, anilines, esters, and nitro groups carry no charge.
    assert!(groups("CC(=O)NC").is_empty());
    assert!(groups("Nc1ccccc1").is_empty());
    assert!(groups("CC(=O)OC").is_empty());
    assert!(groups("C[N+](=O)[O-]").is_empty());
    // A metal bonded to a carboxylate's oxygen leaves it a carboxylate.
    let bound = features("CC(=O)O[Mg+]");
    assert_eq!(
        bound
            .charged_groups()
            .map(|group| group.kind())
            .collect::<Vec<_>>(),
        [K::Carboxylate]
    );
    // Metals are metals, not charged groups.
    assert!(groups("[Zn+2]").is_empty());
    assert_eq!(features("[Zn+2]").metals().len(), 1);

    // An amide nitrogen donates but does not accept; its oxygen accepts.
    let amide = features("CC(=O)NC");
    assert_eq!(
        amide
            .donors()
            .map(|donor| donor.atom().index())
            .collect::<Vec<_>>(),
        [3]
    );
    assert_eq!(
        amide
            .acceptors()
            .iter()
            .map(|atom| atom.index())
            .collect::<Vec<_>>(),
        [2]
    );
    // Pyridine's nitrogen accepts; pyrrole's donates.
    assert_eq!(features("c1ccncc1").acceptors().len(), 1);
    let pyrrole = features("c1cc[nH]c1");
    assert!(pyrrole.acceptors().is_empty());
    assert_eq!(pyrrole.donors().len(), 1);
    // A carboxylic acid's oxygens accept but do not donate.
    let acid = features("CC(=O)O");
    assert_eq!(acid.donors().len(), 0);
    assert_eq!(acid.acceptors().len(), 2);
    // Explicit hydrogens are listed; implicit ones counted.
    let water = features("[H]O[H]");
    let donor = water.donors().next().unwrap();
    assert_eq!(donor.hydrogens().len(), 2);
    assert_eq!(donor.implicit_hydrogen_count(), 0);
    assert_eq!(
        features("O")
            .donors()
            .next()
            .unwrap()
            .implicit_hydrogen_count(),
        2
    );

    // Rings: naphthalene has two; cyclohexane none.
    assert_eq!(features("c1ccc2ccccc2c1").aromatic_rings().len(), 2);
    assert_eq!(features("C1CCCCC1").aromatic_rings().len(), 0);
    // Halogen-bond donors are C–Cl, C–Br, and C–I.
    let halogens = features("FC(Cl)(Br)I.[Cl-]");
    assert_eq!(halogens.halogen_donors().len(), 3);
    // Hydrophobic carbons are bonded only to carbon; the chloride is a
    // metal ligand.
    assert_eq!(halogens.hydrophobic_atoms().len(), 0);
    assert_eq!(halogens.metal_ligands().len(), 1);
    assert_eq!(features("CCO").hydrophobic_atoms().len(), 1);
    assert!(features("CCO").unperceived_definitions().is_empty());
}

#[test]
fn ring_center_and_normal() {
    let ring = hexagon(
        p(1.0, 2.0, 3.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
    );
    let (topology, model) = model("c1ccccc1", &ring);
    let features = InteractionFeatures::new(topology);
    let (center, normal) = features
        .aromatic_rings()
        .next()
        .unwrap()
        .center_and_normal(&model)
        .unwrap();
    let center = center.value_in(NANOMETER).unwrap();
    assert!((center - p(0.1, 0.2, 0.3)).norm() < 1e-12);
    assert!((normal.z.abs() - 1.0).abs() < 1e-12);
}

#[test]
fn hydrogen_bond_distance_without_hydrogens() {
    for (scale, expected) in [(INSIDE, 1), (OUTSIDE, 0)] {
        let points = [p(0.0, 0.0, 0.0), p(3.5 * scale, 0.0, 0.0)];
        let found = find_with("O.O", &points, [0], [1], &InteractionCriteria::default());
        assert_eq!(
            found.of_kind(InteractionKind::HydrogenBond).count(),
            expected
        );
    }
    // Each water could donate to the other: the bond is reported once, the
    // first end donating.
    let points = [p(0.0, 0.0, 0.0), p(2.9, 0.0, 0.0)];
    let found = find_with(
        "O.O",
        &points,
        [0, 1],
        [0, 1],
        &InteractionCriteria::default(),
    );
    let bonds = found
        .of_kind(InteractionKind::HydrogenBond)
        .collect::<Vec<_>>();
    assert_eq!(bonds.len(), 1);
    assert_eq!(bonds[0].first().role(), InteractionRole::Donor);
    assert_eq!(bonds[0].first().atoms(), [TopologyAtomIndex::new(0)]);
    assert!((angstrom(bonds[0].distance()) - 2.9).abs() < 1e-12);
    assert!(bonds[0].hydrogen_distance().is_none());
}

#[test]
fn hydrogen_bond_donor_and_acceptor_angles() {
    // Methanol donates to formaldehyde, which cannot donate back. The
    // X–D…A angle is set at the methanol oxygen, the D…A–Y angle at the
    // formaldehyde oxygen.
    for (degrees_at_donor, expected) in [(90.0 + NUDGE, 1), (90.0 - NUDGE, 0)] {
        let donor = p(0.0, 0.0, 0.0);
        let carbon = polar(donor, 1.43, degrees_at_donor);
        let acceptor = p(3.0, 0.0, 0.0);
        // Atom order: C, O (methanol), C, O (formaldehyde).
        let points = [carbon, donor, p(4.2, 0.0, 0.0), acceptor];
        assert_eq!(
            count("CO.C=O", &points, 2, InteractionKind::HydrogenBond),
            expected
        );
    }
    for (degrees_at_acceptor, expected) in [(90.0 + NUDGE, 1), (90.0 - NUDGE, 0)] {
        let donor = p(0.0, 0.0, 0.0);
        let acceptor = p(3.0, 0.0, 0.0);
        let carbon = polar(acceptor, 1.2, 180.0 - degrees_at_acceptor);
        let points = [p(-1.43, 0.0, 0.0), donor, carbon, acceptor];
        let found = find_with(
            "CO.C=O",
            &points,
            [0, 1],
            [2, 3],
            &InteractionCriteria::default(),
        );
        let bonds = found
            .of_kind(InteractionKind::HydrogenBond)
            .collect::<Vec<_>>();
        assert_eq!(bonds.len(), expected);
        if let Some(bond) = bonds.first() {
            let angles = bond.angles().collect::<Vec<_>>();
            assert_eq!(angles[0].0, InteractionAngle::NeighborDonorAcceptor);
            assert!((degrees(angles[0].1) - 180.0).abs() < 1e-9);
            assert_eq!(angles[1].0, InteractionAngle::DonorAcceptorNeighbor);
            assert!((degrees(angles[1].1) - degrees_at_acceptor).abs() < 1e-6);
        }
    }
}

#[test]
fn hydrogen_bond_with_explicit_hydrogens() {
    // Water with explicit hydrogens donates to formaldehyde along +x.
    // Atom order: H, O, H, C, O.
    let place = |hydrogen_reach: f64, hydrogen_degrees: f64| {
        let hydrogen = p(0.0, 0.0, 0.0);
        let acceptor = p(hydrogen_reach, 0.0, 0.0);
        let donor = polar(hydrogen, 0.96, hydrogen_degrees);
        let other = polar(donor, 0.96, hydrogen_degrees + 75.0);
        [
            hydrogen,
            donor,
            other,
            acceptor + Vector3::new(1.2, 0.0, 0.0),
            acceptor,
        ]
    };
    // The hydrogen–acceptor distance at 2.5 Å, the D–H…A angle 180°.
    for (scale, expected) in [(INSIDE, 1), (OUTSIDE, 0)] {
        let points = place(2.5 * scale, 180.0);
        assert_eq!(
            count("[H]O[H].C=O", &points, 3, InteractionKind::HydrogenBond),
            expected
        );
    }
    // The D–H…A angle at 120°, the hydrogen 1.8 Å from the acceptor.
    for (angle, expected) in [(120.0 + NUDGE, 1), (120.0 - NUDGE, 0)] {
        let points = place(1.8, angle);
        let found = find_with(
            "[H]O[H].C=O",
            &points,
            0..3,
            3..5,
            &InteractionCriteria::default(),
        );
        let bonds = found
            .of_kind(InteractionKind::HydrogenBond)
            .collect::<Vec<_>>();
        assert_eq!(bonds.len(), expected);
        if let Some(bond) = bonds.first() {
            assert_eq!(bond.first().hydrogen(), Some(TopologyAtomIndex::new(0)));
            assert!((angstrom(bond.hydrogen_distance().unwrap()) - 1.8).abs() < 1e-9);
            let (kind, value) = bond.angles().next().unwrap();
            assert_eq!(kind, InteractionAngle::DonorHydrogenAcceptor);
            assert!((degrees(value) - angle).abs() < 1e-6);
        }
    }
}

#[test]
fn salt_bridge_distance() {
    // Acetate's oxygens center at x = 0.6 Å; methylammonium's nitrogen is d
    // beyond it. Atom order: C, C, O, O, C, N.
    for (scale, expected) in [(INSIDE, 1), (OUTSIDE, 0)] {
        let d = 4.0 * scale;
        let points = [
            p(-1.5, 0.0, 0.0),
            p(0.0, 0.0, 0.0),
            p(0.6, 1.0, 0.0),
            p(0.6, -1.0, 0.0),
            p(0.6 + d + 1.47, 0.0, 0.0),
            p(0.6 + d, 0.0, 0.0),
        ];
        let found = find_with(
            "CC(=O)[O-].C[NH3+]",
            &points,
            0..4,
            4..6,
            &InteractionCriteria::default(),
        );
        let bridges = found
            .of_kind(InteractionKind::SaltBridge)
            .collect::<Vec<_>>();
        assert_eq!(bridges.len(), expected);
        if let Some(bridge) = bridges.first() {
            assert_eq!(bridge.first().role(), InteractionRole::Anion);
            assert_eq!(bridge.first().atoms().len(), 3);
            let center = bridge.first().center().value_in(ANGSTROM).unwrap();
            assert!((center - p(0.6, 0.0, 0.0)).norm() < 1e-9);
        }
    }
}

fn two_benzenes(second: Vec<Point3>) -> Vec<Point3> {
    let mut points = hexagon(
        p(0.0, 0.0, 0.0),
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 1.0, 0.0),
    );
    points.extend(second);
    points
}

/// The second ring at `center`, its plane tilted by `tilt` degrees about x.
fn tilted(center: Point3, tilt: f64) -> Vec<Point3> {
    let tilt = tilt.to_radians();
    hexagon(
        center,
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, tilt.cos(), tilt.sin()),
    )
}

#[test]
fn pi_stacking_edges() {
    let benzenes = "c1ccccc1.c1ccccc1";
    let stacking = |points: &[Point3]| {
        let found = find_with(
            benzenes,
            points,
            0..6,
            6..12,
            &InteractionCriteria::default(),
        );
        found
            .interactions()
            .iter()
            .filter(|interaction| {
                matches!(
                    interaction.kind(),
                    InteractionKind::PiStackingParallel | InteractionKind::PiStackingTShaped
                )
            })
            .map(|interaction| (interaction.kind(), interaction.clone()))
            .collect::<Vec<_>>()
    };
    use InteractionKind::{PiStackingParallel as Parallel, PiStackingTShaped as TShaped};
    // Distance, face to face.
    for (scale, expected) in [(INSIDE, 1), (OUTSIDE, 0)] {
        let points = two_benzenes(tilted(p(0.0, 0.0, 5.5 * scale), 0.0));
        assert_eq!(stacking(&points).len(), expected);
    }
    // Offset, 3.5 Å apart.
    for (scale, expected) in [(INSIDE, 1), (OUTSIDE, 0)] {
        let points = two_benzenes(tilted(p(2.0 * scale, 0.0, 3.5), 0.0));
        let found = stacking(&points);
        assert_eq!(found.len(), expected);
        if let Some((kind, interaction)) = found.first() {
            assert_eq!(*kind, Parallel);
            assert!((angstrom(interaction.offset().unwrap()) - 2.0 * scale).abs() < 1e-9);
            assert!(interaction.first().normal().is_some());
        }
    }
    // Angles: parallel up to 30°, T-shaped from 60°, nothing between.
    let at = |tilt: f64| stacking(&two_benzenes(tilted(p(0.0, 0.0, 3.8), tilt)));
    assert_eq!(
        at(30.0 - NUDGE).first().map(|found| found.0),
        Some(Parallel)
    );
    assert!(at(30.0 + NUDGE).is_empty());
    assert!(at(60.0 - NUDGE).is_empty());
    assert_eq!(at(60.0 + NUDGE).first().map(|found| found.0), Some(TShaped));
    let edge = at(90.0);
    assert_eq!(edge[0].0, TShaped);
    let (kind, value) = edge[0].1.angles().next().unwrap();
    assert_eq!(kind, InteractionAngle::RingPlanes);
    assert!((degrees(value) - 90.0).abs() < 1e-9);
}

#[test]
fn cation_pi_edges() {
    // Benzene in the xy plane; methylammonium's nitrogen over it. Atom
    // order: six ring carbons, C, N.
    let place = |offset: f64, height: f64| {
        let mut points = hexagon(
            p(0.0, 0.0, 0.0),
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 1.0, 0.0),
        );
        points.push(p(offset, 0.0, height + 1.47));
        points.push(p(offset, 0.0, height));
        points
    };
    let found = |points: &[Point3]| count("c1ccccc1.C[NH3+]", points, 6, InteractionKind::CationPi);
    for (scale, expected) in [(INSIDE, 1), (OUTSIDE, 0)] {
        assert_eq!(found(&place(0.0, 6.0 * scale)), expected);
        assert_eq!(found(&place(2.0 * scale, 3.0)), expected);
    }
    let points = place(1.0, 3.0);
    let result = find_with(
        "c1ccccc1.C[NH3+]",
        &points,
        0..6,
        6..8,
        &InteractionCriteria::default(),
    );
    let interaction = result.of_kind(InteractionKind::CationPi).next().unwrap();
    // The cation comes before the ring when either could be first; here
    // the ring is in the first set.
    assert_eq!(interaction.first().role(), InteractionRole::Ring);
    let (kind, value) = interaction.angles().next().unwrap();
    assert_eq!(kind, InteractionAngle::RingNormalCation);
    assert!((degrees(value) - (1.0f64 / 3.0).atan().to_degrees()).abs() < 1e-9);
}

#[test]
fn halogen_bond_edges() {
    // Bromomethane's bromine at the origin, its carbon along -x;
    // formaldehyde's oxygen accepts. Atom order: C, Br, C, O.
    let place = |reach: f64, donor_angle: f64, acceptor_angle: f64| {
        let halogen = p(0.0, 0.0, 0.0);
        let acceptor = polar(halogen, reach, 180.0 - donor_angle);
        // The acceptor's carbon at `acceptor_angle` from the line back to
        // the bromine.
        let back = 180.0 + (180.0 - donor_angle);
        let carbon = polar(acceptor, 1.2, back + acceptor_angle);
        [p(-1.94, 0.0, 0.0), halogen, carbon, acceptor]
    };
    let found = |points: &[Point3]| count("CBr.C=O", points, 2, InteractionKind::HalogenBond);
    for (scale, expected) in [(INSIDE, 1), (OUTSIDE, 0)] {
        assert_eq!(found(&place(4.0 * scale, 180.0, 120.0)), expected);
    }
    assert_eq!(found(&place(3.2, 135.0 + NUDGE, 120.0)), 1);
    assert_eq!(found(&place(3.2, 135.0 - NUDGE, 120.0)), 0);
    assert_eq!(found(&place(3.2, 180.0, 90.0 + NUDGE)), 1);
    assert_eq!(found(&place(3.2, 180.0, 90.0 - NUDGE)), 0);
    assert_eq!(found(&place(3.2, 180.0, 150.0 - NUDGE)), 1);
    assert_eq!(found(&place(3.2, 180.0, 150.0 + NUDGE)), 0);

    let points = place(3.2, 170.0, 120.0);
    let result = find_with(
        "CBr.C=O",
        &points,
        0..2,
        2..4,
        &InteractionCriteria::default(),
    );
    let bond = result.of_kind(InteractionKind::HalogenBond).next().unwrap();
    let angles = bond
        .angles()
        .map(|(kind, value)| (kind, degrees(value)))
        .collect::<Vec<_>>();
    assert_eq!(angles[0].0, InteractionAngle::CarbonHalogenAcceptor);
    assert!((angles[0].1 - 170.0).abs() < 1e-6);
    assert_eq!(angles[1].0, InteractionAngle::HalogenAcceptorNeighbor);
    assert!((angles[1].1 - 120.0).abs() < 1e-6);
}

#[test]
fn metal_coordination_and_hydrophobic_edges() {
    for (scale, expected) in [(INSIDE, 1), (OUTSIDE, 0)] {
        let points = [p(0.0, 0.0, 0.0), p(3.0 * scale, 0.0, 0.0)];
        assert_eq!(
            count("[Zn+2].O", &points, 1, InteractionKind::MetalCoordination),
            expected
        );
        // Ethane's second carbon is far from the other ethane.
        let points = [
            p(-1.5, 0.0, 0.0),
            p(0.0, 0.0, 0.0),
            p(4.0 * scale, 0.0, 0.0),
            p(4.0 * scale + 1.5, 0.0, 0.0),
        ];
        assert_eq!(
            count("CC.CC", &points, 2, InteractionKind::HydrophobicContact),
            expected
        );
    }
    // A metal bonded to its ligand still coordinates it.
    let points = [p(-2.0, 0.0, 0.0), p(0.0, 0.0, 0.0), p(2.0, 0.0, 0.0)];
    let found = find_with(
        "O[Zn]O",
        &points,
        [1],
        [0, 2],
        &InteractionCriteria::default(),
    );
    assert_eq!(found.of_kind(InteractionKind::MetalCoordination).count(), 2);
}

#[test]
fn bond_separation_within_one_molecule() {
    // Butane-1,4-diol folded so its hydroxyls are 2.8 Å apart, five bonds
    // apart through the chain. Atom order: O, C, C, C, C, O.
    let points = [
        p(0.0, 0.0, 0.0),
        p(-1.0, -1.0, 0.0),
        p(0.0, -2.4, 0.0),
        p(2.8, -2.4, 0.0),
        p(3.8, -1.0, 0.0),
        p(2.8, 0.0, 0.0),
    ];
    let bonds = |separation: usize| {
        let criteria = InteractionCriteria::default().with_min_bond_separation(separation);
        find_with("OCCCCO", &points, [0], [5], &criteria)
            .of_kind(InteractionKind::HydrogenBond)
            .count()
    };
    assert_eq!(bonds(InteractionCriteria::DEFAULT_MIN_BOND_SEPARATION), 1);
    assert_eq!(bonds(5), 1);
    assert_eq!(bonds(6), 0);
    // Across two structures, no bond joins the ends.
    let (topology, model) = model("OCCCCO", &points);
    let features = InteractionFeatures::new(Arc::clone(&topology));
    let first = select(&topology, [0]);
    let second = select(&topology, [5]);
    let criteria = InteractionCriteria::default().with_min_bond_separation(6);
    let found = find_interactions_between(
        &InteractionSide::new(&features, &model, &first).unwrap(),
        &InteractionSide::new(&features, &model, &second).unwrap(),
        &criteria,
    )
    .unwrap();
    assert_eq!(found.of_kind(InteractionKind::HydrogenBond).count(), 1);
}

#[test]
fn overlapping_sets_report_each_interaction_once() {
    // Two ethanes side by side, all atoms in both sets: each contact once,
    // the lower atom first; the diagonals are 4.09 Å long.
    let points = [
        p(0.0, 0.0, 0.0),
        p(1.5, 0.0, 0.0),
        p(0.0, 3.8, 0.0),
        p(1.5, 3.8, 0.0),
    ];
    let found = find_with(
        "CC.CC",
        &points,
        0..4,
        0..4,
        &InteractionCriteria::default(),
    );
    let pairs = found
        .of_kind(InteractionKind::HydrophobicContact)
        .map(|contact| {
            (
                contact.first().atoms()[0].index(),
                contact.second().atoms()[0].index(),
            )
        })
        .collect::<Vec<_>>();
    assert_eq!(pairs, [(0, 2), (1, 3)]);
}

#[test]
fn within_and_between_agree() {
    // Acetate and methylammonium 3.2 Å apart: one structure, or each in
    // its own.
    let acetate = [
        p(-1.5, 0.0, 0.0),
        p(0.0, 0.0, 0.0),
        p(0.6, 1.0, 0.0),
        p(0.6, -1.0, 0.0),
    ];
    let amine = [p(4.7, 0.0, 0.0), p(3.2, 0.5, 0.0)];
    let mut points = acetate.to_vec();
    points.extend(amine);
    let within = find_with(
        "CC(=O)[O-].C[NH3+]",
        &points,
        0..4,
        4..6,
        &InteractionCriteria::default(),
    );

    let (left, left_model) = model("CC(=O)[O-]", &acetate);
    let (right, right_model) = model("C[NH3+]", &amine);
    let (left_features, right_features) = (
        InteractionFeatures::new(Arc::clone(&left)),
        InteractionFeatures::new(Arc::clone(&right)),
    );
    let (left_atoms, right_atoms) = (AtomSelection::all(left), AtomSelection::all(right));
    let between = find_interactions_between(
        &InteractionSide::new(&left_features, &left_model, &left_atoms).unwrap(),
        &InteractionSide::new(&right_features, &right_model, &right_atoms).unwrap(),
        &InteractionCriteria::default(),
    )
    .unwrap();
    assert!(!within.interactions().is_empty());
    assert_eq!(within.interactions().len(), between.interactions().len());
    for (a, b) in within.interactions().iter().zip(between.interactions()) {
        assert_eq!(a.kind(), b.kind());
        assert_eq!(a.first().atoms(), b.first().atoms());
        let shifted = b
            .second()
            .atoms()
            .iter()
            .map(|atom| TopologyAtomIndex::new(atom.raw() + 4))
            .collect::<Vec<_>>();
        assert_eq!(a.second().atoms(), shifted);
        assert_eq!(a.distance(), b.distance());
        assert_eq!(
            a.angles().collect::<Vec<_>>(),
            b.angles().collect::<Vec<_>>()
        );
    }
    assert!(!between.is_within());
}

#[test]
fn disabled_types_and_invalid_criteria() {
    let points = [p(0.0, 0.0, 0.0), p(2.9, 0.0, 0.0)];
    let criteria = InteractionCriteria::default()
        .with_hydrogen_bond(HydrogenBondCriteria::default().with_enabled(false));
    assert_eq!(
        find_with("O.O", &points, [0], [1], &criteria)
            .interactions()
            .len(),
        0
    );

    let invalid = [
        (
            InteractionCriteria::default().with_salt_bridge(
                SaltBridgeCriteria::default().with_distance(Quantity::new(-1.0, ANGSTROM)),
            ),
            InteractionCriterion::SaltBridgeDistance,
        ),
        (
            InteractionCriteria::default().with_pi_stacking(
                PiStackingCriteria::default().with_parallel_angle(Quantity::new(200.0, DEGREE)),
            ),
            InteractionCriterion::PiStackingParallelAngle,
        ),
        (
            InteractionCriteria::default().with_hydrophobic(
                HydrophobicCriteria::default().with_distance(Quantity::new(4.0, DEGREE)),
            ),
            InteractionCriterion::HydrophobicDistance,
        ),
        (
            InteractionCriteria::default().with_metal_coordination(
                MetalCoordinationCriteria::default()
                    .with_distance(Quantity::new(f64::NAN, ANGSTROM)),
            ),
            InteractionCriterion::MetalCoordinationDistance,
        ),
    ];
    let (topology, model) = model("O.O", &points);
    let features = InteractionFeatures::new(Arc::clone(&topology));
    let all = AtomSelection::all(topology);
    for (criteria, criterion) in invalid {
        assert_eq!(
            criteria.validate(),
            Err(InteractionError::InvalidCriterion(criterion))
        );
        assert_eq!(
            find_interactions_within(&features, &model, &all, &all, &criteria).unwrap_err(),
            InteractionError::InvalidCriterion(criterion)
        );
    }
    // A structure of another layout is refused.
    let (_, other) = self::model("O.O", &points);
    assert_eq!(
        find_interactions_within(
            &features,
            &other,
            &all,
            &all,
            &InteractionCriteria::default()
        )
        .unwrap_err(),
        InteractionError::TopologyMismatch
    );
}

#[test]
fn hydrogen_bond_occupancy_accumulates_and_merges() {
    let (topology, _) = model("O.O", &[Point3::origin(); 2]);
    let features = InteractionFeatures::new(Arc::clone(&topology));
    let all = AtomSelection::all(Arc::clone(&topology));
    let observe = |x: f64| {
        let positions = Positions::new(Quantity::new(
            vec![p(0.0, 0.0, 0.0), p(x, 0.0, 0.0)],
            ANGSTROM,
        ))
        .unwrap();
        let model = Model::new(Arc::clone(&topology), positions).unwrap();
        find_interactions_within(
            &features,
            &model,
            &all,
            &all,
            &InteractionCriteria::default(),
        )
        .unwrap()
    };
    let observations = [(2.9, 1.0), (5.0, 2.0), (3.0, 1.0), (6.0, 4.0)];
    let mut whole = HydrogenBondOccupancyAccumulator::new();
    for (index, (x, weight)) in observations.iter().enumerate() {
        whole.observe(index, &observe(*x), *weight).unwrap();
    }
    let mut left = HydrogenBondOccupancyAccumulator::new();
    let mut right = HydrogenBondOccupancyAccumulator::new();
    for (index, (x, weight)) in observations.iter().enumerate() {
        let target = if index < 2 { &mut left } else { &mut right };
        target.observe(index, &observe(*x), *weight).unwrap();
    }
    left.merge(&right).unwrap();
    let (whole, merged) = (whole.finish().unwrap(), left.finish().unwrap());
    assert_eq!(whole.keys(), merged.keys());
    assert_eq!(whole.weights(), merged.weights());
    assert_eq!(merged.observation_count(), 4);
    let bonds = merged.bonds().collect::<Vec<_>>();
    assert_eq!(bonds.len(), 1);
    assert_eq!(bonds[0].0.donor().atom(), TopologyAtomIndex::new(0));
    assert!(bonds[0].0.donor().is_first() && bonds[0].0.acceptor().is_first());
    assert_eq!(bonds[0].1, 0.25);

    // Errors change nothing.
    let mut accumulator = HydrogenBondOccupancyAccumulator::new();
    let found = observe(2.9);
    assert_eq!(
        accumulator.observe(7, &found, 0.0),
        Err(InteractionError::InvalidWeight { index: 7 })
    );
    assert_eq!(
        accumulator.observe(8, &found, f64::INFINITY),
        Err(InteractionError::InvalidWeight { index: 8 })
    );
    accumulator.observe(0, &found, f64::MAX).unwrap();
    assert_eq!(
        accumulator.observe(1, &found, f64::MAX),
        Err(InteractionError::WeightOverflow { index: 1 })
    );
    assert_eq!(accumulator.observation_count(), 1);
    let other_criteria = InteractionCriteria::default().with_hydrogen_bond(
        HydrogenBondCriteria::default().with_distance(Quantity::new(3.6, ANGSTROM)),
    );
    let positions = Positions::new(Quantity::new(
        vec![p(0.0, 0.0, 0.0), p(2.9, 0.0, 0.0)],
        ANGSTROM,
    ))
    .unwrap();
    let model = Model::new(Arc::clone(&topology), positions).unwrap();
    let other = find_interactions_within(&features, &model, &all, &all, &other_criteria).unwrap();
    let mut fresh = HydrogenBondOccupancyAccumulator::new();
    fresh.observe(0, &found, 1.0).unwrap();
    assert_eq!(
        fresh.observe(1, &other, 1.0),
        Err(InteractionError::ObservationMismatch { index: 1 })
    );
    let mut mismatched = HydrogenBondOccupancyAccumulator::new();
    mismatched.observe(0, &other, 1.0).unwrap();
    assert_eq!(
        fresh.merge(&mismatched),
        Err(InteractionError::MergeMismatch)
    );
    assert_eq!(
        HydrogenBondOccupancyAccumulator::new()
            .finish()
            .unwrap_err(),
        InteractionError::NoObservations
    );
}
