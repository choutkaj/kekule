//! Informational cost of one observation of group contact occupancy in a
//! large structure: a pocket's residues against a ligand, and the residues on
//! both sides of an interface, beside the atom pairs of the pocket and ligand.

use std::error::Error;
use std::hint::black_box;
use std::sync::Arc;
use std::time::{Duration, Instant};

use kekule::geometry::Point3;
use kekule::smiles;
use kekule::statistics::{ContactOccupancyAccumulator, GroupContactOccupancyAccumulator};
use kekule::structure::{Positions, Trajectory, TrajectoryFrame};
use kekule::topology::{AtomSelection, Topology, TopologyAtomIndex, TopologyBuilder};
use kekule::units::{Quantity, NANOMETER};

/// Atoms on a cubic lattice of this many points a side.
const SIDE: usize = 47;
/// About 0.09 atoms per cubic ångström, as in proteins.
const SPACING: f64 = 0.22;
/// Consecutive atoms along z make a residue.
const RESIDUE: usize = 10;
const CUTOFF: f64 = 0.45;
const PASSES: usize = 21;

fn main() -> Result<(), Box<dyn Error>> {
    let atoms = SIDE * SIDE * SIDE;
    let topology = topology(atoms)?;
    let trajectory = frames(&topology, atoms)?;
    let residue_of = |index: usize| index / RESIDUE;
    let residues = atoms / RESIDUE;
    let residue = |r: usize| {
        AtomSelection::from_sorted_indices(
            Arc::clone(&topology),
            (r * RESIDUE..(r + 1) * RESIDUE)
                .map(|index| TopologyAtomIndex::new(index as u32))
                .collect(),
        )
    };
    let point = |index: usize| {
        let (x, rest) = (index / (SIDE * SIDE), index % (SIDE * SIDE));
        Point3::new(
            x as f64 * SPACING,
            (rest / SIDE) as f64 * SPACING,
            (rest % SIDE) as f64 * SPACING,
        )
    };
    let center = point(atoms / 2);
    let distance = |a: Point3, b: Point3| (a - b).norm();

    // A ligand of four residues at the center, and the pocket: the residues
    // with an atom within 1 nm of the ligand's first atom.
    let ligand_residues = (residue_of(atoms / 2)..residue_of(atoms / 2) + 4).collect::<Vec<_>>();
    let ligand = ligand_residues
        .iter()
        .map(|&r| residue(r))
        .try_fold(AtomSelection::empty(Arc::clone(&topology)), |all, r| {
            all.union(&r?).map_err(Box::<dyn Error>::from)
        })?;
    let pocket = (0..residues)
        .filter(|r| !ligand_residues.contains(r))
        .filter(|r| (r * RESIDUE..(r + 1) * RESIDUE).any(|i| distance(point(i), center) <= 1.0))
        .collect::<Vec<_>>();
    // The interface: residues within 0.5 nm of the plane x = center on
    // either side, with their first atom within 2 nm of the center line.
    let near_line = |r: usize| {
        let p = point(r * RESIDUE);
        (p.y - center.y).hypot(p.z - center.z) <= 2.0
    };
    let side = |low: f64, high: f64| {
        (0..residues)
            .filter(|&r| {
                let x = point(r * RESIDUE).x - center.x;
                x >= low && x < high && near_line(r)
            })
            .collect::<Vec<_>>()
    };
    let (left, right) = (side(-0.5, 0.0), side(0.0, 0.5));

    println!(
        "informational only; target={}-{}; atoms={atoms}; cutoff={CUTOFF} nm; median of {PASSES} passes",
        std::env::consts::OS,
        std::env::consts::ARCH,
    );
    println!("| case | groups | pairs | grouped atoms | µs per observation |");
    println!("|---|---:|---:|---:|---:|");
    let cutoff = Quantity::new(CUTOFF, NANOMETER);

    let mut groups = vec![ligand.clone()];
    groups.extend(
        pocket
            .iter()
            .map(|&r| residue(r))
            .collect::<Result<Vec<_>, _>>()?,
    );
    let pairs = (1..groups.len()).map(|g| (g, 0)).collect::<Vec<_>>();
    let grouped = groups.iter().map(AtomSelection::len).sum::<usize>();
    let mut accumulator =
        GroupContactOccupancyAccumulator::new(&topology, groups.clone(), pairs.clone(), cutoff)?;
    let time = observe(&trajectory, |index, frame| {
        accumulator
            .observe(index, frame, 1.0)
            .expect("valid observation")
    });
    row(
        "pocket residues × ligand",
        groups.len(),
        pairs.len(),
        grouped,
        time,
    );

    let ids = topology.atom_ids();
    let atom_pairs = groups[1..]
        .iter()
        .flat_map(|r| r.indices().iter())
        .flat_map(|a| ligand.indices().iter().map(move |b| (*a, *b)))
        .map(|(a, b)| (ids[a.index()], ids[b.index()]))
        .collect::<Vec<_>>();
    let count = atom_pairs.len();
    let mut atom_accumulator = ContactOccupancyAccumulator::new(&topology, atom_pairs, cutoff)?;
    let time = observe(&trajectory, |index, frame| {
        atom_accumulator
            .observe(index, frame, 1.0)
            .expect("valid observation")
    });
    row(
        "their atom pairs (atom-pair accumulator)",
        0,
        count,
        grouped,
        time,
    );

    let mut groups = left
        .iter()
        .map(|&r| residue(r))
        .collect::<Result<Vec<_>, _>>()?;
    groups.extend(
        right
            .iter()
            .map(|&r| residue(r))
            .collect::<Result<Vec<_>, _>>()?,
    );
    let (lefts, rights) = (left.len(), right.len());
    let pairs = (0..lefts)
        .flat_map(|a| (0..rights).map(move |b| (a, lefts + b)))
        .collect::<Vec<_>>();
    let grouped = groups.iter().map(AtomSelection::len).sum::<usize>();
    let mut accumulator =
        GroupContactOccupancyAccumulator::new(&topology, groups.clone(), pairs.clone(), cutoff)?;
    let time = observe(&trajectory, |index, frame| {
        accumulator
            .observe(index, frame, 1.0)
            .expect("valid observation")
    });
    row(
        "interface residues × residues",
        groups.len(),
        pairs.len(),
        grouped,
        time,
    );

    let sides = [
        left.iter()
            .map(|&r| residue(r))
            .try_fold(AtomSelection::empty(Arc::clone(&topology)), |all, r| {
                all.union(&r?)
            })?,
        right
            .iter()
            .map(|&r| residue(r))
            .try_fold(AtomSelection::empty(Arc::clone(&topology)), |all, r| {
                all.union(&r?)
            })?,
    ];
    let grouped = sides.iter().map(AtomSelection::len).sum::<usize>();
    let mut accumulator =
        GroupContactOccupancyAccumulator::new(&topology, sides.clone(), [(0, 1)], cutoff)?;
    let time = observe(&trajectory, |index, frame| {
        accumulator
            .observe(index, frame, 1.0)
            .expect("valid observation")
    });
    row("interface side × side", 2, 1, grouped, time);
    Ok(())
}

fn topology(atoms: usize) -> Result<Arc<Topology>, Box<dyn Error>> {
    let molecule = smiles::to_molecules("C")?.pop().expect("one molecule");
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(molecule)?;
    for _ in 0..atoms {
        builder.add_instance(definition)?;
    }
    Ok(Arc::new(builder.build()?))
}

/// Two frames of the lattice, each atom displaced by up to a third of the
/// spacing along every axis.
fn frames(topology: &Arc<Topology>, atoms: usize) -> Result<Trajectory, Box<dyn Error>> {
    let frame = |seed: u64| {
        let mut state = seed;
        let mut jitter = move || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            ((state >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * SPACING * (2.0 / 3.0)
        };
        let points = (0..atoms)
            .map(|index| {
                let (x, rest) = (index / (SIDE * SIDE), index % (SIDE * SIDE));
                Point3::new(
                    x as f64 * SPACING + jitter(),
                    (rest / SIDE) as f64 * SPACING + jitter(),
                    (rest % SIDE) as f64 * SPACING + jitter(),
                )
            })
            .collect::<Vec<_>>();
        Positions::new(Quantity::new(points, NANOMETER)).map(TrajectoryFrame::new)
    };
    Ok(Trajectory::from_items(
        Arc::clone(topology),
        [frame(1)?, frame(2)?],
    )?)
}

/// The median time of one observation, alternating the frames.
fn observe(
    trajectory: &Trajectory,
    mut observe: impl FnMut(usize, &kekule::structure::TrajectoryFrameView<'_>),
) -> Duration {
    let mut times = (0..PASSES)
        .map(|pass| {
            let frame = trajectory.get(pass % 2).expect("two frames");
            let started = Instant::now();
            observe(pass, black_box(&frame));
            started.elapsed()
        })
        .collect::<Vec<_>>();
    times.sort_unstable();
    times[PASSES / 2]
}

fn row(case: &str, groups: usize, pairs: usize, grouped: usize, time: Duration) {
    println!(
        "| {case} | {groups} | {pairs} | {grouped} | {:.1} |",
        time.as_secs_f64() * 1e6
    );
}
