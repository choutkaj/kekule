//! Informational cost of clustering a set of geometries in memory: the pairwise
//! fitted RMSD matrix, hierarchical clustering on it, and leader clustering.
//! Run: cargo bench -p kekule --bench clustering_fits -- SYSTEM.cif [SAMPLES] [THREADS]
//! The geometries are an external mmCIF model's coordinates, each turned,
//! moved, and perturbed by a fixed pseudo-random amount, and fitted on the
//! model's alpha carbons; this measures fitting cost, not sampling.
//! The checksum line changes when any result bit changes.

use std::error::Error;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::num::NonZeroUsize;
use std::time::Instant;

use kekule::clustering::{ClusteringOptions, ClusteringRun, Cut, Linkage};
use kekule::math::Point3;
use kekule::mmcif;
use kekule::parallel::ScopedThreads;
use kekule::structure::{Geometry, Positions, StructureSet};
use kekule::topology::AtomSelection;
use kekule::units::{Quantity, ANGSTROM};

/// Geometries when none are requested.
const SAMPLES: usize = 400;
/// Runs of each step; the fastest is reported.
const RUNS: usize = 3;

fn main() -> Result<(), Box<dyn Error>> {
    // Cargo passes `--bench` after the arguments.
    let mut arguments = std::env::args()
        .skip(1)
        .filter(|argument| !argument.starts_with("--"));
    let model_path = arguments.next().ok_or("supply an external mmCIF model")?;
    let geometries = arguments
        .next()
        .map_or(Ok(SAMPLES), |count| count.parse())?;
    let threads: NonZeroUsize = match arguments.next() {
        Some(threads) => threads.parse()?,
        None => std::thread::available_parallelism()?,
    };
    let model = mmcif::parse_str(&std::fs::read_to_string(model_path)?, Default::default())?
        .interpret(Default::default())?
        .into_model();
    let topology = model.shared_topology();
    let atoms = AtomSelection::for_label_atom_names(topology.clone(), ["CA"])?;
    let points = model
        .positions()
        .clone()
        .into_values()
        .into_unit(ANGSTROM)?
        .into_value();
    let mut set = StructureSet::new(topology.clone());
    let mut random = Random(0x2545_f491_4f6c_dd1d);
    for index in 0..geometries {
        // Four families of conformations, each bending a quarter of the
        // atoms, so that clusters form at 1 Å.
        let family = index % 4;
        let (angle, offset) = (random.next() * 6.0, random.next() * 10.0);
        let (sin, cos) = angle.sin_cos();
        let moved = points
            .iter()
            .enumerate()
            .map(|(atom, point)| {
                let bend = if atom % 4 == family { 2.5 } else { 0.0 };
                let x = point.x + bend + 0.4 * random.next();
                let y = point.y + 0.4 * random.next();
                let z = point.z + 0.4 * random.next();
                Point3::new(cos * x - sin * y + offset, sin * x + cos * y, z - offset)
            })
            .collect::<Vec<_>>();
        set.push(
            Geometry::new(Positions::new(Quantity::new(moved, ANGSTROM))?),
            None,
        )?;
    }

    let executor = ScopedThreads::new(threads);
    let run = || ClusteringRun::default().with_executor(&executor);
    let radius = Quantity::new(1.0, ANGSTROM);
    let mut checksum = DefaultHasher::new();

    let mut matrix_ms = f64::INFINITY;
    let mut matrix = None;
    for _ in 0..RUNS {
        let start = Instant::now();
        let distances = set.pairwise_fitted_rmsd_with(&atoms, Default::default(), run())?;
        matrix_ms = matrix_ms.min(start.elapsed().as_secs_f64() * 1e3);
        matrix = Some(distances);
    }
    let matrix = matrix.ok_or("no runs")?;
    for value in *matrix.condensed().value() {
        value.to_bits().hash(&mut checksum);
    }

    let mut hierarchical_ms = f64::INFINITY;
    let mut hierarchical_states = 0;
    for _ in 0..RUNS {
        let start = Instant::now();
        let distances = set.pairwise_fitted_rmsd_with(&atoms, Default::default(), run())?;
        let partition =
            set.cluster_hierarchical(&distances, Linkage::Average, Cut::Height(radius))?;
        hierarchical_ms = hierarchical_ms.min(start.elapsed().as_secs_f64() * 1e3);
        hierarchical_states = partition.state_count();
        partition.assignments().hash(&mut checksum);
        partition.representatives().hash(&mut checksum);
    }

    let mut leader_ms = f64::INFINITY;
    let mut leader_states = 0;
    for _ in 0..RUNS {
        let start = Instant::now();
        let partition =
            set.cluster_by_radius_with(radius, &atoms, ClusteringOptions::default(), run())?;
        leader_ms = leader_ms.min(start.elapsed().as_secs_f64() * 1e3);
        leader_states = partition.state_count();
        partition.assignments().hash(&mut checksum);
        partition.representatives().hash(&mut checksum);
    }

    println!(
        "{{\"geometries\":{geometries},\"fitted_atoms\":{},\"threads\":{threads},\"matrix_ms\":{matrix_ms:.1},\"hierarchical_ms\":{hierarchical_ms:.1},\"hierarchical_states\":{hierarchical_states},\"leader_ms\":{leader_ms:.1},\"leader_states\":{leader_states},\"checksum\":\"{:016x}\"}}",
        atoms.len(),
        checksum.finish()
    );
    Ok(())
}

/// A fixed sequence of numbers in [-0.5, 0.5).
struct Random(u64);

impl Random {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1_u64 << 53) as f64 - 0.5
    }
}
