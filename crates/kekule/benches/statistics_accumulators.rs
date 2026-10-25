//! Informational cost of weighted RMSF and contact occupancy over frames of
//! 100,000 atoms: one observation, a read of the same coordinates, and
//! partial accumulators observed in parallel and merged in order.

use std::error::Error;
use std::hint::black_box;
use std::num::NonZeroUsize;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use kekule::geometry::Point3;
use kekule::parallel::{Executor, ScopedThreads};
use kekule::smiles;
use kekule::statistics::{ContactOccupancyAccumulator, RmsfAccumulator};
use kekule::structure::{Model, Positions};
use kekule::topology::{AtomSelection, Topology, TopologyBuilder};
use kekule::units::{Quantity, NANOMETER};

const ATOMS: usize = 100_000;
const FRAMES: usize = 96;
const PASSES: usize = 11;

fn main() -> Result<(), Box<dyn Error>> {
    println!(
        "informational only; target={}-{}; atoms={ATOMS}; frames={FRAMES}; median of {PASSES} passes",
        std::env::consts::OS,
        std::env::consts::ARCH,
    );
    let topology = topology()?;
    let frames = frames(&topology)?;
    let selection = AtomSelection::all(Arc::clone(&topology));
    let atoms = topology.atom_ids();
    let pairs = (0..ATOMS / 2)
        .map(|pair| (atoms[2 * pair], atoms[2 * pair + 1]))
        .collect::<Vec<_>>();
    let cutoff = Quantity::new(0.5, NANOMETER);

    let read = per_frame(|| {
        let started = Instant::now();
        for frame in &frames {
            let mut sum = 0.0;
            for point in frame.positions().values().value().iter() {
                sum += point.x + point.y + point.z;
            }
            black_box(sum);
        }
        started.elapsed()
    });
    let rmsf = per_frame(|| {
        let mut accumulator = RmsfAccumulator::new(&selection).expect("nonempty");
        let started = Instant::now();
        for (index, frame) in frames.iter().enumerate() {
            accumulator.observe(index, frame, 1.0).expect("finite");
        }
        let elapsed = started.elapsed();
        black_box(accumulator.finish().expect("observed"));
        elapsed
    });
    let contacts = per_frame(|| {
        let mut accumulator =
            ContactOccupancyAccumulator::new(&topology, pairs.iter().copied(), cutoff)
                .expect("valid pairs");
        let started = Instant::now();
        for (index, frame) in frames.iter().enumerate() {
            accumulator.observe(index, frame, 1.0).expect("finite");
        }
        let elapsed = started.elapsed();
        black_box(accumulator.finish().expect("observed"));
        elapsed
    });
    println!("| per frame | ms |");
    println!("|---|---:|");
    println!("| read coordinates | {} |", ms(read));
    println!("| RMSF observe, {ATOMS} atoms | {} |", ms(rmsf));
    println!(
        "| contact observe, {} pairs | {} |",
        pairs.len(),
        ms(contacts)
    );

    let merge_rmsf = median(|| {
        let mut left = RmsfAccumulator::new(&selection).expect("nonempty");
        let mut right = RmsfAccumulator::new(&selection).expect("nonempty");
        left.observe(0, &frames[0], 1.0).expect("finite");
        right.observe(1, &frames[1], 1.0).expect("finite");
        let started = Instant::now();
        left.merge(&right).expect("compatible");
        let elapsed = started.elapsed();
        black_box(left);
        elapsed
    });
    println!("| RMSF merge of two partials | {} |", ms(merge_rmsf));

    println!();
    println!("| threads | RMSF frames/s | contact frames/s |");
    println!("|---:|---:|---:|");
    for threads in [1, 2, 4, 6] {
        let executor = ScopedThreads::new(NonZeroUsize::new(threads).expect("nonzero"));
        let parallel_rmsf = median(|| {
            let started = Instant::now();
            let merged = in_parallel(
                &executor,
                &frames,
                || RmsfAccumulator::new(&selection).expect("nonempty"),
                |accumulator, index, frame| accumulator.observe(index, frame, 1.0).expect("finite"),
                |left, right| left.merge(&right).expect("compatible"),
            );
            let elapsed = started.elapsed();
            black_box(merged.finish().expect("observed"));
            elapsed
        });
        let parallel_contacts = median(|| {
            let started = Instant::now();
            let merged = in_parallel(
                &executor,
                &frames,
                || {
                    ContactOccupancyAccumulator::new(&topology, pairs.iter().copied(), cutoff)
                        .expect("valid pairs")
                },
                |accumulator, index, frame| accumulator.observe(index, frame, 1.0).expect("finite"),
                |left, right| left.merge(&right).expect("compatible"),
            );
            let elapsed = started.elapsed();
            black_box(merged.finish().expect("observed"));
            elapsed
        });
        println!(
            "| {threads} | {:.0} | {:.0} |",
            FRAMES as f64 / parallel_rmsf.as_secs_f64(),
            FRAMES as f64 / parallel_contacts.as_secs_f64(),
        );
    }
    Ok(())
}

/// Observes contiguous ranges of `frames`, one per thread, and merges the
/// partial accumulators in range order.
fn in_parallel<A: Send>(
    executor: &dyn Executor,
    frames: &[Model],
    empty: impl Fn() -> A + Sync,
    observe: impl Fn(&mut A, usize, &Model) + Sync,
    merge: impl Fn(&mut A, A),
) -> A {
    let chunks = executor.parallelism().get();
    let partials = (0..chunks).map(|_| Mutex::new(None)).collect::<Vec<_>>();
    executor.run(chunks, &|chunk| {
        let mut accumulator = empty();
        let range = chunk * frames.len() / chunks..(chunk + 1) * frames.len() / chunks;
        for index in range {
            observe(&mut accumulator, index, &frames[index]);
        }
        *partials[chunk].lock().expect("unpoisoned") = Some(accumulator);
    });
    let mut partials = partials
        .into_iter()
        .map(|partial| partial.into_inner().expect("unpoisoned").expect("ran"));
    let mut merged = partials.next().expect("one chunk");
    for partial in partials {
        merge(&mut merged, partial);
    }
    merged
}

fn topology() -> Result<Arc<Topology>, Box<dyn Error>> {
    let molecule = smiles::to_molecules("C")?.pop().expect("one molecule");
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(molecule)?;
    for _ in 0..ATOMS {
        builder.add_instance(definition)?;
    }
    Ok(Arc::new(builder.build()?))
}

/// Frames of a cubic lattice with 0.1 nm spacing, each atom displaced by a
/// deterministic pseudo-random offset of up to 0.05 nm per axis.
fn frames(topology: &Arc<Topology>) -> Result<Vec<Model>, Box<dyn Error>> {
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut jitter = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64 * 0.1 - 0.05
    };
    let side = (ATOMS as f64).cbrt().ceil() as usize;
    (0..FRAMES)
        .map(|_| {
            let points = (0..ATOMS)
                .map(|atom| {
                    let (x, y, z) = (atom % side, atom / side % side, atom / side / side);
                    Point3::new(
                        x as f64 * 0.1 + jitter(),
                        y as f64 * 0.1 + jitter(),
                        z as f64 * 0.1 + jitter(),
                    )
                })
                .collect::<Vec<_>>();
            let positions = Positions::new(Quantity::new(points, NANOMETER))?;
            Ok(Model::new(Arc::clone(topology), positions)?)
        })
        .collect()
}

fn per_frame(pass: impl FnMut() -> Duration) -> Duration {
    median(pass) / FRAMES as u32
}

fn median(mut pass: impl FnMut() -> Duration) -> Duration {
    let mut times = (0..PASSES).map(|_| pass()).collect::<Vec<_>>();
    times.sort_unstable();
    times[PASSES / 2]
}

fn ms(duration: Duration) -> String {
    format!("{:.3}", duration.as_secs_f64() * 1e3)
}
