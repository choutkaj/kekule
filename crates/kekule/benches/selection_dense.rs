//! Informational cost of atom selections built from, and turned into, dense
//! indices and bit masks over a million atoms.

use std::error::Error;
use std::hint::black_box;
use std::sync::Arc;
use std::time::{Duration, Instant};

use kekule::smiles;
use kekule::topology::{AtomSelection, Topology, TopologyAtomIndex, TopologyBuilder};

const ATOMS: usize = 1_000_000;
const PASSES: usize = 21;

fn main() -> Result<(), Box<dyn Error>> {
    println!(
        "informational only; target={}-{}; atoms={ATOMS}; median of {PASSES} passes",
        std::env::consts::OS,
        std::env::consts::ARCH,
    );
    let topology = topology()?;
    println!("| selected | via atom IDs ms | from_indices ms | from_sorted_indices ms | from_mask ms | output floor ms | to_mask ms | ranges ms | runs |");
    println!("|---|---:|---:|---:|---:|---:|---:|---:|---:|");
    for (name, keep) in [
        ("1 %", (|index| index % 100 == 0) as fn(u32) -> bool),
        ("every other", |index| index % 2 == 0),
        ("half, contiguous", |index| (index as usize) < ATOMS / 2),
        ("all", |_| true),
    ] {
        let indices = (0..ATOMS as u32)
            .filter(|index| keep(*index))
            .map(TopologyAtomIndex::new)
            .collect::<Vec<_>>();
        let expected = AtomSelection::from_sorted_indices(Arc::clone(&topology), indices.clone())?;
        let mask = expected.to_mask();
        // The ID round trip `from_indices` made before dense constructors.
        let via_ids = median(|| {
            let started = Instant::now();
            let atoms = indices
                .iter()
                .map(|index| topology.atom_id(*index).expect("valid index"));
            black_box(AtomSelection::from_atoms(Arc::clone(&topology), atoms).expect("valid"));
            started.elapsed()
        });
        let from_indices = median(|| {
            let started = Instant::now();
            black_box(
                AtomSelection::from_indices(Arc::clone(&topology), indices.iter().copied())
                    .expect("valid"),
            );
            started.elapsed()
        });
        let from_sorted = median(|| {
            let owned = indices.clone();
            let started = Instant::now();
            black_box(
                AtomSelection::from_sorted_indices(Arc::clone(&topology), owned).expect("valid"),
            );
            started.elapsed()
        });
        let from_mask = median(|| {
            let started = Instant::now();
            let selection = AtomSelection::from_mask(Arc::clone(&topology), &mask).expect("valid");
            black_box(&selection);
            started.elapsed()
        });
        // Allocating and writing the selection's indices alone, the floor of
        // any constructor that does not reuse its caller's vector.
        let floor = median(|| {
            let started = Instant::now();
            black_box(
                (0..indices.len() as u32)
                    .map(TopologyAtomIndex::new)
                    .collect::<Vec<_>>(),
            );
            started.elapsed()
        });
        assert_eq!(
            AtomSelection::from_mask(Arc::clone(&topology), &mask)?,
            expected
        );
        let to_mask = median(|| {
            let started = Instant::now();
            black_box(expected.to_mask());
            started.elapsed()
        });
        let mut runs = 0;
        let ranges = median(|| {
            let started = Instant::now();
            runs = black_box(expected.ranges().count());
            started.elapsed()
        });
        println!(
            "| {name} | {} | {} | {} | {} | {} | {} | {} | {runs} |",
            ms(via_ids),
            ms(from_indices),
            ms(from_sorted),
            ms(from_mask),
            ms(floor),
            ms(to_mask),
            ms(ranges),
        );
    }
    Ok(())
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

fn median(mut pass: impl FnMut() -> Duration) -> Duration {
    let mut times = (0..PASSES).map(|_| pass()).collect::<Vec<_>>();
    times.sort_unstable();
    times[PASSES / 2]
}

fn ms(duration: Duration) -> String {
    format!("{:.3}", duration.as_secs_f64() * 1e3)
}
