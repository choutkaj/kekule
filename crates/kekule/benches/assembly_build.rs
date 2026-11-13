//! Informational cost of building a biological assembly from an interpreted
//! mmCIF entry: the time `Assembly::to_structure` takes for each of the entry's
//! assemblies, and the atoms it places.
//!
//! The entry is externally supplied, such as PDB 4RHV (a 60-copy virus
//! capsid of 392,520 atoms) or 4V60 (a vault, 2 copies of 483,912 atoms):
//!
//! ```sh
//! curl -sSfLO https://files.rcsb.org/download/4RHV.cif
//! cargo bench -p kekule --bench assembly_build -- 4RHV.cif
//! ```
//!
//! The argument is the mmCIF path; the `KEKULE_ASSEMBLY_CIF` environment
//! variable may supply it instead. Without a path the benchmark reports
//! nothing.

use std::error::Error;
use std::hint::black_box;
use std::time::{Duration, Instant};

use kekule::mmcif;

const ROUNDS: usize = 5;

fn main() -> Result<(), Box<dyn Error>> {
    let Some(path) = std::env::args()
        .skip(1)
        .find(|argument| !argument.starts_with("--"))
        .or_else(|| std::env::var("KEKULE_ASSEMBLY_CIF").ok())
    else {
        println!("assembly_build: no mmCIF path given; nothing measured");
        return Ok(());
    };
    let text = std::fs::read_to_string(&path)?;
    let interpretation =
        mmcif::parse_str(&text, Default::default())?.interpret(Default::default())?;
    let structure = interpretation.structure();
    println!(
        "informational only; target={}-{}; deposited atoms={}; minimum and median of {ROUNDS} \
         rounds",
        std::env::consts::OS,
        std::env::consts::ARCH,
        structure.atom_count(),
    );
    println!("| assembly | copies | atoms | build min ms | build median ms |");
    println!("|---|---:|---:|---:|---:|");
    let milliseconds = |time: Duration| time.as_secs_f64() * 1e3;
    for assembly in interpretation.report().assemblies() {
        let copies: usize = assembly
            .generators()
            .iter()
            .map(|generator| generator.operations().len())
            .sum();
        let mut times = Vec::new();
        let mut atoms = 0;
        for _ in 0..ROUNDS {
            let start = Instant::now();
            let built = assembly.to_structure(structure, None)?;
            times.push(start.elapsed());
            atoms = built.atom_count();
            black_box(built);
        }
        times.sort();
        println!(
            "| {} | {copies} | {atoms} | {:.1} | {:.1} |",
            assembly.id(),
            milliseconds(times[0]),
            milliseconds(times[ROUNDS / 2]),
        );
    }
    Ok(())
}
