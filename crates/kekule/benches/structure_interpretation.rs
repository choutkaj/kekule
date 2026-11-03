//! Informational cost of reading a large mmCIF entry: parsing, interpretation
//! into a model, and default perception of that model.
//!
//! The entry is externally supplied. SeekDelta's request K24 used PDB 6SPB
//! (90,351 atoms) and 7K00 (149,338 atoms, 7,640 molecules):
//!
//! ```sh
//! curl -sSfLO https://files.rcsb.org/download/6SPB.cif
//! cargo bench -p kekule --bench structure_interpretation -- 6SPB.cif
//! ```
//!
//! The argument is the mmCIF path; the `KEKULE_INTERPRETATION_CIF`
//! environment variable may supply it instead. Without a path the benchmark
//! reports nothing. A perception failure is reported with its time, since a
//! failing definition ends perception early.

use std::error::Error;
use std::hint::black_box;
use std::time::{Duration, Instant};

use kekule::mmcif;

const ROUNDS: usize = 7;

fn main() -> Result<(), Box<dyn Error>> {
    let Some(path) = std::env::args()
        .skip(1)
        .find(|argument| !argument.starts_with("--"))
        .or_else(|| std::env::var("KEKULE_INTERPRETATION_CIF").ok())
    else {
        println!("structure_interpretation: no mmCIF path given; nothing measured");
        return Ok(());
    };
    let text = std::fs::read_to_string(&path)?;

    let mut phases = [const { Vec::new() }; 3];
    let mut outcome = None;
    for _ in 0..ROUNDS {
        let start = Instant::now();
        let document = mmcif::parse_str(&text, Default::default())?;
        phases[0].push(start.elapsed());
        let start = Instant::now();
        let interpretation = document.interpret(Default::default())?;
        phases[1].push(start.elapsed());
        let mut model = interpretation.into_model();
        let start = Instant::now();
        let perceived = model.perceive().map(|_| ());
        phases[2].push(start.elapsed());
        outcome = Some((model, perceived));
        black_box(document);
    }
    let (model, perceived) = outcome.expect("at least one round");
    let topology = model.topology();
    println!(
        "informational only; target={}-{}; atoms={}; molecules={}; definitions={}; \
         minimum and median of {ROUNDS} rounds",
        std::env::consts::OS,
        std::env::consts::ARCH,
        model.atom_count(),
        topology.instance_count(),
        topology.definition_count(),
    );
    if let Err(error) = &perceived {
        println!("perception failed: {error}");
    }
    println!("| phase | min ms | median ms |");
    println!("|---|---:|---:|");
    for (name, times) in ["parse", "interpret", "perceive"].iter().zip(&mut phases) {
        times.sort();
        let milliseconds = |time: Duration| time.as_secs_f64() * 1e3;
        println!(
            "| {name} | {:.1} | {:.1} |",
            milliseconds(times[0]),
            milliseconds(times[times.len() / 2]),
        );
    }
    Ok(())
}
