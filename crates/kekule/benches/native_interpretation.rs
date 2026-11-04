//! Informational cost of native persistence of an interpreted mmCIF entry:
//! the size of the native encoding and the time to write and read it, of the
//! model alone and of the model with its interpretation report.
//!
//! The entry is externally supplied, such as PDB 6SPB (90,351 atoms) or 7K00
//! (149,338 atoms):
//!
//! ```sh
//! curl -sSfLO https://files.rcsb.org/download/6SPB.cif
//! cargo bench -p kekule --bench native_interpretation -- 6SPB.cif
//! ```
//!
//! The argument is the mmCIF path; the `KEKULE_INTERPRETATION_CIF`
//! environment variable may supply it instead. Without a path the benchmark
//! reports nothing. Encoding and decoding work in memory.

use std::error::Error;
use std::hint::black_box;
use std::time::{Duration, Instant};

use kekule::mmcif;
use kekule::native::{self, NativeReadOptions, NativeWriteOptions};

const ROUNDS: usize = 7;

fn main() -> Result<(), Box<dyn Error>> {
    let Some(path) = std::env::args()
        .skip(1)
        .find(|argument| !argument.starts_with("--"))
        .or_else(|| std::env::var("KEKULE_INTERPRETATION_CIF").ok())
    else {
        println!("native_interpretation: no mmCIF path given; nothing measured");
        return Ok(());
    };
    let text = std::fs::read_to_string(&path)?;
    let interpretation =
        mmcif::parse_str(&text, Default::default())?.interpret(Default::default())?;
    let report = interpretation.report();
    println!(
        "informational only; target={}-{}; atoms={}; molecules={}; report: {} instances, \
         {} coordination records, {} alternate-location decisions, {} issues, \
         {} polymer entities; minimum and median of {ROUNDS} rounds",
        std::env::consts::OS,
        std::env::consts::ARCH,
        interpretation.model().atom_count(),
        interpretation.topology().instance_count(),
        report.instances().len(),
        report.coordination().len(),
        report.alternate_locations().len(),
        report.issues().len(),
        report.polymer_entities().len(),
    );
    println!("| object | bytes | write min ms | write median ms | read min ms | read median ms |");
    println!("|---|---:|---:|---:|---:|---:|");
    let options = NativeWriteOptions::default();
    let milliseconds = |time: Duration| time.as_secs_f64() * 1e3;
    for (name, with_report) in [("model", false), ("model with report", true)] {
        let (mut writes, mut reads) = (Vec::new(), Vec::new());
        let mut size = 0;
        for _ in 0..ROUNDS {
            let start = Instant::now();
            let bytes = if with_report {
                native::write(&interpretation, options)?
            } else {
                native::write(interpretation.model(), options)?
            };
            writes.push(start.elapsed());
            size = bytes.len();
            let start = Instant::now();
            let object = native::read(bytes.as_slice(), NativeReadOptions::default())?;
            reads.push(start.elapsed());
            black_box(object);
        }
        writes.sort();
        reads.sort();
        println!(
            "| {name} | {size} | {:.1} | {:.1} | {:.1} | {:.1} |",
            milliseconds(writes[0]),
            milliseconds(writes[ROUNDS / 2]),
            milliseconds(reads[0]),
            milliseconds(reads[ROUNDS / 2]),
        );
    }
    Ok(())
}
