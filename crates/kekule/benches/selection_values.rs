//! Informational cost of selections with long value lists over a million
//! atoms, as a viewer's rectangle selection writes them: one object's host
//! term and thousands of index runs.

use std::error::Error;
use std::hint::black_box;
use std::sync::Arc;
use std::time::{Duration, Instant};

use kekule::parallel::Run;
use kekule::select::{
    self, HostRequest, HostTerm, HostTermError, HostValue, SelectEvaluateOptions,
    SelectParseOptions, SelectionResolver,
};
use kekule::smiles;
use kekule::topology::{AtomSelection, Topology, TopologyBuilder};

const ATOMS: usize = 1_000_000;
const RUNS: usize = 5_000;
/// At most this many passes, and at least three.
const PASSES: usize = 11;
/// Passes stop once they have taken this long.
const BUDGET: Duration = Duration::from_secs(2);

/// `shown` is every atom, `hidden` none, as a host's object terms are on
/// the object they name and on every other one.
struct Objects;

impl SelectionResolver for Objects {
    fn resolve(
        &self,
        term: &HostTerm<'_>,
        request: HostRequest<'_>,
    ) -> Result<HostValue, HostTermError> {
        let topology = request.shared_topology();
        match term.text() {
            "shown" => Ok(HostValue::Atoms(AtomSelection::all(topology))),
            "hidden" => Ok(HostValue::Atoms(AtomSelection::empty(topology))),
            other => Err(HostTermError::new(format!("unknown {other}"))),
        }
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    println!(
        "informational only; target={}-{}; atoms={ATOMS}; runs={RUNS}; median of 3 to {PASSES} passes within {BUDGET:?}",
        std::env::consts::OS,
        std::env::consts::ARCH,
    );
    let topology = topology()?;
    // Runs of 100 atoms, 200 apart.
    let step = ATOMS / RUNS;
    let values = (0..RUNS)
        .map(|run| {
            let start = run * step;
            format!("{start}-{}", start + step / 2 - 1)
        })
        .collect::<Vec<_>>()
        .join(" ");
    let options = SelectParseOptions::default()
        .with_max_input_bytes(1 << 20)?
        .with_host_keywords(["shown", "hidden"])?;
    let resolver = Objects;
    let evaluate = SelectEvaluateOptions::default().with_resolver(&resolver);
    println!("| expression | bytes | parses under the default limit | selected | ms |");
    println!("|---|---:|---|---:|---:|");
    for (name, text) in [
        ("index <runs>", format!("index {values}")),
        (
            "shown and index <runs>",
            format!("shown and index {values}"),
        ),
        (
            "hidden and index <runs>",
            format!("hidden and index {values}"),
        ),
        ("instance <runs>", format!("instance {values}")),
    ] {
        let default_limit = select::parse(&text, SelectParseOptions::default()).is_ok();
        let expression = select::parse(&text, options.clone())?;
        let mut selected = 0;
        let time = median(|| {
            let started = Instant::now();
            let atoms = expression
                .evaluate_with(Arc::clone(&topology), evaluate, Run::new())
                .expect("evaluates");
            let elapsed = started.elapsed();
            selected = black_box(atoms).len();
            elapsed
        });
        println!(
            "| `{name}` | {} | {} | {selected} | {} |",
            text.len(),
            if default_limit { "yes" } else { "no" },
            ms(time)
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
    let mut times = Vec::new();
    while times.len() < 3 || (times.len() < PASSES && times.iter().sum::<Duration>() < BUDGET) {
        times.push(pass());
    }
    times.sort_unstable();
    times[times.len() / 2]
}

fn ms(duration: Duration) -> String {
    format!("{:.3}", duration.as_secs_f64() * 1e3)
}
