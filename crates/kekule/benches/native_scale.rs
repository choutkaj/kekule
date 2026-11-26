//! Informational cost of native persistence of a structure of about one
//! million atoms with full per-site metadata: the time to write and read its
//! topology alone and the whole structure, in memory and through a buffered
//! file, the encoding's size and SHA-256, and the heap a read holds at its
//! peak and keeps.
//!
//! The structure is merged copies of an externally supplied mmCIF entry, as
//! a host merges them: each copy after the first takes fresh chain
//! identifiers, on its chains and its atom sites. SeekDelta's project files
//! hold 93 copies of PDB 8QH3 (10,796 atoms), which make 1,004,028 atoms:
//!
//! ```sh
//! curl -sSfLO https://files.rcsb.org/download/8QH3.cif
//! cargo bench -p kekule --bench native_scale -- 8QH3.cif 93
//! ```
//!
//! The arguments are the mmCIF path and the number of copies; the
//! `KEKULE_NATIVE_SCALE_CIF` and `KEKULE_NATIVE_SCALE_COPIES` environment
//! variables may supply them instead. Without a path the benchmark reports
//! nothing. Times come from rounds that count no allocations; one more
//! round of each operation counts the heap. The file rows write to and read
//! from a file in the temporary directory, which is deleted.

use std::alloc::{GlobalAlloc, Layout, System};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::hint::black_box;
use std::io::{BufReader, BufWriter, Write};
use std::path::Path;
use std::sync::atomic::{AtomicBool, AtomicIsize, Ordering};
use std::time::{Duration, Instant};

use kekule::mmcif;
use kekule::native::{self, NativeObject, NativeReadOptions, NativeSource, NativeWriteOptions};
use kekule::structure::{Structure, StructureEditor};
use sha2::{Digest, Sha256};

const COPIES: usize = 93;
const ROUNDS: usize = 7;

/// The system allocator, counting the heap while [`COUNTING`] is set.
struct Counting;

#[global_allocator]
static ALLOCATOR: Counting = Counting;

static COUNTING: AtomicBool = AtomicBool::new(false);
/// Bytes allocated less bytes freed since counting started, now and at the
/// peak.
static CURRENT: AtomicIsize = AtomicIsize::new(0);
static PEAK: AtomicIsize = AtomicIsize::new(0);

fn count(bytes: isize) {
    if COUNTING.load(Ordering::Relaxed) {
        let now = CURRENT.fetch_add(bytes, Ordering::Relaxed) + bytes;
        PEAK.fetch_max(now, Ordering::Relaxed);
    }
}

// SAFETY: every block comes from `System` and returns to it with its layout.
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        // SAFETY: the caller's layout has a nonzero size.
        let block = unsafe { System.alloc(layout) };
        if !block.is_null() {
            count(layout.size() as isize);
        }
        block
    }

    unsafe fn dealloc(&self, block: *mut u8, layout: Layout) {
        count(-(layout.size() as isize));
        // SAFETY: `block` came from `alloc` with this layout.
        unsafe { System.dealloc(block, layout) }
    }

    unsafe fn realloc(&self, block: *mut u8, layout: Layout, size: usize) -> *mut u8 {
        // SAFETY: `block` came from `alloc` with this layout, and `size` is
        // nonzero.
        let moved = unsafe { System.realloc(block, layout, size) };
        if !moved.is_null() {
            count(size as isize - layout.size() as isize);
        }
        moved
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    // Cargo passes `--bench` after the arguments.
    let mut arguments = std::env::args()
        .skip(1)
        .filter(|argument| !argument.starts_with("--"));
    let Some(path) = arguments
        .next()
        .or_else(|| std::env::var("KEKULE_NATIVE_SCALE_CIF").ok())
    else {
        println!("native_scale: no mmCIF path given; nothing measured");
        return Ok(());
    };
    let copies = match arguments
        .next()
        .or_else(|| std::env::var("KEKULE_NATIVE_SCALE_COPIES").ok())
    {
        Some(copies) => copies.parse()?,
        None => COPIES,
    };
    let start = Instant::now();
    let structure = merged_copies(&path, copies)?;
    let built = start.elapsed();
    let topology = structure.topology();
    let bytes = native::write(&structure, NativeWriteOptions::default())?;
    println!(
        "informational only; target={}-{}; {copies} copies built in {:.1} s: atoms={}; \
         sites={}; residues={}; chains={}; definitions={}; instances={}; bytes={}; \
         sha256={}; minimum and median of {ROUNDS} rounds",
        std::env::consts::OS,
        std::env::consts::ARCH,
        built.as_secs_f64(),
        structure.atom_count(),
        topology.atom_sites().len(),
        topology.residues().len(),
        topology.chains().len(),
        topology.definition_count(),
        topology.instance_count(),
        bytes.len(),
        Sha256::digest(&bytes)
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect::<String>(),
    );
    drop(bytes);
    let scratch = std::env::temp_dir().join(format!("kekule-bench-native-{}", std::process::id()));
    std::fs::create_dir_all(&scratch)?;
    let result = benchmark(&structure, &scratch.join("structure.kekule"));
    let _ = std::fs::remove_dir_all(&scratch);
    result
}

/// `copies` copies of the entry at `path`, perceived and merged into one
/// structure; each copy after the first takes fresh chain identifiers.
fn merged_copies(path: &str, copies: usize) -> Result<Structure, Box<dyn Error>> {
    let text = std::fs::read_to_string(path)?;
    let mut copy = mmcif::parse_str(&text, Default::default())?
        .interpret(Default::default())?
        .into_structure();
    copy.perceive()?;
    let mut editor = StructureEditor::new();
    let (mut labels, mut authors) = (0, 0);
    for index in 0..copies {
        let appended = editor.append_structure(&copy)?;
        if index == 0 {
            continue;
        }
        // The chains of one copy that share an author identifier keep
        // sharing it.
        let mut written_authors = HashMap::new();
        for chain in copy.topology().chains() {
            let author = chain.author_id().unwrap_or(chain.label_id()).to_owned();
            let author = written_authors
                .entry(author)
                .or_insert_with(|| {
                    authors += 1;
                    fresh_id(authors)
                })
                .clone();
            labels += 1;
            let label = fresh_id(labels);
            for site in chain.residues().flat_map(|residue| residue.atom_sites()) {
                let old = site.metadata();
                let metadata = old
                    .clone()
                    .with_auth_asym_id(old.auth_asym_id().map(|_| author.clone()))
                    .with_label_asym_id(old.label_asym_id().map(|_| label.clone()));
                let id = appended
                    .atom_site(site.id())
                    .ok_or("a site was not appended")?;
                editor.set_atom_site_metadata(id, metadata)?;
            }
            let id = appended
                .chain(chain.id())
                .ok_or("a chain was not appended")?;
            editor.set_chain_identifiers(id, label, Some(author))?;
        }
    }
    Ok(editor.finish()?)
}

/// The `index`th chain identifier after the source's: `A`..`Z`, `a`..`z`,
/// `0`..`9`, then longer ones in the same order.
fn fresh_id(index: usize) -> String {
    const SYMBOLS: &[u8] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
    let mut index = index;
    let mut id = Vec::new();
    loop {
        id.push(SYMBOLS[index % SYMBOLS.len()]);
        index /= SYMBOLS.len();
        if index == 0 {
            break;
        }
        index -= 1;
    }
    id.reverse();
    String::from_utf8(id).expect("symbols are ASCII")
}

/// How a round measures its operations.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Mode {
    /// The time of each operation, counting no allocations.
    Time,
    /// The heap each operation holds at its peak and keeps.
    Heap,
}

/// One measured operation: its time in each timed round, and the heap above
/// its start, at its peak and kept by its result, in the counting round.
#[derive(Default)]
struct Row {
    times: Vec<Duration>,
    peak: isize,
    kept: isize,
}

impl Row {
    fn measure<T>(
        &mut self,
        mode: Mode,
        operation: impl FnOnce() -> Result<T, Box<dyn Error>>,
    ) -> Result<T, Box<dyn Error>> {
        if mode == Mode::Time {
            let start = Instant::now();
            let value = operation()?;
            self.times.push(start.elapsed());
            return Ok(value);
        }
        CURRENT.store(0, Ordering::Relaxed);
        PEAK.store(0, Ordering::Relaxed);
        COUNTING.store(true, Ordering::Relaxed);
        let value = operation();
        COUNTING.store(false, Ordering::Relaxed);
        self.peak = PEAK.load(Ordering::Relaxed);
        self.kept = CURRENT.load(Ordering::Relaxed);
        value
    }

    fn print(&mut self, name: &str) {
        self.times.sort();
        let milliseconds = |time: Duration| time.as_secs_f64() * 1e3;
        let megabytes = |bytes: isize| bytes as f64 / f64::from(1 << 20);
        println!(
            "| {name} | {:.1} | {:.1} | {:.0} | {:.0} |",
            milliseconds(self.times[0]),
            milliseconds(self.times[self.times.len() / 2]),
            megabytes(self.peak),
            megabytes(self.kept),
        );
    }
}

/// The operations measured, in order; each read reads what the write
/// before it wrote.
const OPERATIONS: [&str; 6] = [
    "write topology",
    "read topology",
    "write structure",
    "read structure",
    "write structure to file",
    "read structure from file",
];

fn benchmark(structure: &Structure, file: &Path) -> Result<(), Box<dyn Error>> {
    let mut rows: Vec<Row> = OPERATIONS.iter().map(|_| Row::default()).collect();
    for _ in 0..ROUNDS {
        round(&mut rows, structure, file, Mode::Time)?;
    }
    round(&mut rows, structure, file, Mode::Heap)?;
    println!("| operation | min ms | median ms | peak heap MiB | kept heap MiB |");
    println!("|---|---:|---:|---:|---:|");
    for (row, name) in rows.iter_mut().zip(OPERATIONS) {
        row.print(name);
    }
    Ok(())
}

/// One round of every operation.
fn round(
    rows: &mut [Row],
    structure: &Structure,
    file: &Path,
    mode: Mode,
) -> Result<(), Box<dyn Error>> {
    let write = NativeWriteOptions::default();
    let read = NativeReadOptions::default();
    let bytes = rows[0].measure(mode, || Ok(native::write(structure.topology(), write)?))?;
    let decoded = rows[1].measure(mode, || Ok(native::read(bytes.as_slice(), read)?))?;
    drop(black_box(decoded));
    drop(bytes);
    let bytes = rows[2].measure(mode, || Ok(native::write(structure, write)?))?;
    let decoded = rows[3].measure(mode, || Ok(native::read(bytes.as_slice(), read)?))?;
    check(&decoded, structure)?;
    drop(black_box(decoded));
    drop(bytes);
    rows[4].measure(mode, || write_file(file, structure.into()))?;
    let decoded = rows[5].measure(mode, || {
        Ok(native::read(BufReader::new(File::open(file)?), read)?)
    })?;
    drop(black_box(decoded));
    Ok(())
}

fn write_file(file: &Path, source: NativeSource<'_>) -> Result<(), Box<dyn Error>> {
    let mut writer = BufWriter::new(File::create(file)?);
    native::write_to(&mut writer, source, NativeWriteOptions::default())?;
    writer.flush()?;
    Ok(())
}

/// Fails unless `decoded` is a structure of as many atoms and sites as
/// `source`.
fn check(decoded: &NativeObject, source: &Structure) -> Result<(), Box<dyn Error>> {
    let NativeObject::Structure(decoded) = decoded else {
        return Err("the input decoded to another kind".into());
    };
    if decoded.atom_count() != source.atom_count()
        || decoded.topology().atom_sites().len() != source.topology().atom_sites().len()
    {
        return Err("the decoded structure differs in size".into());
    }
    Ok(())
}
