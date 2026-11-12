//! Informational cost of reading one sample set over an XTC file in bursts
//! from a fixed pool of threads, as a host's thread pool reads it: each
//! burst reads one frame on every thread and waits for all of them, so the
//! set's pooled readers move between threads from burst to burst.
//!
//! The bench's allocator tags every block with the thread that allocated
//! it and counts the sample arrays (blocks of one frame's positions) that
//! another thread frees. Some allocators, such as mimalloc, do not promptly
//! reuse memory another thread freed, so such frees grow a host's memory.
//! Run: cargo bench -p kekule --bench pooled_reads -- [ATOMS] [FRAMES] [THREADS]
//! The coordinates are synthetic; this measures storage costs, not
//! physical trajectory sampling.

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::error::Error;
use std::hint::black_box;
use std::io::Write;
use std::sync::atomic::{AtomicU32, AtomicUsize, Ordering};
use std::sync::{Arc, Barrier};
use std::time::Instant;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::structure::SampleRead;
use kekule::topology::{Topology, TopologyBuilder};
use kekule::traj::io::open_indexed_sample_set;
use kekule::traj::io::xtc::{XtcWriteOptions, XtcWriter};
use kekule::traj::{FrameBuffer, TrajectoryWriter};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};

/// Atoms when none are requested.
const ATOMS: usize = 100_000;
/// Frames when none are requested.
const FRAMES: usize = 400;
/// Passes over every frame.
const PASSES: usize = 3;

/// Room before each block for the tag of the thread that allocated it.
const HEADER: usize = 16;

/// The system allocator, tagging blocks with their thread.
struct Tagging;

#[global_allocator]
static ALLOCATOR: Tagging = Tagging;

/// Bytes of the blocks counted: one frame's positions.
static COUNTED_BYTES: AtomicUsize = AtomicUsize::new(usize::MAX);
/// Counted blocks freed on another thread than the one that allocated them.
static CROSS_THREAD_FREES: AtomicUsize = AtomicUsize::new(0);
static NEXT_TAG: AtomicU32 = AtomicU32::new(1);

thread_local! {
    static TAG: Cell<u32> = const { Cell::new(0) };
}

/// This thread's tag; 0 while the thread is being torn down.
fn tag() -> u32 {
    TAG.try_with(|tag| {
        if tag.get() == 0 {
            tag.set(NEXT_TAG.fetch_add(1, Ordering::Relaxed));
        }
        tag.get()
    })
    .unwrap_or(0)
}

/// The layout of a block with its header in front, and the header's size.
fn tagged(layout: Layout) -> (Layout, usize) {
    let offset = layout.align().max(HEADER);
    let size = layout
        .size()
        .checked_add(offset)
        .expect("allocations stay below the address range");
    let layout = Layout::from_size_align(size, offset).expect("a valid layout");
    (layout, offset)
}

// SAFETY: blocks come from `System` with a header in front, and are
// returned to it with the same layout.
unsafe impl GlobalAlloc for Tagging {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let (full, offset) = tagged(layout);
        // SAFETY: `full` has a nonzero size.
        let base = unsafe { System.alloc(full) };
        if base.is_null() {
            return base;
        }
        // SAFETY: the block holds `offset >= HEADER` bytes before the
        // caller's, aligned for the tag.
        unsafe {
            let block = base.add(offset);
            block.sub(HEADER).cast::<u32>().write(tag());
            block
        }
    }

    unsafe fn dealloc(&self, block: *mut u8, layout: Layout) {
        let (full, offset) = tagged(layout);
        // SAFETY: `block` came from `alloc` with this layout.
        unsafe {
            let owner = block.sub(HEADER).cast::<u32>().read();
            if layout.size() == COUNTED_BYTES.load(Ordering::Relaxed) && owner != tag() {
                CROSS_THREAD_FREES.fetch_add(1, Ordering::Relaxed);
            }
            System.dealloc(block.sub(offset), full);
        }
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    // Cargo passes `--bench` after the arguments.
    let mut arguments = std::env::args()
        .skip(1)
        .filter(|argument| !argument.starts_with("--"));
    let atoms = arguments.next().map_or(Ok(ATOMS), |atoms| atoms.parse())?;
    let frames = arguments
        .next()
        .map_or(Ok(FRAMES), |frames| frames.parse())?;
    let threads = match arguments.next() {
        Some(threads) => threads.parse()?,
        None => std::thread::available_parallelism()?.get(),
    };
    let topology = topology(atoms)?;
    let scratch = std::env::temp_dir().join(format!("kekule-bench-pooled-{}", std::process::id()));
    std::fs::create_dir_all(&scratch)?;
    let path = scratch.join("synthetic.xtc");
    write(&path, &topology, frames)?;
    let result = benchmark(&path, &topology, threads);
    let _ = std::fs::remove_dir_all(&scratch);
    result
}

fn benchmark(
    path: &std::path::Path,
    topology: &Arc<Topology>,
    threads: usize,
) -> Result<(), Box<dyn Error>> {
    let set = open_indexed_sample_set(path, Arc::clone(topology), Default::default())?;
    let frames = set.len();
    let reads = frames * PASSES;
    COUNTED_BYTES.store(
        topology.atom_count() * std::mem::size_of::<Point3>(),
        Ordering::Relaxed,
    );

    let start = Instant::now();
    for _ in 0..PASSES {
        for index in 0..frames {
            black_box(set.read(index)?);
        }
    }
    let sequential = start.elapsed().as_secs_f64() * 1.0e6 / reads as f64;

    let bursts = frames.div_ceil(threads);
    let start_burst = Barrier::new(threads + 1);
    let end_burst = Barrier::new(threads + 1);
    let before = CROSS_THREAD_FREES.load(Ordering::Relaxed);
    let start = Instant::now();
    std::thread::scope(|scope| {
        for thread in 0..threads {
            let (set, start_burst, end_burst) = (&set, &start_burst, &end_burst);
            scope.spawn(move || {
                for _ in 0..PASSES {
                    for burst in 0..bursts {
                        start_burst.wait();
                        let index = burst * threads + thread;
                        if index < frames {
                            black_box(set.read(index).expect("benchmark read"));
                        }
                        end_burst.wait();
                    }
                }
            });
        }
        for _ in 0..PASSES * bursts {
            start_burst.wait();
            end_burst.wait();
        }
    });
    let pooled = start.elapsed().as_secs_f64() * 1.0e6 / reads as f64;
    let cross = CROSS_THREAD_FREES.load(Ordering::Relaxed) - before;

    println!(
        "{{\"atoms\":{},\"frames\":{frames},\"threads\":{threads},\"sequential_us_per_read\":{sequential:.1},\"bursts_us_per_read\":{pooled:.1},\"arrays_freed_on_another_thread_per_read\":{:.3}}}",
        topology.atom_count(),
        cross as f64 / reads as f64
    );
    Ok(())
}

/// A chain of `atoms` carbons.
fn topology(atoms: usize) -> Result<Arc<Topology>, Box<dyn Error>> {
    let mut molecule = MoleculeEditor::new();
    let carbon = Element::from_symbol("C").ok_or("carbon is a built-in element")?;
    let mut previous = None;
    for _ in 0..atoms {
        let atom = molecule.add_atom(Atom::new(carbon))?;
        if let Some(previous) = previous {
            molecule.add_bond(previous, atom, BondOrder::Single)?;
        }
        previous = Some(atom);
    }
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(molecule.finish()?)?;
    builder.add_instance(definition)?;
    Ok(Arc::new(builder.build()?))
}

/// `frames` frames of the atoms on a grid, each moved by a fixed
/// pseudo-random amount.
fn write(
    path: &std::path::Path,
    topology: &Arc<Topology>,
    frames: usize,
) -> Result<(), Box<dyn Error>> {
    let atoms = topology.atom_count();
    let side = (atoms as f64).cbrt().ceil() as usize;
    let edge = side as f64 * 0.2 + 1.0;
    let file = std::io::BufWriter::new(std::fs::File::create(path)?);
    let mut writer = XtcWriter::new(
        file,
        Arc::clone(topology),
        XtcWriteOptions::default(),
        "bench.xtc",
    )?;
    let mut frame = FrameBuffer::new(Arc::clone(topology));
    frame.sample_mut().set_cell(Some(PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(edge, edge, edge), NANOMETER),
        [true; 3],
    )?));
    let mut state = 0x2545_f491_4f6c_dd1d_u64;
    let mut noise = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        ((state >> 11) as f64 / (1_u64 << 53) as f64 - 0.5) * 0.05
    };
    for index in 0..frames {
        let positions = (0..atoms)
            .map(|atom| {
                let (x, y, z) = (atom % side, (atom / side) % side, atom / (side * side));
                Point3::new(
                    x as f64 * 0.2 + noise(),
                    y as f64 * 0.2 + noise(),
                    z as f64 * 0.2 + noise(),
                )
            })
            .collect::<Vec<_>>();
        frame
            .sample_mut()
            .set_positions(Quantity::new(positions, NANOMETER))?;
        frame.set_step(Some(index as u64));
        frame.set_time(Some(Quantity::new(index as f64 * 0.002, PICOSECOND)))?;
        writer.write_frame(frame.frame_view())?;
    }
    writer.finish()?.flush()?;
    Ok(())
}
