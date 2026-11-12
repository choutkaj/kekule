//! A sample read from a file is allocated on the thread that reads it.
//!
//! A pooled reader serves reads from whichever thread takes it, so a buffer
//! it kept from one read would be allocated on one thread and freed on
//! another, which some allocators do not promptly reuse. This binary's
//! allocator tags every block with the thread that allocated it and counts
//! the blocks of one array's size that another thread frees.

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU32, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

use kekule::diagnostics::UnrepresentablePolicy;
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::structure::{
    Forces, Order, Place, Positions, Sample, SampleRead, SampleSet, TimeAxis, Velocities,
};
use kekule::topology::{Topology, TopologyBuilder};
use kekule::traj::io::dcd::{DcdCellPolicy, DcdWriteOptions};
use kekule::traj::io::trr::{TrrLambdaPolicy, TrrWriteOptions};
use kekule::traj::io::xyz::XyzWriteOptions;
use kekule::traj::io::{open_indexed_sample_set, write_sample_set, TrajectoryWriteOptions};
use kekule::traj::TrajectoryFormat;
use kekule::units::{
    Quantity, KILOJOULE_PER_MOLE_PER_NANOMETER, NANOMETER, NANOMETER_PER_PICOSECOND, PICOSECOND,
};

/// Atoms of the test system: an array of their positions has a size no
/// other block of a read has.
const ATOMS: usize = 1009;
/// Bytes of one array of positions, velocities, or forces.
const ARRAY_BYTES: usize = ATOMS * std::mem::size_of::<Point3>();

/// Room before each block for the tag of the thread that allocated it.
const HEADER: usize = 16;

/// The system allocator, tagging blocks with their thread.
struct Tagging;

#[global_allocator]
static ALLOCATOR: Tagging = Tagging;

/// Blocks of [`ARRAY_BYTES`] that a thread other than their own freed.
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

/// The layout of a block with its header in front.
fn tagged(layout: Layout) -> (Layout, usize) {
    let offset = layout.align().max(HEADER);
    let layout = Layout::from_size_align(layout.size() + offset, offset)
        .expect("test allocations are small");
    (layout, offset)
}

unsafe impl GlobalAlloc for Tagging {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let (full, offset) = tagged(layout);
        // SAFETY: `full` has a nonzero size.
        let base = unsafe { System.alloc(full) };
        if base.is_null() {
            return base;
        }
        // SAFETY: the block holds `offset` bytes before the caller's, and
        // `offset` is a multiple of the tag's alignment.
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
            if layout.size() == ARRAY_BYTES && owner != tag() {
                CROSS_THREAD_FREES.fetch_add(1, Ordering::Relaxed);
            }
            System.dealloc(block.sub(offset), full);
        }
    }
}

/// Tests share the counter, so they run one at a time.
static SERIAL: Mutex<()> = Mutex::new(());

fn topology() -> Arc<Topology> {
    let mut molecule = MoleculeEditor::new();
    let carbon = Element::from_symbol("C").unwrap();
    let mut previous = None;
    for _ in 0..ATOMS {
        let atom = molecule.add_atom(Atom::new(carbon)).unwrap();
        if let Some(previous) = previous {
            molecule
                .add_bond(previous, atom, BondOrder::Single)
                .unwrap();
        }
        previous = Some(atom);
    }
    let mut builder = TopologyBuilder::new();
    let definition = builder
        .add_molecule_definition(molecule.finish().unwrap())
        .unwrap();
    builder.add_instance(definition).unwrap();
    Arc::new(builder.build().unwrap())
}

/// Three samples of a moving chain in a cell, in time order, with
/// velocities and forces when `with_vectors` is set.
fn trajectory(topology: &Arc<Topology>, with_vectors: bool) -> SampleSet {
    let mut set = SampleSet::new(Arc::clone(topology));
    let axis = TimeAxis::new(
        Some(Quantity::new(Vec::new(), PICOSECOND)),
        Some(Vec::new()),
    )
    .unwrap();
    set.set_order(Some(Order::Time(axis))).unwrap();
    for step in 0..3_u32 {
        let shift = f64::from(step);
        let points = (0..ATOMS)
            .map(|atom| Point3::new(0.01 * atom as f64 + 0.1 * shift, 0.2, 0.3))
            .collect::<Vec<_>>();
        let vectors = vec![Vector3::new(0.5, -0.25, shift); ATOMS];
        let cell = PeriodicCell::orthorhombic(
            Quantity::new(Vector3::new(12.0, 12.0, 12.0), NANOMETER),
            [true; 3],
        )
        .unwrap();
        let mut sample =
            Sample::new(Positions::from_vec(Quantity::new(points, NANOMETER)).unwrap())
                .with_cell(cell);
        if with_vectors {
            sample
                .set_velocities(Some(
                    Velocities::new(Quantity::new(vectors.clone(), NANOMETER_PER_PICOSECOND))
                        .unwrap(),
                ))
                .unwrap();
            sample
                .set_forces(Some(
                    Forces::new(Quantity::new(vectors, KILOJOULE_PER_MOLE_PER_NANOMETER)).unwrap(),
                ))
                .unwrap();
        }
        let place = Place::Time {
            time: Some(Quantity::new(shift, PICOSECOND)),
            step: Some(u64::from(step)),
        };
        set.push(sample, Some(place)).unwrap();
    }
    set
}

/// A file of [`trajectory`] in `format`, removed on drop.
struct File(PathBuf);

impl File {
    fn write(topology: &Arc<Topology>, format: TrajectoryFormat) -> Self {
        let extension = format.to_string().to_ascii_lowercase();
        let path = std::env::temp_dir().join(format!(
            "kekule-reading-threads-{}.{extension}",
            std::process::id()
        ));
        let options = TrajectoryWriteOptions::new(format)
            .with_trr_options(
                TrrWriteOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero),
            )
            .with_dcd_options(
                DcdWriteOptions::default()
                    .with_cell_policy(DcdCellPolicy::Write)
                    .with_frame_metadata(UnrepresentablePolicy::Omit),
            )
            .with_xyz_options(
                XyzWriteOptions::default()
                    .with_cell(UnrepresentablePolicy::Omit)
                    .with_frame_metadata(UnrepresentablePolicy::Omit),
            );
        // Only TRR stores velocities and forces.
        let with_vectors = format == TrajectoryFormat::Trr;
        write_sample_set(&path, &trajectory(topology, with_vectors), options).unwrap();
        Self(path)
    }
}

impl Drop for File {
    fn drop(&mut self) {
        let _ = std::fs::remove_file(&self.0);
    }
}

#[test]
fn the_tagging_allocator_counts_an_array_freed_on_another_thread() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let array = vec![Point3::origin(); ATOMS];
    let before = CROSS_THREAD_FREES.load(Ordering::Relaxed);
    std::thread::spawn(move || drop(array)).join().unwrap();
    assert_eq!(CROSS_THREAD_FREES.load(Ordering::Relaxed), before + 1);
}

#[test]
fn a_read_on_another_thread_frees_only_its_own_arrays() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let topology = topology();
    for format in [
        TrajectoryFormat::Xyz,
        TrajectoryFormat::Dcd,
        TrajectoryFormat::Trr,
        TrajectoryFormat::Xtc,
    ] {
        let file = File::write(&topology, format);
        let set =
            open_indexed_sample_set(&file.0, Arc::clone(&topology), Default::default()).unwrap();
        // One thread at a time, so both reads take the pool's one reader.
        std::thread::scope(|scope| {
            scope.spawn(|| drop(set.read(0).unwrap()));
        });
        let before = CROSS_THREAD_FREES.load(Ordering::Relaxed);
        std::thread::scope(|scope| {
            scope.spawn(|| drop(set.read(1).unwrap()));
        });
        assert_eq!(
            CROSS_THREAD_FREES.load(Ordering::Relaxed),
            before,
            "a {format} read freed an array another thread allocated"
        );
    }
}

/// Runs `work` on a thread of its own.
fn on_thread(work: impl FnOnce() + Send) {
    std::thread::scope(|scope| {
        scope.spawn(work);
    });
}

/// Changes the file of `format` at `path` in place so that a read of frame
/// 2, at `offset`, fails after decoding it.
fn spoil_frame_2(path: &std::path::Path, format: TrajectoryFormat, offset: u64) {
    use std::io::{Seek, SeekFrom, Write};
    let mut bytes = std::fs::read(path).unwrap();
    match format {
        // The last atom's x, as text of the same width.
        TrajectoryFormat::Xyz => {
            let line = bytes[..bytes.len() - 1]
                .iter()
                .rposition(|byte| *byte == b'\n')
                .unwrap()
                + 1;
            let x = line + bytes[line..].iter().position(|byte| *byte == b' ').unwrap();
            let x = x + bytes[x..].iter().position(|byte| *byte != b' ').unwrap();
            let end = x + bytes[x..].iter().position(|byte| *byte == b' ').unwrap();
            let nan = format!("{:<width$}", "NaN", width = end - x);
            bytes[x..end].copy_from_slice(nan.as_bytes());
        }
        // The last value of frame 2, before its record's end marker.
        TrajectoryFormat::Dcd => {
            let at = bytes.len() - 8;
            bytes[at..at + 4].copy_from_slice(&f32::NAN.to_le_bytes());
        }
        // The last force value, as single or double precision XDR.
        TrajectoryFormat::Trr => {
            let at = bytes.len() - 8;
            bytes[at..].copy_from_slice(&f64::NAN.to_be_bytes());
        }
        // The step in frame 2's header, which the index recorded.
        TrajectoryFormat::Xtc => {
            let at = usize::try_from(offset).unwrap() + 8;
            bytes[at..at + 4].copy_from_slice(&1000_i32.to_be_bytes());
        }
        _ => unreachable!("the test writes four formats"),
    }
    let mut file = std::fs::OpenOptions::new().write(true).open(path).unwrap();
    file.seek(SeekFrom::Start(0)).unwrap();
    file.write_all(&bytes).unwrap();
}

#[test]
fn a_failed_read_keeps_no_array_for_the_next() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let topology = topology();
    for format in [
        TrajectoryFormat::Xyz,
        TrajectoryFormat::Dcd,
        TrajectoryFormat::Trr,
        TrajectoryFormat::Xtc,
    ] {
        let file = File::write(&topology, format);
        let set =
            open_indexed_sample_set(&file.0, Arc::clone(&topology), Default::default()).unwrap();
        spoil_frame_2(
            &file.0,
            format,
            set.storage().frame_index().frame_offsets()[2],
        );
        on_thread(|| assert!(set.read(2).is_err(), "{format}"));
        let before = CROSS_THREAD_FREES.load(Ordering::Relaxed);
        on_thread(|| drop(set.read(1).unwrap()));
        assert_eq!(
            CROSS_THREAD_FREES.load(Ordering::Relaxed),
            before,
            "a {format} read freed an array a failed read allocated"
        );
    }
}

/// A DCD of [`ATOMS`] atoms in three frames, every atom but each third one
/// fixed, so that frames 1 and 2 hold the free atoms alone.
fn fixed_atom_dcd() -> Vec<u8> {
    let free = (0..ATOMS).filter(|atom| atom % 3 == 0).collect::<Vec<_>>();
    let mut bytes = Vec::new();
    let mut record = |payload: &[u8]| {
        let length = i32::try_from(payload.len()).unwrap().to_le_bytes();
        bytes.extend(length);
        bytes.extend(payload);
        bytes.extend(length);
    };
    let mut header = [0_u8; 84];
    header[..4].copy_from_slice(b"CORD");
    let mut controls = [0_i32; 20];
    controls[0] = 3;
    controls[2] = 1;
    controls[8] = i32::try_from(ATOMS - free.len()).unwrap();
    controls[19] = 24;
    for (chunk, value) in header[4..].as_chunks_mut::<4>().0.iter_mut().zip(controls) {
        chunk.copy_from_slice(&value.to_le_bytes());
    }
    header[40..44].copy_from_slice(&1.0_f32.to_le_bytes());
    record(&header);
    let mut title = [0_u8; 84];
    title[..4].copy_from_slice(&1_i32.to_le_bytes());
    title[4..9].copy_from_slice(b"fixed");
    record(&title);
    record(&i32::try_from(ATOMS).unwrap().to_le_bytes());
    let indices = free
        .iter()
        .flat_map(|atom| i32::try_from(atom + 1).unwrap().to_le_bytes())
        .collect::<Vec<_>>();
    record(&indices);
    for frame in 0..3_u8 {
        let atoms = if frame == 0 {
            (0..ATOMS).collect()
        } else {
            free.clone()
        };
        for axis in 0..3_u8 {
            let values = atoms
                .iter()
                .flat_map(|&atom| {
                    (0.01 * atom as f32 + f32::from(axis) + 0.1 * f32::from(frame)).to_le_bytes()
                })
                .collect::<Vec<_>>();
            record(&values);
        }
    }
    bytes
}

#[test]
fn a_dcd_with_fixed_atoms_frees_only_its_own_arrays() {
    let _serial = SERIAL
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    let topology = topology();
    let path = std::env::temp_dir().join(format!(
        "kekule-reading-threads-fixed-{}.dcd",
        std::process::id()
    ));
    std::fs::write(&path, fixed_atom_dcd()).unwrap();
    let file = File(path);
    let set = open_indexed_sample_set(&file.0, Arc::clone(&topology), Default::default()).unwrap();
    assert_eq!(set.len(), 3);
    // Frame 1 first, then frame 0, then frame 2: each on a thread of its own.
    on_thread(|| drop(set.read(1).unwrap()));
    for index in [0, 2] {
        let before = CROSS_THREAD_FREES.load(Ordering::Relaxed);
        on_thread(|| drop(set.read(index).unwrap()));
        assert_eq!(
            CROSS_THREAD_FREES.load(Ordering::Relaxed),
            before,
            "a read of frame {index} freed an array another thread allocated"
        );
    }
}
