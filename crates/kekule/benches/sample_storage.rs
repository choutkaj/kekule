//! Informational cost of reading samples from an XTC file through a sample
//! set, against the indexed reader it wraps, and of processing them as they
//! are read, against the streaming adapters.
//! Run: cargo bench -p kekule --bench sample_storage -- SYSTEM.cif [TRAJECTORY.xtc]
//! The topology comes from an external mmCIF model, in its cell or a box
//! around its atoms. The frames
//! come from an external XTC file of the same atoms, or else are the model's
//! own coordinates repeated, to isolate storage and processing costs; this
//! measures API overhead, not physical trajectory sampling.

use std::error::Error;
use std::hint::black_box;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::mmcif;
use kekule::structure::{Sample, SampleRead};
use kekule::topology::{AtomSelection, Topology};
use kekule::traj::analysis::FrameSuperposer;
use kekule::traj::io::{
    create_trajectory_stream_writer, open_indexed_sample_set, open_indexed_trajectory,
    TrajectoryWriteOptions,
};
use kekule::traj::periodic::MoleculeImager;
use kekule::traj::{
    FitReference, FrameBuffer, Processed, ProcessingSteps, SeekableTrajectoryReader,
    TrajectoryFormat, TrajectoryReader, TrajectoryWriter,
};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};

/// Frames written when no trajectory is supplied.
const FRAMES: usize = 200;
/// Passes over every frame.
const PASSES: usize = 3;
/// Threads of the concurrent reads.
const THREADS: usize = 4;

fn main() -> Result<(), Box<dyn Error>> {
    // Cargo passes `--bench` after the paths.
    let mut arguments = std::env::args_os()
        .skip(1)
        .filter(|argument| !argument.to_string_lossy().starts_with("--"));
    let model_path = arguments
        .next()
        .ok_or("supply an external mmCIF model, and optionally its XTC trajectory")?;
    let model = mmcif::parse_str(&std::fs::read_to_string(model_path)?, Default::default())?
        .interpret(Default::default())?
        .into_model();
    // A model without a cell, such as a crystal structure of another space
    // group than P1, gets a box 1 nm longer than its extent on every axis.
    let cell = match model.cell() {
        Some(cell) => *cell,
        None => enclosing_cell(model.positions().values().value())?,
    };
    let topology = model.shared_topology();

    let scratch = std::env::temp_dir().join(format!("kekule-bench-samples-{}", std::process::id()));
    let path = match arguments.next() {
        Some(path) => PathBuf::from(path),
        None => {
            std::fs::create_dir_all(&scratch)?;
            let path = scratch.join("repeated.xtc");
            let mut frame = FrameBuffer::from_sample(
                Arc::clone(&topology),
                Sample::new(model.positions().clone()).with_cell(cell),
            )?;
            let file = std::io::BufWriter::new(std::fs::File::create(&path)?);
            let mut writer = create_trajectory_stream_writer(
                file,
                Arc::clone(&topology),
                TrajectoryWriteOptions::new(TrajectoryFormat::Xtc),
                "bench.xtc",
            )?;
            for index in 0..FRAMES {
                frame.set_step(Some(index as u64));
                frame.set_time(Some(Quantity::new(index as f64 * 0.002, PICOSECOND)))?;
                writer.write_frame(frame.frame_view())?;
            }
            writer.finish()?;
            path
        }
    };
    let result = benchmark(&path, &topology);
    let _ = std::fs::remove_dir_all(&scratch);
    result
}

fn benchmark(path: &Path, topology: &Arc<Topology>) -> Result<(), Box<dyn Error>> {
    let mut reader = open_indexed_trajectory(path, Arc::clone(topology), Default::default())?;
    let frames = reader.frame_count().ok_or("the trajectory has no index")?;
    let reads = frames * PASSES;
    let per_frame = |start: Instant| start.elapsed().as_secs_f64() * 1.0e6 / reads as f64;

    let mut buffer = reader.frame_buffer();
    let start = Instant::now();
    for _ in 0..PASSES {
        for index in 0..frames {
            reader.read_frame(index, &mut buffer)?;
            black_box(buffer.positions().len());
        }
    }
    let raw = per_frame(start);

    let set = open_indexed_sample_set(path, Arc::clone(topology), Default::default())?;
    let start = Instant::now();
    for _ in 0..PASSES {
        for index in 0..frames {
            black_box(set.read(index)?.positions().len());
        }
    }
    let read = per_frame(start);

    let start = Instant::now();
    std::thread::scope(|scope| {
        for thread in 0..THREADS {
            let set = &set;
            scope.spawn(move || {
                for _ in 0..PASSES {
                    for index in (thread..frames).step_by(THREADS) {
                        black_box(set.read(index).expect("benchmark read").positions().len());
                    }
                }
            });
        }
    });
    let threaded = per_frame(start);

    let atoms = AtomSelection::all(Arc::clone(topology));
    let imager = MoleculeImager::new(Arc::clone(topology));
    reader.read_frame(0, &mut buffer)?;
    imager.make_whole(0, &mut buffer)?;
    let reference = buffer.clone();
    let reference = reference.frame_view();
    let superposer = FrameSuperposer::new(&reference, &atoms, Default::default());
    let start = Instant::now();
    for _ in 0..PASSES {
        for index in 0..frames {
            reader.read_frame(index, &mut buffer)?;
            imager.make_whole(index, &mut buffer)?;
            superposer.superpose(index, &mut buffer)?;
            black_box(buffer.positions().len());
        }
    }
    let streamed = per_frame(start);

    let steps = ProcessingSteps::default()
        .with_whole(true)
        .with_fit(atoms, FitReference::Sample(0));
    let processed = Processed::new(&set, steps)?;
    let start = Instant::now();
    for _ in 0..PASSES {
        for index in 0..frames {
            black_box(processed.read(index)?.positions().len());
        }
    }
    let lazily = per_frame(start);

    println!(
        "informational only; target={}-{}; parallelism={}",
        std::env::consts::OS,
        std::env::consts::ARCH,
        std::thread::available_parallelism()?.get()
    );
    println!(
        "| atoms | frames | indexed reader us/frame | sample set read us/frame | overhead us/frame | {THREADS} threads us/frame | streaming whole+fit us/frame | processed whole+fit us/frame |"
    );
    println!("|---:|---:|---:|---:|---:|---:|---:|---:|");
    println!(
        "| {} | {frames} | {raw:.2} | {read:.2} | {:.2} | {threaded:.2} | {streamed:.2} | {lazily:.2} |",
        topology.atom_count(),
        read - raw,
    );
    Ok(())
}

/// An orthorhombic box 1 nm longer than the extent of `points` on every
/// axis.
fn enclosing_cell(points: &[Point3]) -> Result<PeriodicCell, Box<dyn Error>> {
    let first = *points.first().ok_or("the model has no atoms")?;
    let (mut low, mut high) = (first, first);
    for point in points {
        low = Point3::new(low.x.min(point.x), low.y.min(point.y), low.z.min(point.z));
        high = Point3::new(
            high.x.max(point.x),
            high.y.max(point.y),
            high.z.max(point.z),
        );
    }
    let lengths = (high - low) + Vector3::new(1.0, 1.0, 1.0);
    Ok(PeriodicCell::orthorhombic(
        Quantity::new(lengths, NANOMETER),
        [true; 3],
    )?)
}
