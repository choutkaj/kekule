//! Informational cost of reading samples from an XTC file through a sample
//! set, against the indexed reader it wraps, and of processing them as they
//! are read, against the streaming adapters.
//! Run: cargo bench -p kekule --bench sample_storage

use std::error::Error;
use std::hint::black_box;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::molecule::{Atom, BondOrder, Element, MoleculeEditor};
use kekule::structure::SampleRead;
use kekule::topology::{AtomSelection, Topology, TopologyBuilder};
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

#[derive(Clone, Copy)]
struct Profile {
    name: &'static str,
    molecules: usize,
    frames: usize,
    passes: usize,
}

fn main() -> Result<(), Box<dyn Error>> {
    println!(
        "informational only; target={}-{}; parallelism={}",
        std::env::consts::OS,
        std::env::consts::ARCH,
        std::thread::available_parallelism()?.get()
    );
    println!(
        "| profile | atoms | frames | indexed reader us/frame | sample set read us/frame | overhead us/frame | 4 threads us/frame | streaming whole+fit us/frame | processed whole+fit us/frame |"
    );
    println!("|---|---:|---:|---:|---:|---:|---:|---:|---:|");
    let directory =
        std::env::temp_dir().join(format!("kekule-bench-samples-{}", std::process::id()));
    std::fs::create_dir_all(&directory)?;
    let result = [
        Profile {
            name: "small",
            molecules: 22,
            frames: 1_000,
            passes: 5,
        },
        Profile {
            name: "large",
            molecules: 4_096,
            frames: 100,
            passes: 3,
        },
    ]
    .into_iter()
    .try_for_each(|profile| benchmark(profile, &directory));
    std::fs::remove_dir_all(&directory)?;
    result
}

fn benchmark(profile: Profile, directory: &Path) -> Result<(), Box<dyn Error>> {
    let topology = topology(profile.molecules)?;
    let path = directory.join(format!("{}.xtc", profile.name));
    write(&path, &topology, profile.frames)?;
    let reads = profile.frames * profile.passes;
    let per_frame = |start: Instant| start.elapsed().as_secs_f64() * 1.0e6 / reads as f64;

    let mut reader = open_indexed_trajectory(&path, Arc::clone(&topology), Default::default())?;
    let mut buffer = reader.frame_buffer();
    let start = Instant::now();
    for _ in 0..profile.passes {
        for index in 0..profile.frames {
            reader.read_frame(index, &mut buffer)?;
            black_box(buffer.positions().len());
        }
    }
    let raw = per_frame(start);

    let set = open_indexed_sample_set(&path, Arc::clone(&topology), Default::default())?;
    let start = Instant::now();
    for _ in 0..profile.passes {
        for index in 0..profile.frames {
            black_box(set.read(index)?.positions().len());
        }
    }
    let read = per_frame(start);

    let threads = 4;
    let start = Instant::now();
    std::thread::scope(|scope| {
        for thread in 0..threads {
            let set = &set;
            scope.spawn(move || {
                for _ in 0..profile.passes {
                    for index in (thread..profile.frames).step_by(threads) {
                        black_box(set.read(index).expect("benchmark read").positions().len());
                    }
                }
            });
        }
    });
    let threaded = per_frame(start);

    let atoms = AtomSelection::all(Arc::clone(&topology));
    let imager = MoleculeImager::new(Arc::clone(&topology));
    reader.read_frame(0, &mut buffer)?;
    imager.make_whole(0, &mut buffer)?;
    let reference = buffer.clone();
    let reference = reference.frame_view();
    let superposer = FrameSuperposer::new(&reference, &atoms, Default::default());
    let start = Instant::now();
    for _ in 0..profile.passes {
        for index in 0..profile.frames {
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
    for _ in 0..profile.passes {
        for index in 0..profile.frames {
            black_box(processed.read(index)?.positions().len());
        }
    }
    let lazily = per_frame(start);

    println!(
        "| {} | {} | {} | {raw:.2} | {read:.2} | {:.2} | {threaded:.2} | {streamed:.2} | {lazily:.2} |",
        profile.name,
        topology.atom_count(),
        profile.frames,
        read - raw,
    );
    Ok(())
}

/// `molecules` three-carbon chains.
fn topology(molecules: usize) -> Result<Arc<Topology>, Box<dyn Error>> {
    let carbon = Element::from_symbol("C").ok_or("carbon is a built-in element")?;
    let mut chain = MoleculeEditor::new();
    let atoms = [(); 3].map(|()| chain.add_atom(Atom::new(carbon)));
    let [first, second, third] = atoms;
    let (first, second, third) = (first?, second?, third?);
    chain.add_bond(first, second, BondOrder::Single)?;
    chain.add_bond(second, third, BondOrder::Single)?;
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(chain.finish()?)?;
    for _ in 0..molecules {
        builder.add_instance(definition)?;
    }
    Ok(Arc::new(builder.build()?))
}

/// Frames of chains on a grid in a 5 nm cell, drifting across its faces.
fn write(path: &Path, topology: &Arc<Topology>, frames: usize) -> Result<(), Box<dyn Error>> {
    let file = std::io::BufWriter::new(std::fs::File::create(path)?);
    let mut writer = create_trajectory_stream_writer(
        file,
        Arc::clone(topology),
        TrajectoryWriteOptions::new(TrajectoryFormat::Xtc),
        "bench.xtc",
    )?;
    let length = 5.0;
    let cell = PeriodicCell::orthorhombic(
        Quantity::new(Vector3::new(length, length, length), NANOMETER),
        [true; 3],
    )?;
    let molecules = topology.atom_count() / 3;
    let side = (molecules as f64).cbrt().ceil() as usize;
    let mut frame = FrameBuffer::new(Arc::clone(topology));
    for index in 0..frames {
        let drift = 0.01 * index as f64;
        let points = (0..topology.atom_count())
            .map(|atom| {
                let molecule = atom / 3;
                let spacing = length / side as f64;
                let base = Vector3::new(
                    (molecule % side) as f64,
                    (molecule / side % side) as f64,
                    (molecule / side / side) as f64,
                ) * spacing;
                let offset = Vector3::new(
                    0.14 * (atom % 3) as f64 + drift,
                    0.05 * (atom % 3) as f64,
                    drift,
                );
                let point = Point3::origin() + base + offset;
                Point3::new(
                    point.x.rem_euclid(length),
                    point.y.rem_euclid(length),
                    point.z.rem_euclid(length),
                )
            })
            .collect::<Vec<_>>();
        frame
            .frame_mut()
            .set_positions(Quantity::new(points, NANOMETER))?;
        frame.frame_mut().set_cell(Some(cell));
        frame.frame_mut().set_step(Some(index as u64));
        frame
            .frame_mut()
            .set_time(Some(Quantity::new(index as f64 * 0.002, PICOSECOND)))?;
        writer.write_frame(frame.frame_view())?;
    }
    writer.finish()?;
    Ok(())
}
