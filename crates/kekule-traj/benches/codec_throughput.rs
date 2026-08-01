//! Informational trajectory codec throughput, indexing, and random-access run.

use std::error::Error;
use std::hint::black_box;
use std::io::Cursor;
use std::time::{Duration, Instant};

use kekule::core::{Atom, Element, Molecule};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::small::SmallMolecule;
use kekule::topology::{MoleculeInstanceMetadata, Topology, TopologyBuilder};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};
use kekule_traj::io::dcd::{DcdReadOptions, DcdReader, DcdWriteOptions, DcdWriter};
use kekule_traj::io::trr::{
    TrrLambdaPolicy, TrrReadOptions, TrrReader, TrrWriteOptions, TrrWriter,
};
use kekule_traj::io::xtc::{XtcReadOptions, XtcReader, XtcWriteOptions, XtcWriter};
use kekule_traj::io::xyz::{XyzReadOptions, XyzReader, XyzWriteOptions, XyzWriter};
use kekule_traj::io::{TrajectoryIoLimits, TrajectoryTopologyBinding};
use kekule_traj::{
    AtomOrderAssertion, FrameBuffer, SeekableTrajectoryReader, TrajectoryReader, TrajectoryWriter,
};

#[derive(Clone, Copy)]
struct Profile {
    name: &'static str,
    atoms: usize,
    frames: usize,
    passes: usize,
    random_reads: usize,
}

fn main() -> Result<(), Box<dyn Error>> {
    println!(
        "informational only; target={}-{}; parallelism={}",
        std::env::consts::OS,
        std::env::consts::ARCH,
        std::thread::available_parallelism()?.get()
    );
    println!(
        "| codec | profile | atoms | frames | MiB | sequential frames/s | sequential MiB/s | index ms | index KiB | random us/frame | writer MiB/s |"
    );
    println!("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|");
    for profile in [
        Profile {
            name: "small",
            atoms: 64,
            frames: 256,
            passes: 5,
            random_reads: 256,
        },
        Profile {
            name: "large",
            atoms: 4096,
            frames: 32,
            passes: 3,
            random_reads: 64,
        },
        Profile {
            name: "high-frame",
            atoms: 1,
            frames: 100_000,
            passes: 1,
            random_reads: 1_000,
        },
    ] {
        benchmark_xyz(profile)?;
        benchmark_dcd(profile)?;
        benchmark_trr(profile)?;
        benchmark_xtc(profile)?;
    }
    Ok(())
}

fn benchmark_xyz(profile: Profile) -> Result<(), Box<dyn Error>> {
    let topology = topology(profile.atoms)?;
    let bytes = encode_xyz(&topology, profile.frames)?;
    let writer_start = Instant::now();
    for _ in 0..profile.passes {
        black_box(encode_xyz(&topology, profile.frames)?);
    }
    let writer_elapsed = writer_start.elapsed();
    benchmark_readers(
        "XYZ",
        profile,
        &topology,
        &bytes,
        writer_elapsed,
        || {
            XyzReader::new(
                Cursor::new(bytes.as_slice()),
                binding(&topology),
                XyzReadOptions::default(),
                TrajectoryIoLimits::default(),
                "bench.xyz",
            )
            .expect("benchmark XYZ reader")
        },
        || {
            XyzReader::new(
                Cursor::new(bytes.as_slice()),
                binding(&topology),
                XyzReadOptions::default(),
                TrajectoryIoLimits::default(),
                "bench.xyz",
            )
            .expect("benchmark XYZ reader")
            .into_indexed()
            .expect("benchmark XYZ index")
        },
    )
}

fn benchmark_dcd(profile: Profile) -> Result<(), Box<dyn Error>> {
    let topology = topology(profile.atoms)?;
    let bytes = encode_dcd(&topology, profile.frames)?;
    let writer_start = Instant::now();
    for _ in 0..profile.passes {
        black_box(encode_dcd(&topology, profile.frames)?);
    }
    let writer_elapsed = writer_start.elapsed();
    benchmark_readers(
        "DCD",
        profile,
        &topology,
        &bytes,
        writer_elapsed,
        || {
            DcdReader::new(
                Cursor::new(bytes.as_slice()),
                binding(&topology),
                DcdReadOptions::default(),
                TrajectoryIoLimits::default(),
                "bench.dcd",
            )
            .expect("benchmark DCD reader")
        },
        || {
            DcdReader::new(
                Cursor::new(bytes.as_slice()),
                binding(&topology),
                DcdReadOptions::default(),
                TrajectoryIoLimits::default(),
                "bench.dcd",
            )
            .expect("benchmark DCD reader")
            .into_indexed()
            .expect("benchmark DCD index")
        },
    )
}

fn benchmark_trr(profile: Profile) -> Result<(), Box<dyn Error>> {
    let topology = topology(profile.atoms)?;
    let bytes = encode_trr(&topology, profile.frames)?;
    let writer_start = Instant::now();
    for _ in 0..profile.passes {
        black_box(encode_trr(&topology, profile.frames)?);
    }
    let writer_elapsed = writer_start.elapsed();
    benchmark_readers(
        "TRR-f32",
        profile,
        &topology,
        &bytes,
        writer_elapsed,
        || {
            TrrReader::new(
                Cursor::new(bytes.as_slice()),
                binding(&topology),
                TrrReadOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero),
                TrajectoryIoLimits::default(),
                "bench.trr",
            )
            .expect("benchmark TRR reader")
        },
        || {
            TrrReader::new(
                Cursor::new(bytes.as_slice()),
                binding(&topology),
                TrrReadOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero),
                TrajectoryIoLimits::default(),
                "bench.trr",
            )
            .expect("benchmark TRR reader")
            .into_indexed()
            .expect("benchmark TRR index")
        },
    )
}

fn benchmark_xtc(profile: Profile) -> Result<(), Box<dyn Error>> {
    let topology = topology(profile.atoms)?;
    let bytes = encode_xtc(&topology, profile.frames)?;
    let writer_start = Instant::now();
    for _ in 0..profile.passes {
        black_box(encode_xtc(&topology, profile.frames)?);
    }
    let writer_elapsed = writer_start.elapsed();
    benchmark_readers(
        "XTC-0.001nm",
        profile,
        &topology,
        &bytes,
        writer_elapsed,
        || {
            XtcReader::new(
                Cursor::new(bytes.as_slice()),
                binding(&topology),
                XtcReadOptions::default(),
                TrajectoryIoLimits::default(),
                "bench.xtc",
            )
            .expect("benchmark XTC reader")
        },
        || {
            XtcReader::new(
                Cursor::new(bytes.as_slice()),
                binding(&topology),
                XtcReadOptions::default(),
                TrajectoryIoLimits::default(),
                "bench.xtc",
            )
            .expect("benchmark XTC reader")
            .into_indexed()
            .expect("benchmark XTC index")
        },
    )
}

#[allow(clippy::too_many_arguments)]
fn benchmark_readers<Sequential, Indexed, MakeSequential, MakeIndexed>(
    codec: &str,
    profile: Profile,
    topology: &Topology,
    bytes: &[u8],
    writer_elapsed: Duration,
    mut make_sequential: MakeSequential,
    mut make_indexed: MakeIndexed,
) -> Result<(), Box<dyn Error>>
where
    Sequential: TrajectoryReader,
    Indexed: SeekableTrajectoryReader,
    MakeSequential: FnMut() -> Sequential,
    MakeIndexed: FnMut() -> Indexed,
{
    let mut destination = FrameBuffer::new(topology.clone());
    let sequential_start = Instant::now();
    let mut decoded_frames = 0_usize;
    for _ in 0..profile.passes {
        let mut reader = make_sequential();
        while reader.read_next(&mut destination)? {
            decoded_frames += 1;
            black_box(destination.frame_view());
        }
    }
    let sequential_elapsed = sequential_start.elapsed();

    let index_start = Instant::now();
    for _ in 0..profile.passes {
        let reader = make_indexed();
        black_box(reader.frame_count());
    }
    let index_elapsed = index_start.elapsed();

    let mut indexed = make_indexed();
    let random_start = Instant::now();
    let mut state = 0x9e37_79b9_u64;
    for _ in 0..profile.random_reads {
        state = state
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1);
        let index = state % profile.frames as u64;
        indexed.read_frame(index, &mut destination)?;
        black_box(destination.frame_view());
    }
    let random_elapsed = random_start.elapsed();

    let mib = bytes.len() as f64 / (1024.0 * 1024.0);
    let sequential_seconds = sequential_elapsed.as_secs_f64();
    let writer_seconds = writer_elapsed.as_secs_f64();
    println!(
        "| {codec} | {} | {} | {} | {:.3} | {:.0} | {:.1} | {:.3} | {:.3} | {:.2} | {:.1} |",
        profile.name,
        profile.atoms,
        profile.frames,
        mib,
        decoded_frames as f64 / sequential_seconds,
        mib * profile.passes as f64 / sequential_seconds,
        index_elapsed.as_secs_f64() * 1000.0 / profile.passes as f64,
        profile.frames as f64 * std::mem::size_of::<u64>() as f64 / 1024.0,
        random_elapsed.as_secs_f64() * 1_000_000.0 / profile.random_reads as f64,
        mib * profile.passes as f64 / writer_seconds,
    );
    Ok(())
}

fn topology(atom_count: usize) -> Result<Topology, Box<dyn Error>> {
    let mut molecule = Molecule::new();
    let carbon = Element::from_symbol("C").expect("carbon is a built-in element");
    for _ in 0..atom_count {
        molecule.add_atom(Atom::new(carbon))?;
    }
    let molecule = SmallMolecule::from_graph(molecule);
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_small_molecule_definition(&molecule)?;
    builder.add_instance(definition, MoleculeInstanceMetadata::default())?;
    Ok(builder.build()?)
}

fn binding(topology: &Topology) -> TrajectoryTopologyBinding {
    TrajectoryTopologyBinding::new(
        topology.clone(),
        AtomOrderAssertion::assert_file_uses_topology_order(topology),
    )
    .expect("benchmark binding")
}

fn positions(topology: &Topology) -> Vec<Point3> {
    (0..topology.atom_count())
        .map(|index| {
            let index = index as f64;
            Point3::new(
                (index % 97.0) * 0.001,
                (index % 89.0) * 0.001,
                (index % 83.0) * 0.001,
            )
        })
        .collect()
}

fn base_frame(topology: &Topology) -> Result<FrameBuffer, Box<dyn Error>> {
    let mut frame = FrameBuffer::new(topology.clone());
    frame.set_positions(Quantity::new(positions(topology), NANOMETER))?;
    Ok(frame)
}

fn dynamic_frame(topology: &Topology) -> Result<FrameBuffer, Box<dyn Error>> {
    let mut frame = base_frame(topology)?;
    frame.set_cell(Some(PeriodicCell::new(
        Quantity::new(
            [
                Vector3::new(5.0, 0.0, 0.0),
                Vector3::new(0.2, 5.1, 0.0),
                Vector3::new(0.1, 0.3, 5.2),
            ],
            NANOMETER,
        ),
        [true; 3],
    )?));
    frame.set_step(Some(0));
    frame.set_time(Some(Quantity::new(0.0, PICOSECOND)))?;
    Ok(frame)
}

fn encode_xyz(topology: &Topology, frames: usize) -> Result<Vec<u8>, Box<dyn Error>> {
    let mut writer = XyzWriter::new(
        Cursor::new(Vec::new()),
        topology.clone(),
        XyzWriteOptions::default(),
        "bench.xyz",
    )?;
    let frame = base_frame(topology)?;
    for _ in 0..frames {
        writer.write_frame(frame.frame_view())?;
    }
    Ok(writer.finish()?.into_inner())
}

fn encode_dcd(topology: &Topology, frames: usize) -> Result<Vec<u8>, Box<dyn Error>> {
    let mut writer = DcdWriter::new(
        Cursor::new(Vec::new()),
        topology.clone(),
        DcdWriteOptions::default().with_step_sequence(0, 1),
        "bench.dcd",
    )?;
    let mut frame = base_frame(topology)?;
    for index in 0..frames {
        frame.set_step(Some(index as u64));
        writer.write_frame(frame.frame_view())?;
    }
    Ok(writer.finish()?.into_inner())
}

fn encode_trr(topology: &Topology, frames: usize) -> Result<Vec<u8>, Box<dyn Error>> {
    let mut writer = TrrWriter::new(
        Cursor::new(Vec::new()),
        topology.clone(),
        TrrWriteOptions::default().with_lambda_policy(TrrLambdaPolicy::RequireZero),
        "bench.trr",
    )?;
    let mut frame = dynamic_frame(topology)?;
    for index in 0..frames {
        frame.set_step(Some(index as u64));
        frame.set_time(Some(Quantity::new(index as f64 * 0.002, PICOSECOND)))?;
        writer.write_frame(frame.frame_view())?;
    }
    Ok(writer.finish()?.into_inner())
}

fn encode_xtc(topology: &Topology, frames: usize) -> Result<Vec<u8>, Box<dyn Error>> {
    let mut writer = XtcWriter::new(
        Cursor::new(Vec::new()),
        topology.clone(),
        XtcWriteOptions::default(),
        "bench.xtc",
    )?;
    let mut frame = dynamic_frame(topology)?;
    for index in 0..frames {
        frame.set_step(Some(index as u64));
        frame.set_time(Some(Quantity::new(index as f64 * 0.002, PICOSECOND)))?;
        writer.write_frame(frame.frame_view())?;
    }
    Ok(writer.finish()?.into_inner())
}
