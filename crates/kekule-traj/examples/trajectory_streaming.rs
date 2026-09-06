//! Align a TRR trajectory to its first frame with bounded frame storage.
//!
//! Run with `SYSTEM.cif INPUT.trr OUTPUT.trr`. Atom order must match the mmCIF.
//! Coordinates are fitted as stored; perform periodic preprocessing first when
//! needed. Explicit f64 output avoids narrowing an f64 TRR source, and the default
//! TRR lambda policy preserves its frame property. Existing output is protected.

use std::{error::Error, fs};

use kekule::{mmcif, topology::AtomSelection};
use kekule_traj::{
    analysis::FrameSuperposer,
    io::{
        create_trajectory_writer, open_trajectory,
        trr::{TrrScalarPrecision, TrrWriteOptions},
        TrajectoryWriteOptions,
    },
    TrajectoryFormat, TrajectoryReader, TrajectoryWriter,
};

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args_os().skip(1);
    let usage = "usage: trajectory_streaming SYSTEM.cif INPUT.trr OUTPUT.trr";
    let topology_path = args.next().ok_or(usage)?;
    let input = args.next().ok_or(usage)?;
    let output = args.next().ok_or(usage)?;
    if args.next().is_some() {
        return Err(usage.into());
    }

    let topology = mmcif::parse_str(&fs::read_to_string(topology_path)?)?
        .interpret()?
        .to_topology();
    let mut reader = open_trajectory(input, topology.clone())?;
    let mut buffer = reader.frame_buffer();
    if !reader.read_next(&mut buffer)? {
        return Err("input has no frames".into());
    }

    // Retain one complete reference before the reusable input buffer advances.
    let reference = buffer.frame_view().to_frame();
    let atoms = AtomSelection::all(&topology);
    let superposer = FrameSuperposer::new(reference.view(&topology)?, &atoms);
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Trr)
        .with_trr_options(TrrWriteOptions::default().with_precision(TrrScalarPrecision::Float64));
    let mut writer = create_trajectory_writer(output, topology.clone(), options)?;
    let mut index = 0;
    loop {
        superposer.superpose_in_place(index, &mut buffer)?;
        writer.write_frame(buffer.frame_view())?;
        index += 1;
        if !reader.read_next(&mut buffer)? {
            break;
        }
    }
    // Publication happens only after clean EOF and successful finalization.
    writer.finish()?;
    println!("Aligned and saved {index} frames");
    Ok(())
}
