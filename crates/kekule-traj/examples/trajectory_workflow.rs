//! Load a trajectory using an mmCIF topology and align it to its first frame.
//!
//! Run with `SYSTEM.cif TRAJECTORY [--make-whole] [--unwrap] [--stride N] [--output PATH]`.
//! The mmCIF must contain one structural block/model, and the resulting topology
//! must match the trajectory's atom order. Periodic preprocessing runs before
//! downsampling. Saving uses the output extension and strict codec defaults;
//! unsupported metadata is an error. Use explicit writer options in your own
//! workflow when changing precision or writing a different format profile.

use std::{error::Error, fs};

use kekule::{
    mmcif,
    topology::AtomSelection,
    units::{ANGSTROM, PICOSECOND},
};
use kekule_traj::io::{read_trajectory, write_trajectory};

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args_os().skip(1);
    let usage = "usage: trajectory_workflow SYSTEM.cif TRAJECTORY [--make-whole] [--unwrap] [--stride N] [--output PATH]";
    let topology_path = args.next().ok_or(usage)?;
    let trajectory_path = args.next().ok_or(usage)?;
    let mut make_whole = false;
    let mut unwrap = false;
    let mut stride = 1;
    let mut output = None;
    while let Some(flag) = args.next() {
        match flag.to_str() {
            Some("--make-whole") => make_whole = true,
            Some("--unwrap") => unwrap = true,
            Some("--stride") => {
                stride = args
                    .next()
                    .ok_or(usage)?
                    .to_str()
                    .ok_or(usage)?
                    .parse::<usize>()?;
                if stride == 0 {
                    return Err("stride must be positive".into());
                }
            }
            Some("--output") => output = Some(args.next().ok_or(usage)?),
            _ => return Err(usage.into()),
        }
    }

    let document = mmcif::parse_str(&fs::read_to_string(topology_path)?)?;
    let topology = document.interpret()?.to_topology();
    let mut trajectory = read_trajectory(trajectory_path, topology.clone())?;

    println!("Frames: {}", trajectory.len());
    println!("Atoms: {}", topology.atom_count());
    println!("Molecules: {}", topology.instance_count());
    println!("Chains: {}", topology.chains().count());
    println!("Residues: {}", topology.residues().count());
    println!(
        "Frames with a periodic cell: {}",
        trajectory
            .frames()
            .filter(|frame| frame.cell().is_some())
            .count()
    );
    trajectory.validate_monotonic_time(false)?;

    if make_whole {
        trajectory.make_molecules_whole_in_place()?;
    }
    if unwrap {
        trajectory.unwrap_in_place()?;
    }
    if stride != 1 {
        trajectory = trajectory.select_frames((0..trajectory.len()).step_by(stride))?;
        println!(
            "Retained {} frames with original times and steps",
            trajectory.len()
        );
    }

    // Fit all atoms; substitute a protein/backbone selection for a solvated system.
    let fit = AtomSelection::all(&topology);
    let aligned = trajectory.superpose_to_frame(0, &fit)?;
    let rmsd = aligned.rmsd_to_frame(0, &fit)?.value_in(ANGSTROM)?;

    for (index, (frame, rmsd)) in aligned.frames().zip(rmsd).enumerate() {
        let time_ps = frame
            .time()
            .map(|time| time.value_in(PICOSECOND))
            .transpose()?;
        println!(
            "Frame {index}: time_ps={time_ps:?}, step={:?}, fitted RMSD={:.3} A",
            frame.step(),
            rmsd
        );
    }
    if let Some(path) = output {
        write_trajectory(path, &aligned)?;
    }
    Ok(())
}
