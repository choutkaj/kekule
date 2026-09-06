//! Load a trajectory using an mmCIF topology and align it to its first frame.
//!
//! Run with `SYSTEM.cif TRAJECTORY [--use-stored-coordinates]`.
//! The mmCIF must contain one structural block/model, and the resulting topology
//! must match the trajectory's atom order. The optional flag allows periodic
//! frames to be fitted as stored; it performs no imaging or unwrapping.

use std::{error::Error, fs};

use kekule::{
    alignment::PeriodicAlignmentPolicy,
    mmcif,
    topology::AtomSelection,
    units::{ANGSTROM, PICOSECOND},
};
use kekule_traj::{analysis::SuperpositionOptions, io::read_trajectory};

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args_os().skip(1);
    let usage = "usage: trajectory_workflow SYSTEM.cif TRAJECTORY [--use-stored-coordinates]";
    let topology_path = args.next().ok_or(usage)?;
    let trajectory_path = args.next().ok_or(usage)?;
    let periodic_policy = match args.next() {
        None => PeriodicAlignmentPolicy::RejectPeriodic,
        Some(flag) if flag == "--use-stored-coordinates" => {
            PeriodicAlignmentPolicy::UseStoredCoordinates
        }
        Some(_) => return Err(usage.into()),
    };
    if args.next().is_some() {
        return Err(usage.into());
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

    // Fit all atoms; substitute a protein/backbone selection for a solvated system.
    let fit = AtomSelection::from_atoms(&topology, topology.atom_ids().iter().copied())?;
    let report = trajectory.superpose_to_frame_with_options(
        0,
        &fit,
        SuperpositionOptions {
            periodic_policy,
            ..Default::default()
        },
    )?;

    for (index, (frame, alignment)) in trajectory.frames().zip(report.alignments()).enumerate() {
        let time_ps = frame
            .time()
            .map(|time| time.value_in(PICOSECOND))
            .transpose()?;
        println!(
            "Frame {index}: time_ps={time_ps:?}, step={:?}, fitted RMSD={:.3} A",
            frame.step(),
            alignment.rmsd().value_in(ANGSTROM)?
        );
    }
    Ok(())
}
