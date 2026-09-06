//! Informational access and superposition timings using an external mmCIF model.
//! Run: cargo bench -p kekule-traj --bench frame_workflow -- SYSTEM.cif
//! Real coordinates and annotations are repeated to isolate container costs;
//! this measures API overhead, not physical trajectory sampling or scaling limits.

use kekule::{mmcif, topology::AtomSelection};
use kekule_traj::{Trajectory, TrajectoryFrame};
use std::{error::Error, hint::black_box, time::Instant};

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args_os()
        .nth(1)
        .ok_or("supply an external mmCIF model")?;
    let model = mmcif::parse_str(&std::fs::read_to_string(path)?)?
        .interpret()?
        .to_model();
    let mut frame = TrajectoryFrame::new(model.positions().clone());
    frame.set_properties(model.properties().clone())?;
    let trajectory =
        Trajectory::from_frames(model.shared_topology(), (0..64).map(|_| frame.clone()))?;
    let atoms = AtomSelection::all(&trajectory.shared_topology());
    let iterations = 1_000_000;
    let start = Instant::now();
    for index in 0..iterations {
        black_box(
            black_box(&trajectory)
                .frame(index % trajectory.len())
                .unwrap()
                .positions()
                .len(),
        );
    }
    let access_ns = start.elapsed().as_nanos() as f64 / iterations as f64;
    let start = Instant::now();
    for _ in 0..iterations / 64 {
        for frame in black_box(&trajectory).frames() {
            black_box(frame.positions().len());
        }
    }
    let iteration_ns = start.elapsed().as_nanos() as f64 / iterations as f64;
    let start = Instant::now();
    for _ in 0..20 {
        let _ = black_box(trajectory.superpose_to_frame(0, &atoms)?);
    }
    let superposition_ms = start.elapsed().as_secs_f64() * 1000.0 / 20.0;
    println!("{{\"atoms\":{},\"frames\":{},\"frame_access_ns\":{access_ns},\"iteration_ns_per_frame\":{iteration_ns},\"superposition_ms\":{superposition_ms}}}", trajectory.topology().atom_count(), trajectory.len());
    Ok(())
}
