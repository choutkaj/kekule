# kekule-traj

`kekule-traj` provides trajectory streaming, periodic preprocessing, analysis,
and pure-Rust XYZ, DCD, TRR, and XTC file I/O for [`kekule`](https://crates.io/crates/kekule).
The in-memory `Trajectory` type itself, including superposition and RMSD, lives in
`kekule::structure`.

## Installation

```sh
cargo add kekule kekule-traj
```

## Basic example

Load topology from an mmCIF file, read a trajectory, align all frames to the first,
and save every tenth frame. The mmCIF must contain one structural block/model and
match the trajectory's atom order.

```rust
use kekule::{mmcif, topology::AtomSelection, units::ANGSTROM};
use kekule_traj::{
    io::{read_trajectory, write_trajectory},
    periodic,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let document = mmcif::parse_str(&std::fs::read_to_string("system.cif")?)?;
    let topology = document.interpret()?.into_topology();
    let mut trajectory = read_trajectory("trajectory.xtc", topology.clone())?;

    println!("Frames: {}", trajectory.len());
    println!("Atoms: {}", topology.atom_count());
    println!("Residues: {}", topology.residues().count());

    // Reconstruct molecules split across periodic boundaries first.
    periodic::make_molecules_whole(&mut trajectory)?;

    // Fit all atoms onto frame 0, or use a protein/backbone selection for a
    // solvated system. Superposition changes the trajectory in place.
    let fit = AtomSelection::all(&topology);
    trajectory.superpose(0, &fit)?;
    let rmsd = trajectory.rmsd(0, &fit)?.value_in(ANGSTROM)?;
    for (index, value) in rmsd.iter().enumerate() {
        println!("Frame {index}: fitted RMSD = {value:.3} A");
    }

    // Frame selection preserves the original times and simulation steps.
    let sampled = trajectory.select((0..trajectory.len()).step_by(10))?;
    write_trajectory("aligned.xtc", &sampled)?;
    Ok(())
}
```

Alignment requires a nonempty trajectory and a non-collinear fitting selection.
Superposition and periodic operations modify the trajectory in place and leave it
unchanged if they fail; clone it first to keep the original. Imaging and temporal
unwrapping are explicit preprocessing steps; unwrap before discarding
intermediate frames.

Saving infers the format from the extension, protects existing files, and rejects
unsupported metadata. Use `read_trajectory_with_options` or
`write_trajectory_with_options` for explicit codec policies and precision.

For more, see the [loaded workflow](examples/trajectory_workflow.rs) and
[streaming workflow](examples/trajectory_streaming.rs) examples. Streaming uses
`open_trajectory` and a reusable frame buffer to process large files one frame at
a time.

Licensed under MIT or Apache-2.0.
