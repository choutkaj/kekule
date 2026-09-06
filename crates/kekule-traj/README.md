# kekule-traj

`kekule-traj` provides fixed-topology trajectory storage, streaming, analysis,
and pure-Rust XYZ, DCD, TRR, and XTC file I/O for [`kekule`](https://crates.io/crates/kekule).

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
use kekule_traj::io::{read_trajectory, write_trajectory};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let document = mmcif::parse_str(&std::fs::read_to_string("system.cif")?)?;
    let topology = document.interpret()?.to_topology();
    let trajectory = read_trajectory("trajectory.xtc", topology.clone())?;

    println!("Frames: {}", trajectory.len());
    println!("Atoms: {}", topology.atom_count());
    println!("Residues: {}", topology.residues().count());

    // If needed, reconstruct molecules split across periodic boundaries first:
    // let trajectory = trajectory.make_molecules_whole()?;

    // Fit all atoms, or use a protein/backbone selection for a solvated system.
    let fit = AtomSelection::all(&topology);
    let aligned = trajectory.superpose_to_frame(0, &fit)?;
    let rmsd = aligned.rmsd_to_frame(0, &fit)?.value_in(ANGSTROM)?;
    for (index, value) in rmsd.iter().enumerate() {
        println!("Frame {index}: fitted RMSD = {value:.3} A");
    }

    // Frame selection preserves the original times and simulation steps.
    let sampled = aligned.select_frames((0..aligned.len()).step_by(10))?;
    write_trajectory("aligned.xtc", &sampled)?;
    Ok(())
}
```

Alignment requires a nonempty trajectory and a non-collinear fitting selection.
Transformations return a new trajectory; use their `_in_place` variants to modify
the original. Imaging and temporal unwrapping are explicit preprocessing steps;
unwrap before discarding intermediate frames.

Saving infers the format from the extension, protects existing files, and rejects
unsupported metadata. Use `read_trajectory_with_options` or
`write_trajectory_with_options` for explicit codec policies and precision.

For more, see the [loaded workflow](examples/trajectory_workflow.rs) and
[streaming workflow](examples/trajectory_streaming.rs) examples. Streaming uses
`open_trajectory` and a reusable frame buffer to process large files one frame at
a time.

Licensed under MIT or Apache-2.0.
