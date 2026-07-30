//! Writes a small canonical DCD file for independent-reader validation.

use std::error::Error;
use std::path::PathBuf;

use kekule::core::{Atom, Element, Molecule};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::small::SmallMolecule;
use kekule::topology::{MoleculeInstanceMetadata, TopologyBuilder};
use kekule::trajectory::{FrameBuffer, TrajectoryFormat, TrajectoryWriter};
use kekule::units::{Quantity, ANGSTROM};
use kekule_trajectory_io::dcd::DcdWriteOptions;
use kekule_trajectory_io::{create_trajectory_writer, OverwritePolicy, TrajectoryWriteOptions};

fn main() -> Result<(), Box<dyn Error>> {
    let output = std::env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .ok_or("usage: dcd_interop_writer OUTPUT.dcd")?;
    let mut molecule = Molecule::new();
    for symbol in ["C", "H", "O"] {
        molecule.add_atom(Atom::new(
            Element::from_symbol(symbol).ok_or("unknown element")?,
        ))?;
    }
    let molecule = SmallMolecule::from_graph(molecule);
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_small_molecule_definition(&molecule)?;
    builder.add_instance(definition, MoleculeInstanceMetadata::default())?;
    let topology = builder.build()?;
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Dcd)
        .with_overwrite_policy(OverwritePolicy::Replace)
        .with_dcd_options(
            DcdWriteOptions::default()
                .with_cells(true)
                .with_step_sequence(0, 1),
        );
    let mut writer = create_trajectory_writer(&output, topology.clone(), options)?;
    let cell = PeriodicCell::new(
        Quantity::new(
            [
                Vector3::new(10.0, 0.0, 0.0),
                Vector3::new(1.0, 11.0, 0.0),
                Vector3::new(2.0, 3.0, 12.0),
            ],
            ANGSTROM,
        ),
        [true; 3],
    )?;
    let mut frame = FrameBuffer::new(topology);
    for (step, shift) in [(0, 0.0), (1, 1.0)] {
        frame.set_positions(Quantity::new(
            [
                Point3::new(0.0 + shift, 1.0, 2.0),
                Point3::new(3.0 + shift, 4.0, 5.0),
                Point3::new(6.0 + shift, 7.0, 8.0),
            ],
            ANGSTROM,
        ))?;
        frame.set_cell(Some(cell));
        frame.set_step(Some(step));
        writer.write_frame(frame.frame_view())?;
    }
    writer.finish()?;
    Ok(())
}
