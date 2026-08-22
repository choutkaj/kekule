//! Writes a compressed XTC file for independent-reader validation.

use std::error::Error;
use std::path::PathBuf;
use std::sync::Arc;

use kekule::core::{Atom, BondOrder, Element};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::topology::{MoleculeInstanceMetadata, TopologyBuilder};
use kekule::units::{Quantity, NANOMETER, PICOSECOND};
use kekule_traj::io::xtc::XtcWriteOptions;
use kekule_traj::io::{create_trajectory_writer, OverwritePolicy, TrajectoryWriteOptions};
use kekule_traj::{FrameBuffer, TrajectoryFormat, TrajectoryWriter};

fn main() -> Result<(), Box<dyn Error>> {
    let output = std::env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .ok_or("usage: xtc_interop_writer OUTPUT.xtc")?;
    let mut molecule = kekule::core::MoleculeEditor::new();
    let mut previous = None;
    for _ in 0..12 {
        let atom = molecule.add_atom(Atom::new(
            Element::from_symbol("C").ok_or("unknown element")?,
        ))?;
        if let Some(parent) = previous {
            molecule.add_bond(parent, atom, BondOrder::Single)?;
        }
        previous = Some(atom);
    }
    let molecule = molecule.finish()?;
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule)?;
    builder.add_instance(definition, MoleculeInstanceMetadata::default())?;
    let topology = Arc::new(builder.build()?);
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xtc)
        .with_overwrite_policy(OverwritePolicy::Replace)
        .with_xtc_options(XtcWriteOptions::default().with_precision(1000.0)?);
    let mut writer = create_trajectory_writer(&output, Arc::clone(&topology), options)?;
    let cell = PeriodicCell::new(
        Quantity::new(
            [
                Vector3::new(2.0, 0.0, 0.0),
                Vector3::new(0.1, 2.1, 0.0),
                Vector3::new(0.2, 0.3, 2.2),
            ],
            NANOMETER,
        ),
        [true; 3],
    )?;
    let mut frame = FrameBuffer::new(topology);
    for (step, shift) in [(0, 0.0), (1, 0.01)] {
        let positions = (0..12)
            .map(|index| {
                let index = index as f64;
                Point3::new(
                    0.1 * index + shift,
                    0.2 * index + shift,
                    0.3 * index + shift,
                )
            })
            .collect::<Vec<_>>();
        frame.set_positions(Quantity::new(positions, NANOMETER))?;
        frame.set_cell(Some(cell));
        frame.set_time(Some(Quantity::new(step as f64 * 0.25, PICOSECOND)))?;
        frame.set_step(Some(step));
        writer.write_frame(frame.frame_view())?;
    }
    writer.finish()?;
    Ok(())
}
