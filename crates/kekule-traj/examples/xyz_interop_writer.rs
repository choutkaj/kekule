//! Writes a small strict XYZ file for independent-reader validation.

use std::error::Error;
use std::path::PathBuf;
use std::sync::Arc;

use kekule::core::{Atom, BondOrder, Element};
use kekule::geometry::Point3;
use kekule::topology::{MoleculeInstanceMetadata, TopologyBuilder};
use kekule::units::{Quantity, ANGSTROM};
use kekule_traj::io::xyz::XyzWriteOptions;
use kekule_traj::io::{create_trajectory_writer, OverwritePolicy, TrajectoryWriteOptions};
use kekule_traj::{FrameBuffer, TrajectoryFormat, TrajectoryWriter};

fn main() -> Result<(), Box<dyn Error>> {
    let output = std::env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .ok_or("usage: xyz_interop_writer OUTPUT.xyz")?;
    let mut molecule = kekule::core::MoleculeEditor::new();
    let mut atoms = Vec::new();
    for symbol in ["O", "H", "H"] {
        atoms.push(molecule.add_atom(Atom::new(
            Element::from_symbol(symbol).ok_or("unknown element")?,
        ))?);
    }
    molecule.add_bond(atoms[0], atoms[1], BondOrder::Single)?;
    molecule.add_bond(atoms[0], atoms[2], BondOrder::Single)?;
    let molecule = molecule.finish()?;
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule)?;
    builder.add_instance(definition, MoleculeInstanceMetadata::default())?;
    let topology = Arc::new(builder.build()?);
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Xyz)
        .with_overwrite_policy(OverwritePolicy::Replace)
        .with_xyz_options(
            XyzWriteOptions::default()
                .with_decimal_places(8)
                .with_comment("Kekule XYZ interoperability fixture"),
        );
    let mut writer = create_trajectory_writer(&output, Arc::clone(&topology), options)?;
    let mut frame = FrameBuffer::new(topology);
    for shift in [0.0, 0.1] {
        frame.set_positions(Quantity::new(
            [
                Point3::new(shift, 0.0, 0.0),
                Point3::new(0.9572 + shift, 0.0, 0.0),
                Point3::new(-0.239_987_2 + shift, 0.927_297, 0.0),
            ],
            ANGSTROM,
        ))?;
        writer.write_frame(frame.frame_view())?;
    }
    writer.finish()?;
    Ok(())
}
