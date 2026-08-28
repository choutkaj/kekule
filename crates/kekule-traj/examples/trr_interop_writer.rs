//! Writes a small full-field TRR file for independent-reader validation.

use std::error::Error;
use std::path::PathBuf;
use std::sync::Arc;

use kekule::core::{Atom, BondOrder, Element};
use kekule::geometry::{PeriodicCell, Point3, Vector3};
use kekule::properties::{PropertyKey, PropertyValue};
use kekule::topology::TopologyBuilder;
use kekule::units::{
    Quantity, CANONICAL_FORCE_UNIT, CANONICAL_VELOCITY_UNIT, DIMENSIONLESS, NANOMETER, PICOSECOND,
};
use kekule_traj::io::trr::{TrrScalarPrecision, TrrWriteOptions};
use kekule_traj::io::{create_trajectory_writer, OverwritePolicy, TrajectoryWriteOptions};
use kekule_traj::{FrameBuffer, TrajectoryFormat, TrajectoryWriter};

fn main() -> Result<(), Box<dyn Error>> {
    let output = std::env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .ok_or("usage: trr_interop_writer OUTPUT.trr")?;
    let mut molecule = kekule::core::MoleculeEditor::new();
    let mut atoms = Vec::new();
    for symbol in ["C", "H", "O"] {
        atoms.push(molecule.add_atom(Atom::new(
            Element::from_symbol(symbol).ok_or("unknown element")?,
        ))?);
    }
    molecule.add_bond(atoms[0], atoms[1], BondOrder::Single)?;
    molecule.add_bond(atoms[0], atoms[2], BondOrder::Single)?;
    let molecule = molecule.finish()?;
    let mut builder = TopologyBuilder::new();
    let definition = builder.add_molecule_definition(&molecule)?;
    builder.add_instance(definition)?;
    let topology = Arc::new(builder.build()?);
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Trr)
        .with_overwrite_policy(OverwritePolicy::Replace)
        .with_trr_options(TrrWriteOptions::default().with_precision(TrrScalarPrecision::Float32));
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
    for (step, shift) in [(0, 0.0), (1, 0.1)] {
        frame.set_positions(Quantity::new(
            [
                Point3::new(0.0 + shift, 0.1, 0.2),
                Point3::new(0.3 + shift, 0.4, 0.5),
                Point3::new(0.6 + shift, 0.7, 0.8),
            ],
            NANOMETER,
        ))?;
        frame.set_cell(Some(cell));
        frame.set_velocities(Some(Quantity::new(
            [
                Vector3::new(1.0, 2.0, 3.0),
                Vector3::new(4.0, 5.0, 6.0),
                Vector3::new(7.0, 8.0, 9.0),
            ],
            CANONICAL_VELOCITY_UNIT,
        )))?;
        frame.set_forces(Some(Quantity::new(
            [
                Vector3::new(10.0, 20.0, 30.0),
                Vector3::new(40.0, 50.0, 60.0),
                Vector3::new(70.0, 80.0, 90.0),
            ],
            CANONICAL_FORCE_UNIT,
        )))?;
        frame.set_time(Some(Quantity::new(step as f64 * 0.25, PICOSECOND)))?;
        frame.set_step(Some(step));
        frame.properties_mut().clear_owner();
        frame.properties_mut().insert(
            PropertyKey::new("gromacs.trr.lambda")?,
            PropertyValue::Real {
                value: 0.125 + step as f64 * 0.125,
                unit: DIMENSIONLESS,
            },
        )?;
        writer.write_frame(frame.frame_view())?;
    }
    writer.finish()?;
    Ok(())
}
