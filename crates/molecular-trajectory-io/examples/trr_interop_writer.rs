//! Writes a small full-field TRR file for independent-reader validation.

use std::error::Error;
use std::path::PathBuf;

use molecular::core::{Atom, Element, Molecule, PropValue};
use molecular::geometry::{PeriodicCell, Point3, Vector3};
use molecular::small::SmallMolecule;
use molecular::topology::{MoleculeInstanceMetadata, TopologyBuilder};
use molecular::trajectory::{FrameBuffer, TrajectoryFormat, TrajectoryWriter};
use molecular::units::{Quantity, MODEL_FORCE_UNIT, MODEL_VELOCITY_UNIT, NANOMETER, PICOSECOND};
use molecular_trajectory_io::trr::{TrrScalarPrecision, TrrWriteOptions};
use molecular_trajectory_io::{create_trajectory_writer, OverwritePolicy, TrajectoryWriteOptions};

fn main() -> Result<(), Box<dyn Error>> {
    let output = std::env::args_os()
        .nth(1)
        .map(PathBuf::from)
        .ok_or("usage: trr_interop_writer OUTPUT.trr")?;
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
    let options = TrajectoryWriteOptions::new(TrajectoryFormat::Trr)
        .with_overwrite_policy(OverwritePolicy::Replace)
        .with_trr_options(TrrWriteOptions::default().with_precision(TrrScalarPrecision::Float32));
    let mut writer = create_trajectory_writer(&output, topology.clone(), options)?;
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
            MODEL_VELOCITY_UNIT,
        )))?;
        frame.set_forces(Some(Quantity::new(
            [
                Vector3::new(10.0, 20.0, 30.0),
                Vector3::new(40.0, 50.0, 60.0),
                Vector3::new(70.0, 80.0, 90.0),
            ],
            MODEL_FORCE_UNIT,
        )))?;
        frame.set_time(Some(Quantity::new(step as f64 * 0.25, PICOSECOND)))?;
        frame.set_step(Some(step));
        frame.props_mut().clear();
        frame.props_mut().insert(
            "gromacs.trr.lambda".into(),
            PropValue::Float(0.125 + step as f64 * 0.125),
        );
        writer.write_frame(frame.frame_view())?;
    }
    writer.finish()?;
    Ok(())
}
