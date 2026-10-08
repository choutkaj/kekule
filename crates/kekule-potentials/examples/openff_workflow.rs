//! Single-point OpenFF energy and minimization of a 3D SDF ligand.
//!
//! ```text
//! cargo run -p kekule-potentials --example openff_workflow -- ligand.sdf path/to/nagl-model
//! ```
//!
//! The first record must contain explicit hydrogens with 3D coordinates. The
//! model directory is an exported NAGL bundle compatible with OpenFF Rosemary.

use std::{env, error::Error, fs};

use kekule::sdf;
use kekule_openff::{ForceField, NaglModel};
use kekule_potentials::{minimize, openff::OpenFfPotential, MinimizeOptions, Potential};

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = env::args().skip(1);
    let (Some(path), Some(nagl)) = (args.next(), args.next()) else {
        return Err("usage: openff_workflow <ligand.sdf> <nagl-model-directory>".into());
    };

    let document = sdf::parse_str(&fs::read_to_string(path)?)?;
    let record = document
        .records()
        .first()
        .ok_or("the SDF document has no records")?;
    let model = record.to_model()?;

    // Parameterization binds the model's exact topology snapshot.
    let nagl = NaglModel::load(nagl)?;
    let parameters = ForceField::rosemary()?.parameterize(model.shared_topology(), &nagl)?;
    let potential = OpenFfPotential::new(&parameters)?;

    let energy = potential.energy(model.view())?;
    println!("initial energy: {:.4} kJ/mol", energy.total().into_value());
    for component in energy.components() {
        println!(
            "  {:<18} {:>12.4}",
            component.kind,
            component.energy.into_value()
        );
    }

    let result = minimize(&potential, model.view(), &MinimizeOptions::default())?;
    let final_energy = result.final_evaluation().energy().total().into_value();
    println!(
        "minimized energy: {final_energy:.4} kJ/mol after {} iterations ({:?})",
        result.iterations(),
        result.status()
    );
    let minimized = result.to_model(model.view())?;
    println!("minimized model has {} atoms", minimized.atom_count());
    Ok(())
}
