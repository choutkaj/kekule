//! Single-point OpenFF energy and minimization of a 3D SDF ligand.
//!
//! ```text
//! cargo run -p kekule-potentials --features openff --example openff_workflow -- ligand.sdf
//! ```
//!
//! The first record must contain explicit hydrogens with 3D coordinates. Charges
//! use the bundled Ash model; an optional second argument names an exported NAGL
//! bundle compatible with OpenFF Rosemary instead.

use std::{env, error::Error, fs};

use kekule::sdf;
use kekule_openff::{ForceField, NaglModel};
use kekule_potentials::{minimize, openff::OpenFfPotential, MinimizeOptions, Potential};

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = env::args().skip(1);
    let Some(path) = args.next() else {
        return Err("usage: openff_workflow <ligand.sdf> [nagl-model-directory]".into());
    };
    let nagl = match args.next() {
        Some(directory) => NaglModel::load(directory)?,
        None => NaglModel::ash()?,
    };

    let document = sdf::parse_str(&fs::read_to_string(path)?)?;
    let record = document
        .records()
        .first()
        .ok_or("the SDF document has no records")?;
    let model = record.to_model()?;

    // Parameterization binds the model's exact topology snapshot.
    let parameters = ForceField::rosemary()?.parameterize(model.shared_topology(), &nagl)?;
    let potential = OpenFfPotential::new(&parameters)?;

    let energy = potential.energy(model.as_model_view())?;
    println!("initial energy: {:.4} kJ/mol", energy.total().into_value());
    for component in energy.components() {
        println!(
            "  {:<18} {:>12.4}",
            component.kind,
            component.energy.into_value()
        );
    }

    let result = minimize(
        &potential,
        model.as_model_view(),
        &MinimizeOptions::default(),
    )?;
    let final_energy = result.final_evaluation().energy().total().into_value();
    println!(
        "minimized energy: {final_energy:.4} kJ/mol after {} iterations ({:?})",
        result.iterations(),
        result.status()
    );
    let minimized = result.to_model(model.as_model_view())?;
    println!("minimized model has {} atoms", minimized.atom_count());
    Ok(())
}
