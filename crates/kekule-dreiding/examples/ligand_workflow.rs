use std::{error::Error, fs};

use kekule::{
    modeling::{minimize, MinimizeOptions},
    sdf::{self, SdfParseOptions, SdfRecord},
    structure::Model,
    units::MODEL_GRADIENT_UNIT,
};
use kekule_dreiding::{DreidingPotential, DreidingPrepareOptions};

fn main() -> Result<(), Box<dyn Error>> {
    // Parse and interpret one SDF record without silently sanitizing it.
    let input = fs::read_to_string("examples/ligand.sdf")?;
    let document = sdf::parse_str(&input, SdfParseOptions::default())?;
    let mut records = sdf::interpret(&document)?.into_records();
    assert_eq!(records.len(), 1, "expected one ligand record");

    // Preserve the record metadata while working on its molecule.
    let record = records.pop().expect("record count was checked");
    let title = record.title().to_owned();
    let data_fields = record.data_fields().to_vec();
    let mut ligand = record.into_molecule();
    ligand.sanitize()?;

    // Inspect the sanitized ligand before modeling it.
    println!("atoms: {}", ligand.atom_count());
    println!("bonds: {}", ligand.bond_count());
    println!("formal charge: {}", ligand.graph().formal_charge());

    // Build a fixed-topology model from the ligand's first conformer.
    let conformer = ligand
        .graph()
        .first_conformer()
        .map(|(id, _)| id)
        .expect("the SDF record has 3D coordinates");
    let mut builder = Model::builder();
    let instance = builder.add_small_molecule(&ligand, conformer)?;
    let model = builder.build()?;

    // DREIDING support is explicitly nonperiodic. This SDF-derived model has no
    // periodic cell, so it is eligible for preparation and minimization.
    assert!(model.cell().is_none());
    let mut potential = DreidingPotential::prepare(
        model.topology(),
        model.view(),
        DreidingPrepareOptions::default(),
    )?;
    let minimized = minimize(
        &model,
        &mut potential,
        MinimizeOptions {
            max_iterations: 10_000,
            gradient_tolerance: 0.05 * MODEL_GRADIENT_UNIT,
            ..MinimizeOptions::default()
        },
    )?;
    println!(
        "{:?} after {} iterations: {} -> {} {}",
        minimized.status,
        minimized.iterations,
        minimized.initial_energy.value(),
        minimized.final_energy.value(),
        minimized.final_energy.unit()
    );

    // Copy the optimized instance positions back to the source conformer.
    minimized
        .model
        .instance_to_conformer(instance, ligand.graph_mut(), conformer)?;

    // Reassemble the original record metadata and write the optimized SDF.
    let output = sdf::write_v2000(&[SdfRecord::new(title, ligand, data_fields)])?;
    fs::write("examples/ligand-minimized.sdf", output)?;
    Ok(())
}
