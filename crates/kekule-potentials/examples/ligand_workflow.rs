use std::{error::Error, fs};

use kekule::{
    geometry::Point3,
    modeling::{minimize, MinimizeOptions},
    sdf::{self, SdfRecordInterpretation},
    structure::{Model, Positions},
    units::{Quantity, ANGSTROM, KILOJOULE_PER_MOLE_PER_ANGSTROM},
};
use kekule_potentials::dreiding::{DreidingPotential, DreidingPrepareOptions};

fn main() -> Result<(), Box<dyn Error>> {
    // Parse and canonically interpret one SDF record without perceiving it.
    let input = fs::read_to_string("examples/ligand.sdf")?;
    let document = sdf::parse_str(&input)?;
    let mut records = sdf::interpret(&document)?.to_records();
    assert_eq!(records.len(), 1, "expected one ligand record");

    // Preserve the record metadata while working on its molecule.
    let record = records.pop().expect("record count was checked");
    let title = record.title().to_owned();
    let data_fields = record.data_fields().to_vec();
    let mut molecules = record.to_molecules();
    assert_eq!(
        molecules.len(),
        1,
        "expected one connected ligand component"
    );
    let mut ligand = molecules.pop().expect("component count was checked");
    ligand.perceive()?;

    // Inspect the canonical, perceived ligand before modeling it.
    println!("atoms: {}", ligand.atom_count());
    println!("bonds: {}", ligand.bond_count());
    println!("formal charge: {}", ligand.formal_charge());

    // Geometry is detached from Molecule. Supply it explicitly when building a model.
    let positions = Positions::new(Quantity::new(
        ligand
            .atom_ids()
            .map(|atom| Point3::new(atom.index() as f64, 0.0, 0.0))
            .collect::<Vec<_>>(),
        ANGSTROM,
    ))?;
    let mut builder = Model::builder();
    builder.add_molecule(&ligand, &positions)?;
    let model = builder.build()?;

    // DREIDING support is explicitly nonperiodic. This SDF-derived model has no
    // periodic cell, so it is eligible for preparation and minimization.
    assert!(model.cell().is_none());
    let mut potential = DreidingPotential::prepare(
        &model.shared_topology(),
        model.view(),
        DreidingPrepareOptions::default(),
    )?;
    let minimized = minimize(
        &model,
        &mut potential,
        MinimizeOptions {
            max_iterations: 10_000,
            gradient_tolerance: 0.05 * KILOJOULE_PER_MOLE_PER_ANGSTROM,
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

    // Reassemble the original record metadata around the minimized geometry-rich model.
    let output = sdf::write_v2000(&[SdfRecordInterpretation::new(
        title,
        minimized.model,
        data_fields,
    )])?;
    fs::write("examples/ligand-minimized.sdf", output)?;
    Ok(())
}
