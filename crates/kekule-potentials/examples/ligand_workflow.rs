use std::{error::Error, fs};

use kekule::{
    core::Conformer,
    geometry::Point3,
    modeling::{minimize, MinimizeOptions},
    sdf::{self, SdfParseOptions, SdfRecordInterpretation},
    structure::Model,
    units::{Quantity, ANGSTROM, MODEL_GRADIENT_UNIT},
};
use kekule_potentials::dreiding::{DreidingPotential, DreidingPrepareOptions};

fn main() -> Result<(), Box<dyn Error>> {
    // Parse and canonically interpret one SDF record without perceiving it.
    let input = fs::read_to_string("examples/ligand.sdf")?;
    let document = sdf::parse_str(&input, SdfParseOptions::default())?;
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
    let mut conformer = Conformer::new(ANGSTROM)?;
    for atom in ligand.atom_ids() {
        conformer.set_position(
            atom,
            Quantity::new(Point3::new(atom.index() as f64, 0.0, 0.0), ANGSTROM),
        )?;
    }
    let mut builder = Model::builder();
    builder.add_molecule(&ligand, &conformer)?;
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

    // Reassemble the original record metadata. Molecule output is geometry-independent.
    let output = sdf::write_v2000(&[SdfRecordInterpretation::new(
        title,
        vec![ligand],
        data_fields,
    )])?;
    fs::write("examples/ligand-minimized.sdf", output)?;
    Ok(())
}
