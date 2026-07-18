use std::{error::Error, fs};

use molecules::{
    modeling::{minimize, InstanceAtomId, MinimizeOptions, Model},
    perception,
    sdf::{self, SdfParseOptions, SdfRecord},
    smiles,
    units::MODEL_GRADIENT_UNIT,
};
use molecules_dreiding::DreidingPotential;

fn main() -> Result<(), Box<dyn Error>> {
    let input = fs::read_to_string("examples/ligand.sdf")?;
    let document = sdf::parse_str(&input, SdfParseOptions::default())?;
    let mut records = sdf::interpret(&document)?.into_records();
    assert_eq!(records.len(), 1, "expected one ligand record");

    let record = records.pop().expect("record count was checked");
    let title = record.title().to_owned();
    let data_fields = record.data_fields().to_vec();
    let mut ligand = record.into_molecule();
    perception::sanitize(&mut ligand)?;

    let charge: i32 = ligand
        .atoms()
        .map(|(_, atom)| i32::from(atom.formal_charge))
        .sum();
    println!("atoms: {}", ligand.atom_count());
    println!("bonds: {}", ligand.bond_count());
    println!("formal charge: {charge}");
    println!("canonical SMILES: {}", smiles::write_canonical(&ligand)?);

    let conformer = ligand
        .graph()
        .first_conformer()
        .map(|(id, _)| id)
        .expect("the SDF record has 3D coordinates");
    let mut builder = Model::builder();
    let instance = builder.add_small_molecule(&ligand, conformer)?;
    let model = builder.build()?;

    let mut potential = DreidingPotential::prepare(&model)?;
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
        "minimization status: {:?} after {} iterations",
        minimized.status, minimized.iterations
    );
    println!(
        "energy: {} -> {} {}",
        minimized.initial_energy.value(),
        minimized.final_energy.value(),
        minimized.final_energy.unit()
    );
    println!(
        "final max gradient: {} {}",
        minimized.final_max_gradient.value(),
        minimized.final_max_gradient.unit()
    );

    let atom_ids = ligand.graph().atom_ids().collect::<Vec<_>>();
    for atom in atom_ids {
        let position = minimized
            .model
            .position(InstanceAtomId::new(instance, atom))?;
        ligand
            .graph_mut()
            .conformer_mut(conformer)?
            .set_position(atom, position)?;
    }

    let output = sdf::write_v2000(&[SdfRecord::new(title, ligand, data_fields)])?;
    fs::write("examples/ligand-minimized.sdf", output)?;
    Ok(())
}
