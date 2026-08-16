use std::error::Error;

use kekule::{small::SmallMolecule, stereo};

fn main() -> Result<(), Box<dyn Error>> {
    // Parse a chiral amino acid, normalize its representation, and perceive chemistry.
    let mut molecule = SmallMolecule::from_smiles("C[C@@H](C(=O)O)N")?;
    molecule.normalize()?;
    molecule.perceive()?;

    // Assign absolute CIP descriptors to the perceived stereo elements.
    let stereochemistry = stereo::assign_cip_descriptors(molecule.graph_mut())?;

    // Inspect basic graph properties and the asserted molecular charge.
    println!("atoms: {}", molecule.atom_count());
    println!("bonds: {}", molecule.bond_count());
    println!("formal charge: {}", molecule.graph().formal_charge());
    for assignment in &stereochemistry.assigned {
        println!(
            "stereo {:?}: {:?}",
            assignment.element, assignment.descriptor
        );
    }

    // Write canonical connectivity and a stereo-preserving SMILES form.
    println!("canonical SMILES: {}", molecule.to_canonical_smiles()?);
    println!("isomeric SMILES: {}", molecule.to_isomeric_smiles()?);
    Ok(())
}
