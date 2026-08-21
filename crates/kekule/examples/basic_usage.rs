use std::error::Error;

use kekule::{core::Molecule, stereo};

fn main() -> Result<(), Box<dyn Error>> {
    // Interpret a chiral amino acid canonically, then perceive chemistry.
    let mut molecules = Molecule::from_smiles("C[C@@H](C(=O)O)N")?;
    let mut molecule = molecules.pop().expect("SMILES contains one component");
    molecule.perceive()?;

    // Assign absolute CIP descriptors to the perceived stereo elements.
    let stereochemistry = stereo::assign_cip_descriptors(&mut molecule)?;

    // Inspect basic graph properties and the asserted molecular charge.
    println!("atoms: {}", molecule.atom_count());
    println!("bonds: {}", molecule.bond_count());
    println!("formal charge: {}", molecule.formal_charge());
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
