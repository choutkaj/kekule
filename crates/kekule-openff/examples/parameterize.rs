use kekule::{hydrogens, smiles};
use kekule_openff::{ForceField, NaglModel};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = std::env::args().skip(1);
    let directory = args
        .next()
        .ok_or("usage: parameterize MODEL_DIRECTORY [SMILES]")?;
    let input = args.next().unwrap_or_else(|| "CCO".into());
    let mut molecules = smiles::to_molecules(&input)?;
    if molecules.len() != 1 {
        return Err("supply one connected molecule".into());
    }
    let mut molecule = molecules.remove(0);
    molecule.perceive()?;
    hydrogens::add_hydrogens(&mut molecule)?;
    let model = NaglModel::load(directory)?;
    let p = ForceField::rosemary()?.parameterize_molecule(&molecule, &model)?;
    println!(
        "{} atoms, {} bonds, {} angles, {} proper torsions, {} improper terms",
        molecule.atom_count(),
        p.bonds().len(),
        p.angles().len(),
        p.proper_torsions().len(),
        p.improper_torsions().len()
    );
    println!("Charge source: {:?}", p.charge_sources());
    println!("Charges (e): {:?}", p.charges().value());
    Ok(())
}
