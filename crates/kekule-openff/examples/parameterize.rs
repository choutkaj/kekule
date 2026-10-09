//! Parameterize one molecule with Rosemary and the bundled Ash model.
//!
//! Run with `[SMILES] [--model MODEL_DIRECTORY]`; the SMILES defaults to
//! ethanol, and `--model` loads an exported bundle instead of Ash.
use kekule::smiles;
use kekule_openff::{ForceField, NaglModel};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let usage = "usage: parameterize [SMILES] [--model MODEL_DIRECTORY]";
    let mut input = "CCO".to_owned();
    let mut directory = None;
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        if arg == "--model" {
            directory = Some(args.next().ok_or(usage)?);
        } else {
            input = arg;
        }
    }
    let mut molecules = smiles::to_molecules(&input)?;
    if molecules.len() != 1 {
        return Err("supply one connected molecule".into());
    }
    let mut molecule = molecules.remove(0);
    molecule.perceive()?;
    molecule.add_hydrogens()?;
    let atoms = molecule.atom_count();
    let model = match directory {
        Some(directory) => NaglModel::load(directory)?,
        None => NaglModel::ash()?,
    };
    let p = ForceField::rosemary()?.parameterize_molecule(molecule, &model)?;
    println!(
        "{atoms} atoms, {} bonds, {} angles, {} proper torsions, {} improper terms",
        p.bonds().len(),
        p.angles().len(),
        p.proper_torsions().len(),
        p.improper_torsions().len()
    );
    println!("Charge source: {:?}", p.charge_sources());
    println!("Charges (e): {:?}", p.charges().value());
    Ok(())
}
