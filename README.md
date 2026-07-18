<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/choutkaj/molecules/main/assets/molecules-logo-dark.svg">
    <img alt="MOLECULES - cheminformatics in Rust" src="https://raw.githubusercontent.com/choutkaj/molecules/main/assets/molecules-logo-light.svg" width="250">
  </picture>
</p>

<p align="center">
  <a href="https://github.com/choutkaj/molecules/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/choutkaj/molecules/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://github.com/choutkaj/molecules/blob/main/Cargo.toml"><img alt="MSRV 1.89" src="https://img.shields.io/badge/MSRV-1.89-blue.svg"></a>
  <a href="https://github.com/choutkaj/molecules/blob/main/LICENSE"><img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-blue.svg"></a>
</p>

`molecules` is an experimental pure-Rust chemistry backend scoped for both small molecules and macromolecules. The capabilities are bundled into features, which are parity-checked against established codebases - RDKit for small molecules and Biopython for macromolecules. This project is human-architected and AI-coded.

For feature overview and parity checks, see the [feature dashboard](https://choutkaj.github.io/molecules/).

> [!NOTE]
> `molecules` is in early development. Breaking API changes will happen without notice.


## Concept

The concept of this project is centered on molecules, which serve as a the most important unit of chemical information. The `Molecule` type is the foundational molecular graph of user-asserted molecular entity. It is usually a fully connected molecular graph of a single covalent structure, although disconnected graphs are allowed (for example for salts or complexes). `SmallMolecule` wraps one `Molecule` with ordinary cheminformatic workflows, while `MacroMolecule` pairs one `Molecule` with an `SmcraHierarchy` for biomolecular labels and structure work.

`Model` is an actual physical model of one or more instances of `SmallMolecule` and/or `MacroMolecule`. It holds an immutable topology and mutable atomic positions. The molecules in `Model` are not flattened into a disorganized bucket of atoms; instead, molecule instances are tracked and can be recognized and deciphered throughout any modeling work. It is the foundational type for molecular modeling. 


```text
Molecule (raw chemical graph)
 ├ SmallMolecule (graph + parameters -> small-molecule cheminformatics)
 ├ MacroMolecule (graph + SMCRA Hierarchy -> macromolecular cheminformatics)   
 │
 └─> Model (immutable topology, mutable positions)
     │
     └─> Optimization ──> Model (with optimized positions)
```

## Basic Usage

Parse and inspect a simple chiral molecule, assign its stereochemistry, and write it back to SMILES:

```rust
use std::error::Error;

use molecules::{perception::stereo, small::SmallMolecule};

fn main() -> Result<(), Box<dyn Error>> {
    // Parse a chiral amino acid and run the explicit sanitization workflow.
    let mut molecule = SmallMolecule::from_smiles_sanitized("C[C@@H](C(=O)O)N")?;

    // Assign absolute CIP descriptors to the perceived stereo elements.
    let stereochemistry = stereo::assign_cip_descriptors(molecule.graph_mut());

    // Inspect basic graph properties and the asserted molecular charge.
    println!("atoms: {}", molecule.atom_count());
    println!("bonds: {}", molecule.bond_count());
    println!("formal charge: {}", molecule.graph().formal_charge());
    for assignment in &stereochemistry.assigned {
        println!("stereo {:?}: {:?}", assignment.element, assignment.descriptor);
    }

    // Write canonical connectivity and a stereo-preserving SMILES form.
    println!("canonical SMILES: {}", molecule.to_canonical_smiles()?);
    println!("isomeric SMILES: {}", molecule.to_isomeric_smiles()?);
    Ok(())
}
```

## Modeling

Load a ligand from SDF, minimize its coordinates with the DREIDING force field, and write the optimized structure back to SDF. Run modeling examples with `--release` for optimized numerical kernels.

```rust
use std::{error::Error, fs};

use molecules::{
    modeling::{minimize, MinimizeOptions, Model},
    sdf::{self, SdfParseOptions, SdfRecord},
    units::MODEL_GRADIENT_UNIT,
};
use molecules_dreiding::DreidingPotential;

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

    // Prepare DREIDING explicitly, then minimize a clone of the model.
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
```

## License

`molecules` is available under the MIT license.
