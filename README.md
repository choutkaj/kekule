<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/choutkaj/kekule/main/assets/kekule-logo-dark.svg">
    <img alt="KEKULE - cheminformatics in Rust" src="https://raw.githubusercontent.com/choutkaj/kekule/main/assets/kekule-logo-light.svg" width="250">
  </picture>
</p>

<p align="center">
  <a href="https://github.com/choutkaj/kekule/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/choutkaj/kekule/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://github.com/choutkaj/kekule/blob/main/Cargo.toml"><img alt="MSRV 1.89" src="https://img.shields.io/badge/MSRV-1.89-blue.svg"></a>
  <a href="#license"><img alt="License: MIT OR Apache-2.0" src="https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg"></a>
</p>

`kekule` is an experimental chemistry backend scoped for both small molecules and macromolecules. The project is intended to cover regular cheminformatic workflows, as well as modeling-oriented tasks. `kekule` is human-architected and AI-coded.

For feature overview and parity checks, see the [feature dashboard](https://choutkaj.github.io/kekule/).

> [!NOTE]
> `kekule` is in early development. Breaking API changes will happen without notice.


## Concept

The `Molecule` type is the raw graph kernel for one user-asserted molecular entity. It is usually a connected molecular graph, although disconnected graphs are allowed (to handle salts, complexes, and such). The `Molecule` is wrapped as either `SmallMolecule` or `MacroMolecule`. `SmallMolecule` handles ordinary cheminformatics workflows; `MacroMolecule` pairs a graph with a `SmcraHierarchy`.

```text
Molecule
┣ SmallMolecule
┗ MacroMolecule (Molecule + SmcraHierarchy)
```

Higher, modeling-based objects are built around `Topology`: an immutable topological construct containing one or multiple molecule instances. 


```text
Topology = reusable definitions + explicit instances + dense ordering
┣ Model      = Topology + one Configuration
┣ Ensemble   = Topology + finite non-temporal members
┗ Trajectory = Topology + ordered frames / reusable streaming buffers
```


## Basic Usage

Parse and inspect a simple chiral molecule, assign its stereochemistry, and write it back to SMILES:

```rust
use std::error::Error;

use kekule::{perception::stereo, small::SmallMolecule};

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

Load a ligand from SDF, minimize its coordinates with the DREIDING force field, and write the optimized structure back to SDF: 

```rust
use std::{error::Error, fs};

use kekule::{
    modeling::{minimize, MinimizeOptions},
    sdf::{self, SdfParseOptions, SdfRecord},
    structure::Model,
    units::MODEL_GRADIENT_UNIT,
};
use kekule_potentials::dreiding::{DreidingPotential, DreidingPrepareOptions};

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
```

## Contributing

Currently not accepting contributions.

## License

`kekule` is available under either the [Apache License 2.0](LICENSE-APACHE) or the [MIT license](LICENSE-MIT), at your option.
