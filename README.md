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

`kekule` is an experimental chemistry backend scoped for both small molecules and macromolecules. The project is intended to cover regular cheminformatics workflows, as well as modeling-oriented tasks. `kekule` is human-architected and AI-coded.

The architectural contract lives in [`ARCHITECTURE.md`](ARCHITECTURE.md). Optional external-reference comparisons are documented in [`benchmarks/README.md`](benchmarks/README.md).

> [!NOTE]
> `kekule` is in early development. Breaking API changes will happen without notice.


## Concept

The `Molecule` type is the canonical represented graph for one connected molecule. `Molecule` is wrapped as either `SmallMolecule` or `MacroMolecule`. `SmallMolecule` handles ordinary cheminformatics workflows, while `MacroMolecule` pairs one connected graph with a `SmcraHierarchy`.

```text
Molecule (one connected graph)
┣ SmallMolecule
┗ MacroMolecule (Molecule + SmcraHierarchy)
```

Higher, modeling-based objects are built around `Topology`: an immutable topological construct containing one or multiple molecule instances. Source-level identities such as an mmCIF entity or chain may map to more than one represented molecule instance when the observed structure contains a genuine unresolved gap.


```text
Topology = reusable definitions + explicit instances + dense ordering
┣ Model      = Topology + Positions + optional cell + AtomData + BondData
┣ Ensemble   = Topology + finite non-temporal members
┗ Trajectory = Topology + ordered frames / reusable streaming buffers
```

### Small-molecule pipeline

Kekule separates the basic chemistry pipeline into clear stages:

- **Parse** reads source syntax into a format-specific `*Document`.
- **Interpret** translates source-asserted chemistry and publishes a canonical, format-independent `Molecule`.
- **Perceive** derives model-dependent chemical interpretation such as valence, rings, and aromaticity.
- **CIP** optionally derives stereochemical descriptors such as `R`/`S` and `E`/`Z`.

Interpretation performs deterministic, meaning-preserving representation canonicalization before a molecule becomes observable. Chemistry-changing standardization such as tautomer or protonation-state selection is a separate future concern. High-level APIs may compose parsing and interpretation while perception remains explicit.

Expert workflows can keep every stage explicit:

```rust
use std::error::Error;

use kekule::{smiles, small::SmallMolecule};

fn load_explicitly(input: &str) -> Result<SmallMolecule, Box<dyn Error>> {
    let document = smiles::parse_str(input)?;
    let mut molecule = smiles::interpret(&document)?.to_molecule()?;
    molecule.perceive()?;
    Ok(molecule)
}
```


## Basic Usage

Parse and inspect a simple chiral molecule, assign its stereochemistry, detect rotatable bonds, and write it back to SMILES:

```rust
use std::error::Error;

use kekule::{
    rotatable_bonds::{self, RotatableBondOptions},
    small::SmallMolecule,
    stereo,
};

fn main() -> Result<(), Box<dyn Error>> {
    // Parse and canonically interpret a chiral amino acid, then perceive it.
    let mut molecule = SmallMolecule::from_smiles("C[C@@H](C(=O)O)N")?;
    molecule.perceive()?;

    // Assign absolute CIP descriptors to the perceived stereo elements.
    let stereochemistry = stereo::assign_cip_descriptors(molecule.as_molecule_mut())?;

    // Detect rotatable bonds using the strict small-molecule definition.
    let rotatable = rotatable_bonds::detect(
        molecule.as_molecule(),
        RotatableBondOptions::STRICT,
    );

    // Inspect basic graph properties and derived chemistry.
    println!("atoms: {}", molecule.atom_count());
    println!("bonds: {}", molecule.bond_count());
    println!("formal charge: {}", molecule.formal_charge());
    println!("strict rotatable bonds: {}", rotatable.len());
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
    // Parse and canonically interpret one SDF record without perceiving it.
    let input = fs::read_to_string("examples/ligand.sdf")?;
    let document = sdf::parse_str(&input, SdfParseOptions::default())?;
    let mut records = sdf::interpret(&document)?.to_records();
    assert_eq!(records.len(), 1, "expected one ligand record");

    // Preserve the record metadata while working on its molecule.
    let record = records.pop().expect("record count was checked");
    let title = record.title().to_owned();
    let data_fields = record.data_fields().to_vec();
    let mut ligand = record.to_molecule();
    ligand.perceive()?;

    // Inspect the canonical, perceived ligand before modeling it.
    println!("atoms: {}", ligand.atom_count());
    println!("bonds: {}", ligand.bond_count());
    println!("formal charge: {}", ligand.formal_charge());

    // Build a fixed-topology model from the ligand's first conformer.
    let conformer = ligand
        .as_molecule()
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
        .instance_to_conformer(instance, ligand.as_molecule_mut(), conformer)?;

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
