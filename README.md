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

`Molecule` is the single universal molecular type: one non-empty, connected,
geometry-independent entity. Its `Graph` owns authoritative chemistry, while
its `Perception` stores reconstructible derived chemistry. Coordinate-independent
residue and chain organization belongs to the system-level `Hierarchy` owned by
`Topology`; the hierarchy may be empty.

```text
Molecule = Graph + Perception

Topology = molecule definitions + instances + dense ordering + Hierarchy

Hierarchy = Chain -> Residue -> AtomSite -> InstanceAtomId

Model      = Topology + Positions + optional cell + AtomData + BondData
Ensemble   = Topology + finite non-temporal realizations
Trajectory = Topology + temporally ordered frames
```

Coordinates are detached from `Molecule`. Dense `Positions` enter at the
`Model`, `Ensemble`, or trajectory layer. Structural construction and editing
use `MoleculeEditor`; only `finish()` can publish a molecule, after validating
that it is non-empty and connected.

Higher, modeling-based objects are built around `Topology`: an immutable
topological system containing one or more molecule instances. Source-level
identities such as an mmCIF entity or chain may map to more than one represented
molecule instance when the observed structure contains a genuine unresolved
gap.

### Molecular pipeline

Kekule separates the basic chemistry pipeline into clear stages:

- **Parse** reads source syntax into a format-specific `*Document`.
- **Interpret** translates source-asserted chemistry and publishes source-ordered, connected `Molecule` components.
- **Perceive** derives model-dependent chemical interpretation such as valence, rings, and aromaticity.
- **CIP** optionally derives stereochemical descriptors such as `R`/`S` and `E`/`Z`.

Interpretation performs deterministic, meaning-preserving representation canonicalization before a molecule becomes observable. Chemistry-changing standardization such as tautomer or protonation-state selection is a separate future concern. High-level APIs may compose parsing and interpretation while perception remains explicit.

Expert workflows can keep every stage explicit:

```rust
use std::error::Error;

use kekule::{core::Molecule, smiles};

fn load_explicitly(input: &str) -> Result<Vec<Molecule>, Box<dyn Error>> {
    let document = smiles::parse_str(input)?;
    let mut molecules = smiles::interpret(&document)?.to_molecules();
    for molecule in &mut molecules {
        molecule.perceive()?;
    }
    Ok(molecules)
}
```


## Basic Usage

Parse and inspect a simple chiral molecule, assign its stereochemistry, detect rotatable bonds, and write it back to SMILES:

```rust
use std::error::Error;

use kekule::{
    core::Molecule,
    rotatable_bonds::{self, RotatableBondOptions},
    stereo,
};

fn main() -> Result<(), Box<dyn Error>> {
    // Parse and canonically interpret a chiral amino acid, then perceive it.
    let mut molecules = Molecule::from_smiles("C[C@@H](C(=O)O)N")?;
    let mut molecule = molecules.pop().expect("SMILES contains one component");
    molecule.perceive()?;

    // Assign absolute CIP descriptors to the perceived stereo elements.
    let stereochemistry = stereo::assign_cip_descriptors(&mut molecule)?;

    // Detect rotatable bonds using the strict small-molecule definition.
    let rotatable = rotatable_bonds::detect(&molecule, RotatableBondOptions::STRICT);

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

Combine a connected molecule with detached coordinates, then minimize the
resulting model with the DREIDING force field:

```rust
use std::error::Error;

use kekule::{
    core::Molecule,
    geometry::Point3,
    modeling::{minimize, MinimizeOptions},
    structure::{Model, Positions},
    units::{Quantity, ANGSTROM, KILOJOULE_PER_MOLE_PER_ANGSTROM},
};
use kekule_potentials::dreiding::{DreidingPotential, DreidingPrepareOptions};

fn main() -> Result<(), Box<dyn Error>> {
    let mut molecules = Molecule::from_smiles("CCO")?;
    let mut ligand = molecules.pop().expect("SMILES contains one component");
    ligand.perceive()?;

    // Molecule carries no geometry. Supply dense coordinates in atom order.
    let positions = Positions::new(Quantity::new(
        ligand
            .atom_ids()
            .enumerate()
            .map(|(index, _)| Point3::new(index as f64, 0.0, 0.0))
            .collect::<Vec<_>>(),
        ANGSTROM,
    ))?;
    let mut builder = Model::builder();
    builder.add_molecule(&ligand, &positions)?;
    let model = builder.build()?;

    // Prepare DREIDING explicitly, then minimize the model.
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

    Ok(())
}
```

## Contributing

Currently not accepting contributions.

## License

`kekule` is available under either the [Apache License 2.0](LICENSE-APACHE) or the [MIT license](LICENSE-MIT), at your option.
