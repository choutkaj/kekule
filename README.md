<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/choutkaj/kekule/main/assets/kekule-logo-dark.svg">
    <img alt="KEKULE - cheminformatics in Rust" src="https://raw.githubusercontent.com/choutkaj/kekule/main/assets/kekule-logo-light.svg" width="250">
  </picture>
</p>

<p align="center">
  <a href="https://github.com/choutkaj/kekule/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/choutkaj/kekule/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://crates.io/crates/kekule"><img alt="crates.io version" src="https://img.shields.io/crates/v/kekule.svg"></a>
  <a href="https://github.com/choutkaj/kekule/blob/main/Cargo.toml"><img alt="MSRV 1.89" src="https://img.shields.io/badge/MSRV-1.89-blue.svg"></a>
  <a href="#license"><img alt="License: MIT OR Apache-2.0" src="https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg"></a>
</p>

`kekule` is an experimental chemistry backend scoped for both small molecules and macromolecules. The project is intended to cover regular cheminformatics workflows, as well as modeling-oriented tasks. `kekule` is human-architected and AI-coded.

The architectural contract lives in [`ARCHITECTURE.md`](ARCHITECTURE.md).

> [!NOTE]
> `kekule` is in early development and should be considered unstable. Breaking API changes will happen without notice.

## Installation

```sh
cargo add kekule
```

Trajectory I/O and selected potentials are available separately through the sister crates
`kekule-traj` and `kekule-potentials`.


## Concepts

`Molecule` is the foundational type storing one molecule without its geometry. `Molecule` must be a connected graph. Its `Graph` owns authoritative chemistry, while its `Perception` stores derived chemical perception such as rings and aromaticity. Topology is a collection of one or more `Molecule`s together with `Hierarchy` (`Chain`, `Residue`, `AtomSite`). Molecules in `Topology` are not stored naively, but as `Definition`s and `Instance`s. For example, a hundred water molecules will be stored as one `Definition` and a hundred `Instances`. Coordinates are detached from `Molecule`: one realization of a system's geometry is a `Conformation`, holding `Positions`, an optional periodic cell, occupancies, B-factors, and per-atom annotations.
```text
Molecule = Graph + Perception + Properties
Topology = collection of Molecules (stored as definitions and instances) + Hierarchy
Hierarchy = Chain -> Residue -> AtomSite
```
Higher, modeling-oriented objects are built around `Topology` and contain actual instances of molecules including their coordinates. `Model` is literally a model of one or more molecules. `Ensemble` and `Trajectory` are distinct collections of realizations of one system: an ensemble is a weighted, unordered sample (conformers, Monte Carlo samples, NMR models), and every member carries a statistical weight; a trajectory is a time-ordered sequence of frames produced by dynamics.
```text
Model           = Topology + Conformation
Ensemble        = Topology + weighted EnsembleMembers (Conformation + weight)
Trajectory      = Topology + time-ordered TrajectoryFrames (Conformation + time, step, velocities, forces)
```

## Examples

### SMILES

Load and inspect a simple chiral molecule, assign its stereochemistry, detect
its rotatable bonds, then write canonical and isomeric SMILES:

```rust
use std::error::Error;

use kekule::{
    rotatable_bonds::{self, RotatableBondOptions},
    smiles::{self, SmilesWriteOptions},
    stereo,
};

fn main() -> Result<(), Box<dyn Error>> {
    // A dot-free SMILES produces one connected molecule.
    let mut molecules = smiles::to_molecules("C[C@@H](C(=O)O)N")?;
    let mut molecule = molecules.pop().expect("SMILES contains one molecule");
    molecule.perceive()?;

    println!("atoms: {}", molecule.atom_count());
    println!("bonds: {}", molecule.bond_count());
    println!("formal charge: {}", molecule.formal_charge());
    
    // Print canonical and isomeric SMILES
    let canonical = smiles::write(&molecule, SmilesWriteOptions::canonical())?;
    let isomeric = smiles::write(&molecule, SmilesWriteOptions::isomeric())?;
    println!("canonical SMILES: {canonical}");
    println!("isomeric SMILES: {isomeric}");
    
    // Assign and print stereochemistry
    let stereochemistry = stereo::assign_cip_descriptors(&mut molecule)?;
    for assignment in &stereochemistry.assigned {
        println!("stereo {}: {:?}", assignment.element, assignment.descriptor);
    }
    
    // Assign and print rotatable bonds
    let rotatable = rotatable_bonds::detect(&molecule, RotatableBondOptions::STRICT)?;
    for &bond_id in rotatable.bond_ids() {
        let bond = molecule.bond(bond_id)?;
        println!("rotatable bond {bond_id}: {}-{}", bond.a(), bond.b());
    }
    Ok(())
}
```

### SDF and mmCIF models

Load one small molecule from an SDF record and another from an mmCIF data block, inspect and save both models, then combine them into a new model and write it as mmCIF:

```rust
use std::{
    error::Error,
    fs::{self, File},
};

use kekule::{
    mmcif::{self, MmcifBlockSource, MmcifInterpretOptions, MmcifWriteOptions},
    sdf::{self, SdfWriteOptions},
    structure::Model,
};

fn print_model(label: &str, model: &Model) {
    println!(
        "{label}: {} molecules, {} atoms, {} bonds",
        model.topology().instance_count(),
        model.atom_count(),
        model.topology().bond_count(),
    );
}

fn main() -> Result<(), Box<dyn Error>> {
    // SDF is record-oriented. This example deliberately loads one record.
    let sdf_document = sdf::parse_str(&fs::read_to_string("ligand.sdf")?)?;
    let sdf_record = sdf_document.records().first().expect("SDF has a record");
    let sdf_model = sdf_record.to_model()?;
    assert_eq!(sdf_model.topology().instance_count(), 1);
    print_model("SDF model", &sdf_model);
    sdf::write_to(
        &mut File::create("ligand-copy.sdf")?,
        [&sdf_model],
        SdfWriteOptions::default(),
    )?;

    // mmCIF is block-oriented. `interpret` requires exactly one atom-site block
    // and selects one coordinate model using the supplied options.
    let cif_document = mmcif::parse_str(&fs::read_to_string("cofactor.cif")?)?;
    let cif_interpretation = mmcif::interpret(&cif_document, MmcifInterpretOptions::default())?;
    let cif_model = cif_interpretation.model();
    assert_eq!(cif_model.topology().instance_count(), 1);
    print_model("mmCIF model", cif_model);
    // Writing with the interpretation report keeps the source entity semantics.
    let block = MmcifBlockSource::model(cif_model)
        .with_reports(std::slice::from_ref(cif_interpretation.report()));
    mmcif::write_to(
        &mut File::create("cofactor-copy.cif")?,
        [block],
        MmcifWriteOptions::default(),
    )?;

    // Initiate Model builder and combine both Models
    let mut builder = Model::builder();
    builder.add_molecule(
        sdf_model.topology().molecules().next().unwrap().molecule().clone(),
        sdf_model.positions(),
    )?;
    builder.add_molecule(
        cif_model.topology().molecules().next().unwrap().molecule().clone(),
        cif_model.positions(),
    )?;
    let combined = builder.build()?;
    print_model("combined model", &combined);

    // Write the combined Model into mmCIF file
    mmcif::write_to(
        &mut File::create("combined.cif")?,
        [&combined],
        MmcifWriteOptions::default(),
    )?;

    Ok(())
}
```
The combined-model example intentionally requires one connected small molecule per input. Multi-record SDF documents, multi-block or multi-model mmCIF documents, and macromolecular systems require explicit source selection.

## Contributing

Currently not accepting contributions.

## License

`kekule` is available under either the [Apache License 2.0](LICENSE-APACHE) or the [MIT license](LICENSE-MIT), at your option.
