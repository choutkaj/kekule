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

The architectural contract lives in [`ARCHITECTURE.md`](ARCHITECTURE.md).

> [!NOTE]
> `kekule` is in early development and should be considered unstable. Breaking API changes will happen without notice.


## Concepts

`Molecule` is the single universal molecular type: one non-empty, connected,
geometry-independent entity. Its `Graph` owns authoritative chemistry, while
its `Perception` stores reconstructible derived chemistry. Coordinate-independent
residue and chain organization belongs to the system-level `Hierarchy` owned by
`Topology`; the hierarchy may be empty.

```text
Molecule = Graph + Perception

Topology = List of Molecules (stored as definitions and instances) + Hierarchy

Hierarchy = Chain -> Residue -> AtomSite

Model      = Topology + Positions
Ensemble   = Topology + one or more EnsembleMembers
Trajectory = Topology + one or more TrajectoryFrames
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

## Examples

### SMILES

Load and inspect a simple chiral molecule, then write canonical and isomeric
SMILES:

```rust
use std::error::Error;

use kekule::smiles;

fn main() -> Result<(), Box<dyn Error>> {
    // A dot-free SMILES produces one connected molecule.
    let mut molecules = smiles::to_molecules("C[C@@H](C(=O)O)N")?;
    let molecule = molecules.pop().expect("SMILES contains one molecule");

    println!("atoms: {}", molecule.atom_count());
    println!("bonds: {}", molecule.bond_count());
    println!("formal charge: {}", molecule.formal_charge());

    println!("canonical SMILES: {}", smiles::write_canonical(&molecule)?);
    println!("isomeric SMILES: {}", smiles::write_isomeric(&molecule)?);
    Ok(())
}
```

### SDF and mmCIF models

Load one small molecule from an SDF record and another from an mmCIF data
block, inspect and save both models, then combine them into a new model and
write it as mmCIF:

```rust
use std::{
    error::Error,
    fs::{self, File},
};

use kekule::{
    mmcif::{
        self, MmcifEntityClassifications, MmcifEntityKind, MmcifInterpretOptions, MmcifWriteOptions,
    },
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
    sdf::write_model_to(
        &mut File::create("ligand-copy.sdf")?,
        &sdf_model,
        SdfWriteOptions::default(),
    )?;

    // mmCIF is block-oriented. `interpret` requires exactly one atom-site block
    // and selects one coordinate model using the supplied options.
    let cif_document = mmcif::parse_str(&fs::read_to_string("cofactor.cif")?)?;
    let cif_interpretation = mmcif::interpret(&cif_document, MmcifInterpretOptions::default())?;
    let cif_model = cif_interpretation.model();
    assert_eq!(cif_model.topology().instance_count(), 1);
    print_model("mmCIF model", cif_model);
    mmcif::write_model_with_report_to(
        &mut File::create("cofactor-copy.cif")?,
        cif_model,
        cif_interpretation.report(),
        MmcifWriteOptions::default(),
    )?;

    // Coordinates are dense in molecule atom order for these one-molecule models.
    let mut builder = Model::builder();
    let sdf_id = builder.add_molecule(
        sdf_model.topology().molecules().next().unwrap().molecule(),
        sdf_model.positions(),
    )?;
    let cif_id = builder.add_molecule(
        cif_model.topology().molecules().next().unwrap().molecule(),
        cif_model.positions(),
    )?;
    let combined = builder.build()?;
    print_model("combined model", &combined);

    // A newly assembled Model has no format-specific entity roles. You have to supply them explicitly.
    let mut classifications = MmcifEntityClassifications::new();
    classifications.insert(sdf_id, MmcifEntityKind::NonPolymer)?;
    classifications.insert(cif_id, MmcifEntityKind::NonPolymer)?;
    mmcif::write_model_with_classifications_to(
        &mut File::create("combined.cif")?,
        &combined,
        &classifications,
        MmcifWriteOptions::default(),
    )?;

    Ok(())
}
```

The combined-model example intentionally requires one connected small molecule
per input. Multi-record SDF documents, multi-block or multi-model mmCIF
documents, and macromolecular systems require an explicit selection and entity
classification policy.

## Contributing

Currently not accepting contributions.

## License

`kekule` is available under either the [Apache License 2.0](LICENSE-APACHE) or the [MIT license](LICENSE-MIT), at your option.
