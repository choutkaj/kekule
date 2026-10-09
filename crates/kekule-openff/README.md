# kekule-openff

SMIRNOFF force-field parameterization with native NAGL partial charges for
[`kekule`](https://crates.io/crates/kekule). It bundles the OpenFF Rosemary force
field and the Ash charge model, so a complete parameterization needs no Python,
external files, or network access.

```rust
use kekule::smiles;
use kekule_openff::{ForceField, NaglModel};

let mut molecule = smiles::to_molecules("CCO")?.remove(0);
molecule.perceive()?;
molecule.add_hydrogens()?;
let model = NaglModel::ash()?;
let parameters = ForceField::rosemary()?.parameterize_molecule(molecule, &model)?;
println!("{} bonds, charges {:?}", parameters.bonds().len(), parameters.charges().value());
# Ok::<(), Box<dyn std::error::Error>>(())
```

Systems are parameterized from a shared topology with
`ForceField::parameterize`, assigning each reusable molecule definition once.
There is no molecule size limit. Custom OFFXML within the supported SMIRNOFF
subset and other exported NAGL bundles are supported; see
[CONTRACT.md](CONTRACT.md) for the supported subset, assignment semantics, error
kinds, and validation record. Energies and gradients live in `kekule-potentials`.

The Rosemary force field and the Ash model are licensed under CC BY 4.0; see
[THIRD_PARTY.md](THIRD_PARTY.md). Building requires a C toolchain for the
official InChI library. The Rust code is MIT OR Apache-2.0.
