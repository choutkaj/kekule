<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="./assets/molecules-logo-dark.svg">
    <img alt="MOLECULES - cheminformatics in Rust" src="./assets/molecules-logo-light.svg" width="250">
  </picture>
</p>

`molecules` is an experimental pure-Rust backend for both small-molecule and macromolecular structure work. This project is human-architected and AI-coded. The cheminformatic capabilities are bundled into features, which are validated by comparison with established codebases - RDkit for small molecules and Biopython for macromolecules.

For already-implemented features, see the [feature dashboard](features/DASHBOARD.html).

> [!NOTE]
> It is early. API may break without warning.


## Basic Usage

### Small molecules

`molecules` uses a small-molecule workflow where parsing, sanitization/perception, and writing are explicit steps. Parse input into a molecular graph, sanitize it when you want chemistry-derived state such as valence or aromaticity, then write or inspect the result. This keeps raw file I/O separate from chemical perception.

```rust
use molecules::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Parse SMILES into a small-molecule graph. Parsing does not sanitize.
    let mut molecule = SmallMolecule::from_smiles("c1ccccc1O")?;

    // Run chemistry perception, including valence, rings, and aromaticity.
    molecule.sanitize()?;

    // Inspect the graph.
    println!("atoms: {}", molecule.atom_count());
    println!("bonds: {}", molecule.bond_count());

    // Write deterministic non-isomeric canonical SMILES.
    let canonical = molecule.to_canonical_smiles()?;
    println!("canonical SMILES: {canonical}");

    Ok(())
}
```

## License

`molecules` is available under the MIT license.
