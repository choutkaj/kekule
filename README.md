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
    let mut molecule = read_smiles_str("c1ccccc1O", SmilesParseOptions)?;

    // Run chemistry perception, including valence, rings, and aromaticity.
    sanitize_small_molecule(&mut molecule, SanitizeOptions::default())?;

    // Inspect the graph.
    println!("atoms: {}", molecule.mol.atom_count());
    println!("bonds: {}", molecule.mol.bond_count());

    // Write deterministic non-isomeric canonical SMILES.
    let canonical = write_canonical_smiles(&molecule, CanonicalSmilesWriteOptions)?;
    println!("canonical SMILES: {canonical}");

    Ok(())
}
```

## License

`molecules` is available under the MIT license.