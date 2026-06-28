<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="./assets/molecules-logo-dark.svg">
    <img alt="MOLECULES - cheminformatics in Rust" src="./assets/molecules-logo-light.svg" width="250">
  </picture>
</p>

`molecules` is an experimental pure-Rust backend for small-molecule and
macromolecular structure work, written in pure Rust. This project is human-architected and AI-coded.
The cheminformatic capabilities are bundled into features, which are validated by comparison with existing codes - DRkit for small molecules, and Biopython for macromolecules.

It has no RDKit or Biopython runtime dependency;
those projects are used only to generate external validation references.

The crate is version `0.0.0`, is not published, and makes no API stability or
chemistry-completeness promise. Successful validation covers only the named
features, fixtures, fields, and reference versions recorded in the generated
[feature dashboard](features/DASHBOARD.html).

> [!NOTE]
> It is early. API may break without warning.



## Validation

`implemented` and `validated` are different states. The dashboard records
successful comparisons, while `cargo xtask validate` recomputes content
addresses and is the authority for whether evidence is current in a checkout.
As of June 21, 2026, tiny, PubChem 100/1000, and PDB 10/100 comparisons are
built and passing locally. PL-REX and Enamine are declared but are not built
validation corpora.

```bash
cargo fmt --all -- --check
cargo check --workspace --all-targets
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
cargo xtask skills --check
cargo xtask corpus check --corpus tiny --require-data
cargo xtask validate --feature all --corpus tiny
```

Large corpus data is intentionally not committed. See
[RELEASE_READINESS.md](RELEASE_READINESS.md) for the clean-runner and release
blockers, and [FUZZING.md](FUZZING.md) for parser robustness campaigns.


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
