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

The public API is organized around a small prelude plus focused format and workflow modules. Parsing, sanitization, validation, and writing are separate steps so callers can choose when interpretation happens.

### Small molecules

Use `SmallMolecule` for the common small-molecule path. Parsing a SMILES string creates the graph without running perception; `sanitize` adds chemistry-derived state such as valence, rings, and aromaticity.

```rust
use molecules::prelude::*;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut mol = SmallMolecule::from_smiles("c1ccccc1O")?;
    mol.sanitize()?;

    assert_eq!(mol.atom_count(), 7);
    assert_eq!(mol.bond_count(), 7);

    let canonical = mol.to_canonical_smiles()?;
    println!("canonical SMILES: {canonical}");

    Ok(())
}
```

### Macromolecules

Use the `mmcif` facade for macromolecular file I/O. Reading mmCIF parses atom-site data into a `MacroMolecule`; validation remains an explicit follow-up step.

```rust
use molecules::mmcif::{self, MmcifParseOptions};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = r#"
data_demo
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_seq_id
C C1 GLY A 1
"#;

    let macro_mol = mmcif::read_str(input, MmcifParseOptions::default())?;
    macro_mol.validate()?;

    println!("atoms: {}", macro_mol.graph().atom_count());

    Ok(())
}
```

## License

`molecules` is available under the MIT license.
