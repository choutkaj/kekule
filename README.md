<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="./assets/molecules-logo-dark.svg">
    <img alt="MOLECULES - cheminformatics in Rust" src="./assets/molecules-logo-light.svg" width="250">
  </picture>
</p>

`molecules` is an experimental pure-Rust backend for both small-molecule and macromolecular structure work. This project is human-architected and AI-coded. The cheminformatic capabilities are bundled into features, which are validated by comparison with established codebases - RDkit for small molecules and Biopython for macromolecules.

For already-implemented features, see the [rendered feature dashboard](https://choutkaj.github.io/molecules/) or inspect the generated [dashboard source](features/DASHBOARD.html).

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

### Molecular modelling

The initial modelling layer creates a fixed-topology `MolecularModel` from one
or more selected small-molecule conformers. Potentials and minimization remain
explicit namespaced operations; modelling types are not part of the prelude.

```rust
use molecules::core::{Atom, BondOrder, Conformer, Element, Molecule, Point3};
use molecules::modeling::potential::{HarmonicBondParameter, HarmonicBondPotential};
use molecules::modeling::{minimize, MinimizeOptions, MolecularModel};
use molecules::small::SmallMolecule;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut graph = Molecule::new();
    let carbon = graph.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
    let oxygen = graph.add_atom(Atom::new(Element::from_symbol("O").unwrap()));
    let source_bond = graph.add_bond(carbon, oxygen, BondOrder::Single)?;

    let mut conformer = Conformer::new();
    conformer.set_position(carbon, Point3::new(0.0, 0.0, 0.0));
    conformer.set_position(oxygen, Point3::new(2.0, 0.0, 0.0));
    let conformer = graph.add_conformer(conformer);
    let molecule = SmallMolecule::from_graph(graph);

    let mut builder = MolecularModel::builder();
    let mapping = builder.add_component(&molecule, conformer)?;
    let model = builder.build()?;
    let bond = mapping.bond(source_bond).unwrap();
    let mut potential = HarmonicBondPotential::new(
        &model,
        [HarmonicBondParameter::new(bond, 1.2, 100.0)],
    )?;
    let minimized = minimize(&model, &mut potential, MinimizeOptions::default())?;

    println!("final energy: {} kJ/mol", minimized.final_energy);
    Ok(())
}
```

### Macromolecules

Use the `mmcif` facade for multi-entity structure input. Parsing preserves the
mmCIF document; a separate interpretation step produces clean macromolecules,
small molecules, ions, and individual solvent molecules.

```rust
use molecules::mmcif::{self, MmcifInterpretOptions, MmcifParseOptions};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let input = r#"
data_demo
loop_
_entity.id
_entity.type
1 polymer
loop_
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.auth_seq_id
C C1 GLY A 1 1
"#;

    let document = mmcif::parse_str(input, MmcifParseOptions::default())?;
    let interpreted = mmcif::interpret(&document, MmcifInterpretOptions::default())?;
    let macro_mol = interpreted
        .contents()
        .macromolecules()
        .next()
        .expect("one polymer molecule");
    macro_mol.validate()?;

    println!("atoms: {}", macro_mol.graph().atom_count());

    Ok(())
}
```

## License

`molecules` is available under the MIT license.
