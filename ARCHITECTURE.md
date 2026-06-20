# Architecture

## Central decision

`molecules` uses one core molecular graph and domain-specific wrappers around it.

```rust
Molecule       // graph, conformers, properties, computed perception state
SmallMolecule  // small-molecule convenience API around Molecule
MacroMolecule  // Molecule plus BioHierarchy
```

The core graph must remain chemically general. Biomolecular hierarchy information belongs in `BioHierarchy`, not in `Atom`, unless the field is chemically general for both small molecules and macromolecules.

## Core graph

The core crate owns typed IDs, atoms, bonds, molecule topology, conformers, arbitrary properties, and perception state. Scoped atom and bond mutation invalidates computed chemistry state only when chemistry-relevant fields change; property and coordinate edits remain state-neutral. Cached perception objects are exposed only while their state is fresh. Parsing and perception are separate operations: a parser should not silently perform full sanitization unless the public API says it does.

The library source follows those ownership boundaries:

```text
crates/molecules/src/
  core/        typed IDs, elements, atom/bond payloads, conformers, graph state
  algorithms/  ring membership/sets, valence, and aromaticity
  chemistry/   transactional small-molecule sanitization
  io/          V2000/SDF, SMILES, and mmCIF parsing/writing
  bio/         macromolecule wrappers and BioHierarchy
  tests/       focused regression modules by domain
```

`lib.rs` is the stable public facade. It deliberately re-exports the same API and prelude while implementation modules remain private and cross-module internals use crate-private visibility.

## Small molecules

`SmallMolecule` is the RDKit-like layer. It should expose operations such as sanitization, valence perception, ring detection, aromaticity detection, stereochemistry assignment, small-molecule file I/O, and eventually substructure/search/canonicalization features.

Small-molecule algorithms should operate on the core `Molecule` where possible so they are reusable for ligands, residues, and extracted fragments from macromolecular structures.

## Macromolecules

`MacroMolecule` is the Biopython-like layer. It owns a `Molecule` plus a `BioHierarchy` that records models, chains, residues, atom-site labels, alternate locations, occupancy, B-factors, and mmCIF author/label identifiers.

The intended hierarchy is:

```text
MacroMolecule
  Molecule
  BioHierarchy
    Model
      Chain
        Residue
          AtomId
```

## File I/O

File readers should be explicit about whether they perform raw parsing only, parsing plus minimal normalization, or parsing plus sanitization/perception. Initial import/export priorities are Molfile/SDF for small molecules and PDB/mmCIF for macromolecules.

## Validation

Reference tools are used only in validation infrastructure: RDKit for small molecules and Biopython for macromolecular parsing/hierarchy behavior. Golden files should be normalized JSON and record the reference tool version used to generate them.

Repository tooling mirrors the validation workflow:

```text
crates/xtask/src/
  cli.rs
  corpus.rs
  features.rs
  dashboard.rs
  skills.rs
  validation/
    manifest.rs
    evidence.rs
    implementation.rs
    compare.rs
    status.rs
```

The command entry point only dispatches; corpus integrity, feature metadata, generated dashboard state, and validation evidence each have a focused owner.

## Feature-driven development

Every nontrivial capability gets a directory under `features/`. `feature.toml` is the machine-readable source of truth, and `feature.md` is the human-readable source of truth. The generated dashboard is derived from feature metadata and should not be hand-edited.
