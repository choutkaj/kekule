# Architecture

## Design goal

`molecules` is a pure-Rust backend for small-molecule cheminformatics and macromolecular structure work. The API should feel small, predictable, and Rust-native at the surface while preserving lower-level access for algorithm development.

The most important architectural rule is that the public API is not allowed to become a flat dump of every internal module. The crate should expose a stable, ergonomic facade first, and keep implementation details behind focused namespaces.

## Central API decision

`molecules` uses one raw molecular graph kernel and domain-specific wrappers around it.

```rust
Molecule        // raw graph kernel: topology, payloads, conformers, properties, caches
SmallMolecule   // chemistry-aware small-molecule wrapper around Molecule
MacroMolecule   // structure-biology wrapper around Molecule plus BioHierarchy
```

`Molecule` is not an RDKit-style fully interpreted chemical object. It is a graph and payload container. It may be incomplete, unsanitized, chemically invalid, or freshly parsed from a lossy file format. Chemical interpretation belongs to explicit perception/sanitization steps and to domain wrappers.

The domain wrappers provide the user-facing meaning:

- `SmallMolecule` owns small-molecule convenience workflows: SMILES, SDF/Molfile, sanitization, valence, ring perception, aromaticity, stereochemistry, canonicalization, descriptors, fingerprints, and eventually substructure search.
- `MacroMolecule` owns macromolecular structure workflows: models, chains, residues, atom-site metadata, alternate locations, occupancy, B-factors, PDB/mmCIF identifiers, and coordinate-heavy structure operations.

Biomolecular hierarchy information must stay in `BioHierarchy`, not in core `Atom`, unless the field is chemically general for both small molecules and macromolecules.

## Public crate shape

The desired public namespace is:

```text
molecules
  core         raw graph kernel: Molecule, Atom, Bond, IDs, Element, Conformer
  small        SmallMolecule and small-molecule convenience API
  bio          MacroMolecule, BioHierarchy, Model, Chain, Residue, AtomSite
  smiles       SMILES parsing/writing facade
  molfile      V2000/V3000 Molfile parsing/writing facade
  sdf          SDF record parsing/writing facade
  perception   sanitization, valence, rings, aromaticity, stereo models
  canon        canonical ranking, identity, canonical graph utilities
  prelude      small set of common user-facing imports
```

Implementation modules may have different internal names while the crate is evolving, but the public surface should move toward these namespaces. In particular, `lib.rs` should not use blanket root re-exports such as `pub use algorithms::*`, `pub use io::*`, or `pub use core::*` as the long-term public API.

Preferred public usage:

```rust
use molecules::prelude::*;

let mut mol = SmallMolecule::from_smiles("c1ccccc1O")?;
mol.sanitize()?;
let canonical = mol.to_canonical_smiles()?;
```

Equivalent explicit namespace usage:

```rust
let mut mol = molecules::smiles::read_str("c1ccccc1O")?;
molecules::perception::sanitize(&mut mol, SanitizeOptions::default())?;
let canonical = molecules::smiles::write_canonical(&mol)?;
```

Low-level graph usage should remain available without pretending the graph is already chemically valid:

```rust
use molecules::core::*;

let mut graph = Molecule::new();
let c = graph.add_atom(Atom::new(Element::from_symbol("C").unwrap()));
let o = graph.add_atom(Atom::new(Element::from_symbol("O").unwrap()));
graph.add_bond(c, o, BondOrder::Double)?;
```

## API tiers

The crate should have three explicit API tiers.

### Tier 1: happy path

This is what most users should discover first.

```rust
let mut mol = SmallMolecule::from_smiles("CC(=O)O")?;
mol.sanitize()?;
let smiles = mol.to_canonical_smiles()?;
```

The happy path should prioritize clarity and common cheminformatics workflows over exposing every option.

### Tier 2: namespaced functional API

This is for users who want explicit staged control.

```rust
let mut mol = molecules::smiles::read_str("CC(=O)O")?;
molecules::perception::sanitize_with_options(&mut mol, options)?;
let smiles = molecules::smiles::write_canonical_with_options(&mol, write_options)?;
```

This tier should be stable enough for normal users.

### Tier 3: expert/algorithm API

This is for low-level algorithm work, validation, and future research.

Examples include custom ring-perception options, aromaticity model choices, canonical atom ranking, resource-limit diagnostics, and internal work counters.

These APIs should live under focused modules such as `perception::rings`, `perception::aromaticity`, and `canon`. They should not all be imported by the default prelude.

## Prelude policy

The prelude should be intentionally small. It should contain the common types and functions needed for typical user code, not every public item in the crate.

Recommended prelude contents:

```rust
pub mod prelude {
    pub use crate::core::{Atom, AtomId, Bond, BondId, BondOrder, Element, Molecule};
    pub use crate::small::{SmallMolecule, SanitizeOptions, SanitizeReport};
    pub use crate::bio::{MacroMolecule, BioHierarchy};
    pub use crate::smiles::{SmilesParseOptions, SmilesWriteOptions, CanonicalSmilesWriteOptions};
}
```

Do not put all parser errors, all algorithm reports, all ring-work diagnostics, all low-level perception functions, or every biomolecular metadata type into the prelude. Users can import those from their modules when needed.

## Core graph

The core graph owns:

- typed IDs: `AtomId`, `BondId`, `ConformerId`
- element and atom payloads
- bond payloads
- graph topology and adjacency
- optional conformers
- arbitrary properties
- cached perception results
- perception freshness state

Core invariants:

- `AtomId`, `BondId`, and `ConformerId` are typed IDs, not plain indices.
- IDs are stable for the lifetime of a molecule object.
- Deleted atoms/bonds leave holes; IDs are not reused.
- `Molecule` may represent chemically invalid or unsanitized input.
- File parsers may create raw `Molecule` values that require later sanitization.
- Cached perception data is only valid when the corresponding perception state is fresh.
- Topology-changing edits invalidate all dependent perception state.
- Chemistry-sensitive payload edits invalidate dependent perception state.
- Property-only and coordinate-only edits should not invalidate topology-derived state.

The perception cache and freshness flags are implementation details. Users should query chemistry through methods such as `rings()`, `ring_count()`, `atom_in_ring()`, `is_aromatic()`, or explicit perception reports, rather than manually editing cache state.

## SmallMolecule

`SmallMolecule` should become a real public abstraction, not just a public field wrapper around `Molecule`.

Preferred shape:

```rust
pub struct SmallMolecule {
    graph: Molecule,
    data: SmallMoleculeData,
}
```

Initially `SmallMoleculeData` may be empty or crate-private, but keeping a slot for it prevents `SmallMolecule` from becoming a vestigial newtype.

Recommended methods:

```rust
impl SmallMolecule {
    pub fn new() -> Self;
    pub fn from_graph(graph: Molecule) -> Self;
    pub fn into_graph(self) -> Molecule;

    pub fn graph(&self) -> &Molecule;
    pub fn graph_mut(&mut self) -> &mut Molecule; // should conservatively invalidate derived chemistry

    pub fn from_smiles(input: &str) -> Result<Self, SmilesParseError>;
    pub fn from_smiles_sanitized(input: &str) -> Result<Self, SmallMoleculeReadError>;

    pub fn sanitize(&mut self) -> Result<SanitizeReport, SanitizeError>;
    pub fn sanitize_with_options(&mut self, options: SanitizeOptions) -> Result<SanitizeReport, SanitizeError>;

    pub fn atom_count(&self) -> usize;
    pub fn bond_count(&self) -> usize;
    pub fn atoms(&self) -> impl Iterator<Item = (AtomId, &Atom)>;
    pub fn bonds(&self) -> impl Iterator<Item = (BondId, &Bond)>;

    pub fn to_smiles(&self) -> Result<String, MolWriteError>;
    pub fn to_canonical_smiles(&self) -> Result<String, MolWriteError>;
}
```

`SmallMolecule::from_smiles` should parse only unless its documentation explicitly says it sanitizes. Sanitizing convenience constructors are allowed, but their names must say so, for example `from_smiles_sanitized` or `smiles::read_sanitized_str`.

The public API should avoid requiring users to write `molecule.mol.atom_count()`. Prefer `molecule.atom_count()` and `molecule.graph().atom_count()`.

Long term, the raw graph field should not be public. Use `graph()`, `graph_mut()`, and `into_graph()` instead.

## MacroMolecule and BioHierarchy

`MacroMolecule` is the macromolecular wrapper:

```rust
pub struct MacroMolecule {
    graph: Molecule,
    hierarchy: BioHierarchy,
}
```

The intended hierarchy is:

```text
MacroMolecule
  Molecule
  BioHierarchy
    Model
      Chain
        Residue
          AtomSite -> AtomId
```

`BioHierarchy` owns structure-biology labels and metadata. It should support mmCIF/PDB author and label identifiers without contaminating the core atom payload with format-specific fields.

`MacroMolecule` should expose:

- `graph()` / `graph_mut()`
- `hierarchy()` / `hierarchy_mut()`
- model, chain, residue, and atom-site iterators
- lookup from `AtomId` to `AtomSite`
- coordinate/conformer access through the graph

Small-molecule algorithms should operate on `Molecule` where reasonable so extracted ligands, residues, and fragments from `MacroMolecule` can reuse them.

## File I/O

Parsing, sanitization, and writing are separate operations.

Default readers should parse raw input and preserve what the format says. They should not silently run full sanitization unless the function name and documentation say so.

Preferred naming pattern:

```rust
molecules::smiles::read_str(input)
molecules::smiles::read_sanitized_str(input)
molecules::smiles::write(&mol)
molecules::smiles::write_canonical(&mol)

molecules::molfile::read_v2000_str(input)
molecules::molfile::write_v2000(&mol)
molecules::molfile::read_v3000_str(input)
molecules::molfile::write_v3000(&mol)

molecules::sdf::read_v2000_str(input)
molecules::sdf::read_v2000_records(input)
molecules::sdf::write_v2000(&records)
```

Compatibility aliases may exist temporarily, but the namespaced API should be the documented target.

Options structs should be future-proof. Avoid public zero-sized options types becoming permanent dead ends. If an options struct is public, prefer a real `Default` implementation and consider `#[non_exhaustive]` before the crate is published.

## Perception and sanitization

Chemical perception is explicit and staged.

The sanitizer is the high-level small-molecule perception pipeline. It should remain transactional: stage changes on a clone or temporary graph, return a useful error if perception fails, and only replace the original molecule after success.

Sanitization stages include, at minimum:

```text
normalize import quirks
perceive valence and implicit hydrogens
perceive ring membership/ring set
perceive aromaticity
perceive stereochemistry
validate consistency
```

Not all stages need to be implemented immediately, but the API should be designed as if they will exist.

Model choices should be explicit:

```rust
ValenceModel::RdkitLike
AromaticityModel::RdkitLike
```

Aromaticity is a toolkit convention, not a single physical truth. The API should allow future models such as Daylight-like, MDL-like, OpenEye-like, or custom research models without rewriting user code.

## Canonicalization and identity

Canonicalization should live under `canon` and be exposed through small-molecule convenience methods where appropriate.

Near-term API:

```rust
molecules::canon::atom_ranking(&mol)
molecules::smiles::write_canonical(&small_molecule)
SmallMolecule::to_canonical_smiles()
```

Future API should distinguish identity modes:

- graph only
- graph plus isotope
- graph plus formal charge
- graph plus stereochemistry
- graph plus atom maps
- explicit-hydrogen-sensitive vs hydrogen-normalized
- aromatic representation-sensitive vs kekule-normalized

Do not bake one identity policy into every equality or canonicalization API.

## Error policy

Prefer domain-specific error types at module boundaries:

- `SmilesParseError`
- `SdfParseError`
- `MolWriteError`
- `SanitizeError`
- `RingPerceptionError`
- `AromaticityError`
- `BioHierarchyError`

A general crate-level error may be useful later for convenience workflows, but the prelude should not export a misleading `Result<T>` alias tied to only one error domain.

Error messages should include enough location/context for user-facing diagnostics: character offset for SMILES, record/line for SDF/Molfile, atom/bond ID for graph/perception errors.

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

## Development rules

Before adding new features, keep the API shape disciplined:

1. Do not add new blanket root re-exports.
2. Prefer focused public modules over a large flat namespace.
3. Keep the prelude small.
4. Keep parsing and sanitization separate.
5. Make convenience methods call the same namespaced functions users can call directly.
6. Keep `Molecule` chemically general and domain-neutral.
7. Keep biomolecular hierarchy out of core atoms and bonds.
8. Keep caches and perception freshness mostly internal.
9. Use explicit model enums for toolkit-convention algorithms.
10. Add tests for public API examples, not only internal helper behavior.

## Feature-driven development

Every nontrivial capability gets a directory under `features/`. `feature.toml` is the machine-readable source of truth, and `feature.md` is the human-readable source of truth. The generated dashboard is derived from feature metadata and should not be hand-edited.

Feature documents should state which API tier they affect: happy path, namespaced functional API, expert API, or internal implementation only.

## Near-term migration plan

The next API-focused milestone should happen before adding large new chemistry features.

Recommended order:

1. Replace blanket root re-exports with explicit public modules.
2. Create a small prelude.
3. Give `SmallMolecule` real methods and stop documenting `.mol` access.
4. Add namespaced facades for `smiles`, `molfile`, `sdf`, `perception`, and `canon`.
5. Keep compatibility aliases only where they are cheap and clearly marked as transitional.
6. Update README examples to use the new happy-path API.
7. Update tests so public examples compile and define the intended user experience.

This migration is allowed to break the API because the crate is still pre-release. The goal is to fix the shape now, before downstream code and additional features make it expensive.
