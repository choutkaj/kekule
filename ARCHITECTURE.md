# Architecture

## Purpose

This document is the normative architecture contract for `kekule`. It defines
the ownership, semantic boundaries, and invariants of the core molecular data
model. Detailed API behavior belongs in Rustdoc and tests.

`kekule` is a pure-Rust foundation for cheminformatics, structural
bioinformatics, molecular structure handling, and molecular modelling.

The foundational chemical object is `Molecule`: one non-empty connected,
geometry-independent chemical graph. The foundational system object is
`Topology`: one coordinate-free molecular system composed from one or more
molecule instances. Geometry belongs above topology.

## Canonical object model

```text
source text / bytes
    -> format-specific parsing / interpretation
    -> Vec<Molecule>
         each Molecule is one connected component
    -> optional Topology construction
         one or more molecule instances
    -> Model      = Topology + Positions + model-level data
       Ensemble   = one Topology + finite non-temporal position sets
       Trajectory = one Topology + ordered position frames
```

The molecule refactor defined here intentionally does not redesign `Topology`,
`Model`, `Ensemble`, or `Trajectory`. Those layers retain their present
geometry/system responsibilities. They may receive only the minimal API
adaptations needed to consume the new universal `Molecule`.

## Semantic layers

Kekule separates three levels:

```text
Molecule
  one connected, geometry-independent molecular entity

Topology
  one geometry-independent system made from one or more Molecule instances

Model
  one geometry-dependent realization of a Topology
```

A salt, noncovalent complex, solvent box, protein-ligand system, or DNA duplex
is therefore not represented by weakening `Molecule` into a disconnected graph.
It is represented by several connected molecules at the system level.

An asserted topological `Bond` contributes to molecular connectedness.
Spatial association, hydrogen bonding, ionic attraction, contact, or any other
non-topological interaction does not.

Examples:

```text
ethanol                         -> 1 Molecule
ubiquitin                       -> 1 Molecule
Na+ + acetate-                  -> 2 Molecules
protein + noncovalent ligand    -> 2 Molecules
covalent protein-ligand adduct  -> 1 Molecule
DNA duplex                      -> typically 2 Molecules
disulfide-linked protein chains -> 1 Molecule
```

## `Molecule`

There is exactly one foundational molecular type. Kekule does not distinguish
`SmallMolecule` and `MacroMolecule` as owning wrappers.

The intended core shape is:

```rust
pub struct Molecule {
    graph: Graph,
    hierarchy: Hierarchy,
    perception: Perception,
}
```

These three fields have deliberately different semantic roles:

```text
Graph
  authoritative represented chemistry
  required
  defines the connected molecular entity

Hierarchy
  authoritative coordinate-independent organization
  may be empty
  residues, chains, polymer annotation, and mappings to graph atoms

Perception
  derived chemical interpretation
  reconstructible from represented chemistry plus an explicit perception model
  does not define molecular identity
```

`Molecule` contains no coordinates, conformers, velocities, unit cell, or other
geometry-dependent state.

### Connectedness invariant

Every published `Molecule` is non-empty and connected through asserted
topological bonds. A single atom is a valid connected molecule.

A disconnected graph is never a valid `Molecule`.

Temporary disconnectedness is permitted only inside construction/edit staging
such as `MoleculeEditor`.

This is a type-level architectural invariant, not merely a convention.
Public construction and editing APIs must make it impossible to publish an
invalid disconnected molecule.

## `Graph`

`Graph` is the authoritative chemical graph of a molecule.

Conceptually it owns:

```text
Graph
  atoms
  bonds
  adjacency/connectivity
  stable local AtomId / BondId identity
  represented stereochemical elements/groups
```

The exact physical storage is an implementation detail, but it should remain
purpose-built for chemistry, compact, deterministic, and efficient.

### Represented atom chemistry

Fundamental asserted atom state belongs in the graph, for example:

```text
element
isotope
formal charge
radical state
represented hydrogen declaration
atom-map identity when retained as canonical represented chemistry
```

### Represented bond chemistry

Fundamental asserted bond state belongs in the graph:

```text
endpoints
localized represented bond order/kind
represented bond stereochemistry where applicable
```

A `Bond` is a topological relation. A noncovalent interaction must not be added
as a `Bond` merely because it relates two atoms spatially or energetically.

### Stereochemistry

Canonical represented stereochemical state belongs to the graph because it is
part of the asserted molecular representation.

Source-format marks such as SMILES directional syntax or molfile wedges are
format/interpreter state. They must be resolved into Kekule's canonical stereo
representation before a molecule is published.

CIP labels are derived and therefore belong to `Perception`, not `Graph`.

### Aromaticity

Aromaticity is perceived chemistry, not a canonical bond order.

Canonical graph bonding is localized. Aromatic source syntax may be accepted by
format readers/interpreters, but a published `Graph` contains ordinary
represented bond orders. Aromatic atom/bond membership belongs to `Perception`.

## `Hierarchy`

`Hierarchy` stores coordinate-independent organization attached to the same
connected molecular entity.

It may be empty. A molecule does not become a different Rust type merely
because hierarchy is present.

Typical hierarchy state includes:

```text
residues
chains
polymer organization
residue/chain identifiers
component names
atom-site annotations
mappings from hierarchical entities to AtomId
```

Hierarchy is orthogonal to chemistry:

```text
Graph answers:
  which atoms exist and how are they chemically bonded?

Hierarchy answers:
  how are those atoms organized into residues, chains, and polymers?
```

Hierarchy does not own independent atoms or bonds and must reference only live
graph atoms.

A small molecule parsed from a biological format may carry hierarchy. A protein
or nucleic acid uses the same `Molecule` type and simply has non-empty
hierarchy.

The presence or absence of hierarchy is not a classifier and must not control
which fundamental molecular type exists.

### Disconnected biological source entities

Biological/source identity may span disconnected graph components, for example
a chain with unresolved missing residues or a noncovalent multi-chain complex.
Kekule must not fabricate bonds to preserve that source grouping.

Each connected component becomes a separate `Molecule`, carrying the subset of
hierarchy that belongs to that component.

Cross-component source provenance or future system-level biological grouping
may be represented above `Molecule`, but this refactor does not redesign
`Topology` to add such grouping.

### Domain-specific APIs

Protein-, nucleic-acid-, or polymer-specific algorithms do not require owning
wrapper types.

They may operate directly on `Molecule`/`Hierarchy`, or future APIs may expose
lightweight borrowed validated views such as `ProteinView<'_>` or
`NucleicAcidView<'_>`.

Such views are interpretations of one `Molecule`; they do not own another
molecular object and are not required for the initial refactor.

## `Perception`

`Perception` is the installed derived interpretation of one exact represented
molecular graph.

It replaces the architectural role currently described by
`PerceptionState`. The shorter name is intentional: within `Molecule`,
`perception: Perception` is idiomatic Rust and unambiguous.

Perception is semantically subordinate to represented chemistry:

```text
Graph
  + perception model/policy
        |
        v
Perception
```

Deleting all perception state must never destroy authoritative information
about the molecule.

Two molecules must not become chemically different merely because one has more
derived perception cached than the other.

Fundamental perception may include sections such as:

```text
valence / implicit hydrogens
ring/cycle information
aromaticity
CIP assignments
```

Not every calculated property belongs in `Perception`. Fingerprints,
descriptors, partial charges, force-field types, pharmacophore features,
rotatable-bond classifications, scoring terms, and other task-specific results
remain separate derived objects unless explicitly promoted by the architecture.

### Structural graph derivations versus chemical perception

Kekule may internally distinguish mathematically unique graph caches from
model-dependent chemical perception.

For example:

```text
graph-derived:
  degree
  connected traversal data
  connected components during editing
  generic cycle membership

chemically perceived:
  valence model
  implicit hydrogens
  aromaticity model
  CIP assignment
```

This distinction may be reflected internally if useful, but the public
`Molecule` architecture remains `graph + hierarchy + perception`.

### Perception installation and invalidation

Perception must always correspond to the current authoritative graph.

Graph-changing edits invalidate affected perception. The initial implementation
should prefer simple, safe invalidation over a complex dependency engine.

Exact reconstruction of externally stored perception may be supported through
checked installation APIs, but installation must validate references and
dimensions and must never rewrite authoritative graph chemistry.

Public convenience APIs may expose perception-backed queries directly through
`Molecule`; callers should not need to duplicate perceived flags into `Atom` or
`Bond`.

## Editing

All structural mutation of `Molecule` happens through `MoleculeEditor` or an
equivalent transactional staging type.

Conceptually:

```text
Molecule
   |
   | edit
   v
MoleculeEditor
   |
   | arbitrary intermediate structural edits
   | temporary disconnection is allowed
   | temporary incomplete chemistry is allowed
   v
finish()
   |
   +-- invalid -> error
   |
   `-- valid -> Molecule
```

`MoleculeEditor` is allowed to violate publication invariants while editing.
The finished `Molecule` is not.

This deliberately avoids APIs such as
`remove_bond_preserving_connectivity(...)`. Editing operations should be simple
graph operations; validity is enforced at publication.

`finish()` must at minimum validate:

```text
non-empty graph
exactly one connected component
valid atom/bond references
valid adjacency
valid stereo references
valid hierarchy references and hierarchy/graph consistency
```

Chemical perception need not be valid during editing. On successful
publication, stale perception must be discarded, recomputed, or explicitly
reinstalled through a checked path.

The same editor concept may be used for construction from scratch; a separate
public `MoleculeBuilder` is not architecturally required unless it provides
clear ergonomic value without duplicating semantics.

Public unrestricted mutable access to graph internals should not bypass the
editor and thereby bypass publication validation.

## Parsing and interpretation

### Component output

The canonical molecule-producing result for one molecular record is:

```rust
Result<Vec<Molecule>>
```

Every returned element is one valid connected molecule.

Disconnected source syntax is partitioned rather than represented as a
disconnected `Molecule`.

Examples:

```text
"CCO"               -> [ethanol]
"[Na+].[Cl-]"       -> [sodium, chloride]
"CC(=O)[O-].[Na+]"  -> [acetate, sodium]
```

Component order follows deterministic source/interpreter order. Element zero
does not carry a semantic guarantee that it is the chemically "main" component.
A caller may choose the first component if that is its desired policy, or apply
an explicit largest/organic/main-component policy separately.

Kekule should not introduce direct parser-to-`Topology` construction as part of
this refactor. A caller that wants a system may subsequently assemble a
`Topology` from the returned molecules.

### Parsing versus source documents

Format parsers may internally or publicly retain format-specific `Document` or
record representations when required for faithful syntax, metadata, coordinates,
or streaming. The architecture requirement is that the canonical
chemistry-producing boundary publishes connected `Vec<Molecule>` values, never
a disconnected molecule.

Multi-record formats should preserve record boundaries, for example as a
stream/iterator whose per-record molecular result is `Vec<Molecule>` rather than
flattening every component of an entire file into one undifferentiated vector.

### Coordinate-bearing formats

Coordinate-bearing formats create an important separation:

- chemistry extracted from a record becomes geometry-independent `Molecule`
  values;
- positions remain geometry-dependent and must not be stored in `Molecule`.

Existing or future APIs that construct `Model`/`Ensemble`/`Trajectory` from
coordinate-bearing formats are outside the scope of this molecule refactor.
The format layer may preserve coordinates until such a geometry-producing path
uses them.

### Interpretation

Parsing recognizes source syntax. Interpretation translates source assertions
into canonical Kekule graph/hierarchy state.

Interpretation may perform deterministic representation rewrites required to
publish a canonical molecule, such as localization of aromatic source bonding
and conversion of source stereo notation into canonical stereo elements.

Interpretation does not run arbitrary chemical standardization, choose a
tautomer/protonation state, or invent bonds merely to force connectedness.

If interpretation yields multiple disconnected components, each component is
published independently as a valid `Molecule`.

## Geometry boundary

`Molecule` is strictly geometry-independent.

Coordinates, conformers, velocities, forces, periodic cells, occupancies,
B-factors, and other coordinate/model state do not belong in `Molecule`.

The existing higher-level separation remains:

```text
Topology + Positions -> Model
```

and:

```text
one Topology + many position sets -> Ensemble / Trajectory
```

This refactor must remove local conformer ownership from `Molecule`.

## `Topology`

`Topology` remains the geometry-independent system layer.

It can contain one or more molecule instances and qualifies their local
identities at system scope.

This refactor is not intended to redesign topology storage, instance semantics,
model binding, ensembles, or trajectories.

However, because the current implementation accepts `SmallMolecule` and
`MacroMolecule`, minimal mechanical changes to `Topology` are expected so it
consumes the new universal `Molecule` directly.

Those changes must preserve existing topology semantics rather than expand the
scope of this refactor.

## Molecular identity and equality

Authoritative molecular identity is defined by represented state, not by derived
cache population.

`Perception` must therefore not make two otherwise identical represented
molecules unequal merely because one has different cache presence.

Whether hierarchy participates in a particular equality/hash/canonicalization
operation must be explicit in that operation's semantics. Chemical graph
identity and full represented-object identity need not be forced into one
ambiguous notion of equality.

## Persistence and reconstruction

Persistence consumers may store graph, hierarchy, and perception separately.

Reconstruction order is:

```text
Graph
  -> validate represented graph
Hierarchy
  -> validate against Graph
Perception
  -> checked install last
Molecule
```

Persisted disconnected graph data must be partitioned into connected molecules
or rejected before publication.

Loading must never weaken the connectedness invariant.

Runtime domain objects are not required to be generic file-format DTOs. Source
metadata that is not canonical represented molecular state should remain in
format records, provenance objects, or other sidecars.

## Mutation and transformations

A normal molecular edit returns one valid connected `Molecule`.

Operations whose semantic purpose is to split a molecule naturally return more
than one molecule, for example a fragmentation transformation may return
`Vec<Molecule>`.

Topology-changing system operations remain system-level transformations.

Coordinate-only operations never mutate `Graph`, `Hierarchy`, or discrete
geometry-independent `Perception`.

## Naming and module style

The intended field/type naming is idiomatic Rust:

```rust
pub struct Molecule {
    graph: Graph,
    hierarchy: Hierarchy,
    perception: Perception,
}
```

Field names use `snake_case`; type names use `UpperCamelCase`. Patterns such as
`graph: Graph`, `hierarchy: Hierarchy`, and `perception: Perception` are normal
Rust style and are preferred over redundant names such as
`molecular_graph: MolecularGraph` unless a real ambiguity appears.

A natural module layout is:

```text
core/
  atom_bond.rs
  graph.rs
  hierarchy.rs
  perception.rs
  molecule.rs
  molecule_edit.rs
  stereo.rs
```

The exact file layout is not normative; the semantic boundaries are.

## Scope of the molecule refactor

The refactor implementing this document should:

- replace the current `Molecule` implementation with the
  `Graph + Hierarchy + Perception` architecture;
- enforce non-empty connected published molecules;
- make `MoleculeEditor` the transactional structural mutation boundary;
- remove owning `SmallMolecule` and `MacroMolecule` wrappers;
- integrate hierarchy directly into universal `Molecule`;
- rename/reframe `PerceptionState` as `Perception` while preserving useful
  perception functionality;
- remove geometry/conformer storage from `Molecule`;
- make molecule-producing parsing/interpretation yield connected
  `Vec<Molecule>` components;
- adapt downstream molecule consumers as required to compile and preserve
  behavior.

The refactor should not redesign `Topology`, `Model`, `Ensemble`, or
`Trajectory`. Touch those layers only where the removal of old wrapper APIs
requires a mechanical compatibility update.

## Design rules

When deciding where new state belongs:

1. Is it authoritative atom/bond/stereo chemistry? Put it in `Graph`.
2. Is it coordinate-independent residue/chain/polymer organization? Put it in
   `Hierarchy`.
3. Is it fundamental chemistry derived from the represented graph? Put it in
   `Perception`.
4. Is it task-specific analysis, typing, scoring, or parameterization? Keep it
   in a separate derived object.
5. Is it coordinate or model state? Keep it outside `Molecule`.
6. Does it combine several connected molecular entities? It belongs at or above
   `Topology`, not in a disconnected `Molecule`.

The core invariant is intentionally simple:

> A Kekule `Molecule` is one connected, geometry-independent molecular entity
> represented by `Graph + Hierarchy + Perception`.
