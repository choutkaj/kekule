# Architecture

## Purpose

This document is the normative architecture contract for `kekule`. It defines
the ownership, semantic boundaries, and invariants of the core molecular data
model. Detailed API behavior belongs in Rustdoc and tests.

`kekule` is a pure-Rust foundation for cheminformatics, structural
bioinformatics, molecular structure handling, and molecular modelling.

The foundational chemical object is `Molecule`: one non-empty connected,
geometry-independent molecular entity. The foundational system object is
`Topology`: one coordinate-free molecular system composed from one or more
molecule instances. Geometry belongs above topology.

## Canonical object model

```text
source text / bytes
    -> format-specific parsing / interpretation
    -> Vec<Molecule>
         each Molecule is one connected component
    -> optional Topology construction
         one or more Molecule instances
    -> Model      = one Topology + one geometric realization
       Ensemble   = one Topology + finite non-temporal realizations
       Trajectory = one Topology + ordered temporal realizations
```

The intended ownership hierarchy is therefore:

```text
Molecule
  one connected geometry-independent molecular entity

Topology
  one geometry-independent system made from one or more Molecule instances

Model
  one geometry-dependent realization of a Topology

Ensemble
  several non-temporal realizations of one Topology

Trajectory
  several temporally ordered realizations of one Topology
```

A salt, noncovalent complex, solvent box, protein-ligand system, or DNA duplex
is not represented by weakening `Molecule` into a disconnected graph. It is
represented by several connected molecules at the `Topology` level.

An asserted topological `Bond` contributes to molecular connectedness. Spatial
association, hydrogen bonding, ionic attraction, contact, or any other
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

These fields have deliberately different semantic roles:

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

`Molecule` contains no coordinates, conformers, velocities, periodic cell, or
other geometry-dependent state.

### Connectedness invariant

Every published `Molecule` is non-empty and connected through asserted
topological bonds. A single atom is a valid connected molecule.

A disconnected graph is never a valid `Molecule`.

Temporary disconnectedness is permitted only inside construction/edit staging
such as `MoleculeEditor`.

This is a type-level architectural invariant, not merely a convention. Public
construction and editing APIs must make it impossible to publish an invalid
disconnected molecule.

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

The current core architecture does not introduce a system-level biological
grouping or provenance framework. Such information may remain in format-layer
records or other external sidecars when needed.

### Domain-specific APIs

Protein-, nucleic-acid-, or polymer-specific algorithms do not require owning
wrapper types.

They may operate directly on `Molecule`/`Hierarchy`, or APIs may expose
lightweight borrowed validated views such as `ProteinView<'_>` or
`NucleicAcidView<'_>`.

Such views are interpretations of one `Molecule`; they do not own another
molecular object.

## `Perception`

`Perception` is the installed derived interpretation of one exact represented
molecular graph.

It replaces the architectural role formerly described by `PerceptionState`.
Within `Molecule`, `perception: Perception` is idiomatic Rust and unambiguous.

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

Graph-changing edits invalidate affected perception. The implementation should
prefer simple, safe invalidation over a complex dependency engine.

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

Editing operations should be simple graph operations; validity is enforced at
publication rather than through specialized mutation APIs.

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

A caller that wants a multi-molecule system subsequently assembles a `Topology`
from the returned molecules.

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
- positions remain geometry-dependent and must not be stored in `Molecule` or
  `Topology`.

Format-specific loaders may construct higher-level geometry objects when that is
the natural API, but the same semantic boundary remains: topology is
coordinate-free and geometry lives above it.

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

## `Topology`

`Topology` is the immutable, geometry-independent representation of one
molecular system.

Its fundamental responsibility is to answer:

> Which molecular entities exist in this system, and how are all of their
> identities laid out at system scope?

A topology contains one or more explicit `Molecule` instances. Because every
`Molecule` is connected and topology introduces no bonds between different
instances, the connected components of the topology's asserted covalent graph
are exactly its molecule instances.

Conceptually:

```text
Topology
  molecule definitions
  molecule instances
  topology-wide atom/bond identity
  canonical dense atom/bond ordering
  identity <-> dense-index mappings
```

Topology contains no positions, velocities, forces, periodic cell, conformers,
frame ordering, or other geometry-dependent state.

### Molecule instances are the public system concept

Scientifically, a topology is a system containing molecules. A caller should
normally think in terms of explicit molecule instances, not storage
normalization.

The primary public abstraction is therefore one instance-qualified molecule.
A borrowed view such as:

```rust
pub struct MoleculeInstanceView<'a> {
    topology: &'a Topology,
    id: MoleculeInstanceId,
}
```

may provide ergonomic access to:

```text
instance identity
underlying Molecule
qualified atoms
qualified bonds
qualified hierarchy
```

The exact type name is not normative, but instance-first navigation is.

Typical APIs should make ordinary iteration natural:

```rust
for molecule in topology.molecules() {
    for (atom_id, atom) in molecule.atoms() {
        // atom_id is instance-qualified
    }
}
```

Lower-level definition/instance APIs remain useful for explicit reuse and
advanced system construction.

### Definitions are a storage/reuse mechanism

Repeated identical molecule instances should not require repeated storage of the
same geometry-independent molecular definition.

Topology may therefore intern reusable `Molecule` values as definitions:

```text
MoleculeDefinition
  owns one Molecule

MoleculeInstance
  has one MoleculeInstanceId
  references one MoleculeDefinitionId
```

For example, a box containing many water molecules may store one water
`MoleculeDefinition` and many `MoleculeInstance`s.

This definition/instance split is part of the storage architecture, but it is
not the primary scientific mental model presented to ordinary callers.

A published topology must not contain unused molecule definitions. Every
`MoleculeDefinition` must be referenced by at least one `MoleculeInstance`.

The minimal core does not attach a generic `MoleculeInstanceMetadata` object to
each instance. Contextual roles, annotations, or source metadata can be added
later only when concrete use cases justify their semantics.

### Instance-qualified identity

`AtomId` and `BondId` are local to one `Molecule` definition. Once a molecule
appears in a topology, system-level identity must qualify the local ID by its
molecule instance.

Conceptually:

```text
InstanceAtomId = (MoleculeInstanceId, AtomId)
InstanceBondId = (MoleculeInstanceId, BondId)
```

The same rule applies to hierarchy-local identities when exposed at topology
scope.

This prevents identity collisions when one definition is instantiated multiple
times.

### Dense topology ordering

Numerical and geometry-bearing state requires a deterministic dense ordering over
the complete system.

Topology therefore owns an authoritative dense atom order and, where useful, a
dense bond order:

```text
InstanceAtomId <-> TopologyAtomIndex
InstanceBondId <-> TopologyBondIndex
```

These concepts have deliberately different roles:

```text
InstanceAtomId / InstanceBondId
  semantic system identity

TopologyAtomIndex / TopologyBondIndex
  dense storage position
```

Dense ordering tells `Model`, `Ensemble`, and `Trajectory` how to interpret their
numerical arrays. The arrays themselves do not own topology identity.

### No topology-level covalent bonds

Topology must not introduce asserted covalent/topological bonds between
molecule instances.

If atoms from two current molecule instances become connected by an asserted
bond, those atoms belong to one connected `Molecule`. The resulting system must
therefore be represented by a new topology containing the newly connected
molecule rather than by adding an inter-instance bond.

Hydrogen bonds, salt bridges, contacts, coordination hypotheses, force-field
interactions, and other spatial or energetic relations are not topology bonds.

### No connected-components API per instance

Because every published `Molecule` is connected, asking for connected components
inside one molecule instance is redundant: the answer is always exactly that
instance.

Topology should therefore expose molecule-instance membership directly rather
than retain an API such as `connected_components(instance)` whose result is
architecturally predetermined.

### Construction and invariants

A topology builder may stage definitions and instances and publish an immutable
`Topology` only after validation.

A published topology must satisfy at least:

```text
at least one molecule instance
every instance references a live definition
every definition is referenced by at least one instance
every referenced Molecule satisfies Molecule invariants
instance-qualified atom/bond identities are valid
dense atom/bond ordering is complete and deterministic
identity/index mappings are mutually consistent
```

Convenience construction may add one fresh definition and one instance in a
single operation. Explicit APIs may separately add a reusable definition and
then instantiate it many times.

### Immutability and topology changes

Published `Topology` is immutable.

Shared exact ownership should use `Arc<Topology>` rather than cloning independent
copies of topology state.

A chemical or structural transformation that changes molecule membership,
connectivity, atom count, bond count, hierarchy identity, or dense layout
produces a new `Topology` rather than mutating an existing topology underneath
geometry-bearing state.

The core architecture does not provide a generic topology-remapping framework.
If a workflow changes topology, geometry or other dense state for the new system
must be constructed explicitly according to that workflow's own semantics.
Topology correspondence, when needed by a specialized algorithm, is a separate
algorithmic result rather than a foundational ownership mechanism.

### Scope of the current Topology design

The current core intentionally remains minimal.

It does not introduce:

```text
generic molecule-instance metadata
system-level biological grouping
system-level provenance hierarchy
geometry-dependent interactions
inter-molecule topology bonds
generic topology remapping
```

These concerns should not be added speculatively. They may be introduced later
only as separate concepts when concrete requirements establish their semantics.

## Geometry boundary

Geometry starts above `Topology`.

The fundamental relationship is:

```text
Topology + Positions -> Model
```

`Positions` and the other dense realization arrays are numerical storage. They
do not own or retain `Arc<Topology>` and they do not resolve semantic atom or
bond identities themselves.

The intended separation is:

```text
Topology
  defines semantic identities and dense atom/bond ordering

Positions / AtomData / BondData / Velocities / Forces
  store dense numerical state
  know their own shape/length and numerical units as appropriate
  do not own Topology
  do not carry topology identity

Model / Ensemble / Trajectory
  own the Topology context exactly once
  validate dense state against that Topology
  provide semantic atom/bond access when topology context is required
```

This means operations such as resolving an `InstanceAtomId` to a coordinate are
operations on `Model`, `Ensemble` member/frame views, or `Trajectory` frame
views, not operations on a detached `Positions` array.

Similarly, topology compatibility is an invariant of the owning aggregate. It
is not established by storing repeated `Arc<Topology>` handles inside every
numerical subobject.

Geometry-dependent quantities such as positions, velocities, forces, periodic
cell, occupancies, B-factors, and other model/frame state do not belong in
`Molecule` or `Topology`.

### Dense numerical containers

`Positions`, `Velocities`, and `Forces` are dense numerical arrays in canonical
model units. They validate numerical shape, units, and finite values as
appropriate, but are otherwise topology-agnostic.

`AtomData` and `BondData` are likewise dense data containers rather than
Topology owners. They may retain a logical item count even when all optional
columns are absent. Whether particular atom/bond data belongs at model, member,
frame, or another level is intentionally not settled by this storage rule and
may be refined separately.

The primitive dense containers must not expose APIs that require a `Topology`
parameter merely to translate semantic IDs. Semantic navigation belongs to the
higher-level object that owns both topology and dense state.

## `Model`

`Model` is one concrete geometry-dependent realization of one topology.

Conceptually:

```text
Model
  shared Topology          <- owned once
  Positions
  optional periodic cell
  model-level AtomData
  model-level BondData
```

A `Model` does not duplicate molecular chemistry. It interprets dense realization
state against its topology's authoritative dense layout.

Construction validates at least:

```text
Positions length == Topology atom count
AtomData logical length == Topology atom count
BondData logical length == Topology bond count
```

After construction, public mutation APIs must preserve those dimensional
invariants.

Semantic operations such as `position(InstanceAtomId)` belong on `Model` or a
model-level borrowed view because only that layer owns both the topology and the
numerical state.

## `Ensemble`

`Ensemble` is a finite collection of non-temporal realizations of one topology.

Conceptually:

```text
Ensemble
  shared Topology          <- owned once
  members[]
    Positions
    optional member-level geometry/data
    optional weight
```

Members do not own or repeat the shared topology. Their dense state is interpreted
in the `Ensemble` topology's authoritative order.

Insertion/construction validates every member's dimensions against the ensemble
topology. Differences between members are geometric or member-level data, not
molecular identity.

An ensemble weight is contextual to membership in that ensemble and therefore
belongs to the member relation rather than to `Topology`.

## `Trajectory`

`Trajectory` is an ordered temporal sequence of realizations of one topology.

Conceptually:

```text
Trajectory
  shared Topology          <- owned once
  ordered frames[]
    Positions
    optional periodic cell
    optional AtomData / BondData
    optional Velocities / Forces
    optional time / step
    optional frame properties
```

Frames do not own or repeat the shared topology. Their dense arrays are
interpreted in the `Trajectory` topology's authoritative order.

Insertion, decoding, and reusable frame-buffer publication validate frame
shapes against the trajectory topology. Streaming infrastructure may use
reusable buffers for allocation efficiency, but those buffers follow the same
ownership rule: topology context is owned once by the buffer/container rather
than repeated inside every numerical array.

A trajectory represents one fixed-topology epoch. Topology-changing chemistry is
not represented by silently mutating one shared topology or by generically
remapping frames. A workflow with changing topology should use separate topology
epochs/objects and explicitly construct the geometry belonging to each epoch.

## Molecular identity and equality

Authoritative molecular identity is defined by represented state, not by derived
cache population.

`Perception` must therefore not make two otherwise identical represented
molecules unequal merely because one has different cache presence.

Whether hierarchy participates in a particular equality/hash/canonicalization
operation must be explicit in that operation's semantics. Chemical graph
identity and full represented-object identity need not be forced into one
ambiguous notion of equality.

Topology layout equality is distinct from graph isomorphism or chemical
identity. Two independently constructed topologies may represent chemically
equivalent systems while still having different instance IDs or dense layouts.

## Persistence and reconstruction

Persistence consumers may store graph, hierarchy, and perception separately.

Molecule reconstruction order is:

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

Topology persistence must reconstruct definitions, instances, qualified
identities, and authoritative dense ordering consistently. Geometry is restored
separately and validated by the owning `Model`, `Ensemble`, or `Trajectory`
against that topology layout.

Runtime domain objects are not required to be generic file-format DTOs. Source
metadata that is not canonical represented molecular or topology state should
remain in format records or other external sidecars.

## Mutation and transformations

A normal molecular edit returns one valid connected `Molecule`.

Operations whose semantic purpose is to split a molecule naturally return more
than one molecule, for example a fragmentation transformation may return
`Vec<Molecule>`.

Topology-changing system operations return a new topology rather than mutating a
published topology in place. They do not automatically remap existing dense
geometry/data into the new topology.

Coordinate-only operations never mutate `Graph`, `Hierarchy`, `Perception`, or
`Topology`.

## Naming and module style

The intended molecule field/type naming is idiomatic Rust:

```rust
pub struct Molecule {
    graph: Graph,
    hierarchy: Hierarchy,
    perception: Perception,
}
```

Field names use `snake_case`; type names use `UpperCamelCase`. Patterns such as
`graph: Graph`, `hierarchy: Hierarchy`, and `perception: Perception` are normal
Rust style and are preferred over redundant names unless a real ambiguity
appears.

The exact file/module layout is not normative; semantic boundaries are.

## Design rules

When deciding where new state belongs:

1. Is it authoritative atom/bond/stereo chemistry? Put it in `Graph`.
2. Is it coordinate-independent residue/chain/polymer organization within one
   connected molecule? Put it in `Hierarchy`.
3. Is it fundamental chemistry derived from the represented graph? Put it in
   `Perception`.
4. Does it identify which connected molecules exist in one coordinate-free
   system or define their topology-wide layout? Put it in `Topology`.
5. Is it task-specific analysis, typing, scoring, or parameterization? Keep it
   in a separate derived object.
6. Is it dense coordinate/model/frame data? Store it in a topology-agnostic
   numerical container above `Topology`.
7. Does an operation need to interpret dense data by semantic atom/bond identity?
   Perform it at the `Model`, `Ensemble`, or `Trajectory` level where topology
   and numerical state meet.
8. Does an asserted new bond connect two current molecule instances? Construct a
   new connected `Molecule` and therefore a new `Topology`.
9. Does a workflow change topology? Construct the new topology and its new dense
   state explicitly; do not rely on a generic remapping layer.

The core invariants are intentionally simple:

> A Kekule `Molecule` is one connected, geometry-independent molecular entity
> represented by `Graph + Hierarchy + Perception`.

> A Kekule `Topology` is one immutable, geometry-independent molecular system
> composed of one or more explicit `Molecule` instances with authoritative
> topology-wide identity and dense layout.

> `Model`, `Ensemble`, and `Trajectory` each own their shared `Topology` exactly
> once. Their dense numerical subobjects do not own topology identity.