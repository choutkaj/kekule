# Architecture

## Purpose

This document is the normative architecture contract for `kekule`. It defines
object ownership, semantic boundaries, and invariants. Detailed API behavior
belongs in Rustdoc and tests; implementation history does not belong here.

`kekule` is a pure-Rust foundation for cheminformatics, structural
bioinformatics, molecular structure handling, and molecular modelling. It does
not make a file-format record, a simulation-engine particle list, or one
flattened graph the universal data model.

The foundational chemical object is `Molecule`: one connected represented
chemical entity. The foundational system object is `Topology`: one immutable,
coordinate-free molecular system composed of one or more connected molecule
instances.

## Canonical object model

```text
source text / bytes
    -> format-specific Document or streaming decoder
    -> explicit interpretation + provenance/report
    -> represented connected molecular objects
       Molecule
         |- SmallMolecule
         `- MacroMolecule
    -> optional normalization / perception / stereo derivation
    -> Topology
         |- Model      = Topology + Positions + optional cell + AtomData + BondData
         |- Ensemble   = Topology + finite non-temporal members
         `- Trajectory = Topology + ordered frames        [kekule-traj]
    -> analysis, topology transformation, or explicit backend preparation
```

`Topology` is shared by models, ensembles, trajectories, compiled selections,
prepared potentials, and analyses. Coordinate updates never modify it.

## Semantic boundaries

These boundaries are architectural requirements:

- **Parsing** recognizes syntax and preserves format information. It does not
  infer chemistry.
- **Interpretation** translates source-asserted semantics into Kekule domain
  objects plus mappings/provenance. It does not run general chemical
  perception or invent bonds merely to force connectedness.
- **Normalization** deterministically rewrites represented chemistry into
  Kekule's canonical equivalent representation. It is idempotent,
  model-independent, and meaning-preserving.
- **Perception** derives chemical meaning from normalized represented chemistry
  and installs a derived view without rewriting primary represented chemistry.
- **CIP assignment** and coordinate-derived stereo are specialized explicit
  derivations, not parsing, interpretation, or normalization.
- **Standardization**, if added, may choose among chemically distinct related
  states such as tautomers or protonation states and is therefore separate from
  normalization and perception.
- **Topology construction** assembles already represented connected molecules
  into one coordinate-free system. It does not invent coordinates or
  force-field state.
- **Coordinate state** is complete and bound to one exact shared topology
  allocation.
- **Analysis** is read-only unless explicitly named as a transformation.
- **Topology-changing operations** return a new topology and an explicit
  mapping. They do not mutate existing models or prepared systems.
- **Backend preparation** creates downstream topology-bound objects and never
  mutates canonical chemistry or coordinate state.
- **Writers** reject unsupported semantics rather than silently dropping or
  coercing them.

High-level conveniences may compose several stages, but composition must not
blur their meanings.

## Molecular chemistry

### `Molecule`

`Molecule` is the canonical boundary for one connected represented chemical
graph. Every non-empty published `Molecule` has exactly one connected atom/bond
component; a single atom is valid. Empty state is permitted as a construction
boundary but is not a valid topology definition.

A `Molecule` owns:

- stable local `AtomId` and `BondId` values;
- represented atoms, bonds, adjacency, stereo elements/groups, and source
  stereo marks;
- optional local source conformers;
- scalar annotations;
- installed derived `PerceptionState`.

Temporary disconnection belongs only to checked construction/editing state.
`MoleculeBuilder` and `MoleculeEditor` may work on incomplete candidates, but
publication is transactional and enforces connectedness. Stable local IDs are
not silently reused after deletion; fixed-width identifier capacity is checked
before mutation.

Disconnected chemistry is represented above this boundary. `[Na+].[Cl-]`, a
protein-ligand complex, solvents, ion pairs, and other noncovalent assemblies
are multiple connected molecules in one topology. An asserted Kekule bond
contributes to connectedness; spatial association does not.

### Represented versus derived chemistry

Primary represented chemistry is stored directly and includes, as applicable:

```text
Atom: element, isotope, formal charge, radical, explicit-H declaration,
      no-implicit-H declaration, atom map
Bond: endpoints, represented BondOrder
Stereo: local StereoElement state, groups, source stereo marks
```

Derived chemistry is a view over that representation:

```text
valence / implicit hydrogens
ring membership and ring basis
aromaticity membership/model
CIP descriptors
```

The initial discrete derived view is `PerceptionState`. It is semantic runtime
state, not an algorithm report. Diagnostics, work counters, candidate lists,
and warnings are sidecars and do not become molecular identity merely because
they were produced during perception.

The architecture does not assume a single possible perception kind. A future
continuous or learned perception state may coexist with the discrete view as a
distinct derived representation; no generic perception registry is required
before such a state exists.

Chemistry- or connectivity-relevant mutation invalidates affected perception.
Coordinate and generic annotation edits are perception-neutral. Operations that
can fail publish new semantic state only after successful validation.

Exact external reconstruction of installed perception is supported. A detached
`PerceptionState` may represent absent or model-neutral valence, complete
implicit-H assignments, ring membership with or without a stored basis,
aromaticity model/membership, and CIP assignments. Checked whole-state
installation validates slot dimensions and live graph/stereo references before
replacement; it does not normalize or re-perceive the molecule.

### `SmallMolecule`

`SmallMolecule` is the ordinary cheminformatics wrapper around one connected
`Molecule`. It may provide ergonomic parse/normalize/perceive workflows, but it
never represents several disconnected salt or mixture components.

### `MacroMolecule` and SMCRA

`MacroMolecule` is one connected `Molecule` plus one validated
`SmcraHierarchy`. The hierarchy stores coordinate-independent structural
identity: chains, residues, atom sites, label/author identifiers, insertion
codes, component names, and mappings to local `AtomId` values.

SMCRA IDs are definition-local. When a reusable macromolecule definition occurs
in a topology, `Topology` qualifies chain/residue/atom-site identity by
`MoleculeInstanceId` and exposes borrowed instance-qualified views.

Coordinate-model IDs, alternate-location choices, occupancy, B-factors, and raw
coordinate text are not hierarchy topology. They belong to interpretation
provenance or model-level data.

Experimental source identity may span several represented molecule instances.
If an observed polymer has a genuine unresolved gap, Kekule does not fabricate
a bond across it: residual connected fragments become separate
`MacroMolecule` instances while shared source entity/chain provenance is
retained.

### Local conformers

Standalone molecule objects may retain local conformers for cheminformatics and
I/O workflows. A molecule definition inserted into `Topology` is
coordinate-free; source conformers are neither copied into topology nor
silently destroyed.

Explicit conversions include:

```text
SmallMolecule + selected conformer   -> Topology + Model
SmallMolecule + selected conformers  -> Topology + Ensemble
several molecules + system positions -> Topology + Model
```

## Chemistry pipeline

The small-molecule semantic pipeline is:

```text
source
  -> parse
  -> format-specific Document
  -> interpret
  -> represented Molecule
  -> normalize
  -> normalized Molecule
  -> perceive
  -> derived PerceptionState
  -> optional coordinate-stereo materialization / CIP assignment
```

Interpretation answers **what the source asserted**. Normalization answers
**how Kekule canonically represents that assertion**. Perception answers
**what chemical interpretation Kekule derives from the normalized graph**.

Normalization has a deliberately narrow contract. A normalized molecule has:

- ordinary localized represented bond orders; no remaining
  `BondOrder::Aromatic` source representation;
- no source `StereoBondMark` state that still needs conversion;
- canonical represented local `StereoElement` state;
- canonical represented charge/hydrogen declarations;
- empty/default derived `PerceptionState`.

Normalization always clears derived perception because it rewrites represented
state. It does not choose an aromaticity model, tautomer, protonation state, or
salt/fragment policy.

The default discrete perception order is:

```text
valence -> rings -> aromaticity
```

CIP assignment is opt-in after the relevant represented stereo and perception
state are available. Geometry-derived stereo is similarly explicit and must not
silently redefine source-asserted represented stereo.

## Persistence and exact reconstruction

Kekule owns invariant-preserving reconstruction of Kekule runtime objects;
persistence consumers own their versioned archive DTOs. Runtime domain objects
are not a generic serialization schema.

Reconstruction builds represented graph/hierarchy/stereo state first and
installs derived perception last. Stable stereo-group slot layout, including
tombstones, is reconstructible without dummy chemistry. Loading must not
normalize, re-perceive, renumber, or silently coerce malformed historical
state.

If persisted data describes several disconnected graph components, it must be
explicitly partitioned or rejected before publication as canonical
`Molecule`s. Persistence never weakens connectedness.

Process-local `Arc<Topology>` sharing is not serialized. Independently rebuilt
topologies may compare equal by complete static layout while remaining distinct
allocations.

## `Topology`

### Definition and ownership

`Topology` is one immutable coordinate-free molecular system. It owns:

- connected molecule definitions;
- molecule instances referencing those definitions;
- instance-qualified atoms and bonds;
- static chemistry, hierarchy, roles, and instance annotations;
- one authoritative dense atom ordering and one dense bond ordering;
- mappings between semantic identities and dense indices.

It does **not** own positions, cells, velocities, forces, time, occupancy,
B-factors, coordinate-derived analyses, force-field parameters, virtual sites,
Drude particles, constraints, backend particles, or execution-engine state.

Raw `Topology` is not cloned as an ordinary value. Shared owners store
`Arc<Topology>`.

### Definitions and instances

Topology separates reusable molecular identity from one occurrence in a system:

```text
MoleculeDefinition
    SmallMolecule | MacroMolecule

MoleculeInstance
    MoleculeInstanceId
    -> MoleculeDefinitionId
    + roles / static instance metadata
```

Definitions are connected and conformer-free. Definition reuse is explicit; a
single water definition may be referenced by many instances. Kekule does not
automatically merge definitions merely because their graphs look equal.

Every topology atom belongs to exactly one molecule instance. Asserted covalent
connectivity belongs within one resulting connected instance. Noncovalent
association does not merge instances.

### Semantic identity and dense order

Local chemical identity remains `AtomId` / `BondId`. System identity is
qualified:

```text
MoleculeDefinitionId
MoleculeInstanceId
InstanceAtomId      = molecule instance + local AtomId
InstanceBondId      = molecule instance + local BondId
InstanceChainId     = molecule instance + local chain ID
InstanceResidueId   = molecule instance + local residue ID
InstanceAtomSiteId  = molecule instance + local atom-site ID
```

Dense storage uses separate numerical indices:

```text
TopologyAtomIndex
TopologyBondIndex
```

Semantic IDs answer **what object is this?** Dense indices answer **where is its
state stored?** The dense order is authoritative and immutable for the lifetime
of a topology. Dense arrays contain live atoms/bonds only; local stable IDs may
retain tombstone positions.

### Compatibility and layout equality

Exact shared-allocation identity is the compatibility criterion for
Topology-bound state:

```rust
Arc::ptr_eq(&a, &b)
```

Clones of one `Arc<Topology>` are compatible. Independently constructed
Topologies are not automatically compatible even when chemically equivalent.

`Topology::same_layout` is a separate explicit static-layout comparison. It may
compare chemistry, hierarchy, definition/instance partitioning, metadata,
semantic IDs, dense order, and index maps, but it does not establish shared
allocation or silently authorize reuse of topology-bound arrays.

Transfer between independent topologies requires an explicit validated mapping.
General order-independent structural equivalence/isomorphism may be added
separately; ambiguous mappings must not be chosen silently.

### Topology transformations

A built topology is immutable. Connectivity- or membership-changing operations
return a new `Arc<Topology>` plus explicit `TopologyMapping` lineage.

Examples include adding/removing hydrogens, deleting atoms or bonds, changing
asserted bond order, merging/splitting molecule instances, solvation, and
chemical reactions. A transformation that disconnects a molecule must either
produce several connected molecule instances or reject the edit.

Existing models, ensembles, trajectories, selections, and prepared systems are
never mutated by a topology edit. Coordinate/data remapping is an explicit
separate operation and cannot invent state for unmapped added atoms.

Any backend workflow that may make or break bonds must likewise publish a new
Topology after explicit connectivity inference rather than mutating the old
one in place.

### Construction, roles, and provenance

`TopologyBuilder` constructs topology without coordinates and supports explicit
definition reuse. It validates connected non-empty definitions, hierarchy
consistency, instance references, identifier capacity, and static metadata.

Instance roles such as `Polymer`, `Branched`, `NonPolymer`, `Solvent`, `Ion`,
`Ligand`, and `Cofactor` are conservative semantic annotations, not hidden
molecular boundaries.

Source provenance belongs in format documents, interpretation reports, or
focused provenance objects. A single source entity may map to several connected
molecule instances.

## Geometry, units, and coordinate state

General 3D primitives (`Point3`, `Vector3`, matrices, `PeriodicCell`,
`RigidTransform`) live in dependency-light geometry code. A periodic cell is
dynamic coordinate state, not topology, and may vary between frames.

Physical values cross public boundaries as `Quantity<T>` with explicit `Unit`
values. Numerical kernels may use canonical raw values only after checked
conversion.

### `Positions`

`Positions` is one complete finite Cartesian array in authoritative
`TopologyAtomIndex` order and retains the exact shared `Arc<Topology>`.
Construction/replacement validates topology compatibility, exact atom count,
units, and finite coordinates.

Matching atom count does not make a `Positions` value compatible with a
different topology.

### `AtomData` and `BondData`

`AtomData` and `BondData` are topology-bound model-level dense data. `AtomData`
has canonical occupancy and isotropic B-factor fields; both types may carry
conservative unit-aware scalar custom columns. Custom properties are
annotations, not topology, represented chemistry, or perception state.

Format provenance such as altloc labels, source model IDs, and raw atom-row IDs
is not generalized into these property containers.

## Structural realizations

### `Model`

`Model` is one concrete realization of one topology:

```text
Model
  Arc<Topology>
  Positions
  Option<PeriodicCell>
  AtomData
  BondData
```

It is also the ordinary application-facing navigation object. Common atom,
bond, molecule-instance, and qualified SMCRA accessors are thin forwards to the
shared topology. Replacing topology-bound state requires the same shared
Topology allocation.

Cloning a model shares topology and copies mutable structural state.

### `Ensemble`

`Ensemble` is one shared topology plus a finite stable-order collection of
non-temporal members. Members carry complete positions and may carry cell,
atom/bond data, weight, and annotations.

Ensemble order has no inherent temporal meaning. Missing/inconsistent atoms are
handled by explicit reconciliation or structured error rather than sparse dense
coordinate arrays.

### `ModelView`

Coordinate-dependent kernels consume `ModelView` (or a narrower equivalent): a
borrowed topology, positions, optional cell, atom data, and bond data. Models,
ensemble members, trajectory frames, and reusable trajectory buffers can expose
the same zero-copy view.

Read-only analyses should operate on this view where practical. A view does not
imply that every algorithm supports every dynamic field; unsupported periodic
or other state is reported explicitly.

### Topology-bound selections

Compiled selections retain the exact shared `Arc<Topology>` and store dense
indices. Selection syntax, semantic resolution, and compiled topology-bound
selection are separate layers. Qualified chain/residue IDs ensure that selecting
one occurrence of a reused definition never implicitly selects another.

## Trajectories

Ordered trajectory state belongs to the one-way `kekule-traj` companion:

```text
kekule <- kekule-traj <- applications
```

A `Trajectory` is an ordered sequence of complete frames sharing one immutable
Topology. It is not an ensemble or a molecule conformer collection. Frames may
vary positions, cell, atom/bond data, velocities, forces, time, step, and frame
metadata as supported.

Trajectory I/O is streaming-first. Reusable `FrameBuffer` storage avoids
materializing entire large files; sequential and seekable reader capabilities
remain distinct. Topology-free formats require an explicit atom-order contract
and never infer identity from atom count alone.

Ordinary trajectories have fixed topology. Reactive data is represented as
fixed-topology segments linked by explicit topology-changing events and
mappings rather than by weakening the fixed-topology invariant.

Codec-specific dialects, byte layouts, indexing rules, and safety limits belong
in `kekule-traj` documentation and tests, not in this core architecture file.

## Format I/O and interpretation

There is no universal text `Document` trait. Each format owns the
loss-preserving representation appropriate to its grammar.

- **SMILES** preserves source syntax and component separators. Component-aware
  interpretation yields one connected `SmallMolecule` per disconnected source
  component.
- **Molfile** preserves V2000/V3000 syntax. The current single-molecule
  interpretation contract requires one connected CTAB graph.
- **SDF** preserves ordered records and raw data fields; record metadata is not
  indiscriminately injected into molecular properties.
- **mmCIF** preserves blocks, items, loops, missing-value markers, unknown
  categories, and source locations.

mmCIF intentionally uses several semantic layers:

```text
MmcifDocument
  -> atom-site/source interpretation
     - explicit model selection
     - explicit altloc selection
     - explicit/special covalent _struct_conn handling
  -> authoritative connectivity completion
     - supplied _chem_comp_bond
     - evidence-backed standard polymer links
     - evidence-backed branched links
  -> residual connected-component partition
  -> Topology + Model / Ensemble + provenance/report
```

Only evidence-backed covalent links establish connectivity. Missing atoms are
not synthesized, sequence gaps are not bridged, and coordinate-distance bond
candidates remain diagnostic only. Residual observed fragments become separate
connected molecule instances while shared source entity/chain identity remains
in provenance.

Multi-model mmCIF interpretation may form an `Ensemble` only after proving
consistent final topology, atom identity, and dense order across selected
models. Alternate-location choice is coordinate-row provenance, not a new atom
identity.

Text writers operate on canonical objects and reject semantics they cannot
represent faithfully.

## Perception, analysis, and transformations

Perception algorithms primarily operate on local connected `Molecule`
definitions. Topology construction preserves installed coordinate-independent
perception and does not normalize or perceive implicitly.

Coordinate-derived analyses return snapshot results; changing positions does
not mutate prior analysis output. Rigid alignment, DSSP, RMSD, contacts, and
similar analyses are read-only unless a separate transformation API explicitly
publishes modified coordinates.

Hydrogen addition/removal and other chemistry-changing operations are explicit
molecule/topology transformations. Topology variants return new topology plus
lineage; newly materialized atoms receive no invented geometry.

## Prepared systems and potentials

Force-field parameters, virtual sites, Drude particles, constraints, mechanical
particles, neighbor-list caches, electronic state, and backend execution state
do not belong in `Topology`, `Model`, `Ensemble`, or `Trajectory`.

A prepared system:

- is created explicitly from one Topology;
- retains that exact shared `Arc<Topology>`;
- maps backend particles to topology identities/indices;
- may contain backend-only particles such as virtual sites;
- does not mutate canonical topology or structural state;
- may evaluate compatible `ModelView` values.

The dependency-light potential contract remains in `kekule`; concrete
implementations may live in one-way companion or adapter crates, for example:

```text
kekule <- kekule-potentials <- applications
```

Preparation may assign implementation-specific parameters, atom types, or
charges. Evaluation must not implicitly normalize, perceive, change topology,
or silently update chemical state.

## Module responsibilities

The public architecture is organized by responsibility rather than by one
universal object:

```text
core         connected Molecule graph, represented chemistry, local IDs,
             perception-state domain model, checked construction/editing
small        SmallMolecule workflows
bio          MacroMolecule and coordinate-independent SMCRA hierarchy
topology     immutable system topology, definitions/instances, qualified IDs,
             dense indices, topology mappings and transforms
structure    Positions, AtomData, BondData, Model, Ensemble, ModelView
geometry     3D primitives, cells, rigid transforms
units        Unit and Quantity
descriptors  read-only molecular descriptors
query        syntax-neutral chemical query representation
substructure query matching algorithms
dssp         read-only secondary-structure assignment
alignment    rigid structural alignment
modeling     potential/prepared-system and numerical workflow contracts
format facades / io
             format documents, interpretation, reports/provenance, writers
```

Implementation modules may be more granular, but dependency direction must
respect these ownership boundaries. Heavy trajectory codecs and potential
implementations belong in one-way companion crates rather than forcing their
dependencies into foundational `kekule`.

## Public API policy

During the `0.x` line, deliberate breaking public API changes require a minor
version increment. The architecture is the contract; pre-1.0 compatibility
shims are not architectural requirements and should not survive without a real
external compatibility need.

Invariant-bearing state is private behind checked constructors and accessors.
Public fields are reserved for deliberate value/options/report payloads.
Extensible public error enums should be `#[non_exhaustive]`.

Parsing, interpretation, normalization, perception, specialized stereo/CIP
derivation, topology construction, coordinate construction, preparation,
analysis, transformation, and writing remain visibly distinct in naming and
documentation even when convenience APIs compose them.

## Architectural invariants

1. Every non-empty canonical `Molecule` is exactly one connected represented
   atom/bond graph.
2. Disconnected chemistry and unresolved experimental fragments are represented
   as multiple molecule instances, with source relationships retained as
   provenance rather than fake bonds.
3. Primary represented chemistry and derived perception are distinct.
4. Interpretation does not infer general chemistry; normalization changes
   representation without changing represented chemical meaning; perception
   derives meaning without rewriting primary representation.
5. Normalization is deterministic, idempotent, model-independent, and separate
   from chemistry-changing standardization.
6. `Topology` is one immutable coordinate-free system of connected molecule
   definitions and instances.
7. Definition identity and instance identity are separate; explicit definition
   reuse is supported.
8. Semantic identities and dense numerical indices are separate, and one
   Topology has one authoritative immutable dense order.
9. Topology-bound state is compatible by exact shared `Arc<Topology>`
   allocation; layout equality alone does not establish compatibility.
10. `Model` is Topology plus one complete structural realization; `Ensemble` is
    finite non-temporal state; `Trajectory` is ordered temporal/sequential state.
11. Periodic cells, atom/bond data, velocities, forces, and time are dynamic
    state, not topology.
12. Topology-changing operations return a new topology and explicit mapping and
    never publish a disconnected molecule definition.
13. Reactive workflows use explicit topology changes rather than mutable
    topology hidden inside coordinate state.
14. Exact persistence reconstruction preserves represented chemistry,
    hierarchy/stereo slot layout, and installed perception without silently
    normalizing or re-perceiving.
15. Only evidence-backed source chemistry establishes topology connectivity;
    geometric connectivity guesses remain explicit diagnostics or explicit
    future perception, never hidden interpretation behavior.
16. Prepared systems and compiled selections retain one exact shared topology
    allocation and do not mutate canonical state.
17. Backend-specific mechanical/execution state never becomes canonical
    chemistry or structure merely because a backend needs it.
