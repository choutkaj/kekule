# Architecture

## Purpose

This document is the normative architecture contract for `kekule`. It defines
object ownership, semantic boundaries, and invariants. Detailed API behavior
belongs in Rustdoc and tests; implementation history does not belong here.

`kekule` is a pure-Rust foundation for cheminformatics, structural
bioinformatics, molecular structure handling, and molecular modelling. It does
not make a file-format record, a simulation-engine particle list, or one
flattened graph the universal data model.

The foundational chemical object is `Molecule`: one connected canonical
represented chemical entity. The foundational system object is `Topology`: one
immutable, coordinate-free molecular system composed of one or more connected
molecule instances.

## Canonical object model

```text
source text / bytes
    -> format-specific Document or streaming decoder
    -> interpretation + canonicalization + provenance/report
    -> canonical represented connected molecular objects
       Molecule
         |- SmallMolecule
         `- MacroMolecule
    -> optional perception / stereo derivation
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
- **Interpretation** translates source-asserted semantics into Kekule's
  canonical represented chemistry. It performs deterministic, model-independent,
  meaning-preserving canonicalization needed to publish a valid `Molecule`,
  including localization of source aromatic bond representations and conversion
  of source stereo marks into canonical stereo elements. It does not run general
  chemical perception or invent bonds merely to force connectedness.
- **Perception** derives chemical meaning from canonical represented chemistry
  and installs a derived view without rewriting primary represented chemistry.
- **CIP assignment** and coordinate-derived stereo are specialized explicit
  derivations, not parsing or interpretation.
- **Standardization**, if added, may choose among chemically distinct related
  states such as tautomers or protonation states and is therefore separate from
  interpretation and perception.
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

A published `Molecule` is already expressed in Kekule's canonical represented
chemistry. In particular, source-only representation artifacts do not survive
this boundary: aromatic source bonds are localized to ordinary bond orders,
source stereo marks are resolved into canonical stereo elements, and any other
supported meaning-preserving representation rewrite required by Kekule has
already been applied.

A `Molecule` owns:

- stable local `AtomId` and `BondId` values;
- represented atoms, localized bonds, adjacency, and canonical stereo
  elements/groups;
- optional local source conformers;
- scalar annotations;
- installed derived `PerceptionState`.

Temporary disconnection or source-specific representation belongs only to
checked construction/interpreter staging. `MoleculeBuilder` and `MoleculeEditor`
may work on incomplete candidates, but publication is transactional, enforces
connectedness, and may only publish states expressible in canonical core types.
Stable local IDs are not silently reused after deletion; fixed-width identifier
capacity is checked before mutation.

Disconnected chemistry is represented above this boundary. `[Na+].[Cl-]`, a
protein-ligand complex, solvents, ion pairs, and other noncovalent assemblies
are multiple connected molecules in one topology. An asserted Kekule bond
contributes to connectedness; spatial association does not.

### Represented versus derived chemistry

Primary represented chemistry is stored directly and includes, as applicable:

```text
Atom: element, isotope, formal charge, radical, represented hydrogen
      declaration, atom map
Bond: endpoints, localized represented BondOrder
Stereo: local StereoElement state and StereoGroup relations
```

Source-only syntax and provenance are not canonical molecular chemistry.
Aromatic source bond kinds, wedge/dash or directional source marks, parser
artifacts, source-format tags, and normalization diagnostics belong to the
format-specific document, interpreter staging, mappings, or reports rather than
the published `Molecule`.

Derived chemistry is a view over the canonical representation:

```text
valence / implicit hydrogens
cycle membership and ring basis
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
replacement; it does not canonicalize or re-perceive the molecule.

### Canonical bond representation

Canonical `BondOrder` describes localized represented bonding only. Aromaticity
is not a canonical bond order.

Aromatic bond syntax may exist in source documents or private interpreter
staging, but interpretation localizes it before publishing `Molecule`. A
canonical aromatic ring is therefore represented by ordinary localized bond
orders plus, after perception, aromatic atom/bond membership in
`PerceptionState`.

The same distinction applies to query semantics: an "aromatic bond" in a query
is a query predicate, not a canonical `BondOrder` value.

### Represented hydrogens versus perceived hydrogens

Hydrogen declarations that are explicitly part of the represented chemical
statement belong to `Atom`. Hydrogens inferred from valence rules belong to
`PerceptionState`.

The exact core API may use a dedicated represented-hydrogen type rather than
multiple partially overlapping fields, but it must preserve this semantic
boundary:

```text
represented hydrogen declaration -> Molecule / Atom
inferred implicit hydrogen count  -> PerceptionState / valence state
```

Materializing or collapsing graph hydrogen atoms is an explicit topology
transformation and invalidates dependent perception.

### Stereo representation versus stereo perception

Canonical stereochemical representation belongs to `Molecule`. Local stereo
focus, carriers, orientation, specifiedness where semantically required, and
stereo-group relationships are represented chemistry.

CIP descriptors are derived chemistry and belong to `PerceptionState`. They are
not stored as fundamental atom/bond/stereo-element identity because descriptors
can change when graph-wide priority relationships change without changing the
local represented stereo orientation.

Source stereo marks such as molfile wedges or SMILES directional bonds are
interpretation inputs, not canonical `Molecule` state. Interpretation converts
them into canonical `StereoElement` state or reports a failure/diagnostic before
publication. Source-format provenance similarly belongs to mappings/reports or
annotations, not to the canonical stereo element unless it is itself chemically
meaningful.

### `PerceptionState`

`PerceptionState` is the installed discrete derived interpretation of one exact
canonical `Molecule` representation. It should remain sectional and narrowly
focused on fundamental derived chemistry:

```text
PerceptionState
  |- ValenceState
  |    |- model/provenance
  |    `- implicit-hydrogen assignments
  |- RingState
  |    |- graph cycle membership
  |    `- optional deterministic ring basis + basis provenance/model
  |- AromaticityState
  |    |- model
  |    |- aromatic atoms
  |    `- aromatic bonds
  `- StereoPerceptionState
       `- CIP descriptor assignments
```

Not every calculated molecular property belongs in `PerceptionState`.
Descriptors, force-field atom types, partial charges, rotatable-bond labels,
functional-group matches, H-bond feature labels, and similar task-specific
results should remain separate derived objects unless they become foundational
requirements shared by core chemistry algorithms.

Graph cycle membership and a chosen ring basis are distinct concepts. Membership
is a graph property; a stored ring basis is an algorithmic choice and should
carry sufficient provenance to identify how it was constructed when that choice
matters downstream.

Perception invalidation should remain simple and dependency-safe. Broad
chemistry edits may clear the complete state. Narrow edits may clear only
provably downstream sections, but Kekule should prefer robust invalidation over
a complex generic dirty-state engine.

Public chemistry accessors should hide physical storage details. Callers may
use molecule-level convenience methods or future molecule-aware atom/bond views
for operations such as aromaticity, ring membership, implicit hydrogens, and CIP
without requiring those values to be duplicated into `Atom` or `Bond`.

### `SmallMolecule`

`SmallMolecule` is the ordinary cheminformatics wrapper around one connected
`Molecule`. It may provide ergonomic parse/interpret/perceive workflows, but it
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
       |- decode source semantics
       |- perform deterministic canonical representation rewrites
       |- localize source aromatic bonds
       `- resolve source stereo into canonical StereoElements
  -> canonical Molecule
  -> perceive
  -> derived PerceptionState
  -> optional coordinate-stereo materialization / CIP assignment
```

Parsing answers **what syntax is present**. Interpretation answers **what
canonical chemical representation Kekule publishes for what the source
asserted**. Perception answers **what chemical interpretation Kekule derives
from that canonical graph**.

There is no separate public "unnormalized Molecule" lifecycle state and no
required `MoleculeDraft` domain object. Format-specific `Document` values and
private interpreter/builder staging already provide the necessary source-side
workspace. Canonicalization may remain factored into reusable internal helper
algorithms, but it is a responsibility of interpretation or checked molecule
publication rather than a distinct public chemistry stage.

A published canonical molecule has:

- ordinary localized represented bond orders; aromaticity is not a core bond
  order;
- no unresolved source `StereoBondMark`-style state;
- canonical represented local `StereoElement` state;
- canonical represented charge/hydrogen declarations;
- no requirement that derived `PerceptionState` be present.

Interpretation does not choose an aromaticity model, tautomer, protonation
state, salt/fragment policy, force-field typing scheme, or other chemically
non-equivalent standardization.

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

Reconstruction builds canonical represented graph/hierarchy/stereo state first
and installs derived perception last. Stable stereo-group slot layout, including
tombstones, is reconstructible without dummy chemistry. Loading must not
re-perceive, renumber, or silently coerce malformed historical state. If a
persistence format stores a legacy noncanonical source representation, its
loader must interpret/canonicalize that representation before publishing a
`Molecule` rather than weakening the runtime invariant.

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

`Topology` does not own coordinates, velocities, forces, cell vectors, energies,
trajectory state, or prepared backend objects.

A `Topology` definition may contain one `Molecule` once and instantiate it many
times. Definition-local IDs remain stable within the definition; system-level
access qualifies them by molecule instance.

### `Model`, `Ensemble`, and `Trajectory`

`Model` is one exact coordinate state over one shared `Topology` allocation.
`Ensemble` is a finite set of non-temporal models sharing one topology.
`Trajectory` is an ordered frame source sharing one topology and may be memory-
or file-backed.

These types do not duplicate chemistry perception per coordinate frame merely
because positions change. Coordinate-dependent analyses or future
position-dependent perception are separate derived objects unless and until an
explicit architecture extends this contract.

## Mutation and transformations

Mutating canonical represented chemistry is allowed only through invariant-
preserving APIs. Any operation that can transiently violate connectedness or
other publication invariants uses checked staging and transactional publication.

Topology-changing operations such as deleting atoms, changing connectivity,
adding/removing materialized hydrogens, reaction transforms, or future
connectivity perception return or publish a new canonical molecular/topological
state and explicit mappings where identity transfer matters.

Coordinate-only operations do not mutate canonical chemistry or its discrete
perception state.

## Design rule

When deciding where new state belongs, apply this order:

1. Is it syntax/source representation? Keep it in the format document or
   interpreter provenance/staging.
2. Is it part of Kekule's canonical represented chemical statement? Store it in
   `Molecule`.
3. Is it fundamental chemistry derived algorithmically from that exact
   representation? Store it in an appropriate `PerceptionState` section.
4. Is it task-specific analysis, typing, scoring, or parameterization? Keep it
   in a separate derived object.
5. Is it coordinate state? Keep it outside `Topology` and canonical molecular
   chemistry.

This boundary is preferred over duplicating derived flags into atoms/bonds or
creating multiple public half-canonical molecule lifecycle states.
