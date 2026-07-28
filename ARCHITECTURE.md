# Architecture

## Status and purpose

This document is the normative architecture contract for `molecular`. The
topology-centered 0.2 design and the initial complete-instance topology-
transformation milestone are implemented. Current work must follow the object
boundaries and invariants defined here.

`molecular` is a pure-Rust foundation for cheminformatics, structural
bioinformatics, molecular structure handling, and molecular modelling. It
serves small molecules and biological macromolecules without making a
file-format record, a simulation-engine particle list, or one flattened graph
the universal data model.

The foundational chemical concept is `Molecule`: one asserted chemical entity.
The foundational system-level concept is `Topology`: one immutable,
coordinate-free molecular system composed of one or more molecule instances.

## Canonical object model

```text
format text or binary data
    -> format-specific Document or streaming decoder
    -> explicit interpretation + report
    -> canonical chemical and structural objects
       Molecule
         |- SmallMolecule
         `- MacroMolecule
       Topology
         |- Model       = Topology + one Configuration
         |- Ensemble    = Topology + finite non-temporal members
         `- Trajectory  = Topology + ordered Frames
    -> explicit perception, validation, transformation, or analysis
    -> downstream prepared System, Potential, or backend object
```

The central ownership relationships are:

```text
Molecule
  local asserted chemical graph
  local AtomId and BondId
  optional source conformers

Topology
  immutable coordinate-free system
  molecule definitions
  molecule instances
  instance-qualified atom and bond identities
  authoritative dense atom and bond orderings

Model
  Topology
  one mutable Configuration

Ensemble
  Topology
  finite collection of EnsembleMember values

Trajectory
  Topology
  ordered collection or stream of TrajectoryFrame values
```

`Topology` is shared by `Model`, `Ensemble`, `Trajectory`, prepared force
fields, compiled selections, and analyses. Coordinate updates never modify it.

## Architectural boundaries

These boundaries are requirements rather than naming conventions:

- Parsing recognizes and preserves format syntax. It does not sanitize
  chemistry, run perception, parameterize a force field, or silently choose
  ambiguous structural records.
- Interpretation applies format semantics and documented policies. It creates
  canonical objects, source-to-canonical mappings, imported annotations, and a
  structured report. It still does not sanitize implicitly.
- Perception derives chemical state from asserted molecular topology.
- Sanitization is an explicit transactional workflow over canonical chemical
  objects.
- Topology construction assembles already asserted molecular entities into one
  coordinate-free system. It does not invent coordinates or force-field state.
- Coordinate construction validates complete state against one exact topology.
- Analysis is read-only unless the operation is explicitly named as a
  transformation.
- Force-field preparation and backend preparation create downstream objects
  bound to a topology. They never mutate canonical chemistry or coordinate
  state.
- Writers reject unsupported semantics rather than silently dropping or
  coercing them.

Physical values cross boundaries as `Quantity<T>` with composable runtime
`Unit` values. Compatible conversion is explicit. Numeric units, source
provenance, and scientific meaning remain separate concepts.

## Geometry and physical state

General three-dimensional geometry belongs in a dependency-light `geometry`
module rather than in a force-field or conformer implementation module. The
common vocabulary includes:

```rust
Point3
Vector3
Matrix3
PeriodicCell
RigidTransform
```

The type system should preserve useful distinctions such as:

```text
Point3 - Point3 -> Vector3
Point3 + Vector3 -> Point3
Vector3 + Vector3 -> Vector3
```

A periodic cell is coordinate state, not topology. It may vary between
trajectory frames, for example under constant-pressure simulation. A
`PeriodicCell` must support orthorhombic and triclinic cells and must validate
finite, non-degenerate vectors and explicit periodic axes.

Coordinates, velocities, forces, gradients, time, and cell dimensions use
explicit units at public boundaries. Numerical kernels may use canonical raw
values after one checked conversion.

## Molecules and domain wrappers

### `Molecule`

`Molecule` is the raw chemical graph kernel and the asserted entity boundary.
It owns:

- stable typed `AtomId` and `BondId` values;
- atoms, bonds, adjacency, graph-adjacent stereo elements, stereo groups, and
  source stereo marks;
- optional local source conformers, each with one explicit length unit;
- arbitrary scalar annotations;
- one internally consistent `PerceptionState`.

Deletion leaves tombstones and stable local identifiers are never reused.
`Molecule` may be disconnected, unsanitized, incomplete, or chemically invalid.
Atom, bond, conformer, stereo-element, and stereo-group insertion checks the
fixed-width identifier slot before mutation and returns a focused structured
capacity error when exhausted. Iteration never reconstructs these identifiers
through unchecked narrowing from platform-sized collection indices.

Graph connectedness is not the definition of a molecular entity. A salt such as
`[Na+].[Cl-]`, a coordination compound such as ferrocene, and an ordinary
connected organic compound may each be represented as one `Molecule`.
Connected-component algorithms report the actual graph components without
overriding the asserted entity boundary.

Topology facts stored directly on atoms, bonds, and stereo elements include
element, isotope, formal charge, radical state, explicit-hydrogen declarations,
atom maps, bond endpoints, `BondOrder`, local stereo, and source stereo marks.

Implicit hydrogens, ring membership, ring sets, aromatic membership,
aromaticity provenance, and derived CIP descriptors are computed state. They
live in private optional sections of `PerceptionState` and are exposed through
read-only queries. One installed perception profile exists at a time.
Alternative-model calculations remain standalone results until explicitly
installed.

Chemistry- or connectivity-relevant mutation invalidates affected perception
state immediately. Coordinate and generic annotation edits are
perception-neutral. Failed transactional operations leave their input
unchanged.

### `SmallMolecule`

`SmallMolecule` is the ordinary cheminformatics wrapper around one `Molecule`.
It provides ergonomic small-molecule workflows while retaining read access and
controlled mutation of the graph.

`SmallMolecule::from_smiles` is an intentional parse-then-interpret convenience
and does not sanitize. Any convenience that also sanitizes must state that
operation in its name.

The physical wrapper type remains dependency-light. I/O, perception,
sanitization, modelling, and workflow facades may depend on it; lower layers do
not depend back on workflow conveniences.

### `MacroMolecule` and `SmcraHierarchy`

`MacroMolecule` is one `Molecule` plus one validated `SmcraHierarchy`.
`SmcraHierarchy` stores coordinate-independent structural identity:

- chains;
- residues;
- atom sites;
- label and author identifiers;
- insertion codes and component names;
- mappings from structural atom sites to local `AtomId` values;
- other static hierarchy annotations.

Coordinate-model membership is not topology. Source model identifiers,
alternate-location choices, occupancy, B-factors, raw Cartesian text, and other
observation-specific values belong to interpretation provenance or
configuration-associated observation data. The hierarchy must not require a
coordinate-model node as the parent of a chain in the final architecture.

Every public `MacroMolecule` is valid: every live graph atom has exactly one
atom site, every atom site references a live graph atom, hierarchy parentage is
consistent, and static identifiers satisfy their documented invariants.
Construction uses checked builders or checked assembly from parts.
Coordinated graph-and-hierarchy mutation is transactional.
Chain, residue, and atom-site capacity failures are structured and occur before
parent lists or lookup maps are changed.

Macromolecule validation is separate from small-molecule sanitization.
Chemically general algorithms operate on `Molecule` where practical.

### Local conformers versus system coordinate state

Standalone `Molecule`, `SmallMolecule`, and `MacroMolecule` values may retain
local conformers because conformer workflows are useful in cheminformatics and
structure I/O.

A molecule definition stored in a `Topology` is coordinate-free. Topology
construction copies or moves only coordinate-independent molecular state.
Source conformers are neither duplicated into topology nor silently discarded
from the source object.

Conversions are explicit:

```text
SmallMolecule + one selected conformer
    -> Topology + Model

SmallMolecule + selected conformers
    -> Topology + Ensemble

several molecules + one coordinate set per instance
    -> Topology + Model
```

## `Topology`

### Definition

`Topology` is the immutable, coordinate-free definition of one fixed molecular
system. It is a first-class public object independent of `Model`.

A `Topology` owns:

- molecule definitions;
- molecule instances referencing those definitions;
- all asserted atoms and covalent bonds available through instance-qualified
  identities;
- system-wide connectivity and adjacency access;
- static chemical and structural metadata;
- molecule roles and instance annotations;
- one authoritative dense atom ordering;
- one authoritative dense bond ordering;
- mappings between semantic identifiers and dense indices;
- exact topology identity used by topology-bound downstream objects.

A `Topology` does not own:

- positions;
- periodic cells;
- velocities;
- forces or gradients;
- time or simulation step;
- coordinate-derived contacts or secondary structure;
- occupancy, B-factors, or raw coordinate text tied to one observation;
- force-field parameters or atom types;
- constraints, virtual sites, Drude particles, or backend particles;
- execution-engine state.

Public `Topology` is a cheap-clone immutable handle, conceptually:

```rust
#[derive(Clone)]
pub struct Topology {
    inner: Arc<TopologyData>,
}
```

Users should work with `Topology`, not with `Arc<Topology>` directly.

### Molecule definitions and molecule instances

Topology distinguishes what a molecule is from one occurrence of that molecule
in a system:

```rust
pub struct MoleculeDefinitionId(u32);

pub enum MoleculeDefinitionPayload {
    Small(SmallMolecule),
    Macro(MacroMolecule),
}

pub struct MoleculeDefinition {
    id: MoleculeDefinitionId,
    payload: MoleculeDefinitionPayload,
}

pub struct MoleculeInstance {
    id: MoleculeInstanceId,
    definition: MoleculeDefinitionId,
    metadata: MoleculeInstanceMetadata,
}
```

A definition stores one conformer-free small- or macromolecular entity. An
instance identifies one occurrence of that definition in the topology.

Definition reuse is explicit. For example, one water definition may be
referenced by many water instances. Each instance still has a unique
`MoleculeInstanceId`, unique instance-qualified atom identities, roles,
annotations, coordinates in associated configurations, and source provenance.

The library must not automatically merge definitions merely because their
graphs compare equal. Equal-looking molecules may differ in isotopes,
perception state, hierarchy data, static annotations, or scientific meaning.

An implementation may initially create one definition per inserted instance,
but the public model and builder must support explicit reuse without a later
semantic redesign.

### Molecular boundaries and connectivity

Every topology atom belongs to exactly one molecule instance.

In the fixed-topology architecture, asserted covalent bonds belong within one
molecule instance. If two structural entities are covalently linked, they are
represented as one asserted molecular entity or are explicitly merged during
interpretation before topology construction.

Molecule instances are not synonyms for connected components. One molecule
definition may itself contain multiple graph components. Downstream operations
must distinguish, as applicable:

- asserted molecule instances;
- connected components;
- charge-assignment groups;
- nonbonded exclusion groups;
- periodic imaging groups;
- user-defined atom groups.

Noncovalent contacts, hydrogen bonds, restraints, spatial neighbors, and
coordinate-derived interactions are not topology bonds.

Topology provides system-wide graph queries using qualified identities while
preserving local molecule graphs. It must support iteration and lookup of:

```text
molecule definitions
molecule instances
instance-qualified atoms
instance-qualified bonds
neighbors and incident bonds
connected components within each instance
instance and definition membership
hierarchy views for macromolecular instances
```

### Semantic identifiers and dense indices

Semantic identity and numerical storage order are separate.

Local chemical identifiers remain:

```rust
AtomId
BondId
```

Topology-level semantic identifiers are:

```rust
MoleculeDefinitionId
MoleculeInstanceId
InstanceAtomId {
    molecule: MoleculeInstanceId,
    atom: AtomId,
}
InstanceBondId {
    molecule: MoleculeInstanceId,
    bond: BondId,
}
```

Dense numerical identifiers are:

```rust
TopologyAtomIndex(u32)
TopologyBondIndex(u32)
```

`InstanceAtomId` answers which local atom of which molecule instance is being
addressed. `TopologyAtomIndex` answers where that atom is stored in complete
position, velocity, force, gradient, and per-atom result arrays.

The dense ordering is authoritative and immutable for the lifetime of one
topology. It need not be identical between independently constructed
topologies representing the same chemistry. Every topology-bound dense array
uses the ordering published by that exact topology.

Local `AtomId` and `BondId` values, including tombstone positions, survive
definition insertion. Qualification adds instance ownership without remapping
local identifiers. Dense indices include only live atoms and bonds.

All conversions from collection lengths to fixed-width public identifiers are
checked. Capacity overflow produces structured errors rather than truncation or
wrapping.

### Identity and layout equality

Exact topology identity is the compatibility criterion for positions, prepared
systems, compiled selections, and frame buffers.

Clones of one `Topology` retain exact identity. Independently constructed
topologies have different identity even when their complete static layouts are
equal.

The API distinguishes:

```rust
topology_a.same_identity(&topology_b)
topology_a.same_layout(&topology_b)
```

`same_layout` compares chemical and hierarchy content, definition and instance
partitioning, instance metadata, semantic identifiers, authoritative dense
order, and the corresponding index maps. It excludes only exact identity.
Therefore it is not order-independent structural equivalence and does not
silently imply compatibility for topology-bound state.

General structural equivalence and validated isomorphism mapping across
different definition, instance, atom, or bond orderings remain planned future
capabilities. Such an operation must account for all relevant chemical,
stereochemical, perception, role, and hierarchy state and must report ambiguous
mappings rather than selecting one silently.

Transferring positions, parameters, or selections between independently
constructed topologies requires an explicit validated mapping.

### Immutability and topology-changing operations

A built topology is immutable. This permits safe sharing by models, ensembles,
trajectories, prepared systems, selections, and caches.

Topology-changing operations return a new topology and explicit lineage
mappings:

```rust
pub struct TopologyEditResult {
    pub topology: Topology,
    pub mapping: TopologyMapping,
}
```

Mappings describe retained, removed, and added molecule definitions, molecule
instances, atoms, bonds, and dense indices as applicable.

Examples include:

```text
add or remove hydrogens
merge or split molecule instances
delete atoms or bonds
change asserted bond order
construct a solvated system
apply a chemical reaction
```

A topology edit never mutates existing models, ensembles, trajectories, or
prepared systems. Remapping coordinate state is a separate explicit operation.
Added atoms may remain without coordinates until a geometry-building operation
supplies them; topology transforms do not invent geometry silently.

Complete molecule instances can be retained or removed through the focused
`topology::transform` namespace. These immutable deletion-only edits preserve
filtered source definition and instance order, explicit definition reuse, and
local atom and bond identifiers while returning complete checked lineage.
Positions, configurations, observations, models, compiled atom selections,
finite ensembles, owned frames, in-memory trajectories, and borrowed frame
state in reusable target buffers provide explicit remapping operations. Every
operation checks exact source and target identity; complete dense arrays reject
unmapped target atoms, and selection loss requires an explicit policy.

### Construction

`TopologyBuilder` constructs topology without requiring coordinates.

It supports:

- adding checked molecule definitions;
- adding multiple instances of one definition;
- attaching instance roles and static annotations;
- reserving capacity for large systems;
- obtaining provisional identifiers where needed;
- transactionally rejecting invalid additions;
- building authoritative dense atom and bond mappings.

The builder must be linear in the amount of appended data. Transactionality is
implemented by validating and staging only the new addition before infallible
append. It must not clone the entire accumulated builder for every inserted
molecule.

Convenience builders may jointly assemble topology and coordinate state, but
topology construction remains independently available.

A built topology rejects:

- no molecule instances;
- empty molecule definitions;
- invalid macromolecule graph/hierarchy pairs;
- invalid instance-to-definition references;
- identifier capacity overflow;
- inconsistent local atom or bond references;
- static metadata violating documented invariants.

Topology construction validates macromolecular graph/hierarchy state without
scanning source conformers. A convenience model builder separately validates
only its explicitly selected conformer while staging positions. Standalone
full macromolecule validation may inspect every retained source conformer.

### Roles, properties, and provenance

A molecule instance may have several roles, including:

```text
Polymer
Branched
NonPolymer
Solvent
Ion
Ligand
Cofactor
```

Roles are conservative semantic annotations and are not inferred from graph
connectedness alone.

Generic property maps are for scalar annotations. Core structural or dynamic
state must use typed fields and validated arrays rather than ad hoc property
keys.

Source provenance remains in format documents, interpretation reports, or
dedicated provenance objects. Exact source record identifiers are not injected
indiscriminately into generic atom and molecule property maps.

## Coordinate state

### `Positions`

`Positions` is one complete finite Cartesian coordinate array in
`TopologyAtomIndex` order.

Construction validates:

- exact topology identity or an explicit assertion against one topology;
- one position for every topology atom;
- compatible length units;
- finite coordinates;
- checked array length.

Internally, positions use the modelling kernel's declared canonical length
unit. Public access returns quantities with explicit units.

A `Positions` value must not be reused with a different topology merely because
the atom counts match.

### `Configuration`

`Configuration` is one complete geometric realization of a topology:

```rust
pub struct Configuration {
    positions: Positions,
    cell: Option<PeriodicCell>,
}
```

The periodic cell is optional and dynamic. Additional coordinate-state fields
are introduced only when they are semantically part of every configuration.

A borrowed `ConfigurationView` exposes topology-compatible arrays without
allocation.

### Observation data

Experimental structure observations may attach coordinate-specific data such
as occupancy, B-factor, alternate-location label, source coordinate-model ID,
or raw source text.

Such data is stored beside a model or ensemble member in a typed
`StructureObservation` or dedicated provenance object. Per-atom observation
arrays use `TopologyAtomIndex` order and validate their lengths. They are not
part of topology identity.

## `Model`

`Model` is one concrete structural realization of one topology:

```text
Model
  Topology
  one Configuration
  optional structure-observation metadata
```

Conceptually, in the common non-periodic case:

```text
Model = Topology + Positions
```

The topology is immutable and cheap to clone. Positions and the periodic cell
may be replaced transactionally while preserving topology identity.

A model rejects incomplete, non-finite, dimensionally incompatible, or
topology-incompatible coordinate state.

Convenience constructors may build a topology and model from one selected local
molecule conformer. Source objects remain unchanged. Explicit operations copy
one model instance's current positions back to a compatible standalone
molecule conformer.

Cloning a model shares topology and copies only mutable coordinate and
observation state.

## `Ensemble`

`Ensemble` represents a finite collection of non-temporal structural
realizations sharing one exact topology:

```rust
pub struct Ensemble {
    topology: Topology,
    members: Vec<EnsembleMember>,
}

pub struct EnsembleMember {
    configuration: Configuration,
    weight: Option<f64>,
    observation: Option<StructureObservation>,
    props: PropMap,
}
```

Typical ensembles include:

- small-molecule conformer ensembles;
- NMR structure ensembles;
- docking poses;
- alternative predicted structures;
- weighted thermodynamic samples;
- crystallographic or structural alternatives after explicit interpretation.

Ensemble order is stable but carries no inherent temporal meaning. A member may
have a finite non-negative weight, score, energy, source label, or other typed
or scalar annotations. Weight normalization is an explicit operation rather
than an insertion side effect.

Every member has complete state for the same exact topology. Formats with
missing or inconsistent atoms require explicit reconciliation or a structured
error; they are not represented by silently sparse dense arrays.

## `Trajectory`

### Semantics

`Trajectory` represents an ordered sequence of frames sharing one exact,
immutable topology:

```text
Trajectory
  Topology
  ordered TrajectoryFrame values
```

Ordering represents time, simulation step, acquisition order, or another
explicit sequential meaning. A trajectory is not merely a molecule conformer
collection and is not interchangeable with an ensemble.

A frame contains:

```rust
pub struct TrajectoryFrame {
    configuration: Configuration,
    velocities: Option<Velocities>,
    forces: Option<Forces>,
    time: Option<Quantity<f64>>,
    step: Option<u64>,
    observation: Option<StructureObservation>,
    props: FrameMetadata,
}
```

Positions are required and complete. Optional per-atom arrays are complete,
finite, unit-aware, and topology-bound. Cell, time, step, velocity, and force
data may vary by frame.

An in-memory `Trajectory` may own a finite collection of frames. Large
trajectory I/O must not require materializing that collection.

### Streaming and random access

Trajectory I/O is streaming-first and uses caller-owned reusable buffers:

```rust
pub trait TrajectoryReader {
    fn topology(&self) -> &Topology;

    fn read_next(
        &mut self,
        destination: &mut FrameBuffer,
    ) -> Result<bool, TrajectoryError>;
}

pub trait SeekableTrajectoryReader: TrajectoryReader {
    fn frame_count(&self) -> Option<u64>;

    fn read_frame(
        &mut self,
        index: u64,
        destination: &mut FrameBuffer,
    ) -> Result<(), TrajectoryError>;
}

pub trait TrajectoryWriter {
    fn topology(&self) -> &Topology;

    fn write_frame(
        &mut self,
        frame: TrajectoryFrameView<'_>,
    ) -> Result<(), TrajectoryError>;
}
```

Sequential and seekable capabilities remain separate because not every
compressed format supports inexpensive random access.

`FrameBuffer` is initialized for one exact topology and reuses allocations
across reads. A frame view allows analyses and potentials to consume decoded
state without constructing or cloning an owned `Model`.

File decoders publish only complete frames through one transactional borrowed-
data operation. That operation validates exact topology identity, complete
array lengths, units, finite values, cell, time, observations, and properties
before destination-visible mutation; it reuses position, velocity, and force
allocations and clears every absent optional field. Clean EOF and any failed
decode leave the caller's buffer unchanged.

Topology-free file readers require an `AtomOrderAssertion` bound to the same
exact topology. The assertion is either constructed from a complete semantic
atom sequence proven equal to authoritative dense order or is an explicit
caller statement that the file uses that order. Neither form is inferred from
atom count, and readers still validate all stronger format metadata.

### Topology-bearing and topology-free formats

Some trajectory formats contain only coordinates and atom count. Their readers
require an externally supplied `Topology` and a documented atom-ordering
contract.

Readers validate all available compatibility information, including atom
count, ordering metadata, units, periodic cell support, and optional array
lengths. A matching atom count alone is not treated as proof of chemical
identity when stronger metadata is available.

Formats that contain topology may construct it explicitly through
`TopologyBuilder` and return interpretation provenance.

### Fixed and changing topology

Ordinary `Trajectory` has one fixed topology. Bond formation, bond breaking,
atom insertion or deletion, and changing molecule-instance boundaries are
outside that contract.

A future reactive trajectory is represented as fixed-topology segments linked
by explicit events and topology mappings:

```text
ReactiveTrajectory
  segment 1: Topology A + frames
  event and TopologyMapping
  segment 2: Topology B + frames
```

The fixed-topology `Trajectory` invariant is not weakened to accommodate
reactions.

## Borrowed structural views

Coordinate-dependent algorithms consume a common borrowed view:

```rust
pub struct ModelView<'a> {
    topology: &'a Topology,
    configuration: ConfigurationView<'a>,
}
```

`Model`, ensemble members, trajectory frames, and reusable frame buffers can
all provide `ModelView` without copying topology or coordinates.

Read-only analyses such as DSSP, RMSD, alignment, contact analysis, and
geometric measurements should accept this view or a narrower equivalent.
Existing convenience functions accepting `&Model` may delegate to the view
kernel.

Rigid alignment lives in the focused `alignment` analysis module. Its first
milestone fits two `ModelView` values sharing one exact topology through one
topology-bound `AtomSelection`. The returned proper `RigidTransform` maps
moving coordinates into reference coordinates, while post-fit weighted RMSD is
reported in the model length unit. Uniform or explicit positive finite
selection-order weights are supported. Periodic configurations are rejected by
default; an explicit stored-coordinate policy ignores cells without imaging or
unwrapping. The analysis is read-only and never materializes replacement
canonical coordinates.

Determined alignment requires at least three selected atoms and rank-two
geometry on both sides. Planar non-collinear selections are valid. Coincident,
collinear, and scale-relatively near-collinear inputs are rejected through
structured errors. Exact topology identity is mandatory; equal layout or atom
count never establishes correspondence.

## Topology-bound selections

Compiled atom and bond selections bind to exact topology identity and store
dense indices:

```rust
pub struct AtomSelection {
    topology: TopologyIdentity,
    indices: Vec<TopologyAtomIndex>,
}
```

Selections may be constructed from:

- explicit instance-qualified identifiers;
- molecule roles;
- molecule definitions or instances;
- chains, residues, and atom names;
- elements and chemical properties;
- connected components;
- chemical substructure matches;
- a future structural selection language.

Selection syntax, selection meaning, and compiled topology-bound results remain
separate layers. A selection is evaluated once against fixed topology and may
then be reused for every model, ensemble member, or trajectory frame sharing
that topology.

## Format documents and interpretation

There is no generic text-format `Document` trait. Each textual format exposes a
loss-preserving type appropriate to its grammar.

- `SmilesDocument` preserves source syntax and interprets to one
  `SmallMolecule`.
- `MolfileDocument` preserves V2000 or V3000 syntax and interprets to one
  `SmallMolecule` plus mappings and a report.
- `SdfDocument` owns ordered records with raw data fields and interprets to
  canonical records. Record metadata is not injected into molecule properties.
- `MmcifDocument` preserves blocks, items, loops, missing-value markers,
  unknown categories, and source locations.

Ordinary mmCIF model interpretation selects exactly one coordinate-model ID.
The default rejects ambiguous multi-model input. Named and first-model
selection are explicit policies. Alternate-location handling is explicit.

Single-model mmCIF interpretation produces:

```text
Topology
Configuration
optional StructureObservation
MmcifInterpretationReport
```

and exposes a convenience `Model`.

Multi-model mmCIF interpretation is a separate ensemble operation. It may
produce one `Ensemble` only after proving a consistent chemical topology,
molecule-instance partition, atom identity mapping, and dense ordering across
members. Coordinate-model identifiers and per-model observation values belong
to ensemble-member metadata. Inconsistent atom presence or topology produces a
structured error unless an explicit reconciliation policy is requested.
Source atom correspondence uses residue sequence and insertion identity,
label/author asymmetry identity, component and atom labels, an explicit
occurrence discriminator when sequence identifiers are absent, and the
selected alternate location where relevant. It does not infer correspondence
from derived molecule insertion order.

Only evidence-backed covalent links establish topology connectivity.
Distance-based candidates, unresolved connections, model selection, altloc
selection, ignored records, and source classification remain visible in
reports and provenance. Interpretation must not assert fabricated bond orders.

Text writers operate on canonical objects and reject unsupported semantics.
Writer-generated one-based numeric serials use widened integer arithmetic;
formatting never increments a `u32` identifier in place.
Trajectory codecs use streaming reader and writer interfaces rather than
forcing a loss-preserving whole-file document abstraction for binary data.

## Perception, transformations, and derived state

Perception algorithms operate primarily on local `Molecule` definitions.
Topology may expose convenient iteration over definitions and instances, but it
does not duplicate perception algorithms.

Topology construction preserves the installed coordinate-independent
perception state of each inserted definition. It does not sanitize or perceive
implicitly.

Hydrogen addition and removal remain explicit topology-changing chemical
transforms. Applied to standalone molecules, they preserve local stable IDs as
documented. Applied to a topology, they return a new topology and lineage
mapping. Newly materialized atoms do not receive invented coordinates.

Coordinate-derived analyses return snapshot results. Updating positions or
advancing a trajectory does not mutate an earlier analysis result. Callers
explicitly rerun analyses for another configuration.

## Downstream prepared systems and potentials

Force-field parameters, virtual sites, Drude particles, constraints,
mechanical particles, electronic state, neighbor-list caches, and
execution-engine objects do not belong in `Topology`, `Model`, `Ensemble`, or
`Trajectory`.

A prepared system:

- is constructed explicitly from one topology;
- binds to exact topology identity;
- provides mappings between backend particles and
  `TopologyAtomIndex`/`InstanceAtomId`;
- may contain particles not in canonical topology, such as virtual sites;
- does not mutate topology or coordinate containers;
- may evaluate supported model views sharing the bound topology.

Potential evaluation consumes `ModelView` or an equivalent borrowed
topology-plus-configuration view. Preparation is performed once and reused
across models, ensemble members, and trajectory frames with the same topology.
Accepting the common view does not imply support for every dynamic
configuration field. Each potential documents capabilities such as
periodic-cell handling and returns a structured error when a compatible view
contains unsupported state.

Policies that group atoms, such as charge equilibration, must state whether
they operate over the whole topology, molecule instances, connected
components, or explicit groups. The architecture never treats those scopes as
interchangeable by accident.

`molecular-dreiding` remains an adapter crate demonstrating this boundary. It
may assign atom types and fixed charges during explicit preparation, but
evaluation does not sanitize, change topology, or update charges implicitly.

## Public module responsibilities

The target public facade is focused:

```text
core
    Molecule graph kernel, local IDs, atom and bond chemistry

geometry
    Point3, Vector3, matrices, periodic cells, rigid transforms

units
    Unit, Quantity, physical constants, modelling units

small
    SmallMolecule and small-molecule workflow conveniences

bio
    MacroMolecule, coordinate-independent SmcraHierarchy, validation

topology
    Topology, TopologyBuilder, molecule definitions and instances,
    instance-qualified IDs, dense topology indices, topology mappings

structure
    Positions, Configuration, Model, Ensemble, borrowed structural views,
    structure-observation state

trajectory
    TrajectoryFrame, FrameBuffer, in-memory Trajectory,
    reader/writer traits and trajectory errors

smiles / molfile / sdf / mmcif
    format-specific documents, interpretation, reports, and writers

perception
    explicit valence, ring, aromaticity, stereo, and sanitization workflows

hydrogens
    explicit hydrogen topology transformations

descriptors
    read-only molecular descriptors with explicit policies

query
    syntax-neutral chemical query graphs and bounded SMARTS parsing

substructure
    query-graph matching algorithms

canon
    canonicalization algorithms

dssp
    read-only DSSP assignment over structural views

alignment
    same-topology selection-based weighted proper rigid fitting

modeling
    prepared-system interfaces, Potential, minimization, numerical workflows
```

Implementation modules may be private and organized differently, but dependency
direction must respect these responsibilities.

The prelude contains only genuinely common domain types. Specialized reports,
format internals, trajectory codecs, modelling objects, and expert algorithms
remain in focused namespaces. Broad root re-exports are not added casually.

Dependency-heavy binary trajectory codecs and force-field adapters should live
in separate crates when required, keeping the foundational `molecular` crate
lightweight.

Production fixed-topology file codecs live in the one-way
`molecular-trajectory-io` workspace companion:

```text
molecular <- molecular-trajectory-io <- applications
```

The companion implements Molecular's streaming traits and uses its topology,
frame-buffer, geometry, unit, format-identity, and typed error contracts. It
does not define duplicate domain objects. Molecular-owned codec source forbids
unsafe code and does not require Chemfiles, a C/C++ compiler, or CMake.

The supported initial companion profile is intentionally explicit:

| Format | Compatibility profile |
|---|---|
| XYZ | strict constant-count multi-frame element/x/y/z text; configured units, angstrom by default |
| DCD | common CHARMM/NAMD/OpenMM 32-bit-record `CORD`, either byte order, common cells, fixed-atom reconstruction, strict `NSET` |
| TRR | GROMACS XDR f32/f64 positions, triclinic cell, optional velocities/forces, time, nonnegative step, explicit lambda policy, cumulative sequential precision metadata |
| XTC | GROMACS 1995/2023 magic, signed nonnegative i32 counts/steps, small uncompressed and ordinary compressed coordinates, explicit lossy precision |

DCD and default XYZ coordinates use angstrom conventions; TRR and XTC use
GROMACS nanometre/picosecond conventions and convert once to Molecular units.
XTC reports nominal lossy resolution as `1 / precision` nanometres. Indexed
opening is O(file size), verifies every complete frame, and stores bounded
checked offsets; indexed frame reads decode one frame. Sequential opening
retains one handle and does not scan the whole file.

All file-controlled lengths and offsets are checked before allocation or seek.
Codecs validate complete finite dense-order state in reusable scratch before
transactional publication. Path writers publish only after consuming finish
flushes, synchronizes, and finalizes metadata; any failed frame write prevents
publication. Exact frame/index limits use a bounded frame-start/EOF probe
without decoding the next frame or growing the index. Indexed reads restore
the stream and all sequential codec state before destination publication.
Unsupported historical dialects and deferred formats remain typed errors
rather than inferred or partially decoded behavior.

## Public API and release policy

The initial release line is `0.x`. Breaking public API changes require a minor
version increment. The topology-centered refactor is a deliberate breaking
transition and targets the next minor release.

Invariant-bearing topology, hierarchy, coordinate, document, provenance, and
error state is private behind accessors and checked constructors. Public direct
fields are reserved for deliberate value, options, and report payloads.

Extensible public error enums are `#[non_exhaustive]`.

Parsing, interpretation, sanitization, perception, topology construction,
coordinate construction, preparation, analysis, and writing remain visibly
separate in names and documentation.

Compatibility shims may exist during the staged refactor, but the final public
surface must use `Topology`, `TopologyAtomIndex`, topology-bound prepared
objects, and shared structural views. Deprecated aliases are removed before the
new minor release unless a deliberate deprecation release is chosen.

## Architectural invariants

The following statements summarize the design:

1. `Molecule` is one asserted local chemical entity; connectedness is not its
   defining invariant.
2. `Topology` is one immutable coordinate-free molecular system.
3. Topology preserves molecule definitions and molecule-instance boundaries.
4. Semantic identities and dense numerical indices are separate.
5. One topology has one immutable authoritative dense ordering.
6. `Model` is one topology plus one complete configuration.
7. `Ensemble` is one topology plus finite non-temporal members.
8. `Trajectory` is one topology plus ordered frames and supports streaming I/O.
9. Periodic cells, velocities, forces, time, and observation data are dynamic
   state, not topology.
10. Prepared systems and compiled selections bind to exact topology identity.
11. Topology-changing operations return a new topology and explicit mappings.
12. Reactive trajectories are segmented by topology rather than weakening
    fixed-topology containers.
13. Parsing, interpretation, perception, sanitization, topology construction,
    coordinate construction, and preparation are explicit distinct operations.
14. Canonical chemistry and structure containers never silently absorb
    backend-specific mechanical state.
