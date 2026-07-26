# Topology-centered refactor plan

## Purpose

This document is the staged implementation plan for migrating `molecular` to
the topology-centered architecture defined in
[`ARCHITECTURE.md`](ARCHITECTURE.md).

The refactor promotes `Topology` into the central immutable system-level object
shared by:

- `Model`;
- `Ensemble`;
- `Trajectory`;
- prepared force fields and backend systems;
- topology-bound selections;
- coordinate-dependent analyses.

The target relationships are:

```text
Model       = Topology + one Configuration
Ensemble    = Topology + finite non-temporal members
Trajectory  = Topology + ordered Frames
```

This is a deliberate pre-1.0 breaking refactor. The final implementation should
prefer a clean public API over preserving obsolete names indefinitely.

## Required outcome

At completion, the repository must provide:

1. A first-class public `Topology` type that is independently constructible
   without coordinates.
2. Explicit molecule definitions and molecule instances inside topology.
3. Instance-qualified semantic atom and bond identifiers.
4. Immutable authoritative dense atom and bond orderings owned by topology.
5. Cheap topology cloning with exact identity retained across clones.
6. Explicit distinction between exact topology identity and structural
   equivalence.
7. A linear-time transactional `TopologyBuilder`.
8. `Positions` and `Configuration` types validated against topology.
9. `Model` implemented as topology plus one configuration.
10. `Ensemble` implemented as topology plus finite non-temporal members.
11. Trajectory frame types and streaming reader/writer contracts over one
    topology.
12. Borrowed structural views usable by models, ensemble members, trajectory
    frames, DSSP, and potentials.
13. Prepared potentials bound to topology identity rather than model identity.
14. Static macromolecular hierarchy separated from coordinate-observation data.
15. Updated mmCIF interpretation and writing consistent with the new object
    model.
16. Updated `molecular-dreiding`, tests, feature contracts, examples,
    documentation, and downstream MolStudio integration.
17. A minor version increment for the breaking public API.

## Non-goals for this refactor

The following are explicitly outside the required initial implementation unless
they become necessary to complete the architecture cleanly:

- production DCD, XTC, TRR, NetCDF, or other binary trajectory codecs;
- a reactive trajectory engine;
- automatic structural interning or deduplication of molecule definitions;
- a complete structural selection language;
- molecular dynamics integration;
- neighbor-list implementation;
- force-field constraints, virtual sites, or Drude particles in canonical
  topology;
- generalized sparse ensembles with inconsistent atom presence;
- automatic coordinate generation for atoms introduced by topology edits;
- broad performance rewrites unrelated to topology and coordinate ownership.

The core trajectory interfaces, in-memory trajectory representation, and
reusable frame buffers are in scope. Dependency-heavy codecs may follow in
adapter crates.

## Global invariants

Every stage must preserve these invariants:

- Parsing, interpretation, perception, sanitization, topology construction,
  coordinate construction, preparation, analysis, and writing remain explicit
  separate operations.
- `Molecule` remains the local asserted chemical graph kernel.
- A `Molecule` may be disconnected.
- Molecule instances are not conflated with connected components.
- `Topology` is coordinate-free and immutable after construction.
- Positions, periodic cells, velocities, forces, time, occupancy, and B-factors
  are not topology state.
- Force-field and backend particle state is not canonical topology state.
- Exact topology identity is required for topology-bound dense arrays,
  selections, frame buffers, and prepared systems.
- Independently constructed structurally equal topologies are not silently
  interchangeable.
- Dense atom and bond orderings never change during a topology's lifetime.
- Topology-changing operations create a new topology and explicit mappings.
- Failed transactional operations leave their input unchanged.
- Public constructors reject malformed, incomplete, non-finite, incompatible,
  or overflowing state with structured errors.
- No public fixed-width identifier is created through unchecked narrowing.
- No step introduces `unsafe` code.
- No implementation stage hides sanitization, perception, or parameterization
  inside a parser or default constructor.
- Writers reject unsupported semantics rather than silently coercing them.

## Current-to-target type mapping

The current modelling layer maps to the target architecture as follows:

| Current | Target |
|---|---|
| `ModelTopology` | public `Topology`/private `TopologyData` |
| `ModelDefinition` | removed; its ownership role moves into `TopologyData` |
| `ModelDefinitionKey` | removed; replaced by `Topology` identity or `TopologyIdentity` |
| `ModelAtomIndex` | `TopologyAtomIndex` |
| no dense bond index | `TopologyBondIndex` |
| `MoleculeInstancePayload` owned by every instance | reusable `MoleculeDefinition` referenced by `MoleculeInstance` |
| `ModelBuilder` creates topology and coordinates together | independent `TopologyBuilder` plus model convenience builders |
| `Model { definition, positions }` | `Model { topology, configuration, ... }` |
| `Quantity<Vec<Point3>>` directly in `Model` | validated `Positions` inside `Configuration` |
| potential bound to `ModelDefinitionKey` | potential bound to exact topology identity |
| DSSP accepts only `&Model` | DSSP kernel accepts borrowed structural view |
| SMCRA model node and coordinate metadata in hierarchy | static hierarchy plus configuration-associated observation/provenance |
| no ensemble | `Ensemble` and `EnsembleMember` |
| no trajectory | frame types, in-memory trajectory, frame buffer, streaming traits |

Temporary aliases may be used inside an intermediate commit if they keep the
workspace buildable, but the final public API must use target names.

## Expected target module layout

The final public surface should approximate:

```text
crates/molecular/src/
  core/
  geometry/
  topology/
    mod.rs
    ids.rs
    definition.rs
    instance.rs
    builder.rs
    mapping.rs
    error.rs
  structure/
    mod.rs
    positions.rs
    configuration.rs
    model.rs
    ensemble.rs
    observation.rs
    view.rs
  trajectory/
    mod.rs
    frame.rs
    buffer.rs
    memory.rs
    reader.rs
    writer.rs
    error.rs
  modeling/
    potential.rs
    minimize.rs
    prepared.rs
```

Exact private file boundaries may differ. Public responsibility boundaries from
`ARCHITECTURE.md` are authoritative.

## Development strategy

Use one dedicated branch, preferably:

```text
codex/topology-refactor
```

Open a draft PR early. Keep commits ordered by the stages below. Each commit
should leave the workspace compiling and should include corresponding tests and
feature-contract updates where practical.

Because the change is broad, a single draft PR with reviewable staged commits
is acceptable. If the patch becomes too large to review, split it into stacked
PRs at the stage boundaries while preserving the dependency order.

Before editing:

1. Read `AGENTS.md`, `ARCHITECTURE.md`, and this document.
2. Inventory all canonical feature IDs under `features/` affected directly or
   through public dependencies.
3. Inspect all uses of:
   - `ModelTopology`;
   - `ModelDefinitionKey`;
   - `ModelAtomIndex`;
   - `ModelBuilder`;
   - `MoleculeInstance`;
   - `Model::topology`;
   - `Potential`;
   - `Dssp`;
   - mmCIF interpretation/writing;
   - `molecular-dreiding`;
   - public examples and generated feature metadata.
4. Record the affected feature IDs in the PR description.
5. Identify MolStudio call sites through the documented sibling-repository
   workflow when `../molstudio` is available.
6. Do not modify generated dashboard output by hand.

## Stage 0: Baseline, contracts, and migration scaffolding

### Goal

Establish a verified baseline and define the feature-contract changes before
moving public types.

### Work

- Run and record the current applicable checks.
- Inventory affected feature directories.
- Add new canonical feature contracts or update existing ones for:
  - topology construction and identity;
  - model coordinate state;
  - ensemble storage;
  - trajectory frames and streaming;
  - topology-bound potential evaluation;
  - mmCIF model/ensemble interpretation as applicable.
- Mark planned pieces accurately rather than claiming implementation.
- Add architecture-focused compile tests demonstrating the intended final API.
  They may begin ignored or gated only if the repository's test policy permits;
  otherwise add them at the stage that makes them compile.
- Decide the public compatibility policy:
  - preferred: hard break within the `0.2.0` branch and update all consumers;
  - acceptable temporary strategy: private/internal aliases only;
  - avoid long-lived public aliases that duplicate the conceptual model.
- Add a changelog entry describing the planned breaking topology transition.

### Required tests

- Existing full workspace suite remains green before functional edits.
- Feature registry and docs generation checks remain green.

### Exit criteria

- Affected feature IDs are known.
- No feature contract contradicts the target architecture.
- The PR documents baseline commands and any unavailable external data.

## Stage 1: Extract shared geometry primitives

### Goal

Move general three-dimensional types out of conformer- or potential-specific
locations so structure and trajectory layers can share them.

### Work

- Create public `geometry` module.
- Move or re-home:
  - `Point3`;
  - `Vector3`.
- Add only the minimum required support for:
  - vector arithmetic;
  - point/vector distinction;
  - finite-value checks;
  - unit scaling.
- Add `Matrix3` and `PeriodicCell` if needed for the configuration API.
- Validate periodic cells:
  - finite vectors;
  - non-zero volume;
  - supported periodic axes;
  - explicit orthorhombic/triclinic representation.
- Update conformers, potentials, minimization, DSSP kernels, and DREIDING to use
  common geometry types.
- Avoid introducing a large third-party linear-algebra dependency unless
  justified by an existing project decision.

### Compatibility

A temporary re-export of `Point3` from `core` may be retained during migration,
but the final canonical path is `geometry::Point3`.

### Required tests

- Point/vector arithmetic semantics.
- Unit conversions for vectors and points.
- Periodic-cell validation.
- Existing minimization and force-field tests.

### Exit criteria

- No potential-specific vector type remains.
- General structure code depends on `geometry`, not on `modeling::potential`.

## Stage 2: Introduce independent `Topology`

### Goal

Extract the coordinate-free system definition from the current model
implementation without yet requiring molecule-definition reuse.

### Work

Create public topology identifiers:

```rust
MoleculeDefinitionId
MoleculeInstanceId
InstanceAtomId
InstanceBondId
TopologyAtomIndex
TopologyBondIndex
```

Create:

```rust
#[derive(Clone)]
pub struct Topology {
    inner: Arc<TopologyData>,
}
```

Initial `TopologyData` should include:

- molecule definitions;
- molecule instances;
- authoritative atom order;
- authoritative bond order;
- semantic-ID-to-dense-index mappings;
- dense-index-to-semantic-ID mappings;
- system-level properties where justified.

Move topology graph access from `ModelTopology` to `Topology`:

- molecule definition and instance lookup;
- atom and bond lookup;
- atom and bond iteration;
- neighbor and incident-bond access;
- dense index conversion;
- molecule and definition membership;
- qualified hierarchy views;
- perception-state queries through definitions.

Provide explicit identity behavior:

```rust
Topology::same_identity
Topology::identity
```

Introduce a private or opaque `TopologyIdentity` suitable for prepared objects,
positions, selections, and buffers.

Do not rely on `PartialEq` to express both identity and structural equality.
Provide an explicit structural-equivalence API if implemented in this stage;
otherwise reserve it clearly for a later stage.

### Initial storage strategy

It is acceptable in this stage for every instance to have a distinct
definition. The public shape must still distinguish definitions from instances
so explicit reuse can be added in Stage 3 without changing semantics.

### Required tests

- Topology can be constructed without coordinates.
- Cloned topology retains exact identity.
- Independently built equal topology has distinct identity.
- Local atom and bond IDs survive insertion.
- Instance qualification disambiguates repeated local IDs.
- Dense mappings are complete, bidirectional, and stable.
- Atom/bond iteration order matches documented topology order.
- Invalid identifiers return structured errors.
- No coordinate field is present in topology.

### Exit criteria

- `ModelTopology` no longer owns the canonical implementation.
- All coordinate-independent system queries are available from `Topology`.
- `Topology` compiles independently of `Model`.

## Stage 3: Add molecule definitions, reusable instances, and scalable building

### Goal

Complete the topology ownership model and eliminate current quadratic builder
behavior.

### Work

Implement explicit molecule definitions:

```rust
pub enum MoleculeDefinitionPayload {
    Small(SmallMolecule),
    Macro(MacroMolecule),
}

pub struct MoleculeDefinition { ... }

pub struct MoleculeInstance {
    id: MoleculeInstanceId,
    definition: MoleculeDefinitionId,
    metadata: MoleculeInstanceMetadata,
}
```

Requirements:

- Definitions are conformer-free.
- Definitions preserve local atom/bond IDs and coordinate-independent
  perception state.
- Instances carry roles and instance-specific scalar annotations.
- One definition may be referenced by many instances.
- No automatic graph-equality interning.
- Accessors allow instance-to-definition and definition-to-instance traversal.
- Graph and hierarchy access through an instance resolve its definition.

Implement `TopologyBuilder`:

```rust
builder.add_small_molecule_definition(...)
builder.add_macro_molecule_definition(...)
builder.add_instance(definition, metadata)
builder.add_small_molecule_instance(...)
builder.add_macro_molecule_instance(...)
builder.reserve_*
builder.build()
```

Use validate-then-commit transactionality:

1. Validate the borrowed definition or instance addition.
2. Construct only the new staged records and mappings.
3. Reserve required capacity.
4. Perform infallible append.
5. Leave the builder unchanged on validation failure.

Remove whole-builder cloning from each addition.

Do not clone all conformers merely to strip them. Add coordinate-free clone or
move helpers on molecular wrappers as needed.

Use checked conversions for all `u32` IDs and indices.

Document authoritative dense ordering. A suitable initial rule is:

1. molecule instance insertion order;
2. live local atom order within the referenced definition;
3. live local bond order within the referenced definition.

The ordering is topology-authoritative, not a claim of canonical chemical
labeling.

### Performance validation

Add focused tests or benchmarks demonstrating:

- many instances of one water definition do not duplicate the molecular graph;
- builder additions do not clone accumulated state;
- construction scales linearly in appended atoms and instances within ordinary
  measurement noise;
- mmCIF systems with many solvent instances avoid quadratic construction.

Performance benchmarks are informational, but implementation structure must
make the asymptotic behavior evident in code review.

### Required tests

- Explicit definition reuse.
- Two instances of one definition have distinct qualified atom IDs.
- Roles and instance metadata remain instance-specific.
- Definition data is shared or stored once as designed.
- Invalid definition references fail transactionally.
- Capacity overflow is structured.
- Empty topology and empty definitions are rejected.
- Invalid macro graph/hierarchy pair is rejected.
- Source molecule and conformers remain unchanged.

### Exit criteria

- Repeated molecules are representable without repeated graph payloads.
- `TopologyBuilder` has no clone-of-self transaction pattern.
- Large solvent-rich topologies are practical.

## Stage 4: Introduce `Positions`, `Configuration`, and the new `Model`

### Goal

Make coordinate state independently validated and implement the central
`Model = Topology + Configuration` relationship.

### Work

Create `Positions`:

- complete dense `Vec<Point3>`;
- exact topology binding or checked construction against topology;
- canonical internal length unit;
- finite-coordinate validation;
- immutable length;
- read-only quantity access;
- transactional full replacement;
- efficient no-allocation or allocation-reusing update path.

Create:

```rust
pub struct Configuration {
    positions: Positions,
    cell: Option<PeriodicCell>,
}
```

Provide `ConfigurationView`.

Refactor `Model` to:

```rust
pub struct Model {
    topology: Topology,
    configuration: Configuration,
    observation: Option<StructureObservation>,
}
```

The observation field may be deferred until Stage 5 if needed, but the model
must not put coordinate-specific source values into topology.

Provide:

```rust
Model::new(topology, configuration)
Model::topology
Model::configuration
Model::positions
Model::set_positions
Model::cell
Model::set_cell
Model::view
```

Retain convenience construction from:

- one `SmallMolecule` conformer;
- one `MacroMolecule` conformer;
- a builder that assembles topology and initial positions together.

The convenience builder must internally produce a standalone `Topology` and
validated `Configuration`; it must not restore the old inseparable ownership.

Replace `ModelAtomIndex` with `TopologyAtomIndex` throughout modelling,
gradients, examples, writers, and adapters.

Remove `ModelDefinition`, `ModelDefinitionKey`, and the old model-owned
topology implementation once all internal call sites migrate.

### Position copying

Update instance-to-conformer copying to use topology instance mappings.
Validate complete local ID compatibility before changing the target conformer.
Leave targets unchanged on failure.

### Required tests

- Model created from independently built topology and positions.
- Wrong topology identity is rejected even when atom count matches, where
  identity-bearing positions are used.
- Wrong length, non-finite coordinates, and incompatible units are rejected.
- Topology remains shared after model clone.
- Coordinate mutation preserves topology identity.
- Periodic cell may change without changing topology.
- Convenience constructors preserve source objects.
- Repeated full updates can reuse allocation.
- Instance-to-conformer mapping remains transactional.

### Exit criteria

- The canonical model contains public `Topology`, not hidden model definition.
- Topology can be reused to construct multiple independent models.
- Prepared objects can begin binding directly to topology.

## Stage 5: Separate static hierarchy from observation state

### Goal

Make macromolecular hierarchy genuinely coordinate-independent so it can be
shared by models, ensembles, and trajectories.

### Work

Refactor `SmcraHierarchy` toward:

```text
chain
  residue
    atom site -> local AtomId
```

Remove coordinate-model nodes as structural parents in the final API. Preserve
source model identifiers in interpretation provenance or
configuration/member-level metadata.

Split `SmcraAtomSiteMetadata` into static and dynamic parts.

Static topology candidates include:

- atom names;
- chain label/author IDs;
- residue label/author IDs;
- component IDs;
- insertion codes;
- element/type assertions where chemically consistent;
- static atom-site identity when it is stable across members.

Dynamic observation candidates include:

- occupancy;
- B-factor;
- alternate-location choice;
- raw Cartesian fields;
- coordinate-model ID;
- other model-specific source values.

Create typed `StructureObservation` with topology-bound per-atom arrays or
records. Validate all lengths and finite numeric values.

Update:

- `MacroMoleculeBuilder`;
- macro validation;
- mmCIF interpretation report/provenance;
- mmCIF writer;
- DSSP hierarchy extraction;
- existing benchmark/reference outputs where semantics intentionally change.

Do not drop source information merely to simplify the canonical hierarchy.
Preserve it in documents, reports, or observation state.

### Migration caution

This stage changes public hierarchy identity and likely benchmark output.
Update feature contracts and versioning explicitly. Do not regenerate goldens
without verifying that the new output reflects the intended contract.

### Required tests

- Topology hierarchy contains no coordinate arrays.
- One hierarchy can be shared by multiple configurations.
- Per-model occupancy and B-factors remain distinct.
- Source model IDs remain available outside topology.
- Macro validation remains transactional.
- mmCIF single-model round trip preserves the documented supported subset.

### Exit criteria

- Static hierarchy is reusable across ensemble members.
- Dynamic structure observations no longer affect topology identity.

## Stage 6: Introduce borrowed structural views

### Goal

Allow coordinate-dependent algorithms to operate on any topology plus
configuration without allocating an owned model.

### Work

Create:

```rust
pub struct ModelView<'a> {
    topology: &'a Topology,
    configuration: ConfigurationView<'a>,
}
```

Provide view creation from:

- `Model`;
- `EnsembleMember`;
- `TrajectoryFrame`;
- `FrameBuffer`.

Consider a narrower trait only if it remains object-safe, simple, and does not
hide topology identity. Prefer concrete borrowed views over an abstract trait
hierarchy unless multiple backends require it.

Move coordinate-dependent read-only kernels to views:

- DSSP;
- geometric utilities;
- potential evaluation;
- future RMSD/alignment/contact kernels.

Keep convenience wrappers accepting `&Model` when useful:

```rust
pub fn assign(model: &Model, options: ...) {
    assign_view(model.view(), options)
}
```

### Required tests

- View results equal owned-model results.
- View construction does not copy coordinates.
- A frame buffer can be analyzed directly.
- Topology identity is preserved through views.
- Snapshot analysis results do not change after source coordinates mutate.

### Exit criteria

- No core coordinate-dependent kernel requires ownership of `Model`.
- Trajectory frames can be evaluated without coordinate copying.

## Stage 7: Rebind potentials and prepared systems to topology

### Goal

Prepare once per topology and evaluate any supported compatible configuration.

### Work

Change `Potential` or its preparation contract so implementations bind to exact
`TopologyIdentity`, not old model-definition identity.

Preferred evaluation shape:

```rust
pub trait Potential {
    fn evaluate(
        &mut self,
        model: ModelView<'_>,
    ) -> Result<PotentialEvaluation, PotentialError>;
}
```

Update `PotentialEvaluation` to associate gradients with
`TopologyAtomIndex`/topology identity as appropriate.

Migrate:

- harmonic bond potential;
- minimization;
- `molecular-dreiding`;
- examples and docs.

Preparation must accept topology when parameters depend only on topology. If a
specific parameterization method legitimately depends on one reference
geometry, make that dependency explicit in its preparation API rather than
pretending it is topology-only.

For DREIDING:

- build atom types and bonded terms from topology;
- make coordinate use explicit only where the external parameterizer truly
  requires it;
- bind the result to topology identity;
- evaluate any supported model/ensemble member/frame with that topology;
- document the adapter as nonperiodic and reject periodic reference and
  evaluation views until a complete periodic policy exists;
- preserve fixed-charge behavior during evaluation.

Clarify QEq scope. Do not call molecule-instance-local calculation
“component-local” unless it actually uses connected components. Introduce an
explicit scope policy if necessary:

```rust
WholeTopology
MoleculeInstances
ConnectedComponents
ExplicitGroups
```

### Minimization

Minimization may continue to clone a `Model` and return a new model. The
potential kernel should evaluate views. Coordinate work buffers should reuse
allocation where practical.

### Required tests

- One prepared potential evaluates two independent models sharing one topology.
- It evaluates ensemble members and trajectory frame views.
- Built-in and adapter potentials document periodic capability and reject
  unsupported periodic model, ensemble, frame, and frame-buffer views.
- It rejects independently constructed topology with equal content.
- Gradient indexing uses `TopologyAtomIndex`.
- DREIDING preparation does not mutate topology or input models.
- QEq grouping semantics are explicitly tested.
- Existing energy/gradient regression tests remain valid.

### Exit criteria

- `ModelDefinitionKey` has no remaining use.
- Prepared systems are reusable across all coordinate containers sharing
  topology.

## Stage 8: Add `Ensemble`

### Goal

Add a finite non-temporal multi-configuration container.

### Work

Create:

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

Requirements:

- exact shared topology identity;
- complete finite coordinates;
- stable member order;
- no implicit temporal semantics;
- finite non-negative optional weights;
- explicit weight normalization;
- topology-compatible observation arrays;
- member access by stable dense member index;
- iterator yielding borrowed model views.

Provide constructors from selected local molecule conformers and from multiple
models sharing topology where useful.

Do not silently merge independently constructed topologies based only on
content. Provide an explicit reconciliation/mapping operation if needed.

### Required tests

- Multiple configurations share one topology.
- Members with wrong topology are rejected.
- Weights validate and normalize explicitly.
- Member order is stable.
- Ensemble views work with DSSP and potentials.
- Conversion from local conformers preserves ordering and source molecules.

### Exit criteria

- NMR models, conformer sets, and docking poses have a canonical non-temporal
  container.
- Ensemble is not implemented as `Vec<Model>` with repeated topology handles
  unless that is a temporary private representation.

## Stage 9: Add trajectory frames, buffers, and in-memory trajectory

### Goal

Add the fixed-topology temporal/sequential data model independently of file
codecs.

### Work

Create topology-bound array types or validated wrappers for:

- velocities;
- forces.

Create:

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

Create reusable:

```rust
pub struct FrameBuffer { ... }
pub struct TrajectoryFrameView<'a> { ... }
```

`FrameBuffer` must:

- bind to one exact topology;
- allocate required position storage once;
- reuse optional array capacity;
- validate state after decoding;
- expose a borrowed model view;
- reject topology mismatch.

Create an in-memory `Trajectory` for deliberately loaded finite data:

```rust
pub struct Trajectory {
    topology: Topology,
    frames: Vec<TrajectoryFrame>,
}
```

Requirements:

- stable ordered frames;
- exact shared topology;
- complete finite positions;
- optional monotonic-time validation only when explicitly requested;
- no claim that every trajectory has time or step metadata;
- iterator yielding frame views/model views.

Do not define trajectory as `Vec<Model>`.

### Required tests

- Frame validation for positions, velocities, forces, time, step, and cell.
- Buffer reuse without reallocation in common sequential reads.
- Trajectory order is stable.
- Wrong topology frames are rejected.
- Potential and DSSP evaluation directly over frames.
- Variable periodic cells do not alter topology.
- Missing optional velocity/force arrays are valid.
- Partial arrays are rejected.

### Exit criteria

- The in-memory fixed-topology trajectory data model is complete.
- No file codec is needed to test the architecture.

## Stage 10: Add streaming trajectory interfaces

### Goal

Support large trajectory I/O without loading all frames.

### Work

Introduce:

```rust
pub trait TrajectoryReader
pub trait SeekableTrajectoryReader: TrajectoryReader
pub trait TrajectoryWriter
```

Use the contracts from `ARCHITECTURE.md`.

Reader requirements:

- expose topology before frames;
- read into caller-owned `FrameBuffer`;
- distinguish end-of-stream from error;
- report structured format/compatibility errors;
- never silently resize a topology-bound buffer to a different atom count;
- document units and precision;
- preserve available time, step, cell, velocity, and force data.

Seekable-reader requirements:

- separate optional frame count from sequential capability;
- define zero-based frame indices;
- report unsupported or unavailable random access structurally.

Writer requirements:

- bind to topology;
- validate each frame;
- reject unsupported fields or document deliberate omission policy;
- never silently reorder atoms.

Provide a simple in-memory/reference reader and writer for tests. A minimal
plain-text or internal test codec may be added, but production binary codecs are
not required.

### Topology-free trajectory sources

Add constructors requiring external topology for coordinate-only sources.
Require explicit atom-order assertion and validate all available metadata.

### Required tests

- Sequential read with one reusable buffer.
- Clean end-of-stream behavior.
- Random-access capability separated from sequential reader.
- Topology mismatch and atom-count mismatch.
- Writer rejection of unsupported state.
- Round trip through the reference codec.
- No per-frame topology clone or coordinate-copy requirement.

### Exit criteria

- Large trajectories can be processed frame-by-frame.
- Core APIs do not require a `Vec<TrajectoryFrame>`.

## Stage 11: Migrate mmCIF and structure I/O

### Goal

Align structural interpretation and writing with independent topology,
configuration, ensemble, and observation state.

### Single-model interpretation

Refactor ordinary mmCIF interpretation to construct:

```text
Topology
Configuration
StructureObservation, when present
MmcifInterpretationReport
```

Expose a convenience `Model`.

Preserve explicit policies for:

- coordinate-model selection;
- alternate-location selection;
- entity classification;
- connection interpretation;
- unresolved or inferred connectivity.

Do not assert every distance-inferred bond as authoritative single bond.
Connectivity inference and bond-order assignment must remain evidence-backed
and visible in reports.

### Multi-model interpretation

Add a separate ensemble interpretation path after the single-model path is
stable.

It must:

1. identify the intended coordinate-model records;
2. interpret static chemistry and molecule boundaries;
3. prove consistent atom identity across members;
4. prove consistent molecule-instance partition and topology;
5. construct one shared topology;
6. construct one ensemble member per coordinate model;
7. attach source model ID and dynamic observation values to members;
8. reject inconsistent atom presence or connectivity with a structured error
   unless an explicit future reconciliation policy is selected.

Do not change the default single-model interpretation into implicit ensemble
loading.

### Writers

Update mmCIF writing to accept the new `Model` and, if implemented, explicit
ensemble writing.

Writers must preserve the documented supported subset and reject unsupported
topology or observation semantics.

### Other formats

Update molfile/SDF conveniences that construct models or copy conformers.
Standalone molecule parsing and writing remain based on local molecule
conformers.

### Required tests

- Existing single-model mmCIF cases.
- Explicit model selection.
- Alternate-location policy.
- Topology/positions/provenance separation.
- Multi-model consistent ensemble.
- Multi-model inconsistent atom set rejection.
- Occupancy/B-factor separation by member.
- Writer supported-subset rejection tests.
- No parser invokes topology perception or force-field preparation implicitly.

### Exit criteria

- mmCIF no longer depends on model-owned topology internals.
- Multi-model structures have an explicit ensemble path.

## Stage 12: Add topology-bound selections and mapping foundations

### Goal

Provide reusable atom selections and explicit topology transfer mechanisms
needed by trajectory analysis and topology edits.

### Work

Create `AtomSelection` and, if useful, `BondSelection`:

- bind to exact topology identity;
- store sorted unique dense indices;
- preserve requested stable order when an API explicitly needs it;
- construct from semantic IDs, roles, definitions, instances, hierarchy labels,
  elements, connected components, and chemical query matches;
- expose iteration over dense and semantic IDs;
- reject incompatible topology use.

Keep chemical SMARTS/query parsing separate from structural selection
semantics.

Implement `TopologyMapping` foundations:

- old/new topology identity;
- definition mapping;
- instance mapping;
- atom mapping;
- bond mapping;
- dense-index mapping;
- retained/removed/added reporting.

It is acceptable for the initial refactor to provide mappings for builder
conversions and selected transforms rather than every possible topology edit.

### Required tests

- Selection reuse over all frames.
- Selection topology mismatch rejection.
- Hierarchy and role selections.
- Query-derived selections.
- Mapping round trips for retained atoms.
- Added/removed atom reporting.
- No implicit position transfer without mapping.

### Exit criteria

- Trajectory analyses can compile selections once.
- Future topology-changing operations have a stable lineage vocabulary.

## Stage 13: Update adapters, examples, and MolStudio

### Goal

Migrate all repository and known consumer call sites to the new public API.

### Work

Update:

- `molecular-dreiding`;
- README examples if owner-approved for any required changes;
- crate-level rustdoc;
- examples and doctests;
- feature documentation;
- changelog;
- roadmap;
- benchmark adapters;
- fuzz targets if public constructors changed;
- MolStudio through the documented sibling workflow.

MolStudio integration should use:

- `Topology` for immutable system structure;
- `TopologyAtomIndex` for render buffers where appropriate;
- molecule-instance and hierarchy identities for semantic selection;
- `Model` for one displayed configuration;
- trajectory frame buffers for streamed playback;
- no flattened synthetic molecule as the application-wide system model.

If `../molstudio` is unavailable, report that consumer validation was not run
and provide the exact expected migration notes.

### Required tests

Run Molecular checks first, then MolStudio checks against the local path patch
when available.

### Exit criteria

- No internal repository call site uses removed public names.
- Known consumer code compiles against the new API.
- Examples describe the architecture accurately.

## Stage 14: Public API cleanup and release transition

### Goal

Remove migration residue and finalize the next minor release contract.

### Work

- Remove public:
  - `ModelTopology`;
  - `ModelDefinition`;
  - `ModelDefinitionKey`;
  - `ModelAtomIndex`;
  - old model-owned topology constructors;
  - obsolete SMCRA coordinate-model hierarchy API;
  - obsolete potential signatures.
- Remove temporary compatibility modules not deliberately retained.
- Ensure public names follow the target module responsibilities.
- Keep the prelude focused.
- Mark extensible error enums `#[non_exhaustive]`.
- Audit direct public fields on invariant-bearing structs.
- Bump workspace version from `0.1.0` to `0.2.0`, unless a later current version
  requires the corresponding next minor.
- Update all inter-workspace dependency versions.
- Update `CHANGELOG.md`.
- Update canonical feature metadata versions for every changed public contract.
- Regenerate procedural dashboard/skills output through `xtask`; never hand-edit
  generated files.
- Audit rustdoc for obsolete terminology:
  - “model topology” where the type is now `Topology`;
  - “model atom index” where the type is now `TopologyAtomIndex`;
  - “component-local” where molecule-instance scope was meant;
  - trajectory claims that imply all frames are in memory.

### Required tests

```bash
cargo fmt --all -- --check
cargo check --workspace --all-targets --all-features --locked
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo test --workspace --all-features --locked
cargo test --workspace --all-features --doc --locked
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --all-features --no-deps --locked
cargo package -p molecular --locked
cargo package -p molecular-dreiding --locked --list
cargo xtask dashboard --check
cargo xtask skills --check
```

Run targeted fuzz smoke tests for parsers affected by constructor or
interpretation changes.

Run targeted external-reference benchmarks only when their data is available
and the changed feature has a registered manifest. Benchmarks remain
informational and are not release gates.

Run MolStudio workspace checks through the sibling patch workflow when
available.

### Exit criteria

- Full workspace checks pass.
- No obsolete type remains in the public API.
- Documentation, feature contracts, implementation, and examples agree.
- The release version reflects the breaking API change.
- Every omitted applicable command is reported with reason.

## Cross-stage test matrix

The final suite should cover at least the following scenarios.

### Topology identity

- clone shares identity;
- independent equal construction does not;
- prepared potential accepts shared identity;
- prepared potential rejects independent identity;
- selection and frame buffer enforce the same rule.

### Molecule definitions and instances

- one definition, many solvent instances;
- distinct instance roles;
- disconnected molecule definition;
- connected-component queries distinct from instance queries;
- macro and small definitions in one topology;
- covalently linked entities represented within one definition.

### Dense indexing

- complete atom and bond inverse mappings;
- stable order after topology clone;
- no tombstones in dense arrays;
- local tombstone IDs still preserved;
- checked overflow behavior.

### Coordinate containers

- complete finite positions;
- unit conversion;
- topology mismatch;
- periodic and non-periodic configurations;
- variable cells across frames;
- no dynamic state in topology.

### Model, ensemble, and trajectory

- several models sharing one topology;
- ensemble with weights and observation data;
- ordered trajectory frames;
- frame streaming through reusable buffer;
- views used by analysis and potential kernels;
- no `Vec<Model>` implementation dependency.

### Hierarchy and mmCIF

- static chain/residue/atom-site identity;
- per-member model ID, occupancy, and B-factor;
- alternate-location selection;
- single-model interpretation;
- consistent multi-model ensemble;
- inconsistent multi-model rejection;
- explicit provenance and issue reports.

### Prepared systems

- DREIDING preparation once, evaluation many configurations;
- explicit nonperiodic DREIDING preparation/evaluation contract;
- explicit charge grouping scope;
- topology mismatch;
- no model mutation during preparation or evaluation.

### Transactionality

- invalid builder addition leaves builder unchanged;
- invalid position update leaves positions unchanged;
- invalid macro edit leaves macro molecule unchanged;
- invalid topology remap leaves source objects unchanged;
- writer errors do not partially mutate canonical inputs.

## Performance requirements

The refactor should not merely rename types. It must correct known scalability
hazards.

Required implementation properties:

- adding a topology instance does not clone the accumulated builder;
- adding one selected conformer does not clone all source conformers;
- repeated molecule definitions can be stored once;
- topology clone is constant-time with respect to atom count;
- model clone copies coordinate state but not topology data;
- frame decoding reuses a caller-owned buffer;
- potential preparation is reusable across frames;
- per-frame analysis does not require creating an owned model;
- dense coordinate lookup is constant-time;
- semantic-to-dense lookup is bounded and suitable for large biomolecular
  systems.

Add targeted benchmarks where useful, but code structure and tests should make
the intended complexity clear even when benchmark corpora are unavailable.

## Error-model requirements

Create focused non-exhaustive public errors rather than one catch-all error:

```text
TopologyBuildError
TopologyError
TopologyMappingError
PositionError
ConfigurationError
ModelError
EnsembleError
TrajectoryError
FrameError
SelectionError
```

Errors should carry semantic IDs when available and distinguish:

- invalid local identifier;
- invalid definition or instance identifier;
- topology identity mismatch;
- count mismatch;
- non-finite values;
- incompatible units;
- invalid periodic cell;
- unsupported format field;
- missing required trajectory topology;
- unsupported random access;
- capacity overflow;
- invalid hierarchy relation;
- inconsistent ensemble member topology.

Do not convert structured internal errors into opaque strings at public
boundaries except for explicitly backend-owned failures.

## Documentation requirements

At each stage:

- update rustdoc for changed public types;
- update affected feature documents;
- include examples of the intended ownership model;
- avoid claiming unimplemented trajectory codecs;
- distinguish ensemble order from trajectory time;
- distinguish molecule instance from connected component;
- distinguish topology identity from structural equivalence;
- distinguish topology from prepared mechanical system;
- describe units and dense ordering explicitly.

The final documentation should include examples equivalent to:

```rust
let topology = Topology::builder()
    .add_small_molecule_definition(...)
    .add_instance(...)
    .build()?;

let model_a = Model::new(topology.clone(), configuration_a)?;
let model_b = Model::new(topology.clone(), configuration_b)?;

let mut potential = DreidingPotential::prepare(&topology, ...)?;
let e_a = potential.evaluate(model_a.view())?;
let e_b = potential.evaluate(model_b.view())?;
```

Each potential documents which configuration state it supports. In 0.2.0 the
built-in harmonic and DREIDING potentials require nonperiodic views and return
structured errors when a periodic cell is present.

The streaming form is:

```rust
let mut reader = open_trajectory(topology.clone(), source)?;
let mut frame = FrameBuffer::new(topology.clone());

while reader.read_next(&mut frame)? {
    let result = potential.evaluate(frame.model_view())?;
    // analyze without creating an owned Model
}
```

Exact API names may evolve during implementation, but examples must express the
same ownership and compatibility rules.

## Review checkpoints

Request focused review after these checkpoints:

1. `Topology` identity, IDs, and dense ordering.
2. Molecule-definition reuse and builder transactionality.
3. `Positions`/`Configuration`/`Model`.
4. Static hierarchy versus observation data.
5. Borrowed view and prepared-potential migration.
6. Ensemble and trajectory contracts.
7. mmCIF migration.
8. Public cleanup and release bump.

At each checkpoint, explicitly ask reviewers to inspect invariants and ownership
rather than only compilation.

## Final acceptance checklist

The refactor is complete only when all statements below are true:

- [x] `Topology` is public and independent of coordinates.
- [x] `Topology` is cheap-clone and immutable.
- [x] Molecule definitions and instances are distinct.
- [x] Explicit definition reuse works.
- [x] Dense atom and bond indices belong to topology.
- [x] Exact identity and structural equivalence are distinct.
- [x] Topology construction is linear and transactional.
- [x] `Model` is topology plus one configuration.
- [x] `Ensemble` shares one topology across members.
- [x] `Trajectory` shares one topology across ordered frames.
- [x] Streaming frame reads reuse a buffer.
- [x] Periodic cells and dynamic arrays are not topology state.
- [x] Static SMCRA hierarchy is coordinate-independent.
- [x] Observation-specific mmCIF data is outside topology.
- [x] DSSP and potentials operate on borrowed structural views.
- [x] Prepared potentials bind to topology identity.
- [x] DREIDING evaluates multiple configurations without re-preparation.
- [x] Molecule instances are not conflated with connected components.
- [x] Topology edits use explicit mappings.
- [x] mmCIF single-model interpretation remains explicit.
- [x] Multi-model mmCIF has a separate ensemble path.
- [x] No obsolete public model-topology types remain.
- [x] Workspace version has the required minor bump.
- [x] Feature contracts and generated dashboard agree.
- [x] Full workspace checks pass.
- [x] MolStudio validation is reported.
- [x] Every command not run is listed with a reason.
