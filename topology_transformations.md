# Topology transformations and state-remapping plan

## Purpose

This document defines the next implementation milestone after the topology-centered
0.2.0 refactor.

The goal is to make immutable `Topology` useful for real system editing without
weakening exact topology identity or silently losing topology-bound state. The
first vertical slice is deliberately narrow:

```text
Topology / Model / Ensemble / Trajectory
                |
                | retain or remove complete molecule instances
                v
new Topology + explicit TopologyMapping
                |
                v
explicitly remapped topology-bound state
```

This milestone must support common operations such as stripping solvent and ions,
extracting a protein-ligand complex, or retaining selected molecule instances. It
must also prove that coordinates, observations, selections, ensemble members, and
trajectory frames can follow an immutable topology edit safely.

The architecture in [`ARCHITECTURE.md`](ARCHITECTURE.md) remains authoritative.
This plan refines the already documented rule that topology-changing operations
return a new topology and explicit lineage mappings, while coordinate transfer is
a separate explicit operation.

## Primary feature and affected contracts

Implementation should start by creating a canonical feature contract:

```text
model.topology-transform
```

Suggested description:

> Create immutable molecule-instance subsets with explicit topology lineage and
> remap compatible topology-bound structure state.

Directly affected existing feature contracts are expected to include:

- `model.topology`;
- `model.system`;
- `model.ensemble`;
- `model.trajectory`;
- `api.public-facade`.

List every actually affected feature ID in the PR. Do not change unrelated
feature versions or generated dashboard output by hand.

## Development strategy

Implement the milestone on one short-lived branch, preferably:

```text
codex/model.topology-transform-instance-subset
```

Open a draft PR early. Keep commits aligned with the stages below and leave the
workspace buildable after each logical commit. If the implementation becomes too
large to review, split at stage boundaries into stacked PRs rather than weakening
tests or transactionality.

## Required outcome

At completion, the repository must provide:

1. Public immutable topology operations for retaining or removing complete
   molecule instances.
2. Deterministic construction of the target topology with explicit reuse of
   retained molecule definitions.
3. A complete checked `TopologyMapping` from the source topology to the target
   topology.
4. Explicit remapping of `Positions`, `Configuration`, `StructureObservation`,
   `Model`, and `AtomSelection`.
5. Explicit remapping of finite `Ensemble` values, owned `TrajectoryFrame`
   values, in-memory `Trajectory`, and reusable `FrameBuffer` destinations.
6. Exact source and target topology-identity validation at every remapping
   boundary.
7. Preservation of all retained dynamic state: positions, periodic cell,
   observations, member weights and properties, velocities, forces, time, step,
   frame properties, and member/frame order.
8. Structured errors for invalid IDs, empty target topologies, incomplete
   mappings, removed selected atoms under strict policy, and unsupported added
   atoms.
9. Transactional behavior: failed transformations or remaps leave source objects
   and caller-owned destinations unchanged.
10. Linear-time behavior in the amount of source or target data being processed,
    without graph isomorphism or repeated cloning of reused definitions.
11. Public API tests, focused regression tests, feature-contract updates, and
    documentation consistent with the final implementation.

## Architectural rules

### Immutability and lineage

A built `Topology` remains immutable. A nontrivial transformation returns:

```rust
TopologyEditResult {
    topology: Topology,
    mapping: TopologyMapping,
}
```

The source topology is never modified. Existing models, ensembles, trajectories,
selections, frame buffers, and prepared potentials remain bound to the source
topology.

`TopologyMapping` in this milestone represents known edit lineage. It must not
perform graph isomorphism, infer correspondence between independently constructed
topologies, or silently choose among ambiguous mappings.

### Exact identity

Every remapping operation validates both ends:

- the mapping source identity must equal the source object's topology identity;
- the mapping target identity must equal the supplied target topology identity.

Equal atom counts or `Topology::same_layout` are not sufficient.

A no-op transformation should return a clone of the original topology, preserving
exact identity, together with a complete identity mapping. Examples are
`remove_instances` with an empty set and `retain_instances` containing every
source instance.

### Deterministic subset ordering

Instance subset inputs are membership requests, not target-order requests.

For a nontrivial subset:

1. validate every requested `MoleculeInstanceId` before construction;
2. ignore duplicate requests;
3. retain source molecule-definition order, filtered to definitions referenced
   by retained instances;
4. retain source molecule-instance order, filtered to retained instances;
5. preserve local `AtomId` and `BondId` values inside retained definitions;
6. rebuild target `MoleculeDefinitionId`, `MoleculeInstanceId`,
   `TopologyAtomIndex`, and `TopologyBondIndex` values deterministically;
7. let the returned mapping record every changed semantic and dense identifier.

A retained definition is stored once even when many retained instances reference
it. Definitions with no retained instances are omitted.

### Complete state, never sentinel state

`Positions`, velocity arrays, force arrays, and observation arrays remain complete
for their target topology. Do not use NaN, zero coordinates, sparse vectors, or
placeholder atoms to represent missing remapped state.

The initial topology transforms are deletion-only and therefore introduce no new
atoms. General remapping utilities should nevertheless detect target atoms that
lack source mappings and return a structured error such as
`AddedAtomsRequireState` or `IncompleteTargetMapping`.

Do not add automatic coordinate generation in this milestone.

### Reused definitions

Whole-instance subsetting does not edit definition content. It only retains or
drops complete instances and the definitions they reference.

Future definition-changing transformations must explicitly distinguish editing
all instances of a reused definition from forking one instance onto a new
definition. Do not introduce an implicit shared-definition mutation path while
implementing this milestone.

### Prepared systems

Prepared potentials and backend systems remain bound to exact source topology
identity. They are not remapped. After a topology edit, callers prepare a new
downstream system for the target topology.

## Target public API shape

Exact names may be improved during implementation, but the public semantics must
remain equivalent.

### Topology subset operations

Prefer a focused namespace such as:

```rust
molecular::topology::transform
```

with operations conceptually equivalent to:

```rust
pub fn retain_instances(
    topology: &Topology,
    instances: impl IntoIterator<Item = MoleculeInstanceId>,
) -> Result<TopologyEditResult, TopologyTransformError>;

pub fn remove_instances(
    topology: &Topology,
    instances: impl IntoIterator<Item = MoleculeInstanceId>,
) -> Result<TopologyEditResult, TopologyTransformError>;
```

Required behavior:

- invalid instance IDs are rejected before target construction;
- retaining no instances or removing all instances returns a focused empty-result
  error;
- duplicates do not duplicate target instances;
- no-op requests preserve source identity;
- nontrivial requests produce a new exact topology identity;
- mappings include definitions, instances, atoms, bonds, atom indices, bond
  indices, and correct added/removed complements.

`remove_instances` may normalize to `retain_instances` internally after validating
the removal set.

### Mapping traversal

Add only the mapping access needed by remapping kernels. Suitable operations
include checked source/target identity queries and stable iterators over validated
pairs:

```rust
mapping.definition_pairs()
mapping.instance_pairs()
mapping.atom_pairs()
mapping.bond_pairs()
mapping.atom_index_pairs()
mapping.bond_index_pairs()
mapping.is_source(topology)
mapping.is_target(topology)
```

Do not expose mutable mapping internals. Pair iteration must follow a documented
stable order, preferably source topology order.

### Structure remapping

Provide explicit checked operations conceptually equivalent to:

```rust
Positions::remap_to(...)
Configuration::remap_to(...)
StructureObservation::remap_to(...)
Model::remap_to(...)
```

The implementation may use shared private kernels rather than duplicating
validation.

For deletion-only remapping:

- allocate or reuse storage sized for the target topology;
- copy source values through dense-index mappings;
- prove every target atom is filled exactly once;
- preserve the periodic cell;
- preserve model-level observation metadata and properties;
- preserve source units through the canonical-unit boundary;
- return target-bound values with exact target identity.

Any mapping with target atoms that lack source state must fail explicitly in this
milestone.

### Selection remapping

Provide:

```rust
pub enum RemovedSelectionPolicy {
    Error,
    Drop,
}
```

or an equivalent explicit policy.

`AtomSelection` remapping must:

- validate mapping and topology identities;
- map selected source atoms or dense indices to target atoms;
- error on removed selected atoms under strict policy;
- drop only those atoms under explicit drop policy;
- restore the target selection invariant of sorted unique dense indices;
- bind the result to exact target topology identity.

No selection may silently lose atoms by default.

### Ensemble remapping

Remap every member against one target topology while preserving:

- member order;
- member configuration and periodic cell;
- member observation;
- weight;
- member properties.

Failure in any member leaves the source ensemble unchanged and returns the member
index in the error context. Weight normalization is not rerun.

### Trajectory remapping

Support owned frame and in-memory trajectory remapping while preserving:

- frame order;
- positions;
- periodic cell;
- optional velocities;
- optional forces;
- time;
- step;
- observation;
- frame properties.

Provide a streaming-friendly operation that remaps a borrowed source frame into a
caller-owned target `FrameBuffer`, conceptually:

```rust
destination.copy_remapped_from(
    source_frame,
    source_topology,
    mapping,
)?;
```

The destination must already be bound to the exact target topology. Validate the
entire operation before changing destination-visible state. Reuse destination
allocations and avoid constructing an owned `Model` per frame.

A generic `MappedTrajectoryReader` wrapper is not required in this milestone.

## Error model

Use focused `#[non_exhaustive]` public errors. Avoid string-only failures and do
not collapse all failures into `TopologyBuildError`.

Expected categories include:

```text
TopologyTransformError
- invalid source instance
- empty target topology
- topology construction failure
- mapping construction failure

TopologyRemapError
- source identity mismatch
- target identity mismatch
- incomplete target atom mapping
- duplicate target assignment
- added atoms require explicit state
- removed selected atom under strict policy
- incompatible destination buffer
- nested position, observation, frame, or unit error
- member/frame index context where applicable
```

Exact enum division may follow module ownership. Errors should expose typed IDs or
indices where useful.

## Staged implementation

## Stage 0: Contracts and baseline

### Work

- Read `AGENTS.md`, `ARCHITECTURE.md`, this plan, and all affected feature
  contracts.
- Create `features/model.topology-transform/feature.toml` and `feature.md` with
  initial planned or experimental status.
- Record baseline repository checks.
- Add public compile tests for the intended subset and remap API before or
  alongside implementation.
- Inventory all current `TopologyMapping`, `TopologyEditResult`, `AtomSelection`,
  structure, ensemble, and trajectory call sites.
- Do not modify `README.md` without separate owner consent.

### Exit criteria

- Scope and affected feature IDs are explicit.
- No feature contract claims support before the implementation and tests exist.
- The final public dependency direction is understood before code movement.

## Stage 1: Mapping and transformation foundations

### Work

- Add stable read-only traversal of validated mapping pairs as required by
  remapping kernels.
- Add source/target topology identity helpers if they improve error handling.
- Introduce focused topology transformation errors.
- Implement internal helpers for deterministic filtered definition and instance
  construction.
- Preserve no-op exact identity.
- Keep mapping validation centralized in `TopologyMapping::from_pairs` or an
  equivalent checked constructor.

### Required tests

- stable pair traversal;
- source/target identity helpers;
- no mutable mapping internals exposed;
- complete identity mapping for a no-op edit;
- invalid mapping relationships remain rejected.

### Exit criteria

- Later stages can consume mappings without reaching into private fields or
  duplicating mapping validation.

## Stage 2: Whole-instance topology subsets

### Work

Implement `retain_instances` and `remove_instances`.

Build the target in one deterministic pass:

1. validate and normalize the membership set;
2. determine referenced definitions;
3. add retained definitions once in filtered source order;
4. add retained instances in filtered source order with cloned metadata;
5. build the target topology;
6. construct complete definition, instance, atom, and bond mappings;
7. validate and return `TopologyEditResult`.

Use local atom and bond IDs to map retained definition content. Never infer atom
correspondence from graph matching.

### Required tests

- retain a ligand and protein while dropping waters and ions;
- remove selected instances;
- repeated water definition remains stored once;
- partial retention of instances sharing one definition;
- roles and instance properties are preserved;
- macro hierarchy content is preserved;
- source definition and instance ordering rules;
- local atom/bond tombstones and IDs remain correct;
- dense atom/bond mappings are complete and inverse-consistent;
- invalid IDs fail before construction;
- duplicate IDs are harmless;
- retain-none and remove-all fail;
- retain-all and remove-none preserve exact identity;
- source topology remains unchanged;
- target mapping added/removed complements are exact.

### Performance expectations

- time is linear in visited definitions, instances, atoms, and bonds;
- retained repeated definitions are cloned once, not once per instance;
- no clone of the accumulated builder occurs per addition;
- a large synthetic solvent-rich topology demonstrates practical scaling.

### Exit criteria

- common system stripping and extraction operations produce a correct immutable
  target topology and complete lineage.

## Stage 3: Base structure and selection remapping

### Work

- Add one shared checked atom-array remapping kernel.
- Remap `Positions`.
- Remap `Configuration` while preserving the cell.
- Remap `StructureObservation` while preserving source model ID and properties.
- Remap `Model`.
- Remap `AtomSelection` under explicit removed-atom policy.
- Reject target atoms without source state.

### Required tests

- nontrivial dense-index reordering;
- exact source and target identity mismatch;
- complete coordinate transfer;
- cell preservation;
- occupancy, B-factor, altloc, raw coordinate text, source model ID, and
  observation properties;
- strict and drop selection policies;
- empty target selection when drop policy removes every selected atom;
- synthetic mapping with added atoms is rejected;
- remap failure leaves every source object unchanged.

### Exit criteria

- one model and its compiled selections can follow a deletion-only topology edit
  without silent data loss.

## Stage 4: Ensemble remapping

### Work

- Remap all ensemble members through the base structure kernel.
- Preserve topology sharing: the target topology is stored once.
- Preserve member order, weights, observations, and properties.
- Add member-index context to failures.

### Required tests

- multiple configurations map to one exact target topology;
- weights and properties are byte-for-byte or value-for-value preserved;
- variable member cells remain distinct;
- one failing member aborts the operation;
- source ensemble remains unchanged;
- target ensemble contains no repeated topology handles per member.

### Exit criteria

- conformer sets, NMR ensembles, and docking poses can be subset consistently.

## Stage 5: Trajectory and reusable-buffer remapping

### Work

- Remap owned `TrajectoryFrame`.
- Remap in-memory `Trajectory`.
- Remap borrowed frame state into an existing target `FrameBuffer`.
- Reuse buffer allocations in repeated frame remaps.
- Preserve optional arrays and all frame metadata.
- Keep fixed-topology trajectory invariants intact.

### Required tests

- positions-only frames;
- velocities and forces;
- variable periodic cells;
- time and step;
- observations and frame properties;
- stable frame order;
- target buffer identity mismatch;
- repeated remapping does not grow position/vector capacities after warm-up;
- transactional destination behavior on validation failure;
- direct potential or analysis use of the remapped target frame view after
  preparing for the target topology.

### Exit criteria

- finite trajectories and streamed frame-by-frame workflows can follow a subset
  edit without constructing an owned model per frame.

## Stage 6: Public cleanup, documentation, and validation

### Work

- Promote `model.topology-transform` to `supported` only after all promised
  behavior and tests exist.
- Update directly affected feature contracts and versions.
- Update `ARCHITECTURE.md` only where the implemented public behavior needs a
  concise availability statement; do not rewrite its object model.
- Update `ROADMAP.md` to mark this milestone complete and identify the next
  topology-editing direction.
- Add rustdoc examples for stripping solvent and remapping a model.
- Keep the prelude focused and avoid broad root re-exports.
- Do not modify `README.md` without separate owner consent.
- Regenerate procedural dashboard/skills output through `xtask`.

### Required validation

Run and report:

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
cargo check --manifest-path fuzz/Cargo.toml --all-targets --locked
```

Run Molecular checks before any optional MolStudio consumer check. This milestone
is additive and does not require MolStudio edits unless an actual public consumer
break is discovered.

### Exit criteria

- no unsupported feature claim remains;
- all topology-bound state promised by this plan is remapped explicitly;
- docs and generated metadata describe the implemented behavior accurately;
- all applicable checks pass or unavailable commands are reported precisely.

## Non-goals

The following are explicitly outside this implementation milestone:

- atom-level or bond-level topology subsetting within a molecule instance;
- splitting or merging molecule instances;
- appending atoms or molecule instances to an existing topology;
- replacing one instance's definition;
- editing a reused definition or forking one instance from it;
- topology-level hydrogen addition/removal;
- residue mutation or chemical reactions;
- coordinate generation for added atoms;
- graph-isomorphism-based topology equivalence;
- inferred correspondence between independently constructed topologies;
- ambiguity resolution for repeated identical definitions or instances;
- automatic remapping of prepared force fields or backend systems;
- reactive trajectories;
- a generic mapped trajectory-reader adapter;
- production DCD, XTC, TRR, NetCDF, or other trajectory codecs;
- a structural selection language;
- unrelated performance or API refactors.

Do not broaden the PR into these areas merely because the new mapping primitives
make them tempting.

## Follow-on topology work

After this milestone is stable, proceed in separate feature branches:

1. topology composition: append definitions and instances with complete lineage;
2. explicit instance-definition replacement;
3. definition-edit scope:
   `AllInstances` versus `ForkInstance(MoleculeInstanceId)`;
4. lifting local hydrogen and chemical transforms into topology-aware operations;
5. atom-level extraction with hierarchy-safe molecule splitting;
6. ambiguity-aware structural correspondence and graph-isomorphism mappings;
7. topology-aware MolStudio editing workflows.

Each follow-on must retain the same identity, transactionality, and explicit-state
principles established here.

## Final acceptance checklist

The milestone is complete only when all statements below are true:

- whole molecule instances can be retained or removed from any valid topology;
- target ordering and identifier remapping are deterministic;
- reused definitions remain reused;
- the source topology and source state are unchanged;
- the returned mapping is complete and validated at every topology layer;
- no topology-bound state is accepted on equal count or equal layout alone;
- positions, configurations, observations, models, selections, ensembles, frames,
  trajectories, and reusable buffers can be remapped explicitly;
- removed selection members are never dropped without explicit policy;
- target atoms without source state are never filled with sentinels;
- prepared systems are not silently reused;
- repeated frame remapping can reuse allocations;
- public errors are structured and typed;
- feature metadata, tests, rustdoc, roadmap, and architecture statements agree;
- all applicable repository checks pass.
