# Molecule Refactor Plan

## Goal

Replace the current molecule/domain-wrapper implementation with the architecture
defined in `ARCHITECTURE.md`:

```rust
pub struct Molecule {
    graph: Graph,
    hierarchy: Hierarchy,
    perception: Perception,
}
```

Every published `Molecule` is non-empty, connected, and geometry-independent.
`SmallMolecule` and `MacroMolecule` cease to be owning foundational types.

This is a molecule-layer refactor. Do not redesign `Topology`, `Model`,
`Ensemble`, or `Trajectory`.

## Scope guard

Allowed changes:

- core molecule/graph/editor/perception code;
- hierarchy code and biological APIs that currently depend on `MacroMolecule`;
- small-molecule APIs that currently depend on `SmallMolecule`;
- parsers/interpreters and writers as required by the new component API;
- tests/examples/docs affected by the new molecule API;
- minimal mechanical updates to `Topology` and downstream consumers required
  because the wrapper types disappear.

Out of scope:

- new topology representation;
- new model/ensemble/trajectory architecture;
- force-field/potential redesign;
- new chemistry perception models;
- broad feature additions unrelated to the refactor.

Prefer deletion and simplification over compatibility shims that preserve the
old wrapper architecture.

## Stage 1 - Introduce the new core ownership model

Create/refactor the core types so `Molecule` owns only:

```text
Graph
Hierarchy
Perception
```

Move current authoritative atom/bond/adjacency/stereo storage into `Graph`.

Requirements:

- stable local `AtomId` and `BondId` semantics remain explicit;
- `Graph` is geometry-independent;
- represented stereo remains authoritative graph state;
- no conformer storage remains in `Molecule`;
- no perception flags are duplicated into `Atom`/`Bond`.

Do not preserve old fields merely to reduce diff size if they violate the new
architecture.

## Stage 2 - Integrate hierarchy

Replace the owning `MacroMolecule` wrapper with `Hierarchy` stored directly in
`Molecule`.

Reuse useful current SMCRA/hierarchy implementation rather than rewriting sound
data structures without reason, but reshape naming and ownership around the new
contract.

Requirements:

- `Hierarchy` may be empty;
- hierarchy references only live `AtomId`s in the same graph;
- hierarchy presence does not classify the molecule into another Rust type;
- hierarchy validation is part of molecule publication;
- a disconnected biological input is split into separate molecules rather than
  forced into one hierarchy-bearing disconnected object.

Remove `MacroMolecule`, `MacroMoleculeBuilder`, `MacroMoleculeEditor`, and
wrapper-only graph-edit APIs once their functionality is represented by
`Molecule`/`MoleculeEditor`.

## Stage 3 - Reframe perception

Rename/rework `PerceptionState` into `Perception`.

Preserve useful existing sectional perception data where compatible:

```text
valence / implicit H
rings
aromaticity
stereo/CIP
```

Requirements:

- perception is derived and reconstructible;
- graph edits cannot leave stale perception published;
- checked exact installation remains possible if currently required by
  persistence/MolStudio;
- perception cache presence does not define molecular identity;
- public molecule-level convenience accessors may expose perceived chemistry.

Avoid a generic cache/dependency framework unless required by existing tests.
Simple robust invalidation is preferred.

## Stage 4 - Make `MoleculeEditor` the mutation boundary

Use one transactional editor for structural construction and editing.

During editing the graph may be:

- empty;
- disconnected;
- temporarily chemically incomplete.

`finish()` publishes only a valid `Molecule`.

At minimum validate:

- non-empty graph;
- exactly one connected component;
- atom/bond/adjacency integrity;
- stereo references;
- hierarchy references and graph consistency.

Remove awkward connectedness-preserving per-operation APIs. Do not require every
intermediate editor state to remain connected.

Do not expose unrestricted mutable graph access that bypasses publication
validation.

A separate builder is optional only if it materially improves ergonomics
without duplicating editor semantics.

## Stage 5 - Remove `SmallMolecule`

Delete the owning `SmallMolecule` wrapper and route generic chemistry workflows
to `Molecule`.

Small-molecule-specific algorithms may remain in logically named modules, but
their input is the universal molecule when no stronger semantic view is needed.

Avoid introducing a replacement wrapper solely to preserve the old API shape.

## Stage 6 - Componentize parsing

For one molecular record, the canonical molecule-producing result is:

```rust
Result<Vec<Molecule>>
```

Requirements:

- one connected input -> one-element vector;
- disconnected input -> one molecule per connected component;
- component order is deterministic and follows source/interpreter order;
- no returned `Molecule` is disconnected;
- no parser directly constructs `Topology` as part of this refactor.

Do not assign semantic "main component" meaning to `vec[0]`; callers may choose
the first component or apply an explicit component-selection policy.

For multi-record formats, preserve record boundaries. Prefer streaming/per-record
results rather than flattening every molecule from an entire file into one
vector.

For coordinate-bearing formats, keep molecular chemistry and positions
separate. Do not reintroduce coordinates into `Molecule` to simplify parsing.

## Stage 7 - Adapt Topology minimally

The current `Topology` imports/accepts `SmallMolecule` and `MacroMolecule`, so a
strict zero-touch topology refactor is impossible.

Make only the minimal changes needed for topology to accept the new universal
`Molecule` and its integrated `Hierarchy`.

Preserve:

- definition/instance semantics;
- instance-qualified atom/bond identity;
- topology immutability;
- dense topology ordering;
- `Model = Topology + Positions`;
- existing ensemble/trajectory behavior.

Do not use this refactor as an opportunity to redesign topology.

## Stage 8 - Update persistence, writers, tests, and examples

Update reconstruction and I/O APIs around the new ownership model.

Important regression cases:

- single atom is a valid molecule;
- empty published molecule is rejected;
- disconnected editor cannot finish;
- editor may disconnect and reconnect before finish;
- `[Na+].[Cl-]` parses to two molecules;
- ordinary connected SMILES parses to one molecule;
- protein/nucleic-acid hierarchy works directly on `Molecule`;
- hierarchy references survive valid edits/remapping or fail transactionally;
- aromaticity/CIP remain derived;
- stale perception is never exposed after graph edits;
- molecule has no geometry/conformers;
- existing Topology/Model/Ensemble/Trajectory behavioral tests still pass.

Remove tests whose only purpose is to preserve `SmallMolecule`/`MacroMolecule`
wrapper behavior.

## Stage 9 - Cleanup

After functionality passes:

- remove dead wrapper modules and exports;
- remove compatibility aliases that no longer serve a real external boundary;
- update crate-level docs and examples;
- ensure terminology consistently uses `Graph`, `Hierarchy`, `Perception`, and
  universal `Molecule`;
- run formatting, clippy, unit tests, integration tests, and doc tests.

## Acceptance contract

The refactor is complete when all of the following are true:

```text
Molecule == Graph + Hierarchy + Perception

Molecule is always:
  non-empty
  connected
  geometry-independent

SmallMolecule:
  removed as owning core/domain wrapper

MacroMolecule:
  removed as owning core/domain wrapper

editing:
  transactional through MoleculeEditor
  temporary invalidity allowed only inside editor
  finish() is the publication gate

parsing one molecular record:
  -> Vec<Molecule>

Topology / Model / Ensemble / Trajectory:
  semantics unchanged apart from minimal adaptation to the new Molecule API
```

Do not optimize for backward compatibility with the discarded molecule
architecture. Optimize for a small, coherent, enforceable core model.
