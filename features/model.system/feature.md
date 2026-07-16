# Fixed-Topology Molecular Model

## Summary

Represent distinct Small or Macro molecule instances as immutable topology plus
one complete mutable coordinate set.

## Behavior/API

- Exposes `MoleculeInstanceId`, `InstanceAtomId`, `InstanceBondId`,
  `ModelAtomIndex`, `MoleculeInstance`, `MoleculeRole`, `ModelTopology`,
  `ModelDefinitionKey`, `Model`, and `ModelBuilder`.
- Builder insertion uses `add_small_molecule[_with_metadata]` and
  `add_macro_molecule[_with_metadata]` and returns a stable instance ID.
- Preserves molecule-local atom and bond IDs, including tombstones; qualification
  adds instance ownership and dense model indices round-trip to qualified IDs.
- Stores typed Small/Macro payloads, multi-valued roles, properties, and a
  qualified read-only `SmcraHierarchy` view.
- Copies one complete finite source conformer into authoritative model positions
  and removes conformers from stored instance topology without mutating sources.
- Converts source coordinates once to `MODEL_LENGTH_UNIT`; model position
  getters and setters expose explicit quantities and accept compatible length units.
- Rejects empty models/molecules, invalid conformers, missing positions, and
  non-finite positions transactionally.
- Validates every `MacroMolecule` graph/hierarchy pair before accepting it as a
  model instance.

## Implementation Notes

- `ModelTopology` owns ordered molecule instances and immutable bidirectional
  `InstanceAtomId`/`ModelAtomIndex` mappings; it never creates a flattened
  multi-entity `Molecule`.
- Built topology and ownership are immutable. Complete positions may change via
  validated setters.
- Clones share one opaque definition key; independently built models receive
  distinct keys even when structurally equal.
- Construction never sanitizes, perceives, prepares, or merges source molecules.

## Validation

- Unit tests cover mixed Small/Macro instances, repeated molecules, stable local
  IDs, tombstones, dense round-trips, qualified hierarchy lookup, roles, source
  immutability, position mutation, and transactional failures.

## Out Of Scope

- Topology mutation, cells, velocities, trajectories, reactions, merging,
  constraints, virtual sites, Drude particles, and backend preparation.

## Revision Notes

- v1: SmallMolecule-only flattened component model.
- v2: Hard break to typed molecule instances, qualified IDs, fixed
  `ModelTopology`, mixed Small/Macro ownership, and authoritative positions.
- v3: Add shared opaque definition identity for binding prepared potentials
  without flattening molecule instances.
- v4: Rename the canonical `MolecularModel`/builder API to `Model` and
  `ModelBuilder`, and rename its qualified hierarchy view to
  `InstanceSmcraHierarchy`.
- v5: Replace implicit model coordinate conventions with quantity-valued
  positions and explicit compatible conversion at model boundaries.
- v6: Make valid macromolecule structure a checked model-construction
  precondition.
