# Topology Construction and Shared Ownership

## Summary

Construct one immutable coordinate-free molecular system from connected,
reusable molecule definitions and explicitly identified instances.

## Behavior/API

- The public `topology` module owns molecule definitions, molecule
  instances, instance-qualified atom and bond identities, and authoritative
  dense atom and bond orderings.
- `Topology` directly owns its private definition, instance, dense-order, and
  index-map collections and deliberately does not implement `Clone`.
  Topology-bound containers retain `Arc<Topology>` values; private
  `Arc::ptr_eq` checks accept clones of one shared allocation and reject
  independently constructed systems.
- `Topology::same_layout` compares complete chemical/static content,
  definition and instance partitioning, semantic IDs, authoritative dense
  order, and index maps. It does not perform order-independent structural
  equivalence or graph isomorphism.
- `TopologyBuilder` is independently usable without coordinates, supports
  explicit definition reuse, uses checked fixed-width identifiers, and commits
  additions only after validation.
- Every small or macro molecule definition must contain one connected graph.
  `TopologyBuildError::DisconnectedMoleculeDefinition` rejects an invalid
  definition transactionally before it can become part of a topology.
- Capacity failures identify definitions, instances, dense atoms, or dense
  bonds through `TopologyIdKind`; shared checked conversion is used for both
  insertion and reservation limits.
- Macro definitions receive static graph/hierarchy validation, including
  connectedness, with coordinate checking disabled. Topology construction
  neither scans nor validates source conformers and strips them only from the
  stored definition clone.
- `TopologyMapping::between_identical_layouts` produces identity maps only
  when complete layouts already match. Explicit mappings validate injectivity,
  definition/instance/atom/bond relationships, mapped bond endpoints, and both
  source and target shared topology allocations.
- Topology-changing operations return a new topology plus checked
  `TopologyMapping` lineage rather than mutating existing coordinate
  containers. Added and removed definitions, instances, atoms, and bonds are
  derived from the validated maps.
- Validated mappings retain the exact source/target `Arc<Topology>` values and stable
  source-order traversal across definition, instance, atom, bond, and dense
  index pairs for explicit remapping kernels.
- `topology::transform::{retain_instances, remove_instances}` creates
  deterministic immutable subsets of complete molecule instances while
  preserving reused definitions and the original source `Arc` for no-op edits.
- Compiled atom selections retain an `Arc<Topology>` and reusable dense
  indices while chemical query syntax remains a separate concern.

## Implementation Notes

- Molecule definitions retain coordinate-independent graph, perception, and
  hierarchy state but do not copy local conformers into topology.
- Dense ordering follows instance insertion order and live local atom or bond
  order and is immutable for the topology lifetime.
- Builder transactionality stages only the new addition and never clones the
  accumulated builder.
- Static macro validation reports zero conformers and coordinates checked, so
  topology construction work is independent of unused source conformer count.

## Tests

- Tests cover shared-allocation compatibility versus equal layout, direct
  topology data ownership, absence of the former public identity API, differing definition,
  instance, and dense insertion order, repeated identical definitions and
  instances, explicit definition reuse, qualified identifiers, dense inverse
  mappings, tombstones, disconnected-definition rejection, checked capacity
  failures, transactional construction, selections, mapping target mismatch,
  duplicate and cross-instance atom maps, mapped-bond endpoint consistency,
  and added/removed reporting at every topology level.
- Macro regressions cover successful coordinate-free construction with many
  unused invalid conformers and verify that stored definitions contain no
  conformers while sources remain unchanged.
- Canonical reconstruction regressions build independent small and macro
  definitions after lossless perception/hierarchy/stereo restoration and
  require `same_layout`.
- Whole-instance transformation regressions cover invalid and duplicate
  requests, empty and no-op results, reused definitions, filtered source order,
  tombstoned local identifiers, roles, properties, hierarchy, and complete
  lineage.

## Out Of Scope

- Automatic definition interning, dynamic coordinates, force-field state,
  reactive trajectories, implicit coordinate transfer, and a structural
  selection language.
- Order-independent structural equivalence, graph-isomorphism mapping,
  ambiguity resolution for indistinguishable definitions or instances, and
  automatic reconciliation of differing dense layouts. These remain planned
  future capabilities.

## Revision Notes

- v1: Track the topology-centered public contract for the initial release.
- v2: Implement independent immutable topology, explicit reusable definitions,
  exact compatibility, dense mappings, selections, and lineage mappings.
- v3: Replace the misleading structural-equivalence name with exact
  `same_layout` semantics and harden mappings across definitions, instances,
  atoms, bond endpoints, added/removed reporting, and edit-result target
  compatibility.
- v4: Make macro topology construction strictly coordinate-independent through
  static-only graph/hierarchy validation that performs no unused-conformer work.
- v5: Replace string-labelled capacity failures with
  `IdentifierCapacityExceeded(TopologyIdKind)` and add synthetic boundary
  regressions for every topology identifier space.
- v6: Add stable checked mapping traversal and deterministic immutable
  whole-instance retain/remove transformations.
- v7: Require independently reconstructed complete molecule definitions,
  including perception, hierarchy enrichment, and stable tombstones, to retain
  exact `same_layout` content.
- v8: Require every topology definition graph to be connected, add structured
  transactional rejection for invalid definitions, and preserve connectedness
  through whole-instance transformations.
- v9: Make `Topology` directly own its data, remove raw cloning and public
  identity machinery, and use retained `Arc<Topology>` values for exact
  compatibility, mappings, selections, and no-op edits.
