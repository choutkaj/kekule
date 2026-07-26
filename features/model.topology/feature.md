# Topology Construction and Identity

## Summary

Construct one immutable coordinate-free molecular system from reusable molecule
definitions and explicitly identified instances.

## Behavior/API

- The public `topology` module owns molecule definitions, molecule
  instances, instance-qualified atom and bond identities, and authoritative
  dense atom and bond orderings.
- Exact topology identity controls compatibility. Clones keep exact identity;
  independently constructed systems do not.
- `Topology::same_layout` compares complete chemical/static content,
  definition and instance partitioning, semantic IDs, authoritative dense
  order, and index maps. It does not perform order-independent structural
  equivalence or graph isomorphism.
- `TopologyBuilder` is independently usable without coordinates, supports
  explicit definition reuse, uses checked fixed-width identifiers, and commits
  additions only after validation.
- Macro definitions receive static graph/hierarchy validation with coordinate
  checking disabled. Topology construction neither scans nor validates source
  conformers and strips them only from the stored definition clone.
- `TopologyMapping::between_identical_layouts` produces identity maps only
  when complete layouts already match. Explicit mappings validate injectivity,
  definition/instance/atom/bond relationships, mapped bond endpoints, and both
  source and target topology identities.
- Topology-changing operations return a new topology plus checked
  `TopologyMapping` lineage rather than mutating existing coordinate
  containers. Added and removed definitions, instances, atoms, and bonds are
  derived from the validated maps.
- Compiled atom selections bind to exact topology identity and reusable dense
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

- Tests cover exact identity versus equal layout, differing definition,
  instance, and dense insertion order, repeated identical definitions and
  instances, explicit definition reuse, qualified identifiers, dense inverse
  mappings, tombstones, checked capacity failures, transactional construction,
  selections, mapping target mismatch, duplicate and cross-instance atom maps,
  mapped-bond endpoint consistency, and added/removed reporting at every
  topology level.
- Macro regressions cover successful coordinate-free construction with many
  unused invalid conformers and verify that stored definitions contain no
  conformers while sources remain unchanged.

## Out Of Scope

- Automatic definition interning, dynamic coordinates, force-field state,
  reactive trajectories, implicit coordinate transfer, and a structural
  selection language.
- Order-independent structural equivalence, graph-isomorphism mapping,
  ambiguity resolution for indistinguishable definitions or instances, and
  automatic reconciliation of differing dense layouts. These remain planned
  future capabilities.

## Revision Notes

- v1: Track the topology-centered public contract for the 0.2.0 transition.
- v2: Implement independent immutable topology, explicit reusable definitions,
  exact identity, dense mappings, selections, and lineage mappings.
- v3: Replace the misleading structural-equivalence name with exact
  `same_layout` semantics and harden mappings across definitions, instances,
  atoms, bond endpoints, added/removed reporting, and edit-result target
  identity.
- v4: Make macro topology construction strictly coordinate-independent through
  static-only graph/hierarchy validation that performs no unused-conformer work.
