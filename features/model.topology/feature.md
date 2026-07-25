# Topology Construction and Identity

## Summary

Construct one immutable coordinate-free molecular system from reusable molecule
definitions and explicitly identified instances.

## Behavior/API

- The planned public `topology` module owns molecule definitions, molecule
  instances, instance-qualified atom and bond identities, and authoritative
  dense atom and bond orderings.
- Exact topology identity is distinct from structural equivalence. Clones keep
  exact identity; independently constructed equal systems do not.
- `TopologyBuilder` is independently usable without coordinates, supports
  explicit definition reuse, uses checked fixed-width identifiers, and commits
  additions only after validation.
- Topology-changing operations return a new topology plus `TopologyMapping`
  lineage rather than mutating existing coordinate containers.
- Compiled atom selections bind to exact topology identity and reusable dense
  indices while chemical query syntax remains a separate concern.

## Implementation Notes

- Molecule definitions retain coordinate-independent graph, perception, and
  hierarchy state but do not copy local conformers into topology.
- Dense ordering follows instance insertion order and live local atom or bond
  order unless the final implementation documents an equivalent immutable rule.
- Builder transactionality stages only the new addition and never clones the
  accumulated builder.

## Tests

- Planned tests cover identity, structural equivalence, explicit definition
  reuse, qualified identifiers, dense inverse mappings, tombstones, checked
  capacity failures, transactional construction, selections, and mappings.

## Out Of Scope

- Automatic definition interning, dynamic coordinates, force-field state,
  reactive trajectories, implicit coordinate transfer, and a structural
  selection language.

## Revision Notes

- v1: Track the topology-centered public contract for the 0.2.0 transition.
