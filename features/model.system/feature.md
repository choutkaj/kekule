# Fixed-Topology Molecular Model

## Summary

Represent one or more small-molecule components as a modelling snapshot with one fixed global topology and exactly one complete coordinate set.

## Behavior/API

- Exposes `MolecularModel`, `ModelDefinitionKey`, `MolecularModelBuilder`, `Component`, `ComponentId`, and `ComponentMapping` under `molecules::modeling`.
- Constructs a single-component model with `MolecularModel::from_conformer` or multiple components through the builder.
- Returns source-to-model atom and bond mappings for each added component.
- Exposes read-only topology and component membership plus validated position access and replacement.
- Provides an opaque immutable-definition key preserved by clones and coordinate updates.
- Rejects empty models/components, invalid conformers, missing coordinates, and non-finite coordinates transactionally.

## Implementation Notes

- Each source `SmallMolecule` is cloned into one fresh, contiguous global `AtomId`/`BondId` space.
- Atom and bond payloads, properties, stored stereo elements/groups, and source bond marks are remapped and preserved.
- Source molecule properties live on the component. Other conformers, derived stereo descriptors, and perception caches are not copied.
- The built topology, components, and atom-membership index form one immutable shared definition; clones share it while owning independent coordinate vectors.
- Independently built models have distinct definition keys even when structurally equal, and atom-to-component lookup is constant time.
- Construction never sanitizes or prepares the source molecule implicitly.

## Validation

- Unit tests cover single- and multi-component construction, exact conformer selection, source tombstones and mappings, source immutability, payload/property/stereo preservation, coordinate mutation, and transactional errors.
- Downstream-style integration tests compile the public modelling namespace.
- Reference molecular goldens are not required for this data-structure feature; `validated` remains false until accepted harness evidence exists.

## Out Of Scope

- `MacroMolecule`, mutable topology, component roles, model merging/transforms, periodic cells, constraints, velocities, serialization, and preparation.

## Revision Notes

- v2: Share immutable model definitions across clones and expose opaque definition identity for prepared-potential binding.
- v1: Add the SmallMolecule-only fixed-topology model and transactional component builder.
