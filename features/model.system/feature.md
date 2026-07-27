# Topology-Bound Molecular Model

## Summary

Represent one concrete molecular structure as an immutable `Topology` plus one
complete mutable `Configuration` and optional observation state.

## Behavior/API

- Exposes `Positions`, `Configuration`, `ConfigurationView`, `Model`,
  `ModelView`, and `ModelBuilder` from `structure`; topology IDs and system
  structure live in the separate `topology` module.
- Builder insertion uses `add_small_molecule[_with_metadata]` and
  `add_macro_molecule[_with_metadata]` and returns a stable instance ID.
- Preserves molecule-local atom and bond IDs, including tombstones; qualification
  adds instance ownership and topology dense indices round-trip to qualified IDs.
- Copies one complete finite source conformer into authoritative model positions
  without copying other source conformers or mutating sources.
- Macro model insertion validates static graph/hierarchy state through topology
  staging and validates only the explicitly selected conformer while staging
  positions. Unselected conformers cannot make an otherwise valid selection
  fail.
- Converts source coordinates once to `MODEL_LENGTH_UNIT`; model position
  getters and setters expose explicit quantities and accept compatible length units.
- `Positions` bind to exact topology identity, reject incomplete or non-finite
  arrays, and reuse their allocation for validated full-coordinate updates.
- `Configuration` owns positions and an optional validated periodic cell.
- `ModelView` borrows topology plus configuration without copying coordinates.
- `StructureObservation` stores topology-bound coordinate-model-specific atom
  values outside topology.
- `Positions`, `Configuration`, `StructureObservation`, and `Model` remap
  explicitly through checked topology lineage. Exact identities are required
  at both ends, complete target atom state is mandatory, cells and observation
  metadata are preserved, and failures leave sources unchanged.
- `Model::instance_to_conformer` maps current instance positions back through
  preserved local atom IDs, converts them to the target conformer unit, and
  commits only after the target live-atom set and all conversions validate.
- Rejects empty topologies/molecules, invalid conformers, missing positions, and
  non-finite positions transactionally.
- Validates every `MacroMolecule` graph/hierarchy pair before accepting it as a
  model instance.

## Implementation Notes

- `Model = Topology + Configuration`; cloning a model shares topology identity
  while copying dynamic configuration and observation state.
- Complete positions and cells may change without changing topology identity.
- Conformer export belongs to the modeling layer, validates exact live local
  atom-ID compatibility, and does not mutate topology or chemistry.
- Construction never sanitizes, perceives, prepares, or merges source molecules.

## Tests

- Unit tests cover independent topology/configuration construction, exact
  identity rejection, shared topology after cloning, unit conversion,
  allocation reuse, periodic cells, source immutability, conformer export, and
  transactional failures.
- Macro construction tests cover one valid selected conformer alongside many
  unrelated invalid conformers, rejection when an invalid conformer is selected,
  and preservation of all source conformers.
- Public transformation tests cover dense-index compaction, equal-layout
  identity rejection, complete coordinate transfer, cells, every observation
  field and property, missing target state, and source immutability.

## Out Of Scope

- Topology mutation, reactions, constraints, virtual sites, Drude particles,
  and backend preparation.

## Revision Notes

- v1: SmallMolecule-only flattened component model.
- v2: Hard break to typed molecule instances, qualified IDs, mixed Small/Macro
  ownership, and authoritative positions.
- v3: Add shared opaque definition identity for binding prepared potentials
  without flattening molecule instances.
- v4: Rename the canonical `MolecularModel`/builder API to `Model` and
  `ModelBuilder`, and rename its qualified hierarchy view to
  `InstanceSmcraHierarchy`.
- v5: Replace implicit model coordinate conventions with quantity-valued
  positions and explicit compatible conversion at model boundaries.
- v6: Make valid macromolecule structure a checked model-construction
  precondition.
- v7: Add transactional `Model::instance_to_conformer` coordinate export using
  preserved molecule-local atom IDs.
- v8: Replace model-owned topology and coordinates with the explicit
  `Model = Topology + Configuration` contract, exact topology-bound positions,
  periodic cells, observation state, and borrowed `ModelView`.
- v9: Validate only the explicitly selected macro conformer during model
  staging, leaving unrelated conformers to optional standalone full validation.
- v10: Add explicit transactional remapping of complete positions,
  configurations, observations, and models through exact topology lineage.
