# Topology-Bound Molecular Model

## Summary

Represent one concrete molecular structure as an immutable `Topology` plus
complete mutable positions, an optional periodic cell, and topology-bound
per-atom and per-bond annotations, with each molecule instance backed by a
connected definition graph.

## Behavior/API

- Exposes `Positions`, `AtomData`, `BondData`, `Model`, `ModelView`, and `ModelBuilder`
  from `structure`; topology IDs and system
  structure live in the separate `topology` module.
- Builder insertion uses `add_small_molecule[_with_metadata]` and
  `add_macro_molecule[_with_metadata]` and returns a stable instance ID.
- Model staging delegates molecule-definition validation to `TopologyBuilder`;
  disconnected input is rejected transactionally through the structured
  topology build error before positions or atom data are committed.
- Preserves molecule-local atom and bond IDs, including tombstones; qualification
  adds instance ownership and topology dense indices round-trip to qualified IDs.
- Copies one complete finite source conformer into authoritative model positions
  without copying other source conformers or mutating sources.
- Macro model insertion validates static graph/hierarchy state through topology
  staging and validates only the explicitly selected conformer while staging
  positions. Unselected conformers cannot make an otherwise valid selection
  fail.
- Converts source coordinates once to `MODEL_LENGTH_UNIT`; `Model::positions`
  and `ModelView::positions` borrow the topology-bound `Positions` directly,
  while `Positions::values` and individual position getters expose explicit
  quantities. Setters accept compatible length units.
- `Positions` retain one `Arc<Topology>`, reject independently allocated,
  incomplete, or non-finite
  arrays, and reuse their allocation for validated full-coordinate updates.
- `Model` directly owns `Arc<Topology>`, `Positions`, `Option<PeriodicCell>`,
  `AtomData`, and `BondData`; there is no intermediate configuration object.
- `AtomData` stores optional dense occupancy and B-factor columns. Wholly
  absent columns allocate no per-atom objects; present columns validate exact
  atom count and finite values and support semantic-ID and dense-index access.
  Occupancy is dimensionless; B-factor APIs use length-squared `Quantity`
  values and normalize storage to square angstroms.
- Canonical occupancy, canonical B-factor, and custom scalar properties reuse
  one private dense unit-aware column implementation. Canonical columns remain
  explicit `AtomData` fields with dedicated APIs and never enter the custom
  property namespace.
- `AtomData::new` starts with no columns. Field-specific setters replace or
  clear occupancy and B-factor columns without a positional all-future-fields
  constructor. `atom_count()` reports topology cardinality and `is_empty()`
  means that no supported metadata column contains data.
- `BondData` is the exact-topology-bound per-bond analogue in authoritative
  `TopologyBondIndex` order and initially has no canonical scientific fields.
- `AtomData` and `BondData` accept conservatively named custom properties only
  as unit-aware scalar `f64` columns with `None` for missing values. Columns
  have exact dense length, reject non-finite values and incompatible units,
  preserve their stored unit across compatible replacements, and normalize
  all-missing values to absence. Atom custom properties cannot shadow
  `occupancy` or `b_factor`.
- Canonical typed atom fields remain separate from custom annotations. Custom
  properties are not authoritative topology, chemical perception, or analysis
  ownership; calculated bond `conformational_entropy` may be attached for
  visualization without making `BondData` the calculation owner.
- `ModelView` directly borrows topology, positions, cell, atom data, and bond data without
  copying coordinates or recreating a configuration wrapper.
- `Model` and `ModelView` provide the same common read-only atom, bond,
  instance, and qualified SMCRA navigation as thin forwards to `Topology`.
  Returned `InstanceChain`, `InstanceResidue`, and `InstanceAtomSite` views
  retain qualification through parent/child navigation, borrow their
  definition-owned metadata, and never reconstruct or copy coordinate arrays.
- `Positions`, `AtomData`, `BondData`, and `Model` remap explicitly through checked topology
  lineage. Exact shared source and target allocations are required, complete
  target atom and bond state is mandatory, cells and data columns are preserved, and
  failures leave sources unchanged.
- `Model::instance_to_conformer` maps current instance positions back through
  preserved local atom IDs, converts them to the target conformer unit, and
  commits only after the target live-atom set and all conversions validate.
- Rejects empty topologies/molecules, disconnected definitions, invalid
  conformers, missing positions, and non-finite positions transactionally.
- Validates every `MacroMolecule` graph/hierarchy pair before accepting it as a
  model instance.

## Implementation Notes

- `Model = Arc<Topology> + Positions + Option<PeriodicCell> + AtomData + BondData`;
  cloning a model shares the topology allocation while copying mutable state.
- One private scalar-column implementation stages complete replacements,
  individual mutation, all-missing normalization, unit conversion, and atom
  remapping for canonical and custom `AtomData` columns. Public errors retain
  canonical field context.
- Complete positions and cells may change without changing the shared topology.
- Conformer export belongs to the modeling layer, validates exact live local
  atom-ID compatibility, and does not mutate topology or chemistry.
- Construction never normalizes, perceives, prepares, or merges source molecules.

## Tests

- Unit tests cover direct model construction without atom metadata, independent
  topology/position allocation rejection, occupancy and B-factor lookup and
  mutation, square-angstrom B-factor round trips, compatible length-squared
  conversion, incompatible-unit and non-finite rejection, dense-column
  validation, shared topology after cloning, coordinate unit conversion,
  allocation reuse, periodic cells, source immutability, conformer export, and
  transactional failures.
- Property tests cover atom and bond insertion, whole-column and individual
  access, missing values, removal and all-missing normalization, conservative
  and reserved names, exact lengths, non-finite rejection, compatible unit
  conversion, incompatible units, transactionality, and a complete
  `conformational_entropy` visualization column.
- Canonical/custom coexistence regressions verify fixed canonical units,
  case-insensitive reserved names, custom-only property iteration, shared dense
  normalization and validation behavior, transactional replacement, remapping,
  and unchanged mmCIF occupancy/B-factor round trips.
- Macro construction tests cover one valid selected conformer alongside many
  unrelated invalid conformers, rejection when an invalid conformer is selected,
  and preservation of all source conformers.
- Public transformation tests cover dense-index compaction, equal-layout
  allocation rejection, complete position and atom-data transfer, cells,
  missing target state, and source immutability.
- Model/view hierarchy tests verify identical qualified identities, qualified
  chained navigation, pointer-identical explicit local-node borrows, and
  unchanged coordinate allocation.

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
- v11: Carry the connected molecule-definition invariant through model staging
  and keep rejection transactional before topology-bound state is published.
- v12: Store `Arc<Topology>` in topology-bound structure state, remove
  identity-specific errors, and use pointer-compatible sharing while preserving
  explicit `same_layout` as a separate comparison.
- v13: Flatten `Model` to direct topology, positions, optional cell, and
  column-oriented `AtomData`; remove configuration and structure-observation
  wrappers and remap positions and atom data together.
- v14: Return borrowed `Positions` directly from model views, make B-factors
  explicitly unitful length-squared quantities, and clarify extensible
  field-specific `AtomData` construction and empty/count semantics.
- v15: Add thin zero-copy model and model-view navigation over topology-owned
  atoms, bonds, instances, and instance-qualified SMCRA hierarchy.
- v16: Forward borrowed instance-qualified hierarchy views so model-level
  chained navigation cannot expose ambiguous definition-local relationship
  IDs accidentally.
- v17: Add exact-topology-bound `BondData`, conservative unit-aware scalar
  custom properties on atom and bond data, zero-copy model views, and complete
  atom/bond property remapping.
- v18: Reuse one private dense unit-aware scalar-column implementation for
  canonical occupancy, canonical B-factor, and custom properties while
  preserving dedicated canonical APIs, units, errors, and format semantics.
