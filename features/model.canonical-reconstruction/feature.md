# Canonical Molecular Reconstruction

## Summary

Let external persistence adapters reconstruct the complete canonical molecular
state that Kekule owns without source reparsing, perception, sanitization,
dummy chemistry, or stable-ID renumbering.

## Behavior/API

- `PerceptionState::builder` constructs detached state with exact optional
  valence, ring, and aromaticity sections plus complete CIP assignments.
- Public immutable section views preserve absent valence versus model-neutral
  valence, named `ValenceModel`, complete implicit-H assignments, membership
  with or without `RingSet`, imported or model-perceived aromaticity, aromatic
  atom/bond sets, and CIP descriptor assignments.
- `RingMembership::from_slot_flags` and `RingSet::from_rings` reconstruct the
  core ring values consumed by the detached state.
- `Molecule::install_perception_state` validates the entire detached state
  against live graph/stereo slots and atomically replaces prior perception.
- `SmcraHierarchy::set_residue_component_ids`, `chain_props_mut`,
  `residue_props_mut`, and `atom_site_props_mut` restore enriched child state
  without exposing parentage, child arrays, lookup maps, or mutable records.
- `Molecule::stereo_group_slot_count`, `stereo_group_slots`, and
  `append_stereo_group_tombstone` expose and replay the exact stable group-slot
  layout. Appending a deleted slot preserves CIP state.
- Stereo-element insertion rejects pre-grouped values transactionally;
  persistence adapters create ungrouped elements and establish membership only
  through `add_stereo_group`. Removed elements are returned ungrouped.
- Live stereo groups are always nonempty. Removing or pruning their final
  member tombstones the group; partial pruning preserves the stable live group.

## Implementation Notes

Perception construction rejects duplicate map/set entries and fixed-width
component overflow. Whole-state installation validates exact ring membership
slot lengths, live atom/bond/stereo references, simple-cycle ring coherence,
basis coverage, and all optional sections before mutation. Algorithm
diagnostics are not reconstructed because they are not semantic perception
state. Installation does not add dependencies between otherwise independent
current sections, so imported aromaticity does not require a perceived ring
basis.

Persistence adapters reconstruct in this order:

1. atoms, bonds, conformers, and their stable slots;
2. stereo elements without group side effects;
3. live stereo-group slots and tombstones in exact order;
4. stereo bond marks;
5. SMCRA hierarchy plus child component IDs and properties;
6. complete perception state last.

Installing perception last ensures ordinary graph and stereo construction
retains its normal invalidation behavior. Later topology/chemistry mutation
clears topology-derived perception, stereo mutation clears CIP, and
coordinate/property-only mutation remains perception-neutral.

## Tests

- Downstream-style tests publicly export and reinstall empty and fully populated
  perception, including model-neutral valence, both ring modes,
  imported/perceived aromaticity, and CIP.
- Malformed/tombstoned references, duplicate entries, slot mismatches,
  malformed cycles, and inconsistent basis membership fail without changing
  prior perception.
- Hierarchy tests preserve distinct residue name/label/author component IDs,
  arbitrary child property maps, atom-site metadata, and invalid-ID rollback.
- Stereo tests replay live/interior-tombstone/live/trailing-tombstone layouts,
  preserve live and next group IDs, tombstone final-member groups, retain
  partially pruned groups, reject empty groups, and preserve CIP on tombstone
  append. They also prove pre-grouped insertion preserves element slots and
  installed perception, while grouped removal returns a reinsertable ungrouped
  element.
- Independently built original/reconstructed small and macro topologies satisfy
  `Topology::same_layout`.
- The sibling MolStudio project adapter is locally tested through the ignored
  path patch with complete perception, enriched hierarchy, and stereo
  tombstone archives; load performs no source reparse or reperception.

## Out Of Scope

- Serde derives or a Kekule-owned archive/wire format.
- General molecule or topology reconstruction DTOs.
- Persistence of process-local `Arc<Topology>` sharing relationships.
- Public incremental perception mutation.
- Mutable SMCRA parentage, child arrays, records, or lookup internals.
- Repairing malformed historical state through coercion or reperception.

## Revision Notes

- v1: Add focused checked reconstruction for installed perception, enriched
  SMCRA child state, and exact stereo-group tombstone layouts, with MolStudio
  consumer proof.
- v2: Enforce ungrouped stereo-element insertion and detached removal so
  canonical replay establishes every relation only through
  `add_stereo_group`.
- v3: Reconstruct semantic ring membership and basis only, excluding algorithm
  diagnostics from canonical molecule state.
