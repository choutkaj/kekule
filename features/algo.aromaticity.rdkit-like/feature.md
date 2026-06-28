# RDKit-like Aromaticity Perception

## Summary

Assign aromatic atom and bond flags for common organic ring systems using the RDKit-like model. This is a perception step and must remain separate from parsing.

## Behavior/API

- Exposes `AromaticityModel::RdkitLike`, `perceive_aromaticity`, and `perceive_aromaticity_with_ring_options`.
- Requires or computes ring perception before assigning aromaticity.
- Marks supported aromatic atoms and bonds and sets aromaticity perception state to fresh.
- Clears prior aromatic flags deterministically before assignment.
- Can be run directly or through the explicit sanitization pipeline.
- Returns `UnsupportedElement` for an explicitly aromatic ring containing an unsupported element instead of silently accepting that representation.
- Returns `InvalidAromaticRepresentation` when imported aromatic bonds cannot be perceived back onto every participating atom.
- Propagates bounded ring-perception failures as `AromaticityError::RingPerception`.

## Implementation Notes

- Operates on the shared core `Molecule` graph.
- Uses per-ring cycle data from `algo.rings.sssr`.
- Integrates with the first-wave valence and ring-set perception stack.
- Applies a 4n+2 electron-count model for common C, N, O, S, Se, Te, and P rings and small fused ring components.
- Computes pi-electron counts from bond order, and uses an atom-contribution path for explicitly imported aromatic-bond rings.
- Uses conservative guards for small rings, hetero fused donors, lactone-like rings, and large macrocycles exposed by external PubChem validation.
- Treats unsupported ring elements as non-aromatic for the current model rather than failing the whole perception pass.
- Leaves unsupported or ambiguous systems non-aromatic rather than claiming full RDKit parity.

## Validation

- Unit tests cover common monocyclic organic rings and stale-state behavior.
- RDKit-generated goldens compare aromatic atom and bond flags for external PubChem fixtures.

## Out Of Scope

- Full RDKit aromaticity parity.
- RDKit-like aromatic bond selection for all fused systems; PubChem-1000 exposes cases where atom aromaticity and bond aromaticity diverge.
- Valence perception, sanitization, kekulization, stereochemistry, and parser behavior.
- Runtime RDKit dependency.

## Revision Notes

- v1: Aromaticity perception for common organic rings.
- v2: Document integration with explicit sanitization and ring/valence perception.
- v3: Per-ring fused aromaticity heuristic passes the RDKit-backed `tiny` corpus; broader required corpora remain pending.
- v4: Add fused-component aromaticity and order-based electron counting for external PubChem fused-ring systems.
- v5: Refine fused heteroaromatic handling and conservative ring-size/electron-count guards to pass PubChem-100.
- v6: Add chalcogen heteroaromatic support and refine fused donor eligibility; PubChem-1000 still exposes fused aromatic bond-selection gaps.
- v7: Reject unsupported elements in explicitly aromatic ring representations and preserve caller state when sanitization propagates the error.
- v8: Count explicitly imported aromatic-bond rings with atom contributions so lowercase aromatic SMILES sanitize without treating every aromatic bond as a localized double bond.
- v9: Propagate configurable structured ring resource limits before mutating aromatic flags.
- v10: Narrow saturated tertiary amine fused-ring clearing to carbon-substituted amines, preserving N-O substituted lactam aromaticity.
