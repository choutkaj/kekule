# RDKit-like Valence Perception

## Summary

Provide conservative valence perception for common organic molecules.

## Behavior/API

- Exposes `ValenceModel`, `ValenceReport`, `ValenceIssue`, and `perceive_valence`.
- Computes explicit valence from bond order and explicit hydrogens.
- Assigns implicit hydrogens when a common allowed valence can be selected.
- Reports unsupported elements or valence excesses instead of silently accepting them.

## Implementation Notes

- The current model covers common organic elements plus selected charged salts, boranes, group-14/group-15 main-group compounds, hypervalent halogens, flexible transition-metal coordination centers, and radical implicit-hydrogen suppression exposed by external PubChem validation.
- Perception state is marked fresh only after the pass completes.

## Validation

- Unit tests cover neutral organics, charged species, and valence error reporting.
- RDKit-generated goldens compare valence status, explicit valence, and implicit hydrogen assignments for external PubChem fixtures.

## Out Of Scope

- Full RDKit valence tables, organometallics, query atoms, valence tautomer handling, and sanitization orchestration.

## Revision Notes

- v1: Conservative valence perception.
- v2: Validation contract narrowed to valence-specific outputs and passes the RDKit-backed `tiny` corpus.
- v3: Add corpus-driven RDKit-compatible valence cases for charged halides, boron anions, alkali counterions, hypervalent halogens, and simple mercury salts.
- v4: Expand corpus-driven RDKit-compatible valence cases for PubChem-100 salts, silicon, phosphonium, and selected metal centers.
- v5: Generalize PubChem-1000-driven valence handling for transition-metal coordination, group-14/group-15 heavy elements, oxonium centers, chalcogens, and radicals; PubChem-1000 still requires further table coverage.
