# RDKit-like Valence Perception

## Summary

Provide conservative valence perception for common organic molecules.

## Behavior/API

- Exposes `ValenceModel`, `ValenceReport`, `ValenceIssue`, and `perceive_valence`.
- Computes explicit valence from bond order and explicit hydrogens.
- Assigns implicit hydrogens when a common allowed valence can be selected.
- Reports unsupported elements or valence excesses instead of silently accepting them.

## Implementation Notes

- The current model covers common organic elements plus selected charged salts, boranes, hypervalent halogens, and organometallic cases exposed by external PubChem validation.
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
