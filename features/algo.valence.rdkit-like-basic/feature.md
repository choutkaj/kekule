# Basic RDKit-like Valence Perception

## Summary

Provide conservative valence perception for common organic molecules.

## Behavior/API

- Exposes `ValenceModel`, `ValenceReport`, `ValenceIssue`, and `perceive_valence`.
- Computes explicit valence from bond order and explicit hydrogens.
- Assigns implicit hydrogens when a common allowed valence can be selected.
- Reports unsupported elements or valence excesses instead of silently accepting them.

## Implementation Notes

- The initial model covers common organic elements and charges used by early fixtures.
- Perception state is marked fresh only after the pass completes.

## Validation

- Unit tests cover neutral organics, charged species, and valence error reporting.
- RDKit goldens are planned but not committed yet.

## Out Of Scope

- Full RDKit valence tables, organometallics, query atoms, valence tautomer handling, and sanitization orchestration.

## Revision Notes

- v1: Initial conservative valence perception.
