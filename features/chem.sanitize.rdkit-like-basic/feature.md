# Basic RDKit-like Sanitization Pipeline

## Summary

Provide an explicit opt-in sanitization pipeline for common small molecules.

## Behavior/API

- Exposes `SanitizeOptions`, `SanitizeReport`, `SanitizeError`, and `sanitize_small_molecule`.
- Runs valence, ring set, and aromaticity perception according to options.
- Does not run automatically from file parsers.

## Implementation Notes

- The pipeline is conservative and returns reports for caller inspection.
- It operates on `SmallMolecule` while using shared core graph algorithms internally.

## Validation

- Unit tests cover parse-without-sanitize behavior and explicit sanitization.
- RDKit goldens are planned but not committed yet.

## Out Of Scope

- Full RDKit sanitization parity, kekulization, stereochemistry assignment, cleanup transforms, and organometallic handling.

## Revision Notes

- v1: Initial explicit sanitization pipeline.
