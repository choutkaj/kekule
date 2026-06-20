# RDKit-like Sanitization Pipeline

## Summary

Provide an explicit opt-in sanitization pipeline for common small molecules.

## Behavior/API

- Exposes `SanitizeOptions`, `SanitizeReport`, `SanitizeError`, and `sanitize_small_molecule`.
- Runs valence, ring set, and aromaticity perception according to options.
- Commits changes only after every requested pass succeeds; any error leaves the input exactly unchanged.
- Marks requested successful passes fresh and ensures skipped passes are not fresh. Aromaticity may compute rings internally, but an unrequested ring result is not retained or exposed.
- Does not run automatically from file parsers.

## Implementation Notes

- The pipeline stages work on a clone and atomically replaces the caller's molecule after success.
- The pipeline is conservative and returns reports for caller inspection.
- It operates on `SmallMolecule` while using shared core graph algorithms internally.
- Applies sanitization-only charge cleanup for hypervalent oxyhalogen patterns before valence perception.
- Its valence, ring, and aromaticity passes are compared together against each required corpus.
- Inherits the current valence and aromaticity improvements, including radical implicit-hydrogen handling, imported aromatic SMILES handling, and conservative unsupported-ring behavior.

## Validation

- Unit tests cover parse-without-sanitize behavior, every option combination, cleanup invalidation, idempotence, and exact rollback after valence or aromaticity failure.
- RDKit-generated goldens compare sanitized atom state for external PubChem fixtures.

## Out Of Scope

- Full RDKit sanitization parity, kekulization, stereochemistry assignment, cleanup transforms, and organometallic handling.

## Revision Notes

- v1: Explicit sanitization pipeline.
- v2: Validated through the corrected valence, ring, and aromaticity passes.
- v3: Add RDKit-like oxyhalogen cleanup and pass PubChem-100 through the corrected valence and aromaticity stack.
- v4: Incorporate broader PubChem-1000-driven valence and aromaticity behavior; PubChem-1000 remains pending on fused aromatic bond selection and remaining valence-table coverage.
- v5: Make sanitization transactional and define fresh/stale state outcomes for every option combination.
- v6: Sanitize imported aromatic SMILES with corrected aromatic valence and atom-contribution aromaticity behavior.
