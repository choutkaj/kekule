# RDKit-like Sanitization Pipeline

## Summary

Provide an explicit opt-in sanitization pipeline for common small molecules.

## Behavior/API

- Exposes `perception::{SanitizeOptions, SanitizeReport, SanitizeError, sanitize, sanitize_with_options, sanitize_with_ring_options}`.
- Runs valence, ring set, aromaticity, and stereo perception according to options.
- The stereo stage validates existing local stereo and assembles explicit source
  marks; coordinate-only stereo assignment remains an explicit
  `stereo.perception` operation.
- Commits changes only after every requested pass succeeds; any error leaves the input exactly unchanged.
- Propagates valence, ring, aromaticity, and stereo failures through
  `SanitizeError` without committing staged mutations.
- Marks requested successful passes fresh and ensures skipped passes are not
  fresh. Aromaticity may compute rings internally, but an unrequested ring
  result is not retained or exposed.
- Normalizes neutral aromatic nitrogen with one perceived donor hydrogen to RDKit-style `nH` atom state: one explicit atom hydrogen, zero implicit hydrogens, and no implicit-hydrogen suppression.
- Does not run automatically from file parsers.

## Implementation Notes

- The pipeline stages work on a clone and atomically replaces the caller's molecule after success.
- The pipeline is conservative and returns reports for caller inspection.
- It operates on `SmallMolecule` while using shared core graph algorithms internally.
- The public facade is `perception`; lower-level sanitizer internals are not root-level API.
- Applies sanitization-only charge cleanup for hypervalent oxyhalogen patterns before valence perception.
- Performs the aromatic nitrogen hydrogen normalization after both valence and aromaticity perception succeed, preserving total hydrogen count while exposing RDKit-like sanitized atom fields.
- Runs stereo perception after aromaticity so local stereo validation and
  source-mark assembly can use sanitized valence and hydrogen semantics.
- Its valence, ring, aromaticity, and stereo passes are compared together against each required corpus.
- Inherits the current valence and aromaticity improvements, including radical implicit-hydrogen handling, imported aromatic SMILES handling, and conservative unsupported-ring behavior.

## Validation

- Unit tests cover parse-without-sanitize behavior, every option combination,
  cleanup invalidation, idempotence, default and skipped stereo perception,
  coordinate-only stereo staying outside sanitization, and exact rollback after
  valence, aromaticity, or stereo failure.
- RDKit-generated goldens compare sanitized atom state for external PubChem fixtures.

## Out Of Scope

- Full RDKit sanitization parity, kekulization, exact CIP assignment,
  coordinate-derived stereo assignment, cleanup transforms, and organometallic
  handling.

## Revision Notes

- v1: Explicit sanitization pipeline.
- v2: Validated through the corrected valence, ring, and aromaticity passes.
- v3: Add RDKit-like oxyhalogen cleanup and pass PubChem-100 through the corrected valence and aromaticity stack.
- v4: Incorporate broader pubchem-1k-driven valence and aromaticity behavior; pubchem-1k remains pending on fused aromatic bond selection and remaining valence-table coverage.
- v5: Make sanitization transactional and define fresh/stale state outcomes for every option combination.
- v6: Sanitize imported aromatic SMILES with corrected aromatic valence and atom-contribution aromaticity behavior.
- v7: Accept explicit ring-work limits and preserve transactional rollback on ring resource errors.
- v8: Move the public small-molecule sanitizer API under the `perception` facade.
- v9: Add PubChem-100k as required broad-corpus validation evidence.
- v10: Normalize pyrrolic aromatic nitrogen donor hydrogens to RDKit-style sanitized `nH` atom state.
- v11: Add stereo perception as an explicit sanitizer stage with options,
  reporting, freshness-state handling, and transactional rollback on stereo
  perception issues.
