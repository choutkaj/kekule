# Agent rules

These rules apply to AI agents working in this repository.

## Mandatory workflow

1. Start from a feature ID under `features/`.
2. Read `ARCHITECTURE.md`, this file, and the feature directory before editing code.
3. Keep changes scoped to the requested feature.
4. Update feature metadata only when the implementation and validation status genuinely changed.
5. Run the relevant checks before handing work back.
6. Commit and push code in logical chunks.

## Anti-slop rules

1. Every nontrivial PR must reference a feature ID.
2. No feature can be marked validated without reference-generated golden data or documented manual validation.
3. No public API can be added without updating the feature spec or architecture docs.
4. The Rust library must not depend on RDKit or Biopython at runtime.
5. Reference code may be used for behavioral comparison. Copied code requires explicit license review and attribution.
6. Mutation must invalidate computed chemistry state.
7. Parsers must distinguish raw parsing from sanitization/perception.
8. Biomolecular labels belong in `BioHierarchy`, not in core `Atom`, unless chemically general.
9. Algorithms must document assumptions and edge cases.
10. The dashboard is generated from feature metadata; do not hand-edit `features/DASHBOARD.md`.
11. Every tracked feature must have schema-v2 `feature.toml` and canonical `feature.md`; do not recreate split feature docs.
