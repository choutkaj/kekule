# Agent rules

These rules apply to AI agents working in this repository.

## Mandatory workflow

1. Start from either a feature ID under `features/` or one numbered remediation stage in `fixes.md`.
2. Read `ARCHITECTURE.md`, this file, and the selected feature directories, before editing code.
3. Keep ordinary work scoped to one feature. Some work may span multiple inter-dependent features, but it must remain logically scoped and feature-oriented. List every affected feature ID.
4. Add or update a regression test that demonstrates the defect or contract being changed.
5. Update feature metadata only when implementation behavior, public API, or the validation contract genuinely changed.
6. Run every applicable check before handoff and report every command that was not run, with the reason.
7. Commit and push code in logical chunks. End every commit message with:

   ```text
   Co-authored-by: codex <codex@openai.com>
   ```

## Anti-slop rules

1. Every nontrivial PR must reference a feature ID.
2. No feature can be marked validated without current reference-generated golden data or documented manual evidence accepted by the validation harness.
3. Do not weaken comparisons, remove asserted fields, delete regression tests, or regenerate goldens merely to make a failure disappear.
4. No public API can be added or changed without updating the feature spec or architecture documentation.
5. The Rust library must not depend on RDKit or Biopython at runtime.
6. Reference code may be used for behavioral comparison. Copied code requires explicit license review and attribution.
7. Mutation must invalidate affected computed chemistry state. Failed transactional operations must not leave partial mutations.
8. Parsers must distinguish raw parsing from sanitization/perception and must return structured errors rather than panic on malformed input.
9. Writers must not silently coerce unsupported chemistry into a different representation; return a structured error instead.
10. Biomolecular labels belong in `BioHierarchy`, not in core `Atom`, unless chemically general.
11. Algorithms must document assumptions, edge cases, and resource limits.
12. The dashboard is generated from feature metadata; do not hand-edit `features/DASHBOARD.html`.
13. Every tracked feature must have `feature.toml` and canonical `feature.md`; do not recreate split feature docs.
14. Molecular validation fixtures must be externally supplied and provenance-pinned; inline toy molecules may be used only for focused unit regressions, not golden validation.
15. Do not claim a check, workflow, branch-protection rule, corpus result, or repository setting was verified unless it was actually inspected or run.
