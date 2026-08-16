# Agent rules

These rules apply to contributors and AI agents working in this repository.

## Workflow

1. Read `ARCHITECTURE.md` and keep its object boundaries and invariants intact.
2. Keep changes scoped; do not mix unrelated cleanup into a functional change.
3. Add or update a regression test for every defect fix or behavior/API contract change.
4. Run the applicable Rust formatting, check, clippy, test, documentation, and package checks before handoff. Report every applicable command not run and why.
5. Use optional external-reference benchmarks only when they are scientifically useful. They are not routine CI or release gates.
6. Do not modify `README.md` without the human's consent.

## Branches and commits

- Do not push feature work directly to `main`; use a short-lived branch based on current `main`.
- Preserve unrelated user changes in a dirty worktree.
- Keep commits reviewable and end every commit message with:

  ```text
  Co-authored-by: codex <codex@openai.com>
  ```

## Scientific tooling

- Keep parsing separate from interpretation, normalization, perception, validation, and preparation as defined by `ARCHITECTURE.md`.
- RDKit, Biopython, DSSP, and similar tools are benchmark references only, never Rust runtime dependencies.
- Benchmark fixtures must be externally supplied and provenance-pinned. Toy molecules belong only in focused unit regressions.
- Do not weaken comparisons, remove asserted fields, or regenerate goldens merely to hide a mismatch.
- Do not claim a check, benchmark, workflow, or repository setting was verified unless it was actually inspected or run.
