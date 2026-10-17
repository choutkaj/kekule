# Agent rules

These rules apply to contributors and AI agents working in this repository.

## Workflow

1. Read `ARCHITECTURE.md` and keep its object boundaries and invariants intact.
2. Keep changes scoped; do not mix unrelated cleanup into a functional change.
3. Run the applicable Rust formatting, check, clippy, test, documentation, and package checks before handoff. Report every applicable command not run and why.
4. Use optional external-reference benchmarks only when they are scientifically useful.
5. Do not modify `README.md` without the human's consent.

## Branches and commits

- Do not push feature work directly to `main`; use a short-lived branch based on current `main`.

## Scientific tooling

- RDKit, gemmi, Biotite, DSSP, and similar tools are benchmark references only, never Rust runtime dependencies.
- Benchmark fixtures must be externally supplied. Toy molecules belong only in focused unit regressions.
- Do not weaken comparisons, remove asserted fields, widen tolerances, or raise known-difference bounds merely to hide a mismatch.