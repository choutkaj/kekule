# Agent rules

These rules apply to contributors and AI agents working in this repository.

## Workflow

1. Read `ARCHITECTURE.md` and keep its object boundaries and invariants intact.
2. Keep changes scoped; do not mix unrelated cleanup into a functional change.
3. Add or update a regression test for every defect fix or behavior/API contract change.
4. Run the applicable Rust formatting, check, clippy, test, documentation, and package checks before handoff. Report every applicable command not run and why.
5. Use optional external-reference benchmarks only when they are scientifically useful.
6. Do not modify `README.md` without the human's consent.


## Tests

- Name a test after the behavior it pins and place it in the module that owns that behavior, not after the defect or audit that prompted it.
- Assert exact expected values; do not settle for `is_ok()`, non-empty, inequality, or substring checks when the exact result is known.
- Add a regression molecule as a row in an existing table or fixture (for example `crates/kekule/tests/fixtures/perception/`) rather than as a new one-off test, and record which independent reference verified its expected values.
- Properties that must hold for every molecule belong in `crates/kekule/tests/invariants`; a known defect there is listed in `KNOWN_DEFECTS`, never hidden by weakening the invariant.

## Branches and commits

- Do not push feature work directly to `main`; use a short-lived branch based on current `main`.

## Scientific tooling

- RDKit, gemmi, Biotite, DSSP, and similar tools are benchmark references only, never Rust runtime dependencies.
- Benchmark fixtures must be externally supplied. Toy molecules belong only in focused unit regressions.
- Do not weaken comparisons, remove asserted fields, widen tolerances, or raise known-difference bounds merely to hide a mismatch. Explain every accepted difference in `benchmarks/known-differences.toml`.
- ChEMBL records are CC BY-SA; never copy them into `crates/*/tests/fixtures`, which ship under the crates' licences. Take regression molecules from PubChem, the CCD or other permissively licensed sources.
