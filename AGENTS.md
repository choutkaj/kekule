# Agent rules

These rules apply to contributors and AI agents working in this repository.

## Workflow

1. Read `SPEC.md` and keep its object boundaries and invariants intact.
2. Keep changes scoped; do not mix unrelated cleanup into a functional change.
3. Write error messages by the "Error messages" rules in `SPEC.md`. Applications show them to their users.
4. Run the gates below before handoff. Report every applicable command not run and why.
5. Use optional external-reference benchmarks only when they are scientifically useful.
6. Do not modify `README.md` without the human's consent.

## Gates

These are the commands CI runs (`.github/workflows/ci.yml`):

```bash
cargo fmt --all -- --check
cargo check --workspace --all-targets --all-features --locked
cargo check --manifest-path fuzz/Cargo.toml --bins --locked
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo clippy -p kekule-potentials --all-targets --no-default-features --locked -- -D warnings
cargo clippy -p kekule-openff --all-targets --no-default-features --locked -- -D warnings
cargo test --workspace --all-features --locked
cargo test -p kekule-potentials --no-default-features --locked
cargo test -p kekule-openff --no-default-features --locked
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --all-features --no-deps --locked
RUSTDOCFLAGS="-D warnings" cargo doc -p kekule-potentials --no-default-features --no-deps --locked
RUSTDOCFLAGS="-D warnings" cargo doc -p kekule-openff --no-default-features --no-deps --locked
cargo package --workspace --exclude kekule-bench --locked
```

When `tools/` or `benchmarks/` change, also run their Python suites:

```bash
python -m unittest discover -s tools/atomic-data -p test_generate.py
python -m unittest discover -s tools/openff -p "test_*.py"
python -m unittest discover -s tools/components -p "test_*.py"
python -m unittest discover -s benchmarks/tests -p "test_*.py"
```

## Branches and commits

- Do not push feature work directly to `main`; use a short-lived branch based on current `main`.

## Scientific tooling

- RDKit, gemmi, Biotite, DSSP, and similar tools are benchmark references only, never Rust runtime dependencies.
- Benchmark fixtures must be externally supplied. Toy molecules belong only in focused unit regressions.
- Do not weaken comparisons, remove asserted fields, widen tolerances, or raise known-difference bounds merely to hide a mismatch.