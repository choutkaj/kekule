# Conjugation and resonance validation

Implementation branch: `codex/conjugation-resonance`. Reference: RDKit
2026.03.3, installed from the repository's pinned reference environment.
No RDKit runtime dependency is introduced.

Strict parity is **not yet complete**. The comparator has not been weakened.
Source-format failures and enumeration disagreements remain visible.

## Scientific comparisons

| Feature and coverage | Agreements | Disagreements | Errors |
|---|---:|---:|---:|
| Conjugation: smoke, full RDKit structures, 100-source PubChem/Enamine/PL-REX samples | 671 | 0 | 4 |
| Connected groups: same coverage | 671 | 0 | 4 |
| Enumeration: complete smoke | 21 | 4 | 0 |
| Enumeration: complete RDKit structures | 32 | 14 | 4 |
| Enumeration: 10-source PubChem sample | 19 | 1 | 0 |
| Enumeration: 10-source Enamine sample | 18 | 2 | 0 |
| Enumeration: 10-source PL-REX sample | 20 | 0 | 0 |

Counts refer to indexed source inputs, including different file representations,
not distinct compounds. Every successful enumeration observation contains all
32 flag masks. The four structure-corpus errors occur before perception: an
unsupported tetrahedral CFG assertion and three unsupported `R`/`R1` element
symbols. They are retained, not counted as agreement.

All observed enumeration disagreements occur at the 1,000-contributor limit.
The smoke differences, repeated across BMS-986142 and ZM374979 representation
variants in the structure corpus, replace one contributor with another tied
under all seven reference ranking criteria. RDKit collects candidates through
`std::unordered_map` and applies `std::sort` with no final tie breaker. Its exact
capped subset therefore depends on container/sort behavior. Kekule currently
uses deterministic indexed ordering for such ties. A separate audit of the
PubChem and both Enamine source representations likewise found exactly one
substituted contributor per differing profile, with equal values under all
seven criteria. All these differences remain strict failures. This audit does
not establish that an arbitrary capped mismatch is a tie.

Accepting only proven cutoff-tie substitutions, or reproducing the pinned
Windows build's container and sorting details, requires a parity-contract
decision. No tolerance or structure field has been removed to hide these cases.

Focused differential checks covered 800 molecule/flag combinations with exact
indexed charge/order collections and no failures. An additional 224 valid
stereochemical and aromatic-anion combinations and all 32 aryne flag masks
matched; one invalid aromatic
test string was rejected by RDKit and was not counted as agreement. Unit tests
also cover aromatic triple bonds, all option masks, absent/empty state,
materialization, reconstruction, invalidation, tombstones, hydrogen expansion,
atom renumbering, atomic failure, and preserved aromaticity/CIP/output.

## Validation commands

The implementation was checked with the repository's commands, using two build
jobs to limit Windows linker memory. Logs and detailed per-case reports are in
the ignored `target/` directory and normal `benchmarks/runs/` history.

- `cargo fmt --all -- --check`
- `cargo check --workspace --all-targets --all-features --locked`
- `cargo +1.89.0 check --workspace --all-targets --all-features --locked`
- `cargo clippy --workspace --all-targets --all-features --locked -- -D warnings`
- `cargo test --workspace --all-features --locked`: 1,503 passed, three existing ignored tests, including the 14-test conjugation/resonance regression suite.
- `cargo test --workspace --all-features --doc --locked`
- `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --all-features --no-deps --locked`
- `cargo test -p kekule-potentials --no-default-features --locked`
- `RUSTDOCFLAGS="-D warnings" cargo doc -p kekule-potentials --no-default-features --no-deps --locked`
- `cargo check --manifest-path fuzz/Cargo.toml --bins --locked`, including the new bounded resonance target.
- `cargo +nightly fuzz run resonance -- -runs=256 -max_len=4096 -seed=13`: passed with supplied PubChem seeds. The first launch lacked the sanitizer DLL; adding the installed MSVC runtime directory to that process's PATH resolved it. Seed-file line endings are trimmed before parsing.
- `cargo package -p kekule --locked --allow-dirty`, including package compilation verification.
- Companion `cargo package --locked --list` checks; `--allow-dirty` is needed for the trajectory crate's updated tests.
- All CI Python unittest suites and `node benchmarks/test_dashboard.cjs`.

The byte-for-byte license-copy check finds pre-existing CRLF differences in the
two OpenFF copies; their text is identical after newline normalization. No
license text was changed to conceal the difference. Linux-only CI execution is
not substituted by the Windows checks. The other registered fuzz targets were
compiled but not run; the runtime smoke exercise is scoped to the new target.
Full companion package verification remains disabled by repository CI until
the foundational crate is published, so their package file lists were checked.
Full bulk-corpus reference generation and comparison were not run; the bounded
samples above are the measured coverage.

The subsequent PR-readiness cleanup completed upstream copyright attribution
and expanded dense benchmark JSON construction without changing its values or
comparison rules. Formatting, diff checks, the resonance observation regression,
scoped benchmark clippy, and core package verification were repeated. The full
workspace test total above precedes this cosmetic cleanup.

## Review fixes

Permutation construction now reserves work for its stored vectors, depth scans,
and sorting before allocation. Checked arithmetic rejects an overflowing work
estimate. Regressions exercise a million-permutation request with a small budget,
sorting-budget exhaustion, overflow, and unchanged contributor priority.
Detached resonance groups are normalized by their smallest bond ID, with a
regression for reversed input order and stable atom/bond group indices. Invalid
empty groups and duplicate members remain installation errors.

Review revalidation passed: 1,506 workspace tests (three existing ignored), all
42 separately run doctests, formatting/diff checks, stable and Rust 1.89
workspace checks, warning-free workspace clippy and documentation, all fuzz
target compilation, and core package verification. Detailed logs are
`target/resonance-review-*.log`.

The optional external comparisons and runtime fuzz smoke were not repeated for
these fixes; focused regressions exercise the two defects and the full workspace
suite checks integration. Python/Node tooling, companion-only minimal-feature
checks and companion package inventories were unchanged and retain their earlier
validation results above. Linux execution remains unavailable on this Windows
host. The previously recorded cutoff-parity differences remain unresolved.
