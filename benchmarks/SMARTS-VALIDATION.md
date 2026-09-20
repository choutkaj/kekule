# Molecular SMARTS validation

Measured on 2026-09-20 against **RDKit 2026.03.3**, using all query rows in the
pinned external fixtures. This is development evidence, not a runtime dependency
or a routine reference-tool CI gate. See [the protocol](SMARTS.md) and
[the public dialect](../crates/kekule/src/query/dialect.md).

## Behavioral comparisons

| Sources | Queries | Targets | Exact agreement | Shared parser rejection | Differences |
| --- | ---: | ---: | ---: | ---: | ---: |
| RDKit functional groups and reactivity tables | 518 | 8 | 4,144 | 0 | 0 |
| OpenFF Sage, water models, and Toolkit fixtures | 521 | 12 | 6,240 | 12 | 0 |

Each successful comparison asserts complete query-to-target mappings, with
chirality enabled and automorphisms retained, ascending-tag projections, query
graph sizes, and target aromatic atom flags. Neither atom order within a mapping
nor duplicate tag tuples is discarded. OpenFF targets have explicit hydrogen
vertices and explicitly selected MDL aromaticity.

The twelve shared rejections are the same malformed upstream
`[N:1](H:2)(H:3)` query on every target. The input is retained unchanged. There
were no resource-limit results or dialect exclusions in these two runs.

Many query/target pairs correctly have no matches. There are 39 nonempty RDKit
pairs and 258 nonempty OpenFF pairs. The latter include Constraints (16), Bonds
(53), Angles (43), ProperTorsions (47), ImproperTorsions (12), vdW (72),
ChargeIncrementModel (1), LibraryCharges (12), and VirtualSites (2).
This small target corpus does not establish universal RDKit equivalence or
validate every parameter's positive environment. Focused Rust regressions cover
the language and matching edge cases separately.

Complete, portable reports retain every input and both observations:

- [RDKit report](reports/smarts-core-rdkit.json.gz), SHA-256
  `16c879b4be3a00a5dbf410d1c0f5a3ff0ba4dce89ba1468e165554222a24a852`.
- [OpenFF report](reports/smarts-openff.json.gz), SHA-256
  `351f552b65e299d0f3c19d639a3efa51ee5cae932c76f8fce010105947e7c00d`.

Report paths are relative to the repository root. Source files, licenses, commit
pins, and byte checksums are in [RDKit provenance](smarts-fixtures/rdkit-queries/PROVENANCE.md)
and [OpenFF provenance](smarts-fixtures/openff-smarts/PROVENANCE.md).
Git text conversion is disabled for these checksum-pinned fixtures.

## Contract coverage

The focused `smarts_core` suite has 16 tests covering recursive anchoring and
independent embeddings, recursive negation and shared budgets, nested syntax,
Boolean bond precedence, hydrogen representation, valence and numeric rings,
stereo Boolean alternatives, atom reordering, tagged automorphisms, enumeration
beyond 1,000 results, streaming stop, topology snapshot identity, repeated
definition reuse, incompatible ring provenance, transactional MDL installation,
and aggregate graph-size errors. The existing query suite has 27 tests, including
tetrahedral carrier permutations and alkene stereo.

The documented Boolean-stereo semantics deliberately differ from RDKit's parser
behavior in some compound expressions. Literal atom/bond Boolean context is
preserved; these are explicit dialect differences, not a claim of reference
agreement. Reaction SMARTS, CXSMARTS, non-tetrahedral stereo, and vendor-specific
extensions remain excluded as documented.

## Repository checks

Validation used isolated target directories under `target/smarts*` to avoid stale
local build artifacts. The following completed successfully:

- `cargo fmt --all -- --check` and `git diff --check`.
- `cargo check --workspace --all-targets --all-features --locked`.
- `cargo +1.89.0 check --workspace --all-targets --all-features --locked`.
- `cargo clippy --workspace --all-targets --all-features --locked -- -D warnings`.
- `cargo test --workspace --all-features --locked`: **1,217 passed**, none failed
  or ignored, including the workspace doctests.
- `cargo test --workspace --all-features --doc --locked` separately.
- `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --all-features --no-deps --locked`.
- `cargo test -p kekule-potentials --no-default-features --locked` and its
  documentation build with warnings denied.
- `cargo check --manifest-path fuzz/Cargo.toml --bins --locked`.
- `cargo package -p kekule --allow-dirty --locked`, including package compilation;
  companion package file lists for `kekule-potentials` and `kekule-traj`.
- Four Python comparison/provenance tests, 21 existing dashboard tests, and the
  JavaScript dashboard checks.

Full companion package publication checks are intentionally not run, matching
the repository's CI policy: their dependency resolves through crates.io. Exact
package file lists were checked instead. `--allow-dirty` permits verification of
the uncommitted implementation; no packages were published.

Windows doctests initially hit PDB/linker errors with the disk nearly full; they
passed after removing task-generated build caches and setting
`RUST_TEST_THREADS=1`. Native Windows fuzz execution did not complete: the ASan
binary stalled during startup, and a sanitizer-free retry failed to link LLVM
coverage symbols. Both targets subsequently passed 256 seeded runs with ASan
under the installed Ubuntu/WSL nightly toolchain:

```text
cargo +nightly fuzz run smarts -- -runs=256 -max_len=4096 -seed=1
cargo +nightly fuzz run smarts_match -- -runs=256 -max_len=4096 -seed=2
```

Each starts with ten externally sourced patterns; see
[seed provenance](../fuzz/SMARTS-SEEDS.md). These short smoke runs check fuzz
integration and are not a substitute for longer fuzz campaigns. Unchanged fuzz
targets were compiled but not re-fuzzed in this implementation pass.

The result supplies the molecular query engine and tested SMIRNOFF matching
primitives. OFFXML interpretation, parameter precedence, charges, fractional bond
orders, improper/virtual-site assignment rules, and parameterized systems remain
the responsibility of the future parameterization layer.
