# Core pipeline audit — 2026-09-12

This audit reviews ring membership and selected rings, valence/hydrogens,
aromaticity, represented stereo/CIP, and SMILES parsing, interpretation and
writing. It starts from `main` at `45c21f6b`, following the existing
[valence/ring audit](valence-rings-audit.md),
[aromaticity audit](aromaticity-audit.md), and
[stereo validation](stereo-validation.md). Those earlier corpus results are
historical evidence, not results rerun by this change.

The changes preserve the `Graph`/`Perception` boundary. No chemical rule,
external fixture, golden output, or README was changed.

## Findings and changes

| Pipeline | Finding | Resolution |
| --- | --- | --- |
| Rings | Input bounds counted live atoms/bonds, but scratch arrays use all stable slots. Deleting atoms could bypass the allocation limits. | Check atom and bond slot counts, including tombstones, before allocating scratch arrays. Include them in initial work accounting. |
| Aromaticity | The six-ring search restriction limits subset size, but does not bound the number of connected subsets or the work needed to construct them. Installed ring sets bypass enumeration limits. | Add an independent, checked aromaticity work budget, covering graph slots, ring entries, fusion discovery, subset visits and workspace copies. Exhaustion is an error and restores the previous perception. |
| Aromaticity | Fusion discovery compared every pair of rings, including independent ring systems. Each component's subset traversal also allocated an exclusion array for the complete molecule's ring set. | Index candidate rings by shared bonds; retain the exactly-one-shared-bond rule and the 24-bond fusion cutoff. Store only exclusions used by the current traversal. |
| SMILES reading | Interpretation rescanned every atom, bond and tetrahedral assertion for each connected component. A record of isolated atoms therefore incurred quadratic work. Completed reports were then copied field by field. | Partition source indexes once by parsed connectivity, retain component-local work, and move completed reports directly into the result. Remove the redundant interpretation wrapper and forwarding function. |
| Valence | Assignments were accumulated in a vector and immediately copied into the installed map. | Build the map directly, retaining transactional installation. |
| Stereo/CIP | Carrier sorting cloned complete expanded ligand trees just to order a handful of carriers. | Sort borrowed signatures and copy only carrier identities. Sequence rules and auxiliary-descriptor semantics remain unchanged. |
| SMILES writing | Reviewed graph projection, hydrogen/isotope handling, traversal and directional-constraint bounds, and supported stereo rejection paths. | No writer algorithm change was warranted by this pass; integrated round-trip coverage exercises the affected upstream pipelines. |

## Public behavior

`RingPerceptionOptions::max_atoms` and `max_bonds` now explicitly bound storage
slots, rather than just live entities. An edited sparse graph can therefore fail
a limit that its live count satisfies. This is intentional: stable identifiers
must not let deleted storage evade the bound.

`perception::aromaticity` adds `AromaticityOptions` and
`perceive_aromaticity_with_options`. The options contain existing ring options
and `max_total_work`, defaulting to 5,000,000 work units. The existing entry
points retain their signatures and use that default. Ring options apply only
when rings must be computed; the aromaticity budget applies in both cases.
`AromaticityError::ResourceLimit` reports the selected budget. Callers needing a
larger search can explicitly select a larger budget and handle failure.

Work units are algorithmic bookkeeping, not elapsed-time or byte limits. They
bound search and scratch-work growth; they do not replace process-level memory
limits or change the supported aromaticity model. Transactional state snapshots
still have cost proportional to existing perception.

## Regression coverage

- Deleted-slot ring limits reject before scratch allocation and preserve existing
  perception for atom, bond and total-work failures.
- Fusion discovery handles 10,000 independent rings with work proportional to
  their ring entries; densely overlapping families and connected-subset searches
  stop at explicit budgets.
- The existing independent bit-mask oracle still checks every graph on five
  vertices and every subset size for exact connected-subset enumeration without
  duplicates. Existing fused chemistry tests retain their full assertions.
- Public aromaticity tests sweep budgets through successive failure stages,
  including newly installed rings and partially replaced aromaticity, and verify
  restoration of valence, rings, aromaticity and CIP. Installed rings skip only
  ring-enumeration limits.
- A 1,280-component SMILES regression preserves complete atom/bond source maps,
  overlapping source spans, source order, represented stereo and created-element
  reports. Its input includes interleaved branch components and ring closures.
- Integrated parse/perceive/CIP/canonical-write tests cover bridged and symmetric
  rings, fused aromatics, charged heterocycles, nitro groups, tetrahedral stereo,
  conjugated E/Z and isotope-bearing stereo. They check descriptors, skipped
  assignments, canonical fixed points and perception idempotence.

## Validation

All required checks passed on the Windows MSVC host. The workspace test command
reported **1,181 passed, zero failed, zero ignored**, including its doctests.
The separate documentation-test command also passed. Documentation ran with
`RUSTDOCFLAGS=-D warnings`.

```text
cargo fmt --all -- --check
cargo check --workspace --all-targets --all-features --locked
cargo +1.89.0 check --workspace --all-targets --all-features --locked
cargo check --manifest-path fuzz/Cargo.toml --bins --locked
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo test --workspace --all-features --locked
cargo test --workspace --all-features --doc --locked
cargo doc --workspace --all-features --no-deps --locked
cargo test -p kekule-potentials --no-default-features --locked
cargo doc -p kekule-potentials --no-default-features --no-deps --locked
cargo package -p kekule --locked --allow-dirty
cargo package -p kekule-potentials --locked --list --allow-dirty
cargo package -p kekule-traj --locked --list --allow-dirty
git diff --check
```

All six packaged license copies match by SHA256. `--allow-dirty` checks the
uncommitted changes; the foundational package was compiled from its packaged
contents. Full companion package builds remain deferred under the repository's
CI policy until the foundational crate is published; their file sets were checked.
Nightly fuzz execution was not run: it is a separate Linux/nightly workflow.
Every registered fuzz binary compiled. No applicable local formatting, check,
Clippy, test, documentation or package gate was skipped.

The optional reference checks used freshly built probe executables and unchanged
external fixtures and goldens:

- All **11** selected smoke targets passed, covering fast/selected rings,
  valence, aromaticity, represented/perceived stereo, CIP, and the four SMILES
  read/write modes. They produced **48 fixture-feature matches**.
- The direct RDKit **2026.03.6** SMILES checker passed **1,000/1,000** distinct
  PubChem-1k inputs across **3,868 encodings**, with zero differences or reference
  errors. It checks both writers, isotope/map/charge/stereo graph identity,
  randomized atom order, and canonical fixed points.
- The standard PubChem-100 SDF runner panicked for each of the four ring,
  valence and aromaticity targets on its existing one-connected-molecule
  assumption (`expected one connected molecule, found 3`). These attempts are
  recorded as failures, not successful comparisons. This is the same benchmark
  harness limitation described in the earlier audits; no production molecule
  invariant or golden was weakened to accommodate it.

The existing component-aware valence/ring diagnostic was copied into the audit
directory and rebuilt against the changed library. It retains complete
source-ordered atom/bond records and compares every expected field; ring atom
lists and the ring family are sorted using the existing comparison convention.
All **144** assessed fixture-feature SHA256 checks match the pinned **2026.03.3**
goldens. Results on original, unmodified SDF input are:

| Corpus | Feature | Exact record matches | Other outcomes |
| --- | --- | ---: | --- |
| PubChem-100 | Each of fast rings, selected rings, valence, aromaticity | 100/100 | None |
| PubChem-1k | Each of fast rings, selected rings, aromaticity | 999/1,000 | CID24959 source-stereo rejection |
| PubChem-1k | Valence | 996/1,000 | CID24959 source-stereo rejection; three raw representation differences |
| PubChem-100k | Aromaticity | 99,991/100,000 | Nine existing valence rejections before aromaticity |

There are no aromatic-mask mismatches and no new work-limit failures. The nine
PubChem-100k failure identities match the preceding valence/ring audit exactly.
CID24959 declares V2000 stereo code 1 on a double bond, which the current source
contract rejects. These remain recorded failures.

The three raw valence differences are CIDs **173868**, **147028** and **423442**,
the same inputs identified in the preceding audit. Publication normalizes their
oxyhalogen representations. A separate live RDKit **2026.03.3** check on their
exact exported canonical graphs matches every valence field for all three. This
does not convert their original-golden differences into exact original-input
matches. Raw expected/actual records and the separate identical-graph comparison
remain in the report directory.

Reproduction: `cargo xtask benchmark --benchmark <feature> --corpus smoke` for
the feature IDs in `target/core-pipeline-audit/benchmarks.json`;
`cargo build -p xtask --examples --locked`, then
`uv --cache-dir target/uv-cache run --offline --python 3.13 benchmarks/reference/rdkit/compare_smiles.py --corpus pubchem-1k --variants 3 --probe target/core-pipeline-audit/smiles_write_probe.exe --output target/core-pipeline-audit/smiles-pubchem-1k.json`.
The probe is copied from the freshly built executable to keep it immutable during
comparison. Use a fresh report path when reproducing the run.

For the SDF diagnostic, build
`target/core-pipeline-audit/parity-probe/Cargo.toml` with `--release --locked
--offline --target-dir target/core-pipeline-audit/parity-build`, then run its
`valence-rings-parity-probe.exe` with repository root, report directory, feature
ID and corpus IDs as positional arguments. Fixture validation is in
`target/core-pipeline-audit/verify_fixtures.py`. The three identical-graph valence
comparisons use `--export-canonical-valence` and the existing
`target/valence-rings-audit/check_canonical_valence.py` with the pinned RDKit runtime.

Command logs, exit statuses, reference reports, input/probe hashes and the
component-aware diagnostic are retained under `target/core-pipeline-audit/`.

## Scope of the production claim

These fixes close the identified defects and reduce avoidable copying. They do
not establish universal chemical correctness or unrestricted input support.
Existing contracts remain relevant: selected rings are an RDKit-like family,
not a promised mathematical minimum basis; zero/dative bonds are excluded from
chemical rings; dative valence has no donor/acceptor semantics; aromaticity uses
the documented bounded fusion model; CIP supports the represented geometries
and explicit ranking limits; and SMILES rejects configurations it cannot preserve.
See [stereo support](stereo-support.md) for the supported geometry/export matrix.

No new production runtime dependency was introduced. RDKit remains an optional
external scientific reference.
