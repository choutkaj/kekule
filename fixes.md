# Molecules Audit Remediation Plan

This file converts the repository audit into a staged execution plan for Codex `/goal` runs.

- **Repository:** `choutkaj/molecules`
- **Audit baseline:** `main` at `323ee995d8b2feb7a5422808e4ad2d4c30d8a5b3`
- **Status:** not started
- **Execution model:** one stage per branch and pull request
- **Ordering rule:** do not begin a stage until all prerequisite acceptance gates pass

The order is intentional. Stage 1 repairs the validation system so later chemistry fixes cannot be declared successful using stale or incomplete evidence.

## How Codex should use this file

For each stage:

1. Create a dedicated branch.
2. Read `AGENTS.md`, `ARCHITECTURE.md`, this file, and every affected feature directory.
3. Paste the stage's `/goal` block into Codex.
4. Add a regression test that demonstrates the audited defect.
5. Implement the smallest coherent fix that satisfies the full stage.
6. Run every required check.
7. Update feature metadata only when behavior, public API, or validation contract changes.
8. Regenerate validation evidence and the dashboard only through repository commands.
9. Open one focused pull request and include all commands run.
10. Mark the stage complete only after merge.

Do not combine unrelated stages. In particular, do not mix chemistry behavior changes with the module-splitting stage.

## Repository-wide constraints

- Keep RDKit and Biopython out of Rust runtime dependencies.
- Preserve raw parsing versus sanitization/perception boundaries.
- Keep `#![forbid(unsafe_code)]`.
- Do not weaken normalized comparisons to hide chemistry differences.
- Do not delete assertions or golden fields merely to obtain a pass.
- Do not silently coerce unsupported chemistry into another representation.
- Return structured errors for unsupported or malformed input.
- Externally validated molecular fixtures must remain provenance-pinned.
- Inline structures are acceptable for focused unit tests, but are not broad validation evidence.
- Never hand-edit `validated = true`, generated corpus status, or the generated dashboard.
- Every commit message must end with:

  ```text
  Co-authored-by: codex <codex@openai.com>
  ```

## Standard verification gate

Run this gate in every implementation stage unless the stage adds stricter commands:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
cargo xtask dashboard --check
cargo xtask skills --check
cargo xtask corpus check --corpus tiny --require-data
cargo xtask validate --feature all --corpus tiny
```

When feature metadata changes:

```bash
cargo xtask dashboard
cargo xtask dashboard --check
```

When non-`tiny` corpus data is available:

```bash
cargo xtask corpus check --corpus all --require-data
cargo xtask validate --feature all --corpus all
```

Every pull request must explicitly list commands that were not run and why.

## Finding closure map

| Finding | Closing stage |
|---|---:|
| Validation evidence can remain green after source, fixture, lock, golden, or generator changes | 1 |
| CI does not run reference validation and large corpora are not reproducible from a clean runner | 1 |
| Perception caches can remain accessible or falsely fresh after mutation | 2 |
| Sanitization can partially mutate a molecule before returning an error | 2 |
| Malformed V2000 and other parser inputs can panic | 3 |
| V2000 stereo and radical representations are not round-trip safe | 4 |
| Lowercase aromatic SMILES does not sanitize correctly | 5 |
| SMILES output can be unreadable by the parser or chemically lossy | 5 |
| mmCIF can merge distinct residues and drops Cartesian coordinates | 6 |
| Ring perception and recursive graph traversal have resource hazards | 7 |
| Large monolithic source files impede safe review and maintenance | 8 |
| Documentation, branch policy, and licensing remain release blockers | 9 |

## Stage checklist

- [ ] Stage 1 — Validation evidence and reproducible CI
- [ ] Stage 2 — Mutation, cache invalidation, and transactional sanitization
- [ ] Stage 3 — Panic-free parsers and fuzzing
- [ ] Stage 4 — V2000/SDF semantic round-trip fidelity
- [ ] Stage 5 — Aromatic SMILES and honest SMILES output
- [ ] Stage 6 — mmCIF residue identity and coordinates
- [ ] Stage 7 — Bounded ring perception and stack-safe graph traversal
- [ ] Stage 8 — Behavior-preserving modularization
- [ ] Stage 9 — Final audit closure and release readiness

---

# Stage 1 — Trustworthy validation evidence and reproducible CI

## Primary feature

- `validation.harness`

Supporting metadata may change for features whose validation contracts need migration.

## Goal

A green validation result must prove that the current implementation, current fixture bytes, current source lock, current golden bytes, current comparison code, and current reference metadata all participated in a successful non-empty comparison.

Old evidence must become invalid automatically when any material input changes.

## Required implementation

### 1. Introduce content-addressed evidence

Create a versioned evidence schema. A conservative whole-workspace implementation digest is acceptable initially and is safer than an incomplete feature-to-source dependency map.

Each passing feature/corpus target must be bound to at least:

- feature/corpus manifest;
- corpus `sources.lock.json`;
- every fixture listed in the manifest;
- every golden used by the comparison;
- Rust implementation source;
- Rust comparison and normalization source;
- applicable Cargo manifests and `Cargo.lock`;
- feature metadata;
- reference generator source;
- reference environment definition or lock;
- comparison mode and evidence schema version.

Requirements:

- Hash deterministic relative paths and file bytes.
- Sort all paths before hashing.
- Reject missing inputs.
- Store component hashes or a canonical evidence document plus a final SHA-256.
- Recompute evidence in `corpus_passed_at`.
- Treat a mismatched or unknown schema as not passed.
- Preserve timestamps on identical repeated `--update` runs.

### 2. Verify golden metadata fully

Extend golden validation to verify:

- `feature_id`;
- `corpus_id`;
- `fixture_path`;
- supported golden schema version;
- `reference.tool` equals manifest tool;
- `reference.version` equals manifest version;
- `input_sha256` equals the current fixture hash;
- generator/environment digest, if present, matches current pinned configuration.

Do not trust status-file copies of reference metadata without checking the golden document.

### 3. Enforce comparison contracts

- Parse and enforce `comparison_mode`.
- Reject unknown modes.
- Implemented features with required corpora must not use zero-fixture manifests.
- A pass requires `fixture_count > 0`.
- A pass requires `compared_count == fixture_count`.
- Missing goldens, skipped fixtures, or partial comparison are failures.
- A parse error may be compared as an expected result only when both reference and implementation explicitly emit the same structured record; it must not disappear through omission.

### 4. Make `--update` failure-safe

For selected targets:

1. Load existing status.
2. Mark selected entries pending in memory.
3. Run every comparison.
4. Insert pass evidence only for fully successful targets.
5. Remove or mark failed selected targets failed.
6. Write status, synchronized feature flags, and dashboard through temporary files and atomic rename.
7. Return nonzero when any selected target fails.

A failed current run must never inherit an old green result.

### 5. Add normal CI validation

Add to ordinary pull-request CI:

```bash
cargo xtask corpus check --corpus tiny --require-data
cargo xtask validate --feature all --corpus tiny
```

The command must execute implementation-versus-golden comparisons; checking stale status files is insufficient.

### 6. Make large-corpus validation reproducible

The manual or scheduled validation workflow must work from a clean checkout although non-`tiny` corpus data is ignored by Git.

It must:

- create the pinned RDKit or Biopython environment;
- build or retrieve selected corpus data from pinned source locks;
- verify corpus integrity with `--require-data`;
- run selected feature/corpus validation;
- upload concise logs on failure;
- avoid automatically committing generated data or statuses;
- use least-privilege workflow permissions.

If a `cargo xtask corpus prepare` command is added, use typed corpus configuration. Do not execute an arbitrary unchecked `build_command` string from repository data.

## Required regression tests

Add tests equivalent to:

- `evidence_invalid_after_rust_source_change`
- `evidence_invalid_after_fixture_change`
- `evidence_invalid_after_source_lock_change`
- `evidence_invalid_after_golden_change`
- `evidence_invalid_after_reference_generator_change`
- `golden_input_hash_must_match_fixture`
- `golden_reference_metadata_must_match_manifest`
- `required_manifest_rejects_zero_fixtures`
- `pass_requires_comparison_for_every_fixture`
- `unknown_comparison_mode_is_rejected`
- `failed_update_clears_old_selected_pass`
- `unchanged_update_is_idempotent`
- `status_and_dashboard_writes_remain_consistent`

Use temporary directories and mutate actual bytes to prove stale evidence is rejected.

## Metadata and documentation

- Increment `validation.harness` once for the contract change.
- Document evidence inputs, update semantics, and local reproduction.
- Regenerate status through the command. Do not preserve old status manually.
- Regenerate the dashboard.

## Acceptance gates

- Changing one implementation byte invalidates old evidence.
- Changing one fixture, golden, lock, or generator byte invalidates old evidence.
- A failed `--update` cannot leave the selected cell green.
- All applicable implemented `tiny` targets perform nonzero comparisons.
- Normal CI runs `tiny` reference validation.
- A clean manual runner can prepare and validate a non-`tiny` corpus.
- No chemistry expectation is weakened in this stage.

## Suggested commits

1. Evidence schema and digest tests.
2. Golden metadata and comparison-mode enforcement.
3. Failure-safe, idempotent status updates.
4. Tiny CI and reproducible large-corpus workflow.
5. Metadata, documentation, regenerated evidence, and dashboard.



---

# Stage 2 — Mutation, cache invalidation, and transactional sanitization

## Primary features

- `core.graph`
- `chem.sanitize.rdkit-like`

Supporting features:

- `algo.valence.rdkit-like`
- `algo.rings.fast`
- `algo.rings.sssr`
- `algo.aromaticity.rdkit-like`

## Goal

No mutation may leave computed state falsely fresh or expose a stale ring cache as current. Any failed sanitization must leave the original molecule unchanged.

## Required implementation

### 1. Centralize invalidation

Create internal mutation/invalidation paths.

At minimum:

- topology edits invalidate valence, rings, aromaticity, and stereo;
- element, formal charge, explicit hydrogen, radical state, aromatic representation, and bond-order edits invalidate affected chemistry;
- ring invalidation clears both `ring_membership` and `ring_set`;
- sanitization cleanup helpers cannot mutate raw storage without invalidation;
- coordinate-only and molecule-property-only changes do not invalidate topology.

A conservative all-chemistry invalidation is acceptable. False freshness is not.

### 2. Restrict perception-state mutation

Remove or restrict APIs that let callers mark state `Fresh` without running an algorithm.

Prefer:

- public read-only state;
- crate-private state setters;
- algorithm-owned cache writes.

Do not add an unsafe or untracked escape hatch.

### 3. Make sanitization transactional

Implement default sanitization against a working clone or explicit change set:

1. Clone/stage the molecule.
2. Apply cleanup transforms.
3. Run requested perception passes.
4. Commit only after all requested work succeeds.
5. On error, return the error and leave the caller's molecule exactly unchanged.

### 4. Define selective option semantics

For every `SanitizeOptions` combination:

- skipped passes remain absent or stale, never fresh;
- cleanup changes invalidate all affected old state;
- requested passes become fresh only after success;
- a failure in a later pass commits nothing;
- repeated successful sanitization is idempotent.

Either document cleanup as always-on or give it an explicit option. Do not hide it behind an unrelated flag.

### 5. Preserve parser boundaries

Raw readers remain raw. Do not call sanitization implicitly from SMILES, Molfile, SDF, or mmCIF readers.

## Required regression tests

- mutation clears current `ring_membership` and `ring_set`;
- atom/bond chemistry mutation cannot leave relevant state fresh;
- cleanup invalidates old state;
- every meaningful option combination has expected final states;
- failed valence sanitization preserves exact input equality;
- failed aromaticity sanitization preserves exact input equality;
- successful sanitization is idempotent;
- success commits cleanup and requested state together;
- raw parse retains absent perception.

## Documentation and metadata

- Update `ARCHITECTURE.md` if public mutation semantics change.
- Increment affected feature versions once for behavior/API changes.
- Document transactional failure semantics.
- Document cache lifetime and invalidation.

## Acceptance gates

- No stale ring cache is returned after mutation.
- No unrequested sanitize pass is marked fresh.
- Every sanitize error preserves the original Rust value exactly.
- Tiny chemistry comparisons still pass.
- No parser begins sanitizing implicitly.



---

# Stage 3 — Panic-free parsers and fuzzing

## Features

- `io.mol.v2000.parse`
- `io.sdf.v2000.parse`
- `io.smiles.parse`
- `io.mmcif.parse`

## Goal

Every public text parser must return a structured error rather than panic for any Rust `&str`, including malformed records, non-ASCII text, zero indices, truncated fixed-width fields, and extreme counts.

## Required implementation

### 1. Harden V2000/SDF arithmetic and slicing

- Replace `a - 1` and similar operations with `checked_sub`.
- Use checked addition/multiplication for block offsets.
- Parse fixed-width structural fields through byte-safe ASCII helpers.
- Never slice a UTF-8 string with unchecked byte ranges.
- Reject non-ASCII fixed-width structural fields where CTfile requires ASCII.
- Validate declared counts before allocation and indexing.
- Validate M-record pair counts and atom indices before mutation.
- Return record and line context in errors.

### 2. Harden SMILES cursor handling

- Keep all offsets at valid character boundaries.
- Reject malformed charges, maps, isotopes, empty brackets, and truncated symbols.
- Reject a pending bond with no following atom/ring endpoint.
- Fully consume supported bracket syntax.
- Return explicit unsupported-syntax errors for unknown bracket decorators.

Aromatic semantics are fixed in Stage 5; this stage fixes safety and malformed-input behavior.

### 3. Harden mmCIF tokenization

- Safely handle quotes, comments, and semicolon blocks.
- Reject ragged loops.
- Reject integer/float overflow.
- Preserve accurate line numbers.
- Document and implement safe default input limits if needed.

### 4. Add fuzz targets

Add non-runtime fuzz targets:

- `mol_v2000`
- `sdf_v2000`
- `smiles`
- `mmcif`

Objectives:

- no panic;
- no abort;
- no pathological allocation from small input;
- successful parse followed by applicable write/read does not panic;
- failures remain structured.

Seed with tiny corpus fixtures and focused malformed cases. Fuzz dependencies must not enter the runtime library dependency graph.

### 5. Add bounded CI smoke coverage

Run deterministic regression seeds or a bounded fuzz smoke job in CI. Keep long campaigns scheduled or manual.

## Required regression cases

- V2000 endpoint `0`;
- endpoint beyond atom count;
- truncated counts line;
- non-ASCII before fixed columns;
- very short atom/bond lines;
- overflowing/inconsistent counts;
- truncated M-record pairs;
- unmatched SMILES branch/ring;
- malformed `%` ring label;
- malformed bracket charge/map/isotope;
- malformed mmCIF quote/semicolon block;
- ragged mmCIF loop;
- numeric overflow.

Use `catch_unwind` only in tests to assert public APIs do not panic.

## Acceptance gates

- Every regression returns `Err` with useful location context.
- Bounded fuzz smoke runs without crashes.
- Unsafe code remains forbidden.
- Runtime dependencies remain clean.
- Tiny parse validation passes without removing comparison fields.



---

# Stage 4 — V2000 and SDF semantic round-trip fidelity

## Features

- `core.atom-bond`
- `io.mol.v2000.parse`
- `io.mol.v2000.write`
- `io.sdf.v2000.parse`
- `io.sdf.v2000.write`

## Goal

Every V2000 feature claimed as supported must survive parse/write/parse without changing meaning. Unsupported representations must fail explicitly.

## Required implementation

### 1. Introduce a lossless radical model

The current electron-count field cannot distinguish all MDL radical multiplicities. Introduce one authoritative representation, such as:

```rust
pub enum AtomRadical {
    Singlet,
    Doublet,
    Triplet,
}
```

Requirements:

- Map every supported `M  RAD` value according to a pinned CTfile specification.
- Make writer mapping the exact inverse.
- Provide algorithm helpers for unpaired-electron behavior.
- Remove or deprecate the old field rather than keep two unsynchronized representations.
- Update valence behavior deliberately.

### 2. Correct bond-stereo mappings

Implement and document a parser/writer mapping table for supported V2000 stereo codes.

Requirements:

- wedge-up and wedge-down round-trip;
- either/unknown round-trips when representable;
- double-bond stereo codes are not conflated with wedge stereo;
- unsupported combinations return a structured write error;
- supported parser/writer tables are inverse-tested.

Use the pinned specification as the source of truth.

### 3. Make metadata symmetric

Verify exact semantic symmetry for:

- formal charge;
- isotope;
- atom map;
- radical state;
- coordinates;
- supported bond order;
- supported bond stereo;
- SDF text fields;
- title/program/comment lines.

Reject unrepresentable values before emitting a record.

### 4. Transactional writing

Construct output only after validating the entire molecule. Do not expose partial records.

### 5. Strengthen validation

- Add external fixtures with supported radical/stereo records when available.
- Keep focused inline records for code-table unit tests.
- Compare all claimed fields.
- Do not normalize away radical or stereo differences.

## Required test matrix

Table-driven parse/write/parse cases for:

- every supported atom-line and `M  CHG` charge;
- isotopes;
- atom maps;
- every radical variant;
- every supported single-bond stereo variant;
- supported double-bond stereo;
- aromatic/dative bonds if claimed;
- unsupported quadruple bonds;
- multiple M-record chunks;
- positive, negative, and zero coordinates;
- multiline SDF fields and multiple records.

Assert semantic equality, not only counts.

## Acceptance gates

- Radical variants round-trip exactly.
- Supported stereo round-trips exactly.
- No supported field changes silently.
- Unsupported cases return structured errors.
- Required V2000/SDF validation passes.
- Public atom-model changes are documented and re-exported correctly.



---

# Stage 5 — Correct aromatic SMILES and honest SMILES output

## Features

- `io.smiles.parse`
- `io.smiles.write`
- `algo.valence.rdkit-like`
- `algo.aromaticity.rdkit-like`
- `chem.sanitize.rdkit-like`

## Goal

Lowercase aromatic SMILES must parse, sanitize, perceive valence/aromaticity, write, and reparse with equivalent chemistry. The writer must never emit unreadable or silently lossy output.

## Required implementation

### 1. Distinguish omitted and explicit bond syntax

Replace the parser's default `BondOrder::Single` pending state with a representation that distinguishes:

- no bond symbol;
- explicit single;
- double;
- triple;
- aromatic;
- unsupported directional stereo.

When no bond symbol connects two aromatic atoms, resolve it as the documented aromatic connection. Apply equivalent logic to ring closures.

Do not treat explicit `-` as identical to an omitted aromatic bond.

### 2. Make ring labels symmetric

- Parse `%nn` labels emitted by the writer.
- Define supported range.
- Ensure the writer cannot exceed parser support.
- Reject ring closures across dot-separated components.
- Resolve bond specifications from both endpoints.
- Reject conflicting explicit specifications.

### 3. Strictly consume bracket syntax

Do not skip unknown characters.

Until implemented, reject with structured unsupported errors:

- `@` and `@@`;
- `/` and `\`;
- wildcard/query atoms;
- SMARTS decorators;
- unsupported atom classes.

A successful bracket parse must consume the full payload.

### 4. Correct pre-aromatic valence and aromaticity

Do not count every imported aromatic bond as a localized double bond.

Choose and document either:

- deterministic kekulization followed by perception; or
- an atom-contribution model that handles aromatic bond representation directly.

Required behavior:

- aromatic carbon has correct implicit H based on substitution;
- pyridine-like `n` has no implicit H;
- `[nH]` preserves explicit H and donor contribution;
- aromatic O/S/P behavior is explicit;
- fused aromatic systems work;
- invalid aromatic systems have an intentional error/non-aromatic outcome;
- successful sanitization produces consistent atom and bond aromatic flags.

### 5. Make the writer non-lossy

The writer must:

- reject `Zero`, `Dative`, and `Quadruple` bonds until faithfully supported;
- reject unsupported atom/bond stereo instead of dropping it;
- emit self-readable ring labels;
- preserve aromatic representation;
- remain deterministic for the same graph order;
- guarantee every successful output reparses with this library.

Do not map unsupported bonds to `-`.

### 6. Add integrated external validation

Add validation for SMILES parse followed by sanitization, not only raw `sanitize=False` import.

Compare:

- graph;
- charge;
- isotope;
- explicit H;
- atom map;
- implicit H after sanitization;
- atom/bond aromatic flags;
- writer output reparsed by this library;
- writer output accepted by RDKit and chemically equivalent.

Exact noncanonical text equality is not a substitute for graph equivalence unless explicitly documented as the contract.

## Required regression set

Supported cases:

```text
c1ccccc1
n1ccccc1
[nH]1cccc1
c1ccoc1
c1ccsc1
c1ccc2ccccc2c1
Cc1ccccc1
c1ccccc1.CC
C%10CCCCC%10
```

Intentional error/unsupported cases:

```text
[C@H](F)Cl
F/C=C/F
C-1CCCCC-1
c1cc1
C1.CCCC1
```

Every case must have a tested, documented outcome.

## Round-trip properties

For every supported representative input:

1. Parse succeeds.
2. Sanitization succeeds when the pinned reference accepts it.
3. Write succeeds.
4. Reparse succeeds.
5. Normalized graph, charge, explicit H, maps, implicit H, and aromaticity are equivalent.
6. Repeated write is deterministic.

## Acceptance gates

- `c1ccccc1` sanitizes as benzene with one implicit H per unsubstituted carbon.
- Pyridine-like and pyrrole-like N match reference behavior.
- Fused aromatic examples match goldens.
- `%` labels round-trip.
- Unsupported stereo/query/bond types return errors.
- Every successful writer output reparses.
- All applicable SMILES, valence, aromaticity, and sanitize validations pass.



---

# Stage 6 — Correct mmCIF residue identity and coordinate preservation

## Features

- `io.mmcif.parse`
- `bio.hierarchy.smcra`
- `core.conformers`

## Goal

Distinct atom-site residues must remain distinct when label sequence IDs are absent, and Cartesian coordinates must be retained.

## Required implementation

### 1. Define a residue-key type

Replace the ad hoc tuple with an internal key policy:

1. When `label_seq_id` exists, use label identity including model, label chain, component, label sequence ID, and insertion code.
2. Otherwise, when `auth_seq_id` exists, use author identity including model, author-or-label chain, component, author sequence ID, and insertion code.
3. When both are absent:
   - strict mode returns an ambiguity error;
   - lenient mode uses a documented conservative occurrence strategy that does not merge unrelated residues.

Preserve both label and author identifiers in hierarchy records.

### 2. Preserve Cartesian coordinates

Parse:

- `_atom_site.Cartn_x`
- `_atom_site.Cartn_y`
- `_atom_site.Cartn_z`

Requirements:

- all three present creates `Point3`;
- all three missing leaves no point;
- partial triplet errors in strict mode;
- invalid numeric input returns a line-specific error;
- define and document multiple-model-to-conformer behavior;
- do not infer bonds or sanitize.

### 3. Correct validation output

Compare actual coordinates rather than unconditional `null`.

Include enough residue identity in expected JSON to prove repeated waters/ligands are not merged.

Do not normalize away meaningful hierarchy identity or order.

### 4. Add external PDB corpus cases

Required:

- two `HOH` residues with absent label sequence IDs and distinct author IDs;
- repeated ligand names with distinct author sequence IDs;
- insertion-code variants;
- alternate locations;
- multiple models;
- complete coordinates;
- missing coordinates;
- malformed partial coordinate triplet.

## Required regression tests

- distinct waters remain separate residues;
- same-residue atoms still group;
- author sequence IDs are preserved;
- insertion codes affect identity;
- coordinates match;
- partial coordinates fail in strict mode;
- raw parse still has zero inferred bonds and absent chemistry perception.

## Acceptance gates

- No distinct regression residue is merged.
- No single residue is accidentally split.
- Coordinates are available through conformer APIs.
- Tiny/PDB validation passes where data is available.
- Strict and lenient ambiguity behavior is documented.



---

# Stage 7 — Bounded ring perception and stack-safe graph traversal

## Features

- `algo.rings.fast`
- `algo.rings.sssr`
- `algo.aromaticity.rdkit-like`
- `chem.sanitize.rdkit-like`
- `io.smiles.write`

## Goal

Ring perception and graph serialization must have documented, bounded behavior on highly symmetric and very large graphs.

## Required implementation

### 1. Instrument first

Measure:

- candidate cycles;
- equivalent shortest paths;
- path expansions;
- queue/stack peaks;
- graph size;
- total ring work;
- behavior on required corpora.

Use measurements to choose limits. Do not introduce unexplained constants.

### 2. Replace or bound all-shortest-path materialization

Replace the strategy that recursively materializes every shortest path for every ring bond.

Acceptable approaches include:

- deterministic bounded minimum-cycle-basis candidates;
- Horton-style candidates;
- another documented polynomial/bounded method preserving required behavior.

When work exceeds limits, return a structured resource-limit error. Do not return a misleading partial ring set.

### 3. Add options and errors

Expose documented limits if needed:

- maximum candidate cycles;
- maximum path expansions;
- maximum cycle size;
- maximum total work.

Propagate ring failures through aromaticity and transactional sanitization.

### 4. Remove recursion hazards

Use explicit heap-backed stacks for graph-size-dependent traversal:

- bridge/ring-membership DFS;
- shortest-path reconstruction;
- SMILES tree collection;
- subtree-size calculation;
- SMILES component emission.

Alternatively enforce a documented safe depth and return an error, but iterative traversal is preferred.

### 5. Add adversarial tests

Build deterministic:

- long chains;
- ladder graphs;
- theta graphs;
- fused/bridged polycycles;
- symmetric cages;
- disconnected mixtures.

Assert correct output or a specific resource-limit error. Prefer work-counter assertions to flaky timing tests.

### 6. Preserve chemical validation

Run ring, aromaticity, and sanitize validation across all available required corpora. Algorithmic speedups are not accepted if they silently change chemical behavior.

## Acceptance gates

- Long chains do not overflow the stack.
- Symmetric adversarial graphs terminate within limits.
- Resource errors leave sanitize input unchanged.
- Required ring/aromaticity goldens still pass.
- Defaults and rationale are documented.
- Ordinary corpus performance does not materially regress.



---

# Stage 8 — Behavior-preserving modularization

## Goal

Split the monolithic library and `xtask` sources into focused modules without changing behavior or public API.

Do this only after Stages 1–7 pass.

## Suggested library layout

```text
crates/molecules/src/
  lib.rs
  core/
    mod.rs
    ids.rs
    element.rs
    props.rs
    atom.rs
    bond.rs
    conformer.rs
    molecule.rs
  algorithms/
    mod.rs
    rings.rs
    valence.rs
    aromaticity.rs
  chemistry/
    mod.rs
    sanitize.rs
  io/
    mod.rs
    mol_v2000.rs
    sdf_v2000.rs
    smiles.rs
    mmcif.rs
  bio/
    mod.rs
    hierarchy.rs
```

## Suggested `xtask` layout

```text
crates/xtask/src/
  main.rs
  cli.rs
  dashboard.rs
  corpus.rs
  validation/
    mod.rs
    manifest.rs
    evidence.rs
    compare.rs
    normalize.rs
    status.rs
```

## Requirements

- Preserve public prelude and intended public paths.
- Use `pub(crate)` rather than widening visibility.
- Move tests next to implementations or focused integration files.
- Keep shared test builders test-only.
- Do not combine behavior changes with moves.
- Do not regenerate chemistry goldens except evidence hashes caused by source movement.
- Do not increment feature versions solely for file moves.
- Keep unsafe code forbidden.
- Keep rustdoc links and examples valid.

## Additional verification

```bash
cargo check --workspace --all-targets
cargo test --workspace --doc
```

Review every prelude re-export before and after.

## Acceptance gates

- No intentional behavior or API change.
- All tests and corpus comparisons pass.
- No feature version bump solely for refactoring.
- No new miscellaneous monolith.
- Module boundaries match `ARCHITECTURE.md`.

## Suggested commits

1. Extract core model.
2. Extract algorithms and sanitization.
3. Extract I/O and bio hierarchy.
4. Extract `xtask` validation/corpus/dashboard modules.
5. Relocate tests and clean re-exports.

Every commit should compile and pass tests.



---

# Stage 9 — Final closure and release readiness

## Goal

Make documentation, CI policy, and release metadata match the corrected implementation. Documentation cannot substitute for unfinished code.

## Required work

### 1. Refresh README

Update it to distinguish:

- implemented capabilities;
- validated capabilities and required corpora;
- unsupported capabilities;
- raw parsing versus sanitization;
- stability expectations for version `0.0.0`;
- local test and validation commands.

Do not describe the repository as only a minimal scaffold after implementation exists.

### 2. Document robustness posture

Document:

- malformed input returns structured errors;
- ring/resource limits;
- fuzz commands;
- crash artifact handling;
- security reporting process if public.

### 3. Enforce repository policy

Owner-verifiable GitHub settings should include:

- required Rust CI;
- required tiny reference validation;
- review requirements for validation generators, locks, goldens, and workflows;
- scheduled long fuzz and large-corpus jobs;
- branch protection appropriate to release policy.

Codex may document required settings but must not claim they are enabled without verification.

### 4. License decision

The owner must select a license. Codex must not choose one autonomously.

After an explicit owner decision:

- add canonical license text;
- update Cargo metadata and README;
- review compatibility of adapted material;
- document redistribution terms for corpus data and goldens.

Public release remains blocked until this is complete.

### 5. Repeat the audit

Re-audit:

- evidence integrity;
- CI execution;
- parser panic safety;
- V2000 round trips;
- aromatic SMILES;
- mmCIF identity/coordinates;
- mutation/cache invariants;
- ring bounds;
- module layout;
- documentation and licensing.

Open explicit issues for intentionally deferred work.

## Final acceptance gates

- Stages 1–8 are complete.
- Full workspace checks pass.
- Tiny validation passes in ordinary CI.
- Large required corpora pass in the designated workflow or have explicit blocking issues.
- No stale evidence is accepted.
- README claims match dashboard evidence.
- License and redistribution decisions are complete before public release.
- An independent audit finds no unresolved original High-severity finding.



---

# Pull request template for each stage

```markdown
## Feature IDs
- `<feature.id>`

## Audit findings closed
- Finding and stage reference

## Behavioral changes
- ...

## Regression tests added
- ...

## Reference validation
- Feature/corpus commands run
- Fixture and comparison counts
- Whether goldens changed and why
- Reference tool/version used for regenerated goldens

## Commands run
- [ ] cargo fmt --all -- --check
- [ ] cargo clippy --workspace --all-targets -- -D warnings
- [ ] cargo test --workspace
- [ ] RUSTDOCFLAGS="-D warnings" cargo doc --workspace --no-deps
- [ ] cargo xtask dashboard --check
- [ ] cargo xtask skills --check
- [ ] cargo xtask corpus check --corpus tiny --require-data
- [ ] cargo xtask validate --feature all --corpus tiny
- [ ] Stage-specific large-corpus checks

## Commands not run
- None, or exact command and reason

## Remaining limitations
- ...
```

# Global definition of done

The audit is closed only when all of these are true:

1. Validation evidence is tied to every material input and rechecked in CI.
2. Failed current validation cannot inherit old green status.
3. Public parsers are panic-free under regression and fuzz testing.
4. Supported V2000 radical, stereo, and metadata representations round-trip exactly.
5. Aromatic SMILES sanitize with reference-compatible valence and aromaticity.
6. SMILES output is self-readable and never silently lossy.
7. mmCIF residue identity and coordinates are preserved.
8. Mutation and sanitization cannot expose falsely fresh or partially committed chemistry state.
9. Ring and graph algorithms have bounded, documented failure behavior.
10. Source modules are architecture-aligned and reviewable.
11. CI, documentation, dashboard, and validation evidence agree.
12. The owner has completed license and release-policy decisions.
