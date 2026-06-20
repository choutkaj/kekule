# Molecules Audit Remediation Plan

This file converts the project audit into staged Codex `/goal` work.

- **Audit baseline:** `main` at `323ee995d8b2feb7a5422808e4ad2d4c30d8a5b3`
- **Status:** not started
- **Execution model:** one numbered stage per branch and pull request
- **Ordering rule:** complete every prerequisite stage before beginning the next

Stage 1 is intentionally first: later chemistry fixes cannot be trusted until validation evidence is current, content-addressed, and executed in CI.

## Working rules

For each stage:

1. Read `AGENTS.md`, `ARCHITECTURE.md`, this file, and every affected feature directory.
2. Create one branch and one PR for the stage.
3. List every affected feature ID in the PR.
4. Add a regression test that demonstrates the audited defect.
5. Implement the smallest coherent fix satisfying every acceptance gate.
6. Update feature versions only for intentional behavior, API, or validation-contract changes.
7. Regenerate evidence and dashboards only through repository commands.
8. Run all stage checks and report every command not run.
9. Perform an independent `feature-review` pass before handoff.
10. Mark a stage complete only after merge.

Do not weaken comparisons, remove asserted fields, delete regression tests, or regenerate goldens merely to obtain a pass. Keep RDKit and Biopython out of Rust runtime dependencies. Preserve raw parsing versus sanitization boundaries. Unsupported chemistry must return a structured error rather than be silently converted.

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

When metadata changes:

```bash
cargo xtask dashboard
cargo xtask dashboard --check
```

When non-`tiny` corpus data is available:

```bash
cargo xtask corpus check --corpus all --require-data
cargo xtask validate --feature all --corpus all
```

## Finding closure map

| Finding | Closing stage |
|---|---:|
| Validation evidence can remain green after source, fixture, lock, golden, or generator changes | 1 |
| CI does not run reference validation and large corpora are not reproducible from a clean runner | 1 |
| Perception caches can remain accessible or falsely fresh after mutation | 2 |
| Sanitization can partially mutate a molecule before returning an error | 2 |
| Malformed V2000 and other parser inputs can panic | 3 |
| V2000 stereo and radical data are not round-trip safe | 4 |
| Lowercase aromatic SMILES does not sanitize correctly | 5 |
| SMILES output can be unreadable by the parser or chemically lossy | 5 |
| mmCIF can merge distinct residues and drops Cartesian coordinates | 6 |
| Ring perception and recursive traversals have resource hazards | 7 |
| Monolithic source files impede safe review and maintenance | 8 |
| Documentation, repository policy, and licensing remain release blockers | 9 |

## Stage checklist

- [ ] Stage 1 — Validation evidence and reproducible CI
- [ ] Stage 2 — Mutation, cache invalidation, and transactional sanitization
- [ ] Stage 3 — Panic-free parsers and fuzzing
- [ ] Stage 4 — V2000/SDF semantic round-trip fidelity
- [ ] Stage 5 — Aromatic SMILES and honest SMILES output
- [ ] Stage 6 — mmCIF residue identity and coordinates
- [ ] Stage 7 — Bounded ring perception and stack-safe traversal
- [ ] Stage 8 — Behavior-preserving modularization
- [ ] Stage 9 — Final audit closure and release readiness

---

# Stage 1 — Trustworthy validation evidence and reproducible CI

## Primary feature

- `validation.harness`

Supporting feature metadata may change when validation contracts are migrated.

## Objective

A green feature/corpus cell must prove that the current implementation, manifest, source lock, fixtures, goldens, comparison code, reference generator, and reference environment all participated in a successful non-empty comparison. Old evidence must become invalid when any material input changes.

## Required implementation

### Content-addressed evidence

Introduce a versioned evidence schema. A conservative whole-workspace implementation digest is acceptable initially and is safer than an incomplete source-to-feature map.

Each pass must bind at least:

- feature/corpus manifest;
- corpus `sources.lock.json`;
- every listed fixture;
- every consumed golden;
- Rust implementation, comparison, and normalization source;
- Cargo manifests and `Cargo.lock`;
- affected feature metadata;
- applicable reference generator source;
- reference environment definition or lock;
- comparison mode and evidence schema version.

Hash deterministic relative paths and bytes in sorted order. Reject missing inputs. Store component hashes or a canonical evidence document plus final SHA-256. Recompute the evidence in `corpus_passed_at`; unknown schemas and mismatches are not passes. Preserve timestamps when an unchanged `--update` repeats identical valid evidence.

### Golden metadata checks

Verify:

- `feature_id`;
- `corpus_id`;
- `fixture_path`;
- supported golden schema version;
- `reference.tool` and `reference.version` against the manifest;
- `input_sha256` against current fixture bytes;
- generator/environment digest when present.

Do not trust copied status metadata without checking the golden document.

### Comparison contract

- Parse and enforce `comparison_mode`.
- Reject unknown modes.
- Implemented features with required corpora must not use zero fixtures.
- A pass requires `fixture_count > 0`.
- A pass requires `compared_count == fixture_count`.
- Missing goldens, skipped fixtures, and partial comparison are failures.
- Structured parse errors may be compared only when both implementations emit the same explicit record; omission is not a pass.

### Failure-safe update

For selected targets:

1. Load existing status.
2. Mark selected entries pending in memory.
3. Run every comparison.
4. Insert pass evidence only for fully successful targets.
5. Remove or explicitly fail selected targets that do not pass.
6. Write status files, synchronized feature flags, and dashboard through temporary files and rename.
7. Return nonzero when any selected target fails.

A failed current run must not inherit an old green result.

### CI

Add to ordinary PR CI:

```bash
cargo xtask corpus check --corpus tiny --require-data
cargo xtask validate --feature all --corpus tiny
```

The comparison must execute; stale status files alone are insufficient.

Make manual/scheduled large-corpus validation reproducible from a clean checkout. It must install the pinned reference environment, build or retrieve data from pinned locks, verify it with `--require-data`, run selected comparisons, upload concise failure logs, avoid committing data/status automatically, and use least-privilege permissions. A new corpus-preparation command must use typed configuration rather than executing unchecked shell text from repository metadata.

## Required regression tests

Add tests equivalent to:

- evidence invalid after Rust source change;
- evidence invalid after fixture change;
- evidence invalid after source-lock change;
- evidence invalid after golden change;
- evidence invalid after reference-generator change;
- fixture hash must match golden metadata;
- reference metadata must match manifest;
- required manifest rejects zero fixtures;
- pass requires comparison for every fixture;
- unknown comparison mode is rejected;
- failed update clears an old selected pass;
- unchanged update is idempotent;
- status/dashboard writes remain synchronized.

Use temporary directories and mutate real bytes.

## Acceptance gates

- One changed byte in implementation, fixture, golden, lock, or generator invalidates old evidence.
- A failed `--update` cannot leave the selected target green.
- Every applicable implemented `tiny` target performs a nonzero comparison.
- Normal CI runs `tiny` reference validation.
- A clean manual runner can prepare and validate one non-`tiny` corpus.
- No chemistry expected payload changes except a schema migration that reproduces equivalent normalized output.

## `/goal`

```text
/goal Implement Stage 1 from fixes.md. Make validation evidence content-addressed over implementation source, manifests, locks, fixtures, goldens, reference generators, and environment metadata; verify golden input/reference metadata; reject empty or partial comparisons; ensure failed --update runs invalidate old selected passes; make generated status/dashboard updates consistent and idempotent; run tiny implementation-versus-golden validation in normal CI; and make manual large-corpus validation reproducible from a clean checkout. Do not change chemistry behavior or weaken comparisons. Complete every Stage 1 acceptance gate.
```

---

# Stage 2 — Mutation, cache invalidation, and transactional sanitization

## Features

- `core.graph`
- `chem.sanitize.rdkit-like`
- `algo.valence.rdkit-like`
- `algo.rings.fast`
- `algo.rings.sssr`
- `algo.aromaticity.rdkit-like`

## Objective

No mutation may leave computed state falsely fresh or expose a cached ring object as current. Failed sanitization must leave the input molecule unchanged.

## Required implementation

- Centralize topology and chemistry invalidation.
- Topology changes invalidate valence, rings, aromaticity, and stereo.
- Charge, explicit-H, radical, element, aromatic representation, and bond-order changes invalidate affected state.
- Ring invalidation clears both cached membership and ring set.
- Coordinate-only and property-only changes must not claim topology changed.
- Remove or restrict public APIs that let callers mark perception fresh without running the algorithm.
- Make `sanitize_small_molecule` operate on a clone or staged change set and commit only after every requested pass succeeds.
- Define every `SanitizeOptions` combination: skipped passes are not fresh; cleanup invalidates affected state; requested passes become fresh only on success; repeated success is idempotent.
- Preserve raw parser behavior; do not introduce implicit sanitization.

## Required tests

- Mutation makes ring cache access return no current value.
- Atom and bond chemistry edits cannot leave relevant state fresh.
- Cleanup invalidates pre-existing perception.
- Every meaningful option combination has the expected final states.
- Valence and aromaticity failures leave the original molecule exactly equal to its pre-call snapshot.
- Successful sanitization is idempotent.
- Raw parsing still leaves perception absent.

## Acceptance gates

- No stale ring cache is returned after mutation.
- No unrequested state is marked fresh.
- Every sanitize error preserves the original molecule.
- Tiny chemistry comparisons still pass.
- Parser/sanitization separation remains explicit.

## `/goal`

```text
/goal Implement Stage 2 from fixes.md. Centralize chemistry/topology invalidation, clear stale ring caches, prevent arbitrary public marking of perception as fresh, and make sanitize_small_molecule atomic so any error leaves the input unchanged. Define and test every SanitizeOptions state transition while preserving raw parser behavior. Complete every Stage 2 acceptance gate.
```

---

# Stage 3 — Panic-free parsers and fuzzing

## Features

- `io.mol.v2000.parse`
- `io.sdf.v2000.parse`
- `io.smiles.parse`
- `io.mmcif.parse`

## Objective

Every public text parser returns a structured error rather than panicking for arbitrary Rust `&str`, including malformed records, non-ASCII text, zero indices, truncated fixed-width fields, and extreme counts.

## Required implementation

### V2000/SDF hardening

- Replace endpoint subtraction with `checked_sub`.
- Use checked arithmetic for block offsets and counts.
- Parse fixed-width structural fields through ASCII byte helpers.
- Never slice UTF-8 strings with unchecked byte ranges.
- Validate declared counts before allocation/indexing.
- Validate M-record counts, pair lengths, and indices before mutation.

### SMILES hardening

- Keep cursor positions on character boundaries.
- Fully validate bracket contents rather than skipping unknown bytes.
- Reject malformed charges, isotopes, maps, empty brackets, truncated element symbols, pending bonds without endpoints, and malformed ring syntax.
- Defer aromatic semantic work to Stage 5, but replace malformed-input acceptance with explicit errors.

### mmCIF hardening

- Preserve line reporting for unterminated quotes and semicolon values.
- Reject ragged loops and numeric overflow without panic.
- Document any input limits through parse options.

### Fuzzing

Add non-runtime fuzz targets for Molfile V2000, SDF V2000, SMILES, and mmCIF. Objectives: no panic/abort, bounded allocation from small input, and no panic on successful parse followed by applicable write/read operations. Seed with tiny corpus records and focused malformed regressions. Add a bounded deterministic CI smoke job; leave longer campaigns scheduled/manual.

## Required cases

- V2000 endpoint zero and endpoint beyond atom count;
- truncated/non-ASCII counts, atom, and bond lines;
- overflowing/inconsistent counts;
- truncated M records;
- unmatched SMILES branches/rings and malformed `%` labels;
- malformed/non-ASCII bracket atoms;
- malformed mmCIF quotes, semicolon blocks, ragged loops, and numeric overflow.

## Acceptance gates

- Every regression returns `Err` with useful location information.
- Bounded fuzz smoke runs without crashes.
- Unsafe code remains forbidden.
- Fuzz-only dependencies do not enter the runtime library graph.
- Tiny parser validation passes without dropping normalized fields.

## `/goal`

```text
/goal Implement Stage 3 from fixes.md. Make Molfile V2000, SDF V2000, SMILES, and mmCIF parsers panic-free for arbitrary &str using checked arithmetic and byte-safe parsing; add malformed-input regressions and non-runtime fuzz targets with bounded CI smoke coverage. Do not change supported chemistry semantics beyond replacing panics or silent malformed-input acceptance with structured errors. Complete every Stage 3 acceptance gate.
```

---

# Stage 4 — V2000 and SDF semantic round-trip fidelity

## Features

- `core.atom-bond`
- `io.mol.v2000.parse`
- `io.mol.v2000.write`
- `io.sdf.v2000.parse`
- `io.sdf.v2000.write`

## Objective

Every V2000 feature claimed as supported must survive parse/write/parse without changing meaning. Unsupported representations must fail explicitly.

## Required implementation

### Radical model

Replace the lossy electron-count field with one authoritative representation capable of distinguishing the supported MDL radical multiplicities, for example `Singlet`, `Doublet`, and `Triplet`. Map supported `M  RAD` codes according to a pinned CTfile specification and make the writer mapping the exact inverse. Provide an algorithm helper for unpaired-electron count. Do not keep two unsynchronized sources of truth.

### Bond stereo

Implement and document a parser/writer mapping table from the pinned specification. Wedge-up, wedge-down, and supported unknown/either forms must round-trip. Do not conflate double-bond stereo codes with wedge stereo. Unsupported combinations return a structured write error.

### Metadata symmetry

Verify formal charge, isotope, atom map, radical, coordinates, supported bond order, supported stereo, SDF fields, and header lines. Validate the complete record before returning output so failures do not expose partial text.

## Required tests

Create a table-driven suite covering every supported charge code, isotope/map, radical variant, stereo variant, supported bond order, unsupported quadruple bond, multiple M-record chunks, signed coordinates, multiline SDF fields, and multiple records. Assert semantic equality, not only counts.

Use external golden records where available and focused inline records for exact code-table regressions. Do not remove stereo/radical fields from normalized comparisons.

## Acceptance gates

- Radical variants round-trip exactly.
- Supported stereo variants round-trip exactly.
- No claimed field changes silently.
- Unsupported representations return errors.
- Required available V2000/SDF corpus comparisons pass.
- Public atom-model changes are documented and re-exported deliberately.

## `/goal`

```text
/goal Implement Stage 4 from fixes.md. Introduce a lossless authoritative atom-radical representation, make V2000 M RAD parsing and writing exact inverses, implement the pinned V2000 bond-stereo mapping without conflating wedge and double-bond semantics, and add complete parse/write/parse tests for every claimed metadata field. Unsupported representations must return structured errors. Complete every Stage 4 acceptance gate.
```

---

# Stage 5 — Correct aromatic SMILES and honest SMILES output

## Features

- `io.smiles.parse`
- `io.smiles.write`
- `algo.valence.rdkit-like`
- `algo.aromaticity.rdkit-like`
- `chem.sanitize.rdkit-like`

## Objective

Lowercase aromatic SMILES must parse, sanitize, perceive valence/aromaticity, write, and reparse with equivalent chemical meaning. Successful writer output must be readable by this parser and must not silently lose unsupported chemistry.

## Required implementation

- Represent “no bond symbol” separately from explicit single bond.
- Resolve omitted bonds between aromatic atoms, including ring closures, as aromatic according to the documented subset.
- Parse the same multi-digit ring-label syntax the writer emits; define a supported range and reject conflicts or cross-component closure.
- Strictly consume bracket syntax. Until implemented, reject `@`, `@@`, `/`, `\`, wildcard/query atoms, SMARTS decorators, and unsupported classes rather than ignoring them.
- Correct pre-aromatic valence and electron handling through deterministic kekulization or a sound atom-contribution model; do not count each aromatic bond as a localized double bond.
- Cover aromatic carbon substitution, pyridine-like `n`, `[nH]`, O/S/P donors, invalid systems, and fused rings.
- Make the writer reject zero, dative, quadruple, and unsupported stereo representations until it can encode them faithfully.
- Guarantee every successful output reparses and is deterministic for the same graph order.

Extend validation to exercise parse followed by sanitization and write/reparse. Compare graph, charge, isotope, explicit/implicit H, maps, aromatic flags, and chemical equivalence. Exact equality to one RDKit noncanonical traversal is not a substitute for graph equivalence unless explicitly part of the contract.

## Required regression set

Supported examples must include:

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

Unsupported/invalid examples must have intentional errors, including atom chirality, directional double-bond notation, conflicting ring bonds, an invalid aromatic ring, and a ring closure crossing a dot component.

## Acceptance gates

- Benzene from lowercase aromatic SMILES has one implicit H on each unsubstituted carbon after sanitize.
- Pyridine-like and pyrrole-like nitrogen behavior matches the pinned reference.
- Fused aromatic examples match required goldens.
- Multi-digit ring labels round-trip.
- Unsupported stereo/query/bond types return errors rather than lose information.
- Every successful writer output is accepted by this parser.
- Applicable available SMILES, valence, aromaticity, and sanitize corpora pass.

## `/goal`

```text
/goal Implement Stage 5 from fixes.md. Correct omitted aromatic-bond and ring-closure parsing, support writer-compatible multi-digit ring labels, strictly reject unsupported bracket/stereo/query syntax, make imported aromatic SMILES sanitize with correct implicit hydrogens and fused-ring aromaticity, and make the writer reject every representation it cannot encode faithfully. Add integrated RDKit-backed parse-sanitize-write-reparse validation and complete every Stage 5 acceptance gate.
```

---

# Stage 6 — Correct mmCIF residue identity and coordinate preservation

## Features

- `io.mmcif.parse`
- `bio.hierarchy.smcra`
- `core.conformers`

## Objective

Distinct residues must remain distinct when label sequence IDs are absent, and parsed Cartesian coordinates must be retained.

## Required implementation

Introduce an internal residue-key type.

Recommended policy:

1. When `label_seq_id` exists, key by model, label chain, component, label sequence ID, and insertion code.
2. Otherwise, when `auth_seq_id` exists, key by model, author-or-label chain, component, author sequence ID, and insertion code.
3. When both are absent, strict mode returns an ambiguity error; lenient mode uses a documented conservative occurrence strategy that does not merge unrelated residues.

Always preserve both label and author identifiers in the hierarchy.

Parse `_atom_site.Cartn_x`, `_atom_site.Cartn_y`, and `_atom_site.Cartn_z`. All three present values create a `Point3`; all missing leaves no coordinate; a partial triplet is an error in strict mode. Define and document how models map to conformers. Do not infer bonds or sanitize.

Update Rust/Biopython normalized output to compare coordinates and enough residue identity to prove repeated waters/ligands are not merged.

## Required cases

- two `HOH` residues with absent label sequence and distinct author sequence IDs;
- repeated ligand names;
- insertion-code variants;
- alternate locations;
- multiple models;
- complete, missing, malformed, and partial coordinates.

Use pinned PDB-derived mmCIF records for goldens.

## Acceptance gates

- Distinct residues are not merged.
- A single residue is not accidentally split.
- Coordinates are available through conformer APIs.
- Raw parse still has no inferred bonds or chemistry perception.
- Tiny and available PDB corpus comparisons pass.
- Strict/lenient ambiguity behavior is documented.

## `/goal`

```text
/goal Implement Stage 6 from fixes.md. Replace the mmCIF residue tuple with a documented label/auth identity policy that never merges distinct residues when label_seq_id is absent, preserve both identifier families, parse Cartn_x/y/z into conformer storage with strict partial-triplet errors, and extend Biopython-backed validation to compare residue identity and coordinates. Do not infer bonds or sanitize. Complete every Stage 6 acceptance gate.
```

---

# Stage 7 — Bounded ring perception and stack-safe graph traversal

## Features

- `algo.rings.fast`
- `algo.rings.sssr`
- `algo.aromaticity.rdkit-like`
- `chem.sanitize.rdkit-like`
- `io.smiles.write`

## Objective

Ring perception and graph serialization must have documented bounded behavior on highly symmetric and very large graphs.

## Required implementation

- Instrument candidate cycles, equivalent shortest paths, path expansions, queue/stack peaks, graph size, total work, and corpus behavior before choosing limits.
- Replace or strictly bound recursive all-shortest-path materialization. A deterministic bounded minimum-cycle-basis or Horton-style candidate approach is acceptable if required behavior is preserved.
- Return a structured resource-limit error rather than partial misleading rings or unbounded allocation.
- Expose documented limits where needed: candidates, path expansions, cycle size, and total work.
- Propagate ring errors through aromaticity and transactional sanitization.
- Convert graph-size-dependent recursion to explicit stacks for bridge DFS, shortest-path reconstruction, SMILES tree collection, subtree sizing, and component emission.

Add deterministic adversarial graph tests: long chain, ladder, theta graphs, fused/bridged polycycles, symmetric cages, and disconnected mixtures. Assert output or a specific limit error using work counters rather than flaky timing.

## Acceptance gates

- Long chains do not overflow the call stack.
- Symmetric adversarial graphs terminate within configured limits.
- Resource errors leave sanitize input unchanged.
- Required ring/aromaticity goldens still pass.
- Defaults and their rationale are documented.
- Ordinary corpus performance does not materially regress.

## `/goal`

```text
/goal Implement Stage 7 from fixes.md. Instrument ring work, replace or strictly bound all-shortest-path candidate enumeration, add structured ring resource-limit errors and options, propagate failures transactionally through aromaticity and sanitization, and convert graph-size-dependent recursive traversals to iterative forms. Preserve required RDKit ring/aromaticity behavior and complete every Stage 7 acceptance gate.
```

---

# Stage 8 — Behavior-preserving modularization

## Objective

Split the monolithic library and `xtask` sources into focused architecture-aligned modules without changing behavior or public API. Begin only after Stages 1–7 pass.

## Suggested library layout

```text
crates/molecules/src/
  lib.rs
  core/{mod,ids,element,props,atom,bond,conformer,molecule}.rs
  algorithms/{mod,rings,valence,aromaticity}.rs
  chemistry/{mod,sanitize}.rs
  io/{mod,mol_v2000,sdf_v2000,smiles,mmcif}.rs
  bio/{mod,hierarchy}.rs
```

## Suggested `xtask` layout

```text
crates/xtask/src/
  main.rs
  cli.rs
  dashboard.rs
  corpus.rs
  validation/{mod,manifest,evidence,compare,normalize,status}.rs
```

## Requirements

- Preserve public prelude and intended public paths.
- Use `pub(crate)` rather than widening visibility.
- Move tests beside implementation or into focused integration files.
- Keep shared test builders test-only.
- Do not combine behavior changes with moves.
- Do not regenerate chemistry goldens except evidence hashes caused by source movement.
- Do not increment feature versions solely for moving code.
- Keep unsafe code forbidden and rustdoc links valid.

Additional verification:

```bash
cargo check --workspace --all-targets
cargo test --workspace --doc
```

## Acceptance gates

- No intentional behavior or API change.
- All tests and corpus comparisons pass.
- No version bump solely for refactoring.
- No new miscellaneous monolith.
- Module boundaries match `ARCHITECTURE.md`.

## `/goal`

```text
/goal Implement Stage 8 from fixes.md as a behavior-preserving refactor. Split crates/molecules/src/lib.rs and crates/xtask/src/main.rs into architecture-aligned modules, preserve public exports and prelude paths, keep internals crate-private, relocate tests coherently, and prove all unit, doc, dashboard, skills, corpus, and validation checks still pass. Do not mix functional fixes or feature-version bumps into this stage.
```

---

# Stage 9 — Final closure and release readiness

## Objective

Make documentation, CI policy, and release metadata match the corrected implementation. Documentation cannot substitute for unfinished code.

## Required work

- Refresh README to distinguish implemented, validated, unsupported, raw-parse, and sanitize behavior and to state `0.0.0` stability expectations.
- Document parser error behavior, ring/resource limits, fuzz commands, crash artifacts, and a security-reporting path if public.
- Verify owner-managed repository policy: required Rust CI, required tiny validation, review controls for generators/locks/goldens/workflows, and scheduled long fuzz/large-corpus jobs. Codex must not claim these settings without inspection.
- Owner selects the license. Codex must not choose one autonomously. After an explicit decision, add canonical text, Cargo/README metadata, compatibility review, and corpus/golden redistribution terms.
- Repeat the original audit and open explicit issues for intentional deferrals.

## Acceptance gates

- Stages 1–8 are complete.
- Full workspace checks pass.
- Tiny validation passes in ordinary CI.
- Large required corpora pass in the designated workflow or have explicit blockers.
- No stale validation evidence is accepted.
- README claims match dashboard evidence.
- License and redistribution decisions are complete before public release.
- An independent audit finds no unresolved original High-severity finding.

## `/goal`

```text
/goal Implement the code-and-documentation portions of Stage 9 from fixes.md after Stages 1-8 are complete. Refresh README and robustness documentation so claims match generated evidence, add final audit checks, and document required repository settings. Do not select a license or claim branch-protection settings without an explicit owner decision and verification. Complete every non-owner-blocked Stage 9 acceptance gate.
```

---

# Pull-request template for each stage

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
- Reference tool/version for regenerated goldens

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

The audit is closed only when:

1. Validation evidence is tied to every material input and rechecked in CI.
2. Failed current validation cannot inherit old green status.
3. Public parsers are panic-free under regression and fuzz testing.
4. Supported V2000 radical, stereo, and metadata representations round-trip exactly.
5. Aromatic SMILES sanitize with reference-compatible valence and aromaticity.
6. SMILES output is self-readable and never silently lossy.
7. mmCIF residue identity and coordinates are preserved.
8. Mutation and sanitization cannot expose falsely fresh or partially committed state.
9. Ring and graph algorithms have bounded documented failure behavior.
10. Source modules are architecture-aligned and reviewable.
11. CI, documentation, dashboard, and validation evidence agree.
12. The owner has completed license and release-policy decisions.
