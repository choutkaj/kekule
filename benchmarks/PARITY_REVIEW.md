# Benchmark parity review — 2026-09-17

Work proceeds one feature per turn on `codex/benchmark-parity`.
The baseline is the completed local run `run-1789652366683-27088.json`
(started 2026-09-17), with 25 features across five datasets. Its recorded revision
is `5c3324698daff5b69a076d6e3de431966f4fe336` with a dirty worktree; this review
started from clean `main` at `1fa0fe15`. Source locks remain unchanged; independently
justified reference corrections are documented per feature below.
Raw cases and analysis stay local under ignored `target/benchmark-parity/`.
The baseline case archive is 78.6 GB; `case-index.json` records byte ranges by
feature/dataset so later reviews need not repeatedly scan unrelated features.
See the [consolidated handoff](#consolidated-handoff) for the remaining issues
and the completion evidence across all features.

## Feature progress

| Feature | State |
| --- | --- |
| `io.smiles.parse` | Reviewed; adapter fix verified; residuals below |
| `io.smiles.write` | Reviewed; setup/diagnostic fixes verified on full corpus |
| `io.smiles.canonical` | Reviewed; full corpus agrees where supported; no new fix needed |
| `io.smiles.isomeric` | Reviewed; full corpus agrees where supported; no new fix needed |
| `io.mol.parse` | Reviewed; missing reference SDF metadata restored; stereo/policy residuals below |
| `io.mol.v2000.write` | Reviewed; complete-record reader validation added; residual stereo/H-policy differences |
| `io.mol.v3000.write` | Reviewed; coordinate precision and fixed-H encoding corrected; query-output validation added |
| `io.sdf.parse` | Reviewed; SDF reader routing corrected and private reference fields restored; chemistry residuals documented |
| `io.sdf.v2000.write` | Reviewed; record/field validation repaired and private reference fields restored; residual chemistry/format limits documented |
| `io.mmcif.parse` | Reviewed; embedded quote parsing fixed; all corpus differences are permitted multiline whitespace policy |
| `algo.rings.fast` | Reviewed; corrected reference API; all membership observations agree, with CXSMILES input limits |
| `algo.rings.sssr` | Reviewed; cycle-path comparison strengthened; one bond-order-sensitive alternative remains |
| `algo.valence.rdkit-like` | Reviewed; reference cleanup corrected; all supported corpus observations agree |
| `algo.aromaticity.rdkit-like` | Reviewed; one SMILES radical-inference gap remains; explicit-radical SDF agrees |
| `algo.canonical-ranking` | Reviewed; hydrogen-equivalence invariant corrected; one upstream radical-inference gap remains |
| `algo.substructure.vf2` | Reviewed; SMARTS bond typing and explicit ring-bond syntax corrected; one upstream radical-inference gap remains |
| `query.smarts` | Reviewed; prior ring-bond fix resolves 37,030 parser failures; remaining 20,904 errors require stereo queries |
| `chem.perception.default` | Reviewed; four prior aromatic-triple fixes verified; hydrogen/radical policies and source-stereo differences remain |
| `chem.hydrogen-transforms` | Reviewed; generated-H declarations corrected in 187,615 cases; removal/source policies and stereo residuals remain |
| `descriptor.molecular` | Reviewed; atomic-data generator path repaired; formulas match throughout, mass constants/conventions remain distinct |
| `descriptor.rotatable-bonds.rdkit-strict` | Reviewed; charged resonance exclusions corrected; hydrogen-representation differences and 86 aromaticity-approximation cases remain |
| `stereo.representation` | Reviewed; 2,881 private reference fields restored; full stereo/source audit completed and broader interpretation limits documented |
| `stereo.perception` | Reviewed; repeated-H candidates and overcoordinated coordinate inference repaired; reference metadata restored; broader symmetry/stereo-family limits remain |
| `stereo.cip` | Reviewed; benchmark pseudoasymmetric descriptor spelling aligned; remaining source-stereo differences and explicit ranking limits documented |
| `bio.secondary-structure.dssp` | Reviewed; float precision, nullable-label validation and input batching repaired; complete independent-reference audit verified |

## 1. SMILES parsing

The adapter previously replaced every aromatic bond type with `AROMATIC`, losing
triple bond order. It now applies that representation only to single/double
bonds, retaining higher bond orders and the separate aromatic flag. This follows
RDKit's [markAtomsBondsArom implementation](https://github.com/rdkit/rdkit/blob/Release_2026_03/Code/GraphMol/Aromaticity.cpp).
Rust and independent RDKit regressions cover arynes and assert that changing the
reported triple bond to aromatic still fails comparison. No fields, tolerances,
reference values or fixture memberships were removed or weakened.

Full rerun: `cargo benchmark --feature io.smiles.parse --dataset all --output target/benchmark-parity/io-smiles-parse.json`.
The run completed and correctly exited 1 for the remaining differences.

| Dataset | Cases | Agree before → after | Disagree after | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 100,000 | 92,557 → 92,559 | 7,432 | 9 |
| Enamine | 50,240 | 37,065 → 37,065 | 10,294 | 2,881 |
| Smoke | 8 | 8 → 8 | 0 | 0 |

PL-Rex and PDB have no supplied SMILES; their 164 and 1,000 selected IDs remain
not applicable, as do 12 smoke IDs. The two corrected external cases are PubChem
141133 and 141134. No runtime code changed.

### Unresolved distinctions and defects

- **Hydrogen representation:** 6,520 PubChem cases and all 10,294 Enamine
  disagreements differ in explicit H, inferred H and represented valence.
  Both engines give identical per-atom explicit-plus-implicit H totals in all
  these cases (437,070 atom comparisons). RDKit sanitization moves some inferred
  aromatic hydrogens into its explicit count; Kekule preserves declarations.
  The guide now explains why these native-API comparisons are not a pure
  chemical-equivalence score. Raw differences remain asserted.
- **Radicals:** 933 PubChem cases contain radical differences (911 radical-only,
  21 overlapping hydrogen differences, and one overlapping aromaticity).
  Kekule deliberately does not infer radical state from bracket atoms, as tested
  by `bracket_atoms_do_not_infer_radicals_from_a_valence_model`. RDKit uses
  electron-deficit rules, including parity conventions for isolated metals.
  A trial interpretation change exposed conflicts with the existing SMILES
  writer and round-trip contracts and was fully reverted. Defer a coordinated
  parser/writer policy change; do not manufacture radical values in the adapter.
  PubChem 181201 also differs in aromaticity at a bracket carbon radical.
- **Stereo cleanup:** PubChem 39368 retains represented double-bond stereo whose
  nitrogen endpoint has two equivalent oxygen ligands; RDKit cleanup removes
  it. Kekule's represented-stereo policy preserves source assertions. Leave for
  the stereo feature review rather than silently cleaning the observation.
- **CXSMILES:** all 2,881 Enamine errors explicitly reject source extensions.
  Supporting these requires grammar and represented-stereo/group integration;
  discarding extensions would create false agreements.
- **Invalid valence:** the same nine PubChem cases fail in both engines. Matching
  failures remain errors, never agreements.

### Validation

All final checks passed:

- `cargo fmt --all -- --check`
- `cargo check --workspace --all-targets --all-features --locked --offline`
- `cargo clippy --workspace --all-targets --all-features --locked --offline -- -D warnings`
- `cargo test --workspace --all-features --locked --offline` (includes doctests)
- `cargo test -p kekule-bench --locked --offline` (42 unit and four integration tests)
- `cargo doc --workspace --all-features --no-deps --locked --offline`, with `RUSTDOCFLAGS=-D warnings`
- `cargo package -p kekule-bench --locked --allow-dirty --offline --list`
- Pinned RDKit Python: `python -m unittest discover -s benchmarks/reference/rdkit -p 'test_*.py'` (17 tests)
- `git diff --check`

Logs are `target/benchmark-parity/check-*.log`, `bench-tests.log` and the benchmark
report/log. An earlier trial radical implementation failed existing regressions;
all final workspace tests passed after reverting the trial.

Not run for this benchmark-only change: runtime crate package builds/file lists,
MSRV and no-default-feature matrices, fuzz builds/runs, other reference suites
and external trajectory checks. Their code/configuration is unchanged. Full
packaging of the unpublished benchmark is replaced by the package file-list
check, consistent with the existing benchmark workflow.

The original next feature was `io.smiles.write`. Every writer case in the baseline errored;
inspect the recorded writer-validation cause and reference Python setup first.
The installed pinned interpreter is
`C:/Users/chout/AppData/Roaming/mamba/envs/molecular-rdkit-reference/python.exe`.
Its `Library/bin` directory may need to precede PATH. Plain `python` is absent
from this shell's PATH.

## 2. Plain SMILES writing

All 150,248 applicable baseline cases errored because the selected interpreter
could not import `rdkit.Chem`. This was reference setup failure, not evidence that
all molecules failed to serialize. The original emitted outcomes remain in the
baseline archive's `written` field.

With the pinned RDKit 2026.03.3 environment, both complete pre-fix and final reruns produced:

| Dataset | Applicable cases | Agree | Disagree | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 100,000 | 87,455 | 0 | 12,545 |
| Enamine | 50,240 | 41,874 | 0 | 8,366 |
| Smoke | 8 | 6 | 0 | 2 |

The 1,176 selected IDs without supplied SMILES remain not applicable.
All 129,335 successful outputs preserve the reference's complete canonical
isomeric CXSMILES identity. There is no observed identity disagreement to repair
in the plain writer. This is coverage of successfully emitted cases, not a claim
that unsupported inputs are handled.

Remaining errors are 18,023 explicit stereo rejections (12,536 PubChem, 5,485
Enamine, two smoke), 2,881 unsupported Enamine CXSMILES extensions, and nine
PubChem valence errors also present in the reference. The ordinary writer
intentionally rejects represented stereo; separate isomeric and canonical
features exercise stereo-preserving modes. Do not strip stereo from expectations
or switch modes under this feature's name to hide those errors.

### Benchmark corrections

- Before the first applicable input for each writer feature/dataset, verify the
  selected interpreter can load the reference and matches the stored version.
  Setup failure produces an incomplete run with a clear error and no fabricated
  per-molecule measurements. Datasets without applicable writer inputs do not
  require an interpreter.
- Preserve original implementation errors through both a failed reader process
  and individual reader results. Independent validation cannot replace a writer
  error or turn failed emission into successful output.
- Per-case diagnostics distinguish `kekule_failed` and
  `writer_validation_failed`, matching the already separate summary counters.
  The previous per-case flag incorrectly attributed reader failures to Kekule.

Regression coverage includes unavailable-reader CLI behavior and mixed batches
containing both a successful emission and an implementation error, with reader
process failure and per-case reader failure. Existing version-drift validation
continues to pass. No runtime code, source data, goldens or scientific comparison
fields/tolerances were changed.

### Validation

The final implementation passed formatting, workspace all-target/all-feature
check and clippy with warnings denied, workspace tests including doctests,
Rustdoc with warnings denied, and the benchmark package file-list check. Exact
commands match section 1; logs are in
`target/benchmark-parity/smiles-write-checks/`. The focused benchmark suite has
43 unit and five integration tests. The shared Python reference runner's eight
tests, dashboard's 21 Python tests and JavaScript live-update checks also passed.
Runtime packaging, MSRV/no-default-feature matrices, fuzzing, Biopython/DSSP and
trajectory suites were not rerun because their implementations/configurations
are unchanged. The RDKit adapter was unchanged this turn; its 17 tests passed in
the preceding feature review.

Pre-fix report: `target/benchmark-parity/io-smiles-write-before.json`.
The final full rerun completed (`complete: true`) and correctly exited 1 for
unsupported inputs and invalid valence. Report:
`target/benchmark-parity/io-smiles-write.json`, using:

```text
cargo benchmark --feature io.smiles.write --dataset all --writer-python C:/Users/chout/AppData/Roaming/mamba/envs/molecular-rdkit-reference/python.exe --output target/benchmark-parity/io-smiles-write.json
```

A record-by-record comparison of all 151,424 rows verified unchanged source
identities, emitted text, reference observations and scientific outcomes. All
20,913 original implementation errors are preserved with correct per-case flags.
Evidence: `target/benchmark-parity/smiles-write-verification.json`.

Next feature: `io.smiles.canonical`.

## 3. Canonical SMILES writing

The baseline's reference interpreter failed, as for plain SMILES writing. Its
saved emitted outcomes nevertheless showed that all 147,358 supported inputs
completed canonical writing, reverse-numbering invariance and read/write
fixed-point checks. Re-ran the entire feature with the pinned RDKit interpreter
to independently validate the emitted molecular identities:

```text
cargo benchmark --feature io.smiles.canonical --dataset all --writer-python C:/Users/chout/AppData/Roaming/mamba/envs/molecular-rdkit-reference/python.exe --output target/benchmark-parity/io-smiles-canonical-before.json
```

| Dataset | Applicable cases | Agree | Disagree | Errors |
| --- | ---: | ---: | ---: | ---: |
| pubchem-100k | 100,000 | 99,991 | 0 | 9 |
| enamine-diversity | 50,240 | 47,359 | 0 | 2,881 |
| smoke | 8 | 8 | 0 | 0 |

The run completed and correctly exited 1 for 2,881 unsupported Enamine CXSMILES
inputs and nine PubChem valence errors that also fail in RDKit. The 1,176 IDs
without supplied SMILES remain not applicable. No independent reader failures,
observation-schema errors, canonicalization contract failures or identity
disagreements occurred. No additional writer or adapter fix was warranted.

The benchmark compares complete canonical isomeric CXSMILES identity after RDKit
independently reads emitted bytes, including isotope, charge, radical, atom-map
and stereo information. It does not require Kekule's traversal string to match
RDKit's. The separate Kekule fixed-point and reverse-numbering probes are active;
these are useful but do not exhaust all atom/bond permutations. Existing unit
regressions additionally cover permutations, stereo carriers and bounded labeling.
No comparison fields or assertions were weakened, and no goldens were changed.

All 151,424 report rows were checked against the summary counts. All successful
observations equal the stored expectations exactly; all 2,890 implementation
errors retain their original messages and correct attribution. Local evidence:
`target/benchmark-parity/canonical-verification.json` and
`io-smiles-canonical-before-analysis.json`.

This turn changed only this review log. Formatting and `git diff --check` passed.
Rust check, clippy, tests, Rustdoc, package checks, reference-adapter tests and
other feature/fuzz/MSRV matrices were not repeated: no runtime, adapter, schema,
test or build-configuration code changed. The executable, contract and reference
adapter hashes exactly match the validated final plain-writer run from section
2; its successful workspace and benchmark checks apply to this same code.

Next feature: `io.smiles.isomeric`.

## 4. Isomeric SMILES writing

The original run's universal reader errors concealed successful emission for
147,358 cases. A full rerun with the pinned RDKit reference environment completed:

```text
cargo benchmark --feature io.smiles.isomeric --dataset all --writer-python C:/Users/chout/AppData/Roaming/mamba/envs/molecular-rdkit-reference/python.exe --output target/benchmark-parity/io-smiles-isomeric.json
```

| Dataset | Applicable cases | Agree | Disagree | Errors |
| --- | ---: | ---: | ---: | ---: |
| pubchem-100k | 100,000 | 99,991 | 0 | 9 |
| enamine-diversity | 50,240 | 47,359 | 0 | 2,881 |
| smoke | 8 | 8 | 0 | 0 |

All 147,358 successful outputs match the stored complete canonical isomeric
CXSMILES identity exactly when independently read by RDKit. All 18,023 inputs
rejected by the plain writer for represented stereochemistry pass here; a
case-by-case comparison with the plain-writer run verified their identities and
coverage. There are no identity disagreements, reader failures or schema errors.

The 2,881 Enamine CXSMILES inputs remain explicitly unsupported, and the same
nine PubChem valence errors occur in both engines. The 1,176 IDs without supplied
SMILES remain not applicable. The completed run correctly returns exit 1 for
these errors; they are not counted as agreement.

The adapter selects the actual public isomeric mode, writes each component and
passes only emitted bytes to the independent reader. The runtime emits source
order with localized bonds and adjusts stereo orientation to its emitted
neighbor order. The benchmark appropriately requires molecular identity without
imposing canonical traversal invariance on this noncanonical writer. Full
isomeric identity remains asserted; enhanced groups and other unsupported stereo
representations are not silently discarded. No additional writer or adapter fix
was justified by this feature's supplied corpus.

All 151,424 report rows were checked against the summary, including exact
successful expectations and error attribution. Local evidence:
`target/benchmark-parity/isomeric-verification.json` and
`io-smiles-isomeric-analysis.json`. No code, goldens, fixtures, schemas or
comparison thresholds changed this turn; only this review log was extended.

Formatting and `git diff --check` passed. Rust check, clippy, tests, Rustdoc,
package checks, Python reference tests and other feature/MSRV/fuzz matrices were
not repeated because their code and configuration remain unchanged. The tested
executable, reference-adapter and contract hashes match the final plain-writer
run in section 2, whose workspace and benchmark checks passed.

Next feature: `io.mol.parse`.

## 5. MOL parsing

The reference adapter silently omitted private SDF fields: RDKit's default
`GetPropNames()` excludes names beginning with an underscore. The adapter now
uses the ordered source headers to select values from RDKit. This preserves
source fields without accidentally including internal RDKit properties. A
regression covers private/public ordering, multiline values, exclusion of
internal metadata and duplicate private names, through both MOL and SDF parsing.

All 50,240 Enamine reference records were independently regenerated with the
pinned RDKit 2026.03.3 environment into a new directory. An exhaustive old/new
comparison proved that the only changes restore `_CXSMILES_Data` in 2,881
records. All chemical observations, previously asserted fields, source digests
and case identities are unchanged. The original archive/manifest remain under
`target/benchmark-parity/mol-original-goldens/`; the corrected pair was promoted
only after that audit. The manifest records the actual adapter fingerprint.
No contract, tolerance, source input, comparison field or assertion was weakened.
Other feature archives retain their existing provenance and await their reviews.

Commands:

```text
cargo benchmark --feature io.mol.parse --dataset all --output target/benchmark-parity/io-mol-parse-before.json
cargo benchmark generate --feature io.mol.parse --dataset enamine-diversity --python C:/Users/chout/AppData/Roaming/mamba/envs/molecular-rdkit-reference/python.exe --goldens target/benchmark-parity/mol-corrected-goldens --output target/benchmark-parity/mol-reference-generation.json
cargo benchmark --feature io.mol.parse --dataset all --output target/benchmark-parity/io-mol-parse-after.json
```

| Dataset | Applicable cases | Agree in baseline → after | Disagree after | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 100,000 | 93,949 → 93,951 | 6,040 | 9 |
| Enamine | 50,240 | 32,625 → 32,634 | 17,606 | 0 |
| PL-Rex | 328 | 290 → 290 | 38 | 0 |
| Smoke | 17 | 12 → 12 | 5 | 0 |

The two PubChem improvements come from section 1's aromatic triple-bond adapter
fix. Restoring metadata removes a false field mismatch in all 2,881 affected
Enamine records; nine now agree completely, while 2,872 retain independent
hydrogen/stereo distinctions. All 1,000 PDB IDs and three smoke IDs remain not
applicable. The completed run correctly exits 1 for residual disagreements and
errors. No runtime chemistry code changed in this feature review.

### Unresolved distinctions and defects

- Hydrogen counts: every affected atom in the MOL observations has matching
  total hydrogens (25,803 Enamine atoms and seven smoke atoms). RDKit moves
  aromatic N hydrogens into the explicit count. Kekule resolves omitted wedge
  hydrogens into `Fixed(1)`, whereas RDKit keeps explicit H with inference enabled.
  These are representation and future-inference-policy distinctions, not evidence
  of lost H atoms. Their individual fields remain asserted.
- Stereo cleanup: many extra elements involve equivalent substituents; PubChem
  32 has a double bond with equivalent ester arms, and PL-Rex JAK1/4E5W includes
  a center with equivalent ring paths. Matching RDKit's sanitization requires
  a general symmetry-aware cleanup policy, including stereogenic dependencies.
  Source interpretation must not silently install general perception. Deferred.
- Crowded drawing orientation: 58 PubChem cases differ solely in parity, with
  further parity differences among cases whose stereo arrays differ structurally.
  PubChem 10524 reproduces inversion at two tetrahedral centers. The pinned
  [RDKit wedge interpreter](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/Chirality.cpp)
  handles angular ordering, collinear bonds and conflicting pseudo-3D wedges;
  Kekule uses displaced-coordinate volume. A reliable replacement needs those
  drawing ambiguities addressed together, rather than a parity flip. Deferred.
- Three PubChem records are assigned an axis instead of a pyramidal tetrahedral
  element; 146091 is a sulfoxide. Kekule approximates SP2 eligibility with ring
  membership or a double bond. The pinned [RDKit axis detector](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/Atropisomers.cpp)
  explicitly uses conjugation-aware hybridization to avoid this defect. A shared
  chemically complete endpoint predicate is needed; an element-specific ban
  would be insufficient. Deferred.
- Five PubChem records retain radical differences. For example, 24755 supplies
  divalent carbon valence declarations without RAD records; RDKit infers two
  radical electrons at each deficient carbon, while Kekule preserves the source
  declaration. This requires the explicit radical-inference policy discussed in
  section 1, not treating every valence deficit as a radical.
- Across all substantive stereo differences, the audit also records 64 missing
  PubChem double-bond elements. These remain unresolved along with the extra and
  differently oriented elements; no claim of complete stereo parity is made.

Coordinate deviations in these rows are within the existing source-precision
threshold. The analysis excludes tolerated differences when categorizing defects.
Evidence: `mol-substantive-audit.json`, `mol-inspection.txt`,
`mol-golden-verification.json` and the original externally supplied record excerpts
under `target/benchmark-parity/`.

The final before/after audit checked all 151,588 report rows and matched their
counts to the completed summary. Every actual Kekule observation is identical;
only the independently verified reference metadata changed. Exactly nine rows
move from disagreement to agreement, with no regressions or altered errors.
Evidence: `target/benchmark-parity/mol-rerun-verification.json`.

### Validation

Passed: workspace formatting; all-target/all-feature check and clippy with
warnings denied; all workspace tests and doctests; workspace Rustdoc with
warnings denied; benchmark package file-list check; 18 pinned RDKit tests;
eight shared reference-runner tests; 21 dashboard Python tests and the Node
dashboard regression; `git diff --check`. Rust commands used `--locked --offline`
where applicable. Logs are under `target/benchmark-parity/mol-checks/`.

Runtime package builds, the MSRV/optional-feature matrices, Linux fuzzing,
Biopython/DSSP-specific tests and unrelated optional external benchmarks were
not rerun: this change affects only the RDKit reference adapter and its metadata
observations, with no Rust implementation or build-configuration changes. The
unpublished benchmark package uses the documented file-list packaging check.

Next feature: `io.mol.v2000.write`.

## 6. V2000 MOL writing

The baseline classified every applicable case as an error because its Python
reader could not import RDKit. With the pinned reader, the full feature exposes
real writer behavior. The existing preflight and error-attribution fixes from
section 2 apply here as well.

The independent reader also accepted incomplete output coverage: appending an
SDF separator, an SDF field, or a second complete molecule to one `.mol` payload
still returned only the first molecule. RDKit stops at `M  END`; successful parsing
alone did not prove that the entire emitted file was valid. The adapter now
requires each emitted MOL file to contain one complete record, with only
whitespace after its terminator. RDKit still independently validates the chemistry.

The new regression failed in all six malformed-output cases before the fix and
passes afterwards: all three trailing-content forms are tested for both MOL
versions, with trailing whitespace accepted. This strengthens output validation
without changing any molecular observations, goldens, source fixtures or tolerances.

Full before/after commands:

```text
cargo benchmark --feature io.mol.v2000.write --dataset all --writer-python C:/Users/chout/AppData/Roaming/mamba/envs/molecular-rdkit-reference/python.exe --output target/benchmark-parity/io-mol-v2000-before.json
cargo benchmark --feature io.mol.v2000.write --dataset all --writer-python C:/Users/chout/AppData/Roaming/mamba/envs/molecular-rdkit-reference/python.exe --output target/benchmark-parity/io-mol-v2000-after.json
```

| Dataset | Applicable cases | Agree | Disagree | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 100,000 | 99,833 | 158 | 9 |
| Enamine | 50,240 | 41,987 | 5,372 | 2,881 |
| PL-Rex | 328 | 328 | 0 | 0 |
| Smoke | 17 | 13 | 4 | 0 |

Both full runs completed with these counts and correctly exited 1 for the
remaining discrepancies. All 151,588 case rows are byte-for-byte identical
before and after stricter output validation, including actual observations,
emitted bytes, expectations, differences and error attribution. Their summary
counts were also cross-checked. All 1,000 PDB IDs and three smoke IDs remain
not applicable. Evidence: `target/benchmark-parity/v2000-rerun-verification.json`.

### Remaining behavior

- PubChem: 140 cases differ in tetrahedral parity; 18 gain a tetrahedral
  assignment absent from the reference. Tracing each failing focus to the native
  MOL interpretation identifies 68 cases with inherited orientation defects and
  73 with orientation defects introduced during wedge projection; one case
  belongs to both categories. PubChem 30535 demonstrates a correct internal
  orientation being emitted incorrectly. All 18 added assignments are already
  present after interpretation (example 159021).
- The wedge interpreter and writer's self-check use the same geometric rule,
  so a successful self-round-trip does not establish external stereo correctness.
  A focused bond-stretching reproducer demonstrates an orientation flip without
  changing angular order. Simply normalizing bond vectors fixes 10524 but fails
  another independently checked RDKit drawing. No partial geometry fix was kept.
  The reproducer and probe remain local in `wedge-length-reproducer.rs` and
  `probe_wedge.py`; a coherent angular/ambiguous-wedge implementation is deferred.
- Enamine: 5,371 cases differ only in `no_implicit_hydrogens`; one additional
  case, Z4834848811, also differs in stereo parity. The writer faithfully encodes
  Kekule's fixed wedge-H declaration as a valence field, which RDKit reads with
  inference disabled. Counts, connectivity and other chemical fields match.
  Changing that policy belongs with the interpretation distinction in section 5.
- All 2,881 Enamine implementation errors explicitly reject enhanced stereo
  groups, which V2000 cannot express. These inputs remain counted, without
  silently stripping groups or emitting V3000 under a V2000 feature label.
- Four smoke cases differ only in the hydrogen-inference flag. All 328 PL-Rex
  cases agree. Nine PubChem source-valence failures also fail the independent
  emitted-output reader; they remain errors, not agreement.

The benchmark compares complete independently decoded indexed graphs and source
coordinates, with only the declared V2000 quantization allowance. It observes
the public model writer's empty-title/no-SDF-data contract, retaining unmodified
source observations in reports. Private SDF fields restored in section 5 do not
change MOL writer expectations because MOL has no data section.

Evidence: `io-mol-v2000-before-analysis.json`, `io-mol-v2000-before-examples.json`
and `v2000-stereo-trace.json` under `target/benchmark-parity/`.

### Validation

Passed: formatting; workspace all-target/all-feature check and clippy with
warnings denied; workspace tests and doctests; Rustdoc with warnings denied;
benchmark package file-list check; 19 pinned RDKit tests; eight shared
reference-runner tests; 21 dashboard Python tests and the Node dashboard
regression; `git diff --check`. Rust commands used `--locked --offline` where
applicable. Logs are under `target/benchmark-parity/v2000-checks/` and the
`v2000-reference-regression-before/after.log` files.

Runtime package builds, MSRV/optional-feature matrices, Linux fuzzing,
Biopython/DSSP tests and unrelated optional benchmarks were not rerun because
the retained change affects only independent RDKit output validation. The
temporary geometry regression was removed after its proposed fix failed
independent adjudication; no runtime code or build configuration changed.

Next feature: `io.mol.v3000.write`.

## 7. MOL V3000 writing

The original fresh run could not import its RDKit writer reader, as described
in section 2. A new full before-run with the pinned reader exposed two runtime
defects and one benchmark blind spot.

### Fixes

- Shared model serialization rounded coordinates to four decimal places before
  stereo projection, and the V3000 renderer also printed only four places.
  V3000 now preserves the converted floating-point coordinates using round-trip
  decimal formatting. V2000 retains its required fixed-width rounding.
  Stereo validation uses the coordinates actually emitted by the selected
  version. Auto fallback rebuilds the V3000 projection rather than reusing a
  rounded V2000 record or failing before it can try V3000.
- Positive fixed hydrogen declarations were emitted as `HCOUNT`, which RDKit
  reads as a query atom. The writer now uses total molecular `VAL`, including
  attached fixed hydrogens; zero valence remains `VAL=-1`. This follows the
  same general valence encoding used by the V2000 writer, with no molecule-
  specific exceptions.
- Independent MOL/SDF output validation now rejects any query atom or bond
  before sanitization. Its ordinary graph observations had omitted query
  constraints, allowing chemically different output to appear to agree.
  Input reference generation remains unchanged.

Focused Rust regressions cover high-precision coordinates through direct,
generic, sink, SDF and Auto model writers; an alkene whose drawing collapses
under four-place rounding; and fixed-H encoding for neutral, charged, aromatic
and stereochemical atoms. Independent RDKit regressions demonstrate that
`HCOUNT=4` creates a query while `VAL=4` preserves methane without a query,
and reject both query atoms and query bonds. The regressions failed before
their respective fixes. Two existing geometry tests now modify atom fields
without depending on an exact decimal spelling; all original assertions remain.

Full before/after commands:

```text
cargo benchmark --feature io.mol.v3000.write --dataset all --writer-python C:/Users/chout/AppData/Roaming/mamba/envs/molecular-rdkit-reference/python.exe --output target/benchmark-parity/io-mol-v3000-before.json
cargo benchmark --feature io.mol.v3000.write --dataset all --writer-python C:/Users/chout/AppData/Roaming/mamba/envs/molecular-rdkit-reference/python.exe --output target/benchmark-parity/io-mol-v3000-after.json
```

| Dataset | Applicable cases | Agree before → after | Disagree after | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 100,000 | 99,833 → 99,833 | 158 | 9 |
| Enamine | 50,240 | 47,358 → 41,999 | 8,241 | 0 |
| PL-Rex | 328 | 328 → 328 | 0 | 0 |
| Smoke | 17 | 15 → 13 | 4 | 0 |

Both full runs completed and correctly exited 1 for the remaining discrepancies.
All 1,000 PDB IDs and three smoke IDs remain not applicable. The lower agreement
count reflects the replacement of query encoding with molecular valence, which
exposes source-vs-interpreted hydrogen-inference policy differences. It is not
a reason to retain chemically incorrect query output or remove asserted fields.

The paired audit verifies all 151,588 rows, including unchanged source
expectations and case identities. All 2,883 coordinate-disagreement cases are
resolved: 2,881 Enamine and two smoke cases. The remaining atom-observation
changes are exclusively `no_implicit_hydrogens`; hydrogen counts, valence,
charges, radicals, bonds, stereo assertions and stereo groups are unchanged.
Small coordinate differences elsewhere reflect unit-conversion roundoff rather
than four-place rounding and pass the existing numerical comparison.

Previously, 8,245 outputs contained `HCOUNT` query atoms, including 5,375 cases
reported as agreement. No rerun output contains `HCOUNT`, and the independent
reader's stronger check finds no emitted query atoms or bonds. Twelve Enamine
cases and both precision-related smoke cases now fully agree. Other corrected
coordinate cases still expose the separate hydrogen-policy/stereo differences.

### Remaining behavior

- PubChem has the same 140 tetrahedral-parity and 18 added-assignment cases
  documented for V2000 in section 6. Nine independent source-valence failures
  also fail the emitted-output reader. No PubChem case changed status.
- All 8,241 Enamine disagreements include a hydrogen-inference flag difference;
  8,238 differ only in that flag, and three also have unchanged tetrahedral-
  parity differences. This is the same source-vs-represented fixed-H distinction
  described in sections 5–6. Correct molecular valence makes it visible instead
  of obscuring it through query semantics.
- Four smoke cases differ only in the hydrogen-inference flag. All PL-Rex cases
  agree. V3000 preserves the enhanced stereo groups that V2000 cannot represent;
  the paired audit finds no group changes from these fixes.

Evidence: `v3000-rerun-verification.json`, `verify_v3000.py`,
`io-mol-v3000-before-analysis.json` and `io-mol-v3000-after-examples.json` under
`target/benchmark-parity/`. No goldens, fixtures, comparison fields or tolerances
were changed for this feature.

### Validation

Passed on Windows: `cargo fmt --all -- --check`; workspace all-target/all-feature
check and clippy with warnings denied; workspace all-feature tests and explicit
doctests; Rustdoc with warnings denied; Rust 1.89 workspace all-target/all-feature
check; potentials tests and documentation without default features; full
`cargo package -p kekule --allow-dirty` including package verification; package
file-list checks for benchmarks, potentials and trajectory; license-copy hash
comparisons; 20 pinned RDKit tests; eight shared reference-runner tests; 21
dashboard Python tests and the Node regression; `git diff --check`.
Rust commands used `--locked --offline` where applicable. Logs are in
`target/benchmark-parity/v3000-checks/`; `tests-final.log` and
`clippy-final.log` supersede the intermediate logs with the old formatting-
dependent fixture failure. Before/after focused failures and passes are also
retained in `v3000-regression-*.log` and `v3000-query-regression-*.log`.

Linux fuzz-target build and nightly fuzz smoke were not run on this Windows
host. Full companion package builds were not run because they require a
published foundational crate; their file lists were checked as prescribed by
CI. Biopython/DSSP tests and unrelated external benchmarks were not rerun
because this feature changes molecular serialization and RDKit output checking.

Next feature: `io.sdf.parse`.

## 8. SDF parsing

The fresh before-run reproduces the original SDF disagreements except for the
two PubChem aromatic-triple observations already corrected in section 1. Its
complete disagreement patterns match the MOL interpretation review in section
5, including the missing private Enamine data fields.

### Benchmark corrections

The SDF feature previously chose a reader solely from the filename suffix.
Both adapters therefore used their MOL parser on `.mol`/`.mdl` inputs instead
of exercising the SDF reader. Focused regressions demonstrate that a second
record and the first record's private field were silently omitted. Both
adapters now use their SDF reader for the SDF feature. Kekule uses its public
EOF-termination option for standalone MOL inputs and retains default delimiter
checks for `.sdf` inputs. The common SDF path uses each record's rich public
interpretation, and carries fields from the same parse instead of reparsing
the document just to recover metadata.

The new Rust and independent RDKit tests require both record titles, the full
record count, field placement and equal observations for SDF/MOL/MDL suffixes.
They fail before the routing fix and pass after it. The earlier private-field
regression also explicitly covers this feature. No runtime chemistry changed.

Fresh independent references are generated in a separate directory and checked
against every old outcome before promotion. The first full Enamine comparison
verified exactly 2,881 `_CXSMILES_Data` additions and no changes to any existing
field, chemical observation, case identity or source digest.

```text
cargo benchmark --feature io.sdf.parse --dataset all --output target/benchmark-parity/io-sdf-parse-before.json
cargo benchmark generate --feature io.sdf.parse --dataset all --python C:/Users/chout/AppData/Roaming/mamba/envs/molecular-rdkit-reference/python.exe --goldens target/benchmark-parity/sdf-final-goldens --output target/benchmark-parity/sdf-final-reference-generation.json
```

The complete fresh reference audit passed across all 151,588 rows. Outside the
2,881 Enamine field additions, every observation is unchanged, including all
nine PubChem reference errors. The new SDF-reader routing therefore preserves
the reference chemistry on the full supplied corpus. Original Enamine files
were backed up under `target/benchmark-parity/sdf-original-goldens/`, then only
that corrected archive and its generation manifest were promoted. The other
historical archives remain untouched; the fresh audit archives remain local.
Evidence: `sdf-all-golden-verification.json` and `audit_sdf_goldens.py`.

The final comparison uses the fully audited fresh references:

```text
cargo benchmark --feature io.sdf.parse --dataset all --goldens target/benchmark-parity/sdf-final-goldens --output target/benchmark-parity/io-sdf-parse-after.json
```

| Dataset | Applicable cases | Agree before → after | Disagree after | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 100,000 | 93,951 → 93,951 | 6,040 | 9 |
| Enamine | 50,240 | 32,625 → 32,634 | 17,606 | 0 |
| PL-Rex | 328 | 290 → 290 | 38 | 0 |
| Smoke | 17 | 12 → 12 | 5 | 0 |

The run completed and correctly exited 1 for the remaining disagreements and
errors. All 1,000 PDB IDs and three smoke IDs remain not applicable. Restoring
metadata removes a false field mismatch in every affected Enamine record; nine
now agree completely, while 2,872 retain separate chemistry/policy differences.

### Remaining behavior

The residuals match section 5's MOL interpretation issues. PubChem has 5,974
cases with different stereo arrays, 58 parity-only cases, three differing
axis/tetrahedral assignments and five radical-only cases. These require general
symmetry-aware stereo cleanup, robust drawing interpretation, conjugation-aware
axis eligibility or an explicit radical-inference policy. Nine valence failures
remain errors in both engines. Enamine and smoke primarily differ in declared
versus inferred hydrogen representation and inference policy, with some Enamine
stereo differences; PL-Rex's 38 disagreements concern stereo arrays. None is
addressed by molecule-specific exceptions or by removing asserted fields.

The paired audit checks all 151,588 rows against both the SDF before-run and
the reviewed MOL after-run. Every actual Kekule observation is unchanged, and
every final actual/expected outcome matches the reviewed MOL results. Only the
2,881 independently verified reference-field additions differ from the SDF
before-run; exactly nine cases improve, with no regressions or altered errors.
All affected hydrogen totals agree across engines (25,803 Enamine atoms and
seven smoke atoms), despite their declaration/inference differences. No
coordinate, connectivity, isotope, charge or additional chemistry discrepancy
was introduced. Evidence: `sdf-rerun-verification.json` and
`verify_sdf_rerun.py` under `target/benchmark-parity/`.

### Validation

Passed: Rust formatting; workspace all-target/all-feature check and clippy with
warnings denied; workspace all-feature tests including doctests; Rustdoc with
warnings denied; benchmark package file-list check; 21 pinned RDKit tests;
eight shared reference-runner tests; 21 dashboard Python tests and the Node
regression; `git diff --check`. Rust commands used `--locked --offline` where
applicable. Logs are in `target/benchmark-parity/sdf-checks/` and
`sdf-dispatch-regression-before.log`, `sdf-bench-regression-after.log`, and
`sdf-reference-regression-before/after.log`.

Runtime package builds, MSRV and optional-feature matrices were not rerun:
this turn changes benchmark adapters and reference data, with no runtime crate
or build-configuration changes. Linux fuzz checks were not run on this Windows
host. Biopython/DSSP and unrelated external benchmarks were not rerun because
their adapters and data are unaffected.

Next feature: `io.sdf.v2000.write`.

## 9. SDF V2000 writing

The original fresh run could not import its independent writer reader, as in
section 2. A new full before-run uses the pinned RDKit interpreter to expose
the actual molecular and metadata behavior.

### Benchmark corrections

The output-format check split text on every `$$$$` substring. This falsely
rejected valid titles, field names and values containing those characters,
while accepting an inline substring as the required final record delimiter.
Record splitting now recognizes complete delimiter lines, preserves blank
titles, checks every record's format and requires termination of every emitted
record. Trailing whitespace is accepted; empty records are rejected.

The source-field enumeration also scanned the entire file for header-shaped
text, including titles and multiline field values. It now scans only after
the CTAB terminator, consumes each value through its blank-line boundary, and
selects the resulting names from RDKit's native property values. Duplicate
actual field names still produce explicit reference errors. For emitted SDF,
the same scanner additionally rejects missing CTAB terminators and unexpected
text outside fields, which RDKit can warn about and silently ignore.

Three focused independent regressions failed before these fixes: literal
delimiter text in metadata, header-like titles/value lines, and malformed record
framing. They cover multiple records, blank titles, private fields, trailing
whitespace, missing terminators, junk outside fields and empty records. A Rust
regression confirms that the actual writer already preserves delimiter-like
metadata and record boundaries correctly. No runtime writer change was needed.

The SDF writer's Enamine references also retained the old private-field omission.
All 50,240 Enamine reference outcomes were independently regenerated, with no
errors. A full old/new audit proves exactly 2,881 `_CXSMILES_Data` additions,
with every previously asserted chemical value, data field and case identity
unchanged. A second generation using the corrected scanner produced a
byte-identical compressed archive. Only the audited Enamine pair was promoted,
with its final generation fingerprint; the original archive and manifest remain
under `target/benchmark-parity/sdf-writer-original-goldens/`. Evidence:
`sdf-writer-golden-verification.json` and
`sdf-writer-final-golden-verification.json`.

Full before/after commands:

```text
cargo benchmark --feature io.sdf.v2000.write --dataset all --writer-python C:/Users/chout/AppData/Roaming/mamba/envs/molecular-rdkit-reference/python.exe --output target/benchmark-parity/io-sdf-v2000-before.json
cargo benchmark generate --feature io.sdf.v2000.write --dataset enamine-diversity --python C:/Users/chout/AppData/Roaming/mamba/envs/molecular-rdkit-reference/python.exe --goldens target/benchmark-parity/sdf-writer-final-goldens --output target/benchmark-parity/sdf-writer-final-reference-generation.json
cargo benchmark --feature io.sdf.v2000.write --dataset all --writer-python C:/Users/chout/AppData/Roaming/mamba/envs/molecular-rdkit-reference/python.exe --output target/benchmark-parity/io-sdf-v2000-after.json
```

| Dataset | Applicable cases | Agree | Disagree | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 100,000 | 99,833 | 158 | 9 |
| Enamine | 50,240 | 41,987 | 5,372 | 2,881 |
| PL-Rex | 328 | 328 | 0 | 0 |
| Smoke | 17 | 13 | 4 | 0 |

Both full runs completed with these counts and correctly exited 1 for residual
disagreements and unsupported/error inputs. All 1,000 PDB IDs and three smoke
IDs remain not applicable.

An exhaustive cross-check of all 151,588 rows against the reviewed MOL V2000
run confirms identical molecular observations for all 147,695 decoded outputs,
identical errors for 2,890 applicable cases and identical non-applicability for
1,003 IDs. Every successful SDF output preserves its source title and ordered
data fields. Evidence: `sdf-mol-writer-equivalence.json` and
`compare_sdf_mol_writers.py`. The latter's final accounting keeps non-applicable
input diagnostics separate from actual errors.

### Remaining behavior

- PubChem has the same 140 parity differences and 18 added stereo assignments
  documented in section 6. Nine source-valence errors also fail independent
  output reading.
- Enamine has 5,371 hydrogen-inference-flag-only disagreements; Z4834848811 also
  has a parity disagreement. Its 2,881 enhanced-group inputs are explicitly
  rejected because V2000 cannot encode those groups. Restoring their reference
  metadata does not make unsupported serialization pass.
- Four smoke cases differ only in the hydrogen-inference flag. All PL-Rex
  records pass. No title, field-value, field-order or record-count discrepancy
  remains in successfully decoded corpus outputs.

The chemistry issues remain deferred for the reasons in sections 5–6. No
comparison fields, tolerances, source fixtures or molecular goldens were
weakened. The final paired audit checks all 151,588 rows and confirms identical
emitted text, independently read observations and statuses before and after the
validator changes. Every decoded output retains its source title, field values,
field order and record count. The only reference changes are the independently
verified private-field additions in 2,881 inputs already rejected for enhanced
stereo; no error was concealed or changed to agreement. Summary counts match
the complete case archive. Evidence: `sdf-writer-rerun-verification.json`,
`verify_sdf_writer_rerun.py` and `io-sdf-v2000-after-examples.json`.

### Validation

Passed: formatting; workspace all-target/all-feature check and clippy with
warnings denied; all-feature workspace tests including doctests; Rustdoc with
warnings denied; benchmark package file-list check. Rust commands used
`--locked --offline` where applicable. Logs are under
`target/benchmark-parity/sdf-writer-checks/`. The isolated candidate passes all
21 preexisting RDKit tests plus the three new framing regressions. All 24
tests also pass against the installed final adapter. Eight shared reference-
runner tests, 21 dashboard Python tests, the Node regression and
`git diff --check` pass. The before-fix failures and isolated candidate probes
are retained in the `sdf-framing-*.log` files.

Runtime package builds, MSRV and optional-feature matrices were not rerun
because this turn changes benchmark validation, reference data and focused
tests, with no runtime crate or build-configuration change. Linux fuzz checks
were not run on this Windows host. Biopython/DSSP and unrelated external
benchmarks were not rerun because their adapters and data are unaffected.

Next feature: `io.mmcif.parse`.

## 10. mmCIF parsing

The original PDB run's 79 disagreements involve exactly 103 values. Every
difference is trailing spaces or tabs on lines inside semicolon-delimited text;
block names, tag membership, column lengths and all other text agree. The
independent audit preserves the actual and reference strings in
`target/benchmark-parity/mmcif-baseline-audit.json`.

[CIF 1.1 syntax, paragraphs 17–20](https://www.iucr.org/what-we-do/digital-standards/cif/cif1/file-syntax)
permits eliding trailing whitespace in multiline text. Biopython does this;
Kekule preserves the source text. Both conform here. Keeping more source text
fits the document layer's preservation contract, so neither runtime trimming
nor comparison normalization is warranted. All strict mismatches remain visible;
no goldens, asserted fields or tolerances change.

### Parser correction

The review found a separate general lexical defect: a matching quote inside a
quoted CIF 1.1 value ended the token even when followed by a non-whitespace
character. For example, `'don't stop'` was rejected. The tokenizer now closes
such a token only before whitespace or end of line, following the CIF delimiter
rule and Biopython's implementation. This is a small change in the document
tokenizer; interpretation and chemical state are untouched.

Focused regressions cover both quote styles, adjacent embedded quotes, literal
hash characters, empty values, loop rows, space/tab/newline separators, and
unterminated values. The valid-quote regression fails before the fix. Existing
multiline preservation/limit tests remain intact. Independent Biopython tests
confirm quote behavior and its trailing-whitespace policy, while the benchmark
regression asserts that strict observations retain the entire source text.

### Benchmark scope and remaining limitations

This feature compares decoded values from every category, not biomolecular
interpretation or complete CIF conformance. `MMCIF2Dict` loses quoted-versus-bare
missing-token distinctions, accepts some malformed syntax (including duplicate
tags), and cannot faithfully represent multiple data blocks. The adapter's
post-parse multi-block guard is incomplete for blocks encountered inside a loop.
These remain reference limitations outside the supplied single-block corpus;
implementing an independent general CIF validator is not a small parity fix.
The guide now makes this boundary explicit. These source-level distinctions
require focused syntax tests rather than dictionary comparison.

### Validation

Passed: `cargo fmt --all -- --check`; workspace all-target/all-feature check
and clippy with warnings denied; workspace all-feature tests including doctests;
Rustdoc with warnings denied; Rust 1.89 workspace all-target/all-feature check;
potentials tests and documentation without default features; full
`cargo package -p kekule --allow-dirty` including package verification; package
file lists for benchmarks, potentials and trajectory; license-copy hash checks;
all eight pinned Biopython reference tests; and `git diff --check`.
Rust validation commands used `--locked --offline` where applicable.

The existing debug target directory initially linked an old library into the
benchmark regression even though the library's own new tests passed. A fresh
`target/mmcif-validation` directory resolved this. Workspace check, clippy,
tests, docs, optional-feature checks and the final release benchmark use that
fresh directory. `mmcif-*-fresh.log`, `mmcif-docs.log`, `mmcif-optional-*.log`
and `io-mmcif-parse-final.*` are authoritative; the redundant intermediate
`io-mmcif-parse-after` process was stopped and is not a completed validation.

Linux fuzz-target checks and nightly fuzz smoke were not run on this Windows
host. Full companion package verification was replaced by the CI-prescribed
file-list checks because it resolves the foundational crate through crates.io.
RDKit/shared-runner/dashboard tests and other feature benchmarks were not rerun:
no changes to those adapters, runner, comparison code or dashboard were made
for this feature. All corpus goldens and source locks remain unchanged.

Full rerun:
`cargo run --release --locked --offline --target-dir target/mmcif-validation -p kekule-bench -- --feature io.mmcif.parse --dataset all --output target/benchmark-parity/io-mmcif-parse-final.json`.
It completed and correctly exited 1 for the retained exact-text differences.

| Dataset | Applicable | Agree before → after | Disagree | Errors |
| --- | ---: | ---: | ---: | ---: |
| PDB | 1,000 | 921 → 921 | 79 | 0 |
| Smoke | 1 | 1 → 1 | 0 | 0 |
| Total | 1,001 | 922 → 922 | 79 | 0 |

Another 150,423 selected IDs have no mmCIF source and remain not applicable.
The valid supplied corpus did not exercise the embedded-quote defect; focused
regressions provide the before/after evidence for that fix.

The paired audit verifies all 151,424 selected rows are identical before and
after, including every actual/reference observation, status and difference.
Evidence: `target/benchmark-parity/mmcif-rerun-verification.json` and
`verify-mmcif-rerun.py`. The 79 retained discrepancies are therefore exactly
the fully audited whitespace-policy cases, with no new differences.

Next feature: `algo.rings.fast`.

## 11. Fast ring membership

The original run and a full current rerun have no membership disagreements.
The 2,881 errors are Enamine SMILES inputs with CXSMILES extensions, rejected
before ring perception rather than silently discarding their source semantics.
The corresponding SDF inputs are independently included and succeed.

The Rust implementation marks non-bridge edges and their endpoints using an
iterative linear-time traversal, excluding zero-order and dative bonds. It
handles all connected components of the eligible-edge graph. The adapter
observes every atom flag and every bond flag with its source endpoints;
comparison retains both membership and correspondence. No runtime fix was
supported by this review.

### Reference adapter correction

The fast-membership adapter invoked `Chem.GetSymmSSSR`, doing selected-ring
perception for a feature that only requires cycle membership. It now uses
`Chem.FastFindRings`, the API documented by the
[RDKit Book](https://github.com/rdkit/rdkit/blob/master/Docs/Book/RDKit_Book.rst)
for this purpose. The selected-ring feature keeps its independent `GetSymmSSSR`
call. This changes which reference operation is exercised, without weakening
the asserted observations.

The regression fails on the old adapter when selected-ring perception is
unavailable. With the corrected adapter, it verifies two rings joined by a
bridge (cyclic endpoints, noncyclic connecting bond) and cycles interrupted by
zero-order or dative bonds. All 25 pinned RDKit tests pass.

The corrected reference independently regenerated all five datasets. A full
old/new comparison verifies every one of the 301,834 golden rows is unchanged,
including all observations, input hashes, identities, reference versions and
missing-format entries. All five pairs were promoted with their actual
generation manifests; predecessors are backed up under
`target/benchmark-parity/rings-fast-original-goldens/`. Different compressed
bytes do not indicate changed observations. No fields or tolerances changed.

Full rerun against the corrected reference:
`target/mmcif-validation/release/kekule-bench.exe --feature algo.rings.fast --dataset all --goldens target/benchmark-parity/rings-fast-goldens --output target/benchmark-parity/algo-rings-fast-after.json`.
It completed and correctly exited 1 for the unsupported CXSMILES inputs.

| Dataset | Applicable | Agree before → after | Disagree | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 200,000 | 200,000 → 200,000 | 0 | 0 |
| Enamine | 100,480 | 97,599 → 97,599 | 0 | 2,881 |
| PL-Rex | 328 | 328 → 328 | 0 | 0 |
| Smoke | 25 | 25 → 25 | 0 | 0 |
| Total | 300,833 | 297,952 → 297,952 | 0 | 2,881 |

Both supplied molecular formats are evaluated separately. Another 1,001
selected cases have no applicable input (1,000 PDB and one smoke case).
The complete native error audit confirms all 2,881 errors have the explicit
CXSMILES rejection diagnostic. Adding full CXSMILES support remains outside
this scoped ring review; the matching SDF cases stay included.

Evidence: `rings-fast-golden-verification.json`, `verify-fast-ring-goldens.py`,
`rings-fast-generation.log` and `algo-rings-fast-after.*` in
`target/benchmark-parity/`. The paired native audit also verifies all 301,834
before/after rows are identical, including observations, errors and comparison
details; see `rings-fast-rerun-verification.json` and
`verify-fast-ring-rerun.py`.

### Validation

Passed: `cargo fmt --all -- --check`; workspace all-target/all-feature check
and clippy with warnings denied; workspace all-feature tests including doctests;
Rustdoc with warnings denied; benchmark package file-list check; 25 pinned RDKit
tests, eight shared reference-runner tests and `git diff --check`. Rust commands used
`--locked --offline` where applicable and `--target-dir target/mmcif-validation`
for compilation, following the clean-build validation in section 10. Logs are
`target/benchmark-parity/rings-fast-*.log`.

Runtime packaging, MSRV, optional-feature matrices, Linux fuzzing, Biopython/DSSP
and dashboard tests were not repeated: this feature changes one Python
reference operation, its test and documentation, with no Rust runtime,
runner, schema or dashboard change. The unpublished benchmark uses the
package file-list check, as in the preceding benchmark-only reviews.

Next feature: `algo.rings.sssr`.

## 12. Selected ring sets

The original and current full runs contain one disagreement: PubChem 296819,
SMILES record 620 in `data/packs/pack_061.smi`. Kekule selects five rings of
sizes 4, 6, 6, 6 and 7; RDKit selects those plus an alternative seven-membered
ring. The corresponding MOL source agrees. All other successful cases agree.

### Remaining selection limitation

The candidate search is sensitive to bond iteration order. An independent
experiment on the original supplied SMILES removes and reinserts one bond in
RDKit, changing its iteration position while preserving the complete indexed
chemical graph. RDKit 2026.03.3 then returns exactly Kekule's five cyclic paths.
The original RDKit ordering reproduces all six stored reference paths.
Both results contain valid cycles, cover all 27 cyclic bonds, and span the
graph's rank-five cycle space. Kekule still lacks the additional alternative
required for exact agreement with this reference ordering.

This is a candidate-selection limitation rather than a missing cyclic bond or
invalid cycle. The trace and the pinned
[RDKit implementation](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_3/Code/GraphMol/FindRings.cpp)
support the traversal-order explanation. No bounded change was established
that resolves it while retaining the current selected-ring policy. The runtime
algorithm and the visible disagreement remain unchanged. There is no
case-specific ring insertion, source rewrite or replacement of the comparison
with cycle-space equivalence.

A general pruning trial excluded noncyclic bonds before candidate search using
the existing membership flags. It fixed this case and passed the ring regressions,
but a full corpus rerun introduced eight new differences: PubChem 361336 and
492392, plus both formats of Enamine Z6200214296, Z6200212988 and Z6200174517.
Some alternatives were added and another was lost, so this did not establish
a uniformly better selection or improve reference parity. The trial was
reverted completely. Evidence is retained in `sssr-pruning-trial.patch`,
`algo-rings-sssr-pruned.*` and `sssr-pruning-disagreements.json` under the local
analysis directory. The final runtime files are unchanged from before this
feature review.

Independent evidence: `target/benchmark-parity/audit-sssr-order.py` and
`sssr-order-audit.json`. Temporary diagnostic tracing was removed; no runtime
or runtime-test change remains from this review.

### Benchmark correction

The comparator previously sorted all atoms within a ring, which could equate
different cycles passing through the same vertices. It now canonicalizes each
cyclic path by rotation and reversal, then sorts the ring list. Consecutive
atom pairs, ring multiplicity and source correspondence remain significant.
This strengthens comparison of information already present in the stored
observations; no observation schema, goldens or tolerances changed.

The regression fails before the correction: paths `[0,1,2,3]` and `[0,1,3,2]`
were treated as equal. It now distinguishes them while accepting rotated or
reversed paths and reordered ring lists, and rejecting an omitted ring.

### Full rerun

`cargo run --release --locked --offline --target-dir target/mmcif-validation -p kekule-bench -- --feature algo.rings.sssr --dataset all --output target/benchmark-parity/algo-rings-sssr-after.json`
completed and correctly exited 1 for the retained disagreement and unsupported
CXSMILES inputs.

| Dataset | Applicable | Agree before → after | Disagree | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 200,000 | 199,999 → 199,999 | 1 | 0 |
| Enamine | 100,480 | 97,599 → 97,599 | 0 | 2,881 |
| PL-Rex | 328 | 328 → 328 | 0 | 0 |
| Smoke | 25 | 25 → 25 | 0 | 0 |
| Total | 300,833 | 297,951 → 297,951 | 1 | 2,881 |

Another 1,001 selected cases are not applicable. The paired audit verifies all
301,834 rows retain identical raw actual/reference observations and statuses.
Only the original disagreement's rendered details change, now retaining cyclic
paths. All 2,881 errors retain the explicit CXSMILES rejection diagnostic.
Evidence: `verify-sssr-rerun.py`, `sssr-rerun-verification.json` and the before/
after reports under `target/benchmark-parity/`. The final full run after restoring
the runtime trial is `algo-rings-sssr-final.*`; its case archive is byte-identical
to the verified after-run (`sssr-restoration-verification.json`). The table above
therefore also describes the final restored state.

### Validation

Passed: workspace formatting, all-target/all-feature check, clippy with warnings
denied, all-feature tests including doctests, Rustdoc with warnings denied,
benchmark package file-list check, 25 pinned RDKit tests and eight shared
reference-runner tests, Rust 1.89 workspace check, 21 Python dashboard tests
and the Node dashboard regression. After restoring the trial, all 96 focused
ring tests and the comparator regression passed again; the three trial runtime
files have no remaining diff. `git diff --check` passes. Rust commands used `--locked --offline`
where applicable and the clean `target/mmcif-validation` build directory
(the MSRV check used the existing target directory). Logs are
`target/benchmark-parity/sssr-*.log`.

Runtime package builds and the potentials optional-feature matrix were not
repeated: the final change is confined to the benchmark comparator, its
regression and documentation. Linux fuzz checks were not run on this Windows
host. Biopython/DSSP and unrelated feature benchmarks were not rerun because
they do not observe selected molecular rings. The unpublished benchmark uses
its package file-list check.

Next feature: `algo.valence.rdkit-like`.

## 13. RDKit-like valence

The fresh baseline has 316 disagreements, all in PubChem: 158 source IDs in
both SDF and SMILES. Every changed atom differs only in formal charge and
explicit valence. Explicit and implicit hydrogen counts already agree.
The representative input `OCl(=O)=O` compares Kekule's published
charge-separated chlorate representation against RDKit's unnormalized source.

### Reference preparation correction

Kekule's molecule publication already implements the general oxohalogen rule
in the pinned RDKit
[halogenCleanup implementation](https://raw.githubusercontent.com/rdkit/rdkit/Release_2026_03_3/Code/GraphMol/MolOps.cpp).
For a neutral Cl, Br or I with oxygen-only neighbors and explicit valence
3, 5 or 7, double bonds to oxygen become single bonds with corresponding
charge separation. Ester oxygen is allowed; a carbon neighbor prevents this
conversion. This is represented-chemistry normalization, not a defect in the
subsequent implicit-hydrogen calculation.

The reference adapter previously called only non-strict `UpdatePropertyCache`
on unsanitized input. It now calls RDKit `Cleanup` on its private copy first.
The reference operation is independently implemented by RDKit and includes
its nitrogen and phosphorus cleanup too; it is not a hand-coded expected-value
rewrite restricted to observed failing inputs. Full sanitization is deliberately
absent because it would additionally infer radicals and perceive aromaticity.
All atom fields, source correspondence and exact comparison remain intact.
No runtime code or observation schema changed in this feature review.

Three regressions cover all three halogens at the three relevant valences,
ester oxygen, non-oxygen neighbors, precharged halogens, nitrogen/phosphorus
cleanup, source-copy preservation and exclusion of full sanitization. They
produce 12 failing subtests before the adapter correction and pass afterward.

### Full reference audit

All five dataset/feature pairs were regenerated independently using RDKit
2026.03.3. The paired audit checked all 301,834 rows: exactly the original
316 PubChem discrepancies change, with only formal charge and explicit valence
affected. Every changed record preserves total charge and both hydrogen counts;
all source IDs, paths, record indices, input digests, titles, atom correspondence
and reference versions remain identical. The other 301,518 rows are unchanged,
including all missing-format cases. No reference error was introduced.

The full operation also exercises RDKit's nitrogen/phosphorus cleanup without
restricting it to halogens; no additional change was observed in this corpus.
This does not establish normalization parity for every possible input outside
the supplied corpus.

After this audit, the five generated archive/manifest pairs were promoted and
their copied file hashes verified. The prior pairs remain in
`target/benchmark-parity/valence-original-goldens/`. Evidence:
`valence-generation.*`, `verify-valence-goldens.py` and
`valence-golden-verification.json`. No expected value was derived from Kekule.

The full baseline selected 301,834 rows: 297,636 agree, 316 disagree, 2,881
error and 1,001 are not applicable. Every error is the explicit unsupported
CXSMILES-extension diagnostic. Missing-format cases remain separate from
these implementation errors. Evidence: `algo-valence-before.*`,
`valence-baseline-summary.json` and `valence-original-disagreements.json`
under `target/benchmark-parity/`.

### Validation

Passed: workspace formatting, all-target/all-feature check, clippy with warnings
denied, all-feature tests including doctests, Rustdoc with warnings denied,
benchmark package file-list check, 28 pinned RDKit tests and eight shared
reference-runner tests. Rust commands used `--locked --offline` where applicable
and `--target-dir target/mmcif-validation`. Logs are
`target/benchmark-parity/valence-*.log`.

Runtime package builds, the Rust 1.89 check and the potentials optional-feature
matrix were not repeated: this feature changes only the benchmark's Python
reference adapter, its tests, reference data and documentation. Linux fuzz
checks were not run on this Windows host. Biopython/DSSP, dashboard tests and
unrelated feature benchmarks were not rerun because their code and observations
are unaffected. The unpublished benchmark uses its package file-list check.

### Full native rerun

`target/mmcif-validation/release/kekule-bench.exe --feature algo.valence.rdkit-like --dataset all --output target/benchmark-parity/algo-valence-after.json`
completed against the promoted references. It correctly exits 1 because the
unsupported CXSMILES cases remain errors.

| Dataset | Applicable | Agree before → after | Disagree before → after | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 200,000 | 199,684 → 200,000 | 316 → 0 | 0 |
| Enamine | 100,480 | 97,599 → 97,599 | 0 → 0 | 2,881 |
| PL-Rex | 328 | 328 → 328 | 0 → 0 | 0 |
| Smoke | 25 | 25 → 25 | 0 → 0 | 0 |
| Total | 300,833 | 297,636 → 297,952 | 316 → 0 | 2,881 |

Another 1,001 selected cases are not applicable. There is no remaining valence
disagreement in the supported supplied corpus. CXSMILES input support remains
deferred as documented in the earlier feature reviews.

The final paired audit checks all 301,834 before/after rows. Native observations
are identical throughout. Exactly 316 statuses move from disagreement to
agreement against the independently corrected reference; all other complete
rows, including error diagnostics and missing-format results, remain identical.
Evidence: `verify-valence-rerun.py`, `valence-rerun-verification.json` and the
before/after reports under `target/benchmark-parity/`. `git diff --check` passes.

Next feature: `algo.aromaticity.rdkit-like`.

## 14. RDKit-like aromaticity — 2026-09-18

The complete feature rerun reproduces the original result exactly: one
disagreement, in PubChem 181201's supplied SMILES
`C1=CC(=CC=[C]1)N`. All six ring atom flags and all six ring bond flags differ.
The supplied SDF for the same ID agrees. This is an upstream interpretation
difference, not evidence of a different aromaticity assignment on the same
represented radical state.

### Deferred radical-inference gap

RDKit's full sanitization infers one unpaired electron at SMILES atom 5 before
aromaticity assignment. Kekule's source interpretation preserves the bracket
hydrogen declaration without inferring a radical; its default perception never
rewrites represented chemistry. The independent RDKit 2026.03.3 experiment
turns off only `SANITIZE_FINDRADICALS`: every resulting atom and bond flag then
matches the native observation exactly. With full sanitization, all flags
reproduce the stored reference. The supplied SDF contains an explicit doublet;
both engines identify its six aromatic ring atoms and six aromatic ring bonds.

This follows the documented RDKit
[sanitization order and neutral-carbon radical eligibility](https://www.rdkit.org/docs/RDKit_Book.html).
The existing runtime regression
`neutral_carbon_radical_can_complete_an_aromatic_sextet` also exercises this
structure with an explicitly installed doublet and verifies its one-electron
donation and aromatic sextet. The separate
`bracket_atoms_do_not_infer_radicals_from_a_valence_model` regression asserts
the current interpretation policy across carbon, nitrogen, oxygen and aromatic
bracket forms.

The missing SMILES radical inference is a real reference-parity gap; the review
does not claim that the native interpretation is scientifically superior.
Resolving it requires an explicit general policy for deriving represented
radicals from valence, including its effects on parsing, perception and writers.
That changes an established architectural/API contract beyond a local
aromaticity correction. No case-specific electron-count adjustment was added,
and radical inference was not disabled in the benchmark reference.
Evidence: `target/benchmark-parity/audit-aromaticity.py` and
`aromaticity-source-audit.json`.

### Benchmark audit and errors

The adapters retain every atom flag and every bond flag, keyed by source atom
correspondence. The comparator checks those booleans exactly and does not infer
aromatic bonds from aromatic endpoints. The guide now states explicitly that
this feature includes each engine's default preparation and its upstream
interpretation differences; it is not an isolated comparison on equalized
radical states. No executable code, golden, schema or tolerance changed.

Eighteen PubChem errors are nine IDs in both supplied formats, rejected during
valence preparation by both engines. Independent RDKit reevaluation of the nine
SMILES reproduces `AtomValenceException` for each: 24594, 61654, 77880, 118990,
139908, 139910, 139911, 141144 and 167661. They involve halogen/interhalogen or
metal-bound halide valences outside the strict model. This is not a claim that
these compounds cannot exist chemically; both configured models reject the
supplied representations. The errors stay visible and are not counted as
agreement. The remaining 2,881 errors are the established CXSMILES rejection.

### Full rerun

`target/mmcif-validation/release/kekule-bench.exe --feature algo.aromaticity.rdkit-like --dataset all --output target/benchmark-parity/algo-aromaticity-before.json`
completed and correctly exited 1 for the retained disagreement and errors.
No second full run was necessary because this review changes documentation
only. The paired audit compares every raw result to the original run, including
actual/reference observations, source correspondence, diagnostics and statuses:
all 301,834 rows are identical.

| Dataset | Applicable | Agree | Disagree | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 200,000 | 199,981 | 1 | 18 |
| Enamine | 100,480 | 97,599 | 0 | 2,881 |
| PL-Rex | 328 | 328 | 0 | 0 |
| Smoke | 25 | 25 | 0 | 0 |
| Total | 300,833 | 297,933 | 1 | 2,899 |

Another 1,001 selected cases are not applicable. The full evidence is
`verify-aromaticity-rerun.py`, `aromaticity-rerun-verification.json` and
`algo-aromaticity-before.*` under `target/benchmark-parity/`.

### Validation

Passed: workspace formatting, all-target/all-feature check, clippy with warnings
denied, all-feature tests including doctests and the two named radical-policy
regressions, Rustdoc with warnings denied, benchmark package file-list check,
28 pinned RDKit tests and eight shared reference-runner tests. `git diff --check`
passes. Rust commands used `--locked --offline` where applicable and
`--target-dir target/mmcif-validation`; logs are `aromaticity-*.log` in the
local analysis directory.

Runtime package builds, the Rust 1.89 check, the potentials optional-feature
matrix and dashboard/Biopython/DSSP tests were not repeated: this feature review
changes documentation only. Linux fuzz checks were not run on this Windows
host. No golden regeneration was needed because reference preparation and
observations were retained. The unpublished benchmark uses its package
file-list check. No new regression test was needed because there is no defect
fix or behavior/API change; the relevant existing regressions passed again.

Next feature: `algo.canonical-ranking`.

## 15. Canonical ranking — 2026-09-19

The full baseline has 83 disagreements: 81 PubChem SMILES inputs and both
formats of Enamine Z4453364263. Every difference concerns an atom-equivalence
partition. No title, component correspondence or other observation field differs.

### Chemical hydrogen-equivalence correction

The ranking invariant included both the total hydrogen count and its storage
policy: the declared count and whether implicit hydrogens were allowed. This
split otherwise equivalent atoms, for example the terminal carbons of
`[CH3]CC`, and propagated the distinction throughout symmetric structures.
It also distinguished oxygens in normalized chlorate/perchlorate ions based on
whether their original spelling was bracketed.

Ranking now uses the current total hydrogen count, independent of fixed versus
inferred storage. The two redundant declaration-policy fields were removed
from the initial ranking signature; graph refinement and the ring-topology
algorithm are unchanged. Different totals, isotopes, charges, atom maps and
connectivity still contribute. The input graph, hydrogen declarations and
installed perception are never rewritten. This is a bounded correction to
derived chemical equivalence, unlike changing radical interpretation in feature
14. The public API documentation now states that installed implicit counts
contribute to ranking and that only declared counts contribute before valence
perception.

The old declaration-sensitive ranking regression was updated to the corrected
contract while retaining its hydrogen-count assertions and adding checks that
the stored declarations and molecule remain unchanged. A second regression
checks several hydrogen spellings of propane and rejects equivalence when a
terminal hydrogen count actually differs. Both fail before the correction and
pass afterward. Existing isotope/map, atom-order and aromatic-bond ranking
regressions remain intact.

The canonical writer already projects hydrogen representation on a private copy
and includes the serialized atom label in its complete labeling certificate.
It does not need the removed distinction in general chemical ranking to
preserve represented output. Its projection comment was clarified; no writer
algorithm changed.

### Independent audit and benchmark contract

The reference invokes RDKit's
[canonical ranking API](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.CanonicalRankAtoms)
with tie breaking and chirality disabled, and isotope/atom-map distinctions
enabled. The benchmark compares partitions of source atom indices, so arbitrary
numeric rank labels do not affect agreement. Every class and member remains
asserted; no reference, golden, schema, tolerance or comparator changed.

The independent source audit recomputes the original RDKit 2026.03.3 partitions
for all 83 discrepant inputs. An analysis-only refinement of those partitions
by source hydrogen declarations and ordinary neighbor propagation reproduces
the previous native partition exactly in 82 cases. That diagnostic does not
rewrite a molecule or supply expected benchmark values. It isolates the two
extra native invariants as the cause. The remaining case is PubChem 181201's
SMILES radical-inference gap already established in feature 14. Source and
diagnostic evidence are `audit-ranking-hydrogens.py` and
`ranking-hydrogen-audit.json` under `target/benchmark-parity/`.

### Full rerun

`cargo run --release --locked --offline --target-dir target/mmcif-validation -p kekule-bench -- --feature algo.canonical-ranking --dataset all --output target/benchmark-parity/algo-ranking-after.json`
completed against the unchanged stored references. It correctly exits 1 for
the retained radical-inference disagreement and input/preparation errors.

| Dataset | Applicable | Agree before → after | Disagree before → after | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 200,000 | 199,901 → 199,981 | 81 → 1 | 18 |
| Enamine | 100,480 | 97,597 → 97,599 | 2 → 0 | 2,881 |
| PL-Rex | 328 | 328 → 328 | 0 → 0 | 0 |
| Smoke | 25 | 25 → 25 | 0 → 0 | 0 |
| Total | 300,833 | 297,851 → 297,933 | 83 → 1 | 2,899 |

Another 1,001 selected cases are not applicable. The 18 shared valence failures
and 2,881 unsupported CXSMILES errors retain their causes described in feature
14. No new disagreement or error appears. The only deferred ranking case is
the upstream radical-inference policy; no unresolved hydrogen-equivalence
difference remains in the supplied corpus.

The final paired audit checks all 301,834 before/after rows. Exactly the 82
independently explained cases change from disagreement to agreement, with only
their native class partitions changing. All reference observations and source
identities remain identical. Every other complete row, including the remaining
disagreement, all errors and all missing-format results, is unchanged. Evidence:
`verify-ranking-rerun.py`, `ranking-rerun-verification.json` and the before/after
reports in the local analysis directory.

### Validation

Passed: `cargo fmt --all -- --check`; workspace all-target/all-feature check;
clippy with warnings denied; workspace all-feature tests including doctests;
Rustdoc with warnings denied; Rust 1.89 workspace all-target/all-feature check;
potentials tests and docs with default features disabled; full `kekule` package
build and verification; companion and benchmark package file-list checks; and
packaged-license hash comparisons. The six focused ranking tests pass, including
the two regressions that failed before the fix. Also passed: 28 pinned RDKit
tests, eight shared reference-runner tests, 21 Python dashboard tests and the
Node dashboard regression. `git diff --check` passes.

Canonical and isomeric SMILES writer smoke benchmarks each have eight applicable
cases, all agreeing with zero errors. Canonical writer verification includes
its fixed-point and atom/bond renumbering checks. Existing full workspace
serialization tests also pass. Their complete external writer corpora were not
repeated; the complete external corpus run in this turn is canonical ranking.

Rust commands used `--locked --offline` where applicable, with
`--target-dir target/mmcif-validation` except the separate MSRV check; package
commands used `--allow-dirty`. Logs and writer reports are `ranking-*.log` and
`ranking-*-writer-smoke.*` under `target/benchmark-parity/`. Full companion
package builds were not run because the CI contract uses file-list checks
until the foundational crate is published. Linux fuzz build/smoke checks were
not run on this Windows host. Biopython/DSSP and unrelated feature benchmarks
were not repeated. No golden regeneration was needed or performed.

Next feature: `algo.substructure.vf2`.

## 16. VF2 substructure matching — 2026-09-19

The full baseline has 20,512 disagreements: 7,743 PubChem, 12,757 Enamine,
10 PL-Rex and two smoke cases. The query-by-query audit isolates incorrect
carbon–oxygen single-bond matches in aromatic heterocycles and four incorrect
benzene-pattern matches in aryne inputs. The remaining source is PubChem 181201's
known SMILES radical-inference gap. Counts by query can exceed case counts
because a source can contain several components.

### SMARTS bond semantics

The query IR correctly defines `BondPredicate::Order` as a localized represented
order and `BondPredicate::Aromatic` as perceived membership. The SMARTS parser
previously used either primitive alone for SMARTS bond types. That let `-` and
`=` match localized bonds inside aromatic rings, and let `:` and omitted bonds
match a triple bond carrying an aromaticity flag.

The parser now expresses single/double types as the requested localized order
AND nonaromatic membership. Aromatic type combines aromatic membership with
localized single-or-double order. Triple and quadruple types retain their
represented order. Omitted bonds combine the corrected single and aromatic
types. This follows the independent RDKit behavior and the documented
[SMARTS bond types](https://www.rdkit.org/docs/RDKit_Book.html#smarts-support-and-extensions).
It applies uniformly across atoms, elements, rings and source formats.
The VF2 search, atom predicates, concrete molecule representation and low-level
programmatic predicates are unchanged.

Matching explicit SMARTS single/double bonds now requires installed aromaticity,
as other aromaticity-dependent predicates already do. It does not perceive a
target implicitly. The public parser documentation and query example were
updated accordingly. Direct programmatic order queries still work on an
unperceived graph and still inspect localized orders.

The new ring-closure regression also exposed a syntax guard that rejected any
explicit bond before a ring label, despite the subsequent code supporting it.
The guard now accepts the already-validated bond token. Tests cover bond
declarations at either end of a ring closure and retain rejection of conflicting
declarations. No SMARTS grammar workaround or molecule-specific rule was added.

### Regressions and independent checks

The two new regressions fail before the bond correction. They cover furan,
biphenyl, an aryne, a nonaromatic alkene ring and an acyclic alkyne across
single, double, triple, aromatic, omitted and any-bond queries, including both
mapping orientations. The counts were independently checked with pinned RDKit
2026.03.3. They also verify that an aromatic-flagged triple retains its triple
type, that the single link between aromatic biphenyl rings still matches, and
that programmatic order/membership predicates retain their separate meanings.
All 16 query tests pass after the correction, including parser mutation,
precedence, hydrogen, uniqueness, non-induced matching and hard-limit tests.

### Benchmark audit

The reference uses `uniquify=False, maxMatches=0`; native matching uses
`uniquify: false, max_matches: usize::MAX`. The existing reference regression
checks more than 1,000 matches and both orientations of a carbon–carbon edge.
Both adapters preserve query-atom order within a mapping and sort only the
complete mapping list. All 18 queries and every source component remain
asserted. Search/candidate exhaustion is reported as an error, rather than a
partial successful result. No match was removed from a reference or normalized
away during comparison, and no reference/golden/schema/tolerance changed.
The guide now makes these conventions explicit.

Source evidence: `vf2-original-disagreements.json`, `vf2-original-summary.json`
and `vf2-original-scan.log` under `target/benchmark-parity/`. The independent
bond-type probe is preserved in `probe-vf2-bond-types.py` and
`vf2-bond-type-probe.json` in the same directory.

### Full rerun

`cargo run --release --locked --offline --target-dir target/mmcif-validation -p kekule-bench -- --feature algo.substructure.vf2 --dataset all --output target/benchmark-parity/algo-vf2-after.json`
completed against the unchanged stored references. It correctly exits 1 for
the remaining radical-inference disagreement and preparation/input errors.

| Dataset | Applicable | Agree before → after | Disagree before → after | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 200,000 | 192,239 → 199,981 | 7,743 → 1 | 18 |
| Enamine | 100,480 | 84,842 → 97,599 | 12,757 → 0 | 2,881 |
| PL-Rex | 328 | 318 → 328 | 10 → 0 | 0 |
| Smoke | 25 | 23 → 25 | 2 → 0 | 0 |
| Total | 300,833 | 277,422 → 297,933 | 20,512 → 1 | 2,899 |

Another 1,001 selected cases are not applicable. The sole remaining disagreement
is the SMILES radical-inference gap of PubChem 181201 from feature 14; atom
aromaticity differences propagate into several queries. The same 18 shared
strict-valence failures and 2,881 CXSMILES rejections remain errors. No new
disagreement or error appears.

The final paired audit checks all 301,834 raw before/after rows. Exactly 20,511
cases move from disagreement to agreement. Their only changed observations are
the removal of 36,401 false carbon–oxygen single-bond mappings and 48 false
benzene mappings across four aryne cases. All reference observations, other
queries, source identities and component metadata are unchanged. Every other
complete row, including the remaining disagreement and every error or missing
format, is identical. Evidence: `verify-vf2-rerun.py`,
`vf2-rerun-verification.json` and `algo-vf2-before.*` / `algo-vf2-after.*` under
`target/benchmark-parity/`.

### Validation

Passed: workspace formatting, all-target/all-feature check, clippy with warnings
denied, all-feature tests including doctests, Rustdoc with warnings denied,
Rust 1.89 workspace check, potentials tests/docs without default features,
full `kekule` package build and verification, companion/benchmark package
file-list checks and packaged-license hash comparisons. All 16 focused query
tests pass. Also passed: 28 pinned RDKit tests, eight shared reference-runner
tests, 21 Python dashboard tests and the Node dashboard regression. The final
parser-documentation clarification passed another warnings-denied Rustdoc
build. `git diff --check` passes.

Rust commands used `--locked --offline` where applicable and
`--target-dir target/mmcif-validation` except the separate MSRV check; package
commands used `--allow-dirty`. Logs are `vf2-*.log` under the local analysis
directory. Full companion package builds were not run because CI uses file-list
checks until the foundational crate is published. Linux fuzz build/smoke checks
were not run on this Windows host. Biopython/DSSP and unrelated feature
benchmarks were not repeated. No golden regeneration was needed or performed.

Next feature: `query.smarts`.

## 17. SMARTS parsing

The explicit ring-bond syntax correction from feature 16 resolves 37,030 cases
in this feature: three PubChem records and 37,027 Enamine records. The parser
already supported pending bond expressions on ring closures, but its entry
guard rejected a bond token before a ring label. The general guard correction
accepts these valid queries, with no source-specific handling. Existing new
regressions exercise explicit bonds at either end and conflicting declarations.
No additional runtime or reference change was needed in this turn.

### Full corpus and transition audit

`target/mmcif-validation/release/kekule-bench.exe --feature query.smarts --dataset all --output target/benchmark-parity/query-smarts-before.json`
completed against the unchanged stored references. The `before` filename denotes
the start of this feature's review, after the feature-16 correction. It correctly
exits 1 for unsupported queries.

| Dataset | Applicable | Original agrees → fresh agrees | Original errors → fresh errors | Disagreements |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 100,000 | 87,461 → 87,464 | 12,539 → 12,536 | 0 |
| Enamine | 50,240 | 4,847 → 41,874 | 45,393 → 8,366 | 0 |
| Smoke | 8 | 6 → 6 | 2 → 2 | 0 |
| Total | 150,248 | 92,314 → 129,344 | 57,934 → 20,904 | 0 |

PL-Rex and PDB have no applicable SMARTS text files. Together with 12 smoke
records, 1,176 selected records are not applicable.

The paired audit compares all 151,424 original/fresh rows, including source
identity, reference observations, failures and missing formats. Exactly 37,030
rows change from the ring-label syntax error to exact agreement. Another 3,001
previous ring-label errors now reach a later unsupported-stereochemistry error.
Every other complete row is unchanged. All expected observations and source
identities are unchanged, including in changed rows. There are no regressions.
Evidence: `verify-smarts-rerun.py`, `smarts-rerun-verification.json`,
`smarts-original-summary.json` and `query-smarts-before.*` under
`target/benchmark-parity/`.

### Remaining scope and benchmark interpretation

Every remaining error explicitly rejects unsupported query stereochemistry:
20,826 have stereochemical atom syntax and 78 have directional/stereochemical
bond syntax. These counts classify the first reported unsupported construct;
a query can contain both. They are feature limitations, not reference failures.
The query representation currently contains local atom/bond predicates and
connectivity, with no stereo constraints or mapping-relative carrier ordering.
Correct support would extend the query representation and matching semantics,
including checking parity under query-to-target permutations. Simply accepting
or stripping the syntax would lose the requested chemical constraint. That
larger implementation is deferred under this review's scope. The existing
unsupported-semantics regression asserts structured errors for both categories.

Both adapters interpret the first whitespace-delimited token of each text
record as SMARTS and the rest as its title. The current external datasets supply
SMILES strings, exercising the overlapping grammar. The result compares parse
acceptance and atom/bond counts; agreement does not establish equivalent query
predicates, molecular interpretation, or handling of CXSMILES extensions. The
18 independent behavioral substructure queries from feature 16 and focused
runtime tests provide complementary coverage, not a complete SMARTS conformance
suite. The guide now states these boundaries explicitly. Unsupported queries
remain errors in the full denominator, and all source records remain selected.

A fresh independent RDKit 2026.03.3 audit reparses all 150,248 stored query
strings and verifies successful reference parsing and both graph counts. All
reference observations agree with that audit; no golden regeneration or
comparison weakening is needed. Evidence: `audit-smarts-reference.py`,
`smarts-reference-audit.json` and its log in the local analysis directory.

### Validation

Passed: workspace formatting, all-target/all-feature check, clippy with warnings
denied, all-feature tests including query regressions and doctests, Rustdoc with
warnings denied, benchmark package file-list check, 28 pinned RDKit reference
tests and eight shared runner tests. Rust commands used `--locked --offline`
where applicable and `--target-dir target/mmcif-validation`; packaging used
`--allow-dirty`. Logs are `smarts-*.log` under `target/benchmark-parity/`.

No new runtime, reference or comparator code changed in this turn. The MSRV,
optional-feature tests/docs, full runtime package verification, companion
package lists, license hashes and dashboard tests passed in feature 16 and were
not repeated for this guide/review-only change. Full companion package builds
remain deferred by CI until the foundational crate is published. Linux fuzz
checks were not run on this Windows host. Other feature benchmarks and unrelated
Biopython/DSSP tests were not repeated.

Next feature: `chem.perception.default`.

## 18. Default chemical perception

The complete rerun finds no new isolated default-perception defect. It verifies
the earlier aromatic triple-bond observation fix in four cases: PubChem 141133
and 141134 in both SMILES and SDF. The remaining differences are hydrogen
representation/inference policy, radical inference and source stereo or stereo
cleanup. No new runtime, reference, comparator or golden change was made here.

### Full corpus verification

`target/mmcif-validation/release/kekule-bench.exe --feature chem.perception.default --dataset all --output target/benchmark-parity/chem-perception-before.json`
completed against unchanged references. The `before` filename denotes the start
of this feature review, after earlier fixes. Exit 1 correctly reports remaining
disagreements and errors.

| Dataset | Applicable | Original agrees → fresh agrees | Fresh disagreements | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 200,000 | 186,506 → 186,510 | 13,472 | 18 |
| Enamine | 100,480 | 69,699 → 69,699 | 27,900 | 2,881 |
| PL-Rex | 328 | 290 → 290 | 38 | 0 |
| Smoke | 25 | 20 → 20 | 5 | 0 |
| Total | 300,833 | 256,515 → 256,519 | 41,415 | 2,899 |

Another 1,001 selected cases are not applicable. The complete before/after audit
checks all 301,834 rows. Exactly four change from disagreement to exact
agreement; each only corrects one aromatic-flagged triple bond's reported type
from `AROMATIC` to `TRIPLE`. Every other observation and reference value in those
rows is unchanged. All other 301,830 rows are byte-for-byte identical to the
original run. There are no regressions. Evidence: `verify-perception-rerun.py`,
`verify-perception-changes.py`, `perception-rerun-verification.json`,
`perception-final-verification.json` and `chem-perception-before.*` under
`target/benchmark-parity/`.

### Remaining disagreements

Classification uses the comparator's complete differences after documented bond
ordering normalization, rather than mistaking raw bond-list order for changed
chemistry. Every remaining difference belongs to the following categories;
counts overlap where a case has multiple causes.

| Category | Cases | Finding |
| --- | ---: | --- |
| Explicit/implicit hydrogen storage | 27,593 | All 30,313 affected atoms have equal total H; explicit-valence differences exactly track declared-H differences |
| Future implicit-H inference policy | 8,245 | 14,148 atoms differ in `no_implicit_hydrogens`, principally omitted wedge-H interpretation |
| Represented radicals | 938 | 933 SMILES and five SDF observations retain the previously identified radical-inference gap |
| Source stereo / cleanup | 6,113 | 6,036 PubChem, 39 Enamine and 38 PL-Rex cases need stereo interpretation/cleanup review |
| Aromaticity downstream of radicals | 1 | PubChem 181201 SMILES, also counted in the radical category |

The disjoint classification has 34,364 hydrogen-only cases, 6,100 stereo-only
cases, 916 radical-only cases, 21 hydrogen/radical cases, 13 hydrogen/stereo
cases and the one radical/aromaticity case. These sum to all 41,415 remaining
disagreements. The hydrogen audit checks each changed graph atom's explicit-plus-
implicit H difference is zero and its represented-valence difference equals its
declared-H difference. No disagreement in element, charge, isotope, atom map,
atom count, bond count or connectivity remains unexplained by this audit.

Hydrogen storage and inference-policy distinctions are not evidence of missing
hydrogen atoms. Changing them during default perception would rewrite represented
chemistry, contrary to the architecture. A coordinated interpretation/policy
change belongs outside this bounded perception review. The earlier radical
investigation also found parser/writer contract dependencies; default perception
must not silently install represented radical state to force equality.

Stereo differences remain fully asserted and will be revisited under the stereo
features. Earlier source investigations include non-stereogenic assertions that
RDKit cleans up, ambiguous/crowded wedge interpretation, missing double-bond
elements and three mistaken axis assignments. Default perception does not
materialize stereo or rewrite source assertions, so these are not fixes to its
derived-state pipeline. No blanket stereo filtering or parity flip was applied.

The 18 PubChem errors are the same nine source IDs in both formats whose strict
valence preparation fails in both engines. The 2,881 Enamine errors reject
unsupported CXSMILES extensions. Matching failures remain errors, not agreements.
The earlier independent source audits remain applicable because these rows are
byte-identical.

Evidence: `scan-perception.py`, `perception-original-summary.json`,
`perception-original-disagreements.json`, `classify-perception.py` and
`perception-classification.json` in the local analysis directory. The initial
raw-order diagnostic was superseded by `perception-normalized-scan.log` and the
comparator-based classification; raw bond-list order is not a chemical defect.

### Benchmark audit and scope

The native adapter calls the public `Molecule::perceive` on each connected
component, then reads the full indexed graph, basic atom properties and valence.
The reference separately parses the source, preserves explicit H vertices,
splits all components, removes conformers, and runs full RDKit sanitization and
stereo cleanup on a copy. Its broader preparation is intentional for the
end-to-end prepared-state comparison, and its consequences remain visible.
RDKit documents radical assignment, hydrogen adjustment and chirality cleanup
as part of [molecular sanitization](https://www.rdkit.org/docs/RDKit_Book.html#molecular-sanitization).
The guide now distinguishes this comparison from perception on identical
represented graphs. Coordinates, SDF metadata and selected ring lists are not
claimed as observations of this feature; their separate features cover them.
All graph and valence fields, component identities and source membership remain
asserted. No reference values, fields or tolerances were changed.

Fresh pinned RDKit 2026.03.3 source probes reproduce the reported reference
differences for one external example from every category: PubChem 32, 117,
2831 and 181201, plus Enamine Z9077138344. They also verify that preparation
leaves each source fragment unchanged. These are representative source probes,
not a regeneration of all reference records. Evidence:
`audit-perception-sources.py`, `perception-source-audit.json` and its log.

### Validation

Passed: workspace formatting, all-target/all-feature check, clippy with warnings
denied, all-feature tests including doctests, Rustdoc with warnings denied,
benchmark package file-list check, 28 pinned RDKit reference tests and eight
shared runner tests. Existing regressions cover separation of interpretation
from perception, represented stereo preservation, absence of coordinate-only
stereo assignment, idempotence and transactional failure. Rust commands used
`--locked --offline` where applicable and `--target-dir target/mmcif-validation`;
packaging used `--allow-dirty`. Logs are `perception-*.log` under the analysis
directory. `git diff --check` passes.

No runtime/reference/comparator code changed this turn. MSRV, optional-feature
tests/docs, full runtime package verification, companion package lists, license
hashes and dashboard tests passed in feature 16 and were not repeated for this
guide/review-only change. Full companion package builds remain deferred by CI
until the foundational crate is published. Linux fuzz checks were not run on
this Windows host. Other feature benchmarks and unrelated Biopython/DSSP tests
were not repeated; no golden regeneration was necessary.

Next feature: `chem.hydrogen-transforms`.

## 19. Hydrogen transformations

### Generated hydrogen declarations

New graph hydrogens were explicitly assigned `HydrogenDeclaration::Fixed(0)`,
although the ordinary neutral-hydrogen constructor defaults to inferred
hydrogens. RDKit `AddHs` uses the ordinary default for these atoms. This created
a widespread mismatch in `added_graph.atoms[].no_implicit_hydrogens` despite
equal current hydrogen counts and connectivity.

The addition loop now inserts `Atom::new(hydrogen)` directly. Its single bond
already satisfies hydrogen's valence, so subsequent perception infers zero
additional H. This removes the special declaration override for every generated
hydrogen, irrespective of parent element, charge, aromaticity or source format.
It preserves every original atom's inference policy and consumes only the
hydrogen counts being materialized. The removal algorithm is unchanged.

The new regression fails on the old code and passes after the correction.
It covers methane, ammonia, water, ammonium, aromatic nitrogen and bracket-fixed
carbon in both full and explicit-only addition modes. It asserts normal
generated-atom declarations, single attachment, preservation of parent policies,
zero inferred H on generated atoms, and no graph growth on repeated addition.
Independent RDKit 2026.03.3 probes confirm its generated-atom default and
idempotence across those inputs. These are focused unit examples, not additions
to the externally supplied benchmark corpus.

### Removal policies and reference limitations

The original corpus has 781 cases, spanning 785 component observations, where
RDKit retains 897 more explicit hydrogen vertices after removal. The complete
atom arrays in every such comparison have identical heavy-atom sequences
(element, isotope, charge, radical and map) and equal total hydrogens, including
both encoded and graph H. This audit does not claim full stereo equivalence in
every source. Source-stereo differences remain independently visible.

RDKit's default policy preserves certain hydrogens defining double-bond stereo,
as described by its [hydrogen-removal API](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RemoveHs).
Kekule can retain their role with typed implicit-hydrogen carriers. The external
PubChem 586 probe confirms that enabling RDKit's `removeDefiningBondStereo`
option produces the same collapsed atom count as Kekule; the benchmark keeps
RDKit's unmodified default. Existing regressions cover collapse/re-addition of
double-bond hydrogen carriers and preservation of tetrahedral carriers and groups.

Two additional SDF cases expose an actual reference loss of hydrogen count:
PubChem 88146, an iridium complex, changes from 31 total H to 30 under RDKit's
default removal; PubChem 136981, a molybdenum compound, changes from two to zero.
The native collapsed observations retain those hydrogens in the encoded count.
Independent pinned RDKit probes reproduce the losses from the original source
records. Kekule's count-preserving behavior is preferable here and is retained.
A regression now protects collapse and re-addition for parents whose valence
model does not infer replacement H. No reference hydrogen count was patched,
and these disagreements remain reported rather than being hidden.

Evidence under `target/benchmark-parity/`: `probe-added-hydrogens.py`,
`added-hydrogens-reference-probe.json`, `probe-hydrogen-removal.py`,
`hydrogens-removal-reference-probe.json`, `audit-hydrogen-retention.py` and
`hydrogens-retention-audit.json`. The complete original differences and normalized
field classification are in `hydrogens-original-disagreements.jsonl`,
`hydrogens-original-summary.json` and `hydrogens-classification.json`.

### Benchmark interpretation

Both adapters expand the prepared molecule, retain the complete expanded graph,
count newly added H by original parent, remove all eligible H, and retain the
complete collapsed graph. The label `round_trip` denotes this last state; it
does not promise restoration of the input's explicit-H representation or dense
atom indices. Original H vertices are excluded from the newly-added count.
All component-local atom indices, graph fields, stereo and parent counts remain
asserted. The guide now documents this contract and the reference limitations.
No reference implementation, golden, schema, comparison field or tolerance was
changed.

### Full rerun

Before:
`target/mmcif-validation/release/kekule-bench.exe --feature chem.hydrogen-transforms --dataset all --output target/benchmark-parity/chem-hydrogens-before.json`.
After:
`cargo run --release --locked --offline --target-dir target/mmcif-validation -p kekule-bench -- --feature chem.hydrogen-transforms --dataset all --output target/benchmark-parity/chem-hydrogens-after.json`.
Both completed against unchanged references and correctly exited 1 for remaining
differences and errors.

| Dataset | Applicable | Agree before → after | Disagree before → after | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 200,000 | 93,855 → 192,131 | 106,127 → 7,851 | 18 |
| Enamine | 100,480 | 0 → 89,324 | 97,599 → 8,275 | 2,881 |
| PL-Rex | 328 | 282 → 282 | 46 → 46 | 0 |
| Smoke | 25 | 6 → 21 | 19 → 4 | 0 |
| Total | 300,833 | 94,143 → 281,758 | 203,791 → 16,176 | 2,899 |

Another 1,001 selected cases are not applicable. The new correction resolves
187,615 complete case disagreements. Preparation/input errors are unchanged:
18 shared strict-valence failures and 2,881 unsupported CXSMILES inputs.

The separate original-to-before audit verifies all 301,834 rows. Only the four
known PubChem aryne observations (141133/141134, SMILES/SDF) changed upstream;
each corrects the triple-bond type in both expanded and collapsed graphs.
The two SDF cases become exact agreements, while the SMILES cases initially
retain generated-H declaration differences addressed by this turn's fix.
All other 301,830 rows are byte-identical to the original run. Evidence:
`verify-hydrogens-original-rerun.py`, `verify-hydrogens-upstream-changes.py` and
`hydrogens-upstream-verification.json` in the local analysis directory.

The final paired audit checks all 301,834 before/after rows. It verifies exactly
4,110,948 corrected declaration fields on newly generated hydrogen atoms in
the expanded graphs. One generated hydrogen is also retained in a dihydrogen
component, so its corrected declaration persists in the collapsed observation;
its atom index is unchanged because that component loses no atoms. This is the
same general constructor correction, not a special molecular rule.

Exactly 187,615 cases move to exact agreement, while 9,036 improved cases retain
independent disagreements. The other 105,183 complete rows are byte-identical.
All hydrogen counts, parent assignments, connectivity, stereo, other atom
properties, reference observations and source identities remain unchanged.
No new disagreement or error appears. The report audit also verifies complete
runs, matching source selection and golden identities, and the new executable
fingerprint. Evidence: `verify-hydrogens-rerun.py`,
`hydrogens-rerun-verification.json`, `verify-hydrogens-reports.py` and
`hydrogens-report-verification.json` in the local analysis directory.

The remaining 16,176 disagreements have the following overlapping categories:

| Category | Cases | Remaining scope |
| --- | ---: | --- |
| Original parent inference policy | 8,245 | Source declarations, principally omitted wedge-H interpretation; preserved by the transform |
| Stereo / carrier representation | 6,852 | Source cleanup/orientation differences plus explicit-versus-implicit stereo carriers; retained for stereo review |
| Collapsed encoded-H counts | 3,535 | Storage differences with equal total H, except the two documented reference metal-hydride losses |
| Explicit-H retention | 781 | Equal heavy-atom sequences and total H, with different retained graph H |
| Radical inference | 938 | The previously documented source-preparation gap |
| Aromaticity downstream of radicals | 1 | PubChem 181201, also counted above |

The final comparator differences contain no additional unexplained category.
The original full scan also confirms that added-H counts and parent assignments
already agreed throughout the successful corpus; this correction changes the
generated declarations rather than the number of hydrogens. Larger source-stereo
and declaration-policy changes remain for their owning interpretation features;
no removal-policy override or reference-value adjustment was introduced.

### Validation

Passed: workspace formatting, all-target/all-feature check, clippy with warnings
denied, all-feature tests including doctests, Rustdoc with warnings denied,
Rust 1.89 workspace check, potentials tests/docs without default features,
full `kekule` package build and verification, companion/benchmark package
file-list checks, and packaged-license hash comparisons. All 16 focused
hydrogen tests pass, including both new regressions. Also passed: 28 pinned
RDKit tests, eight shared runner tests, 21 Python dashboard tests and the Node
dashboard regression. Final formatting and clippy checks include both new tests.
`git diff --check` passes.

Rust commands used `--locked --offline` where applicable and
`--target-dir target/mmcif-validation` except the separate MSRV check; package
commands used `--allow-dirty`. Logs are `hydrogens-*.log` under
`target/benchmark-parity/`. Full companion package builds were not run because
CI uses file-list checks until the foundational crate is published. Linux fuzz
build/smoke checks were not run on this Windows host. Biopython/DSSP and unrelated
feature benchmarks were not repeated. No golden regeneration was performed.

Next feature: `descriptor.molecular`.

## 20. Molecular descriptors

### Fix and source verification

The atomic-data generator's default output used the wrong repository ancestor:
`parents[3]` pointed to `C:/Users/chout/repos/crates/...`, outside this checkout.
It now uses `parents[2]`, placing the table in this repository's descriptor module.
Two dependency-free regressions exercise the real entry point with mocked source
acquisition and writes: the default destination must be inside the repository,
and an explicit output must be respected. The default-path regression failed
before the fix; both pass afterward. CI now runs them. No runtime mass formula,
constant, benchmark comparator or golden was changed.

The three published source files were downloaded into the ignored analysis
directory and all existing pinned SHA-256 checks passed. Regenerating into a
separate file produces exactly the checked-in table text, including all 3,557
isotope entries. Offline regeneration from that verified cache also succeeds.
The data sources are CIAAW's 2024 abridged standard weights and natural isotope
compositions, and AME2020 atomic masses. Public descriptor documentation and
the benchmark guide now state the mass and charge conventions explicitly.

### Independent corpus audit

All 338,837 successful component records have identical complete formulas,
including isotope labels and formal charges. Every other non-mass observation
also matches. There are 35,790 charged components and 1,030 components with
explicit isotope labels. All numerical differences are in the two mass fields:
156,199 average-mass fields and 338,833 monoisotopic-mass fields differ in their
raw values. These field counts include differences within the comparator's
existing floating-point allowance; they are not counts of failed cases.

An independent reconstruction from each formula verifies every successful
native and reference mass. The native reconstruction reads the checksum-verified
published sources through the generator's parsers; the reference reconstruction
reads the pinned RDKit 2026.03.3 periodic-table API. This verifies the aggregation
and conventions independently of the molecular traversal in each descriptor.
Changing summation order produces only bounded floating-point roundoff, with
maximum residuals below 1.6e-11 Da. The analysis uses a constituent-count roundoff
bound for this reconstruction only; the benchmark tolerance is unchanged.

Kekule uses CIAAW abridged standard weights for unlabeled atoms and AME2020 masses
for explicitly labeled isotopes. Its monoisotopic mass uses the most abundant
natural isotope, selected from CIAAW compositions, for an unlabeled atom. Both
methods subtract formal charge times the CODATA 2022 electron mass,
0.0005485799090441 Da. RDKit's average mass uses its own weights without an
electron correction; its exact mass uses its own isotope masses and an electron
correction of approximately 0.00054857991 Da, independently probed with a proton.
For example, native protium is 1.007825031898 Da and the pinned RDKit value is
1.007825032 Da. These differences explain the entire successful corpus.

The published data and consistent ion-mass correction are scientifically
defensible and are retained. Abridged standard weights are conventional point
estimates with uncertainties, so this does not claim that every native average
weight is more accurate for an arbitrary real sample. Numerical compatibility
with RDKit would require choosing its constants and average-mass convention;
there is no unexplained descriptor arithmetic defect to fix. No compatibility
mode or corpus-specific correction was added.

The primary sources are [CIAAW abridged weights](https://ciaaw.org/abridged-atomic-weights.htm),
[CIAAW isotope compositions](https://ciaaw.org/isotopic-abundances.htm),
[AME2020 atomic masses](https://amdc.impcas.ac.cn/masstables/Ame2020/mass_1.mas20),
and [CODATA 2022 constants](https://physics.nist.gov/cuu/pdf/wall_2022.pdf).
Evidence under `target/benchmark-parity/`: `audit-descriptor-results.py`,
`descriptor-results-audit.json`, `audit-descriptor-masses.py`,
`descriptor-mass-audit.json`, the downloaded `atomic-data-sources/`, and
`regenerated-descriptor-data.rs`. The mass audit records every atomic constant
used in the corpus and representative charged/isotope-labeled observations.

### Full rerun and retained limitations

Before and after commands:
`target/mmcif-validation/release/kekule-bench.exe --feature descriptor.molecular --dataset all --output target/benchmark-parity/descriptor-molecular-{before,after}.json`
(each suffix was run separately). Both completed and correctly exited 1 for
remaining numerical disagreements and errors. The same executable is appropriate
because this turn changes generator tooling and documentation, not runtime code.

| Dataset | Applicable | Agree | Disagree | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 200,000 | 4 | 199,942 | 54 |
| Enamine | 100,480 | 0 | 97,599 | 2,881 |
| PL-Rex | 328 | 0 | 328 | 0 |
| Smoke | 25 | 0 | 25 | 0 |
| Total | 300,833 | 4 | 297,894 | 2,935 |

Another 1,001 selected cases are not applicable. All 301,834 original-to-before
rows are byte-identical, so previous feature fixes did not alter these results.

Errors comprise 2,881 unsupported CXSMILES inputs, 18 previously documented
strict-valence failures, and 36 missing standard-weight cases: Pm (2), Ac (7),
Fm (2), Pu (6), Tc (6), Np (2), Cf (4), Lr (4), and No (3). CIAAW does not supply
standard atomic weights for these elements. The native API deliberately reports
that absence instead of treating a representative isotope mass number as an
average atomic weight. Explicit isotope labels can supply a defined isotope
mass; existing regressions cover that distinction and unavailable isotope data.
These limitations remain visible in the benchmark.

The final paired audit confirms that the complete before/after case files are
byte-identical (SHA-256 `7027b5609845a9445393af23862af14c66dc6a0eda200520052b17b78dc13a6a`).
Both reports retain identical source selection, reference identities, contract
and executable fingerprints. Evidence: `verify-descriptor-original-rerun.py`,
`descriptor-original-rerun-verification.json`, `verify-descriptor-rerun.py` and
`descriptor-rerun-verification.json` in the local analysis directory.

### Validation

Passed: `cargo fmt --all -- --check`; workspace all-target/all-feature check;
clippy with warnings denied; workspace all-feature tests including doctests;
workspace Rustdoc with warnings denied; and the `kekule-bench` package file-list
check. Rust commands used `--locked --offline` where applicable and
`--target-dir target/mmcif-validation`; packaging used `--allow-dirty`.
Also passed: both new generator regressions, source checksum verification,
offline regeneration with exact table-text comparison, all 28 pinned RDKit
regressions, eight shared runner regressions, the two full descriptor runs and
their observation audits. Benchmark exit 1 represents the documented remaining
differences, not a failed runner. `git diff --check` passes.

Logs are `descriptor-*.log` under `target/benchmark-parity/`. Rust 1.89 checks,
potentials no-default-feature tests/docs, full `kekule` package verification,
companion package file-list checks, packaged-license comparisons and dashboard
tests were not repeated: they passed in feature 19, and this turn changes only
the Python generator, its CI regression and documentation. Full companion
package builds remain deferred as in CI until the foundational crate is
published. Linux fuzz checks were not run on this Windows host. Biopython/DSSP
and unrelated feature benchmarks were not repeated. No golden was regenerated.

Next feature: `descriptor.rotatable-bonds.rdkit-strict`.

## 21. Strict rotatable bonds

### Fix

The resonance filter incorrectly required a neutral doubly bonded N/O/S atom,
with a separate narrower exception for cationic nitrogen. RDKit's strict
`[N,O,S]` predicate constrains element and aliphaticity, not formal charge.
Charged imidates, amidines and related conjugated linkages therefore lost their
resonance exclusion in the native detector. The filter now accepts the charged
double-bond partner under the same general rule and removes the redundant
cation-specific branch.

The regression fails before the fix on `CC(=[NH2+])OC` and passes afterward.
It covers neutral, positive and negative nitrogen, O/S/N linkages, both encoded
and materialized hydrogen representations, the option to include restricted
bonds, and preservation of the following unrestricted alkyl bond. All 15 focused
rotatable-bond tests pass. No special source IDs, element combinations beyond
the existing strict definition, or fixture-specific exclusions were introduced.

The reference definition was checked against RDKit's
[Lipinski.cpp](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Descriptors/Lipinski.cpp)
and independent probes of the pinned RDKit 2026.03.3 query and descriptor APIs.
The native detector remains read-only and uses represented connectivity without
mutating or installing perception.

### Benchmark audit and hydrogen policy

Both the complete endpoint set and descriptor count remain asserted. The
reference explicitly selects `Strict`, enumerates its source SMARTS without a
match cap, and fails if that set's count differs from the native RDKit count.
It is not the separate `StrictLinkages` option. All 44,196 component observations
in the 36,919 initially disagreeing cases were regenerated from the original
external fixtures with pinned RDKit and exactly reproduced their stored
reference records, including titles and component indices. The reference
count-versus-endpoint invariant also passed throughout this audit.

Of those cases, 36,797 initially match Kekule exactly after graph hydrogens are
removed from an analysis-only copy of the reference molecule. Original heavy
atom indices are retained as atom properties and mapped back before comparing
endpoints; the audit does not replace atom correspondence with bond counts.
This experiment is not a benchmark-adapter or golden change. RDKit's query uses
graph degree and hydrogen-count predicates, so materializing H can change its
terminal, symmetric-group and resonance classification. Kekule deliberately
uses heavy-atom connectivity. That representation invariance is appropriate
for a molecular flexibility descriptor and is retained rather than copying
the reference's hydrogen-representation dependence.

The guide now explains these policies and the native resonance approximation.
No assertion, reference observation, schema or comparison rule was weakened,
and no golden was regenerated. Evidence under `target/benchmark-parity/`:
`scan-rotatable.py`, `rotatable-summary.json`, `rotatable-disagreements.json`,
`audit-rotatable.py` and `rotatable-source-audit.json`.

### Full rerun

Before:
`target/mmcif-validation/release/kekule-bench.exe --feature descriptor.rotatable-bonds.rdkit-strict --dataset all --output target/benchmark-parity/rotatable-before.json`.
After:
`cargo run --release --locked --offline --target-dir target/mmcif-validation -p kekule-bench -- --feature descriptor.rotatable-bonds.rdkit-strict --dataset all --output target/benchmark-parity/rotatable-after.json`.
Both runs completed and correctly exited 1 for remaining disagreements/errors.

| Dataset | Applicable | Agree before → after | Disagree before → after | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 200,000 | 163,298 → 163,332 | 36,684 → 36,650 | 18 |
| Enamine | 100,480 | 97,599 → 97,599 | 0 → 0 | 2,881 |
| PL-Rex | 328 | 94 → 94 | 234 → 234 | 0 |
| Smoke | 25 | 24 → 24 | 1 → 1 | 0 |
| Total | 300,833 | 261,015 → 261,049 | 36,919 → 36,885 | 2,899 |

Another 1,001 selected cases are not applicable. The 18 PubChem errors are
reference sanitization failures on the previously identified valence-invalid
records; this structural detector itself does not run strict valence perception.
The 2,881 Enamine errors are unsupported CXSMILES inputs.

The original-to-before audit verifies all 301,834 rows are byte-identical.
The before/after audit verifies 301,798 rows remain byte-identical. Exactly 36
cases lose 48 incorrect axes; every removed axis is absent from the unchanged
reference endpoint set. Thirty-four cases become exact agreements, and two
retain independent hydrogen-representation differences. No axes were added,
no error or new disagreement appeared, and all other observation fields,
source identities and reference values remain unchanged. Report audits verify
complete runs and the rebuilt executable fingerprint.

Evidence: `verify-rotatable-original-rerun.py`,
`rotatable-original-rerun-verification.json`, `verify-rotatable-rerun.py` and
`rotatable-rerun-verification.json` in the local analysis directory.

### Deferred aromaticity integration

After the charge fix, 36,799 of the remaining disagreements are explained
entirely by hydrogen representation. The other 86 cases retain a difference
even against the analysis-only hydrogen-collapsed reference:

| Cause | Cases | Axes |
| --- | ---: | ---: |
| A nonaromatic conjugated ring is treated as aromatic by the localized approximation | 84 | 90 incorrectly included |
| A selenium-containing aromatic ring is not recognized by that approximation | 2 | 2 incorrectly excluded |

Every extra axis has both endpoints excluded by the reference resonance query
and a participating nonaromatic cyclic C=N/O/S double bond. The missing axes
occur in the SMILES and SDF representations of PubChem 356034, with a
selenium-containing aromatic ring. These classifications are checked against
the original structures in `classify-rotatable-residuals.py` and
`rotatable-residual-classification.json`.

The current detector approximates aromaticity using local conjugation in five-
or six-member cycles. Conjugation alone does not establish aromaticity, and the
limited element/ring rules do not cover the full model. Adding selenium or
individual exocyclic-bond exceptions would leave the underlying defect in place.
A clean correction should reuse the existing aromaticity implementation and
define how its fallible, resource-bounded perception fits this currently
infallible, read-only detector API. It should also preserve results regardless
of whether perception was previously installed. That integration is deferred
as a broader API/design change; Rustdoc now explicitly describes the limitation.

### Validation

Passed: `cargo fmt --all -- --check`; workspace all-target/all-feature check;
clippy with warnings denied; workspace all-feature tests including doctests;
workspace Rustdoc with warnings denied; Rust 1.89 workspace all-target/all-feature
check; potentials tests and docs without default features; full `kekule` package
build and verification; companion and benchmark package file-list checks; and
packaged-license hash comparisons. All 15 focused rotatable-bond tests pass,
including the new regression that failed before the fix. The 28 pinned RDKit
tests and eight shared runner tests also pass. Full benchmark reruns and paired
source/observation audits passed their invariants; benchmark exit 1 retains the
documented differences. `git diff --check` passes.

Rust commands used `--locked --offline` where applicable and
`--target-dir target/mmcif-validation`, except the separate Rust 1.89 check.
Package commands used `--allow-dirty`. Logs are `rotatable-*.log` under
`target/benchmark-parity/`. Full companion package builds were not run because
CI checks their file lists until the foundational crate is published. Linux
fuzz build/smoke checks were not run on this Windows host. Dashboard tests,
atomic-data generator tests, Biopython/DSSP and unrelated feature benchmarks
were not repeated because their implementations were unchanged in this turn.
No golden regeneration was performed.

Next feature: `stereo.representation`.

## 22. Stereo representation

### Reference correction

This feature's Enamine reference archive still contained the private-property
omission repaired in sections 5 and 8. All 100,480 Enamine observations, spanning
both SDF and SMILES, were independently regenerated with RDKit 2026.03.3 into
`target/benchmark-parity/stereo-representation-corrected-goldens/`.
The exhaustive comparison proved that exactly 2,881 SDF records regain their
source `_CXSMILES_Data` field. All chemical observations, previously asserted
fields, source digests, component indices, titles and SMILES records are unchanged.
The corrected pair was promoted only after that verification; the original
archive and manifest remain under `stereo-representation-original-goldens/` in
the local analysis directory. The existing private-property regression now
explicitly exercises `stereo.representation` as well as MOL and SDF parsing.

The generation command used the validated release benchmark executable with
`generate --feature stereo.representation --dataset enamine-diversity`, the pinned
RDKit Python path, and the separate golden directory. It completed all 100,480
outcomes without a reference error. No runtime chemistry code, reference
algorithm, comparison field, tolerance or schema was changed this turn.

### Stereo and source audit

The feature compares the complete prepared graph, not an isolated stereo score.
The guide now makes explicit that hydrogen declarations, inferred radicals,
coordinates and source properties accompany stereo focuses, normalized carriers,
orientations and enhanced groups. RDKit performs its default stereo cleanup;
Kekule's default perception does not rewrite represented stereo. Source readers
also differ in how drawing coordinates are interpreted. RDKit's
[documented Molfile policy](https://www.rdkit.org/docs/RDKit_Book.html#interpretation-of-the-2d-3d-flag-and-stereochemistry)
includes automatic coordinate-based interpretation; Kekule's source normalization
decodes format-local marks and leaves general coordinate materialization to the
separate operation measured by `stereo.perception`.

The original run contains 6,113 cases with stereo-field disagreements. A
focus-keyed comparison identifies the following differences without confusing
array reordering with changed stereocenters:

| Difference | Elements |
| --- | ---: |
| Additional specified double-bond elements in Kekule | 7,542 |
| Reference double-bond elements absent in Kekule | 64 |
| Different tetrahedral orientation at the same focus and carriers | 102 |
| Additional specified tetrahedral elements in Kekule | 80 |
| Axis instead of a tetrahedral element | 3 axes / 3 missing tetrahedra |

These occur in 5,930 cases with only extra double-bond elements, 31 with only
missing double-bond elements, 58 with only tetrahedral orientation differences,
71 with only extra tetrahedral elements, 12 with extra double-bond plus orientation
differences, eight with extra double-bond plus extra tetrahedral elements, and
three axis/tetrahedral substitutions. No additional carrier-convention or
enhanced-group mismatch category remains after aligning by type and focus.
Only one extra double-bond element comes from a SMILES input; all other listed
stereo differences come from SDF.

An independent source audit reparses every one of these 6,113 external inputs
and reproduces all stored reference stereo/group observations across 7,259
components. It also interrogates RDKit's separate `FindPotentialStereo` operation:
7,421 extra double-bond focuses and five extra tetrahedral focuses are absent
from its potential set, while 121 extra double-bond focuses and 75 extra
tetrahedral focuses are retained as potential centers. Thus absence from the
default cleaned representation does not uniformly mean that the native focus
is chemically impossible. These counts are diagnostic and do not replace any
benchmark assertion with candidate-count agreement.

### Remaining scope

- Stereo cleanup remains the main issue. Many extra double bonds have equivalent
  substituents, as in PubChem 32. Removing redundant or invalid assertions safely
  requires a symmetry-aware policy that handles stereo-dependent equivalence;
  blanket deletion or a local neighbor comparison is insufficient. The differing
  results of RDKit's cleanup and candidate APIs are retained in the audit.
- The 102 orientation differences occur in crowded wedge drawings; PubChem
  10524 has two examples. The existing decoder uses displaced-coordinate volume,
  whereas the reference handles angular ordering and ambiguous/collinear drawing
  configurations together. No global parity inversion is justified.
- Three sulfoxide examples (146091, 461502 and 461520) are interpreted as axes
  instead of pyramidal tetrahedral centers. Correcting this needs a reliable
  shared endpoint-hybridization/conjugation rule, not a sulfur-specific ban.
- Sixty-three missing reference double-bond elements are endocyclic C=N or N=N
  bonds in rings of at least eight atoms. The source interpreter explicitly
  excludes ring double bonds with a noncarbon endpoint, and the candidate layer
  has the same exclusion. Existing tests encode that boundary. Extending it
  requires reconciling source assertions, candidate detection and downstream
  behavior; the restriction was not removed in just one layer to improve this
  feature's score.
- The remaining missing double-bond observation is PubChem 372730's P=N bond.
  RDKit reports `STEREOANY`, not a specified E/Z configuration. The P endpoint
  has three other neighbors and is SP3 in RDKit; the N-side drawing is collinear.
  This exceeds the native double-bond carrier invariant, and RDKit's candidate
  finder does not identify that bond as potentially stereogenic. The native
  representation was not weakened to create this unsupported unknown element.

Evidence under `target/benchmark-parity/`: `scan-stereo-representation.py`,
`stereo-representation-original-summary.json`,
`stereo-representation-original-disagreements.json`,
`classify-stereo-representation.py`, `stereo-representation-classification.json`,
`audit-stereo-representation-sources.py`, `stereo-representation-source-audit.json`,
`stereo-missing-double-audit.json`, `stereo-missing-double-ring-audit.json`,
and the original P=N source/probe in
`stereo-pn-source.sdf` and `stereo-pn-audit.json`.

### Rerun

The original-to-before audit verifies all 301,834 rows: 301,830 are byte-identical,
and the other four contain only the already verified aromatic triple-bond
observation correction for PubChem 141133/141134 in both source formats.
Stereo fields and every stored reference are unchanged by those upstream fixes.
The final rerun against the corrected metadata reference has completed:

| Dataset | Applicable | Agree before → after | Disagree before → after | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 200,000 | 186,510 → 186,510 | 13,472 → 13,472 | 18 |
| Enamine | 100,480 | 69,690 → 69,699 | 27,909 → 27,900 | 2,881 |
| PL-Rex | 328 | 290 → 290 | 38 → 38 | 0 |
| Smoke | 25 | 20 → 20 | 5 → 5 | 0 |
| Total | 300,833 | 256,510 → 256,519 | 41,424 → 41,415 | 2,899 |

Another 1,001 selected cases are not applicable. Both before and after runs used
`target/mmcif-validation/release/kekule-bench.exe --feature stereo.representation --dataset all`
with separate `stereo-representation-before.json` / `stereo-representation-after.json`
output paths under `target/benchmark-parity/`. Both completed and correctly
exited 1 for remaining disagreements/errors. The 18 strict-valence failures
and 2,881 unsupported CXSMILES cases remain unchanged.

The final paired audit checks all 301,834 rows and proves that every native
observation is unchanged. Only the 2,881 verified reference property lists change;
nine cases move to complete agreement and no case regresses. Both reports retain
identical executable, adapter and contract fingerprints, input selection and
error counts. Only the audited Enamine golden identity changes.

The remaining 41,415 disagreements have these overlapping causes: stereo in
6,113 cases, hydrogen storage in 27,593, hydrogen-inference policy in 8,245,
radicals in 938, and the previously documented radical/aromaticity dependency
in one. Every hydrogen-storage delta preserves total H and the corresponding
represented-valence relationship. All 674,454 coordinate differences recorded
within these disagreeing cases are inside the existing numerical allowance;
none is an additional unexplained geometry mismatch. No other difference
category remains. Evidence: `verify-stereo-representation-rerun.py`,
`stereo-representation-rerun-verification.json`,
`verify-stereo-representation-goldens.py`, and
`stereo-representation-golden-verification.json` in the local analysis directory.

### Validation

Passed: workspace formatting, all-target/all-feature check, clippy with warnings
denied, all-feature workspace tests including doctests, Rustdoc with warnings
denied, and the benchmark package file-list check. All 28 pinned RDKit tests,
including the expanded private-property regression, and eight shared runner
tests pass. The independent golden and source-stereo audits described above
also pass. `git diff --check` passes.

Rust commands used `--locked --offline` where applicable and
`--target-dir target/mmcif-validation`; the package list used `--allow-dirty`.
Logs are `stereo-representation-*.log` under `target/benchmark-parity/`.
Rust 1.89, potentials no-default-feature tests/docs, full runtime packaging,
companion package lists and license-copy checks were not repeated: they passed
in feature 21 and no Rust implementation, dependency or build configuration
changed here. Full companion package builds remain deferred as in CI until the
foundational crate is published. Linux fuzz checks were not run on this Windows
host. Dashboard and atomic-data tests, Biopython/DSSP and unrelated feature
benchmarks were not repeated because their code is unchanged in this turn.

Next feature: `stereo.perception`.

## 23. Stereo perception

### Scope and independently checked reference behavior

This feature compares the complete prepared graph, materialized coordinate
stereo and candidate identities. It inherits the hydrogen, radical and source
interpretation differences described in features 18 and 22. Native candidate
detection is a local eligibility test; RDKit's `FindPotentialStereo` additionally
resolves ligand equivalence and dependencies between stereo sites. Neither a
candidate count nor absence of a CIP label establishes equivalent behavior.
The public candidate API and benchmark guide now state this distinction.

The original report contains 169,088 cases with candidate or stereo differences.
Every one was reread from its externally supplied source and evaluated with
pinned RDKit 2026.03.3. All 175,100 affected components reproduce their stored
reference candidate sets and stereo observations. The audit matches both fields
together, including multiplicity for repeated component signatures. Evidence:
`audit-stereo-perception-sources.py`, `stereo-perception-source-audit.json`,
`scan-stereo-perception.py` and `stereo-perception-classification.json` under
`target/benchmark-parity/`.

RDKit's coordinate call keeps its default replacement policy for existing tags
on 3D conformers; native materialization preserves represented assertions.
A pinned probe confirms that the reference call leaves tags unchanged on a 2D
conformer. The benchmark does not mistake that no-op for coordinate inference.
These preparation and supported-family differences remain explicit in the
comparison rather than being normalized away.

### General runtime fixes

Tetrahedral and double-bond candidates now reject repeated terminal hydrogen
ligands. Previously, a center with two encoded hydrogens was excluded, but the
same center became a candidate after those hydrogens were expanded into graph
atoms. Even methane carbons could then be proposed as tetrahedral stereo sites.
The same representation defect affected terminal alkene endpoints. The rule
compares isotope, charge and radical state, and treats an implicit ordinary H
like an ordinary terminal graph H; atom-map labels do not distinguish ligands.
Bridging hydrogens are not assumed equivalent to terminal hydrogen ligands.

Distinct isotope ligands remain eligible. Focused regressions cover repeated H,
repeated D, H/D and D/T tetrahedral sites, mixed graph/implicit H and terminal
alkenes. The pinned reference agrees on the tetrahedral cases and repeated-H
alkenes but omits the H/D alkene `[H]C([2H])=CF`. Native eligibility is retained
there: both alkene endpoints have chemically distinct ligands, so suppressing
the site merely because two ligands are hydrogen would lose isotope-dependent
geometric stereochemistry. This is a local capability check, not a claim of
complete isotope-sensitive stereo perception. `stereo-perception-probes.json`
records the independent reference observations.

Double-bond candidates also enforce at most two substituents at either endpoint,
counting all graph neighbors and checking the resulting carrier-list size.
Previously an overcoordinated P=N or related center could become a coordinate
proposal which failed the core stereo invariant when materialized, failing the
entire case. Candidate detection now respects that invariant before inference.
A focused P=N regression checks both eligibility and successful transactional
materialization. Both new regression tests fail against the original runtime
and pass with the fixes. Existing coordinate, source-stereo and CIP regressions
continue to pass. No represented assertion is deleted by these changes.

### Reference metadata correction

The Enamine reference archive was independently regenerated for all 100,480
cases. An exhaustive before/after comparison found exactly 2,881 additions of
the private source field `_CXSMILES_Data`; all other fields, candidate sets,
stereo, chemistry and source identities are unchanged. The remaining 97,599
observations are identical. The original golden and manifest are retained in
`stereo-perception-original-goldens/enamine-diversity/`, and the corrected pair
was promoted only after this audit. The private-property regression now covers
`stereo.perception` explicitly. No comparison field or tolerance was weakened.

The corrected compressed archive is 36,682,593 bytes, SHA-256
`b1fd0284f9ac584c73213b151581e26548f8bdb3806a41736b38123caae093e5`.
Its manifest retains the locked source identity and records RDKit 2026.03.3 and
the current adapter fingerprint. `verify-stereo-perception-goldens.py` and
`stereo-perception-golden-verification.json` preserve the full metadata audit.

### Full rerun and remaining work

The rebuilt release benchmark ran `--feature stereo.perception --dataset all`
before and after the fixes, producing `stereo-perception-before.json` and
`stereo-perception-after.json` with their full case archives. The completed
comparison exits 1 for the remaining disagreements, as expected.

| Dataset | Cases | Agree before → after | Disagree after | Errors before → after |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 200,000 | 72,717 → 130,808 | 69,174 | 85 → 18 |
| Enamine | 100,480 | 33,769 → 33,788 | 63,811 | 2,921 → 2,881 |
| PL-Rex | 328 | 56 → 198 | 130 | 0 → 0 |
| Smoke | 25 | 11 → 14 | 11 | 0 → 0 |

Total agreement rises from 106,553 to 164,808, with 133,126 disagreements and
2,899 errors remaining. All 107 coordinate-materialization errors are resolved;
the remaining errors are the previously documented 18 strict-valence and 2,881
unsupported CXSMILES cases. The 1,001 not-applicable IDs remain unchanged.

The paired audit checks every one of the 301,834 rows. All 106,553 previously
agreeing cases still agree; 58,205 disagreements and 50 errors become complete
agreements, and the other 57 repaired errors become reportable disagreements.
Among previously successful native outcomes, 89,597 change: exactly 538,586
ineligible candidates and 1,748 inferred stereo elements are removed. The audit
independently applies the hydrogen/coordination rules to the observed graphs,
checks every remaining candidate and stereo element, and verifies all other
native graph fields are unchanged. The only reference changes are the 2,881
audited private fields. No source identity or not-applicable status changes.
Evidence: `verify-stereo-perception-rerun.py` and
`stereo-perception-rerun-verification.json`.

Remaining differences have these overlapping causes: candidate detection in
108,493 cases, stereo observations in 6,205, hydrogen storage in 27,593,
hydrogen-inference policy in 8,245, radicals in 938 and the known radical-related
aromaticity dependency in one. All hydrogen-storage differences preserve total
H and the expected represented-valence relationship. Every reported coordinate
difference is within the existing numerical tolerance; no unexplained field
category remains.

The original candidate audit includes 533,364 extra carbon tetrahedral sites
with repeated hydrogen ligands. General equivalence remains substantial even
after this bounded correction: the original audit also found 124,750 extra
carbon tetrahedral sites without repeated H, alongside nitrogen eligibility,
heteroatom lone-pair and coordination questions. Original missing tetrahedral
sites include 4,187 phosphorus and 245 arsenic centers without repeated H.
Additional gaps include endocyclic heteroatom double bonds, non-tetrahedral
coordination stereo and atropisomer candidates. Completing these requires a
consistent symmetry/dependency and stereo-family policy shared with coordinate
inference and CIP; adding a few element-specific exceptions would leave that
contract incomplete. These are deferred, together with the source orientation,
cleanup and hydrogen/radical differences already documented in earlier reviews.

### Validation

Passed: 68 focused stereo tests; workspace formatting, all-target/all-feature
check, clippy with warnings denied, all-feature workspace tests including
doctests, Rustdoc with warnings denied, Rust 1.89 all-target/all-feature check,
potentials no-default-feature tests and docs, full `kekule` package verification,
companion and benchmark package file lists, and all packaged license comparisons.
All 28 pinned RDKit adapter tests and eight shared runner tests pass.
`git diff --check` passes.

Rust commands used `--locked --offline` where applicable and
`--target-dir target/mmcif-validation`, except the MSRV check used its established
default target. Package checks used `--allow-dirty`; documentation used
`RUSTDOCFLAGS=-D warnings`. Logs are `stereo-perception-*.log` in the local
analysis directory. Full companion package builds remain deferred as in CI
until the foundational crate is published. Linux fuzz compilation/smoke runs
were not run on this Windows host. Dashboard, atomic-data and Biopython/DSSP
tests and unrelated feature benchmarks were not repeated because their code
is unchanged in this turn.

Next feature: `stereo.cip`.

## 24. CIP assignment

### Benchmark correction

The native adapter serialized `StereoDescriptor::SeqCis` and `SeqTrans` using
their Rust variant names, while RDKit's CIPLabeler renders the corresponding
descriptors as lowercase `z` and `e`. Its
[descriptor conversion](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/CIPLabeler/Descriptor.h)
explicitly distinguishes these from ordinary uppercase `Z` and `E`.
The native benchmark adapter now uses the same spelling. The runtime enum,
sequence rules, represented stereo and descriptor assignments are unchanged.
The comparator still distinguishes lowercase from uppercase labels; no label,
focus, count, observation or tolerance is removed or weakened.

A regression exercises both configurations of the existing sequence-cis/trans
validation scaffold and explicitly rejects replacing the lowercase result with
its uppercase counterpart. It fails against the original adapter and passes
after the correction. This resolves PubChem SDF cases 385445 (`seqCis` versus
`z`) and 403847 (`seqTrans` versus `e`) by a general descriptor conversion,
without source-ID checks or golden regeneration.

### Source and ranking audit

The fresh before-run is byte-for-byte identical to the user's original CIP
results in all 301,834 rows. Changes made during earlier feature reviews did
not alter this feature. Pinned RDKit 2026.03.3 independently reproduces every
reference observation in the 302 disagreeing cases and the 104 native
CIP-ranking error cases. Clearing existing atom and bond `_CIPCode` properties
before rerunning `AssignCIPLabels` produces the same labels throughout these
406 cases; stale legacy labels do not explain the differences. CIP observations
contain no source-property fields, so the private-SDF metadata correction from
other features does not require a CIP golden update.

After the spelling correction, all 300 disagreeing cases have corresponding
source-stereo differences already audited in feature 22. Their disjoint case
profiles are: 126 with extra bond descriptors (451 descriptor observations),
30 with missing bond descriptors (63 observations), 70 with opposite atom
descriptors (102 observations), 71 with extra atom descriptors (72 observations),
and three with a native axis label in place of the reference sulfoxide atom
descriptor. The 102 opposite atom descriptors have opposite represented
tetrahedral parity with the same carriers, rather than an isolated ranking
disagreement. Missing bonds are the previously documented endocyclic heteroatom
stereo exclusions. This audit preserves disconnected-component offsets when
matching source stereo to global descriptor focuses.

An additional analysis-only experiment supplies the native represented stereo
to the corresponding prepared RDKit graphs without changing production adapters
or goldens. RDKit then reproduces the native labels in 199 of the 300 cases.
The other 101 contain assertions that cannot be transferred directly: 391
localized aromatic double-bond tags and three sulfoxide axis tags. RDKit's
labeler rejects a double-bond stereo assertion on a bond whose prepared type
is aromatic. Omitting these unsupported diagnostic assertions leaves identical
atom labels and identical labels for every common bond in all 101 cases; only
native extra bond labels remain. This is evidence about the layer responsible,
not permission to suppress those disagreements in the actual benchmark.

The remaining source-wedge interpretation, symmetry cleanup, ring eligibility,
aromatic represented-stereo and sulfoxide family issues therefore remain as
described in features 22–23. They need consistent source/perception policy,
not local inversions or descriptor substitutions in the CIP ranker.

Native defaults limit rooted expansion to depth 32 and 100,000 nodes; all 104
CIP assignment errors report depth exhaustion, not a proven ligand tie.
RDKit uses its existing 1,000,000-recursive-iteration bound. The guide now states
both limits and explains that the units are not equivalent computational
budgets. Default resource exhaustion stays visible as an error; no unsupported
assignment is silently skipped and no default limit is raised to improve the
score.

An optional standalone diagnostic raised only depth to 64, keeping the
100,000-node bound. It completed 56 of the 104 error cases before being
intentionally stopped after roughly 416 CPU seconds: 51 then agree, three
produce descriptor differences at already differing represented stereo focuses,
and two (PubChem 158374 and 163705) hit the node limit instead. The other 48
cases have no raised-depth result; this partial experiment is not a full
coverage or timing comparison. It demonstrates that a larger cap changes
computational cost and does not eliminate the underlying source/resource
questions. Production defaults and benchmark outcomes remain unchanged.
`stereo-cip-depth64-summary.json` records the explicit stop and incomplete
coverage; `stereo-cip-depth64-disagreement-audit.json` verifies the three exposed
source-stereo differences. The probe uses supplied corpus records, not new
benchmark fixtures. Improving expansion efficiency or changing public resource
defaults is deferred rather than treated as a small parity repair.

Audit evidence under `target/benchmark-parity/` includes
`verify-stereo-cip-original.py`, `stereo-cip-original-verification.json`,
`classify-stereo-cip.py`, `stereo-cip-classification.json`,
`audit-stereo-cip-sources.py`, `stereo-cip-source-audit.json`,
`audit-stereo-cip-represented.py` and `stereo-cip-represented-audit.json`.
The supplied error sources and the analysis-only depth probe also remain local.

### Full rerun

The rebuilt release executable ran `--feature stereo.cip --dataset all` after
the correction. Both full comparisons completed normally and correctly exited
1 for the remaining differences.

| Dataset | Cases | Agree before → after | Disagree after | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 200,000 | 199,576 → 199,578 | 300 | 122 |
| Enamine | 100,480 | 97,599 → 97,599 | 0 | 2,881 |
| PL-Rex | 328 | 328 → 328 | 0 | 0 |
| Smoke | 25 | 25 → 25 | 0 | 0 |

Total agreement rises from 297,528 to 297,530. The 3,003 errors comprise the
104 explicit depth-limit failures, 18 previously documented strict-valence
failures and 2,881 unsupported CXSMILES records. The 1,001 not-applicable IDs
remain unchanged. The before/after verifier checks all 301,834 rows: exactly
the two descriptor spellings change, every reference observation and golden
identity is unchanged, every other complete case row is identical, and all
previous agreements remain agreements. Reports and evidence are
`stereo-cip-before.json`, `stereo-cip-after.json`,
`verify-stereo-cip-rerun.py` and `stereo-cip-rerun-verification.json`.

### Validation

Passed: the new focused adapter regression, workspace formatting,
all-target/all-feature check, clippy with warnings denied, all-feature workspace
tests including doctests, Rustdoc with warnings denied, benchmark package file
list, 28 pinned RDKit tests, eight shared runner tests and `git diff --check`.
Rust commands used `--locked --offline` and
`--target-dir target/mmcif-validation`; the package list used `--allow-dirty`.
Logs are `stereo-cip-*.log` under the local analysis directory.

Rust 1.89, potentials no-default-feature tests/docs, full runtime packaging,
companion package lists and license-copy checks were not repeated: they passed
in feature 23 and no runtime, dependency or build configuration changed here.
Full companion package builds remain deferred as in CI until the foundational
crate is published. Linux fuzz compilation/smoke runs were not run on this
Windows host. Dashboard, atomic-data and Biopython/DSSP tests and unrelated
feature benchmarks were not repeated because their code is unchanged.

Next feature: `bio.secondary-structure.dssp`.

## 25. DSSP secondary structure

The numerical audit follows the pinned
[DSSP 4.6.1 implementation](https://github.com/PDB-REDO/dssp/blob/v4.6.1/libdssp/src/dssp.cpp).
Its private points, vector operations and angular calculations use `float`.
Kekule previously used double precision for torsions, TCO and bend angles while
already using single precision for distances and reconstructed hydrogens. This
crossed reference output rounding boundaries. The private DSSP geometry now
uses single precision consistently, without changing model coordinates or the
public `f64` result types. Undefined or overflowing angular intermediates return
no angle instead of a non-finite observation. Hydrogen-bond energy arithmetic,
conformer selection and secondary-structure rules are unchanged.

A regression from supplied PDB 1A2L, chain B residue 67, fails before the change:
the old phi rounds to -51.9, while DSSP 4.6.1 reports -52.0. The corrected
implementation passes, as do the sign-convention and overflow regressions.

The reference adapter computes omega independently using Biopython's
double-precision vectors constructed from single-precision parsed coordinates.
The comparator now allows `16 * f32::EPSILON * 180` degrees for that arithmetic
boundary, in addition to its existing double-precision conversion allowance.
This is based on precision and angular range, not a fit to corpus errors. The
regression failed before the change and now accepts arithmetic roundoff while
still rejecting a 0.01-degree geometry change and applying no such allowance to
other features. All raw values remain in the reports. No other tolerance changes.

The benchmark schema also incorrectly required a numeric polymer sequence label
on every residue. Thirteen supplied structures contain 39 native eligible
residues with an author sequence number but no polymer label. A null polymer
label is now valid; omitting the field remains an error, and null versus a
number remains a structural disagreement. The regression demonstrates the old
schema error and verifies all three distinctions. Residues and partner identities
are retained, not filtered to match the reference.

### Residuals and reference scope

The numerical rerun has 183 disagreements. Of those, 181 supplied structures
contain alternative backbone coordinates. Inspection of pinned DSSP source
shows that `addAtom` overwrites backbone points as alternate atom rows arrive,
whereas Kekule analyzes its model's selected conformer. Biopython's independently
computed omega uses selected atoms as well, so the reference observation can
combine measurements from different conformers. The guide now makes this
limitation explicit: these end-to-end observations do not demonstrate algorithm
parity on identical coordinates. Backbone alternatives are evidence of differing
inputs, not proof that they explain every remaining field. A controlled common-
conformer benchmark needs explicit fixture/interpretation policy; changing model
selection or stripping alternatives just to improve this score is deferred.

The two structures without alternate backbone coordinates are 1HTR (chain-break
classification) and 6IY2 (weak hydrogen-bond slots, energies and partner identities).
DSSP determines breaks by peptide distance, even across label-chain changes;
Kekule preserves chain boundaries. The weak-bond differences require a focused
residue-order/reconstructed-hydrogen audit before an algorithm change is justified.
Neither receives a case-specific adjustment.

The thirteen schema-error cases expose a separate eligibility difference:
Kekule's documented analysis accepts complete backbones without polymer labels,
including free amino-acid residues. DSSP builds its residue table from
`_pdbx_poly_seq_scheme`. These observations should disagree explicitly rather
than fail JSON validation. Four other inputs (5J3E, 6C1U, 6C1V and 6D1T) contain
unknown atom element `X` rejected by model interpretation; supporting unknown
elements is a broader model decision. Reference process/metadata failures and
inputs without analyzable protein remain errors, never agreements.

Local evidence includes `dssp-original-summary.json`, `dssp-classification.json`,
`dssp-after-audit.json`, `dssp-conformer-audit.json`,
`dssp-observation-errors.json` and their scripts. Difference paths refer to the
comparator's normalized residue order; the early structural-detail helper's
raw-array examples must not be used to identify residues by that index.

The first default-concurrency run aborted on memory allocation. A later serial
native run and concurrent reference generation also aborted under memory pressure;
neither incomplete report is counted as a full result. A subsequent run was
explicitly stopped when superseded by the batching fix. The runner previously
accumulated 256 structures without considering their source sizes. It now flushes
at 256 cases or 8 MiB of accumulated input, whichever comes first, with regression
coverage for large inputs and the retained small-case limit. An individual large
input can exceed this threshold; it is not a hard process-memory cap. This applies
equally to generation and comparison, without changing observations or timing
boundaries. Full native reruns use `--jobs 1`, now documented for constrained hosts.

The original stored DSSP golden payloads have not been replaced. The incomplete
independent reference generation is retained in
`target/benchmark-parity/dssp-audit-goldens-verified`.
The first 246 complete records reproduce all structural observations; differences
are omega roundoff of at most 2.85e-14 degrees and temporary paths in two error
messages. This partial check is not full reference coverage.

The retry with bounded batches subsequently completed all 1,000 PDB inputs and
all 20 smoke selections (one applicable structure and 19 not-applicable IDs).
The exhaustive independent audit finds 297 PDB outcomes identical and 695
differing only in 11,408 omega values by at most 2.842170943040401e-14 degrees.
Every such difference fits the original double-precision conversion allowance,
without using the new single-precision allowance. Eight error messages differ
only in temporary input filenames; all other fields, identities and outcomes
match. All 20 smoke records are identical. The same 244 PDB reference failures
remain; generation correctly exits 1 for them, while smoke generation exits 0.
No stored golden payload was replaced. Evidence is
`dssp-bounded-reference-audit.json`, `dssp-bounded-smoke-reference-audit.json`
and the two `dssp-bounded-generation-*.json` reports.

The shared `contract.json` now records the omega precision boundary and nullable
polymer labels. All 125 golden manifests were rebound after independently checking
every compressed archive hash and source-lock hash. This changes only the
comparison-contract digest: no observation, case membership, tool version,
generation fingerprint or provenance changed. `dssp-contract-rebinding.json`
records all verified hashes. Existing reports retain their historical contract
identities; the final DSSP and smoke reruns use the updated contract.

### Full rerun

The rebuilt release executable completed `--feature bio.secondary-structure.dssp
--dataset all --jobs 1 --output target/benchmark-parity/dssp-bounded-final.json`.
It correctly exited 1 for the remaining disagreements and errors.

| Dataset | Applicable cases | Agree before → after | Disagree after | Errors after |
| --- | ---: | ---: | ---: | ---: |
| PDB | 1,000 | 0 → 556 | 196 | 248 |
| Smoke | 1 | 0 → 1 | 0 | 0 |

The other 150,423 selected IDs remain not applicable. All 557 new agreements
come from the numerical correction; the nullable-label repair moves 13 schema
errors into explicit structural disagreements. The 248 remaining errors include
233 inputs with no analyzable residues in either engine, eight other reference
failures with native results, 8OLZ with failures in both engines, two reference
no-analyzable-residue outcomes with native results (3OXM and 9U3I), and the four
unknown-element interpretation failures listed above.

`dssp-before-serial.json` and `dssp-after.json` are the complete numerical
before/after reports. Their exhaustive 151,424-row audit checks unchanged source
identities and reference observations, with 761 changed native observations and
no loss of prior agreement. `dssp-bounded-final.json` is authoritative after the
nullable-label, finite-intermediate and batching corrections. Its exhaustive
151,424-row verification passed: every native and reference observation is
unchanged from the numerical rerun, exactly the 13 schema errors become
disagreements, and every other complete case row is identical.
`dssp-bounded-verification.json` records this separate audit.

The final cross-feature smoke run, `final-smoke-bounded.json`, completed all 25
features under the updated contract. All 578 case rows are byte-identical to
`final-smoke.json` from before the batching/schema changes. Known residuals remain
visible and the run correctly exits 1. The failed attempt to reuse an existing
report filename was rejected by the runner before execution; the successful run
uses its own unique output path.

### Validation

Passed after the final runtime change: workspace formatting, all-target/all-feature
check, clippy with warnings denied, all-feature workspace tests including doctests,
Rustdoc with warnings denied, Rust 1.89 check and full `kekule` package verification.
After the subsequent benchmark-only changes, formatting, workspace check/clippy,
all 50 benchmark unit tests and integration tests, Rustdoc and Rust 1.89 were
rerun and passed. Commands use `--locked --offline`; the authoritative current
target directory is `target/mmcif-validation` (the MSRV check uses its own default
target). Package commands use `--allow-dirty`.

Also passed: potentials no-default-feature tests/docs, companion and benchmark
package file lists, license-copy checks, eight Biopython tests, 28 pinned RDKit
tests, eight shared reference-runner tests, two atomic-data tests, 21 Python
dashboard tests, Node dashboard checks and `git diff --check`. The first Biopython
test invocation failed while loading a numerical dependency without the runtime
DLL directory; rerunning with the pinned environment's `Library/bin` on Python's
actual `PATH` passed all eight tests. An initial shared-test filename matched
zero tests; the corrected `test_runner.py` invocation ran and passed eight.

Full companion package builds remain deferred as in CI until the foundational
crate is published. Linux fuzz compilation/smoke runs were not run on this
Windows host. Optional external benchmarks are evidence, not release gates;
failed/aborted attempts above are not reported as passing. Root `README.md`
remains unchanged. Detailed logs are `dssp-*.log` and `final-*.log` locally.

All 25 feature reviews and feature-specific full reruns are recorded, and the
bounded independent DSSP reference audit is complete. The consolidated handoff
below covers the full review.

## Consolidated handoff

The review covers all 25 original features, with one feature addressed per turn
and a complete rerun across all five datasets for each feature. A completed
comparison is not necessarily a passing comparison: residual differences,
unsupported inputs and resource failures remain visible. The per-feature sections
above distinguish implementation defects, reference-adapter defects, supported
policy differences and limits that need broader work.

The completion audit reopens the 25 final feature reports, verifies each report's
`complete` state and outcome accounting, and checks that every selected case in
each golden archive was covered. It verifies all 125 current archive hashes.
Four non-Enamine SDF parsing archives were retained after an earlier independent
regeneration proved their observations unchanged; the completion audit again
compares every record in those differently compressed archives (101,348 rows).
All observations and identities match the archives used in the feature reports.
The final 25-feature smoke comparison retains all 578 case rows byte for byte.
Machine-readable evidence is `target/benchmark-parity/completion-audit.json`;
`final-feature-reports.json` identifies the full report for every feature.

These reports were produced at their respective feature-review stages. They are
not represented as one newly executed full cross-product benchmark after the
last change. The last runtime change is confined to DSSP; its full rerun and
exhaustive row audit are recorded above. The subsequent batching change preserves
all native/reference observations in the full DSSP rerun and every smoke row.

### Remaining work and deliberate differences

| Area | Unresolved issue or retained policy |
| --- | --- |
| Input/query coverage | CXSMILES and stereochemical SMARTS need grammar and stereo/group integration. Invalid-valence inputs remain explicit errors. |
| Hydrogen and radical policy | Native declared/inferred H storage differs from RDKit normalization. Bracket-radical inference needs coordinated parser/writer policy; one such case propagates into aromaticity, ranking and substructure. Native metal-hydrogen preservation is retained where RDKit removal loses chemically represented H. |
| Stereo | Broader symmetry cleanup, source wedge interpretation, endocyclic/heteroatom double-bond eligibility and sulfoxide/axis-family handling remain. Local candidate detection is not a complete potential-stereochemistry implementation. |
| CIP resources | 104 supplied cases exhaust the explicit default depth bound. Raising it is not an algorithmic repair; the partial depth-64 probe and its costs are recorded in feature 24. |
| Rotatable bonds | Hydrogen representation dominates residuals; 86 cases require replacing the descriptor's aromaticity approximation through a consistent, fallible perception path. |
| Molecular masses | All audited formulas agree. Modern native atomic constants and charged-mass conventions differ from RDKit's values; replacing them simply to match would discard independently supported data. Missing standard weights remain explicit errors. |
| Ring/mmCIF conventions | One valid alternative cycle basis and 79 permitted multiline-whitespace differences remain. Ring paths and all decoded values are still asserted. |
| DSSP interpretation | Alternate-conformer selection, nonpolymer eligibility, chain-boundary policy, one weak-hydrogen-bond case and unknown element `X` remain as detailed in feature 25. Matching source bytes do not guarantee matching selected coordinates. |
| Resource/coverage limits | Batching no longer accumulates 256 large inputs blindly, but one large structure can still require several GB. Reference deadlines, explicit computation limits and incomplete-run reporting remain necessary. Corpus preselection and the fixed SMARTS query set limit generalization. |

No molecule IDs or corpus-specific branches were added to runtime algorithms.
Corrections use chemical rules, parser/serializer contracts, reference API
semantics or documented numerical precision. Focused regression tests may use
supplied examples; they are not special cases in the implementations.

### Final validation scope

The final runtime state passed workspace all-feature tests and doctests, formatting,
all-target/all-feature check, clippy with warnings denied, Rustdoc with warnings
denied, Rust 1.89, potentials without default features, full foundational-crate
packaging, companion/benchmark package lists and license-copy checks. Benchmark
tests, workspace checks, documentation and MSRV were repeated after the last
benchmark changes. Python reference, atomic-data and dashboard tests and Node
dashboard checks also passed. Feature 25 gives exact log names and the corrected
environment/test invocations; earlier sections preserve their own validation.

Linux fuzz compilation and fuzz smoke execution were not performed on this
Windows host. Full companion package builds remain deferred under the repository's
CI policy until `kekule` is published. These omissions are not presented as passes.
External reference tools remain benchmark-only dependencies. `ARCHITECTURE.md`
ownership boundaries and root `README.md` are unchanged. Changes remain uncommitted
on `codex/benchmark-parity`; no feature work was pushed to `main`.

The bounded independent DSSP regeneration and its exhaustive audit completed,
with no structural changes to the stored expectations. Its total error count
also exposed a final reporting defect: generation incremented `errors` for
reference failures but left `reference_errors` at zero. The generator now
increments both counters. The existing mixed-success/failure generation
regression was extended and failed before this correction. This changes only
generation accounting, not observations, comparison behavior or golden values.
The extended regression now passes. A real three-case DSSP generation using
1XHJ, 3RR8 and 6UFJ produces two successes and one retained reference failure:
`errors = reference_errors = 1`, `kekule_errors = 0`, with all three observations
identical to the completed independent audit. Its expected exit code is 1.
The initial one-case probe selected a successful structure because PDB limits
follow source-lock order; that probe assumption was corrected without changing
the dataset. Evidence is `generation-counter-mixed-probe.json` and its log.

After this one-line accounting correction, formatting, workspace check/clippy,
all 50 benchmark unit tests and integration tests, Rustdoc, Rust 1.89 and the
benchmark package list passed again (`counter-final-*.log`). The final all-feature
smoke run, `final-smoke-counter.json`, completed with all 578 case rows identical
to the preceding smoke runs. Full reference generation was not repeated after
the accounting-only change: the mixed-outcome regression and real reference
probe verify the changed counter, while the complete preceding audit verifies
the unchanged observations. No runtime or comparison algorithm changed here.

The final completion audit passed: 25 complete full-feature reports, all five
datasets per feature, all golden memberships and archive hashes checked,
equivalent reference observations across retained/recompressed archives,
1,020 independently regenerated DSSP reference entries verified, the generation
error counter verified, and all 578 final smoke case rows unchanged. The
requested scoped review is complete as of 2026-09-19; unresolved issues and
validation exclusions are explicitly listed above.
