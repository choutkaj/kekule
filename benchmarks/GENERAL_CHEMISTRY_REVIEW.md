# General chemistry follow-up

This work follows the scoped parity review at commit `13a6f93d`. The goal is
general chemical correctness, with RDKit conventions adopted where they fit
Kekule's represented-chemistry and perception boundaries. Public API changes
are authorized when needed for a clean implementation. No case-specific runtime
rules, discarded benchmark assertions, or silent resource-limit successes are
acceptable.

## Required work

- Hydrogen/radical semantics: retain graph, declared and perceived hydrogen
  distinctions; make consumers use their intended counts; preserve known total
  hydrogen content in SMILES export; implement coordinated bracket-radical
  interpretation and writing without inferring physical spin from electron count.
- Rotatable bonds: remove the local aromaticity approximation, use compatible
  perception or calculate it on temporary state, expose failures, and preserve
  the documented heavy-atom hydrogen convention.
- Stereo: shared chemical eligibility, symmetry-aware cleanup, drawing-aware
  wedge interpretation, and coherent ring/heteroatom/axis handling. Represented
  cleanup remains explicit and transactional.
- CXSMILES: preserve extension syntax and source mapping; interpret supported
  chemical fields; keep metadata in the format sidecar; distinguish explicit
  base-SMILES projection from full chemical interpretation.
- Stereo SMARTS: canonical query stereo constraints checked under atom mappings,
  initially tetrahedral and double-bond constraints, with explicit unsupported
  syntax errors.
- CIP: retain ligand expansion state, expand unresolved comparisons on demand,
  preserve path-dependent auxiliary semantics, and measure every previously
  exhausted case before changing resource defaults.
- Scientific diagnostics: retain raw mass, hydrogen, ring-path and mmCIF value
  comparisons while explaining supported convention differences. Preserve
  current verified CIAAW/AME mass data and mmCIF source whitespace.
- Rings: add permutation and independent cycle-validity regressions; retain
  the remaining selected-ring difference rather than insert a specific cycle.
- Benchmark coverage/resources: preserve running inputs and outputs; record
  execution provenance and explicit resource failures; broaden externally
  sourced query/structure coverage and add focused invariance regressions.
- DSSP interpretation remains deferred as agreed; existing differences stay visible.

## Validation and handoff

Each behavior change requires regression tests and affected full feature reruns.
Complete with one benchmark on a fixed final revision, applicable Rust formatting,
check, clippy, tests, documentation, MSRV and packaging checks, plus changed Python
or dashboard tests. Report platform and publication-related exclusions explicitly.
Reference changes must be independently justified before any golden adoption.

Use `target/general-chemistry-validation` for builds and
`target/general-chemistry-review` for diagnostic outputs. Benchmark history for
this work uses a separate `KEKULE_BENCHMARK_RUNS_DIR`. Existing runs and golden
archives are not modified as a side effect of a comparison.

## Progress

The starting worktree was clean. No live `kekule*` process was found during the
initial process check; outputs remain isolated regardless of that snapshot.
This is a chronological evidence ledger. Later completed audits supersede earlier
in-progress handoffs; the completion audit at the end records the final scope.

### Hydrogen export

Removed the blanket rejection of `Infer { explicit: n }` from SMILES writer
validation. Existing emission already projects the known declared-plus-inferred
total into bracket syntax and rejects missing perception. The source molecule
and its perception remain unchanged. The public writer documentation now states
that total hydrogen content, rather than the storage/inference policy, is exported.

The revised missing-perception regression and the new carbon/nitrogen/oxygen
total-count regression failed before the fix and passed afterward. All 776 core
unit tests, integration tests and 20 doctests passed with:

```text
cargo test -p kekule --all-features --locked --offline --target-dir target/general-chemistry-validation
cargo fmt --all -- --check
git diff --check
```

The isolated release benchmark build passed. Its binary hash is
`d78ec23c27350e3b60eb3360fd85d09d283202e664c9eb7d07cd8e47feb8e288`.
The full `io.smiles.write` comparison against all datasets completed with four
jobs and pinned RDKit 2026.03.3: 129,335 agreements, zero disagreements,
20,913 retained errors and 1,176 not-applicable rows. Output:
`target/general-chemistry-review/hydrogen-smiles-write.json`, with the adjacent
log and case JSONL. An exhaustive parsed-record comparison verified all 151,424
case rows identical to the previous full run (`hydrogen-smiles-write-audit.json`).
The new API regression covers a programmatically constructed declaration that
does not occur in the supplied SMILES inputs. The isomeric full rerun completed
from the retained `hydrogen-writer-bench.exe` with the same verified binary hash:
147,358 agreements, zero disagreements, 2,890 errors and 1,176 not-applicable
rows. Its exhaustive audit found all 151,424 parsed rows unchanged from baseline
(`hydrogen-smiles-isomeric-audit.json`). The canonical rerun completed after the
coordinated bracket changes described below.
Workspace all-target/all-feature checking passed, logged in `hydrogen-check.log`.
Remaining workspace-test/clippy/docs/MSRV/package gates are pending, not passed
or waived.

### Radical-contract evidence for the next change

The pinned reference probe assigns bracket radical-electron counts of four,
three, two, one and zero to `[C]`, `[CH]`, `[CH2]`, `[CH3]` and `[CH4]`.
The current native parser assigns none. A separate direct MOL parser probe
shows that `M RAD` singlet/doublet/triplet codes produce RDKit electron counts
two/one/two. Native `AtomRadical::Singlet.unpaired_electron_count()` returns zero
and that method currently feeds valence and normalization. Spin multiplicity,
radical-electron occupancy and unspecified spin therefore need explicit distinct
semantics before bracket-radical inference is installed. This is not permission
to guess spin multiplicity from a bracket atom's valence deficit.

### Radical representation and occupancy

Replaced the multiplicity-only `AtomRadical` enum with a validated immutable
electron count and optional local spin multiplicity. `new(count, multiplicity)`
rejects zero electron counts and incompatible multiplicities. Unknown spin is
preserved, and the type can represent lower-spin coupling of several electrons.
This is not a molecular-spin or ground-state prediction.

Valence, aromaticity and represented normalization now use `electron_count()`.
The regression for a singlet two-electron carbon center failed before the fix
(four inferred hydrogens instead of two); singlet, triplet and unspecified-spin
two-electron states now each infer two hydrogens. Older charge-only regressions
used the old singlet variant as a zero-electron placeholder. Their inputs now
use no radical, preserving the original tested electron occupancy and assertions.

MOL interpretation preserves the source singlet/doublet/triplet assertions with
two/one/two electrons. V2000/V3000 writers reject unspecified spin and electron
counts not representable by their radical codes, rather than asserting a spin
or changing electron occupancy. Regression coverage includes a four-electron
triplet, which cannot be represented by the two-electron triplet MOL code.

Workspace all-target/all-feature checking, all-feature tests (including 779 core
tests and all integration/doctests), formatting, and clippy with warnings denied
passed (`radical-check.log`, `radical-workspace-tests.log`, `radical-clippy.log`).
The first formatting check caught a line-wrap change from the constructor
migration; formatting was applied and the check passed before clippy.
Workspace Rustdoc with warnings denied passed (`radical-doc.log`). A release
build and the full `algo.valence.rdkit-like` benchmark completed
(`radical-build.log`, `radical-valence.json` and its log). All 301,834 parsed rows
were identical to baseline: 297,952 agreements, zero disagreements, 2,881 errors
and 1,001 not-applicable rows (`radical-valence-audit.json`). This was the radical
representation stage, before bracket inference. The complete scientifically
corrected benchmark radical-observation contract remains pending. The interim native
adapter retains all existing fields and reports unspecified/higher spin explicitly;
its historical `unpaired_electrons` field carries electron occupancy and needs
coordinated reference/schema documentation work. No reference files, contracts
or goldens were changed during the existing writer runs.

Further consumer audit also found that canonical atom-rank signatures currently
omit radical state entirely. The coordinated radical work must test and correct
that distinction before claiming complete downstream support.

### Bracket radicals and canonical ranking

SMILES interpretation now assigns bracket radical-electron occupancy after
aromatic-bond localization, using one octet/duet and allowed-valence rule shared
with writer validation. The model follows RDKit 2026.03.3 `assignRadicals` in
`Code/GraphMol/MolOps.cpp`; its periodic-table outer-electron values belong to
this model, not to a new ground-state prediction. Ordinary organic-subset atoms
retain their hydrogen-inference policy. No spin multiplicity is inferred.

All three SMILES writers require radical atoms to use brackets and reject
brackets whose implied occupancy would differ from represented chemistry.
Explicit spin assertions still require a format that can represent them.
Canonical atom-rank signatures now include electron count and optional spin;
otherwise distinct radical assertions were assigned identical symmetry ranks.
Both bracket inference and radical-ranking regressions failed before their fixes.
Additional tests cover main-group radicals, hypervalence, isolated and bonded
metals, charges, isotope preservation, all three writer round trips, and rejected
lossy output. RDKit probes independently checked the tested electron counts.

The complete workspace suite passed, including 782 core tests, integration tests
and doctests (`bracket-workspace-tests-3.log`). Earlier runs caught integration
fixtures that depended on the old meaning of `[C]`: hydrogen-storage tests now
explicitly construct their intended radical-free atom states, and the stereo
bond-deletion export test explicitly caps the severed fixed-H carbon. Original
hydrogen and stereo assertions remain in place. The additional element/charge
cases passed afterward (`bracket-elements-tests.log`). Workspace all-target,
all-feature clippy with warnings denied passed (`bracket-clippy.log`).

A direct MOL probe confirmed that RDKit does not retain singlet/triplet source
assertions in its atom properties: both become two radical electrons. The old
reference adapter fabricates a triplet label from that count. Correcting this
requires a coordinated observation contract and independent source-metadata
handling, with new reference output audited before adoption. Existing reference
files and golden archives remain unchanged; that work is not complete.

An additional exploratory notation-invariance probe exposed a pre-existing
aromatic-localization limitation: `[n+]1ccccc1` is rejected, while the localized
`[N+]1=CC=CC=C1` can represent the one-electron nitrogen center. RDKit accepts
both and perceives this ring as nonaromatic. The localizer requires a valence
deficit of exactly zero or one before radical inference; its neutral-carbon
exception does not generalize to charged heteroatoms. No new element-specific
exception was added. A general source-aromatic localization review is required;
this charged-heteroatom case is not claimed as supported by the new inference.

Neutral carbon-radical notation invariance is covered for both carbocyclic and
nitrogen-containing rings, including all three writer round trips
(`bracket-notation-tests-2.log`). The earlier exploratory log retains the
charged-nitrogen failure. Formatting and diff whitespace checks passed.
Warnings-denied Rustdoc, Rust 1.89 all-target/all-feature workspace checking and
verified core packaging passed (`bracket-doc.log`, `bracket-msrv.log`,
`bracket-package.log`). Packaging used `--allow-dirty` to validate current work.

The full SMILES parser comparison completed against unchanged goldens with
129,632 agreements, 17,726 disagreements, 2,890 errors and 1,176 not-applicable
rows. Its retained binary SHA-256 is
`7d84b112c90f94b6f17afa26e40f1f04906c7c964da4de7625b3f3267073931a`
(`bracket-radical-bench.exe`). The initial exhaustive audit found 933 changed
records, 1,466 changed radical-electron observations, and one six-membered ring's
aromaticity changes. No case changed comparison status. The canonical writer
full rerun completed from this same binary with reference code and inputs frozen
throughout the run. This is intermediate verification, not completion of
the coordinated hydrogen/radical work or the final fixed-revision benchmark.

The completed parser audit checked all 151,424 rows against baseline and
confirmed that all 1,466 changed electron counts match the unchanged independent
reference. No changed electron count disagrees. The aromaticity changes are
confined to PubChem 181201, the previously documented carbon-radical gap.
Remaining differences include the old reference's inferred spin labels; neither
those assertions nor any golden payload was discarded.

Companion package file-list checks and the no-default-feature potentials tests
and warnings-denied documentation passed (`bracket-package-potentials.log`,
`bracket-package-traj.log`, `bracket-no-default-tests.log`,
`bracket-no-default-doc.log`). Linux/fuzz checks have not been run in this Windows
environment. The other affected full feature reruns and final fixed-revision
verification remain pending the coordinated reference-contract changes and
remaining implementation work.

The full canonical-writer rerun retained 147,358 agreements, zero disagreements,
2,890 errors and 1,176 not-applicable rows (`bracket-io.smiles.canonical.json`).
An exhaustive audit of all 151,424 rows found 71 changed emitted strings but no
changed decoded observations, comparisons or statuses (`bracket-canonical-audit.json`).

### Rotatable bonds: shared perception

Removed the local five-/six-member conjugation approximation. Strict resonance
classification now uses the existing aromatic atom membership, including the
query's distinction between elemental nitrogen and aliphatic oxygen/sulfur.
Installed default valence, Figueras ring selection and RDKit-like aromaticity
are reused when compatible. Otherwise default perception runs on a temporary
copy. The source chemistry and perception stay unchanged. Disabling resonance
exclusions does not invoke chemical perception.

`rotatable_bonds::detect` now returns `Result<RotatableBondSet, PerceptionError>`;
the benchmark propagates failures. Strict valence is revalidated even on the
cached path, because model provenance alone does not distinguish a permissive
valence installation. The existing valence assignment kernel is reused for this
check; matching cached assignments avoid recomputing rings and aromaticity.
Hydrogen-invariant terminal and symmetric-group classification is unchanged.

A focused nonaromatic conjugated-ring regression failed before the fix. It and
the sulfur/selenium aromatic-ring regressions now pass without element-specific
aromaticity exceptions. Additional tests check absent, installed and model-neutral
perception; unchanged source state; permissively cached invalid valence; and the
perception-free general mode. The unsupported-bond-order fixture now uses sulfur
at its focus to represent quadruple-plus-single valence without invalid carbon.
The first new cache test incorrectly assumed an internal H setter erased model
provenance; it now constructs a model-neutral state through the public builder.

The full workspace suite passed, including 786 core tests and all integrations
and doctests (`rotatable-workspace-tests-2.log`). Formatting, whitespace checks,
all-target/all-feature clippy with warnings denied, warnings-denied Rustdoc,
Rust 1.89 all-target/all-feature checking and verified core packaging passed
(`rotatable-clippy.log`, `rotatable-doc.log`, `rotatable-msrv.log`,
`rotatable-package.log`). Packaging used `--allow-dirty`. Linux/fuzz and the
final whole-goal verification remain pending as noted above.

The full five-dataset feature rerun completed from retained binary
`rotatable-perception-bench.exe`, SHA-256
`5eba724cd098885e8b6afbcd23f95635c560af17bbbf408bb7da3361ef2c874b`,
against unchanged reference snapshots. Results: 261,129 agreements, 36,805
disagreements, 2,899 errors and 1,001 not-applicable rows.
The exhaustive 301,834-row audit verified every expected observation unchanged.
All 86 previously classified aromaticity-approximation cases now have the exact
bond sets from the independently computed hydrogen-collapsed reference. Eighty
became full agreements; six retain raw hydrogen-representation differences.
Eighteen already-failing cases now also report the native strict-valence failure.
No other case changed (`rotatable-perception-audit.json`). The collapsed-reference
comparison is diagnostic evidence only; raw benchmark assertions remain intact.

### Benchmark radical contract

Observation contract 3 replaces the misleading `unpaired_electrons` field and
guessed `radical` labels with required `radical_electrons` and required nullable
`spin_multiplicity`. Both measurements remain compared exactly; zero spin and
missing measurements are rejected by the schema. Unknown spin is not a triplet
or doublet assertion. Report/archive framing remains at schema 2.

`reference/rdkit/source_radicals.py` independently retains the explicit spin
metadata RDKit loses. It reads V2000 atom-block/property precedence, V3000 atom
order and continued/quoted records, and CX radical codes while respecting label
and coordinate fields. Its rules follow BIOVIA CTfile Formats 2020 and Chemaxon
CXSMILES documentation, linked in the source and guide. RDKit still owns electron
counts and chemistry. If RDKit changes an explicitly supplied count, observation
fails instead of fabricating agreement. A direct probe found that RDKit ignores
the legacy V2000 atom-block doublet code; that limitation remains a reference
failure. Ordinary metadata cannot inject the reader's internal assertions.

Fifty Python reference/runner tests pass (`radical-contract-reference-tests-5.log`),
including singlet/triplet distinction, unknown bracket spin, all seven CX codes,
literal metadata, fragment/H transformations, and internal-property collision
protection. Native schema tests assert occupancy and spin independently. All 53
benchmark unit tests passed (`radical-contract-unit-tests-2.log`). Integration
tests initially stopped on the intentionally stale smoke contracts; fresh
independent smoke generation and audits are being adopted before rerunning them.

The runner fingerprint now includes the new metadata reader, with a regression
that changing this source changes the fingerprint. A real Windows test failure
also exposed timestamp collisions in parallel temporary fixtures; an atomic
sequence and concurrent-fixture regression resolve it. These are benchmark
infrastructure changes, not chemical comparison exceptions.

All 23 molecular smoke features were independently regenerated under contract 3.
Their audits found no changed observations beyond the measured-field migration;
the first parser audit covers 61 renamed atom observations and the SDF audit 447.
The two biological smoke archives have no radical fields and unchanged comparison
rules: their payloads remain byte-for-byte identical, with original reference
provenance retained and the contract migration recorded in their manifests.
Evidence: `radicals-3-all-smoke-adoption-audit.json` and preceding smoke audits.

The active full-corpus candidate is
`target/general-chemistry-review/reference-candidates/radicals-3`. Earlier isolated
candidates were stopped when missing fingerprint coverage and metadata lexical
issues were found; their full-corpus output is not adopted. The active generator
uses retained binary `radical-contract-bench.exe`, SHA-256
`60bee1982c10179658cf837f5e3a8386f3cbf2cf1cad4e89315b67efb133ab3c`.
PubChem SMILES generation completed with 100,000 outcomes and the same nine
reference failures. Generation reports failure status for retained case errors,
so subsequent datasets continue only after checking the report is complete.
The full reference audit and native rerun completed. Independent regeneration
preserved all electron counts and all unrelated observations across 100,000
PubChem rows (2,280,748 atom observations), 50,240 Enamine rows (1,296,159 atom
observations), and the 1,164 inapplicable biological rows. The old PubChem
reference had guessed 1,466 spin labels from electron counts; none of these
inputs explicitly asserted spin. The new reference records unknown spin while
continuing to compare electron counts exactly. Evidence:
`radical-reference-audit-{dataset}.json` and `radical-reference-audit-all.log`.

The full native SMILES rerun produced 130,544 agreements, 16,814 disagreements,
2,890 errors and 1,176 not-applicable rows (`radical-contract-smiles.json`). An
exhaustive audit of all 151,424 rows found exactly 912 disagreement-to-agreement
transitions, no other status transitions, and no native observation changes
beyond the field migration (`radical-contract-cases-audit.json`). These resolved
comparisons reflect the corrected spin contract, not additional native chemistry
changes during this stage. Remaining failures stay visible.

All benchmark unit and integration tests pass after smoke adoption
(`radical-contract-bench-tests-7.log`). Warnings-denied workspace documentation
and clippy passed; the final runner diagnostic change also passed benchmark
clippy and Rust 1.89 all-target checking (`radical-contract-clippy-2.log`,
`radical-contract-msrv-2.log`). The runner now reports an entirely inapplicable
selection accurately instead of describing it as reference-generation errors;
an integration regression verifies that such a selection still does not pass.
Formatting and whitespace checks passed. This stage changed no runtime crate
code beyond the previously validated stages. Full workspace tests and package
checks were not repeated here; final whole-goal verification and Linux/fuzz
checks remain pending.

Other full feature archives still require audited contract migration or
regeneration. Full non-smoke candidates have not replaced the original archives;
the completed SMILES comparison explicitly used the audited candidate directory.
This stage does not establish final benchmark or goal completion.

### CTfile charge/radical property precedence

The BIOVIA CTfile Formats 2020 specification (V2000 properties, printed page 49;
V3000 atom fields, printed page 8) confirms that either V2000 `M CHG` or `M RAD`
replaces all legacy atom-block charge/radical values. Both formats support radical
code zero. Native parsing previously retained unlisted legacy values and rejected
zero. The parser now clears the legacy values once before applying the first
charge/radical property, so repeated lines and either property ordering accumulate
correctly. Zero maps to the same nonradical state as the format default.
Source: https://discover.3ds.com/sites/default/files/2020-08/biovia_ctfileformats_2020.pdf

Focused regressions reproduced both failures before the fixes and now pass. They
cover unlisted legacy charges/radicals, property order, multiple property lines,
explicit zero and V3000 round trips. Full workspace tests (788 core unit tests,
integrations and doctests), warnings-denied clippy and Rustdoc, Rust 1.89
all-target/all-feature checking, and verified core packaging passed. Logs:
`ctab-workspace-tests.log`, `ctab-clippy.log`, `ctab-doc.log`, `ctab-msrv.log`,
`ctab-package.log`. Formatting and whitespace checks passed. Packaging validates
the current worktree with `--allow-dirty`; no publication was performed.

Retained native binary: `ctab-properties-bench.exe`, SHA-256
`d07cbc06a7a1c4d8d1cf0a17e8cf51920f195dd7b9f60ac87928bf73149879d5`.
Independent full MOL and SDF reference generation is running in
`reference-candidates/ctab-properties` with unchanged reference source. Both
native revisions will be compared against the same candidates. The original
goldens remain untouched. A suspected stall was disproved: Windows directory
listings showed stale sizes for open compressed files, while direct file reads
confirmed progress. Single-record and batch process probes also passed; the
original runs were not interrupted or restarted.

Both reference generations completed across all five datasets, with nine
retained reference errors per feature. Exhaustive audits found no changed
electron counts or unrelated observations. Each PubChem feature contains
4,338,472 atom observations: 378 explicit doublet assertions are retained, and
ten formerly guessed triplet labels become unknown spin. Enamine contributes
1,296,159 atom observations per feature, and PL-REX 16,932. Evidence:
`ctab-reference-audit-pubchem.log`, `ctab-reference-audit-other.log` and their
per-dataset JSON audits. Original archives have not been replaced.

The first native comparisons exhausted disk space while writing raw case JSONL.
All four processes terminated and their reports are incomplete; they provide no
full-corpus validation. Their outputs remain in place. Rebuildable incremental
compilation state under this task's build directory (6,329,877,969 bytes) was
removed, and transparent NTFS compression was enabled for the task output
directory. Compression of the four failed case files is checked against SHA-256
hashes before fresh comparisons. This local storage measure leaves logical file
contents and paths intact. A portable compressed-case output option remains an
identified benchmark resource improvement for the resource/coverage stage.

Compression finished successfully: 4,803,524,485 logical bytes now occupy
1,342,775,296 bytes. All four SHA-256 checksums are unchanged
(`ctab-storage-compression-audit.json`). The completed debug-build cache was
also removed (9,816,840,116 logical bytes); test logs and retained executables
remain. Fresh comparison workflows are running against all five datasets,
first the previous native binary and then the fixed binary for each feature.
Their reports use `ctab-{mol,sdf}-verified-{before,after}.json`, avoiding the
failed report paths. Newly created case files have the inherited compressed
attribute. `audit_ctab_cases.py` is prepared to verify all expected observations,
native changes and status transitions once both workflows finish.

The complete baseline runs each contain 126,887 agreements, 23,689 disagreements,
nine errors and 1,003 not-applicable rows. Their plain case files are each
6,155,244,426 bytes. The subsequent fixed-parser runs again exhausted space and
remain incomplete (`ctab-{mol,sdf}-verified-after.json`); they are not validation
evidence. Three completed case files from earlier stages were transparently
compressed in place after their reports were verified complete. Their paths and
logical contents remain unchanged. This freed additional space while a portable
output fix was implemented below.

### Streaming compressed benchmark case records

New comparison runs now write `.cases.jsonl.gz` using streaming fast gzip. Every
case, raw observation and difference remains in the JSONL payload; report schema
and chemical comparisons are unchanged. The `cases` path in each report identifies
its artifact, so existing plain JSONL reports are unaffected. Compression is
explicitly finalized and flushed before marking a run complete. Setup/evaluation
failures also finalize the partial stream where possible, and output failures
keep the report incomplete. No reference source or golden changed in this stage.

The new CLI regression failed on the old plain-file behavior and now verifies
the compressed path, checksum/trailer, all case counts and retained actual and
expected observations. Existing setup-error regressions now read compressed
streams to EOF and require zero decoded records. A unit regression checks errors
from both finalization writes and flushing. All 54 benchmark unit tests and seven
integration tests pass (`compressed-cases-tests.log`). All 21 dashboard tests pass.
Benchmark all-target/all-feature check, warnings-denied clippy and Rustdoc, Rust
1.89 checking, package file listing, formatting and whitespace checks passed
(`compressed-cases-{check,clippy,doc,msrv,package-list}.log`). Runtime workspace
tests and verified core packaging were not repeated because this stage changes
only the benchmark layer; their CTfile-stage results above remain applicable.
The benchmark crate is not publishable. Linux/fuzz and final whole-goal checks
remain pending.

Retained executable `compressed-cases-bench.exe` has SHA-256
`ea6a4dfe01f412e37cc3f3f7428fa0ab7803ddac0ac5d5be2ecec3103d98ccc7`.
Both full fixed-parser comparisons completed from it against the audited CTfile
reference candidates, with reports `ctab-{mol,sdf}-compressed-after.json`.
Each retains 126,887 agreements, 23,689 disagreements, nine errors and 1,003
not-applicable rows. The exhaustive audit checked all 151,588 rows per feature
(303,176 total), finding identical expected observations, native observations
and comparison statuses (`ctab-cases-audit-summary.json`). The corrected source
rules do not change these supplied corpus records; their focused regressions
cover the formerly incorrect behavior. Compression preserves complete case
content while reducing the two 6,155,244,426-byte files to 668,628,175 and
668,305,121 bytes. No incomplete report has been treated as a passing run.

### Atomic report publication

Disk exhaustion exposed a second output defect: truncating reports before
rewriting them destroyed the last valid progress snapshot, and archive copies
could leave empty JSON files. `runner/report.rs` now writes and syncs a temporary
file in the destination directory before publication. First publication uses a
non-overwriting hard link; updates atomically rename the completed snapshot over
the report owned by this run. Failed writes or publication retain the previous
snapshot and clean up the temporary file. Dashboard summary archives use the
same publication boundary. Existing reports from other runs remain protected.

Regressions cover complete successive snapshots, refusal to replace another
run's file, injected partial-write failures before and after first publication,
temporary-file cleanup, and Windows replacement failure under an open reader
that denies deletion. The previous report remains readable and a later retry
succeeds once that reader closes. All 57 benchmark unit tests and seven
integration tests pass (`atomic-report-tests.log`). Benchmark all-target checking,
warnings-denied clippy and documentation, Rust 1.89 checking, package file listing,
formatting and whitespace checks passed (`atomic-report-{check,clippy,doc,msrv,
package-list}.log`). Runtime source did not change; workspace tests and core
package verification were not repeated. The prior 21 dashboard tests remain
applicable; the Rust history integrations exercise the changed archive writer.

Retained executable `atomic-report-bench.exe` has SHA-256
`fc96687f363709da205d9d5b2192ac4e236caa17cd4b197aa9d55cc563ec6d9d`.
Full smoke reruns across all 25 features completed under both the previous and
atomic report writers. All 578 case records are identical except timing fields
(`atomic-report-smoke-audit.json`): 375 agreements, 73 disagreements, four errors
and 126 not-applicable rows. Existing chemical differences remain visible.
Reports: `atomic-report-smoke-baseline.json`, `atomic-report-all-smoke.json`.
The empty archives from the earlier disk failures remain historical failed
artifacts, which the dashboard explicitly skips; no historical result was
fabricated or overwritten. Remaining chemistry and final whole-goal verification
are still pending.

### CXSMILES record preservation and explicit base projection

The SMILES document now separates base syntax, an optional opaque CX extension,
and a record name while preserving the complete source and original byte spans.
Source atom mappings expose record-global indices, independently of connected
component partitioning and component-local atom IDs. Names and raw extensions
remain format metadata on the interpretation; they do not enter the molecule.

`SmilesDocument::interpret_base` and its namespace helper explicitly omit extension
semantics and report the omitted source span. Default interpretation rejects every
nonempty extension until chemical field support is implemented. Empty extensions
carry no omitted chemistry. The benchmark calls the same full library pipeline;
it neither strips extensions nor uses the lossy projection. Record limits include
all metadata bytes, and grammar errors remain located in the base source.

This is input-preservation groundwork, not complete CX chemical support. The
2,881 external Enamine CX records contain enhanced stereo groups and the `r`
flag. The Chemaxon specification assigns chemical meaning to these fields; `r`
must not be discarded solely because the pinned RDKit reader ignores it alone.
Coordinates can also override base stereo, and labels can encode query atoms.
These fields require interpretation or an explicit unsupported error, rather
than blanket treatment as harmless metadata. Source:
https://docs.chemaxon.com/latest/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.html

Five focused core regressions cover metadata, original byte offsets, interleaved
components, ring closures across dots, global atom numbering, malformed framing,
input limits, and strict versus explicitly lossy interpretation. A benchmark
regression covers names and record order through the shared parser. All workspace
tests pass, including 793 core and 58 benchmark unit tests, integrations and
doctests (`cx-record-workspace-tests-3.log`). The first run caught an incorrect
new test expectation: benchmark observations number parsed records, not physical
lines. The corrected test retains that contract. A subsequent integration run
hit Windows access denied while publishing a report; its isolated retry and the
full workspace rerun passed. The original failure remains in
`cx-record-workspace-tests-2.log`; its exact external cause is not established.

The retained release executable `cx-record-bench.exe` has SHA-256
`e4f3f66d9330aabdbb32f0289749db105543f72e0076491ceb9e9c8c8da38560`.
The full `io.smiles.parse` rerun uses the unchanged audited `radicals-3`
reference candidates. It retains 130,544 agreements, 16,814 disagreements,
2,890 errors and 1,176 not-applicable rows (`cx-record-smiles.json`). Unsupported
CX errors now originate in the library and include the extension's source offset.
The exhaustive audit verifies all 151,424 rows: exactly 2,881 unsupported-CX
diagnostic messages changed, with no other non-timing field changes
(`cx-record-cases-audit.json`). All 25 smoke features were also rerun because
they share input adapters. All 578 rows are unchanged, including their existing
73 disagreements and four errors (`cx-record-smoke-audit.json`). The reference
archives and observations were not modified.

Workspace all-target/all-feature checking, warnings-denied clippy and Rustdoc,
Rust 1.89 checking, formatting, whitespace checks, verified core packaging and
companion/benchmark package file lists passed. Logs use `cx-record-{check,clippy,
doc,msrv-2,package,package-list,potentials-package-list,traj-package-list}.log`.
No-default-features potentials tests and warnings-denied documentation also passed
(`cx-record-no-default-{tests,doc}.log`).
Python reference/data-generator and dashboard tests were not repeated: their
sources did not change in this stage, and the CLI integrations exercise dashboard
publication. Linux CI and `cargo check --manifest-path fuzz/Cargo.toml --bins
--locked` / `cargo +nightly fuzz run ...` remain unrun in this Windows session;
they require the supported CI/fuzz environment before final whole-goal handoff.
Full CX chemistry, the other required feature work, and the fixed-final-revision
benchmark remain pending.

### CXSMILES radical and enhanced stereo interpretation

Added a bounded, source-indexed CX field reader for radical codes `^1` through
`^7`, tetrahedral `a`, `oN` and `&N` groups, and the record-level `r` flag. It
checks integer overflow, index ranges, field delimiters, duplicate assertions,
group membership and ownership. Repeated clauses for the same group merge by
type and source group number; numbers are format labels, not canonical IDs.
Absolute memberships partition across components. Relative relationships across
components return an error because the molecular owner cannot represent them.
No successful interpretation silently skips an unsupported field.

Explicit radicals override bracket inference after aromatic localization and
before represented-stereo publication. Spin is retained only when a CX code
asserts it. Both bracketed and unbracketed atoms use the existing hydrogen policy;
no perception is installed implicitly. Coordinates, coordinate bonds, labels,
query annotations and other unsupported fields remain available in the document
but require the explicit base projection to omit their semantics.

The `r` flag becomes a `Relative` group on otherwise ungrouped tetrahedral
assertions; explicit enhanced groups take precedence. This records a common
relative configuration without claiming a pure sample or a racemic composition.
It does not invert alkene geometry. Chemaxon's current query guide explicitly
distinguishes the legacy chiral flag from enhanced groups and describes an
unflagged structure as either enantiomer or their mixture:
https://docs.chemaxon.com/latest/jchem-base_stereochemistry.html
The CX field definitions are documented at:
https://docs.chemaxon.com/latest/formats_chemaxon-extended-smiles-and-smarts-cxsmiles-and-cxsmarts.html
RDKit's omission of standalone `r` is not copied.

Five core regressions cover every radical code, source-to-component mapping,
aromatic radicals, source assertion preservation, mixed and repeated enhanced
groups, relative-flag precedence, and malformed/unsupported fields. The benchmark
regression verifies radical/spin and stereo-group observations and proves that
omitting a group still fails comparison. All workspace tests pass: 798 core and
59 benchmark unit tests, integrations and doctests (`cx-chemistry-workspace-tests.log`).
Workspace check, warnings-denied clippy/documentation and Rust 1.89 checks passed
(`cx-chemistry-{check,clippy,doc,msrv}.log`).

A new test initially assumed strict perception rejected `[CH4] |^1:0|`. Inspection
showed that the existing documented RDKit-like model deliberately bypasses its
radical occupancy check for fixed-H atoms. This stage does not change that model
contract. Import and perception retain both explicit assertions; neither silently
removes the radical as RDKit's sanitization does. General electronic-state
consistency validation remains a limitation to report, not a CX-specific patch.

Retained binary `cx-chemistry-bench.exe` has SHA-256
`4cc82282e1e59c935a57c55911967546cde7d6075b0ccb630c3fde554508432a`.
The full SMILES comparison against unchanged audited `radicals-3` references
completed (`cx-chemistry-smiles.json`): 132,796 agreements, 17,443 disagreements,
nine errors and 1,176 not-applicable rows. All 2,881 formerly unsupported CX
records now interpret: 2,252 agree and 629 disagree. Exhaustive comparison of
all 151,424 rows proves no changed reference observation or non-CX non-timing
field (`cx-chemistry-cases-audit.json`).

Every explicit enhanced group agrees with the reference. Of the CX records,
2,711 have identical complete group observations; the remaining 170 add only the
relative relationship of otherwise ungrouped tetrahedral centers. An independent
RDKit 2026.03.3 source/index audit checks every one of those 170, including salts:
all explicit groups are unchanged and every relative member is precisely an
ungrouped tetrahedral source center (`cx-relative-groups-audit.json`). This is a
documented interpretation difference: Kekule honors the supplied relative flag,
whereas RDKit assumes absolute configuration for those centers. The remaining
raw differences are hydrogen declaration/inference and their explicit-valence
observations. They are retained, not normalized away.

All 25 smoke features were rerun from the retained executable. All 578 case rows
are unchanged except timing (`cx-chemistry-smoke-audit.json`), retaining 375
agreements, 73 disagreements, four errors and 126 not-applicable rows. Verified
core packaging, companion/benchmark package lists, no-default-features potentials
tests/documentation, formatting and whitespace checks also passed; logs use the
`cx-chemistry-` prefix. Python reference/data-generator and standalone dashboard
tests were not repeated because their code did not change. Linux/fuzz checks
remain unrun for the Windows-environment reason recorded in the previous stage.
Unsupported CX field families remain explicit limitations; stereo SMARTS, shared
stereo rules, CIP expansion, scientific diagnostics, broader external coverage
and fixed-final-revision verification remain required goal work.

### Stereochemical SMARTS and mapping constraints

`QueryGraph` now owns checked, syntax-independent local tetrahedral and
double-bond constraints. The builder validates live focuses, distinct adjacent
carriers and duplicate assertions, canonicalizes tetrahedral carrier order with
its parity, and revalidates constraints when publishing the graph. These are
query predicates, not represented molecular chemistry or installed perception.
The matcher checks them under a complete atom mapping before deduplication and
the match cap. It neither assigns CIP nor perceives stereo implicitly. An
achiral query remains unconstrained; a stereo query requires specified target
stereo. Two or more omitted tetrahedral carriers can complete either handedness;
one omitted carrier has a unique remaining position.

The SMARTS frontend retains source neighbor order, including ring placeholders,
and supports `@`, `@@`, `@TH1`, `@TH2` and paired directional bonds. Inline H
at a root has the appropriate source-order parity when one carrier is omitted.
Hydrogen-count predicates remain separate from graph neighbors, including when
a query matches an explicit graph hydrogen. Negated or disjunctive stereo,
unspecified stereo alternatives and other stereo classes still fail explicitly.
Contradictory directions are rejected rather than copied from RDKit's permissive
first-direction behavior. Enhanced stereo group matching is outside this local
constraint API; target group relationships are not interpreted.

Reference semantics were checked against the Daylight SMARTS theory and pinned
RDKit 2026.03.3 probes:
https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
https://www.rdkit.org/docs/RDKit_Book.html#smarts-support-and-extensions
Directional bonds use RDKit's single-or-aromatic predicate. A lone direction
does not define a cis/trans relationship. This is deliberately distinct from
an explicit `-` query, which excludes aromatic bonds.

External testing also found that the unbracketed atom reader could consume `Cn`
as copernicium instead of carbon followed by aromatic nitrogen. Its lexical
rule now admits only `Cl` and `Br` as two-letter unbracketed elements. Bracketed
elements still use the complete periodic table. Regressions cover both forms,
carrier permutations, explicit/omitted H, partial environments, tetrahedral and
directional ring closures, reversed endpoints, conjugated double bonds,
uniqueness/limits, builder validation and deterministic malformed mutations.

The RDKit substructure adapter now explicitly sets `useChirality=True`, with a
reference-side positive/opposite/unspecified regression for both supported stereo
types. Its existing 18 queries are achiral; their results do not establish stereo
coverage. Likewise `query.smarts` checks acceptance and graph counts, not predicate
equivalence. The guide states these limits rather than advertising complete
SMARTS coverage from the current molecular inputs.

Retained executable `stereo-smarts-verified-bench.exe` has SHA-256
`b61be8cba0467dbc9fa98be14ccbf516d7ea7c9c3abfc5f12a8bc96222112e51`.
The full parser feature rerun (`stereo-smarts-verified.json`) agrees on all
150,248 applicable records, with zero differences/errors and 1,176 not-applicable
rows. Before this stage, 20,904 records failed: 20,826 at atom stereo and 78 at
directional bonds. Exhaustive paired-case auditing proves every old successful
observation and every reference observation is unchanged
(`stereo-smarts-cases-audit.json`).

The substructure rerun selected up to 1,000 source IDs per dataset and covered
both supplied molecular formats: all 4,353 applicable records agree; 1,001 remain
not applicable. All 5,354 case rows are identical except timing to the retained
pre-stage binary (`stereo-smarts-matching-audit.json`). Full query and selected
matching goldens were independently regenerated into separate candidate dirs
under the current contract, then compared to the original stored rows. All
151,424 parser and 5,354 selected matching reference rows are unchanged, including
their source provenance (`stereo-query-reference-audit.json`). No original full
golden was replaced and no comparison was weakened.

A separate differential matching audit derives reordered, unspecified and
single-center/bond-inverted variants from externally supplied connected molecules
of at most 60 atoms. It selects tetrahedral and double-bond strata independently,
up to 100 of each per dataset, retaining source paths, file hashes, record indices,
reference version and input/executable hashes. It found 281 source molecules
(165 PubChem, 116 Enamine; 214 tetrahedral and 81 double-bond-bearing, overlapping).
Across 2,276 complete mapping comparisons, 2,274 agree and two remain errors, with
no mapping disagreement (`stereo-query-audit.json`). Both errors concern the
same independently generated variant of PubChem CID 446377: source-aromatic
localization rejects an aromatic atom whose incident written bonds are explicit
directional/single or double bonds. This input-interpretation limitation is
retained for the shared stereo/aromatic review, not bypassed in query matching.
The audit is supplementary local evidence, not a new checked-in toy corpus.

Validation for this stage: all workspace tests and doctests pass, including 804
core and 59 benchmark unit tests (`smarts-verified-workspace-tests.log`). The final
focused query run passes 23 selected tests. The 43 RDKit reference tests pass
(`smarts-python.log`). `cargo fmt --all -- --check`, workspace/all-target/all-feature
`cargo check` and warnings-denied `cargo clippy`, warnings-denied workspace
documentation, Rust 1.89 workspace checks, verified core packaging, companion and
benchmark package-file lists, no-default-features potentials tests/documentation,
and `git diff --check` pass. Logs use the `smarts-` prefix under the review output
directory; builds disable incremental compilation and use the validation/MSRV
directories noted above.

One intermediate workspace run failed in `writer_setup`: its report retained the
initial incomplete snapshot with `error: null`, so the test could not read the
expected writer-setup error. The isolated retry and final complete workspace run
passed. The exact transient cause was not established; the failed test log is
retained as `smarts-final-workspace-tests.log`. No retry or comparison rule was
weakened to make it pass. An early parser-comparison dashboard refresh also used
an incorrectly typed Python path; the final reruns use the correct interpreter
and refreshed history successfully. The four previously documented empty history
archives remain untouched and skipped by the dashboard.

Linux CI, `cargo check --manifest-path fuzz/Cargo.toml --bins --locked`, and
`cargo +nightly fuzz run ...` remain unrun on this Windows host and require the
supported Linux/fuzz environment. Unchanged atomic-data generator and standalone
dashboard tests were not repeated. Full companion package verification remains
deferred until the foundational crate is published, following the CI policy;
their required file-list checks passed. Full-corpus substructure comparison and
the fixed-final-revision whole-project benchmark remain for the final validation
stage; this stage's full corpus run is `query.smarts` and its matching evidence
has the explicitly recorded sample/stratification limits above.

The remaining implementation work is shared stereo eligibility/cleanup and
drawing interpretation, CIP expansion, scientific mass/hydrogen/ring/mmCIF
diagnostics, independent ring/permutation tests, and broader persistent external
query/structure coverage. Unsupported CX fields, the source-aromatic variant
limitation and general electronic-state consistency validation remain explicit
limitations. DSSP interpretation stays deferred. The overall goal is not complete.

## Shared local stereo eligibility

The next stage repairs local candidate geometry without changing the represented
versus perceived stereo contract. Candidate discovery and coordinate inference
now share tetrahedral carrier construction and a conservative terminal-ligand
equivalence proof. Three-direction phosphorus and arsenic environments include
the lone-pair carrier, including phosphines/arsines with one hydrogen ligand.
Four ligand directions can include a double bond, covering phosphoryl centers.
Existing sulfur/selenium valence and charge restrictions remain in place.

Terminal equivalence compares element, isotope, charge, radical occupancy/spin,
attachment bond and all terminal hydrogen identities. Atom maps do not distinguish
chemical ligands. Implicit and explicit hydrogen representations give the same
answer when hydrogen assignments are complete. Unknown hydrogen counts do not
constitute a proof. A terminal group with all-distinct hydrogen isotopes can
itself be stereogenic and is deliberately left to the general stereo dependency
analysis; the regression with opposite CHDT configurations preserves the central
pseudoasymmetric candidate. This is a local equivalence proof, not a replacement
for general branch or ring symmetry analysis.

Source directional-bond decoding, MOL drawing interpretation and candidate
perception share the represented double-bond geometry bounds: at most three
graph neighbors per endpoint and no containing ring shorter than eight atoms.
The shortest-ring
test is independent of a selected cycle basis. Heteroatom endpoints no longer
exclude larger-ring imines/azo groups, and an atom's aromatic flag does not by
itself exclude a nonaromatic double bond. Perception still excludes aromatic
bonds. Source interpretation does not install perception, default perception
does not clean/materialize represented stereo, and CIP continues to rank
represented assertions rather than silently applying candidate filters.

The local rules were checked against pinned RDKit 2026.03.3
[`FindStereo.cpp`](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/FindStereo.cpp),
[`QueryOps.cpp`](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/QueryOps.cpp),
and [`ConjugHybrid.cpp`](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/ConjugHybrid.cpp),
with independent Python probes. RDKit's iterative stereochemical ranking and
dependent-ring treatment are not equivalent to ordinary canonical color ranks.
Nitrogen inversion eligibility also depends on conjugation and bridgehead/ring
geometry; neither a count of ring bonds nor a short list of molecule-specific
exceptions supplies that model. Those general algorithms remain the next stereo
work, together with explicit transactional cleanup and wedge interpretation.

Four new regressions cover pnictogen/multiple-bond tetrahedral carriers,
terminal-equivalence isotope/map/hydrogen invariance, ring-size/heteroatom
eligibility, and large-ring imine stereo through both V2000 and V3000 model
round trips. The final workspace run passes all tests and doctests, including
808 core and 59 benchmark unit tests (`stereo-eligibility-workspace-tests.log`).
Formatting, all-target/all-feature check and warnings-denied clippy, warnings-denied
workspace documentation, Rust 1.89 checks, verified core packaging, companion and
benchmark package-file lists, package license checks, no-default-features
potentials tests/documentation, and `git diff --check` pass. Validation uses the
same offline/nonincremental build directories as the previous stage. Two initial
clippy attempts identified obsolete imports after centralizing the geometry rule;
both imports were removed before the final clean run.

Linux/fuzz commands remain unrun on this Windows host. Unchanged atomic-data,
standalone dashboard and Python reference tests were not repeated in this stage.
Companion package verification remains subject to the foundational-crate
publication policy documented above. No README or reference adapter changed.

The first full comparison exposed a real intermediate regression: requiring two
graph neighbors excluded an imine endpoint carrying an implicit hydrogen. The
lower degree bound was removed, and implicit/explicit hydrogen regressions for
`CC=N` and `CC=[NH]` were added before repeating the checks and full comparison.
The retained final executable is `stereo-eligibility-verified-bench.exe`, SHA-256
`64cd9058653371b7bf81664e9de7be4e1b9e9562de01bd11a6897dc0965dee3a`.
The earlier executable and its reports remain diagnostic evidence, not the final
result. Formatting, checks, tests, clippy, docs, MSRV and verified core packaging
were repeated after that correction.

Independent RDKit generation into `reference-candidates/stereo-eligibility`
completed for all datasets under the current contract. All reference outcomes,
including 18 reference errors, were retained. The pre-stage retained SMARTS
executable and final executable were compared against exactly those same
observations; no original full golden was replaced. The full feature was run
per dataset, with these results:

| Dataset | Before agreement | Final agreement | Final disagreement | Errors | Not applicable |
| --- | ---: | ---: | ---: | ---: | ---: |
| PubChem | 131,558 | 162,188 | 37,794 | 18 | 0 |
| Enamine | 35,060 | 49,174 | 51,306 | 0 | 0 |
| PL-REX | 198 | 272 | 56 | 0 | 0 |
| PDB | 0 | 0 | 0 | 0 | 1,000 |
| Smoke | 14 | 14 | 11 | 0 | 1 |
| Total | 166,830 | 211,648 | 89,167 | 18 | 1,001 |

Across all 301,834 paired rows, 46,931 formerly disagreeing rows now agree and
2,113 formerly agreeing rows now disagree: a net improvement of 44,818. The
exhaustive audit proves identical references, source identities and every native
observation outside candidates and stereo. It independently checks terminal
equivalence for all 83,722 removed candidates, and element/coordination or shortest
ring bounds for all 8,234 additions. Surviving stereo assertions are unchanged.
There are 92 removed inferred assertions and 131 newly interpreted assertions;
68 of the latter have no eligible perceived candidate and expose the outstanding
source cleanup issue. Reports use `stereo-eligibility-before-*` and
`stereo-eligibility-verified-*`; paired audits and the aggregate
`stereo-eligibility-verified-audit-summary.json` are retained beside them.

The new disagreements remain explicitly recorded, rather than obscured by the
net gain. PubChem has 2,103: 2,082 involve newly supported local candidates that
RDKit excludes, with general symmetry analysis still missing; 18 involve newly
interpreted large-ring heterocycle drawing assertions absent from the perceived
candidate set, including assertions on aromatic bonds; three
remove candidates with provably identical terminal ligands. The last three are
RDKit representation inconsistencies: guanidine, an isopropylidene carbanion
fragment and a boron/carbon double-bond fragment lose those RDKit candidates on
hydrogen expansion, while the repeated substituents are unchanged. These probes
are recorded in `stereo-eligibility-reference-terminal-symmetry.json`; native
hydrogen-invariant exclusion is retained as the scientifically correct result.
Enamine has ten new disagreements, representing five molecules in two supplied
formats. Exact chirality-respecting graph self-matches prove odd ligand
permutations fixing their new phosphorus/sulfur centers
(`stereo-eligibility-enamine-symmetry.json`): these need general ring symmetry,
not element-specific exclusions. Thus 2,110 newly disagreeing rows remain native
symmetry/cleanup work; the three justified reference differences stay visible.

The final all-feature smoke rerun covers 578 rows with unchanged outcomes and
only the two justified candidate removals; its paired audit proves no other
observational change (`stereo-eligibility-verified-smoke-retry-audit.json`). The
first final smoke attempt stopped after the first feature with Windows access
denied (`os error 5`), leaving the initial incomplete report with `error: null`.
Its log/report are retained under `stereo-eligibility-verified-smoke-all`; retrying
to a new output completed. This recurring publication problem is not considered
fixed and requires a separate benchmark-layer investigation. No comparison,
golden, or test assertion was weakened to overcome it. The four older empty
history archives remain untouched and skipped by the dashboard.

Local eligibility is only this stage of the agreed work. General branch/ring
stereo dependencies, nitrogen inversion/conjugation, explicit transactional
cleanup of represented stereo, drawing interpretation, CIP expansion, scientific
mass/hydrogen/ring/mmCIF diagnostics, independent ring/permutation checks and
broader persistent query coverage remain. The overall goal remains active.

### General stereo symmetry and explicit cleanup

Stereo candidate detection now uses bounded exact graph automorphisms after
local eligibility. Color refinement only rejects impossible mappings; matching
colors alone never establish equivalent ligands. The detached analysis graph
materializes known non-graph hydrogens and compares element, isotope, charge,
radical state and typed bonds, including directional dative bonds. Atom maps and
hydrogen storage declarations do not change ligand identity. Iterative search
avoids recursion proportional to input size and reports resource exhaustion.

An orientation-reversing mapping must preserve other active stereo constraints.
Specified absolute centers may exchange when their mapped configurations agree;
unspecified sites are held fixed with their local orientation to retain dependent
ring stereo. Removing a site releases its constraints in subsequent rounds.
Enhanced groups and represented axes are conservatively fixed pending their full
dependency model. This is a documented limitation, not complete group/axis support.

The public `detect_stereo_candidates` API is now fallible. The options-bearing
variant bounds analysis vertices (including detached hydrogens) and total trial
mappings. Coordinate stereo inference uses the same prepared perception and
candidate set. Existing compatible, complete valence/aromaticity perception is
reused; otherwise preparation occurs on a temporary clone. Analysis is read-only.

`cleanup_stereo` explicitly edits a `MoleculeEditor` transactionally and reports
removed and unclassified assertions. It removes chemically excluded tetrahedral
and double-bond assertions while the editor maintains stereo-group membership.
Unsupported tetrahedral geometries and axes are preserved and reported, with
their reference frames protected in symmetry analysis. Nitrogen inversion and
conjugation are not yet modeled fully; treating every unsupported three-coordinate
center as invalid would erase potentially valid bridgehead assertions. Ordinary
perception continues to preserve represented stereo as required by architecture.

The `stereo.perception` benchmark explicitly invokes cleanup before and after
coordinate materialization. The reference already cleans represented stereo
during preparation. Other feature adapters and the independent reference remain
unchanged. Raw graph, surviving stereo, group, candidate and coordinate comparisons
remain asserted. No golden values were regenerated for this stage.

Regressions cover equivalent branches, ring symmetry, specified and unspecified
stereo dependencies, homomorphic/enantiomorphic ligands, explicit hydrogen and
atom-map invariance, atom/bond renumbering, transactional resource failure, group
pruning, aromatic double-bond cleanup and preservation of unclassified geometry.
The old coordinate-axis unit fixture used manually assigned aromatic flags on an
acyclic saturated graph. It now uses actual sp2 endpoint chemistry with the same
geometry, priority relationships and orientation assertions.

The initial full comparison used `stereo-symmetry-verified-bench.exe`, SHA-256
`f3ede0e6e54be945640dcc6e1bb61e032e94d2366d4c83654d02f302f1eb6088`.
Across all 301,834 rows it yielded 262,665 agreements, 38,142 disagreements,
26 errors and 1,001 not-applicable rows, improving agreement by 51,017 versus
the previous local-eligibility stage. Eight new search exhaustion errors (four
structures in both SDF and SMILES) prevented considering this implementation
finished. Exhaustive paired audits allowed only candidate/assertion removals and
corresponding group pruning, proving other observations and all expectations
unchanged (`stereo-symmetry-*-audit.json`).

There were 27 newly disagreeing PubChem rows and one Enamine row, apart from the
resource errors. Independent RDKit self-substructure mappings now provide an
orientation-reversing graph-automorphism witness for every removed candidate in
those 27 PubChem rows, preserving the other stereo constraints. Source graphs
were verified against the complete reference observations before the proof search.
Absolute centers may exchange with the correct mapped parity; ordinary terminal
hydrogens outside stereo frames were distinguished solely to restrict redundant
proof enumeration. All 27 searches finished below their enumeration cap. Evidence:
`prove_stereo_symmetry_regressions.py` and
`stereo-symmetry-pubchem-regression-proofs.{json,log}` in the review directory.
The Enamine witness is in `stereo-symmetry-enamine-regression-proof.json` and its
general ring/exocyclic-equivalent-ligand rule has a focused unit regression.
These raw benchmark differences are retained rather than copying reference
false-positive candidates into the native implementation.

#### Removing redundant symmetry search

The exhausted structures exposed permutations of equivalent terminal hydrogens
unrelated to stereo. Search now prioritizes active stereo frames and retains one
unused target representative among same-colored vertices having identical typed
neighborhoods. Every stereo focus, carrier and unsupported reference anchor is
excluded from this equivalence reduction. Swapping any two retained-equivalent
vertices is itself an automorphism fixing all stereo constraints and existing
assignments, so this removes redundant search without changing its answer.

A family of small substituted bridged-cation unit fixtures exhausted 20,000
trial mappings before the change and passes at that same bound afterward. The
default million-state bound was not raised. All six stereo symmetry tests and
the full workspace suite, including 814 core tests, passed. The retained release
executable is `stereo-symmetry-pruned-bench.exe`, SHA-256
`ccf887c3bb7c19cb5f86ed694a63bc4be79c292b5c913c7e2bdf54889a736c74`.

Full reruns from this executable completed for Enamine (71,840 agree, 28,640
disagree), PL-REX (328 agree), PDB (1,000 not applicable), and stereo smoke
(14 agree, 11 disagree, one not applicable), with no errors. Every case row is
byte-identical to the preceding symmetry run
(`stereo-symmetry-pruning-unchanged-audit.json`). The full PubChem rerun completed
with 190,491 agreements, 9,491 disagreements and the 18 pre-existing errors.
All eight search-exhausted rows now agree with the reference. Exhaustive paired
comparison verified the other 199,992 PubChem rows byte for byte and verified
unchanged inputs and expectations for the eight repaired rows
(`stereo-symmetry-pruning-pubchem-audit.json`). Final totals for this stage are
262,673 agreements, 38,142 disagreements, 18 errors and 1,001 not-applicable rows:
51,025 additional agreements compared with the local-eligibility stage.

The all-feature smoke run completed all 578 rows. Its paired audit against the
previous local-eligibility stage found three changed stereo-perception cases
(four candidate removals), unchanged outcomes, and no other observational change
(`stereo-symmetry-pruned-all-smoke-audit.json`).

Workspace all-target/all-feature checking, warnings-denied clippy and Rustdoc,
Rust 1.89 checking, verified core packaging, companion/benchmark package file
lists, packaged license copies, no-default-feature potentials tests/docs, formatting and diff whitespace
checks passed. Logs use `stereo-symmetry-pruning-*`. Earlier workspace tests,
documentation and packaging attempts failed when the disk filled. Old generated
build directories were removed without touching fixtures or benchmark evidence.
Rustdoc then exposed a truncated generated `src-files.js`; rebuilding the generated
documentation directory passed (`stereo-symmetry-pruning-doc-retry.log`).
Linux/fuzz checks remain unavailable on this Windows host. Full companion package
verification remains deferred until foundational-crate publication, per CI policy.
Unchanged Python reference, dashboard and atomic-data tests were not repeated in
this stereo-only stage; their preceding successful checks remain recorded above.

The overall goal remains active. Remaining work includes the
nitrogen/group/axis/drawing models, CIP
expansion, scientific diagnostics, independent ring checks, broader external query
coverage and final fixed-revision validation. Recurring Windows benchmark report
publication failures also remain unresolved; DSSP stays deferred.

### Closed-shell nitrogen stereo eligibility

Three-coordinate neutral, nonradical nitrogen now participates in the shared
tetrahedral model when its three single-bond ligand directions are unconjugated
and it belongs to a three-membered ring or satisfies the RDKit-like bridgehead
criterion. Three-membered rings are recognized directly from adjacency. The
bridgehead convention requires at least three incident ring bonds and requires
every selected ring containing the atom to share at least two bonds with another
such ring. Counting three ring bonds alone would also admit ordinary fused-ring
junctions. This is a discrete inversion-eligibility convention, not a calculation
of configurational lifetime or an inversion barrier.

Conjugation checks reuse the existing RDKit-like outer-electron and pi-electron
rules, accounting for neutral hypervalence, the model's higher-row donor limits,
available electrons, total coordination and an adjacent multiple/aromatic bond.
There is no list of amide, enamine, aromatic-amine or corpus-specific patterns.
The pinned reference sources are `FindStereo-2026.03.3.cpp`,
`ConjugHybrid-2026.03.3.cpp` and `QueryOps-2026.03.3.cpp` in the review directory.
Their public upstream locations are under
`https://github.com/rdkit/rdkit/tree/Release_2026_03_3/Code/GraphMol`.

Prepared stereo analysis now requires a complete installed ring set as well as
valence and aromaticity; incomplete perception is prepared on a temporary copy.
Source parsing and ordinary perception still preserve represented stereo.
Explicit cleanup can now distinguish unsupported geometry from neutral nitrogen
that is classified as inversion-labile or conjugated, and removes only the latter
when chemically ineligible. Symmetry analysis no longer treats those invalid
assertions as immutable reference frames. Supported nitrogen sites participate in
the existing general stereo-dependency analysis, including dependent carbon sites.

Hydrogen ligands have identical meaning in graph, declared and inferred forms.
An independent RDKit probe demonstrates a remaining convention difference:
`N1CC1C` acquires a potential nitrogen stereocenter after explicit hydrogen
expansion, despite unchanged chemical composition. Native perception includes it
consistently in both forms. Charged/radical three-coordinate nitrogen is outside
this closed-shell lone-pair model. Existing assertions remain unclassified and
preserved: the singly occupied orbital of a three-coordinate radical cation must
not silently be represented as a lone pair. Four-coordinate nitrogen remains
handled by the existing four-ligand geometry.

Regressions cover ordinary amines, aziridines with carbon or hydrogen ligands,
larger monocyclic amines, conjugated carbonyl/alkenyl/cyano/aryl attachments,
higher-row multiple-bond neighbors, bridged and fused rings, ammonium geometry,
explicit cleanup and unclassified radicals. Hydrogen expansion and atom/bond
renumbering are checked. The bridged nitrogen/carbon dependency has an explicit
two-site regression. Independent reference probes for 16 inputs in both hydrogen
representations are saved in `stereo-nitrogen-reference-probe.json`.

The retained release executable is `stereo-nitrogen-bench.exe`, SHA-256
`206e6ada02408c1925fb0f2582040e85e40d284da6c4f3ea0247c95ad11ba364`.
Full `stereo.perception` reruns against unchanged reference files completed:

| Dataset | Agree | Disagree | Errors | Not applicable |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 191,296 | 8,686 | 18 | 0 |
| Enamine | 71,955 | 28,525 | 0 | 0 |
| PL-REX | 328 | 0 | 0 | 0 |
| PDB | 0 | 0 | 0 | 1,000 |
| Smoke | 14 | 11 | 0 | 1 |

This adds 920 agreements with no additional errors. The exhaustive Enamine audit
verified all 100,480 rows: 115 disagreements became agreements, no agreement
regressed, and 144 changed cases added 152 nitrogen and 102 dependent carbon
candidates. Other observations and expectations are unchanged. PL-REX, PDB and
stereo smoke are unchanged. The all-feature smoke audit also found all 578 rows
unchanged. The full PubChem audit verified all 200,000 rows: 844 disagreements
became agreements and 39 agreements became disagreements, with the same 18
errors. All 950 changed cases contain only added tetrahedral candidates: 1,322
nitrogen, 180 carbon, 38 silicon and 14 phosphorus sites. The non-nitrogen sites
are stereo dependencies affected by recognition of the nitrogen sites. Surviving
candidates, represented stereo, groups, all other observations and expectations
are unchanged (`stereo-nitrogen-pubchem-audit.json`).

Every one of the 39 new disagreements was independently reproduced from its
original SMILES with pinned RDKit. Its complete reference component graph was
verified before comparison. After `Chem.AddHs`, RDKit's entire candidate list
matches the native list exactly in all 39 cases, including three dependent carbon
sites and 44 N-H nitrogen sites. These are the documented hydrogen-representation
convention difference, not unexplained regressions. Evidence:
`prove_nitrogen_hydrogen_invariance.py` and
`stereo-nitrogen-hydrogen-invariance-proofs.json`. The native representation-invariant
behavior and all raw disagreements are retained. Final stage totals are 263,593
agreements, 37,222 disagreements, 18 errors and 1,001 not-applicable rows.

Workspace all-feature tests and doctests passed on retry, including 816 core and
59 benchmark unit tests. The first attempt reached the missing-goldens integration
test and failed with Windows access denied during benchmark execution; the
isolated retry passed. Both logs are retained. This recurrence is not considered
a repair of the previously documented benchmark-publication defect and must be
investigated independently. The later focused dependency regression also passed
all eight stereo symmetry tests (`stereo-nitrogen-dependency-tests.log`).

Formatting, whitespace, workspace all-target/all-feature checking, warnings-denied
clippy and Rustdoc, Rust 1.89 checking, verified foundational-crate packaging,
companion/benchmark package lists and no-default-feature potentials tests/docs
passed (`stereo-nitrogen-*` logs). Unchanged Python reference, dashboard and atomic
data checks were not repeated. Linux/fuzz checks remain unrun on this Windows
host, and full companion packaging remains deferred until foundational-crate
publication under the CI policy. No runtime reference dependency, changed golden,
comparison exception, corpus-specific rule or README edit was introduced.

The goal remains active. Group/axis/drawing stereo, charged/radical nitrogen
geometry, CIP resource work, scientific diagnostics, ring validation, broader
external query coverage, reliable report publication and final fixed-revision
validation remain. DSSP stays deferred.

### Windows report publication under transient reader locks

Reproduced the historical Windows access-denied failure with an ordinary reader
that omits `FILE_SHARE_DELETE`. The new regression failed before the fix with
the same OS error 5 as the intermittent missing-goldens integration failure.
Windows requires delete sharing for concurrent rename/replacement; see the
[CreateFile sharing contract](https://learn.microsoft.com/en-us/windows/win32/api/fileapi/nf-fileapi-createfilew).
The identity of the process behind earlier incidental locks is not established.

Publication now retries Windows access/sharing/lock errors (5, 32 and 33) with
bounded backoff for up to two seconds. Only the final hard-link/rename operation
is retried, using the already serialized and synchronized staging file. Chemistry,
reference evaluation, case output and serialization never rerun. Other errors
return immediately. First publication still cannot replace another run's report,
and a persistently blocked update preserves the previous complete snapshot bytes.
No new dependency or platform-specific unsafe code was introduced.

Report creation, synchronization and publication errors now identify their
operation and destination path. If an execution failure is followed by failure
to publish its diagnostic report, the CLI preserves both errors rather than
replacing the original missing-goldens or computation error with an opaque I/O
failure. The regression for a persistent reader lock checks this combined
diagnostic and verifies that the old snapshot survives. The temporary-lock
regression verifies eventual publication and exactly one serialization call.

All benchmark unit/integration tests passed, including the 60 benchmark unit
tests and the existing collision, partial-write, missing-reference and history
tests (`report-retry-tests.log`). The formerly flaky missing-goldens CLI test
also passed 20 complete repetitions; each still exercises both missing and
corrupt reference files (`report-retry-missing-goldens-*.log`). These observations
verify recovery for the reproduced transient-lock mechanism, not immunity to
permanent filesystem permissions or locks exceeding the bounded retry period.

The retained executable is `report-retry-bench.exe`, SHA-256
`bcfacba15dd29546d3031b85ee17f639886aea56d1ff51797681d5451d391f44`.
The real 25-feature smoke run completed while a separate Windows reader acquired
44 short-lived, replacement-denying locks on its report. Its scientific exit
status remains failure for the existing comparison disagreements, while the
report is complete with no execution error. All 578 case rows are byte-identical
to the preceding nitrogen-stage smoke run. Evidence is
`run_locked_report_smoke.py`, `report-retry-locked-smoke*` and
`report-retry-smoke-audit.json` in the review directory. Existing empty history
archives were retained; the dashboard continues to identify and skip them.

Workspace all-target/all-feature check, warnings-denied clippy, all-feature tests
and doctests (including 816 core tests), warnings-denied Rustdoc, Rust 1.89 check,
formatting, whitespace checks and benchmark package file listing passed
(`report-retry-{check,clippy,workspace-tests,doc,msrv,package}.log`). Runtime crates,
reference adapters, atomic data and dashboard code are unchanged, so their
separate packaging/no-default-feature/Python/Node checks were not repeated; the
preceding validations remain recorded above. Linux/fuzz checks remain unrun on
this Windows host. Full external chemistry corpora were not rerun for this
publication-only change: the complete smoke comparison plus controlled filesystem
regressions exercise the changed behavior without recomputing unchanged chemistry.

The overall goal remains active. Group/axis/drawing stereo, charged/radical
nitrogen geometry, CIP expansion, scientific convention diagnostics, independent
ring validation, broader external query coverage and final fixed-revision
validation remain. DSSP stays deferred.

### CIP retained constitutional expansion and budget boundaries

Constitutional ranking now retains its rooted ligand trees while increasing
depth. Only the unexpanded boundary retains molecular paths; each occurrence
keeps its own root/path/duplicate semantics. Existing nodes and priorities are
not reconstructed or counted again. Auxiliary descriptor construction and its
path-dependent semantics are unchanged. Comparisons use increasing depth
intervals to amortize work on large tied trees, with an additional comparison
before either the depth or node bound prevents further expansion.

This fixes a general premature-exhaustion defect: depth doubling could pass a
distinguishing shell and then fail the node budget before attempting to compare
that shell. The new regression uses two atom orderings of a tetrahedral center
with oxygen- and carbon-terminated chains. Its ten-node assignment failed before
the fix and succeeds afterward; a nine-node bound still fails transactionally.
Existing tests retain the prohibition against resolving truncated constitutional
ties through isotope/auxiliary rules or declaring them nonstereogenic. No default
limit changed: depth 32 and 100,000 nodes remain.

The final executable is `target/general-chemistry-review/cip-retained-bench.exe`,
SHA-256 `4420b69b6aecc53cf7f7260a6a8e312288345cf0fb054caba66dcf505a625acb`.
Fresh independent RDKit 2026.03.3 CIP expectations are frozen in
`reference-candidates/cip-expansion`; tracked goldens were not replaced. The
preceding `report-retry-bench.exe` and final executable were compared against
those same expectations. All 301,834 case observations are identical:

| Dataset | Agree | Disagree | Error | Not applicable |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 199,590 | 288 | 122 | 0 |
| Enamine | 100,480 | 0 | 0 | 0 |
| PL-REX | 328 | 0 | 0 | 0 |
| PDB | 0 | 0 | 0 | 1,000 |
| Smoke | 25 | 0 | 0 | 1 |

The final 25-feature smoke rerun also preserves all 578 observations. Reports
are `cip-retained-all.json` and `cip-retained-all-smoke.json`; exhaustive evidence
is `audit_cip_retained.py` and `cip-retained-audit.json`. The initial implementation
that compared every shell also preserved the full corpus; its results remain as
`cip-expansion-after-*`, with its executable retained in `cip-every-shell`.
Concurrent diagnostic runs make these timings unsuitable for claiming a speedup.

All 104 historically exhausted inputs were remeasured. Their expectations were
independently verified identical to the fresh reference rows. At depth 32, all
still return explicit depth exhaustion. An initial serial depth-64 diagnostic
completed 56 inputs before being stopped during an expensive subsequent input;
its complete preceding rows remain in `cip-expansion-depth64.jsonl`. It did not
establish an outcome for the interrupted input.

The final depth-64 feasibility probe then measured **every** input in its own
process, retaining the 100,000-node bound and applying a diagnostic-only 30-second
wall limit with four concurrent processes. Results: 92 agreements, four
disagreements, two node-limit failures, and six timeouts. Timeouts are unresolved
measurements under that concurrency/time policy, not chemical disagreements or
proofs that unlimited evaluation cannot finish. This bound is not installed in
the library or benchmark. Evidence includes `run_cip_limits_bounded.py`,
`cip-retained-depth64-bounded{.jsonl,.meta.json,.log}`, the depth-32 probe, and
`cip-retained-limits-audit.json`.

The two node-limit inputs are PubChem 158374 and 163705. Three deeper-run
disagreements (122322, 123741, 195709) contain a carbon assertion that RDKit drops
with an overlapping-neighbor drawing warning. The fourth (445597) contains four
native Z labels on bonds RDKit represents as aromatic without bond stereo. All
other atom and bond labels match for these four cases; the observations remain
unchanged. `inspect_cip_limit_disagreements.py` and
`cip-retained-depth64-differences.json` record those distinctions. The timeout
inputs are 4732, 73148, 131910, 157414, 176534 and 441915. These results do not
justify raising the production defaults in this stage.

Final-source formatting, workspace all-target/all-feature check, warnings-denied
clippy, all-feature workspace tests/doctests (817 core tests), warnings-denied
Rustdoc, Rust 1.89 check, full dirty-tree core package verification, no-default
potential tests/docs and whitespace checks passed (`cip-retained-*.log`).
Companion/benchmark package file listings passed earlier in this stage and were
not repeated for the subsequent private comparison-scheduling change; their
package contents are unchanged. Production Python, dashboard, atomic data and
license files are unchanged, so their separate Python/Node/generator/license
checks were not repeated. Linux/fuzz checks remain unrun on this Windows host.

CIP work is not complete: selecting only unresolved comparisons for expansion
and reducing repeated work in tied/auxiliary ranking still require a correctness
argument that preserves path-dependent descriptors. The goal also retains
group/axis/drawing stereo, charged/radical nitrogen geometry, scientific
diagnostics, independent ring validation, broader external coverage and final
fixed-revision validation. DSSP stays deferred.

### CIP descendant-order reuse within immutable comparisons

The sequence-rule comparator now reuses descendant sorting within one immutable
comparison. Cache keys identify individual ligand-tree occurrences and the
sequence rule; the selected Rule 6 reference belongs to the comparison context.
Consequently, different paths to the same molecular atom cannot share their
duplicates, descriptors or orderings. Cached references cannot survive a tree
expansion or molecular edit. Pointer identities are keys only and are never
dereferenced. Storage grows with visited branching occurrences and the fixed
number of sequence rules, rather than arbitrary pairs of ligand nodes.

The sequence rules and breadth-first comparison order are unchanged. Unreferenced
Rules 4b and 5 still do not independently choose references while sorting
descendants. New regressions cover constitutional ties in branched carbon
ligands, a distinguishing isotope in one repeated atom occurrence, reversal of
child storage order, and different Rule 6 references. All existing path-dependent
auxiliary, ring-duplicate, fractional-number and pseudoasymmetry regressions pass.

The retained executable is `cip-comparison-bench.exe`, SHA-256
`31a27c75081458fadd7286af87f8ca884a31c282c402e3b9098a218fa53b7b9a`.
The full CIP run uses the same frozen `reference-candidates/cip-expansion`
expectations. `audit_cip_comparison.py` verifies all 301,834 rows unchanged, and
all 578 observations in the complete 25-feature smoke rerun unchanged. Counts
remain 300,423 agreements, 288 disagreements, 122 errors and 1,001 not-applicable
rows. Reports and exhaustive evidence are `cip-comparison-all.json`,
`cip-comparison-all-smoke.json`, and `cip-comparison-audit.json`.

Every historical exhaustion input was measured again at depths 32 and 64.
Depth 32 still exhausts on all 104. The same diagnostic-only depth-64 policy
(100,000 nodes, 30 seconds per process, four concurrent processes) gives 91
agreements, four disagreements, two node-limit failures and seven timeouts.
All 96 cases completed in both bounded runs have identical labels/errors.
73148 now completes; 441921 and 461697 now time out. The other five timeout
inputs remain unchanged. This variation is a resource measurement, not changed
chemistry, and does not establish a performance regression or improvement.
Evidence is `cip-comparison-depth64-bounded*`, `cip-comparison-limits-audit.json`
and `cip-comparison-depth64-transition-audit.json`.

Paired sequential timings then ran both implementations twice, alternating their
order, on five external inputs. All assignments match, including 441921 and
461697 completing with both executables. Observed median assignment seconds:

| PubChem input | Retained expansion | With comparison reuse |
| --- | ---: | ---: |
| 43217 | 0.538 | 0.532 |
| 60348 | 0.452 | 0.426 |
| 114830 | 7.844 | 7.700 |
| 441921 | 10.552 | 10.100 |
| 461697 | 13.729 | 13.430 |

These modest reductions on a small measured subset do not support a broad
speedup claim or solve the expensive auxiliary cases. The comparison script,
binary hashes, complete outcomes and individual timings are retained in
`compare_cip_runtime.py` and `cip-comparison-runtime.{json,log}`. Resource defaults
remain unchanged.

Formatting, workspace all-target/all-feature check, warnings-denied clippy,
workspace tests/doctests (819 core tests), warnings-denied Rustdoc, Rust 1.89,
full dirty-tree core package verification, companion/benchmark package listings,
no-default potential tests/docs and whitespace checks passed
(`cip-comparison-*.log`). Python/Node/atomic-data/license checks were not repeated
because those production files are unchanged. Linux/fuzz checks remain unrun on
this Windows host.

The next concrete CIP investigation is auxiliary traversal: the current code
reconstructs the same root-to-ancestor path at every visited node, and eagerly
materializes trees for auxiliary ranking. Reusing a traversal's orientation and
expanding only unresolved comparisons must preserve the original occurrence
paths and auxiliary dependency order. This stage does not complete that work or
the remaining stereo, scientific-diagnostic, ring, coverage and final-revision
requirements. The overall goal remains active; DSSP stays deferred.

### CIP auxiliary traversal orientation

Auxiliary ranking now computes its original ancestor path once per ranking root.
An immutable traversal records the edges directed toward that root and excludes
them while walking away from it. Neighbor iteration borrows the existing child
arrays rather than allocating and reconstructing the ancestor path at every
node. Tetrahedral carriers and each bond/axis endpoint share their traversal.
Original ligand occurrences, duplicate nodes, stored paths, neighbor order and
the ordering of auxiliary descriptor batches are unchanged.

The new regression builds a branched cyclic ligand with ring duplicates, roots
its auxiliary tree at every occurrence, and compares each outgoing sequence
against an independent undirected traversal that excludes the incoming edge.
Every occurrence is reached exactly once from every root. Existing CIP tests
also pass, including path-dependent auxiliary labels and equivalent-carrier
ordering invariance. No default limits or chemical eligibility rules changed.

The retained benchmark executable is `cip-traversal-bench.exe`, SHA-256
`8e62455c3b147f51c20199cec8fe39dd379254e1762cbc8be3da6ac657d72166`.
The standalone diagnostic executable is `cip-traversal-probe.exe`, SHA-256
`efd19ca3aa6adc664940354a098edcb946dfcf56ae4046cc195145f3da96bd79`.
The successful full comparison is `cip-traversal-recheck-all.json`, using the
unchanged frozen `reference-candidates/cip-expansion` expectations. All 301,834
case observations and all 578 complete smoke observations are unchanged.
Totals remain 300,423 agreements, 288 disagreements, 122 errors and 1,001
not-applicable inputs. `audit_cip_traversal_recheck.py` and
`cip-traversal-recheck-audit.json` verify the complete rows, not just totals.

The depth-32 probe still reports depth exhaustion for all 104 historical inputs.
The completed depth-64 feasibility run retains the same 100,000-node, 30-second,
four-process policy: 91 agreements, four disagreements, two node-limit failures
and seven timeouts. All 94 cases completed in both this and the preceding
comparison-cache run have identical labels/errors. Timeout membership changed
while compilation and filesystem compression were active, so those wall limits
are not comparative performance evidence. Full outcomes and the transition audit
are `cip-traversal-recheck-depth64-bounded*`,
`cip-traversal-recheck-limits-audit.json` and
`cip-traversal-depth64-transition-audit.json`.

After builds, corpus runs and compression finished, both executables ran twice
on the same five external inputs, alternating order. All assignments matched.
Observed median assignment seconds were:

| PubChem input | Comparison cache | Reused auxiliary traversal |
| --- | ---: | ---: |
| 43217 | 0.533 | 0.447 |
| 60348 | 0.470 | 0.418 |
| 114830 | 8.175 | 7.135 |
| 441921 | 10.627 | 9.775 |
| 461697 | 13.957 | 11.972 |

These are 8–16% reductions on this measured subset, not a universal speedup
claim. The script, executable hashes, complete assignments and individual times
are `compare_cip_traversal_runtime.py` and `cip-traversal-runtime.{json,log}`.

Disk exhaustion interrupted the first validation attempt: test compilation and
package verification failed, and the original corpus/diagnostic files are
incomplete. They remain under the original `cip-traversal-*` names and must not
be treated as passed runs. NTFS compression of the completed historical
`target/benchmark-parity/chem-hydrogens-before.cases.jsonl` preserved its
14,672,150,560 logical bytes in 3,189,518,336 allocated bytes. No fixtures,
reports or retained executables were deleted. The successful reruns use new
`cip-traversal-recheck-*` paths. Validation wrappers now stop on command-launch
errors rather than reading a stale successful exit code. Compression completed
successfully; its transcript is `cip-traversal-evidence-compression.log`.

Formatting, all-target/all-feature workspace check, warnings-denied clippy,
workspace tests/doctests (820 core tests), warnings-denied Rustdoc, Rust 1.89,
full dirty-tree core package verification, companion/benchmark package listings,
no-default potential tests/docs and whitespace checks passed. Successful rerun
logs for the affected gates are `cip-traversal-recheck-{test,doc,package}.log`;
the other completed checks use `cip-traversal-*.log`. Production Python,
dashboard, atomic-data and license files are unchanged, so their separate
Python/Node/generator/license checks were not repeated. Linux/fuzz checks remain
unrun on this Windows host.

Auxiliary ranking still eagerly materializes whole ligand trees. The next CIP
step is to reuse bounded incremental expansion for this graph view, stopping
only when a sequence-rule comparison is conclusive and preserving the original
occurrence paths and descriptor dependency order. The broader stereo,
scientific-diagnostic, ring, coverage and fixed-revision requirements remain;
the overall goal is active and DSSP remains deferred.

### Incremental auxiliary CIP expansion and measured default depth

Constitutional and auxiliary ranking now share one private bounded expansion
driver. Auxiliary ranking grows views of the existing occurrence graph one shell
at a time, retaining the frontier, and stops once every carrier pair has a proven
order. Original occurrence identities, paths, duplicate nodes, root orientation,
auxiliary descriptor batches and Rule 6 handling are preserved. An unfinished
constitutional tie cannot advance to isotope or auxiliary rules or establish
nonstereogenicity. Bounds are checked before adding a shell, after first trying
to resolve its current signatures. This eliminates eager reconstruction of
complete auxiliary ligand trees; it introduces no molecule-specific rules.

The regression `auxiliary_ranking_stops_at_a_proven_order_but_never_promotes_a_truncated_tie`
first failed against the eager implementation: distinct carrier elements could
not be ranked with one node per ligand (`cip-lazy-aux-regression-before.log`).
It now agrees with complete expansion under that bound, while a separate
isotopically distinguished pair with an unfinished constitutional tie still
returns a depth error. The shared driver expands an unresolved carrier set
together; this is not individual pair freezing, and construction of the main
auxiliary graph remains bounded but eager.

Before changing defaults, the optimized depth-32 executable was run against the
complete frozen CIP corpus. All 301,834 observations were unchanged. Its retained
executable is `cip-lazy-aux-bench.exe`, SHA-256
`43edb5d8b9bf3b6f4d2dfca71afb47402d195c63e2e87ae14603b722989fc2a2`.
The intermediate report is `cip-lazy-aux-depth32-all.json`.

All 104 historical depth failures were then measured at depth 64 with the
unchanged 100,000-node bound and previous diagnostic policy of four processes
and 30 seconds per process. Every process completed: 98 agreements, four
disagreements, two explicit node-limit errors, no timeouts. The slowest measured
assignment was 0.799 seconds. Evidence is
`cip-lazy-aux-depth64-bounded.{jsonl,meta.json,log}`. This supplied-corpus evidence
justifies raising the default depth from 32 to 64 without raising the node cap.
A new public regression resolves a remote constitutional difference at default
depth, verifies transactional failure at depth 32, and independently matches
RDKit's R assignment. Existing large-equivalent-ligand regressions scale their
chain lengths with the default so they retain their original boundary coverage.

The final retained executable is `cip-depth64-bench.exe`, SHA-256
`4ee91cc75519542ed18b24501b6794b3f4a7425918172c7021ba427771cfe78d`.
The final diagnostic executable is `cip-depth64-probe.exe`, SHA-256
`1af18d64da8fcf7747681ba5d4525753966817ea599a4a18ab0be8fec1d3bd35`.
Its sequential run over all 104 historical inputs, without the diagnostic wall
cutoff, reproduced every bounded-run assignment and error
(`cip-depth64-all-historical.{jsonl,log}`). Frozen reference expectations were
independently rechecked against all 104 source manifest entries.

The final full benchmark `cip-depth64-all.json` is complete and retains the same
RDKit 2026.03.3 expectations. Results are:

| Dataset | Agreements | Disagreements | Errors | Not applicable |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 199,688 | 292 | 20 | 0 |
| Enamine diversity | 100,480 | 0 | 0 | 0 |
| PL-REX | 328 | 0 | 0 | 0 |
| PDB | 0 | 0 | 0 | 1,000 |
| Smoke | 25 | 0 | 0 | 1 |
| Total | 300,521 | 292 | 20 | 1,001 |

The exhaustive three-stage audit verifies that optimization alone leaves all
301,834 rows unchanged. Raising the depth changes exactly the historical 104
error rows: 98 become agreements, four become disagreements, and two remain
errors with explicit node exhaustion. All other 301,730 complete observations
are unchanged, including their expectations. All 578 all-feature smoke
observations are also unchanged. The audit is `audit_cip_depth64.py`, with
`cip-depth64-audit.{json,log}`. The first audit attempt had an assertion mistake
that counted not-applicable inputs inside the report's applicable `cases` field;
the corrected audit verifies 300,833 applicable plus 1,001 not-applicable rows.
No benchmark comparison or expected result was changed.

PubChem inputs 158374 and 163705 still reach the independent node bound. The
four newly visible disagreements remain source-stereo interpretation issues:
122322, 123741 and 195709 contain extra native S carbon labels where RDKit
rejects overlapping-neighbor drawings; 445597 has four native Z aromatic-bond
labels omitted by RDKit. Other descriptors agree. These remain visible for the
general drawing/aromatic stereo investigation. The other 18 reported errors
are pre-existing reference/parse failures.

After builds and full benchmarks finished, both retained probe executables ran
twice on each of five supplied difficult inputs, alternating order. Every
assignment matched. Median assignment seconds were:

| PubChem input | Previous auxiliary traversal | Incremental auxiliary ranking |
| --- | ---: | ---: |
| 43217 | 0.424 | 0.00843 |
| 60348 | 0.385 | 0.01361 |
| 114830 | 6.711 | 0.05187 |
| 441921 | 9.721 | 0.09936 |
| 461697 | 11.956 | 0.07017 |

These timings characterize this subset only. Individual measurements, outputs
and executable hashes are retained in `compare_cip_lazy_aux_runtime.py` and
`cip-lazy-aux-runtime.{json,log}`.

Formatting, workspace all-target/all-feature check, warnings-denied clippy,
workspace tests/doctests (822 core tests), warnings-denied Rustdoc, Rust 1.89,
full dirty-tree core package verification, companion/benchmark package listings,
no-default potential tests/docs and whitespace checks passed. Logs use the
`cip-depth64-*` prefix. Production Python, dashboard, atomic-data and license
files did not change in this stage, so their separate Python/Node/generator/
license checks were not repeated. Linux/fuzz checks remain unrun on this
Windows host. Full benchmark and smoke comparisons intentionally retain exit
status 1 for the scientific differences described above; both reports are
complete with no run-level error.

This completes the measured auxiliary-expansion/default-depth stage, not the
overall goal. Remaining node-bound cases, drawing/group/axis and charged/radical
nitrogen stereo, scientific mass/hydrogen/ring/mmCIF diagnostics, broader supplied
coverage and final fixed-revision validation remain to be addressed. DSSP stays
deferred.

### Ambiguous tetrahedral source drawings (validation in progress)

The preceding turn completed the incremental auxiliary CIP/depth stage and was
progress. This stage investigates remaining source stereo differences. The four
aromatic-bond labels on PubChem 445597 are preserved: CIP ranks represented
assertions, whereas aromaticity belongs to perception. Suppressing them inside
CIP would violate that ownership boundary. General cleanup/convention work
remains separate.

Source wedge decoding had a concrete defect: missing or degenerate coordinates
could fall back to atom-list parity and invent a configuration. That fallback is
removed. Single and redundant wedges now use the same checked drawing path;
failed geometry produces the existing source warning even for one mark, without
creating a specified or explicitly-unknown placeholder. Wavy marks still retain
their explicit unknown semantics independently of geometry.

Drawing calculations translate to the center and scale uniformly before adding
the wedge displacement. Unit-vector separation detects overlapping bond
directions independently of coordinate units and radial bond lengths. The
0.001 squared-separation convention follows RDKit's
[drawing interpretation](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Chirality.cpp).
The wedge's out-of-plane direction is applied before overlap checking, so a
projected overlap legitimately disambiguated by wedging remains supported.
Zero-length unmarked directions and zero chiral volumes cannot create labels.
No supplied molecule identifiers enter the implementation.

Coordinate-free Molfile writers previously depended on the same fallback.
They now reject specified stereo and require a model with suitable coordinates.
Models with degenerate emitted drawings also reject lossy export. Existing
stereo round-trip tests now retain/provide their drawings, preserving their
carrier, parity, hydrogen and group assertions. Public writer documentation and
the benchmark guide describe this behavior. Graph ownership and ordinary
perception remain unchanged; this is not a coordinate-generation feature.

The new ambiguous-drawing regression failed before the change
(`stereo-drawing-regression-before.log`). Native tests cover 60 invalid and 48
valid V2000/V3000 drawings across reflections, translations and uniform scales,
plus missing-geometry normalization and coordinate-free/degenerate export.
Pinned RDKit 2026.03.3 independently confirms all 108 drawing outcomes
(`prove_stereo_drawing.py`, `stereo-drawing-reference.{json,log}`). Ordinary
geometric round trips and explicit-unknown behavior remain covered. This does
not claim to implement all RDKit drawing heuristics or resolve every wedge,
axis, group or nitrogen issue.

Workspace tests/doctests passed, including 826 core tests
(`stereo-drawing-workspace-tests-3.log`). Earlier validation caught obsolete
coordinate-free test fixtures and error-message expectations; those were
corrected without removing stereo assertions. Formatting, workspace all-target/
all-feature check, warnings-denied clippy and Rustdoc, Rust 1.89 checking, core
package verification, companion/benchmark package listings, no-default potential
tests/docs and whitespace checks passed (`stereo-drawing-*` gate logs).
Production Python/dashboard/atomic-data/license code did not change in this
stage, so their separate Python/Node/generator/license checks were not repeated.
Linux/fuzz checks remain unrun on this Windows host.

The retained executable is `stereo-drawing-bench.exe`, SHA-256
`b30cb01e2acb43a8bf85691e60061de4664243b447d13cc849e9ddf19786d7ab`.
The complete full CIP comparison (`stereo-drawing-cip-all.json`) contains
300,595 agreements, 218 disagreements, 20 errors and 1,001 not-applicable rows.
Exactly 74 disagreements became agreements; all other 301,760 observations and
all expected values are unchanged. No new native error or regression occurred.
The exhaustive audit is `audit_stereo_drawing.py cip` and
`stereo-drawing-cip-audit.{json,log}`. The 20 errors retain the previous 18
reference/parse failures and two node-limit failures.

Every changed case was also read directly from its supplied SDF record with
RDKit. Its full atom-label list matches the new output; surviving atom labels
and all bond labels are unchanged. RDKit emits an ambiguous-drawing warning for
each of the 75 removed atom labels across those 74 cases. This includes the
three historical depth-recovered inputs 122322, 123741 and 195709. Evidence is
`prove_stereo_drawing_corpus.py`, `stereo-drawing-corpus-proof.json` and
`stereo-drawing-corpus-proof-2.log`; the initial diagnostic was corrected to map
source indices into component-contiguous output order for disconnected inputs.

The stage is not complete yet. Full stereo perception, all-feature smoke,
fresh independent stereo-representation generation, and full before/after
representation comparisons are still running. Expected values are generated
only by the unchanged reference adapter into a new `reference-candidates/
stereo-drawing` directory; previous goldens and reports remain untouched.
Windows directory enumeration reported a stale zero-byte size for its open
temporary output; opening the file confirmed growing compressed reference data.

Continuation handles at this checkpoint:

- Session 8493: retained final executable; CIP completed, then full
  `stereo.perception` and all-feature smoke. Outputs use
  `stereo-drawing-perception-all` and `stereo-drawing-all-smoke`.
- Session 70605: independent full `stereo.representation` generation followed by
  its baseline comparison with retained `cip-depth64-bench.exe`. Outputs use
  `stereo-drawing-reference-all` and `stereo-drawing-before-all`.
- Session 52791: waits for the reference report's explicit complete state,
  runs the final full representation comparison (`stereo-drawing-after-all`),
  then audits perception, smoke and representation after their complete reports
  exist. Logs are `stereo-drawing-{perception,smoke,representation}-audit.log`.

Re-poll these existing handles before scheduling any rerun. Inspect all complete
reports, exhaustive row transitions, changed rows and any regressions before
closing this stage. Full writer-corpus reruns are not yet performed at this
checkpoint; model-writer regression coverage and the pending all-feature smoke
run do not substitute for the final fixed-revision all-feature benchmark.
The full goal remains active and DSSP remains deferred.

### Independent ring validity and numbering checks

The preceding drawing turn was progress: it removed invented source stereo,
validated the complete CIP corpus and launched the broader comparisons. Its
session 8493 has now completed. The full stereo-perception report is complete
with no run-level error: 263,666 agreements, 37,149 disagreements, 18 errors and
1,001 not-applicable rows. This is a net gain of 73 agreements, all in PubChem;
its exhaustive changed-row audit remains queued in session 52791. The complete
all-feature smoke comparison has 578 observations, each parsed row identical to
`cip-depth64-all-smoke.cases.jsonl.gz`. Its formal audit is queued in that same
session. Do not restart these completed benchmarks.

While those existing runs continued, the ring work added independent validity
regressions in `crates/kekule/src/algorithms/rings/tests.rs`. The oracle enumerates
edge subsets and recognizes a simple cycle by connectedness and degree two;
it does not call the production bridge or ring-selection algorithms.

All 772 connected simple graphs on one through five vertices are checked under
three atom orders and two bond insertion orders, also reversing stored bond
endpoints. For every variant the tests verify exact atom and bond cycle
membership, distinct real edges and vertices in every selected ring, agreement
between its cyclic atom path and bond set, absence of duplicate cycles, and
coverage of all cyclic edges. The selected ring list is deliberately not
required to be invariant: the documented model promises valid selected cycles
and coverage, not a unique minimum or full-rank cycle basis. The 4,632 variant
checks took 0.34 seconds in the first focused run.

A second regression tests all 729 assignments of single, zero and dative orders
on a four-vertex complete molecular graph. The independent oracle excludes zero
and dative edges, including combinations that disconnect the eligible subgraph,
and verifies native membership and selected-cycle coverage. Short-ring boundary
tests now cover double, triple and quadruple orders as well. No runtime
algorithm, reference adapter, comparison field, golden, or chemical rule was
changed to obtain these results. The historical PubChem 296819 alternative-ring
difference stays subject to the raw path comparison described in
`PARITY_REVIEW.md`; no particular cycle was inserted.

All 828 core tests passed (`ring-invariants-core-tests.log`). All-target/all-
feature clippy for the core crate with warnings denied, Rust 1.89 test-target
checking, formatting and whitespace checks passed (`ring-invariants-{clippy,
msrv,fmt,diff}.log`). The focused run is `ring-invariants-tests.log`.
Workspace integration/doctest, documentation and package checks were not
repeated for this test-only change; the immediately preceding drawing stage's
complete gates remain applicable to its unchanged production sources. The
unchanged Python/dashboard/atomic-data/license checks and unavailable Linux/fuzz
checks retain the exclusions documented above.

Both full ring benchmarks are now being rerun against freshly generated pinned
RDKit expectations in a separate `reference-candidates/ring-invariants`
directory. The retained `stereo-drawing-bench.exe` is used because runtime code
has not changed; its SHA-256 remains
`b30cb01e2acb43a8bf85691e60061de4664243b447d13cc849e9ddf19786d7ab`.
These runs validate supplied-corpus behavior beyond the exhaustive small-graph
regressions. They are pending, not passed.

Current continuation handles:

- 70605 remains live: stereo-representation reference generation followed by
  its baseline comparison. PubChem's 200,000 reference observations are already
  published; subsequent datasets continue.
- 52791 remains live: final stereo-representation comparison after explicit
  reference completion, then perception/smoke/representation audits.
- 95897 remains live: full `algo.rings.fast` generation and comparison followed
  by full `algo.rings.sssr` generation and comparison. Outputs are
  `ring-invariants-{fast,sssr}-{reference,all}.{json,log}`.

Sessions 66461, 72442 and 63919 completed successfully. Once the ring reports
complete, run `audit_ring_invariants.py` to compare every observation and frozen
expected value against the preserved preceding full runs. Inspect any assertion
failure or changed row rather than weakening that audit. Re-poll live handles
before scheduling reruns. Complete the queued drawing audit, ring corpus audit,
remaining scientific diagnostics/coverage/stereo work and final fixed-revision
validation before any overall completion claim. DSSP remains deferred.

### Selected-ring comparison diagnostics

The preceding ring turn was progress: its independent validity regressions
passed, while full reference generation remained in progress. This turn adds
an explanatory diagnostic to `algo.rings.sssr` differences. Raw cyclic paths,
structural/numerical differences, exactness and agreement are unchanged.

`benchmarks/src/compare/rings.rs` constructs undirected edge sets from each
reported simple path and performs sparse Gaussian elimination over GF(2).
It reports selected-cycle counts, covered-edge counts and cycle-span ranks,
plus equality of the covered edge sets and spanned spaces. It has no fixed
vertex/edge bit-mask limit. Diagnostics require successful records with matching
record indices and titles and reject malformed paths. They describe the
reported selections only; they do not establish validity in the source graph,
a complete graph cycle basis, or minimality. The benchmark still compares the
complete raw path lists and still reports alternative selections as disagreements.

Four focused regressions check alternative bases and redundant cycles, equal
edge coverage with unequal spans, missing cycles, cyclic path normalization,
record identity and malformed paths. They also assert that diagnostics cannot
turn disagreement into agreement or mutate the observed values. In particular,
K4's Hamilton cycles cover every edge but span only the even-cycle subspace;
three independent triangles span its full cycle space. This prevents edge
coverage from being mistaken for equal cycle spans.

Validation passed: all 64 benchmark unit tests and seven integration tests,
workspace/all-target/all-feature check and warnings-denied clippy,
warnings-denied workspace documentation, benchmark package-file listing,
Rust 1.89 benchmark/all-target check, formatting and whitespace checks.
Logs are `ring-diagnostics-{focused-tests,check,clippy,test,doc,package,msrv,
fmt,diff}.log`. The release build passed and was retained as
`target/general-chemistry-review/ring-diagnostics-bench.exe`, SHA-256
`9b92235bbbd1ecfacb139dcb3ff5cf3c5c2b29ebad65fe1a131ebc4a350e9967`.
The core runtime and public API did not change. Core/companion tests, doctests,
verified core packaging, and optional-feature documentation were not repeated:
their preceding drawing and ring gates cover unchanged sources. Python,
dashboard, atomic-data and license checks were not repeated because their
sources did not change; Linux/fuzz checks remain unavailable on this Windows
host. No README or reference expectations were changed this turn.

The complete all-feature smoke benchmark has 578 rows: 377 agreements,
73 disagreements, two errors and 126 not-applicable rows. Every parsed row is
identical to `stereo-drawing-all-smoke.cases.jsonl.gz`; no smoke selection needed
the new diagnostic. Evidence is `ring-diagnostics-all-smoke.{json,log}` and
`ring-diagnostics-smoke-audit.json`. Scientific disagreements produce exit 1;
the report is complete with no run-level error. This is a regression check,
not a claim that existing smoke disagreements have been resolved.

Current continuation handles:

- 70605 completed: fresh full stereo-representation expectations and the
  baseline comparison are published. The baseline comparison exited 1 for
  retained scientific disagreements.
- 52791 remains live: the final full stereo-representation comparison has also
  finished (`stereo-drawing-after-all.json`, complete, no run-level error).
  It is now running the queued exhaustive perception/smoke/representation
  audits. Inspect those audits before claiming absence of regressions.
- 95897 remains live: full fast-ring reference generation continues, followed
  by its comparison and full selected-ring generation/comparison. The completed
  PubChem fast-ring reference has 200,000 rows and no reference errors.
- 16937 waits for an explicitly complete/error-free
  `ring-invariants-sssr-all.json`, then runs the full selected-ring comparison
  using the new retained diagnostic executable and the same freshly generated
  goldens. It subsequently runs `audit_ring_diagnostics.py all`. The audit
  removes only the additive diagnostics and requires every remaining complete
  observation, expected value and status to equal the retained baseline.
  Outputs are `ring-diagnostics-sssr-all.{json,log}` and
  `ring-diagnostics-all-audit.{json,log}`. Its reproducible wrapper is
  `target/general-chemistry-review/run_ring_diagnostics.ps1`.

Sessions 86120, 64843 and 30441 finished; the MSRV check also finished.
Do not restart existing jobs. Once the original ring runs finish, inspect the
historical `audit_ring_invariants.py` separately from the diagnostic audit:
it checks an older chemistry/reference stage and may reveal changes from the
intervening general fixes. Investigate any assertion rather than suppressing it.
The full ring rerun and drawing audits remain pending at this checkpoint.
Remaining scientific diagnostics, coverage/stereo work and final fixed-revision
validation remain part of the active goal. DSSP remains deferred.

Checkpoint update before handoff: session 52791 completed the perception and
smoke audits successfully and continues with representation. Perception has
301,760 unchanged rows and 74 changed rows: 73 disagreements become agreements,
one remains a disagreement, with no agreement regressions or new errors.
The changed input identities are exactly the same 74 supplied records already
independently checked in the drawing CIP corpus proof. The smoke audit confirms
all 578 rows unchanged. Fast-ring reference generation is also complete with
no run-level error and no reference errors; session 95897 now continues its
scheduled comparisons and selected-ring generation. The selected-ring
diagnostic rerun remains queued in session 16937.

### mmCIF text-policy diagnostics

The preceding goal turn made progress by adding selected-ring diagnostics and
completing its smoke audit. This turn adds a per-value mmCIF diagnostic for
multiline strings differing only in ASCII spaces/tabs at line ends. The parser
continues preserving source text; all raw values, exactness and disagreement
statuses remain asserted. Leading whitespace, changed text, line-break counts,
missing-value tokens and non-ASCII whitespace are not normalized. Single-line
values are excluded because the decoded comparison lacks their quoting syntax.
Diagnostics require corresponding block names, tags and value positions and
do not classify an entire input from one explained value.

The independent basis is CIF 1.1 paragraph 17, rechecked directly at
https://www.iucr.org/what-we-do/digital-standards/cif/cif1/file-syntax,
and the installed Biopython 1.87 `MMCIF2Dict._tokenize` implementation, which
strips each semicolon-text line before joining it. All 103 historically audited
value differences are multiline. No input-specific tags or molecule IDs enter
the implementation in `benchmarks/src/compare/mmcif.rs`.

Four regressions cover the whitespace relation in both directions, preservation
of original values and disagreement, exclusion of other textual differences,
coexisting unexplained changes, and matching block/column identities. All 68
benchmark unit tests and seven integration tests passed. Workspace/all-target/
all-feature check and warnings-denied clippy, warnings-denied workspace docs,
Rust 1.89 benchmark/all-target check, benchmark package-file listing, formatting
and whitespace checks passed. Logs use `mmcif-diagnostics-{focused-tests,check,
clippy,test,doc,msrv,package,fmt,diff}.log`.

The release executable is retained as `mmcif-diagnostics-bench.exe`, SHA-256
`2e264f4abd294e966563d17ff8ec764f2a3fab0a46feaac54970487978492289`.
All 578 complete smoke observations match the preceding executable exactly
(`mmcif-diagnostics-all-smoke.{json,log}` and
`mmcif-diagnostics-smoke-audit.json`). As before, existing disagreements give
exit 1; the run itself is complete without an error. Core/companion tests,
doctests, core package verification and optional-feature docs were not repeated
for this benchmark-only change; prior gates cover unchanged runtime sources.
Python/dashboard/atomic-data/license checks were not repeated because those
sources did not change. Linux/fuzz gates remain unavailable on this Windows host.

Full independent mmCIF validation is in progress, using the existing pinned
`molecular-biopython-reference` environment (Biopython 1.87) and separate new
goldens under `reference-candidates/mmcif-diagnostics`. This is not golden
adoption or replacement of existing reference values.

- Session 29562 generates the full reference, then compares the retained
  `ring-diagnostics-bench.exe` baseline. Outputs use
  `mmcif-diagnostics-{reference,before}-all.{json,log}`.
- Session 28893 waits for the complete/error-free baseline and retained new
  executable, runs the final full comparison, then audits all rows. Outputs use
  `mmcif-diagnostics-after-all.{json,log}` and
  `mmcif-diagnostics-audit.{json,log}`. Its reproducible wrapper is
  `run_mmcif_diagnostics.ps1`; `audit_mmcif_diagnostics.py` requires unchanged
  historical expected/actual values and statuses, exact baseline/final row
  identity after removing only diagnostics, and independently checks the
  whitespace relation for every added diagnostic. These checks remain pending.

The earlier drawing representation audit is now complete: 301,760 unchanged
rows, 74 changed, 68 disagreements resolved, no agreement regressions or new
errors. Its changed identities are exactly the same 74 supplied records in the
independent drawing CIP proof. Session 52791 finished successfully; do not rerun.

The full fast-ring benchmark also completed successfully: 300,833 agreements,
zero disagreements/errors and 1,001 not-applicable observations. Its historical
audit preserved every expected value, found 298,953 identical rows and 2,881
changed rows (`ring-invariants-fast-audit.{json,log}`).
All 2,881 changes are Enamine inputs previously rejected as unsupported
CXSMILES; each is now an agreement. No previous agreement regressed.
The first combined audit
reached the not-yet-created selected-ring report and stopped; the script now
accepts an explicit feature selection without relaxing its assertions, and the
fast-only audit completed in session 90417. Session 95897 continues selected-ring
reference generation/comparison; session 16937 still queues the diagnostic
rerun. Their PubChem selected-ring expectations (200,000, zero errors) are now
published. Re-poll the existing handles before scheduling further runs.

Sessions 12500, 22667 and 90417 are terminal; the MSRV and smoke commands also
finished. Remaining mass/hydrogen diagnostics, external coverage/stereo work
and final fixed-revision validation remain in the full objective. DSSP remains
deferred. The goal is active, and pending full runs are not counted as passed.

### Independent molecular-mass analysis

The preceding goal turn made progress by implementing mmCIF diagnostics and
launching their full validation. This turn preserves the verified mass data and
adds `benchmarks/analyze_masses.py`, a reusable optional analysis of completed
descriptor reports. It develops the earlier one-off independent corpus audit
into maintained tooling; it never rewrites benchmark values or tolerances.

The analyzer requires the existing checksum-pinned CIAAW/AME source cache,
regenerates and checks the current native atomic table, and requires installed
RDKit to match each applicable result's recorded reference version. It reconstructs
each engine's two masses independently from the complete isotope-resolved
formula. Native average and monoisotopic masses include the CODATA electron
correction; the reference average does not, and its exact-mass electron correction
is independently probed. Differences in atomic data and charge convention are
reported separately, with both observed values, reconstructed values, residuals
and a constituent-count summation-roundoff bound. The bound is not experimental
uncertainty or an added benchmark tolerance.

Matching formulas alone cannot yield a verified classification: each mass must
also reconstruct within the numerical bound. Formula/identity differences,
invalid observations, unavailable data and unexplained masses remain explicit.
Error and not-applicable cases stay in the summary counts and are not claimed
as reconstructed. The streaming compressed detail output identifies every
analyzed component. Report/case/analyzer/source/parser/table hashes and reference
version are retained. Incomplete reports, wrong reference versions, truncated
case files and existing output paths are rejected. No network downloads occur.

Nine dependency-free tests cover isotope/ion decomposition, corrupted masses
despite matching formulas, roundoff, formula/identity/status differences,
invalid or missing masses, invalid counts, duplicate terms, retained error/NA
counts, input immutability, reference-version checks and truncation. They pass
in the Biopython environment without RDKit. CI now runs them, and the benchmark
guide documents the optional command and its limits. Atomic-data generator
regressions, formatting, benchmark package-file listing (including the new
analyzer/tests), and whitespace checks also pass. Logs use
`mass-diagnostics-{tests-final,tests-no-rdkit-final,generator-tests,fmt,package,
final-diff}.log`.

The final analyzer independently verifies all 25 descriptor components in the
current smoke report (`mass-diagnostics-smoke-final-analysis.{json,log}` and
its compressed details). The historical full report also completed analysis:
all 338,837 successful components reconstruct, using 222 distinct element/isotope
terms, with no unexplained masses. Its 2,935 errors and 1,001 not-applicable cases
remain explicitly counted (`mass-diagnostics-historical-analysis.{json,log}`).
That historical analysis preceded the final invalid-record-status guard and
added parser fingerprint; its recorded analyzer hash identifies that version.
The old report is independent evidence, not the current full rerun.

CIAAW's published abridged weights/compositions still identify the 2024 tables,
and the AMDC evaluation page identifies AME2020. These primary pages were
rechecked this turn; cached source hashes and exact regeneration passed. No
runtime mass table, charge convention, formula, or README was changed.
Workspace check/clippy/tests/doctests/docs/MSRV, optional-feature tests/docs and
core/companion packaging were not repeated for this Python/CI/documentation-only
change; preceding gates cover unchanged Rust sources. Reference-adapter and
dashboard suites were not repeated because those sources did not change. Linux/
fuzz checks retain the previously documented Windows exclusion.

Fresh full descriptor validation remains in progress:

- Session 65429 generates pinned RDKit expectations in
  `reference-candidates/mass-diagnostics`, then runs the complete descriptor
  comparison with the retained `mmcif-diagnostics-bench.exe` (unchanged runtime).
  Outputs are `mass-diagnostics-{reference-all,all}.{json,log}`. PubChem reference
  generation is complete: 200,000 outcomes, including 18 explicit reference errors.
- Session 82250 waits for that comparison's complete/error-free run report,
  then executes the maintained analyzer. Outputs are
  `mass-diagnostics-analysis.{json,log}` and its compressed details. The wrapper
  is `target/general-chemistry-review/run_mass_diagnostics.ps1`. A nonzero
  analysis exit may identify retained unexplained records; inspect its
  classifications before drawing a conclusion. The full rerun is pending.

### Completed ring and mmCIF corpus audits

Sessions 95897 and 16937 are terminal. The final selected-ring comparison has
300,832 agreements, one raw path disagreement, no errors and 1,001 NA cases.
The diagnostic audit checked all 301,834 rows, permitting only the additive
diagnostic field; every remaining field and status was identical. For the one
retained difference, both selections cover the same 27 edges and span the same
rank-five cycle space, with five native versus six reference selected cycles.
No cycle was inserted or assertion relaxed. The historical selected-ring audit
also passed: 298,953 unchanged rows and 2,881 former CXSMILES errors now agreeing,
with all expected values preserved and no regressions. Evidence:
`ring-diagnostics-all-audit.json`, `ring-invariants-sssr-audit.{json,log}`.

Sessions 29562 and 28893 are terminal. Full mmCIF comparison has 922 agreements,
79 retained disagreements, no errors and 150,423 NA cases. Its audit checked all
151,424 rows against the historical observations and current baseline, and
proved all raw values/statuses unchanged. Exactly 103 multiline values in those
79 inputs receive diagnostics; each reference is the native text with ASCII
line-end spaces/tabs removed. No unexplained differing values remain in this
corpus. Evidence is `mmcif-diagnostics-after-all.json` and
`mmcif-diagnostics-audit.{json,log}`. Benchmark exit 1 is the retained strict
disagreement, not a run-level failure. Do not restart these completed jobs.

Session 79803 (historical mass analysis) also completed. Continue the fresh mass
run, remaining hydrogen diagnostics, supplied coverage/stereo work and final
fixed-revision validation. All requirements remain in the active goal; DSSP is
still deferred.

### Hydrogen representation diagnostics

The preceding goal turn made progress by adding maintained independent mass
analysis and launching its fresh full run. This turn adds conservative graph-
scoped diagnostics to `chem.hydrogen-transforms`, preserving every raw field,
comparison count and status. It does not change hydrogen chemistry or removal
policy. The implementation is `benchmarks/src/compare/hydrogens.rs`.

A graph receives `hydrogen_representation_only` only when its indexed atom
correspondence, bonds, stereo, groups and every other non-hydrogen-state field
are identical. Every corresponding atom must have equal declared-plus-inferred
hydrogen counts and equal explicit valence after subtracting declared non-graph
hydrogens. Graph H atoms, isotope labels, formal charge, radical occupancy/spin,
coordinates and atom maps therefore stay asserted. Counts must be valid and
cannot claim inferred H while implicit inference is disabled. Overflow and
inconsistent valence are not explained away. This describes the current graph;
it does not make different future inference policies or entire records agree.

Diagnostics list changed atoms and each side's declared/inferred counts and
policy flag. They apply separately to the expanded and collapsed graphs.
Differences in retained explicit graph hydrogens or atom correspondence are
deliberately left unexplained by this classification. No isotope stripping,
stereo-insensitive matching, or input-specific chemical rule was introduced.

Four regressions cover declared/inferred redistribution, a policy difference
with zero H, unchanged raw disagreement/values, rejection of changed totals,
charge/radical/isotope/stereo/bond/atom states, and malformed counts or record
correspondence. All 72 benchmark unit tests and seven integration tests passed,
as did workspace/all-target/all-feature check and warnings-denied clippy,
warnings-denied workspace documentation, Rust 1.89 benchmark/all-target check,
benchmark package-file listing, formatting and whitespace checks. Logs use
`hydrogen-diagnostics-{focused-tests,check,clippy,test,doc,msrv,package,fmt,diff}.log`.

The final executable is retained as `hydrogen-diagnostics-bench.exe`, SHA-256
`d24aff8c6027cc20f923449c882379e63bd660d80c5510eea6aacf93d5f7775b`.
The complete smoke rerun has all 578 raw rows and statuses identical to the
previous run after removing only the added hydrogen diagnostics. Four rows
receive eight graph diagnostics for eight atom observations. The independent
audit rechecks complete non-H graph equality, per-atom totals, policy consistency
and valence arithmetic for every diagnostic, rather than trusting its label.
Evidence is `hydrogen-diagnostics-all-smoke.{json,log}`,
`hydrogen-diagnostics-smoke-audit.{json,log}` and the compressed annotations.

Core/companion tests, doctests, core package verification and optional-feature
docs were not repeated for this benchmark-only change; preceding gates cover
unchanged runtime sources. Python/reference/dashboard/atomic-data/license tests
were not repeated because those sources did not change. Linux/fuzz checks remain
unavailable on this Windows host. README and existing reference archives are
unchanged.

Current full-run continuation handles:

- Session 59876 generates fresh pinned RDKit expectations for all hydrogen
  transformation inputs, then compares the retained pre-diagnostic
  `mmcif-diagnostics-bench.exe`. The separate reference directory is
  `reference-candidates/hydrogen-diagnostics`; outputs are
  `hydrogen-diagnostics-{reference,before}-all.{json,log}`.
- Session 18492 waits for the explicitly complete/error-free baseline, runs the
  final full comparison using `hydrogen-diagnostics-bench.exe`, then invokes
  `audit_hydrogen_diagnostics.py all`. Outputs are
  `hydrogen-diagnostics-after-all.{json,log}` and
  `hydrogen-diagnostics-all-audit.{json,log}`. The audit requires every raw
  observation/status unchanged and independently proves each diagnostic.
  Its wrapper is `run_hydrogen_diagnostics.ps1`. These full runs remain pending.

### Completed fresh molecular-mass validation

Sessions 65429 and 82250 completed. The full descriptor comparison has four
agreements, 300,775 raw mass disagreements, 54 errors and 1,001 NA cases.
The maintained analyzer verified all 341,775 successful component records,
including complete formula equality and all four reconstructed mass values,
with no formula differences or unexplained masses. Its source and executable
fingerprints are retained in `mass-diagnostics-analysis.json`; detailed residuals
and convention contributions are in the adjacent compressed details file.

The historical row audit also passed: 298,953 rows unchanged, 2,881 prior
CXSMILES errors now successfully compared, every expected observation preserved,
and no agreement regressions. Those newly supported inputs remain raw mass
disagreements as required by the unchanged constants/conventions. Evidence is
`mass-diagnostics-run-audit.{json,log}` and `mass-diagnostics-changes.jsonl.gz`.
The 54 native errors comprise 18 strict-valence failures and 36 missing CIAAW
standard weights, with no substitute isotope weights invented
(`mass-diagnostics-errors.json`). No successful case was excluded from the
reconstruction audit. Session 28300 (historical row audit) is terminal.

Sessions 52105 and 77714, the MSRV command and smoke run also finished. Re-poll
59876/18492 before scheduling any rerun. Complete their full audit, remaining
supplied coverage/stereo work and final fixed-revision validation before an
overall completion claim. The full goal remains active; DSSP remains deferred.

### Externally supplied SMARTS coverage

The preceding goal turn made progress by adding hydrogen representation
diagnostics and completing the fresh mass audit. This turn broadens persistent
query coverage with the independent `rdkit-queries` corpus. Existing corpus
locks, inputs and expectations remain unchanged while their jobs run.

The new corpus contains all 518 non-comment query rows from the installed,
pinned RDKit 2026.03.3 distribution: 38 from `FunctionalGroups.txt`, 52 from
`Functional_Group_Hierarchy.txt`, and 428 from `SmartsLib/RLewis_smarts.txt`.
Each source file matched its conda package manifest's per-file SHA-256. The
archive URL/hash, installed paths and source hashes are retained in the source
lock; original bytes, copyright headers and the BSD 3-Clause license are bundled.
No query was selected based on acceptance, rewritten, deduplicated or omitted.
Source IDs retain the table and original line; tests prove complete mechanical
extraction. `PROVENANCE.md` describes the source and scope. Git attributes keep
upstream and packed fixture checksums stable across platforms.

The dataset registry, Git/Cargo inclusion rules and dashboard support the new
small bundled corpus. Its `.smarts` files are applicable only to `query.smarts`,
not molecular algorithms or writers. Source-only completeness/applicability
tests and the staging regression cover this boundary. Dashboard format coverage
now reports SMARTS explicitly. The source-provenance test runs in CI without
RDKit. No runtime query chemistry, README or existing golden was changed.

Independent generation succeeded on every query, with zero reference errors.
The candidate query strings/titles, source membership and all archive checksums
were checked before adopting the new corpus's 25 feature goldens. The other
24 features correctly contain only not-applicable outcomes. All 12,950 rows in
the candidate and adopted comparisons are identical. The initial comparison's
dashboard refresh encountered the temporarily incomplete new catalogue; after
verified adoption the final comparison refreshed it successfully. Both reports
retain their original complete scientific results.

The new corpus has 384 agreements, zero graph-size disagreements and 134 explicit
native unsupported-syntax errors. The errors are 92 recursive/grouped atom
queries; 39 connectivity, ring-bond-count, valence or hybridization primitives;
one exact ring-membership count; one ring-size predicate; and one composite
bond expression. They remain measured errors, not exclusions or empty successes.
These results identify further general grammar work; they do not justify
query-specific substitutions. This benchmark measures syntax and atom/bond
counts, not predicate equivalence, as the guide and provenance now state.

The full `query.smarts` rerun across all six datasets completed with 150,632
agreements, zero disagreements, 134 errors and 1,176 not-applicable observations.
It uses byte-identical independently generated goldens from the earlier complete
stereo-SMARTS run plus the new corpus; copy provenance is retained under
`reference-candidates/query-corpus-full`. An exhaustive audit proved all 151,424
old-dataset observations unchanged. All 578 smoke observations also remain
identical, including the preceding hydrogen diagnostics. Evidence is
`query-corpus-{reference-all,first-all,final-all,query-all,all-smoke}.{json,log}`,
`query-corpus-adoption-audit.json` and `query-corpus-regression-audit.{json,log}`.

The retained executable is `query-corpus-bench.exe`, SHA-256
`45ea988db90d0f0287cca4d445e08ab0d41fa278ae7f6e4f7b02ce645635405d`.
All 73 benchmark unit tests and seven integration tests passed. Workspace/all-
target/all-feature check, warnings-denied clippy, warnings-denied workspace docs,
Rust 1.89 benchmark/all-target check, formatting, package-file listing, provenance
tests, Python dashboard tests, Node live-update checks and whitespace checks
passed. Logs use `query-corpus-{check,clippy,test,doc,msrv,fmt,package,diff}.log`,
`query-corpus-test_query_corpus.py.log`, `query-corpus-test_dashboard.py.log` and
`query-corpus-dashboard-js.log`. Package contents include the new data, upstream
sources/license, provenance and all reference files. Core/companion tests,
doctests, verified core packaging and optional-feature docs were not repeated
for unchanged runtime sources; preceding gates cover them. Unchanged reference-
adapter/atomic-data/license-copy checks and unavailable Linux/fuzz gates retain
their preceding exclusions.

Sessions 96139 and 80436 are terminal; their scientific comparison exit 1 is
the explicit unsupported-query coverage, not a run-level failure. The full query
and smoke commands and regression audit also finished. Sessions 59876 and 18492
were re-polled and remain live for the separate full hydrogen reference/baseline/
final comparisons and diagnostic audit. Do not restart them.

Continue by assessing the newly exposed general SMARTS primitives and remaining
stereo/externally supplied structure coverage, then complete the pending hydrogen
audit and final fixed-revision validation. Recursive queries need a real query
model/matcher design if pursued; no parser-only acceptance should be added.
The full goal remains active and DSSP remains deferred.

### General SMARTS connectivity predicates

This turn adds two syntax-independent atom predicates and their actual matcher
semantics: `TotalConnectivity` (`X`) and `RingBondCount` (`xN`). `X` counts graph
neighbors plus declared and inferred nongraph hydrogens, so materializing a
hydrogen does not change the parent's connectivity. Radical electrons do not
count as connections. `xN` counts incident bonds flagged cyclic by installed ring
membership, independently of any selected ring basis. Bare `X` means exactly
one; bare `x` means any ring membership. Boolean expressions, negation, zero and
bounded numeric counts use the existing expression/parser machinery. Two-letter
element recognition remains intact, including xenon. The matcher requires
valence or ring perception even when these predicates occur under negation; it
never installs perception implicitly.

The definitions were checked against Daylight SMARTS theory
<https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html> and the RDKit Book
<https://www.rdkit.org/docs/RDKit_Book.html>, then tested against the pinned
RDKit 2026.03.3 executable. `v` remains explicitly unsupported: RDKit assigns a
dative bond's valence at its acceptor, whereas the current Kekule valence model
excludes dative orders at both endpoints. Adding a superficially compatible
valence predicate without resolving that model contract would be misleading.
Hybridization, recursive queries and selected-ring count/size predicates also
remain explicit unsupported syntax; no query-specific rewrite was introduced.

Five new focused regressions cover implicit/declared/isotopic/graph hydrogens,
radicals, charge, hydrogen materialization, zero/dative graph neighbors, missing
perception, overflow, and simple/fused/bridged/spiro/cage ring connectivity without
a selected ring basis. The dative regression constructs the supported graph
API directly, because arrow SMILES is not supported. All 27 query tests and all
833 core unit tests pass, along with the workspace integration/companion tests.
The independent focused probe proves identical atom mappings for 43 cases,
including the numeric ring-count matrix, against RDKit.

The full six-dataset `query.smarts` rerun uses the unchanged frozen reference
archives at `reference-candidates/query-corpus-full`. It contains 150,663
agreements, zero disagreements, 103 explicit native unsupported errors and
1,176 not-applicable observations. All 151,942 rows were audited: 31 previous
errors became agreements, eight progressed to their still-unsupported recursive
syntax, and 151,903 rows were unchanged. Every reference observation is identical.
The external query corpus now accepts 415/518 queries; its remaining errors are
100 recursive queries, one exact selected-ring count, one ring-size predicate
and one composite bond expression. These remain measured errors, not exclusions.
Evidence: `query-connectivity-all.{json,log,cases.jsonl.gz}` and
`query-connectivity-run-audit.json`.

Actual matching was independently checked for all 31 newly accepted external
queries against the first 128 supplied SMILES records from each of PubChem and
Enamine. All 7,750 query/target comparisons on 250 connected molecules agree on
every query-ordered atom mapping, including 312 positive comparisons covering
21 distinct queries. Six disconnected records lie outside this connected-molecule
probe and are explicitly retained in its coverage/provenance record. Ten queries
had no positive in this fixed sample, so this probe does not claim exhaustive
positive coverage. Every source line/path/hash, request, expected mapping and
actual result is retained in `query-connectivity-matching-inputs.jsonl`,
`query-connectivity-matching-details.jsonl.gz` and
`query-connectivity-matching-audit.json`. The runner is
`probe_query_connectivity.py`; focused proof is
`query-connectivity-focused-reference.json` and its probe script. No synthetic
molecule was added to the benchmark corpus.

The existing substructure feature also passes all 4,353 applicable cases with
its original 1,000-source limit. An exhaustive comparison proves all 5,354 prior
observations unchanged; the added query-only dataset contributes 518 correct
not-applicable records. Evidence:
`query-connectivity-substructure-verified.{json,log,cases.jsonl.gz}` and
`query-connectivity-substructure-audit.json`. The initial unrestricted before/
after invocations used references generated for that smaller sample, so both
reported 296,480 missing-reference errors beyond it. These are retained as
`query-connectivity-substructure-{before,after}.*`; they are not a successful
full-reference validation. This exposed a separate dashboard defect:
`dashboard.py:load_report` rejects valid coverage-gap reports by requiring all
selected observations to fit inside `golden.cases`, instead of allowing and
showing their explicit missing-reference errors. Address this in the next
benchmark-layer turn with a focused regression; preserve all raw outcomes.

All 578 smoke observations across 25 features remain identical in
`query-connectivity-verified-smoke.json` and `query-connectivity-smoke-audit.json`.
The initial smoke attempt omitted the pinned writer interpreter and stopped
explicitly; the verified rerun supplies `--writer-python` with the pinned RDKit
Python. No golden was regenerated or changed for this implementation.

The retained executable is `query-connectivity-bench.exe`, SHA-256
`c5d392830b86ccd660c6cb3872db54130c028a576ce68f5188252451218ee2e1`.
Workspace/all-target/all-feature check and warnings-denied clippy, complete
workspace tests and doctests, warnings-denied workspace documentation, optional-
feature companion tests/docs, verified core packaging, all three remaining
package-file listings, Rust 1.89 workspace/all-target/all-feature check, format
and whitespace checks passed. Logs use `query-connectivity-*.log`; final fmt and
clippy were rerun after the regression edits. The package-list wrapper initially
had a PowerShell variable-shadowing error; direct commands for each real package
succeeded. Unchanged Python reference/dashboard/atomic-data and Node tests were
not repeated for this runtime-only change. Full companion packaging still needs
the unpublished foundational dependency; package listings were checked instead.
Linux/fuzz gates remain unavailable in this Windows environment as previously
recorded. README and reference adapters are unchanged in this turn.

The pending hydrogen reference/baseline session 59876 is now terminal. Reference
generation completed with 18 explicit reference errors; the complete baseline
contains 285,482 agreements, 15,333 raw disagreements, 18 errors and 1,001
not-applicable observations. Session 18492 remains live for the final hydrogen
comparison and independent diagnostic audit. Do not restart it. Query build,
comparison, focused-test, gate and package/MSRV sessions finished. The separate
partial-reference before/after audit (18246) also finished: all 302,352 raw rows,
including the missing-reference errors, are unchanged except for timing.

Continue with the discovered dashboard coverage-gap issue, then finish the
hydrogen diagnostic audit, remaining coherent stereo/externally supplied
structure coverage and final fixed-revision requirement audit. The full goal
remains active. DSSP remains deferred by agreement.

### Dashboard coverage and incomplete-report validation

The preceding query turn made concrete progress by implementing and validating
`X`/`x`. This turn repairs the dashboard defect exposed by its unrestricted
comparison against a sampled reference archive. `load_report` no longer treats
the number of stored reference records as a ceiling on all selected observations.
Missing-reference and input errors may exceed that archive, and not-applicable
inputs do not require reference lookup. Other applicable cases remain bounded
by stored coverage. Overlapping error-origin counters are capped at the number
of error cases when calculating that necessary bound; this does not classify
all reference errors as missing references or change any observation.

Full coverage continues to mean selected source membership. Completed full
selections without input failures must account for at least all stored records.
A failed input can stand for multiple unreadable records, and an interrupted
run may not have reached every selected input. Neither is discarded for having
fewer observations. Outcome sums, exact-agreement bounds, error-origin counts,
source membership, provenance, finite values and pass-status checks remain.
Missing references stay errors in the denominator; failed or interrupted runs
cannot become successful through dashboard rendering. GUIDE documents these
coverage meanings. No runtime chemistry, reference adapter, golden or raw
comparison was changed.

Five regressions cover missing-reference overrun and rendered-data preservation,
additional not-applicable inputs, an unreadable multi-record input, interrupted
full selection, and rejection of fabricated reference coverage using native,
writer or observation errors. Four fail against the old loader and all pass
with the fix. Existing rejection of an unexplained short completed full selection
still passes. All 26 dashboard tests, all 36 top-level benchmark Python tests,
the Node live-update checks, and all 73 Rust benchmark unit tests plus seven
integration tests pass. Formatting, whitespace and benchmark package-file checks
also pass. Logs: `dashboard-coverage-{before-tests,tests,python-all,node,
rust-tests,fmt,package,diff}.log`.

A fresh end-to-end rerun selects 1,001 PubChem sources against the unchanged
1,000-source substructure references. It produces exactly 2,000 agreements and
two missing-reference errors, remains failed, and successfully refreshes a
separate live dashboard history. Both errored rows have successful native
observations and the explicit missing-golden reference outcome. Evidence:
`dashboard-coverage-rerun.{json,log,cases.jsonl.gz}` and
`dashboard-coverage-runs/`. It uses retained `query-connectivity-bench.exe` because
no compiled source changed; the executable invokes the current dashboard script.
The no-RDKit Biopython environment's Python renders the page using only the
standard library.

Both previously rejected unrestricted reports also load and render successfully,
each preserving all 296,480 errors. The independent audit verifies every
aggregate counter, original report hashes, fresh per-case outcomes and the live
feed's failed status. No input report is rewritten. Evidence:
`audit_dashboard_coverage.py`, `dashboard-coverage-audit.{json,log}` and
`dashboard-coverage-verified.html`.

Workspace check/clippy/core tests/docs/MSRV/verified core packaging were not
repeated for unchanged Rust sources; the preceding query turn's complete gates
cover the exact current compiled sources. Unchanged reference-adapter and
atomic-data tests likewise retain their preceding evidence. Linux/fuzz and full
companion publication exclusions are unchanged. README is unchanged.

Hydrogen session 18492 was re-polled and remains live. Its final comparison has
published a complete report and the wrapper is now running the independent
full-row diagnostic audit. Do not restart it. Dashboard test/rerun sessions are
terminal. Remaining goal work is the hydrogen audit, coherent handling or
explicitly justified deferral of the residual stereo models, broader externally
supplied structure coverage and the fixed-final-revision completion audit.

### Complete hydrogen diagnostic audit and residual classification

The preceding dashboard turn made concrete progress. This turn finishes the
previously live full hydrogen validation and investigates its residuals. Session
18492 is now terminal with exit zero: both comparisons are complete, and the
independent audit verifies all 301,834 rows. Removing only the new additive
annotations leaves every original field, reference observation and status
identical. Final counts remain 285,482 agreements, 15,333 disagreements,
18 errors and 1,001 not-applicable observations.

All 8,386 annotated rows were checked independently, covering 16,668 graph
annotations and 28,465 changed atom-storage states. Each proof checks every
non-hydrogen-storage atom field, all other graph fields, current total H content,
explicit-valence accounting and fixed-H consistency; no annotation asserts a
chemical equivalence based only on matching aggregate counts. Evidence is
`hydrogen-diagnostics-all-audit.{json,log}` and
`hydrogen-diagnostics-all-annotations.jsonl.gz`, produced by
`audit_hydrogen_diagnostics.py all` from the retained before/after executables.

A separate exhaustive residual classification verifies the full report row and
disagreement counts. Exactly 8,384 rows have every raw difference covered by
these representation-only graph proofs. Another 6,949 have residual differences;
two of those also contain a proved annotation on a different graph. Across all
disagreements, 59,124 raw field differences remain asserted, of which 30,273 lie
outside annotated graphs. The leading residual paths concern stereo, followed
by hydrogen storage in graphs whose stereo also differs. No diagnosis is
extended across such a stereo difference. Complete residual observations and
path/example counts are retained in `hydrogen-diagnostics-residuals.jsonl.gz`
and `hydrogen-diagnostics-residual-summary.{json,log}`; the script is
`classify_hydrogen_diagnostics.py` (session 90611 completed successfully).

All 785 round-trip graph-count differences were reproduced independently from
the supplied files with pinned RDKit 2026.03.3. Changing only RDKit's diagnostic
`RemoveHsParameters.removeDefiningBondStereo` from its default false to true
makes every atom count match Kekule. The difference totals 897 graph H atoms.
Every original default reference graph was reproduced exactly, and the native
collapsed and expanded graphs preserve the reference's complete isotope-resolved
composition. This proves a suppression-policy explanation for those counts;
it does not claim full stereo equivalence where other assertions differ. The
option is documented in the official RDKit API:
<https://www.rdkit.org/docs/cppapi/structRDKit_1_1MolOps_1_1RemoveHsParameters.html>.
No production reference option or golden was changed. Proofs include original
source hashes, row/component identities and all counts in
`hydrogen-diagnostics-removal-policy-proof.{json,log}`, from
`prove_hydrogen_removal_policy.py` (session 33611 completed successfully).
Existing native regressions already verify double-bond carrier replacement and
the required parity change during lossless hydrogen collapse.

The residual radical differences are confined to ten atoms in five supplied
SDF rows. Every affected source carbon explicitly declares valence two. RDKit
has zero radical electrons before sanitization and infers two afterward; Kekule
retains zero. The same discrepancy is present after addition and removal. Thus
it is an upstream CTAB/electron-deficit interpretation difference, not evidence
of hydrogen loss. An independent source replay records every atom line, file
checksum, original/component index and before/after occupancy in
`hydrogen-diagnostics-radical-source-proof.{json,log}`, from
`prove_hydrogen_radical_sources.py`. No carbon-specific or valence-two-specific
rule was added. A future repair would need one general, explicitly scoped
source-inference policy and coordinated CTAB reader/writer semantics; ordinary
perception must continue to preserve represented radical state. Keep this
limitation in the final issue report rather than guessing spin or changing a
hydrogen transform to mask it.

No production code, golden, README or adapter changed in this verification turn.
The full rerun uses retained `hydrogen-diagnostics-bench.exe`, SHA-256
`d24aff8c6027cc20f923449c882379e63bd660d80c5510eea6aacf93d5f7775b`.
The feature's four regressions and all applicable compiled-source gates passed
in the implementation and subsequent query turns. They were not repeated for
unchanged source in this audit-only turn; the completed full feature run and
independent proofs were the missing validation. Platform and publication-related
check exclusions remain as previously recorded.

There are no outstanding hydrogen diagnostic sessions: 59876, 18492, 90611 and
33611 are terminal. Remaining work is the residual stereo model assessment,
broader externally supplied structure coverage and fixed-final-revision full
validation/completion audit. The CTAB radical-inference difference above must be
explicitly assessed or reported as deferred, not silently marked repaired. DSSP
remains deferred. The complete goal remains active.

### External structure coverage and first complete comparison

Added the bundled `rdkit-structures` corpus: 50 byte-identical input files from
RDKit `Release_2026_03_3`, commit
`e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`. The selection includes all 47 input
variants of six named medicinal compounds under the upstream atropisomer test
directory and all three root `chebi_*.mol` files. The selection was fixed before
native evaluation; upstream generated `.expected` outputs are excluded, while
all matching `Bad`, 3D and alternative-drawing inputs remain. This is targeted
representation coverage for a small set of compounds, not 50 independent
molecules or a population estimate. ChEBI's R-group inputs measure unsupported
format coverage; their reference descriptor outputs are not evidence of a
physically defined molecular mass for an unspecified substituent.

The exact upstream commit, complete directory tree and BSD 3-Clause license are
bundled. Every source was checked against its Git blob SHA-1, SHA-256 and byte
length. The offline provenance regression verifies the complete selection,
source identities, exact bytes and absence of additional fixtures. Source IDs
retain the original upstream filenames. The new dataset regression verifies
that all 50 inputs are used and query/SMILES/mmCIF/DSSP features remain
not-applicable. Git staging, Cargo package rules, dashboard descriptions, GUIDE
and the CI provenance check include this corpus. README remains unchanged.

Before native comparison, inspection established that all 47 upstream `.sdf`
files actually contain one MOL block ending at `M  END`, with no SDF delimiter
or data fields. Their local extensions are consequently `.mol`; every byte is
unchanged, and original paths and the format-routing explanation remain in the
lock and provenance. This avoids turning an upstream naming convention into
47 strict SDF-delimiter failures. The initial original-suffix reference run and
source lock are retained under the evidence directory but were never adopted.

All 25 features were generated independently with pinned RDKit 2026.03.3 using
the final source lock, yielding 900 successful applicable reference observations
and 350 not-applicable entries, with zero reference errors. Adoption checked
all source/record identities, input/archive hashes, contract and reference-code
hashes, versions and report totals before copying the new corpus archives into
`benchmarks/goldens/rdkit-structures`. Existing goldens were not changed in this
turn. Evidence: `structure-corpus-mol-reference-all.{json,log}`,
`reference-candidates/structure-corpus-mol`,
`adopt_structure_goldens.py` and `structure-corpus-adoption-audit.json`.

The first complete native comparison covers all 1,250 feature/input entries:
413 agreements, 287 disagreements, 200 errors and 350 not-applicable entries.
Exit one is the expected scientific failure result; the report is complete with
no run-level error. The independent audit checks every row against its unchanged
reference observation, fixture bytes and report counters. It also verifies that
all 105 new corpus/provenance/golden files appear in the Cargo package listing.
Evidence: `structure-corpus-final-all.{json,log,cases.jsonl.gz}`,
`audit_structure_corpus.py` and `structure-corpus-audit.json`.

The new findings are retained for the residual stereo assessment:

- Eleven inputs fail molecular interpretation across the 18 applicable features:
  eight report conflicting atropisomeric wedge marks, and three contain
  unsupported `R` or `R1` placeholders. Conflicting assertions must be assessed
  through a general source-stereo policy, not accepted by filename. Placeholder
  atoms require a coherent query/document representation rather than a fabricated
  chemical element.
- Two additional error observations come from one V2000 coordinate-width
  rejection, repeated by the MOL and SDF writers. V3000 preserves its separate
  result. No output coordinate has been clipped to satisfy the reference.
- CIP has 37 agreements, two disagreements and 11 import errors. The two
  disagreements are missing tetrahedral assertions in the BMS-986142 and
  Sotorasib 3D inputs. The axial candidate comparison also exposes represented
  axes absent from native candidate reporting. These provide external evidence
  for the next stereo work; they are not repaired by this coverage addition.
- The five rotatable-bond disagreements all reproduce the existing intentional
  heavy-atom convention. An independent replay first reproduces the complete
  original RDKit observation, then removes graph hydrogens only for diagnosis.
  All five resulting RDKit counts and source-indexed selected bond sets equal
  Kekule exactly. The production comparison and raw disagreements stay intact.
  Evidence: `prove_structure_rotatable_policy.py` and
  `structure-corpus-rotatable-policy-proof.json`.

The all-feature smoke rerun preserves all 578 previous raw case observations
exactly after excluding timing fields: 377 agreements, 73 disagreements, two
errors and 126 not-applicable entries. Evidence:
`structure-corpus-smoke.{json,log,cases.jsonl.gz}` and
`structure-corpus-smoke-audit.json`. Both runs use retained
`structure-corpus-bench.exe`, SHA-256
`6a56b2be5919e5f1d0f5b5bc1ec47b48f78445a029a5dd6e4c69b45da657f55e`.

All 37 top-level benchmark Python tests, the Node dashboard checks, 74 benchmark
unit tests and seven integration tests pass. Workspace all-target/all-feature
check and clippy with warnings denied, workspace docs with warnings denied,
Rust 1.89 workspace/all-target/all-feature check, formatting, benchmark package
listing and whitespace checks pass. Logs use `structure-corpus-{python,fmt,
check,clippy,test,doc,package,msrv,dashboard-js,diff}.log`.

Core and companion runtime sources did not change in this coverage turn, so
their full test/doctest, optional-feature and verified core-package checks were
not repeated; the preceding query turn's complete gates cover them. Reference
adapter and atomic-data tests likewise retain their preceding evidence.
Linux/fuzz checks are unavailable in this Windows environment. Full companion
publication packaging remains dependent on publication of the foundational
crate; no publication was attempted. All new generation, comparison and gate
sessions are terminal. Residual stereo assessment and the fixed-final-revision
full validation/completion audit remain outstanding; the goal stays active.

### Molfile interpretation of unmarked 3D tetrahedral stereo

The external structure addition was concrete progress. This turn implements the
next general defect repair: Molfile interpretation now materializes unmarked
tetrahedral configurations supplied by nondegenerate 3D coordinates. It reuses
the existing chemical eligibility, bounded symmetry and coordinate-orientation
algorithms; there are no compound, atom-index or benchmark-case rules. Ordinary
perception is unchanged and does not materialize stereo. Temporary perception
used during source interpretation is not installed on the published graph.
Inferred omitted hydrogen carriers are resolved by the existing Molfile source
hydrogen policy before publication. Existing explicit unknown assertions block
inference, and this change does not infer unmarked atropisomeric axes.

This is a format-level chemical interpretation, consistent with the documented
RDKit rule that 3D Molfile coordinates carry tetrahedral configuration even in
the absence of wedging:
<https://www.rdkit.org/docs/RDKit_Book.html#interpretation-of-the-2d-3d-flag-and-stereochemistry>.
Existing source assertions and V3000 atom-CFG conflict handling are retained;
this change does not claim every RDKit source-precedence convention is identical.

Import and export are coordinated. A writer uses the same detached inference
to find tetrahedral configurations that the graph leaves unasserted, then emits
unknown bond annotations in temporary writer state. It neither mutates the
source model nor allows a supplied conformation to silently assert new stereo.
Both V2000 and V3000 projections still validate their complete set of marks.

Full tests caught a contract regression in the first implementation: strict
valence preparation rejected a previously readable raw V3000 document. The
final implementation retains such represented graphs and returns a typed
`CoordinateStereoValenceUnsupported` warning, including the original valence
issues, without guessing a configuration. Resource exhaustion, malformed
geometry/state and other inference failures remain errors. Writers reject a
3D projection when stereo preservation cannot be verified. The original raw
metadata/charge/isotope/radical regression still checks all fields and now also
checks the diagnostic and absence of installed perception or fabricated stereo.

Four new regressions cover three and four graph ligands, reflected configurations,
omitted hydrogen resolution, unknown annotations, equivalent ligands, obliquely
planar degeneracy, both Molfile versions, unchanged source models, and V3000
rotation/translation/scales from 1e-100 to 1e100. Reference replay independently
checks all four base/reflected R/S assignments. A separate native writer probe
produces four outputs, which RDKit 2026.03.3 confirms contain no specified
tetrahedral configuration. Removing only their unknown annotation makes RDKit
infer one configuration in each case, proving why that annotation is necessary.
Evidence: `probe_3d_writer.rs`, `probe-3d-writer-v2.exe`, `prove_3d_writer.py`,
the four `molfile-3d-{implicit,explicit}-{v2000,v3000}.mol` files, and
`molfile-3d-writer-reference-proof.{json,log}`. The initial probe incorrectly
expected RDKit to retain a bond direction after interpretation; its corrected
audit checks the consumed `_UnknownStereo` source property and the independent
unmarked control. No benchmark assertion or reference policy was changed.

The final retained executable is `molfile-3d-v2-bench.exe`, SHA-256
`26d66abe6f6b2e1e36f4324bcc982e6b0cbf13e30892259171c5fb9ee8ac6ff1`.
All 25 new-corpus features complete with 420 agreements, 280 disagreements,
200 retained errors and 350 not-applicable entries. Seven disagreements become
agreements; none become errors. Exactly 18 raw observations change across two
source structures and nine features. Both formerly missing 3D CIP assignments
now agree, giving 39 agreements, zero CIP disagreements and 11 import errors.
The complete original reference observations remain asserted. The final results
also match all 1,250 observations from the initial implementation, proving the
valence-policy correction did not alter this corpus. All 578 smoke observations
remain unchanged. Evidence: `molfile-3d-v2-{rdkit-structures,smoke}.{json,log}`,
their case archives, `molfile-3d-structure-transitions.json`,
`molfile-3d-v2-structures-audit.json` and `molfile-3d-v2-smoke-audit.json`.

The final executable's complete all-dataset CIP run has 300,634 agreements,
218 disagreements, 31 errors and 1,519 not-applicable entries, including both
new corpora. The original implementation's full representation and perception
runs also completed; their totals are respectively 259,772 / 41,082 / 29 /
1,519 and 263,667 / 37,187 / 29 / 1,519. These totals are not by themselves
proof that every earlier observation is unchanged. The exhaustive raw audits
and final-executable representation/perception reruns remain in progress as
described below; do not mark their validation complete yet.

All 837 core tests and all workspace/integration tests and doctests pass.
Workspace check and clippy with warnings denied, workspace docs with warnings
denied, optional-feature potentials tests/docs, Rust 1.89 workspace/all-target/
all-feature check, verified core packaging, all companion/benchmark package-file
listings, formatting and whitespace checks pass. Logs use `molfile-3d-v2-`
with the gate name; `run_molfile_3d_v2_gates.ps1` records exact commands. An
initial clippy failure on a test-only unnecessary vector was fixed. The earlier
strict-valence test failure above was fixed before the final complete gate run.
GUIDE and public Rustdoc describe the new source policy; GUIDE's stale CIP depth
description was also corrected to the already-implemented default of 64.

Unchanged Python/reference-adapter/atomic-data and Node checks retain their
preceding evidence; no code on those surfaces changed. Linux/fuzz gates remain
unavailable here. Full companion publication packaging still requires the
foundational crate to be published; exact package listings passed instead.
README is unchanged. This turn does not complete the overall goal.

Live validation handoff: session 2181 runs the final executable's full stereo
comparisons sequentially. CIP and representation are complete; perception is
currently active. Session 46736 runs `audit_molfile_3d.py` against the
completed first-implementation reports. Both handles were re-polled and confirmed
live; do not restart them. Final gate session 92902, final bundled comparisons
60214, original full-comparison session 31105 and core-test session 46959 are
terminal. After session 2181 completes, run `audit_molfile_3d.py molfile-3d-v2`
and record the exhaustive final results here. The final smoke audit already
proves all 578 observations unchanged. The first audit attempt encountered an incomplete
in-progress report and correctly failed its completeness assertion; it did not
publish a partial success. Residual axial/group/heteroatom assessment and the
fixed-final-revision full validation/completion audit remain outstanding.

### Complete 3D audits and conflicting axis source marks

The previous turn made concrete runtime progress. Its outstanding full-data
verification is now complete: sessions 2181, 46736 and 10590 are terminal.
`molfile-3d-v2-audit.{json,log}` proves all 301,834 original observations unchanged
in each of CIP, representation and perception, while retaining the additional
568 query/structure entries per feature. The complete final 3D implementation
reports have the same totals recorded above. Its 18 changed new-corpus rows
are confined to the two expected 3D structures; seven become agreements and
all original reference observations remain unchanged. This supersedes the
preceding live-validation handoff.

This turn makes conflicting axial drawing marks nonfatal in the same way as
ambiguous tetrahedral drawing marks. The source-normalization kernel checks
every mark, records `ConflictingAtropisomericWedgeMarks` in the normalization
warnings, and creates no axis assertion when those marks contradict one another.
Molfile reports map the warning back to the source axis-bond line and retain
the number of marks. Connectivity and unrelated stereo remain intact. The
format document retains the original marks. Explicitly unknown assertions and
consistent redundant marks keep their previous behavior. Writer projection
still rejects any set of marks that produces a normalization warning.

This does not accept conflicting configurations as a specified stereoisomer.
It preserves the readable chemical graph and exposes the failed interpretation
of the drawing. The existing source error variant was replaced by the diagnostic;
no filename, source ID, compound family or atom-number exception was introduced.
The source-selection, reference adapters, contract and goldens are unchanged.
RDKit's documented multiple-wedge consistency rules and the complete independent
new-corpus observations support this distinction:
<https://www.rdkit.org/docs/RDKit_Book.html#defining-atropisomers>.

The regression checks two conflicting wedges and two conflicting hashes in
both V2000 and V3000. The complete interpreted molecule equals its unmarked
counterpart, with no installed perception or invented stereo. The source
document remains unchanged and the diagnostic maps to the actual source bond
line. Existing opposing-mark P/M assignments and all explicitly unknown-mark
combinations remain asserted. An initial test compile failure from accessing
private mapping fields was corrected to use the public accessors.

All 25 new-corpus feature comparisons complete with 502 agreements, 324
disagreements, 74 errors and 350 not-applicable entries. Seven previously
rejected inputs now produce meaningful observations: 82 feature outcomes agree
and 44 disagree, replacing 126 errors. Another 18 observations expose the next
independent failure in the eighth input: V3000 atom CFG on an unsupported
nitrogen center. That annotation remains an explicit error; it is not reinterpreted
as an axial assertion to satisfy a reference. Three R-group inputs remain
unsupported, and the existing V2000 coordinate-width failure remains visible.
The corpus provenance now explicitly explains that R-group descriptor output
cannot establish the physical mass of an unspecified substituent.

The new-corpus audit verifies all 1,250 rows. Only the eight previously failing
source inputs change; all reference observations remain identical. For each of
the seven recovered inputs it independently compares the complete native and
reference stereo arrays, proving that unrelated tetrahedral assertions survive
and no conflicting axis is fabricated. CIP now has 46 agreements, no
disagreements and four import errors. Evidence: `axis-conflicts-rdkit-structures`
report/log/case archive, `audit_axis_conflicts.py`,
`axis-conflicts-structures-audit.json` and `axis-conflicts-changed-rows.jsonl.gz`.

The full all-dataset `stereo.representation` rerun completes with 259,774
agreements, 41,087 disagreements, 22 errors and 1,519 not-applicable entries.
An exhaustive audit of 302,402 observations proves 302,394 byte-identical rows.
Only the eight new-corpus inputs change: two errors become agreements, five
become disagreements, and one retains the independent CFG error. Every original
corpus observation and every reference observation is unchanged. Evidence:
`axis-conflicts-representation-all.{json,log,cases.jsonl.gz}`,
`audit_axis_conflicts_full.py` and `axis-conflicts-full-audit.{json,log}`.
The all-feature smoke run likewise preserves all 578 previous raw observations
(`axis-conflicts-smoke-audit.json`). The retained executable is
`axis-conflicts-bench.exe`, SHA-256
`2f5d7af079bc8a2bbac46992cc5a93a5bb3b6187063cb27ae2e05624364e6c23`.

All 837 core tests, workspace/integration tests and doctests, workspace
check/clippy/docs, optional-feature potentials tests/docs, Rust 1.89 checks,
verified core packaging, companion/benchmark package listings, formatting and
whitespace checks pass. Exact commands are in `run_axis_conflicts_gates.ps1`;
logs use `axis-conflicts-` plus the gate name. The focused regression and final
format/clippy checks were repeated after extending the test to V3000. Unchanged
Python/reference/atomic-data/Node checks retain their earlier evidence. Linux/fuzz
and companion-publication exclusions are unchanged. README was not modified.

The remaining stereo-model scope has now been explicitly assessed. Across all
successfully interpreted new-corpus structures, candidate differences consist
of exactly 33 missing axial candidates. Each axis is already present in the
native represented stereo; there is no loss of its configuration. Evidence:
`axis-model-residuals.json`. The current candidate engine proves eligibility
and orientation-reversing symmetry for tetrahedral/double-bond sites; axes and
non-absolute enhanced groups constrain those proofs conservatively. Cleanup
reports axial and unsupported heteroatom assertions as unclassified and preserves
them. Simply copying source axes into the candidate list would conceal this
capability gap without proving eligibility or stereo dependencies.

Complete axial/group candidate cleanup requires extending the shared chemical
eligibility model and the automorphism constraints to axial orientation and
coupled group inversions, with corresponding permutations and mixed-stereo
regressions. Charged/radical three-coordinate nitrogen additionally needs an
explicit supported electronic-geometry model; an unpaired electron must not be
silently treated as a lone pair. These are deferred under the user's limited-
reengineering scope, rather than patched in the benchmark. Source assertions,
unknown states and all raw disagreements remain visible. GUIDE states this
limitation and distinguishes a represented/geometric axis from demonstrated
configurational stability. The remaining invalid atom-CFG and R-group coverage
limits must also appear in the final issue report.

All sessions from the 3D and axis-conflict work are terminal, including full
comparison 3836, exhaustive audit 12514 and gate 7137. There are no remaining
validation jobs to restart. The next required work is the requirement-by-
requirement completion audit and fixed-final-implementation full validation;
the overall goal remains active. DSSP and the other explicitly documented
unsupported models remain deferred, not claimed as fixed.

### Completion audit: implementation scope

The architecture and the required-work list were rechecked against the final
implementation and the feature-specific evidence above. This audit distinguishes
implemented general rules from explicitly unsupported models.

| Requirement | Implementation and evidence | Remaining scope |
| --- | --- | --- |
| Hydrogen/radical semantics | Total-H SMILES export; radical electron count separated from spin; bracket interpretation/writing and canonical ranking; full hydrogen raw-observation audit | CTAB carbon VAL=2 electron-deficit inference needs coordinated reader/writer rules; stereo-defining graph-H removal policies remain different |
| Rotatable bonds | Fallible descriptor uses compatible installed perception or temporary default perception; H-invariant heavy-atom convention; full feature rerun and independent explicit-H policy proof | Raw RDKit differences under different graph-H conventions remain visible |
| Stereo | Shared tetrahedral/double-bond eligibility, bounded exact symmetry, transactional cleanup, neutral closed-shell N geometry, coordinate-aware source wedges, unmarked 3D import, conflict diagnostics | Axial/enhanced-group candidate symmetry and charged/radical heteroatom geometry need broader model work; source assertions are preserved and reported unclassified |
| CXSMILES | Source-preserving sidecar, explicit base projection, supported radicals/enhanced groups, source mapping and regression coverage | Unsupported chemical extensions are explicit rather than silently discarded |
| Stereo SMARTS | Canonical tetrahedral/double-bond constraints checked under mappings; X/x predicates; external query corpus and independent exact-mapping checks | Recursive/ring-size/ring-count/composite-bond grammar and v/^ semantics remain unsupported |
| CIP resources | Retained and incremental ligand expansion, comparison/order reuse, path-dependent auxiliary handling; measured depth change to 64, node bound remains 100,000 | Two full-corpus node exhaustions and residual chemical convention differences remain explicit |
| Scientific diagnostics | Raw H/mass/ring/mmCIF comparisons retained; guarded explanatory diagnostics; independent mass reconstruction against verified CIAAW/AME/CODATA data | Native mass tables and source whitespace retained; no parity-driven constant or text changes |
| Rings | Independent cycle oracle and atom-permutation regressions; full corpus comparisons | One selected-ring basis difference remains, with equal cycle-space diagnostics |
| Benchmark integrity/coverage | Contract 3 independently retains electron and spin measurements; streaming compressed cases; atomic report publication; coverage/error validation; pinned external query/structure corpora | Historic full archives with old contracts are not silently adopted; unsupported inputs/resource failures remain errors |
| DSSP | Existing assertions and raw discrepancies retained | Interpretation work explicitly deferred |

All implementation gates recorded in the axis-conflict section are complete.
The fixed-revision cross-feature benchmark is the remaining handoff check. It
uses a deterministic 1,000-source selection per dataset (complete smaller
corpora), matching available substructure coverage, and complements the full
feature-specific runs above. It is not a claim of a new all-source run for every
feature. Fresh references, where needed for contract 3, must be generated without
Kekule evaluation; retained references keep their original provenance and hashes.


### Final fixed-revision validation and handoff

The scoped implementation is committed and pushed on `codex/benchmark-parity`
as `11dbce145092d5f21ac2c2ec3bd94442ece4b313`. The final benchmark ran from that
clean revision with the retained `axis-conflicts-bench.exe` (SHA-256
`2f5d7af079bc8a2bbac46992cc5a93a5bb3b6187063cb27ae2e05624364e6c23`).
All later changes are documentation only; runtime, adapters, contract and input
locks remain identical to the benchmark revision. README was not modified.

The run covers all 25 features and all seven datasets, with a deterministic
1,000-source limit: 1,000 PubChem, 1,000 Enamine, 164 PL-REx, all 1,000 PDB,
20 smoke, all 518 RDKit query and all 50 RDKit structure source IDs. This is
3,752 selected source IDs and 175 dataset/feature pairs. It complements the
full-corpus feature reruns above; it does not replace them or claim a new
all-source PubChem/Enamine evaluation of every feature. Formats and repeated
compound representations are separate observations, not independent molecules.

The frozen final reference directory is
`target/general-chemistry-review/final-reference-11dbce14-v3`.
It combines 131 checksum-verified retained pairs with 44 fresh independent
reference generations under contract 3. Four retained DSSP manifests record
an explicit contract-only migration: DSSP has no radical measurements, its
comparison rules and archive bytes are unchanged, and its original generator
provenance remains intact. No reference was changed using Kekule's answers.
The 44 fresh generations were independently checked for archive checksums,
source-lock hashes, selected-source membership, row counts, duplicate identity
and clean-revision generation provenance. Existing archives were not overwritten.

`run_final_validation.py` records the exact generation and comparison commands.
The comparison uses `--feature all --dataset all --limit 1000 --jobs 4`, the
frozen reference directory and the pinned RDKit 2026.03.3 writer interpreter.
`final-reference-consolidation.json` and `final-fresh-reference-audit.json`
record reference provenance and the independent archive audit. These and all
following named outputs are under `target/general-chemistry-review`.

The report `final-11dbce14-all-features.json` is complete with no execution error.
It correctly returns scientific failure (exit 1): 71,326 agreements, 9,525
raw disagreements, 923 error observations and 41,056 not-applicable observations
remain. Native and reference error counts overlap (913 and 244 respectively).
There are no input, writer-validation or observation-layer errors and no
missing-reference errors. An error never contributes to agreement.

| Feature | Agree | Disagree | Error | Not applicable |
| --- | ---: | ---: | ---: | ---: |
| `io.smiles.parse` | 1,718 | 290 | 0 | 1,744 |
| `io.smiles.write` | 1,740 | 0 | 268 | 1,744 |
| `io.smiles.canonical` | 1,951 | 0 | 57 | 1,744 |
| `io.smiles.isomeric` | 1,951 | 0 | 57 | 1,744 |
| `io.mol.parse` | 1,910 | 481 | 4 | 1,521 |
| `io.mol.v2000.write` | 2,213 | 120 | 62 | 1,521 |
| `io.mol.v3000.write` | 2,213 | 178 | 4 | 1,521 |
| `io.sdf.parse` | 1,910 | 481 | 4 | 1,521 |
| `io.sdf.v2000.write` | 2,213 | 120 | 62 | 1,521 |
| `io.mmcif.parse` | 922 | 79 | 0 | 2,751 |
| `algo.rings.fast` | 4,399 | 0 | 4 | 1,519 |
| `algo.rings.sssr` | 4,399 | 0 | 4 | 1,519 |
| `algo.valence.rdkit-like` | 4,399 | 0 | 4 | 1,519 |
| `algo.aromaticity.rdkit-like` | 4,399 | 0 | 4 | 1,519 |
| `algo.canonical-ranking` | 4,399 | 0 | 4 | 1,519 |
| `algo.substructure.vf2` | 4,399 | 0 | 4 | 1,519 |
| `query.smarts` | 2,423 | 0 | 103 | 1,226 |
| `chem.perception.default` | 3,629 | 770 | 4 | 1,519 |
| `chem.hydrogen-transforms` | 4,096 | 303 | 4 | 1,519 |
| `descriptor.molecular` | 0 | 4,397 | 6 | 1,519 |
| `descriptor.rotatable-bonds.rdkit-strict` | 3,764 | 635 | 4 | 1,519 |
| `stereo.representation` | 3,628 | 771 | 4 | 1,519 |
| `stereo.perception` | 3,699 | 700 | 4 | 1,519 |
| `stereo.cip` | 4,395 | 4 | 4 | 1,519 |
| `bio.secondary-structure.dssp` | 557 | 196 | 248 | 2,751 |

The independent final audit (`audit_final_validation.py`,
`final-11dbce14-audit.json` and its log) reads all 122,830 case records, verifies
identity uniqueness, exact selected-source coverage and every summary count,
and rejects any alleged agreement containing a failed native or reference
observation. All 1,828 smoke and structure-corpus rows are identical to the
preceding validated implementation's rows. The finalized case archive is
`final-11dbce14-all-features.cases.jsonl.gz`, SHA-256
`97032a41d1b4d7630b9be9c5fb7cef4b4c4b477e4e85ee94a8a5387fea7a21db`.
No successful comparison was substituted for an unsupported operation.

The final error audit distinguishes intentional format contracts from missing
capabilities. The 268 plain-SMILES export errors reject configured stereo;
canonical/isomeric export is tested separately. Each of those latter features
has 57 enhanced-group export errors. V2000 likewise rejects enhanced groups
(57 cases in each writer feature) and one out-of-range fixed-width coordinate;
V3000 is the supported richer format. No group or configuration is discarded
to make these comparisons pass. Two additional mass errors are missing CIAAW
standard weights for actinium and nobelium, not fabricated numerical values.

The larger work deliberately left for future model extensions remains:

- Axial and enhanced-group candidate symmetry/cleanup, and charged or radical
  heteroatom geometry. Existing represented axes remain preserved; missing
  eligibility proofs are not replaced with a copy of the source assertions.
- Recursive SMARTS, ring-size/ring-count and composite-bond grammar; explicit
  v/^ semantics and enhanced-group SMILES export. Unsupported syntax stays
  explicit, including all 103 query-corpus errors.
- General CTAB electron-deficit radical inference with coordinated reading and
  writing; arbitrary R/R1 chemical placeholders and unsupported atom-CFG sites.
- The two full-corpus CIP node-limit cases and residual CIP conventions already
  documented above. The final sample's four CIP disagreements are retained;
  it does not establish zero disagreement outside its selection.
- DSSP interpretation, as agreed: the full PDB run retains 196 disagreements
  and 248 error observations, including 244 reference failures. Matching
  no-analyzable-residue failures remain errors.

The verified CIAAW/AME/CODATA mass policy, explicit/implicit hydrogen policy,
heavy-atom rotatable-bond convention, valid alternative selected-ring basis and
source-preserving mmCIF whitespace policy remain documented scientific choices.
Their raw comparisons and guarded diagnostics stay visible. They are not
silently normalized into benchmark agreement.

Applicable engineering gates passed: 837 core unit tests and all workspace
integration tests/doctests; workspace all-target/all-feature check and clippy
with warnings denied; warnings-denied workspace docs; potentials tests/docs
without default features; Rust 1.89 all-target/all-feature checking; verified
core packaging; companion and benchmark package-file listings; formatting and
whitespace checks. `run_axis_conflicts_gates.ps1` and the `axis-conflicts-*`
logs retain exact commands. The final V3000 regression extension was separately
rerun along with clippy and formatting. No runtime code changed afterward.

At the clean implementation revision, all 98 Python tests passed: 37 benchmark,
43 RDKit adapter, eight Biopython adapter, eight shared reference-runner and two
atomic-data generator tests. Node dashboard checks and all six license-copy
comparisons passed (`final-python-*.log`, `final-atomic-data.log`,
`final-node.log`, `final-source-state.json`). No engineering check was changed
to accept a chemistry mismatch.

Not run: Linux CI/fuzz compilation and fuzz smoke, because this execution host
is Windows; full companion publication packaging, because it resolves the
unpublished foundational crate through crates.io. The companion package file sets were
checked instead, as prescribed by CI. External trajectory reference runs were
not repeated because trajectory behavior was unchanged. Documentation-only
handoff edits do not require repeating the unchanged runtime suite.

The requirement matrix above is complete within the user's explicitly limited
reengineering scope. Remaining unsupported models and convention differences
are reported rather than claimed fixed. Benchmark session 25308 and audit
session 64990 are terminal; no validation process remains outstanding.
