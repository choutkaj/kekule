# SMILES text and identity, 2026-09-28

This change moves native canonical SMILES ordering toward RDKit conventions and
separates exact text agreement from the existing molecular-identity comparison.
It does not promise complete RDKit compatibility or add a SMARTS writer.

The subsequent [four-rule follow-up](SMILES-RDKIT-FOLLOWUP.md) raises full-corpus
exact agreement from the 85.38% recorded here to 96.29%. This document retains
the initial implementation's measurements and provenance.

## Contracts

The dashboard shows **SMILES** for `io.smiles.text.write`,
`io.smiles.text.isomeric` and `io.smiles.text.canonical`, and **SMILES identity**
for the existing `io.smiles.write`, `io.smiles.isomeric` and
`io.smiles.canonical` features. Existing feature IDs, identity observations,
goldens and historical result meanings are preserved.

The text features observe `records[].smiles` directly from native emission, with
titles in a separate asserted field. They never invoke the independent reader,
normalize the emitted string, or remove CX fields. `CCO` versus `OCC` is a text
disagreement even when the identity check agrees. Ordinary stereo rejections
remain errors. Text comparison runs using stored references without RDKit;
identity comparison still requires the pinned RDKit reader.

[smiles-text.json](smiles-text.json) defines the independent RDKit writer options
and contributes to feature-specific contract hashes. It also contributes to the
reference-adapter fingerprint. Canonical text uses canonical isomeric CXSMILES;
ordinary text uses noncanonical aromatic CXSMILES; isomeric text uses
noncanonical Kekule CXSMILES. All use `CX_ALL`. The underlying molecular source
preparation is the same as the existing benchmark. This explicitly defines the
comparison instead of treating RDKit's default flags as interchangeable with
every Kekule writer mode.

## Native ordering

The writer now builds its initial partition using atom map, degree, atomic
number, isotope, hydrogen count and formal charge priorities, with represented
radical state and emitted atom spelling retained as additional distinctions.
Neighborhood refinement compares descending bond/class lists. Complete bounded
individualization/refinement still resolves remaining ties using the full graph
and stereo certificate. The public molecular symmetry-class algorithm is
unchanged.

Emission starts at the least-ranked atom instead of trying every root and
choosing by branch count and lexicographic string. Traversal visits acyclic
substituents before ring continuations and emits the last child as the main
continuation. Ring membership is computed once per planning/emission pass. The
existing resource bounds and failure behavior are retained conservatively.

The full-corpus run also exposed eight substituted conjugated chains whose
directional-bond carrier selection depended on source bond numbering. The
directional constraint solver now seeds components in canonical endpoint order
and couples redundant marks at shared alkene endpoints. Their slash/backslash
phases remain compatible during emission and independent of numbering. All
eight PubChem cases have permutation and read/write fixed-point regressions;
the full-corpus checks retain those assertions too.

These rules were checked against RDKit release `Release_2026_03_3`:
[atom comparison](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/new_canon.h),
[SMILES writing](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/SmilesParse/SmilesWrite.cpp),
and [DFS traversal](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/Canon.cpp).
This is a native implementation of selected ordering conventions, not a port of
all RDKit ranking, stereo cleanup or CX policies. RDKit remains a benchmark-only
dependency.

## Measured improvement

Compared the same deterministic 1,000-ID selections from each of PubChem and
Enamine and all eight supplied smoke SMILES records, using RDKit **2026.03.3**.
The baseline is clean revision `f036b0b643d93616239e9861b39f3ee13707e291`.
Baseline raw strings and reference observations are retained under
`target/string-parity-audit/`. The new measurements are under
`target/smiles-rdkit-rules/`. Every sampled new canonical-text expectation was
checked against the corresponding baseline RDKit identity string; all 2,008
expectations are unchanged apart from the observation field name.

| Corpus | Cases | Exact strings before | Exact strings after |
| --- | ---: | ---: | ---: |
| PubChem | 1,000 | 127 | 846 |
| Enamine | 1,000 | 33 | 882 |
| Smoke | 8 | 3 | 7 |
| Total | 2,008 | 163 (8.12%) | 1,735 (86.40%) |

There are no writer errors in these canonical selections. Of the paired cases,
1,576 changed from different to equal, 159 stayed equal, 269 stayed different,
and four changed from equal to different. No fuzzy similarity measure is used.

The selections are samples of the bulk corpora, not full-corpus agreement rates.
The original PubChem selection is biased by its upstream selection criteria.
Smoke is not necessarily independent of the bulk corpora. Other corpus IDs have
no applicable SMILES source format for this feature and remain not applicable;
MOL/SDF files were not silently converted into additional SMILES cases.

Focused regressions now assert the RDKit strings for terminal-root selection,
carboxyl branches, furan, aspirin, naphthalene and a chain with unlike terminal
halogens, using alternate source traversals and read/write fixed points. The
existing randomized atom/bond permutation regressions were extended with these
aromatic and acyclic cases. Stereo, isotope, hydrogen projection, symmetry and
resource-limit regressions remain asserted.

The final full-corpus canonical text run (`accepted-full-text-canonical.json`)
also checks numbering invariance and read/write fixed points:

| Corpus | Cases | Exact strings | Different strings | Errors |
| --- | ---: | ---: | ---: | ---: |
| PubChem | 100,000 | 83,614 | 16,377 | 9 |
| Enamine | 50,240 | 44,663 | 5,577 | 0 |
| Smoke | 8 | 7 | 1 | 0 |
| Total | 150,248 | 128,284 (85.38%) | 21,955 | 9 |

All nine errors are valence rejections by both engines; matching failures do
not count as agreement. There are no remaining native-only errors in this run.
The eight conjugated-chain cases above now pass the harness checks and were
also independently reread with RDKit to confirm molecular identity. No
full-corpus baseline rate was measured; the before/after comparison is the
paired 2,008-case sample above.

For completeness, the same 2,008-case selection gives 121 exact matches and
1,887 differences in noncanonical isomeric mode; ordinary mode gives 339 exact
matches, 1,401 differences and 268 errors. All ordinary-mode errors are its
explicit rejection of stereo-bearing inputs. These modes preserve their
distinct writer behavior and are not canonical-equivalence measurements.

Independent full text references were generated for all three modes and all
seven datasets: 21 new manifests and nine bundled small-corpus payloads. Bulk
payloads remain local. Their compressed hashes, source locks, complete source-ID
coverage, contract hashes, reference versions and generator fingerprints were
verified before adoption. Each mode retains the same nine RDKit valence
failures. Existing identity references were not replaced.

## Remaining differences

Identity agrees on 2,007/2,008 sampled canonical cases. This must not be described
as an absence of individual identity regressions: baseline failures
`Z6781268208` and `Z9466031432` now agree, while **Z7371431020** now disagrees.
All three contain explicit enhanced groups plus the legacy CX `r` flag. The
ordering change selects a different representative for Kekule's relative stereo
group. RDKit treats the otherwise ungrouped tetrahedron differently, so its
reread identifies a different absolute stereo representative. Native relative
group semantics are preserved; this known interpretation incompatibility is
neither hidden nor reclassified as agreement.

Remaining text differences include refined atom ordering, ring-closure labels
and their order, stereo-group representatives, and CX radical field spelling.
For example, smoke caffeine still differs:

```text
Kekule: Cn1cnc2c1c(=O)n(C)c(=O)n2C
RDKit:  Cn1c(=O)c2c(ncn2C)n(C)c1=O
```

Canonical output spelling is intentionally changed by this work; callers using
it as a persisted key need to account for the version of the canonicalizer.

## Reproduction and validation

```text
cargo benchmark generate --feature io.smiles.text.canonical --dataset all --python PATH --goldens NEW_DIRECTORY
cargo benchmark --feature io.smiles.text.canonical --dataset all --limit 1000 --jobs 4 --goldens NEW_DIRECTORY
cargo benchmark --feature io.smiles.canonical --dataset all --limit 1000 --jobs 4 --goldens IDENTITY_REFERENCES --writer-python PATH
```

The text feature is expected to exit unsuccessfully while exact differences
remain; identity likewise retains the documented CX disagreement. Generation
also reports failure when retained reference errors are present, while preserving
their outcomes and finalized manifests. Generation does not run Kekule. New
references must be published under the text feature
names, not substituted for identity observations.

Validation passed: workspace/all-target/all-feature checking and Clippy with
warnings denied; workspace tests including doctests (1,475 passed, three
pre-existing ignored tests); warnings-denied workspace documentation; Rust 1.89
workspace/all-target/all-feature checking; fuzz-target compilation; potentials
tests and documentation without default features; foundational package build
and verification; companion and benchmark package-file listings (including all
21 new manifests and exactly nine new bundled payloads); formatting and
diff whitespace checks. Python RDKit reference tests (48), reference runner
tests (8), dashboard tests (28), and JavaScript dashboard tests passed.

Linux CI and fuzz execution were not run on this Windows host. Full companion
package verification remains a release-stage check requiring the foundational
crate on crates.io; exact package file sets were checked as in CI. Scientific
agreement numbers above describe the explicitly declared samples or full corpora.
