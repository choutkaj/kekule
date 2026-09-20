# Standalone MDL aromaticity benchmark

`algo.aromaticity.mdl` independently compares Kekule's explicit
`AromaticityModel::Mdl` with RDKit **2026.03.3**
`Chem.AromaticityModel.AROMATICITY_MDL`. Both adapters perform their ordinary
input interpretation and perception/sanitization, then explicitly install MDL.
The reference calls `Chem.Kekulize(mol, clearAromaticFlags=True)` first so that
neither atom nor bond flags from the default aromaticity model survive.

Every atom flag and every endpoint-qualified bond flag is compared exactly.
There is no tolerance, atom-set reduction, or omission of difficult molecules.
Disconnected source records retain all their components. Hydrogen representation
is preserved. The measured workflow includes parsing, ordinary perception, MDL
application, and observation construction, rather than isolated algorithm timing.

The feature appears in `cargo benchmark --list`, `--feature all`, stored-reference
validation, reports, and the benchmark dashboard. It has its own reference files;
the existing `algo.aromaticity.rdkit-like` references remain unchanged.

```text
cargo benchmark --feature algo.aromaticity.mdl --dataset smoke
cargo benchmark --feature algo.aromaticity.mdl --dataset all
```

## Reference provenance

New reference observations were generated through the normal reference-only
runner from every locked input in the five existing corpora. Each adjacent
`goldens/<dataset>/algo.aromaticity.mdl.jsonl.meta.json` records source-lock,
compressed-file, contract, reference-version, and adapter-code hashes. Smoke
goldens are checked in; large payloads remain local under the existing policy.
Reference errors remain stored outcomes, not exclusions or successful cases.
PDB-only inputs have no applicable molecular format and remain explicitly
not applicable.

The smoke corpus contains 11 cases whose reference MDL assignment differs from
the default aromaticity feature. All 25 applicable smoke cases agree with Kekule,
so the smoke validation exercises the model distinction as well as unchanged
aromatic and aliphatic environments.

## Measured comparisons

The 2026-09-20 run uses complete source selections, without a `--limit` filter.
Counts are source/format cases, not distinct chemical structures.

| Corpus | Exact agreements | Aromaticity mismatches | Errors | Not applicable |
| --- | ---: | ---: | ---: | ---: |
| Smoke | 25 | 0 | 0 | 1 |
| PubChem 100k | 199,981 | 1 | 18 | 0 |
| Enamine diversity | 97,599 | 0 | 2,881 | 0 |
| PL-REX | 328 | 0 | 0 | 0 |
| PDB | 0 | 0 | 0 | 1,000 |
| Total | 297,933 | 1 | 2,899 | 1,001 |

The single disagreement is PubChem CID **181201**, SMILES
`C1=CC(=CC=[C]1)N`, at `data/packs/pack_037.smi`, zero-based record 786.
RDKit marks the six carbon atoms and six ring bonds aromatic; Kekule marks them
nonaromatic. The same difference occurs under the default aromaticity model, so
it is not isolated to MDL selection. RDKit represents atom 5 as a carbon radical
with one radical electron. This observation identifies a case for a subsequent
chemistry investigation; it does not establish the root cause. The discrepancy
is retained in [the complete failing case](reports/mdl-discrepancy-181201.json).
The runtime implementation was not changed as part of adding this benchmark.

The 18 PubChem errors are simultaneous failures in both engines: RDKit reports
`normalization_or_perception_error`, and Kekule reports valence-perception issues
(16 cases with one issue and two with two issues). They are errors, not
agreements. The PubChem comparison therefore exits nonzero too.

All 2,881 Enamine errors are Kekule's explicit CXSMILES rejection, before MDL
perception. For example, `Z9078819835` in `data/packs/pack_001.smi`, record 8,
reports `CXSMILES extensions are not supported by Kekule; refusing to discard
them`. RDKit evaluated those inputs successfully. They remain failed cases in
the denominator; the Enamine comparison exits nonzero. No molecular input or
asserted field was removed to obtain agreement.

Local summaries and complete observations are retained under `runs/mdl-*.json`
and the corresponding `.cases.jsonl` files. A portable copy of all five original
run summaries is in [MDL results](reports/mdl-results.json). The PDB reference archive records
all 1,000 not-applicable entries; its standalone generation command exits
nonzero because no applicable case was measured, rather than because RDKit
failed. The existing dashboard displays not-applicable counts separately.

## Regression checks

Focused Rust and RDKit-adapter regressions distinguish pyrrole from biphenyl:
MDL clears pyrrole's atom and bond flags while retaining biphenyl's aromatic
rings and its nonaromatic connecting bond. The Rust regression proves that
changing either an atom flag or a bond flag fails comparison, and that omitting
bond observations fails schema validation. The reference regression also checks
that model selection leaves the source molecule unchanged.

Passed repository checks: `cargo fmt --all -- --check`, workspace all-target and
all-feature `cargo check` and clippy with warnings denied, Rust 1.89 checking of
all benchmark targets, all 46 benchmark unit/integration tests, and benchmark
documentation with warnings denied. All Cargo checks used `--locked`.
The 17 RDKit-adapter tests, eight reference-protocol tests, 21 Python dashboard
tests, JavaScript dashboard checks, and final dashboard refresh also passed.
Runtime crate tests, fuzz targets, and package checks are not rerun for this
benchmark-only change: no runtime source, public API, package metadata, or fuzz
target changed. README files are unchanged.
