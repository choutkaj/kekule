# General SMARTS behavioral benchmark

`query.smarts` version 2 retains each input's SMARTS text, source title, atom count
and bond count, and adds complete query-to-target atom mappings. The independent
reference is RDKit **2026.03.3**. Normal benchmark comparisons use stored reference
values and require neither RDKit nor OpenFF.

## Contract

Every source query is matched against all 16 molecules embedded in
[`query-smarts.json`](query-smarts.json). Eight are the existing PubChem smoke
inputs; eight come from the pinned PubChem-100k corpus, selected by the smallest
SHA256(`query-smarts-targets-v2:ID`) values before matching. Source paths, record
indices and checksums are recorded. The general RDKit query corpus contains all
518 rows from three upstream tables, including duplicates at distinct source
locations. Tests verify the complete original extraction and its copied bytes.

Both engines use their default perception/aromaticity and preserve the supplied
hydrogen vertices. Kekule prepares immutable topology targets explicitly, so
disconnected molecules retain a single source atom-index space. Matching does
not change chemistry. MDL aromaticity is measured separately by
`algo.aromaticity.mdl`; this benchmark adds no force-field semantics.

The comparison asserts:

- Every complete mapping, with tuple positions in query-atom order and values in
  source target-atom order. Only the list of tuples is sorted.
- Stereo-aware matches, with automorphisms retained (`uniquify=false`).
- All outer-query tags, sorted by tag number then query-atom index, and every
  corresponding tag projection. Repeated projections are preserved. General
  SMARTS tags may be zero or repeated; force-field-specific rules are out of scope.
- Explicit empty results for non-matches; failure is never converted to an empty
  successful result.

One source query is one scored case. All 16 target observations must agree for
that case to agree. Existing molecular SMILES inputs continue to be interpreted
as queries, preserving their previous parser checks. Missing formats remain
not-applicable; no input is excluded because an engine fails.

Kekule uses complete enumeration, with 100,000 mappings, 1,000,000 search states
and 1,000,000 candidate pairs allowed per query/target. Exceeding a limit fails
the case. RDKit disables both outer and recursive enumeration caps and rejects
results above 100,000 mappings; its worker retains the existing 300-second
process deadline. These are different resource mechanisms, not a speed comparison.
RDKit also tests each disconnected query fragment using RDKit itself: if any
fragment cannot match, the complete query cannot match. Otherwise the original
full query produces the mappings. This necessary-condition precheck avoids
factorial water-fragment searches for hydrated salts with an absent fragment.
It neither assembles mappings from fragments nor changes query predicates.

## Coverage and interpretation

The full run on 2026-09-20 retained 150,766 applicable cases: **150,765 exact
agreements, zero mapping disagreements, and one resource error**. There were
1,176 missing-format cases, counted separately. The report fails overall because
errors do not count as agreement. See the [recorded summary](reports/query-smarts-results.json).

| Corpus | Exact agreements | Errors | Not applicable |
| --- | ---: | ---: | ---: |
| PubChem-100k | 99,999 | 1 | 0 |
| Enamine diversity | 50,240 | 0 | 0 |
| RDKit queries | 518 | 0 | 0 |
| Smoke | 8 | 0 | 12 |
| PL-REX | 0 | 0 | 164 |
| PDB-1000 | 0 | 0 | 1,000 |

PubChem CID **447702**, `C.C.C.C.O`, exceeds the mapping limit against target
`pubchem:174291` in both engines. It is retained as an error, not truncated or
excluded. The [case record](reports/query-smarts-resource-limit-447702.json)
preserves both outcomes. Its old parser-only reference succeeded; adding matching
exposes the resource failure. All other successful old parser observations were
checked against the new values and preserved exactly.

The 518 RDKit patterns agree on all 8,288 query/target pairs. Of these, 136 pairs
are positive, covering 58 distinct query rows and 265 full mappings; the other
8,152 pairs are non-matches. This is a small molecular panel, not evidence of
universal SMARTS conformance. In particular, 460 query rows have no positive
example in the panel. Expanding externally supplied target coverage remains
valuable. PubChem's original selection is also biased by format, size and RDKit
success, as described in its provenance.

Focused regressions additionally check bond-predicate differences with identical
graph counts, nested predicate matching, negation, disconnected queries, both
stereo kinds, ordered tags, benzene automorphisms, enumeration beyond 1,000
matches, explicit exhaustion, missing observation fields and stale contracts.
Toy structures appear only in these regression tests.

## Reference migration and reproduction

The feature-specific contract hash includes the target panel and limits. Old
count-only goldens are rejected; historical dashboard reports remain visible
as stale. Other feature contracts, including MDL, retain their existing hashes.
The previous five datasets' query references are preserved under
`goldens/legacy/query-smarts-v1/`. The migration adds assertions rather than
replacing the original parser observations. Each manifest records the adapter
digest actually used; the later disconnected-fragment optimization preserves
the same observation contract.

```text
cargo benchmark --feature query.smarts --dataset rdkit-queries
cargo benchmark --feature query.smarts --dataset all
```

The second command requires the existing local bulk inputs and golden archives.
Independent regeneration is optional development evidence and must use a fresh
output directory:

```text
cargo benchmark generate --feature query.smarts --dataset rdkit-queries --goldens target/new-query-references --python PATH_TO_PINNED_RDKIT_PYTHON
```

## Implementation checks

### Integration with current main

After merging the general chemistry follow-up (`e4f8054b`), the complete query
comparison was rerun across all seven registered corpora. It retains 150,765
exact agreements, no mapping disagreements and the same CID 447702 resource
error. The new structure corpus contributes 50 additional not-applicable query
cases (1,226 total). See the [post-merge report](reports/query-smarts-merged-results.json).

Contract 3 changes radical/spin and DSSP observations. SMARTS mappings and MDL
atom/bond flags contain none of those changed fields, so their compatibility
metadata was explicitly migrated. All 12 preexisting SMARTS/MDL manifests retain
exactly their payload digest, generator digest, source lock, reference identity
and case count. Each migrated manifest records its prior contract in `origin`.
Historical reports retain their original hashes and remain distinguishable.

Full merged-tree runtime, packaging and fuzz-build checks are recorded in
[core SMARTS integration validation](CORE-SMARTS-MERGE-VALIDATION.md).

### Earlier benchmark-only checks

Passed: benchmark-package formatting, check, clippy with warnings denied,
46 Rust unit tests and four integration tests, documentation with warnings
denied, Rust 1.89 check, and the nonpublished benchmark package's file-list check.
Python checks cover 20 RDKit adapter tests, eight reference-protocol tests and
28 dashboard/provenance/conformance tests; the dashboard JavaScript checks pass.

Published-library package builds, workspace runtime tests and fuzz builds were
not rerun: this change touches the benchmark package, its reference adapter,
fixtures and dashboard, and changes no runtime library or fuzz target. RDKit and
OpenFF remain absent from Rust runtime dependencies. `README.md` is unchanged.
