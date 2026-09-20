# Stored reference coverage

The registry now contains 182 dataset/feature pairs: 26 features across seven
datasets. The 78 bundled archives cover smoke, the 518-query RDKit corpus and
the 50-file RDKit structure corpus. The latter supplies several representations
of a small set of compounds; it is not a random sample of 50 independent molecules.
Bulk inputs, full reference archives, reports and diagnostic outputs remain local.

The general chemistry follow-up uses observation contract 3, which asserts
radical electron occupancy and explicitly supplied spin independently. Compatible
references and their provenance are recorded in `GENERAL_CHEMISTRY_REVIEW.md`.
Historical full-dataset archives below are not silently migrated or regenerated
by comparison. An archive with an old contract is rejected; generate a fresh
independent reference into a new directory or select the audited compatible
archive explicitly with `--goldens` as described in `GUIDE.md`.

## Original five-dataset snapshot

The following counts describe the original five-dataset reference inventory,
before the general chemistry follow-up. They are reference-availability counts,
not current Kekule scores. Every locked ID was accounted for, including missing
formats and reference failures.

| Dataset | Source IDs | Reference values | Reference errors | Missing-format cases |
| --- | ---: | ---: | ---: | ---: |
| pubchem-100k | 100,000 | 3,599,739 | 261 | 200,000 |
| enamine-diversity | 50,240 | 1,808,640 | 0 | 100,480 |
| pl-rex | 164 | 5,904 | 0 | 1,148 |
| pdb-1000 | 1,000 | 1,756 | 244 | 23,000 |
| smoke | 20 | 452 | 0 | 126 |

An ID contributes to several features and sometimes through both SDF and SMILES.
Missing formats are reported separately from applicable cases; they never count
as agreement. Reference failures remain errors. Inspect each report's coverage
before interpreting its agreement counts.

At that snapshot, the local compressed goldens occupied 1,845,149,606 bytes.
The 25 smoke archives (141,639 bytes) were the only tracked payloads; the
other 100 archives stayed local. All 125 small manifests remain tracked to
preserve their expected bytes and provenance. The two later bundled corpora
add 50 archives and manifests. Full golden archives are not published by the
benchmark workflow.

The goldens use RDKit 2026.03.3,
Biopython 1.87, and mkdssp 4.6.1. Adjacent schema-2 manifests pin the compressed
file, source lock, comparison contract, reference version and case count. See
[GUIDE.md](GUIDE.md) for the comparison and provenance rules.

## Fast-ring reference operation correction, 2026-09-17

The `algo.rings.fast` reference now calls RDKit `FastFindRings` instead of
computing a symmetrized selected-ring set. All five dataset archives were
independently regenerated with RDKit 2026.03.3 and compared in full: every
one of the 301,834 rows retains exactly the same observation, source identity,
input digest, reference version and missing-format status. No asserted value
changed. The regenerated pairs carry the actual generator fingerprint;
their predecessors remain backed up locally. Compression sizes changed,
but decompressed row equality was verified independently of compression.

## Valence reference preparation correction, 2026-09-17

The `algo.valence.rdkit-like` reference now invokes RDKit `Cleanup` before
its non-strict property-cache update, matching the represented-chemistry
normalization stage of the native input path. All five dataset archives were
independently regenerated with RDKit 2026.03.3 and audited in full. Exactly
316 PubChem rows change: only oxohalogen/oxygen formal charges and explicit
valences change under RDKit's cleanup rule. Total charge, explicit and implicit
hydrogen counts, every other atom field and all source provenance remain
unchanged. The other 301,518 rows are identical, including missing-format cases.
There are no new reference errors. The old archives and manifests are backed
up locally; the promoted pairs record the actual generator fingerprint.
No comparison field or tolerance was removed or weakened.

## Independent corrections, 2026-09-16

The former adapter silently read 2,881 supplied Enamine CXSMILES extensions as
titles. Those cases were independently reevaluated with RDKit across the 17
molecular features that consume them. Their corrected titles and molecular
observations retain the full extension meaning. Unaffected golden lines were
copied verbatim. All 2,881 source records remain in the benchmark, even where
Kekule rejects an extension.

DSSP references now include explicit beta-partner identities and independently
calculated CA-C-N-CA omega angles. Regeneration used original supplied mmCIF
bytes. Every previously asserted value in every successful DSSP case was checked
against the preceding golden and preserved. The same 244 PDB reference failures
remain, with their actual causes preserved instead of a generic empty-output
message. Of these inputs, 230 contain only nucleic-acid polymers; the remaining
14 previously failed through no analyzable residues, missing DSSP output, or
missing input metadata.

The current values outside these corrections were retained. Older observations
lack the original generator fingerprint; their manifests say so explicitly.
Correction provenance identifies the scope and adapter digest. This historical
limitation cannot be repaired by inventing provenance or substituting Kekule
results. Fresh generation records its actual adapter fingerprint.

The unused historical corpus goldens, indexes, projection fixture and support
scripts were removed from the benchmark tree. They remain recoverable from Git
history and are not part of the current distribution. Current goldens retain all
asserted fields; recompression of updated files saved about 15 MB.

## SDF metadata correction, 2026-09-17

The RDKit adapter's default property enumeration omitted source names beginning
with an underscore. It now retrieves values through RDKit in source-header order,
preserving private source fields without adding internal RDKit properties.

All 50,240 Enamine observations for each of `io.mol.parse`, `io.sdf.parse`
and `io.sdf.v2000.write`
were independently regenerated with RDKit 2026.03.3 in separate directories.
Complete old/new comparisons verified that exactly 2,881 records per feature
gained their supplied `_CXSMILES_Data` field. Every previously asserted value,
chemical observation, case identity and source digest was unchanged. Corrected
archives and their generation manifests replace only these three dataset/feature
pairs; the original pairs remain backed up locally. The adapter fingerprints
are recorded for these complete generations.

The SDF writer archive was independently generated both before and after the
record/field scanner correction. The two compressed archives are byte-identical,
confirming that the scanner fix changes no supplied reference observations.
The promoted manifest records the final scanner's actual fingerprint.

The SDF review also regenerated all other datasets through RDKit's SDF reader,
including standalone MOL inputs. Every observation outside the Enamine metadata
additions was unchanged, including the same nine PubChem reference errors. This
full audit covers 151,588 rows, including not-applicable IDs. The fresh audit
archives remain local; unchanged historical archives retain their provenance.

The `stereo.representation` review on 2026-09-19 independently regenerated all
100,480 Enamine observations, including both source formats. The exhaustive
comparison again found only the 2,881 missing `_CXSMILES_Data` fields in SDF
records; all chemistry, SMILES observations and source identities were unchanged.
Its corrected archive and generation manifest were promoted after this audit,
with the original pair retained locally. The metadata regression now explicitly
covers stereo representation as well.

The `stereo.perception` review repeated this independent 100,480-observation
audit. Only the same 2,881 private SDF fields changed; candidate sets, stereo,
all other chemistry and source identities were unchanged. Its corrected archive
and manifest were promoted with the original pair retained locally, and this
feature is now included in the private-property regression.

Other feature archives have not been regenerated for this correction. Their
SDF-property coverage will be audited in the corresponding feature reviews.

## Scope

The 2026-09-19 DSSP review updates the shared comparison contract to record the
single-precision omega allowance and nullable polymer sequence labels. All 125
manifests were rebound to that contract only after verifying their compressed
archive and source-lock hashes. No payload, case membership, reference version,
generator fingerprint or generation provenance changed in this rebinding.
`PARITY_REVIEW.md` documents the numerical derivation and schema regression.
The exhaustive local audit is `target/benchmark-parity/dssp-contract-rebinding.json`.

All 50,240 Enamine SDF records, all 1,000 PDB entries, and the existing source
membership of the other datasets remain included. PubChem's historical
preselection by format, size and RDKit success remains a sampling limitation;
these files cannot demonstrate correctness on excluded chemistry.

The 18 fixed SMARTS searches provide limited behavioral coverage. The separate
SMARTS feature checks acceptance and graph size. The mmCIF feature checks decoded
syntax, while DSSP exercises only part of biomolecular interpretation. Agreement
with one reference is evidence for the measured fields on these inputs, not a
proof of general correctness.

## Core SMARTS and MDL additions

The query feature now asserts full stereo-aware mappings against 16 pinned
external targets. Its previous parser-only references remain archived, and the
new contract is feature-specific. MDL aromaticity is an independent feature.
See [SMARTS validation](QUERY-SMARTS-VALIDATION.md) and
[MDL validation](MDL-VALIDATION.md) for recorded observations and limitations.
These reports predate the merge of the general chemistry follow-up and remain
historical evidence with their original contract and implementation hashes.
