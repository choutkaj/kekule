# Stored reference coverage

The current benchmark has 156 dataset/feature pairs: 26 features across six
datasets. Every locked ID is accounted for, including missing formats and
reference failures. These are reference-availability counts, not Kekule scores.

| Dataset | Source IDs | Reference values | Reference errors | Missing-format cases |
| --- | ---: | ---: | ---: | ---: |
| pubchem-100k | 100,000 | 3,799,720 | 280 | 200,000 |
| enamine-diversity | 50,240 | 1,909,120 | 0 | 100,480 |
| pl-rex | 164 | 6,232 | 0 | 1,148 |
| pdb-1000 | 1,000 | 1,756 | 244 | 24,000 |
| smoke | 20 | 477 | 0 | 127 |
| rdkit-queries | 518 | 518 | 0 | 12,950 |

An ID contributes to several features and sometimes through both SDF and SMILES.
Missing formats are reported separately from applicable cases; they never count
as agreement. Reference failures remain errors. Inspect each report's coverage
before interpreting its agreement counts.

The 52 smoke and RDKit-query archives are tracked. The other 104 current
archives stay local; all 156 current manifests remain tracked so their expected
bytes and provenance can be verified. Archived version-1 SMARTS manifests and
their smoke payload are retained separately. No bulk golden archive is published
by the benchmark workflow.

The goldens use RDKit 2026.03.3,
Biopython 1.87, and mkdssp 4.6.1. Adjacent schema-2 manifests pin the compressed
file, source lock, comparison contract, reference version and case count. See
[GUIDE.md](GUIDE.md) for the comparison and provenance rules.

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

## Scope

All 50,240 Enamine SDF records, all 1,000 PDB entries, and the existing source
membership of the other datasets remain included. PubChem's historical
preselection by format, size and RDKit success remains a sampling limitation;
these files cannot demonstrate correctness on excluded chemistry.

The 18 fixed substructure searches provide limited behavioral coverage. The
`query.smarts` feature now additionally compares full mappings for every supplied
query on 16 pinned targets; see [SMARTS validation](QUERY-SMARTS-VALIDATION.md).
Its version-1 parser observations are archived rather than silently reused.
The separate MDL feature is described in [MDL validation](MDL-VALIDATION.md).
The mmCIF feature checks decoded
syntax, while DSSP exercises only part of biomolecular interpretation. Agreement
with one reference is evidence for the measured fields on these inputs, not a
proof of general correctness.
