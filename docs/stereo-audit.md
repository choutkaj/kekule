# Stereo and CIP audit

Audit date: 2026-09-11. Current reference: RDKit **2026.03.6**.

This audit covers represented stereo, source interpretation, coordinate inference,
hydrogen transforms, CIP assignment, perception installation, and interchange.
The implementation remains pure Rust. See the [operation support contract](stereo-support.md)
and [reproduction instructions](stereo-validation.md).

## Method

The parity target is RDKit's accurate `AssignCIPLabels` implementation, using
ordered sequence rules and rooted digraphs. It is separate from RDKit's legacy
approximate ranking and from stereo candidate perception.

- [IUPAC Blue Book P-9](https://iupac.qmul.ac.uk/BlueBook/P9.html).
- [Hanson et al., 2018](https://doi.org/10.1021/acs.jcim.8b00324), including the
  proposed Rule 6 refinement used by RDKit.
- [Pinned RDKit implementation](https://github.com/rdkit/rdkit/tree/Release_2026_03_6/Code/GraphMol/CIPLabeler).
- [Accurate CIP API](https://www.rdkit.org/docs/source/rdkit.Chem.rdCIPLabeler.html).
- [Published CIP Validation Suite](https://cipvalidationsuite.github.io/ValidationSuite/),
  revision `6b9f9db46dadc6749da8234b05164e1e0fb413b9`.

The durable runner records parsing, sanitization, and prepared atom tags separately.
Its `sanitized` mode reproduces the original broad-corpus procedure;
`assertions` excludes `SANITIZE_CLEANUPCHIRALITY`. Both retain explicit hydrogen
vertices, clear previous labels, and invoke the accurate labeler. Neither procedure
can undo a configuration change made while parsing SMILES.

## Changes

CIP ranking now:

- Applies sequence rules cumulatively while sorting descendants.
- Uses isotope masses at natural-weight boundaries, including the VS176 oxygen
  case, with an attributed compact mass table.
- Computes complete mancude component averages before installing fractional atomic
  numbers, removing an atom-order dependency.
- Assigns auxiliary tetrahedral, double-bond, and axial descriptors to occurrences
  in the rooted digraph, including the appropriate ring duplicates.
- Handles balanced enantiomorphic descriptor pairs and fixed reference ordering.
- Uses Rule 6 ligand partitions and reference parity instead of atom-ID, ring-path,
  or remote-descriptor counting shortcuts; propagates the pseudoasymmetric
  comparison metadata corrected in RDKit 2026.03.6.
- Expands actual atropisomer digraphs and cancels two pseudoasymmetric endpoint
  inversions correctly.
- Expands constitutional comparisons progressively, reports depth/node exhaustion,
  and installs results transactionally.
- Ranks explicit double-bond assertions without reapplying candidate-perception
  exclusions for small rings, aromatic Kekule bonds, or heteroatoms.

Global assignment rounds, auxiliary descriptor seeding, and associated ordering
heuristics were removed. Incomplete ligand comparisons cannot become chemical
ties or pass to lower sequence rules.

Representation and source handling now:

- Reject malformed carriers, duplicate focuses, invalid axis endpoints, and
  descriptor families incompatible with their represented geometry.
- Preserve valid stereo IDs, groups, and previous perception on failed edits or
  assignment.
- Remap hydrogen carriers before deleting each explicit hydrogen.
- Handle implicit-H and lone-pair tetrahedra, fully substituted alkenes, and
  coordinate-scale independence.
- Decode consistent redundant wedges, retain explicit unknown marks, and reject
  contradictory atropisomer wedges.
- Interpret V3000 atom CFG 1/2/3, including CTfile's hydrogen-last carrier order.
- Read and write enhanced atropisomer groups using RDKit's endpoint-`ATOMS`
  convention, with ambiguity and cross-component checks.
- Write specified Model E/Z only when the rounded emitted coordinates encode it.
  Unasserted alkenes that would acquire E/Z from the drawing receive either/crossed
  syntax; rereading may introduce an explicit unknown element.

These are deliberate validation/API changes. New error variants may require
updates to exhaustive matches. Canonical SMILES behavior is unchanged.

## Adjudicated reference differences

**VS132, Troger's base.** The published structure and Kekule give S/S. Default
RDKit sanitization removes both nitrogen tags; preserving those tags instead gives
R/S from the published SMILES. The remaining discrepancy occurs in RDKit's SMILES
chirality normalization, before ranking. Independent Rule 1a ordering and signed
volumes from the published 3D SDF give S/S. Transporting those local configurations
onto the original RDKit SMILES graph also gives S/S with unchanged ligand ranks.

The published SDF, exact source hashes, local-carrier regression, and
`vs132_reproducer.py` preserve this evidence. This is not an omitted CIP label
or an unexplained accepted mismatch.

**V3000 explicit-H CFG.** CTfile requires hydrogen to be the highest-numbered
carrier even when its atom row occurs earlier. RDKit 2026.03.6's optional
`AssignAtomChiralTagsFromMolParity` helper ignores this exception. Kekule follows
[CTfile Appendix A](https://www.wincept.eu/toxlab/pdf/ctfile.pdf), with a regression
for hydrogen in each row position. Ordinary RDKit Molfile reading retains atom CFG
as parity metadata without applying that helper; normal output interoperability is
verified through wedge syntax.

**SMILES bracket policy.** RDKit brackets atoms adjacent to metals even when the
input used ordinary organic-subset syntax. That changes explicit/implicit H
declarations and explicit valence without changing the molecular graph. These
fields remain asserted against a documented source-preservation projection;
RDKit's full raw emission and decoded fields remain separate reference evidence.
When chemical normalization creates a formal charge, brackets are mandatory and
the expected declaration is the minimum legal SMILES projection, with original
source fields retained.

The eight Rule 6 reordering differences observed in RDKit 2026.03.5 are resolved
in 2026.03.6. All five base structures and twenty reordered encodings now agree.
The original .5 outputs remain historical evidence, not current discrepancies.

## External validation

Comparisons retain complete descriptor maps, label absence, lowercase labels,
graph counts, disconnected records, explicit failure statuses, and input hashes.
Identical input SMILES are deduplicated within each differential run.

| Check | Result |
| --- | --- |
| Unchanged RDKit 2026.03.3 `stereo.cip` snapshots | 178/178 fixtures |
| Enamine stereo-bearing SMILES, RDKit 2026.03.6 | 8,366/8,366 exact |
| PubChem 100k stereo-bearing SMILES, RDKit 2026.03.6 | 12,536/12,536 exact |
| Combined published/RDKit regression encodings, assertion mode | 1,517/1,571 exact; every other case classified below |
| Published suite, assertion mode, distinct SMILES | 280/298 exact; same classified boundaries |
| Molfile cross-tool validation | 134 requests / 265 outputs; complete chemical graphs, labels, and groups match |
| Schema-v2 representation | 166/166 fixtures; 151,506 records |
| Schema-v2 perception | 167/167 fixtures; 151,507 records |
| Schema-v2 isomeric output, including stereo-free controls | 14/14 fixtures; 1,103 records |
| Independent SMILES declaration-projection checks | 11/11 external structures |

Neither broad corpus lost atom tags during reference preparation. The 54 nonmatches
in the combined adversarial run comprise 51 unsupported inputs (15 allene/cumulene
and 36 wildcard/query encodings), VS132, and two unavailable references.
In assertion mode, ten of the unsupported axial inputs also cause RDKit's
`Received an invalid Atom Descriptor` error. RDKit exceeds one million recursive
iterations on VS009 and VS226. These are not successful comparisons; the reports
retain both implementations' results, including one concurrent debug-probe timeout
on VS009. A separate run completed that input.

The unsupported geometry cases are: allene VS078/079/120/144/231/232/243/287;
longer axial cumulenes VS141/166; extended cumulene E/Z VS063/118/135/154/164.
Helical, planar, and general coordination stereo also require new geometry types.
These feature additions and the canonical SMILES stereo/isotope issue remain
separate work.

The benchmark schema migration preserves the original gzip snapshots and hashes
under `benchmarks/reference/stereo/schema-v1`. The independent migration reads
published source bytes, prior assertions, RDKit's graph, and defined carrier/geometry
transformations; it never uses current Rust output as an expected result. Source
marks now belong to document evidence, and SDF perception uses its real coordinates.
All-record isomeric coverage restores previously filtered stereo-free controls.
All 89 previously overlapping isomeric assertions remain unchanged.
Independent verification confirmed all 347 archived snapshots exactly match their
original Git bytes, all replacement and source hashes match the migration manifest,
and all 178 CIP snapshots and their input hashes remain unchanged.

## Verification

All 13 applicable Rust validation commands passed: formatting; workspace and
Rust 1.89 checks; fuzz-binary check; Clippy with warnings denied; workspace and
documentation tests; documentation builds; no-default-feature potential tests and
documentation; foundational package build/verification; and the two companion
package-content checks prescribed by CI. Documentation used
`RUSTDOCFLAGS=-D warnings`. Package commands used `--allow-dirty` for these
uncommitted changes. Full companion package dependency verification is not enabled
by CI until the foundational crate is published; its exact file-list checks passed.

The workspace passed **1,156 tests across 53 suites**, with no failures or ignored
tests. The three external adapters passed their four unit tests. All ten registered
Linux fuzz targets passed **256 seeded runs each**, with AddressSanitizer, a
4,096-byte input limit, and the checked-in seed corpora. Generated cases were written
under `target`, preserving the original seeds.

The independent reference checkers and migration passed all 18 Python regression
tests, including failure-record retention, chemical graph loss, mandatory bracket
projection, preserved historical fields, and read-only migration behavior.

All 242 Rust source files remained unchanged during final verification. All six
packaged license copies match their originals, and `git diff --check` passed.
The final compiled CIP adapter again matched all 20,902 broad-corpus structures.
Results and command logs are retained under
`target/stereo-production/`; original audit evidence remains under
`target/stereo-audit/`. The durable tools, source fixtures, migration provenance,
regressions, and contracts are checked-in source files.

No external reference is a runtime dependency or an added CI release gate.
`README.md` is unchanged.
