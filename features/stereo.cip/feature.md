# CIP Stereochemistry

## Summary

Assign bounded CIP descriptors as a derived cache over validated local stereo
elements.

## Behavior/API

`molecules::perception::stereo` exposes:

- `assign_cip_descriptors`
- `assign_cip_descriptors_with_options`
- `CipAssignmentOptions`
- `CipAssignmentReport`
- `CipAssignment`
- `CipSkipped`
- `CipSkippedReason`
- `CipAssignmentIssue`

Assignment mutates `StereoElement.descriptor` on the underlying `Molecule`.
Descriptors are derived cache, not graph truth: assignment clears existing
descriptors first, and topology or stereo-invalidating mutations clear them
again.

The implemented contract assigns isotope- and E/Z-descriptor-sensitive `R`/`S`
for specified tetrahedral elements and `E`/`Z` for specified double-bond
elements when the local stereo is valid and carrier priorities are unique under
the implemented bounded ranking rules. Unspecified, unknown, or
invalid-cleared elements are skipped. Axis elements, invalid local stereo,
unresolved priorities, and resource-limit exhaustion are reported without
assigning lossy descriptors.

## Implementation Notes

The implementation is a descriptor layer over `stereo.representation` and
`stereo.perception`; it does not create stereo elements. Callers should sanitize
or otherwise ensure explicit valence and hydrogen semantics before assignment.

Carrier ranking uses bounded graph expansion with explicit node limits
(defaulting to a 100,000-node budget), sanitized implicit and explicit
hydrogens, implicit lone-pair carriers, cycle guards, and terminal duplicate
nodes for higher-order bonds and ring closures. Ligands are compared as
branch-preserving paired digraphs using RDKit-like recursive sequence-rule
ordering: Rule 1a atomic number comparison is exhausted through the digraph
before Rule 1b duplicate-node priority is considered; Rule 2 isotope mass is
considered only after both earlier rules remain tied; and Rule 3 orders
embedded `Z` descriptors before embedded `E` descriptors before unlabeled
double bonds. Duplicate nodes do not carry isotope mass, and duplicate nodes
for higher-order bonds back to the original stereocenter are suppressed.

Assignment is descriptor-aware and iterative. Descriptors that are unique under
constitutional rules are assigned first, then previously unresolved elements
are retried so embedded E/Z labels can break otherwise equal ligand ties without
storing descriptor flags on atoms or bonds.

The layer validates existing stereo by default and returns structured issues
instead of guessing when the current graph cannot support the stored local
stereo or when the implemented ranking rules cannot distinguish carriers.

## Validation

Unit tests cover tetrahedral descriptors, double-bond descriptors, recursive
Rule 1a/1b/2 ordering, Rule 3 embedded E/Z ordering, isotope priority, Rule 1b
duplicate-node ordering, implicit lone-pair carriers, unsupported double-bond
stereo exclusions, unresolved equivalent ligands, bounded resource failures,
and descriptor invalidation after mutation.

Smoke, PubChem 100, PubChem 1k, and PubChem 100k validation use externally
supplied PubChem isomeric SMILES fixtures. CIP goldens are generated with RDKit
and compare atom and bond descriptor maps, not bytewise SMILES spelling or
internal stereo element IDs. Validation records include molecules where RDKit or
the implementation assigns at least one CIP descriptor; no-descriptor molecules
are filtered out so broad CIP validation is not dominated by unrelated parser or
sanitizer coverage for structures with no stereochemical labels. Bond
descriptors are keyed by endpoint atom indexes and descriptor instead of
parser-local bond IDs, because SMILES bond insertion order is not a portable
chemical identity. Molecules validation maps removable plain explicit
hydrogens out of descriptor records to match RDKit default atom indexing.
PubChem 100k is enabled as a broad RDKit parity gate for current
descriptor-bearing coverage.

## Out Of Scope

Full exact machine-oriented CIP coverage remains out of scope for this version:
Rule 4/5 handling, pseudoasymmetric `r`/`s`, `seqCis`/`seqTrans`, mancude and
fractional atomic numbers, axial `M`/`P`, non-tetrahedral geometries, enhanced
stereo relation semantics, parity beyond the current descriptor-bearing
validation corpora, isomeric SMILES emission, and stereo enumeration.

## Revision Notes

- v1: Feature contract reserved.
- v2: Reframe CIP as a derived-cache layer over representation and perception,
  with deterministic ranking and sanitized chemistry as dependencies.
- v3: Add bounded descriptor assignment for validated local tetrahedral and
  double-bond stereo elements, with explicit skip and issue reporting.
- v4: Correct isotope priority so isotope mass refines equal current-sphere
  atoms before deeper substituent atoms are considered.
- v5: Add RDKit-aligned terminal duplicate nodes, Rule 1b ring-duplicate
  priority before isotope priority, and zero isotope mass for duplicate nodes.
- v6: Switch CIP validation to RDKit-backed descriptor maps and require smoke,
  PubChem 100, and PubChem 1k parity corpora.
- v7: Use branch-preserving paired breadth-first ligand comparison, normalize
  SMILES directional double-bond marks through stereo perception, and key CIP
  bond validation by endpoint atoms rather than parser-local bond IDs. Raise
  the default CIP node budget to cover larger fused-ring PubChem parity cases
  while preserving explicit resource-limit failures.
- v8: Compare descriptor-bearing records in CIP validation and add PubChem 100k
  as a broad RDKit parity gate. The large gate exposes remaining exact-CIP
  ligand-ordering mismatches after unrelated no-descriptor parse/sanitize noise
  is filtered out.
- v9: Apply recursive RDKit-like Rule 1a, then Rule 1b, then Rule 2 comparison;
  add implicit lone-pair carrier support for supported heteroatom centers;
  suppress root-adjacent multiple-bond duplicates; skip unsupported aromatic
  and endocyclic hetero double-bond stereo; and align validation output with
  RDKit default explicit-hydrogen indexing.
- v10: Add descriptor-aware iterative assignment and RDKit-like Rule 3 ordering
  for embedded `Z` versus `E` double-bond descriptors.
