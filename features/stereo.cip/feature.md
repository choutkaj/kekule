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

The implemented contract assigns isotope-sensitive `R`/`S` for specified
tetrahedral elements and `E`/`Z` for specified double-bond elements when the
local stereo is valid and carrier priorities are unique under the implemented
bounded ranking rules. Unspecified, unknown, or invalid-cleared elements are
skipped. Axis elements, invalid local stereo, unresolved priorities, and
resource-limit exhaustion are reported without assigning lossy descriptors.

## Implementation Notes

The implementation is a descriptor layer over `stereo.representation` and
`stereo.perception`; it does not create stereo elements. Callers should sanitize
or otherwise ensure explicit valence and hydrogen semantics before assignment.

Carrier ranking uses bounded graph expansion with explicit node limits
(defaulting to a 100,000-node budget),
sanitized implicit and explicit hydrogens, cycle guards, and terminal duplicate
nodes for higher-order bonds and ring closures. Ligands are compared as
branch-preserving paired digraphs: immediate child priorities are compared
before recursive branch tie-breaks, and tied child pairs are processed
breadth-first. Per-node priority applies atomic number first, then
RDKit-aligned Rule 1b ring-duplicate priority, then isotope mass. Duplicate
nodes do not carry isotope mass.

The layer validates existing stereo by default and returns structured issues
instead of guessing when the current graph cannot support the stored local
stereo or when the implemented ranking rules cannot distinguish carriers.

## Validation

Unit tests cover tetrahedral descriptors, double-bond descriptors, isotope
priority, Rule 1b duplicate-node ordering, unresolved equivalent ligands,
bounded resource failures, and descriptor invalidation after mutation.

Smoke, PubChem 100, and PubChem 1k validation use externally supplied PubChem
isomeric SMILES fixtures. CIP goldens are generated with RDKit and compare atom
and bond descriptor maps, not bytewise SMILES spelling or internal stereo
element IDs. Bond descriptors are keyed by endpoint atom indexes and descriptor
instead of parser-local bond IDs, because SMILES bond insertion order is not a
portable chemical identity. PubChem 1k is enabled as a parity gate and currently
records non-passing cases for the documented exact-CIP follow-on rules.

## Out Of Scope

Full exact machine-oriented CIP parity remains out of scope for this version:
Rule 3/4/5 handling, pseudoasymmetric `r`/`s`, full duplicate-node Rule 1b
parity across all digraph and mancude cases, mancude and fractional atomic
numbers, axial `M`/`P`, non-tetrahedral geometries, enhanced stereo relation
semantics, broad RDKit parity, isomeric SMILES emission, and stereo enumeration.

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
