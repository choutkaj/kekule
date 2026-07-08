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

Carrier ranking uses bounded graph expansion with explicit node limits,
sanitized implicit and explicit hydrogens, cycle guards, and terminal duplicate
nodes for higher-order bonds and ring closures. At each expansion depth,
atomic-number spheres are compared first, then RDKit-aligned Rule 1b
ring-duplicate priority, then isotope-mass spheres. Duplicate nodes do not carry
isotope mass.

The layer validates existing stereo by default and returns structured issues
instead of guessing when the current graph cannot support the stored local
stereo or when the implemented ranking rules cannot distinguish carriers.

## Validation

Unit tests cover tetrahedral descriptors, double-bond descriptors, isotope
priority, Rule 1b duplicate-node ordering, unresolved equivalent ligands,
bounded resource failures, and descriptor invalidation after mutation.

Smoke validation uses externally supplied PubChem SMILES fixtures and compares
semantic JSON containing the CIP assignment report plus stereo elements with
derived descriptors.

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
