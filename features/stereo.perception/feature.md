# Stereochemistry Perception

## Summary

Validate graph-local stereo elements, detect candidate stereochemical units,
and assign local stereo elements from supported source marks and coordinates.

## Behavior/API

- Exposes `perception::stereo::{StereoPerceptionOptions,
  StereoPerceptionReport, StereoCandidate, StereoPerceptionIssue,
  validate_stereo, validate_stereo_with_options, perceive_stereo,
  perceive_stereo_with_options}`.
- `validate_stereo` is read-only. It reports invalid local stereo elements,
  potential tetrahedral atom candidates, potential double-bond candidates, and
  source-mark assembly diagnostics without mutating the graph.
- `perceive_stereo` is mutating. It runs the same checks, assembles supported
  source marks and supported coordinate geometry into first-class stereo
  elements, records created stereo element IDs, and marks stereo perception
  fresh.
- Candidate detection uses current graph and hydrogen state. Sanitization or
  valence perception should run first when implicit hydrogens matter.
- Validation accepts stored implicit lone-pair tetrahedral carriers for
  supported no-implicit heteroatom centers, but candidate detection does not
  broadly perceive lone-pair stereocenters.
- Tetrahedral candidates are local geometry candidates only; no exact ligand
  ranking, symmetry pruning, or CIP descriptor assignment is performed.
- Stored axis elements are validated structurally: the axis bond must exist and
  the local reference carriers must be atom carriers adjacent to opposite
  endpoints of the axis bond. Broad axis candidate detection remains out of
  scope, but the Molfile wedge subset below can assemble supported stored
  axes.
- Double-bond candidates include atom and implicit-hydrogen carriers. Paired
  directional `/` and `\` source marks are normalized relative to each alkene
  endpoint and can create specified double-bond stereo elements when compatible
  marked single bonds are available on both ends. A substituted alkene endpoint
  may carry two redundant directional marks when they cover both atom carriers
  with opposite endpoint-normalized directions. Aromatic-focus bonds,
  double bonds between aromatic atoms, double bonds in rings smaller than
  eight atoms, and endocyclic double bonds with a non-carbon endpoint are
  skipped until those unsupported stereo cases have explicit
  format/perception semantics.
- Molfile wedge up/down marks can create specified tetrahedral stereo elements.
  When conformer coordinates are available, wedge up/down orientation is
  assembled from those coordinates with the marked bond direction treated as
  the local out-of-plane sense. If the local center has exactly one implicit
  hydrogen carrier, perception synthesizes a virtual implicit-H coordinate from
  the three explicit substituent vectors after applying the wedge z offset
  rather than falling back to a fixed up/down orientation.
  Molfile wedge/either marks can create explicit unknown tetrahedral stereo
  elements. In both cases the marked bond's first endpoint is treated as the
  local stereo center and the marked carrier is placed first in carrier order.
- Molfile wedge up/down marks that cannot support tetrahedral stereo can
  conservatively create specified atropisomeric axis elements when the marked
  bond is adjacent to a single-bond axis whose endpoints are SP2-like, have
  total degree two or three, and have one or two explicit reference bonds each.
  Ring membership, aromatic flags, aromatic-order bonds, and local double bonds
  provide the current SP2-like evidence. When a Molfile source mark can support
  both non-ring and ring-internal axis candidates, non-ring candidates are
  preferred; ring-internal macrocyclic axes are allowed when no non-ring
  candidate remains. Redundant wedge marks that define the same axis are
  consumed as the same local axis evidence instead of becoming hard tetrahedral
  ambiguity. The stored local reference carriers are the lowest-ID explicit
  carriers on the two axis endpoints, matching RDKit's documented internal
  convention. Coordinates determine the stored axis handedness; for 2D
  Molfiles, the wedge/hash mark supplies a virtual out-of-plane sign used to
  lift the reference carriers consistently.
- Coordinate-derived assignment uses the first conformer conservatively. It
  assigns tetrahedral stereo only when all four carriers are explicit atoms
  with nondegenerate 3D coordinates, and assigns double-bond stereo only when
  each side has exactly one explicit atom carrier with nondegenerate 2D or 3D
  geometry. It does not infer coordinates for implicit hydrogens.

## Implementation Notes

This feature identifies candidate tetrahedral atoms and double bonds, validates
existing local stereo elements against current topology and hydrogen semantics,
including supported implicit lone-pair tetrahedral carriers and structurally
valid stored axis carriers, assembles SMILES-style paired directional bond
marks with endpoint-relative normalization, and assembles supported Molfile
tetrahedral wedge/either source marks plus supported Molfile atropisomeric
wedge subsets. Coordinate-bearing Molfile tetrahedral
wedges with one implicit hydrogen use a virtual carrier coordinate derived
from explicit substituent geometry. Molfile atrop axes store lowest-neighbor
endpoint references and compute local handedness from coordinates so equivalent
wedge placements around the same axis collapse to the same local axis state.
Ring-internal single bonds adjacent to the same wedged endpoint are ignored
when a non-ring candidate is available so alternate wedged substituent
placement does not create ambiguous axes; otherwise ring-internal macrocyclic
axis candidates are accepted. Exocyclic axes may have either two ring endpoints
or one ring endpoint plus one locally SP2 endpoint, and macrocyclic
ring-internal axes follow the same endpoint and reference-carrier structural
checks, matching the official RDKit atropisomer fixtures covered by smoke CIP
parity validation.
Coordinate-derived assignment is local and conservative; exact CIP descriptors
belong to `stereo.cip`. Small-ring double-bond exclusion uses a bounded
shortest-path check around the candidate bond so direct perception and
sanitizer-driven perception follow the same RDKit-like stereogenic-bond
boundary.

Small-molecule perception should run as an explicit staged workflow and may be
run from the small-molecule sanitizer when `SanitizeOptions::perceive_stereo`
is enabled. The sanitizer uses the source-mark assembly subset and leaves
coordinate-derived assignment to explicit stereo perception calls. It should not
run over whole `MacroMolecule` structures by default.

## Validation

- Unit tests cover read-only validation, candidate detection after sanitization,
  directional double-bond assembly, unsupported double-bond exclusions including
  the small-ring alkene boundary, Molfile wedge/either assembly, unsupported
  source-mark diagnostics, structural validation for stored axis elements,
  Molfile atropisomeric axis assembly from official RDKit fixtures, including
  alternate wedged substituent placement around the same exocyclic axis,
  redundant wedge marks around the same atrop axis, one-ring-endpoint SP2
  exocyclic axes, ring-internal macrocyclic axes,
  virtual implicit-H Molfile tetrahedral orientation,
  coordinate-derived tetrahedral and double-bond assignment, sanitizer
  integration, transactional rollback, and preservation of explicit unknown
  versus absent stereo.
- Smoke, PubChem 100, PubChem 1k, PubChem 100k, Enamine diversity, and PL-REX
  validation record semantic perception JSON for externally pinned isomeric
  SMILES fixtures covering absent stereo, stored tetrahedral stereo, and
  directional double-bond source-mark assembly. PL-REX adds coordinate-bearing
  ligand SDF packs for Molfile wedge/either, conformer-aware tetrahedral
  orientation, and supported atropisomeric source-mark assembly coverage. The
  broader PubChem, Enamine, and PL-REX tiers are implementation-golden semantic
  regression gates for perception stability, while exact RDKit descriptor
  parity belongs to `stereo.cip`. Broad semantic validation records sanitize
  failures per record when unrelated unsupported valence chemistry prevents
  stereo perception.

## Out Of Scope

Exact CIP descriptors, broad axis candidate perception beyond the supported
Molfile wedge subsets, coordinate-only axis assignment, CXSMILES
atropisomeric syntax, isomeric SMILES writing, enhanced stereo serialization,
implicit-hydrogen coordinate reconstruction, stereo enumeration, and reaction
stereo transfer.

## Revision Notes

- v1: Feature contract reserved for stereo candidate detection and local
  validation.
- v2: Add public stereo perception API, local element validation, conservative
  tetrahedral/double-bond candidate detection, paired directional source-mark
  assembly, unit coverage, and smoke semantic validation.
- v3: Assemble supported Molfile wedge up/down source marks into specified
  tetrahedral elements and wedge/either source marks into explicit unknown
  tetrahedral elements.
- v4: Integrate stereo perception into the explicit small-molecule sanitization
  pipeline with opt-out options, reporting, freshness-state handling, and
  transactional rollback on stereo issues.
- v5: Add conservative coordinate-derived local assignment for explicit-atom
  tetrahedral centers and double bonds using the first conformer.
- v6: Normalize SMILES directional bond marks relative to alkene endpoints and
  accept redundant two-mark substituted endpoints when the marks cover both atom
  carriers with opposite normalized directions.
- v7: Validate supported implicit lone-pair tetrahedral carriers and skip
  unsupported aromatic or endocyclic hetero double-bond stereo candidates before
  source-mark assembly.
- v8: Add PubChem 100 and PubChem 1k semantic regression requirements for
  stereo perception over externally supplied isomeric SMILES.
- v9: Add PubChem 100k and Enamine diversity semantic regression requirements
  and preserve unsupported-sanitization records as per-record validation output
  instead of aborting broad perception fixtures.
- v10: Exclude double-bond stereo candidates in rings smaller than eight atoms
  using the RDKit-like stereogenic-bond boundary while preserving cyclooctene
  and larger ring alkene candidates.
- v11: Derive Molfile wedge up/down tetrahedral orientation from conformer
  coordinates when present so coordinate-bearing V2000 records preserve
  RDKit-like local stereo sense.
- v12: Validate stored axis elements structurally instead of treating every
  axis element as unsupported, enabling the CIP layer to assign descriptors for
  explicitly stored axes.
- v13: Assemble a conservative Molfile atropisomeric wedge subset into stored
  axis elements using ring-atom axis eligibility and 2D opposite-side carrier
  selection.
- v14: Restrict Molfile wedge-derived atrop axis candidates to non-ring axis
  bonds, matching the exocyclic subset covered by official RDKit RP-6306
  variants and avoiding ambiguous ring-internal candidates.
- v15: Use virtual implicit-H geometry for coordinate-bearing Molfile
  tetrahedral wedges and store Molfile atrop axes with RDKit-style
  lowest-neighbor endpoint references plus coordinate-derived handedness.
- v16: Broaden Molfile atrop axis assembly to RDKit-like SP2 endpoint
  eligibility, consume redundant same-axis wedge marks, and cover official
  BMS/ZM atropisomer regressions.
- v17: Allow ring-internal macrocyclic Molfile atrop axes when no non-ring
  candidate is available from the same source mark, with official RDKit
  macrocycle regressions.
- v18: Add PL-REX ligand SDF packs to the perception validation contract for
  coordinate-bearing Molfile stereo and source-mark assembly regression
  coverage.
