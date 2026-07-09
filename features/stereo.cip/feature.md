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

The implemented contract assigns `R`/`S`/`r`/`s` for specified tetrahedral
elements, `E`/`Z` for specified double-bond elements, and `M`/`P` for
specified stored axis elements when the local stereo is valid and carrier
priorities are unique under the implemented bounded ranking rules, including
natural-vs-indicated isotope priority, E/Z descriptor priority, descriptor
class and pair priority, and RDKit-like mancude fractional priority for bond
duplicate nodes. Rule helpers also understand stored
`seqCis`/`seqTrans`
descriptor values for RDKit-like sequence-rule ordering and use a
parity-stable symmetric S4-style Rule 6 retry for fully equivalent
tetrahedral carrier sets. Unspecified, unknown, or invalid-cleared elements
are skipped. Invalid local stereo, unresolved priorities, and resource-limit
exhaustion are reported without assigning lossy descriptors. Double-bond
elements in rings smaller than eight atoms do not receive E/Z descriptors,
matching the RDKit-style stereogenic-bond boundary.

## Implementation Notes

The implementation is a descriptor layer over `stereo.representation` and
`stereo.perception`; it does not create stereo elements. Callers should sanitize
or otherwise ensure explicit valence and hydrogen semantics before assignment.

Carrier ranking uses bounded graph expansion with explicit node limits
(defaulting to a 100,000-node budget), sanitized implicit and explicit
hydrogens, implicit lone-pair carriers, cycle guards, and terminal duplicate
nodes for higher-order bonds and ring closures. Ligands are compared as
branch-preserving paired digraphs using RDKit-like recursive sequence-rule
ordering: Rule 1a compares integer atomic numbers for normal nodes and
mancude neighbor-averaged atomic-number fractions for affected bond-duplicate
nodes, and is exhausted through the digraph before Rule 1b duplicate-node
priority is considered; Rule 2 compares indicated isotope masses against
natural atomic weights only after both earlier rules remain tied; and Rule 3
orders embedded `Z` descriptors before embedded `E` descriptors before
unlabeled double bonds. Rule 4a orders
uppercase sequence descriptors (`R`/`S`/`M`/`P`) and `seqTrans`/`seqCis` before
pseudo or geometric descriptors (`r`/`s`/`E`/`Z`) before unlabeled nodes; Rule
4b derives a reference descriptor from the first
descriptor-bearing equivalent ligand level and compares like versus unlike
descriptor families relative to that reference; Rule 4c orders
pseudoasymmetric `r` before `s`; and Rule 5 uses descriptor-pair lists over
the `R`/`M`/`seqCis` versus `S`/`P`/`seqTrans` descriptor families so like
descriptor pairs outrank unlike pairs. Rule 6 is a contextual tetrahedral retry
that selects a reference atom from unresolved equivalent carrier partitions and
gives priority to ligand nodes that point back to that reference. Duplicate
nodes do not carry isotope mass, duplicate nodes for higher-order bonds back to
the original stereocenter are suppressed, bond duplicates carry the parent
atom's mancude fraction when it is fractional, and negative atoms with
fractional mancude atomic numbers receive one terminal duplicate child
following RDKit's digraph expansion rule. For fully equivalent tetrahedral
carrier sets, Rule 6 retries all atom carriers that can produce a complete
ranking and accepts the result only when every successful reference choice
preserves the same carrier permutation parity. Natural atomic weights are
stored as integer-scaled ranking values matching RDKit reference validation,
while duplicate nodes keep zero Rule 2 mass.

Assignment is descriptor-aware and iterative. Each assignment round evaluates
pending elements against a descriptor snapshot and applies successful
assignments only after the round completes, so a center cannot commit using
only a prefix of newly assigned neighboring descriptors. Descriptors that are
unique under constitutional rules are assigned first, then previously
unresolved elements are retried so embedded E/Z labels and descriptor-pair
labels can break otherwise equal ligand ties without storing descriptor flags
on atoms or bonds. After ordinary descriptor propagation stalls, deferred
tetrahedral assignment builds a root-centered auxiliary occurrence graph for
the current center. Tetrahedral occurrences in that graph are labeled by
virtually re-rooting their local ligand view within the same occurrence graph,
so ring-duplicate paths and repeated visits keep their local identity instead
of collapsing to a molecule-level atom descriptor. These auxiliary labels are
scoped to the current comparison context, cached to avoid recursive explosion,
and skipped for the center currently being assigned, matching the RDKit-style
separation between primary descriptors and local digraph auxiliary
descriptors. Precomputed auxiliary lookups fall back to already assigned
primary tetrahedral descriptors when no local occurrence label exists, allowing
absolute neighboring centers to resolve later coupled ties. Deferred Rule 6
still handles selected coupled ring-center ties in the same retry phase, and
when a deferred batch contains any absolute tetrahedral assignments those are
committed before pseudoasymmetric assignments are retried against the stronger
descriptor snapshot. When Rule 5 supplies the unique ordering for a
tetrahedral center, assignment emits pseudoasymmetric `r`/`s` descriptors.
Double-bond descriptor assignment also applies the small-ring alkene exclusion
directly from topology, so manually inserted or cache-stale local double-bond
elements cannot receive lossy E/Z labels.

Stored axis descriptor assignment treats the axis bond endpoints as the two
ranking roots. The stored axis carriers are local reference carriers, one
adjacent to each axis endpoint; assignment ranks all carriers at each endpoint
and inverts the stored clockwise/counterclockwise handedness whenever the
stored reference carrier is not the highest-priority carrier at that endpoint.
Counterclockwise top-anchor handedness maps to `M`, and clockwise top-anchor
handedness maps to `P`, matching RDKit's atropisomeric bond convention.
Molfile atropisomeric wedge perception, including the conservative exocyclic
axis subset, remains in `stereo.perception`; this layer only consumes the
resulting stored axis element. For Molfile-derived axes, perception stores the
lowest-ID explicit carrier at each endpoint as the local reference pair, so CIP
priority flips are applied from a stable RDKit-like local axis convention
instead of from whichever wedge carrier happened to trigger perception.
The consumed perception subset includes redundant same-axis Molfile atrop
wedges and exocyclic axes with one ring endpoint plus one locally SP2 endpoint.
For axis endpoint ranking only, all-carbon aromatic bond components use a
uniform aromatic duplicate-node count so retained Molfile single/double
spelling in phenyl-like ligands does not change the `M`/`P` descriptor. This
normalization is deliberately scoped away from tetrahedral and ordinary
double-bond descriptor assignment, where existing mancude and source-order
regressions remain the controlling behavior until a full RDKit-style canonical
kekulizer is implemented.

The layer validates existing stereo by default and returns structured issues
instead of guessing when the current graph cannot support the stored local
stereo or when the implemented ranking rules cannot distinguish carriers.

## Validation

Unit tests cover tetrahedral descriptors, double-bond descriptors, recursive
Rule 1a/1b/2 ordering, mancude fractional atomic-number ordering, Rule 3
embedded E/Z ordering, Rule 4a descriptor-class ordering, Rule 4b
reference-descriptor and like/unlike pairing, Rule 4c pseudo-descriptor
ordering, Rule 5 descriptor-pair ordering, pseudoasymmetric tetrahedral
`r`/`s` assignment, sequence cis/trans descriptor-family ordering, Rule 6
reference-atom tie breaking for tetrahedral retry ranking,
natural-vs-indicated isotope priority, Rule 1b duplicate-node ordering,
negative-fraction duplicate expansion, implicit lone-pair carriers, unsupported
double-bond stereo exclusions including the ring-size boundary, unresolved
equivalent ligands, bounded resource failures, and descriptor invalidation
after mutation. Axis regressions cover stored local reference carriers,
endpoint priority flips, and `M`/`P` descriptor assignment. Targeted Rule 6
regressions cover both parity-stable and parity-unstable symmetric S4-style
reference retries. Auxiliary occurrence-graph regressions cover coupled
pseudoasymmetric cyclobutane, fused-ring, spiro-fused, and absolute-neighbor
bicyclic centers from the Enamine diversity corpus. Rule 2 regressions cover
natural atoms versus indicated isotopes and duplicate-node zero mass.
Equivalent-ring regressions cover isolated unsubstituted ring bridges that
must not become stereogenic through atom-id tie breaking.

Smoke, PubChem 100, PubChem 1k, PubChem 100k, and Enamine diversity validation
use externally supplied isomeric SMILES fixtures. Smoke validation also includes
official RDKit atropisomer Molfile fixtures with bond-centered `P`
descriptors, including alternate wedged substituent placement around the same
exocyclic axis, all-carbon aromatic source-Kekule variants for RP-6306 and
BMS-986142, redundant same-axis wedge marks, a JDQ443 heteromancude guardrail,
Mrtx1719 and ZM374979 axes, one-ring-endpoint SP2 axes, plus BMS, Sotorasib,
and ZM374979 fixtures that combine exocyclic Molfile atrop axes with
implicit-H tetrahedral centers adjacent to those axes.
PL-REX validation uses externally supplied ligand SDF packs to cover Molfile
and coordinate-bearing records. CIP goldens are generated with RDKit and
compare atom and bond
descriptor maps, not bytewise SMILES spelling or internal stereo element IDs.
Validation records include molecules where RDKit or the implementation assigns
at least one CIP descriptor; no-descriptor molecules are filtered out so broad
CIP validation is not dominated by unrelated parser or sanitizer coverage for
structures with no stereochemical labels. Bond descriptors are keyed by
endpoint atom indexes and descriptor instead of parser-local bond IDs, because
SMILES bond insertion order is not a portable chemical identity. Molecules
validation maps removable plain explicit hydrogens out of descriptor records
to match RDKit default atom indexing. PubChem 100k, the Enamine Discovery
Diversity Set, and PL-REX are enabled as broad RDKit parity gates for current
descriptor-bearing coverage.

## Out Of Scope

Full exact machine-oriented CIP coverage remains out of scope for this version:
perception or assignment of sequence cis/trans descriptors, kekulization of
aromatic-only inputs for mancude parity, remaining exact Rule 6 edge cases
beyond the parity-stable tetrahedral fallback, broad axis perception beyond
the conservative exocyclic Molfile atropisomeric wedge subset, full
RDKit-canonical aromatic kekulization for every alternate aromatic Molfile
layout, coordinate-only axis assignment, non-tetrahedral geometries beyond
stored axis descriptors, enhanced
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
- v11: Add RDKit-like Rule 4a descriptor-class ordering for uppercase sequence
  descriptors versus pseudo/geometric descriptors.
- v12: Add Rule 4c pseudo-descriptor ordering, Rule 5 descriptor-pair ordering,
  and pseudoasymmetric tetrahedral `r`/`s` assignment.
- v13: Add Rule 4b reference-descriptor selection and like/unlike descriptor
  family comparison.
- v14: Add contextual Rule 6 reference-atom tie breaking for unresolved
  tetrahedral carrier partitions.
- v15: Add `seqCis`/`seqTrans` descriptor vocabulary to sequence-rule ordering,
  including Rule 4a class and Rule 4b/5 descriptor-family handling.
- v16: Add RDKit-like mancude fractional atomic-number comparison for Rule 1a
  and negative-fraction duplicate expansion in the ligand digraph.
- v17: Generalize symmetric tetrahedral Rule 6 retry to test all successful
  atom-reference rankings and reject parity-unstable S4 outcomes.
- v18: Add RDKit-like Rule 2 mass ranking for natural atoms versus indicated
  isotopes while preserving zero mass for duplicate nodes.
- v19: Add staged descriptor-assignment rounds and deferred Rule 6 support for
  coupled pseudoasymmetric ring centers, including fused and small-ring
  descriptor timing cases.
- v20: Add deferred path-aware auxiliary tetrahedral descriptors so coupled
  pseudoasymmetric ligand paths can use local digraph labels instead of only
  molecule-level primary descriptors.
- v21: Promote Enamine diversity SMILES packs into the CIP validation contract
  as a broad RDKit parity gate for drug-like descriptor-bearing molecules.
- v22: Replace path-only auxiliary descriptor lookup with a root-centered
  auxiliary occurrence graph, virtual local re-rooting, absolute-before-pseudo
  deferred assignment, and Enamine regressions for fused, spiro, and coupled
  absolute/pseudo tetrahedral centers.
- v23: Apply the RDKit-like rule that double bonds in rings smaller than eight
  atoms are not E/Z stereogenic while preserving cyclooctene and larger
  endocyclic alkene assignment.
- v24: Add PL-REX ligand SDF packs to the CIP validation contract and compare
  every descriptor-bearing SDF record against RDKit-backed atom and bond
  descriptor maps.
- v25: Assign `M`/`P` CIP descriptors for structurally valid stored axis
  elements by ranking endpoint anchors and applying RDKit-like atropisomeric
  clockwise/counterclockwise handedness.
- v26: Validate Molfile wedge-derived atropisomeric axes against RDKit smoke
  goldens using the official RP-6306 atropisomer fixture.
- v27: Add an official RP-6306 alternate wedged substituent fixture to smoke
  parity validation, covering exocyclic-axis selection when ring-internal
  single bonds are also adjacent to the marked endpoint.
- v28: Use virtual implicit-H geometry for Molfile wedge-derived tetrahedral
  centers, store Molfile atrop axes with RDKit-style lowest-neighbor endpoint
  references, and add official BMS/Sotorasib atrop fixtures to smoke parity
  validation.
- v29: Normalize all-carbon aromatic duplicate counts during stored-axis
  endpoint ranking and add RP-6306/BMS source-Kekule atrop regressions plus a
  JDQ443 heteromancude guardrail.
- v30: Expand Molfile atrop parity coverage with official BMS/JDQ/Mrtx/RP,
  Sotorasib, and ZM374979 variants, including redundant same-axis wedges and
  one-ring-endpoint SP2 axes from the perception layer.
