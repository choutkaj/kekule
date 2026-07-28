# Stereochemistry Representation

## Summary

Store stereochemistry as a first-class layer adjacent to the molecular graph.
The stored truth is local stereo state plus relation groups, not CIP labels and
not atom/bond payload flags.

## Behavior/API

- `core::stereo` defines stable stereo element IDs, stereo group IDs,
  tetrahedral atom elements, double-bond elements, reserved axis elements,
  atom, implicit-hydrogen, or implicit-lone-pair carriers, local orientations,
  specifiedness, source metadata, optional derived descriptors, source bond
  marks, stereo groups, and group kinds.
- `Molecule` stores stereo elements, stereo groups, and source bond marks with
  focused insertion, lookup, iteration, removal, and topology-aware pruning
  methods.
- Stereo-element and stereo-group insertion checks the fixed-width ID space
  before commit and returns the focused molecule capacity error when exhausted.
- Stereo-element insertion accepts only ungrouped elements. A pre-set relation
  group is rejected transactionally before reference validation, slot
  allocation, mutation, or perception invalidation; membership is established
  only through `add_stereo_group`.
- Stereo-element replacement validates every carrier and central reference
  before commit. An element already in a relation group cannot change group
  membership implicitly; callers must use the group operations explicitly.
- Stereo-element removal prunes or tombstones its prior live group as required
  and returns the detached element with `group = None`, ready for explicit
  re-insertion and regrouping.
- Stereo groups must contain at least one unique live stereo element.
- Complete stereo-group slot iteration includes interior and trailing
  tombstones. A checked append operation reserves one deleted stable ID without
  creating dummy chemistry or invalidating CIP.
- Local stereo is the authoritative representation. `R/S`, `E/Z`, `M/P`,
  sequence cis/trans, and pseudoasymmetric descriptors are optional derived
  descriptors and must be treated as cacheable views over local stereo.
- Unknown, unspecified, invalid-cleared, and specified stereo are distinct
  states. Missing stereo elements mean absent stereo, not explicit unknown
  stereo.
- Stereo groups model relation semantics separately from local parity, including
  absolute, relative, racemic, AND, and OR group kinds.
- SMILES `@`/`@@` markers are preserved as tetrahedral elements using SMILES
  local orientation and carrier order. Carrier order follows the SMILES-local
  sequence, including the incoming atom, bracket hydrogens, branches, ring
  digits, and supported implicit lone-pair placeholders for three-neighbor
  heteroatom centers. SMILES `/` and `\` markers are preserved as source bond
  marks without double-bond perception.
- Supported V2000 and V3000 bond stereo fields are preserved as source bond
  marks. Atom `CFG`/parity and enhanced stereo collections remain unsupported
  until explicit format features are implemented.
- Writers read the new stereo model. Non-isomeric SMILES rejects stored stereo
  unless canonical non-isomeric output explicitly opts to ignore it; Molfile
  writers emit supported source bond marks and reject unsupported stereo
  elements or incompatible bond-mark/order combinations.

## Implementation Notes

- `Atom::chiral`, `AtomStereo`, `Bond::stereo`, and `BondStereo` are no longer
  the authoritative public model. Core atom and bond payloads remain chemically
  general graph payloads; stereo lives on `Molecule`.
- Topology deletion prunes stereo elements and source bond marks that reference
  deleted atoms or bonds, removes pruned stereo elements from relation groups,
  and tombstones a group when its final member is removed.
- Topology or stereo mutation invalidates the stereo perception cache state.
- Source bond marks intentionally preserve parser syntax or Molfile wedge/either
  fields even before perception can assemble them into validated stereo
  elements.
- Implicit lone-pair carriers are local stereo placeholders only. They preserve
  supported imported tetrahedral syntax for heteroatom centers without
  converting lone pairs into graph atoms.
- Macromolecules may carry stereo metadata through the shared graph, but
  small-molecule stereo perception is a later explicit workflow and should not
  run over whole `MacroMolecule` structures by default.

## Tests

- Unit tests cover stereo element, group, and source bond mark CRUD; invalid
  references; synthetic ID-capacity boundaries; mutation invalidation;
  topology-aware pruning; exact live/tombstone slot replay; stable next IDs;
  CIP-neutral tombstone append; transactional rejection of pre-grouped
  insertion; detached removal/re-insertion; and parser/writer adapter behavior.

## Benchmarks

- Smoke, PubChem 100, PubChem 1k, PubChem 100k, Enamine diversity, and PL-REX
  validation record semantic stereo JSON for externally supplied isomeric
  SMILES fixtures, including `stereo_elements`, `stereo_groups`,
  `stereo_bond_marks`, source marks, and specifiedness. PL-REX adds
  coordinate-bearing ligand SDF packs so Molfile wedge/either and source-mark
  preservation stay covered outside SMILES-only fixtures. The broader PubChem,
  Enamine, and PL-REX tiers are implementation-golden semantic regression gates
  for representation stability, while exact RDKit descriptor parity belongs to
  `stereo.cip`.
- Optional external-reference manifests are available for `pubchem-1k`, `pubchem-100k`, `pl-rex`, `enamine-diversity`.
- Benchmark observations are informational and never determine this feature's release status or repository health.

## Out Of Scope

- Candidate stereo perception, coordinate/wedge assignment, local stereo
  validation, exact CIP assignment, isomeric SMILES writing, enhanced
  V3000/CXSMILES stereo, stereo enumeration, and reaction stereo transfer.

## Revision Notes

- v1: Feature contract reserved.
- v2: Add first-class `core::stereo` representation, graph-adjacent storage on
  `Molecule`, parser preservation adapters, writer rejection/mark handling, and
  smoke semantic stereo validation.
- v3: Generalize double-bond stereo carriers from atom-only IDs to
  `StereoCarrier` so alkene perception can represent implicit-hydrogen
  substituents.
- v4: Preserve SMILES-local tetrahedral carrier order for bracket hydrogens and
  ring-digit closures in smoke semantic validation.
- v5: Add implicit lone-pair stereo carriers so supported three-neighbor
  heteroatom tetrahedral markers can be represented without adding graph atoms.
- v6: Add sequence cis/trans entries to the derived descriptor vocabulary.
- v7: Add PubChem 100 and PubChem 1k semantic regression requirements for
  stereo representation over externally supplied isomeric SMILES.
- v8: Add PubChem 100k and Enamine diversity semantic regression requirements
  for broader drug-like stereo representation preservation coverage.
- v9: Add PL-REX ligand SDF packs to the representation benchmark contract for
  coordinate- and Molfile-stereo source-mark regression coverage.
- v10: Keep every ignored non-smoke corpus as explicit local-only validation
  instead of repository-wide required evidence.
- v11: Replace unchecked mutable stereo-element access with validated
  transactional replacement and reject empty or duplicate-member groups.
- v12: Use PubChem-1k as the required baseline benchmark corpus after retiring the former smoke corpus from public validation.
- v13: Reclassify external-reference parity from a required gate to optional benchmarking without changing implementation behavior or golden expectations.
- v14: Check stereo-element and stereo-group identifier capacity before
  transactional insertion and iterate their stable slots without narrowing.
- v15: Expose exact stereo-group slots plus checked tombstone append, and
  tombstone groups whose final live member is pruned.
- v16: Reject pre-grouped stereo-element insertion before mutation or
  invalidation, and return removed stereo elements detached from their prior
  relation group.
