# Stereochemistry Representation

## Summary

Store stereochemistry as a first-class layer adjacent to the molecular graph.
The stored truth is local stereo state plus relation groups, not CIP labels and
not atom/bond payload flags.

## Behavior/API

- `core::stereo` defines stable stereo element IDs, stereo group IDs,
  tetrahedral atom elements, double-bond elements, reserved axis elements,
  atom or implicit-hydrogen carriers, local orientations, specifiedness, source
  metadata, optional derived descriptors, source bond marks, stereo groups, and
  group kinds.
- `Molecule` stores stereo elements, stereo groups, and source bond marks with
  focused insertion, lookup, iteration, removal, and topology-aware pruning
  methods.
- Local stereo is the authoritative representation. `R/S`, `E/Z`, `M/P`, and
  pseudoasymmetric descriptors are optional derived descriptors and must be
  treated as cacheable views until the CIP feature is implemented.
- Unknown, unspecified, invalid-cleared, and specified stereo are distinct
  states. Missing stereo elements mean absent stereo, not explicit unknown
  stereo.
- Stereo groups model relation semantics separately from local parity, including
  absolute, relative, racemic, AND, and OR group kinds.
- SMILES `@`/`@@` markers are preserved as tetrahedral elements using SMILES
  local orientation and carrier order. Carrier order follows the SMILES-local
  sequence, including the incoming atom, bracket hydrogens, branches, and ring
  digits. SMILES `/` and `\` markers are preserved as source bond marks without
  double-bond perception.
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
  deleted atoms or bonds, and removes pruned stereo elements from relation
  groups.
- Topology or stereo mutation invalidates the stereo perception cache state.
- Source bond marks intentionally preserve parser syntax or Molfile wedge/either
  fields even before perception can assemble them into validated stereo
  elements.
- Macromolecules may carry stereo metadata through the shared graph, but
  small-molecule stereo perception is a later explicit workflow and should not
  run over whole `MacroMolecule` structures by default.

## Validation

- Unit tests cover stereo element, group, and source bond mark CRUD; invalid
  references; mutation invalidation; topology-aware pruning; and parser/writer
  adapter behavior.
- Smoke validation records semantic stereo JSON for externally supplied PubChem
  isomeric SMILES fixtures, including `stereo_elements`, `stereo_groups`,
  `stereo_bond_marks`, source marks, and specifiedness.

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
