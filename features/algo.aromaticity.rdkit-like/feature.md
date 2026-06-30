# RDKit-like Aromaticity Perception

## Summary

Assign aromatic atom and bond flags for common organic ring systems using the RDKit-like model. This is a perception step and must remain separate from parsing.

## Behavior/API

- Exposes `AromaticityModel::RdkitLike`, `perceive_aromaticity`, and `perceive_aromaticity_with_ring_options`.
- Requires or computes ring perception before assigning aromaticity.
- Marks supported aromatic atoms and bonds and sets aromaticity perception state to fresh.
- Clears prior aromatic flags deterministically before assignment.
- Can be run directly or through the explicit sanitization pipeline.
- Returns `UnsupportedElement` for an explicitly aromatic ring containing an unsupported element instead of silently accepting that representation.
- Returns `InvalidAromaticRepresentation` when imported aromatic bonds cannot be perceived back onto every participating atom.
- Propagates bounded ring-perception failures as `AromaticityError::RingPerception`.

## Implementation Notes

- Operates on the shared core `Molecule` graph.
- Uses per-ring cycle data from `algo.rings.sssr`.
- Integrates with the first-wave valence and ring-set perception stack.
- Applies a 4n+2 electron-count model for common C, N, O, S, Se, Te, and P rings and small fused ring components.
- Computes pi-electron counts from bond order, and uses an atom-contribution path for explicitly imported aromatic-bond rings.
- Accepts imported aromatic-bond rings when RDKit-like variable atom-contribution ranges contain a valid 4n+2 count.
- Routes imported aromatic-order atom contributions through a shared RDKit-style donor classifier with vacant, one-electron, two-electron, and variable donor ranges.
- Routes localized monocyclic and fused-component pi-electron counts through the same donor classifier, so ring-pi atoms, lone-pair donors, anionic carbon donors, vacant exocyclic centers, and non-donors are handled by one path.
- Uses the same donor classifier for aromatic candidate-path and active hetero-donor checks instead of duplicating separate symbol-based rules.
- Reuses one localized per-ring donor analysis for candidate-path and active hetero-donor predicates, mirroring RDKit's table-driven ring candidate flow.
- Caches explicit per-ring aromaticity analysis objects in the top-level perception pass so ring aromaticity and adjacent protection decisions share the same classified donor state.
- Reuses cached per-ring donor analysis when selecting fused-ring candidates, avoiding a second donor classification pass.
- Reuses cached per-ring donor analysis when selecting fused single-exocyclic-carbon rings.
- Reuses cached per-ring donor analysis in fused lactone and saturated ether cleanup decisions.
- Reuses cached per-ring donor analysis in saturated tertiary amine cleanup and fused carbon/nitrogen fallback decisions.
- Feeds initial ring aromaticity gates from the same localized donor table used for Huckel counting, reducing raw hetero-atom presence checks.
- Uses cached heavy-chalcogen candidate state for the five-member terminal chalcogen exocyclic-pi gate rather than raw terminal chalcogen shape.
- Uses cached active donor state for fused fallback admission and terminal-exocyclic atom retention instead of raw nitrogen/chalcogen presence.
- Uses cached active neighbor donor state for exocyclic alkene chalcogen cleanup rather than raw nitrogen/chalcogen adjacency.
- Clears only cached active ring hetero-donor neighbors during exocyclic alkene chalcogen cleanup, avoiding raw hetero-atom adjacency outside the ring context.
- Uses cached inactive chalcogen donor state for ring-oxo chalcogen cleanup instead of raw chalcogen/terminal-oxo matching alone.
- Reuses RDKit-like radical eligibility when admitting inactive terminal and ring-oxo chalcogen cleanup atoms, so radical heteroatoms do not enter cleanup solely by raw element shape.
- Uses cached active oxygen donor state for saturated fused ether cleanup rather than raw hetero-atom presence.
- Uses ring-local donor gates for fused lactam/enone and saturated fused ether cleanup rather than global heavier-chalcogen molecule skips.
- Uses cached active chalcogen donor state for fused lactone bridge cleanup rather than raw chalcogen presence.
- Reuses RDKit-like radical eligibility when recognizing saturated chalcogen bridge atoms, so radical heteroatoms do not enter bridge cleanup by raw bridge shape alone.
- Uses cached inactive chalcogen donor state for terminal chalcogen-oxo cleanup rather than raw terminal chalcogen exocyclic-pi matching alone.
- Uses cached active nitrogen and oxygen donor state for fused lactam-bridge cleanup rather than raw hetero-atom presence.
- Uses cached candidate nitrogen state for imide carbonyl cleanup rather than raw nitrogen presence, preserving vacant cationic candidates while rejecting non-candidate nitrogens.
- Uses cached active nitrogen donor state for fused lactam/enone cleanup rather than raw nitrogen presence.
- Uses cached active nitrogen donor state for terminal chalcogen-oxo cleanup rather than raw nitrogen presence.
- Uses cached active nitrogen donor state for saturated fused lactam-carbonyl cleanup rather than raw nitrogen presence.
- Uses cached active chalcogen donor state for the saturated tertiary amine guard rather than raw saturated chalcogen presence.
- Uses cached active nitrogen donor state when selecting fused lactam/enone carbon neighbors, so inactive nitrogen atoms do not satisfy the per-carbon cleanup selector by raw neutral nitrogen shape alone.
- Reuses RDKit-like radical eligibility when classifying saturated tertiary amine guard nitrogens, so radical heteroatoms do not enter cleanup by raw saturated-amine shape alone.
- Reuses RDKit-like radical eligibility for terminal chalcogen substituents in the saturated tertiary amine guard, so radical heteroatom substituents do not satisfy the guard by raw terminal-pi shape alone.
- Reuses RDKit-like radical eligibility for carbon substituents in the saturated tertiary amine guard, so charged carbon radicals do not satisfy the guard by raw carbon shape alone.
- Uses cached chalcogen candidate state for saturated fused chalcogen bridge cleanup rather than raw saturated bridge shape alone.
- Uses cached candidate nitrogen state when admitting carbon/nitrogen fused fallback components rather than raw nitrogen presence.
- Uses RDKit-like nitrogen candidate eligibility for aromatic amidine cleanup rather than raw neighboring nitrogen presence.
- Uses RDKit-like carbon candidate eligibility for aromatic amidine cleanup rather than raw aromatic carbon shape.
- Uses RDKit-like carbon and nitrogen candidate eligibility for terminal aromatic imine cleanup rather than raw aromatic `N=C` symbol shape.
- Uses electronegativity-aware exocyclic pi-bond checks for aromatic-order carbon electron stealing rather than treating every hetero exocyclic atom as equivalent.
- Reuses localized active donor state for imported aromatic-order exocyclic-carbon electron stealing instead of raw nitrogen/chalcogen presence.
- Handles single electronegative exocyclic pi bonds in imported six-member chalcogen aromatic-order rings through the donor classifier instead of a fixed six-electron override.
- Uses cached localized donor state, rather than raw hetero-element presence, when deciding which non-aromatic fused single bonds should remain protected.
- Reuses one localized donor analysis inside saturated tertiary amine fused cleanup instead of recalculating candidate-path and active-donor facts separately.
- Routes nitrogen lone-pair donor checks through the same per-ring donor analysis instead of separate hydrogen/charge symbol logic.
- Evaluates simple and fused Huckel aromaticity through the shared donor-analysis object, keeping all-candidate and electron-count checks together.
- Handles fused all-carbon components with exocyclic pi bonds through classified donor state instead of a fixed six-electron override.
- Admits low-unsaturation hetero fused candidates from active donor analysis rather than hetero atom presence alone.
- Applies a bounded connected fused-subset Huckel search before older fused fallback heuristics, marking accepted subsets additively.
- Iterates accepted fused ring subsets additively with an RDKit-like done-bond stop condition instead of returning only the first accepted subset.
- Uses RDKit's 24-atom fused-ring candidate size limit for large carbon/nitrogen fused candidates instead of an older smaller local cap.
- Counts fused-system Huckel donor atoms using RDKit-style fused-subsystem atom multiplicity, excluding atoms buried in more than two selected rings.
- Marks accepted all-carbon fused Huckel systems with RDKit-style selected-subsystem perimeter bonds, preserving non-aromatic shared bonds for azulene-like systems.
- Uses the same RDKit-style selected-subsystem perimeter-bond definition for fused-system marking and fused-combination completion tracking.
- Rejects aromatic candidate atoms with more than one explicit double/triple bond through the shared donor classifier.
- Applies RDKit-like radical candidate eligibility in the shared donor classifier: radical heteroatoms and charged radical carbons are non-candidates, while neutral carbon radicals remain eligible.
- Applies RDKit-like candidate coordination eligibility, rejecting atoms whose bond plus hydrogen degree exceeds three before donor typing.
- Applies RDKit-like default-valence candidate eligibility, rejecting atoms whose total valence exceeds the charge-adjusted default valence before donor typing.
- Gates RDKit-like candidate admission through the classified donor type plus candidate options, including the exocyclic multiple-bond toggle used by RDKit's model helpers.
- Counts localized saturated, vacant, and lone-pair atom donors with an RDKit `countAtomElec`-style helper using default valence, outer-shell electrons, charge, radical electrons, effective hydrogens, and exocyclic electronegativity.
- Allows localized two-electron Huckel rings, preserving RDKit-like cyclopropenyl cation aromaticity.
- Evaluates localized simple rings of arbitrary size through the shared Hückel donor-count path instead of applying a small-ring-only cutoff.
- Treats terminal hetero exocyclic pi carbons as non-donating in imported six-member nitrogen/chalcogen aromatic-order rings with multiple terminal exocyclic pi bonds.
- Treats exocyclic hetero pi carbons as non-donating in imported five-member nitrogen/chalcogen aromatic-order rings when needed for RDKit-like Huckel counts.
- Clears terminal aromatic imine fragments and orphan aromatic atoms left outside any aromatic bond path after fused-subsystem cleanup.
- Clears five-member neutral imide carbonyl ring atoms when a saturated ring nitrogen is flanked by two terminal ring carbonyls.
- Clears saturated fused all-carbon ring atoms that are not retained by neighboring aromatic rings, saturated aromatic carbon centers, and localized cyclic amidine centers.
- Reuses RDKit-like radical eligibility for saturated ring carbon cleanup, so charged carbon radicals do not satisfy cleanup gates by raw saturated-carbon shape alone.
- Clears localized fused lactam/enone bridge carbons and saturated fused oxygen bridge atoms when they are not part of a conjugated RDKit-like aromatic path.
- Clears saturated fused nitrogen carbonyl ring atoms that are not retained by a neighboring aromatic ring.
- Uses conservative guards for small rings, hetero fused donors, lactone-like rings, and large macrocycles exposed by external PubChem validation.
- Treats unsupported ring elements as non-aromatic for the current model rather than failing the whole perception pass.
- Leaves unsupported or ambiguous systems non-aromatic rather than claiming full RDKit parity.

## Validation

- Unit tests cover common monocyclic organic rings and stale-state behavior.
- RDKit-generated goldens compare aromatic atom and bond flags for external PubChem fixtures.

## Out Of Scope

- Full RDKit aromaticity parity.
- RDKit-like aromatic bond selection for all fused systems; PubChem-1000 exposes cases where atom aromaticity and bond aromaticity diverge.
- Valence perception, sanitization, kekulization, stereochemistry, and parser behavior.
- Runtime RDKit dependency.

## Revision Notes

- v1: Aromaticity perception for common organic rings.
- v2: Document integration with explicit sanitization and ring/valence perception.
- v3: Per-ring fused aromaticity heuristic passes the RDKit-backed `tiny` corpus; broader required corpora remain pending.
- v4: Add fused-component aromaticity and order-based electron counting for external PubChem fused-ring systems.
- v5: Refine fused heteroaromatic handling and conservative ring-size/electron-count guards to pass PubChem-100.
- v6: Add chalcogen heteroaromatic support and refine fused donor eligibility; PubChem-1000 still exposes fused aromatic bond-selection gaps.
- v7: Reject unsupported elements in explicitly aromatic ring representations and preserve caller state when sanitization propagates the error.
- v8: Count explicitly imported aromatic-bond rings with atom contributions so lowercase aromatic SMILES sanitize without treating every aromatic bond as a localized double bond.
- v9: Propagate configurable structured ring resource limits before mutating aromatic flags.
- v10: Narrow saturated tertiary amine fused-ring clearing to carbon-substituted amines, preserving N-O substituted lactam aromaticity.
- v11: Preserve valid imported aromatic SMILES components while clearing saturated fused thioether and ring-oxo chalcogen bridges.
- v12: Refine imported aromatic nitrogen, pyrone, fused lactone, saturated fused carbon, and fluorenone-like carbonyl bridge handling exposed by canonical SMILES PubChem validation.
- v13: Broaden fused carbonyl bridge and cationic imide cleanup using exocyclic-pi and saturated-bridge criteria; PubChem-1000 still exposes fused subsystem selection gaps.
- v14: Add RDKit-like aromatic-order electron ranges plus final aromatic consistency cleanup for terminal imine and orphan aromatic fragments.
- v15: Generalize imide carbonyl cleanup from cationic systems to neutral five-member imides with saturated nitrogen between two ring carbonyls.
- v16: Refine fused saturated carbon and cyclic amidine cleanup so sp3/enone fused atoms do not remain aromatic while valid aromatic iminium and carbonyl systems are preserved.
- v17: Refine saturated tertiary amine cleanup to treat oxidized chalcogen substituents as non-donor sulfone-like groups, preserving aliphatic sulfonamide ring nitrogens.
- v18: Clear RDKit-like exocyclic alkene ring carbons between nitrogen and chalcogen donors so fused thiazine/thiazolium systems do not over-aromatize deactivated atoms.
- v19: Count terminal hetero exocyclic pi carbons as non-donating in imported nitrogen/chalcogen aromatic-order rings with multiple terminal exocyclic pi bonds, preserving RDKit-like thione-rich heterocycles without over-aromatizing singly carbonylated fused systems.
- v20: Clear localized fused lactam/enone bridge carbons and saturated fused oxygen bridge atoms, preserving RDKit-like canonical reparse semantics for oxygen/nitrogen polycyclic lactam systems while leaving heavier-chalcogen fused aromaticity on its existing path.
- v21: Clear saturated fused nitrogen carbonyl ring atoms outside neighboring aromatic rings, preserving RDKit-like benzodiazepinone lactam canonical reparse semantics.
- v22: Count exocyclic hetero pi carbons as non-donating in imported five-member nitrogen/chalcogen aromatic-order rings, preserving cationic imine thiadiazolium canonical reparse semantics.
- v23: Centralize imported aromatic-order atom contributions behind a RDKit-style donor classifier and cover a PubChem macrocycle anionic-nitrogen canonical round trip.
- v24: Use the shared donor classifier for localized monocyclic and fused-component electron counts, preserving anionic carbon donors while rejecting neutral saturated carbon non-donors.
- v25: Route conjugated candidate-path and active hetero-donor checks through the shared donor classifier, with a phosphorus lone-pair donor regression.
- v26: Add bounded connected fused-subset Huckel search ahead of legacy fused fallback marking, preserving additive subsystem aromaticity exposed by PubChem.
- v27: Count fused-system donor atoms with RDKit-style atom multiplicity and reject candidate atoms with more than one explicit pi bond.
- v28: Add RDKit-like radical candidate eligibility to the shared donor classifier, preserving neutral carbon radicals while rejecting radical heteroatoms and charged radical carbons.
- v29: Add RDKit-like candidate coordination eligibility, rejecting over-coordinated ring atoms before aromatic donor typing.
- v30: Count localized saturated and vacant donors through an RDKit-style atom-electron helper and accept localized two-electron Huckel rings.
- v31: Mark accepted all-carbon fused Huckel systems with RDKit-style perimeter-bond selection, keeping azulene-like shared bonds non-aromatic.
- v32: Iterate accepted fused ring subsets additively with a done-bond stop condition, moving the fused search closer to RDKit's combination loop.
- v33: Remove the seven-atom simple-ring cutoff so localized larger rings are evaluated by the RDKit-like Hückel donor-count path.
- v34: Reject aromatic candidate atoms whose total valence exceeds RDKit's charge-adjusted default valence, preserving substituted hypervalent phosphorus rings as non-aromatic.
- v35: Reuse RDKit-style selected-subsystem perimeter bonds for both fused-system marking and fused-combination done-bond tracking.
- v36: Align the large fused carbon/nitrogen candidate size limit with RDKit's 24-atom fused-ring cap.
- v37: Use active donor analysis, not hetero atom presence alone, for low-unsaturation hetero fused candidate admission.
- v38: Protect non-aromatic fused single bonds using cached localized active-donor state instead of raw hetero-element counts.
- v39: Reuse localized donor analysis as the input to initial ring aromaticity gates and the simple Huckel fallback.
- v40: Use active donor state for fused fallback admission and fused-context terminal-exocyclic atom retention.
- v41: Use active donor state for exocyclic alkene chalcogen cleanup seed and neighbor checks.
- v42: Use active oxygen donor state for saturated fused ether cleanup instead of raw oxygen/chalcogen presence.
- v43: Use active chalcogen donor state for fused lactone bridge cleanup instead of raw chalcogen presence.
- v44: Use active nitrogen donor state for fused lactam/enone cleanup instead of raw nitrogen presence.
- v45: Use active localized donor state for imported aromatic-order exocyclic-carbon electron stealing instead of raw nitrogen/chalcogen presence.
- v46: Use active nitrogen donor state for terminal chalcogen-oxo cleanup instead of raw nitrogen presence.
- v47: Use active chalcogen donor state for the saturated tertiary amine guard instead of raw saturated chalcogen presence.
- v48: Use active nitrogen donor state for saturated fused lactam-carbonyl cleanup instead of raw nitrogen presence.
- v49: Use active nitrogen and oxygen donor state for fused lactam-bridge cleanup instead of raw hetero-atom presence.
- v50: Use candidate nitrogen state for imide carbonyl cleanup instead of raw nitrogen presence.
- v51: Use candidate nitrogen state for carbon/nitrogen fused fallback admission instead of raw nitrogen presence.
- v52: Clear only active ring hetero-donor neighbors during exocyclic alkene chalcogen cleanup instead of raw hetero-atom neighbors.
- v53: Replace global heavier-chalcogen molecule skips in lactam/enone and saturated ether cleanup with ring-local donor gates.
- v54: Route ring-oxo chalcogen cleanup through cached inactive chalcogen donor state.
- v55: Route terminal chalcogen-oxo cleanup through cached inactive chalcogen donor state.
- v56: Route saturated fused chalcogen bridge cleanup through cached chalcogen candidate state.
- v57: Route aromatic amidine cleanup through RDKit-like nitrogen candidate eligibility.
- v58: Route exocyclic pi-bond electron stealing gates through electronegativity-aware checks.
- v59: Reuse RDKit-like radical eligibility for inactive terminal and ring-oxo chalcogen cleanup, avoiding raw terminal chalcogen shape checks for radical heteroatoms.
- v60: Reuse RDKit-like radical eligibility for saturated tertiary amine guard nitrogens, avoiding raw saturated-amine shape checks for radical heteroatoms.
- v61: Reuse RDKit-like radical eligibility for saturated chalcogen bridge recognition, avoiding raw bridge-shape checks for radical heteroatoms.
- v62: Reuse RDKit-like radical eligibility for terminal chalcogen substituents in the saturated tertiary amine guard.
- v63: Reuse RDKit-like radical eligibility for carbon substituents in the saturated tertiary amine guard.
- v64: Use localized heavy-chalcogen candidate state for the five-member terminal chalcogen exocyclic-pi gate.
- v65: Reuse RDKit-like radical eligibility for saturated ring carbon cleanup.
- v66: Use active neighbor nitrogen donor state when selecting fused lactam/enone cleanup carbons.
- v67: Use RDKit-like candidate eligibility for terminal aromatic imine cleanup.
- v68: Use RDKit-like carbon candidate eligibility for aromatic amidine cleanup.
- v69: Handle single electronegative exocyclic pi bonds in imported six-member chalcogen aromatic-order rings through donor classification instead of a fixed electron-count override.
- v70: Handle fused all-carbon components with exocyclic pi bonds through donor classification instead of a fixed six-electron override.
