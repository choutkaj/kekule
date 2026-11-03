# Conjugation and resonance

Kekule keeps four distinct layers:

- Aromaticity assigns aromatic atoms and bonds using the existing aromaticity models.
- Conjugation marks bonds using RDKit 2026.03.3's local electronic rules, including aromatic bonds. Default perception installs these flags after aromaticity.
- Resonance delocalization marks the non-aromatic bonds whose order changes when a formal charge moves, in connected groups. Default perception installs it after conjugation.
- Explicit enumeration generates localized contributors without changing the source molecule.

A bond is delocalized when aromaticity or resonance delocalizes it:
`aromaticity.bond_is_aromatic(id) || resonance.bond_is_delocalized(id)`.

## Conjugation

`perception::conjugation::perceive_conjugation` requires installed aromaticity and known hydrogen counts. `ConjugationModel::RdkitLike` reuses the low-level valence and electron helpers; it does not modify aromatic donor classification or ring evaluation. Conjugation is broader than delocalization: a neutral diene, enone or amide is conjugated, yet its bonds keep their orders.

## Resonance delocalization

`perception::resonance::perceive_resonance(&mut molecule, ResonanceModel::ChargeMoves)` requires installed aromaticity. `ResonancePerception` records the model and answers per bond whether resonance delocalizes it (`bond_is_delocalized`) and which group holds it (`bond_group`, `atom_group`). A group is a connected set of delocalized bonds with their atoms; groups are ordered by their smallest bond ID.

`ResonanceModel::ChargeMoves` moves formal charges:

- A mobile charge is an atom with formal charge ±1 and no radical. An anion is C, N, O, S or Se at −1 that can form one more bond as a neutral atom, such as carboxylate oxygen. An onium is N, O, P, S or Se at +1 with a complete octet, which can give up one bond as a neutral atom, such as iminium nitrogen. A vacant cation is C at +1.
- A move alternately raises and lowers bond orders along a simple path. An anion or a vacant cation first raises; an onium first lowers. Raising takes a single bond to double, or a double bond outside rings to triple. Lowering takes a double or triple bond to one order less. Aromatic, zero-order, dative and quadruple bonds never change.
- The path passes through B, C, N, O, Si, P, S, As and Se. Atoms inside the path keep their charge and valence.
- The path ends at a neutral atom that takes over the same kind of charge within its allowed valences: an anion after lowering onto C, N, O, S or Se; an onium after raising onto N, O, P, S or Se; a vacant cation after lowering onto C.
- Every bond of such a path is delocalized.

| Delocalized | Localized |
| --- | --- |
| carboxylate C–O bonds | neutral carboxylic acid |
| guanidinium and amidinium C–N bonds | amide, also drawn as a zwitterion |
| nitro and nitrate N–O bonds | diene, enone, styrene |
| terminal P–O of phosphates and phosphodiesters, S–O of sulfonates | P–OH, P–OR, sulfoxide, sulfone, N-oxide |
| enolate, 1,3-diketonate, amidate, azide, allyl ions | phenolate C–O, ammonium, iminium |

A charge moves only to an atom that then carries the same kind of charge, so moves keep the number of charges and complete octets. Neutral push-pull systems stay localized: an amide's second contributor separates charge, and so does an enone's or an enamine's. An amide drawn as a zwitterion is localized too, since its charges would cancel rather than move. An iminium's carbocation contributor gives up an octet, so a lone C=N⁺ bond stays localized.

Aromatic bonds belong to aromaticity, and a charge does not travel through an aromatic ring. Phenolate's C–O bond is localized, and p-nitrophenolate delocalizes its nitro group only. Protonated histidine (imidazolium) is aromatic and has no resonance group.

### Hypervalent phosphorus and sulfur

Phosphorus, arsenic, sulfur and selenium bonded to terminal chalcogens can be drawn with double bonds (P=O) or as ylides (P⁺–O⁻). The search reads each ylidic pair, an X⁺ single-bonded to a terminal O⁻, S⁻ or Se⁻ without hydrogens, as X=Y with both atoms neutral when that is a valid valence. Both drawings give the same groups: the terminal P–O bonds of a phosphate dianion or phosphodiester, the S–O bonds of a sulfonate. A neutral sulfoxide, sulfone or phosphine oxide has no group in either drawing. The molecule itself keeps its drawing.

### Protonation

Charges drive delocalization, so the result follows the represented protonation state. Structures read from mmCIF take the Chemical Component Dictionary's states: neutral carboxylic acids (COOH) and P–OH, charged arginine, lysine and histidine. A protein's carboxylates and phosphates delocalize only when they are represented as charged, for example after protonation assignment. Arginine's guanidinium is delocalized as read. Lysine's ammonium has no path, and histidine's imidazolium is aromatic.

### Cost

Only charged atoms start a search, and the search reads only the atoms and bonds its paths reach. A neutral molecule costs one scan of its formal charges. The path search is bounded by work that grows with the graph; exhaustion returns `ResonancePerceptionError::LimitExceeded` and leaves all state intact. Cancellable owner perception polls its cancellation during the search.

## Installed state

Detached construction uses `Perception::builder().with_resonance(model, groups)`. Installation requires aromaticity and checks that each group is non-empty, connected, made of live non-aromatic bonds not shared with another group, and that its atoms are exactly the bonds' endpoints. It does not rerun chemistry.

Valence, ring or aromaticity replacement invalidates conjugation and resonance. Resonance does not depend on conjugation, so replacing conjugation keeps it. Standalone conjugation, resonance perception and enumeration preserve installed CIP. Default owner workflows perceive each reusable definition transactionally; coordinates and trajectory frames remain separate.

## Contributor enumeration

`enumerate_resonance(&molecule, ResonanceOptions::default())` requires installed conjugation and returns a result borrowing its exact source molecule. Contributors vary over the connected conjugated groups, which `groups()` returns; these follow RDKit and differ from installed resonance groups. Contributors retain atom and bond IDs and expose complete formal-charge and localized-bond-order assignments. `to_molecule(index)` explicitly materializes a contributor, fixes known hydrogen counts, and invalidates derived state. Bond changes pass through the normal stereo-aware molecule editor.

All five RDKit options can be combined with `|`: `ALLOW_INCOMPLETE_OCTETS`, `ALLOW_CHARGE_SEPARATION`, `KEKULE_ALL`, `UNCONSTRAINED_CATIONS`, and `UNCONSTRAINED_ANIONS`. The unconstrained-cation option implies incomplete octets and charge separation; unconstrained anions imply charge separation. Default enumeration returns at most 1,000 contributors. `limit_reached()` means completeness is not claimed. The independent work budget returns an error on exhaustion, never a success containing an unfinished search. A zero structure limit returns no contributors. RDKit 2026.03.3 cannot materialize conjugated results at a structure limit of one; Kekule reports `InvalidStructureLimit` for that case.

No contributor enumeration runs during ordinary molecule, topology, model, ensemble or trajectory perception. No runtime RDKit dependency or SMARTS extension is introduced. Stereo and rotatable-bond algorithms retain their existing rules.

Earlier external comparisons against pinned RDKit 2026.03.3 covered conjugation, connected conjugated groups and contributor enumeration under all 32 option masks, with complete indexed charges, localized orders and contributor multiplicity. Focused synthetic examples are unit regressions. Resonance delocalization has no external reference; focused regressions pin its rules.

Exact parity is not yet complete at contributor cutoffs: tied choices can differ from the reference's C++ container and sorting behavior. See the [validation report](https://github.com/choutkaj/kekule/blob/ead90ca94f66e6b4f3082e5db49d6ddc042cc32c/benchmarks/RESONANCE-VALIDATION.md) for coverage and differences measured at its revision.
