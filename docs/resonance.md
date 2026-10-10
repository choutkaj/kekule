# Conjugation and resonance

Kekule keeps three distinct layers:

- Aromaticity assigns aromatic atoms and bonds using the existing aromaticity models.
- Conjugation marks bonds using RDKit 2026.03.3's local electronic rules, including aromatic bonds. Default perception installs these flags after aromaticity.
- Resonance preparation partitions conjugated bonds into connected groups. Explicit enumeration generates localized contributors without changing the source molecule.

`perception::conjugation::perceive_conjugation` requires installed aromaticity and known hydrogen counts. `ConjugationModel::RdkitLike` reuses the low-level valence and electron helpers; it does not modify aromatic donor classification or ring evaluation. Group preparation, `perception::resonance::perceive_resonance`, requires installed conjugation and installs only the connected partition. Groups can include aromatic and non-aromatic bonds. Membership does not imply equivalent bonds or equal contributor weights.

`Perception::has_conjugation()` and `has_resonance()` distinguish uncomputed state from a successfully computed empty result. The corresponding state accessors expose model provenance, membership and groups. Detached construction uses `Perception::builder()`; installation validates live references, endpoint sets, uniqueness, connectivity, and the complete partition without rerunning chemistry.

`enumerate_resonance(&molecule, ResonanceOptions::default())` returns a result borrowing its exact source molecule. Contributors retain atom and bond IDs and expose complete formal-charge and localized-bond-order assignments. `to_molecule(index)` explicitly materializes a contributor, fixes known hydrogen counts, and invalidates derived state. Bond changes pass through the normal stereo-aware molecule editor.

All five RDKit options can be combined with `|`: `ALLOW_INCOMPLETE_OCTETS`, `ALLOW_CHARGE_SEPARATION`, `KEKULE_ALL`, `UNCONSTRAINED_CATIONS`, and `UNCONSTRAINED_ANIONS`. The unconstrained-cation option implies incomplete octets and charge separation; unconstrained anions imply charge separation. Default enumeration returns at most 1,000 contributors. `limit_reached()` means completeness is not claimed. The independent work budget returns an error on exhaustion, never a success containing an unfinished search. A zero structure limit returns no contributors. RDKit 2026.03.3 cannot materialize conjugated results at a structure limit of one; Kekule reports `InvalidStructureLimit` for that case.

Valence, ring or aromaticity replacement invalidates conjugation and prepared groups. Replacing conjugation invalidates groups. Standalone conjugation, group preparation and enumeration preserve installed CIP. Default owner workflows perceive each reusable definition transactionally; coordinates and trajectory frames remain separate.

No contributor enumeration runs during ordinary molecule, topology, model, ensemble or trajectory perception. No runtime RDKit dependency, normalization, SMARTS extension, or file-format representation is introduced. Stereo and rotatable-bond algorithms retain their existing rules.

Earlier external comparisons against pinned RDKit 2026.03.3 covered conjugation, group preparation and contributor enumeration under all 32 option masks, with complete indexed charges, localized orders and contributor multiplicity. Focused synthetic examples are unit regressions.

Exact parity is not yet complete at contributor cutoffs: tied choices can differ from the reference's C++ container and sorting behavior. See the [validation report](https://github.com/choutkaj/kekule/blob/ead90ca94f66e6b4f3082e5db49d6ddc042cc32c/benchmarks/RESONANCE-VALIDATION.md) for coverage and differences measured at its revision.
