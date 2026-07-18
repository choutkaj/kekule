# Molecular Descriptors

## Summary

Provide read-only, explicit-policy molecular formula, average mass, and
monoisotopic mass descriptors for `SmallMolecule` values.

## Behavior/API

- The planned `descriptors` facade exposes `MolecularFormula`,
  `HydrogenCountPolicy`, `MolecularDescriptorError`, `molecular_formula`,
  `average_mass`, and `monoisotopic_mass`; it is not re-exported from the
  crate root or prelude.
- Every operation accepts `&SmallMolecule`, never mutates it, never sanitizes
  or perceives chemistry implicitly, and requires the caller to select one of
  two hydrogen policies:
  - `StoredOnly` counts live hydrogen atoms and stored explicit-hydrogen
    declarations but ignores perceived implicit hydrogens.
  - `IncludePerceived` additionally counts installed implicit hydrogens and
    returns a structured atom-qualified error when the count is absent for a
    live atom whose asserted state does not suppress implicit hydrogens.
- Formula construction retains separate counts for unlabeled elements and
  explicit isotopes, plus the aggregate asserted formal charge. It counts all
  live atoms across every graph component because one `Molecule` is one
  caller-asserted entity.
- `MolecularFormula` has private invariant-bearing state with read-only count,
  isotope, charge, and iteration accessors. Its stable text form uses Hill
  order: carbon, hydrogen, then remaining elements alphabetically when carbon
  is present; otherwise all elements alphabetically. Unlabeled atoms precede
  ascending explicit isotope mass numbers for the same element, isotopes use
  `[13C]`-style notation without `D`/`T` aliases, counts of one are omitted,
  and total charge is appended as `+`, `-`, `+2`, `-2`, and so on.
- `average_mass` returns `Quantity<f64>` in `units::DALTON`. An explicitly
  labeled atom uses its isotope mass; an unlabeled atom uses the scalar 2024
  CIAAW abridged standard atomic weight. It returns a structured error for an
  unknown isotope or an unlabeled element without a standard atomic weight.
- `monoisotopic_mass` returns `Quantity<f64>` in `units::DALTON`. An explicitly
  labeled atom uses its isotope mass; an unlabeled atom uses the exact mass of
  the most abundant naturally occurring isotope from the pinned reference
  composition. It returns a structured error when either value is unavailable
  and never substitutes an integer mass number.
- Both mass functions apply the total formal-charge correction
  `mass -= charge * electron_mass`, using the pinned electron rest mass in
  daltons. Radical state needs no separate correction because formal charge
  determines the electron-count difference from the neutral constituent
  atoms.
- Empty molecules produce an empty zero-charge formula and zero-dalton masses.
  Count accumulation is checked and reports overflow rather than wrapping.

## Implementation Notes

- Formula and mass calculations are one pass over live atoms, with time linear
  in live atoms plus stored hydrogen declarations and memory linear in the
  number of distinct element/isotope terms. They do not traverse bonds,
  conformers, properties, or connected components and use no recursion.
- Stored explicit-hydrogen declarations represent unlabeled hydrogen. An
  isotopic hydrogen must be a live isotope-labeled hydrogen atom.
- The implementation vendors immutable, source-pinned numeric tables rather
  than adding a Rust runtime or network dependency:
  - scalar average weights from the
    [CIAAW 2024 abridged standard atomic weights](https://ciaaw.org/abridged-atomic-weights.htm);
  - natural-isotope selection from the
    [CIAAW 2024 isotopic compositions](https://ciaaw.org/isotopic-abundances.htm);
  - exact isotope masses from
    [AME 2020](https://amdc.impcas.ac.cn/web/masseval.html);
  - electron mass from the
    [2022 CODATA recommended constants](https://physics.nist.gov/cuu/Constants/).
- Standard atomic weights describe normal materials and some are intervals.
  The average-mass contract deliberately uses CIAAW's abridged scalar values;
  uncertainty propagation and material-specific isotope distributions remain
  outside this feature.
- Public errors are non-exhaustive and identify the atom and unavailable
  element/isotope or missing implicit-hydrogen state where applicable.

## Validation

- Before promotion to `experimental`, add `pubchem-1k` as required validation
  and generate provenance-pinned reference rows with the repository's pinned
  RDKit environment.
- Compare formula composition, isotope counts, total charge, average mass, and
  monoisotopic mass independently. Formula comparison is structured rather
  than relaxed string comparison; mass comparison uses explicit tolerances
  narrow enough to detect table or hydrogen-count drift.
- RDKit `CalcMolFormula` with separated, non-abbreviated isotopes, `MolWt`, and
  `ExactMolWt` are the common-organic reference behavior. Differences caused
  by newer pinned standards, charge-mass correction, or deliberate rejection
  of elements without defined reference values must be explicit expected
  outcomes, not silently excluded records.
- Focused unit regressions cover Hill ordering with and without carbon,
  disconnected salts, explicit and perceived hydrogens, isotope labels,
  positive and negative ions, radicals, empty molecules, missing perception,
  unavailable standard weights or isotopes, and checked count overflow.

## Out Of Scope

- LogP, polar surface area, hydrogen-bond counts, rotatable-bond counts,
  fragment descriptors, fingerprints, and 3D descriptors.
- Automatic sanitization, valence perception, protonation, neutralization,
  fragment selection, isotope-distribution envelopes, abundance-weighted
  isotopologue probabilities, uncertainty propagation, and material-specific
  atomic weights.
- Macromolecule- or model-level aggregation and component-separated formula
  rendering.
- Ionization-energy, molecular binding-energy, and other high-precision mass
  corrections beyond the documented additive atomic-mass and electron-mass
  model.

## Revision Notes

- v1: Feature contract reserved.
- v2: Define the explicit hydrogen policy, structured formula representation,
  average and monoisotopic mass semantics, authoritative data provenance,
  charge correction, resource bounds, validation gate, and focused scope.
