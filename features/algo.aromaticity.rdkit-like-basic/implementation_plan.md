# Implementation Plan

## Feature ID

`algo.aromaticity.rdkit-like-basic`

## Goal

Assign aromatic atom and bond flags for common organic rings using an initial RDKit-like model. This feature should depend on ring membership and should remain a perception step explicitly separate from parsing.

The goal is practical compatibility for common small molecules, not full RDKit parity for every edge case.

## Public API

Add an API equivalent to:

- `AromaticityModel`
- `AromaticityModel::RdkitLikeBasic`
- `perceive_aromaticity(mol: &mut Molecule, model: AromaticityModel) -> Result<(), AromaticityError>`
- Optional convenience on `SmallMolecule`: `sanitize_aromaticity_basic()` only if docs make clear it runs perception.

Expected behavior:

- Requires or computes ring membership.
- Sets `Atom::aromatic` and `Bond::aromatic` for perceived aromatic systems.
- Sets `PerceptionState::aromaticity` to fresh after successful perception.
- Leaves parsing and raw bond order loading separate.

## Internal Modules Touched

Expected scope:

- Aromaticity algorithm module or the current single crate file while small.
- Ring membership access from `algo.rings.fast`.
- `Atom` and `Bond` aromatic flags.
- Unit tests and optional validation fixtures.
- Feature docs under `features/algo.aromaticity.rdkit-like-basic/`.

Do not implement SDF parsing, valence sanitization, stereochemistry, or full kekulization in this feature.

## Data Model

Use existing atom and bond fields:

- `Atom::aromatic`
- `Bond::aromatic`
- `BondOrder::Aromatic`
- `PerceptionState::aromaticity`

The algorithm may need temporary per-atom electron contribution classifications:

- contributes 0
- contributes 1
- contributes 2
- unsupported or ambiguous

Temporary classifications should not be public API until the model is stable.

## Algorithm Outline

Initial RDKit-like basic model:

1. Ensure ring membership is available.
2. Build candidate components from ring atoms and ring bonds.
3. Restrict candidates to common organic aromatic elements: C, N, O, S, P, and selected charged variants.
4. Exclude atoms with unsupported valence, radicals, or ambiguous valence state.
5. Estimate pi-electron contribution from element, charge, explicit hydrogens, and bond order pattern.
6. For each candidate cyclic component, apply Huckel `4n + 2` electron rule.
7. Mark atoms and bonds in accepted components aromatic.
8. Leave unsupported candidate components non-aromatic with no hard parse failure unless the API chooses a perception warning.
9. Invalidate or overwrite previous aromatic flags deterministically before assigning new ones.

This feature should document every assumption and known divergence from RDKit.

## Tests

Add tests for:

- Benzene aromatic atoms and bonds.
- Cyclohexane non-aromatic.
- Cyclobutadiene non-aromatic.
- Pyridine aromatic.
- Pyrrole or furan aromatic if supported by the basic model.
- Acyclic conjugated chain non-aromatic.
- Ring membership must be considered; aromatic-looking acyclic bonds do not count.
- Existing aromatic flags are cleared or overwritten deterministically.
- Unsupported heteroatom ring stays non-aromatic or returns documented warning.
- Aromaticity perception sets aromaticity state to fresh.
- Topology mutation after perception invalidates aromaticity state.

## Reference Validation

Use RDKit as a reference tool through `validation.harness`.

Golden JSON should include:

- RDKit version.
- Input SMILES or SDF fixture.
- Atom aromatic flags.
- Bond aromatic flags.
- Notes for fixtures intentionally excluded from the basic model.

The comparison should focus on common organic rings supported by this feature. Do not claim full RDKit aromaticity parity.

## Required Checks

Run:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
cargo xtask validate --feature algo.aromaticity.rdkit-like-basic
```

## Risks

- RDKit aromaticity behavior includes many special cases; claiming broad parity too early would be misleading.
- Aromatic flags can be confused with raw parser-loaded aromatic bond orders.
- Valence perception is out of scope, but aromaticity often relies on valence assumptions.
- Heteroatom contribution rules can be subtle.
- Fused ring systems can require component-level reasoning beyond simple single-ring checks.

## Edge Cases

- Charged heteroatoms.
- Exocyclic double bonds.
- Fused aromatic systems.
- Five-member heteroaromatics.
- Aromatic bonds loaded from input before perception.
- Molecules with stale ring membership.
- Disconnected molecules.
- Unsupported elements in ring systems.

## Explicitly Out of Scope

- Full RDKit aromaticity parity.
- Valence perception or sanitization.
- Kekulization.
- Stereochemistry.
- Ring enumeration beyond membership needed for this model.
- Parser implementation.
- Runtime RDKit dependency.
