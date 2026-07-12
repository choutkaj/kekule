# DREIDING Force-Field Adapter

## Summary

Provide an explicit SmallMolecule modelling adapter that prepares DREIDING atom types,
fixed QEq charges, bonded terms, nonbonded terms, and complete Cartesian gradients.

## Behavior/API

- Exposes `DreidingPotential` and `DreidingPrepareError` from the separate
  `molecules-dreiding` crate.
- Prepares a potential with `DreidingPotential::prepare(&MolecularModel)` and implements
  the core `Potential` evaluation contract.
- Binds preparation to the source model definition, accepting coordinate-modified clones
  and rejecting independently built models.
- Exposes read-only per-atom type and partial-charge diagnostics.
- Rejects unresolved implicit-hydrogen state; every atom must carry an explicit zero
  implicit-hydrogen count or a no-implicit-hydrogens assertion.
- Uses angstrom coordinates, kJ/mol energies, and kJ/mol/angstrom gradients.

## Implementation Notes

- Uses pinned `dreid-forge` and matching `dreid-kernel` releases; upstream types do not
  cross the adapter's public API.
- Uses the shared model-definition key rather than rebuilding an adapter-specific topology
  signature during each evaluation.
- Maps aromatic-flagged localized single and double bonds to DREIDING aromatic bonds
  without changing the bond orders stored by Molecules.
- Runs QEq separately for each model component using its formal-charge sum and keeps the
  resulting charges fixed during evaluation and minimization.
- Evaluates harmonic bonds, cosine angles, torsions, inversions, Lennard-Jones,
  electrostatic, and directional hydrogen-bond terms.
- Excludes 1-2 and 1-3 nonbonded pairs and includes full-strength 1-4 and inter-component
  pairs. Nonbonded work is all-pairs and therefore O(N^2).
- Preparation never sanitizes, adds hydrogens, or mutates the source model.

## Validation

- Unit tests compare Cartesian gradients with central finite differences and cover
  component charge isolation, exclusions, topology binding, singular geometry, and
  minimization integration.
- No external force-field golden corpus is currently accepted; `validated` remains false.

## Out Of Scope

- Macromolecules, periodic cells, cutoffs and neighbor lists, constraints, dynamics,
  charge updates during optimization, custom DREIDING parameters, and scientific accuracy
  claims beyond analytic regression coverage.

## Revision Notes

- v2: Bind to shared model identity and report structured evaluation geometry errors.
- v1: Add explicit DREIDING preparation and fixed-topology energy/gradient evaluation.
