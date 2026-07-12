# Model Potential Interface

## Summary

Provide a minimal object-safe energy-and-gradient contract for fixed-topology molecular models and a transparent caller-parameterized harmonic bond potential.

## Behavior/API

- Exposes `Potential`, `PotentialEvaluation`, `PotentialError`, `PotentialGeometryError`, and `Vector3` under `molecules::modeling::potential`.
- Requires one finite Cartesian gradient vector per model atom and rejects non-finite energy or gradients.
- Exposes `HarmonicBondParameter` and `HarmonicBondPotential` for explicit model `BondId` parameters.
- Uses angstroms, kJ/mol energies, kJ/mol/angstrom gradients, and kJ/mol/angstrom-squared force constants.
- Distinguishes incompatible models, coordinate singularities, malformed outputs, and backend failures.

## Implementation Notes

- `Potential::evaluate` takes `&mut self` so implementations may retain caches while remaining object-safe.
- Prepared potentials bind to the model's opaque definition key and remain compatible with coordinate-modified clones.
- Harmonic terms use `0.5 * k * (r - r0)^2` and validate positive finite parameters, unique bond terms, and the topology observed at construction.
- Coincident bonded atoms return a structured coordinate-geometry failure because a nonzero-rest-length harmonic gradient has no defined Cartesian direction there.
- The built-in potential performs no parameter inference and contains no angle, torsion, or nonbonded interactions.

## Validation

- Unit tests compare analytic harmonic gradients against central finite differences in arbitrary orientations.
- Tests cover invalid bonds, duplicate or invalid parameters, malformed evaluations, topology mismatch, additive terms, and coincident atoms.
- Reference molecular goldens are not required for this analytic infrastructure; `validated` remains false until accepted harness evidence exists.

## Out Of Scope

- Automatic force-field typing, prepared backend lifecycles, energy-only evaluation, QM backends, nonbonded interactions, and runtime reference-tool dependencies.

## Revision Notes

- v2: Add shared model-definition binding and structured evaluation-failure categories.
- v1: Add the potential contract, validated evaluation container, and explicit harmonic bond potential.
