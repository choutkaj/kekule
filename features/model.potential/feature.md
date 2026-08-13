# Model Potential Interface

## Summary

Provide a minimal object-safe energy-and-gradient contract over borrowed
topology-plus-model-state views and a transparent caller-parameterized
harmonic bond potential.

## Behavior/API

- Exposes `Potential`, `PotentialEvaluation`, `PotentialError`, and
  `PotentialGeometryError` under `kekule::modeling::potential`; shared
  `Vector3` lives in `kekule::geometry`.
- Requires one finite Cartesian gradient vector per model atom and rejects non-finite energy or gradients.
- Exposes `HarmonicBondParameter` and `HarmonicBondPotential` for explicit
  `InstanceBondId` parameters; atom errors use `InstanceAtomId` and gradients
  remain dense in `TopologyAtomIndex` order.
- Accepts explicit compatible quantities for harmonic parameters and potential
  outputs, then converts once to the modelling kernel's declared length,
  energy, gradient, and force-constant units.
- Distinguishes incompatible topology allocations, coordinate singularities,
  unsupported periodic configurations, malformed outputs, and backend failures.
- `ModelView` is a common transport for topology plus direct model state, not
  a promise that every potential supports every field. Each implementation
  documents its periodic-cell capability.
- `HarmonicBondPotential` is nonperiodic and returns
  `PotentialError::UnsupportedPeriodicCell` whenever an evaluated view carries
  a periodic cell.

## Implementation Notes

- `Potential::evaluate(ModelView)` takes `&mut self` so implementations may
  retain caches while remaining object-safe.
- Prepared potentials retain one `Arc<Topology>` and remain compatible
  with supported model and ensemble-member views in `kekule`, plus
  `kekule-traj` frame and buffer views sharing that topology.
- Harmonic terms use `0.5 * k * (r - r0)^2` and validate positive finite parameters, unique bond terms, and the topology observed at construction.
- Coincident bonded atoms return a structured coordinate-geometry failure because a nonzero-rest-length harmonic gradient has no defined Cartesian direction there.
- The built-in potential performs no parameter inference and contains no angle, torsion, or nonbonded interactions.

## Tests

- Unit tests compare analytic harmonic gradients against central finite differences in arbitrary orientations.
- Tests cover invalid bonds, duplicate or invalid parameters, malformed evaluations, topology mismatch, additive terms, and coincident atoms.
- Periodic-policy tests place bonded atoms on opposite box faces and verify
  explicit rejection through core model and ensemble views plus cross-crate
  `kekule-traj` frame and frame-buffer views, while the same nonperiodic
  coordinates remain evaluable.
- Reference molecular goldens are not currently defined for this analytic
  infrastructure, so no external parity result is recorded.

## Out Of Scope

- Automatic force-field typing, prepared backend lifecycles, energy-only evaluation, QM backends, nonbonded interactions, and runtime reference-tool dependencies.

## Revision Notes

- v1: Add the potential contract, validated evaluation container, and explicit harmonic bond potential.
- v2: Qualify topology references by molecule instance and include ownership in
  mismatch detection.
- v3: Bind prepared potentials to shared model-definition identity and add
  structured coordinate-geometry and backend failures.
- v4: Migrate potential signatures to the renamed canonical `Model` API.
- v5: Make energies, Cartesian gradients, harmonic lengths, and force constants
  quantity-valued instead of relying on documentation-only units.
- v6: Evaluate borrowed `ModelView` values, bind preparation and gradients to
  one exact topology/dense order, and move vector geometry to the common
  `geometry` module.
- v7: Add a structured unsupported-periodic-cell error and make the built-in
  harmonic potential's nonperiodic capability explicit across every structural
  view.
- v8: Move trajectory-owned view regressions to `kekule-traj` while preserving
  the dependency-light `ModelView` potential contract in `kekule`.
- v9: Use retained `Arc<Topology>` values and
  pointer-compatible evaluation across model, ensemble, trajectory, and buffer
  views.
