# DREIDING Force-Field Adapter

## Summary

Provide an explicit topology-bound adapter that prepares DREIDING atom types,
fixed QEq charges, bonded terms, nonbonded terms, and complete Cartesian
gradients once for supported nonperiodic structural views.

## Behavior/API

- Exposes `DreidingPotential`, `DreidingPrepareOptions`, `QeqGrouping`, and
  `DreidingPrepareError` from the `kekule-potentials` companion package under
  the focused `kekule_potentials::dreiding` namespace. The companion crate
  does not re-export DREIDING types from its root.
- The package's default `dreiding` feature owns the optional `dreid-forge` and
  `dreid-kernel` dependencies, allowing future independent implementations to
  be selected without forcing DREIDING dependencies.
- Prepares with an explicit `&Arc<Topology>`, reference `ModelView`, and QEq grouping
  policy, then implements the core borrowed-view `Potential` contract.
- The adapter is nonperiodic. Preparation rejects a reference view carrying a
  periodic cell with `DreidingPrepareError::UnsupportedPeriodicCell`, and
  evaluation rejects periodic views with
  `PotentialError::UnsupportedPeriodicCell`.
- Binds preparation to one shared topology allocation, accepting models, ensemble
  members, and trajectory frames sharing it while rejecting independently
  constructed equal topology.
- Exposes read-only per-atom type diagnostics and quantity-valued partial charges.
- Rejects unresolved implicit-hydrogen state; every atom must carry an explicit zero
  implicit-hydrogen count or a no-implicit-hydrogens assertion.
- Consumes the view's declared coordinate quantity and returns explicit
  kJ/mol energy and kJ/mol/angstrom gradient quantities.

## Implementation Notes

- Uses pinned `dreid-forge` and matching `dreid-kernel` releases; upstream types do not
  cross the adapter's public API.
- Keeps the dependency-light `Potential` evaluation contract in `kekule`; the
  companion package owns concrete implementations and namespaces them by
  preparation model.
- Retains the shared `Arc<Topology>` rather than rebuilding an
  adapter-specific signature during each evaluation.
- Maps aromatic-flagged localized single and double bonds to DREIDING aromatic bonds
  without changing the bond orders stored by Kekule.
- Makes QEq grouping explicit as whole topology, molecule instances, or actual
  connected components; molecule-instance grouping is the default. Charges
  remain fixed during evaluation and minimization.
- Evaluates harmonic bonds, cosine angles, torsions, inversions, Lennard-Jones,
  electrostatic, and directional hydrogen-bond terms. Eligible Small and Macro
  instances use the same chemistry requirements.
- Excludes 1-2 and 1-3 nonbonded pairs and includes full-strength 1-4 and inter-instance
  pairs. Nonbonded work is all-pairs and therefore O(N^2).
- Preparation never normalizes, perceives, adds hydrogens, or mutates topology or the
  reference model.
- No periodic cell is ignored and no orthorhombic-only minimum-image shortcut
  is applied.

## Tests

- Unit tests compare Cartesian gradients with central finite differences and cover
  molecule-instance charge isolation, exclusions, topology binding, singular geometry, and
  minimization integration.
- Tests prepare once and evaluate model, ensemble-member, and trajectory-frame
  views, reject independently allocated equal-layout topologies, and exercise every QEq grouping
  policy.
- Periodic-policy tests use atoms on opposite box faces and verify structured
  preparation failure plus consistent model, ensemble, trajectory-frame, and
  frame-buffer evaluation failure; the same coordinates without a cell remain
  evaluable.
- No external force-field golden corpus is currently accepted, so no parity
  result is recorded.

## Out Of Scope

- Periodic cells, cutoffs and neighbor lists, constraints, dynamics,
  charge updates during optimization, custom DREIDING parameters, and scientific accuracy
  claims beyond analytic regression coverage.

## Revision Notes

- v1: Add explicit DREIDING preparation and fixed-topology energy/gradient evaluation.
- v2: Migrate to molecule-qualified IDs, per-instance QEq, mixed Small/Macro
  models, and instance-boundary topology signatures.
- v3: Replace adapter-specific topology signatures with shared model-definition
  identity and report structured evaluation geometry errors.
- v4: Build adjacency and nonbonded exclusions through dense model indexes so
  repeated instances and tombstoned molecule-local atom IDs remain isolated.
- v5: Migrate preparation and evaluation signatures to the renamed canonical
  `Model` API.
- v6: Integrate explicit coordinate, energy, gradient, and charge quantities at
  the adapter boundary while retaining raw numeric inner kernels.
- v7: Bind preparation to the exact topology used at preparation, evaluate `ModelView`, make
  reference-geometry use explicit, and distinguish whole-topology,
  molecule-instance, and connected-component QEq grouping.
- v8: Declare DREIDING nonperiodic and reject periodic reference and evaluation
  views structurally instead of silently applying direct Cartesian geometry.
- v9: Hard-rename the adapter package and Rust import root to
  `kekule-dreiding` and `kekule_dreiding`.
- v10: Generalize the companion package to `kekule-potentials` and move the
  DREIDING public surface under `kekule_potentials::dreiding` without changing
  preparation or evaluation behavior.
- v11: Set the publishable `kekule-potentials` package and its Kekule companion
  dependencies to the shared initial `0.1.0` release line.
- v12: Retain `Arc<Topology>` during preparation and use shared-allocation
  compatibility without a separate token.
