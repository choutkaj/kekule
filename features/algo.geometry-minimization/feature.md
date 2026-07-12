# Geometry Minimization

## Summary

Minimize a fixed-topology molecular model without mutating the input, using deterministic normalized steepest descent and Armijo backtracking.

## Behavior/API

- Exposes `modeling::minimize`, `MinimizeOptions`, `MinimizationResult`, `MinimizationStatus`, and `MinimizationError`.
- Reports initial/final energy, final maximum-gradient norm, accepted iterations, potential evaluations, and terminal status.
- Treats convergence, maximum iterations, and line-search stalling as explicit non-error statuses.
- Rejects invalid optimizer options and propagates structured potential or position errors.
- Backtracks trial evaluations that report coordinate-only invalid geometry while propagating incompatible-model and backend failures.

## Implementation Notes

- Search directions are negative gradients normalized so the largest atom-vector norm is one; the line-search step therefore expresses maximum atom displacement in angstroms.
- Armijo backtracking guarantees accepted steps provide sufficient energy decrease.
- Every attempted potential call, including rejected invalid-geometry trials, contributes to the evaluation count.
- Default limits are 1000 iterations, `1e-4` kJ/mol/angstrom gradient tolerance, `0.1` angstrom initial step, `1e-8` angstrom minimum step, factor `0.5`, Armijo coefficient `1e-4`, and 24 backtracks.
- The input model remains unchanged on success, non-convergence, or error.

## Validation

- Unit tests minimize distorted one- and multi-component harmonic systems and verify energy decrease, equilibrium distances, convergence reports, source immutability, iteration limits, and line-search stalling.
- The implementation uses focused analytic regressions rather than external molecular fixtures; `validated` remains false until accepted harness evidence exists.

## Out Of Scope

- L-BFGS, FIRE, constraints, frozen atoms, topology changes, dynamics, trajectories, and production force-field validation.

## Revision Notes

- v2: Backtrack recoverable coordinate singularities and preserve fatal evaluation errors.
- v1: Add normalized steepest descent with Armijo backtracking and structured convergence results.
