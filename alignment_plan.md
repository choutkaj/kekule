# Rigid alignment implementation plan

> **Status:** Active implementation plan for the first rigid-alignment milestone.
> `ARCHITECTURE.md` remains authoritative. After the feature is implemented and
> its canonical feature contracts describe the shipped behavior, remove this
> plan rather than retaining it as historical documentation.

## Purpose

Implement a minimal, scientifically sound Kabsch/orthogonal-Procrustes alignment
capability in `molecular` so clients can align compatible molecular coordinate
snapshots without owning a duplicate geometry or chemistry implementation.

The immediate upstream consumer is MolStudio ensemble and conformer viewing. The
capability is nevertheless a general read-only structural analysis belonging in
`molecular`, not a viewer-specific helper.

The first milestone fits one moving coordinate snapshot onto one reference
snapshot over an explicit topology-bound atom selection and returns:

- the proper rigid transform that maps moving coordinates into the reference
  coordinate system;
- the post-fit root-mean-square deviation;
- basic diagnostics needed by downstream callers.

The operation derives a transform only. It never modifies either input model,
configuration, position array, topology, selection, or trajectory frame.

## Architectural placement

The public feature should live in a focused module:

```text
molecular::alignment
```

This is preferable to placing the complete API in `geometry`:

- `geometry` owns dependency-light mathematical value types such as `Point3`,
  `Vector3`, `Matrix3`, and `RigidTransform`;
- alignment depends on topology identity, topology-bound selections, coordinate
  snapshots, units, and structural views;
- the numerical kernel may use private geometry helpers, but the structural
  analysis boundary belongs above the primitive geometry layer.

Use the existing `geometry::RigidTransform`. Do not introduce another public
matrix, quaternion, pose, affine-transform, or rigid-transform type. Do not
expose a third-party linear-algebra type through the public API.

The primary canonical feature ID should be:

```text
algo.rigid-alignment
```

Expected direct feature dependencies are:

```text
model.topology
model.system
core.units
```

The implementation must also update `api.public-facade` and every other feature
contract actually affected by the final public surface.

## Required outcome

At completion, the repository must provide:

1. A public read-only Kabsch alignment operation over two `ModelView` values.
2. Correspondence defined by one exact-topology `AtomSelection`.
3. Uniform weighting as the default.
4. Explicit positive finite per-selected-atom weights as an advanced option.
5. An explicit periodic-coordinate policy.
6. A proper right-handed rigid transform mapping moving coordinates to reference
   coordinates.
7. Post-fit weighted RMSD in the model length unit.
8. Structured errors for incompatible topology, invalid selection, invalid
   weights, underdetermined geometry, unsupported periodic use, and numerical
   failure.
9. Deterministic, scale-aware numerical behavior in `f64`.
10. Focused unit, integration, and MolStudio consumer validation.

## Non-goals for the first milestone

The following are explicitly out of scope:

- graph isomorphism or atom-correspondence inference;
- alignment between independently constructed topology identities;
- symmetry-aware atom permutation;
- sequence or structural alignment;
- consensus or iterative mean-structure construction;
- automatic protein-backbone, ligand, chain, or residue selection;
- atomic-mass weighting until Molecular owns an authoritative mass source and a
  separate convenience contract;
- periodic imaging, molecule unwrapping, or minimum-image fitting;
- mutation of `Model`, `Configuration`, `Positions`, `Ensemble`, or trajectory
  coordinates;
- a batch ensemble-alignment container or cache;
- GPU acceleration, threading policy, or viewer-owned render transforms;
- non-rigid, flexible, affine, or scale fitting.

These may be added later as focused capabilities without weakening this
milestone's semantics.

## Target public API

Exact ergonomic details may be refined during implementation, but the public
semantics should remain equivalent to the following shape:

```rust
pub fn kabsch(
    moving: ModelView<'_>,
    reference: ModelView<'_>,
    selection: &AtomSelection,
) -> Result<RigidAlignment, AlignmentError>;

pub fn kabsch_with_options(
    moving: ModelView<'_>,
    reference: ModelView<'_>,
    selection: &AtomSelection,
    options: KabschOptions<'_>,
) -> Result<RigidAlignment, AlignmentError>;

#[derive(Debug, Clone, Copy)]
pub struct KabschOptions<'a> {
    pub weighting: AlignmentWeighting<'a>,
    pub periodic_policy: PeriodicAlignmentPolicy,
}

#[derive(Debug, Clone, Copy, Default)]
pub enum AlignmentWeighting<'a> {
    #[default]
    Uniform,
    Explicit(&'a [f64]),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
#[non_exhaustive]
pub enum PeriodicAlignmentPolicy {
    #[default]
    RejectPeriodic,
    UseStoredCoordinates,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RigidAlignment {
    transform: RigidTransform,
    rmsd: Quantity<f64>,
    selected_atom_count: usize,
}
```

`RigidAlignment` should expose read-only getters. Its fields should remain
private so future diagnostics can be added without allowing inconsistent manual
construction.

`kabsch` is the convenience form using default options:

```text
uniform weights + reject periodic coordinates
```

Do not add this API to the broad prelude merely for convenience. The focused
`alignment` module is the canonical public path.

## Transform direction and units

The most important directional invariant is:

```text
aligned_point = result.transform().transform_point(moving_point)
```

The result maps **moving -> reference**. This direction must be stated in
rustdoc, feature documentation, and asymmetric regression tests. Do not expose
an ambiguously named transform whose direction callers must infer.

Coordinates are already stored internally in `MODEL_LENGTH_UNIT`. The alignment
kernel operates on those canonical raw values. The returned `RigidTransform`
therefore uses the existing geometry convention, while RMSD is returned as a
`Quantity<f64>` in `MODEL_LENGTH_UNIT`.

The fitted transform is derived from the selected atoms but may be applied by a
caller to every point in the moving structure. Molecular does not need to
materialize transformed full-coordinate arrays for this milestone.

## Topology and correspondence contract

The moving view, reference view, and selection must all belong to the same exact
`Topology` identity.

Equal atom counts and `Topology::same_layout` are not sufficient. This first
milestone is designed for models, ensemble members, conformer adaptations, and
trajectory frames that already share one exact topology.

Correspondence is positional in the selection's validated dense-index order:

```text
moving[position at selected index i] <-> reference[position at selected index i]
```

`AtomSelection` already guarantees sorted unique dense indices. The alignment
implementation must call its compatibility check rather than reaching into or
reconstructing selection identity.

Alignment across different topology identities requires a future explicit
validated correspondence or topology mapping and must not be guessed here.

## Weighting contract

Uniform weighting assigns equal weight to every selected atom.

For `AlignmentWeighting::Explicit`:

- the slice length must equal `selection.indices().len()`;
- weights correspond to the selection's dense-index order, not complete topology
  order;
- every weight must be finite and strictly positive;
- multiplying every weight by one positive constant must not alter the fitted
  transform or RMSD;
- the implementation may normalize weights internally but must not mutate the
  caller's slice.

Zero weights are rejected in this milestone rather than silently changing the
active correspondence set. Callers that want to exclude atoms should construct
a different selection.

## Periodic-coordinate contract

Periodic coordinates are ambiguous unless molecules have been imaged or
unwrapped according to an explicit scientific policy.

The default `RejectPeriodic` policy returns a structured error if either input
configuration contains a periodic cell.

`UseStoredCoordinates` means exactly:

- ignore the cell for fitting;
- treat the stored Cartesian coordinates as already prepared by the caller;
- perform no wrapping, unwrapping, centering, minimum-image correction, or
  molecule reconstruction.

This explicit escape hatch supports known-prepared data without pretending that
raw periodic alignment is generally well defined.

## Mathematical contract

For moving points `x_i`, reference points `y_i`, and positive weights `w_i`, the
implementation must minimize:

```text
sum_i w_i * ||R x_i + t - y_i||^2
```

subject to:

```text
R^T R = I
det(R) = +1
```

The required steps are mathematically equivalent to:

1. compute weighted moving and reference centroids;
2. center both selected point sets;
3. accumulate the weighted cross-covariance matrix;
4. solve the proper orthogonal-Procrustes/Kabsch problem;
5. correct or exclude reflection so the result is right-handed;
6. compute `t = reference_centroid - R * moving_centroid`;
7. construct and validate the existing `RigidTransform`;
8. compute post-fit weighted RMSD from the final transform.

The weighted RMSD is:

```text
sqrt(sum_i w_i * ||R x_i + t - y_i||^2 / sum_i w_i)
```

A mirrored structure must never be reported as a zero-RMSD rigid fit through an
improper rotation.

## Determinacy and degeneracy

The first milestone requires at least three selected atoms with positive
weights.

Point count alone is not sufficient. Both centered point sets must span at least
two geometric dimensions under a scale-aware numerical rank test.

Required behavior:

- three or more planar, non-collinear points are valid;
- one point is underdetermined;
- two points are underdetermined;
- coincident points are underdetermined;
- collinear points are underdetermined;
- near-degenerate inputs must be classified using a documented scale-relative
  tolerance, not one absolute coordinate threshold.

Do not reject all planar selections. Protein backbones, aromatic groups, and
other scientifically useful selections may be approximately or exactly planar.

The implementation should return a structured degeneracy error rather than an
arbitrary rotation for rank-deficient input.

## Numerical implementation requirements

The kernel must use `f64` and be deterministic for a fixed input.

Preferred implementation properties:

- O(n) time in selected atom count;
- O(1) additional numerical storage beyond fixed-size matrices and accumulators;
- no cloning of complete coordinate arrays;
- no owned `Model` construction;
- stable centroid and covariance accumulation for coordinates far from the
  origin;
- scale-aware rank and residual checks;
- a final proper orthonormal rotation accepted by `RigidTransform::new`;
- no `unsafe` code.

A private fixed-size SVD, polar-decomposition method, or mathematically
equivalent quaternion eigensolver is acceptable. The implementation choice must
be documented and tested. Avoid a broad linear-algebra dependency unless a
focused dependency is clearly justified, minimally configured, and does not
leak into the public API.

Do not weaken the global `RigidTransform` validation merely to accommodate a
numerically imprecise solver. Correct or re-orthonormalize the private numerical
result when necessary, then validate through the existing constructor.

## Error model

Use a focused public `#[non_exhaustive]` error type. Exact variant names may be
refined, but callers must be able to distinguish at least:

```text
AlignmentError
- moving/reference topology identity mismatch
- selection topology identity mismatch
- insufficient selected points
- degenerate selected geometry
- explicit weight-count mismatch
- non-finite explicit weight
- non-positive explicit weight
- periodic coordinates rejected by policy
- numerical solution failure
- invalid resulting rigid transform
```

Errors should expose relevant counts or selected-weight positions where useful.
Do not collapse failures into `String`, `PositionError`, or one generic
`NumericalFailure` when a stable caller-relevant category is known.

Because public positions are already finite and topology-bound, the kernel may
rely on those constructor invariants rather than duplicating unrelated position
validation.

## Input immutability

Alignment is an analysis. Success and every error path must leave unchanged:

- moving and reference topology identity;
- moving and reference coordinates;
- periodic cells;
- observations and properties;
- the atom selection;
- any ensemble or trajectory container from which the views were borrowed.

Downstream applications apply the returned transform as derived view or render
state. Canonical molecular coordinates are not rewritten by alignment.

## Required tests

Add focused tests covering at least the following cases.

### Correct fits

- identical coordinates produce identity transform and zero RMSD;
- translation-only coordinates recover the correct moving-to-reference
  translation;
- a known asymmetric non-axis-aligned rotation plus translation is recovered;
- applying the result to every selected moving point reproduces the reference;
- a subset selection determines the transform while unselected atoms do not
  affect the fit;
- planar non-collinear points succeed;
- noisy coordinates produce the expected nonzero post-fit RMSD and improve over
  the unaligned RMSD;
- explicit weights change the optimum in a deliberately asymmetric fixture;
- uniformly rescaling explicit weights leaves the result unchanged.

### Proper rotation and direction

- mirrored coordinates yield a proper transform with determinant `+1` and
  nonzero RMSD;
- an asymmetric fixture proves the API maps moving -> reference rather than the
  inverse direction;
- returned rotations satisfy the existing `RigidTransform` validation.

### Validation failures

- empty selection;
- one selected atom;
- two selected atoms;
- coincident selected points;
- collinear selected points;
- deliberately near-degenerate geometry on both sides of the documented rank
  threshold;
- different moving/reference topology identities, including equal-layout
  topologies;
- stale or foreign atom selection;
- explicit weight-count mismatch;
- zero, negative, NaN, and infinite weights;
- periodic input under default policy.

### Periodic and stability behavior

- `UseStoredCoordinates` permits periodic configurations and ignores their cell;
- no imaging or cell-dependent coordinate modification occurs;
- large absolute coordinate offsets with small internal geometry remain
  accurate;
- neither input model nor selection changes on success or failure.

### Public contract

- rustdoc examples compile;
- result RMSD uses `MODEL_LENGTH_UNIT`;
- the focused module is accessible through the public facade without an
  unnecessary prelude expansion.

## Optional external validation

An informational benchmark against an independent NumPy/SciPy, RDKit, or other
well-defined Kabsch reference is useful, especially for randomized proper
transforms, weighted cases, reflections, planar inputs, and large coordinate
offsets.

Such a benchmark must follow the repository's provenance-pinned external
benchmark rules. It is not a runtime dependency and should not become a release
or CI requirement merely to mark the feature supported.

## Staged implementation

## Stage 0: Contracts and baseline

- Read `AGENTS.md`, `ARCHITECTURE.md`, this plan, and the affected topology,
  model, units, and public-facade feature contracts.
- Create `features/algo.rigid-alignment/feature.toml` and `feature.md` with honest
  initial status.
- Inventory existing geometry, `RigidTransform`, `ModelView`, `AtomSelection`,
  and MolStudio capability-error call sites.
- Record the four locked Molecular baseline gates.

Exit criterion: public semantics and affected feature IDs are explicit before
numerical implementation begins.

## Stage 1: Public types and validation boundary

- Add the focused `alignment` module.
- Add options, weighting, periodic policy, result, getters, and structured
  errors.
- Implement exact topology/selection compatibility and option validation.
- Add compile tests and validation-error regressions.

Exit criterion: the final public API shape compiles and rejects invalid inputs
without mutating state.

## Stage 2: Numerical kernel

- Implement weighted centroids and covariance accumulation.
- Implement the private proper orthogonal-Procrustes solver.
- Add scale-aware rank detection that accepts planar non-collinear input.
- Construct the moving-to-reference `RigidTransform`.
- Compute post-fit weighted RMSD.
- Add analytic correctness, reflection, degeneracy, and numerical-stability
  tests.

Exit criterion: all focused alignment tests pass deterministically.

## Stage 3: Documentation and feature integration

- Complete rustdoc with transform direction, weighting, periodic semantics,
  degeneracy, units, and examples.
- Mark the feature supported only after behavior and tests are complete.
- Update `api.public-facade`, `ARCHITECTURE.md`, `CHANGELOG.md`, and other actually
  affected feature contracts.
- Regenerate feature dashboard artifacts through repository tooling; do not
  hand-edit generated HTML.

Exit criterion: canonical architecture and feature documentation describe the
implemented behavior without relying on this plan.

## Stage 4: Consumer validation and cleanup

- Validate the sibling MolStudio checkout through the untracked local Cargo
  patch.
- Replace or exercise the alignment capability boundary sufficiently to prove
  the public API satisfies the blocked same-topology ensemble/conformer use
  case; do not commit a path dependency or local lockfile rewrite.
- Run all required locked Molecular gates and the relevant locked MolStudio
  consumer gates.
- Remove `alignment_plan.md` in the implementation PR once the canonical feature
  contracts and architecture documentation are complete.

Exit criterion: the feature is available at a reproducibly pinned public Git
revision and MolStudio can consume it without a local substitute.

## Completion criteria

The milestone is complete only when:

- the public API and transform direction are unambiguous;
- proper weighted Kabsch behavior is implemented and tested;
- planar non-collinear selections work;
- underdetermined geometry fails structurally;
- exact topology identity and explicit periodic policy are enforced;
- inputs remain immutable;
- all four locked Molecular gates pass;
- the sibling MolStudio consumer check passes against the candidate revision;
- feature metadata and generated artifacts are current;
- the implementation is published at a public Git revision;
- this temporary plan has been removed in favor of canonical architecture,
  rustdoc, and feature contracts.
