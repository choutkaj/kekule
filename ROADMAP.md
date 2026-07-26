# Roadmap

The canonical project status is the feature registry under `features/`. Its
statuses and dependency graph are rendered procedurally into
`features/DASHBOARD.html`; this document records release direction rather than
duplicating that inventory.

## 0.2 release line

The 0.2 line establishes immutable coordinate-free `Topology`, explicit
definition reuse and exact identity, topology-bound configurations and borrowed
views, finite ensembles, streaming-first fixed-topology trajectories,
coordinate-independent SMCRA hierarchy, separate single-/multi-model mmCIF
interpretation, and topology-bound prepared potentials. Supported features and
optional benchmark availability/observations are listed separately in the
generated dashboard.

The complete-instance topology-transformation milestone is implemented:
callers can retain or remove whole molecule instances with complete checked
lineage and explicitly remap models, selections, ensembles, frames,
trajectories, and reusable buffers without weakening exact topology identity.
The next topology-editing direction is explicit topology composition, followed
by instance-definition replacement and definition-edit scope, each as a
separate feature contract.

## Next tracked capabilities

Feature contracts currently reserved with `planned` status include:

- `descriptor.molecular`: explicit-policy molecular formula, mass, and related
  descriptors.
- `fp.morgan`: a defined-shape Morgan-style circular fingerprint with explicit
  perception dependencies.

Additional work should begin as a feature contract with explicit dependencies,
resource limits, and tests before it is treated as a release commitment.
External-reference benchmarks may be added where useful, but they are not a
release-status requirement.
