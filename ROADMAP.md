# Roadmap

The canonical project status is the feature registry under `features/`. Its
statuses and dependency graph are rendered procedurally into
`features/DASHBOARD.html`; this document records release direction rather than
duplicating that inventory.

## 0.1 release line

The first release establishes the graph kernel, typed small- and
macromolecules, fixed-coordinate `Model` workflows, staged structure I/O,
sanitization and perception, stereochemistry, bounded query matching,
macromolecular hierarchy and secondary-structure analysis, and the DREIDING
adapter. Supported features and their current validation evidence are listed in
the generated dashboard.

## Next tracked capabilities

Two feature contracts are currently reserved with `planned` status:

- `descriptor.molecular`: explicit-policy molecular formula, mass, and related
  descriptors.
- `fp.morgan`: a defined-shape Morgan-style circular fingerprint with explicit
  perception dependencies.

Additional work should begin as a feature contract with explicit dependencies,
resource limits, and validation requirements before it is treated as a release
commitment.
