# Public API Facade

## Summary

Expose the architecture-defined public facade instead of a flat root namespace.

## Behavior/API

- Public modules are focused around `core`, `small`, `bio`, `smiles`, `molfile`, `sdf`, `mmcif`, `perception`, `canon`, and `modeling`.
- The crate root no longer blanket re-exports implementation modules.
- The prelude is intentionally small and limited to common user-facing types.
- `SmallMolecule` owns small-molecule convenience methods and hides its raw graph field behind `graph()`, `graph_mut()`, and `into_graph()`.
- `MacroMolecule` hides its raw graph and hierarchy fields behind `graph()`, `graph_mut()`, `hierarchy()`, and `hierarchy_mut()`.
- `MacroMolecule` exposes direct hierarchy iterators, atom-site lookup, and separate macro validation/sanitization APIs.
- SMILES, Molfile, SDF, and mmCIF expose format-specific Documents and explicit
  interpretation; superseded direct reader APIs are absent.
- `Molecule` is one asserted entity and may have disconnected graph topology.
- `mmcif::interpret` returns a selected-coordinate `MolecularModel` plus report;
  `MolecularContents` and `Solvent` are removed.
- Expert perception functions live under focused modules such as `perception::rings`, `perception::aromaticity`, and `perception::valence`.
- Fixed-topology modelling types, potentials, and minimization live under `modeling` and are not added to the prelude.

## Implementation Notes

- Existing algorithm and I/O internals remain available through focused facade modules rather than root aliases.
- `SmallMolecule::from_smiles` orchestrates parse/interpret without sanitizing;
  `from_smiles_sanitized` names the additional operation explicitly.
- `graph_mut()` itself is state-neutral; chemistry and topology mutators on the
  returned graph perform their own targeted invalidation, allowing perception
  operations to consume already-installed prerequisite state.
- Macro validation is read-only; macro sanitization is conservative and rejects unsupported preparation-like options instead of silently guessing.
- Internal validation tooling uses the same public namespaces as user code.

## Validation

- External integration tests compile public happy-path, namespaced, low-level graph, and macro-molecule API examples as downstream user code.
- Workspace tests exercise the migrated validation tooling and existing chemistry/IO behavior through the new wrapper accessors.

## Out Of Scope

- Implementing new chemistry perception, stereochemistry, preparation, or invasive macromolecule sanitization behavior.
- Keeping root-level compatibility aliases for the previous pre-release API.

## Revision Notes

- v1: Introduce architecture-aligned facade modules, a small prelude, and non-public wrapper graph fields.
- v2: Move expert perception APIs under focused facade modules and add separate macro validation/sanitization surface.
- v3: Add downstream-style integration tests for the architecture-level public API.
- v4: Add the focused SmallMolecule modelling, potential, and minimization namespace without expanding the prelude.
- v5: Add staged mmCIF document interpretation and molecular-content containers without expanding the prelude.
- v6: Hard-break the historical direct mmCIF reader and all compatibility re-exports.
- v7: Molecule-first hard break: format Documents, private `PerceptionState`,
  instance-based `ModelTopology`, mmCIF model output, and deletion of all
  superseded readers/components/content containers.
- v8: Make wrapper mutable graph access state-neutral and rely on concrete graph
  mutators for invalidation, preventing perception prerequisites from being
  erased before stereo and CIP operations.
- v9: Expose opaque shared model-definition identity and instance-qualified
  structured potential failures through the focused modelling namespace.
