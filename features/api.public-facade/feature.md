# Public API Facade

## Summary

Expose the architecture-defined public facade instead of a flat root namespace.

## Behavior/API

- Public modules are focused around `core`, `small`, `bio`, `smiles`, `molfile`, `sdf`, `perception`, and `canon`.
- The crate root no longer blanket re-exports implementation modules.
- The prelude is intentionally small and limited to common user-facing types.
- `SmallMolecule` owns small-molecule convenience methods and hides its raw graph field behind `graph()`, `graph_mut()`, and `into_graph()`.
- `MacroMolecule` hides its raw graph and hierarchy fields behind `graph()`, `graph_mut()`, `hierarchy()`, and `hierarchy_mut()`.
- `MacroMolecule` exposes direct hierarchy iterators, atom-site lookup, and separate macro validation/sanitization APIs.
- SMILES, Molfile, SDF, perception, and canonicalization expose namespaced functions matching the staged architecture.
- Expert perception functions live under focused modules such as `perception::rings`, `perception::aromaticity`, and `perception::valence`.

## Implementation Notes

- Existing algorithm and I/O internals remain available through focused facade modules rather than root aliases.
- `SmallMolecule::from_smiles` parses without sanitizing; `from_smiles_sanitized` and `smiles::read_sanitized_str` make sanitization explicit in the name.
- `graph_mut()` conservatively invalidates topology-derived perception state before handing out mutable graph access.
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
