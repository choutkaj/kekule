# Public API Facade

## Summary

Expose the architecture-defined public facade instead of a flat root namespace.

## Behavior/API

- Public modules are focused around `core`, `small`, `bio`, `smiles`, `molfile`, `sdf`, `perception`, and `canon`.
- The crate root no longer blanket re-exports implementation modules.
- The prelude is intentionally small and limited to common user-facing types.
- `SmallMolecule` owns small-molecule convenience methods and hides its raw graph field behind `graph()`, `graph_mut()`, and `into_graph()`.
- `MacroMolecule` hides its raw graph and hierarchy fields behind `graph()`, `graph_mut()`, `hierarchy()`, and `hierarchy_mut()`.
- SMILES, Molfile, SDF, perception, and canonicalization expose namespaced functions matching the staged architecture.

## Implementation Notes

- Existing algorithm and I/O internals remain available through focused facade modules rather than root aliases.
- `SmallMolecule::from_smiles` parses without sanitizing; `from_smiles_sanitized` and `smiles::read_sanitized_str` make sanitization explicit in the name.
- `graph_mut()` conservatively invalidates topology-derived perception state before handing out mutable graph access.
- Internal validation tooling uses the same public namespaces as user code.

## Validation

- Unit tests compile public happy-path and namespaced API examples.
- Workspace tests exercise the migrated validation tooling and existing chemistry/IO behavior through the new wrapper accessors.

## Out Of Scope

- Implementing new chemistry perception, stereochemistry, preparation, or macromolecule sanitization behavior.
- Keeping root-level compatibility aliases for the previous pre-release API.

## Revision Notes

- v1: Introduce architecture-aligned facade modules, a small prelude, and non-public wrapper graph fields.
