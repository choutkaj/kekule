# Public API Facade

## Summary

Expose the architecture-defined public facade instead of a flat root namespace.

## Behavior/API

- Public modules are focused around `core`, `units`, `small`, `bio`, `smiles`,
  `molfile`, `sdf`, `mmcif`, `perception`, `hydrogens`, `query`,
  `substructure`, `canon`, and `modeling`.
- The crate root no longer blanket re-exports implementation modules.
- The prelude is intentionally small and limited to common user-facing types.
- `SmallMolecule` owns small-molecule convenience methods and hides its raw graph field behind `graph()`, `graph_mut()`, and `into_graph()`.
- `MacroMolecule` exposes read-only graph/hierarchy access plus checked
  construction and transactional coordinated editing; completed values cannot
  be independently mutated into an invalid graph/hierarchy pair.
- `MacroMolecule` exposes direct hierarchy iterators, atom-site lookup, and
  read-only validation. The placeholder macro sanitization surface is absent.
- SMILES, Molfile, SDF, and mmCIF expose format-specific Documents and explicit
  interpretation results with reports/mappings; superseded direct reader APIs
  are absent.
- SMILES and Molfile retain simple default-bounded `parse_str` entry points and
  expose focused parse-options overloads; SDF and mmCIF accept their parse
  options directly.
- `mmcif::write` exposes explicit supported `Model` serialization with
  format-specific options and structured rejection errors.
- `Molecule` is one asserted entity and may have disconnected graph topology.
- `mmcif::interpret` returns a selected-coordinate `Model` plus report;
  `MolecularContents` and `Solvent` are removed.
- Expert perception functions live under focused modules such as `perception::rings`, `perception::aromaticity`, and `perception::valence`.
- Fixed-topology modelling types, potentials, and minimization live under `modeling` and are not added to the prelude.
- Explicit small-molecule hydrogen topology transforms live under `hydrogens`
  and as `SmallMolecule` convenience methods; they are not hidden in parsing or
  sanitization.
- Syntax-independent query graphs and bounded SMARTS translation live under
  `query`; matching lives under `substructure`, preserving one-way dependency
  on the query IR. Neither namespace is added to the prelude.

## Implementation Notes

- Existing algorithm and I/O internals remain available through focused facade modules rather than root aliases.
- `SmallMolecule::from_smiles` orchestrates parse/interpret without sanitizing;
  `from_smiles_sanitized` names the additional operation explicitly.
- `graph_mut()` itself is state-neutral; chemistry and topology mutators on the
  returned graph perform their own targeted invalidation, allowing perception
  operations to consume already-installed prerequisite state.
- Internal validation tooling uses the same public namespaces as user code.
- Invariant-bearing hierarchy, provenance, document, model, and structured
  error state is private behind accessors or checked constructors.
- Extensible public error enums are non-exhaustive. Deliberate value, options,
  and report payloads may retain direct public fields.
- Published crates start at `0.1.0`; breaking changes in the `0.x` line require
  a minor version increment.

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
- v10: Add the focused `hydrogens` namespace and `SmallMolecule`
  conveniences for transactional explicit/implicit hydrogen normalization.
- v11: Add focused `query` and `substructure` namespaces for syntax-neutral
  query graphs, bounded SMARTS parsing, and matching without expanding the prelude.
- v12: Add the foundational `mmcif::write` model-serialization surface without
  expanding the crate root or prelude.
- v13: Hard-break the modelling facade to `Model`/`ModelBuilder` and the full
  biomolecular hierarchy vocabulary to `Smcra*` names without compatibility
  aliases.
- v14: Add the focused `units` namespace and migrate coordinate and modelling
  boundaries to explicit quantities without expanding the prelude.
- v15: Establish the hard-break release facade: real format interpretation
  results/reports, checked macromolecule lifecycle, private invariant-bearing
  hierarchy/provenance/error state, non-exhaustive extensible errors, and the
  published `0.1.0` contract.
- v16: Expose configurable resource-bounded SMILES, Molfile, and SDF parsing
  without widening the crate root or prelude.
