# Public API Facade

## Summary

Expose the architecture-defined public facade instead of a flat root namespace.

## Behavior/API

- Public modules are focused around `core`, `units`, `small`, `bio`, `smiles`,
  `molfile`, `sdf`, `mmcif`, `perception`, `hydrogens`, `query`,
  `substructure`, `canon`, `descriptors`, `geometry`, `topology`, `structure`,
  `trajectory`, and `modeling`.
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
- `Molecule::add_atom` is fallible. Graph, conformer, stereo, SMCRA hierarchy,
  and topology insertion surfaces report focused fixed-width identifier
  capacity errors before mutation rather than truncating or panicking.
- `Molecule::formal_charge` exposes the asserted live-atom charge aggregate as
  an `i64` without hiding sanitization or perception.
- `mmcif::interpret` returns a selected-coordinate `Model` plus report;
  `MolecularContents` and `Solvent` are removed.
- `mmcif::interpret_ensemble` is a separate shared-topology multi-model path
  that rejects inconsistent atom identity or topology.
- `MmcifInterpretation::into_model` consumes an interpretation when callers do
  not need to retain its report.
- Expert perception functions live under focused modules such as `perception::rings`, `perception::aromaticity`, and `perception::valence`.
- Immutable system structure, coordinate containers, and frames live under
  `topology`, `structure`, and `trajectory`; potentials and minimization remain
  under `modeling`. These focused types are not added to the prelude.
- The topology facade names complete static equality `Topology::same_layout`
  and restricts `TopologyMapping::between_identical_layouts` to identity maps
  over that exact layout. Checked explicit mappings and topology-edit results
  expose structured consistency failures.
- Immutable whole-instance edits live under the focused
  `topology::transform` namespace. Structure and trajectory containers expose
  explicit remapping methods and typed errors without expanding the crate root
  or prelude.
- `Model::instance_to_conformer` provides an explicit transactional path from
  instance-qualified model positions back to a compatible local conformer.
- Explicit small-molecule hydrogen topology transforms live under `hydrogens`
  and as `SmallMolecule` convenience methods; they are not hidden in parsing or
  sanitization.
- Syntax-independent query graphs and bounded SMARTS translation live under
  `query`; matching lives under `substructure`, preserving one-way dependency
  on the query IR. Neither namespace is added to the prelude.
- Read-only molecular formula and mass calculation lives under `descriptors`,
  requires an explicit hydrogen-count policy, and is not added to the prelude.

## Implementation Notes

- Existing algorithm and I/O internals remain available through focused facade modules rather than root aliases.
- `SmallMolecule::from_smiles` orchestrates parse/interpret without sanitizing;
  `from_smiles_sanitized` names the additional operation explicitly.
- `graph_mut()` itself is state-neutral; chemistry and topology mutators on the
  returned graph perform their own targeted invalidation, allowing perception
  operations to consume already-installed prerequisite state.
- Internal benchmark tooling uses the same public namespaces as user code.
- Invariant-bearing hierarchy, provenance, document, model, and structured
  error state is private behind accessors or checked constructors.
- Extensible public error enums are non-exhaustive. Deliberate value, options,
  and report payloads may retain direct public fields.
- The topology-centered hard break is released as `0.2.0`; later breaking
  changes in the `0.x` line likewise require a minor version increment.

## Tests

- External integration tests compile public happy-path, namespaced, low-level
  graph, macro-molecule, topology/configuration/model, mmCIF, and borrowed-view
  examples as downstream user code.
- Public topology tests compile the exact-identity, same-layout, checked
  identity-mapping, and checked edit-result surface.
- Downstream transformation tests compile instance retain/remove, stable
  mapping traversal, model and selection remapping, ensemble remapping, owned
  trajectory remapping, and reusable target-buffer remapping.
- Workspace tests exercise the benchmark tooling and existing chemistry/IO behavior through the new wrapper accessors.

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
  instance-qualified system structure, mmCIF model output, and deletion of all
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
- v17: Add the core formal-charge aggregate, transactional model-instance
  conformer export, and concise basic versus modeling public examples.
- v18: Add consuming `MmcifInterpretation::into_model` access through the
  focused mmCIF facade.
- v19: Add the focused `descriptors` facade for explicit-policy molecular
  formula and mass calculation without expanding the prelude.
- v20: Add the 0.2 `geometry`, `topology`, `structure`, and `trajectory`
  modules; remove obsolete model-owned topology concepts; expose separate
  mmCIF ensemble interpretation; and migrate DSSP/potentials to borrowed views.
- v21: Replace the misleading topology structural-equivalence names with
  complete static-layout semantics and expose checked mapping/edit-result
  construction.
- v22: Make core atom insertion fallible, expose structured identifier-capacity
  error kinds across graph, conformer, hierarchy, and topology APIs, and remove
  unchecked ID reconstruction and one-based writer arithmetic.
- v23: Add the focused immutable topology-transform and explicit
  topology-bound state-remapping public surface without expanding the prelude.
