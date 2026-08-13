# Public API Facade

## Summary

Expose the architecture-defined public facade instead of a flat root namespace.

## Behavior/API

- The foundation package and Rust import root are both named `kekule`; there is
  no compatibility package or alias under a previous project name.
- Public modules are focused around `core`, `units`, `small`, `bio`, `smiles`,
  `molfile`, `sdf`, `mmcif`, `perception`, `hydrogens`, `query`,
  `substructure`, `canon`, `descriptors`, `geometry`, `topology`, `structure`,
  `alignment`, and `modeling`. Ordered trajectory state, codecs, and focused
  trajectory workflows live in the one-way `kekule-traj` companion, whose
  specialized options, reports, and errors remain under `kekule_traj::analysis`
  and `kekule_traj::io`.
- The crate root no longer blanket re-exports implementation modules.
- The prelude is intentionally small and limited to common user-facing types.
- `SmallMolecule` owns small-molecule convenience methods and hides its raw graph field behind `graph()`, `graph_mut()`, and `into_graph()`.
- `MacroMolecule` exposes read-only graph/hierarchy access plus checked
  construction and transactional coordinated editing; completed values cannot
  be independently mutated into an invalid graph/hierarchy pair.
- Macro graph growth uses atomic bonded-atom insertion after the first atom and
  reports `MacroGraphEditError` rather than exposing an unattached public atom.
- `MacroMolecule` exposes direct hierarchy iterators, atom-site lookup, and
  read-only validation. The placeholder macro sanitization surface is absent.
- SMILES, Molfile, SDF, and mmCIF expose format-specific Documents and explicit
  interpretation results with reports/mappings; superseded direct reader APIs
  are absent.
- SMILES and Molfile retain simple default-bounded `parse_str` entry points and
  expose focused parse-options overloads; SDF and mmCIF accept their parse
  options directly.
- Dot-separated SMILES interpretation exposes one connected component result
  per component. Single-molecule accessors and convenience readers return a
  structured component-count error when the document contains more than one.
- Molfile and SDF interpretation reject CTAB records whose graph contains more
  than one component. mmCIF completes declared connectivity and then partitions
  any remaining components into separate connected molecule instances.
- `mmcif::write` exposes explicit supported `Model` serialization with
  format-specific options and structured rejection errors.
- `Molecule` is one asserted connected entity; empty and single-atom values are
  valid boundary cases. `Molecule::builder()` and transactional `edit()` are the
  public topology-construction surfaces. Their staging graphs are not exposed;
  only `build()` and `commit()` can publish a completed graph, and both reject a
  disconnected result without exposing partial mutation.
- Graph staging, conformer, stereo, SMCRA hierarchy, and topology insertion
  surfaces report focused fixed-width identifier capacity errors before
  mutation rather than truncating or panicking.
- `Molecule::formal_charge` exposes the asserted live-atom charge aggregate as
  an `i64` without hiding sanitization or perception.
- `mmcif::interpret` returns a selected-coordinate `Model` plus report;
  `MolecularContents` and `Solvent` are removed.
- `mmcif::interpret_ensemble` is a separate shared-topology multi-model path
  that rejects inconsistent atom identity or topology.
- `MmcifInterpretation::into_model` consumes an interpretation when callers do
  not need to retain its report.
- Expert perception functions live under focused modules such as `perception::rings`, `perception::aromaticity`, and `perception::valence`.
- Focused canonical reconstruction remains under `core` and `bio`: detached
  exact perception construction plus whole-state installation, targeted SMCRA
  child enrichment, and exact stereo-group tombstone replay do not widen the
  crate root or prelude.
- Immutable system topology and direct model-state containers live under
  `topology` and `structure`; potentials and minimization remain under
  `modeling`. These focused types are not added to the prelude.
- Definition-local SMCRA IDs remain under `bio`; `topology` supplies canonical
  instance-qualified chain, residue, and atom-site IDs and navigation. `Model`
  and `ModelView` forward the common read-only surface without wrapper objects
  or duplicated hierarchy state.
- The topology facade names complete static equality `Topology::same_layout`
  and restricts `TopologyMapping::between_identical_layouts` to identity maps
  over that exact layout. Checked explicit mappings and topology-edit results
  expose structured consistency failures.
- Raw `Topology` directly owns its private layout data and is not cloneable;
  topology-bound containers use `Arc<Topology>` for cheap sharing and exact
  compatibility without exposing a public identity type.
- Immutable whole-instance edits live under the focused
  `topology::transform` namespace. Core structure containers expose explicit
  remapping methods and typed errors; `kekule-traj` builds frame and trajectory
  remapping on that public lineage contract.
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
- Same-topology selection-based proper rigid fitting lives under `alignment`,
  returns the existing geometry transform plus unit-bearing RMSD, and is not
  added to the prelude.

## Implementation Notes

- Existing algorithm and I/O internals remain available through focused facade modules rather than root aliases.
- `SmallMolecule::from_smiles` orchestrates parse/interpret without sanitizing
  and accepts exactly one dot component; `from_smiles_sanitized` names the
  additional operation explicitly and has the same component-count boundary.
- `graph_mut()` itself is state-neutral; chemistry and topology mutators on the
  returned graph perform their own targeted invalidation, allowing perception
  operations to consume already-installed prerequisite state.
- Internal benchmark tooling uses the same public namespaces as user code.
- Invariant-bearing hierarchy, provenance, document, model, and structured
  error state is private behind accessors or checked constructors.
- Builder/editor staging `Molecule` state is likewise private; public callers
  receive focused topology operations and perform ordinary graph operations
  only after successful finalization.
- Extensible public error enums are non-exhaustive. Deliberate value, options,
  and report payloads may retain direct public fields.
- The topology-centered API and Kekule package names form the public contract;
  breaking changes in the `0.x` line require a minor version increment.

## Tests

- External integration tests compile public happy-path, namespaced, low-level
  graph, macro-molecule, topology/model, mmCIF, and borrowed-view
  examples as downstream user code.
- Downstream reconstruction tests compile and exercise every new checked
  perception, hierarchy enrichment, stereo-slot, and topology-layout surface.
- Public topology tests compile the shared-Arc, same-layout, checked
  identical-layout mapping, and checked edit-result surface, while compile-fail
  coverage keeps the removed raw-clone and identity APIs unavailable.
- Public and unit hierarchy tests cover reused macro definitions, qualified
  parent navigation, small-molecule absence, zero-copy model views, and
  instance-specific chain/residue selections.
- Downstream transformation tests compile instance retain/remove, stable
  mapping traversal, model and selection remapping, and ensemble remapping;
  `kekule-traj` separately covers owned trajectory and reusable-buffer remaps.
- A downstream alignment test compiles the focused options, result, structured
  error, topology-selection, transform-direction, and RMSD-unit surface.
- A downstream `kekule-traj` analysis test compiles explicit direct RMSD,
  transactional superposition, and fused aligned RMSD without adding companion
  analysis types to the foundational Kekule facade or prelude.
- Workspace tests exercise the benchmark tooling and existing chemistry/IO behavior through the new wrapper accessors.
- Downstream tests verify connected builder/editor behavior, structured
  multi-component SMILES handling, disconnected Molfile/SDF rejection, and
  connected mmCIF topology instances.
- Downstream compile-fail tests verify builder/editor staging cannot be cloned
  or taken out as a public `Molecule`.

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
  initial `0.1.0` contract.
- v16: Expose configurable resource-bounded SMILES, Molfile, and SDF parsing
  without widening the crate root or prelude.
- v17: Add the core formal-charge aggregate, transactional model-instance
  conformer export, and concise basic versus modeling public examples.
- v18: Add consuming `MmcifInterpretation::into_model` access through the
  focused mmCIF facade.
- v19: Add the focused `descriptors` facade for explicit-policy molecular
  formula and mass calculation without expanding the prelude.
- v20: Add the `geometry`, `topology`, `structure`, and `trajectory`
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
- v24: Add the focused same-topology weighted rigid-alignment analysis surface
  without expanding the crate root aliases or prelude.
- v25: Expose focused canonical persistence reconstruction through existing
  `core` and `bio` modules without adding Serde, root aliases, or prelude items.
- v26: Hard-rename the published foundation package and Rust import root to
  `kekule` for the initial release without compatibility aliases.
- v27: Keep the foundational facade focused by moving ordered frames,
  trajectory storage, streaming traits, and codecs to `kekule-traj` without a
  compatibility module in `kekule`.
- v28: Add the companion's focused `analysis` namespace for explicit direct
  RMSD, transactional superposition, and fused aligned RMSD while keeping the
  foundational crate and prelude unchanged.
- v29: Set the foundation and companion packages to the shared initial
  `0.1.0` release line and synchronize their internal dependency requirements.
- v30: Align the facade with connected `Molecule` boundaries, component-aware
  SMILES interpretation, disconnected CTAB rejection, and connected mmCIF
  instance partitioning without widening the crate root or prelude.
- v31: Seal builder/editor staging graphs behind focused operations so `build`
  and `commit` are the only public routes to a completed `Molecule`.
- v32: Make raw `Topology` a directly owned non-cloneable value and move exact
  compatibility and cheap sharing to `Arc<Topology>` throughout the focused
  structure, selection, mapping, alignment, potential, and companion APIs.
- v33: Replace the configuration and structure-observation exports with
  topology-bound `AtomData` and the flattened `Model`/`ModelView` API.
- v34: Reconcile the flattened API around directly borrowed `Positions`,
  field-specific `AtomData` construction, and quantity-valued B-factors.
- v35: Add canonical instance-qualified SMCRA identities and expose coherent
  zero-copy hierarchy navigation through `Topology`, `Model`, and `ModelView`.
