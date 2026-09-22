# Selection audit — 2026-09-22

Audited baseline: `dab25109` (`origin/main`). The implementation added with this
audit extends the existing selection contract without changing ownership or
structural-subset semantics. Current API inventories and examples live in the
[atom selection Rustdoc](../crates/kekule/src/topology/selection.rs) and
[bond selection Rustdoc](../crates/kekule/src/topology/selection/bonds.rs).

## Baseline capabilities

`kekule::topology::AtomSelection` owned an `Arc<Topology>` and a sorted, unique
`Vec<TopologyAtomIndex>`. Its authoritative ordering was topology dense order,
not input order. Empty selections were supported through the checked constructors.

| Selection source or operation | Baseline API |
| --- | --- |
| All atoms, explicit qualified atom IDs, dense indices | `all`, `from_atoms`, `from_indices` |
| Complete molecule occurrences; every occurrence of a reusable definition | `for_instances`, `for_definitions` |
| Element and intrinsic molecule classification | `for_elements`, `for_molecule_classes` |
| Chain, residue, atom site, residue class, chain label | `for_chains`, `for_residues`, `for_atom_sites`, `for_residue_classes`, `for_chain_label` |
| Exact atom names in separate identifier namespaces | `for_label_atom_names`, `for_author_atom_names` |
| Molecule-local substructure results | `from_query_matches(topology, instance, matches)` |
| Combining sets | `union`, `intersection`, `difference` |
| Whole-residue expansion | `expand_to_residues` |
| Access and validation | `atom_ids`, `indices`, `topology`, `semantic_ids`, `ensure_compatible` |

Additional selection-related capabilities already existed outside that type:

- [`measure::within`](../crates/kekule/src/structure/measure.rs) returns candidate
  atoms within an inclusive Cartesian cutoff of any reference atom. It consumes
  a borrowed realization and two selections on its exact snapshot. It is
  nonperiodic and produces a static result for that realization.
- [SMARTS topology matching](../crates/kekule/src/algorithms/substructure/topology.rs)
  returns `TopologyQueryMatch` objects with qualified atom IDs and a retained
  topology handle. Matching and tagged output are separate from selection sets.
- [Hierarchy lookup](../crates/kekule/src/topology/lookup.rs) supports exact
  label/author addressing, insertion codes, and explicit ambiguity errors.
- `Topology::subset`, `Model::slice`, `Ensemble::slice`, and `Trajectory::slice`
  accept atom selections. They construct induced structural subsets and
  repartition cut molecules, transferring supported entity data.
- Editors already accept source atom/bond IDs through `atom_handle` and
  `bond_handle`, then operate on draft handles. This is distinct from selection.

## Findings addressed

The audit found missing interaction capabilities, rather than evidence that
the existing sorted-set and snapshot contracts were incorrect.

1. There was no independently usable bond selection. Added `BondSelection` with
   validated ID/index constructors, complete-instance/definition constructors,
   predicates, semantic iteration, and the same exact-snapshot contract.
2. Basic UI actions required reconstructing sets manually. Both types now expose
   `empty`, `len`, `is_empty`, `contains`, `contains_index`, `insert`, `remove`,
   `toggle`, `clear`, `complement`, `symmetric_difference`, `is_subset`, and
   `is_disjoint`, alongside union/intersection/difference.
3. Atom-to-bond conversion had no explicit semantics. Added
   `to_bonds(BondSelectionMode::{Internal, Incident, Boundary})` for both,
   either, or exactly one selected endpoint. `BondSelection::to_atoms` selects
   both endpoints. The round trip may add other bonds between those atoms.
4. Whole-molecule, whole-chain, and bounded graph-neighborhood expansion were
   absent. Added `expand_to_instances`, `expand_to_chains`, and `expand_bonded`.
   Hierarchy expansion retains selected atoms without hierarchy assignments;
   covalent traversal never crosses instance boundaries.
5. Custom property/chemistry selection needed manual ID collection. Added
   `from_predicate` and `filter` for atoms and bonds. These accept qualified IDs
   and borrowed chemical records, so callers can inspect existing topology
   properties and perception without a second expression language.
6. Snapshot-bound topology matches lacked a checked selection constructor.
   Added `from_topology_query_matches`, which checks every match's source before
   combining atoms. The older local-match constructor remains available and now
   documents that its caller is responsible for target provenance.
7. Applications could borrow the topology from a selection but could not clone
   its retained handle publicly. `shared_topology` is now public for atoms and
   bonds, allowing new selections to retain the same snapshot.

## Viewer/editor integration boundaries

Keep an atom selection and a bond selection when the application allows
independent picking. Derive highlighted bonds from selected atoms only when
that is the application's intended interaction. Chain/residue constructors
produce atom sets; they do not retain independently selected hierarchy nodes.
Empty hierarchy nodes therefore cannot be represented by an atom selection.

Use union for adding a group, difference for removing a group, symmetric
difference for toggling a group, and complement for inversion. Single-member
mutations validate before changing state. `insert`/`remove` return whether the
set changed; `toggle` returns the resulting selected state. Invalid membership
queries return false; invalid mutations return an error.

Selections retain exact snapshot identity, not just chemically equivalent
contents or compatible array lengths. Geometry-only changes preserve selection
bindings. Structural edits and owning perception can publish new snapshots;
old selections remain attached to the old snapshot. Rebuild deliberately when
publishing a new one. General editor publication does not currently expose a
source-to-result selection correspondence. Subsets have their own explicit
correspondence; this audit does not introduce universal remapping.

Bare semantic IDs and dense indices do not encode their source snapshot. An
in-range ID copied from another snapshot cannot be recognized as foreign by a
constructor. Associate renderer pick IDs with the displayed topology and check
selection compatibility before applying results to an owner/editor.

Selection sets intentionally discard click order, query order, duplicate hits,
tags, and match boundaries. Retain an active pick or ordered measurement list
separately. Implicit hydrogens are counts, not independently selectable atoms.

Screen-space picking, lasso/frustum selection, named selection storage, undo,
and live query reevaluation are application concerns built on these sets.
Spatial selection currently supplies Cartesian proximity, without periodic
imaging, a spatial index, or additional shape constructors. Reevaluate explicitly
when positions change. Molecular-local selections without a topology and an
independent hierarchy-node selection type are not provided.

The storage remains a sorted vector: dense-index membership is logarithmic,
single insertion/removal can shift a linear number of members, and binary set
algebra merges two sorted vectors linearly. Large group operations should use
constructors and set algebra rather than repeated individual insertion. This
audit makes no benchmark claim for large interactive scenes.

## Regression coverage

[Public editor-selection regressions](../crates/kekule/tests/selection_editor_public_api.rs)
cover an independent bit-mask oracle for set algebra, invalid mutation atomicity,
sorting/deduplication, equal-layout foreign snapshots, reused definitions,
endpoint conversion rules, ring round trips, isolated atoms, graph expansion,
predicates, and query-match provenance.
[Scientific selection regressions](../crates/kekule/tests/scientific_selection_public_api.rs)
exercise hierarchy expansion across molecule boundaries and retention of
unassigned atoms. Rustdoc examples compile interactive selection and selected
bond deletion workflows. Existing spatial, subset, and trajectory tests remain
the regression coverage for their unchanged contracts.

## Validation

Executed on Windows; all commands below passed. The workspace test run reported
1,383 passing tests including doc tests. A Clippy warning in the new test was
corrected and the full Clippy command rerun successfully.

```text
cargo fmt --all -- --check
cargo check --workspace --all-targets --all-features --locked
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo test --workspace --all-features --locked
cargo test --workspace --all-features --doc --locked
cargo doc --workspace --all-features --no-deps --locked
cargo test -p kekule-potentials --no-default-features --locked
cargo doc -p kekule-potentials --no-default-features --no-deps --locked
cargo package -p kekule --locked --allow-dirty
cargo package -p kekule-potentials --locked --allow-dirty --list
cargo package -p kekule-traj --locked --allow-dirty --list
cargo +1.89.0 check --workspace --all-targets --all-features --locked
cargo check --manifest-path fuzz/Cargo.toml --bins --locked
git diff --check
```

Documentation builds used `RUSTDOCFLAGS="-D warnings"`. Packaging used
`--allow-dirty` to validate the uncommitted working-tree changes. Companion
packages used the file-list checks prescribed by CI, since full packaging
depends on registry publication of the foundational crate.

Nightly fuzz execution and external-reference benchmarks were not run: this
change adds deterministic selection operations, with no format parser or
scientific interpretation change. Fuzz targets were compilation-checked.
Unrelated Python/JavaScript maintenance tests and Linux CI were not run locally.
