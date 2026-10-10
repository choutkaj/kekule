# Hydrogen representation

Explicit hydrogens are hydrogen atoms in the molecular graph. Implicit hydrogens
are attached hydrogens represented without separate atoms.

`Molecule`, `MoleculeEditor`, and `Topology` expose the same counting contract:

| Method | Counts | Return value |
| --- | --- | --- |
| `explicit_hydrogens(atom)` | Hydrogen neighbors in the graph, including isotopic H | `Result<usize>` |
| `implicit_hydrogens(atom)` | Non-graph hydrogens | `Result<Option<usize>>` |
| `total_hydrogens(atom)` | Explicit plus implicit hydrogens | `Result<Option<usize>>` |

The error type depends on the owner. Invalid atom IDs return an error. `None`
means an inferred count is unresolved; it never means zero. These getters
neither mutate the graph nor run perception.

## Fixed and inferred counts

`Atom::hydrogens` is an `ImplicitHydrogens`, and an atom's implicit count is
either wholly fixed or wholly inferred:

- `Fixed(n)` represents exactly `n` implicit hydrogens. This count is available
  without perception, and valence perception never adds to it.
- `Inferred` lets the valence model supply the complete count. It is unknown
  until valence perception is installed.

For example, `C` infers its count, `[CH4]` fixes four implicit hydrogens, and
`[C]` fixes zero. `[H][CH3]` has one explicit hydrogen and three fixed implicit
ones. SMILES bracket notation fixes a count; it does not create hydrogen atoms.
Molfile atoms fix their count when they declare `HCOUNT` or a valence.

Chemical edits invalidate inferred counts and keep fixed ones. After
reperception, an inferred count can change with bond order or charge. A fixed
count that conflicts with an edit produces a valence error under strict
perception, rather than silently removing hydrogens. Property-only edits do not
invalidate chemistry.

`Perception::inferred_hydrogens` and `ValencePerception::inferred_hydrogens`
expose the valence model's raw assignments, and `PerceptionBuilder::with_valence`
accepts them. An assignment is zero for every fixed atom: installing a detached
perception that infers hydrogens on a fixed count fails with
`PerceptionInstallError::InferredHydrogensOnFixedAtom`. Applications counting
chemical hydrogens use the molecule-level getters instead.

## Conversions

`add_hydrogens` materializes resolved implicit hydrogens as graph atoms and bonds.
Fixed counts need no perception and become `Fixed(0)`. Inferred atoms require an
installed assignment; an absent or incomplete assignment fails before mutation.
`AddHydrogensOptions::fixed_only` materializes only fixed counts without
requiring inference. Origins in the report are `Fixed` or `Inferred`.

`remove_hydrogens` suppresses eligible explicit hydrogens while preserving the
total count and stereo relationships. A parent keeps `Inferred` only when
inference on the collapsed graph reproduces its count; fixed parents, stereo
parents, aromatic `[nH]`, and parents whose count inference would not reproduce
(for example the hydrogens of `SH4` or `PH5`, or hydrogens on metals) receive
the count as `Fixed`. Each report adjustment records the resulting declaration
and complete implicit count. Removal retains hydrogens with individual
information that cannot be represented by a count, including isotope labels,
maps, charges, properties, and unsupported stereo roles.

Conversions change representation, not protonation. They invalidate perception
after graph edits; recompute it before reading inferred counts. Conversion does
not promise to reconstruct the original SMILES spelling.

## Migrating callers

Kekule 0.3.0 changes this API and its behavior from 0.2.1:

| 0.2.1 API | Replacement |
| --- | --- |
| `HydrogenDeclaration` | `ImplicitHydrogens` |
| `Infer { explicit: 0 }` | `Inferred` |
| `Infer { explicit: n }`, `n > 0` | `Fixed(n)` for a fixed count; graph hydrogens to keep inference |
| `explicit_count()` | `fixed_count()` (`None` when inferred) or `represented_count()` (zero when inferred) |
| `allows_implicit()` | `is_inferred()` |
| `with_explicit_count(n)` | Assign `Fixed(n)` |
| `Molecule::implicit_hydrogens(atom)` and `Topology::implicit_hydrogens(atom)`, which returned only the perceived part | `implicit_hydrogens` now returns the complete non-graph count (on the topology, through the atom view); `perception().inferred_hydrogens(atom)` returns the raw assignment |
| `Perception::implicit_hydrogens(atom)`, `ValencePerception::implicit_hydrogens()` | `inferred_hydrogens` |
| `AddHydrogensOptions::explicit_only` | `fixed_only` |
| `AddedHydrogenOrigin::ExplicitCount` / `Implicit` | `Fixed` / `Inferred` |
| Adjustment `explicit_hydrogens` / `implicit_hydrogens` (`u8`) | `usize` counts, plus `hydrogens` (the resulting declaration) |
| `HydrogenTransformError::HydrogenCountNotPreserved` | Removed; collapse preserves counts by construction |

`HydrogenCountPolicy::StoredOnly` now counts graph hydrogens and fixed counts.
