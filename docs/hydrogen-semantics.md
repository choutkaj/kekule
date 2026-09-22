# Hydrogen representation and API migration

Explicit hydrogens are hydrogen atoms in the molecular graph. Implicit hydrogens
are attached hydrogens represented without separate atoms, including both a
specified count and any additional count inferred by the valence model.

`Molecule`, `MoleculeEditor`, and `Topology` expose the same counting contract:

| Method | Counts | Return value |
| --- | --- | --- |
| `explicit_hydrogens(atom)` | Hydrogen neighbors in the graph, including isotopic H | `Result<usize>` |
| `implicit_hydrogens(atom)` | Specified plus inferred non-graph hydrogens | `Result<Option<usize>>` |
| `total_hydrogens(atom)` | Explicit plus implicit hydrogens | `Result<Option<usize>>` |

The error type depends on the owner. Invalid atom IDs return an error. `None`
means inference is unresolved; it never means zero. Counts are widened before
adding contributions. These getters neither mutate the graph nor run perception.

## Specified counts and inference

`Atom::hydrogens` remains a `HydrogenDeclaration`:

- `Fixed(n)` specifies exactly `n` implicit hydrogens. This count is available
  without perception. Separate graph hydrogen neighbors are additional.
- `Infer { specified: n }` preserves `n` implicit hydrogens and permits the
  valence model to infer more. The complete implicit count is unknown until that
  additional contribution has been perceived, even when `n` is nonzero.

For example, `C` permits inference, `[CH4]` fixes four implicit hydrogens, and
`[C]` fixes zero. `[H][CH3]` has one explicit hydrogen and three implicit ones.
SMILES bracket notation specifies a count; it does not create hydrogen atoms.

Chemical edits invalidate inferred counts while retaining the declaration.
After reperception, an inference-enabled count can change with bond order or
charge. Fixed counts remain fixed: a conflicting edit produces a valence error
under strict perception, rather than silently removing specified hydrogens.
Property-only edits do not invalidate chemistry. The architecture's separation
between represented declarations and derived perception remains intact.

`inferred_hydrogens` on molecule, editor, topology, and perception surfaces is an
expert diagnostic for the installed inferred contribution alone. Likewise,
`ValencePerception::inferred_hydrogens` iterates those assignments, and
`PerceptionBuilder::with_valence` accepts inferred contributions. Applications
counting chemical hydrogens should normally use the combined getters instead.

## Conversions

`add_hydrogens` materializes resolved implicit hydrogens as graph atoms and bonds.
Fixed counts need no perception. Inference-enabled atoms require an installed
assignment; an absent or incomplete assignment fails before mutation.
`AddHydrogensOptions::specified_only` materializes only the specified contribution
without requiring inference. Origins in the report are `Specified` or `Inferred`.

`remove_hydrogens` suppresses eligible explicit hydrogens while preserving the
total count and stereo relationships. It retains hydrogens with individual
information that cannot be represented by a count, including isotope labels,
maps, charges, properties, and unsupported stereo roles. The report distinguishes
remaining explicit neighbors from the complete resulting implicit count; the
specified/inferred split remains available as diagnostic detail.

Conversions change representation, not protonation. They invalidate perception
after graph edits; recompute it before reading counts that require inference.
Suppression may use a temporary valence calculation to verify that its chosen
declaration preserves composition. Conversion does not promise to reconstruct
the original SMILES spelling or the original specified/inferred split.

## Migrating callers

This is an API and behavior change:

| Previous API | Replacement |
| --- | --- |
| `implicit_hydrogens` for the inferred contribution only | `inferred_hydrogens` |
| Manually adding declaration and inferred counts | `implicit_hydrogens` |
| `HydrogenDeclaration::Infer { explicit: n }` | `Infer { specified: n }` |
| `explicit_count()` / `with_explicit_count(n)` | `specified_count()` / `with_specified_count(n)` |
| `allows_implicit()` | `allows_inference()` |
| `AddHydrogensOptions::explicit_only` | `specified_only` |
| `AddedHydrogenOrigin::ExplicitCount` / `Implicit` | `Specified` / `Inferred` |
| Adjustment report's old `explicit_hydrogens` | `specified_hydrogens` |
| Adjustment report's old `implicit_hydrogens` | `inferred_hydrogens` |
| `PerceptionComponent::ImplicitHydrogens` | `InferredHydrogens` |
| `PerceptionBuildError::DuplicateImplicitHydrogen` | `DuplicateInferredHydrogen` |

In particular, do not retain a manual addition of the specified count around a
call to the new `implicit_hydrogens`: that would count the contribution twice.
DREIDING's `CountedHydrogens` error now reports a single combined `implicit` count.

The reference benchmark's serialized fields retain RDKit's external terminology:
`explicit_hydrogens` is its stored non-graph count and `implicit_hydrogens` is its
inferred count. The adapter deliberately uses the diagnostic APIs. Existing
goldens and strict comparisons remain unchanged, including storage disagreements.
