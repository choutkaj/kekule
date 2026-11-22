# Hydrogen representation

Explicit hydrogens are hydrogen atoms in the molecular graph. Implicit hydrogens
are attached hydrogens represented without separate atoms.

`Molecule`, `MoleculeEditor`, and `Topology` expose the same counting contract:

| Method | Counts | Return value |
| --- | --- | --- |
| `explicit_hydrogens(atom)` | Hydrogen neighbors in the graph joined by bonds of order one or more, including isotopic H | `Result<usize>` |
| `implicit_hydrogens(atom)` | Non-graph hydrogens | `Result<Option<usize>>` |
| `total_hydrogens(atom)` | Explicit plus implicit hydrogens | `Result<Option<usize>>` |

The error type depends on the owner. Invalid atom IDs return an error. `None`
means an inferred count is unresolved; it never means zero. These getters
neither mutate the graph nor run perception. A hydrogen joined only by a
zero-order or dative bond is not counted, as valence does not count it.
`descriptors::molecular_formula` and the masses count every hydrogen and fail
on an unresolved count.

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

The valence model infers no hydrogens on alkali, alkaline-earth, and p-block
metals: an unbonded neutral sodium is Na, not NaH, and a tin atom with three
carbons is not a stannane. Nor does a hydrogen atom joined only by zero-order
or dative bonds take a hydrogen of its own. Molfile and PDB or mmCIF
interpretation report a metal whose valence its element does not allow.

Molfile atoms fix their count only when they declare a valence (the V2000
valence field or V3000 `VAL`). The V2000 hydrogen-count field and V3000
`HCOUNT` are query constraints ("at least n hydrogens"), not hydrogen
declarations: a nonzero value is unsupported content, rejected by default or
reported and omitted on request.

A Molfile atom without a declared valence whose bonds imply a charge takes
the charge, not hydrogens. An onium centre, a neutral atom with more bonded
neighbors than its default valence such as the sulfur of trimethylsulfonium,
takes +1. A neutral atom with more bonds than any neutral valence allows, such
as a quaternary nitrogen, takes +1 or -1. A singly bonded chalcogen with no
other non-metal neighbor, or a boron or aluminium neighbor, takes -1 with its
+1, so DMSO drawn with a single S-O bond is S+ with O-. A non-metal atom
bonded to metals whose bond valence its element allows only with some of those
bonds not counted, such as the phosphorus of a phosphine on palladium, is a
donor: it keeps its charge and takes no hydrogens. The Molfile writers declare
a valence wherever this reader or the default-valence rule of other readers
would infer a different count, and write dative bonds toward their metal end.

Chemical edits invalidate inferred counts. After reperception, an inferred
count can change with bond order or charge. Topology and structure edits keep a
fixed count whole: an atom gains a hydrogen for each unit of bond valence a
deletion, a lower order, or a subset cut takes from it, and gives one up for
each unit a new bond or a higher order adds. A deleted explicit hydrogen of a
fixed atom stays as an implicit one. A gain takes a declared radical electron
first, then a hydrogen, then fills open valence; past the largest allowed
valence it is owed: later losses pay it back first, and publication fails with
`TopologyEditError::HydrogenDeficit` while any is owed, rather than silently
dropping valence. A loss adds hydrogens only for valence the radical electrons
and hydrogens leave unfilled, so deleting a hydrogen of an atom declared a
radical leaves the radical. A metal's count never changes. A non-metal's bond
to a metal changes no count when the non-metal donates to it, as an ammine
does; otherwise it is covalent, so a chloride or a methyl cut from its metal
takes a hydrogen.
To remove a hydrogen for good, as for hydroxide from water, delete it and set
the atom's charge and count with `TopologyEditor::replace_atom`. A molecule draft (`MoleculeEditor`) edits the
graph as stated and keeps counts as they are;
`perception::valence::open_valences` names fixed atoms short of their valence.
Property-only edits do not invalidate chemistry.

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
`AddHydrogensScope::FixedOnly` materializes only fixed counts without
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

`StructureEditor::add_hydrogens` and `StructureEditor::remove_hydrogens` apply the
same conversions to chosen atoms of a structure draft, placing new hydrogens at
ideal positions; inferred counts come from the valence model on the draft's
current graph, so they need no perception. See
[placement](../crates/kekule/src/structure/placement.rs).

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
| `Molecule::implicit_hydrogens(atom)` and `Topology::implicit_hydrogens(atom)`, which returned only the perceived part | `implicit_hydrogens` now returns the complete non-graph count (on the topology, through the atom view); `perception().valence()` and `ValencePerception::inferred_hydrogen_count(atom)` return the raw assignment |
| `Perception::implicit_hydrogens(atom)`, `ValencePerception::implicit_hydrogens()` | `inferred_hydrogens` |
| `AddHydrogensOptions::explicit_only` | `fixed_only` |
| `AddedHydrogenOrigin::ExplicitCount` / `Implicit` | `Fixed` / `Inferred` |
| Adjustment `explicit_hydrogens` / `implicit_hydrogens` (`u8`) | `usize` counts, plus `hydrogens` (the resulting declaration) |
| `HydrogenTransformError::HydrogenCountNotPreserved` | Removed; collapse preserves counts by construction |

`HydrogenCountPolicy` is removed: formula and masses take no policy and count
every hydrogen.
