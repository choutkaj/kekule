
## Molecular SMARTS dialect

`QueryGraph` is the reusable compiled query: parse once and match many targets.
`substructure::PreparedTarget` caches target facts behind an immutable borrow.
`PreparedTopologyTarget` also reuses connected, component-local searches once per
molecule definition and expands their mappings to occurrences. Queries with
disconnected outer or recursive graphs search the full topology instead.

The practical reference is RDKit **2026.03.3**, with OpenFF Toolkit **0.18.0**
tag projection conventions. This is a declared molecular dialect, not a claim
to implement every Daylight, RDKit, or CXSMARTS extension.

| Feature | Meaning / prerequisite | Regression coverage |
| --- | --- | --- |
| `*`, symbols, `#n`, isotope, charge | Represented element/isotope/charge; element symbols impose aromatic or aliphatic state | `query`, `smarts_core` |
| `a`, `A`, aromatic symbols | Installed aromaticity; a named target model is selected explicitly | `mdl_aromaticity_is_explicit_and_does_not_rewrite_chemistry` |
| `D` | Graph neighbors, default 1 | `query` |
| `X` | Graph neighbors plus nongraph hydrogens, default 1; valence required | `query` |
| `H` | Total attached graph, declared, and inferred hydrogens, default 1; valence required | `hydrogen_valence_and_ring_primitives_are_distinct` |
| `[H]`, `[2H+]`, `[H:1]` | Hydrogen element when the bracket content is a hydrogen atom declaration | `query`, `smarts_core` |
| `h`, `hN` | Declared plus inferred nongraph hydrogens; bare `h` means at least one | `hydrogen_valence_and_ring_primitives_are_distinct` |
| `v` | Localized bond-order sum plus nongraph hydrogens, default 1; valence required | `hydrogen_valence_and_ring_primitives_are_distinct` |
| `R`, `R0`, `r`, `r0`, `x` | True cyclic membership, independent of ring basis; ring perception required | `query`, `smarts_core` |
| `RN`, `rN` | Selected symmetrized ring count and smallest selected ring size; identified Figueras SSSR-like basis required | `hydrogen_valence_and_ring_primitives_are_distinct`, `incompatible_ring_models_fail_and_mdl_failure_is_transactional` |
| `xN` | Count of incident cyclic bonds | `query` |
| Bonds `- = # $ : ~ @` | Nonaromatic single/double, triple, quadruple, aromatic, any, cyclic; appropriate perception required | `bond_boolean_operators_and_precedence`, `query` |
| Omitted bond | Single or aromatic | `query` |
| `!`, `&`, implicit AND, `,`, `;` | Atom and bond Boolean logic, in descending precedence | `bond_boolean_operators_and_precedence`, `query` |
| `$()` | Nested existential subquery, first atom anchored to the tested atom; independent of the outer embedding | `recursive_queries_are_anchored_nested_and_independent_of_outer_embedding` |
| Branches / ring closures | Non-induced injective graph matching | `query` |
| `.` | Query fragments may match in the same or different target components | `topology_disconnected_queries_keep_occurrences_and_snapshot` |
| `:n` | Query label, never a test of target atom maps | `tags_are_outputs_and_do_not_match_target_map_values` |
| `@`, `@@`, `@TH1`, `@TH2` | Local represented tetrahedral configuration, checked under atom mapping | `query`, `boolean_stereo_keeps_atom_expression_context` |
| `/`, `\` | Relative bond directions, jointly evaluated under a complete mapping; lone directions impose no stereo relationship | `query`, `directional_alternatives_keep_boolean_context` |

### Explicit boundaries

- Reactions, component grouping such as `(C).(C)`, CXSMARTS, `@?`, `/?`,
  non-tetrahedral classes, hybridization predicates, numeric ranges, and directed
  dative extensions are outside this dialect. Valid recognized unsupported
  syntax reports `Unsupported`; malformed syntax reports `InvalidSyntax`.
- Numeric hydrogen, valence, degree, and ring predicates use `u8`; isotopes use
  `u16`, charge uses `i8`, and tags use `u32`. Out-of-range values are errors,
  never wrapped. Recursive nesting has an absolute ceiling of 32.
- Molecular elements are 1 through 118. The core has no dummy target atoms.
- Numeric `R/r` reject fallback or unidentified ring bases. The default ring
  algorithm labels deterministic depth-first recovery as `DepthFirstFallback`;
  Boolean cyclic membership remains usable. This prevents an approximate ring
  set from silently answering an exact ring-count question.
- `v` follows represented localized valence. Kekule's existing zero/dative
  bonds contribute zero at both endpoints; directional dative chemistry is
  outside this compatibility claim.
- Boolean stereo is interpreted literally in its expression, rather than
  copying RDKit's incidental parser behavior. For example,
  `N[C@,C@@](F)(Cl)Br` accepts both specified configurations here, and
  `F/,\C=C/Cl` accepts either specified alkene configuration. RDKit 2026.03.3
  effectively retains only one configuration in these examples. Negation is
  logical complement, so negated atom stereo can accept unspecified targets.
  An unnegated stereo primitive requires specified target stereo.
- Stereo alternatives on one atom must share the same inline-hydrogen carrier
  convention. Conflicting frames and contradictory literal directions are
  rejected. CIP descriptors and enhanced stereo groups do not affect matching.

### Enumeration and errors

Existing `find_substructure_matches` convenience behavior remains bounded:
1,000 results and target-atom-set deduplication by default. Use
`find_substructure_matches_complete` for complete enumeration; its match limit
is a hard bound and exceeding it discards the collected result with an error.
Set `uniquify: false` to retain automorphisms. Streaming callbacks report
`Complete` or `Stopped`; callbacks delivered before an error are provisional.

`TaggedQuery` validates positive unique outer tags and returns tuples in ascending
tag order. Tags on recursive subqueries are local metadata. Tagged search always
retains tag permutations; optional deduplication removes only identical ordered
tagged tuples. Its match bound counts full embeddings before tag projection.

Parsing bounds aggregate atoms, bonds, and expression nodes across recursive
queries. Matching shares candidate and search budgets across recursion and
topology definitions. A failed recursive search is an error, including underneath
negation. Results are deterministic for one graph; traversal order is not a
canonical chemical ordering. Topology matches retain their exact `Arc<Topology>`.

### SMIRNOFF preparation

SMARTS support supplies tagged environments, not OFFXML parsing or force-field
assignment. SMIRNOFF uses MDL aromaticity and explicit hydrogen vertices:

```rust
use kekule::{core::AromaticityModel, hydrogens, perception::aromaticity,
             query::parse_smarts, smiles, substructure::*};
let mut molecule = smiles::to_molecules("CO")?.pop().unwrap();
molecule.perceive()?;
hydrogens::add_hydrogens(&mut molecule)?;
molecule.perceive()?;
aromaticity::perceive_aromaticity(&mut molecule, AromaticityModel::Mdl)?;
let target = PreparedTarget::new(&molecule);
let query = parse_smarts("[#6:1]-[#8:2]")?;
let tagged = TaggedQuery::new(&query)?;
let matches = tagged.find_matches(&target, SubstructureMatchOptions {
    max_matches: 100_000, uniquify: false, ..Default::default()
}, true)?;
assert_eq!(matches.len(), 1);
# Ok::<(), Box<dyn std::error::Error>>(())
```

The `Mdl` model follows the pinned RDKit MDL implementation: eligible C/N
one-electron donors, no exocyclic multiple bonds or triple bonds, and fused
Hückel evaluation with a minimum of six contributing ring atoms. It preserves
localized chemistry, validates ring provenance, and installs transactionally.
The ordinary perception default remains `RdkitLike`.

Parameter precedence, improper-torsion permutations, virtual-site orientation,
charge generation, and fractional bond orders belong to the future parameterizer.

References: [Daylight SMARTS](https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html),
[RDKit SMARTS](https://www.rdkit.org/docs/RDKit_Book.html#smarts-support-and-extensions),
[SMIRNOFF](https://openforcefield.github.io/standards/standards/smirnoff/), and
[OpenFF matching adapter](https://docs.openforcefield.org/projects/toolkit/en/0.18.0/_modules/openff/toolkit/utils/rdkit_wrapper.html).
