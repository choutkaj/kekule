# Selection language

Kekule's selection language names atoms with short text expressions such as
`protein and name CA` or `within 5 of resname LIG and not water`. It follows
the conventions of VMD, PyMOL and MDAnalysis where they agree, so people and
language models can write it from what they already know. An expression is
parsed once into a `select::SelectionExpr` and evaluated against a topology or a
structure into a `topology::AtomSelection`.

The API lives in [`kekule::select`](../crates/kekule/src/lib.rs); the
implementation is in [`selection_language`](../crates/kekule/src/selection_language/mod.rs).

```rust
use kekule::select::{self, SelectEvaluateOptions, SelectParseOptions};

let expr = select::parse("protein and name CA", SelectParseOptions::default())?;
let atoms = expr.evaluate(topology.clone(), SelectEvaluateOptions::default())?;
// Spatial terms need coordinates:
let near = select::parse("within 5 of resname LIG", SelectParseOptions::default())?;
let pocket = near.evaluate_structure(&model, SelectEvaluateOptions::default())?;
```

## Grammar

```text
expression := or
or         := and ("or" and)*
and        := unary ("and" unary)*
unary      := "not" unary
            | "within" number [unit] "of" unary
            | ("byres" | "bychain" | "bymolecule") unary
            | "same" ("residue" | "chain" | "molecule") "as" unary
            | "bonded" "to" unary
            | primary
primary    := "(" expression ")" | term
term       := flag | keyword values | keyword comparison value | host
values     := value+
value      := word | quoted | number | range
range      := number "-" number | number "to" number
comparison := "<" | "<=" | ">" | ">=" | "==" | "=" | "!="
```

- **Precedence.** `not` and the other prefix operators bind tightest, then
  `and`, then `or`: `not water and ion or protein` means
  `((not water) and ion) or protein`. A prefix operator applies to the one
  term or parenthesized group after it, so
  `within 5 of resname LIG and not water` selects atoms near the ligand that
  are not water. Write `within 5 of (resname LIG and not water)` for the
  other reading.
- **Lists.** Several values after a keyword match any of them:
  `resname ALA GLY` is `resname ALA or resname GLY`. A list ends at the next
  operator, keyword, parenthesis, or host term. Juxtaposed terms without an
  operator are an error, never an implicit `and`.
- **Words and quotes.** A bare word is any run of characters other than
  whitespace, parentheses, quotes, and the comparison characters `<>=!`, and
  may contain `'` and `*` after its first character, as atom names such as
  `C1'` do. Quote a value with `"..."` or `'...'` when it is a reserved word,
  contains other characters, or is empty; `\"`, `\'` and `\\` escape inside
  quotes. Keywords are lower case and case-sensitive, as in VMD, so upper-case
  values such as the residue names `ION` or `ALL` never read as keywords.
  Values are case-sensitive except element symbols.
- **Numbers and ranges.** Ranges are inclusive: `resid 10-50`, `resid 10 to 50`,
  `charge -2 to -1` (or `-2--1`). Residue numbers with insertion codes order by
  number, then code, with no code first: `resid 50-52` includes `51A` but not
  `52A`, and `resid 52` matches only residue 52 without an insertion code.
- **Units.** Distances are in ångströms unless a unit follows the number:
  `A` or `Å`, `nm`, `pm`, or `bohr`. B factors compare in Å². A real property
  value compares in its column's unit unless a unit follows the number: one
  of the units above or `Da`, `e`, `K`, `ps`, `fs`, `kJ/mol`, `kcal/mol`,
  `deg`, or `rad`. Comparisons are exact; a value converted between units
  may differ from a stored value in its last bits.

## Terms

Identifier terms follow **author identifiers by default**, as PDB files,
viewers and papers do. An entity without an author identifier, such as a chain
of a MOL/SDF synthetic hierarchy, matches on its label identifier. The
`label_*` terms match mmCIF label identifiers only.

| Term | Matches | Example |
| --- | --- | --- |
| `all`, `none` | Every atom, no atom | `all` |
| `index` | Dense atom index | `index 0-99` |
| `element` (`elem`) | Element symbol | `element C N O` |
| `name` | Atom name: author atom ID, else label atom ID | `name CA "C1'"` |
| `label_name` | Label atom ID only | `label_name OXT` |
| `charge` | Formal charge | `charge < 0` |
| `isotope` | Mass number; atoms without one never match | `isotope 13` |
| `bfactor` | Isotropic B factor in Å² (structure) | `bfactor > 50` |
| `occupancy` | Occupancy (structure) | `occupancy < 1` |
| `property` | Atom property by key; a value list, a comparison, or a bare key for `true` booleans | `property partial_charge < -0.5` |
| `resname` | Residue name | `resname ALA GLY` |
| `resid` | Author residue number with insertion code, else label number | `resid 10-50 52A` |
| `label_resid` | Label sequence number | `label_resid 1-20` |
| `resclass` | Residue class | `resclass amino_acid` |
| `chain` | Author chain ID, else label chain ID | `chain A B` |
| `label_chain` | Label chain ID | `label_chain C` |
| `molclass` | Molecule class | `molclass small_molecule` |
| `instance` | Molecule instance ID | `instance 0 3` |
| `definition` | Molecule definition ID | `definition 1` |
| `aromatic` | Aromatic atoms (perception required) | `aromatic` |
| `ring` | Ring atoms (perception required) | `ring and element N` |
| `smarts` | Atoms of any match of a SMARTS query (perception required) | `smarts "[OX2H]"` |
| `hydrogen` | Hydrogen atoms of any isotope | `hydrogen` |
| `heavy` | Non-hydrogen atoms | `heavy` |
| `polar_hydrogen` | Hydrogens bonded to N or O | `polar_hydrogen` |
| `protein` | Atoms of protein molecules | `protein` |
| `nucleic` | Atoms of DNA and RNA molecules | `nucleic` |
| `water`, `ion`, `ligand` | Atoms of water, ion, and small-molecule molecules | `ligand` |
| `backbone` | `N CA C O` of amino-acid residues in proteins; `P OP1 OP2 OP3 O1P O2P O5' C5' C4' C3' O3'` of nucleotides | `backbone` |
| `sidechain` | Atoms of amino-acid residues in proteins that are neither backbone nor `OXT`, nor hydrogens bonded to those | `sidechain and heavy` |

Residue classes are `amino_acid`, `dna`, `rna`, `carbohydrate`, `water`,
`ion`, and `other`. Molecule classes are `protein`, `dna`, `rna`,
`carbohydrate`, `water`, `ion`, `small_molecule`, and `other`. After `resclass`
and `molclass`, class names are values even where they are also keywords, so
`molclass water ion` lists two classes. The macros use
Kekule's intrinsic classification, not contextual roles: `ligand` is every
small molecule, cofactors included, and never water or an ion. Implicit
hydrogens are counts, not atoms, so no term selects them.

`property` looks the key up in the structure's atom properties (when a
structure is evaluated), the topology's atom properties, and the molecule
definitions' atom properties, and fails when the key is in more than one
scope or in none. Atoms without a value never match.

Terms marked *structure* and `within` need coordinates and fail with
`SelectEvaluateError::RequiresStructure` when evaluated against a topology
alone. Chemistry terms read installed perception and fail with
`SelectEvaluateError::MissingPerception` when it is absent; evaluation never
perceives.

## Operators

| Operator | Selects |
| --- | --- |
| `a and b`, `a or b`, `not a` | Intersection, union, complement |
| `within d of a` | Atoms at Cartesian distance at most `d` from any atom of `a`, `a` included; no periodic imaging |
| `byres a` (`same residue as a`) | Every atom of a residue containing an atom of `a`, plus `a`'s atoms without a residue |
| `bychain a` (`same chain as a`) | The same for chains |
| `bymolecule a` (`same molecule as a`) | Every atom of a molecule instance containing an atom of `a` |
| `bonded to a` | Atoms bonded to an atom of `a` |

`within` uses a uniform spatial grid with cells of the cutoff size, so its
cost grows with the atoms near the reference, not with the product of the
two sets. Points a host supplies come as `select::ReferencePoints`, which
keep their grid and the neighbor lists it gathers, so a host evaluating one
selection over many structures builds them once and grids them once for
each cutoff. Repeated reference
points, as of copies stacked in one place, are measured once.
`SelectEvaluateOptions::with_executor` measures ranges of the evaluated atoms
as separate work items on a `parallel::Executor`, such as a host's thread
pool; every executor selects the same atoms.

## Host extensions

Applications extend the language without forking it:

- `#...` (a reference such as `#3`, `#2-5`, or `#1.2`) and `$name` (a named
  selection) always parse, and evaluate through the
  `select::SelectionResolver` the caller supplies.
- Bare keywords listed in `SelectParseOptions::with_host_keywords`, such as
  `displayed` or `selected`, parse as host terms too.

A resolver returns atoms of the evaluated topology, or points, for spatial
queries across objects: `within 4 of #2` evaluated on one object can measure
against another object's coordinates. Points are accepted only as the
reference of `within`. Evaluating a host term without a resolver fails with
`SelectEvaluateError::NoResolver`.

A `within` reference that is not a host term, such as `ligand` in
`within 5 of ligand`, selects atoms of the evaluated structure, and a resolver
may add points to it through `SelectionResolver::reference_points`. A host
that keeps a receptor and its ligand as two structures evaluates the
reference, itself a `SelectionExpr`, on the ligand structure and returns the
positions of the atoms it selects there, so `within 5 of ligand` evaluated on
the receptor selects the pocket. References nested in it reach the resolver
again during that evaluation. The points may include the atoms the reference
selects in the evaluated structure itself, so a host evaluating one selection
over many structures answers each reference with one `ReferencePoints` of
every structure, shared by all of them. Evaluation still reads one structure: others
contribute only the points the host supplies, and without them a reference
reads the evaluated structure alone.

## Canonical form

`SelectionExpr` displays its canonical text: lower-case keywords with aliases
replaced (`same residue as` becomes `byres`, `elem` becomes `element`, `=`
becomes `==`), single spaces, flattened `and`/`or` chains, parentheses only
where precedence requires them, and values quoted only when necessary.
Parsing the canonical text gives an equal expression, and printing that again
gives the same text.

## Errors

Parse errors are `diagnostics::FormatError` values of format
`Format::Selection`, kind `Syntax`, `InvalidValue`, or `Limit`, carrying the
byte span of the offending text (`FormatError::span`) and its line and
column. `SelectParseOptions` bounds the text length and the nesting depth.
Evaluation errors are `select::SelectEvaluateError`, which names the failing
term and its span.

## Conformance

The conformance suite in
[`tests/fixtures/select/conformance.txt`](../crates/kekule/tests/fixtures/select/conformance.txt)
lists expressions with their expected dense atom indices on reference
structures. Its expectations come from an independent per-atom reference
implementation in the test, not from the evaluator.
