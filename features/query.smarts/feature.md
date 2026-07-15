# Bounded SMARTS Query Parsing

## Summary

Translate a deliberately bounded SMARTS subset into `query::QueryGraph` without sharing implementation or dependency direction with the SMILES parser.

## Behavior/API

- `query::parse_smarts` parses with conservative defaults;
  `parse_smarts_with_options` exposes input, topology, branch, ring-closure,
  expression-node, and expression-depth limits.
- Supported atom syntax includes `*`, atomic-number queries, isotopes, all
  bracketed element symbols, organic-subset unbracketed elements, aromatic
  `b/c/n/o/p/s/as/se`, `A`/`a`, formal charge, `D`, `H`, `R`/`R0`, bare `r`,
  negation, high-precedence `&`, `,`, low-precedence `;`, and implicit AND.
- Supported topology syntax includes branches, one- and two-digit ring labels,
  and ungrouped dot-disconnected query components.
- Supported bonds are omitted single-or-aromatic, `-`, `=`, `#`, `$`, `:`,
  `~`, `@`, and `!@`.
- `SmartsParseError` reports `Empty`, `InvalidSyntax`, `Unsupported`, or
  `ResourceLimit` with a byte span and message.

## Implementation Notes

- SMARTS is a frontend only. It constructs `QueryGraph` and has no dependency
  on substructure matching or `io.smiles.parse`.
- Operator precedence follows Daylight SMARTS: `!`, high-precedence AND,
  comma OR, then low-precedence `;` AND. Omitted bonds mean single or aromatic.
- Unsupported constructs are rejected rather than approximated. These include
  recursive SMARTS, atom-expression parentheses, stereochemistry, atom maps,
  exact ring counts other than `R0`, ring sizes, connectivity, ring-bond-count,
  valence, and hybridization primitives, bond-expression boolean algebra, and
  component-level grouping.
- Default limits are 16,384 input bytes, 256 atoms, 512 bonds, branch depth 64,
  128 ring closures, 512 nodes per atom expression, and expression depth 32.

## Validation

- Regression tests cover supported topology and logic, precedence, structured
  spans, unsupported constructs, malformed inputs, and every resource class.
- RDKit goldens parse externally supplied non-stereochemical PubChem SMILES as
  valid SMARTS and compare parse status and query atom/bond topology.

## Out Of Scope

Full SMARTS compatibility, recursive queries, stereochemical matching,
component grouping, CXSMARTS, reaction SMARTS, query serialization, and silent
fallback for unsupported primitives.

## Revision Notes

- v1: Feature contract reserved.
- v2: Implement the bounded parser as a one-way frontend for `query.graph`.
