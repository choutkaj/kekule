# Syntax-Independent Query Graph

## Summary

Represent molecular graph queries independently of any input syntax and independently of the concrete `Molecule` kernel.

## Behavior/API

- `query::QueryExpression<P>` represents constants, primitives, negation,
  conjunction, and disjunction. `AtomExpression` and `BondExpression` bind it
  to the supported atom and bond primitive enums.
- `evaluate_with` lets syntax-independent downstream consumers supply primitive
  semantics; `contains_predicate` supports capability and prerequisite discovery.
- Atom primitives cover element, isotope, formal charge, aromaticity, explicit
  graph degree, total attached hydrogens, and ring membership.
- Bond primitives cover concrete bond order, aromaticity, and ring membership.
- `QueryGraphBuilder` creates a non-empty immutable `QueryGraph` with dense
  `QueryAtomId` and `QueryBondId` values. Self-bonds, duplicate bonds, and
  invalid endpoints are rejected.
- A query graph may be disconnected. It carries no source text, SMARTS token,
  parser span, or matcher state.

## Implementation Notes

- The representation is separate from `Molecule`; query-only alternatives and
  negations never become optional fields on concrete `Atom` or `Bond` values.
- Boolean constructors flatten same-kind operators, remove identity constants,
  and preserve a maximum of 4,096 nodes and 64 levels. Private expression-tree
  fields prevent construction of values that violate those bounds.
- Query topology is immutable after `QueryGraphBuilder::build`, which keeps
  adjacency and endpoint invariants local and makes one graph reusable across
  match operations.
- The IR deliberately describes chemical facts, not how a frontend spells
  them. Other query syntaxes can target the same representation.

## Validation

- Unit regressions cover expression normalization and hard depth bounds, query
  IDs and adjacency, duplicate/self-bond rejection, and disconnected graphs.
- Syntax translation and chemical matching receive separate RDKit-backed
  validation under `query.smarts` and `algo.substructure.vf2`.

## Out Of Scope

Recursive subqueries, stereochemical predicates, atom-map capture semantics,
property predicates, mutable query graphs, and a general query optimizer.

## Revision Notes

- v1: Introduce the bounded syntax-neutral expression tree and immutable query graph.
