# Bounded VF2-Style Substructure Search

## Summary

Match a syntax-independent `QueryGraph` against a concrete `Molecule` with deterministic traversal, explicit perception prerequisites, uniqueness control, and hard work bounds.

## Behavior/API

- `substructure::find_substructure_match` returns the first embedding;
  `find_substructure_matches` returns unique target atom sets; and
  `find_substructure_matches_with_options` exposes limits and uniqueness.
- `QueryMatch` stores target `AtomId` values in query-atom order and supports
  lookup by `QueryAtomId`.
- Matching is injective and non-induced: every query edge must match, while
  additional target edges are allowed. Disconnected query components may match
  within one target component, consistent with ungrouped SMARTS dots.
- Unique mode removes query-automorphism duplicates by target atom set.
- Element, isotope, charge, graph degree, hydrogen count, aromaticity, ring
  membership, bond order, aromatic bond, and ring-bond predicates are supported.

## Implementation Notes

- The matcher depends on `query.graph`, never on SMARTS. Programmatic queries
  and future parser frontends use the same path.
- Candidate compatibility is precomputed, the next query atom maximizes mapped
  neighbors and constraint tightness, mapped edges are checked immediately,
  and unmapped neighbors receive forward feasibility checks.
- Matching never runs perception. Queries using total hydrogen count,
  ring membership, or aromaticity require current valence, ring, or aromaticity
  state respectively and return `MissingPerception` otherwise.
- Default limits are 256 query atoms, 2,000,000 candidate pairs, 1,000,000
  visited search states, and 1,000 returned matches. `max_matches` intentionally
  caps enumeration in the same style as RDKit; search-work limits are errors.
- Because backtracking depth is one frame per query atom, 256 is also an
  absolute stack-safety ceiling; callers may lower but not raise it.
- Target and candidate iteration follow stable IDs, giving deterministic output.

## Validation

- Unit regressions cover programmatic IR matching, SMARTS-produced queries,
  logical precedence, non-induced and disconnected matching, automorphism
  uniqueness, perception prerequisites, and hard candidate/search-state limits.
- RDKit goldens compare normalized unique target atom sets for a fixed bounded
  query suite over externally supplied PubChem SDF molecules.

## Out Of Scope

Induced matching, maximum common substructure, recursive queries,
stereochemical matching, tautomer or resonance equivalence, fingerprint
prefilters, parallel search, and motif-specific chemical exceptions.

## Revision Notes

- v1: Feature contract reserved.
- v2: Implement deterministic bounded matching on `query.graph` and remove the
  dependency on the SMARTS frontend.
