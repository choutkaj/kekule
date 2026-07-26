# Noncanonical SMILES Writer

## Summary

Write small molecules as deterministic noncanonical SMILES for round-trip workflows.

## Behavior/API

- Exposes `smiles::{SmilesWriteOptions, write, write_with_options}`.
- Emits graph-order-based noncanonical SMILES with branches, ring closures, dot fragments, common bond symbols, and bracket atoms when needed.
- Emits `[nH]` for sanitized aromatic donor nitrogen when the perceived hydrogen must survive reparse.
- Preserves bracket-only no-implicit-hydrogen semantics.
- Rejects zero, dative, quadruple, stored stereo elements, source bond stereo
  marks, and graphs requiring more than 99 ring labels instead of silently
  coercing them. Radical atoms are writeable only when reparsing the emitted
  bracket atom will infer the same multiplicity from its valence; arbitrary
  stored radical/valence combinations are rejected.
- Does not canonicalize or sanitize before writing.

## Implementation Notes

- The writer targets readability and deterministic output, not canonical ranking.
- A deterministic DFS tree is rendered with preassigned ring closures at both endpoints and branch children before the selected continuation path.
- One-based ring-closure labels are assigned in `u64` after enforcing the
  supported maximum of 99 simultaneous closures.
- Tree collection, subtree sizing, and component emission use explicit stacks so graph depth does not consume the Rust call stack.
- Unsupported stereo/query details are read from the first-class stereo representation and return structured write errors until isomeric SMILES support can encode them faithfully.
- Representable radicals need no nonstandard token: bracket atom syntax carries
  the valence state and the parser reconstructs doublet through quintet
  multiplicity deterministically.

## Tests

- Unit tests cover parse/write/parse round trips for branches, rings, brackets, fragments, aromatic examples, and unsupported lossy bond/stereo cases from the graph-adjacent stereo model.

## Benchmarks

- RDKit-generated goldens compare sanitize/write/reparse atom identity, labeled-neighbor topology, bond order/aromaticity, charge, isotope, hydrogen, map, and valence records for external PubChem SMILES fixtures rather than exact RDKit noncanonical traversal strings.
- Optional external-reference manifests are available for `pubchem-1k`, `pubchem-100k`, `enamine-diversity`.
- Benchmark observations are informational and never determine this feature's release status or repository health.

## Out Of Scope

- Canonical SMILES, isomeric SMILES parity, SMARTS, reactions, and full stereochemical output.

## Revision Notes

- v1: Noncanonical writer.
- v2: Deterministic ring-closure and branch emission passes the RDKit-backed `smoke` corpus.
- v3: Make writer output self-readable for aromatic SMILES, preserve aromatic donor `[nH]`, and reject unencoded lossy bond/stereo representations.
- v4: Make graph-size-dependent writer traversals iterative while preserving deterministic output.
- v5: Move the public noncanonical writer API under the `smiles` facade.
- v6: Add PubChem-100k as required broad-corpus external-parity evidence.
- v7: Reject first-class stereo elements and source bond marks instead of reading removed atom/bond payload flags.
- v8: Round-trip valence-implied bracket radicals through quintet while
  rejecting radical multiplicities that the emitted atom valence cannot encode.
- v9: Keep every ignored non-smoke corpus as explicit local-only validation
  instead of repository-wide required evidence.
- v10: Use PubChem-1k as the required baseline benchmark corpus after retiring the former smoke corpus from public validation.
- v11: Reclassify external-reference parity from a required gate to optional benchmarking without changing implementation behavior or golden expectations.
- v12: Generate one-based ring labels in `u64` so formatting never relies on
  overflow-prone collection-index arithmetic.
