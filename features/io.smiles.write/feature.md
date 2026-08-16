# Noncanonical SMILES Writer

## Summary

Write one connected small molecule as deterministic noncanonical SMILES for
round-trip workflows.

## Behavior/API

- Exposes the single noncanonical entry point `smiles::write` and the
  `SmallMolecule::to_smiles()` convenience; there is no empty options type or
  redundant options overload.
- Emits graph-order-based noncanonical SMILES with branches, ring closures,
  common bond symbols, and bracket atoms when needed. Multi-component documents
  are interpreted and written component-by-component by callers.
- Emits `[nH]` for aromatic donor nitrogen from either a represented explicit
  hydrogen or one perceived implicit hydrogen, without requiring perception to
  rewrite the atom payload.
- Preserves bracket-only no-implicit-hydrogen semantics.
- Rejects zero, dative, quadruple, stored stereo elements, source bond stereo
  marks, radical atoms, and graphs requiring more than 99 ring labels instead
  of silently coercing them. The supported SMILES grammar has no explicit
  radical token, so no emitted bracket spelling can promise a lossless radical
  round trip.
- Does not canonicalize, normalize, or perceive before writing.

## Implementation Notes

- The writer targets readability and deterministic output, not canonical ranking.
- A deterministic DFS tree is rendered with preassigned ring closures at both endpoints and branch children before the selected continuation path.
- One-based ring-closure labels are assigned in `u64` after enforcing the
  supported maximum of 99 simultaneous closures.
- Tree collection, subtree sizing, and molecule emission use explicit stacks so
  graph depth does not consume the Rust call stack.
- Unsupported stereo/query details are read from the first-class stereo representation and return structured write errors until isomeric SMILES support can encode them faithfully.
- Radical output remains unsupported until the parser and writer share an
  explicit source-semantic radical representation.

## Tests

- Unit tests cover parse/write/parse round trips for branches, rings, brackets,
  individually interpreted components, aromatic examples, and unsupported
  lossy bond/stereo cases from the graph-adjacent stereo model.
- Pyrrole regressions cover both represented `[nH]` and an aromatic donor whose
  hydrogen remains derived in `PerceptionState` after default perception.

## Benchmarks

- RDKit-generated goldens compare normalize/perceive/write/reparse atom
  identity, labeled-neighbor topology, bond order/aromaticity, charge, isotope,
  hydrogen, map, and valence records for external PubChem SMILES fixtures
  rather than exact RDKit noncanonical traversal strings.
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
- v13: Align writer input and multi-component workflows with the connected
  `SmallMolecule` boundary.
- v14: Preserve aromatic donor `[nH]` output from the explicit represented-H
  plus perceived implicit-H split after removing sanitizer representation
  feedback.
- v15: Reject radical atoms explicitly after removing model-driven bracket
  radical inference from SMILES interpretation.
- v16: Remove the empty writer-options type and redundant options overload.
