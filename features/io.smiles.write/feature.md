# Noncanonical SMILES Writer

## Summary

Write small molecules as deterministic noncanonical SMILES for round-trip workflows.

## Behavior/API

- Exposes `SmilesWriteOptions` and `write_smiles`.
- Emits graph-order-based noncanonical SMILES with branches, ring closures, dot fragments, common bond symbols, and bracket atoms when needed.
- Emits `[nH]` for sanitized aromatic donor nitrogen when the perceived hydrogen must survive reparse.
- Preserves bracket-only no-implicit-hydrogen semantics.
- Rejects zero, dative, quadruple, atom/bond stereochemistry, radicals, and graphs requiring more than 99 ring labels instead of silently coercing them.
- Does not canonicalize or sanitize before writing.

## Implementation Notes

- The writer targets readability and deterministic output, not canonical ranking.
- A deterministic DFS tree is rendered with preassigned ring closures at both endpoints and branch children before the selected continuation path.
- Tree collection, subtree sizing, and component emission use explicit stacks so graph depth does not consume the Rust call stack.
- Unsupported advanced stereo/query details return structured write errors until later feature work can encode them faithfully.

## Validation

- Unit tests cover parse/write/parse round trips for branches, rings, brackets, fragments, aromatic examples, and unsupported lossy bond/stereo cases.
- RDKit-generated goldens compare sanitize/write/reparse atom identity, labeled-neighbor topology, bond order/aromaticity, charge, isotope, hydrogen, map, and valence records for external PubChem SMILES fixtures rather than exact RDKit noncanonical traversal strings.

## Out Of Scope

- Canonical SMILES, isomeric SMILES parity, SMARTS, reactions, and full stereochemical output.

## Revision Notes

- v1: Noncanonical writer.
- v2: Deterministic ring-closure and branch emission passes the RDKit-backed `tiny` corpus.
- v3: Make writer output self-readable for aromatic SMILES, preserve aromatic donor `[nH]`, and reject unencoded lossy bond/stereo representations.
- v4: Make graph-size-dependent writer traversals iterative while preserving deterministic output.
