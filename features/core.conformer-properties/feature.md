# Conformer Properties

## Summary

Attach arbitrary metadata to a conformer without placing coordinate-set identity on the molecular graph.

## Behavior/API

- `Conformer::props` exposes read-only arbitrary properties.
- `Conformer::props_mut` supports explicit metadata updates.
- New conformers start with an empty property map.

## Implementation Notes

- Properties reuse the chemically general core `PropMap` and `PropValue` types.
- Format-specific coordinate-model identifiers and provenance remain in
  interpretation reports rather than this generic property map.

## Validation

- Unit tests verify generic conformer property storage independently of format
  provenance.
- No external validation evidence is required for the generic metadata container.

## Out Of Scope

- A controlled property schema, serialization, provenance graphs, and automatic property invalidation.

## Revision Notes

- v1: Add arbitrary conformer properties.
- v2: Remove mmCIF coordinate-model provenance from generic conformer
  properties and keep it in qualified interpretation reports.
