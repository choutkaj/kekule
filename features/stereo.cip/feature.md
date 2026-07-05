# CIP Stereochemistry

## Summary

Plan CIP assignment for later stereochemical workflows.

## Behavior/API

No public API is implemented yet. CIP labels will be derived descriptors over
validated stereo elements, not stored graph truth.

## Implementation Notes

This feature should depend on explicit stereo representation, local stereo
perception/validation, sanitized valence and hydrogen semantics, and
deterministic ranking helpers. Exact assignment should follow the
machine-oriented CIP model with bounded exploration and explicit resource
limits for highly symmetric graphs.

## Validation

Future validation should compare manually reviewed and RDKit-backed stereochemical fixtures.

## Out Of Scope

Current first-wave implementation.

## Revision Notes

- v1: Feature contract reserved.
- v2: Reframe CIP as a derived-cache layer over representation and perception,
  with deterministic ranking and sanitized chemistry as dependencies.
