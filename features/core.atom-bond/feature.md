# Core Atom and Bond Model

## Summary

Represent chemically general atom and bond data shared by small-molecule and macromolecular workflows.

## Behavior/API

- Provides `Atom`, `Bond`, `Element`, `BondOrder`, property maps, and simple chemistry annotations.
- Atom and bond endpoint mutation stays controlled by `Molecule` topology operations.
- Mutable chemistry-relevant payload access conservatively invalidates perception state.

## Implementation Notes

- Element handling covers periodic-table symbols and atomic numbers.
- Atom fields stay chemically general and do not contain biomolecular hierarchy labels.
- Bond order storage is descriptive and does not imply valence validation or sanitization.

## Validation

- Current coverage is unit-test based.
- No RDKit or Biopython golden data is required for the initial representation feature.

## Out Of Scope

- Valence validation, aromaticity assignment, stereochemistry perception, parsing, and validation generation.
- Runtime RDKit or Biopython dependency.

## Revision Notes

- v1: Initial chemically general atom and bond payload model.
