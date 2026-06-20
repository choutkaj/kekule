# Core Atom and Bond Model

## Summary

Represent chemically general atom and bond data shared by small-molecule and macromolecular workflows.

## Behavior/API

- Provides `Atom`, `AtomRadical`, `Bond`, `Element`, `BondOrder`, `BondStereo`, property maps, and simple chemistry annotations.
- Stores radicals as one authoritative multiplicity (`Singlet`, `Doublet`, or `Triplet`) with a helper for unpaired-electron count.
- Atom and bond endpoint mutation stays controlled by `Molecule` topology operations.
- Mutable chemistry-relevant payload access conservatively invalidates perception state.

## Implementation Notes

- Element handling covers periodic-table symbols and atomic numbers.
- Atom fields stay chemically general and do not contain biomolecular hierarchy labels.
- Bond order storage is descriptive and does not imply valence validation or sanitization.

## Validation

- Current coverage is unit-test based.
- No RDKit or Biopython golden data is required for the representation feature.

## Out Of Scope

- Valence validation, aromaticity assignment, stereochemistry perception, parsing, and validation generation.
- Runtime RDKit or Biopython dependency.

## Revision Notes

- v1: Chemically general atom and bond payload model.
- v2: Add an explicit `BondStereo::Any` value for raw file formats that preserve unspecified double-bond stereochemistry.
- v3: Replace lossy radical electron counts with authoritative radical multiplicity and re-export `AtomRadical`.
