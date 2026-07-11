# Core Atom and Bond Model

## Summary

Represent chemically general atom and bond data shared by small-molecule and macromolecular workflows.

## Behavior/API

- Provides `Atom`, `AtomRadical`, `Bond`, `Element`, `BondOrder`, property maps, and simple chemistry annotations.
- Stores radicals as one authoritative multiplicity (`Singlet`, `Doublet`,
  `Triplet`, `Quartet`, or `Quintet`) with a helper for unpaired-electron count.
- Stores whether an atom suppresses implicit hydrogen assignment, matching bracket-atom and explicit-valence workflows.
- Does not store authoritative atom or bond stereochemistry payload flags; stereo state lives in the graph-adjacent `core::stereo` model on `Molecule`.
- Atom and bond endpoint mutation stays controlled by `Molecule` topology operations.
- Mutable chemistry-relevant payload access conservatively invalidates perception state.

## Implementation Notes

- Element handling covers periodic-table symbols and atomic numbers.
- Atom fields stay chemically general and do not contain biomolecular hierarchy labels.
- Bond order storage is descriptive and does not imply valence validation, sanitization, or stereochemical assignment.

## Validation

- Current coverage is unit-test based.
- No RDKit or Biopython golden data is required for the representation feature.

## Out Of Scope

- Valence validation, aromaticity assignment, stereochemistry representation or perception, parsing, and validation generation.
- Runtime RDKit or Biopython dependency.

## Revision Notes

- v1: Chemically general atom and bond payload model.
- v2: Add an explicit `BondStereo::Any` value for raw file formats that preserve unspecified double-bond stereochemistry.
- v3: Replace lossy radical electron counts with authoritative radical multiplicity and re-export `AtomRadical`.
- v4: Add `Atom::no_implicit_hydrogens` so parsers can preserve explicit no-implicit-hydrogen atom semantics.
- v5: Remove atom/bond stereo payloads from the authoritative public model in favor of graph-adjacent stereo elements and source marks.
- v6: Extend authoritative radical multiplicity through quartet and quintet so
  bracket-SMILES high-spin states are not collapsed to triplets.
