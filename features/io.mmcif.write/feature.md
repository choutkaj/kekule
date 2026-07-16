# mmCIF Molecular Model Writer

## Summary

Write a supported canonical `Model` as one structural PDBx/mmCIF data
block without hiding validation, perception, sanitization, or preparation.

## Behavior/API

- Exposes `mmcif::write(model, MmcifWriteOptions)` and returns structured
  `MmcifWriteError` failures.
- Emits deterministic `_entity`, `_struct_asym`, `_atom_site`, and optional
  `_struct_conn` loops in model topology order.
- Reads coordinates only from the authoritative `Model` position array;
  `coordinate_precision` controls fixed-decimal output and defaults to three
  decimal places.
- Explicitly converts model length quantities to the mmCIF Cartesian angstrom
  convention before formatting.
- Preserves supported macro hierarchy labels, residue identifiers, selected
  alternate-location metadata, occupancy, B-factor, element, formal charge, and
  author identifiers. Small-molecule instances receive deterministic structural
  labels required by mmCIF.
- Emits explicit single, double, triple, and quadruple bonds through
  `_struct_conn.pdbx_value_order`.
- Preserves `Polymer`, `Branched`, `NonPolymer`, and `Solvent` as mmCIF entity
  kinds. The `Ion` role must agree with the interpreter's single-charged-atom
  rule; unsupported `Ligand` and `Cofactor` roles fail explicitly.
- Rejects invalid/incomplete hierarchies, duplicate structural-instance IDs,
  ambiguous connection selectors, unrepresentable disconnected multi-asym
  molecule boundaries, unsupported atom/stereo fields, formal charges outside
  the PDBx range, and zero/aromatic/dative bond orders.
- Does not mutate, sanitize, canonicalize, perceive, or prepare the input model.

## Implementation Notes

- The writer targets `Model`, not `MacroMolecule`, because model
  positions are complete and authoritative and one structure may contain
  multiple Small and Macro molecule instances.
- Macro atom sites remain a `SmcraHierarchy` sidecar over local `AtomId`s; writer
  rows qualify them only while resolving model positions and connectivity.
- Values are emitted as single CIF tokens. Source formatting, comments,
  unknown categories, and original atom-site row IDs belong to `MmcifDocument`
  rather than canonical-model writing.
- Work and temporary memory are linear in the number of atoms, bonds, and
  hierarchy records. Output is currently accumulated in one `String`.

## Validation

- Unit tests cover public-facade writing, writer-output parsing, supported model
  round trips, all four supported bond orders, explicit unknown-order rejection,
  unsupported aromatic connectivity, missing atom-site rejection, duplicate atom
  identities, and unencodable model roles.
- No external writer golden is currently accepted by the validation harness, so
  no external writer parity result is recorded despite targeted unit regression
  coverage.

## Out Of Scope

- Loss-preserving `MmcifDocument` serialization and source-format round trips.
- Multiple coordinate models, unselected alternate locations, assemblies,
  crystallographic cell/symmetry data, anisotropic displacement, sequences,
  secondary structure, and chemical-component dictionaries.
- Aromatic/dative/zero-order connectivity, atom stereo, isotopes, radicals,
  atom maps, implicit-hydrogen declarations, ligand/cofactor role encoding, and
  streaming output.

## Revision Notes

- v1: Add foundational canonical-model writing with explicit support boundaries
  and order-preserving covalent connectivity.
- v2: Accept the renamed canonical `Model` and `SmcraHierarchy` APIs without
  changing emitted mmCIF semantics.
- v3: Convert explicit model length quantities to the mmCIF Cartesian angstrom
  convention before serialization.
