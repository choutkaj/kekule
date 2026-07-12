# mmCIF Molecular Interpretation

## Summary

Decipher a parsed mmCIF structure document into clean molecule objects without treating the entire experimental structure as one graph.

## Behavior/API

- Exposes `mmcif::interpret` with explicit options, structured errors, and an interpretation report.
- Selects exactly one coordinate-containing data block per call.
- Uses `_entity`, `_struct_asym`, and `_atom_site` metadata to classify and partition structural instances.
- Produces polymer and branched instances as `MacroMolecule`, non-polymers and ions as `SmallMolecule`, and each water as an individual solvent molecule.
- Deduplicates atom identity across coordinate models and stores coordinates as conformers carrying the source model ID.
- Applies an explicit alternate-location policy; the default selects highest occupancy deterministically.
- Declared covalent `_struct_conn` links merge molecular instances and add single bonds. Non-covalent, unresolved, and unsupported links are reported without topology changes.
- Strict mode rejects missing entity classification; default mode makes conservative inferences and reports them.
- Never sanitizes, prepares, adds hydrogens, or infers missing template bonds.

## Implementation Notes

- Polymer/branched structural instances establish provisional molecule boundaries while CCD and polymer template connectivity remain unavailable.
- Non-polymer and water occurrences are separated by residue identity, including conservative occurrence grouping when sequence identifiers are absent.
- Source atom/component/asym identifiers remain on atom and graph properties; macromolecular hierarchy metadata remains in `BioHierarchy`.
- The report counts coordinate models, output categories, applied connections, inferred classifications, ignored/unresolved connections, and multi-atom graphs still awaiting template bonds.
- The parse-then-interpret pipeline is the only mmCIF molecular import path; no direct whole-file molecule reader or compatibility alias exists.

## Validation

- Unit tests cover mixed polymers, ligands, ions, repeated waters, multiple coordinate models, alternate locations, entity inference, strict metadata, and declared covalent links. Successful fuzzed documents are also passed through the public interpretation boundary.
- A downstream-style integration test compiles the public parse-then-interpret workflow.
- No corpus-backed molecule-boundary evidence exists yet, so the feature remains unvalidated.

## Out Of Scope

- CCD/template bond lookup, polymer bond construction, ligand sanitization, PDB input, coordination chemistry, assembly generation, model conversion, preparation, and serialization.

## Revision Notes

- v1: Add conservative staged interpretation of structural mmCIF documents into clean molecular contents.
- v2: Remove the historical direct reader and make parse-then-interpret the exclusive mmCIF molecular import path.
