# Roadmap

## P0: repository and core foundation

- Cargo workspace and CI.
- Minimal `molecules` crate skeleton.
- Feature registry, dashboard, and templates.
- Codex skills for research, planning, implementation, and review.
- Core atom/bond/graph/property/conformer data structures.

## P0: first chemistry features

- Molfile V2000 parse.
- SDF V2000 parse and data field handling.
- Basic valence and sanitization pipeline.
- Fast ring membership detection.
- SSSR ring basis.
- Basic RDKit-like aromaticity model.

## P0: first macromolecule features

- PDB parse.
- mmCIF parse.
- SMCRA-like hierarchy.
- Atom-site labels and alternate locations.
- Basic sequence extraction and selection API.

## P1: practical usability

- Molfile/SDF writing.
- SMILES parse and noncanonical write.
- Tetrahedral and double-bond stereochemistry representation.
- Kabsch superposition and geometry helpers.
- Component dictionary support for common biomolecular residue connectivity.

## P2: deeper cheminformatics

- V3000 parse/write.
- Canonical ranking and canonical SMILES.
- CIP stereochemistry.
- SMARTS-like substructure search.
- Fingerprints and basic descriptors.

## P3: advanced macromolecular support

- Biological assemblies.
- Symmetry records.
- Secondary structure annotations.
- Nucleic acid sequence support.
- Ligand extraction and polymer connectivity refinement.
