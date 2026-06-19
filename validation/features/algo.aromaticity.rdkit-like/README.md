# algo.aromaticity.rdkit-like Validation Fixtures

Validation fixtures for this feature are externally supplied PubChem SDF records declared in
`tiny.toml`.

Future RDKit goldens should record atom and bond aromatic flags and explicitly mark cases outside
the RDKit-like model. The goal is stable pressure on common organic aromaticity, not full RDKit parity.
