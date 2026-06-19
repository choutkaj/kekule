# algo.rings.fast Validation Fixtures

Validation fixtures for this feature are externally supplied PubChem SDF records declared in
`tiny.toml`.

Future RDKit goldens should compare only atom and bond membership in at least one cycle. They should
not compare SSSR, ring-family, or ring-ordering behavior.
