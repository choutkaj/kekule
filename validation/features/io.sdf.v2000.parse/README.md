# io.sdf.v2000.parse Validation Fixtures

Validation fixtures for this feature are externally supplied PubChem SDF records declared in
`tiny.toml`.

Malformed parser inputs belong in Rust unit tests. Reference validation uses external molecule
records plus RDKit-generated normalized goldens.
