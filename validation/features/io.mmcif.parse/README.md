# io.mmcif.parse Validation Fixtures

Validation fixtures for this feature are externally supplied RCSB PDB mmCIF records declared in
`tiny.toml`.

Malformed parser inputs belong in Rust unit tests. Reference validation uses external structure
records plus Biopython-generated normalized goldens.
