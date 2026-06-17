# Implementation Plan

## Feature ID

`io.sdf.v2000.parse`

## Goal

Parse multi-record SDF V2000 input into `SmallMolecule` or `Molecule` values using the shared core graph. The parser should perform raw parsing and minimal structural loading only; sanitization, valence perception, ring detection, aromaticity perception, stereochemistry assignment, and canonicalization remain separate features.

## Public API

Add an SDF parsing API equivalent to:

- `SdfParseOptions`
- `SdfRecord`
- `SdfParseError`
- `read_sdf_v2000_str(input: &str, options: SdfParseOptions) -> Result<Vec<SmallMolecule>, SdfParseError>`
- Optional streaming API: `SdfV2000Reader<R: BufRead>`

Expected behavior:

- Parse multiple records separated by `$$$$`.
- Parse V2000 counts line.
- Parse atom block into `Atom` records.
- Parse bond block into core graph bonds.
- Attach SDF data fields to molecule properties.
- Preserve raw record title or header fields as molecule properties if the API documents names.
- Return record index and line context in parse errors.

Avoid adding APIs that imply sanitization or perception.

## Internal Modules Touched

Expected scope:

- `crates/molecules/src/lib.rs` initially, or an `io::sdf` module if module splitting is introduced.
- Unit tests under the parser module.
- Test fixtures under a local test fixture directory if useful.
- Feature docs under `features/io.sdf.v2000.parse/`.

Do not touch mmCIF parsing, ring detection, aromaticity, or validation generation except for validation fixture planning.

## Data Model

Parser-loaded molecules should use:

- Core `Molecule` graph for atoms and bonds.
- `Atom` fields for element, isotope, formal charge, explicit hydrogens, atom map, and raw carried annotations where V2000 fields support them.
- `Bond` fields for endpoints, `BondOrder`, stereo annotation, and properties.
- Molecule `PropMap` for SDF data fields.

Recommended property keys:

- `sdf.title`
- `sdf.program`
- `sdf.comment`
- `sdf.field.<name>` or the field name directly if the API chooses raw SDF field names.

The parser should preserve unknown or unsupported atom/bond annotations as properties rather than silently performing perception.

## Algorithm Outline

1. Split input into records at lines containing `$$$$`.
2. For each non-empty record, read the three header lines.
3. Read and validate the counts line.
4. Require V2000 format for this feature; reject V3000 with an explicit unsupported error.
5. Parse exactly the declared atom count of atom block lines.
6. Convert element symbols to `Element`; reject unknown symbols.
7. Parse exactly the declared bond count of bond block lines.
8. Convert atom indices from one-based SDF indices to `AtomId`s created during atom insertion.
9. Convert bond order codes to `BondOrder` storage values.
10. Read property and metadata lines until `M  END`.
11. Parse data fields after `M  END` until record end.
12. Return parsed molecules in input order.

All parsing should be deterministic and should distinguish raw parsing errors from future sanitization errors.

## Tests

Add tests for:

- Empty input returns no molecules or a documented empty-input error.
- Single minimal V2000 molecule parses.
- Multiple records parse in order.
- Atom count and bond count match counts line.
- One-based bond endpoints map to stable `AtomId`s.
- Single, double, triple, aromatic, zero, and dative-supported mappings behave as documented.
- Unknown element symbol fails.
- Bad counts line fails with record and line context.
- Bond endpoint outside atom range fails.
- V3000 input fails as unsupported.
- SDF data fields are attached as molecule properties.
- Blank data field values are preserved.
- Parser does not set perceived rings, aromaticity, valence, or stereo state to fresh.
- Parser does not call sanitization.

## Reference Validation

Use RDKit only as a reference tool in validation infrastructure.

Golden JSON should include:

- RDKit version.
- Input fixture name and checksum.
- Molecule count.
- Per-record atom count and bond count.
- Element atomic numbers.
- Bond endpoint index pairs and bond order labels.
- SDF data fields.
- Unsupported or intentionally ignored fields.

Do not copy RDKit parsing code. Use RDKit only to generate comparable expected data.

## Required Checks

Run:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
cargo xtask validate --feature io.sdf.v2000.parse
```

## Risks

- V2000 fixed-width fields are easy to parse incorrectly with simple splitting.
- Some SDF files contain malformed but common variants.
- Mapping MDL charge and isotope records can blur parsing and sanitization if not documented.
- Aromatic bond code handling can be confused with aromaticity perception.
- Data field names may contain characters awkward for property keys.

## Edge Cases

- Missing `M  END`.
- Missing final `$$$$`.
- Empty records between delimiters.
- Counts line declares more atoms or bonds than present.
- Atom coordinates are present but conformer storage may not exist yet.
- Multi-line data fields.
- Duplicate data field names.
- Bond order codes outside supported values.
- Atom list/query features unsupported in this first parser.

## Explicitly Out of Scope

- SDF V3000 parsing.
- Molfile writing.
- Sanitization.
- Valence perception.
- Ring detection.
- Aromaticity perception.
- Stereochemistry assignment.
- Canonicalization.
- Runtime RDKit dependency.
