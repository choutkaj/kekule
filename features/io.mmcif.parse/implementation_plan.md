# Implementation Plan

## Feature ID

`io.mmcif.parse`

## Goal

Parse mmCIF atom-site data into `MacroMolecule`, preserving biomolecular hierarchy labels in `BioHierarchy` and chemistry topology in the core `Molecule`. Parsing must remain separate from sanitization and chemical perception.

The first parser should prioritize `_atom_site` records and enough metadata to populate SMCRA hierarchy.

## Public API

Add an API equivalent to:

- `MmcifParseOptions`
- `MmcifParseError`
- `read_mmcif_str(input: &str, options: MmcifParseOptions) -> Result<MacroMolecule, MmcifParseError>`
- Optional streaming or file helper later.

Expected options:

- strict vs permissive missing-field handling
- whether to include hydrogens if present
- whether to create core bonds from `_struct_conn` when available; default can be no bond perception

Expected behavior:

- Parse mmCIF loops for `_atom_site`.
- Create core atoms from element symbols.
- Populate `BioHierarchy` models, chains, residues, and atom sites.
- Preserve label and author identifiers separately.
- Preserve alternate location, occupancy, B-factor, and coordinates if conformer storage exists; otherwise store coordinates as planned metadata or leave documented TODO.

## Internal Modules Touched

Expected scope:

- mmCIF parser module or `crates/molecules/src/lib.rs` while small.
- `MacroMolecule` and `BioHierarchy` APIs from `bio.hierarchy.smcra`.
- Optional parser fixtures.
- Feature docs under `features/io.mmcif.parse/`.

Do not implement aromaticity, ring detection, valence perception, or full polymer bond perception.

## Data Model

Input data should map as follows:

- `_atom_site.type_symbol` -> core `Atom.element`
- `_atom_site.label_atom_id` -> `BioHierarchy` atom-site label atom ID
- `_atom_site.auth_atom_id` -> `BioHierarchy` atom-site author atom ID
- `_atom_site.label_alt_id` -> alternate location
- `_atom_site.label_comp_id` and `_atom_site.auth_comp_id` -> residue labels
- `_atom_site.label_asym_id` and `_atom_site.auth_asym_id` -> chain labels
- `_atom_site.label_seq_id` and `_atom_site.auth_seq_id` -> residue sequence labels
- `_atom_site.pdbx_PDB_ins_code` -> insertion code
- `_atom_site.pdbx_PDB_model_num` -> model
- `_atom_site.occupancy` -> atom-site occupancy
- `_atom_site.B_iso_or_equiv` -> atom-site B-factor
- coordinates -> conformer storage when available, otherwise documented holding field or out-of-scope note

Unknown categories should be ignored or preserved as properties only if documented.

## Algorithm Outline

1. Tokenize mmCIF respecting quoted strings, semicolon text blocks, comments, `?`, and `.` missing values.
2. Parse data blocks and loops.
3. Locate `_atom_site` loop.
4. Validate required atom-site columns for the selected strictness.
5. Iterate atom-site rows in file order.
6. Create or reuse model, chain, and residue hierarchy records based on label/auth IDs.
7. Create a core atom from `type_symbol`.
8. Attach the new `AtomId` to the residue with atom-site metadata.
9. Store coordinates if supported by the current data model; otherwise document that coordinate attachment is deferred.
10. Optionally parse `_struct_conn` only if explicitly in scope for a later increment.
11. Return `MacroMolecule` with populated core graph and hierarchy.

The parser should not infer missing bonds by distance or run perception.

## Tests

Add tests for:

- Minimal mmCIF atom-site loop parses.
- Multiple models are preserved.
- Label and author chain IDs are preserved separately.
- Label and author residue IDs are preserved separately.
- Insertion code is preserved.
- Alternate location is preserved.
- Occupancy and B-factor parse.
- Missing optional fields are represented as `None`.
- Missing required fields fail in strict mode.
- `?` and `.` missing values are handled.
- Quoted values and semicolon text blocks tokenize correctly.
- Unknown element fails with parse error.
- Parser does not run sanitization or perception.
- Core `Atom` does not receive biomolecular hierarchy labels.

## Reference Validation

Use Biopython through `validation.harness`.

Golden JSON should include:

- Biopython version.
- Fixture name and checksum.
- Model, chain, residue, and atom-site hierarchy.
- Label and author IDs.
- Atom element symbols.
- Coordinates if implemented.
- Occupancy, B-factor, altloc, and insertion code.

Validation should focus on hierarchy preservation, not chemical bond perception.

## Required Checks

Run:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
cargo xtask validate --feature io.mmcif.parse
```

## Risks

- mmCIF tokenization is deceptively complex; ad hoc whitespace splitting will fail.
- Label and author identifiers can be confused.
- Some files omit columns that are common but not guaranteed.
- Coordinates need a clear storage destination.
- Inferring bonds from macromolecular data can become a large chemistry feature if not deferred.

## Edge Cases

- Multiple data blocks.
- Loop columns in unusual order.
- Missing values represented by `?` or `.`.
- Quoted strings with spaces.
- Semicolon-delimited multiline values.
- Alternate locations.
- Multiple models.
- Non-polymer ligands and waters.
- Negative residue numbers or insertion codes.
- Unknown or nonstandard element symbols.

## Explicitly Out of Scope

- Full mmCIF category coverage.
- PDB parsing.
- Distance-based bond perception.
- Polymer template bond assignment.
- Ring detection.
- Aromaticity.
- Valence perception.
- Stereochemistry.
- Runtime Biopython dependency.
