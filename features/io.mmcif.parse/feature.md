# mmCIF Document Parser

## Summary

Parse a structural mmCIF file into a format-level document before chemical interpretation.

## Behavior/API

- Exposes only `mmcif::parse_str` as the mmCIF reader, returning `MmcifDocument`.
- Preserves every data block, scalar item, loop table, unknown category, missing-value marker, and value source line.
- Offers case-insensitive lookup while preserving original data names and values.
- Distinguishes quoted values from bare syntax controls, preserves `#` inside
  bare values, recognizes comments only where a token may begin, and supports
  explicit `stop_` loop terminators.
- Rejects content outside a data block, unnamed blocks, missing scalar values, duplicate data names, duplicate loop tags, and ragged loop rows with structured parse errors.
- Bounds input bytes, token count, token bytes, and atom-site rows through `MmcifParseOptions`.
- Does not expose a direct whole-file `MacroMolecule` reader or any compatibility alias.

## Implementation Notes

- The document layer owns strings and has no live dependency on the input buffer.
- Parsing does not infer entities, molecule boundaries, bonds, coordinate models, hierarchy, or chemistry.
- Molecular meaning is assigned only by the separate `mmcif::interpret` stage.

## Validation

- Unit tests cover unknown scalar/loop content, quoted controls, multiline values, multiple blocks, lookup, malformed structure, resource limits, and deterministic mutation safety. The parser fuzz target exercises the public document API and bounded parse options.
- Biopython `MMCIF2Dict` goldens compare the nineteen canonical `_atom_site`
  columns as format-level strings or missing values without asserting molecular
  interpretation. PDB-100 is required baseline evidence and PDB-1000 is the
  broad manifest-backed tier.
- The former Biopython structure summaries described the removed direct
  `MacroMolecule` reader and are intentionally excluded from this document contract.

## Out Of Scope

- Dictionary save frames, DDL validation, category schemas, binary CIF, writing, and chemical interpretation.

## Revision Notes

- v1-v7: Historical direct atom-site-to-`MacroMolecule` reader.
- v8: Hard-break the historical reader and redefine the canonical parser as format-level `MmcifDocument` construction only.
- v9: Add Biopython-backed PDB-100 atom-site document parity as required baseline validation.
- v10: Preserve embedded hash characters in bare values while retaining token-boundary comments.
