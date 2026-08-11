# SMILES Document Parsing and Interpretation

## Summary

Parse SMILES syntax into a loss-preserving `SmilesDocument`, then explicitly
interpret each dot-separated component as one connected `SmallMolecule` with
component-qualified source mappings.

## Behavior/API

- Exposes `smiles::parse_str(&str) -> SmilesDocument`,
  `smiles::parse_str_with_options(&str, SmilesParseOptions)`, and
  `smiles::interpret(&SmilesDocument) -> SmilesInterpretation` with distinct
  errors.
- Preserves source text, tokens and UTF-8 spans, branch/ring/stereo marks, and
  dot-component token ranges. Parsing also builds and validates the private
  syntax program consumed by interpretation.
- Each `SmilesComponentInterpretation` contains one connected molecule and a
  `SmilesInterpretationReport` that maps the component's parsed atom and bond
  records, including source spans, to stable canonical IDs.
- Interpretation supports the established atom, bond, branch, ring, charge,
  isotope, map, radical, aromatic-token, and local-stereo subset.
- Dots delimit distinct asserted molecules. `[Na+].[Cl-]` produces two ordered
  component interpretations rather than one disconnected graph.
- `components()` and `into_molecules()` expose all component results.
  `molecule()`, `report()`, `into_molecule()`, and `into_parts()` require
  exactly one component and return `SmilesComponentCountError` otherwise; they
  never panic on valid multi-component input.
- Interpretation may install imported aromatic annotations with input
  provenance. It never sanitizes or runs full perception.
- `SmilesParseOptions` bounds input bytes plus parsed atom and bond counts.
  Defaults permit large molecules while preventing unbounded parser allocation;
  callers handling intentionally larger records can opt into higher limits.
- `SmallMolecule::from_smiles` and `from_smiles_sanitized` remain deliberate
  single-component convenience orchestrators, use the default parser limits,
  and report `SmallMoleculeReadError::ComponentCount` for dot-separated input.

## Tests

- Focused unit and regression tests cover the behavior, error, and resource-limit contracts described above.

## Benchmarks

- Unit/fuzz tests cover raw document spans, ordered component interpretation,
  structured single-component rejection, malformed syntax, parse/interpret
  separation, local stereo, and parse-write-parse behavior. Existing external
  benchmark corpora remain the semantic reference.
- Optional external-reference manifests are available for `pubchem-1k`, `pubchem-100k`, `enamine-diversity`.
- Benchmark observations are informational and never determine this feature's release status or repository health.

## Out Of Scope

- SMARTS, reactions, sanitization during parsing, and full grammar parity.

## Revision Notes

- v1-v12: Practical raw SMILES reader and expanded chemistry/stereo coverage.
- v13: Hard break to `SmilesDocument` parse/interpret and remove direct readers
  and `SmilesParseOptions`.
- v14: Keep every ignored non-smoke corpus as explicit local-only validation
  instead of repository-wide required evidence.
- v15: Make the parse/interpret boundary real: parsing validates the complete
  supported grammar into a private syntax program, interpretation consumes it
  without reparsing source text, and returns canonical atom/bond mappings.
- v16: Add configurable input, atom, and bond resource limits while retaining
  `parse_str` as the default-bounded convenience entry point.
- v17: Use PubChem-1k as the required baseline benchmark corpus after retiring the former smoke corpus from public validation.
- v18: Reclassify external-reference parity from a required gate to optional benchmarking without changing implementation behavior or golden expectations.
- v19: Interpret each dot component as a separate connected molecule, expose
  ordered component results, and replace panicking single-molecule access with
  structured component-count errors.
