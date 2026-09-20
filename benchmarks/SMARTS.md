# SMARTS behavioral validation

This is an optional scientific comparison, not a routine CI or release gate.
RDKit/OpenFF are reference sources, never Rust runtime dependencies.

Build the observer with:

```text
cargo build -p kekule-bench --bin smarts_conformance --locked
```

Run `smarts_conformance.py` in the existing RDKit 2026.03.3 environment, passing
`--binary`, `--queries` (one or more complete `.smarts` or `.offxml` files),
`--molecules` (external whitespace-delimited SMILES files), and `--output`.
The script checks the RDKit version. Every input row is retained, with source
checksums, query origin, target origin, both complete observations, and an explicit
comparison outcome. It returns nonzero on differing successful results or
unmatched implementation/reference errors.

The Rust observer exchanges one JSON object per line. A request has `smarts`,
`smiles`, `explicit_hydrogens`, and `mdl`. Omit `smiles` for syntax-only inspection.
Success includes query sizes, every full mapping, every ordered tagged mapping,
tag labels, and every target atom's aromaticity. The mapping lists are sorted;
atom order inside each mapping is never sorted. Query automorphisms and repeated
tagged tuples are retained. RDKit uses `useChirality=True`, `uniquify=False`, and
a checked one-million-result bound. Rust exhaustion is an error, not partial
success. The observer is intentionally confined to the unpublished benchmark crate.

`both_rejected` means the source query is rejected by both parsers; it is not a
successful matching case. `dialect_exclusion` requires a successful reference
observation and a recognized unsupported construct in kekule. The fixture and
comparison-contract tests run without RDKit or network access:

```text
python -m unittest discover -s benchmarks -p test_smarts_conformance.py
```

The conformance data sets are:

- All 518 rows of the existing pinned RDKit functional-group/reactivity corpus,
  crossed with all eight existing PubChem SMILES smoke targets.
- All 521 rows in `smarts-fixtures/openff-smarts`, crossed with those eight targets and
  its four externally sourced PubChem targets. These use explicit hydrogens and
  MDL aromaticity.

See the corpus provenance files and `SMARTS-VALIDATION.md` for measured results.
Focused Rust regressions separately cover syntax and scientific edge cases that
these external targets do not exercise, including recursive negation, numeric
rings, graph permutations, Boolean stereo, repeated topology definitions, and
match-budget exhaustion. Deliberate Boolean-stereo differences from RDKit are
documented in the public query dialect; corpus agreement is not evidence of
universal RDKit equivalence.
