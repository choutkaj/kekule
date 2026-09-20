# RDKit query tables

This corpus contains every non-comment query row from three tables distributed
with RDKit 2026.03.3: `FunctionalGroups.txt` (38),
`Functional_Group_Hierarchy.txt` (52), and `SmartsLib/RLewis_smarts.txt` (428).
They describe functional groups and reactivity filters, not molecular fixtures.

The supplied files came from the installed conda-forge `librdkit` package.
Every original file matched its package manifest's SHA-256. The archive URL,
archive checksum, installed paths and individual checksums are recorded in
`sources.lock.json`. Original bytes and the BSD 3-Clause license are in
`upstream/`. Source headers and author attributions are preserved.

Extraction excludes only blank/comment lines. The two functional-group tables
provide label then SMARTS as their first two nonempty tab-delimited fields.
The Lewis table provides SMARTS then the remaining text as a label. Surrounding
field whitespace is removed; query text is otherwise unchanged. Packed records
contain SMARTS, a tab, and the source label. Source IDs retain table names and
one-based original line numbers. Duplicate queries at different source rows
remain distinct observations. Tests verify complete, unfiltered extraction.

All 518 queries are registered in the standard `query.smarts` benchmark. Each is
parsed and matched against every target in `benchmarks/query-smarts.json`.
Unsupported grammar and resource failures remain errors; no query is removed
based on either engine's result. Original extraction is identical to the earlier
`smarts-fixtures/rdkit-queries` copy, which remains available for historical runs.
See [SMARTS benchmark validation](../../QUERY-SMARTS-VALIDATION.md).
