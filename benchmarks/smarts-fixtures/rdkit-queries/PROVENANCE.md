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

All 518 queries parse in the pinned independent reference. Native unsupported
grammar remains an error, never an exclusion or an inferred empty query. This
corpus measures parsing and graph size; predicate semantics require behavioral
matching tests. The initial native result is 384 agreements and 134 explicit
unsupported-grammar errors, with no count disagreements.

The broad SMARTS implementation additionally uses every row in full-mapping
comparisons; see [SMARTS validation](../../SMARTS.md) and its measured report.
