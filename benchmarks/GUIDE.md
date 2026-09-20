# Correctness benchmarks

The benchmark measures agreement with independent references on externally
supplied, hash-locked inputs. Generation runs only the reference software;
comparison runs Kekule against stored observations. Matching failures never
establish correctness. Scientific runs are optional, outside CI/release gates.

## Run

```text
cargo benchmark --list
cargo benchmark --feature io.smiles.parse --dataset smoke
cargo benchmark --feature io.sdf.parse --dataset enamine-diversity
cargo benchmark --feature all --dataset smoke --writer-python PATH_TO_RDKIT_PYTHON
cargo benchmark generate --feature FEATURE --dataset DATASET --python PATH --goldens NEW_DIRECTORY
```

A fresh checkout includes the smoke inputs and goldens. Full dataset inputs and
golden archives stay local and are excluded from Git and Cargo packages; no
download or upload is automatic. Their source locks and golden manifests remain
tracked. To run a full dataset on another machine, supply the original inputs
and matching golden archives locally at the paths below. Comparison verifies
their hashes against the tracked provenance. Missing files stop the run.

The default selects all locked source IDs. `--limit N` selects a deterministic
subset before either engine runs; molecular subsets are nested and PDB subsets
retain lock-file order. `--jobs N` controls Kekule concurrency. Reports and per-case
JSONL are saved under `benchmarks/runs/`, with unique run names. `--output FILE`
selects another report location; its summary is also archived in `benchmarks/runs/`
without duplicating its case records. Existing reports are never overwritten.

Normal parser and algorithm checks need Rust and the stored goldens. Writers
also need RDKit to read the emitted text. Select it with `--writer-python`,
`KEKULE_WRITER_PYTHON`, or the activated environment's `python`.
Create the pinned reference environment from `reference/environment.yml`;
engine-specific environment files are available for narrower installations.
DSSP requires mkdssp, its shared libraries, and its CIF dictionary on the
activated environment's paths.

Generation is explicit and never evaluates Kekule. Write regenerated goldens to
a new directory, investigate changes, and review them before adopting them.
Missing, stale, malformed or duplicated goldens stop comparison; they do not
trigger regeneration. A partial generation does not cover omitted cases.

## Dashboard

Every comparison refreshes `benchmarks/runs/index.html`, including comparisons that
fail or stop with an error. Open that page once: it checks the adjacent local data
file every five seconds and selects the newest run when history changes. Unchanged
history preserves your selected run and search. Reports are ordered by their recorded
start time, not file modification time. Reference generation is excluded.

Python 3.11+ and its standard library are needed to render the dashboard. The command
uses `KEKULE_DASHBOARD_PYTHON`, an executable path in the ignored file
`benchmarks/.dashboard-python`, the selected writer Python, or Python on `PATH`.
If rendering fails, a warning identifies the problem; reports and the scientific
exit status are preserved. `KEKULE_BENCHMARK_RUNS_DIR` can override the local history
directory. There is no server or upload.

To refresh history manually, including checkpoints left by an aborted process:

```text
python benchmarks/dashboard.py
```

Invalid or incompatible reports are skipped with a warning; identical reports are
deduplicated. Incomplete comparison reports remain labeled in provenance. Original
reports are retained. Old runner filenames carry a sortable start time; other
untimestamped imports sort after timestamped runs, without inferring dates from mtime.

You can also build a fixed snapshot from explicit schema-2 comparison reports:

```text
python benchmarks/dashboard.py target/smoke.json target/full-corpus-sample.json
```

Open `target/benchmark-dashboard/index.html` in a browser, or choose another
destination with `--output PATH.html`. Python 3.11+ and its standard library are
sufficient. This explicit-report mode produces a single file without live polling,
a server, network requests or external assets. It can later be hosted unchanged.
Rerun the command to refresh this derived page; input reports are never modified.

The dashboard shows a searchable feature/dataset comparison table and an overview
of the data corpora: source ID counts (N), supplied formats
and selection notes. Plots group exact and within-precision agreement into one
"Agrees" category; report counts and comparison tolerances are unchanged. Runs
remain separate. Cells distinguish full selections, deterministic samples, stale
provenance and absent results; run completion is recorded in the provenance panel.
No overall correctness score or cross-library speed ratio is computed. Counts
repeat inputs across features and formats; unavailable formats stay outside the
measured denominator.

Only aggregate counts, source membership sizes and hashes/tool versions are
embedded. Inputs, golden payloads, per-case observations, local paths and raw
error messages stay local. The page can be generated without the full datasets
or golden archives. Run history and local Python configuration are excluded from
Git and Cargo packages; fixed snapshots under `target/` are likewise local.
Publication is a separate, deliberate step; this command uploads nothing.
"Download plotted data" saves a sanitized dashboard export, not the original
runner report. Rebuild the page from the original comparison summaries.

```text
python -m unittest discover -s benchmarks -p test_dashboard.py
node benchmarks/test_dashboard.cjs
```

## Observations and comparisons

The versioned contract is [contract.json](contract.json). Successful observations
must satisfy the feature's schema: required fields, types and collection presence
are checked, and unknown fields are rejected. Nonfinite measurements fail before
serialization. Nullable measurements distinguish an unavailable value from a
missing required field. Coordinates occur only when a conformer exists.

Every structural difference and numerical difference is retained, with its path
and both raw values. Discrete quantities compare exactly. For floating-point
values, the allowance is `16 * f64::EPSILON * max(abs(expected), abs(actual))`,
covering the fixed decimal parsing/unit conversion chain. Additional allowances
follow output precision, rather than an empirical fit to results:

| Observation | Additional absolute allowance |
| --- | ---: |
| DSSP legacy energies, phi/psi, reported alpha/kappa | 0.05 kcal/mol or degrees |
| DSSP TCO | 0.0005 |
| V2000 writer coordinates | 0.00005 angstrom |
| All other numerical observations | None |

Angles use circular distance. Values are never rounded or replaced. Each case
reports both exact equality and agreement within the declared precision; mass
constant differences and other larger numerical differences remain failures.
Biopython's independently calculated omega uses its parsed coordinate precision;
its deviations are retained, without a corpus-fitted tolerance.

Only representation order is normalized: undirected bond endpoints/list order,
ring membership sets, and arbitrary DSSP sheet/strand/ladder labels and residue
order. Dative direction, atom correspondence, stereo parity and group membership
remain significant. Complete indexed molecular graphs include atoms, bonds,
charges, isotopes, radicals, maps, declared and inferred hydrogens, represented
valence, aromaticity, stereo and relation groups. Both engines apply their
ordinary perception/sanitization before observing this state. SDF properties
remain ordered; duplicate names that RDKit cannot preserve are reference errors.

SMILES writers compare complete canonical isomeric CXSMILES identity, without
requiring the same traversal string as RDKit. Canonical mode additionally checks
a read/write fixed point and invariance under reversing atom and bond numbering.
These probes do not exhaust all permutations. CXSMILES source extensions are
parsed by RDKit; unsupported extensions cause an explicit Kekule error.

Versioned MOL/SDF writers must emit the requested format. RDKit reads the emitted
bytes independently. The MOL model API has no title/property input, so its
expected serialized title is empty and its property list is empty; the original
source observation remains intact in the report. SDF document writers retain
titles and properties. Coordinates and all chemical fields remain asserted.

Valence and descriptors read each library's direct results. No adapter adjusts
aromatic nitrogen, masses, formal charges, or hydrogen-removal policy to force
agreement. Substructure retains every query-to-target mapping without a match
cap. The shared 18-query input is [queries.smarts](queries.smarts).
`query.smarts` retains syntax and graph-size checks and compares complete ordered
query-to-target mappings on the 16 externally supplied molecules pinned in
[query-smarts.json](query-smarts.json). Stereo is enabled; automorphisms are
retained. Tags and their ordered projections are asserted without applying
force-field-specific tag rules. Each source query is one case, and every target
must agree. The `rdkit-queries` corpus contains all 518 query rows from three
RDKit tables. Molecular SMILES inputs remain query cases as before. See
[SMARTS benchmark validation](QUERY-SMARTS-VALIDATION.md) for coverage and limits.

```text
cargo benchmark --feature query.smarts --dataset rdkit-queries
```

`algo.aromaticity.rdkit-like` and `algo.aromaticity.mdl` are independent features.
Both compare all indexed atom flags and all endpoint-qualified bond flags,
exactly. The MDL adapter explicitly selects `AromaticityModel::Mdl`; the RDKit
reference kekulizes with `clearAromaticFlags=True` before calling
`SetAromaticity(..., AROMATICITY_MDL)`. Default perception/sanitization and the
explicit model application are included in the timed workflow. Molecules keep
their supplied hydrogen representation. Invalid inputs and perception failures
remain reported errors. See [MDL results](MDL-VALIDATION.md) for measured coverage.

```text
cargo benchmark --feature algo.aromaticity.mdl --dataset smoke
cargo benchmark --feature algo.aromaticity.mdl --dataset all
```

mmCIF compares all decoded tags/values, including non-atom categories and distinct
`.`/`?` tokens. It does not claim full biomolecular topology interpretation
coverage. Biopython's multi-block limitation remains a reference error. DSSP uses
the first model and original source bytes, including archive metadata. It
compares secondary structure, backbone measurements, hydrogen-bond energies and
partner identities, helix/sheet/ladder information, explicit beta partners, and
CA-C-N-CA omega. Alternate-location/model-policy differences remain visible.

## Integrity and reports

Goldens live in `goldens/<dataset>/<feature>.jsonl.gz`, with adjacent
`.jsonl.meta.json` manifests. Manifests bind the compressed bytes, input lock,
contract/schema, reference version and case count. Newly generated files also
record the adapter source digest. Imported older values explicitly mark the
unrecorded original generator fingerprint; it has not been fabricated. The
manifest is published only after complete generation, as the commit marker.
Repository text fingerprints normalize CRLF to LF, so checkout line endings do
not invalidate a run. External input bytes and compressed goldens use exact hashes.

SMARTS version 2 additionally binds its target panel and limits into a
feature-specific contract hash: SHA256 of the base contract's hex digest, a
newline, and the LF-normalized `query-smarts.json`. Reports carry this override
in `implementation.feature_contracts`. Count-only references fail validation;
historical reports remain visible as stale. Other features retain the base
contract. Previous SMARTS references are archived in
`goldens/legacy/query-smarts-v1/`; the standard catalogue ignores that archive.

Goldens are streamed with one lookahead record, ordered by fixture, record index
and source ID. Successful records require reference and input identities; selected
cases verify their exact input bytes. The whole file is checked, even with a
case limit, so a small selection can still incur substantial validation I/O.
Memory does not grow with corpus size. Writer readers must match the stored
reference tool/version.

Reports record their start time in `started_at_unix_ms`, the Git revision/worktree
state, executable hash, current adapter
and contract hashes, and each feature's golden manifest. A summary is written
before evaluation and checkpointed after each feature. Incomplete runs carry
`complete: false` and the error; per-case JSONL preserves completed observations
and emitted writer text. Reference processes have a 300-second deadline; DSSP
invocations have a 120-second deadline. Failed reference batches are isolated
case by case; Kekule panics become errors. These are resource failures, not
agreements. A process abort/OOM can still terminate a run; its incomplete report
must not be interpreted as a complete score.

`cases`, `agrees`, `disagrees` and `errors` cover applicable selected inputs.
`not_applicable` accounts separately for selected IDs with no supplied format.
Input, reference, Kekule, writer-validation and schema errors are reported
separately. Error categories can overlap when both engines fail. An independent
writer-read failure is not automatically attributed to Kekule. No absent format is
counted as agreement. A run passes only if it measured at least one applicable
case and all applicable cases agree. Always inspect coverage alongside agreement;
there is no meaningful universal correctness percentage across these features.

`kekule_ms` measures parallel public-API work plus observation construction and
writer contract checks. `reference_ms` measures generation, or writer reading in
comparison mode. Loading, process startup, transport and comparison are excluded.
These are workflow timings, not evidence of a cross-library speed ratio.

## Inputs and maintenance

`corpora/<dataset>/sources.lock.json` defines membership and supplied file hashes.
All matching formats are tested: one ID can contribute both an SDF and a SMILES
case. Source/feature applicability is independent of either engine's success.
All 50,240 Enamine SDF records and all 1,000 PDB entries remain selected. The
historical PubChem sample was preselected for V2000, size and RDKit success; it
cannot establish performance on compounds that were excluded. Do not describe
that sample as unbiased. Current counts are in [GOLDENS.md](GOLDENS.md).

Large inputs belong in `corpora/<dataset>/data/`; full golden archives belong in
`goldens/<dataset>/`. Keep local bundles, alternate generations and reports under
the ignored repository `target/` directory, or outside the checkout. Smoke and
the small RDKit query corpus, their reference payloads, and provenance belong in
commits. `data.py verify`, `data.py pack`, and
`data.py unpack` operate on supplied inputs only (see `--help`).
Do not synthesize benchmark molecules. Toy inputs belong in focused regressions.

```text
cargo test -p kekule-bench --locked
python -m unittest discover -s benchmarks/reference -p test_runner.py
python -m unittest discover -s benchmarks/reference/rdkit -p "test_*.py"
python -m unittest discover -s benchmarks/reference/biopython -p "test_*.py"
```

Fix a demonstrated adapter defect or investigate Kekule when values differ.
Never replace a reference with Kekule output or adjust a threshold to conceal a
mismatch. The separate optional trajectory workflow is documented in
[trajectory validation](reference/trajectory/VALIDATION.md); it is not included
in this executable's feature coverage.
