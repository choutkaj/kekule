# Correctness and timing benchmarks

The benchmark executable is the unpublished `kekule-bench` workspace package.
It calls Kekule's public API. External scientific tools are reference tools only.
Benchmark execution remains optional, outside ordinary CI and release gates.

This guide describes the replacement interface. The earlier benchmark READMEs
describe the removed `xtask` interface and were left unchanged under the repository's
README editing policy.

## Run

```text
cargo benchmark --list
cargo benchmark --feature io.smiles.parse --dataset pubchem-100k
cargo benchmark --feature io.sdf.v2000.parse --dataset enamine-diversity
cargo benchmark --feature stereo.cip --dataset pl-rex
cargo benchmark --feature bio.secondary-structure.dssp --dataset pdb-1000
cargo benchmark --feature all --dataset smoke
```

`cargo benchmark` builds in release mode. `--feature all` and `--dataset all`
select the available pinned evidence. A concrete unsupported combination errors.
Each input record is evaluated once per feature. The same output is used for correctness;
there is no warmup or repeated timing pass.
Evaluation uses all available CPUs by default. `--jobs N` limits concurrency;
`--jobs 1` runs serially. Small batches span input files, including individual PDB
structures, and results retain source order regardless of completion order.
Each feature shows a progress bar labeled with its dataset and feature, completed
inputs, total selected inputs, and percentage. Errors and unsupported inputs count
toward completion. Bars update between evaluation batches, outside measured Time,
and remain above each feature's summary. Bars are hidden when stderr is redirected
or the terminal declares `TERM=dumb`.
`--output FILE` chooses a new JSON result file. The default is a unique file in
`target/benchmarks`. Existing reports are never overwritten.

The old `cargo xtask`, corpus builders, feature TOML manifests and golden-acceptance
command have been removed. There are no compatibility aliases.

## Inputs and selections

| Dataset | Purpose | N |
| --- | --- | ---: |
| `pubchem-100k` | Broad cheminformatics, supplied SDF and SMILES | 100,000 |
| `enamine-diversity` | Compound-library workload, supplied SDF and SMILES | 50,240 |
| `pl-rex` | Primary ligands, including supplied 3D coordinates | 164 |
| `pdb-1000` | Proteins, nucleic acids, complexes, heterogens, multiple models | 1,000 |
| `smoke` | Small checked-in external edge fixtures | 20 |

N is the number of source records in the pinned dataset. Every feature selects
the whole dataset by default, including expensive features. There are no automatic
dataset or feature limits. An explicit `--limit N` selects a smaller run;
selecting more than a dataset contains selects the whole dataset.

Molecular IDs are ordered by SHA-256 of `kekule-benchmark-v1:<id>`, with the ID
as a tie-breaker. Limits are nested selections of the same source membership,
independent of feature and input format. PDB uses its existing locked order,
preserving the 10/100/1000 subsets. Limits count source records, not connected
components, coordinate models, file packs, or successfully processed molecules.

Features use the original pinned input format recorded in their evidence. They
do not turn every dataset into RDKit-generated SMILES. A selected ID with no
input/reference for that feature is reported as unsupported. This also makes
the intentionally narrow applicability of individual smoke fixtures visible.

Feature adapters do not exclude SMILES containing stereochemistry or wildcard
atoms. They attempt every supplied record and report actual outcomes. CIP also
evaluates molecules without stereo: an empty descriptor list is a valid result,
while parsing, perception, and assignment failures remain visible.

PubChem membership is unchanged: the existing selection uses the CID 1–500,000
shard, V2000 records, bounded molecule sizes, and RDKit-based selection checks.
It is not an unbiased sample of all chemical inputs. Keep rare external fixtures
visible separately. PubChem's 2D depictions are not a SASA input; future geometry
features should use PL-REX's supplied 3D coordinates or PDB structures.

Subset datasets have been removed. Use `pubchem-100k --limit N` or
`pdb-1000 --limit N` for smaller runs.

## Data storage

`corpora/<dataset>/sources.lock.json` pins source IDs, file membership and hashes.
`references.json` indexes existing evidence: tool/version, input paths and
preparation notes. It does not configure algorithms, tolerances or plugin loading.
Changing packaging preserved the compressed `golden/` outputs. Removing the old
SMILES and CIP filters subsequently extended the RDKit 2026.03.3 evidence for
previously skipped records. Existing asserted records and source hashes were
preserved; only `unsupported` placeholders and missing CIP records were filled
with independent RDKit results.

Inputs other than the checked-in smoke fixtures remain ignored local files.
The runner verifies each used input against both its source lock and its reference.
Missing files, changed bytes, bad reference versions and malformed outputs fail
explicitly and appear in the JSON report.

Use a standard Python 3.11+ interpreter to verify or transport pinned inputs:

```text
python benchmarks/data.py verify --dataset pubchem-100k
python benchmarks/data.py pack --dataset pubchem-100k --archive target/pubchem-100k.tar.gz
python benchmarks/data.py unpack --dataset pubchem-100k --archive DOWNLOADED.tar.gz --sha256 EXPECTED_SHA256
```

`pack` prints the archive hash. Store that archive on ordinary artifact storage
and distribute its hash with it. No hosted archive is assumed or automatically
published. `unpack` stages and verifies all members before installation, rejects
unexpected paths/links/duplicates, and never replaces different existing input
bytes. `verify` checks all pinned source files, including raw provenance files
that a benchmark may not need when a complete pack is available.

## Correctness

The existing asserted JSON fields, ordering normalizations, preparation choices
and numeric tolerances are retained in the feature code and comparator. Each
selected record is evaluated independently, restoring its source record index
before comparison. One bad record cannot abort all other records in a pack.

Reports separate agreement, disagreement, unsupported input, record errors, and
execution/provenance errors. Matching parse failures do not count as successful
scientific agreement. A missing CIP record is not a successful comparison; an
explicit empty descriptor list must agree with the reference. A mismatch retains source ID, fixture path, record
index and the first differing field for that record. Every failing record is kept.

Numeric outputs additionally report mean and maximum absolute errors by field
path, without mixing different physical quantities. Agreement uses the existing
field-specific tolerances. Floating-point metrics describe the normalized values
actually compared. Correlation is not used as proof of numerical agreement.

Historical `*-manual-semantic` outputs are implementation snapshots. They remain
available and explicitly labeled, but are not independent scientific validation.
There is no command to accept implementation output as a reference.

The biomolecular parsing feature compares all retained atom-site fields.
DSSP uses the first model and the existing highest-occupancy reference preparation;
its exact policy and compared fields remain in `src/features/bio.rs` and
`reference/biopython/run_feature.py`. Source identity orders residue comparisons.

## Timing and live references

`Time` measures elapsed time for parallel evaluation batches: parsing, the selected
feature, and materializing comparison outputs. Each batch holds at most 256 records.
The timer stops before outputs are compared with the reference and dropped.
Inputs are preloaded. Thread-pool creation, file reading,
reference decompression, checksum validation, correctness comparisons, JSON
encoding and report writing occur outside measurement. These are adapter workflow
timings, not isolated kernel microbenchmarks. Scheduling and failed evaluations
are included. Each fixture lists all evaluated record IDs, including failures.

Results use schema version 3 and contain `time_ms` per feature/dataset,
obtained by summing elapsed batch times, not individual worker durations. Batches
can span fixtures, so times are not attributed to individual files. Worker counts
are not recorded. There are no timing
repetitions or `--samples` option. These single-pass times can vary between runs.
Results also record source IDs, dataset hashes, reference versions,
Git revision/dirty status, Rust version,
build mode, OS, architecture and available machine information. Compare revisions
on the same machine, with the same selection and measurement settings.

The default uses pinned reference outputs and reports Kekule timing. To also run
and time an installed independent reference, supply its Python interpreter:

```text
cargo benchmark --feature io.smiles.parse --dataset pubchem-100k --limit 1000 --reference-python PATH_TO_REFERENCE_PYTHON
```

Use the supplied RDKit or Biopython/DSSP environment definitions. Stereo SMILES
schema 2 requires RDKit 2026.03.6; older feature evidence includes RDKit 2026.03.3.
Choose a suitable environment for the requested feature. Live results record the
actual version, evaluate each input once, compare against the already computed
Kekule output, and never replace pinned expectations.
Differences against either the pinned or live reference return a nonzero exit.

Optional live reference adapters currently evaluate serially, once per input,
after the corresponding Kekule batch finishes. Their results and `time_ms` values
are stored in `live_references`, with the input IDs for each batch.
Reference timing excludes interpreter startup, imports and JSON transport. The
Biopython/DSSP adapter includes staged-file reading and DSSP subprocess execution.
Some historical feature adapters perform different preparation work, so the
report deliberately does not compute a cross-library speed ratio. Use the labeled
timings for each workflow; only compare speed ratios after matching their scopes.

Generate a standalone candidate for independent review with:

```text
python benchmarks/reference/run.py --feature io.smiles.parse --input INPUT.smi --output target/candidate.json
```

The output path must be new. Existing reference-specific diagnostic probes remain
available under `reference/`, with their Rust examples in this benchmark package.
Trajectory reference checks remain documented in `reference/trajectory/VALIDATION.md`;
static PDB models do not replace real trajectory inputs.

## Maintenance

Add an ordinary function in `src/features/`, register the feature, and provide
independent reference evidence and focused comparison regressions. Input
preparation and comparison policies belong beside their feature. The shared
runner owns only selection, execution, timing and reporting. Avoid new manifest
languages, generic chemistry adapters and framework layers.

```text
cargo test -p kekule-bench --locked
python -m unittest discover -s benchmarks/reference -p test_runner.py
python -m unittest discover -s benchmarks/reference/rdkit -p "test_*.py"
```

The last command requires RDKit. The runner and comparison regressions use
ordinary Rust tests and do not execute the broad external benchmarks in CI.
