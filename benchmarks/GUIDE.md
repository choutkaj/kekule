# Correctness comparisons

Reference values are calculated once in a separate generation step and stored.
Normal benchmarks calculate Kekule values from the same supplied inputs and
compare against those stored goldens. They never regenerate or modify them.
A successful run
requires every selected case to agree. Reference errors, Kekule errors, missing
input formats and missing reference adapters all fail the run. Matching failures
are never scientific agreement.

## Run

Normal parser and algorithm benchmarks need only Rust and stored goldens:

```text
cargo benchmark --list
cargo benchmark --feature io.smiles.parse --dataset smoke
cargo benchmark --feature io.sdf.parse --dataset enamine-diversity
cargo benchmark --feature stereo.cip --dataset pl-rex
```

Progress bars show the current dataset/feature, completed cases, total and phase.
Agreement, disagreement and error counts follow each completed feature. Progress
includes multiple source formats and cases that fail.

Goldens live in `benchmarks/goldens/<dataset>/<feature>.jsonl.gz`. Complete files
are prepared for all 25 features across all five datasets, without case limits.
See [GOLDENS.md](GOLDENS.md) for coverage, reference-error counts and reference
versions. Missing-format cases and reference failures remain stored errors;
complete files do not mean every comparison passes.

To prepare a new dataset/feature, activate the reference environment in
`reference/environment.yml` and explicitly generate:

```text
cargo benchmark generate --feature FEATURE --dataset DATASET
```

The generation command runs the reference software only; it does not run Kekule.
Use `--python PATH` on this command to select a different reference environment.
RDKit supplies molecular goldens; Biopython and mkdssp supply biomolecular ones.
The activated environment must expose mkdssp, its libraries and CIF dictionary.

Existing goldens are never overwritten. To regenerate deliberately, write to a
new directory with `--goldens DIR`, review the new results, then use that same
directory for comparison. Incomplete generation is not published. Reference
errors are stored and remain failures; successful cases do not disappear because
another case failed.

The default selects every source ID. `--limit N` explicitly selects a smaller,
deterministic subset, including during generation. A partial golden set causes
errors for any selected cases it does not contain. Missing, malformed, duplicated
or stale goldens never trigger automatic regeneration. An absent or unreadable
golden file stops the run with its cause before that feature evaluates Kekule;
it is not reported as thousands of individual comparison errors. `--jobs N`
controls Kekule concurrency; `--output FILE` selects a new report.

Writer benchmarks retain independent validation: RDKit reads newly emitted
Kekule text and compares its meaning against **stored** expectations. It never
recalculates those expectations. Use an activated reference environment, set
`KEKULE_WRITER_PYTHON`, or select it explicitly:

```text
cargo benchmark --feature io.smiles.isomeric --dataset smoke --writer-python PATH
```

This is the only normal benchmark path that uses Python. The former required
`--reference-python` option is gone.

`io.mol.parse` and `io.sdf.parse` read supplied V2000 and V3000 records directly.
There is no parser benchmark that generates its own input with Kekule's writer.
The versioned writer features specify the requested output version. A format
that cannot encode the input must report an error, not quietly discard fields.

## Outputs

A normal run writes a summary `FILE.json` and per-case comparisons
`FILE.cases.jsonl`. Writer cases retain the original emitted text. Stored golden
files are read only. Generation writes compressed reference JSONL files plus its
own summary; it does not produce Kekule observations.

Each evaluation is either `{"status":"ok","value":...}` or
`{"status":"error","message":"..."}`. The transport rejects unknown fields,
missing outputs and nonfinite numbers. Stored entries identify the source ID,
path, record index, input hash, reference tool/version and measured values.
Loading verifies dataset/feature identity and rejects duplicate cases or mixed
reference versions; comparison verifies the current input hash.

There are no implementation snapshot sources, acceptance commands, evidence
migrations, hardware fingerprints or per-case comparison overrides. The old
`corpora/*/references.json` indexes and compressed `golden/` files are historical
artifacts and are never read by this runner. Their earlier scores do not validate
the new comparisons. Former standalone SMILES/CIP/MOL comparison commands have
been retired; use the features above.

## Comparison contract

Keys, types, values, cardinality and ordering are checked. Numbers compare
exactly: no rounding, tolerances or averaging can turn a difference into a pass.
Small floating-point differences and large discrepancies both remain visible;
a disagreement alone does not establish the cause or its scientific magnitude.

Only representation order is normalized: undirected bond endpoints and bond
list order, ring membership-set order, and arbitrary DSSP sheet/strand/ladder
labels and residue order. Dative direction is retained. Stereo parity is expressed
relative to explicit atom indices and a fixed carrier order. Distinct global
graphs, components and stereo placements cannot collapse to local atom hashes.

Molecular reader checks compare complete indexed component graphs: atoms, bonds,
charges, isotopes, radicals, atom maps, declared and inferred hydrogens, represented
valence, aromaticity, stereo and relation groups. Both adapters run ordinary
perception/sanitization before observing this state. Coordinates are in angstroms;
SDF data fields remain an ordered list. Kekule preserves duplicate names;
RDKit cannot represent them, so such cases report a reference error. Component
and atom order follow the input correspondence. This is deliberately stronger
than a count or formula comparison.

SMILES writer outputs are read by RDKit and compared through full canonical
isomeric CXSMILES identity. MOL/SDF writer outputs are also read by RDKit, with
coordinates and properties retained. Writers do not validate themselves by
re-reading with Kekule. Canonical SMILES additionally checks a read/write fixed
point and invariance after reversing atom and bond numbering while retaining all
stereo and relation groups. These deterministic probes do not exhaust all
permutations. Molecule output components are sorted in canonical mode.

Valence reads Kekule's actual `represented_valence` and `implicit_hydrogens`
results. No benchmark nitrogen/aromaticity rules reconstruct or adjust them.
Dative valence therefore exposes any difference between the libraries' models.
Hydrogen transformations compare the complete added and collapsed graphs, using
RDKit's default removal policy without matching it to Kekule's policy.
Reference masses are the direct RDKit descriptor values, without charge-based
corrections. Rotatable-bond references retain explicit hydrogen vertices.
Substructure comparisons retain every query-to-target mapping, without the old
1000-match cap or atom-set projection. Search resource errors remain failures.
The fixed query list is in the two adapters. `query.smarts` separately compares
syntax acceptance and graph size; it is not a complete proof of predicate
semantics. Behavioral predicate checks come from the substructure feature.

mmCIF compares every decoded tag/value, including non-atom categories, entity
IDs, charges and distinct `.`/`?` values. Biopython's single-block limitation is
reported as a reference error. DSSP consumes the original supplied mmCIF and
uses the first model through the reference APIs. No alternate-location snapshot
or archive-category deletion edits the reference input. Differences in the
libraries' model/alternate-location policies remain visible.

The stereo adapter's index/parity conversion follows public Kekule definitions
and RDKit conventions, with geometry-based regression checks. In particular,
RDKit's atropisomer convention is defined in
[Atropisomers.cpp](https://github.com/rdkit/rdkit/blob/master/Code/GraphMol/Atropisomers.cpp).

## Inputs and coverage

`corpora/<dataset>/sources.lock.json` supplies membership and input hashes.
Selection never consults reference success or implementation support. Every
available matching-format source is evaluated, so an ID with both SDF and
SMILES can contribute two cases to an algorithm feature. Multi-record source
files are split before either engine runs. One failed record cannot remove the
remaining records. Missing formats are recorded as errors for the selected ID;
`all` includes every feature, including combinations lacking source material.
This means heterogeneous datasets such as `smoke` can contain coverage errors.

Enamine includes all 50,240 supplied SDF records. Packs 049–051 restore the 2,881
V3000 records that the old builder filtered out, without changing their bytes.
All 1,000 PDB entries reach DSSP, including the five previously absent from its
reference index. Nothing filters stereo, disconnected records, wildcards or
reference failures.

The existing PubChem membership remains a historical, preselected 100,000-ID
sample: its old builder imposed V2000, size and RDKit-success requirements.
These source files cannot establish performance on the compounds that builder
did not retain. The new runner does not repeat those filters, but it cannot
recover missing PubChem inputs. Do not describe this sample as unbiased.

Input files outside `smoke` remain local/ignored. Transport them with `data.py`:

```text
python benchmarks/data.py verify --dataset enamine-diversity
python benchmarks/data.py pack --dataset enamine-diversity --archive target/enamine.tar.gz
python benchmarks/data.py unpack --dataset enamine-diversity --archive ARCHIVE --sha256 HASH
```

## Measurement and maintenance

Timing is supplementary. `kekule_ms` measures parallel batches of public API
work and result construction, including writer contract checks. `reference_ms`
measures reference generation in generation reports, and only independent writer
validation in normal reports (zero for other features). Golden loading, startup,
JSON transport and comparison are excluded. These are workflow timings; no
cross-library speed ratio is claimed.

```text
cargo test -p kekule-bench --locked
python -m unittest discover -s benchmarks/reference -p test_runner.py
python -m unittest discover -s benchmarks/reference/rdkit -p "test_*.py"
python -m unittest discover -s benchmarks/reference/biopython -p "test_*.py"
```

Scientific runs remain optional, outside routine CI/release gates. Fix the
implementation or a demonstrably incorrect adapter when values differ; never
substitute Kekule output for a golden. New features need independent reference
code and regressions that prove meaningful mutations fail comparison.

Trajectory development checks remain a separate workflow documented in
`reference/trajectory/VALIDATION.md`; this executable does not claim coverage
of their trajectory operations or historical results.
