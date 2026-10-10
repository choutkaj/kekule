# Scientific benchmark

Measures where Kekule's results differ from established toolkits, on two
purpose-built datasets, and keeps every known difference in one checked-in list.
A run passes when every difference is known and no known difference has grown.

| Dataset | Contents | Reference toolkits |
| --- | --- | --- |
| `small` | ~4,800 small molecules: ChEMBL clinical, natural-product, large, isotope, uncommon-element, charged, multi-component, stereo-rich and uniform strata; wwPDB CCD ligands, including metal and cluster components; RDKit's atropisomer and V3000 test inputs; 519 RDKit SMARTS table rows | RDKit |
| `bio` | 200 wwPDB entries across 11 strata: ultra-high-resolution, typical and cryo-EM structures, NMR ensembles, nucleic acids and complexes, glycans, modified residues, cofactors, antibodies, neutron structures | gemmi, Biotite, mkdssp |

Records are chosen only from upstream metadata, never from whether either
toolkit can read them. `datasets/<name>/selection.toml` states the rules,
`selection.lock.json` what they resolved to, `manifest.jsonl` every file with
its upstream URL and SHA-256, and `LICENSES.md` the attributions. The curated
tier (`curated.txt`, `curated/`) is committed; the full tier is a hash-locked
release archive.

## Run

A run needs only Rust and Python's standard library:

```sh
python benchmarks/bench.py run small --tier curated
python benchmarks/bench.py fetch bio && python benchmarks/bench.py run bio --tier full
```

`run` builds the release observer, observes every case with Kekule, compares
with the stored reference observations, triages, and writes `report.md`,
`report.json` and `cases.jsonl` under `target/kekule-bench/reports/`. It exits
non-zero on a new difference, a stale known difference or a bound exceeded.
`--task` selects tasks by glob, `--limit` the first records, `--jobs` the
observer threads. `KEKULE_BENCH_CACHE` moves the data and report caches, which
otherwise live in `target/kekule-bench`.

Reference observations are computed once, by maintainers, in the pinned
environment, and stored: `datasets/<name>/references-curated.jsonl.gz` in git,
the full set inside the dataset's release archive. Regenerate them after
changing a dataset or a reference observer it uses: each file records the
dataset version, a digest of the tier's input files and a fingerprint of the
observer source, and a run refuses a file whose record does not match. Generation takes about two minutes for `small`
and four for `bio`; RDKit's resonance enumeration is capped at a fixed number of steps
(`RESONANCE_WORK_LIMIT`) so that a few porphyrins cannot stall it.

```sh
micromamba create -y -p target/bench-reference -f benchmarks/conda-linux-64.lock   # or conda-win-64.lock
micromamba run -p target/bench-reference python benchmarks/bench.py reference small
micromamba run -p target/bench-reference python benchmarks/bench.py pack small
```

## How a comparison works

`observer/` (the unpublished `kekule-bench` crate) reports what Kekule's public
API computes and `kbench/references/` what the reference computes, as the same
facts `[key..., value]`. Small-molecule atoms are keyed by source position and
mmCIF atoms by `_atom_site.id`, so storage choices never count as differences:
the parse compares elements, isotopes, charges, radicals, total hydrogens, maps,
bonds, components and stereo parity, not how hydrogens are stored. Molfile and
SDF writers are checked by Kekule reading its own output back against the
reference's reading of the original, coordinates included; SMILES round trips,
which reorder atoms, are pinned by the invariant suite. Tasks that
build on the parse are *blocked* when the parse already differs in what they
depend on, so one cause is counted once. `kbench/tasks.py` is the only routing
table.

Each difference becomes a signature `path|context|transition`, for example
`atom.hydrogens|N|ar|1→0`. Context comes from the reference's view of the
input, so it does not move when Kekule changes. Numeric tolerances live in one
table in `kbench/compare.py`, each half the printed precision of the compared
value; values are never rounded.

## Triage

`known-differences.toml` lists every accepted difference with a verdict
(`kekule-bug`, `kekule-gap`, `reference-limitation`, `intended-policy`), a
reason, example records and per-tier bounds. When a run reports a new
signature, either fix Kekule, or add an entry that explains it, with the
smallest example records. When an entry's count falls, rerun with
`--update-known` to lower its bound; bounds never rise. A stale entry, whose
examples no longer show the difference, is deleted. Never widen a tolerance or
drop a compared field to make a run pass.

## Datasets

`bench.py fetch` downloads and verifies the full-tier archives. Maintainers
change a dataset by editing `selection.toml`, bumping its version and running
`bench.py select <name>` (which contacts the upstream databases),
`bench.py reference <name>` and `bench.py pack <name>`, then uploading the
archive to the release named in `selection.lock.json`. ChEMBL records are CC BY-SA: never copy them into
`crates/*/tests/fixtures`.
