# Benchmark inputs

Each dataset contains `sources.lock.json` and externally supplied files in
`data/`. The lock defines membership, source provenance and expected file hashes.
The tracked smoke corpus is small; larger source files are local and ignored.

The current datasets are `smoke`, `pubchem-100k`, `enamine-diversity`, `pdb-1000`
and `pl-rex`. The runner selects from their locks, independently of reference or
Kekule success. Available matching formats determine the cases for each feature.

Current independent goldens live separately in `../goldens/`. See
[GUIDE.md](../GUIDE.md) for selection, integrity checks, generation, comparison
and the historical PubChem selection limitation. Use `../data.py` to verify,
pack or unpack supplied inputs without changing their bytes or membership.
