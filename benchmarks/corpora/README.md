# Benchmark Corpora

Each corpus is self-contained:

```text
<corpus-id>/
  corpus.toml
  sources.lock.json
  data/
  features/
  golden/
  results.toml
```

`data/` is generated locally and ignored. Source locks, feature manifests,
deterministic goldens, and result observations are tracked. Missing manifests
are normal: availability exists only for feature/corpus pairs with a manifest.

PubChem-1k and PDB-100 are the small default selection. PubChem-100k,
Enamine Diversity, PL-REX, and PDB-1000 are deliberate broad selections
included by explicit `--corpus all`.

Checked-in source-lock entry IDs and categories are canonical membership.
Normal builders reconstruct those records. Creating a different corpus version
requires explicit reselection and a new `selection_id`.

PubChem-100 is an exact prefix of PubChem-1k. PDB-10 and PDB-100 are exact
prefixes of PDB-1000. Historical smoke, PubChem-100, and PDB-10 directories
remain useful internal fixtures but are not registered dashboard columns.
