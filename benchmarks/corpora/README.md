# Benchmark Corpora

Each corpus is self-contained:

```text
<corpus-id>/
  corpus.toml
  sources.lock.json
  data/
  features/
  golden/
```

`data/` is generated locally and ignored. Source locks, benchmark manifests, and deterministic goldens are tracked. The historical `features/` directory name and `feature_id` manifest field are retained only to avoid rewriting manifests and compressed goldens; they are benchmark-schema vocabulary.

The runner discovers every corpus descriptor and manifest directly from this tree and sorts selected `(benchmark ID, corpus ID)` pairs deterministically. Missing benchmark/corpus combinations are normal unless requested concretely.

Checked-in source-lock entry IDs and categories are canonical membership. Normal builders reconstruct those records. Creating a different corpus version requires explicit reselection and a new `selection_id`.

PubChem-100 is an exact prefix of PubChem-1k. PDB-10 and PDB-100 are exact prefixes of PDB-1000. The smoke corpus provides small tracked fixtures; broader corpora may require locally built ignored data.
