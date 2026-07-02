# Reference Generators

Reference generators produce normalized JSON golden data from external tools using externally supplied
fixtures declared in each feature's validation manifest.

These scripts are validation infrastructure only:

- RDKit is used for small-molecule parser, ring, and aromaticity reference output.
- Biopython is used for mmCIF atom-site and SMCRA hierarchy reference output.
- Neither tool is a Rust runtime dependency.

Create the matching micromamba environments from the repository root:

```bash
micromamba create -f validation/reference/rdkit/environment.yml
micromamba create -f validation/reference/biopython/environment.yml
```

Run dependency checks through those environments:

```bash
micromamba run -n molecules-rdkit-reference python validation/reference/rdkit/run_feature.py --feature io.sdf.v2000.parse --corpus smoke --check-deps
micromamba run -n molecules-biopython-reference python validation/reference/biopython/run_feature.py --feature io.mmcif.parse --corpus smoke --check-deps
```

Generate goldens for a feature:

```bash
micromamba run -n molecules-rdkit-reference python validation/reference/rdkit/run_feature.py --feature io.sdf.v2000.parse --corpus smoke
micromamba run -n molecules-rdkit-reference python validation/reference/rdkit/run_feature.py --feature algo.rings.fast --corpus smoke
micromamba run -n molecules-rdkit-reference python validation/reference/rdkit/run_feature.py --feature algo.aromaticity.rdkit-like --corpus smoke
micromamba run -n molecules-biopython-reference python validation/reference/biopython/run_feature.py --feature io.mmcif.parse --corpus smoke
micromamba run -n molecules-biopython-reference python validation/reference/biopython/run_feature.py --feature bio.hierarchy.smcra --corpus smoke
```

By default, output goes to `validation/corpora/<corpus-id>/golden/<feature-id>/`. Use `--fixture`
to limit generation to a listed fixture, and `--output-dir` to write elsewhere for review.

Golden files should be reviewed before committing. Creating a golden file does not automatically
record validation evidence; run the Rust comparison with `--update` only after the corpus is ready.

Do not create molecule fixtures by hand for reference validation. Add compact records from external
sources under the corpus `data/` directory, record their source URL and SHA-256 in
`sources.lock.json`, then generate goldens with RDKit or Biopython.

The PubChem builder uses the official `CID-SMILES.gz` snapshot and the first `CURRENT-Full` SDF
shard. Selection remains seeded and deterministic; the shard constraint and source checksums are
recorded in corpus metadata and locks.
