# Reference Generators

Reference generators produce normalized JSON goldens from external tools using externally supplied fixtures declared in feature validation manifests.

- RDKit supplies small-molecule parser, writer, perception, ring, stereo, and aromaticity reference output.
- Biopython supplies format-level mmCIF atom-site rows and, together with `mkdssp`, DSSP reference output.
- Neither tool is a Rust runtime dependency.

Create the matching micromamba environments from the repository root:

```bash
micromamba create -f validation/reference/rdkit/environment.yml
micromamba create -f validation/reference/biopython/environment.yml
```

Run dependency checks through those environments:

```bash
micromamba run -n molecular-rdkit-reference python validation/reference/rdkit/run_feature.py --feature io.sdf.v2000.parse --corpus pubchem-1k --check-deps
micromamba run -n molecular-biopython-reference python validation/reference/biopython/run_feature.py --feature io.mmcif.parse --corpus pdb-100 --check-deps
```

Generate feature goldens:

```bash
micromamba run -n molecular-rdkit-reference python validation/reference/rdkit/run_feature.py --feature io.sdf.v2000.parse --corpus pubchem-1k
micromamba run -n molecular-rdkit-reference python validation/reference/rdkit/run_feature.py --feature algo.aromaticity.rdkit-like --corpus pubchem-1k
micromamba run -n molecular-biopython-reference python validation/reference/biopython/run_feature.py --feature io.mmcif.parse --corpus pdb-100
micromamba run -n molecular-biopython-reference python validation/reference/biopython/build_dssp_validation.py --corpus pdb-100 --jobs 4
```

Construct or refresh the nested macromolecular tiers with:

```bash
micromamba run -n molecular-biopython-reference python validation/reference/biopython/build_corpus.py
```

The PDB builder intersects the current RCSB holdings with typed RCSB Search API candidate pools, ranks each category deterministically from the corpus seed, verifies the downloaded mmCIF records, and preserves PDB-100 as the exact prefix of PDB-1000.

DSSP generation defaults to four independent process workers. Provisioned hosts may override the bound explicitly with `--jobs N`.

By default, output goes to `validation/corpora/<corpus-id>/golden/<feature-id>/`. Use `--fixture` to limit generation to a listed fixture and `--output-dir` to write elsewhere for review.

Golden files should be reviewed before committing. Creating a golden does not record validation evidence; run the Rust comparison with `--update` only after the corpus is ready.

Do not create molecular reference fixtures by hand. Add externally sourced records under the corpus `data/` directory, pin their source URL and SHA-256 in `sources.lock.json`, then generate goldens with the matching reference environment.