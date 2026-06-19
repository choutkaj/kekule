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
micromamba run -n molecules-rdkit-reference python validation/reference/rdkit/run_feature.py --feature io.sdf.v2000.parse --check-deps
micromamba run -n molecules-biopython-reference python validation/reference/biopython/run_feature.py --feature io.mmcif.parse --check-deps
```

Generate goldens for a feature:

```bash
micromamba run -n molecules-rdkit-reference python validation/reference/rdkit/run_feature.py --feature io.sdf.v2000.parse
micromamba run -n molecules-rdkit-reference python validation/reference/rdkit/run_feature.py --feature algo.rings.fast
micromamba run -n molecules-rdkit-reference python validation/reference/rdkit/run_feature.py --feature algo.aromaticity.rdkit-like
micromamba run -n molecules-biopython-reference python validation/reference/biopython/run_feature.py --feature io.mmcif.parse
micromamba run -n molecules-biopython-reference python validation/reference/biopython/run_feature.py --feature bio.hierarchy.smcra
```

By default, output goes to `validation/features/<feature-id>/golden/`. Use `--fixture` to limit
generation to a listed fixture, and `--output-dir` to write elsewhere for review.

Golden files should be reviewed before committing. Creating a golden file does not automatically make
a feature validated; update feature metadata only after the validation criteria are actually met.

Do not create molecule fixtures by hand for reference validation. Add compact records from external
sources under `validation/external_sources/`, record their source URL and SHA-256 in
`fixture_sources`, then generate goldens with RDKit or Biopython.
