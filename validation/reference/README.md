# Reference Generators

Reference generators produce normalized JSON golden data from external tools.

These scripts are validation infrastructure only:

- RDKit is used for small-molecule parser, ring, and aromaticity reference output.
- Biopython is used for mmCIF atom-site and SMCRA hierarchy reference output.
- Neither tool is a Rust runtime dependency.

Run dependency checks from the repository root after creating the matching conda environment:

```bash
python validation/reference/rdkit/run_feature.py --feature io.sdf.v2000.parse --check-deps
python validation/reference/biopython/run_feature.py --feature io.mmcif.parse --check-deps
```

Generate goldens for a feature:

```bash
python validation/reference/rdkit/run_feature.py --feature io.sdf.v2000.parse
python validation/reference/rdkit/run_feature.py --feature algo.rings.fast
python validation/reference/rdkit/run_feature.py --feature algo.aromaticity.rdkit-like
python validation/reference/biopython/run_feature.py --feature io.mmcif.parse
python validation/reference/biopython/run_feature.py --feature bio.hierarchy.smcra
```

By default, output goes to `validation/features/<feature-id>/golden/`. Use `--fixture` to limit
generation to a listed fixture, and `--output-dir` to write elsewhere for review.

Golden files should be reviewed before committing. Creating a golden file does not automatically make
a feature validated; update feature metadata only after the validation criteria are actually met.
