# Benchmark Reference Generators

Reference generators reproduce normalized JSON goldens from external tools and
externally supplied, provenance-pinned fixtures. They are development tools,
not Rust runtime dependencies.

Create the pinned environments:

```bash
micromamba create -f benchmarks/reference/rdkit/environment.yml
micromamba create -f benchmarks/reference/biopython/environment.yml
```

Check dependencies or generate selected goldens:

```bash
micromamba run -n rdkit-benchmark-reference python benchmarks/reference/rdkit/run_feature.py --feature io.sdf.v2000.parse --corpus pubchem-1k --check-deps
micromamba run -n biopython-dssp-benchmark-reference python benchmarks/reference/biopython/run_feature.py --feature io.mmcif.parse --corpus pdb-100 --check-deps
micromamba run -n rdkit-benchmark-reference python benchmarks/reference/rdkit/run_feature.py --feature algo.aromaticity.rdkit-like --corpus pubchem-1k
micromamba run -n biopython-dssp-benchmark-reference python benchmarks/reference/biopython/build_dssp_benchmark.py --corpus pdb-100 --jobs 4
```

Build normal corpora from their checked-in locked membership:

```bash
micromamba run -n rdkit-benchmark-reference python benchmarks/reference/rdkit/build_corpus.py
micromamba run -n biopython-dssp-benchmark-reference python benchmarks/reference/biopython/build_corpus.py
```

The builders use source-lock entry order and categories by default. Their
explicit reselection options require new repository-neutral corpus version
identifiers.

Output defaults to
`benchmarks/corpora/<corpus-id>/golden/<benchmark-id>/`. Use `--fixture` to limit
generation or `--output-dir` for separate review. Golden changes require
independent review; running a generator does not claim repository health.
