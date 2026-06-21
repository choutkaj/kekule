<p align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="./assets/molecules-logo-dark.svg">
    <img alt="MOLECULES - cheminformatics in Rust" src="./assets/molecules-logo-light.svg" width="250">
  </picture>
</p>

`molecules` is an experimental pure-Rust backend for small-molecule and
macromolecular structure work. It has no RDKit or Biopython runtime dependency;
those projects are used only to generate external validation references.

The crate is version `0.0.0`, is not published, and makes no API stability or
chemistry-completeness promise. Successful validation covers only the named
features, fixtures, fields, and reference versions recorded in the generated
[feature dashboard](features/DASHBOARD.html).

## Implemented

- Typed atom, bond, conformer, graph, property, and computed-state models.
- Raw V2000 Molfile and SDF parsing, including supported charge, isotope,
  radical, map, coordinate, and bond-stereo fields.
- V2000 Molfile and SDF writing with structured errors for unsupported
  representations.
- Raw SMILES parsing for the documented organic/bracket/aromatic subset,
  branches, fragments, and ring labels `0` through `99`.
- Deterministic noncanonical SMILES writing for the supported subset.
- Explicit valence, ring, aromaticity, and transactional sanitization passes.
- Raw mmCIF atom-site parsing into `MacroMolecule` and `BioHierarchy`, including
  models, label/author residue identity, alternate locations, and coordinates.

Parsing and sanitization are separate. The V2000, SDF, SMILES, and mmCIF
readers return raw parsed state and do not silently run chemistry perception.
Call `sanitize_small_molecule` explicitly when its current RDKit-like subset is
appropriate. A failed sanitize operation leaves the input molecule unchanged.

## Unsupported

The roadmap includes, but the implementation does not currently provide,
V3000, canonical SMILES/ranking, SMARTS/substructure search, CIP assignment,
general stereo perception, fingerprints, or descriptors. SMILES query syntax,
wildcards, atom chirality, and directional bond notation are rejected.
Writers reject unsupported bond orders or stereo instead of coercing them.

The aromaticity and valence implementations intentionally cover a tested
subset, not full RDKit behavior. A parser success is not a claim that the
result is chemically valid; sanitization can return a structured error.

## Validation

`implemented` and `validated` are different states. The dashboard records
successful comparisons, while `cargo xtask validate` recomputes content
addresses and is the authority for whether evidence is current in a checkout.
As of June 21, 2026, tiny, PubChem 100/1000, and PDB 10/100 comparisons are
built and passing locally. PL-REX and Enamine are declared but are not built
validation corpora.

```bash
cargo fmt --all -- --check
cargo check --workspace --all-targets
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
cargo xtask skills --check
cargo xtask corpus check --corpus tiny --require-data
cargo xtask validate --feature all --corpus tiny
```

Large corpus data is intentionally not committed. See
[RELEASE_READINESS.md](RELEASE_READINESS.md) for the clean-runner and release
blockers, and [FUZZING.md](FUZZING.md) for parser robustness campaigns.

## Limits And Errors

Public text parsers return structured errors for malformed input rather than
panicking. V2000 is limited by the format subset to 999 atoms and 999 bonds.
mmCIF defaults cap input at 256 MiB, tokens at 10 million, a token at 16 MiB,
and atom-site rows at 5 million. Ring perception defaults cap graph size,
candidate cycles, shortest-path work, cycle size, and total work; callers can
override `RingPerceptionOptions`.

These are defensive implementation limits, not a resource-usage guarantee for
hostile workloads. Security reports should follow [SECURITY.md](SECURITY.md).

## License

No license has been selected. The repository is private and `publish = false`.
No public release or redistribution permission should be assumed until the
owner selects a license and reviews dependency, corpus, and golden-data terms.
