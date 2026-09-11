# Stereo validation

RDKit is an optional scientific reference, never a Rust dependency or release
gate. Differential checks pin **2026.03.6** and use accurate `AssignCIPLabels`,
with modern perception and a one-million-iteration bound. The implementation
follows [IUPAC P-9](https://iupac.qmul.ac.uk/BlueBook/P9.html) and the
[Hanson et al. refinements](https://doi.org/10.1021/acs.jcim.8b00324) used by
[RDKit's labeler](https://github.com/rdkit/rdkit/tree/Release_2026_03_6/Code/GraphMol/CIPLabeler).
See the [support contract](stereo-support.md) and
[benchmark contract](../benchmarks/reference/stereo/SCHEMA.md).

## Running the checks

From the repository root, with external corpus data installed:

```text
cargo build -p xtask --examples --locked
cargo test -p xtask --examples --locked
uv run --with rdkit==2026.3.6 --python 3.13 python -m unittest discover -s benchmarks/reference/rdkit -p "test_*.py"
uv run --python 3.13 benchmarks/reference/rdkit/compare_cip.py --corpus enamine-diversity --stereo-only --mode sanitized --probe target/debug/examples/cip_probe --output target/stereo-validation/enamine.json
uv run --python 3.13 benchmarks/reference/rdkit/compare_cip.py --corpus pubchem-100k --stereo-only --mode sanitized --probe target/debug/examples/cip_probe --output target/stereo-validation/pubchem.json
uv run --with rdkit==2026.3.6 --python 3.13 python benchmarks/reference/rdkit/compare_molfile_interchange.py --probe target/debug/examples/molfile_interchange_probe --output target/stereo-validation/molfile.json
uv run --python 3.13 benchmarks/reference/rdkit/isomeric_projection_reproducer.py --probe target/debug/examples/smiles_write_probe --output target/stereo-validation/isomeric.json
```

Append `.exe` to probe paths on Windows. Choose fresh output paths for each run;
reports preserve input/probe hashes and refuse overwrites. Rust, documentation,
package, and bounded Linux fuzz checks follow `.github/workflows/ci.yml`.
Use a separate `CARGO_TARGET_DIR` for WSL builds alongside Windows builds.

## What is compared

The CIP runner asserts complete descriptor maps, including absent and lowercase
labels, plus atom/bond counts. Indices retain explicit hydrogen vertices and all
components; bond keys use sorted endpoint indices. Identical input SMILES are
deduplicated. Failures and timeouts remain records, and any error or difference
causes a nonzero exit. Omitting `--stereo-only` includes stereo-free controls.
`--input <file>` accepts external SMILES; `--input-format json` takes input SMILES
from a report and recomputes expectations. `--input-format suite` also retains
published labels as independent evidence.

Both reference modes parse unsanitized SMILES, retain explicit H, perceive stereo
with `cleanIt=False`, clear previous CIP properties, and run the accurate labeler.
`sanitized` applies `SANITIZE_ALL`; `assertions` excludes
`SANITIZE_CLEANUPCHIRALITY` to retain supplied tags. Parsed, sanitized, and prepared
tags are recorded separately. Neither mode undoes changes made by the parser.

The Molfile checker compares complete indexed chemical graphs, CIP labels, and
enhanced groups after RDKit rereads V2000/V3000 output. External PubChem alkenes and
an RDKit atropisomer fixture also exercise cleared E/Z and V3000 atom CFG. Every
request, output, and explicit unknown assertion is checked.

The isomeric checker independently decodes eleven PubChem structures, asserting
the complete chemical graph and source hydrogen declarations. It covers nine
optional metal-bracketing cases and two charge-normalization cases. Source URLs,
CIDs, and hashes live beside the fixture. This check does not call Kekule's
canonical writer or require identical emitted strings.

## Known reference differences and limits

- **VS132 (Troger's base):** the published 3D structure and Kekule give S/S.
  RDKit's default sanitization removes the nitrogen tags; retaining them gives
  R/S from the published SMILES. Independent ligand ordering and signed volumes,
  then transporting those configurations onto the same RDKit graph, give S/S.
  Run `uv run --with rdkit==2026.3.6 --python 3.13 python benchmarks/reference/rdkit/vs132_reproducer.py`;
  its sourced fixture is in `crates/kekule/tests/fixtures/cip/`.
- **V3000 explicit-H CFG:** [CTfile Appendix A](https://www.wincept.eu/toxlab/pdf/ctfile.pdf)
  puts hydrogen last in carrier order regardless of atom-row position. RDKit's
  optional `AssignAtomChiralTagsFromMolParity` ignores that exception. Kekule
  follows CTfile. The checker explicitly invokes the helper for CFG cases;
  ordinary wedge output uses RDKit's normal reader.
- **SMILES hydrogen declarations:** RDKit optionally brackets metal neighbors,
  changing explicit/implicit H policy without changing the molecular graph.
  Kekule preserves legal source declarations. Charge introduced by normalization
  requires brackets, so the expected declaration fixes the perceived total H on
  those atoms. Both the original source and RDKit's full output remain evidence.
- **Unavailable/unsupported cases:** wildcard/query atoms and allene/cumulene
  geometries remain explicit errors. RDKit also rejects ten axial suite inputs
  and exceeds the iteration bound on VS009/VS226 in assertion mode. Such results
  must remain visible and cannot count as matches.

External suite inputs are pinned to
[revision 6b9f9db](https://github.com/CIPValidationSuite/ValidationSuite/tree/6b9f9db46dadc6749da8234b05164e1e0fb413b9):
[SMILES](https://raw.githubusercontent.com/CIPValidationSuite/ValidationSuite/6b9f9db46dadc6749da8234b05164e1e0fb413b9/compounds.smi)
and [3D SDF](https://raw.githubusercontent.com/CIPValidationSuite/ValidationSuite/6b9f9db46dadc6749da8234b05164e1e0fb413b9/compounds_3d.sdf).
Finite comparisons establish tested coverage, not universal chemical completeness.
