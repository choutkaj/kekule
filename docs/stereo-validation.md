# Reproducing stereo validation

RDKit is an optional scientific reference, never a Rust runtime dependency. The
current comparison pins RDKit **2026.03.6** and calls the accurate
`AssignCIPLabels` implementation. It does not use the legacy `_CIPRank` algorithm.
The existing `stereo.cip` corpus snapshots retain their original 2026.03.3 target.

## Differential CIP runner

Build the adapter and check its record-preservation contract:

```text
cargo build -p xtask --example cip_probe --locked
cargo test -p xtask --example cip_probe --locked
uv run --with rdkit==2026.3.6 --python 3.13 python -m unittest discover -s benchmarks/reference/rdkit -p test_compare_cip.py
```

With the external corpus data already installed, run from the repository root:

```text
uv run --python 3.13 benchmarks/reference/rdkit/compare_cip.py --corpus enamine-diversity --stereo-only --mode sanitized --probe target/debug/examples/cip_probe --output target/stereo-production/enamine-2026.03.6.json
uv run --python 3.13 benchmarks/reference/rdkit/compare_cip.py --corpus pubchem-100k --stereo-only --mode sanitized --probe target/debug/examples/cip_probe --output target/stereo-production/pubchem-2026.03.6.json
```

On Windows append `.exe` to the probe path. Outputs must be new files: the runner
refuses to overwrite previous evidence. Use `--input <file>` for other external
SMILES files. `--input-format json` reads the `cases[].smiles` input field of a
previous audit and computes fresh reference results; it never imports expected
descriptors. `--input-format suite` reads the published validation suite's tabular
SMILES, case IDs, and published labels. Published labels are retained as separate
evidence, not substituted for the RDKit result.

Each input SMILES string appears once per run, including disconnected structures.
Without `--stereo-only`, stereo-free controls are included. Syntax, interpretation,
perception, CIP, process, protocol, timeout, and reference errors remain visible
records. A reference failure is never a successful comparison. Empty selections
are errors. The process returns nonzero for any error or descriptor difference.

The report stores source and binary SHA-256 checksums, versions, limits, and both
complete descriptor maps. Atom counts, bond counts, label absence, and lowercase
descriptors are asserted. Atom indices include explicit hydrogen vertices and all
components. Bonds are identified by their sorted endpoint indices, avoiding
toolkit-specific bond insertion order. Reference-only runs are explicitly counted
as `reference_only`, never `match`.

## Reference preparation is part of the experiment

Both modes parse unsanitized SMILES, retain explicit hydrogen vertices, run modern
stereo perception with `cleanIt=False`, clear previous CIP properties, and invoke
`AssignCIPLabels(maxRecursiveIterations=1000000)`.

| Mode | Sanitization | Purpose |
| --- | --- | --- |
| `sanitized` | `SANITIZE_ALL` | Reproduce the original broad-corpus pipeline |
| `assertions` | Exclude `SANITIZE_CLEANUPCHIRALITY` | Keep supplied tetrahedral tags for a closer comparison of asserted stereo |

Parsed, sanitized, and prepared atom tags are stored separately. Neither mode can
undo a tag transformation that RDKit already applied while parsing SMILES. VS132
demonstrates why a disagreement must be traced through interpretation before it is
attributed to CIP ranking. These modes describe observable procedures; neither is
a blanket assertion that every input format has identical semantics in both tools.

`benchmarks/reference/rdkit/vs132_reproducer.py` independently ranks the ligands
in the published 3D structure and transports its local configurations onto the
original SMILES graph. The fixture and its source hashes are retained in
`crates/kekule/tests/fixtures/cip/`. Run the reproducer with
`uv run --with rdkit==2026.3.6 --python 3.13 python benchmarks/reference/rdkit/vs132_reproducer.py`.

The [published CIP Validation Suite](https://cipvalidationsuite.github.io/ValidationSuite/)
is pinned to revision `6b9f9db46dadc6749da8234b05164e1e0fb413b9`.
Its [SMILES source](https://raw.githubusercontent.com/CIPValidationSuite/ValidationSuite/6b9f9db46dadc6749da8234b05164e1e0fb413b9/compounds.smi)
and [3D SDF source](https://raw.githubusercontent.com/CIPValidationSuite/ValidationSuite/6b9f9db46dadc6749da8234b05164e1e0fb413b9/compounds_3d.sdf)
are external fixtures. Retain the downloaded bytes and the runner's input hashes.
Allene/cumulene geometries and wildcard atoms remain explicit unsupported inputs;
they must not be deleted from a differential report to obtain a passing total.

## Molfile interchange

The cross-tool check uses externally supplied PubChem alkenes and a sourced RDKit
atropisomer fixture. It compares complete indexed chemical graphs, accurate CIP
labels, and enhanced groups after RDKit rereads both output versions. It includes
cleared E/Z assertions and V3000 atom CFG. Source and binary hashes plus every
request and response are retained; an existing report is never overwritten.

```text
cargo build -p xtask --example molfile_interchange_probe --locked
uv run --with rdkit==2026.3.6 --python 3.13 python benchmarks/reference/rdkit/compare_molfile_interchange.py --probe target/debug/examples/molfile_interchange_probe --output target/stereo-production/molfile-interchange.json
```

On Windows append `.exe` to the probe path. Atom CFG comparisons explicitly call
RDKit's parity-conversion helper; ordinary wedge output is tested through the
normal Molfile reader. Explicit-H CFG ordering follows CTfile and has separate
Rust regressions because that RDKit helper ignores its hydrogen-last exception.

## Rust and Linux validation

The applicable Rust checks are listed in `.github/workflows/ci.yml`. In addition
to the workspace checks, run the adapter's example tests above. Documentation is
checked with `RUSTDOCFLAGS=-D warnings`. Dirty package validation uses
`--allow-dirty`; this does not waive compilation or package-content checks.

Linux fuzz smoke follows the CI configuration with nightly Rust and
`cargo-fuzz 0.13.2`:

```bash
set -euo pipefail
seed=0
while read -r target; do
    seed=$((seed + 1))
    cargo +nightly fuzz run "$target" -- -runs=256 -max_len=4096 -seed="$seed"
done < <(cargo +nightly fuzz list)
```

When using WSL alongside Windows builds, set `CARGO_TARGET_DIR` to a separate
workspace directory such as `target/linux-fuzz`. Fuzz smoke is bounded execution,
not proof that no malformed input can fail.

## Isomeric writer representation checks

RDKit's isomeric writer brackets every atom bonded to a metal and writes its
total hydrogen count inside those brackets. This changes the explicit/implicit
hydrogen split, `noImplicit`, and explicit valence even when the chemical graph
is unchanged. For example, it writes `O[Fe]=O` as `[OH][Fe]=[O]`. Kekule preserves
the source hydrogen declaration when the SMILES grammar permits it. The policy
is visible in the pinned [RDKit 2026.03.6 writer source](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/SmilesParse/SmilesWrite.cpp#L138).

The targeted external check retains all source declaration fields, RDKit's
literal emitted SMILES and complete decoded projection, and Kekule's emitted
SMILES independently decoded by RDKit. It checks the whole chemical graph using
RDKit's canonical isomeric representation and separately asserts the complete
minimum legal source projection. Normalization can introduce a formal charge on
an atom that originally allowed implicit hydrogens. SMILES then requires brackets;
the target fixes the independently perceived total hydrogen count only for those
atoms. The report retains the original fields and lists every required change.
The check does not call Kekule's canonical writer or require the two writers to
emit identical strings.

```text
cargo build -p xtask --example smiles_write_probe --locked
cargo test -p xtask --example smiles_write_probe --locked
uv run --python 3.13 benchmarks/reference/rdkit/isomeric_projection_reproducer.py --probe target/debug/examples/smiles_write_probe --output target/stereo-production/metal-writer-projection.json
```

On Windows append `.exe` to the probe path. Eleven PubChem input records cover
nine optional metal-bracketing cases and two chlorine-normalization cases. Their
source URLs, CIDs, hashes, and corpus locations are retained in
`benchmarks/reference/rdkit/fixtures/metal-writer-projection.json`. The check pins
RDKit 2026.03.6, refuses to overwrite an existing report, and records the fixture
and probe hashes. Omitting the probe performs a reference-only check and reports
that distinction explicitly.
