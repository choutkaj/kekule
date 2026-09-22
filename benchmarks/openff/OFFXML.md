# OFFXML generalization validation, 2026-09-22

`ForceField::from_file(path)` and `ForceField::from_offxml(xml)` compile custom
rules within the [supported subset](../../crates/kekule-openff/CONTRACT.md#custom-offxml-subset).
`ForceField::rosemary()` remains a preset using that same compiler. Ash model
loading, fingerprints, features, lookup semantics, charge precedence, and the
parameterization API are unchanged.

This pass generalizes file loading, handler defaults, optional sections, legacy
headers, rule metadata, and unit expressions. It does not add alternate charge
methods, OFFXML merging, virtual sites, fractional-bond-order parameters, or new
potential energy functions.

## Independent parser reference

[reference_offxml.py](reference_offxml.py) loads the externally supplied,
checksum-pinned Rosemary document using OpenFF Toolkit **0.19.0**, applies the
explicit transformations below, and freezes Toolkit's interpretation. No native
executable participates in reference generation. The fixture contains each XML
variant, every ordered parameter, and every exposed nonbonded setting in
canonical units. Source and fixture hashes are recorded in
[offxml-generalization.lock.json](offxml-generalization.lock.json).

| Reference variant | What it exercises |
| --- | --- |
| Original Rosemary | Every original rule and nonbonded setting remains compatible |
| Specification defaults | Omitted functional forms, mixing rule, scales, cutoffs, switch widths and methods |
| Legacy 0.3, Coulomb | Bonds/proper-torsion versions and old nonbonded method attributes |
| Legacy 0.3, default PME | Default electrostatic method converts to Ewald3D-ConductingBoundary |
| Optional sections absent | No Constraints, ImproperTorsions, or LibraryCharges rules |
| Compact units and anonymous override | Whitespace-free units, division, parentheses, scientific notation, plural units, an explicit bond-value change and an omitted rule ID |

All six variants pass the Rust comparison of every parsed rule and setting.
Parameter values use the existing `atol=1e-10`, `rtol=1e-12` policy; structure,
rule order, SMIRKS, IDs, counts, and method names compare exactly. Automatic
proper-torsion division retains the engine's deferred-resolution marker; the
reference stores that explicit convention rather than inventing a numeric
division factor before a molecule is assigned.

The omitted-section variant also assigns all molecules from the original
23-input external audit. Its bonds, angles, proper torsions and vdW labels match
the preset, while its constraint and improper label sets are empty. Regressions
reject malformed expressions, dimension mismatches, nonfinite results,
excessive nesting/exponents, unsupported charge models, unknown attributes and
versions, and mixtures of legacy and current method syntax. File loading checks
paths containing spaces, I/O errors, malformed XML, invalid UTF-8 and contextual
diagnostics.

The defaults follow the [SMIRNOFF specification](https://openforcefield.github.io/standards/standards/smirnoff/)
and the installed Toolkit handler implementation. Omitted electrostatic
`scale14` means **0.833333**, preserving the reference's decimal value; the
explicit Rosemary value remains **0.8333333333**. Defaults do not silently
substitute a different charging model.

## Reproduction

Generate a new independent fixture in the pinned reference environment. The
generator refuses to overwrite an existing output:

```text
micromamba run -p target/openff-reference python benchmarks/openff/reference_offxml.py --output target/new-offxml-reference.json.gz
cargo test -p kekule-openff --locked
python -m unittest discover -s benchmarks/openff -p "test_*.py"
```

The archived fixture is
[offxml-generalization.json.gz](../../crates/kekule-openff/tests/fixtures/offxml-generalization.json.gz).
Rust tests compare the archived interpretation without Python or OpenFF installed.
The Python provenance test verifies byte fingerprints. External reference runs
remain optional scientific validation, rather than a routine CI requirement.

The main [validation report](VALIDATION.md) retains its original observations,
figures and timing measurements. This parser extension does not regenerate
historical reference data.

## Completed checks

The rebuilt engine also passed all **66 live Rosemary reference cases** in the
original small-molecule suite, including parameters, atom-order permutations,
features, lookup charges, and forced inference. The new observations are saved
separately in [offxml-rosemary-recheck.json.gz](offxml-rosemary-recheck.json.gz)
and fingerprinted in the lock file. The maximum final charge difference remains
`2.9996037483215332e-5 e`, within the unchanged `5e-5 e` tolerance.

The following checks passed on Windows:

```text
cargo fmt --all -- --check
cargo check --workspace --all-targets --all-features --locked
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo test --workspace --all-features --locked
cargo test -p kekule-openff --test contracts --locked -- --ignored
cargo +1.89.0 check --workspace --all-targets --all-features --locked
cargo doc --workspace --all-features --no-deps --locked
cargo package -p kekule-openff --allow-dirty --locked --list
cargo package -p kekule-bench --allow-dirty --locked --list
python -m unittest discover -s benchmarks/openff -p "test_*.py"
git diff --check
```

Workspace tests report **1,415 passed**, zero failed, and one model-dependent
test ignored; that test passed separately with `KEKULE_OPENFF_MODEL` pointing to
the exported bundle. All **23 Python tests** pass. Rustdoc used
`RUSTDOCFLAGS=-D warnings`.

The initial workspace run failed the unrelated trajectory test
`path_writer_failure_poisoning_prevents_partial_publication` at its output-file
existence assertion. Both its entire test binary and the subsequent full
workspace rerun passed, without trajectory code changes. Initial and rerun logs
remain in `target/openff-generalization-*.log`.

Not repeated: the 110-molecule energy panel (its immutable observations are
retained; this pass checks parser interpretations and replays the original 66
cases), a separate workspace doctest command (included by workspace tests),
unchanged core/other-companion package checks and no-default-feature checks,
fuzz execution, or Linux/macOS runs. Full companion package verification remains
deferred until the core 0.2.1 dependency is published; package file lists were
checked, matching the existing CI policy.
