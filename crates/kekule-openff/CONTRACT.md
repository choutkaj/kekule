# SMIRNOFF parameterization and the Rosemary preset

`ForceField::rosemary()` compiles the packaged
`openff_no_water-3.0.0-alpha2b.offxml` (2026-09-11). This is the pinned Rosemary
prerelease, not a promise to track future releases. The supported SMIRNOFF
handlers are Bonds, Angles, ProperTorsions, ImproperTorsions, Constraints, vdW,
Electrostatics, LibraryCharges and NAGLCharges. It does not include a water model.

`ForceField::from_file(path)` reads a user-supplied UTF-8 OFFXML file;
`ForceField::from_offxml(xml)` compiles a string. Paths are literal, without
force-field registry resolution or downloads. File-loading errors include the
path. Both entry points use the same compiler as the Rosemary preset.

```rust,no_run
let force_field = kekule_openff::ForceField::from_file("my-forcefield.offxml")?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Custom OFFXML subset

The root must declare SMIRNOFF `0.3` and `OEAroModel_MDL`. Charge-model decisions
are unchanged: `NAGLCharges` version `0.3` must declare the pinned Ash filename
and hash, and parameterization still takes `&NaglModel`. Other charge methods,
model loading, and user charge overrides remain unsupported.

| Handler | Accepted versions | Presence |
| --- | --- | --- |
| Bonds, ProperTorsions | 0.3, 0.4 | Required |
| Angles | 0.3 | Required |
| vdW, Electrostatics | 0.3, 0.4 | Required |
| Constraints, ImproperTorsions, LibraryCharges | 0.3 | Optional |
| NAGLCharges | 0.3 | Required; pinned Ash |

Absent optional handlers produce no corresponding rules. Missing required
handlers remain errors, including for partial force-field fragments. Files are
loaded individually; merging multiple OFFXML files is not implemented.

Supported attributes may omit their specification defaults: harmonic bonds and
angles, periodic torsions with automatic division factors, Lennard–Jones 12-6
with Lorentz–Berthelot mixing, nonbonded scales, cutoffs, switch widths and
methods. Bonds also accept the equivalent `(k/2)*(r-length)^2` spelling.
The default electrostatic 1–4 scale is exactly `0.833333`; an explicit value is
preserved. Default cutoffs are 0.9 nm, vdW switch width 0.1 nm, and electrostatic
switch width zero. An explicit positive switch width is retained as settings;
this does not claim support for evaluating electrostatic switching.

Version 0.3 vdW accepts only `method="cutoff"` (also the default), mapped to
periodic cutoff and nonperiodic no-cutoff. Version 0.3 electrostatics accepts
`method="PME"` (default) or `"Coulomb"`; PME maps to
`Ewald3D-ConductingBoundary`. The version 0.4 method attributes remain separate;
mixing old and new attribute spellings is rejected. vdW supports cutoff and
no-cutoff settings; electrostatics supports Ewald/Coulomb periodic settings and
Coulomb nonperiodic/exception settings. These describe backend policy; the
validated benchmark energy consumer is still vacuum NoCutoff.

Parameter `id` is optional and represented as an empty source ID when absent.
SMIRKS and the physical parameter values remain required. Cosmetic `name`,
`parent_id`, and `description` attributes are accepted without affecting rule
precedence; unknown attributes are rejected.

Unit expressions support whitespace-independent multiplication/division,
parentheses, signed integer powers, scientific notation, and the existing
energy/length/angle/charge units plus common plural and `nm`, `kcal`, `kJ`, `mol`
aliases. For example, `2*kcal/(mol*angstrom**2)` is accepted. Expressions are
parsed without evaluation of code, limited to 4,096 ASCII bytes, 32 nesting
levels and exponents from −32 through +32. Unknown units, dimensional
mismatches, nonfinite results and unsupported operators fail explicitly.
Fractional-bond-order indexed parameters remain unsupported; valid but unused
`none`/`AM1-Wiberg` and linear-interpolation headers do not trigger QM work.

Six independently generated Toolkit fixtures check every loaded parameter and
nonbonded setting, including legacy headers, defaults, anonymous rules, and
omitted optional sections. See [OFFXML validation](../../benchmarks/openff/OFFXML.md).

## API and ownership

```rust,no_run
use kekule::{hydrogens, smiles};
use kekule_openff::{ForceField, NaglModel};

let mut molecule = smiles::to_molecules("CCO")?.remove(0);
molecule.perceive()?;
hydrogens::add_hydrogens(&mut molecule)?;
let model = NaglModel::load("target/openff-ash")?;
let parameters = ForceField::rosemary()?.parameterize_molecule(&molecule, &model)?;
# Ok::<(), Box<dyn std::error::Error>>(())
```

For systems, call `parameterize(Arc<Topology>, &model)`. Results retain that exact
snapshot; reusable definitions are assigned once and expanded to qualified
instance atoms. Dense charge and vdW arrays follow `topology.atom_ids()` order.
Input molecules must already contain all hydrogen atoms. The engine neither
adds atoms nor changes the input topology or its properties.

Private copies receive RDKit-compatible neutral valence-five nitrogen and
phosphorus cleanup, preserving total charge, followed by MDL aromaticity for
SMIRKS assignment. NAGL uses its own pinned normalization, resonance averaging,
and selected SSSR ring features. This is not a general RDKit sanitization port.

Complete parameters live in a typed `ParameterizedTopology`, with units and
rule IDs/SMIRKS. This avoids encoding angle/torsion identities and multi-term
parameters into unrelated property columns. A future property projection or
potential adapter can consume this result without changing its ownership.

## Assignment semantics

Rules apply in file order; the last matching rule wins. Matching enumerates
tagged tuples, including automorphisms, and canonicalizes handler-specific
symmetry. Every atom, bond, angle and proper torsion must receive parameters.
Improper assignments are optional. Exhausted SMARTS searches fail explicitly.

Lengths, energies, angles and charges use nm, kJ/mol, radians and e. Bond and
angle energies include the harmonic one-half factor. A Fourier term means
`k / idivf * (1 + cos(periodicity * theta - phase))`. Proper `auto` divisors use
the number of proper paths about the central bond. Each improper produces
three cyclic terms with `auto` divisor 3; **emitted dihedral tuples put the
central atom first**, matching Interchange. Label keys and SMIRKS instead put
it second. Supported improper phases are integer multiples of pi.

Constraints without a distance use the assigned equilibrium bond length.
Constrained bond/angle terms remain available; a future evaluation backend must
choose how to treat constrained degrees of freedom. Explicit pair exceptions
cover shortest graph distances 1, 2 and 3, with independent vdW/electrostatic
scales. All other pairs have scale 1, including pairs between instances.
Nonbonded cutoffs, switching widths, methods and combining rules are retained.

Charge precedence is complete LibraryCharges coverage, then the Ash lookup
table, then neural inference. A partial library assignment falls back for the
whole molecule. Lookup uses the full fixed-H InChI string and bounded graph
mapping, relaxing formal charge/bond order and finally stereo as upstream does.
It retains the table's small asymmetries rather than averaging equivalent atoms.
NAGL results receive the toolkit's uniform total-charge correction.

## Model data and native dependencies

Inference is Rust CPU code: six mean-aggregation GraphSAGE/ReLU layers, the
sigmoid readout, and charge-conserving pooling. The exact Ash 1.0.0 checkpoint,
exported tensor layout, manifest and weights are checksum-pinned. No Python,
PyTorch, RDKit, network access or pickle loading occurs in the Rust runtime.

Prepare the bundle once in the separately installed reference environment:

```text
micromamba create --override-channels -c conda-forge -p target/openff-reference -f benchmarks/openff/environment.yml
micromamba run -p target/openff-reference python benchmarks/openff/export_model.py target/openff-ash
cargo run -p kekule-openff --release --example parameterize -- target/openff-ash CCO
```

The exporter accepts `--checkpoint PATH` and verifies its hash before loading
its Python serialization. Otherwise the OpenFF resolver may download the model.
The bundle is about 13 MB (`model.json` and little-endian float32 `weights.bin`)
and is not checked into Git. Loading checks exact file sizes before allocation,
caps reads against file growth, and verifies both fingerprints. Re-export in
the pinned environment produced byte-identical files.

The `inchi` and `inchi-sys` dependencies are pinned to 0.1.4, using the official
InChI 1.07.5 C implementation. This companion therefore requires a C toolchain
at build time; core `kekule` has no InChI dependency. The adapter handles isotope
labels, tetrahedral parity and alkene stereo, including absolute stereo groups.
Relative/mixture groups and axial stereo currently fail explicitly.

## Boundaries and validation

This implements the Rosemary functional forms, not every SMIRNOFF extension or
NAGL architecture. Unsupported handlers, section versions, parameter attributes,
unit expressions, fractional-bond-order interpolation and virtual sites fail.
The OFFXML's unused AM1-Wiberg defaults do not require AM1 calculations.
Custom OFFXML must satisfy the supported subset above and use the same Ash model.
The runtime supplies parameter assignments. An optional benchmark consumer uses
`potentials` to compare vacuum energies against OpenMM; it is not a dynamics
backend or a public force-evaluation API.

Native regression tests use an unchanged, externally sourced 23-molecule audit
fixture: all rule labels, fixed-H identifiers and all 22 NAGL feature columns.
Focused tests cover precedence, constraints, library charge fallback, repeated
instances, topology identity and improper expansion. An optional model test
also checks charges against the archived official observations:

```powershell
$env:KEKULE_OPENFF_MODEL = (Resolve-Path target/openff-ash).Path
cargo test -p kekule-openff --test contracts -- --ignored
```

The extended live comparison contains those 23 molecules plus ten independently
retrieved PubChem cases, each in forward and reversed atom order: **66/66 pass**.
It checks fixed-H identifiers, all feature columns, forced neural inference,
lookup/library/inference system charges, every bonded/vdW parameter, torsion
multiplicity and improper energies on three deterministic coordinate sets.
The largest forced-inference error was `2.54e-7 e`; the largest system charge
error was `3.00e-5 e`. The user-approved absolute charge cutoff is `5e-5 e`.
Parameter labels remain exact; numerical parameters use `atol=1e-10, rtol=1e-12`.

Exhaustive identity checks against the 13,944 stored lookup entries yield
13,235 identical stored keys. Of 709 disagreements, the current official toolkit
rejects 689 inputs; 19 native identifiers agree with the current toolkit but
differ from the stored key; one input is accepted by the toolkit but rejected
by native InChI and is outside the Ash inference domain. Raw failures are retained
in the implementation reports. These are identity checks, not exhaustive charge
remapping or whole-domain validation. Broader molecules, all normalizations and
large-system performance still warrant independent reference coverage.

The [validation report](../../benchmarks/openff/VALIDATION.md) presents parameter
and energy parity figures, numerical differences, and measured CPU timings.
The subsequent [robustness panel](../../benchmarks/openff/ROBUSTNESS.md) contains
100 independently selected PubChem molecules and ten prepared PDB protein chains
(up to 2,940 atoms). All 220 atom-order parameterization cases and 660 geometry
energy comparisons pass. Charge error is at most `2.23e-7 e`; total energy error
with identical charges is at most `6.18e-9 kJ/mol`. Strict reports still flag six
standalone identifier failures on three proteins above the InChI adapter's
1,024-atom limit. Lookup safely bypasses InChI for inputs larger than the frozen
table's largest entry (11 atoms), so those proteins parameterize successfully.
Charge lookup ignores isotopic masses on a private copy, matching the reference
toolkit while preserving the input graph. The panel does not establish protein
force-field accuracy, periodic electrostatics, forces, or trajectory stability.

Limits fail rather than silently truncate: 256 atoms/one million search states
for lookup mapping, 4096 atoms for features, 200 normalization applications per
rule, one million resonance path visits, and bounded resonance state/product
queues. Higher limits or broader chemical support need explicit validation.

See [reference execution](../../benchmarks/openff/IMPLEMENTATION.md),
[attribution](THIRD_PARTY.md), and the earlier
[prerequisite audit](../../benchmarks/openff/AUDIT.md).
