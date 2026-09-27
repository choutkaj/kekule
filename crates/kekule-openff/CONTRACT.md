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

The root must declare SMIRNOFF `0.3` and `OEAroModel_MDL`. `NAGLCharges`
version `0.3` declares a nonempty `model_file` and a 64-digit SHA-256
`model_file_hash`. `ForceField::charge_model()` exposes that identity.
Parameterization takes `&NaglModel` and requires an exact model identifier and
checkpoint digest match before assignment, including for library-only inputs.
Hexadecimal digests are case-insensitive. Model identifiers are literal strings,
not paths to resolve or download. Other charge methods and implicit overrides
remain unsupported. The Rosemary preset continues to require Ash 1.0.0.

| Handler | Accepted versions | Presence |
| --- | --- | --- |
| Bonds, ProperTorsions | 0.3, 0.4 | Required |
| Angles | 0.3 | Required |
| vdW, Electrostatics | 0.3, 0.4 | Required |
| Constraints, ImproperTorsions, LibraryCharges | 0.3 | Optional |
| NAGLCharges | 0.3 | Required; matching supported NAGL bundle |

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
omitted optional sections. See [OFFXML validation](../../benchmarks/openff/VALIDATION.md).

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

Charge precedence is complete LibraryCharges coverage, then the model's optional lookup
table, then neural inference. A partial library assignment falls back for the
whole molecule. Lookup uses the full fixed-H InChI string and bounded graph
mapping, relaxing formal charge/bond order and finally stereo as upstream does.
It retains the table's small asymmetries rather than averaging equivalent atoms.
NAGL results receive the toolkit's uniform total-charge correction. Both lookup
and inference provenance retain the model identifier and original checkpoint
SHA-256. Library provenance retains the matched parameter IDs.

## Model data and native dependencies

Inference is Rust CPU code. Bundle loading, feature configuration, GraphSAGE
evaluation, and charge assignment are separate internal modules with one public
`NaglModel` type. No Python, PyTorch, RDKit, network access or pickle loading
occurs in the Rust runtime. Direct `.pt` import is not implemented.

Schema-2 bundles contain `model.json` and little-endian float32 `weights.bin`.
The manifest declares original checkpoint identity, weights checksum, complete
tensor layout, NAGL configuration, chemical domain, optional lookup tables and
preparation profile. The original checksum-pinned schema-1 Ash export remains
accepted unchanged. New exports use schema 2.

| Configuration | Supported subset |
| --- | --- |
| Preparation | `openff-nagl-0.6.1`: existing normalization, resonance and SSSR semantics |
| NAGL config | Version `0.1`, no bond features, exactly one charge readout |
| Atom features | Element and connectivity one-hot categories in declared order; average formal charge; ring sizes 3–6 |
| Categories | Unique canonical element symbols; unique connectivity values in 0–6; missing input categories fail |
| Convolution | SAGEConv with mean aggregation, 1–16 layers, widths 1–2048 |
| Activations | ReLU, sigmoid, identity |
| Readout | Atom pooling, 0–16 hidden layers, widths 1–2048, regularized charge equilibration with three outputs |
| Dropout | Training probability 0–1 accepted; always disabled during evaluation |
| Lookup | Zero or one table matching the readout name; duplicate keys rejected; map/charge correspondence checked on a hit |

Feature lists may contain up to 64 entries and produce up to 256 columns. The
loader caps metadata at 64 MiB and weights at 256 MiB, including file-growth
checks. It validates tensor shapes, offsets, complete nonoverlapping storage,
finite weights and checksums. Unknown fields, schemas, profiles, features or
operations fail explicitly. An empty domain element list means unrestricted
elements at the domain-check stage; feature categories still must cover inputs.

Prepare an Ash bundle once in the reference environment (output must be new):

```text
micromamba run -p target/openff-reference python benchmarks/openff/scripts/export_model.py target/ash-bundle
cargo run -p kekule-openff --release --example parameterize -- target/ash-bundle CCO
```

For another trusted checkpoint, supply `--checkpoint PATH` and
`--checkpoint-sha256 SHA256`; optionally supply its `--license PATH`. The
exporter verifies the checkpoint before Python deserialization and records its
configuration. The Rust loader determines whether that configuration is
supported. The default Ash resolver may download its model; custom checkpoint
paths are local. See [model validation and reproduction](../../benchmarks/openff/VALIDATION.md).

Bundle checksums detect corruption; the declared original checkpoint hash is
provenance from the exporter, not proof that arbitrary edited bundles reproduce
that checkpoint. Distribute bundles with independently pinned fingerprints and
retain source licenses. The two validation bundles are externally supplied and
are not checked into Git. Supporting their configuration does not establish
scientific accuracy for arbitrary user-trained weights.

The `inchi` and `inchi-sys` dependencies are pinned to 0.1.4, using the official
InChI 1.07.5 C implementation. This companion therefore requires a C toolchain
at build time; core `kekule` has no InChI dependency. The adapter handles isotope
labels, tetrahedral parity and alkene stereo, including absolute stereo groups.
Relative/mixture groups and axial stereo currently fail explicitly.

## Boundaries and validation

The subsequent two-model end-to-end suite checks Ash 1.0.0 (22 features,
13,944 lookup entries) and OpenFF `0.1.0-rc.2` (21 features, no lookup table).
All **132 cases** pass: 33 externally sourced molecules in both atom orders for
each model. It checks complete parameterization, every valence/vdW value and
multiplicity, all feature columns, direct/assigned charges, model provenance,
input immutability and rejection of mismatched models or malformed bundles.
Forced-inference tolerance is `1e-6 e`; final-charge tolerance remains `5e-5 e`.
Measured maxima and commands are recorded in [VALIDATION.md](../../benchmarks/openff/VALIDATION.md).

This implements the Rosemary functional forms, not every SMIRNOFF extension or
NAGL architecture. Unsupported handlers, section versions, parameter attributes,
unit expressions, fractional-bond-order interpolation and virtual sites fail.
The OFFXML's unused AM1-Wiberg defaults do not require AM1 calculations.
Custom OFFXML must satisfy the supported subset above and declare the supplied
model identity. New charge methods and new network architectures require their
own implementations and reference validation.
The runtime supplies parameter assignments. An optional benchmark consumer uses
`potentials` to compare vacuum energies against OpenMM; it is not a dynamics
backend or a public force-evaluation API.
That consumer also validates Cartesian energy gradients, numerical derivatives
and rigid transformations. The 110-molecule panel retains one near-linear
geometry with a documented OpenMM force-regularization discrepancy and a strict
native torsion rotation failure; see the validation report before interpreting
energy parity as evidence of general force-evaluation robustness.

Native regression tests use an unchanged, externally sourced 23-molecule audit
fixture: all rule labels, fixed-H identifiers and all 22 NAGL feature columns.
Focused tests cover precedence, constraints, library charge fallback, repeated
instances, topology identity and improper expansion. An optional model test
also checks charges against the archived official observations:

```powershell
$env:KEKULE_OPENFF_MODEL = (Resolve-Path target/openff-ash).Path
cargo test -p kekule-openff --test contracts -- --ignored
```

The original Ash live comparison contains those 23 molecules plus ten independently
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
The panel contains
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

See the report for reference execution commands and the historical prerequisite
audit, and [THIRD_PARTY.md](THIRD_PARTY.md) for attribution.
