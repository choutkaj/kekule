# kekule-openff: numerical validation

Experimental Rust implementation of the OpenFF parameterization engine.
The Rosemary parameterization port reproduces parameter assignments and charges
for **100 small molecules and 10 protein chains**, with energy agreement checked
at three geometries per molecule and in two atom orders.
A separate panel of [four larger protein chains](#large-protein-chains-2026-10-09),
with 4,615–13,818 atoms, validates parameterization beyond the former size limit.

**Revalidated 2026-09-27** after the configurable NAGL model changes. Fresh OpenFF
and OpenMM observations reproduce the previous reference parameters, labels,
features, charges, exceptions and energies exactly. Native comparisons retain
the same numerical maxima and six standalone identifier diagnostics. CPU timing
now covers all 110 inputs. [Raw results and fingerprints](results/manifest.json).

| Validation measure | Result |
| --- | ---: |
| Parameterization cases passing all parameter, feature, charge and energy comparisons | **220 / 220** |
| Geometry comparisons | **660 / 660** |
| Paired atomic charges and vdW parameters | **29,252** |
| Maximum partial-charge difference | **2.229 × 10⁻⁷ e** |
| Maximum total-energy difference, identical charges | **6.170 × 10⁻⁹ kJ/mol** |
| Maximum total-energy difference, independently assigned charges | **6.011 × 10⁻⁴ kJ/mol** |
| Strict cases, including standalone InChI diagnostics | **214 / 220** |
| Strict gradient cases, including rigid transformations | **218 / 220** |
| Numerical directional energy-derivative checks | **25,410 / 25,410** |


## Parameter and charge parity

![Parity and residual plots for Lennard–Jones sigma, epsilon, and atomic partial charges.](figures/parameter-parity.png)

**Figure 1.** Kekule versus OpenFF Lennard–Jones parameters and final partial
charges. Dashed lines in panels a–c indicate exact agreement; panels d–f show
signed residuals, **Δ = Kekule − OpenFF**, at expanded scales. Every paired atom
from both atom orders is included. Repeated parameter values overlap; points are
not jittered or subsampled. Blue denotes small molecules and orange protein
chains. Exact epsilon agreement is shown on a zero-only residual axis.
[PDF](figures/parameter-parity.pdf) · [SVG](figures/parameter-parity.svg)

| Quantity | Unit | Paired values | Maximum absolute difference | RMSE |
| --- | --- | ---: | ---: | ---: |
| Partial charge | e | 29,252 | 2.229 × 10⁻⁷ | 2.691 × 10⁻⁸ |
| vdW sigma | nm | 29,252 | 5.551 × 10⁻¹⁷ | 1.261 × 10⁻¹⁷ |
| vdW epsilon | kJ/mol | 29,252 | 0 | 0 |
| Bond equilibrium length | nm | 29,712 | 0 | 0 |
| Bond force constant | kJ mol⁻¹ nm⁻² | 29,712 | 2.328 × 10⁻¹⁰ | 4.233 × 10⁻¹¹ |
| Angle equilibrium value | rad | 53,358 | 0 | 0 |
| Angle force constant | kJ mol⁻¹ rad⁻² | 53,358 | 0 | 0 |
| Proper-torsion coefficient | kJ/mol | 108,504 | 0 | 0 |
| Improper-torsion coefficient | kJ/mol | 18,624 | 0 | 0 |
| Constraint distance | nm | 14,408 | 0 | 0 |
| Atom feature value | feature-specific | 643,544 | 9.934 × 10⁻⁹ | 1.918 × 10⁻¹⁰ |

Torsion phases, periodicities, and division factors also agree exactly. Torsion
counts include all Fourier terms and improper multiplicities. Rule IDs, SMIRKS,
interaction coverage, and exception pairs are checked separately from numerical
values. The maximum deviation of summed atomic charge from formal molecular
charge is **5.596 × 10⁻¹⁴ e**.

Numerical parameters use absolute tolerance `1e-10` and relative tolerance
`1e-12`, through Python `math.isclose`; the relative term accommodates the tiny
rounding difference in the large bond force constants. Charges use absolute
tolerance **5e-5 e**, molecular charge sums `1e-6 e`, and features `1e-6`.
RMSE is computed over all paired values, so larger molecules contribute more
atoms. The two atom orders test invariance and are not independent molecules.

The separate model suite covers **132 cases**: 33 molecules in two atom orders
with each of Ash 1.0.0 and `openff-gnn-am1bcc-0.1.0-rc.2`. All pass. The latter
uses 21 features and no lookup table; Ash uses 22 features and 13,944 lookup
entries. It exercises additional lookup/library-charge paths. Ash's maximum final charge
difference is **2.9996 × 10⁻⁵ e**, caused by assigning slightly asymmetric stored
lookup charges to symmetry-equivalent atoms. The maximum forced-network
difference is **2.5332 × 10⁻⁷ e** for both models. Those results are kept separate from the
110-molecule statistics above; the original values and mappings are preserved
in the [model comparison](results/models.json.gz) and the historical archive.
Six independently generated OFFXML variants also test supported defaults, units,
method versions and optional handlers through the Rust file-loading API. These
are compatibility tests of the [supported subset](../../crates/kekule-openff/CONTRACT.md),
not validation of arbitrary force fields or arbitrary neural architectures.

## Energy parity

![Total-energy parity and signed residuals with identical charges and independently assigned charges.](figures/energy-parity.png)

**Figure 2.** Rust energy evaluation (recorded with the then benchmark-local
evaluator built on `potentials`) versus OpenMM, for all
660 geometry/atom-order comparisons. Left: both engines use the OpenFF charges,
isolating parameter assignment and energy evaluation. Right: Kekule uses its
own assigned charges, exposing their propagated electrostatic effect. The
energy axes use a symmetric logarithmic scale, linear between −100 and
+100 kJ/mol; the residual axes are linear and have different labeled scales.
All observations, including large energies from strained geometries, are shown.
[PDF](figures/energy-parity.pdf) · [SVG](figures/energy-parity.svg)

| Energy component | Maximum absolute difference (kJ/mol) | RMSE (kJ/mol) |
| --- | ---: | ---: |
| Bonds | 3.027 × 10⁻⁹ | 1.356 × 10⁻¹⁰ |
| Angles | 7.640 × 10⁻¹¹ | 6.068 × 10⁻¹² |
| Proper torsions | 4.184 × 10⁻¹¹ | 2.824 × 10⁻¹² |
| Improper torsions | 3.467 × 10⁻¹² | 1.993 × 10⁻¹³ |
| vdW | 3.024 × 10⁻¹⁰ | 2.355 × 10⁻¹¹ |
| Electrostatics, identical charges | 6.054 × 10⁻⁹ | 3.407 × 10⁻¹⁰ |
| Total, identical charges | 6.170 × 10⁻⁹ | 3.345 × 10⁻¹⁰ |
| Electrostatics, independently assigned charges | 6.011 × 10⁻⁴ | 9.790 × 10⁻⁵ |
| Total, independently assigned charges | 6.011 × 10⁻⁴ | 9.790 × 10⁻⁵ |

Each row contains 660 comparisons. The reference Hamiltonian is **vacuum,
nonperiodic NoCutoff**, with constrained bond terms retained in the energy.
Both engines use identical coordinates, Lorentz–Berthelot mixing, and matching
graph-based nonbonded exceptions. The three geometries are the prepared input
and seeded displacements of 0.002 and 0.01 nm standard deviation. These are
unminimized stress cases; large positive energies do not indicate stable conformers.

Identical-charge energies use absolute tolerance `1e-7 kJ/mol` and relative
tolerance `1e-10`. Independently assigned charges are evaluated separately
against a Coulomb error bound calculated from actual charge differences,
distances, and exception scales. This preserves the strict identical-charge
comparison while making the effect of charge roundoff explicit.

## Cartesian gradient validation

The benchmark now evaluates analytical **gradients, dE/dx**, in kJ mol⁻¹ nm⁻¹.
OpenMM returns forces; these are negated once at the reference boundary. Charges
stay fixed during differentiation: NAGL in this model depends on the molecular
graph, not coordinates. The recorded snapshot predates the public evaluator; see
[the public evaluator rerun](#public-evaluator-rerun-2026-10-08).

The same 110 molecules, three geometries and two atom orders give **660 geometry
comparisons**. Every Cartesian component is checked separately for bonds, angles,
proper torsions, impropers, vdW, electrostatics and the total: **1,842,876 scalar
comparisons per charge policy**. Identical-charge comparisons use absolute
tolerance `1e-5 kJ mol⁻¹ nm⁻¹` and relative tolerance `1e-10`. Independent-charge
electrostatic and total comparisons additionally use a per-atom, per-axis Coulomb
error bound derived from the actual charge-product differences, pair distances
and exception scales. The same bound handles charge roundoff between atom orders;
the identical-charge tolerance is unchanged.

| Gradient component | Maximum absolute difference, all 660 comparisons | Maximum in the other 658 comparisons |
| --- | ---: | ---: |
| Bonds | 2.750 × 10⁻⁹ | 2.750 × 10⁻⁹ |
| Angles | **2.021 × 10³** | 1.039 × 10⁻⁹ |
| Proper torsions | 9.313 × 10⁻¹⁰ | 2.440 × 10⁻¹⁰ |
| Improper torsions | 1.406 × 10⁻¹⁰ | 1.406 × 10⁻¹⁰ |
| vdW | 1.030 × 10⁻⁸ | 1.030 × 10⁻⁸ |
| Electrostatics | 4.241 × 10⁻¹¹ | 4.241 × 10⁻¹¹ |
| Total | **2.021 × 10³** | 1.019 × 10⁻⁸ |

All values are kJ mol⁻¹ nm⁻¹, with identical charges. The second numeric column
explicitly excludes **only the original geometry of PubChem 443915 in both atom
orders**; those observations remain in the first column and in the failures.
All ten proteins pass. The independently assigned-charge total gradient error
outside that geometry is at most `4.451e-4 kJ mol⁻¹ nm⁻¹`, within its calculated
charge-propagation bound.

Each of the 330 distinct geometries is also checked in eleven directions: nine
coordinate axes on the first, middle and last atoms, plus two seeded dense unit
directions involving every atom. Central energy differences at `1e-5` and
`5e-6 nm` are Richardson-extrapolated, separately for all seven components.
Near-linear geometries additionally use `1e-6/5e-7` and `1e-7/5e-8 nm` pairs;
the last two extrapolations must agree. **All 25,410 checks pass**, with maximum
absolute residual `4.254e-4 kJ mol⁻¹ nm⁻¹` against the analytical derivatives.
The finite-difference tolerance is `atol=2e-3, rtol=2e-6`; every step and energy
observation is retained. These are sampled directional derivatives, not exhaustive
finite differences of every protein coordinate.

Net gradient and net torque pass for every component and charge policy, with
maximum residuals `6.058e-10 kJ mol⁻¹ nm⁻¹` and `5.139e-10 kJ/mol` respectively.
Every geometry is also rigidly rotated and translated, and its gradients are
checked as vectors. Reversed atom orders are explicitly mapped back before
comparison. Constrained bond terms remain present; these are potential-energy
gradients, not constraint forces or gradients projected onto a constraint surface.

### Remaining near-linear case

PubChem 443915 has a nearly straight angle at mapped atoms 16–17–18 in its first
geometry: `sin(theta)=5.22349e-5`, with a bond-vector cross-product norm of
`8.76465e-7 nm²`. The assigned equilibrium angle is about 112 degrees, so this
is a strongly strained geometry with a large angular restoring force.

OpenMM's Reference angle implementation floors that cross-product norm at
`1e-6 nm²`. This reduces the returned angular force without making the same
change to the energy. Applying that documented cap predicts the observed OpenMM
angle gradient within `4.42e-8 kJ mol⁻¹ nm⁻¹`.
[Reference backend implementation](https://github.com/openmm/openmm/blob/8.4.0/platforms/reference/src/SimTKReference/ReferenceAngleBondIxn.cpp#L117-L130).
The linked source shows the convention; the stored observations were measured
with the pinned OpenMM **8.6.1** environment.

Numerically differentiating OpenMM's own energies independently confirms the
discrepancy: two dense directions differ from its reported angle gradients by
about **232 and 137 kJ mol⁻¹ nm⁻¹**. The same numerical estimates are within
0.05 of the uncapped native derivatives, with visible small-step roundoff from
OpenMM's angle energy. Native angle evaluation now uses stable `atan2` geometry
and its exact Cartesian derivative, avoiding both `acos` cancellation and the
different sine floor in `potentials::angle::Harm::derivative`.

The nearly singular torsions in this geometry also retain a strict rotation
failure: the native maximum residual is about `2.21e-4 kJ mol⁻¹ nm⁻¹`; OpenMM
shows a smaller residual of `6.66e-5` and passes the mixed tolerance. Native
torsion gradients at the original coordinates agree with OpenMM within `1e-9`.
This sensitivity is recorded, not relabeled as a passing comparison. Both
atom-order cases remain failures, and the full gradient command exits nonzero.

The public evaluator adopts an explicit near-singular policy: exact, uncapped
derivatives, rejecting only undefined quantities. This validation does not
establish periodic forces, minimizer behavior, constraint integration or MD stability.
The [independent reference](data/gradients.json.gz),
[full native results](results/gradients.json.gz), and
[initial failing observations](archive/gradient-before-fix.json.gz) retain the
evidence; the earlier energy references and tolerances were not changed.

```powershell
micromamba run -p target/openff-reference python benchmarks/openff/run.py gradients reference --output target/openff-gradients-reference-new.json.gz
micromamba run -p target/openff-reference python benchmarks/openff/run.py gradients compare --reference target/openff-gradients-reference-new.json.gz --output target/openff-gradients-native-new.json.gz
# Reproduce the diagnosed case alone:
micromamba run -p target/openff-reference python benchmarks/openff/run.py gradients compare --case pubchem-443915 --output target/openff-gradients-edge-new.json.gz
```

### Public evaluator rerun, 2026-10-08

The benchmark now calls the public `kekule_potentials::openff::OpenFfPotential`
instead of its own evaluator, so these comparisons cover shipped code. Both
commands were rerun against the unchanged frozen references
(`data/reference.json.gz` and `data/gradients.json.gz`) and unchanged tolerances;
the fingerprinted snapshot in `results/` and the figures above were not replaced.

- Energy: all 660 geometry comparisons pass for both charge policies, and every
  case keeps its recorded status: 214 passed, with the same six InChI-limit
  diagnostic failures. The identical-charge total maximum is unchanged at
  `6.170e-9 kJ/mol`. The angle maximum rose from `7.640e-11` to `4.280e-9 kJ/mol`,
  entirely in the first PubChem 443915 geometry: the recorded value used `acos`
  angles, and the public evaluator uses the `atan2` form adopted for gradients.
  Recomputing that molecule's angle energy both ways from the frozen reference
  gives exactly this `-4.2799e-9 kJ/mol` difference. All other component maxima
  changed by less than `3e-11 kJ/mol`.
- Gradients: all 220 cases keep their recorded status, 218 passed. The two
  failures are the same PubChem 443915 near-linear observations, with the same
  failed checks, components and counts, including the `2.021e3 kJ mol⁻¹ nm⁻¹`
  angle difference from OpenMM's capped derivative. Only their total-gradient
  rotation residuals, about `2.205e-4 kJ mol⁻¹ nm⁻¹`, moved, by at most
  `4e-12`: the total is now accumulated in one array instead of summed from
  separately evaluated components.

## CPU inference and parameterization time

![CPU execution times versus atom count, with median points and observed minimum-to-maximum error bars.](figures/cpu-timings.png)

**Figure 3.** Warm CPU execution times on an **Intel Core i5-11400F**, Windows
x86-64, Rust 1.94.1, release build, default CPU target and `ndarray` backend.
Points show medians; bars show the observed minimum and maximum, not confidence
intervals. No scaling fit is implied. The eight-atom molecule uses a charge
lookup in the full pipeline, making that path faster than forced inference.
[PDF](figures/cpu-timings.pdf) · [SVG](figures/cpu-timings.svg)

| Input | Atoms including H | Forced inference + preparation (ms) | Normal charge assignment (ms) | Full parameterization (ms) | Previous full (ms) |
| --- | ---: | ---: | ---: | ---: | ---: |
| pubchem-11642 | 8 | 2.15 | 0.088 | 0.74 | 0.63 |
| pubchem-111320 | 27 | 4.36 | 3.855 | 5.38 | 6.22 |
| pubchem-446684 | 38 | 4.44 | 4.774 | 7.05 | 7.59 |
| pubchem-476793 | 102 | 11.01 | 12.109 | 15.53 | 16.95 |
| pdb-6EE9-X | 378 | 40.00 | 39.518 | 55.72 | 56.32 |
| pdb-1W1F-A | 974 | 100.81 | 103.236 | 152.80 | 139.76 |
| pdb-3ABD-A | 2,940 | 324.83 | 328.676 | 462.13 | 447.83 |

The figure includes **all 110 inputs**; this table shows the seven inputs also
measured in the earlier run. Molecules were parsed before timing; model and force
field stayed loaded. Each operation received one warm-up, then five charge/feature
measurements or three full parameterization measurements. Forced inference includes
chemical preparation, domain checks, feature construction, the network, and charge
normalization. Full parameterization includes normal charge selection and force-field
assignment, but no energy calculation. The columns are independent measurements,
not timings to be added together. Model loading took **33.0 ms** in one separately
recorded observation and is excluded from the table.

Full parameterization changed by −13.5% to +16.5% across the seven shared inputs;
the 974-atom protein was 9.3% slower and the largest protein 3.2% slower.
Largest-protein forced inference remained about 325 ms. These are separate runs,
not a controlled before/after statistical experiment, so they do not establish
either a systematic regression or identical performance.

These measurements describe single-molecule latency on one machine. They do not
establish batch throughput, cross-platform performance, or GPU acceleration;
the current engine runs on CPU. The [timing results](results/timings.json) preserve
every repetition, input-panel hash, build information, executable fingerprint, and
source hashes. The reusable probe lives in
[`openff_timing.rs`](../src/bin/openff_timing.rs).

## Large protein chains, 2026-10-09

Removing the former 4,096-atom Ash limit was validated on four additional
protein chains with **4,615–13,818 atoms** including hydrogens:
`pdb-8D3R-C`, `pdb-7P9L-A`, `pdb-8ROT-B` and `pdb-3D6N-A`. They come from the
same locked PDB corpus and preparation as the main panel, selected in the same
ranked order with 300–1,000 residues and at least 4,097 atoms. The panel was
frozen before any native observation. Fresh OpenFF/OpenMM references and native
comparisons use the same script, tolerances and two atom orders as the main
panel. [Inputs](data/large-inputs.json.gz) ·
[reference](data/large-reference.json.gz) ·
[native results](results/large-native.json.gz)

| Validation measure | Result |
| --- | ---: |
| Parameterization cases passing all parameter, feature, charge and energy comparisons | **8 / 8** |
| Geometry comparisons | **24 / 24** |
| Paired atomic charges and vdW parameters | **63,010** |
| Maximum partial-charge difference | **1.860 × 10⁻⁷ e** |
| Maximum feature difference | **9.934 × 10⁻⁹** |
| Maximum total-energy difference, identical charges | **2.058 × 10⁻⁷ kJ/mol** |
| Maximum total-energy difference, independently assigned charges | **3.165 × 10⁻³ kJ/mol** |
| Strict cases, including standalone InChI diagnostics | **0 / 8** |

The identical-charge maximum is on `pdb-8ROT-B`, whose total energies are about
−4.3 × 10⁴ kJ/mol; it passes the relative tolerance of `1e-10`. Every case
passes the Coulomb error bound for independently assigned charges. Each
chain's native charges sum to its formal charge within 4 × 10⁻¹³ e.

All eight strict failures are the standalone fixed-H InChI diagnostic, which
exceeds the InChI adapter's 1,024-atom limit, as for the three largest
main-panel proteins. Charge assignment does not compute this identifier for
molecules larger than Ash's 11-atom largest lookup entry.

The OpenMM reference's split-versus-unsplit self-consistency check now uses
the comparison tolerance `ENERGY_RTOL` (`1e-10`) rather than `1e-12`. On
`pdb-8ROT-B` the two OpenMM evaluations differed by 7.05 × 10⁻⁸ kJ/mol, a
relative 1.65 × 10⁻¹², which is OpenMM summation roundoff at that size.

One single-run native observation took 2.9–10.2 s per chain for full
parameterization, including JSON transfer, on the machine above. This is not a
controlled timing measurement.

## Dataset, reference versions, and coverage

Small molecules were selected deterministically from the Kekule `pubchem-100k`
benchmark, using chemical diversity descriptors without consulting Kekule
success. They contain 3–40 heavy atoms, or 8–102 atoms with hydrogens.
Protein chains were sourced from `pdb-1000`: ten independently preparable chains
from distinct entries, with 25–180 residues and 378–2,940 atoms including H.
Hydrogen preparation used OpenMM at pH 7; OpenFF supplied bond orders and formal
charges. Waters, other chains, and ligands were excluded. There were 66 recorded
read/preparation failures before selection completed; these remain archived and
are outside the 110 prepared-input denominator. Full source hashes, preparation
details and all preparation attempts are in [the frozen inputs](data/inputs.json.gz).
Protein identities, residue counts and formal charges are also listed in
[the machine-readable summary](results/summary.json).

| Reference component | Pinned version |
| --- | --- |
| Rosemary OFFXML | `openff_no_water-3.0.0-alpha2b` |
| Ash model | `openff-gnn-am1bcc-1.0.0.pt`, checksum-pinned |
| OpenFF Toolkit / NAGL / Interchange | 0.19.0 / 0.6.1 / 0.5.5 |
| RDKit / PyTorch / NumPy | 2026.03.3 / 2.10.0 / 2.5.3 |
| OpenMM reference / Rust `potentials` consumer | 8.6.1 / 0.1.0 |

| SMIRNOFF handler | Rules exercised in the 110-molecule panel | Rules in Rosemary |
| --- | ---: | ---: |
| Bonds | 66 | 90 |
| Angles | 35 | 42 |
| Proper torsions | 116 | 187 |
| Improper torsions | 7 | 7 |
| vdW | 23 | 35 |
| Constraints | 1 | 1 |

The panel has 218 inference and two lookup observations. Coverage is substantial
but incomplete: additional rule coverage, all lookup mappings, periodic PME and
switching, constraint enforcement, solvated complexes, and MD
stability require separate validation. Protein graph compatibility does not
establish Rosemary's empirical suitability for protein simulations. Scientific
parity here was measured on Windows; Linux/macOS numerical parity remains to be
established. Ash inference has no molecule size limit; see the
[large protein chains](#large-protein-chains-2026-10-09).

The six strict diagnostic failures are the original and reversed forms of
`pdb-3ABD-A`, `pdb-9B3P-B` and `pdb-8J90-E`: each exceeds the InChI adapter's
1,024-atom limit. Ash's largest lookup entry has 11 atoms, so charge assignment
bypasses lookup for these proteins and their parameterizations pass. The energy
comparison command deliberately retains a nonzero exit status for these diagnostic
failures; they are not removed from the report.

For small eligible molecules, a remaining compatibility difference is that
native InChI generation errors propagate, whereas upstream treats an empty
InChI as a lookup miss and attempts inference. The supported runtime boundaries
are documented in the [crate contract](../../crates/kekule-openff/CONTRACT.md).


## Reproduce and inspect

Run from the repository root. Scientific tools are optional benchmark dependencies,
never Rust runtime dependencies. The environment includes the pinned reference
packages and Matplotlib used for these figures. Commands write fresh results under
`target/` and refuse to replace existing raw outputs; choose unused output names
for another run.

```powershell
micromamba create -y -p target/openff-reference -f benchmarks/openff/environment.yml
cargo build -p kekule-bench --release --locked --bin openff_parameterize --bin openff_timing
micromamba run -p target/openff-reference python benchmarks/openff/run.py export target/openff-models/openff-gnn-am1bcc-1.0.0

# Generate fresh independent OpenFF/OpenMM observations and freeze the JSONL stream.
micromamba run -p target/openff-reference python benchmarks/openff/run.py energy reference --output target/openff-reference-new.jsonl
micromamba run -p target/openff-reference python benchmarks/openff/run.py energy freeze --reference target/openff-reference-new.jsonl --output target/openff-reference-new.json.gz
micromamba run -p target/openff-reference python benchmarks/openff/run.py energy compare --reference target/openff-reference-new.json.gz --output target/openff-native-new.json.gz

# Large protein chains: the same comparison on the frozen large-protein panel.
micromamba run -p target/openff-reference python benchmarks/openff/scripts/robustness.py reference --inputs benchmarks/openff/data/large-inputs.json.gz --output target/large-reference-new.jsonl
micromamba run -p target/openff-reference python benchmarks/openff/scripts/robustness.py freeze --inputs benchmarks/openff/data/large-inputs.json.gz --reference target/large-reference-new.jsonl --output target/large-reference-new.json.gz
micromamba run -p target/openff-reference python benchmarks/openff/scripts/robustness.py compare --inputs benchmarks/openff/data/large-inputs.json.gz --reference target/large-reference-new.json.gz --output target/large-native-new.json.gz

# Full 110-input performance panel; supply the actual CPU description.
python benchmarks/openff/run.py timings --cpu "Intel Core i5-11400F" --output target/openff-timings-new.json

# Rebuild figures and summary from the checked-in, fingerprinted validation snapshot.
micromamba run -p target/openff-reference python benchmarks/openff/run.py plot
python -m unittest discover -s benchmarks/openff/scripts -p 'test_*.py'
```

The default binaries have Windows `.exe` paths; on other platforms pass `--binary`
with the corresponding executable. Model identity and tensor metadata are checked
when loading the bundle. The separate two-model end-to-end workflow is:

```powershell
micromamba run -p target/openff-reference python benchmarks/openff/run.py models-reference --bundles target/openff-models-new --output target/openff-models-reference-new.json.gz
python benchmarks/openff/run.py models-compare --binary target/release/openff_parameterize.exe --bundles target/openff-models-new --reference target/openff-models-reference-new.json.gz --output target/openff-models-native-new.json.gz
$env:KEKULE_OPENFF_MODELS = (Resolve-Path target/openff-models-new).Path
cargo test -p kekule-openff --test models -- --ignored
```

Use `python benchmarks/openff/run.py --help` for the command list and append
`--help` to any command for its options. `prepare` rebuilds the panel from locally
available PubChem/PDB corpora; `offxml-reference` regenerates the six parser
reference cases; `audit` retains the prerequisite perception/SMARTS investigation.
Frozen runtime reference fixtures stay in `crates/kekule-openff/tests/fixtures`.

| Location | Contents |
| --- | --- |
| `data/` | Frozen 110-input panel, fresh OpenFF/OpenMM reference and model provenance |
| `results/` | Native observations, all timing repetitions, summary and fingerprints |
| `figures/` | Three Matplotlib figures in PNG, PDF and SVG |
| `fixtures/` | External source extracts, licenses and source locks |
| `scripts/` | Active reference, comparison, export and offline verification code |
| `archive/previous-validation.zip` | Historical observations, reports and retired scripts |

The archive preserves 57 previous files byte for byte, including pre-fix failures,
the exhaustive lookup-identity study and original timing probe. Its `index.json`
records every entry's SHA-256; offline tests verify the archive and its original
report locks. Retired standalone scripts and overlapping audit/implementation
documents were consolidated here, rather than leaving multiple apparent current
reports. No comparison tolerances or stored failure observations were removed.

## Repository checks for this rerun

The gradient addition reran formatting, all-feature workspace check/Clippy/tests,
the Rust 1.89 check, documentation with warnings denied, Python checks, and the
benchmark package file listing. It passed **1,416 Rust tests and 29 Python tests**.
The three optional exported-model tests were not repeated for this benchmark-only
addition; they passed in the preceding run recorded below. Standalone OpenFF
package verification was not repeated; its existing registry-dependency failure
is also retained below. The scientific gradient run itself exits nonzero for
the two documented near-linear cases.

| Check | Result |
| --- | --- |
| `cargo test --workspace --locked` | 1,415 passed; three external-model tests run separately |
| `cargo test -p kekule-openff --locked -- --ignored` with both exported bundles | Three passed |
| Python offline provenance/comparison tests | 26 passed |
| All nine original `run.py` command help paths and prerequisite source inventory | Passed |
| `cargo fmt --all -- --check` | Passed |
| `cargo check --workspace --all-targets --locked` | Passed |
| `cargo +1.89.0 check --workspace --all-targets --locked` | Passed |
| `cargo clippy --workspace --all-targets --locked -- -D warnings` | Passed |
| `cargo doc --workspace --no-deps --locked`, `RUSTDOCFLAGS=-D warnings` | Passed |
| `cargo package --allow-dirty --list --locked` for `kekule-openff` and `kekule-bench` | Passed |
| `cargo package -p kekule-openff --allow-dirty --locked --offline` | Tarball created; dependency verification failed |

The standalone package resolves the registry's `kekule 0.2.1`, which lacks the
workspace's MDL aromaticity and prepared/tagged SMARTS APIs, among other additions.
A release must publish the updated core under a new version and update the
companion dependency before the OpenFF package can verify against the registry.
This limitation does not affect the workspace benchmark results above.
