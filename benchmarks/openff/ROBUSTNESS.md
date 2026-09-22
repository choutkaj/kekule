# Rosemary robustness validation, 2026-09-22

All **110 molecules parameterize successfully** after this pass: 100 PubChem
small molecules and ten PDB protein chains. Both original and reversed atom
orders pass parameter, feature, charge, exception, and energy comparisons:
**220/220 parameterization cases, 660 geometry evaluations, 29,252 paired atoms**.

The stricter result including the standalone InChI diagnostic is **214/220**.
The remaining six failures are the two atom orders of three proteins above the
InChI Rust adapter's 1,024-atom limit. They remain failures in the archive and
the comparison command returns nonzero. Their parameterization and energies
pass; charge assignment no longer depends on generating those identifiers.

## Frozen sources and preparation

`robustness-inputs.json.gz` is the self-contained input panel. It includes source
corpus/pack hashes, source SMILES or PDB entry/chain, mapped explicit-H SMILES,
coordinates, preparation versions, protein residue sequences and hydrogen
variants, and the prepared protein PDB text. `robustness-reports.lock.json`
fingerprints the inputs, independent reference, baseline failures, final native
observations, and summary. The bulk corpora are needed only to rebuild selection.

Small molecules come from the existing `pubchem-100k` source lock. All SMILES
packs are verified before selecting a hash-ranked pool of 4,000 entries. The
eligible pool contains connected, nonradical molecules with 3–40 heavy atoms
and Ash-supported elements. A deterministic greedy selection covers element,
corpus category, size, charge, ring size/count, and atom-environment descriptors;
ties retain hash order. Selection does not use native or reference success.
RDKit ETKDGv3 embeds explicit hydrogens with seed 20260922. All 100 selected
inputs prepared successfully; they contain 8–102 atoms including hydrogens.
Undefined source stereochemistry remains undefined for parameter assignment.

Protein entries come from `pdb-1000`, in deterministic hash order. Eligible chains
have 25–180 canonical amino-acid residues. The first ten independently preparable
chains from distinct entries are retained. OpenMM reads the first coordinate
model, existing hydrogens are removed, and `Modeller.addHydrogens` uses pH 7,
`amber14-all.xml`, seed 20260922, and the Reference platform. Amber is used only
for hydrogen preparation. OpenFF `Topology.from_pdb` assigns the prepared graph's
bond orders, formal charges, and stereochemistry. No missing heavy atoms are
invented, chains with inter-chain covalent bonds are rejected, and disconnected
chains are rejected. Waters, other chains, and ligands are excluded.

There are 76 recorded read/preparation attempts: ten selections and 66 failures,
including missing atoms, incomplete termini, and inter-chain bonds. These
preparation failures are retained separately from the 110-input validation
denominator. Selection finishes before any native parameterization runs.

| PDB chain | Residues | Atoms including H | Formal charge |
| --- | ---: | ---: | ---: |
| 1PJV A | 32 | 507 | +6 |
| 2JUC A | 55 | 937 | −6 |
| 1IWC A | 34 | 504 | +2 |
| 3ABD A | 180 | 2,940 | −7 |
| 2N4H A | 50 | 700 | +1 |
| 9B3P B | 78 | 1,281 | +10 |
| 2F52 A | 67 | 1,014 | −6 |
| 8J90 E | 77 | 1,274 | +3 |
| 6EE9 X | 25 | 378 | 0 |
| 1W1F A | 60 | 974 | 0 |

## Independent reference and energy protocol

The pinned `environment.yml` uses Toolkit 0.19.0, NAGL 0.6.1, Interchange 0.5.5,
RDKit 2026.03.3, OpenMM 8.6.1, PyTorch 2.10.0, and NumPy 2.5.3. The same
fingerprint-pinned Rosemary OFFXML and Ash checkpoint used in the initial audit
are retained. Reference generation never invokes Kekule.

Each reference stores exact rule IDs/SMIRKS, every numerical parameter and
torsion multiplicity, all 22 feature columns, final system charges, nonbonded
exceptions, mapped coordinates, energy components, and the exported OpenMM
system XML. No reference values were changed in response to native failures.

The reference Hamiltonian is explicitly vacuum, nonperiodic **NoCutoff**, with
all valence forces retained, including constrained bonds. Interchange exports
combined nonbonded forces; the observer splits the exported forces into six
energy groups without changing their parameters. The split total is checked
against the original exported system. This avoids the rounded Coulomb constant
used in Interchange's separate 1–4 custom-force export path.

The Rust benchmark consumer uses [`potentials` 0.1.0](https://docs.rs/potentials/0.1.0/potentials/):
harmonic bonds/angles receive `k/2`, periodic torsions receive `k/idivf`, sigma
mixes arithmetically, epsilon geometrically, and explicit graph exceptions apply
the emitted scales. The Coulomb prefactor is OpenMM's
`138.93545764438198 kJ mol^-1 nm e^-2`. `potentials` is a benchmark dependency;
the parameterization crate retains its existing ownership boundary.

Three geometries are evaluated: the prepared coordinates, then seeded Gaussian
displacements of standard deviation 0.002 and 0.01 nm. Coordinates are identical
in both engines and keyed by atom map, including reversed native atom order.
These unminimized geometries exercise strained terms and can contain clashes;
their energies are numerical validation data, not predictions of stability.

Two energy observations are retained: native parameters with reference charges
isolate assignment/evaluation error, while native parameters with native charges
measure the full result. The latter is checked against an explicit Coulomb error
bound derived from the actual per-atom charge differences and distances. This
does not relax the charge tolerance or the energy comparison with shared charges.

## Numerical results

| Quantity | Maximum absolute difference |
| --- | ---: |
| Partial charge | 2.229e−7 e |
| Molecular charge conservation | 5.596e−14 e |
| Feature column | 9.934e−9 |
| Bond energy, shared charges | 3.027e−9 kJ/mol |
| Angle energy, shared charges | 7.640e−11 kJ/mol |
| Proper-torsion energy, shared charges | 4.184e−11 kJ/mol |
| Improper-torsion energy, shared charges | 3.467e−12 kJ/mol |
| vdW energy, shared charges | 3.025e−10 kJ/mol |
| Electrostatic energy, shared charges | 6.054e−9 kJ/mol |
| Total energy, shared charges | 6.171e−9 kJ/mol |
| Total energy, native charges | 6.012e−4 kJ/mol |

Exact label/exception coverage and multiplicities are asserted. Numerical
parameters use `atol=1e-10, rtol=1e-12`; charges retain the approved `5e-5 e`
cutoff; features use `1e-6`; shared-charge energies use `atol=1e-7 kJ/mol,
rtol=1e-10`. The largest observed energy residual is far below that threshold.
There are 218 inference and two lookup observations; the original 66-case suite
remains necessary for broader lookup/library coverage and was rerun successfully.

| Handler | Rules exercised | Rosemary rules |
| --- | ---: | ---: |
| Bonds | 66 | 90 |
| Angles | 35 | 42 |
| Proper torsions | 116 | 187 |
| Improper torsions | 7 | 7 |
| vdW | 23 | 35 |
| Constraints | 1 | 1 |

## Hardening changes

* Connected SMARTS extensions use the adjacency of an already matched neighbor,
  preserving candidate order and complete-match semantics. This fixes quadratic
  scans that exhausted the existing ten-million-state limit on 3ABD A. No search
  limit was raised. A 1,500-atom focused regression asserts all expected matches
  under a 30,000-state budget; existing recursive, disconnected, stereo, and
  exhaustion tests remain in place.
* Ash lookup uses an isotope-free private copy to match OpenFF Toolkit's charge
  representation. PubChem CID 143783 exposed the discrepancy. Kekule's input
  isotopes and the general internal InChI adapter remain unchanged. A regression
  checks the external case and input immutability; the model test also checks
  isotope/non-isotope lookup charge equality.
* Inputs larger than the frozen table's largest entry (11 atoms) proceed directly
  to inference. A full-molecule lookup match cannot have a different charge-array
  length. This removes an unnecessary dependency on InChI's 1,024-atom limit.
  The public identifier diagnostic retains that explicit limitation and its six
  failures remain visible in this panel.
* Manifest and weights are checked for their exact expected sizes before
  allocation, and reads are capped against file growth. Truncated, oversized,
  missing, and nonfinite inputs have focused tests. The 4,096-atom feature limit
  is checked before expensive preparation/inference work.
* The optional model regression covers shared immutable engines on four threads,
  serial/concurrent charge equality, input immutability, the large-molecule
  lookup bypass, and early rejection above the inference size limit.

The baseline archive contains eight failures: two isotope-key differences, two
SMARTS exhaustion failures, and four oversized identifier failures. The final
archive retains six identifier failures while completing all parameter/energy
checks, including the largest protein that previously failed before charging.

## Reproduction

Use the environment and model-export commands in [IMPLEMENTATION.md](IMPLEMENTATION.md).
Outputs must be new paths; reference journals resume only the same input hash.

```powershell
# Optional: requires locally supplied full benchmark corpora.
micromamba run -p target/openff-reference python benchmarks/openff/prepare_robustness.py --output target/new-inputs.json.gz

# Generate an independent reference, then freeze its append-only journal.
micromamba run -p target/openff-reference python benchmarks/openff/robustness.py reference --output target/new-reference.jsonl
python benchmarks/openff/robustness.py freeze --reference target/new-reference.jsonl --output target/new-reference.json.gz

# Replay the archived reference. Exit 1 includes the documented InChI diagnostics.
cargo build -p kekule-bench --bin openff_parameterize --release --locked
micromamba run -p target/openff-reference python benchmarks/openff/robustness.py compare --output target/new-native.json.gz
python benchmarks/openff/summarize_robustness.py --output target/new-summary.json
python -m unittest discover -s benchmarks/openff -p test_robustness.py
```

Pass `--inputs` and `--reference` to use newly generated artifacts; defaults
replay the checked-in panel. A `.jsonl` progress journal accompanies each native
run, and a per-case timeout preserves explicit failures rather than hanging the
entire panel. The original 66-case recheck is also archived.

## Remaining scope

This establishes substantially broader implementation parity, including ten
protein-sized graphs. Protein-chain compatibility does not establish empirical
protein force-field accuracy. Periodic PME, cutoff/switching behavior, dynamics
and constraint enforcement, analytic force parity, full complexes/solvation,
all lookup permutations, the remaining parameter rules, and cross-platform
floating-point behavior still need separate validation. Large standalone InChI
identifiers require an upstream adapter improvement. The runtime remains bounded
to 4,096 atoms for Ash inference.

## Maintenance checks

Passed on Windows with the final changes:

```text
cargo fmt --all -- --check
cargo check --workspace --all-targets --all-features --locked
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo test --workspace --all-features --locked
cargo test -p kekule-openff --test contracts --locked -- --ignored
cargo +1.89.0 check --workspace --all-targets --all-features --locked
cargo check --manifest-path fuzz/Cargo.toml --bins --locked
cargo doc --workspace --all-features --no-deps --locked
cargo package -p kekule --allow-dirty --locked
cargo package -p kekule-openff --allow-dirty --locked --list
cargo package -p kekule-bench --allow-dirty --locked --list
python -m unittest discover -s benchmarks/openff -p "test_*.py"
git diff --check
```

Workspace tests report **1,406 passed**, zero failed, one model-dependent test
ignored. That model test was separately run successfully with
`KEKULE_OPENFF_MODEL` set to the absolute exported-bundle path, including its
large-input, isotope, and concurrency regressions. The 18 standard-library
Python contract/provenance tests pass. Rustdoc used `RUSTDOCFLAGS=-D warnings`.
The original 66-case live OpenFF comparison also passes in full.

Not run: a separate `cargo test --workspace --all-features --doc --locked`
command, since workspace tests already executed the doctests; full companion
package verification/publication, because the workspace's core 0.2.1 is not
published and CI therefore checks companion file lists; unchanged companion
package checks and the unchanged `kekule-potentials --no-default-features`
checks. Nightly fuzz execution and Linux/macOS runs remain CI/platform work;
all registered fuzz targets were checked for compilation locally. Other optional
external corpora were not rerun. The new scientific panel and original 66 cases
are the external comparisons relevant to this pass.

External reference runs remain optional scientific benchmarks, not routine
CI/release gates. CI runs the offline harness/provenance tests.
