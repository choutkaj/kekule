# OpenFF prerequisite audit

This document preserves the initial audit. The subsequent implementation and
its additional reference checks are described in [IMPLEMENTATION.md](IMPLEMENTATION.md)
and the [crate contract](../../crates/kekule-openff/CONTRACT.md).

This is a development audit for the proposed `kekule-openff` crate, not a force
field implementation or a claim of complete SMIRNOFF/NAGL compatibility.
The initial audit targets Rosemary `openff_no_water-3.0.0-alpha2b.offxml`
(2026-09-11) and Ash `openff-gnn-am1bcc-1.0.0.pt`. It changes no core APIs.

## Reproducible inputs and executable cases

`sources.lock.json` pins the full OFFXML, upstream NAGL source/test files, model
card, and licenses by immutable commit and SHA-256. Files under `fixtures/` are
unmodified external source bytes, not code imported by the audit. The model is
not vendored: its SHA-256 must equal the value specified by the OFFXML.
The copied force-field/model materials have their upstream licenses; NAGL source
and tests have the included MIT license. Python/OpenFF/RDKit remain reference-only.

The 23 molecules comprise 11 literal selections from the upstream NAGL lookup,
resonance and shared fixtures, all eight existing PubChem SMILES smoke inputs,
and all four existing PubChem SMARTS auxiliary inputs. Every selected literal
is checked against its source AST. Reused PubChem bytes are verified against
their existing source locks. Selection is intentional and small; it does not
establish whole-domain coverage. No benchmark molecules were synthesized.

The audit has three independently executable stages:

```text
python benchmarks/openff/audit.py inventory --output target/openff-inventory.json
cargo build -p kekule-bench --bin openff_prerequisites --bin smarts_conformance --locked
python benchmarks/openff/audit.py primitives --graph-binary target/debug/openff_prerequisites --smarts-binary target/debug/smarts_conformance --output target/openff-primitives.json.gz
python benchmarks/openff/audit.py reference --graph-binary target/debug/openff_prerequisites --smarts-binary target/debug/smarts_conformance --output target/openff-reference.json.gz
```

On Windows use the `.exe` executable suffixes. Inventory needs only Python's
standard library. Primitives needs RDKit 2026.03.3. The reference stage uses the
exact top-level versions in `environment.yml`; create a separate environment:

```text
micromamba create --override-channels -c conda-forge -p target/openff-reference -f benchmarks/openff/environment.yml
micromamba run -p target/openff-reference python benchmarks/openff/audit.py reference --graph-binary target/debug/openff_prerequisites.exe --smarts-binary target/debug/smarts_conformance.exe --output target/openff-reference.json.gz
```

`--model PATH` selects an already downloaded checkpoint. Without it the OpenFF
model resolver may download the published model. Package versions and relevant
runtime source hashes are recorded separately from the source-review commit
pins. The CPU reference run should not be interpreted as a speed benchmark.

Offline maintenance tests do not load any scientific tools:

```text
python -m unittest discover -s benchmarks/openff -p test_audit.py -v
cargo test -p kekule-bench --bin openff_prerequisites --locked
```

## Observations and comparison contract

The native graph observer accepts one fully explicit, uniquely mapped connected
molecule per JSONL request. It reports element, formal charge, graph degree,
selected-ring membership at sizes 3/4/5/6, and MDL atom aromaticity, keyed by atom
map rather than traversal index. Missing/duplicate labels, implicit hydrogens,
multiple components, and fallback ring bases are errors. It uses public kekule
APIs and installs perception only on its private input molecule.

The SMARTS stage crosses every one of the **371** Rosemary patterns with every
case. It uses the existing complete-mapping adapter, retaining full mappings,
tag permutations and duplicate ordered tag tuples, query sizes, and all target
aromatic atom flags. Both engines receive the same hydrogen-suppressed unmapped
serialization, recorded beside the original external input. This explicit
representation conversion avoids confusing RDKit's input-H removal/append order
with Kekule's graph-H preservation. It is not a test of atom-order preservation
by the input parser; the separate mapped graph observer checks correspondence.

The full reference stage additionally records fixed-H and standard InChI,
normalized molecular representation, resonance-averaged formal charges, every
configured NAGL input feature, domain decisions, lookup hit/miss and charges,
raw model and Toolkit-normalized charges, parameter IDs and tagged atom tuples,
and charge invariance under reversed atom order (absolute tolerance 5e-5 e,
raised from the original diagnostic threshold of 2e-6 e at the user's request).
It executes upstream lookup assertions for reordered atoms, alternate represented
charges/bonds, and a miss. Those small regression charges are upstream test data,
not substituted production-model charges.

Errors remain attached to cases. An equal pair of error objects is never an
agreement. Missing native output fails the stream. Nonfinite JSON is rejected.
Reports are replaced atomically after each full case, with `complete: false`
until the selected run finishes. A completed audit can still contain explicit
unsupported-domain charge results; inspect coverage, not just its completion flag.

## Measured primitive results (2026-09-22)

`primitives.json.gz` records **8,533/8,533 exact SMARTS comparisons**, and
**23/23 exact mapped atom-feature comparisons**, with no errors or exclusions.
There are 425 nonempty query/target pairs: Constraints 21, Bonds 93, Angles 84,
ProperTorsions 72, ImproperTorsions 18, vdW 135, LibraryCharges 2. Empty matches
are retained. This is not positive coverage of all 371 parameters.

The primitive report has zero charge cases and a null lookup-contract outcome;
it must not be used as evidence of NAGL inference or lookup compatibility.

## Measured full reference results

The pinned reference is Toolkit 0.19.0, NAGL 0.6.1, NAGL Models 2026.09.0,
Interchange 0.5.5, RDKit 2026.03.3, PyTorch 2.10.0 and NumPy 2.5.3.
`reference.json.gz` retains the complete observations, not just summary counts.
`reports.lock.json` binds both compressed report files and their generator
fingerprints. The offline tests verify these archived bytes and keep the
reference disagreements visible. New experiments should write under `target/`;
do not replace these observations just to obtain a passing result.

- All 23 Rosemary Interchange systems are constructed successfully. This is
  reference-only evidence; no native parameterizer exists yet.
- The four lookup-contract assertions pass, including the same-InChI / different
  represented-chemistry assertion and the explicit upstream charge remapping.
- 22 cases produce NAGL charges. Sodium is rejected by the NAGL domain and feature
  encoder but receives +1 e from Rosemary's LibraryCharges in the complete system.
  The rejection is retained; it is not counted as successful NAGL inference.
- Nine cases hit the production lookup table. The upstream small-test-table miss
  for ethane is intentionally distinct from a hit in the production table.
- The checksum-pinned checkpoint contains **13,944** lookup entries. Its upstream
  model card says 13,939; the executable inventory retains the actual count.
- Reversing atom order exposes **three charge-permutation disagreements**: the
  two hydrogen-sulfide representations and ethane. The sulfide discrepancy is
  approximately 4.992e-6 e; the ethane discrepancy is approximately 2.0005e-5 e.
  The production lookup charges are slightly asymmetric between equivalent
  atoms, and graph remapping can exchange them. The 2e-6 e threshold is unchanged;
  raw arrays and failed assertions are preserved. No symmetrization is applied.

The archived report records the original **exit code 1** under its 2e-6 e
diagnostic threshold. The current user-approved cutoff is **5e-5 e**; all three
observations pass under this policy. The original report bytes, failed flags,
raw arrays, and measured differences remain unchanged. Current classification
uses the recorded numerical error and the current cutoff. Equivalent-atom tie
breaking is not an implementation blocker; lookup charges are not averaged.

There is no positive 3-/4-membered-ring coverage in this small panel. Additional
external fused/bridged/small-ring cases, production-table collision cases,
normalization/reaction coverage, and larger peptide/protein performance cases
remain follow-up work. The full NAGL feature vectors are reference observations;
only the explicitly listed graph/MDL features are compared to native Rust.

## Prerequisite decisions

| Area | Finding | Next action / owner |
| --- | --- | --- |
| Explicit H and identity | Existing molecular conversion and qualified topology identities are suitable | Keep atom addition explicit; bind parameterization to a final snapshot |
| MDL and SMARTS | This Rosemary panel agrees; full tagged enumeration already exists | Reuse in `kekule-openff`; expand positive environments before claiming coverage |
| Ring-size features | Can be extracted from the public selected-ring set without a core API change | Expand fused/bridged/small-ring feature cases; do not substitute smallest-ring SMARTS |
| Owning perception options | Topology perception currently installs its fixed default profile | Optional core enhancement for named profiles; temporary per-definition copies also work |
| Normalization | Core source normalization is not NAGL's chemical normalization policy | Port the pinned NAGL transformations on temporary chemistry; do not silently rewrite the input |
| Resonance average | No equivalent public NAGL/Gilson enumerator was found | Implement bounded donor/acceptor fragment enumeration and exact selection rules; report exhaustion |
| Lookup identity | Fixed-H InChI equivalence is broader than represented graph equality | Resolve native-vs-FFI identity implementation before claiming exact lookup parity |
| GNN | No runtime inference implementation exists | Export pinned model weights/configuration/table, port forward pass and postprocessing |
| Assignment | No OFFXML/SMIRNOFF assignment layer exists | New crate: ordered rules, handler-specific symmetry, precedence, completeness checks |
| Properties | Typed/unit-aware scalar atom/bond columns already suffice for projections | Keep complete parameterization in a separate typed result, export selected annotations |
| Evaluation | Existing potential contract binds an exact topology | Later adapter in `kekule-potentials`; this audit does not implement energy/forces |

The force field contains 90 bond, 42 angle, 187 proper, 7 improper, 35 vdW,
1 constraint and 9 library-charge patterns. Its Bonds/ProperTorsions sections
declare AM1-Wiberg defaults but contain no interpolated bond-order parameters.
Do not infer a need for an AM1 runtime from those defaults alone. It has no water
model or virtual sites. Constraint distance inference, torsion `idivf`/trefoils,
charge precedence, units, exclusions and nonbonded settings still need assignment
tests; successful SMARTS comparisons do not validate these behaviors.

## InChI decision

Do not begin a complete native InChI implementation as an incidental prerequisite
of this audit. First freeze the exact identity and remapping contract.
Upstream uses the **full fixed-hydrogen InChI string**, not an InChIKey, then
strict isomorphism followed by relaxed bond-order/formal-charge matching and,
if necessary, relaxed stereo. The upstream nitromethane case supplies two
represented graphs with the same InChI but different charges and bond orders.
Ordinary kekule equality or canonical SMILES cannot replace that key faithfully.

General InChI support would be useful independently. The smallest established
implementation route is an optional adapter to the official IUPAC InChI C
library, isolated from core chemistry. That would introduce a native C runtime
dependency and requires an explicit architectural choice: it does not fulfill
a strict all-Rust implementation requirement. A pure Rust InChI implementation
is a separate substantial project with its own conformance corpus.

Keep lookup identity behind a narrow boundary in the planned charge engine.
Conversion of the checkpoint's finite lookup table to native graph keys is a
candidate experiment, not an established equivalent replacement: checking every
stored entry does not prove equivalent recognition of alternate incoming forms.
Do not silently disable the table or return a lookup miss because InChI is absent.
The full lookup contract and its source are supplied here to inform that choice.

## Maintenance validation

The audit was developed on branch `codex/openff-prerequisite-audit`, based on
`e7883a80`. The working tree was initially clean. No core/runtime crate, README,
or dependency manifest was changed.

Passed:

- `cargo fmt --all -- --check` and `git diff --check`.
- `cargo check --workspace --all-targets --all-features --locked`.
- `cargo clippy --workspace --all-targets --all-features --locked -- -D warnings`.
- `cargo test --workspace --all-features --locked`: 1,385 passed, including
  doctests and the two new observer tests; `RUST_TEST_THREADS=1` on Windows.
- `RUSTDOCFLAGS="-D warnings" cargo doc --workspace --all-features --no-deps --locked`.
- `cargo +1.89.0 check -p kekule-bench --bin openff_prerequisites --locked`.
- `cargo package -p kekule-bench --allow-dirty --locked --list` (the benchmark
  crate is unpublished; nothing was published).
- Twelve offline audit/provenance tests and the existing four SMARTS contract tests.
- Inventory and primitive reference commands above. Full reference execution
  completed with the three preserved permutation disagreements described above.

Not repeated: the separate workspace doctest command (already included in
`cargo test`), runtime crate package builds and companion package lists (no
runtime/package contents changed), no-default-feature potential checks, fuzz
checks/runs, and unrelated external benchmarks (no affected implementations).
The MSRV check targets the new observer; the full workspace MSRV matrix was not
rerun because other Rust sources and dependency versions are unchanged.
