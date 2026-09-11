# Valence and rings audit — 2026-09-11

The reference is **RDKit 2026.03.3**, pinned by the existing benchmark
environment. The audit covers represented bond valence, hydrogen inference,
strict and permissive validation, cycle membership, ring selection and
symmetrization, resource limits, and their integration with aromaticity.

Sources were checked against pinned [Atom.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/Atom.cpp),
[atomic_data.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/atomic_data.cpp),
[FindRings.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/FindRings.cpp),
and [MolOps.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/MolOps.cpp),
with executable confirmation from that same RDKit release.

## Valence findings

The previous calculation combined explicit-valence validation with hydrogen
inference. That conflated rules with different applicability and led to several
edge-case errors. The corrected calculation follows separate stages:

1. Read represented bond contributions and declared hydrogens.
2. Select original and charge-adjusted element rules, retaining unrestricted
   entries and clamping the effective atomic number as the reference does.
3. Validate explicit valence in strict mode, including the distinct hydride and
   hypervalent-anion rules.
4. Return zero immediately when hydrogen inference is disabled.
5. Handle isolated hydrogen charges and unrestricted target elements.
6. Infer hydrogens from occupied valence, radical electrons and any hypervalent
   charge offset; reject excess occupancy only under the applicable strict rules.

Confirmed corrections include radical checks when inference is disabled,
radical/charge-only excess hidden by saturating subtraction, original-element
unrestricted and zero-valence exceptions, isolated hydrogen charge behavior,
charge clamping, and the ordering of hypervalent-anion inference. The duplicate
default-valence lookup was removed; default and allowed valences now derive from
one table.

Diagnostics distinguish ordinary explicit-valence overflow, occupied-valence
overflow, and unreasonable isolated-hydrogen charge. Occupancy errors retain
the explicit contribution, radical count, charge offset and unsubtracted target
limit, avoiding contradictory reports such as explicit valence zero exceeding
an allowed value of zero.

Assignments are collected before installation. Failure preserves all installed
perception. Success retains ring information and invalidates dependent
aromaticity/CIP. Represented graph chemistry is unchanged by valence perception.

## Ring findings

The audit compared complete selected ring families, not only counts or atom/bond
membership. It found differences in degree-three search ordering, when extra
rings were pruned, the set of witnesses used for symmetrization, and recovery
when the initial search did not supply the reference's expected ring count.

The revised implementation follows the reference's connecting-cycle recovery
and iterative depth-first fallback. Symmetrization compares extras only against
the original selected rings. It does not let previously accepted extras become
new replacement witnesses. The PubChem CID125634 regression protects a missing
six-membered ring: the original implementation returned eight selected rings;
the reference and corrected implementation return nine.

The short-ring helper now uses the same edge policy as normal ring perception:
zero-order and dative edges do not form chemical rings. Search setup, recovery
copies/scans and ring comparisons are charged to the existing total-work budget,
and resource failures leave previous perception intact. Degree-two root selection
now scans each fragment once instead of repeatedly scanning its prefix. Obsolete
bespoke recovery code was removed.

### Confirmed limitations in the reference

RDKit's selected ring family is not guaranteed to have full cycle-space rank.
For example, the connected six-vertex graph with edges
`(0,2),(0,3),(0,4),(0,5),(1,2),(1,3),(1,4),(1,5),(2,5),(3,4)`
has cycle rank `10 - 6 + 1 = 5`. Both `GetSSSR` and `GetSymmSSSR` in the pinned
reference return only four triangles, with binary cycle rank four. Kekule does
not invent a stronger mathematical guarantee for an RDKit-compatible ring model.

Four other generated graphs expose a more serious reference defect: its selected
ring family omits a bond that is cyclic according to an independent alternate
path and RDKit's own `FastFindRings`. Kekule preserves its existing agreement
between true cycle membership and installed ring coverage. Those cases use a
bounded depth-first fallback instead of copying incomplete reference output.
This is an explicit, regression-tested parity exception.

## Oxyhalogen normalization

The existing publication cleanup used a terminal-oxygen heuristic. That missed
perchlorate esters, including the prior audit's CID80794 failure, while rewriting
some iodine compounds that the reference leaves unchanged.

The corrected predicate requires a neutral Cl/Br/I center, represented explicit
valence three, five or seven, and oxygen-only neighbors. Bridging ester oxygen
is permitted. Double-bonded oxygen is converted to single-bonded O- and the
halogen charge records the number of converted bonds. Focused tests cover the
ester, declared-H and already-charged oxygen cases, alongside carbon-neighbor
and even-valence negative controls.

This rewrite happens before canonical publication, not inside perception.
Stereo assertions focused on changed bonds are pruned and stale perception is
cleared. Graphs outside the modeled cleanup domain remain unchanged for the
valence validator to assess.

The predicate bounds the result to at most +3. Consequently the former
`FormalChargeOutOfRange` variants on `MoleculePublicationError` and
`NormalizationError`, and their internal fallible wrapper, were removed.
This is a source-API removal. The new `InvalidFormalCharge` and
`ValenceOccupancyExceeded` valence issue variants are source-API additions.
Existing error consumers may need to update exhaustive matches.

## Validation methodology

The valence golden generator calls `UpdatePropertyCache(strict=False)` on raw
input; it does not first run complete sanitization. Kekule's canonical
publication can change resonance representations and materialize a hydrogen
needed to preserve source stereo. Raw-golden differences therefore remain
reported separately from algorithm checks on identical represented chemistry.

No expected fields were removed. Canonical-graph checks preserve element,
isotope, charge, radical state, explicit-H declaration, inference policy,
atom-map identity, ordered endpoints and localized bond order before asking
RDKit to calculate the complete valence record. Corpus ring comparison uses the
existing comparator's ordering normalization only: sort atom IDs within each
ring, then sort the ring list. The diagnostic graph sweep additionally compares
complete bond-cycle families. Membership comparisons retain every atom and
bond boolean.

Generated graphs and atom states are isolated diagnostic probes and focused unit
regressions, not replacements for externally supplied benchmark fixtures.
Tracked fixtures and goldens are unchanged. Reproduction scripts, pinned source
files, hashes and detailed results are retained under
`target/valence-rings-audit/`.

### Diagnostic sweeps

- **322,140 atom states**, covering all 118 supported elements, 21 charges
  including both `i8` endpoints, five radical counts, 13 explicit-H counts and
  both hydrogen-inference policies. Strict and permissive success/failure and
  implicit-H results give **644,280 comparisons with zero mismatches**. The
  original implementation disagreed on 28,364 of these states. Bonded-atom unit
  regressions additionally cover localized bond orders, including quadruple bonds.
- **35,432 connected graphs**, exhausting labeled graphs with five and six
  vertices and adding seeded samples with seven, eight, ten and twelve vertices
  (sampled maximum degree four; seed `20260911`). **35,428 complete atom and
  bond-cycle families match RDKit exactly**. The remaining four are the same
  cycle-coverage exceptions described above; their fallback families match
  RDKit's `FastFindRings`. **All 35,432 atom/bond membership masks match**, with
  zero perception errors. The original implementation differed on 352 graphs.

### External corpora

Smoke, PL-REX, PubChem-100, PubChem-1k, PubChem-100k and Enamine diversity supply
**594,506 feature-record comparisons** across valence, fast rings, selected rings
and dependent aromaticity. These are record occurrences, not unique compounds;
the corpora overlap. All 493 valence/ring feature-fixture SHA256 checks match the
unchanged reference inputs. Final original-input results are:

| Feature | Exact matches | Source-stereo interpretation failures | Other differences |
| --- | ---: | ---: | --- |
| Valence | 143,092 | 4,212 | 1,322 raw representation differences |
| Fast ring membership | 144,415 | 4,212 | None |
| Selected ring families | 144,414 | 4,212 | None |
| Aromaticity masks | 144,406 | 4,212 | 9 records RDKit also rejects |

All **1,322 raw valence differences match every RDKit valence field on the same
represented graph**. Of these, 1,161 Enamine records differ in explicit/implicit
H storage due to source-stereo materialization; 161 PubChem record occurrences
differ in oxyhalogen charge/bond representation. The corrected normalization adds
CID80794 and CID159402 to that raw-representation group; both pass the same-graph
check. Original valence-golden differences remain visible.

The 4,212 source-stereo failures are 4,211 Enamine records and one PubChem-1k
record. In a separate experiment, temporary copies clear only V2000 stereo
fields. All 4,212 projected records match the original ring/aromaticity goldens,
with pinned RDKit independently confirming those masks are unchanged. The
projected canonical graphs also give **4,212 exact full-field valence matches**.
Clearing stereo changes RDKit's explicit-H materialization on the Enamine
records, so those projected valence results are not relabeled as original-golden
successes.

CID125634 now matches all nine selected rings on original input. CID80794 now
passes the full original-input pipeline. The nine remaining aromaticity
perception failures are exactly the source records rejected by the pinned
reference's normalization/perception pipeline; no additional mask mismatch or
perception failure remains among successfully interpreted records.

The native benchmark runner still aborts whole packs on source-stereo errors
and assumes one connected component per SDF record. The isolated diagnostic
probe instead records each failure and recombines components by original
atom/bond identity, preserving every asserted field. This audit does not change
that runner. Detailed corpus tables, exact expected/actual records and
reproduction commands are retained in
`target/valence-rings-audit/PARITY-VALIDATION.md` and `final-summary.json`.

### Representation boundaries

Canonical `Element` excludes dummy/query atoms, and aromatic source bonding is
localized before perception. These APIs do not implement every RDKit sanitizer
or query-graph feature. Dative bonds currently have no defined donor/acceptor
valence contract across graph canonicalization, source interpretation and CIP;
their existing zero contribution at both endpoints differs from RDKit's
directional acceptor contribution. Resolving that requires a coordinated graph
contract change rather than changing one valence sum in isolation.

## Rust validation

All required checks passed on the Windows MSVC host:

```text
cargo fmt --all -- --check
cargo check --workspace --all-targets --all-features --locked
cargo +1.89.0 check --workspace --all-targets --all-features --locked
cargo check --manifest-path fuzz/Cargo.toml --bins --locked
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo test --workspace --all-features --locked
cargo test --workspace --all-features --doc --locked
cargo doc --workspace --all-features --no-deps --locked
cargo test -p kekule-potentials --no-default-features --locked
cargo doc -p kekule-potentials --no-default-features --no-deps --locked
cargo package -p kekule --locked --allow-dirty
cargo package -p kekule-potentials --locked --list --allow-dirty
cargo package -p kekule-traj --locked --list --allow-dirty
cargo build -p kekule --lib --locked
git diff --check
```

Documentation ran with `RUSTDOCFLAGS=-D warnings`. All six packaged license
copies match the repository licenses by SHA256. `--allow-dirty` permits checking
the uncommitted changes; the foundational package was compiled from its packaged
contents. The final workspace run includes all 695 foundational library unit
tests, its integration tests, companion-crate tests and documentation tests.

No required formatting, check, clippy, test, documentation or package gate was
skipped. Full companion package builds remain deferred by the repository's CI
policy until the foundational crate is published; their exact file sets were
checked. Nightly fuzz execution belongs to the separate Linux fuzz workflow and
was not run during this audit; all registered fuzz binaries compiled successfully.

Changes are on `codex/valence-rings-audit`, based on `main` at `7bba7c08`.
`README.md`, external fixtures and golden outputs were not modified.
