# Aromaticity audit — 2026-09-11

The parity target is **RDKit 2026.03.3**, default/RDKIT aromaticity, matching the
version pinned by `benchmarks/reference/rdkit/environment.yml`. This audit covers
the complete localized-graph aromaticity path, its valence and ring inputs,
publication/invalidation behavior, and the existing external-reference tests.

## Findings and changes

1. **Incomplete element data.** The old donor lookup omitted Be, Mg, Al and Si.
   It also reused that restricted lookup for exocyclic neighbors, so halogens
   and transition metals could fail to withdraw electrons. The lookup now uses
   RDKit's complete outer-electron data. All 119 entries, including the unused
   zero slot, were checked against the pinned executable reference. Regressions
   assert complete atom and bond masks for eligible and ineligible elements,
   electron-withdrawing neighbors, and equal/lower-electronegativity controls.
2. **Incorrect zero-electron donor classification.** An atom with zero available
   electrons and an exocyclic multiple bond was rejected. It must remain a
   vacant donor. A represented carbon-radical regression checks both the donor
   classification and the resulting complete aromatic masks.
3. **Different unsaturation arithmetic.** Counting only extra bond orders did
   not reproduce the reference's explicit-valence-minus-raw-degree calculation.
   Declared hydrogens and zero-order neighbors matter to that calculation. The
   corrected calculation is shared by donor capping and candidate exclusion;
   a cyclic-allene/zero-order-neighbor regression exercises the difference.
   Multiple-bond detection also includes represented quadruple bonds.
4. **Standalone/full-pipeline disagreement.** Aromaticity contained a separate
   implicit-hydrogen target table. A five-membered sulfur cation with hydrogen
   inference enabled produced no aromatic flags standalone, but five aromatic
   atoms and bonds through `Molecule::perceive()` and RDKit. Both paths now use
   one valence calculation. Standalone perception still respects installed
   hydrogen overrides and does not install valence state. The regression was
   observed failing before this fix; a second test protects partial overrides.
5. **Combinatorial fused-ring traversal.** The old search generated combinations
   and materialized connected subsets before evaluating them, then repeated
   connectivity checks with a different overlap predicate. The replacement
   computes one fusion adjacency map, visits connected subsets directly, and
   stops immediately when all component bonds are marked. A 300-ring chain has
   only 295 connected six-ring subsets but over 962 billion unrestricted
   six-ring combinations. A regression covers that case and another verifies
   immediate cancellation. Exhaustive comparison against an independent mask
   and reachability oracle covers every undirected five-node graph, checking
   both completeness and absence of duplicate subsets.

Unused private options for unimplemented alternative aromaticity models, the
duplicate hydrogen table, and orphaned union-find helpers were removed. The
public API remains unchanged. Aromaticity entry points now document cached ring
and hydrogen reuse and transactional failure behavior.

## Algorithm review

| Stage | Verified behavior |
| --- | --- |
| Input | Localized graph bonds; source aromatic syntax is interpreted before publication. |
| Valence | Installed H counts take precedence; otherwise use the same atom-level calculation as valence perception. |
| Rings | Reuse installed ring basis, or compute the symmetrized SSSR-like basis with explicit resource limits. |
| Donors | Default-valence, degree, lone-pair, charge and radical arithmetic; cyclic/exocyclic multiple-bond handling; 0/1/2-electron contribution. |
| Candidates | Element, coordination, charge-adjusted valence, radical and multiple-unsaturation restrictions. |
| Fusion | Exactly one shared bond; rings of at most 24 bonds may fuse. Larger rings remain eligible individually. |
| Search | Connected subsets through six rings; components exceeding 300 rings stop at pairs. |
| Electron count | Count atoms occurring in one or two selected rings; exclude buried atoms occurring in three or more. Accept 2 or 4n+2 electrons. |
| Assignment | Mark bonds occurring once in a successful subset; mark endpoints through single/double bonds. Fusion bonds need independent support. |
| State | Replace derived aromaticity, invalidate dependent CIP, preserve localized chemistry, and roll back perception on errors. |

These reference rules were reviewed against the pinned
[aromaticity implementation](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/Aromaticity.cpp),
[periodic-table comparison](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/PeriodicTable.h),
and [atomic data](https://github.com/rdkit/rdkit/blob/Release_2026_03_3/Code/GraphMol/atomic_data.cpp).
The [RDKit Book](https://www.rdkit.org/docs/RDKit_Book.html#aromaticity) provides
the conceptual model; the pinned implementation and executable settle edge cases.

The new traversal can visit subsets in a different order from RDKit. Donors are
fixed before traversal, and successful subsets only add flags. Therefore order
does not alter the union of aromatic flags, including when the search ends after
all component bonds have been marked. An independent review checked this and
the recursive exclusion/unwind invariants.

Existing regressions also cover cached-ring behavior, resource-limit rollback,
neutral carbon radicals, charged carbons, aromatic fusion singles, localization,
and perception invalidation. No ownership or publication invariant was changed.

## Scope of parity

Parity means aromatic membership on equivalent represented localized chemistry.
It does not mean that every RDKit input, sanitizer operation or aromaticity model
is supported:

- Canonical `Element` excludes dummy/query atoms. RDKit's variable-electron dummy
  handling is outside this graph model.
- `BondOrder::Dative` currently has no directional valence semantics. Kekule
  counts it as zero at both endpoints; RDKit counts the acceptor contribution.
  Changing this requires a separate graph/valence contract change.
- Source localization deliberately accepts some forms, such as `cc`, without
  asserting perceived aromaticity. This is different from RDKit's full SMILES
  sanitization contract.
- Installed custom ring bases and hydrogen overrides are trusted inputs. Ring
  enumeration remains Kekule's own deterministic implementation, not a copied
  RDKit ring finder. Corpus agreement is evidence, not proof for every graph.
- The inherited six-ring, 24-atom fusion and 300-ring search restrictions remain.
  The search does not establish unrestricted aromaticity for arbitrary graphs.

RDKit remains an optional external reference, with no Rust runtime dependency.
Externally supplied fixtures and golden files were not changed.

## External-reference validation

All 165 external fixture files were SHA256-checked against their tracked
references. Every golden identifies RDKit 2026.03.3. The comparison retains
complete source-ordered atom and bond boolean arrays, record identity and status;
it does not compare just aromatic counts. Multi-component records are compared
by recombining their components through source mappings.

| Corpus | Records | Exact original-record matches | Upstream failures before aromatic comparison |
| --- | ---: | ---: | --- |
| Smoke | 4 | 4 | 0 |
| PL-REX | 164 | 164 | 0 |
| PubChem 100 | 100 | 100 | 0 |
| PubChem 1k | 1,000 | 999 | 1 source stereo code |
| PubChem 100k | 100,000 | 99,990 | 10 valence failures |
| Enamine diversity | 47,359 | 43,148 | 4,211 source stereo/publication failures |
| **Total record occurrences** | **148,627** | **144,405** | **4,222** |

There were **zero aromatic-mask mismatches** among successfully processed
original records. Counts are record occurrences; the corpora are not asserted
to be mutually disjoint. The full assay was repeated against the final source
and produced identical results.

The standard `cargo xtask benchmark --benchmark algo.aromaticity.rdkit-like
--corpus all` cannot complete: source-stereo interpretation aborts Enamine packs,
and its one-molecule SDF reader panics on disconnected PubChem records. Smoke
passes 4/4 fixtures and PL-REX 2/2. A separate diagnostic probe under
`target/aromaticity-audit/parity-probe/` permits per-record failure accounting
and all-component comparison without changing the production harness or goldens.

Two further experiments isolate upstream failures; these are **not** claimed as
successful original-input interpretation:

- Clearing only V2000 stereo fields on the 4,212 stereo-failing records yields
  **4,212 exact aromatic-mask matches**. The pinned live RDKit was run on both
  original and projected inputs and confirmed that both full masks equal the
  original golden for every record. Elements, charges, radicals, hydrogens,
  connectivity and bond orders were preserved.
- Of the ten valence failures, nine are also rejected by RDKit. The remaining
  record, **PubChem CID80794**, is acyclic tert-butyl perchlorate. RDKit cleanup
  changes Cl to +3, three O atoms to -1 and three Cl=O bonds to singles. Kekule's
  existing oxyhalogen-cleanup guard requires a terminal oxygen and skips the
  bridging ester oxygen. Applying the reference normalization separately yields
  the same all-false masks for all 18 atoms and 17 bonds. This is an independent
  represented-chemistry normalization gap, left outside this aromaticity change.

Thus every one of the 148,618 RDKit-successful record occurrences has matching
aromatic masks either directly or after the explicitly separated upstream
projections. The original-input success count remains 144,405.

Raw logs, source hashes, projection scripts, the diagnostic probe and validation
commands are retained locally under `target/aromaticity-audit/`.

## Rust validation

All of the following passed on the Windows MSVC host:

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
git diff --check
```

Documentation ran with `RUSTDOCFLAGS=-D warnings`. All six packaged license
copies match the repository licenses by SHA256. `--allow-dirty` permits checking
the reviewable, uncommitted changes; the foundational package was compiled from
its packaged contents. Focused aromaticity tests passed (13 unit tests and two
integration tests), and the shared-valence extraction also passed 27 focused
valence tests before the full workspace run.

No applicable formatting, check, clippy, test, documentation or package gate was
skipped. Full companion package builds remain deferred by the repository's CI
policy until the foundational crate is published; their exact file sets were
checked. Nightly fuzz execution belongs to the separate fuzz workflow and was
not run during this audit; all registered fuzz binaries compiled successfully.

Changes are on `codex/aromaticity-audit`, based on current `main` at audit start.
`README.md`, external fixtures and golden outputs were not modified.
