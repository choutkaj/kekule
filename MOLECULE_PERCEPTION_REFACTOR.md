# Molecule / Perception Refactor Plan

This plan implements the canonical `Molecule` architecture defined in
`ARCHITECTURE.md`.

The target invariant is simple:

> A published `Molecule` is always one connected canonical Kekule
> representation. Source-only representation is resolved during interpretation;
> derived chemistry lives in `PerceptionState`.

Keep the refactor incremental. Each stage should compile, pass tests, and avoid
unrelated API redesign.

## Stage 1 — Remove aromaticity from canonical `BondOrder`

- Remove `BondOrder::Aromatic` from the canonical core enum.
- Introduce source-/interpreter-local aromatic bond representation where needed
  by SMILES, molfile, or other readers.
- Reuse the existing deterministic aromatic-localization logic during
  interpretation before a `Molecule` is published.
- Update writers, query code, tests, and algorithms so canonical `Molecule`
  bonds are always localized.
- Preserve perceived aromaticity exclusively through `PerceptionState`.
- Do not redesign ring or aromaticity algorithms beyond what this migration
  requires.

Exit condition: no published `Molecule` can contain an aromatic bond order, and
all existing aromatic input/output behavior remains covered by tests.

## Stage 2 — Fold public normalization into interpretation/publication

- Remove the ordinary public `normalize()` lifecycle from `SmallMolecule`.
- Make format interpretation return canonical molecules directly.
- Keep reusable canonicalization helpers internally where useful; deleting
  working algorithms merely to remove the word "normalization" is not a goal.
- Ensure checked manual construction/editing cannot publish source-only or
  otherwise noncanonical core state.
- Replace `normalize_and_perceive()` with the simpler canonical
  interpret/read-then-`perceive()` workflow as appropriate.

Exit condition: the public chemistry pipeline is parse -> interpret -> perceive,
with no public unnormalized `Molecule` state.

## Stage 3 — Remove source stereo marks from canonical `Molecule`

- Move `StereoBondMark` and equivalent wedge/directional source state into
  format-specific documents or private interpreter staging.
- Resolve source stereo into canonical `StereoElement` / `StereoGroup` state
  during interpretation.
- Remove source-mark storage and mutation APIs from canonical `Molecule`.
- Preserve source diagnostics/provenance in interpretation reports rather than
  chemical identity.

Exit condition: a published `Molecule` contains canonical stereo representation
only; no unresolved source stereo marks remain.

## Stage 4 — Audit canonical stereo metadata

Review `StereoSource` and `StereoSpecifiedness` against the represented-vs-source
boundary.

- Move pure source provenance out of canonical stereo objects.
- Keep only specifiedness distinctions that carry chemical semantics after
  interpretation.
- Treat invalid/cleared source-state diagnostics as interpretation output rather
  than canonical molecular state where possible.
- Preserve exact persistence/reconstruction requirements intentionally rather
  than accidentally.

Exit condition: canonical stereo state contains chemistry, not parser history.

## Stage 5 — Clean up represented hydrogen declarations

Evaluate replacing the current `explicit_hydrogens` +
`no_implicit_hydrogens` pair with one explicit represented-hydrogen abstraction
if it materially simplifies invariants.

Requirements:

- preserve source-represented hydrogen semantics;
- keep inferred implicit-H counts exclusively in valence perception;
- preserve SMILES/molfile behavior and hydrogen add/remove transforms;
- avoid changing chemistry merely for API aesthetics.

This stage may conclude that the current representation is already preferable;
if so, document that decision and make no churn.

## Stage 6 — Regularize `PerceptionState`

- Keep valence, rings, aromaticity, and stereo-derived information sectional.
- Introduce a `StereoPerceptionState` wrapper for CIP descriptors if it improves
  symmetry without needless API churn.
- Distinguish graph cycle membership from an algorithm-selected ring basis.
- Add ring-basis model/provenance only where a downstream algorithm can actually
  depend on that choice.
- Keep invalidation simple and dependency-safe; do not build a generic cache or
  dirty-bit framework.

Exit condition: `PerceptionState` is a small, coherent container of fundamental
derived chemistry with no duplicate atom/bond flags.

## Stage 7 — Ergonomic molecule-aware views (optional)

If useful after the core refactor, add lightweight atom/bond views so callers can
access represented and perceived chemistry through one ergonomic API without
physically duplicating derived state into `Atom` or `Bond`.

Examples:

```text
atom.element()
atom.formal_charge()
atom.implicit_hydrogens()
atom.is_aromatic()
atom.is_in_ring()
```

This is an API convenience stage, not a prerequisite for the architectural
invariants above.

## Non-goals

- no generic perception registry;
- no continuous/learned perception implementation yet;
- no force-field typing or descriptor migration into `PerceptionState`;
- no broad rewrite of working ring/aromaticity/CIP algorithms;
- no unrelated `Topology`, `Model`, or trajectory refactor.
