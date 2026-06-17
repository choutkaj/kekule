# Implementation Plan

## Feature ID

`core.atom-bond`

## Goal

Define the chemically general atom and bond record model used by the shared core `Molecule` graph. This feature should make atom and bond data explicit, typed, and suitable for both small-molecule and macromolecular workflows without adding perception, sanitization, parsing, or biomolecular hierarchy fields.

`core.atom-bond` depends on `core.graph`; graph topology, stable IDs, mutation, deletion, and adjacency behavior remain owned by `Molecule`.

## Public API

Refine or complete the existing public model around:

- `Element` as a chemically general element identifier.
- `Atom` as the chemically general atom payload stored in `Molecule`.
- `Bond` as the chemically general bond payload stored in `Molecule`.
- `BondOrder` for raw or parsed bond order values.
- `AtomStereo` and `BondStereo` as carried annotations only, not assigned perception results.
- `PropMap` and `PropValue` for arbitrary molecule, atom, and bond properties.

Expected constructors and accessors:

- `Element::from_atomic_number(u8) -> Option<Element>`.
- `Element::from_symbol(&str) -> Option<Element>`.
- `Element::atomic_number() -> u8`.
- `Element::symbol() -> &'static str`.
- `Atom::new(element: Element) -> Atom`.
- `Bond::new(a: AtomId, b: AtomId, order: BondOrder) -> Bond`.
- `Bond::a() -> AtomId`, `Bond::b() -> AtomId`, and `Bond::endpoints() -> (AtomId, AtomId)`.
- Read/write access to chemically general atom and bond fields, either through public fields or narrow setters where invariants require them.

Do not add public APIs for valence perception, ring detection, aromaticity perception, stereochemistry assignment, canonicalization, file parsing, or biomolecular hierarchy labels.

## Internal Modules Touched

Expected scope:

- `crates/molecules/src/lib.rs` initially, unless the model becomes large enough to justify a small module split such as `atom.rs`, `bond.rs`, or `chem.rs`.
- `prelude` exports for the finalized public types.
- Unit tests colocated with the model.
- Feature documentation under `features/core.atom-bond/` if the public model changes materially.

Do not touch parsers, validation generators, or downstream algorithm features except where tests need to compile against the finalized atom and bond API.

## Data Model

### Element

`Element` should store atomic number as the canonical internal representation.

Required behavior:

- Accept atomic numbers `1..=118`.
- Reject `0` and values above `118`.
- Support common symbols needed by initial tests and fixtures.
- Use canonical case-sensitive symbols, for example `Cl`, not `CL`.
- Keep symbol parsing deterministic and independent of sanitization.

The implementation may start with a hand-written match table and later move to a generated static table if coverage grows.

### Atom

`Atom` fields should remain chemically general:

- `element: Element`.
- `isotope: Option<u16>`.
- `formal_charge: i8`.
- `radical_electrons: u8`.
- `explicit_hydrogens: u8`.
- `implicit_hydrogens: Option<u8>`.
- `aromatic: bool`.
- `chiral: Option<AtomStereo>`.
- `atom_map: Option<u32>`.
- `props: PropMap`.

These fields are allowed because they are useful for both small molecules and macromolecular extracted chemistry. They do not include model, chain, residue, author atom name, alternate location, occupancy, B-factor, or mmCIF label fields; those belong in `BioHierarchy`.

`Atom::new(element)` should initialize neutral, non-aromatic, non-radical, non-stereochemical defaults with an empty property map.

### Bond

`Bond` should store:

- Private endpoints `AtomId` and `AtomId`.
- `order: BondOrder`.
- `aromatic: bool`.
- `stereo: Option<BondStereo>`.
- `props: PropMap`.

Endpoints should remain private so graph topology cannot be rewritten behind `Molecule` adjacency storage. Endpoint mutation, if ever needed, must be a `Molecule` topology operation that updates adjacency and invalidates perception state.

`Bond::new(a, b, order)` should set endpoints, store the order, initialize an empty property map, and set `aromatic = true` only when `order == BondOrder::Aromatic`.

### BondOrder

Initial variants:

- `Zero`.
- `Single`.
- `Double`.
- `Triple`.
- `Quadruple`.
- `Aromatic`.
- `Dative`.

Bond order storage is descriptive only. This feature must not validate valence, infer aromaticity, or normalize parsed bond orders.

### Stereo Annotations

`AtomStereo` and `BondStereo` represent stored annotations, not assigned perception results.

Initial variants can include:

- `AtomStereo::TetrahedralClockwise`.
- `AtomStereo::TetrahedralCounterClockwise`.
- `AtomStereo::Unspecified`.
- `BondStereo::E`.
- `BondStereo::Z`.
- `BondStereo::Up`.
- `BondStereo::Down`.
- `BondStereo::Unspecified`.

Future stereochemistry work may replace or extend this model, but this feature should only preserve annotations.

### Properties

Use `BTreeMap<String, PropValue>` for deterministic property ordering.

Initial `PropValue` variants:

- `String(String)`.
- `Int(i64)`.
- `Float(f64)`.
- `Bool(bool)`.

Property maps should be available on atoms and bonds and remain opaque to chemistry perception unless a later feature explicitly documents a chemically interpreted property.

## Algorithm Outline

This feature is mostly data modeling rather than algorithmic chemistry.

1. Audit existing `Atom`, `Bond`, `Element`, stereo, bond-order, and property APIs against this plan.
2. Fill missing accessors and constructors without changing graph topology semantics.
3. Keep bond endpoint mutation unavailable outside `Molecule`.
4. Ensure defaults are explicit and deterministic.
5. Ensure changing atom or bond payload data through `Molecule` APIs continues to invalidate computed state conservatively where mutable access can affect chemistry.
6. Keep parsing and perception separate: parsed input may populate fields, but this feature does not sanitize, infer, or validate chemistry.

## Unit Tests

Add or update unit tests for:

- `Element::from_atomic_number` accepts `1` and `118`.
- `Element::from_atomic_number` rejects `0` and `119`.
- `Element::from_symbol` accepts representative one-letter and two-letter symbols.
- `Element::from_symbol` rejects unknown or incorrectly cased symbols.
- `Element::symbol` and `Display` return canonical symbols.
- `Atom::new` initializes all defaults.
- Atom isotope, formal charge, radical, hydrogen, aromatic, stereo, atom-map, and property fields can be set and read.
- `Bond::new` stores endpoints and order.
- `Bond::new` sets `aromatic` only for `BondOrder::Aromatic`.
- Bond endpoint accessors return the original `AtomId`s.
- Bond properties and stereo annotations can be set and read.
- `PropValue` equality works for string, int, float, and bool values.
- `SmallMolecule` and `MacroMolecule` continue to use `Molecule` atom and bond payloads without duplicating atom or bond models.

If mutable atom or bond access remains on `Molecule`, include tests that chemistry-relevant payload mutation invalidates fresh perception state.

## Reference Validation

No RDKit or Biopython golden data is required for the first implementation. This feature defines local data representation, not reference chemistry behavior.

If future validation data is added, use normalized JSON with:

- Reference tool name and version, when applicable.
- Input fixture identifier.
- Atom field values.
- Bond field values.
- Property values.
- Notes describing any normalization performed before comparison.

Do not mark `validated = true` unless there is reference-generated golden data or documented manual validation.

## Required Checks

Run:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
cargo xtask validate --feature core.atom-bond
```

Regenerate the dashboard only if feature metadata changes. This planning-only update should not change `implemented` or `validated`.

## Risks

- Overfitting atom fields to small molecules could leak biomolecular labels into core `Atom`.
- Making bond endpoints public would allow topology changes without adjacency updates.
- Treating `aromatic` or stereo fields as perceived truth too early could blur raw parsing and perception.
- Expanding `Element` symbol coverage manually can introduce table mistakes; tests should cover representative boundary and common cases.
- `PropValue::Float` equality is exact; future serialized validation may need explicit float normalization rules.

## Edge Cases

- Unknown element symbols should be rejected rather than mapped to a placeholder element.
- `Element::symbol()` should only return `"?"` for impossible internal states, if such states remain constructible internally.
- Isotopes should allow mass numbers above the natural range because parsers may carry unusual labels; detailed isotope validation is out of scope.
- Formal charge range is bounded by `i8`, but chemically unreasonable charges are not rejected here.
- `BondOrder::Zero` and `BondOrder::Dative` are storage values only and do not imply valence behavior.
- `AtomStereo::Unspecified` and `BondStereo::Unspecified` should preserve explicit source annotations distinct from `None`.

## Explicitly Out of Scope

- Ring detection.
- Aromaticity perception.
- Valence perception or sanitization.
- Stereochemistry assignment.
- Canonicalization.
- SDF, PDB, or mmCIF parsing.
- RDKit or Biopython runtime dependencies.
- Biomolecular hierarchy fields on core `Atom`.
- Any change to `core.graph` stable ID or adjacency semantics.
