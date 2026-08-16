# Small-Molecule Pipeline Refactor Plan

## Purpose

This document is the implementation plan for the small-molecule chemistry pipeline defined in `ARCHITECTURE.md`.

The target pipeline is:

```text
source text / bytes
    -> parse
    -> format-specific *Document
    -> interpret
    -> represented Molecule
    -> normalize
    -> normalized Molecule
    -> perceive(model)
    -> Molecule + derived perception state
    -> optional CIP assignment or other specialized derivation
```

The stages have fixed meanings:

- **Parse** reads source syntax only.
- **Interpret** maps source-asserted chemistry into format-independent represented chemistry without general chemical inference.
- **Normalize** deterministically rewrites equivalent representation into Kekule's canonical internal form without changing represented chemical meaning.
- **Perceive** derives model-dependent chemical interpretation without rewriting primary representation.
- **CIP** derives stereochemical descriptors as an explicit specialized step.
- **Standardization** is a future separate layer for chemistry-changing choices such as protonation, tautomer selection, salt handling, or fragment-parent selection.

This plan is intentionally implementation-oriented. `ARCHITECTURE.md` remains the normative source for the conceptual design.

## Target invariants

After **interpretation**:

- the result is format-independent represented chemistry;
- only source-asserted chemistry has been added;
- no general valence, ring, aromaticity, or coordinate-stereo inference has run;
- source aromatic syntax and source stereo marks may still be present;
- ordinary interpretation does not install a chemical `PerceptionState`.

After **normalization**:

- represented chemical meaning is unchanged;
- normalization is deterministic and idempotent;
- normalization is model-independent;
- `BondOrder::Aromatic` is absent from normalized primary representation;
- resolvable source stereo marks have become canonical represented `StereoElement` state;
- unresolved or invalid source representation is reported structurally rather than guessed through a perception model;
- successful normalization consumes source stereo marks that have been represented canonically;
- installed perception is cleared: normalization publishes represented chemistry, not derived interpretation.

After **perception**:

- primary represented atoms, bonds, bond orders, charges, explicit-H declarations, stereo elements, stereo groups, source-independent properties, and conformers are unchanged;
- derived valence/implicit-H, rings, and aromaticity are installed in `PerceptionState`;
- different perception models may legitimately produce different derived views of the same normalized molecule.

After **CIP assignment**:

- primary representation remains unchanged;
- descriptor assignment is transactional and remains an explicit specialized operation.

## Classification and migration map

### Interpret

Keep or move here:

- elements and isotopes explicitly encoded by the source;
- formal charges explicitly encoded by the source;
- radical state when directly encoded by the source;
- atom maps and source properties;
- asserted bond endpoints and represented bond order;
- explicit hydrogen declarations and source no-implicit-H policy;
- coordinates into conformers;
- explicit stereo groups;
- SMILES `@`/`@@` local stereo when directly resolvable from source semantics;
- SMILES `/` and `\\`, Molfile wedge/hash/either, and similar source stereo syntax into format-independent `StereoBondMark` state when not already canonical local stereo;
- source aromatic syntax as represented source chemistry, not installed semantic aromaticity.

Refactor later:

- SMILES parsing currently builds a private `SmilesProgram` using core chemistry types. Parsing should ultimately retain format-specific syntax/AST state and interpretation should own conversion into `Atom`, `BondOrder`, stereo source state, etc.
- V2000/V3000 parsing currently constructs core chemistry values early. The eventual target is the same parse/interpret separation.

Remove from interpretation:

- `AromaticityProvenance::Imported` installation from SMILES interpretation;
- RDKit-like radical inference performed as a general chemistry guess during SMILES interpretation;
- Molfile hydrogen inference through general allowed-valence rules.

### Normalize

Move or implement here:

- deterministic imported-aromatic localization / Kekulization into ordinary localized bond orders;
- canonical meaning-preserving charge/bond normalization such as the current supported hypervalent oxo-halide cleanup;
- source stereo mark assembly into canonical `StereoElement` representation;
- structural validation needed to prove that the represented source form can be converted without inventing chemistry;
- complete invalidation/removal of installed perception state before publishing the normalized representation.

Normalization must not depend on:

- installed implicit-H assignments;
- a named valence perception model;
- perceived aromaticity;
- perceived ring state except for purely graph-theoretical temporary work that is not installed as semantic perception;
- coordinate-derived chemical inference.

Specific rewrites required:

1. The current imported-aromatic matching/Kekulization code in `algorithms/aromaticity.rs` mutates bond orders but consults perceived implicit H and RDKit-like valence rules. Do not simply move it. Rebuild the normalization kernel so it operates from represented source state and fixed normalization rules only.
2. The source-stereo half of current `perceive_stereo()` depends on perceived H/ring/aromatic state in several helpers. Split and rewrite those helpers so decoding source assertions does not require a chemical perception profile.
3. `normalize_aromatic_nitrogen_hydrogens()` should be removed rather than moved. Perceived aromaticity/implicit-H must not feed back into primary representation merely to imitate another toolkit's sanitized atom state.

### Perceive

Keep here:

- `perceive_valence*` and implicit-H assignment;
- ring membership;
- deterministic ring basis;
- aromaticity perception;
- read-only stereo-candidate detection;
- coordinate-derived stereochemical inference as a derived result.

Target default discrete order:

```text
valence / implicit H
    -> rings
    -> aromaticity
```

Important constraints:

- aromaticity perception must operate on normalized localized bond orders and must no longer Kekulize or otherwise rewrite `BondOrder`;
- perception may install or replace `PerceptionState`, but must not mutate primary represented chemistry;
- the default `perceive()` operation should be one well-defined discrete profile, not a new bag of independent boolean flags;
- experts continue to have focused valence/ring/aromaticity operations for custom workflows.

### Coordinate-derived stereo

Current coordinate stereo inference is genuinely perception, but current publication semantics violate the new boundary because the mutating stereo perception function installs new `StereoElement`s.

Target:

- coordinate stereo inference is read-only and returns derived candidates/results;
- applying/materializing inferred coordinate stereo into represented `StereoElement` state is a separately named explicit transform if needed;
- do not redesign this until the source-stereo normalization split is complete.

### CIP

Keep the current transactional assignment model.

CIP remains separate from the default general perception pipeline and may depend on normalized represented stereo plus the relevant installed perception state.

### Explicit transforms

Keep outside normalization:

- add hydrogens;
- remove hydrogens;
- any graph-atom materialization/collapse operation;
- future chemistry-changing standardization.

The current feature/documentation term `Hydrogen Topology Normalization` should eventually be renamed to avoid confusion with canonical pipeline normalization.

### Delete / retire

Retire once replacements exist:

- `sanitize_*` orchestration functions;
- `SanitizeOptions`;
- `SanitizeReport`;
- `SanitizeError`;
- `SmallMolecule::sanitize*`;
- `SmallMolecule::from_smiles_sanitized`;
- the `chem.sanitize.rdkit-like` feature identity and sanitization terminology;
- `AromaticityProvenance::Imported`;
- `normalize_aromatic_nitrogen_hydrogens()`;
- transitional source-stereo assembly inside general stereo perception.

Do not keep long-lived compatibility wrappers. Kekule is pre-1.0 and this refactor may hard-break the current API.

## Target public API shape

Low-level format use remains explicit:

```rust
let document = smiles::parse_str(input)?;
let mut molecule = smiles::interpret(&document)?.into_molecule();

normalization::normalize(molecule.graph_mut())?;
perception::perceive(molecule.graph_mut())?;
stereo::assign_cip_descriptors(molecule.graph_mut())?;
```

Focused expert perception remains available:

```rust
perception::valence::perceive(...)
perception::rings::perceive(...)
perception::aromaticity::perceive(...)
```

A later ordinary fast-forward convenience should compose normalization + default perception transactionally, for example conceptually:

```rust
molecule.prepare()?;
```

The exact final convenience naming should be decided after the low-level stages are clean. Do not introduce `prepare()` prematurely just to preserve the old sanitizer shape.

## Staged implementation

Each stage should be one focused branch/PR unless a very small follow-up correction is required. Do not implement later stages early.

### Stage 1 — Introduce normalization as a first-class layer

Goal: establish the new architectural seam without yet rewriting imported aromaticity or stereo internals.

Implement:

- [x] add a focused `normalization` module/facade for small-molecule representation normalization;
- [x] define `NormalizationError` and, only if genuinely useful, a minimal success report for representation warnings;
- [x] move the current representation-only hypervalent oxo-halide cleanup out of `sanitize.rs` into normalization;
- [x] make normalization transactional and idempotent;
- [x] successful normalization clears installed `PerceptionState` before publishing the normalized representation;
- [x] add `SmallMolecule::normalize()` convenience if this remains a thin forward to the canonical operation;
- [x] add focused tests for meaning-preserving rewrite, idempotence, transactional failure where applicable, and perception clearing;
- [x] update public facade/docs/features to name normalization explicitly.

Temporary state after Stage 1:

- old sanitization still exists and may call normalization plus its existing perception pipeline;
- imported aromatic localization still lives in aromaticity perception;
- source stereo assembly still lives in stereo perception;
- no attempt yet to make normalized molecules globally `BondOrder::Aromatic`-free.

Do not:

- rewrite aromaticity;
- rewrite stereo;
- change parser/interpretation boundaries;
- delete sanitize APIs yet;
- add a generic pipeline framework.

Acceptance:

- normalization is independently callable and transactional;
- the existing sanitizer delegates its representation cleanup to normalization rather than owning that cleanup itself;
- normalizing twice is identical to normalizing once;
- no perception algorithm has moved yet.

### Stage 2 — Separate imported aromatic localization from aromaticity perception

Goal: make aromatic bond localization representation normalization and make aromaticity purely perceptual.

Implement:

- [x] extract/rewrite imported-aromatic localization into normalization;
- [x] normalization converts accepted `BondOrder::Aromatic` source representation into deterministic localized ordinary bond orders;
- [x] keep the localization kernel independent of installed perception state and selectable chemical perception models;
- [x] move invalid imported-aromatic representation and normalization resource-limit errors into `NormalizationError`;
- [x] make aromaticity perception assume normalized localized representation and never rewrite primary bond orders;
- [x] remove `AromaticityProvenance::Imported`; installed aromaticity is model-perceived semantic state only;
- [x] update SMILES interpretation so source aromatic syntax is represented without installing semantic aromaticity;
- [x] add regressions for aromatic SMILES interpretation -> normalize -> perceive and for normalization idempotence.

Temporary compatibility is allowed only inside this branch while moving callers; do not leave a permanent dual path where aromaticity silently normalizes input.

Acceptance:

- after successful normalization, no live bond has `BondOrder::Aromatic`;
- aromaticity perception does not mutate atoms/bonds/stereo representation;
- interpreted aromatic source has empty semantic aromaticity until perception runs.

### Stage 3 — Simplify valence around normalized localized graphs

Goal: remove historical coupling needed only because valence could receive imported aromatic representation.

Implement:

- [x] audit and simplify RDKit-like valence logic assuming normalized ordinary bond orders;
- [x] remove any valence dependence on pre-installed semantic aromaticity that is no longer chemically necessary;
- [x] preserve correct implicit-H behavior for aromatic systems through localized valence structure;
- [x] retain transactional valence error semantics;
- [x] expand focused regressions for benzene, pyridine, pyrrole-like explicit-H cases, charged aromatics, radicals, and relevant existing corpus fixtures.

Do not redesign `PerceptionState` broadly.

Acceptance:

- default valence can run first on a normalized molecule without requiring aromaticity state;
- resulting implicit H supports subsequent aromaticity perception correctly.

### Stage 4 — Move source-declared stereo assembly into normalization

Goal: make represented stereo complete before chemical perception.

Implement:

- [x] split current stereo source-mark assembly from coordinate-derived stereo inference;
- [x] move wedge/directional/either/source-axis decoding into normalization;
- [x] source-declared stereo becomes canonical `StereoElement` state transactionally;
- [x] successful normalization consumes/removes source marks that were canonically represented;
- [x] normalization warnings/errors replace source-mark warnings/errors currently exposed by stereo perception;
- [x] source-stereo normalization must not require installed valence/aromaticity/ring perception;
- [x] narrow structural `validate_stereo()` so represented stereo validity can be checked without `PerceptionState`; move semantic availability checks to appropriate perception/CIP helpers;
- [x] preserve direct SMILES `@`/`@@` stereo elements where already canonical.

Do not redesign coordinate-stereo publication yet.

Acceptance:

- normalized represented molecules contain canonical source-declared `StereoElement`s and no unresolved consumed source marks;
- default chemical perception need not assemble source stereo.

### Stage 5 — Remove perception-to-representation feedback

Goal: enforce the one-way boundary from represented molecule to perception state.

Implement:

- [x] delete `normalize_aromatic_nitrogen_hydrogens()` and adapt tests/writers/algorithms to read represented explicit H plus perceived implicit H correctly;
- [x] audit valence/rings/aromaticity operations for primary-representation mutation and remove any remaining representation writes;
- [x] add a strong regression that snapshots primary represented state before default perception and proves it is unchanged afterward except for `PerceptionState`.

Acceptance:

- general perception mutates only derived perception state.

### Stage 6 — Introduce the default perception pipeline and retire sanitization

Goal: replace the broad sanitizer concept with explicit normalization and default perception.

Implement:

- [x] add one transactional default discrete `perception::perceive()` pipeline using the agreed dependency order;
- [x] define focused `PerceptionError` wrapping valence/ring/aromaticity failures as needed;
- [x] do not recreate `SanitizeOptions` as a boolean-filled `PerceptionOptions`;
- [x] experts use individual focused algorithms for custom subsets/models;
- [x] remove `SanitizeOptions`, `SanitizeReport`, `SanitizeError`, sanitizer functions, sanitizer facade exports, and `SmallMolecule::sanitize*`;
- [x] remove/replace `from_smiles_sanitized`; do not finalize the ordinary fast-forward constructor yet unless the new low-level API has already made its semantics obvious;
- [x] retire the sanitization feature/docs terminology and update examples/benchmarks/tests.

Acceptance:

- [x] the central public vocabulary is normalize + perceive, not sanitize;
- [x] default perception is transactional and leaves representation unchanged on success and failure.

### Stage 7 — Clean parse / interpret boundaries

Goal: make format `*Document` values genuinely format-specific and interpretation the sole owner of source-to-chemistry mapping.

Implement incrementally by format, starting with SMILES and then Molfile/SDF:

- parsers produce typed format syntax/AST state rather than core `Atom`/`BondOrder` chemistry objects where practical;
- interpretation converts source syntax to core represented chemistry;
- remove general chemical inference from interpretation;
- replace SMILES bracket-radical inference with exact format semantics if the radical is semantically determined by SMILES; otherwise leave it to perception rather than guessing;
- remove Molfile `preserve_molfile_tetrahedral_hydrogens()` valence inference and encode only actual source hydrogen/valence semantics;
- keep source-to-canonical mappings and interpretation reports.

Do not force all formats through a generic `Document` trait.

Acceptance:

- parse can represent syntactically valid chemistry that interpretation or normalization may later reject;
- interpretation returns represented chemistry with empty general `PerceptionState`.

### Stage 8 — Coordinate stereo and ordinary fast-forward API

Goal: finish the user-facing workflow after the low-level contracts are clean.

Implement:

- make coordinate-derived stereo inference a read-only perception/result operation;
- provide an explicit transform to materialize inferred coordinate stereo into represented `StereoElement` state if useful;
- decide and implement the ordinary fast-forward convenience API that composes parse/interpret/normalize/default-perceive for common small-molecule use;
- update README examples to use the ordinary path while retaining explicit expert examples;
- consider whether CIP belongs in the ordinary convenience or remains opt-in; default preference is opt-in unless there is a strong usability reason otherwise.

Acceptance:

- ordinary users have one obvious path to a ready-to-use molecule;
- expert users retain direct access to every semantic stage;
- convenience composition does not weaken the underlying stage boundaries.

### Stage 9 — Final desloppification and terminology cleanup

Goal: remove migration residue.

Audit and remove:

- sanitize naming and compatibility shims;
- unused perception helpers that existed only for source representation handling;
- `AromaticityProvenance::Imported` residue;
- old source-stereo perception report types or fields no longer meaningful;
- redundant test-only compatibility paths where feasible;
- `Hydrogen Topology Normalization` terminology if it conflicts with the canonical normalization concept;
- outdated feature docs, dashboard entries, README examples, benchmark plumbing, and comments.

Run full repository formatting, clippy, tests, doc tests, feature/dashboard consistency, and relevant benchmark smoke checks.

## Refactor discipline

Throughout the refactor:

- prefer focused functions and explicit domain types over generic pipeline abstractions;
- do not add state freshness generations, caches, registries, or universal report frameworks unless a concrete need emerges;
- preserve transactional publication for operations that can fail;
- never silently choose a tautomer, protonation state, fragment parent, or other chemically distinct state under normalization;
- do not preserve obsolete public APIs merely for compatibility unless a specific external consumer requires it;
- keep alternative future perception states possible without building a generic multi-state framework now;
- keep representation transforms and derived perception visibly separate in code ownership and public naming.

## Completion criterion

The refactor is complete when an ordinary normalized molecule has one stable represented graph independent of perception model, default discrete perception can be rerun or replaced without changing that representation, source-format interpretation no longer installs semantic chemical perception, and the old sanitization abstraction has disappeared from the central API and code organization.
