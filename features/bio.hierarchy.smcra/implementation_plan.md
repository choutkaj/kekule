# Implementation Plan

## Feature ID

`bio.hierarchy.smcra`

## Goal

Represent biomolecular model, chain, residue, and atom-site hierarchy as a sidecar over the shared core `Molecule` graph. The hierarchy should preserve labels needed by PDB/mmCIF-style workflows without putting biomolecular identifiers into core `Atom`.

SMCRA means structure, model, chain, residue, atom. In this crate, `MacroMolecule` is the structure-level owner, `Molecule` is the chemistry graph, and `BioHierarchy` stores biomolecular organization.

## Public API

Add or refine APIs equivalent to:

- `BioHierarchy`
- `ModelId`
- `ChainId`
- `ResidueId`
- `AtomSiteId` if atom-site records need identity beyond `AtomId`
- `Model`
- `Chain`
- `Residue`
- `AtomSite`
- `MacroMolecule { mol: Molecule, hierarchy: BioHierarchy }`

Expected operations:

- Add model.
- Add chain under model.
- Add residue under chain.
- Attach an existing `AtomId` to a residue with atom-site metadata.
- Iterate models, chains, residues, and atom sites in input order.
- Lookup hierarchy placement by `AtomId`.
- Remove or invalidate hierarchy entries when the referenced atom is deleted, if deletion integration exists.

## Internal Modules Touched

Expected scope:

- `crates/molecules/src/lib.rs` initially, or a `bio` module if splitting becomes warranted.
- `BioHierarchy` data structures.
- `MacroMolecule` wrapper.
- Unit tests for hierarchy insertion and lookup.
- Feature docs under `features/bio.hierarchy.smcra/`.

Do not add mmCIF parsing in this feature. The parser should consume this API later.

## Data Model

Hierarchy sidecar:

- `BioHierarchy`
  - ordered models
  - ordered chains per model
  - ordered residues per chain
  - atom-site entries referencing core `AtomId`
  - lookup map from `AtomId` to atom-site placement
  - `props: PropMap`

Chemically non-general labels belong here:

- model number or ID
- chain author ID
- chain label ID
- residue name
- residue sequence ID
- insertion code
- atom name
- alternate location ID
- occupancy
- B-factor
- mmCIF label identifiers
- mmCIF author identifiers

Core `Atom` should not grow these fields.

## Algorithm Outline

1. Define typed IDs for model, chain, residue, and optional atom-site records.
2. Store hierarchy records in stable slot or append-only vectors.
3. Enforce parent-child relationships when adding records.
4. When attaching an atom, validate that the `AtomId` exists in `MacroMolecule.mol`.
5. Record atom placement and atom-site metadata in `BioHierarchy`.
6. Maintain `AtomId -> placement` lookup for efficient access.
7. Decide and document duplicate placement policy. Initial recommendation: one primary hierarchy placement per `AtomId`; alternate locations can be separate atom sites if they correspond to separate core atoms.
8. Keep ordering deterministic and matching input order.

## Tests

Add tests for:

- Empty hierarchy.
- Add model, chain, residue, and atom site.
- Iteration order is insertion order.
- Lookup placement by `AtomId`.
- Reject attaching an invalid atom ID.
- Reject chain with missing model.
- Reject residue with missing chain.
- Reject atom site with missing residue.
- Preserve author and label IDs separately.
- Preserve residue insertion code.
- Preserve alternate location, occupancy, and B-factor on atom-site metadata.
- `MacroMolecule` uses one core `Molecule` plus `BioHierarchy`.
- Core `Atom` remains free of biomolecular labels.

## Reference Validation

Use Biopython through `validation.harness` for hierarchy behavior.

Golden JSON should include:

- Biopython version.
- Input fixture name.
- Models.
- Chains.
- Residues with author and label IDs where available.
- Atom-site labels and corresponding atom ordering.
- Alternate locations, occupancy, and B-factors.

Validation should compare hierarchy structure and labels, not chemistry perception.

## Required Checks

Run:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace
cargo xtask dashboard --check
cargo xtask validate --feature bio.hierarchy.smcra
```

## Risks

- Accidentally duplicating molecule graph topology in `BioHierarchy`.
- Leaking biomolecular labels into core `Atom`.
- Confusing author IDs and label IDs from mmCIF.
- Alternate-location handling can imply multiple coordinates or atom identities.
- Atom deletion can leave stale hierarchy references unless policy is explicit.

## Edge Cases

- Multiple models.
- Empty chains.
- Residues with insertion codes.
- Non-polymer residues and ligands.
- Water residues.
- Alternate locations.
- Missing occupancy or B-factor.
- Duplicate atom names within a residue due to altlocs.
- mmCIF label IDs and author IDs disagree.

## Explicitly Out of Scope

- mmCIF parsing.
- PDB parsing.
- Sequence alignment.
- Secondary structure assignment.
- Polymer chemistry perception.
- Coordinate/conformer algorithms.
- RDKit or Biopython runtime dependencies.
