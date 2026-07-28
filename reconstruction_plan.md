# Canonical reconstruction unblock plan

> **Status:** Active implementation plan for the MolStudio persistence blocker.
> This is intentionally a narrow API milestone, not a general serialization
> framework. Remove this file in the implementation PR once canonical feature
> contracts, rustdoc, and architecture documentation own the shipped behavior.

## Purpose

MolStudio persists canonical edited Molecular state rather than relying on source
reparsing. Its versioned DTO layer can currently serialize most Molecular
objects, but three valid Molecular-owned states cannot be reconstructed exactly
through public checked APIs:

1. installed `PerceptionState`;
2. enriched SMCRA child metadata and child property maps;
3. the complete stable stereo-group slot/tombstone layout.

MolStudio correctly rejects these states instead of silently re-running
perception, simplifying hierarchy data, or renumbering stereo groups. This
milestone supplies the smallest safe public APIs needed to remove those
rejections.

The primary feature ID should be:

```text
model.canonical-reconstruction
```

Update every directly affected existing feature contract, especially the core
molecule/perception, stereochemistry, biological hierarchy, and public-facade
contracts identified during implementation.

## Design principles

- Molecular owns validation and invariant-preserving reconstruction of Molecular
  state.
- MolStudio continues to own its versioned archive DTOs and Serde representation.
- Do not derive or require Serde on Molecular runtime types.
- Do not add a general `MoleculeReconstructionState`, project archive, wire
  format, or second chemistry model.
- Do not expose unrestricted mutable access to `PerceptionState` internals.
- Restore installed perception as one complete transaction, never through a
  sequence of public partially mutating installation calls.
- Add only targeted SMCRA mutation where the fields are metadata or properties
  and do not affect parentage or atom placement.
- Preserve stable stereo-group identifiers, including interior and trailing
  tombstones, without dummy stereo elements or groups.
- Loading must never sanitize, re-perceive, renumber, or silently discard state.

## Scope

The milestone consists of three focused API additions and consumer validation.

### 1. Public exact `PerceptionState` construction and installation

External callers must be able to:

- inspect every installed perception section without information loss;
- construct an equivalent detached `PerceptionState` through public APIs;
- validate and atomically install it into one existing `Molecule`;
- receive structured errors for invalid references or inconsistent state;
- leave the molecule completely unchanged on failure.

The exact ergonomic design may use immutable public component value types or a
`PerceptionStateBuilder`. It should remain close to the existing internal model
and expose semantics equivalent to:

```rust
impl PerceptionState {
    pub fn builder() -> PerceptionStateBuilder;

    pub fn valence_state(&self) -> Option<ValencePerceptionStateRef<'_>>;
    pub fn ring_state(&self) -> Option<RingPerceptionStateRef<'_>>;
    pub fn aromaticity_state(&self) -> Option<AromaticityPerceptionStateRef<'_>>;
    pub fn cip_descriptors(
        &self,
    ) -> impl ExactSizeIterator<Item = (StereoElementId, StereoDescriptor)> + '_;
}

impl PerceptionStateBuilder {
    pub fn with_valence(
        self,
        model: Option<ValenceModel>,
        implicit_hydrogens: impl IntoIterator<Item = (AtomId, u8)>,
    ) -> Result<Self, PerceptionStateBuildError>;

    pub fn with_rings(
        self,
        ring_atoms: impl IntoIterator<Item = AtomId>,
        ring_bonds: impl IntoIterator<Item = BondId>,
        ring_set: Option<InstalledRingSet>,
    ) -> Result<Self, PerceptionStateBuildError>;

    pub fn with_aromaticity(
        self,
        provenance: AromaticityProvenance,
        atoms: impl IntoIterator<Item = AtomId>,
        bonds: impl IntoIterator<Item = BondId>,
    ) -> Result<Self, PerceptionStateBuildError>;

    pub fn with_cip_descriptors(
        self,
        descriptors: impl IntoIterator<Item = (StereoElementId, StereoDescriptor)>,
    ) -> Result<Self, PerceptionStateBuildError>;

    pub fn build_for(
        self,
        molecule: &Molecule,
    ) -> Result<PerceptionState, PerceptionStateInstallError>;
}

impl Molecule {
    pub fn install_perception_state(
        &mut self,
        state: PerceptionState,
    ) -> Result<(), PerceptionStateInstallError>;
}
```

This is an illustrative contract, not mandatory naming. A small set of public
immutable component structs with constructors/getters is equally acceptable.
The final design must avoid exposing raw mutable maps or the existing
crate-private incremental `install_*` methods as public workflow APIs.

#### Perception information that must round-trip

The public surface must preserve the distinction between an absent section and a
present section with empty or partially identified content. In particular:

- no valence section;
- a valence section with `model: None` and installed implicit-H assignments;
- a valence section with a named `ValenceModel`;
- ring membership without an installed ring basis;
- ring membership plus the complete installed deterministic `RingSet`, including
  `RingWork`;
- no aromaticity section;
- aromaticity provenance `Imported` versus `Perceived(model)`;
- complete aromatic atom and bond membership;
- complete CIP descriptor assignments.

Existing convenience queries such as `has_valence()` must not be used as the
only serialization signal because they do not encode all section-presence
states.

#### Perception installation validation

Before replacing the current state, validate at least:

- every implicit-H assignment references a live atom;
- ring atom and bond membership references only live slots;
- ring membership dimensions or generated flags match the molecule's complete
  stable atom/bond slot layout;
- every installed ring references live atoms and bonds and is internally
  coherent with the installed membership;
- aromaticity atoms and bonds are live;
- every CIP descriptor references a live stereo element;
- duplicate input entries are rejected rather than last-write-wins;
- fixed-width and allocation bounds are checked before allocation;
- no stronger dependency is invented than Molecular currently requires (for
  example, imported aromaticity must not be rejected merely because a perceived
  ring basis is absent).

Installation is all-or-nothing. A failure must preserve the prior perception
state exactly. A successful installation must subsequently follow ordinary
invalidation rules: graph chemistry/connectivity mutation clears affected
perception, stereo mutation clears CIP descriptors, and coordinate/property-only
changes remain perception-neutral.

### 2. Targeted enriched SMCRA restoration

The normal hierarchy builder already reconstructs chain/residue/atom-site
parentage, IDs, ordering, and atom placement. Do not replace it with a parallel
hierarchy DTO inside Molecular.

Add only the missing safe metadata/property operations, approximately:

```rust
impl SmcraHierarchy {
    pub fn set_residue_component_ids(
        &mut self,
        residue: SmcraResidueId,
        label_comp_id: Option<String>,
        author_comp_id: Option<String>,
    ) -> Result<(), SmcraHierarchyError>;

    pub fn chain_props_mut(
        &mut self,
        chain: SmcraChainId,
    ) -> Result<&mut PropMap, SmcraHierarchyError>;

    pub fn residue_props_mut(
        &mut self,
        residue: SmcraResidueId,
    ) -> Result<&mut PropMap, SmcraHierarchyError>;

    pub fn atom_site_props_mut(
        &mut self,
        site: SmcraAtomSiteId,
    ) -> Result<&mut PropMap, SmcraHierarchyError>;
}
```

Equivalent whole-map replacement methods are acceptable. These operations must:

- validate the requested child ID before mutation;
- leave hierarchy parentage, child ordering, IDs, and atom lookup unchanged;
- preserve arbitrary supported `PropMap` values exactly;
- permit distinct `name`, `label_comp_id`, and `author_comp_id` values;
- remain usable through `MacroMoleculeBuilder` and `MacroMoleculeEditor` without
  bypassing final graph/hierarchy validation.

Do not add broad mutable access to complete `SmcraChain`, `SmcraResidue`, or
`SmcraAtomSite` records.

### 3. Stereo-group stable slot/tombstone reconstruction

Expose the full stereo-group slot layout and provide one capacity-checked append
operation for a deleted slot. The public semantics should be equivalent to:

```rust
impl Molecule {
    pub fn stereo_group_slot_count(&self) -> usize;

    pub fn stereo_group_slots(
        &self,
    ) -> impl ExactSizeIterator<Item = (StereoGroupId, Option<&StereoGroup>)>;

    pub fn append_stereo_group_tombstone(
        &mut self,
    ) -> Result<StereoGroupId, MoleculeError>;
}
```

A slot count plus the existing live-group iterator is acceptable instead of a
slot iterator, provided interior and trailing tombstones are unambiguous.

`append_stereo_group_tombstone` must:

- append exactly one `None` slot;
- preserve every existing group and stereo-element reference;
- return the stable ID of the appended deleted slot;
- check fixed-width identifier capacity before mutation;
- not create a dummy group or claim a live chemical object exists;
- not invalidate CIP descriptors solely because an empty deleted slot was
  appended.

MolStudio can then reconstruct serialized slots in order:

```text
Some(group) -> add_stereo_group(group)
None        -> append_stereo_group_tombstone()
```

The resulting group IDs and the next future group ID must match the source
molecule exactly.

#### Empty-group invariant

Current public insertion rejects empty stereo groups, but removing or pruning a
stereo element can leave a live group with no members. Make the invariant
consistent:

- a live stereo group must contain at least one live stereo element;
- when element removal/pruning empties a group, that group becomes a tombstone;
- remaining nonempty groups retain their stable IDs and membership;
- reconstruction rejects live empty groups.

This is a focused invariant correction required so every valid state produced by
current APIs remains reconstructible through the same public contract.

## Loading order for downstream DTO adapters

Document the safe order for MolStudio and other persistence consumers:

```text
1. reconstruct atoms, bonds, conformers, and their existing stable layouts
2. reconstruct stereo elements with no pre-installed group side effects
3. replay stereo-group live slots and tombstones in slot order
4. restore stereo bond marks
5. reconstruct SMCRA hierarchy and enrich child metadata/properties
6. construct and install the complete PerceptionState last
```

Perception is installed last because ordinary graph and stereo construction
correctly invalidates computed state.

## Non-goals

This milestone does not add:

- Serde derives or a stable binary/JSON format in Molecular;
- a MolStudio project archive implementation;
- a general `MoleculeReconstructionState` or `TopologyReconstructionState`;
- serialization of process-local `TopologyIdentity`;
- public incremental perception mutation;
- public mutable access to hierarchy parentage or child arrays;
- dummy chemistry for ID padding;
- re-perception, sanitization, or source reparsing during load;
- reconstruction of malformed historical states by silently coercing them.

## Error model

Use focused public `#[non_exhaustive]` errors. Exact names may follow repository
conventions, but callers must be able to distinguish at least:

```text
perception reconstruction
- invalid atom reference
- invalid bond reference
- invalid stereo-element reference
- duplicate assignment/reference
- ring membership slot mismatch
- malformed installed ring
- inconsistent installed ring membership
- identifier/allocation capacity failure

SMCRA enrichment
- invalid chain ID
- invalid residue ID
- invalid atom-site ID

stereo slot reconstruction
- stereo-group identifier capacity exceeded
- invalid or empty live stereo group through existing group construction
```

Do not return free-form strings where an existing or focused structured variant
can represent the failure.

## Required tests

### Perception round-trip

- Export and publicly reconstruct an empty `PerceptionState`.
- Round-trip installed valence with a named model and atom-wise implicit H.
- Round-trip a present valence section with `model: None`.
- Round-trip ring membership without a ring basis.
- Round-trip ring membership plus `RingSet` and `RingWork`.
- Round-trip imported aromaticity and perceived aromaticity distinctly.
- Round-trip CIP assignments.
- Round-trip one molecule containing every installed section and assert exact
  `PerceptionState` equality.
- Reject references to atom, bond, and stereo-element tombstones.
- Reject malformed ring data and inconsistent membership.
- Prove failed installation preserves the previous state exactly.
- Prove graph mutation after installation invalidates topology-derived state and
  stereo mutation invalidates CIP as before.

### SMCRA enrichment

- Reconstruct chains, residues, and atom sites with nonempty property maps at
  every hierarchy level.
- Preserve distinct residue `name`, `label_comp_id`, and `author_comp_id`.
- Preserve all existing atom-site metadata plus atom-site properties.
- Reject invalid child IDs without changing any hierarchy state.
- Build a `MacroMolecule` from the reconstructed hierarchy and assert complete
  equality with the source hierarchy.

### Stereo groups

- Reconstruct a slot sequence containing live, interior tombstone, live, and
  trailing tombstone entries.
- Assert live `StereoGroupId` values remain identical.
- Append the same next group to source and reconstructed molecules and assert the
  returned ID is identical.
- Remove the last member of a group and assert the group slot becomes a
  tombstone.
- Remove one member from a multi-member group and assert the group remains live.
- Reject reconstruction of a live empty group.
- Verify appending a tombstone alone does not clear installed CIP state.

### Topology and consumer validation

- Build topologies from original and reconstructed small/macro definitions and
  assert `same_layout`.
- Through the ignored sibling Cargo patch, replace MolStudio's structured
  rejection for these three states with the new public APIs.
- Add MolStudio project round-trips containing installed perception, enriched
  SMCRA child state, and an interior/trailing stereo-group tombstone layout.
- Verify no source reparse or re-perception occurs during load.
- Keep MolStudio changes separate from the Molecular PR unless explicitly
  requested; do not commit path dependencies or path-source lockfile changes.

## Staged implementation

### Stage 0: Audit and contract lock

- Read `AGENTS.md`, `ARCHITECTURE.md`, this plan, and all affected feature
  contracts.
- Inspect the MolStudio DTO rejection and record the exact three unsupported
  cases.
- Create `features/model.canonical-reconstruction/feature.toml` and `feature.md`
  with honest initial status.
- Confirm no additional state is blocked before broadening scope.

Exit criterion: the public semantics and affected feature IDs are explicit.

### Stage 1: Perception state

- Add lossless public read access for every installed section.
- Add minimal public construction API.
- Add checked atomic installation on `Molecule`.
- Add round-trip, malformed-state, and invalidation tests.

Exit criterion: MolStudio can encode and reconstruct all current installed
perception states without invoking perception algorithms.

### Stage 2: SMCRA enrichment

- Add residue component-ID restoration.
- Add checked child property-map mutation/replacement.
- Add hierarchy and `MacroMolecule` round-trip tests.

Exit criterion: all current enriched SMCRA child fields and properties round-trip
through public APIs.

### Stage 3: Stereo-group slots

- Expose slot count/layout.
- Add capacity-checked tombstone append.
- Enforce the nonempty-live-group invariant after member pruning.
- Add stable-ID and next-ID regressions.

Exit criterion: every valid stereo-group slot layout produced by Molecular can be
reconstructed exactly.

### Stage 4: Canonical documentation and consumer proof

- Update rustdoc, `ARCHITECTURE.md`, `CHANGELOG.md`, the new feature contract, and
  every directly affected feature contract.
- Regenerate dashboard artifacts through repository tooling; never hand-edit
  generated HTML.
- Validate the sibling MolStudio DTO path using the untracked local patch and
  remove the three structured capability errors locally.
- Delete `reconstruction_plan.md` once canonical docs own the contract.

Exit criterion: Molecular is available at a public Git revision that lets
MolStudio round-trip the three formerly blocked canonical states losslessly.

## Validation

Run at minimum:

```bash
cargo fmt --all -- --check
cargo check --workspace --all-targets --all-features --locked
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo test --workspace --all-features --locked
cargo test --workspace --all-features --doc --locked
RUSTDOCFLAGS="-D warnings" cargo doc --workspace --all-features --no-deps --locked
cargo xtask dashboard --check
cargo xtask skills --check
```

Also run focused core, perception, stereo, bio/hierarchy, topology-layout, and
MolStudio project-persistence tests. Report every omitted check and reason.

## Completion criteria

The milestone is complete only when:

- installed perception is publicly inspectable, constructible, checked, and
  atomically installable;
- present-versus-absent perception sections retain their exact meaning;
- enriched SMCRA child component IDs and properties round-trip;
- interior and trailing stereo-group tombstones preserve stable and next IDs;
- live empty stereo groups can no longer be produced by normal deletion/pruning;
- malformed state fails structurally without partial mutation;
- original and reconstructed topology definitions satisfy `same_layout`;
- MolStudio removes its three capability rejections under local consumer
  validation;
- all locked Molecular gates pass;
- canonical docs replace this temporary plan.
