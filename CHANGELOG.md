# Changelog

All notable changes to this project will be documented in this file.

The project follows Cargo semantic-versioning conventions. During the `0.x`
series, breaking public API changes increment the minor version.

## Unreleased

- Hard-rename the project, foundational crate, DREIDING adapter, trajectory I/O
  companion, Rust import roots, repository URLs, generated Molfile provenance,
  and visual identity to Kekule for the breaking 0.3 release. Benchmark input
  digests now label core sources through a package-neutral namespace.
- Add focused canonical molecular reconstruction for external persistence:
  exact detached `PerceptionState` inspection/build/install with transactional
  graph, stereo, ring, slot, duplicate, and capacity validation; targeted
  enriched SMCRA component-ID and child-property restoration; and exact
  stereo-group tombstone replay with stable next IDs, CIP-neutral deleted-slot
  append, automatic final-member group tombstoning, transactional rejection of
  pre-grouped stereo-element insertion, and detached grouped-element removal.
  Add public
  original/reconstructed `same_layout` regressions and validate MolStudio
  project round-trips through the ignored sibling path patch without Serde on
  Kekule runtime objects, source reparsing, or reperception.
- Add the core contracts required by external trajectory codecs: complete
  transactional allocation-reusing `FrameBuffer` publication, explicit
  exact-topology atom-order assertions, non-exhaustive trajectory format
  identity, typed file/codec error context, and downstream trait-implementation
  tests. Add the supported one-way `kekule-trajectory-io` companion with
  bounded signature-plus-extension detection, metadata/reports/limits,
  one-handle sequential and verified indexed access, strict atomic path
  writers that cannot publish after a failed frame, and strict multi-frame XYZ
  read/write with explicit units,
  element-order validation, and a provenance-pinned ASE fixture. Add
  supported common-profile DCD read/write with both byte orders, cells,
  fixed-atom reconstruction, checked records and counts, indexed access,
  explicit step/time policy, and a provenance-pinned MDAnalysis fixture. Add
  supported pure-Rust TRR/XDR f32/f64 read/write with triclinic cells,
  optional velocities and forces, time, step, explicit lambda preservation,
  mixed-precision indexing, and reciprocal MDAnalysis interoperability. Add a
  supported checked XTC decoder for the 1995/2023 magic profiles, small and
  compressed coordinates, explicit lossy precision and box policies, bounded
  full-decode indexing, and a private panic-contained writer adapter over
  audited pure-Rust `molly` 0.6.1 with reciprocal MDAnalysis interoperability.
  Harden the supported profiles with pre-decode frame/index-limit probes,
  restoration-before-publication random reads, strict DCD `NSET` consistency,
  cumulative sequential TRR precision metadata, and signed-i32 XTC counts and
  steps. Bound index allocation with capped geometric growth, require every
  production writer to finish at least one frame, and remove ordinary
  per-frame DCD EOF read/seek probes.
- Add focused `kekule::alignment` weighted Kabsch fitting over exact-topology
  `ModelView` and `AtomSelection` inputs. Results map moving coordinates into
  reference coordinates with a proper `RigidTransform`, unit-bearing post-fit
  RMSD, structured rank/weight/topology/periodic failures, and no coordinate
  mutation or periodic imaging.
- Add immutable whole-instance topology subset transformations with complete
  checked lineage plus explicit transactional remapping for topology-bound
  positions, models, selections, ensembles, trajectory frames, trajectories,
  and reusable frame buffers.

## 0.2.0 - 2026-07-25

- Complete the topology-centered hard break: add immutable coordinate-free
  `Topology` with reusable molecule definitions, explicit instances,
  authoritative dense order, exact identity, explicit `same_layout` comparison,
  selections, and checked lineage mappings. Identity maps use
  `between_identical_layouts`; general order-independent structural equivalence
  and isomorphism mapping remain planned.
- Define `Model = Topology + Configuration`, add topology-bound `Positions`,
  validated periodic cells, borrowed `ModelView`, and coordinate-model-specific
  `StructureObservation`.
- Add finite shared-topology `Ensemble` members and fixed-topology trajectory
  frames, reusable `FrameBuffer`, in-memory storage, sequential/seekable reader
  contracts, writer contracts, and an atom-order-asserted coordinate source.
- Make SMCRA hierarchy coordinate-independent and move source model ID,
  alternate location, occupancy, B-factor, source atom-site ID, and raw
  Cartesian values into observation state.
- Keep ordinary mmCIF interpretation explicitly single-model; add a separate
  consistency-proving multi-model ensemble path. Distance heuristics now report
  connectivity candidates without asserting unsupported single bonds.
- Strengthen ensemble source-atom correspondence with sequence, insertion,
  asymmetry, component, occurrence, atom, and selected-altloc provenance rather
  than derived molecule insertion order.
- Keep topology construction coordinate-independent by validating only static
  macromolecular graph/hierarchy state; model construction validates only its
  selected source conformer, while explicit standalone macro validation retains
  full-conformer checking.
- Make every collection-backed public ID insertion checked and transactional:
  `Molecule::add_atom` is now fallible, graph/conformer/stereo/SMCRA/topology
  capacity failures are structured by identifier kind, iterators avoid
  narrowing slot indices, and writer-generated one-based serials use `u64`.
- Resolve declared mmCIF `_struct_conn` partners from all supplied label,
  author, insertion, and alternate-location identity fields. Ambiguous,
  conflicting, selected-out, and missing partners are source-aware report
  issues and never bind by source-row order or spatial proximity.
- Move DSSP and potential evaluation to borrowed structural views. Bind
  harmonic and DREIDING preparation to exact topology identity; DREIDING now
  exposes explicit whole-topology, molecule-instance, and connected-component
  QEq grouping.
- Remove the obsolete public model-owned topology types and move general
  geometry to `kekule::geometry`.
- Reclassify RDKit, Biopython, and DSSP external parity from a required
  repository gate to optional benchmarking; preserve corpora and goldens, add
  neutral result observations, and remove repository identity from benchmark
  input digests and corpus selection identifiers.
- Consolidate the then-current repository and publishable-crate branding,
  including the DREIDING adapter, Rust import paths, benchmark tooling,
  generated writer provenance, and project assets.

## 0.1.0 - 2026-07-16

Initial release.

- Stable-ID molecular graph kernel with explicit perception state.
- Separate small-molecule and macromolecule domain boundaries.
- Staged SMILES, Molfile, SDF, and mmCIF parsing and interpretation.
- Configurable parser resource limits and structured rejection of malformed or
  unsafe record boundaries.
- Qualified biomolecular hierarchy and mmCIF provenance.
- Fixed-topology molecular models and the DREIDING adapter.
- Explicit sanitization, hydrogen normalization, query, substructure,
  canonicalization, and modelling workflows.
- Bounded parser fuzz smoke tests in CI and longer scheduled campaigns.
- Generated feature/corpus benchmark matrix with PubChem-1k and PDB-100
  defaults, broader PubChem/Enamine/PL-REX/PDB-1000 tiers, and reproducible
  RDKit/Biopython/mkdssp reference artifacts.
