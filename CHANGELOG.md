# Changelog

All notable changes to this project will be documented in this file.

The project follows Cargo semantic-versioning conventions. During the `0.x`
series, breaking public API changes increment the minor version.

## Unreleased

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
- Move DSSP and potential evaluation to borrowed structural views. Bind
  harmonic and DREIDING preparation to exact topology identity; DREIDING now
  exposes explicit whole-topology, molecule-instance, and connected-component
  QEq grouping.
- Remove the obsolete public model-owned topology types and move general
  geometry to `molecular::geometry`.
- Reclassify RDKit, Biopython, and DSSP external parity from a required
  repository gate to optional benchmarking; preserve corpora and goldens, add
  neutral result observations, and remove repository identity from benchmark
  input digests and corpus selection identifiers.
- Hard-rename the repository and publishable crates from `molecules` to
  `molecular`, including the DREIDING adapter, Rust import paths, benchmark
  tooling, generated writer provenance, and project branding.

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
