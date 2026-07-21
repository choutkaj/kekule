# DSSP Protein Secondary Structure

## Summary

Assign protein secondary structure from one explicit three-dimensional `Model`
coordinate snapshot using the DSSP hydrogen-bond and geometric definitions. The
feature is structural analysis, not sequence-based secondary-structure
prediction.

The implementation targets DSSP 4.6.1 behavior, including polyproline-II
(kappa) helices. Analysis is read-only and returns a derived result; it never
installs secondary-structure labels into `Molecule`, `MacroMolecule`,
`SmcraHierarchy`, or `Model`.

## Behavior/API

- The public namespace is `molecular::dssp`.
- The primary entry point is
  `dssp::assign(&Model, DsspOptions) -> Result<DsspResult, DsspError>`.
- Input is one already-constructed `Model`, whose positions are the sole
  authoritative coordinate set. DSSP does not parse files, choose an mmCIF
  model or alternate location, infer molecular boundaries, sanitize chemistry,
  or repair hierarchy data.
- All macro-molecule instances are considered together so hydrogen bonds and
  beta-sheet topology may cross chain and molecule-instance boundaries. Small
  molecules and non-peptide residues are ignored and counted in the report.
- `DsspResidueKey` qualifies each local `SmcraResidueId` with its
  `MoleculeInstanceId`. Results never flatten or renumber model topology.
- `DsspSecondaryStructure` represents the complete DSSP 4 summary alphabet:
  loop (` `), alpha helix (`H`), isolated beta bridge (`B`), extended beta
  strand (`E`), 3-10 helix (`G`), pi helix (`I`), polyproline-II/kappa helix
  (`P`), hydrogen-bonded turn (`T`), and bend (`S`). Conversion to the canonical
  one-character DSSP code is explicit.
- `DsspResult` provides deterministic hierarchy-order iteration and lookup by
  `DsspResidueKey`. Each analyzed residue exposes its summary assignment,
  chain-break status, optional backbone/geometric values (`phi`, `psi`,
  `omega`, `alpha`, `kappa`, and `TCO`), up to two strongest donor and acceptor
  hydrogen-bond partners with energies, and beta bridge/ladder/sheet
  relationships when present.
- `DsspStatistics` reports analyzed residues and chain segments, hydrogen bonds,
  secondary-structure counts, and intra-chain versus inter-chain beta bridges.
- `DsspReport` records ignored instances, non-peptide residues, incomplete or
  skipped residues, detected gaps, reconstructed amide hydrogens, and consumed
  work. Missing backbone atoms create an explicit skipped-residue issue and a
  chain break rather than silently assigning loop.
- Ambiguous backbone atom mapping, non-finite required geometry, absence of any
  analyzable protein residue, and exceeded resource limits return structured
  `DsspError` values. A failure leaves the input unchanged.
- `DsspOptions::default()` follows the pinned DSSP 4 reference behavior,
  including polyproline-II assignment. `DsspLimits` bounds analyzed residues,
  candidate hydrogen-bond pairs, and generated bridge/ladder topology.
- A result is a snapshot of the coordinates used during assignment. Coordinate
  changes do not update it; callers must recompute explicitly.

## Implementation Notes

- The algorithm lives inside the `molecular` crate behind the focused
  `molecular::dssp` facade. Do not create a separate DSSP crate or a runtime
  binding to `libdssp`.
- Use `SmcraHierarchy` chain and residue order plus atom-site labels to identify
  peptide backbone atoms. Preserve author and label identifiers only as source
  metadata; result identity is the qualified canonical residue key.
- Treat a polymer-sequence residue with the required `N`, `CA`, `C`, and `O`
  backbone positions as analyzable. A label sequence ID is the canonical
  polymer-membership signal; covalently attached non-polymer components remain
  ignored even when their atom names resemble a peptide backbone. Determine
  peptide continuity with DSSP-compatible chain and
  coordinate checks because mmCIF interpretation does not infer template
  polymer bonds. Missing or discontinuous residues split local pattern
  recognition.
- Backbone amide hydrogens are reconstructed deterministically according to the
  DSSP convention; proline is not treated as an amide-hydrogen donor. DSSP's
  compatibility behavior of deriving a non-initial residue's hydrogen from the
  preceding complete residue table row is retained even across a reported
  chain break, while continuity still gates torsions and local patterns.
- The Kabsch-Sander electrostatic hydrogen-bond calculation, assignment
  precedence, beta bridge/ladder construction, and DSSP 4 polyproline-II rule
  are implemented directly. A deterministic 9 angstrom C-alpha spatial index
  bounds candidate enumeration instead of using an unbounded all-pairs pass.
- DSSP 4.6.1 stores coordinate points and evaluates its distance kernel as
  single-precision floats. Only that reference-compatibility boundary is
  reproduced so threshold and top-two tie behavior matches `mkdssp`; public
  coordinates, angles, energies, and result fields remain `f64`.
- Sheet and ladder identifiers must be deterministic under atom insertion-order
  changes while retaining hierarchy chain/residue order in the public result.
- This repository has no standalone units dependency. Consistent with the live
  `Model` and potential APIs, physical values are documented `f64` values:
  coordinates and distance limits are angstroms, angles are degrees, and
  hydrogen-bond energies are kcal/mol. Unit-bearing names prevent ambiguity.
- The scientific anchors are Kabsch and Sander's original definition
  ([DOI 10.1002/bip.360221211](https://doi.org/10.1002/bip.360221211)) and
  DSSP 4 ([DOI 10.1002/pro.70208](https://doi.org/10.1002/pro.70208)). The
  official [PDB-REDO DSSP implementation](https://github.com/PDB-REDO/dssp)
  is the validation reference, not a Rust runtime dependency.

## Validation

- Golden generation is pinned to Biopython 1.87 and `mkdssp version 4.6.1`.
  The Windows reference executable used for the committed evidence has SHA256
  `963f7e3bfc46818817639430485ad698faee3fd4d26a75d25d895af8925b3d1f`.
  `validation/reference/biopython/environment.yml` recreates the version-pinned
  reference environment, and each manifest and golden records the command,
  executable checksum, fixture checksum, and `runtime_dependency = false`.
- The reference runner first constructs the same highest-occupancy coordinate
  snapshot used by normal mmCIF interpretation. Biopython invokes legacy DSSP
  on that explicit snapshot; a second `mkdssp --output-format=mmcif` call on the
  same snapshot supplies DSSP4-only geometry and topology fields. Neither tool
  is a Rust runtime dependency.
- Exact implementation-golden comparison covers residue author and label
  identity, residue inclusion/order, chain breaks, the nine-state summary,
  helix-position flags, sheet/strand/ladder topology and orientation, and both
  retained donor/acceptor slots. Phi, psi, kappa, and alpha use a 0.15 degree
  tolerance; TCO uses 0.0015; one-decimal legacy hydrogen-bond energies use
  0.051 kcal/mol. Tolerances were fixed before broad validation.
- Required baseline evidence uses all 100 fixtures in the provenance-pinned
  `pdb-100` corpus. The nested `pdb-1000` corpus supplies deliberate broad
  macromolecular coverage across the same five structural categories. Both
  source datasets are local-only and must be built before validation.
- Focused regressions directly cover the complete nine-state alphabet,
  alpha/3-10/pi/polyproline-II construction, parallel and antiparallel bridge
  formulas, deterministic top-two retention, snapshot immutability, and every
  resource-limit class. The broad external corpora additionally exercise chain
  gaps, termini, proline donors, incomplete and non-standard residues,
  inter-chain interactions, turns/bends, and multi-sheet topology.
- Exact categorical fields must match the pinned reference. Floating-point
  tolerances may not be widened merely to make mismatches pass.

## Out Of Scope

- Sequence-based secondary-structure prediction.
- Nucleic-acid secondary structure.
- Parsing or writing legacy DSSP files, or writing DSSP annotations into mmCIF.
- Solvent-accessible surface-area calculation.
- Mutating or caching derived assignments in `SmcraHierarchy`, `Model`, atom
  properties, or residue properties.
- Automatic model/alternate-location selection, hierarchy repair, missing-atom
  construction, template connectivity, sanitization, or preparation.
- Trajectory-wide, ensemble, or consensus secondary-structure analysis.
- Runtime dependencies on `mkdssp`, `libdssp`, Biopython, or other reference
  implementations.

## Revision Notes

- v2: Implement the read-only DSSP 4.6.1 kernel and public API, pin Biopython
  1.87 plus DSSP 4.6.1 reference generation, require checked-in smoke evidence,
  and add explicit PDB-10/PDB-100 supplemental validation.
- v1: Establish the planned DSSP 4-compatible, read-only analysis contract.
- v3: Use PDB-100 as the required macromolecular baseline and retire smoke and PDB-10 as validation corpora.
- v4: Ignore backbone-like non-polymer components, compare validation records by source residue identity, and exclude explicit reference-tool failures from generated manifests.
