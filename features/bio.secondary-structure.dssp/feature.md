# DSSP Protein Secondary Structure

## Summary

Assign protein secondary structure from one explicit three-dimensional `Model`
coordinate snapshot using the DSSP hydrogen-bond and geometric definitions. The
feature is structural analysis, not sequence-based secondary-structure
prediction.

The target is DSSP 4-compatible behavior, including polyproline-II (kappa)
helices. Analysis is read-only and returns a derived result; it never installs
secondary-structure labels into `Molecule`, `MacroMolecule`, `SmcraHierarchy`,
or `Model`.

## Behavior/API

- The planned public namespace is `molecules::dssp`.
- The primary entry point is conceptually
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
- `DsspStatistics` reports analyzed residues and chains, hydrogen bonds,
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

- Implement the algorithm inside the `molecules` crate behind the focused
  `molecules::dssp` facade. Do not create a separate DSSP crate or a runtime
  binding to `libdssp`.
- Use `SmcraHierarchy` chain and residue order plus atom-site labels to identify
  peptide backbone atoms. Preserve author and label identifiers only as source
  metadata; result identity is the qualified canonical residue key.
- Treat a residue with the required `N`, `CA`, `C`, and `O` backbone positions
  as analyzable. Determine peptide continuity with DSSP-compatible chain and
  coordinate checks because mmCIF interpretation does not infer template
  polymer bonds. Missing or discontinuous residues split local pattern
  recognition.
- Reconstruct backbone amide-hydrogen positions deterministically according to
  the DSSP convention; proline is not treated as an amide hydrogen donor.
- Reproduce the Kabsch-Sander electrostatic hydrogen-bond calculation and DSSP
  assignment precedence, then add DSSP 4 polyproline-II recognition. Use
  spatial pruning for candidate pairs rather than an unbounded all-pairs pass.
- Sheet and ladder identifiers must be deterministic under atom insertion-order
  changes while retaining hierarchy chain/residue order in the public result.
- Physical values must use the canonical quantity types supplied to
  `molecules` by the standalone units dependency. This feature must not define
  a second units system or expose undocumented raw-unit `f64` values.
- The scientific anchors are Kabsch and Sander's original definition
  ([DOI 10.1002/bip.360221211](https://doi.org/10.1002/bip.360221211)) and
  DSSP 4 ([DOI 10.1002/pro.70208](https://doi.org/10.1002/pro.70208)). The
  official [PDB-REDO DSSP implementation](https://github.com/PDB-REDO/dssp)
  is the validation reference, not a Rust runtime dependency.

## Validation

- Before setting `implemented = true`, pin an exact DSSP 4 release, executable
  checksum, command line, and reference environment for golden generation.
- Externally supplied, provenance-pinned PDB/mmCIF fixtures must compare
  residue inclusion, chain breaks, nine-state summary assignments, the two
  strongest donor/acceptor hydrogen bonds and their energies, beta partners,
  ladders/sheets, and available backbone/geometric values.
- Start with targeted `pdb-10` coverage and require broader `pdb-100` comparison
  before claiming general compatibility. Add those corpora to
  `validation_required` only when their DSSP manifests and current evidence
  exist.
- Focused unit regressions should cover alpha, 3-10, pi, and polyproline-II
  helices; parallel and antiparallel beta topology; turns and bends; chain
  gaps; termini; proline donors; incomplete and non-standard residues;
  inter-chain interactions; degenerate geometry; determinism; and every
  resource limit.
- Exact categorical fields must match the pinned reference. Floating-point
  energies and angles use documented tolerances chosen before golden
  generation; tolerances may not be widened merely to make mismatches pass.
- The feature remains `implemented = false`, `validated = false`, and has no
  required validation corpus while this is only a contract.

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

- v1: Establish the planned DSSP 4-compatible, read-only analysis contract.
