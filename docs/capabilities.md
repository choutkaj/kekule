# Current capabilities

This is the entry point for current support boundaries. Detailed contracts live
beside their implementations; dated audits and benchmark reports record results
at their stated revisions and are not a current list of missing features.

| Area | Current support and boundaries | Authoritative contract |
| --- | --- | --- |
| Chemistry ownership | Connected molecules; systems, coordinates, and derived perception have separate owners | [Architecture](../ARCHITECTURE.md) |
| Conjugation and resonance | RDKit 2026.03.3 conjugation in default perception; explicit connected-group preparation and bounded contributor enumeration with all five RDKit flags; source chemistry remains unchanged | [Guide](resonance.md), [conjugation](../crates/kekule/src/algorithms/conjugation.rs), [resonance](../crates/kekule/src/algorithms/resonance.rs) |
| Measurements | Cartesian distances, angles, and signed dihedrals; optional consecutive-bond validation; deterministic bond dihedrals use CIP priority and atom-ID ties with fixed references and optional values; periodic preprocessing is explicit | [Measurement contract](../crates/kekule/src/structure/measure.rs), [bond dihedrals](../crates/kekule/src/structure/measure/bond_dihedrals.rs) |
| Geometry editing | Atomic Model distance, angle, and dihedral setters; reusable topology-bound edits with automatic fragments or explicit moving selections; selection translation, axis rotation, and rigid transforms; stored Cartesian coordinates only, without relaxation or ring-closure solving | [Geometry editing contract](../crates/kekule/src/structure/geometry_edit.rs) |
| Selections | Snapshot-bound atom and bond sets; checked picking, set algebra, hierarchy and graph expansion, predicates, query matches, and Cartesian proximity | [Atom selection contract](../crates/kekule/src/topology/selection.rs), [bond selection contract](../crates/kekule/src/topology/selection/bonds.rs), [spatial selection](../crates/kekule/src/structure/measure.rs) |
| Hydrogens | Explicit graph atoms and implicit counts; fixed counts survive reperception; conversions preserve information and stereo | [Hydrogen semantics](hydrogen-semantics.md) |
| SMILES | Parsing, interpretation, plain/isomeric/canonical writing, and supported CX radical/group input and tetrahedral group output; unsupported content fails explicitly | [Public format API](../crates/kekule/src/lib.rs), [stereo support](stereo-support.md) |
| SMARTS | Recursive queries, Boolean atom/bond logic, connectivity/valence/ring predicates, tetrahedral and alkene stereo, enhanced tetrahedral groups, tagged and topology matching | [Dialect and explicit exclusions](../crates/kekule/src/query/dialect.md) |
| Stereo/CIP | Represented tetrahedral, double-bond and atropisomer stereo; bounded CIP assignment; explicit cleanup and candidate detection remain separate operations | [Stereo support](stereo-support.md), [CIP options](../crates/kekule/src/algorithms/cip/mod.rs) |
| MOL/SDF | Independent records, V2000/V3000 interpretation and writing; automatic version selection and explicit unsupported-content errors | [MOL/SDF writer](../crates/kekule/src/io/molfile_write.rs) |
| mmCIF | Independent blocks, coherent residue altloc selection and overrides, compatible coordinate/conformer ensembles, hierarchy/source reports, and model/ensemble output | [Alternate locations](mmcif-alternate-locations.md), [interpretation](../crates/kekule/src/io/mmcif_interpret/mod.rs), [writer](../crates/kekule/src/io/mmcif_write.rs) |
| Trajectories | Fixed topology, XYZ/DCD/TRR/XTC, streaming, periodic transformations and analyses | [Trajectory I/O](../crates/kekule-traj/src/io/mod.rs), [analysis](../crates/kekule-traj/src/analysis.rs) |
| DSSP | Read-only analysis of selected model coordinates; conformer, residue eligibility and chain policies matter to reference comparisons | [Analysis contract](../crates/kekule/src/dssp/mod.rs) |
| Solvation | Periodic box sizing, TIP3P water packing, monovalent counterions and salt; no force-field assignment or equilibration | [Solvation contract](../crates/kekule/src/structure/solvation/mod.rs), [box sizing](../crates/kekule/src/structure/solvation/boxes.rs) |
| Potentials | Topology-bound energy/gradient contract; OpenFF evaluation of a `ParameterizedTopology` in vacuum without cutoffs (periodic cells rejected, constraints not applied); L-BFGS minimization of any model, ensemble member, or frame | [Potential contract and singular-geometry policy](../crates/kekule-potentials/src/lib.rs), [OpenFF capabilities](../crates/kekule-potentials/src/openff.rs), [minimizer](../crates/kekule-potentials/src/minimize.rs) |
| OpenFF | Bundled Rosemary and Ash, or compatible composed OFFXML rule sets and exported NAGL bundles; typed topology-bound parameters, explicit intramolecular distance constraints, complete library-charge assignment without a model, exact-molecule NAGL lookup and native NAGL CPU inference with declared model identity checks, no molecule size limit, and typed error kinds; explicit H required | [OpenFF contract](../crates/kekule-openff/CONTRACT.md) |

The DSSP investigation and force-field scalability/periodicity extensions remain
deferred. Raw reference disagreements are not automatically defects: comparisons
must account for input selection and documented scientific conventions without
discarding asserted fields or relabeling errors as agreements.

For maintenance validation, see [CI](../.github/workflows/ci.yml) and
[fuzzing](../FUZZING.md). External scientific comparisons are optional and never
release gates.

The retired OpenFF [validation report](https://github.com/choutkaj/kekule/blob/ead90ca94f66e6b4f3082e5db49d6ddc042cc32c/benchmarks/openff/VALIDATION.md) recorded
parity plots, numerical differences and CPU timings for the 110-molecule panel and
two-model NAGL checks at its revision. The crate contract defines current support.
