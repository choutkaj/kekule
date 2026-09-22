# Current capabilities

This is the entry point for current support boundaries. Detailed contracts live
beside their implementations; dated audits and benchmark reports record results
at their stated revisions and are not a current list of missing features.

| Area | Current support and boundaries | Authoritative contract |
| --- | --- | --- |
| Chemistry ownership | Connected molecules; systems, coordinates, and derived perception have separate owners | [Architecture](../ARCHITECTURE.md) |
| Hydrogens | Explicit graph atoms and implicit counts; specified counts survive reperception; conversions preserve information and stereo | [Hydrogen semantics](hydrogen-semantics.md) |
| SMILES | Parsing, interpretation, plain/isomeric/canonical writing, and supported CX radical/group input and tetrahedral group output; unsupported content fails explicitly | [Public format API](../crates/kekule/src/lib.rs), [stereo support](stereo-support.md) |
| SMARTS | Recursive queries, Boolean atom/bond logic, connectivity/valence/ring predicates, tetrahedral and alkene stereo, enhanced tetrahedral groups, tagged and topology matching | [Dialect and explicit exclusions](../crates/kekule/src/query/dialect.md) |
| Stereo/CIP | Represented tetrahedral, double-bond and atropisomer stereo; bounded CIP assignment; explicit cleanup and candidate detection remain separate operations | [Stereo support](stereo-support.md), [CIP options](../crates/kekule/src/algorithms/cip/mod.rs) |
| MOL/SDF | Independent records, V2000/V3000 interpretation and writing; automatic version selection and explicit unsupported-content errors | [MOL/SDF writer](../crates/kekule/src/io/molfile_write.rs) |
| mmCIF | Independent blocks, compatible coordinate ensembles, hierarchy/source reports, and model/ensemble output | [Interpretation](../crates/kekule/src/io/mmcif_interpret/mod.rs), [writer](../crates/kekule/src/io/mmcif_write.rs) |
| Trajectories | Fixed topology, XYZ/DCD/TRR/XTC, streaming, periodic transformations and analyses | [Trajectory I/O](../crates/kekule-traj/src/io/mod.rs), [analysis](../crates/kekule-traj/src/analysis.rs) |
| DSSP | Read-only analysis of selected model coordinates; conformer, residue eligibility and chain policies matter to reference comparisons | [Analysis contract](../crates/kekule/src/dssp/mod.rs) |
| Potentials | Explicit DREIDING preparation and evaluation on the same topology; current adapter is nonperiodic | [DREIDING contract](../crates/kekule-potentials/src/dreiding/mod.rs) |

The DSSP investigation and force-field scalability/periodicity extensions remain
deferred. Raw reference disagreements are not automatically defects: comparisons
must account for input selection and documented scientific conventions without
discarding asserted fields or relabeling errors as agreements.

For maintenance validation, see [CI](../.github/workflows/ci.yml) and
[fuzzing](../FUZZING.md). External scientific comparisons remain optional; their
execution and provenance rules are in the [benchmark guide](../benchmarks/GUIDE.md).
The manually dispatched [reference-adapter tests](../.github/workflows/reference-tests.yml)
exercise pinned RDKit and Biopython integrations separately from ordinary CI.
